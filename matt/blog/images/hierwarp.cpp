
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// sampling/hierwarp.cpp*
#include <pbrt/sampling/hierwarp.h>
#include <pbrt/sampling/sampling.h>
#include <pbrt/util/math.h>
#include <pbrt/util/print.h>

#include <algorithm>

namespace pbrt {

Hierarchical2DWarp::Hierarchical2DWarp(absl::Span<const Float> values, int nx, int ny,
                                       const Bounds2f &domain)
    : domain(domain) {
    levels.push_back(Array2D<Float>(nx, ny));
    for (int y = 0; y < ny; ++y)
        for (int x = 0; x < nx; ++x)
            levels.back()(x, y) = values[y * nx + x];

    nx = RoundUpPow2(nx);
    ny = RoundUpPow2(ny);
    while (nx > 2 || ny > 2) {
        nx = std::max(1, nx / 2);
        ny = std::max(1, ny / 2);
        size_t prev = levels.size() - 1;
        levels.push_back(Array2D<Float>(nx, ny));
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                levels.back()(x, y) += (Lookup(prev, 2 * x,     2 * y) +
                                        Lookup(prev, 2 * x + 1, 2 * y) +
                                        Lookup(prev, 2 * x,     2 * y + 1) +
                                        Lookup(prev, 2 * x + 1, 2 * y + 1));
    }
    std::reverse(levels.begin(), levels.end());
}

Point2i Hierarchical2DWarp::SampleDiscrete(Point2f u, Float *pdfOut) const {
    Float pdf = 1;
    Point2i p(0, 0);

    for (size_t i = 0; i < levels.size(); ++i) {
        if (i > 0 && Resolution(i).x > Resolution(i-1).x) p.x *= 2;
        if (i > 0 && Resolution(i).y > Resolution(i-1).y) p.y *= 2;

        Float wx[2] = { Lookup(i, p.x, p.y) + Lookup(i, p.x, p.y + 1),
                        Lookup(i, p.x + 1, p.y) + Lookup(i, p.x + 1, p.y + 1) };
        Float sampPdf;
        p.x += pbrt::SampleDiscrete(wx, u.x, &sampPdf, &u.x);
        pdf *= sampPdf;

        Float wy[2] = { Lookup(i, p.x, p.y), Lookup(i, p.x, p.y + 1) };
        p.y += pbrt::SampleDiscrete(wy, u.y, &sampPdf, &u.y);
        pdf *= sampPdf;
    }

    if (pdfOut != nullptr) *pdfOut = pdf;
    return p;
}

Float Hierarchical2DWarp::DiscretePdf(Point2i p) const {
    Float pdf = 1;

    for (int i = levels.size() - 1; i >= 0; --i) {
        const Array2D<Float> &level = levels[i];

        Point2i pe(p.x & ~1, p.y & ~1);  // even coordinates, rounded down

        Float v = Lookup(i, p.x, p.y);
        if (v == 0) return 0;
        pdf *= v / (Lookup(i, pe.x,     pe.y) +
                    Lookup(i, pe.x + 1, pe.y) +
                    Lookup(i, pe.x,     pe.y + 1) +
                    Lookup(i, pe.x + 1, pe.y + 1));

        if (i > 0 && Resolution(i-1).x < Resolution(i).x) p.x /= 2;
        if (i > 0 && Resolution(i-1).y < Resolution(i).y) p.y /= 2;
    }

    return pdf;
}

Point2f Hierarchical2DWarp::SampleContinuous(Point2f u, Float *pdfOut) const {
    Float pdf = 1;
    Point2i p(0, 0);

    for (size_t i = 0; i < levels.size(); ++i) {
        if (i > 0 && Resolution(i).x > Resolution(i-1).x) p.x *= 2;
        if (i > 0 && Resolution(i).y > Resolution(i-1).y) p.y *= 2;

        Float wx[2] = { Lookup(i, p.x, p.y) + Lookup(i, p.x, p.y + 1),
                        Lookup(i, p.x + 1, p.y) + Lookup(i, p.x + 1, p.y + 1) };
        Float sampPdf;
        p.x += pbrt::SampleDiscrete(wx, u.x, &sampPdf, &u.x);
        pdf *= sampPdf;

        Float wy[2] = { Lookup(i, p.x, p.y), Lookup(i, p.x, p.y + 1) };
        p.y += pbrt::SampleDiscrete(wy, u.y, &sampPdf, &u.y);
        pdf *= sampPdf;
    }

    pdf *= (Resolution().x * Resolution().y) / domain.Area();
    if (pdfOut != nullptr) *pdfOut = pdf;

    return domain.Lerp(Point2f((p.x + u.x) / levels.back().xSize(),
                               (p.y + u.y) / levels.back().ySize()));
}

Float Hierarchical2DWarp::ContinuousPdf(const Point2f &p) const {
    DCHECK(Inside(p, domain));
    Vector2f o = domain.Offset(p);
    return DiscretePdf(Point2i(o.x * Resolution().x, o.y * Resolution().y)) *
        (Resolution().x * Resolution().y) / domain.Area();
}

absl::optional<Point2f> Hierarchical2DWarp::Inverse(const Point2f &p) const {
    DCHECK(Inside(p, domain));

    Bounds2f b;
    b.pMin = domain.pMin;
    b.pMax = domain.Lerp({Float(RoundUpPow2(Resolution().x)) / Float(Resolution().x),
                          Float(RoundUpPow2(Resolution().y)) / Float(Resolution().y)});
    Bounds2f u(Point2f(0, 0), Point2f(1, 1));
    Point2i pi(0, 0);

    for (size_t i = 0; i < levels.size(); ++i) {
        DCHECK(Inside(p, b));
        if (i > 0 && Resolution(i).x > Resolution(i-1).x) pi.x *= 2;
        if (i > 0 && Resolution(i).y > Resolution(i-1).y) pi.y *= 2;

        // update u.x
        if (Resolution(i).x > 1) {
            Float wx[2] = { Lookup(i, pi.x, pi.y) + Lookup(i, pi.x, pi.y + 1),
                            Lookup(i, pi.x + 1, pi.y) + Lookup(i, pi.x + 1, pi.y + 1) };
            Float xMid = (b.pMin.x + b.pMax.x) * 0.5f;
            if (p.x >= xMid) {
                if (wx[1] == 0) return {};
                u.pMin.x = Lerp(wx[0] / (wx[0] + wx[1]), u.pMin.x, u.pMax.x);
                ++pi.x;
                b.pMin.x = xMid;
            } else {
                if (wx[0] == 0) return {};
                u.pMax.x = Lerp(wx[0] / (wx[0] + wx[1]), u.pMin.x, u.pMax.x);
                b.pMax.x = xMid;
            }
        }

        // update u.y
        if (Resolution(i).y > 1) {
            Float wy[2] = { Lookup(i, pi.x, pi.y), Lookup(i, pi.x, pi.y + 1) };
            Float yMid = (b.pMin.y + b.pMax.y) * 0.5f;
            if (p.y >= yMid) {
                if (wy[1] == 0) return {};
                u.pMin.y = Lerp(wy[0] / (wy[0] + wy[1]), u.pMin.y, u.pMax.y);
                ++pi.y;
                b.pMin.y = yMid;
            } else {
                if (wy[0] == 0) return {};
                u.pMax.y = Lerp(wy[0] / (wy[0] + wy[1]), u.pMin.y, u.pMax.y);
                b.pMax.y = yMid;
            }
        }
    }

    Vector2f delta = b.Offset(p);
    return u.Lerp({delta.x, delta.y});
}

std::string Hierarchical2DWarp::ToString() const {
    std::string s = StringPrintf("[ Hierarchical2DWarp nLevels: %d [ ", levels.size());
    for (const Array2D<Float> &level : levels) {
        s += StringPrintf(" nx: %d ny: %d values: [ ", level.xSize(), level.ySize());
        for (int y = 0; y < level.ySize(); ++y)
            for (int x = 0; x < level.xSize(); ++x)
                s += StringPrintf("%f ", level(x, y));
        s += "] ";
    }
    s += "] ";
    return s;
}

}  // namespace pbrt
