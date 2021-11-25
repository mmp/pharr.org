
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SAMPLING_HIERWARP_H
#define PBRT_SAMPLING_HIERWARP_H

// sampling/hierwarp.h*
#include <pbrt/core/types.h>
#include <pbrt/geometry/vecmath.h>
#include <pbrt/geometry/bounds.h>
#include <pbrt/util/array2d.h>
#include <pbrt/util/check.h>

#include <absl/types/optional.h>
#include <absl/types/span.h>
#include <string>
#include <ostream>

namespace pbrt {

class Hierarchical2DWarp {
 public:
    // Take optional Bounds2f domain, use it.
    Hierarchical2DWarp() = default;
    explicit Hierarchical2DWarp(const Array2D<Float> &values)
        : Hierarchical2DWarp(values, values.xSize(), values.ySize()) { }
    Hierarchical2DWarp(const Array2D<Float> &values, const Bounds2f &domain)
        : Hierarchical2DWarp(values, values.xSize(), values.ySize(), domain) { }
    Hierarchical2DWarp(absl::Span<const Float> values, int nx, int ny,
                       const Bounds2f &domain = { Point2f(0, 0), Point2f(1, 1) });

    Point2i SampleDiscrete(Point2f u, Float *pdf = nullptr) const;
    Float DiscretePdf(Point2i p) const;

    Point2f SampleContinuous(Point2f u, Float *pdf = nullptr) const;
    Float ContinuousPdf(const Point2f &p) const;

    absl::optional<Point2f> Inverse(const Point2f &p) const;

    std::string ToString() const;

    Bounds2f Domain() const { return domain; }
    Point2i Resolution() const { return Resolution(levels.size() - 1); }

 private:
    Point2i Resolution(int level) const {
        DCHECK(level >= 0 && level < levels.size());
        return {levels[level].xSize(), levels[level].ySize()};
    }

    Float Lookup(int level, int x, int y) const {
        DCHECK(level >= 0 && level < levels.size());
        DCHECK(x >= 0 && y >= 0);

        if (x >= levels[level].xSize() || y >= levels[level].ySize())
            return 0;
        return levels[level](x, y);
    }

    Bounds2f domain;
    std::vector<Array2D<Float>> levels;
};

inline std::ostream &operator<<(std::ostream &os, const Hierarchical2DWarp &w) {
    return os << w.ToString();
}

}  // namespace pbrt

#endif  // PBRT_SAMPLING_HIERWARP_H

