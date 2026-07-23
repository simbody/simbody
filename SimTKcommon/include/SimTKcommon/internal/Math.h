#ifndef SimTK_SIMTKCOMMON_MATH_H_
#define SimTK_SIMTKCOMMON_MATH_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors:                                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
 * Low-level constexpr math utilities used internally by SimTK scalar types.
 * These are implementation details; do not rely on them in user code.
 */

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace SimTK {

// Constexpr square root. Uses std::sqrt when __cpp_lib_constexpr_cmath is
// available (C++23), otherwise falls back to Newton-Raphson seeded with the
// IEEE 754 biased-exponent halving trick. The fallback handles all IEEE 754
// special cases but may be off by 1 ulp from the correctly-rounded value.
constexpr double constexpr_sqrt(double x) noexcept {
#if defined(__cpp_lib_constexpr_cmath) && __cpp_lib_constexpr_cmath >= 202202L
    return std::sqrt(x);
#else
    if (x != x)  return x;                                        // NaN
    if (x < 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return x;                                       // preserves -0
    if (x > std::numeric_limits<double>::max()) return x;         // +inf

    // Scale subnormals into the normal range so the biased-exponent trick works.
    double scale = 1.0;
    if (x < std::numeric_limits<double>::min()) {
        x     *= 4503599627370496.0;  // 2^52
        scale  = 1.0 / 67108864.0;   // 2^-26 = 1/sqrt(2^52)
    }

    auto b = std::bit_cast<std::uint64_t>(x);
    b = (b + (1023ULL << 52)) >> 1;        // halve the biased exponent
    auto g = std::bit_cast<double>(b);
    for (;;) {
        const auto next = (g + x / g) * 0.5;
        if (next >= g) return g * scale;
        g = next;
    }
#endif
}

constexpr float constexpr_sqrt(float x) noexcept {
#if defined(__cpp_lib_constexpr_cmath) && __cpp_lib_constexpr_cmath >= 202202L
    return std::sqrt(x);
#else
    if (x != x)  return x;
    if (x < 0.f) return std::numeric_limits<float>::quiet_NaN();
    if (x == 0.f) return x;
    if (x > std::numeric_limits<float>::max()) return x;

    float scale = 1.0f;
    if (x < std::numeric_limits<float>::min()) {
        x     *= 16777216.0f;    // 2^24
        scale  = 1.0f / 4096.0f; // 2^-12 = 1/sqrt(2^24)
    }

    auto b = std::bit_cast<std::uint32_t>(x);
    b = (b + (127U << 23)) >> 1;
    auto g = std::bit_cast<float>(b);
    for (;;) {
        const auto next = (g + x / g) * 0.5f;
        if (next >= g) return g * scale;
        g = next;
    }
#endif
}

} // namespace SimTK

#endif // SimTK_SIMTKCOMMON_MATH_H_
