/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2026 Stanford University and the Authors.           *
 * Authors: Adam Kewley                                                       *
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

#include "SimTKcommon/internal/NTraits.h"

#include "SimTKcommon/Testing.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

using namespace SimTK;

namespace
{
    // Returns the ULP distance between two `double`s.
    uint64_t ulp_distance(double a, double b) noexcept
    {
        static_assert(std::numeric_limits<double>::is_iec559);

        if (std::isnan(a) || std::isnan(b)) {
            return std::numeric_limits<uint64_t>::max();
        }

        if (a == b) {
            return 0;
        }

        if (std::signbit(a) != std::signbit(b)) {
            return std::numeric_limits<uint64_t>::max();
        }

        uint64_t ua = std::bit_cast<uint64_t>(a);
        uint64_t ub = std::bit_cast<uint64_t>(b);
        return (ua >= ub) ? (ua - ub) : (ub - ua);
    }
}

// Perform some very rudimentary tests that `consteval_sqrt` effectively
// returns the same as `std::sqrt` on the target platform, but at
// compile-time.
//
// It may be _slightly_ off (1 ULP), though.
void testConstevalSqrtReturnsSameAsSqrt()
{
    static_assert(detail::consteval_sqrt(0.0) == 0.0);
    static_assert(detail::consteval_sqrt(1.0) == 1.0);
    static_assert(detail::consteval_sqrt(4.0) == 2.0);

    {
        constexpr double result = detail::consteval_sqrt(5.0);
        SimTK_TEST(ulp_distance(std::sqrt(5.0), result) <= 1);
    }

    {
        constexpr double result = detail::consteval_sqrt(-10.0);
        SimTK_TEST(std::isnan(result));
    }

    {
        constexpr double result = detail::consteval_sqrt(
            std::numeric_limits<double>::quiet_NaN()
        );
        SimTK_TEST(std::isnan(result));
    }

    {
        constexpr double result = detail::consteval_sqrt(2.0);  // irregular
        SimTK_TEST(ulp_distance(std::sqrt(2.0), result) <= 1);
    }

    {
        constexpr double result = detail::consteval_sqrt(0.5);
        SimTK_TEST(ulp_distance(std::sqrt(0.5), result) <= 1);
    }
}

int main() {

    SimTK_START_TEST("TestNTraits");

        SimTK_SUBTEST(testConstevalSqrtReturnsSameAsSqrt);

    SimTK_END_TEST();
}

