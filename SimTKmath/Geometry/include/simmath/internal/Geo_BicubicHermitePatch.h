#ifndef SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_
#define SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
Provides primitive operations for a single bicubic Hermite patch using either
single or double precision. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {


//==============================================================================
//                         GEO BICUBIC HERMITE PATCH
//==============================================================================
/** A primitive useful for computations involving a single bicubic Hermite
patch. Note that a bicubic Hermite spline surface would not necessarily be
composed of these, but could use the static methods here for patch 
computations. **/
template <class P>
class Geo::BicubicHermitePatch_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized patch; control points will be garbage. **/
BicubicHermitePatch_() {}
/** Construct a bicubic Hermite patch using the given geometry matrix B. **/
explicit BicubicHermitePatch_(const Mat<4,4,Vec3P>& geometry) 
: B(geometry) {} 


/**@name                 Utility methods
These static methods work with given control points. **/
/**@{**/

/**@}**/

private:
Mat<4,4,Vec3P> B;
};



} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_BICUBIC_HERMITE_PATCH_H_
