#ifndef SimTK_SIMMATH_GEO_H_
#define SimTK_SIMMATH_GEO_H_

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
Defines geometric primitive shapes and algorthms. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {

//==============================================================================
//                                    GEO
//==============================================================================
/** The Geo class collects geometric primitives intended to deal with raw, 
fixed-size geometric shapes occupying minimal memory and providing
maximum performance through small inline methods and larger high performance
algorithms. Subclasses collect algorithms relevant to particular shapes. 
There are no virtual methods or class hierarchies here; each subclass is a
"POD" (plain old data) class. The general idea is to make it so that these
common methods are implemented in only one place in Simbody.

The Geo class itself is dataless and provides only static methods. It is also
used as a namespace for geometric primitives to allow these names to be used
elsewhere for more significant objects. **/
class SimTK_SIMMATH_EXPORT Geo {
public:
template <class P> class Point_;
template <class P> class Sphere_;
template <class P> class LineSeg_;
template <class P> class Line_;
template <class P> class Plane_;
template <class P> class Circle_;
template <class P> class Box_;
template <class P> class AlignedBox_;
template <class P> class OrientedBox_;
template <class P> class Triangle_;
template <class P> class CubicHermiteCurve_;
template <class P> class BicubicHermitePatch_;
template <class P> class CubicBezierCurve_;
template <class P> class BicubicBezierPatch_;

typedef Point_<Real>        Point;
typedef Sphere_<Real>       Sphere;
typedef LineSeg_<Real>      LineSeg;
typedef Line_<Real>         Line;
typedef Plane_<Real>        Plane;
typedef Circle_<Real>       Circle;
typedef Box_<Real>          Box;
typedef Triangle_<Real>     Triangle;
typedef CubicHermiteCurve_<Real>    CubicHermiteCurve;
typedef BicubicHermitePatch_<Real>  BicubicHermitePatch;
typedef CubicBezierCurve_<Real>     CubicBezierCurve;
typedef BicubicBezierPatch_<Real>   BicubicBezierPatch;

/** Return the default tolerance to use for degeneracy tests and other tests
for "too small" or "near enough" that arise in dealing with geometry primitives.
The value depends on the precision being used; we use the SimTK constant
SignificantReal which is eps^(7/8) where eps is the resolution of P. That
makes this tolerance around 2e-14 in double precision and 9e-7 in float. **/
template <class P> static P getDefaultTol() 
{   return NTraits<P>::getSignificant(); }
/** Return machine precision for floating point calculations at precision P. **/
template <class P> static P getEps() 
{   return NTraits<P>::getEps(); }
/** Return a NaN (not a number) at precision P. **/
template <class P> static P getNaN() 
{   return NTraits<P>::getNaN(); }
/** Return Infinity at precision P. **/
template <class P> static P getInfinity() 
{   return NTraits<P>::getInfinity(); }

/** Stretch a dimension by a given tolerance amount. The result is the
given \a length increased by at least an absolute amount \a tol, or by a
relative amount length*tol if length > 1. Cost is 3 flops. **/
template <class P> static P stretchBy(P length, P tol)
{   assert(tol >= getEps<P>()); 
    return length + std::max(length*tol, tol); }

/** Stretch a dimension using the default tolerance for this precision as
the tolerance in stretchBy(). Cost is 3 flops. **/
template <class P> static P stretch(P length)
{   return stretchBy(length, getDefaultTol<P>()); }

};


} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_H_
