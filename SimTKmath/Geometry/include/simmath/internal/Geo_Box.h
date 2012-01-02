#ifndef SimTK_SIMMATH_GEO_BOX_H_
#define SimTK_SIMMATH_GEO_BOX_H_

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
Defines primitive operations involving 3d rectangular boxes. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"

#include <cassert>
#include <cmath>
#include <algorithm>

namespace SimTK {

//==============================================================================
//                                  GEO BOX
//==============================================================================
/** TODO: A 3d box aligned with an unspecified frame and centered at that 
frame's origin. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::Box_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized Box object; the dimensions will be garbage. **/
Box_() {}
/** Construct a Box with the given half-dimensions. **/
Box_(const Vec3P& halfLengths) {setHalfLengths(halfLengths);} 

/** Change the half-dimensions of this box. **/
Box_& setHalfLengths(const Vec3P& halfLengths) 
{   assert(halfLengths >= 0);
    h = halfLengths; return *this; }

/** Increase the half-dimensions of this box. **/
Box_& addToHalfLengths(const Vec3P& incr) 
{   assert(incr >= 0);
    h += incr; return *this; }

/** Get the half-dimensions of this box. **/
const Vec3P& getHalfLengths() const {return h;}

private:
Vec3P   h;  // half-dimensions of the box
};



//==============================================================================
//                              GEO ALIGNED BOX
//==============================================================================
/** TODO: A 3d box aligned with an unspecified frame and centered at a given 
point measured from that frame's origin. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::AlignedBox_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized AlignedBox object; the dimensions and location
will be garbage. **/
AlignedBox_() {}
/** Construct an AlignedBox with the given box shape with the center located
as given. **/
AlignedBox_(const Vec3P& center, const Geo::Box_<P>& box) 
: c(center), box(box) {} 
/** Construct an AlignedBox with the given center location and 
half-dimensions. **/
AlignedBox_(const Vec3P& center, const Vec3P& halfLengths) 
: c(center), box(halfLengths) {} 

/** Change the center location of this box. **/
AlignedBox_& setCenter(const Vec3P& center) 
{   c=center; return *this; }

/** Change the dimensions of this box. **/
AlignedBox_& setHalfLengths(const Vec3P& halfLengths) 
{   box.setHalfLengths(halfLengths); return *this; }

const Vec3P& getCenter() const {return c;}
Vec3P& updCenter() {return c;}
const Vec3P& getHalfLengths() const {return box.getHalfLengths();}
// no updHalfLengths()
const Box_<P>& getBox() const {return box;}
Box_<P>& updBox() {return box;}

/** Stretch this box in place by a small amount to ensure that there will 
be no roundoff problems if this is used as a bounding box. The amount to 
stretch depends on the default tolerance for this precision, the dimensions, 
and the position of the box in space. A very large box, or a box that is 
very far from the origin, must be stretched more than a small one at the 
origin. Cost is 6 flops.
@see Geo class for tolerance information. **/
AlignedBox_& stretchBoundary() {
    const RealP tol    = Geo::getDefaultTol<P>();
    const RealP maxdim = max(getCenter().abs());
    const RealP maxrad = max(getHalfLengths());
    const RealP scale  = std::max(maxdim, maxrad);
    const RealP incr   = std::max(scale*Geo::getEps<P>(), tol);
    box.addToHalfLengths(Vec3P(incr));
    return *this; 
}

private:
Vec3P           c;
Geo::Box_<P>    box;
};


//==============================================================================
//                              GEO ORIENTED BOX
//==============================================================================
/** TODO: A 3d box oriented and positioned with respect to an unspecified 
frame F. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::OrientedBox_ {
typedef P               RealP;
typedef Vec<3,P>        Vec3P;
typedef Rotation_<P>    RotationP;
typedef Transform_<P>   TransformP;

public:
/** Construct an uninitialized OrientedBox object; the dimensions and pose
will be garbage. **/
OrientedBox_() {}
/** Construct an OrientedBox with the given box shape with positioned and
oriented according to the given Transform X_FB which gives the box local
frame B (at the box center) in an unspecifed frame F. **/
OrientedBox_(const TransformP& X_FB, const Geo::Box_<P>& box) 
:   X_FB(X_FB), box(box) {} 
/** Construct an OrientedBox with the given location and 
half-dimensions. **/
OrientedBox_(const TransformP& X_FB, const Vec3P& halfLengths) 
:   X_FB(X_FB), box(halfLengths) {} 


/** Change the pose of this box. **/
OrientedBox_& setTransform(const TransformP& newX_FB) 
{   X_FB=newX_FB; return *this; }

/** Change the dimensions of this box. **/
OrientedBox_& setHalfLengths(const Vec3P& halfLengths) 
{   box.setHalfLengths(halfLengths); return *this; }

const Vec3P& getCenter() const {return X_FB.p();}
Vec3P& updCenter() {return X_FB.updP();}
const RotationP& getOrientation() const {return X_FB.R();}
RotationP& updOrientation() {return X_FB.updR();}
const TransformP& getTransform() const {return X_FB;}
TransformP& updTransform() {return X_FB;}
const Vec3P& getHalfLengths() const {return box.getHalfLengths();}
// no updHalfLengths()
const Box_<P>& getBox() const {return box;}
Box_<P>& updBox() {return box;}

/** Stretch this box in place by a small amount to ensure that there will 
be no roundoff problems if this is used as a bounding box. The amount to 
stretch depends on the default tolerance for this precision, the dimensions, 
and the position of the box in space. A very large box, or a box that is 
very far from the origin, must be stretched more than a small one at the 
origin. Cost is 6 flops.
@see Geo class for tolerance information. **/
OrientedBox_& stretchBoundary() {
    const RealP tol    = Geo::getDefaultTol<P>();
    const RealP maxdim = max(getCenter().abs());
    const RealP maxrad = max(getHalfLengths());
    const RealP scale  = std::max(maxdim, maxrad);
    const RealP incr   = std::max(scale*Geo::getEps<P>(), tol);
    box.addToHalfLengths(Vec3P(incr));
    return *this; 
}

private:
TransformP      X_FB;
Geo::Box_<P>    box;
};


} // namespace SimTK

#endif // SimTK_SIMMATH_GEO_BOX_H_
