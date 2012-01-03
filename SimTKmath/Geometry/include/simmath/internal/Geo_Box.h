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
/** A 3d rectangular box aligned with an unspecified frame F and centered at 
that frame's origin. The box has a local frame B, centered at the box center 
and oriented along the box edges, and B==F. **/
template <class P>
class SimTK_SIMMATH_EXPORT Geo::Box_ {
typedef P               RealP;
typedef Vec<3,RealP>    Vec3P;

public:
/** Construct an uninitialized Box object; the dimensions will be garbage. **/
Box_() {}
/** Construct a Box with the given nonnegative half-dimensions. **/
Box_(const Vec3P& halfLengths) {setHalfLengths(halfLengths);} 

/** Change the half-dimensions of this box. Dimensions must be nonnegative. **/
Box_& setHalfLengths(const Vec3P& halfLengths) 
{   assert(halfLengths >= 0);
    h = halfLengths; return *this; }

/** Change the half-dimensions of this box by adding the given vector. The
result must be nonnegative. Cost is 3 flops. **/
Box_& addToHalfLengths(const Vec3P& incr) 
{   h += incr; assert(h >= 0); return *this; }

/** Return the half-lengths of this box as a Vec3 from the center to the
first quadrant vertex. **/
const Vec3P& getHalfLengths() const {return h;}

// Don't allow update access to the half lengths in case we want to add
// some precalculations later.

/** Calculate the volume of this box. Cost is 4 flops. **/
RealP findVolume() const {return 8*h[0]*h[1]*h[2];}
/** Calculate the surface area of this box, ignoring degeneracy (meaning that
all pairs of sides are counted even if coincident). Cost is 6 flops. **/
RealP findArea() const {return 8*(h[0]*h[1] + h[0]*h[2] + h[1]*h[2]);}

/** Given a point measured and expressed in the box frame, determine whether
it is strictly contained in the box (just touching doesn't count). The point
must be measured from the box center. Cost is about 5 flops. **/
bool containsPoint(const Vec3P& pt) const {
    const Vec3P absPt = pt.abs(); // reflect to first quadrant
    return absPt < h;
}

private:
Vec3P   h;  // half-dimensions of the box
};



//==============================================================================
//                              GEO ALIGNED BOX
//==============================================================================
/** A 3d box aligned with an unspecified frame F and centered at a given 
point measured from that frame's origin. The box frame B is aligned with F
but the origin Bo is shifted from Fo. **/
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
:   center(center), box(box) {} 
/** Construct an AlignedBox with the given center location and 
half-dimensions. **/
AlignedBox_(const Vec3P& center, const Vec3P& halfLengths) 
:   center(center), box(halfLengths) {} 

/** Change the center location of this box. **/
AlignedBox_& setCenter(const Vec3P& center) 
{   this->center=center; return *this; }

/** Change the dimensions of this box. **/
AlignedBox_& setHalfLengths(const Vec3P& halfLengths) 
{   box.setHalfLengths(halfLengths); return *this; }

/** Return the location of the center of this box (box frame origin). **/
const Vec3P& getCenter() const {return center;}
/** Return a writable reference to the center location of this box. **/
Vec3P& updCenter() {return center;}
/** Return the half-lengths of this box as a Vec3 from the center to the
first quadrant vertex. **/
const Vec3P& getHalfLengths() const {return box.getHalfLengths();}
// no updHalfLengths()
const Box_<P>& getBox() const {return box;}
Box_<P>& updBox() {return box;}

/** Given a point measured and expressed in the base frame F, determine whether
it is strictly contained in the box (just touching doesn't count). Cost is 
about 8 flops. **/
bool containsPoint(const Vec3P& pt_F) const
{   return box.containsPoint(pt_F - center); } // shift to box frame B

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
Vec3P           center;
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


/** Given a point measured and expressed in the base frame F, determine whether
it is strictly contained in the box (just touching doesn't count). Cost is 
about 23 flops. **/
bool containsPoint(const Vec3P& pt_F) const
{   return box.containsPoint(~X_FB*pt_F); } // shift to box frame B

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
