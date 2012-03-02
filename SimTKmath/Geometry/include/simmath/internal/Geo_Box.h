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
and oriented along the box edges, and B==F. We keep track of the relative
edge lengths to facilitate short-to-long processing. **/
template <class P>
class Geo::Box_ {
typedef P               RealP;
typedef Vec<3,P>        Vec3P;
typedef Mat<3,3,P>      Mat33P;
typedef Rotation_<P>    RotationP;
typedef Transform_<P>    TransformP;

public:
/** Construct an uninitialized Box object; the dimensions will be garbage. **/
Box_() {}
/** Construct a Box with the given nonnegative half-dimensions. Cost is 4
flops to sort the edges. **/
Box_(const Vec3P& halfLengths) {setHalfLengths(halfLengths);} 

/** Change the half-dimensions of this box. Dimensions must be nonnegative. 
Cost is 4 flops to sort the edges. **/
Box_& setHalfLengths(const Vec3P& halfLengths) {
    SimTK_ERRCHK3(halfLengths >= 0, "Geo::Box_::setHalfLengths()",
        "Half lengths must be nonnegative; got %g,%g,%g.",
        (double)halfLengths[0],(double)halfLengths[1],(double)halfLengths[2]);
    h = halfLengths; 
    sortEdges();
    return *this; 
}

/** Change the half-dimensions of this box by adding the given vector. The
result must be nonnegative. Cost is 7 flops, including resorting the edges. **/
Box_& addToHalfLengths(const Vec3P& incr) {
    h += incr; 
    SimTK_ERRCHK3(h >= 0, "Geo::Box_::addToHalfLengths()",
        "Half lengths must be nonnegative but were %g,%g,%g after change.",
        (double)h[0],(double)h[1],(double)h[2]);
    sortEdges();
    return *this; 
}

/** Return the half-lengths of this box as a Vec3 from the center to the
first quadrant vertex. **/
const Vec3P& getHalfLengths() const {return h;}

/** Get lengths in order shortest to longest; 0 is shortest, 2 is longest. **/
RealP getOrderedHalfLength(int i) const {
    SimTK_INDEXCHECK(i, 3, "Geo::Box_::getOrderedHalfLength()");
    return h[order[i]];
}

/** Get axes in order shortest to longest; 0 is shortest, 2 is longest. **/
CoordinateAxis getOrderedAxis(int i) const {
    SimTK_INDEXCHECK(i, 3, "Geo::Box_::getOrderedAxis()");
    return CoordinateAxis(order[i]);
}

/** Calculate the volume of this box. Cost is 4 flops. **/
RealP findVolume() const {return 8*h[0]*h[1]*h[2];}
/** Calculate the surface area of this box, ignoring degeneracy (meaning that
all pairs of sides are counted even if coincident). Cost is 6 flops. **/
RealP findArea() const {return 8*(h[0]*h[1] + h[0]*h[2] + h[1]*h[2]);}

/** Given a point measured and expressed in the box frame, determine whether
it is inside the box (we count touching the surface as inside). The point
must be measured from the box center. Cost is about 5 flops. **/
bool containsPoint(const Vec3P& pt) const {
    const Vec3P absPt = pt.abs(); // reflect to first quadrant
    return absPt <= h;
}

/** Return the square of the distance from this box to a given point whose
location is measured from and expressed in the box frame (at the box center). 
If the point is on or inside the box the returned distance is zero. Cost is 
about 14 flops. **/
RealP findDistanceSqrToPoint(const Vec3P& pt) const {
    const Vec3P absPt = pt.abs(); // reflect to first quadrant
    RealP d2 = 0;
    if (absPt[0] > h[0]) d2 += square(absPt[0]-h[0]);
    if (absPt[1] > h[1]) d2 += square(absPt[1]-h[1]);
    if (absPt[2] > h[2]) d2 += square(absPt[2]-h[2]);
    return d2;
}

/** Return the square of the distance from this box to a given sphere whose
center location is measured from and expressed in the box frame (at the box 
center). If the sphere intersects the box the returned distance is zero. Cost 
is about 17 flops. **/
RealP findDistanceSqrToSphere(const Geo::Sphere_<P>& sphere) const {
    const Vec3P absCtr = sphere.getCenter().abs(); // reflect to first quadrant
    const Vec3P grow = h + sphere.getRadius(); // 3 flops
    RealP d2 = 0;
    if (absCtr[0] > grow[0]) d2 += square(absCtr[0]-grow[0]);
    if (absCtr[1] > grow[1]) d2 += square(absCtr[1]-grow[1]);
    if (absCtr[2] > grow[2]) d2 += square(absCtr[2]-grow[2]);
    return d2;
}

/** Return the square of the distance from this box to an axis-aligned box whose
center location is measured from and expressed in this box frame (at the box 
center). If the boxes intersect the returned distance is zero. Cost 
is about 17 flops. **/
RealP findDistanceSqrToAlignedBox(const Geo::AlignedBox_<P>& aab) const {
    const Vec3P absCtr = aab.getCenter().abs(); // reflect to first quadrant
    const Vec3P grow = h + aab.getHalfLengths();
    RealP d2 = 0;
    if (absCtr[0] > grow[0]) d2 += square(absCtr[0]-grow[0]);
    if (absCtr[1] > grow[1]) d2 += square(absCtr[1]-grow[1]);
    if (absCtr[2] > grow[2]) d2 += square(absCtr[2]-grow[2]);
    return d2;
}

/** Given a sphere with center measured and expressed in the box frame, return
true if the box and sphere intersect. We are treating both objects as solids,
so we'll say yes even if one object completely contains the other. We also 
return true if they are just touching. Cost is about 8 flops. **/
bool intersectsSphere(const Geo::Sphere_<P>& sphere) const {
    const Vec3P absCtr = sphere.getCenter().abs(); // reflect to first quadrant
    const RealP r = sphere.getRadius();
    if (absCtr[0] > h[0]+r) return false;
    if (absCtr[1] > h[1]+r) return false;
    if (absCtr[2] > h[2]+r) return false;
    return true;
}

/** Given an aligned box with center measured and expressed in the from of
this box, return true if the two boxes intersect. We are treating both objects 
as solids, so we'll say yes even if one box completely contains the other. We 
also return true if they are just touching. Cost is about 8 flops. **/
bool intersectsAlignedBox(const Geo::AlignedBox_<P>& aab) const {
    const Vec3P absCtr = aab.getCenter().abs(); // reflect to first quadrant
    const Vec3P& aabh = aab.getHalfLengths();
    if (absCtr[0] > h[0]+aabh[0]) return false;
    if (absCtr[1] > h[1]+aabh[1]) return false;
    if (absCtr[2] > h[2]+aabh[2]) return false;
    return true;
}

/** Given an oriented box whose pose is measured and expressed in the frame
of this box, return true if the two boxes intersect. We are treating both 
objects as solids, so we'll say yes even if one box completely contains the 
other. We also return true if they are just touching. This is an exact but
fairly expensive test if the boxes are separated; if you don't mind some
false positives, use mayIntersectOrientedBox() instead. Cost is about 200
flops worst case (when boxes are intersecting) although it can return 
\c false in as few as 16 flops. **/
SimTK_SIMMATH_EXPORT bool 
intersectsOrientedBox(const Geo::OrientedBox_<P>& ob) const;

/** Given an oriented box whose pose is measured and expressed in the frame
of this box, return true if the two boxes may be intersecting. Only relatively
cheap operations are performed at the expense of returning false positives
sometimes (allegedly less than 10% of the time). If you need an exact 
determination, use intersectsOrientedBox(). Cost is about 75 flops worst 
case (when boxes appear to be intersecting) but can return \c false in as few 
as 16 flops. **/
SimTK_SIMMATH_EXPORT bool 
mayIntersectOrientedBox(const Geo::OrientedBox_<P>& ob) const;


private:
// Call this whenever an edge length changes. Each axis will appear once.
void sortEdges() {
    CoordinateAxis shortest = XAxis, longest = ZAxis; 
    if (h[YAxis] < h[shortest]) shortest=YAxis;
    if (h[ZAxis] < h[shortest]) shortest=ZAxis;
    if (h[XAxis] > h[longest])  longest=XAxis;
    if (h[YAxis] > h[longest])  longest=YAxis;
    order[0] = shortest; order[2] = longest; 
    order[1] = shortest.getThirdAxis(longest); // not shortest or longest
}

int intersectsOrientedBoxHelper(const OrientedBox_<P>& O,
                                Mat33P&  absR_BO, 
                                Vec3P&   absP_BO) const;

Vec3P           h;         // half-dimensions of the box
unsigned char   order[3];  // 0,1,2 reordered short to long
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
