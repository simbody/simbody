#ifndef SimTK_SIMMATH_GEO_BOX_H_
#define SimTK_SIMMATH_GEO_BOX_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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
typedef Vec<2,P>        Vec2P;
typedef Vec<3,P>        Vec3P;
typedef UnitVec<P,1>    UnitVec3P;
typedef Mat<3,3,P>      Mat33P;
typedef Rotation_<P>    RotationP;
typedef Transform_<P>   TransformP;

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

/** Given a point location in the box frame, return the closest point of
the solid box, and a flag saying whether the given point was inside the
box, using the same definition of "inside" as the containsPoint() method.
Here we define the closest point to an inside point to be the point itself.
Cost is about 9 flops. **/
Vec3P findClosestPointOfSolidBox(const Vec3P& pt, bool& ptWasInside) const {
    Vec3P c(pt);
    ptWasInside = true; // tentatively
    for (int i=0; i<3; ++i) {
        if      (c[i] < -h[i]) {c[i]=-h[i]; ptWasInside=false;}
        else if (c[i] >  h[i]) {c[i]= h[i]; ptWasInside=false;}
    }
    return c;
}

/** Given a point location in the box frame, return the closest point on
the box surface, and a flag saying whether the given point was inside the
box, using the same definition of "inside" as the containsPoint() method.
Cost is about 9 flops for outside points, 18 for inside points. **/
Vec3P findClosestPointOnSurface(const Vec3P& pt, bool& ptWasInside) const {
    Vec3P c = findClosestPointOfSolidBox(pt, ptWasInside);

    if (ptWasInside) { // every |c[i]| <= h[i]
        RealP dToSide = h[0]-std::abs(c[0]); // distance to closest x-face
        int which=0; RealP minDist=dToSide;
        dToSide = h[1]-std::abs(c[1]);
        if (dToSide < minDist) {which=1; minDist=dToSide;}
        dToSide = h[2]-std::abs(c[2]);
        if (dToSide < minDist) {which=2; minDist=dToSide;}
        // Now project the point to the nearest side.
        c[which] = c[which] < 0 ? -h[which] : h[which];
    }
    return c;
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


/** Find a supporting point on the surface of the box in the given direction,
which must be expressed in the box frame. The direction vector does not have
to be a unit vector. The returned point will always be one of the eight
vertices; we treat zeroes here as positive. Consequently if the input vector is
exactly zero, the vertex in the positive orthant is returned as it would be
if the input direction were (1,1,1). Cost is about 5 flops. **/
Vec3P findSupportPoint(const Vec3& d) const {
    // Basically transferring the sign from d to h, but with 0 treated as 1.
    return Vec3P(d[0]<0 ? -h[0]:h[0], d[1]<0 ? -h[1]:h[1], d[2]<0 ? -h[2]:h[2]);
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


/** @name                  Box mesh methods
Methods to use if you want to think of the box as a convex mesh.
**/
static int getNumVertices() {return 8;}
static int getNumEdges()    {return 12;}
static int getNumFaces()    {return 6;}

/** Use bits in the vertex number to pick the signs, with 0=negative,
1=positive:
<pre>
    000  -hx -hy -hz
    001  -hx -hy  hz
    ...
    111   hx  hy  hz
</pre>
**/
Vec3P getVertexPos(int vx) const {
    SimTK_INDEXCHECK(vx,8,"Geo::Box::getVertexPos()");
    return Vec3P(vx&0x4 ? h[0]:-h[0], vx&0x2 ? h[1]:-h[1], vx&0x1 ? h[2]:-h[2]);
}

/** Find the vertex (0-7) that is furthest in the direction d, which is given
in the box frame. Zero coordinates in d are treated as though positive. **/
int findSupportVertex(const Vec3P& d) const {
    int vx = 0;
    if (d[0] >= 0) vx += 0x4; // see table in getVertexPos().
    if (d[1] >= 0) vx += 0x2;
    if (d[2] >= 0) vx += 0x1;
    return vx;
}

/** Vertex normals point diagonally outwards from the box corners. These are
unit vectors with each coordinate +/- 1/sqrt(3). **/
UnitVec3P getVertexNormal(int vx) const {
    SimTK_INDEXCHECK(vx,8,"Geo::Box::getVertexNormal()");
    const RealP c = 1/(RealP)SimTK_SQRT3, nc = -c;
    return UnitVec3P(Vec3P(vx&0x4 ? c:nc, vx&0x2 ? c:nc, vx&0x1 ? c:nc),true);
}

/** Return the center point of the specified edge, in the box frame. **/
Vec3P getEdgeCenter(int ex) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeCenter()");
    int faces[2], which[2];
    getEdgeFaces(ex, faces, which);
    Vec3P c(0);
    for (int i=0; i<2; ++i) {
        const CoordinateDirection fd = getFaceCoordinateDirection(faces[i]);
        c[fd.getAxis()] = fd.getDirection() * h[fd.getAxis()];
    }
    return c;
}

/** Edge normals point diagonally outwards from the edge. These are unit
vectors with one coordinate zero and the others +/- 1/sqrt(2). **/
UnitVec3P getEdgeNormal(int ex) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeNormal()");
    const RealP oosqrt2 = (RealP)SimTK_OOSQRT2;
    int faces[2], which[2];
    getEdgeFaces(ex, faces, which);
    Vec3P n(0);
    for (int i=0; i<2; ++i) {
        const CoordinateDirection fd = getFaceCoordinateDirection(faces[i]);
        n[fd.getAxis()] = fd.getDirection() * oosqrt2;
    }
    return UnitVec3P(n, true);
}



/** Return the direction of an edge, going from its first vertex to its second
vertex, as a CoordinateDirection. **/
CoordinateDirection getEdgeCoordinateDirection(int ex) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeDirection()");
    static const int axis[12] = {2,1,2,1,0,2,0,0,1,1,2,0};
    return CoordinateDirection(CoordinateAxis(axis[ex]), 1); // all in + dir
}


/** Return a unit vector aligned with the selected edge, pointing in the
direction from the first vertex towards the second vertex, in the box
frame. For a box, all edges are aligned with the coordinate system axes so
the returned vector will have only one non-zero component which will be 1
or -1. **/
UnitVec3P getEdgeDirection(int ex) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeUnitVec()");
    return UnitVec3P(getEdgeCoordinateDirection(ex));
}

/** Return the outward normal for the given face as a CoordinateDirection. **/
CoordinateDirection getFaceCoordinateDirection(int fx) const {
    SimTK_INDEXCHECK(fx,6,"Geo::Box::getFaceDirection()");
    static const int axis[6] = {0,  1, 2, 0, 1, 2};
    static const int dir[6]  = {-1,-1,-1, 1, 1, 1};
    return CoordinateDirection(CoordinateAxis(axis[fx]), dir[fx]);
}

/** Return the center point position for the given face. This will be a
vector with one non-zero component of magnitude equal to one of the box
half-dimensions, with the appropriate sign. **/
Vec3P getFaceCenter(int fx) const {
    SimTK_INDEXCHECK(fx,6,"Geo::Box::getFaceCenter()");
    const CoordinateDirection fd = getFaceCoordinateDirection(fx);
    Vec3P c(0);
    c[fd.getAxis()] = fd.getDirection() * h[fd.getAxis()];
    return c;
}

/** Return the outward normal for the given face as a unit vector in the box
frame. This will have only one non-zero component which will be 1 or -1. **/
UnitVec3P getFaceNormal(int fx) const {
    SimTK_INDEXCHECK(fx,6,"Geo::Box::getFaceNormal()");
    return UnitVec3P(getFaceCoordinateDirection(fx));
}

/** A face has four vertices ordered counterclockwise about the face normal. **/
void getFaceVertices(int fx, int v[4]) const {
    SimTK_INDEXCHECK(fx,6,"Geo::Box::getFaceVertices()");
    static const int verts[6][4] = {{0,1,3,2},{0,4,5,1},{0,2,6,4},
                                    {7,5,4,6},{7,6,2,3},{7,3,1,5}};
    for (int i=0; i<4; ++i) v[i] = verts[fx][i];
}

/** Each vertex has three incident faces. Return the face numbers (0-5) and
which of the four face vertices is this one (0-3). **/
void getVertexFaces(int vx, int f[3], int w[3]) const {
    SimTK_INDEXCHECK(vx,8,"Geo::Box::getVertexFaces()");
    static const int faces[8][3] = {{0,1,2},{0,1,5},{0,2,4},{0,4,5},
                                    {1,2,3},{1,3,5},{2,3,4},{3,4,5}};
    static const int which[8][3] = {{0,0,0},{1,3,2},{3,1,2},{2,3,1},
                                    {1,3,2},{2,1,3},{2,3,1},{0,0,0}};
    for (int i=0; i<3; ++i) {f[i]=faces[vx][i]; w[i]=which[vx][i];}
}

/** An edge connects two vertices. **/
void getEdgeVertices(int ex, int v[2]) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeVertices()");
    static const int verts[12][2] = {{0,1},{1,3},{2,3},{0,2},   // 0-3
                                     {0,4},{4,5},{1,5},{2,6},   // 4-7
                                     {4,6},{5,7},{6,7},{3,7}};  // 8-11
    for (int i=0; i<2; ++i) v[i] = verts[ex][i];
}
/** Each vertex has three incident edges. Return the edge numbers (0-11) and
which of the two edge vertices is this one (0-1). **/
void getVertexEdges(int vx, int e[3], int w[3]) const {
    SimTK_INDEXCHECK(vx,8,"Geo::Box::getVertexEdges()");
    static const int edges[8][3] = {{0,3,4},{0,1,6},{2,3,7},{1,2,11},
                                    {4,5,8},{5,6,9},{7,8,10},{9,10,11}};
    static const int which[8][3] = {{0,0,0},{1,0,0},{0,1,0},{1,1,0},
                                    {1,0,0},{1,1,0},{1,1,0},{1,1,1}};
    for (int i=0; i<3; ++i) {e[i]=edges[vx][i]; w[i]=which[vx][i];}
}


/** A face has four edges, ordered by the vertex ordering: v0-v1, v1-v2,
v2-v3, v3-v1. **/
void getFaceEdges(int fx, int e[4]) const {
    SimTK_INDEXCHECK(fx,6,"Geo::Box::getFaceEdges()");
    static const int edges[6][4] = {{0,1,2, 3},{ 4,5,6, 0},{ 3,7,8,4},
                                    {9,5,8,10},{10,7,2,11},{11,1,6,9}};
    for (int i=0; i<4; ++i) e[i] = edges[fx][i];
}
/** An edge is between two faces. Return the face numbers (0-5) and which
one of the four edges on each face is this edge (0-3). **/
void getEdgeFaces(int ex, int f[2], int w[2]) const {
    SimTK_INDEXCHECK(ex,12,"Geo::Box::getEdgeFaces()");
    static const int faces[12][2] = {{0,1},{0,5},{0,4},{0,2},   // 0-3
                                     {1,2},{1,3},{1,5},{2,4},   // 4-7
                                     {2,3},{3,5},{3,4},{4,5}};  // 8-11
    static const int which[12][2] = {{0,3},{1,1},{2,2},{3,0},   // 0-3
                                     {0,3},{1,1},{2,2},{1,1},   // 4-7
                                     {2,2},{0,3},{3,0},{3,0}};  // 8-11
    for (int i=0; i<2; ++i) {f[i] = faces[ex][i]; w[i] = which[ex][i];}
}



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
