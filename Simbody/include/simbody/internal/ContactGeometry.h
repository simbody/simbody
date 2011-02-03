#ifndef SimTK_SIMBODY_CONTACT_GEOMETRY_H_
#define SimTK_SIMBODY_CONTACT_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-11 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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
Defines the ContactGeometry class and its API-visible local subclasses for
individual contact shapes. **/

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/OrientedBoundingBox.h"

#include <cassert>

namespace SimTK {

/** @class SimTK::ContactGeometryTypeId
This is a unique integer type for quickly identifying specific types of 
contact geometry for fast lookup purposes. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactGeometryTypeId);

class ContactGeometryImpl;
class OBBTreeNodeImpl;



//==============================================================================
//                             CONTACT GEOMETRY
//==============================================================================
/** A ContactGeometry object describes the physical shape of contact surface.
It is used with GeneralContactSubsystem or ContactTrackerSubsystem for doing
collision detection and contact modeling. This is the base class for the
geometry handles; user code will typically reference one of the local classes
it defines instead for specific shapes. **/
class SimTK_SIMBODY_EXPORT ContactGeometry {
public:
class HalfSpace;
class Sphere;
class Ellipsoid;
class TriangleMesh;
class HalfSpaceImpl;
class SphereImpl;
class EllipsoidImpl;
class TriangleMeshImpl;

/** Base class default constructor creates an empty handle. **/
ContactGeometry() : impl(0) {}
/** Copy constructor makes a deep copy. **/
ContactGeometry(const ContactGeometry& src);
/** Copy assignment makes a deep copy. **/
ContactGeometry& operator=(const ContactGeometry& src);
/** Base class destructor deletes the implementation object.\ Note that this
is not virtual; handles should consist of just a pointer to the 
implementation. **/
~ContactGeometry();

/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param[in]  position    The point in question.
@param[out] inside      On exit, this is set to true if the specified point is 
                        inside this object, false otherwise.
@param[out] normal      On exit, this contains the surface normal at the 
                        returned point.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;

/** Determine whether this object intersects a ray, and if so, find the 
intersection point.
@param[in]  origin      The position at which the ray begins.
@param[in]  direction   The ray direction.
@param[out] distance    If an intersection is found, the distance from the ray 
                        origin to the intersection point is stored in this. 
                        Otherwise, it is left unchanged.
@param[out] normal      If an intersection is found, the surface normal of the
                        intersection point is stored in this. Otherwise, it is 
                        left unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                   Real& distance, UnitVec3& normal) const;

/** Get a bounding sphere which completely encloses this object.
@param[out] center  On exit, this contains the location of the center of the 
                    bounding sphere.
@param[out] radius  On exit, this contains the radius of the bounding 
                    sphere. **/
void getBoundingSphere(Vec3& center, Real& radius) const;

/** This utility method is useful for characterizing the relative geometry of
two locally-smooth surfaces in contact, in a way that is useful for later
application of Hertz compliant contact theory for generating forces. We assume
that contact points Q1 on surface1 and Q2 on surface2 have been determined with
the following properties:
    - the surface normals are aligned but opposite
    - points Q1 and Q2 are separated only along the normal (no tangential 
      separation)

Then the local regions near Q1 and Q2 may be fit with paraboloids P1 and P2
that have their origins at Q1 and Q2, and have the same normals and curvatures
at the origins as do the original surfaces. We will behave here as though
Q1 and Q2 are coincident in space at a point Q; imagine sliding them along
the normal until that happens. Now we define the equations of P1 and P2 in
terms of the maximum and minimum curvatures of surface1 and surface2 at Q:<pre>
    P1: -2z = kmax1 x1^2 + kmin1 y1^2
    P2:  2z = kmax2 x2^2 + kmin2 y2^2  
</pre>
Although the origin Q and z direction are shared, the x,y directions for the 
two paraboloids, though in the same plane z=0, are relatively rotated. Note
that the kmins might be negative; the surfaces do not have to be convex. 

For Hertz contact, we need to know the difference (relative) surface
between the two paraboloids. The difference is a paraboloid P with equation
<pre>
    P: -2z = kmax x^2 + kmin y^2
</pre>
It shares the origin Q and z direction (oriented as for P1), but has its
own principal directions x,y which are coplanar with x1,y1 and x2,y2 but
rotated into some unknown composite orientation. The purpose of this method
is to calculate kmax and kmin, and optionally (depending which signature you
call), x and y, the directions of maximum and minimum curvature (resp.). The
curvature directions are also the principal axes of the contact ellipse formed
by the deformed surfaces, so are necessary (for example) if you want to draw
that ellipse.

Cost is about 220 flops. If you don't need the curvature directions, call the
other overloaded signature which returns only kmax and kmin and takes only 
about 1/3 as long. 

@param[in]          R_SP1  
    The orientation of the P1 paraboloid's frame, expressed in some frame S 
    (typically the frame of the surface to which P1 is fixed). R_SP1.x() is 
    the direction of maximum curvature; y() is minimum curvature; z is the
    contact normal pointing away from surface 1.
@param[in]          k1      
    The maximum (k1[0]) and minimum (k1[1]) curvatures for P1 at the contact 
    point Q1 on surface1. Negative curvatures are handled correctly here but
    may cause trouble for your force model if the resulting contact is 
    conforming.
@param[in]          x2      
    The direction of maximum curvature for paraboloid P2. \a x2 must be in the 
    x1,y1 plane provided in \a R_SP1 and expressed in the S frame.
@param[in]          k2      
    The maximum (k2[0]) and minimum (k2[1]) curvatures for P2 at the contact 
    point Q2 on surface2. Negative curvatures are handled correctly here but
    may cause trouble for your force model if the resulting contact is 
    conforming.
@param[out]         R_SP    
    The orientation of the difference paraboloid P's frame, expressed in the 
    same S frame as was used for P1. R_SP.x() is the direction of maximum 
    curvature of P at the contact point; y() is the minimum curvature 
    direction; z() is the unchanged contact normal pointing away from surface1.
@param[out]         k       
    The maximum (k[0]) and minimum(k[1]) curvatures for the difference 
    paraboloid P at the contact point Q. If either of these is negative or
    zero then the surfaces are conforming and you can't use a point contact
    force model. Note that if k1>0 and k2>0 (i.e. surfaces are convex at Q)
    then k>0 too. If some of the surface curvatures are concave, it is still
    possible that k>0, depending on the relative curvatures.

@see The other signature for combineParaboloids() that is much cheaper if
you just need the curvatures \a k but not the directions \a R_SP. **/
static void combineParaboloids(const Rotation& R_SP1, const Vec2& k1,
                               const UnitVec3& x2, const Vec2& k2,
                               Rotation& R_SP, Vec2& k);

/** This is a much faster version of combineParaboloids() for when you just
need the curvatures of the difference paraboloid, but not the directions 
of those curvatures. Cost is about 70 flops. See the other overload of
this method for details. **/
static void combineParaboloids(const Rotation& R_SP1, const Vec2& k1,
                               const UnitVec3& x2, const Vec2& k2,
                               Vec2& k);

/** Get a string which uniquely identifies the type of geometry this object 
represents. Typically each subclass of ContactGeometry defines its own
value. **/
const std::string& getType() const;
/** Get an integer which uniquely identifies the type of geometry this object
represents. A unique index is generated automatically for each unique type 
value as returned by getType(). **/
int getTypeIndex() const;

/** ContactTrackerSubsystem uses this id for fast identification of specific
surface shapes. **/
ContactGeometryTypeId getTypeId() const;

explicit ContactGeometry(ContactGeometryImpl* impl); /**< Internal use only. **/
bool isOwnerHandle() const;                          /**< Internal use only. **/
bool isEmptyHandle() const;                          /**< Internal use only. **/
bool hasImpl() const {return impl != 0;}             /**< Internal use only. **/
/** Internal use only. **/
const ContactGeometryImpl& getImpl() const {assert(impl); return *impl;}
/** Internal use only. **/
ContactGeometryImpl& updImpl() {assert(impl); return *impl; }

protected:
ContactGeometryImpl* impl; /**< Internal use only. **/
};



//==============================================================================
//                                 HALF SPACE
//==============================================================================
/** This ContactGeometry subclass represents an object that occupies the 
entire half-space x>0. This is useful for representing walls and floors. **/
class SimTK_SIMBODY_EXPORT ContactGeometry::HalfSpace : public ContactGeometry {
public:
HalfSpace();

/** Return true if the supplied ContactGeometry object is a halfspace. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const halfspace. **/
static const HalfSpace& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const HalfSpace&>(geo); }
/** Cast the supplied ContactGeometry object to a writable halfspace. **/
static HalfSpace& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<HalfSpace&>(geo); }

/** Obtain the unique id for HalfSpace contact geometry. **/
static ContactGeometryTypeId classTypeId();
};



//==============================================================================
//                                  SPHERE
//==============================================================================
/** This ContactGeometry subclass represents a sphere centered at the 
origin. **/
class SimTK_SIMBODY_EXPORT ContactGeometry::Sphere : public ContactGeometry {
public:
explicit Sphere(Real radius);
Real getRadius() const;
void setRadius(Real radius);

/** Return true if the supplied ContactGeometry object is a sphere. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const sphere. **/
static const Sphere& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const Sphere&>(geo); }
/** Cast the supplied ContactGeometry object to a writable sphere. **/
static Sphere& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<Sphere&>(geo); }

/** Obtain the unique id for Sphere contact geometry. **/
static ContactGeometryTypeId classTypeId();

const SphereImpl& getImpl() const;
SphereImpl& updImpl();
};



//==============================================================================
//                                  ELLIPSOID
//==============================================================================
/** This ContactGeometry subclass represents an ellipsoid centered at the 
origin, with its principal axes pointing along the x, y, and z axes and
half dimensions a,b, and c (all > 0) along those axes, respectively. The
implicit equation f(x,y,z)=0 of the ellipsoid surface is <pre>
    f(x,y,z) = Ax^2+By^2+Cz^2 - 1
    where A=1/a^2, B=1/b^2, C=1/c^2     
</pre>
A,B, and C are the squares of the principal curvatures ka=1/a, kb=1/b, and
kc=1/c.

The interior of the ellipsoid consists of all points such that f(x,y,z)<0 and
points exterior satisfy f(x,y,z)>0. The region around any point (x,y,z) on an 
ellipsoid surface is locally an elliptic paraboloid with equation <pre>
    -2 z' = kmax x'^2 + kmin y'^2   
</pre>
where z' is measured along the the outward unit normal n at (x,y,z), x' is
measured along the the unit direction u of maximum curvature, and y' is
measured along the unit direction v of minimum curvature. kmax,kmin are the 
curvatures with kmax >= kmin > 0. The signs of the mutually perpendicular
vectors u and v are chosen so that (u,v,n) forms a right-handed coordinate 
system for the paraboloid. **/
class SimTK_SIMBODY_EXPORT ContactGeometry::Ellipsoid : public ContactGeometry {
public:
/** Construct an Ellipsoid given its three principal half-axis dimensions a,b,c
(all positive) along the local x,y,z directions respectively. The curvatures 
(reciprocals of radii) are precalculated here at a cost of about 50 flops. **/
explicit Ellipsoid(const Vec3& radii);
/** Obtain the three half-axis dimensions a,b,c used to define this
ellipsoid. **/
const Vec3& getRadii() const;
/** Set the three half-axis dimensions a,b,c (all positive) used to define this
ellipsoid, overriding the current radii and recalculating the principal 
curvatures at a cost of about 50 flops. 
@param[in] radii    The three half-dimensions of the ellipsoid, in the 
                    ellipsoid's local x, y, and z directions respectively. **/
void setRadii(const Vec3& radii);

/** For efficiency we precalculate the principal curvatures whenever the 
ellipsoid radii are set; this avoids having to repeatedly perform these three
expensive divisions at runtime. The curvatures are ka=1/a, kb=1/b, and kc=1/c 
so that the ellipsoid's implicit equation can be written Ax^2+By^2+Cz^2=1, 
with A=ka^2, etc. **/
const Vec3& getCurvatures() const;

/** Given a point \a P =(x,y,z) on the ellipsoid surface, return the unique unit
outward normal to the ellipsoid at that point. If \a P is not on the surface, 
the result is the same as for the point obtained by scaling the vector 
\a P - O until it just touches the surface. That is, we compute 
P'=findPointInThisDirection(P) and then return the normal at P'. Cost is about
50 flops regardless of whether P was initially on the surface. 
@param[in] P    A point on the ellipsoid surface, measured and expressed in the
                ellipsoid's local frame. See text for what happens if \a P is
                not actually on the ellipsoid surface.
@return The outward-facing unit normal at point \a P (or at the surface point
pointed to by \a P).
@see findPointInSameDirection() **/
UnitVec3 findUnitNormalAtPoint(const Vec3& P) const;

/** Given a unit direction \a n, find the unique point P on the ellipsoid 
surface at which the outward-facing normal is \a n. Cost is about 50 flops. 
@param[in] n    The unit vector for which we want to find a match on the
                ellipsoid surface, expressed in the ellipsoid's local frame. 
@return The point on the ellipsoid's surface at which the outward-facing
normal is the same as \a n. The point is measured and expressed in the 
ellipsoid's local frame. **/
Vec3 findPointWithThisUnitNormal(const UnitVec3& n) const;

/** Given a direction d defined by the vector Q-O for an arbitrary point in 
space Q=(x,y,z)!=O, find the unique point P on the ellipsoid surface that is 
in direction d from the ellipsoid origin O. That is, P=s*d for some scalar 
s > 0 such that f(P)=0. Cost is about 50 flops. 
@param[in] Q    A point in space measured from the ellipsoid origin but not the
                origin.
@return P, the intersection of the ray in the direction Q-O with the ellipsoid
surface **/
Vec3 findPointInSameDirection(const Vec3& Q) const;

/** Given a point Q on the surface of the ellipsoid, find the approximating
paraboloid at Q in a frame P where OP=Q, Pz is the outward-facing unit
normal to the ellipsoid at Q, Px is the direction of maximum curvature
and Py is the direction of minimum curvature. k=(kmax,kmin) are the returned
curvatures with kmax >= kmin > 0. The equation of the resulting paraboloid 
in the P frame is -2z = kmax*x^2 + kmin*y^2. Cost is about 270 flops; you can
save a little time if you already know the normal at Q by using the other
overloaded signature for this method.
 
@warning It is up to you to make sure that Q is actually on the ellipsoid
surface. If it is not you will quietly get a meaningless result.

@param[in]  Q       A point on the surface of this ellipsoid, measured and
                    expressed in the ellipsoid's local frame.
@param[out] X_EP    The frame of the paraboloid P, measured and expressed in
                    the ellipsoid local frame E. X_EP.p() is \a Q, X_EP.x() 
                    is the calculated direction of maximum curvature kmax; y() 
                    is the direction of minimum curvature kmin; z is the 
                    outward facing normal at \a Q.
@param[out] k       The maximum (k[0]) and minimum (k[1]) curvatures of the
                    ellipsoid (and paraboloid P) at point \a Q.
@see findParaboloidAtPointWithNormal() **/
void findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const;

/** If you already have both a point and the unit normal at that point, this 
will save about 50 flops by trusting that you have provided the correct normal;
be careful -- no one is going to check that you got this right. The results are
meaningless if the point and normal are not consistent. Cost is about 220 flops.
@see findParaboloidAtPoint() for details **/
void findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& n,
    Transform& X_EP, Vec2& k) const;

/** Return true if the supplied ContactGeometry object is an Ellipsoid. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const Ellipsoid. **/
static const Ellipsoid& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const Ellipsoid&>(geo); }
/** Cast the supplied ContactGeometry object to a writable Ellipsoid. **/
static Ellipsoid& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<Ellipsoid&>(geo); }

/** Obtain the unique id for Ellipsoid contact geometry. **/
static ContactGeometryTypeId classTypeId();

/** Internal use only. **/
const EllipsoidImpl& getImpl() const;
/** Internal use only. **/
EllipsoidImpl& updImpl();
};



//==============================================================================
//                              TRIANGLE MESH
//==============================================================================
/** This ContactGeometry subclass represents an arbitrary shape described by a 
mesh of triangular faces. The mesh surface must satisfy the following 
requirements:
  - It must be closed, so that any point can unambiguously be classified as 
    either inside or outside.
  - It may not intersect itself anywhere, even at a single point.
  - It must be an oriented manifold.
  - The vertices for each face must be ordered counter-clockwise when viewed
    from the outside. That is, if v0, v1, and v2 are the locations of the 
    three vertices for a face, the cross product (v1-v0)%(v2-v0) must point 
    outward.
  - The length of every edge must be non-zero.

It is your responsibility to ensure that any mesh you create meets these 
requirements. The constructor will detect many incorrect meshes and signal them
by throwing an exception, but it is not guaranteed to detect all possible 
problems.  If a mesh fails to satisfy any of these requirements, the results of
calculations performed with it are undefined. For example, collisions involving
it might fail to be detected, or contact forces on it might be calculated 
incorrectly. **/
class SimTK_SIMBODY_EXPORT ContactGeometry::TriangleMesh 
:   public ContactGeometry {
public:
class OBBTreeNode;
/** Create a TriangleMesh.
@param vertices     The positions of all vertices in the mesh.
@param faceIndices  The indices of the vertices that make up each face. The 
                    first three elements are the vertices in the first face, 
                    the next three elements are the vertices in the second 
                    face, etc.
@param smooth       If true, the mesh will be treated as a smooth surface, and 
                    normal vectors will be smoothly interpolated between 
                    vertices. If false, it will be treated as a faceted mesh 
                    with a constant normal vector over each face. **/
TriangleMesh(const ArrayViewConst_<Vec3>& vertices, const ArrayViewConst_<int>& faceIndices, bool smooth=false);
/** Create a TriangleMesh based on a PolygonalMesh object. If any faces of the 
PolygonalMesh have more than three vertices, they are automatically 
triangulated.
@param mesh      The PolygonalMesh from which to construct a triangle mesh.
@param smooth    If true, the mesh will be treated as a smooth surface, and 
                 normal vectors will be smoothly interpolated between vertices.
                 If false, it will be treated as a faceted mesh with a constant
                 normal vector over each face. **/
explicit TriangleMesh(const PolygonalMesh& mesh, bool smooth=false);
/** Get the number of edges in the mesh. **/
int getNumEdges() const;
/** Get the number of faces in the mesh. **/
int getNumFaces() const;
/** Get the number of vertices in the mesh. **/
int getNumVertices() const;
/** Get the position of a vertex in the mesh.
@param index  The index of the vertex to get.
@return The position of the specified vertex. **/
const Vec3& getVertexPosition(int index) const;
/** Get the index of one of the edges of a face. Edge 0 connects vertices 0 
and 1. Edge 1 connects vertices 1 and 2. Edge 2 connects vertices 0 and 2.
@param face    The index of the face.
@param edge    The index of the edge within the face (0, 1, or 2).
@return The index of the specified edge. **/
int getFaceEdge(int face, int edge) const;
/** Get the index of one of the vertices of a face.
@param face    The index of the face.
@param vertex  The index of the vertex within the face (0, 1, or 2).
@return The index of the specified vertex. **/
int getFaceVertex(int face, int vertex) const;
/** Get the index of one of the faces shared by an edge
@param edge    The index of the edge.
@param face    The index of the face within the edge (0 or 1).
@return The index of the specified face. **/
int getEdgeFace(int edge, int face) const;
/** Get the index of one of the vertices shared by an edge.
@param edge    The index of the edge.
@param vertex  The index of the vertex within the edge (0 or 1).
@return The index of the specified vertex. **/
int getEdgeVertex(int edge, int vertex) const;
/** Find all edges that intersect a vertex.
@param vertex  The index of the vertex.
@param edges   The indices of all edges intersecting the vertex will be added
               to this. **/
void findVertexEdges(int vertex, Array_<int>& edges) const;
/** Get the normal vector for a face. This points outward from the mesh.
@param face    The index of the face. **/
const UnitVec3& getFaceNormal(int face) const;
/** Get the area of a face.
@param face    The index of the face. **/
Real getFaceArea(int face) const;
/** Calculate the location of a point on the surface, in the local frame of
the TriangleMesh. Cost is 11 flops. 
@param face    The index of the face containing the point.
@param uv      The point within the face, specified by its barycentric uv 
               coordinates. **/
Vec3 findPoint(int face, const Vec2& uv) const;
/** Calculate the location of a face's centroid, that is, the point uv=(1/3,1/3)
which is the average of the three vertex locations. This is a common special 
case of findPoint() that can be calculated more quickly (7 flops).
@param face    The index of the face whose centroid is of interest. **/
Vec3 findCentroid(int face) const;
/** Calculate the normal vector at a point on the surface.
@param face    The index of the face containing the point.
@param uv      The point within the face, specified by its barycentric uv 
               coordinates. **/
UnitVec3 findNormalAtPoint(int face, const Vec2& uv) const;
/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param position    The point in question.
@param inside      On exit, this is set to true if the specified point is 
                   inside this object, false otherwise.
@param normal      On exit, this contains the surface normal at the returned 
                   point.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;
/** Given a point, find the nearest point on the surface of this object. If 
multiple points on the surface are equally close to the specified point, this 
may return any of them.
@param position    The point in question.
@param inside      On exit, this is set to true if the specified point is 
                   inside this object, false otherwise.
@param face        On exit, this contains the index of the face containing the 
                   returned point.
@param uv          On exit, this contains the barycentric coordinates (u and v)
                   of the returned point within its face.
@return The point on the surface of the object which is closest to the 
specified point. **/
Vec3 findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const;
/** Determine whether this mesh intersects a ray, and if so, find the 
intersection point.
@param origin     The position at which the ray begins.
@param direction  The ray direction.
@param distance   If an intersection is found, the distance from the ray origin
                  to the intersection point is stored in this. Otherwise, it is
                  left unchanged.
@param normal     If an intersection is found, the surface normal of the 
                  intersection point is stored in this. Otherwise, it is left 
                  unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const;
/** Determine whether this mesh intersects a ray, and if so, find what face it 
hit.
@param origin     The position at which the ray begins.
@param direction  The ray direction.
@param distance   If an intersection is found, the distance from the ray origin
                  to the intersection point is stored in this. Otherwise, it is
                  left unchanged.
@param face       If an intersection is found, the index of the face hit by the
                  ray is stored in this. Otherwise, it is left unchanged.
@param uv         If an intersection is found, the barycentric coordinates (u 
                  and v) of the intersection point within the hit face are 
                  stored in this. Otherwise, it is left unchanged.
@return \c true if an intersection is found, \c false otherwise. **/
bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, Vec2& uv) const;
/** Get the OBBTreeNode which forms the root of this mesh's Oriented Bounding 
Box Tree. **/
OBBTreeNode getOBBTreeNode() const;

/** Generate a PolygonalMesh from this TriangleMesh; useful mostly for debugging
because you can create a DecorativeMesh from this and then look at it. **/
PolygonalMesh createPolygonalMesh() const;

/** Return true if the supplied ContactGeometry object is a triangle mesh. **/
static bool isInstance(const ContactGeometry& geo)
{   return geo.getTypeId()==classTypeId(); }
/** Cast the supplied ContactGeometry object to a const triangle mesh. **/
static const TriangleMesh& getAs(const ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<const TriangleMesh&>(geo); }
/** Cast the supplied ContactGeometry object to a writable triangle mesh. **/
static TriangleMesh& updAs(ContactGeometry& geo)
{   assert(isInstance(geo)); return static_cast<TriangleMesh&>(geo); }

/** Obtain the unique id for TriangleMesh contact geometry. **/
static ContactGeometryTypeId classTypeId();

/** Internal use only. **/
const TriangleMeshImpl& getImpl() const;
/** Internal use only. **/
TriangleMeshImpl& updImpl();
};



//==============================================================================
//                       TRIANGLE MESH :: OBB TREE NODE
//==============================================================================
/** This class represents a node in the Oriented Bounding Box Tree for a 
TriangleMesh. Each node has an OrientedBoundingBox that fully encloses all 
triangles contained within it or its  children. This is a binary tree: each 
non-leaf node has two children. Triangles are stored only in the leaf nodes. **/
class SimTK_SIMBODY_EXPORT ContactGeometry::TriangleMesh::OBBTreeNode {
public:
OBBTreeNode(const OBBTreeNodeImpl& impl);
/** Get the OrientedBoundingBox which encloses all triangles in this node or 
its children. **/
const OrientedBoundingBox& getBounds() const;
/** Get whether this is a leaf node. **/
bool isLeafNode() const;
/** Get the first child node. Calling this on a leaf node will produce an 
exception. **/
const OBBTreeNode getFirstChildNode() const;
/** Get the second child node. Calling this on a leaf node will produce an 
exception. **/
const OBBTreeNode getSecondChildNode() const;
/** Get the indices of all triangles contained in this node. Calling this on a
non-leaf node will produce an exception. **/
const Array_<int>& getTriangles() const;
/** Get the number of triangles inside this node. If this is not a leaf node,
this is the total number of triangles contained by all children of this
node. **/
int getNumTriangles() const;

private:
const OBBTreeNodeImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_GEOMETRY_H_
