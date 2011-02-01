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


#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/OrientedBoundingBox.h"

#include <cassert>

namespace SimTK {

SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactGeometryTypeId);

class ContactGeometryImpl;
class OBBTreeNodeImpl;

/**
 * A ContactGeometry object describes the physical shape of a body.  It is used with GeneralContactSubsystem
 * for doing collision detection and contact modeling.
 */
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
    ContactGeometry() : impl(0) {
    }
    ContactGeometry(const ContactGeometry& src);
    explicit ContactGeometry(ContactGeometryImpl* impl);
    virtual ~ContactGeometry();
    /**
     * Given a point, find the nearest point on the surface of this object.  If multiple points on the
     * surface are equally close to the specified point, this may return any of them.
     *
     * @param position    the point in question
     * @param inside      on exit, this is set to true if the specified point is inside this object, false otherwise
     * @param normal      on exit, this contains the surface normal at the returned point
     * @return the point on the surface of the object which is closest to the specified point
     */
    Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;
    /**
     * Determine whether this object intersects a ray, and if so, find the intersection point.
     *
     * @param origin     the position at which the ray begins
     * @param direction  the ray direction
     * @param distance   if an intersection is found, the distance from the ray origin to the intersection point
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @param normal     if an intersection is found, the surface normal of the intersection point
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @return true if an intersection is found, false otherwise
     */
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const;
    /**
     * Get a bounding sphere which completely encloses this object.
     *
     * @param center     on exit, this contains the location of the center of the bounding sphere
     * @param radius     on exit, this contains the radius of the bounding sphere
     */
    void getBoundingSphere(Vec3& center, Real& radius) const;
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    /**
     * Copy assignment makes a deep copy.
     */
    ContactGeometry& operator=(const ContactGeometry& src);
    bool hasImpl() const {
        return impl != 0;
    }
    const ContactGeometryImpl& getImpl() const {
        assert(impl);
        return *impl;
    }
    ContactGeometryImpl& updImpl() {
        assert(impl);
        return *impl;
    }
    /**
     * Get a string which uniquely identifies the type of geometry this object represents.
     * Typically each subclass of ContactGeometry defines its own value.
     */
    const std::string& getType() const;
    /**
     * Get an integer which uniquely identifies the type of geometry this object represents.
     * A unique index is generated automatically for each unique type value as returned by getType().
     */
    int getTypeIndex() const;

    ContactGeometryTypeId getTypeId() const;

   
protected:
    ContactGeometryImpl* impl;
};

/**
 * This ContactGeometry subclass represents an object that occupies the entire 
 * half-space x>0. This is useful for representing walls and floors.
 */
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

/**
 * This ContactGeometry subclass represents a sphere centered at the origin.
 */
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
    /** Construct an Ellipsoid given its three principal half-axis dimensions
    a,b,c (all positive) along the local x,y,z directions respectively. 
    The curvatures (reciprocals of radii) are precalculated here at a cost
    of about 50 flops. **/
    explicit Ellipsoid(const Vec3& radii);
    /** Obtain the three half-axis dimensions a,b,c used to define this 
    ellipsoid. **/
    const Vec3& getRadii() const;
    /** Set the three half-axis dimensions a,b,c (all positive) used to define
    this ellipsoid, overriding the current radii and recalculating the
    principal curvatures at a cost of about 50 flops. **/
    void setRadii(const Vec3& radii);

    /** For efficiency we precalculate the principal curvatures 
    whenever the ellipsoid radii are set; this avoids having to repeatedly
    perform these three expensive divisions at runtime. The curvatures are 
    ka=1/a, kb=1/b, and kc=1/c so that the ellipsoid's implicit equation can 
    be written Ax^2+By^2+Cz^2=1, with A=ka^2, etc. **/
    const Vec3& getCurvatures() const;

    /** Given a point \a P=(x,y,z) on the ellipsoid surface, return the unique
    unit outward normal to the ellipsoid at that point. If \a P is not
    on the surface, the result is the same as for the point obtained by 
    scaling the vector \a P - O until it just touches the surface. That is, we
    compute P'=findPointInThisDirection(P) and then return the normal at P'.
    Cost is about 50 flops regardless of whether P was initially on the
    surface. @see findPointInSameDirection() **/
    UnitVec3 findUnitNormalAtPoint(const Vec3& P) const;

    /** Given a unit direction \a n, find the unique point P on the ellipsoid 
    surface at which the outward-facing normal is \a n. Cost is about 
    50 flops. **/
    Vec3 findPointWithThisUnitNormal(const UnitVec3& n) const;

    /** Given a direction d defined by the vector Q-O for an arbitrary point
    in space Q=(x,y,z)!=O, find the unique point P on the ellipsoid surface 
    that is in direction d from the ellipsoid origin O. That is, P=s*d for some
    scalar s > 0 such that f(P)=0. Cost is about 50 flops. 
    @param[in]  Q   a point in space measured from the ellipsoid origin but 
                    not the origin
    @return     P, the intersection of the ray in the direction Q-O with the 
                ellipsoid surface **/
    Vec3 findPointInSameDirection(const Vec3& Q) const;

    /** Given a point Q on the surface of the ellipsoid, find the approximating
    paraboloid at Q in a frame P where OP=Q, Pz is the outward-facing unit
    normal to the ellipsoid at Q, Px is the direction of maximum curvature
    and Py is the direction of minimum curvature. k=(kmax,kmin) are the returned
    curvatures with kmax >= kmin > 0. The equation of the resulting paraboloid 
    in the P frame is -2z = kmax*x^2 + kmin*y^2. Cost is about 270 flops. 
    @warning It is up to you to make sure that Q is actually on the ellipsoid
    surface. If it is not you will quietly get a meaningless result.
    @see findParaboloidAtPointWithNormal() **/
    void findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const;

    /** If you already have both a point and the unit normal at that point,
    this will save about 50 flops by trusting that you have provided the 
    correct normal; be careful -- no one is going to check that you got this
    right. The results are meaningless if the point and normal are not 
    consistent. Cost is about 220 flops.
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

    const EllipsoidImpl& getImpl() const;
    EllipsoidImpl& updImpl();
};

/**
 * This ContactGeometry subclass represents an arbitrary shape described by a mesh of triangular faces.
 * The mesh surface must satisfy the following requirements:
 *
 * <ol>
 * <li>It must be closed, so that any point can unambiguously be classified as either inside or outside.</li>
 * <li>It may not intersect itself anywhere, even at a single point.</li>
 * <li>It must be an oriented manifold.</li>
 * <li>The vertices for each face must be ordered counter-clockwise when viewed from the outside.  That is, if
 * v0, v1, and v2 are the locations of the three vertices for a face, the cross product (v1-v0)%(v2-v0) must
 * point outward.</li>
 * <li>The length of every edge must be non-zero.</li>
 * </ol>
 *
 * It is your responsibility to ensure that any mesh you create meets these requirements.  The constructor
 * will detect many incorrect meshes and signal them by throwing an exception, but it is not guaranteed to
 * detect all possible problems.  If a mesh fails to satisfy any of these requirements, the results of calculations
 * performed with it are undefined.  For example, collisions involving it might fail to be detected, or contact
 * forces on it might be calculated incorrectly.
 */
class SimTK_SIMBODY_EXPORT ContactGeometry::TriangleMesh : public ContactGeometry {
public:
    class OBBTreeNode;
    /**
     * Create a TriangleMesh.
     *
     * @param vertices     the positions of all vertices in the mesh
     * @param faceIndices  the indices of the vertices that make up each face.  The first three
     *                     elements are the vertices in the first face, the next three elements are
     *                     the vertices in the second face, etc.
     * @param smooth       if true, the mesh will be treated as a smooth surface, and normal vectors
     *                     will be smoothly interpolated between vertices.  If false, it will be treated
     *                     as a faceted mesh with a constant normal vector over each face.
     */
    TriangleMesh(const ArrayViewConst_<Vec3>& vertices, const ArrayViewConst_<int>& faceIndices, bool smooth=false);
    /**
     * Create a TriangleMesh based on a PolygonalMesh object.  If any faces of the PolygonalMesh
     * have more than three vertices, they are automatically triangulated.
     *
     * @param mesh      the PolygonalMesh from which to construct a triangle mesh
     * @param smooth    if true, the mesh will be treated as a smooth surface, and normal vectors
     *                  will be smoothly interpolated between vertices.  If false, it will be treated
     *                  as a faceted mesh with a constant normal vector over each face.
     */
    explicit TriangleMesh(const PolygonalMesh& mesh, bool smooth=false);
    /**
     * Get the number of edges in the mesh.
     */
    int getNumEdges() const;
    /**
     * Get the number of faces in the mesh.
     */
    int getNumFaces() const;
    /**
     * Get the number of vertices in the mesh.
     */
    int getNumVertices() const;
    /**
     * Get the position of a vertex in the mesh.
     *
     * @param index  the index of the vertex to get
     * @return the position of the specified vertex
     */
    const Vec3& getVertexPosition(int index) const;
    /**
     * Get the index of one of the edges of a face.  Edge 0 connects vertices 0 and 1.
     * Edge 1 connects vertices 1 and 2.  Edge 2 connects vertices 0 and 2.
     *
     * @param face    the index of the face
     * @param edge    the index of the edge within the face (0, 1, or 2)
     * @return the index of the specified edge
     */
    int getFaceEdge(int face, int edge) const;
    /**
     * Get the index of one of the vertices of a face.
     *
     * @param face    the index of the face
     * @param vertex  the index of the vertex within the face (0, 1, or 2)
     * @return the index of the specified vertex
     */
    int getFaceVertex(int face, int vertex) const;
    /**
     * Get the index of one of the faces shared by an edge
     *
     * @param edge    the index of the edge
     * @param face    the index of the face within the edge (0 or 1)
     * @return the index of the specified face
     */
    int getEdgeFace(int edge, int face) const;
    /**
     * Get the index of one of the vertices shared by an edge
     *
     * @param edge    the index of the edge
     * @param vertex  the index of the vertex within the edge (0 or 1)
     * @return the index of the specified vertex
     */
    int getEdgeVertex(int edge, int vertex) const;
    /**
     * Find all edges that intersect a vertex.
     *
     * @param vertex  the index of the vertex
     * @param edges   the indices of all edges intersecting the vertex will be added to this
     */
    void findVertexEdges(int vertex, Array_<int>& edges) const;
    /**
     * Get the normal vector for a face.  This points outward from the mesh.
     *
     * @param face    the index of the face
     */
    const UnitVec3& getFaceNormal(int face) const;
    /**
     * Get the area of a face.
     *
     * @param face    the index of the face
     */
    Real getFaceArea(int face) const;
    /**
     * Calculate the location of a point on the surface, in the local frame of
     * the TriangleMesh. Cost is 11 flops.
     *
     * @param face    the index of the face containing the point
     * @param uv      the point within the face, specified by its barycentric uv coordinates
     */
    Vec3 findPoint(int face, const Vec2& uv) const;
    /**
     * Calculate the location of a face's centroid, that is, the point
     * uv=(1/3,1/3) which is the average of the three vertex locations. This is 
     * a common special case of findPoint() that can be calculated more quickly 
     * (7 flops).
     *
     * @param face    the index of the face whose centroid is of interest
     */
    Vec3 findCentroid(int face) const;
    /**
     * Calculate the normal vector at a point on the surface.
     *
     * @param face    the index of the face containing the point
     * @param uv      the point within the face, specified by its barycentric uv coordinates
     */
    UnitVec3 findNormalAtPoint(int face, const Vec2& uv) const;
    /**
     * Given a point, find the nearest point on the surface of this object.  If multiple points on the
     * surface are equally close to the specified point, this may return any of them.
     *
     * @param position    the point in question
     * @param inside      on exit, this is set to true if the specified point is inside this object, false otherwise
     * @param normal      on exit, this contains the surface normal at the returned point
     * @return the point on the surface of the object which is closest to the specified point
     */
    Vec3 findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const;
    /**
     * Given a point, find the nearest point on the surface of this object.  If multiple points on the
     * surface are equally close to the specified point, this may return any of them.
     *
     * @param position    the point in question
     * @param inside      on exit, this is set to true if the specified point is inside this object, false otherwise
     * @param face        on exit, this contains the index of the face containing the returned point
     * @param uv          on exit, this contains the barycentric coordinates (u and v) of the returned point
     *                    within its face
     * @return the point on the surface of the object which is closest to the specified point
     */
    Vec3 findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const;
    /**
     * Determine whether this mesh intersects a ray, and if so, find the intersection point.
     *
     * @param origin     the position at which the ray begins
     * @param direction  the ray direction
     * @param distance   if an intersection is found, the distance from the ray origin to the intersection point
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @param normal     if an intersection is found, the surface normal of the intersection point
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @return true if an intersection is found, false otherwise
     */
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const;
    /**
     * Determine whether this mesh intersects a ray, and if so, find what face it hit.
     *
     * @param origin     the position at which the ray begins
     * @param direction  the ray direction
     * @param distance   if an intersection is found, the distance from the ray origin to the intersection point
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @param face       if an intersection is found, the index of the face hit by the ray
     *                   is stored in this.  Otherwise, it is left unchanged.
     * @param uv         if an intersection is found, the barycentric coordinates (u and v) of the intersection point
     *                   within the hit face are stored in this.  Otherwise, it is left unchanged.
     * @return true if an intersection is found, false otherwise
     */
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, Vec2& uv) const;
    /**
     * Get the OBBTreeNode which forms the root of this mesh's Oriented Bounding Box Tree.
     */
    OBBTreeNode getOBBTreeNode() const;

    /**
     * Generate a PolygonalMesh from this TriangleMesh; useful mostly for
     * debugging because you can create a DecorativeMesh from this and then
     * look at it.
     */
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

    const TriangleMeshImpl& getImpl() const;
    TriangleMeshImpl& updImpl();
};

/**
 * This class represents a node in the Oriented Bounding Box Tree for a TriangleMesh.  Each node
 * has an OrientedBoundingBox that fully encloses all triangles contained within it or its
 * children.  This is a binary tree: each non-leaf node has two children.  Triangles
 * are stored only in the leaf nodes.
 */

class SimTK_SIMBODY_EXPORT ContactGeometry::TriangleMesh::OBBTreeNode {
public:
    OBBTreeNode(const OBBTreeNodeImpl& impl);
    /**
     * Get the OrientedBoundingBox which encloses all triangles in this node or its children.
     */
    const OrientedBoundingBox& getBounds() const;
    /**
     * Get whether this is a leaf node.
     */
    bool isLeafNode() const;
    /**
     * Get the first child node.  Calling this on a leaf node will produce an exception.
     */
    const OBBTreeNode getFirstChildNode() const;
    /**
     * Get the second child node.  Calling this on a leaf node will produce an exception.
     */
    const OBBTreeNode getSecondChildNode() const;
    /**
     * Get the indices of all triangles contained in this node.  Calling this on a non-leaf node will produce an exception.
     */
    const Array_<int>& getTriangles() const;
    /**
     * Get the number of triangles inside this node.  If this is not a leaf node, this is the total number
     * of triangles contained by all children of this node.
     */
    int getNumTriangles() const;
private:
    const OBBTreeNodeImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_GEOMETRY_H_
