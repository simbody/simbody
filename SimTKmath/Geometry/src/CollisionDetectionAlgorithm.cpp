/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKmath.h"

#include <set>

using std::map;
using std::pair;
using std::set;

namespace SimTK {

//==============================================================================
//                        COLLISION DETECTION ALGORITHM
//==============================================================================

CollisionDetectionAlgorithm::AlgorithmMap::~AlgorithmMap() {
    // Clean up algorithms to satisfy valgrind
    for (iterator i = begin(); i != end(); ++i)
        delete i->second;
}

CollisionDetectionAlgorithm::AlgorithmMap
    CollisionDetectionAlgorithm::algorithmMap;

static int registerStandardAlgorithms() {
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::HalfSpace::classTypeId(),
        ContactGeometry::Sphere::classTypeId(),
        new CollisionDetectionAlgorithm::HalfSpaceSphere());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::Sphere::classTypeId(),
        ContactGeometry::Sphere::classTypeId(),
        new CollisionDetectionAlgorithm::SphereSphere());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::HalfSpace::classTypeId(),
        ContactGeometry::Ellipsoid::classTypeId(),
        new CollisionDetectionAlgorithm::HalfSpaceEllipsoid());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::Ellipsoid::classTypeId(),
        ContactGeometry::Sphere::classTypeId(),
        new CollisionDetectionAlgorithm::ConvexConvex());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::Ellipsoid::classTypeId(),
        ContactGeometry::Ellipsoid::classTypeId(),
        new CollisionDetectionAlgorithm::ConvexConvex());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::HalfSpace::classTypeId(),
        ContactGeometry::TriangleMesh::classTypeId(),
        new CollisionDetectionAlgorithm::HalfSpaceTriangleMesh());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::Sphere::classTypeId(),
        ContactGeometry::TriangleMesh::classTypeId(),
        new CollisionDetectionAlgorithm::SphereTriangleMesh());
    CollisionDetectionAlgorithm::registerAlgorithm
       (ContactGeometry::TriangleMesh::classTypeId(),
        ContactGeometry::TriangleMesh::classTypeId(),
        new CollisionDetectionAlgorithm::TriangleMeshTriangleMesh());
    return 1;
}

static int staticInitializer = registerStandardAlgorithms();

void CollisionDetectionAlgorithm::registerAlgorithm
   (ContactGeometryTypeId type1, ContactGeometryTypeId type2,
    CollisionDetectionAlgorithm* algorithm) {
    algorithmMap[std::make_pair(type1,type2)] = algorithm;
}

CollisionDetectionAlgorithm* CollisionDetectionAlgorithm::
getAlgorithm(ContactGeometryTypeId type1, ContactGeometryTypeId type2) {
    AlgorithmMap::const_iterator
        iter = algorithmMap.find(std::make_pair(type1, type2));
    if (iter == algorithmMap.end())
        return NULL;
    return iter->second;

}



//==============================================================================
//                             HALF SPACE - SPHERE
//==============================================================================
void CollisionDetectionAlgorithm::HalfSpaceSphere::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1,
    const Transform& transform1,
    ContactSurfaceIndex index2, const ContactGeometry& object2,
    const Transform& transform2,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::Sphere& sphere =
        ContactGeometry::Sphere::getAs(object2);
    // Location of the sphere in the half-space's coordinate frame
    Vec3 location = (~transform1)*transform2.p();
    Real r = sphere.getRadius();
    Real depth = r+location[0];
    if (depth > 0) {
        // They are overlapping.

        Vec3 normal = transform1.R()*Vec3(-1, 0, 0);
        Vec3 contactLocation =
            transform1*Vec3(depth/2, location[1], location[2]);
        contacts.push_back(PointContact(index1, index2, contactLocation,
                                        normal, r, depth));
    }
}



//==============================================================================
//                              SPHERE - SPHERE
//==============================================================================
void CollisionDetectionAlgorithm::SphereSphere::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1,
    const Transform& transform1,
    ContactSurfaceIndex index2, const ContactGeometry& object2,
    const Transform& transform2,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::Sphere& sphere1 =
        ContactGeometry::Sphere::getAs(object1);
    const ContactGeometry::Sphere& sphere2 =
        ContactGeometry::Sphere::getAs(object2);
    Vec3 delta = transform2.p()-transform1.p();
    Real dist = delta.norm();
    if (dist == 0)
        return; // No sensible way to deal with this.
    Real r1 = sphere1.getRadius();
    Real r2 = sphere2.getRadius();
    Real depth = r1+r2-dist;
    if (depth > 0) {
        // They are overlapping.

        Real radius = r1*r2/(r1+r2);
        Vec3 normal = delta/dist;
        Vec3 location = transform1.p()+(r1-depth/2)*normal;
        contacts.push_back(PointContact(index1, index2, location, normal, radius, depth));
    }
}



//==============================================================================
//                            HALF SPACE - ELLIPSOID
//==============================================================================
void CollisionDetectionAlgorithm::HalfSpaceEllipsoid::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1,
    const Transform& transform1,
    ContactSurfaceIndex index2, const ContactGeometry& object2,
    const Transform& transform2,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::Ellipsoid& ellipsoid =
        ContactGeometry::Ellipsoid::getAs(object2);

    const Vec3& radii = ellipsoid.getRadii();
    const Vec3 r2(radii[0]*radii[0], radii[1]*radii[1], radii[2]*radii[2]);
    const Vec3 ri2(1/r2[0], 1/r2[1], 1/r2[2]);
    // Transform giving half space frame in the ellipsoid's frame.
    const Transform trans = (~transform1)*transform2;
    Vec3 normal = ~trans.R()*Vec3(-1, 0, 0);
    Vec3 location(normal[0]*r2[0], normal[1]*r2[1], normal[2]*r2[2]);
    location /= -std::sqrt(  normal[0]*location[0]
                           + normal[1]*location[1]
                           + normal[2]*location[2]);
    Real depth = (trans*location)[0];
    if (depth > 0) {
        // They are overlapping. We need to calculate the principal radii of
        // curvature.

        Vec3 v1 = ~trans.R()*Vec3(0, 1, 0);
        Vec3 v2 = ~trans.R()*Vec3(0, 0, 1);
        Real dxx = v1[0]*v1[0]*ri2[0] + v1[1]*v1[1]*ri2[1] + v1[2]*v1[2]*ri2[2];
        Real dyy = v2[0]*v2[0]*ri2[0] + v2[1]*v2[1]*ri2[1] + v2[2]*v2[2]*ri2[2];
        Real dxy = v1[0]*v2[0]*ri2[0] + v1[1]*v2[1]*ri2[1] + v1[2]*v2[2]*ri2[2];
        Vec<2, complex<Real> > eigenvalues;
        PolynomialRootFinder::findRoots(Vec3(1, -dxx-dyy, dxx*dyy-dxy*dxy),
                                        eigenvalues);
        Vec3 contactNormal = transform2.R()*normal;
        Vec3 contactPoint = transform2*(location+(depth/2)*normal);
        contacts.push_back(PointContact(index1, index2, contactPoint,
                                        contactNormal,
                                        1/std::sqrt(eigenvalues[0].real()),
                                        1/std::sqrt(eigenvalues[1].real()),
                                        depth));
    }
}



//==============================================================================
//                          HALF SPACE - TRIANGLE MESH
//==============================================================================
void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1,
    const Transform& X_GH,
    ContactSurfaceIndex index2, const ContactGeometry& object2,
    const Transform& X_GM,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::TriangleMesh& mesh =
        ContactGeometry::TriangleMesh::getAs(object2);
    // Transform giving mesh (S2) frame in the halfspace (S1) frame.
    const Transform X_HM = (~X_GH)*X_GM;
    set<int> insideFaces;
    Vec3 axisDir = ~X_HM.R()*Vec3(-1, 0, 0);
    Real xoffset = ~axisDir*(~X_HM*Vec3(0));
    processBox(mesh, mesh.getOBBTreeNode(), X_HM, axisDir, xoffset, insideFaces);
    if (insideFaces.size() > 0)
        contacts.push_back(TriangleMeshContact(index1, index2, X_HM,
                                               set<int>(), insideFaces));
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processBox
   (const ContactGeometry::TriangleMesh&              mesh,
    const ContactGeometry::TriangleMesh::OBBTreeNode& node,
    const Transform& X_HM, const Vec3& axisDir, Real xoffset,
    set<int>& insideFaces) const
{   // First check against the node's bounding box.

    OrientedBoundingBox bounds = node.getBounds();
    const Vec3 b = bounds.getSize() / 2;
    Vec3 boxCenter = bounds.getTransform()*b;
    Real radius = ~b*(~bounds.getTransform().R()*axisDir).abs();
    Real dist = ~axisDir*boxCenter-xoffset;
    if (dist > radius)
        return;
    if (dist < -radius) {
        addAllTriangles(node, insideFaces);
        return;
    }

    // If it is not a leaf node, check its children.

    if (!node.isLeafNode()) {
        processBox(mesh, node.getFirstChildNode(), X_HM, axisDir, xoffset, insideFaces);
        processBox(mesh, node.getSecondChildNode(), X_HM, axisDir, xoffset, insideFaces);
        return;
    }

    // Check the triangles.

    const Array_<int>& triangles = node.getTriangles();
    const Row3 xdir = X_HM.R().row(0);
    const Real tx = X_HM.p()[0];
    for (int i = 0; i < (int) triangles.size(); i++) {
        if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 0))+tx > 0)
            insideFaces.insert(triangles[i]);
        else if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 1))+tx > 0)
            insideFaces.insert(triangles[i]);
        else if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 2))+tx > 0)
            insideFaces.insert(triangles[i]);
    }
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::addAllTriangles
   (const ContactGeometry::TriangleMesh::OBBTreeNode& node,
    std::set<int>& insideFaces) const
{
    if (node.isLeafNode()) {
        const Array_<int>& triangles = node.getTriangles();
        for (int i = 0; i < (int) triangles.size(); i++)
            insideFaces.insert(triangles[i]);
    }
    else {
        addAllTriangles(node.getFirstChildNode(), insideFaces);
        addAllTriangles(node.getSecondChildNode(), insideFaces);
    }
}



//==============================================================================
//                           SPHERE - TRIANGLE MESH
//==============================================================================
void CollisionDetectionAlgorithm::SphereTriangleMesh::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1, // sphere
    const Transform& X_GS,
    ContactSurfaceIndex index2, const ContactGeometry& object2, // mesh
    const Transform& X_GM,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::Sphere& sphere =
        ContactGeometry::Sphere::getAs(object1);
    const ContactGeometry::TriangleMesh& mesh =
        ContactGeometry::TriangleMesh::getAs(object2);

    const Transform X_SM = ~X_GS*X_GM;
    // Want the sphere center measured and expressed in the mesh frame.
    const Vec3 p_MC = (~X_SM).p();
    set<int> insideFaces;
    processBox(p_MC, sphere.getRadius()*sphere.getRadius(),
               mesh, mesh.getOBBTreeNode(), insideFaces);
    if (insideFaces.size() > 0)
        contacts.push_back(TriangleMeshContact(index1, index2, X_SM,
                                               set<int>(), insideFaces));
}

void CollisionDetectionAlgorithm::SphereTriangleMesh::processBox
   (const Vec3& center, Real radius2,
    const ContactGeometry::TriangleMesh& mesh,
    const ContactGeometry::TriangleMesh::OBBTreeNode& node,
    set<int>& insideFaces) const
{   // First check against the node's bounding box.

    Vec3 nearestPoint = node.getBounds().findNearestPoint(center);
    if ((nearestPoint-center).normSqr() > radius2)
        return;

    // If it is not a leaf node, check its children.

    if (!node.isLeafNode()) {
        processBox(center, radius2, mesh, node.getFirstChildNode(), insideFaces);
        processBox(center, radius2, mesh, node.getSecondChildNode(), insideFaces);
        return;
    }

    // Check the triangles.

    const Array_<int>& triangles = node.getTriangles();
    for (int i = 0; i < (int) triangles.size(); i++) {
        Vec2 uv;
        Vec3 nearestPoint = mesh.findNearestPointToFace
                                            (center, triangles[i], uv);
        if ((nearestPoint-center).normSqr() < radius2)
            insideFaces.insert(triangles[i]);
    }
}



//==============================================================================
//                        TRIANGLE MESH - TRIANGLE MESH
//==============================================================================
void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::
processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& object1,
    const Transform& X_GM1,
    ContactSurfaceIndex index2, const ContactGeometry& object2,
    const Transform& X_GM2,
    Array_<Contact>& contacts) const
{
    const ContactGeometry::TriangleMesh& mesh1 =
        ContactGeometry::TriangleMesh::getAs(object1);
    const ContactGeometry::TriangleMesh& mesh2 =
        ContactGeometry::TriangleMesh::getAs(object2);

    // Get mesh2's frame measured and expressed in mesh1's frame.
    const Transform X_M1M2 = (~X_GM1)*X_GM2;
    set<int> triangles1;
    set<int> triangles2;
    OrientedBoundingBox mesh2Bounds = X_M1M2*mesh2.getOBBTreeNode().getBounds();
    processNodes(mesh1, mesh2, mesh1.getOBBTreeNode(), mesh2.getOBBTreeNode(),
                 mesh2Bounds, X_M1M2, triangles1, triangles2);
    if (triangles1.size() == 0)
        return; // No intersection.

    // There was an intersection.  We now need to identify every triangle and vertex of each mesh that is inside the other mesh.

    findInsideTriangles(mesh1, mesh2, ~X_M1M2, triangles1);
    findInsideTriangles(mesh2, mesh1,  X_M1M2, triangles2);
    contacts.push_back(TriangleMeshContact(index1, index2, X_M1M2,
                                           triangles1, triangles2));
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::
processNodes
   (const ContactGeometry::TriangleMesh&                mesh1,
    const ContactGeometry::TriangleMesh&                mesh2,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node1,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node2,
    const OrientedBoundingBox&                          node2Bounds,
    const Transform&                                    X_M1M2,
    set<int>&                                           triangles1,
    set<int>&                                           triangles2) const
{   // See if the bounding boxes intersect.

    if (!node1.getBounds().intersectsBox(node2Bounds))
        return;

    // If either node is not a leaf node, process the children recursively.

    if (!node1.isLeafNode()) {
        if (!node2.isLeafNode()) {
            OrientedBoundingBox firstChildBounds = X_M1M2*node2.getFirstChildNode().getBounds();
            OrientedBoundingBox secondChildBounds = X_M1M2*node2.getSecondChildNode().getBounds();
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        }
        else {
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2, node2Bounds, X_M1M2, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2, node2Bounds, X_M1M2, triangles1, triangles2);
        }
        return;
    }
    else if (!node2.isLeafNode()) {
        OrientedBoundingBox firstChildBounds = X_M1M2*node2.getFirstChildNode().getBounds();
        OrientedBoundingBox secondChildBounds = X_M1M2*node2.getSecondChildNode().getBounds();
        processNodes(mesh1, mesh2, node1, node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
        processNodes(mesh1, mesh2, node1, node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        return;
    }

    // These are both leaf nodes, so check triangles for intersections.

    const Array_<int>& node1triangles = node1.getTriangles();
    const Array_<int>& node2triangles = node2.getTriangles();
    for (int i = 0; i < (int) node2triangles.size(); i++) {
        Vec3 a1 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 0));
        Vec3 a2 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 1));
        Vec3 a3 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 2));
        Geo::Triangle A(a1,a2,a3);
        for (int j = 0; j < (int) node1triangles.size(); j++) {
            const Vec3& b1 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 0));
            const Vec3& b2 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 1));
            const Vec3& b3 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 2));
            Geo::Triangle B(b1,b2,b3);
            if (A.overlapsTriangle(B)) {
                // The triangles intersect.
                triangles1.insert(node1triangles[j]);
                triangles2.insert(node2triangles[i]);
            }
        }
    }
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::
findInsideTriangles(const ContactGeometry::TriangleMesh&    mesh,       // M
                    const ContactGeometry::TriangleMesh&    otherMesh,  // O
                    const Transform&                        X_OM,
                    set<int>&                               triangles) const
{   // Find which triangles are inside.
    const int Unknown = UNKNOWN; // work around gcc bug
    Array_<int> faceType(mesh.getNumFaces(), Unknown);
    for (set<int>::iterator iter = triangles.begin();
                            iter != triangles.end(); ++iter)
        faceType[*iter] = BOUNDARY;

    for (int i = 0; i < (int) faceType.size(); i++) {
        if (faceType[i] == UNKNOWN) {
            // Trace a ray from its center to determine whether it is inside.
            const Vec3     origin_O    = X_OM    * mesh.findCentroid(i);
            const UnitVec3 direction_O = X_OM.R()* mesh.getFaceNormal(i);
            Real distance;
            int face;
            Vec2 uv;
            if (   otherMesh.intersectsRay(origin_O, direction_O, distance,
                                           face, uv)
                && ~direction_O*otherMesh.getFaceNormal(face) > 0)
            {
                faceType[i] = INSIDE;
                triangles.insert(i);
            } else
                faceType[i] = OUTSIDE;

            // Recursively mark adjacent triangles.
            tagFaces(mesh, faceType, triangles, i, 0);
        }
    }
}

//TODO: the following method uses depth-first recursion to iterate through
//unmarked faces. For a large mesh this was observed to produce a stack
//overflow in OpenSim. Here we limit the recursion depth; after we get that
//deep we'll pop back out and do another expensive intersectsRay() test in
//the method above.
static const int MaxRecursionDepth = 500;

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::
tagFaces(const ContactGeometry::TriangleMesh&   mesh,
         Array_<int>&                           faceType,
         set<int>&                              triangles,
         int                                    index,
         int                                    depth) const
{
    for (int i = 0; i < 3; i++) {
        int edge = mesh.getFaceEdge(index, i);
        int face = (mesh.getEdgeFace(edge, 0) == index ? mesh.getEdgeFace(edge, 1) : mesh.getEdgeFace(edge, 0));
        if (faceType[face] == UNKNOWN) {
            faceType[face] = faceType[index];
            if (faceType[index] > 0)
                triangles.insert(face);
            if (depth < MaxRecursionDepth)
                tagFaces(mesh, faceType, triangles, face, depth+1);
        }
    }
}



//==============================================================================
//                              CONVEX - CONVEX
//==============================================================================
// This is an implementation based on the Minkowski Portal Refinement method
// used by XenoCollide. See G. Snethen, “Xenocollide: Complex collision made
// simple,” in Game Programming Gems 7, 2008. MPR is used to obtain an initial
// guess for the contact location which is then refined to numerical precision
// by using the smooth representation of the colliding objects.

Vec3 CollisionDetectionAlgorithm::ConvexConvex::
computeSupport(const ContactGeometry& obj1,
               const ContactGeometry& obj2,
               const Transform& transform, UnitVec3 direction) {
    return   obj1.calcSupportPoint(direction)
           - transform*obj2.calcSupportPoint(~transform.R()*-direction);
}

void CollisionDetectionAlgorithm::ConvexConvex::processObjects
   (ContactSurfaceIndex index1, const ContactGeometry& obj1,
    const Transform& transform1,
    ContactSurfaceIndex index2, const ContactGeometry& obj2,
    const Transform& transform2,
    Array_<Contact>& contacts) const
{
    Transform transform = ~transform1*transform2;

    // Compute a point that is known to be inside the Minkowski difference, and
    // a ray directed from that point to the origin.

    Vec3 v0 = (  computeSupport(obj1, obj2, transform, UnitVec3(1, 0, 0))
               + computeSupport(obj1, obj2, transform, UnitVec3(-1, 0, 0))) / 2;
    if (v0 == 0.0) {
        // This is a pathological case: the two objects are directly on top of
        // each other with their centers at exactly the same place. Just
        // return *some* vaguely plausible contact.

        Vec3 point1 = obj1.calcSupportPoint(UnitVec3(1, 0, 0));
        Vec3 point2 = obj2.calcSupportPoint(~transform.R()*UnitVec3(-1, 0, 0));
        addContact(index1, index2, obj1, obj2, transform1, transform2,
                   transform, point1, point2, contacts);
        return;
    }

    // Select three points that define the initial portal.

    UnitVec3 dir1 = UnitVec3(-v0);
    Vec3 v1 = computeSupport(obj1, obj2, transform, dir1);
    if (~v1*dir1 <= 0.0)
        return;
    if (v1%v0 == 0.0) {
        Vec3 point1 = obj1.calcSupportPoint(dir1);
        Vec3 point2 = obj2.calcSupportPoint(~transform.R()*-dir1);
        addContact(index1, index2, obj1, obj2, transform1, transform2,
                   transform, point1, point2, contacts);
        return;
    }
    UnitVec3 dir2 = UnitVec3(v1%v0);
    Vec3 v2 = computeSupport(obj1, obj2, transform, dir2);
    if (~v2*dir2 <= 0.0)
        return;
    UnitVec3 dir3 = UnitVec3((v1-v0)%(v2-v0));
    if (~dir3*v0 > 0) {
        UnitVec3 swap1 = dir1;
        Vec3 swap2 = v1;
        dir1 = dir2;
        v1 = v2;
        dir2 = swap1;
        v2 = swap2;
        dir3 = -dir3;
    }
    Vec3 v3 = computeSupport(obj1, obj2, transform, dir3);
    if (~v3*dir3 <= 0.0)
        return;
    while (true) {
        if (~v0*(v1%v3) < -1e-14) {
            dir2 = dir3;
            v2 = v3;
        }
        else if (~v0*(v3%v2) < -1e-14) {
            dir1 = dir3;
            v1 = v3;
        }
        else
            break;
        dir3 = UnitVec3((v1-v0)%(v2-v0));
        v3 = computeSupport(obj1, obj2, transform, dir3);
    }

    // We have a portal that the origin ray passes through. Now we need to
    // refine it.

    while (true) {
        UnitVec3 portalDir = UnitVec3((v2-v1)%(v3-v1));
        if (~portalDir*v0 > 0)
            portalDir = -portalDir;
        Real dist1 = ~portalDir*v1;
        Vec3 v4 = computeSupport(obj1, obj2, transform, portalDir);
        Real dist4 = ~portalDir*v4;
        if (dist1 >= 0.0) {
            // The origin is inside the portal, so we have an intersection.
            // Compute the barycentric coordinates of the origin in the outer
            // face of the portal.

            Vec3 origin = v0+v0*(~portalDir*(v1-v0)/(~portalDir*v0));
            Real totalArea = ((v2-v1)%(v3-v1)).norm();
            Real area1 = ~portalDir*((v2-origin)%(v3-origin));
            Real area2 = ~portalDir*((v3-origin)%(v1-origin));
            Real u = area1/totalArea;
            Real v = area2/totalArea;
            Real w = 1-u-v;

            // Compute the contact properties.

            Vec3 point1 =  u*obj1.calcSupportPoint(dir1)
                         + v*obj1.calcSupportPoint(dir2)
                         + w*obj1.calcSupportPoint(dir3);
            Vec3 point2 =   u*obj2.calcSupportPoint(~transform.R()*-dir1)
                          + v*obj2.calcSupportPoint(~transform.R()*-dir2)
                          + w*obj2.calcSupportPoint(~transform.R()*-dir3);
            addContact(index1, index2, obj1, obj2,
                       transform1, transform2, transform,
                       point1, point2, contacts);
            return;
        }
        if (dist4 <= 0.0)
            return;
        Vec3 cross = v4%v0;
        if (~v1*cross > 0.0) {
            if (~v2*cross > 0.0) {
                dir1 = portalDir;
                v1 = v4;
            }
            else {
                dir3 = portalDir;
                v3 = v4;
            }
        }
        else {
            if (~v3*cross > 0.0) {
                dir2 = portalDir;
                v2 = v4;
            }
            else {
                dir1 = portalDir;
                v1 = v4;
            }
        }
    }
}

void CollisionDetectionAlgorithm::ConvexConvex::addContact
   (ContactSurfaceIndex index1, ContactSurfaceIndex index2,
    const ContactGeometry& object1,
    const ContactGeometry& object2,
    const Transform& transform1, const Transform& transform2,
    const Transform& transform12,
    Vec3 point1, Vec3 point2, Array_<Contact>& contacts) {
    // We have a rough estimate of the contact points. Use Newton iteration to
    // refine them.

    Vec6 err = computeErrorVector(object1, object2, point1, point2, transform12);
    while (err.norm() > Real(1e-12)) {
        Mat66 J = computeJacobian(object1, object2, point1, point2, transform12);
        FactorQTZ qtz;
        qtz.factor(Matrix(J), Real(1e-6));
        Vector deltaVec(6);
        qtz.solve(Vector(err), deltaVec);
        Vec6 delta(&deltaVec[0]);

        // Line search for safety in case starting guess bad.

        Real f = 2; // scale back factor
        Vec3 point1old = point1, point2old = point2;
        Vec6 errold = err;
        do {
            f /= 2;
            point1 = point1old - f*delta.getSubVec<3>(0);
            point2 = point2old - f*delta.getSubVec<3>(3);
            err = computeErrorVector(object1, object2, point1, point2,
                                     transform12);
        } while (err.norm() > errold.norm());
        if (f < 0.1) {
            // We're clearly outside the region where Newton iteration is going
            // to work properly. Just project the points onto the surfaces and
            // then exit.

            bool inside;
            UnitVec3 normal;
            point1 = object1.findNearestPoint(point1, inside, normal);
            point2 = object2.findNearestPoint(point2, inside, normal);
            break;
        }
    }

    // Compute the curvature of the two surfaces.

    Vec2 curvature1, curvature2;
    Rotation orientation1, orientation2;
    object1.calcCurvature(point1, curvature1, orientation1);
    object2.calcCurvature(point2, curvature2, orientation2);
    Vec2 curvature;
    UnitVec3 maxDir2(transform12.R()*orientation2(0));
    ContactGeometry::combineParaboloids(orientation1, curvature1,
                                        maxDir2, curvature2, curvature);

    // Record the contact.

    Vec3 p1 = transform1*point1;
    Vec3 p2 = transform2*point2;
    Vec3 position = (p1+p2)/2;
    UnitVec3 normal(p1-p2);
    contacts.push_back(PointContact(index1, index2, position, normal,
                                    1/curvature[0], 1/curvature[1],
                                    (p1-p2).norm()));
}

Vec6 CollisionDetectionAlgorithm::ConvexConvex::computeErrorVector
   (const ContactGeometry& object1,
    const ContactGeometry& object2,
    Vec3 pos1, Vec3 pos2, const Transform& transform12) {
    // Compute the function value and normal vector for each object.

    const Function& f1 = object1.getImplicitFunction();
    const Function& f2 = object2.getImplicitFunction();
    Vector x(3);
    Array_<int> components(1);
    Vec3 grad1, grad2;
    Vec3::updAs(&x[0]) = pos1;
    for (int i = 0; i < 3; i++) {
        components[0] = i;
        grad1[i] = f1.calcDerivative(components, x);
    }
    Real error1 = f1.calcValue(x);
    Vec3::updAs(&x[0]) = pos2;
    for (int i = 0; i < 3; i++) {
        components[0] = i;
        grad2[i] = f2.calcDerivative(components, x);
    }
    Real error2 = f2.calcValue(x);

    // Construct a coordinate frame for each object.

    UnitVec3 n1(-grad1);
    UnitVec3 n2(-transform12.R()*grad2);
    UnitVec3 u1(fabs(n1[0]) > 0.5 ? n1%Vec3(0, 1, 0) : n1%Vec3(1, 0, 0));
    UnitVec3 u2(fabs(n2[0]) > 0.5 ? n2%Vec3(0, 1, 0) : n2%Vec3(1, 0, 0));
    Vec3 v1 = n1%u1; // Already a unit vector, so we don't need to normalize it.
    Vec3 v2 = n2%u2;

    // Compute the error vector. The components indicate, in order, that n1
    // must be perpendicular to both tangents of object 2, that the separation
    // vector should be zero or perpendicular to the tangents of object 1, and
    // that both points should be on their respective surfaces.

    Vec3 delta = pos1-transform12*pos2;
    return Vec6(~n1*u2, ~n1*v2, ~delta*u1, ~delta*v1, error1, error2);
}

Mat66 CollisionDetectionAlgorithm::ConvexConvex::computeJacobian
   (const ContactGeometry& object1,
    const ContactGeometry& object2,
    Vec3 pos1, Vec3 pos2, const Transform& transform12) {
    Real dt = Real(1e-7);
    Vec6 err0 = computeErrorVector(object1, object2, pos1, pos2, transform12);
    Vec3 d1 = dt*Vec3(1, 0, 0);
    Vec3 d2 = dt*Vec3(0, 1, 0);
    Vec3 d3 = dt*Vec3(0, 0, 1);
    Vec6 err1 = computeErrorVector(object1, object2, pos1+d1, pos2, transform12)
                - err0;
    Vec6 err2 = computeErrorVector(object1, object2, pos1+d2, pos2, transform12)
                - err0;
    Vec6 err3 = computeErrorVector(object1, object2, pos1+d3, pos2, transform12)
                - err0;
    Vec6 err4 = computeErrorVector(object1, object2, pos1, pos2+d1, transform12)
                - err0;
    Vec6 err5 = computeErrorVector(object1, object2, pos1, pos2+d2, transform12)
                - err0;
    Vec6 err6 = computeErrorVector(object1, object2, pos1, pos2+d3, transform12)
                - err0;
    return Mat66(err1, err2, err3, err4, err5, err6) / dt;
}

} // namespace SimTK

