/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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


#include "simbody/internal/CollisionDetectionAlgorithm.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometryImpl.h"
#include <vector>
#include <set>

using std::map;
using std::pair;
using std::set;
using std::vector;

namespace SimTK {

static int registerStandardAlgorithms() {
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::HalfSpaceImpl::Type(), ContactGeometry::SphereImpl::Type(), new CollisionDetectionAlgorithm::HalfSpaceSphere());
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::SphereImpl::Type(), ContactGeometry::SphereImpl::Type(), new CollisionDetectionAlgorithm::SphereSphere());
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::HalfSpaceImpl::Type(), ContactGeometry::TriangleMeshImpl::Type(), new CollisionDetectionAlgorithm::HalfSpaceTriangleMesh());
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::SphereImpl::Type(), ContactGeometry::TriangleMeshImpl::Type(), new CollisionDetectionAlgorithm::SphereTriangleMesh());
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::TriangleMeshImpl::Type(), ContactGeometry::TriangleMeshImpl::Type(), new CollisionDetectionAlgorithm::TriangleMeshTriangleMesh());
    return 1;
}

static int staticInitializer = registerStandardAlgorithms();

map<pair<int, int>, CollisionDetectionAlgorithm*>& CollisionDetectionAlgorithm::getAlgorithmMap() {
    static map<pair<int, int>, CollisionDetectionAlgorithm*> algorithmMap;
    return algorithmMap;
}

void CollisionDetectionAlgorithm::registerAlgorithm(const std::string& type1, const std::string& type2, CollisionDetectionAlgorithm* algorithm) {
    int typeIndex1 = ContactGeometryImpl::getIndexForType(type1);
    int typeIndex2 = ContactGeometryImpl::getIndexForType(type2);
    getAlgorithmMap()[pair<int, int>(typeIndex1, typeIndex2)] = algorithm;
}

CollisionDetectionAlgorithm* CollisionDetectionAlgorithm::getAlgorithm(int typeIndex1, int typeIndex2) {
    const map<pair<int, int>, CollisionDetectionAlgorithm*>& algorithmMap = getAlgorithmMap();
    map<pair<int, int>, CollisionDetectionAlgorithm*>::const_iterator iter = algorithmMap.find(pair<int, int>(typeIndex1, typeIndex2));
    if (iter == algorithmMap.end())
        return NULL;
    return iter->second;

}

void CollisionDetectionAlgorithm::HalfSpaceSphere::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::SphereImpl& sphere = dynamic_cast<const ContactGeometry::SphereImpl&>(object2.getImpl());
    Vec3 location = (~transform1)*transform2.T(); // Location of the sphere in the half-space's coordinate frame
    Real r = sphere.getRadius();
    Real depth = r+location[0];
    if (depth > 0) {
        // They are overlapping.

        Real contactRadius = std::sqrt(depth*r);
        Vec3 normal = transform1.R()*Vec3(-1, 0, 0);
        Vec3 contactLocation = transform1*Vec3(0.5*depth, location[1], location[2]);
        contacts.push_back(PointContact(index1, index2, contactLocation, normal, contactRadius, depth));
    }
}

void CollisionDetectionAlgorithm::SphereSphere::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::SphereImpl& sphere1 = dynamic_cast<const ContactGeometry::SphereImpl&>(object1.getImpl());
    const ContactGeometry::SphereImpl& sphere2 = dynamic_cast<const ContactGeometry::SphereImpl&>(object2.getImpl());
    Vec3 delta = transform2.T()-transform1.T();
    Real dist = delta.norm();
    if (dist == 0)
        return; // No sensible way to deal with this.
    Real r1 = sphere1.getRadius();
    Real r2 = sphere2.getRadius();
    Real depth = r1+r2-dist;
    if (depth > 0) {
        // They are overlapping.
        
        Real curvature = r1*r2/(r1+r2);
        Real contactRadius = std::sqrt(depth*curvature);
        Vec3 normal = delta/dist;
        Vec3 location = transform1.T()+(r1-0.5*depth)*normal;
        contacts.push_back(PointContact(index1, index2, location, normal, contactRadius, depth));
    }
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::TriangleMesh& mesh = static_cast<const ContactGeometry::TriangleMesh&>(object2);
    Transform transform = (~transform1)*transform2; // Transform from the mesh's coordinate frame to the half-space's coordinate frame
    set<int> insideFaces;
    Vec3 axisDir = ~transform.R()*Vec3(-1, 0, 0);
    Real xoffset = ~axisDir*(~transform*Vec3(0));
    processBox(mesh, mesh.getOBBTreeNode(), transform, axisDir, xoffset, insideFaces);
    if (insideFaces.size() > 0)
        contacts.push_back(TriangleMeshContact(index1, index2, set<int>(), insideFaces));
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processBox(const ContactGeometry::TriangleMesh& mesh, const ContactGeometry::TriangleMesh::OBBTreeNode& node, const Transform& transform, const Vec3& axisDir, Real xoffset, set<int>& insideFaces) const {
    // First check against the node's bounding box.
    
    OrientedBoundingBox bounds = node.getBounds();
    const Vec3 b = 0.5*bounds.getSize();
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
        processBox(mesh, node.getFirstChildNode(), transform, axisDir, xoffset, insideFaces);
        processBox(mesh, node.getSecondChildNode(), transform, axisDir, xoffset, insideFaces);
        return;
    }
    
    // Check the triangles.
    
    const vector<int>& triangles = node.getTriangles();
    const Row3 xdir = transform.R().row(0);
    const Real tx = transform.T()[0];
    for (int i = 0; i < (int) triangles.size(); i++) {
        if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 0))+tx > 0)
            insideFaces.insert(triangles[i]);
        else if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 1))+tx > 0)
            insideFaces.insert(triangles[i]);
        else if (xdir*mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], 2))+tx > 0)
            insideFaces.insert(triangles[i]);
    }
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::addAllTriangles(const ContactGeometry::TriangleMesh::OBBTreeNode& node, std::set<int>& insideFaces) const {
    if (node.isLeafNode()) {
        const vector<int>& triangles = node.getTriangles();
        for (int i = 0; i < (int) triangles.size(); i++)
            insideFaces.insert(triangles[i]);
    }
    else {
        addAllTriangles(node.getFirstChildNode(), insideFaces);
        addAllTriangles(node.getSecondChildNode(), insideFaces);
    }
}

void CollisionDetectionAlgorithm::SphereTriangleMesh::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::Sphere& sphere = static_cast<const ContactGeometry::Sphere&>(object1);
    const ContactGeometry::TriangleMesh& mesh = static_cast<const ContactGeometry::TriangleMesh&>(object2);
    Vec3 center = ~transform2*transform1.T();
    set<int> insideFaces;
    processBox(center, sphere.getRadius()*sphere.getRadius(), mesh, mesh.getOBBTreeNode(), insideFaces);
    if (insideFaces.size() > 0)
        contacts.push_back(TriangleMeshContact(index1, index2, set<int>(), insideFaces));
}

void CollisionDetectionAlgorithm::SphereTriangleMesh::processBox(const Vec3& center, Real radius2, const ContactGeometry::TriangleMesh& mesh, const ContactGeometry::TriangleMesh::OBBTreeNode& node, set<int>& insideFaces) const {
    // First check against the node's bounding box.

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
    
    const vector<int>& triangles = node.getTriangles();
    const ContactGeometry::TriangleMeshImpl& impl = mesh.getImpl();
    for (int i = 0; i < (int) triangles.size(); i++) {
        Vec2 uv;
        Vec3 nearestPoint = impl.findNearestPointToFace(center, triangles[i], uv);
        if ((nearestPoint-center).normSqr() < radius2)
            insideFaces.insert(triangles[i]);
    }
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::TriangleMesh& mesh1 = static_cast<const ContactGeometry::TriangleMesh&>(object1);
    const ContactGeometry::TriangleMesh& mesh2 = static_cast<const ContactGeometry::TriangleMesh&>(object2);
    Transform transform = (~transform1)*transform2; // Transform from mesh2's coordinate frame to mesh1's coordinate frame
    set<int> triangles1;
    set<int> triangles2;
    OrientedBoundingBox mesh2Bounds = transform*mesh2.getOBBTreeNode().getBounds();
    processNodes(mesh1, mesh2, mesh1.getOBBTreeNode(), mesh2.getOBBTreeNode(), mesh2Bounds, transform, triangles1, triangles2);
    if (triangles1.size() == 0)
        return; // No intersection.
    
    // There was an intersection.  We now need to identify every triangle and vertex of each mesh that is inside the other mesh.
    
    findInsideTriangles(mesh1, mesh2, ~transform, triangles1);
    findInsideTriangles(mesh2, mesh1, transform, triangles2);
    contacts.push_back(TriangleMeshContact(index1, index2, triangles1, triangles2));
}

extern "C" int tri_tri_overlap_test_3d(const Real p1[3], const Real q1[3], const Real r1[3], 
			    const Real p2[3], const Real q2[3], const Real r2[3]);

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::processNodes(const ContactGeometry::TriangleMesh& mesh1, const ContactGeometry::TriangleMesh& mesh2,
        const ContactGeometry::TriangleMesh::OBBTreeNode& node1, const ContactGeometry::TriangleMesh::OBBTreeNode& node2, const OrientedBoundingBox& node2Bounds,
        const Transform& transform, set<int>& triangles1, set<int>& triangles2) const {
    // See if the bounding boxes intersect.
    
    if (!node1.getBounds().intersectsBox(node2Bounds))
        return;
    
    // If either node is not a leaf node, process the children recursively.
    
    if (!node1.isLeafNode()) {
        if (!node2.isLeafNode()) {
            OrientedBoundingBox firstChildBounds = transform*node2.getFirstChildNode().getBounds();
            OrientedBoundingBox secondChildBounds = transform*node2.getSecondChildNode().getBounds();
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getFirstChildNode(), firstChildBounds, transform, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getSecondChildNode(), secondChildBounds, transform, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getFirstChildNode(), firstChildBounds, transform, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getSecondChildNode(), secondChildBounds, transform, triangles1, triangles2);
        }
        else {
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2, node2Bounds, transform, triangles1, triangles2);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2, node2Bounds, transform, triangles1, triangles2);
        }
        return;
    }
    else if (!node2.isLeafNode()) {
        OrientedBoundingBox firstChildBounds = transform*node2.getFirstChildNode().getBounds();
        OrientedBoundingBox secondChildBounds = transform*node2.getSecondChildNode().getBounds();
        processNodes(mesh1, mesh2, node1, node2.getFirstChildNode(), firstChildBounds, transform, triangles1, triangles2);
        processNodes(mesh1, mesh2, node1, node2.getSecondChildNode(), secondChildBounds, transform, triangles1, triangles2);
        return;
    }
    
    // These are both leaf nodes, so check triangles for intersections.
    
    const vector<int>& node1triangles = node1.getTriangles();
    const vector<int>& node2triangles = node2.getTriangles();
    for (int i = 0; i < (int) node2triangles.size(); i++) {
        Vec3 a1 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 0));
        Vec3 a2 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 1));
        Vec3 a3 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 2));
        for (int j = 0; j < (int) node1triangles.size(); j++) {
            const Vec3& b1 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 0));
            const Vec3& b2 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 1));
            const Vec3& b3 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 2));
            if (tri_tri_overlap_test_3d(&a1[0], &a2[0], &a3[0], &b1[0], &b2[0], &b3[0])) {
                // The triangles intersect.
                
                triangles1.insert(node1triangles[j]);
                triangles2.insert(node2triangles[i]);
            }
        }
    }
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::findInsideTriangles(const ContactGeometry::TriangleMesh& mesh, const ContactGeometry::TriangleMesh& otherMesh,
        const Transform& transform, set<int>& triangles) const {
    // Find which triangles are inside.
    
    vector<int> faceType(mesh.getNumFaces(), UNKNOWN);
    for (set<int>::iterator iter = triangles.begin(); iter != triangles.end(); ++iter)
        faceType[*iter] = BOUNDARY;
    for (int i = 0; i < (int) faceType.size(); i++) {
        if (faceType[i] == UNKNOWN) {
            // Trace a ray from its center to determine whether it is inside.
            
            Vec3 origin = mesh.getVertexPosition(mesh.getFaceVertex(i, 0))+mesh.getVertexPosition(mesh.getFaceVertex(i, 1))+mesh.getVertexPosition(mesh.getFaceVertex(i, 2));
            origin = transform*(origin/3.0);
            UnitVec3 direction = transform.R()*mesh.getFaceNormal(i);
            Real distance;
            int face;
            Vec2 uv;
            if (otherMesh.intersectsRay(origin, direction, distance, face, uv) && ~direction*otherMesh.getFaceNormal(face) > 0) {
                faceType[i] = INSIDE;
                triangles.insert(i);
            }
            else
                faceType[i] = OUTSIDE;
            
            // Recursively mark adjacent triangles.
            
            tagFaces(mesh, faceType, triangles, i);
        }
    }
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::tagFaces(const ContactGeometry::TriangleMesh& mesh, vector<int>& faceType,
        set<int>& triangles, int index) const {
    for (int i = 0; i < 3; i++) {
        int edge = mesh.getFaceEdge(index, i);
        int face = (mesh.getEdgeFace(edge, 0) == index ? mesh.getEdgeFace(edge, 1) : mesh.getEdgeFace(edge, 0));
        if (faceType[face] == UNKNOWN) {
            faceType[face] = faceType[index];
            if (faceType[index] > 0)
                triangles.insert(face);
            tagFaces(mesh, faceType, triangles, face);
        }
    }
}

} // namespace SimTK

