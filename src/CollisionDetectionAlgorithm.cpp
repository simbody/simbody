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
    map<pair<int, int>, CollisionDetectionAlgorithm*> algorithmMap = getAlgorithmMap();
    map<pair<int, int>, CollisionDetectionAlgorithm*>::iterator iter = algorithmMap.find(pair<int, int>(typeIndex1, typeIndex2));
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
        contacts.push_back(Contact(index1, index2, contactLocation, normal, contactRadius, depth));
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
        contacts.push_back(Contact(index1, index2, location, normal, contactRadius, depth));
    }
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::TriangleMesh& mesh = static_cast<const ContactGeometry::TriangleMesh&>(object2);
    Transform transform = (~transform1)*transform2; // Transform from the mesh's coordinate frame to the half-space's coordinate frame
    
    // First check against the mesh's bounding box.
    
    OrientedBoundingBox bounds = transform*mesh.getOBBTreeNode().getBounds();
    const Mat33& r = bounds.getTransform().R().asMat33();
    const Vec3 b = 0.5*bounds.getSize();
    const Vec3 meshCenter = bounds.getTransform()*b;
    Real radius = ~b*r.col(0).abs();
    if (meshCenter[0] < -radius)
        return;
    
    // Find the location of each vertex in the half-space's coordinate frame, and record which ones overlap it.
    
    vector<Vec3> vertexPos(mesh.getNumVertices());
    set<int> insideVertices;
    Real maxDepth = 0;
    for (int i = 0; i < vertexPos.size(); i++) {
        vertexPos[i] = transform*mesh.getVertexPosition(i);
        if (vertexPos[i][0] > 0) {
            insideVertices.insert(i);
            maxDepth = std::max(maxDepth, vertexPos[i][0]);
        }
    }
    if (insideVertices.size() == 0)
        return; // No intersection.
    
    // Find connected groups of vertices and record the location and size of each one.
    
    vector<Vec2> points;
    vector<bool> processed(mesh.getNumVertices(), false);
    Vec3 normal = transform1.R()*Vec3(-1, 0, 0);
    while (insideVertices.size() > 0) {
        // Record a set of locations over the surface of the intersection.  This includes the vertices as well as the locations
        // where edges intersect the plane x == 0.
        
        points.clear();
        int firstVertex = *insideVertices.begin();
        processed[firstVertex] = true;
        set<int> insideFaces;
        processVertex(mesh, firstVertex, vertexPos, points, processed, insideVertices, insideFaces);
        
        // Calculate a center position and radius from the points.
        
        Vec2 center(0);
        for (int i = 0; i < points.size(); i++)
            center += points[i];
        center /= points.size();
        Real radius = 0;
        for (int i = 0; i < points.size(); i++)
            radius += (points[i]-center).normSqr();
        radius = std::sqrt(radius/points.size());
        Vec3 contactPoint = transform1*Vec3(0.5*maxDepth, center[0], center[1]);
        contacts.push_back(TriangleMeshContact(index1, index2, contactPoint, normal, radius, maxDepth, set<int>(), insideFaces));
    }
}

void CollisionDetectionAlgorithm::HalfSpaceTriangleMesh::processVertex(const ContactGeometry::TriangleMesh& mesh, int vertex, const vector<Vec3>& vertexPositions, vector<Vec2>& points, vector<bool>& processed, set<int>& insideVertices, std::set<int>& insideFaces) const {
    insideVertices.erase(vertex);
    const Vec3& pos = vertexPositions[vertex];
    points.push_back(Vec2(pos[1], pos[2]));
    vector<int> edges;
    mesh.findVertexEdges(vertex, edges);
    vector<int> needToProcess;
    for (int i = 0; i < edges.size(); i++) {
        insideFaces.insert(mesh.getEdgeFace(edges[i], 0));
        insideFaces.insert(mesh.getEdgeFace(edges[i], 1));
        int otherVertex = (mesh.getEdgeVertex(edges[i], 0) == vertex ? mesh.getEdgeVertex(edges[i], 1) : mesh.getEdgeVertex(edges[i], 0));
        const Vec3& otherVertexPos = vertexPositions[otherVertex];
        if (otherVertexPos[0] < 0) {
            // Add the point where the edge intersects the plane.
            
            Real weight1 = pos[0]/(pos[0]-otherVertexPos[0]);
            Real weight2 = 1.0-weight1;
            points.push_back(Vec2(weight2*pos[1]+weight1*otherVertexPos[1], weight2*pos[2]+weight1*otherVertexPos[2]));
        }
        else if (!processed[otherVertex]) {
            processed[otherVertex] = true;
            needToProcess.push_back(otherVertex);
            
        }
    }
    
    // Recursively process other vertices.
    
    for (int i = 0; i < needToProcess.size(); i++)
        processVertex(mesh, needToProcess[i], vertexPositions, points, processed, insideVertices, insideFaces);
}


void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::processObjects(int index1, const ContactGeometry& object1, const Transform& transform1,
        int index2, const ContactGeometry& object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::TriangleMesh& mesh1 = static_cast<const ContactGeometry::TriangleMesh&>(object1);
    const ContactGeometry::TriangleMesh& mesh2 = static_cast<const ContactGeometry::TriangleMesh&>(object2);
    Transform transform = (~transform1)*transform2; // Transform from mesh2's coordinate frame to mesh1's coordinate frame
    set<int> triangles1;
    set<int> triangles2;
    vector<Vec3> intersectionPoints;
    processNodes(mesh1, mesh2, mesh1.getOBBTreeNode(), mesh2.getOBBTreeNode(), transform, triangles1, triangles2, intersectionPoints);
    if (intersectionPoints.size() == 0)
        return; // No intersection.
    
    // There was an intersection.  We now need to identify every triangle and vertex of each mesh that is inside the other mesh.
    
    vector<int> vertices1;
    vector<int> vertices2;
    findInsideTrianglesAndVertices(mesh1, mesh2, ~transform, triangles1, vertices1);
    findInsideTrianglesAndVertices(mesh2, mesh1, transform, triangles2, vertices2);
    
    // Estimate the contact normal by averaging the face normals.
    
    Vec3 normal(0);
    for (set<int>::iterator iter = triangles1.begin(); iter != triangles1.end(); ++iter)
        normal += mesh1.getFaceNormal(*iter)*mesh1.getFaceArea(*iter);
    for (set<int>::iterator iter = triangles2.begin(); iter != triangles2.end(); ++iter)
        normal -= transform.R()*mesh2.getFaceNormal(*iter)*mesh2.getFaceArea(*iter);
    normal = normal.normalize();
    
    // Estimate the penetration depth.
    
    Real minPosition = MostPositiveReal;
    Real maxPosition = MostNegativeReal;
    for (int i = 0; i < vertices1.size(); i++) {
        Real position = ~normal*mesh1.getVertexPosition(vertices1[i]);
        maxPosition = std::max(maxPosition, position);
    }
    for (int i = 0; i < vertices2.size(); i++) {
        Real position = ~normal*(transform*mesh2.getVertexPosition(vertices2[i]));
        minPosition = std::min(minPosition, position);
    }
    for (int i = 0; i < intersectionPoints.size(); i++) {
        Real position = ~normal*intersectionPoints[i];
        minPosition = std::min(minPosition, position);
        maxPosition = std::max(maxPosition, position);
    }
    Real depth = maxPosition-minPosition;
    
    // Estimate the center and radius from the intersection points.

    Vec3 center(0);
    for (int i = 0; i < intersectionPoints.size(); i++)
        center += intersectionPoints[i];
    center /= intersectionPoints.size();
    Real centerPosition = ~normal*center;
    center += normal*(0.5*(minPosition+maxPosition)-centerPosition);
    Real radius = 0;
    for (int i = 0; i < intersectionPoints.size(); i++) {
        Vec3 delta = intersectionPoints[i]-center;
        delta -= normal*(~delta*normal);
        radius += delta.normSqr();
    }
    radius = std::sqrt(radius/intersectionPoints.size());
    Vec3 contactPoint = transform1*center;
    Vec3 contactNormal = transform1.R()*normal;
    contacts.push_back(TriangleMeshContact(index1, index2, contactPoint, contactNormal, radius, depth, triangles1, triangles2));
}

extern "C" int tri_tri_intersection_test_3d(const Real p1[3], const Real q1[3], const Real r1[3], 
				 const Real p2[3], const Real q2[3], const Real r2[3],
				 int * coplanar, 
				 Real source[3],Real target[3]);

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::processNodes(const ContactGeometry::TriangleMesh& mesh1, const ContactGeometry::TriangleMesh& mesh2,
        const ContactGeometry::TriangleMesh::OBBTreeNode& node1, const ContactGeometry::TriangleMesh::OBBTreeNode& node2,
        const Transform& transform, set<int>& triangles1, set<int>& triangles2, vector<Vec3>& intersectionPoints) const {
    // See if the bounding boxes intersect.
    
    if (!node1.getBounds().intersectsBox(transform*node2.getBounds()))
        return;
    
    // If either node is not a leaf node, process the children recursively.
    
    if (!node1.isLeafNode()) {
        if (!node2.isLeafNode()) {
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getFirstChildNode(), transform, triangles1, triangles2, intersectionPoints);
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2.getSecondChildNode(), transform, triangles1, triangles2, intersectionPoints);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getFirstChildNode(), transform, triangles1, triangles2, intersectionPoints);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2.getSecondChildNode(), transform, triangles1, triangles2, intersectionPoints);
        }
        else {
            processNodes(mesh1, mesh2, node1.getFirstChildNode(), node2, transform, triangles1, triangles2, intersectionPoints);
            processNodes(mesh1, mesh2, node1.getSecondChildNode(), node2, transform, triangles1, triangles2, intersectionPoints);
        }
        return;
    }
    else if (!node2.isLeafNode()) {
        processNodes(mesh1, mesh2, node1, node2.getFirstChildNode(), transform, triangles1, triangles2, intersectionPoints);
        processNodes(mesh1, mesh2, node1, node2.getSecondChildNode(), transform, triangles1, triangles2, intersectionPoints);
        return;
    }
    
    // These are both leaf nodes, so check triangles for intersections.
    
    const vector<int>& node1triangles = node1.getTriangles();
    const vector<int>& node2triangles = node2.getTriangles();
    for (int i = 0; i < node2triangles.size(); i++) {
        Vec3 a1 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 0));
        Vec3 a2 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 1));
        Vec3 a3 = transform*mesh2.getVertexPosition(mesh2.getFaceVertex(node2triangles[i], 2));
        for (int j = 0; j < node1triangles.size(); j++) {
            const Vec3& b1 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 0));
            const Vec3& b2 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 1));
            const Vec3& b3 = mesh1.getVertexPosition(mesh1.getFaceVertex(node1triangles[j], 2));
            int coplanar = 0;
            Vec3 source, target;
            if (tri_tri_intersection_test_3d(&a1[0], &a2[0], &a3[0], &b1[0], &b2[0], &b3[0], &coplanar, &source[0], &target[0])) {
                // The triangles intersect.
                
                triangles1.insert(node1triangles[j]);
                triangles2.insert(node2triangles[i]);
                if (!coplanar) {
                    intersectionPoints.push_back(source);
                    intersectionPoints.push_back(target);
                }
            }
        }
    }
}

void CollisionDetectionAlgorithm::TriangleMeshTriangleMesh::findInsideTrianglesAndVertices(const ContactGeometry::TriangleMesh& mesh, const ContactGeometry::TriangleMesh& otherMesh,
        const Transform& transform, set<int>& triangles, vector<int>& vertices) const {
    // Find which triangles are inside.
    
    vector<int> faceType(mesh.getNumFaces(), UNKNOWN);
    for (set<int>::iterator iter = triangles.begin(); iter != triangles.end(); ++iter)
        faceType[*iter] = BOUNDARY;
    for (int i = 0; i < faceType.size(); i++) {
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
    
    // Now find which vertices are inside.  In most cases, this can be determined from the faces.
    
    vector<int> vertexType(mesh.getNumVertices(), UNKNOWN);
    for (int i = 0; i < faceType.size(); i++) {
        if (faceType[i] == INSIDE || faceType[i] == OUTSIDE) {
            vertexType[mesh.getFaceVertex(i, 0)] = faceType[i];
            vertexType[mesh.getFaceVertex(i, 1)] = faceType[i];
            vertexType[mesh.getFaceVertex(i, 2)] = faceType[i];
        }
    }
    for (int i = 0; i < vertexType.size(); i++) {
        if (vertexType[i] == UNKNOWN) {
            // Trace a ray to find out whether it is inside.
            
            Vec3 origin = transform*mesh.getVertexPosition(i);
            UnitVec3 direction = UnitVec3(1, 0, 0);
            Real distance;
            int face;
            Vec2 uv;
            if (otherMesh.intersectsRay(origin, direction, distance, face, uv) && ~direction*otherMesh.getFaceNormal(face) > 0)
                vertices.push_back(i);
        }
        else if (faceType[i] == INSIDE)
            vertices.push_back(i);
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

