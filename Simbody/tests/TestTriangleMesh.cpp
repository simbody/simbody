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

#include "SimTKsimbody.h"
#include <vector>
#include <exception>

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

template <class P, int N>
void assertEqual(UnitVec<P,N> val1, UnitVec<P,N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

void testTriangleMesh() {
    // Create a mesh representing a tetrahedron (4 vertices, 4 faces, 6 edges).
    
    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(0, 1, 0));
    vertices.push_back(Vec3(0, 0, 1));
    int faces[4][3] = {{0, 2, 1}, {0, 3, 2}, {0, 1, 3}, {1, 2, 3}};
    vector<int> faceIndices;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);
    ASSERT(mesh.getNumVertices() == 4);
    ASSERT(mesh.getNumFaces() == 4);
    ASSERT(mesh.getNumEdges() == 6);
    
    // Verify that all faces and vertices are correct.
    
    for (int i = 0; i < (int) vertices.size(); i++)
        assertEqual(vertices[i], mesh.getVertexPosition(i));
    for (int i = 0; i < 4; i++) {
        ASSERT(faces[i][0] == mesh.getFaceVertex(i, 0));
        ASSERT(faces[i][1] == mesh.getFaceVertex(i, 1));
        ASSERT(faces[i][2] == mesh.getFaceVertex(i, 2));
    }
    
    // Verify that all indices are consistent.
    
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        for (int j = 0; j < 3; j++) {
            int edge = mesh.getFaceEdge(i, j);
            ASSERT(mesh.getEdgeFace(edge, 0) == i || mesh.getEdgeFace(edge, 1) == i);
            for (int k = 0; k < 2; k++)
                ASSERT(mesh.getEdgeVertex(edge, k) == mesh.getFaceVertex(i, 0) ||
                        mesh.getEdgeVertex(edge, k) == mesh.getFaceVertex(i, 1) ||
                        mesh.getEdgeVertex(edge, k) == mesh.getFaceVertex(i, 2));
        }
    }
    for (int i = 0; i < mesh.getNumEdges(); i++) {
        for (int j = 0; j < 2; j++) {
            int face = mesh.getEdgeFace(i, j);
            ASSERT(mesh.getFaceEdge(face, 0) == i || mesh.getFaceEdge(face, 1) == i || mesh.getFaceEdge(face, 2) == i);
        }
    }
    for (int i = 0; i < mesh.getNumVertices(); i++) {
        Array_<int> edges;
        mesh.findVertexEdges(i, edges);
        ASSERT(edges.size() == 3);
        for (int j = 0; j < (int) edges.size(); j++)
            ASSERT(mesh.getEdgeVertex(edges[j], 0) == i || mesh.getEdgeVertex(edges[j], 1) == i);
    }
    
    // Check the face normals and areas.
    
    assertEqual(mesh.getFaceArea(0), 0.5);
    assertEqual(mesh.getFaceArea(1), 0.5);
    assertEqual(mesh.getFaceArea(2), 0.5);
    assertEqual(mesh.getFaceArea(3), std::sin(Pi/3.0));
    assertEqual(mesh.getFaceNormal(0), Vec3(0, 0, -1));
    assertEqual(mesh.getFaceNormal(1), Vec3(-1, 0, 0));
    assertEqual(mesh.getFaceNormal(2), Vec3(0, -1, 0));
    assertEqual(mesh.getFaceNormal(3), Vec3(1, 1, 1)/Sqrt3);
}

/**
 * Given an invalid mesh, verify that the constructor throws an exception.
 */

void verifyMeshThrowsException(vector<Vec3> vertices, int faces[][3], int numFaces) {
    vector<int> faceIndices;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);
    try {
        ContactGeometry::TriangleMesh mesh(vertices, faceIndices);
        ASSERT(false);
    }
    catch (std::exception&) {
    }
}

void testIncorrectMeshes() {
    vector<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(0, 1, 0));
    vertices.push_back(Vec3(0, 0, 1));
    {
        // The last face has its vertices ordered incorrectly.
        
        int faces[4][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 2, 3}};
        verifyMeshThrowsException(vertices, faces, 4);
    }
    {
        // The last face repeats a vertex.
        
        int faces[4][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 3}};
        verifyMeshThrowsException(vertices, faces, 4);
    }
    {
        // The surface is not closed.
        
        int faces[3][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}};
        verifyMeshThrowsException(vertices, faces, 3);
    }
    {
        // The last face a vertex that is out of range.
        
        int faces[4][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 2, 4}};
        verifyMeshThrowsException(vertices, faces, 4);
    }
}

void addOctohedron(vector<Vec3>& vertices, vector<int>& faceIndices, Vec3 offset) {
    int start = (int)vertices.size();
    vertices.push_back(Vec3(0, 1, 0)+offset);
    vertices.push_back(Vec3(1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, 1)+offset);
    vertices.push_back(Vec3(-1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, -1)+offset);
    vertices.push_back(Vec3(0, -1, 0)+offset);
    int faces[8][3] = {{0, 2, 1}, {0, 3, 2}, {0, 4, 3}, {0, 1, 4}, {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]+start);
}

void validateOBBTree(const ContactGeometry::TriangleMesh& mesh, ContactGeometry::TriangleMesh::OBBTreeNode node, ContactGeometry::TriangleMesh::OBBTreeNode parent, vector<int>& faceReferenceCount) {
    if (node.isLeafNode()) {
        const Array_<int>& triangles = node.getTriangles();
        ASSERT(triangles.size() > 0);
        ASSERT(triangles.size() == node.getNumTriangles());
        for (int i = 0; i < (int) triangles.size(); i++) {
            faceReferenceCount[triangles[i]]++;
            for (int j = 0; j < 3; j++) {
                ASSERT(node.getBounds().containsPoint(mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], j))));
                ASSERT(parent.getBounds().containsPoint(mesh.getVertexPosition(mesh.getFaceVertex(triangles[i], j))));
            }
        }
        try {
            node.getFirstChildNode();
            ASSERT(false); // This should have produced an exception.
        }
        catch (std::exception ex) {
        }
        try {
            node.getFirstChildNode();
            ASSERT(false); // This should have produced an exception.
        }
        catch (std::exception ex) {
        }
    }
    else {
        try {
            node.getTriangles();
            ASSERT(false); // This should have produced an exception.
        }
        catch (std::exception ex) {
        }
        validateOBBTree(mesh, node.getFirstChildNode(), node, faceReferenceCount);
        validateOBBTree(mesh, node.getSecondChildNode(), node, faceReferenceCount);
        ASSERT(node.getNumTriangles() == node.getFirstChildNode().getNumTriangles()+node.getSecondChildNode().getNumTriangles());
    }
}

void testOBBTree() {
    // Create a mesh consisting of a bunch of octohedra.
    
    vector<Vec3> vertices;
    vector<int> faceIndices;
    addOctohedron(vertices, faceIndices, Vec3(0, 0, 0));
    addOctohedron(vertices, faceIndices, Vec3(2.5, 0, 0));
    addOctohedron(vertices, faceIndices, Vec3(0, 2.5, 0));
    addOctohedron(vertices, faceIndices, Vec3(2.5, 2.5, 0));
    addOctohedron(vertices, faceIndices, Vec3(0, 0, 2.5));
    addOctohedron(vertices, faceIndices, Vec3(2.5, 0, 2.5));
    addOctohedron(vertices, faceIndices, Vec3(0, 2.5, 2.5));
    addOctohedron(vertices, faceIndices, Vec3(2.5, 2.5, 2.5));
    addOctohedron(vertices, faceIndices, Vec3(1.25, 1.25, 1.25));
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

    // Validate the OBBTree.
    
    vector<int> faceReferenceCount(mesh.getNumFaces(), 0);
    validateOBBTree(mesh, mesh.getOBBTreeNode(), mesh.getOBBTreeNode(), faceReferenceCount);
    for (int i = 0; i < (int) faceReferenceCount.size(); i++)
        ASSERT(faceReferenceCount[i] == 1);
}

void testRayIntersection() {
    // Create an octrohedral mesh.
    
    vector<Vec3> vertices;
    vector<int> faceIndices;
    addOctohedron(vertices, faceIndices, Vec3(0, 0, 0));
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);
    
    // Check various rays.
    
    Real distance;
    UnitVec3 normal;
    ASSERT(!mesh.intersectsRay(Vec3(2, 0, 0), UnitVec3(1, 0, 0), distance, normal));
    ASSERT(mesh.intersectsRay(Vec3(2, 0, 0), UnitVec3(-1, 0, 0), distance, normal));
    assertEqual(1.0, distance);
    ASSERT(mesh.intersectsRay(Vec3(0, 0, 0), UnitVec3(1, 1, 1), distance, normal));
    assertEqual(1.0/Sqrt3, distance);
    assertEqual(normal, Vec3(1, 1, 1)/Sqrt3);
    ASSERT(mesh.intersectsRay(Vec3(0.1, 0.1, 0.1), UnitVec3(-1, -1, -1), distance, normal));
    assertEqual(std::sqrt(3*0.1*0.1)+1.0/Sqrt3, distance);
    assertEqual(normal, Vec3(-1, -1, -1)/Sqrt3);
    ASSERT(!mesh.intersectsRay(Vec3(-1, -1, -1), UnitVec3(-1, -1, -1), distance, normal));
}

void testSmoothMesh() {
    // Create two octrohedral meshes: one smooth and one not.
    
    vector<Vec3> vertices;
    vector<int> faceIndices;
    addOctohedron(vertices, faceIndices, Vec3(0, 0, 0));
    ContactGeometry::TriangleMesh mesh1(vertices, faceIndices);
    ContactGeometry::TriangleMesh mesh2(vertices, faceIndices, true);
    
    // At the center of every face, the normals should be identical.
 
    for (int i = 0; i < 8; i++)
        assertEqual(mesh1.findNormalAtPoint(i, Vec2(1.0/3.0, 1.0/3.0)), mesh2.findNormalAtPoint(i, Vec2(1.0/3.0, 1.0/3.0)));
    
    // At the vertices, the smooth mesh normal should point directly outward, while the faceted mesh should
    // be the same as the face normal.
    
    for (int i = 0; i < 8; i++) {
        assertEqual(mesh1.findNormalAtPoint(i, Vec2(1, 0)), mesh1.getFaceNormal(i));
        assertEqual(mesh1.findNormalAtPoint(i, Vec2(0, 1)), mesh1.getFaceNormal(i));
        assertEqual(mesh1.findNormalAtPoint(i, Vec2(0, 0)), mesh1.getFaceNormal(i));
        assertEqual(mesh2.findNormalAtPoint(i, Vec2(1, 0)), mesh2.getVertexPosition(mesh2.getFaceVertex(i, 0)));
        assertEqual(mesh2.findNormalAtPoint(i, Vec2(0, 1)), mesh2.getVertexPosition(mesh2.getFaceVertex(i, 1)));
        assertEqual(mesh2.findNormalAtPoint(i, Vec2(0, 0)), mesh2.getVertexPosition(mesh2.getFaceVertex(i, 2)));
    }
}

void testFindNearestPoint() {
    // Create an octrohedral mesh.
    
    vector<Vec3> vertices;
    vector<int> faceIndices;
    addOctohedron(vertices, faceIndices, Vec3(0, 0, 0));
    ContactGeometry::TriangleMesh mesh(vertices, faceIndices);
    
    // Test some points.
    
    Random::Gaussian random(0, 1);
    for (int i = 0; i < 100; i++) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        bool inside;
        UnitVec3 normal;
        Vec3 nearest = mesh.findNearestPoint(pos, inside, normal);
        ASSERT(inside == (~pos*normal < 1/Sqrt3));
        assertEqual(~nearest*normal, 1/Sqrt3);
        for (int j = 0; j < 3; j++) {
            ASSERT(pos[j]*nearest[j] >= 0 || std::abs(nearest[j]) < 100*Eps);
            ASSERT(pos[j]*normal[j] >= 0);
        }
    }
}

void testBoundingSphere() {
    Random::Uniform random(0, 10);
    for (int i = 0; i < 100; i++) {
        // Create a mesh consisting of a random number of octohedra at random places.
        
        vector<Vec3> vertices;
        vector<int> faceIndices;
        int numOctohedra = random.getIntValue()+1;
        for (int i = 0; i < numOctohedra; i++)
            addOctohedron(vertices, faceIndices, 
            Vec3(random.getValue(), random.getValue(), random.getValue()));
        ContactGeometry::TriangleMesh mesh(vertices, faceIndices);

        // Verify that all points are inside the bounding sphere.
        
        Vec3 center;
        Real radius;
        mesh.getBoundingSphere(center, radius);
        for (int i = 0; i < mesh.getNumVertices(); i++) {
            Real dist = (center-mesh.getVertexPosition(i)).norm();
            ASSERT(dist <= radius);
        }
        
        // Make sure the bounding sphere is reasonably compact.
        
        Vec3 boxRadius = 0.5*mesh.getOBBTreeNode().getBounds().getSize();
        ASSERT(radius <= boxRadius.norm());
    }
}

int main() {
    try {
        testTriangleMesh();
        testIncorrectMeshes();
        testOBBTree();
        testRayIntersection();
        testSmoothMesh();
        testFindNearestPoint();
        testBoundingSphere();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
