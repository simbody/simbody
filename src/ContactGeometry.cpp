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

#include "simbody/internal/ContactGeometryImpl.h"
#include <pthread.h>
#include <map>
#include <set>
#include <vector>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

ContactGeometry::ContactGeometry(ContactGeometryImpl* impl) : impl(impl) {
    assert(impl);
    impl->setMyHandle(*this);
}

ContactGeometry::~ContactGeometry() {
    if (isOwnerHandle())
        delete impl;
    impl = 0;
}

bool ContactGeometry::isOwnerHandle() const {
    return (impl == 0 || impl->getMyHandle() == this);
}

bool ContactGeometry::isEmptyHandle() const {
    return (impl == 0);
}

ContactGeometry::ContactGeometry(const ContactGeometry& src) : impl(0) {
    if (src.impl) {
        impl = src.impl->clone();
        impl->setMyHandle(*this);
    }
}

ContactGeometry& ContactGeometry::operator=(const ContactGeometry& src) {
    if (&src != this) {
        if (isOwnerHandle())
            delete impl;
        impl = 0;
        if (src.impl) {
            impl = src.impl->clone();
            impl->setMyHandle(*this);
        }
    }
    return *this;
}

const string& ContactGeometry::getType() const {
    return impl->getType();
}

int ContactGeometry::getTypeIndex() const {
    return impl->getTypeIndex();
}

ContactGeometryImpl::ContactGeometryImpl(const string& type) : myHandle(0), type(type), typeIndex(getIndexForType(type)) {
}

int ContactGeometryImpl::getIndexForType(std::string type) {
    static map<string, int> typeMap;
    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);
    map<string, int>::iterator index = typeMap.find(type);
    int indexForType;
    if (index == typeMap.end()) {
        indexForType = typeMap.size();
        typeMap[type] = indexForType;
    }
    else
        indexForType = index->second;
    pthread_mutex_unlock(&lock);
    return indexForType;
}

ContactGeometry::HalfSpace::HalfSpace() : ContactGeometry(new HalfSpaceImpl()) {
}

ContactGeometry::Sphere::Sphere(Real radius) : ContactGeometry(new SphereImpl(radius)) {
}

Real ContactGeometry::Sphere::getRadius() const {
    return getImpl().getRadius();
}

void ContactGeometry::Sphere::setRadius(Real radius) {
    updImpl().setRadius(radius);
}

const ContactGeometry::SphereImpl& ContactGeometry::Sphere::getImpl() const {
    assert(impl);
    return static_cast<const SphereImpl&>(*impl);
}

ContactGeometry::SphereImpl& ContactGeometry::Sphere::updImpl() {
    assert(impl);
    return static_cast<SphereImpl&>(*impl);
}

ContactGeometry::TriangleMesh::TriangleMesh(const vector<Vec3>& vertices, const vector<int>& faceIndices) : ContactGeometry(new TriangleMeshImpl(vertices, faceIndices)) {
}

int ContactGeometry::TriangleMesh::getNumEdges() const {
    return getImpl().edges.size();
}

int ContactGeometry::TriangleMesh::getNumFaces() const {
    return getImpl().faces.size();
}

int ContactGeometry::TriangleMesh::getNumVertices() const {
    return getImpl().vertices.size();
}

const Vec3& ContactGeometry::TriangleMesh::getVertexPosition(int index) const {
    assert(index >= 0 && index < getNumVertices());
    return getImpl().vertices[index].pos;
}

int ContactGeometry::TriangleMesh::getFaceEdge(int face, int edge) const {
    assert(face >= 0 && face < getNumFaces());
    assert(edge >= 0 && edge < 3);
    return getImpl().faces[face].edges[edge];
}

int ContactGeometry::TriangleMesh::getFaceVertex(int face, int vertex) const {
    assert(face >= 0 && face < getNumFaces());
    assert(vertex >= 0 && vertex < 3);
    return getImpl().faces[face].vertices[vertex];
}

int ContactGeometry::TriangleMesh::getEdgeFace(int edge, int face) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(face >= 0 && face < 2);
    return getImpl().edges[edge].faces[face];
}

int ContactGeometry::TriangleMesh::getEdgeVertex(int edge, int vertex) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(vertex >= 0 && vertex < 2);
    return getImpl().edges[edge].vertices[vertex];
}

void ContactGeometry::TriangleMesh::findVertexEdges(int vertex, std::vector<int>& edges) const {
    // Begin at an arbitrary edge which intersects the vertex.
    
    int firstEdge = getImpl().vertices[vertex].firstEdge;
    int previousEdge = firstEdge;
    int previousFace = getImpl().edges[firstEdge].faces[0];
    
    // Walk around the vertex, using each edge to find the next face and each face to find the next edge.
    
    do {
        edges.push_back(previousEdge);
        const ContactGeometry::TriangleMeshImpl::Edge& edge = getImpl().edges[previousEdge];
        int nextFace = (edge.faces[0] == previousFace ? edge.faces[1] : edge.faces[0]);
        const ContactGeometry::TriangleMeshImpl::Face& face = getImpl().faces[nextFace];
        int nextEdge;
        if (face.edges[0] != previousEdge && (face.vertices[0] == vertex || face.vertices[1] == vertex))
            nextEdge = face.edges[0];
        else if (face.edges[1] != previousEdge && (face.vertices[1] == vertex || face.vertices[2] == vertex))
            nextEdge = face.edges[1];
        else
            nextEdge = face.edges[2];
        previousEdge = nextEdge;
        previousFace = nextFace;
    } while (previousEdge != firstEdge);
}

ContactGeometry::TriangleMesh::OBBTreeNode ContactGeometry::TriangleMesh::getOBBTreeNode() const {
    return OBBTreeNode(getImpl().obb);
}

const ContactGeometry::TriangleMeshImpl& ContactGeometry::TriangleMesh::getImpl() const {
    assert(impl);
    return static_cast<const TriangleMeshImpl&>(*impl);
}

ContactGeometry::TriangleMeshImpl& ContactGeometry::TriangleMesh::updImpl() {
    assert(impl);
    return static_cast<TriangleMeshImpl&>(*impl);
}

ContactGeometry::TriangleMeshImpl::TriangleMeshImpl(const std::vector<Vec3>& vertexPositions, const vector<int>& faceIndices) : ContactGeometryImpl(Type()) {
    SimTK_APIARGCHECK_ALWAYS(faceIndices.size()%3 == 0, "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl", "The number of indices must be a multiple of 3.");
    int numFaces = faceIndices.size()/3;
    
    // Create the vertices.
    
    for (int i = 0; i < vertexPositions.size(); i++)
        vertices.push_back(Vertex(vertexPositions[i]));
    
    // Create the faces and build lists of all the edges.
    
    map<pair<int, int>, int> forwardEdges;
    map<pair<int, int>, int> backwardEdges;
    for (int i = 0; i < numFaces; i++) {
        int start = i*3;
        int v1 = faceIndices[start], v2 = faceIndices[start+1], v3 = faceIndices[start+2];
        SimTK_APIARGCHECK1_ALWAYS(v1 >= 0 && v1 < vertices.size() && v2 >= 0 && v2 < vertices.size() && v3 >= 0 && v3 < vertices.size(),
                "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                "Face %d contains a vertex with an illegal index.", i);
        faces.push_back(Face(v1, v2, v3));
        int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};
        for (int j = 0; j < 3; j++) {
            SimTK_APIARGCHECK1_ALWAYS(edges[j][0] != edges[j][1], "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                    "Vertices %d appears twice in a single face.", edges[j][0]);
            if (edges[j][0] < edges[j][1]) {
                SimTK_APIARGCHECK2_ALWAYS(forwardEdges.find(pair<int, int>(edges[j][0], edges[j][1])) == forwardEdges.end(),
                        "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                        "Multiple faces have an edge between vertices %d and %d in the same order.", edges[j][0], edges[j][1]);
                forwardEdges[pair<int, int>(edges[j][0], edges[j][1])] = i;
            }
            else {
                SimTK_APIARGCHECK2_ALWAYS(backwardEdges.find(pair<int, int>(edges[j][1], edges[j][0])) == backwardEdges.end(),
                        "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                        "Multiple faces have an edge between vertices %d and %d in the same order.", edges[j][1], edges[j][0]);
                backwardEdges[pair<int, int>(edges[j][1], edges[j][0])] = i;
            }
        }
    }
    
    // Create the edges.
    
    SimTK_APIARGCHECK_ALWAYS(forwardEdges.size() == backwardEdges.size(), "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
            "Each edge must be shared by exactly two faces.");
    for (map<pair<int, int>, int>::iterator iter = forwardEdges.begin(); iter != forwardEdges.end(); ++iter) {
        int vert1 = iter->first.first;
        int vert2 = iter->first.second;
        int face1 = iter->second;
        map<pair<int, int>, int>::iterator iter2 = backwardEdges.find(pair<int, int>(vert1, vert2));
        SimTK_APIARGCHECK_ALWAYS(iter2 != backwardEdges.end(), "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                "Each edge must be shared by exactly two faces.");
        int face2 = iter2->second;
        edges.push_back(Edge(vert1, vert2, face1, face2));
    }
    
    // Record the edges for each face.
    
    for (int i = 0; i < edges.size(); i++) {
        Edge& edge = edges[i];
        int f[2] = {edge.faces[0], edge.faces[1]};
        for (int j = 0; j < 2; j++) {
            Face& face = faces[f[j]];
            if ((edge.vertices[0] == face.vertices[0] || edge.vertices[0] == face.vertices[1]) &&
                    (edge.vertices[1] == face.vertices[0] || edge.vertices[1] == face.vertices[1]))
                face.edges[0] = i;
            else if ((edge.vertices[0] == face.vertices[1] || edge.vertices[0] == face.vertices[2]) &&
                    (edge.vertices[1] == face.vertices[1] || edge.vertices[1] == face.vertices[2]))
                face.edges[1] = i;
            else if ((edge.vertices[0] == face.vertices[2] || edge.vertices[0] == face.vertices[0]) &&
                    (edge.vertices[1] == face.vertices[2] || edge.vertices[1] == face.vertices[0]))
                face.edges[2] = i;
            else
                SimTK_ASSERT_ALWAYS(false, "Face and edge vertices are inconsistent.");
        }
    }
    
    // Record a single edge for each vertex.
    
    for (int i = 0; i < edges.size(); i++) {
        vertices[edges[i].vertices[0]].firstEdge = i;
        vertices[edges[i].vertices[1]].firstEdge = i;
    }
    for (int i = 0; i < vertices.size(); i++)
        SimTK_APIARGCHECK1_ALWAYS(vertices[i].firstEdge >= 0, "ContactGeometry::TriangleMeshImpl", "TriangleMeshImpl",
                "Vertex %d is not part of any face.", i);
    
    // Create the OBBTree.
    
    vector<int> allFaces(faces.size());
    for (int i = 0; i < allFaces.size(); i++)
        allFaces[i] = i;
    createObbTree(obb, allFaces);
}

void ContactGeometry::TriangleMeshImpl::createObbTree(OBBTreeNodeImpl& node, const vector<int>& faceIndices) {
    // Find all vertices in the node and build the OrientedBoundingBox.

    set<int> vertexIndices;
    for (int i = 0; i < faceIndices.size(); i++) 
        for (int j = 0; j < 3; j++)
            vertexIndices.insert(faces[faceIndices[i]].vertices[j]);
    Vector_<Vec3> points(vertexIndices.size());
    int index = 0;
    for (set<int>::iterator iter = vertexIndices.begin(); iter != vertexIndices.end(); ++iter)
        points[index++] = vertices[*iter].pos;
    node.bounds = OrientedBoundingBox(points);
    if (faceIndices.size() > 3) {

        // Order the axes by size.

        int axisOrder[3];
        const Vec3& size = node.bounds.getSize();
        if (size[0] > size[1]) {
            if (size[0] > size[2]) {
                axisOrder[0] = 0;
                if (size[1] > size[2]) {
                    axisOrder[1] = 1;
                    axisOrder[2] = 2;
                }
                else {
                    axisOrder[1] = 2;
                    axisOrder[2] = 1;
                }
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 0;
                axisOrder[2] = 1;
            }
        }
        else if (size[0] > size[2]) {
            axisOrder[0] = 1;
            axisOrder[1] = 0;
            axisOrder[2] = 2;
        }
        else {
            if (size[1] > size[2]) {
                axisOrder[0] = 1;
                axisOrder[1] = 2;
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 1;
            }
            axisOrder[2] = 0;
        }

        // Try splitting along each axis.

        for (int i = 0; i < 3; i++) {
            vector<int> child1Indices;
            vector<int> child2Indices;
            splitObbAxis(faceIndices, child1Indices, child2Indices, axisOrder[i]);
            if (child1Indices.size() > 0 && child2Indices.size() > 0) {
                // It was successfully split, so create the child nodes.

                node.child1 = new OBBTreeNodeImpl();
                node.child2 = new OBBTreeNodeImpl();
                createObbTree(*node.child1, child1Indices);
                createObbTree(*node.child2, child2Indices);
                return;
            }
        }
    }
    
    // This is a leaf node.
    
    node.triangles.insert(node.triangles.begin(), faceIndices.begin(), faceIndices.end());
}

void ContactGeometry::TriangleMeshImpl::splitObbAxis(const vector<int>& parentIndices, vector<int>& child1Indices, vector<int>& child2Indices, int axis) {
    // For each face, find its minimum and maximum extent along the axis.
    
    Vector minExtent(parentIndices.size());
    Vector maxExtent(parentIndices.size());
    for (int i = 0; i < parentIndices.size(); i++) {
        int* vertexIndices = faces[parentIndices[i]].vertices;
        Real minVal = vertices[vertexIndices[0]].pos[axis];
        Real maxVal = vertices[vertexIndices[0]].pos[axis];
        minVal = std::min(minVal, vertices[vertexIndices[1]].pos[axis]);
        maxVal = std::max(maxVal, vertices[vertexIndices[1]].pos[axis]);
        minExtent[i] = std::min(minVal, vertices[vertexIndices[2]].pos[axis]);
        maxExtent[i] = std::max(maxVal, vertices[vertexIndices[2]].pos[axis]);
    }
    
    // Select a split point that tries to put as many faces as possible entirely on one side or the other.
    
    Real split = 0.5*(median(minExtent)+median(maxExtent));
    
    // Choose a side for each face.
    
    for (int i = 0; i < parentIndices.size(); i++) {
        if (maxExtent[i] <= split)
            child1Indices.push_back(parentIndices[i]);
        else if (minExtent[i] >= split)
            child2Indices.push_back(parentIndices[i]);
        else if (0.5*(minExtent[i]+maxExtent[i]) <= split)
            child1Indices.push_back(parentIndices[i]);
        else
            child2Indices.push_back(parentIndices[i]);
    }
}

OBBTreeNodeImpl::OBBTreeNodeImpl(const OBBTreeNodeImpl& copy) : bounds(copy.bounds), triangles(copy.triangles) {
    if (copy.child1 == NULL) {
        child1 = NULL;
        child2 = NULL;
    }
    else {
        child1 = new OBBTreeNodeImpl(*copy.child1);
        child2 = new OBBTreeNodeImpl(*copy.child2);
    }
}

OBBTreeNodeImpl::~OBBTreeNodeImpl() {
    if (child1 != NULL)
        delete child1;
    if (child2 != NULL)
        delete child2;
}

ContactGeometry::TriangleMesh::OBBTreeNode::OBBTreeNode(const OBBTreeNodeImpl& impl) : impl(&impl) {
}

const OrientedBoundingBox& ContactGeometry::TriangleMesh::OBBTreeNode::getBounds() const {
    return impl->bounds;
}

bool ContactGeometry::TriangleMesh::OBBTreeNode::isLeafNode() const {
    return (impl->child1 == NULL);
}

const ContactGeometry::TriangleMesh::OBBTreeNode ContactGeometry::TriangleMesh::OBBTreeNode::getFirstChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child1, "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child1);
}

const ContactGeometry::TriangleMesh::OBBTreeNode ContactGeometry::TriangleMesh::OBBTreeNode::getSecondChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child2, "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child2);
}

const std::vector<int>& ContactGeometry::TriangleMesh::OBBTreeNode::getTriangles() const {
    SimTK_ASSERT_ALWAYS(impl->child2 == NULL, "Called getTriangles() on a non-leaf node");
    return impl->triangles;
}
