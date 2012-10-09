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
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <iostream>
#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::cout; using std::endl;

//==============================================================================
//                  CONTACT GEOMETRY :: TRIANGLE MESH
//==============================================================================

ContactGeometry::TriangleMesh::TriangleMesh
   (const ArrayViewConst_<Vec3>& vertices, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(vertices, faceIndices, smooth)) {}

ContactGeometry::TriangleMesh::TriangleMesh
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(mesh, smooth)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::TriangleMesh::classTypeId() 
{   return ContactGeometry::TriangleMesh::Impl::classTypeId(); }


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

const UnitVec3& ContactGeometry::TriangleMesh::getFaceNormal(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].normal;
}

Real ContactGeometry::TriangleMesh::getFaceArea(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].area;
}

void ContactGeometry::TriangleMesh::
findVertexEdges(int vertex, Array_<int>& edges) const {
    // Begin at an arbitrary edge which intersects the vertex.
    
    int firstEdge = getImpl().vertices[vertex].firstEdge;
    int previousEdge = firstEdge;
    int previousFace = getImpl().edges[firstEdge].faces[0];
    
    // Walk around the vertex, using each edge to find the next face and each 
    // face to find the next edge.
    
    do {
        edges.push_back(previousEdge);
        const ContactGeometry::TriangleMesh::Impl::Edge& 
            edge = getImpl().edges[previousEdge];
        int nextFace = (edge.faces[0] == previousFace ? edge.faces[1] 
                                                      : edge.faces[0]);
        const ContactGeometry::TriangleMesh::Impl::Face& 
            face = getImpl().faces[nextFace];
        int nextEdge;
        if (    face.edges[0] != previousEdge
            && (face.vertices[0] == vertex || face.vertices[1] == vertex))
            nextEdge = face.edges[0];
        else if (   face.edges[1] != previousEdge 
                 && (face.vertices[1] == vertex || face.vertices[2] == vertex))
            nextEdge = face.edges[1];
        else
            nextEdge = face.edges[2];
        previousEdge = nextEdge;
        previousFace = nextFace;
    } while (previousEdge != firstEdge);
}

Vec3 ContactGeometry::TriangleMesh::findPoint(int face, const Vec2& uv) const {
    return getImpl().findPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findCentroid(int face) const {
    return getImpl().findCentroid(face);
}

UnitVec3 ContactGeometry::TriangleMesh::
findNormalAtPoint(int face, const Vec2& uv) const {
    return getImpl().findNormalAtPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, int& face, Vec2& uv) const {
    return getImpl().findNearestPoint(position, inside, face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    return getImpl().findNearestPointToFace(position, face, uv);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, 
    UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, 
    Vec2& uv) const {
    return getImpl().intersectsRay(origin, direction, distance, face, uv);
}

ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::getOBBTreeNode() const {
    return OBBTreeNode(getImpl().obb);
}

PolygonalMesh ContactGeometry::TriangleMesh::createPolygonalMesh() const {
    PolygonalMesh mesh;
    getImpl().createPolygonalMesh(mesh);
    return mesh;
}

const ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::getImpl() const {
    assert(impl);
    return static_cast<const TriangleMesh::Impl&>(*impl);
}

ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::updImpl() {
    assert(impl);
    return static_cast<TriangleMesh::Impl&>(*impl);
}



//==============================================================================
//                CONTACT GEOMETRY :: TRIANGLE MESH :: IMPL
//==============================================================================

DecorativeGeometry ContactGeometry::TriangleMesh::Impl::createDecorativeGeometry() const {
    PolygonalMesh mesh;
    createPolygonalMesh(mesh);
    return DecorativeMesh(mesh);
}

Vec3 ContactGeometry::TriangleMesh::Impl::findPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    return             uv[0] * vertices[f.vertices[0]].pos
           +           uv[1] * vertices[f.vertices[1]].pos
           +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].pos;
}

// same as findPoint(face, (1/3,1/3)) but faster
Vec3 ContactGeometry::TriangleMesh::Impl::findCentroid(int face) const {
    const Face& f = faces[face];
    return (  vertices[f.vertices[0]].pos
            + vertices[f.vertices[1]].pos
            + vertices[f.vertices[2]].pos) / 3;
}

UnitVec3 ContactGeometry::TriangleMesh::Impl::findNormalAtPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    if (smooth)
        return UnitVec3(            uv[0] * vertices[f.vertices[0]].normal
                        +           uv[1] * vertices[f.vertices[1]].normal
                        +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].normal);
    return f.normal;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    int face;
    Vec2 uv;
    Vec3 nearestPoint = findNearestPoint(position, inside, face, uv);
    normal = findNormalAtPoint(face, uv);
    return nearestPoint;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const 
{
    Real distance2;
    Vec3 nearestPoint = obb.findNearestPoint(*this, position, MostPositiveReal, distance2, face, uv);
    Vec3 delta = position-nearestPoint;
    inside = (~delta*faces[face].normal < 0);
    return nearestPoint;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              UnitVec3& normal) const {
    int face;
    Vec2 uv;
    if (!intersectsRay(origin, direction, distance, face, uv))
        return false;
    normal = findNormalAtPoint(face, uv);
    return true;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    Real boundsDistance;
    if (!obb.bounds.intersectsRay(origin, direction, boundsDistance))
        return false;
    return obb.intersectsRay(*this, origin, direction, distance, face, uv);
}

void ContactGeometry::TriangleMesh::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = boundingSphereCenter;
    radius = boundingSphereRadius;
}

void ContactGeometry::TriangleMesh::Impl::
createPolygonalMesh(PolygonalMesh& mesh) const {
    for (unsigned vx=0; vx < vertices.size(); ++vx)
        mesh.addVertex(vertices[vx].pos);
    for (unsigned fx=0; fx < faces.size(); ++fx) {
        const Face& face = faces[fx];
        const ArrayViewConst_<int> verts(face.vertices, face.vertices+3);
        mesh.addFace(verts);
    }
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const ArrayViewConst_<Vec3>& vertexPositions, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) {
    init(vertexPositions, faceIndices);
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) 
{   // Create the mesh, triangulating faces as necessary.
    Array_<Vec3>    vertexPositions;
    Array_<int>     faceIndices;
    for (int i = 0; i < mesh.getNumVertices(); i++)
        vertexPositions.push_back(mesh.getVertexPosition(i));
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        int numVert = mesh.getNumVerticesForFace(i);
        if (numVert < 3)
            continue; // Ignore it.
        if (numVert == 3) {
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
        }
        else if (numVert == 4) {
            // Split it into two triangles.
            
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 3));
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
        }
        else {
            // Add a vertex at the center, then split it into triangles.
            
            Vec3 center(0);
            for (int j = 0; j < numVert; j++)
                center += vertexPositions[mesh.getFaceVertex(i, j)];
            center /= numVert;
            vertexPositions.push_back(center);
            int newIndex = vertexPositions.size()-1;
            for (int j = 0; j < numVert-1; j++) {
                faceIndices.push_back(mesh.getFaceVertex(i, j));
                faceIndices.push_back(mesh.getFaceVertex(i, j+1));
                faceIndices.push_back(newIndex);
            }
        }
    }
    init(vertexPositions, faceIndices);
    
    // Make sure the mesh normals are oriented correctly.
    
    Vec3 origin(0);
    for (int i = 0; i < 3; i++)
        origin += vertices[faces[0].vertices[i]].pos;
    origin /= 3; // this is the face centroid

    const UnitVec3 direction = -faces[0].normal;
    // Calculate a ray origin that is guaranteed to be outside the
    // mesh. If the topology is right (face 0 normal points outward), we'll be
    // outside on the side containing face 0. If it is wrong, we'll be outside
    // on the opposite side of the mesh. Then we'll shoot a ray back along the
    // direction we came from (that is, towards the interior of the mesh from
    // outside). We'll hit *some* face. If the topology is right, the hit 
    // face's normal will be pointing back at us. If it is wrong, the face 
    // normal will also be pointing inwards, in roughly the same direction as 
    // the ray.
    origin -= max(obb.bounds.getSize())*direction;
    Real distance;
    int face;
    Vec2 uv;
    bool intersects = intersectsRay(origin, direction, distance, face, uv);
    assert(intersects);
    // Now dot the hit face normal with the ray direction; correct topology
    // will have them pointing in more-or-less opposite directions.
    if (dot(faces[face].normal, direction) > 0) {
        // We need to invert the mesh topology.
        
        for (int i = 0; i < (int) faces.size(); i++) {
            Face& f = faces[i];
            int temp = f.vertices[0];
            f.vertices[0] = f.vertices[1];
            f.vertices[1] = temp;
            temp = f.edges[1];
            f.edges[1] = f.edges[2];
            f.edges[2] = temp;
            f.normal *= -1;
        }
        for (int i = 0; i < (int) vertices.size(); i++)
            vertices[i].normal *= -1;
    }
}

void ContactGeometry::TriangleMesh::Impl::init
   (const Array_<Vec3>& vertexPositions, const Array_<int>& faceIndices) 
{   SimTK_APIARGCHECK_ALWAYS(faceIndices.size()%3 == 0, 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl", 
        "The number of indices must be a multiple of 3.");
    int numFaces = faceIndices.size()/3;
    
    // Create the vertices.
    
    for (int i = 0; i < (int) vertexPositions.size(); i++)
        vertices.push_back(Vertex(vertexPositions[i]));
    
    // Create the faces and build lists of all the edges.
    
    map<pair<int, int>, int> forwardEdges;
    map<pair<int, int>, int> backwardEdges;
    for (int i = 0; i < numFaces; i++) {
        int start = i*3;
        int v1 = faceIndices[start], v2 = faceIndices[start+1], 
            v3 = faceIndices[start+2];
        SimTK_APIARGCHECK1_ALWAYS
           (   v1 >= 0 && v1 < (int) vertices.size() 
            && v2 >= 0 && v2 < (int) vertices.size() 
            && v3 >= 0 && v3 < (int) vertices.size(),
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d contains a vertex with an illegal index.", i);
        Vec3 cross =   (vertexPositions[v2]-vertexPositions[v1])
                     % (vertexPositions[v3]-vertexPositions[v1]);
        Real norm = cross.norm();
        cross *= 1.0/norm;
        SimTK_APIARGCHECK1_ALWAYS(norm > 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d is degenerate.", i);
        faces.push_back(Face(v1, v2, v3, cross, 0.5*norm));
        int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};
        for (int j = 0; j < 3; j++) {
            SimTK_APIARGCHECK1_ALWAYS(edges[j][0] != edges[j][1], 
                "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                "Vertices %d appears twice in a single face.", edges[j][0]);
            if (edges[j][0] < edges[j][1]) {
                SimTK_APIARGCHECK2_ALWAYS
                   (forwardEdges.find(pair<int, int>(edges[j][0], edges[j][1])) 
                    == forwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][0], edges[j][1]);
                forwardEdges[pair<int, int>(edges[j][0], edges[j][1])] = i;
            }
            else {
                SimTK_APIARGCHECK2_ALWAYS
                   (backwardEdges.find(pair<int, int>(edges[j][1], edges[j][0]))
                    == backwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][1], edges[j][0]);
                backwardEdges[pair<int, int>(edges[j][1], edges[j][0])] = i;
            }
        }
    }
    
    // Create the edges.
    
    SimTK_APIARGCHECK_ALWAYS(forwardEdges.size() == backwardEdges.size(), 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
        "Each edge must be shared by exactly two faces.");
    for (map<pair<int, int>, int>::iterator iter = forwardEdges.begin(); 
         iter != forwardEdges.end(); ++iter) {
        int vert1 = iter->first.first;
        int vert2 = iter->first.second;
        int face1 = iter->second;
        map<pair<int, int>, int>::iterator iter2 = 
            backwardEdges.find(pair<int, int>(vert1, vert2));
        SimTK_APIARGCHECK_ALWAYS(iter2 != backwardEdges.end(), 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Each edge must be shared by exactly two faces.");
        int face2 = iter2->second;
        edges.push_back(Edge(vert1, vert2, face1, face2));
    }
    
    // Record the edges for each face.
    
    for (int i = 0; i < (int) edges.size(); i++) {
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
                SimTK_ASSERT_ALWAYS(false, 
                    "Face and edge vertices are inconsistent.");
        }
    }
    
    // Record a single edge for each vertex.
    
    for (int i = 0; i < (int) edges.size(); i++) {
        vertices[edges[i].vertices[0]].firstEdge = i;
        vertices[edges[i].vertices[1]].firstEdge = i;
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        SimTK_APIARGCHECK1_ALWAYS(vertices[i].firstEdge >= 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Vertex %d is not part of any face.", i);
    
    // Calculate a normal for each vertex.
    
    Vector_<Vec3> vertNorm(vertices.size(), Vec3(0));
    for (int i = 0; i < (int) faces.size(); i++) {
        const Face& f = faces[i];
        UnitVec3 edgeDir[3];
        for (int j = 0; j < 3; j++) {
            edgeDir[j] = UnitVec3(  vertices[f.vertices[(j+1)%3]].pos
                                  - vertices[f.vertices[j]].pos);
        }
        for (int j = 0; j < 3; j++) {
            Real angle = std::acos(~edgeDir[j]*edgeDir[(j+2)%3]);
            vertNorm[f.vertices[j]] += f.normal*angle;
        }
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        vertices[i].normal = UnitVec3(vertNorm[i]);
    
    // Create the OBBTree.
    
    Array_<int> allFaces(faces.size());
    for (int i = 0; i < (int) allFaces.size(); i++)
        allFaces[i] = i;
    createObbTree(obb, allFaces);
    
    // Find the bounding sphere.
    Array_<const Vec3*> points(vertices.size());
    for (int i = 0; i < (int) vertices.size(); i++)
        points[i] = &vertices[i].pos;
    const Geo::Sphere bnd = Geo::Point::calcBoundingSphereIndirect(points);
    boundingSphereCenter = bnd.getCenter();
    boundingSphereRadius = bnd.getRadius();
}

void ContactGeometry::TriangleMesh::Impl::createObbTree
   (OBBTreeNodeImpl& node, const Array_<int>& faceIndices) 
{   // Find all vertices in the node and build the OrientedBoundingBox.
    node.numTriangles = faceIndices.size();
    set<int> vertexIndices;
    for (int i = 0; i < (int) faceIndices.size(); i++) 
        for (int j = 0; j < 3; j++)
            vertexIndices.insert(faces[faceIndices[i]].vertices[j]);
    Vector_<Vec3> points((int)vertexIndices.size());
    int index = 0;
    for (set<int>::iterator iter = vertexIndices.begin(); 
                            iter != vertexIndices.end(); ++iter)
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
            Array_<int> child1Indices, child2Indices;
            splitObbAxis(faceIndices, child1Indices, child2Indices, 
                         axisOrder[i]);
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
    
    node.triangles.insert(node.triangles.begin(), faceIndices.begin(), 
                          faceIndices.end());
}

void ContactGeometry::TriangleMesh::Impl::splitObbAxis
   (const Array_<int>& parentIndices, Array_<int>& child1Indices, 
    Array_<int>& child2Indices, int axis) 
{   // For each face, find its minimum and maximum extent along the axis.
    Vector minExtent(parentIndices.size());
    Vector maxExtent(parentIndices.size());
    for (int i = 0; i < (int) parentIndices.size(); i++) {
        int* vertexIndices = faces[parentIndices[i]].vertices;
        Real minVal = vertices[vertexIndices[0]].pos[axis];
        Real maxVal = vertices[vertexIndices[0]].pos[axis];
        minVal = std::min(minVal, vertices[vertexIndices[1]].pos[axis]);
        maxVal = std::max(maxVal, vertices[vertexIndices[1]].pos[axis]);
        minExtent[i] = std::min(minVal, vertices[vertexIndices[2]].pos[axis]);
        maxExtent[i] = std::max(maxVal, vertices[vertexIndices[2]].pos[axis]);
    }
    
    // Select a split point that tries to put as many faces as possible 
    // entirely on one side or the other.
    
    Real split = 0.5*(median(minExtent)+median(maxExtent));
    
    // Choose a side for each face.
    
    for (int i = 0; i < (int) parentIndices.size(); i++) {
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

Vec3 ContactGeometry::TriangleMesh::Impl::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    // Calculate the distance between a point in space and a face of the mesh.
    // This algorithm is based on a description by David Eberly found at 
    // http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf.
    
    const ContactGeometry::TriangleMesh::Impl::Face& fc = faces[face];
    const Vec3& vert1 = vertices[fc.vertices[0]].pos;
    const Vec3& vert2 = vertices[fc.vertices[1]].pos;
    const Vec3& vert3 = vertices[fc.vertices[2]].pos;
    const Vec3 e0 = vert2-vert1;
    const Vec3 e1 = vert3-vert1;
    const Vec3 delta = vert1-position;
    const Real a = e0.normSqr();
    const Real b = ~e0*e1;
    const Real c = e1.normSqr();
    const Real d = ~e0*delta;
    const Real e = ~e1*delta;
    const Real f = delta.normSqr();
    const Real det = a*c-b*b;
    Real s = b*e-c*d;
    Real t = b*d-a*e;
    if (s+t <= det) {
        if (s < 0) {
            if (t < 0) {
                // Region 4

                if (d < 0) {
                    s = (-d >= a ? 1 : -d/a);
                    t = 0;
                }
                else {
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            }
            else {
                // Region 3

                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 5

            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
            t = 0;
        }
        else {
            // Region 0

            const Real invDet = 1.0/det;
            s *= invDet;
            t *= invDet;
        }
    }
    else {
        if (s < 0) {
            // Region 2

            Real temp0 = b+d;
            Real temp1 = c+e;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1-s;
            }
            else {
                s = 0;
                t = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 6

            Real temp0 = b+e;
            Real temp1 = a+d;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                t = (numer >= denom ? 1 : numer/denom);
                s = 1-t;
            }
            else {
                s = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -d/a));
                t = 0;
            }
        }
        else {
            // Region 1

            const Real numer = c+e-b-d;
            if (numer <= 0)
                s = 0;
            else {
                const Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
            }
            t = 1-s;
        }
    }
    uv = Vec2(1-s-t, s);
    return vert1 + s*e0 + t*e1;
}


//==============================================================================
//                            OBB TREE NODE IMPL
//==============================================================================

OBBTreeNodeImpl::OBBTreeNodeImpl(const OBBTreeNodeImpl& copy) 
:   bounds(copy.bounds), triangles(copy.triangles), 
    numTriangles(copy.numTriangles) {
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

Vec3 OBBTreeNodeImpl::findNearestPoint
   (const ContactGeometry::TriangleMesh::Impl& mesh, 
    const Vec3& position, Real cutoff2, 
    Real& distance2, int& face, Vec2& uv) const 
{
    Real tol = 100*Eps;
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance2 = MostPositiveReal, 
             child2distance2 = MostPositiveReal;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        Vec3 child1point, child2point;
        Real child1BoundsDist2 = 
            (child1->bounds.findNearestPoint(position)-position).normSqr();
        Real child2BoundsDist2 = 
            (child2->bounds.findNearestPoint(position)-position).normSqr();
        if (child1BoundsDist2 < child2BoundsDist2) {
            if (child1BoundsDist2 < cutoff2) {
                child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
                if (child2BoundsDist2 < child1distance2 && child2BoundsDist2 < cutoff2)
                    child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
            }
        }
        else {
            if (child2BoundsDist2 < cutoff2) {
                child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
                if (child1BoundsDist2 < child2distance2 && child1BoundsDist2 < cutoff2)
                    child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
            }
        }
        if (   child1distance2 <= child2distance2*(1+tol) 
            && child2distance2 <= child1distance2*(1+tol)) {
            // Decide based on angle which one to use.
            
            if (  std::abs(~(child1point-position)*mesh.faces[child1face].normal) 
                > std::abs(~(child2point-position)*mesh.faces[child2face].normal))
                child2distance2 = MostPositiveReal;
            else
                child1distance2 = MostPositiveReal;
        }
        if (child1distance2 < child2distance2) {
            distance2 = child1distance2;
            face = child1face;
            uv = child1uv;
            return child1point;
        }
        else {
            distance2 = child2distance2;
            face = child2face;
            uv = child2uv;
            return child2point;
        }
    }    
    // This is a leaf node, so check each triangle for its distance to the point.
    
    distance2 = MostPositiveReal;
    Vec3 nearestPoint;
    for (int i = 0; i < (int) triangles.size(); i++) {
        Vec2 triangleUV;
        Vec3 p = mesh.findNearestPointToFace(position, triangles[i], triangleUV);
        Vec3 offset = p-position;
        // TODO: volatile to work around compiler bug
        volatile Real d2 = offset.normSqr(); 
        if (d2 < distance2 || (d2 < distance2*(1+tol) && std::abs(~offset*mesh.faces[triangles[i]].normal) > std::abs(~offset*mesh.faces[face].normal))) {
            nearestPoint = p;
            distance2 = d2;
            face = triangles[i];
            uv = triangleUV;
        }
    }
    return nearestPoint;
}

bool OBBTreeNodeImpl::
intersectsRay(const ContactGeometry::TriangleMesh::Impl& mesh,
              const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance, child2distance;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        bool child1intersects = child1->bounds.intersectsRay(origin, direction, child1distance);
        bool child2intersects = child2->bounds.intersectsRay(origin, direction, child2distance);
        if (child1intersects) {
            if (child2intersects) {
                // The ray intersects both child nodes.  First check the closer one.
                
                if (child1distance < child2distance) {
                    child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                    if (!child1intersects || child2distance < child1distance)
                        child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                }
                else {
                    child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                    if (!child2intersects || child1distance < child2distance)
                        child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                }
            }
            else
                child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
        }
        else if (child2intersects)
            child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
        
        // If either one had an intersection, return the closer one.
        
        if (child1intersects && (!child2intersects || child1distance < child2distance)) {
            distance = child1distance;
            face = child1face;
            uv = child1uv;
            return true;
        }
        if (child2intersects) {
            distance = child2distance;
            face = child2face;
            uv = child2uv;
            return true;
        }
        return false;
    }
    
    // This is a leaf node, so check each triangle for an intersection with the 
    // ray.
    
    bool foundIntersection = false;
    for (int i = 0; i < (int) triangles.size(); i++) {
        const UnitVec3& faceNormal = mesh.faces[triangles[i]].normal;
        double vd = ~faceNormal*direction;
        if (vd == 0.0)
            continue; // The ray is parallel to the plane.
        const Vec3& vert1 = mesh.vertices[mesh.faces[triangles[i]].vertices[0]].pos;
        double v0 = ~faceNormal*(vert1-origin);
        double t = v0/vd;
        if (t < 0.0)
            continue; // Ray points away from plane of triangle.
        if (foundIntersection && t >= distance)
            continue; // We already have a closer intersection.

        // Determine whether the intersection point is inside the triangle by projecting onto
        // a plane and computing the barycentric coordinates.

        Vec3 ri = origin+direction*t;
        const Vec3& vert2 = mesh.vertices[mesh.faces[triangles[i]].vertices[1]].pos;
        const Vec3& vert3 = mesh.vertices[mesh.faces[triangles[i]].vertices[2]].pos;
        int axis1, axis2;
        if (std::abs(faceNormal[1]) > std::abs(faceNormal[0])) {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[1])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 0;
                axis2 = 2;
            }
        }
        else {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[0])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 1;
                axis2 = 2;
            }
        }
        Vec2 pos(ri[axis1]-vert1[axis1], ri[axis2]-vert1[axis2]);
        Vec2 edge1(vert1[axis1]-vert2[axis1], vert1[axis2]-vert2[axis2]);
        Vec2 edge2(vert1[axis1]-vert3[axis1], vert1[axis2]-vert3[axis2]);
        double denom = 1.0/(edge1%edge2);
        edge2 *= denom;
        double v = edge2%pos;
        if (v < 0.0 || v > 1.0)
            continue;
        edge1 *= denom;
        double w = pos%edge1;
        if (w < 0.0 || w > 1.0)
            continue;
        double u = 1.0-v-w;
        if (u < 0.0 || u > 1.0)
            continue;
        
        // It intersects.
        
        distance = t;
        face = triangles[i];
        uv = Vec2(u, v);
        foundIntersection = true;
    }
    return foundIntersection;
}




//==============================================================================
//            CONTACT GEOMETRY :: TRIANGLE MESH :: OBB TREE NODE
//==============================================================================

ContactGeometry::TriangleMesh::OBBTreeNode::
OBBTreeNode(const OBBTreeNodeImpl& impl) : impl(&impl) {}

const OrientedBoundingBox& 
ContactGeometry::TriangleMesh::OBBTreeNode::getBounds() const {
    return impl->bounds;
}

bool ContactGeometry::TriangleMesh::OBBTreeNode::isLeafNode() const {
    return (impl->child1 == NULL);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getFirstChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child1, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child1);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getSecondChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child2, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child2);
}

const Array_<int>& ContactGeometry::TriangleMesh::OBBTreeNode::
getTriangles() const {
    SimTK_ASSERT_ALWAYS(impl->child2 == NULL, 
        "Called getTriangles() on a non-leaf node");
    return impl->triangles;
}

int ContactGeometry::TriangleMesh::OBBTreeNode::getNumTriangles() const {
    return impl->numTriangles;
}

