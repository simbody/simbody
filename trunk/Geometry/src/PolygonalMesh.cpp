/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
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

#include "PolygonalMeshImpl.h"
#include <cassert>
#include <sstream>
#include <string>

using std::string;
using std::stringstream;

namespace SimTK {

PolygonalMeshImpl::PolygonalMeshImpl() {
    faceVertexStart.push_back(0);
}

PolygonalMesh::PolygonalMesh() : impl(new PolygonalMeshImpl()) {
}

PolygonalMesh::PolygonalMesh(const PolygonalMesh& copy) : impl(new PolygonalMeshImpl(*copy.impl)) {
}

PolygonalMesh& PolygonalMesh::operator=(const PolygonalMesh& copy) {
    delete impl;
    impl = new PolygonalMeshImpl(*copy.impl);
    return *this;
}

PolygonalMesh::~PolygonalMesh() {
    delete impl;
}

int PolygonalMesh::getNumFaces() const {
    return getImpl().faceVertexStart.size()-1;
}

int PolygonalMesh::getNumVertices() const {
    return getImpl().vertices.size();
}

const Vec3& PolygonalMesh::getVertexPosition(int vertex) const {
    assert(vertex >= 0 && vertex < getNumVertices());
    return getImpl().vertices[vertex];
}

int PolygonalMesh::getNumVerticesForFace(int face) const {
    assert(face >= 0 && face < getNumFaces());
    const Array_<int>& faceVertexStart = getImpl().faceVertexStart;
    return faceVertexStart[face+1]-faceVertexStart[face];
}

int PolygonalMesh::getFaceVertex(int face, int vertex) const {
    assert(face >= 0 && face < getNumFaces());
    assert(vertex >= 0 && vertex < getNumVerticesForFace(face));
    return getImpl().faceVertexIndex[getImpl().faceVertexStart[face]+vertex];
}

int PolygonalMesh::addVertex(const Vec3& position) {
    updImpl().vertices.push_back(position);
    return getImpl().vertices.size()-1;
}

int PolygonalMesh::addFace(const Array_<int>& vertices) {
    for (int i = 0; i < (int) vertices.size(); i++)
        updImpl().faceVertexIndex.push_back(vertices[i]);
    updImpl().faceVertexStart.push_back(getImpl().faceVertexIndex.size());
    return getImpl().faceVertexStart.size()-2;
}

void PolygonalMesh::scaleMesh(Real scale) {
    Array_<Vec3>& vertices = updImpl().vertices;
    for (int i = 0; i < (int) vertices.size(); i++)
        vertices[i] *= scale;
}

void PolygonalMesh::transformMesh(const Transform& transform) {
    Array_<Vec3>& vertices = updImpl().vertices;
    for (int i = 0; i < (int) vertices.size(); i++)
        vertices[i] = transform*vertices[i];
}

const PolygonalMeshImpl& PolygonalMesh::getImpl() const {
    return *impl;
}

PolygonalMeshImpl& PolygonalMesh::updImpl() {
    return *impl;
}

void PolygonalMesh::loadObjFile(std::istream& file) {
    string line;
    Array_<int> indices;
    int initialVertices = getNumVertices();
    while (!file.eof()) {
        getline(file, line);
        while (line.size() > 0 && line[line.size()-1] == '\\') {
            line[line.size()-1] = ' ';
            string continuation;
            getline(file, continuation);
            line += continuation;
        }
        stringstream s(line);
        string command;
        s >> command;
        if (command == "v") {
            // A vertex
            
            Real x, y, z;
            s >> x;
            s >> y;
            SimTK_ASSERT1_ALWAYS(s >> z, "Found invalid vertex description: %s", line.c_str());
            addVertex(Vec3(x, y, z));
        }
        else if (command == "f") {
            // A face
            
            indices.clear();
            int index;
            while (s >> index) {
                s.ignore(line.size(), ' ');
                if (index < 0)
                    index += getNumVertices()-initialVertices;
                else
                    index--;
                indices.push_back(index);
            }
            addFace(indices);
        }
    }
}

} // namespace SimTK
