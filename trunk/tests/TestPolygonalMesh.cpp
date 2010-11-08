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

#include "SimTKcommon.h"

#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

void testCreateMesh() {
    PolygonalMesh mesh;
    ASSERT(mesh.getNumFaces() == 0);
    ASSERT(mesh.getNumVertices() == 0);
    ASSERT(mesh.addVertex(Vec3(0)) == 0);
    ASSERT(mesh.addVertex(Vec3(0, 1, 0)) == 1);
    ASSERT(mesh.addVertex(Vec3(0, 0, 1)) == 2);
    Array_<int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(0);
    ASSERT(mesh.addFace(v) == 0);
    ASSERT(mesh.getNumFaces() == 1);
    ASSERT(mesh.getNumVertices() == 3);
    ASSERT(mesh.getFaceVertex(0, 0) == 1);
    ASSERT(mesh.getFaceVertex(0, 1) == 2);
    ASSERT(mesh.getFaceVertex(0, 2) == 0);
    ASSERT(mesh.getVertexPosition(0) == Vec3(0));
    ASSERT(mesh.getVertexPosition(1) == Vec3(0, 1, 0));
    ASSERT(mesh.getVertexPosition(2) == Vec3(0, 0, 1));

    // Make sure copy and assignment are shallow.
    PolygonalMesh mesh2(mesh); // shallow copy construction
    ASSERT(&mesh2.getImpl() == &mesh.getImpl());
    ASSERT(mesh.getImplHandleCount()==2);

    PolygonalMesh mesh3;
    mesh3 = mesh; // shallow assignment
    ASSERT(&mesh3.getImpl() == &mesh.getImpl());
    ASSERT(mesh.getImplHandleCount()==3);

    PolygonalMesh mesh4;
    mesh4.copyAssign(mesh); // deep copy
    ASSERT(mesh4.getNumVertices() == mesh.getNumVertices());
    ASSERT(&mesh4.getImpl() != &mesh.getImpl());
    ASSERT(mesh4.getImplHandleCount()==1);
    ASSERT(mesh.getImplHandleCount()==3);
}

void testLoadObjFile() {
    string file;
    file += "# This is a comment\n";
    file += "v -1.0 1.0 2.0\n";
    file += "v -2.0 2.0 3.0\n";
    file += "v -3.0 3.0 \\\n";
    file += "4.0\n";
    file += "v -4.0 4.0 5.0\n";
    file += "v -5.0 5.0 6.0\n";
    file += "f 1 2 3\n";
    file += "f 2// 3/4/5 4//2\n";
    file += "f -3 -2/3/4 -1\n";
    file += "f 1 3\\\n";
    file += "5 2\n";
    PolygonalMesh mesh;
    stringstream stream(file);
    mesh.loadObjFile(stream);
    ASSERT(mesh.getNumVertices() == 5);
    ASSERT(mesh.getNumFaces() == 4);
    for (int i = 0; i < mesh.getNumVertices(); i++) {
        const Vec3& pos = mesh.getVertexPosition(i);
        ASSERT(pos[0] == -(i+1));
        ASSERT(pos[1] == i+1);
        ASSERT(pos[2] == i+2);
    }
    for (int i = 0; i < 3; i++) {
        ASSERT(mesh.getNumVerticesForFace(i) == 3);
        ASSERT(mesh.getFaceVertex(i, 0) == i);
        ASSERT(mesh.getFaceVertex(i, 1) == i+1);
        ASSERT(mesh.getFaceVertex(i, 2) == i+2);
    }
    ASSERT(mesh.getNumVerticesForFace(3) == 4);
    ASSERT(mesh.getFaceVertex(3, 0) == 0);
    ASSERT(mesh.getFaceVertex(3, 1) == 2);
    ASSERT(mesh.getFaceVertex(3, 2) == 4);
    ASSERT(mesh.getFaceVertex(3, 3) == 1);
}

int main() {
    try {
        testCreateMesh();
        testLoadObjFile();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
