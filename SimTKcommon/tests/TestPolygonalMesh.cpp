/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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
