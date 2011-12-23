/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

/* Tests for low-level geometric primitives and algorithms. */

#include "SimTKmath.h"
#include <vector>
#include <exception>

using namespace SimTK;
using namespace std;

static void addOctohedron(vector<Vec3>& vertices, vector<int>& faceIndices, 
                          Vec3 offset) {
    int start = (int)vertices.size();
    vertices.push_back(Vec3(0, 1, 0)+offset);
    vertices.push_back(Vec3(1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, 1)+offset);
    vertices.push_back(Vec3(-1, 0, 0)+offset);
    vertices.push_back(Vec3(0, 0, -1)+offset);
    vertices.push_back(Vec3(0, -1, 0)+offset);
    int faces[8][3] = {{0, 2, 1}, {0, 3, 2}, {0, 4, 3}, {0, 1, 4}, 
                       {5, 1, 2}, {5, 2, 3}, {5, 3, 4}, {5, 4, 1}};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]+start);
}


void testTriMeshBoundingSphere() {
    Random::Uniform random(0, 10);
    for (int i = 0; i < 100; i++) {
        // Create a mesh consisting of a random number of octohedra at random 
        // places.
        
        vector<Vec3> vertices;
        vector<int> faceIndices;
        int numOctohedra = random.getIntValue()+1;
        for (int j = 0; j < numOctohedra; j++)
            addOctohedron(vertices, faceIndices, 
            Vec3(random.getValue(), random.getValue(), random.getValue()));
        ContactGeometry::TriangleMesh mesh(vertices, faceIndices);


        //if (i==97) {
        //    Geo::Sphere bs = Geo::Sphere::calcBoundingSphere
        //        (Array_<Vec3>(vertices));
        //    for (int vx=0; vx<(int)vertices.size(); ++vx)
        //        if (bs.isPointOutside(vertices[vx]))
        //            printf("point %d outside by %g\n", vx,
        //                (vertices[vx]-bs.getCenter()).norm()-bs.getRadius());
        //    Array_<int> which;
        //    const Vec3& p21 = vertices[21];
        //    Geo::Sphere bs3 = Geo::Sphere::calcBoundingSphere
        //        (vertices[22],vertices[1],vertices[26], which);
        //    cout << "(22,1,26)->" << which << "\n";
        //    cout << "p21 in by " << (p21-bs3.getCenter()).norm()
        //        - bs3.getRadius() << "\n";
        //    Geo::Sphere bs4 = Geo::Sphere::calcBoundingSphere
        //        (vertices[30],vertices[22],vertices[26],vertices[1], which);
        //    cout << "(30,22,26,1)->" << which << "\n";
        //    cout << "p21 in by " << (p21-bs4.getCenter()).norm()
        //        - bs4.getRadius() << "\n";
        //}

        // Verify that all points are inside the bounding sphere.
        
        Vec3 center;
        Real radius;
        mesh.getBoundingSphere(center, radius);
        for (int j = 0; j < mesh.getNumVertices(); j++) {
            Real dist = (center-mesh.getVertexPosition(j)).norm();
            SimTK_TEST(dist <= radius);
        }

        
        // Make sure the bounding sphere is reasonably compact.
        
        Vec3 boxRadius = 0.5*mesh.getOBBTreeNode().getBounds().getSize();
        SimTK_TEST(radius <= boxRadius.norm());
    }
}

int main() {
    SimTK_START_TEST("TestGeo");
        SimTK_SUBTEST(testTriMeshBoundingSphere);
    SimTK_END_TEST();
}
