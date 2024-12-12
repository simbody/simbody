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
#include <fstream>
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

    ASSERT(!mesh.hasNormals())
    ASSERT(!mesh.hasTextureCoordinates())
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
    ASSERT(!mesh.hasNormals())
    ASSERT(!mesh.hasTextureCoordinates())
}

void testLoadObjFileWithNormalsTexture() {
    string file;
    file += "# This is a comment\n";
    file += "v 1.000000 1.000000 -1.000000\n";
    file += "v 1.000000 -1.000000 -1.00000\n";
    file += "v 1.000000 1.000000 1.000000\n";
    file += "v 1.000000 -1.000000 1.000000\n";
    file += "v -1.000000 1.000000 -1.00000\n";
    file += "v -1.000000 -1.000000 -1.000\n";
    file += "v -1.000000 1.000000 1.000000\n";
    file += "v -1.000000 -1.000000 1.000000\n";
    file += "vn -0.0000 1.0000 -0.0000\n";
    file += "vn -0.0000 -0.0000 1.0000\n";
    file += "vn -1.0000 -0.0000 -0.0000\n";
    file += "vn -0.0000 -1.0000 -0.0000\n";
    file += "vn 1.0000 -0.0000 -0.0000\n";
    file += "vn -0.0000 -0.0000 -1.0000\n";
    file += "vt 0.625000 0.500000\n";
    file += "vt 0.875000 0.500000\n";
    file += "vt 0.875000 0.750000\n";
    file += "vt 0.625000 0.750000\n";
    file += "vt 0.375000 0.750000\n";
    file += "vt 0.625000 1.000000\n";
    file += "vt 0.375000 1.000000\n";
    file += "vt 0.375000 0.000000\n";
    file += "vt 0.625000 0.000000\n";
    file += "vt 0.625000 0.250000\n";
    file += "vt 0.375000 0.250000\n";
    file += "vt 0.125000 0.500000\n";
    file += "vt 0.375000 0.500000\n";
    file += "vt 0.125000 0.750000\n";
    file += "s 0\n";
    file += "usemtl Material\n";
    file += "f 1/1/1 5/2/1 7/3/1 3/4/1\n";
    file += "f 4/5/2 3/4/2 7/6/2 8/7/2\n";
    file += "f 8/8/3 7/9/3 5/10/3 6/11/3\n";
    file += "f 6/12/4 2/13/4 4/5/4 8/14/4\n";
    file += "f 2/13/5 1/1/5 3/4/5 4/5/5\n";
    file += "f 6/11/6 5/10/6 1/1/6 2/13/6\n";

    PolygonalMesh mesh;
    stringstream stream(file);
    mesh.loadObjFile(stream);
    ASSERT(mesh.getNumVertices() == 8);
    ASSERT(mesh.getNumFaces() == 6);
    const UnitVec3 norms[] = {
        SimTK::CoordinateAxis::YCoordinateAxis(),
        SimTK::CoordinateAxis::ZCoordinateAxis(),
        -SimTK::CoordinateAxis::XCoordinateAxis(),
        -SimTK::CoordinateAxis::YCoordinateAxis(),
        SimTK::CoordinateAxis::XCoordinateAxis(),
        -SimTK::CoordinateAxis::ZCoordinateAxis()
    };
    const Vec2 textureCoords[] = {
        Vec2{0.625, 0.5},
        Vec2{0.875, 0.5},
        Vec2{0.875, 0.75},
        Vec2{0.625, 0.75}
    };
    for (int face = 0; face < mesh.getNumFaces(); face++) {
        ASSERT(mesh.getNumVerticesForFace(face) == 4);
        int numVerts = mesh.getNumVerticesForFace(face);
        UnitVec3 nextNorm = norms[face];
        for (int v = 0; v < numVerts; v++) {
            int idx = mesh.getFaceVertex(face, v);
            const Vec3& pos = mesh.getVertexPosition(idx);
            ASSERT(pos[0] == 1.0 && idx < 4 || pos[0] == -1.0 && idx >= 4);
            ASSERT(pos[1] == 1.0 && idx % 2==0 || pos[1] == -1.0 && idx % 2==1);
            // Get normals through the face/vertx indices
            UnitVec3 norm = mesh.getVertexNormal(face, v);
            ASSERT(norm == nextNorm);
            Vec2 tc = mesh.getVertexTextureCoordinate(face, v);
            if (face == 0) {
                ASSERT(textureCoords[v]==tc);
            }
        }
     }
    ASSERT(mesh.hasNormals())
    ASSERT(!mesh.hasNormalsAtVertices())
    ASSERT(mesh.hasNormalsAtFaces())
    ASSERT(mesh.hasTextureCoordinates())
    ASSERT(mesh.hasTextureCoordinatesAtFaces())
    ASSERT(!mesh.hasTextureCoordinatesAtVertices())
}

void testConvertObjFileToVisualizationFormat() {
    string file;
    file += "# This is a comment\n";
    file += "v 1.000000 1.000000 -1.000000\n";
    file += "v 1.000000 -1.000000 -1.00000\n";
    file += "v 1.000000 1.000000 1.000000\n";
    file += "v 1.000000 -1.000000 1.000000\n";
    file += "v -1.000000 1.000000 -1.00000\n";
    file += "v -1.000000 -1.000000 -1.000\n";
    file += "v -1.000000 1.000000 1.000000\n";
    file += "v -1.000000 -1.000000 1.000000\n";
    file += "vn -0.0000 1.0000 -0.0000\n";
    file += "vn -0.0000 -0.0000 1.0000\n";
    file += "vn -1.0000 -0.0000 -0.0000\n";
    file += "vn -0.0000 -1.0000 -0.0000\n";
    file += "vn 1.0000 -0.0000 -0.0000\n";
    file += "vn -0.0000 -0.0000 -1.0000\n";
    file += "vt 0.625000 0.500000\n";
    file += "vt 0.875000 0.500000\n";
    file += "vt 0.875000 0.750000\n";
    file += "vt 0.625000 0.750000\n";
    file += "vt 0.375000 0.750000\n";
    file += "vt 0.625000 1.000000\n";
    file += "vt 0.375000 1.000000\n";
    file += "vt 0.375000 0.000000\n";
    file += "vt 0.625000 0.000000\n";
    file += "vt 0.625000 0.250000\n";
    file += "vt 0.375000 0.250000\n";
    file += "vt 0.125000 0.500000\n";
    file += "vt 0.375000 0.500000\n";
    file += "vt 0.125000 0.750000\n";
    file += "s 0\n";
    file += "usemtl Material\n";
    file += "f 1/1/1 5/2/1 7/3/1 3/4/1\n";
    file += "f 4/5/2 3/4/2 7/6/2 8/7/2\n";
    file += "f 8/8/3 7/9/3 5/10/3 6/11/3\n";
    file += "f 6/12/4 2/13/4 4/5/4 8/14/4\n";
    file += "f 2/13/5 1/1/5 3/4/5 4/5/5\n";
    file += "f 6/11/6 5/10/6 1/1/6 2/13/6\n";

    PolygonalMesh mesh;
    stringstream stream(file);
    mesh.loadObjFile(stream);
    // Convert mesh to triangles, create parallel arrays of indices, normals,
    // texture-coords
    // Results below were identical to those produced using VTK's OBJReader that 
    // visualize/render correctly in viewers
    std::vector<std::array<int, 3>> triangles;
    std::vector<int> vertexIndices;
    std::vector<UnitVec3> normals;
    std::vector<Vec2> textures;
    bool hasNormals = true;  // By construction for obj files
    bool hasTextureCoordinatesAtFaces = mesh.hasTextureCoordinatesAtFaces();
    for (int face = 0; face < mesh.getNumFaces(); face++) {
        int numVerts = mesh.getNumVerticesForFace(face);
        // First triangle is 0, 1, 2, then 0, 2, 3, up to ... 0, n-2, n-1
        for (int tri = 0; tri < numVerts - 2; tri++) {
            std::array<int, 3> indices = {0, tri + 1, tri + 2};
            // Make triangle of vertices 0, tri+1, tri+2
            triangles.push_back(indices);
            if (tri == 0) {
                vertexIndices.push_back(mesh.getFaceVertex(face, indices[0]));
                vertexIndices.push_back(mesh.getFaceVertex(face, indices[1]));
            }
            vertexIndices.push_back(mesh.getFaceVertex(face, indices[2]));
            if (mesh.hasNormalsAtFaces()) {
                if (tri == 0) {
                    normals.push_back(mesh.getVertexNormal(face, indices[0]));
                    normals.push_back(mesh.getVertexNormal(face, indices[1]));
                }
                normals.push_back(mesh.getVertexNormal(face, indices[2]));
            }
            if (mesh.hasTextureCoordinatesAtFaces()) {
                if (tri == 0) {
                    textures.push_back(
                        mesh.getVertexTextureCoordinate(face, indices[0]));
                    textures.push_back(
                        mesh.getVertexTextureCoordinate(face, indices[1]));
                }
                textures.push_back(
                    mesh.getVertexTextureCoordinate(face, indices[2]));
            }
        }
    }
    std::vector<UnitVec3> expectedNormals = std::vector<UnitVec3>{
        UnitVec3{ 0., 1., 0. },  UnitVec3{ 0., 1., 0. },  UnitVec3{ 0., 1., 0. },  UnitVec3{ 0., 1., 0. },
        UnitVec3{ 0., 0., 1. },  UnitVec3{ 0., 0., 1. },  UnitVec3{ 0., 0., 1. },  UnitVec3{ 0., 0., 1. },
        UnitVec3{ -1., 0., 0. }, UnitVec3{ -1., 0., 0. }, UnitVec3{ -1., 0., 0. }, UnitVec3{ -1., 0., 0. },
        UnitVec3{ 0., -1., 0. }, UnitVec3{ 0., -1., 0. }, UnitVec3{ 0., -1., 0. }, UnitVec3{ 0., -1., 0. },
        UnitVec3{ 1., 0., 0. },  UnitVec3{ 1., 0., 0. },  UnitVec3{ 1., 0., 0. },  UnitVec3{ 1., 0., 0. },
        UnitVec3{ 0., 0., -1. }, UnitVec3{ 0., 0., -1. }, UnitVec3{ 0., 0., -1. }, UnitVec3{ 0., 0., -1. }
    };
    std::vector<Vec2> expectedTextures = std::vector<Vec2>{
        Vec2{0.625, 0.5},    Vec2{ 0.875, 0.5 },  Vec2{ 0.875, 0.75 }, Vec2{ 0.625, 0.75 },
        Vec2{ 0.375, 0.75 }, Vec2{ 0.625, 0.75 }, Vec2{ 0.625, 1. },   Vec2{ 0.375, 1. },
        Vec2{ 0.375, 0. },   Vec2{ 0.625, 0. },   Vec2{ 0.625, 0.25 }, Vec2{ 0.375, 0.25 },
        Vec2{ 0.125, 0.5 },  Vec2{ 0.375, 0.5 },  Vec2{ 0.375, 0.75 }, Vec2{ 0.125, 0.75 },
        Vec2{ 0.375, 0.5 },  Vec2{ 0.625, 0.5 },  Vec2{ 0.625, 0.75 }, Vec2{ 0.375, 0.75 },
        Vec2{ 0.375, 0.25 }, Vec2{ 0.625, 0.25 }, Vec2{ 0.625, 0.5 },  Vec2{ 0.375, 0.5 }};

    assert(normals == expectedNormals);
    assert(textures == expectedTextures);
}
void testLoadVtpFile() {
    PolygonalMesh mesh;
    string fileContent;
    fileContent += "<?xml version='1.0'?>";
    fileContent +=
        "<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian' "
        "compressor='vtkZLibDataCompressor'>";
    fileContent += "<PolyData>";
    fileContent +=
        "<Piece NumberOfPoints='4' NumberOfVerts='0' NumberOfLines='0' "
        "NumberOfStrips='0' NumberOfPolys='1'>";
    fileContent += "<PointData Normals='Normals' TCoords='TextureCoordinates'>";
    fileContent +=
        "<DataArray type='Float32' Name='Normals' NumberOfComponents='3' "
        "format='ascii' RangeMin='1' RangeMax='1'>";
    fileContent += "    0 0 1 0 0 1";
    fileContent += "    0 0 1 0 0 1";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Float32' Name='TextureCoordinates' "
        "NumberOfComponents='2' format='ascii' RangeMin='0' "
        "RangeMax='1.4142135624'>";
    fileContent += "    0 0 1 0 0 1";
    fileContent += "    1 1";
    fileContent += "</DataArray>";
    fileContent += "</PointData>";
    fileContent += "<CellData>";
    fileContent += "</CellData>";
    fileContent += "<Points>";
    fileContent +=
        "<DataArray type='Float32' Name='Array 0476C968' "
        "NumberOfComponents='3' format='ascii' RangeMin='0.70710678119' "
        "RangeMax='0.70710678119'>";
    fileContent += "    -0.5 -0.5 0 0.5 -0.5 0";
    fileContent += "    -0.5 0.5 0 0.5 0.5 0";
    fileContent += "</DataArray>";
    fileContent += "</Points>";
    fileContent += "<Verts>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Verts>";
    fileContent += "<Lines>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Lines>";
    fileContent += "<Strips>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Strips>";
    fileContent += "<Polys>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "    0 1 3 2";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "    4";
    fileContent += "</DataArray>";
    fileContent += "</Polys>";
    fileContent += "</Piece>";
    fileContent += "</PolyData>";
    fileContent += "</VTKFile>";

    std::ofstream filePtr("plane.vtp", std::ios::out);
    filePtr << fileContent;
    filePtr.close();

    mesh.loadVtpFile("plane.vtp");
    // verts = -0.5 -0.5 0, 0.5 -0.5 0, -0.5 0.5 0,  0.5 0.5 0
    // Normals = 0 0 1, 0 0 1,0 0 1, 0 0 1
    std::array<Vec3, 4> verts={Vec3{-.5, -.5, .0}, Vec3{.5, -.5, .0},
                                 Vec3{-.5, .5, .0}, Vec3{.5, .5, .0}};
    ASSERT(mesh.getNumVertices() == 4);
    ASSERT(mesh.getNumFaces() == 1);
    for (int v = 0; v < verts.size(); ++v) {
        ASSERT(mesh.getVertexNormal(v) == UnitVec3(0, 0, 1));
        ASSERT(mesh.getVertexPosition(v) == verts[v]);
    }
    ASSERT(mesh.hasNormals())
    ASSERT(mesh.hasNormalsAtVertices())
    ASSERT(!mesh.hasNormalsAtFaces())
    ASSERT(mesh.hasTextureCoordinatesAtVertices())
}

void testLoadVtpFileNoNormals() {
    PolygonalMesh mesh;
    string fileContent;
    fileContent += "<?xml version='1.0'?>";
    fileContent +=
        "<VTKFile type='PolyData' version='0.1' byte_order='LittleEndian' "
        "compressor='vtkZLibDataCompressor'>";
    fileContent += "<PolyData>";
    fileContent +=
        "<Piece NumberOfPoints='4' NumberOfVerts='0' NumberOfLines='0' "
        "NumberOfStrips='0' NumberOfPolys='1'>";
    fileContent += "<PointData TCoords='TextureCoordinates'>";
    fileContent +=
        "<DataArray type='Float32' Name='TextureCoordinates' "
        "NumberOfComponents='2' format='ascii' RangeMin='0' "
        "RangeMax='1.4142135624'>";
    fileContent += "    0 0 1 0 0 1";
    fileContent += "    1 1";
    fileContent += "</DataArray>";
    fileContent += "</PointData>";
    fileContent += "<CellData>";
    fileContent += "</CellData>";
    fileContent += "<Points>";
    fileContent +=
        "<DataArray type='Float32' Name='Array 0476C968' "
        "NumberOfComponents='3' format='ascii' RangeMin='0.70710678119' "
        "RangeMax='0.70710678119'>";
    fileContent += "    -0.5 -0.5 0 0.5 -0.5 0";
    fileContent += "    -0.5 0.5 0 0.5 0.5 0";
    fileContent += "</DataArray>";
    fileContent += "</Points>";
    fileContent += "<Verts>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Verts>";
    fileContent += "<Lines>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Lines>";
    fileContent += "<Strips>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "</DataArray>";
    fileContent += "</Strips>";
    fileContent += "<Polys>";
    fileContent +=
        "<DataArray type='Int32' Name='connectivity' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "    0 1 3 2";
    fileContent += "</DataArray>";
    fileContent +=
        "<DataArray type='Int32' Name='offsets' format='ascii' "
        "RangeMin='1e+299' RangeMax='-1e+299'>";
    fileContent += "    4";
    fileContent += "</DataArray>";
    fileContent += "</Polys>";
    fileContent += "</Piece>";
    fileContent += "</PolyData>";
    fileContent += "</VTKFile>";

    std::ofstream filePtr("planeNoNormals.vtp", std::ios::out);
    filePtr << fileContent;
    filePtr.close();

    mesh.loadVtpFile("planeNoNormals.vtp");
    // verts = -0.5 -0.5 0, 0.5 -0.5 0, -0.5 0.5 0,  0.5 0.5 0
    std::array<Vec3, 4> verts = {Vec3{-.5, -.5, .0}, Vec3{.5, -.5, .0},
                                 Vec3{-.5, .5, .0}, Vec3{.5, .5, .0}};
    ASSERT(mesh.getNumVertices() == 4);
    ASSERT(mesh.getNumFaces() == 1);
    for (int v = 0; v < verts.size(); ++v) {
        ASSERT(mesh.getVertexPosition(v) == verts[v]);
    }
    ASSERT(!mesh.hasNormals())
    ASSERT(mesh.hasTextureCoordinatesAtVertices())
}

int main() {
    try {
        testCreateMesh();
        testLoadObjFile();
        testLoadObjFileWithNormalsTexture();
        testConvertObjFileToVisualizationFormat();
        testLoadVtpFile();
        testLoadVtpFileNoNormals();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
