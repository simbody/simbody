/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
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

#include "PolygonalMeshImpl.h"
#include "SimTKcommon/internal/Xml.h"

#include <cassert>
#include <sstream>
#include <string>
#include <set>
#include <map>

using namespace SimTK;

//==============================================================================
//                            POLYGONAL MESH
//==============================================================================

// If the handle is empty, reconstruct it to be an owner handle whose
// implementation is present but contains no vertices.
void PolygonalMesh::initializeHandleIfEmpty() {
    if (isEmptyHandle())
        new(this) PolygonalMesh(new PolygonalMeshImpl());
}

// default (shallow) copy constructor, copy assignment, destructor

void PolygonalMesh::clear() {
    if (!isEmptyHandle()) updImpl().clear();
}

int PolygonalMesh::getNumFaces() const {
    return isEmptyHandle() ? 0 : getImpl().faceVertexStart.size()-1;
}

int PolygonalMesh::getNumVertices() const {
    return isEmptyHandle() ? 0 : getImpl().vertices.size();
}

const Vec3& PolygonalMesh::getVertexPosition(int vertex) const {
    assert(0 <= vertex && vertex < getNumVertices());
    return getImpl().vertices[vertex];
}

int PolygonalMesh::getNumVerticesForFace(int face) const {
    assert(0 <= face && face < getNumFaces());
    const Array_<int>& faceVertexStart = getImpl().faceVertexStart;
    return faceVertexStart[face+1]-faceVertexStart[face];
}

int PolygonalMesh::getFaceVertex(int face, int vertex) const {
    assert(face >= 0 && face < getNumFaces());
    assert(vertex >= 0 && vertex < getNumVerticesForFace(face));
    return getImpl().faceVertexIndex[getImpl().faceVertexStart[face]+vertex];
}

int PolygonalMesh::addVertex(const Vec3& position) {
    initializeHandleIfEmpty();
    updImpl().vertices.push_back(position);
    return getImpl().vertices.size()-1;
}

int PolygonalMesh::addFace(const Array_<int>& vertices) {
    initializeHandleIfEmpty();
    for (int i = 0; i < (int) vertices.size(); i++)
        updImpl().faceVertexIndex.push_back(vertices[i]);

    // faceVertexStart is preloaded to have its first element 0 before any
    // faces have been added. So the back() element of faceVertexStart is
    // already the starting entry for the face we're adding.
    // This is where the *next* face will begin.
    updImpl().faceVertexStart.push_back(getImpl().faceVertexIndex.size());
    // The current face start is now at end()-2 (back()-1).
    return getImpl().faceVertexStart.size()-2;
}

PolygonalMesh& PolygonalMesh::scaleMesh(Real scale) {
    if (!isEmptyHandle()) {
        Array_<Vec3>& vertices = updImpl().vertices;
        for (int i = 0; i < (int) vertices.size(); i++)
            vertices[i] *= scale;
    }
    return *this;
}

PolygonalMesh& PolygonalMesh::transformMesh(const Transform& X_AM) {
    if (!isEmptyHandle()) {
        Array_<Vec3>& vertices = updImpl().vertices;
        for (int i = 0; i < (int) vertices.size(); i++)
            vertices[i] = X_AM*vertices[i];
    }
    return *this;
}

void PolygonalMesh::loadObjFile(std::istream& file) {
    const char* methodName = "PolygonalMesh::loadObjFile()";
    SimTK_ERRCHK_ALWAYS(file.good(), methodName,
        "The supplied std::istream object was not in good condition"
        " on entrance -- did you check whether it opened successfully?");

    std::string line;
    Array_<int> indices;
    int initialVertices = getNumVertices();
    while (!file.eof()) {
        SimTK_ERRCHK_ALWAYS(file.good(), methodName,
            "An error occurred while reading the input file.");

        std::getline(file, line);
        while (line.size() > 0 && line[line.size()-1] == '\\') {
            line[line.size()-1] = ' ';
            std::string continuation;
            std::getline(file, continuation);
            line += continuation;
        }
        std::stringstream s(line);
        std::string command;
        s >> command;
        if (command == "v") {
            // A vertex
            
            Real x, y, z;
            s >> x;
            s >> y;
            s >> z;
            SimTK_ERRCHK1_ALWAYS(!s.fail(), methodName,
                "Found invalid vertex description: %s", line.c_str());
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

/* Use our XML reader to parse VTK's PolyData file format and add the polygons
found there to whatever is currently in this PolygonalMesh object. OpenSim uses
this format for its geometric objects. 

Here is a somewhat stripped down and annotated version of Kitware's description
from vtk.org:

All the metadata is case sensitive.

PolyData -- Each PolyData piece specifies a set of points and cells 
independently from the other pieces. [SimTK Note: we will read in only the
first Piece element.] The points are described explicitly by the
Points element. The cells are described explicitly by the Verts, Lines, Strips,
and Polys elements.
    <VTKFile type="PolyData" ...>
        <PolyData>
            <Piece NumberOfPoints="#" NumberOfVerts="#" NumberOfLines="#"
                   NumberOfStrips="#" NumberOfPolys="#">
                <PointData>...</PointData>
                <CellData>...</CellData>
                <Points>...</Points>
                <Verts>...</Verts>
                <Lines>...</Lines>
                <Strips>...</Strips>
                <Polys>...</Polys>
            </Piece>
        </PolyData>
    </VTKFile>

PointData and CellData -- Every dataset describes the data associated with 
its points and cells with PointData and CellData XML elements as follows:
    <PointData Scalars="Temperature" Vectors="Velocity">
        <DataArray Name="Velocity" .../>
        <DataArray Name="Temperature" .../>
        <DataArray Name="Pressure" .../>
    </PointData>

VTK allows an arbitrary number of data arrays to be associated with the points 
and cells of a dataset. Each data array is described by a DataArray element 
which, among other things, gives each array a name. The following attributes 
of PointData and CellData are used to specify the active arrays by name:
    Scalars � The name of the active scalars array, if any.
    Vectors � The name of the active vectors array, if any.
    Normals � The name of the active normals array, if any.
    Tensors � The name of the active tensors array, if any.
    TCoords � The name of the active texture coordinates array, if any.
That is, for each attribute of the form Sometype="Somename" there must be a 
DataArray element with attribute Name="Somename" whose text contains 
NumberOfPoints values each of type Sometype.

Points -- The Points element explicitly defines coordinates for each point 
individually. It contains one DataArray element describing an array with 
three components per value, each specifying the coordinates of one point.
    <Points>
        <DataArray NumberOfComponents="3" .../>
    </Points>

Verts, Lines, Strips, and Polys -- The Verts, Lines, Strips, and Polys elements
define cells explicitly by specifying point connectivity. Cell types are 
implicitly known by the type of element in which they are specified. Each 
element contains two DataArray elements. The first array specifies the point 
connectivity. All the cells� point lists are concatenated together. The second
array specifies the offset into the connectivity array for the end of each
cell.
    <Polys>
        <DataArray type="Int32" Name="connectivity" .../>
        <DataArray type="Int32" Name="offsets" .../>
    </Polys>
(same format for Verts, Lines, and Strips)

DataArray -- All of the data and geometry specifications use DataArray elements
to describe their actual content as follows:

The DataArray element stores a sequence of values of one type. There may be 
one or more components per value. [SimTK Note: there are also "binary" and
"appended" formats which we do not support -- be sure to check that the
format attribute for every DataArray is "ascii".]
    <DataArray type="Int32" Name="offsets" format="ascii">
    10 20 30 ... </DataArray>

The attributes of the DataArray elements are described as follows:
    type -- The data type of a single component of the array. This is one of 
        Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Float32, 
        Float64. 
    Name -- The name of the array. This is usually a brief description of the
        data stored in the array. [Note that the PolyData element uses the 
        DataArray Name attribute to figure out what's being provided.]
    NumberOfComponents -- The number of components per value in the array.
    format -- The means by which the data values themselves are stored in the
        file. This is "ascii", "binary", or "appended". [SimTK only supports
        "ascii".]
    format="ascii" -- The data are listed in ASCII directly inside the 
        DataArray element. Whitespace is used for separation.
*/
void PolygonalMesh::loadVtpFile(const String& pathname) {
  try
  { const char* method = "PolygonalMesh::loadVtpFile()";
    Xml vtp(pathname);
    // The file has been read in and parsed into memory by the Xml system.

    SimTK_ERRCHK1_ALWAYS(vtp.getRootTag() == "VTKFile", method,
        "Expected to see document tag <VTKFile> but saw <%s> instead.",
        vtp.getRootTag().c_str());
    // This is a VTKFile document.

    Xml::Element root = vtp.getRootElement();
    SimTK_ERRCHK1_ALWAYS(root.getRequiredAttributeValue("type") == "PolyData",
        method, "Expected VTK file type='PolyData' but got type='%s'.",
        root.getRequiredAttributeValue("type").c_str());
    // This is a VTK PolyData document.

    Xml::Element polydata = root.getRequiredElement("PolyData");
    Xml::Element piece    = polydata.getRequiredElement("Piece");
    Xml::Element points   = piece.getRequiredElement("Points");
    const int numPoints = 
        piece.getRequiredAttributeValueAs<int>("NumberOfPoints");
    const int numPolys  = 
        piece.getRequiredAttributeValueAs<int>("NumberOfPolys");

    // Remember this because we'll have to use it to adjust the indices we use 
    // when referencing the vertices we're about to read in. This number is
    // the index that our first vertex will be assigned.
    const int firstVertex = getNumVertices();

    // The lone DataArray element in the Points element contains the points'
    // coordinates. Read it in as a Vector of Vec3s.
    Xml::Element pointData = points.getRequiredElement("DataArray");
    SimTK_ERRCHK1_ALWAYS(pointData.getRequiredAttributeValue("format")
                         == "ascii", method, 
        "Only format=\"ascii\" is supported for .vtp file DataArray elements,"
        " got format=\"%s\" for Points DataArray.",
        pointData.getRequiredAttributeValue("format").c_str());

    Vector_<Vec3> coords = 
        points.getRequiredElementValueAs< Vector_<Vec3> >("DataArray");

    SimTK_ERRCHK2_ALWAYS(coords.size() == numPoints, method,
        "Expected coordinates for %d points but got %d.",
        numPoints, coords.size());

    // Now that we have the point coordinates, use them to create the vertices
    // in our mesh.
    for (int i=0; i < numPoints; ++i)
        addVertex(coords[i]);

    // Polys are given by a connectivity array which lists the points forming
    // each polygon in a long unstructured list, then an offsets array, one per
    // polygon, which gives the index+1 of the *last* connectivity entry for
    // each polygon.
    Xml::Element polys = piece.getRequiredElement("Polys");

    // Find the connectivity and offset DataArrays.
    Xml::Element econnectivity, eoffsets;
    for (Xml::element_iterator p = polys.element_begin("DataArray");
         p != polys.element_end(); ++p) 
    {       
        const String& name = p->getRequiredAttributeValue("Name");
        SimTK_ERRCHK2_ALWAYS(p->getRequiredAttributeValue("format")
                             == "ascii", method, 
            "Only format=\"ascii\" is supported for .vtp file DataArray"
            " elements, but format=\"%s\" for DataArray '%s'.",
            p->getRequiredAttributeValue("format").c_str(), name.c_str());

        if (name == "connectivity") econnectivity = *p;
        else if (name == "offsets") eoffsets = *p; 
    }

    SimTK_ERRCHK_ALWAYS(econnectivity.isValid() && eoffsets.isValid(), method, 
        "Expected to find a DataArray with name='connectivity' and one with"
        " name='offsets' in the VTK PolyData file's <Polys> element but at"
        " least one of them was missing.");

    // Read in the arrays.
    Array_<int> offsets = eoffsets.getValueAs< Array_<int> >();
    // Size may have changed if file is bad.
    SimTK_ERRCHK2_ALWAYS(offsets.size() == numPolys, method,
        "The number of offsets (%d) should have matched the stated "
        " NumberOfPolys value (%d).", offsets.size(), numPolys);

    // We expect that the last entry in the offsets array is one past the
    // end of the last polygon described in the connectivity array and hence
    // is the size of the connectivity array.
    const int expectedSize = numPolys ? offsets.back() : 0;
    Array_<int> connectivity = econnectivity.getValueAs< Array_<int> >();

    SimTK_ERRCHK2_ALWAYS(connectivity.size()==expectedSize, method,
        "The connectivity array was the wrong size (%d). It should"
        " match the last entry in the offsets array which was %d.",
        connectivity.size(), expectedSize);

    int startPoly = 0;
    for (int i=0; i < numPolys; ++i) {
        // Now read in the face in [startOffs,endOffs]
        addFace(connectivity(startPoly, offsets[i]-startPoly));
        startPoly = offsets[i]; // move to the next poly
    }

  } catch (const std::exception& e) {
      // This will throw a new exception with an enhanced message that
      // includes the original one.
      SimTK_ERRCHK2_ALWAYS(!"failed", "PolygonalMesh::loadVtpFile()",
          "Attempt to load a VTK PolyData (.vtp) file from file name"
          " '%s' failed with message:\n  %s", pathname.c_str(), e.what());
  }
}

//------------------------------------------------------------------------------
//                            CREATE SPHERE MESH
//------------------------------------------------------------------------------
///*static*/ PolygonalMesh PolygonalMesh::
//createSphereMesh(Real radius, int resolution) {
//
//}


//------------------------------------------------------------------------------
//                            CREATE BRICK MESH
//------------------------------------------------------------------------------
// Resolution 0 is just the four corners and six quads.
// Resolution n means divide longest edge with n extra vertices; then make the
// other faces about the same size.

// Given the number of vertices along the x,y,z edges and the three 
// (i,j,k) coordinates of a particular vertex, return a packed representation
// of the (i,j,k) location that we can use as an index. Using a long long here
// just to avoid trouble, an int probably would have been fine.
static long long locCode(int nv[3], int i, int j, int k)
{   return (long long)i*nv[1]*nv[2] + (long long)j*nv[2] + (long long)k; }

/*static*/ PolygonalMesh PolygonalMesh::
createBrickMesh(const Vec3& halfDims, int resolution) {
    SimTK_ERRCHK3_ALWAYS(halfDims > 0, "PolygonalMesh::createBrickMesh()",
        "Bad brick dimensions %g %g %g.", halfDims[0], halfDims[1], halfDims[2]);
    SimTK_ERRCHK1_ALWAYS(resolution >= 0, "PolygonalMesh::createBrickMesh()",
        "Resolution %d illegal.", resolution);

    const Vec3 dims(2.*halfDims);

    PolygonalMesh brick;
    const Real longest = max(dims);
    const Real edgeLengthTarget = longest/(resolution+1);
    int nv[3]; // number of vertices along each edge
    for (int i=0; i<3; ++i)
        nv[i] = 1 + std::max((int)(dims[i]/edgeLengthTarget + 0.49), 1);
    const Vec3 edgeLengths(dims[0]/(nv[0]-1), dims[1]/(nv[1]-1), 
                           dims[2]/(nv[2]-1)); 

    // Add regularly-spaced vertices on the surfaces.
    std::map<long long,int> vertLoc2Vert; // map i,j,k -> vertex number
    for (int i=0; i < nv[0]; ++i) {
        bool xface = (i==0 || i==nv[0]-1);
        const Real iloc = i*edgeLengths[0] - halfDims[0];
        for (int j=0; j < nv[1]; ++j) {
            bool yface = (j==0 || j==nv[1]-1);
            const Real jloc = j*edgeLengths[1] - halfDims[1];
            for (int k=0; k < nv[2]; ++k) {
                bool zface = (k==0 || k==nv[2]-1);
                if (!(xface||yface||zface)) 
                    continue; // skip interior vertices
                const Real kloc = k*edgeLengths[2] - halfDims[2];
                const int vnum = brick.addVertex(Vec3(iloc,jloc,kloc));
                vertLoc2Vert[locCode(nv,i,j,k)] = vnum;
            }
        }
    }

    // Add quad faces oriented with normal outwards.
    Array_<int> face(4);

    // This is the -x surface (y,z rectangle).
    for (int j=0; j < nv[1]-1; ++j)
        for (int k=0; k < nv[2]-1; ++k) {
            face[0] = vertLoc2Vert[locCode(nv,0,j,  k+1)];
            face[1] = vertLoc2Vert[locCode(nv,0,j+1,k+1)];
            face[2] = vertLoc2Vert[locCode(nv,0,j+1,k)];
            face[3] = vertLoc2Vert[locCode(nv,0,j,  k)];
            brick.addFace(face);
        }
    // This is the +x surface (y,z rectangle).
    for (int j=0; j < nv[1]-1; ++j)
        for (int k=0; k < nv[2]-1; ++k) {
            face[3] = vertLoc2Vert[locCode(nv,nv[0]-1,j,  k+1)];
            face[2] = vertLoc2Vert[locCode(nv,nv[0]-1,j+1,k+1)];
            face[1] = vertLoc2Vert[locCode(nv,nv[0]-1,j+1,k)];
            face[0] = vertLoc2Vert[locCode(nv,nv[0]-1,j,  k)];
            brick.addFace(face);
        }
    // This is the -y surface (x,z rectangle).
    for (int i=0; i < nv[0]-1; ++i)
        for (int k=0; k < nv[2]-1; ++k) {
            face[0] = vertLoc2Vert[locCode(nv,i,  0,k)];
            face[1] = vertLoc2Vert[locCode(nv,i+1,0,k)];
            face[2] = vertLoc2Vert[locCode(nv,i+1,0,k+1)];
            face[3] = vertLoc2Vert[locCode(nv,i,  0,k+1)];
            brick.addFace(face);
        }
    // This is the +y surface (x,z rectangle).
    for (int i=0; i < nv[0]-1; ++i)
        for (int k=0; k < nv[2]-1; ++k) {
            face[3] = vertLoc2Vert[locCode(nv,i,  nv[1]-1,k)];
            face[2] = vertLoc2Vert[locCode(nv,i+1,nv[1]-1,k)];
            face[1] = vertLoc2Vert[locCode(nv,i+1,nv[1]-1,k+1)];
            face[0] = vertLoc2Vert[locCode(nv,i,  nv[1]-1,k+1)];
            brick.addFace(face);
        }
    // This is the -z surface (x,y rectangle).
    for (int i=0; i < nv[0]-1; ++i)
        for (int j=0; j < nv[1]-1; ++j) {
            face[0] = vertLoc2Vert[locCode(nv,i,  j+1,0)];
            face[1] = vertLoc2Vert[locCode(nv,i+1,j+1,0)];
            face[2] = vertLoc2Vert[locCode(nv,i+1,j,  0)];
            face[3] = vertLoc2Vert[locCode(nv,i,  j,  0)];
            brick.addFace(face);
        }
    // This is the +z surface (x,y rectangle).
    for (int i=0; i < nv[0]-1; ++i)
        for (int j=0; j < nv[1]-1; ++j) {
            face[3] = vertLoc2Vert[locCode(nv,i,  j+1,nv[2]-1)];
            face[2] = vertLoc2Vert[locCode(nv,i+1,j+1,nv[2]-1)];
            face[1] = vertLoc2Vert[locCode(nv,i+1,j,  nv[2]-1)];
            face[0] = vertLoc2Vert[locCode(nv,i,  j,  nv[2]-1)];
            brick.addFace(face);
        }

    return brick; // just a shallow copy
}


//------------------------------------------------------------------------------
//                           CREATE CYLINDER MESH
//------------------------------------------------------------------------------
// Minimum end surface is a hexagon (resolution 0).
// Think of the axis as the +z axis, with the end caps in x-y planes at
// height -z and +z.
/*static*/ PolygonalMesh PolygonalMesh::
createCylinderMesh(const UnitVec3& axis, Real radius, Real halfLength, 
                   int resolution) 
{
    SimTK_ERRCHK1_ALWAYS(radius > 0, "PolygonalMesh::createCylinderMesh()",
        "Bad radius %g.", radius);
    SimTK_ERRCHK1_ALWAYS(halfLength > 0, "PolygonalMesh::createCylinderMesh()",
        "Bad half length %g.", halfLength);
    SimTK_ERRCHK1_ALWAYS(resolution >= 0, "PolygonalMesh::createCylinderMesh()",
        "Resolution %d illegal.", resolution);

    int rezAround = 6*(resolution+1);

    Real angle = 2*Pi/rezAround;
    Real chordLen = 2*radius*std::sin(angle/2);
    int rezRadial = (int)(radius/chordLen+0.5);
    Real edgeLenRad = radius/rezRadial;

    int rezAlong  = 1 + std::max((int)(halfLength/edgeLenRad + 0.5), 1);
    Real edgeLenAlong = 2*halfLength/(rezAlong-1);

    PolygonalMesh cyl;
    Rotation R_ZG = Rotation(axis, ZAxis);

    int nv[3] = {rezRadial,rezAround,rezAlong};
    // Do the tube.
    std::map<long long, int> rak2Vert;
    for (int k=0; k < rezAlong; ++k) {
        bool isEndCap = (k==0 || k==rezAlong-1);
        Real z = -halfLength + k*edgeLenAlong;
        for (int a=0; a < rezAround; ++a) {
            Real x = edgeLenRad*std::sin(a*angle);
            Real y = edgeLenRad*std::cos(a*angle);
            for (int r=1; r <= rezRadial; ++r) {
                if (r < rezRadial && !isEndCap)
                    continue; // skip interior vertices
                int vnum = cyl.addVertex(R_ZG*Vec3(r*x,r*y,z));
                rak2Vert[locCode(nv,r,a,k)] = vnum;
            }
        }
    }

    Array_<int> qface(4), tface(3);
    for (int a=0; a < rezAround; ++a) {
        int ap = (a+1)%rezAround;
        for (int k=0; k < rezAlong-1; ++k) {
            int r = rezRadial;
            qface[0] = rak2Vert[locCode(nv,r, a,k)];
            qface[1] = rak2Vert[locCode(nv,r, a,k+1)];
            qface[2] = rak2Vert[locCode(nv,r,ap,k+1)];
            qface[3] = rak2Vert[locCode(nv,r,ap,k)];
            cyl.addFace(qface);
        }
    }

    // Add central face vertices.
    rak2Vert[locCode(nv,0,0,0)] = cyl.addVertex(R_ZG*Vec3(0,0,-halfLength));
    rak2Vert[locCode(nv,0,0,rezAlong-1)] = cyl.addVertex(R_ZG*Vec3(0,0,halfLength));

    // Tri faces from center to first ring.
    for (int a=0; a < rezAround; ++a) {
        int ap = (a+1)%rezAround; 
        tface[0] = rak2Vert[locCode(nv,0,0, 0)];
        tface[1] = rak2Vert[locCode(nv,1,a, 0)];
        tface[2] = rak2Vert[locCode(nv,1,ap,0)];
        cyl.addFace(tface);
        tface[0] = rak2Vert[locCode(nv,0,0, rezAlong-1)];
        tface[1] = rak2Vert[locCode(nv,1,ap, rezAlong-1)];
        tface[2] = rak2Vert[locCode(nv,1,a,rezAlong-1)];
        cyl.addFace(tface);
    }


    // Quad faces from first ring out.
    for (int a=0; a < rezAround; ++a) {
        int ap = (a+1)%rezAround; 
        for (int r=1; r <= rezRadial-1; ++r) {
        qface[0] = rak2Vert[locCode(nv,r,  a, 0)];
        qface[1] = rak2Vert[locCode(nv,r+1,a, 0)];
        qface[2] = rak2Vert[locCode(nv,r+1,ap,0)];
        qface[3] = rak2Vert[locCode(nv,r,  ap,0)];
        cyl.addFace(qface);
        qface[3] = rak2Vert[locCode(nv,r,  a, rezAlong-1)];
        qface[2] = rak2Vert[locCode(nv,r+1,a, rezAlong-1)];
        qface[1] = rak2Vert[locCode(nv,r+1,ap,rezAlong-1)];
        qface[0] = rak2Vert[locCode(nv,r,  ap,rezAlong-1)];
        cyl.addFace(qface);
        }
    }

    return cyl; // just a shallow copy
}

