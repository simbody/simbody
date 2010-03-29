/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-10 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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
#include "SimTKcommon/internal/Xml.h"

#include <cassert>
#include <sstream>
#include <string>

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
    Scalars — The name of the active scalars array, if any.
    Vectors — The name of the active vectors array, if any.
    Normals — The name of the active normals array, if any.
    Tensors — The name of the active tensors array, if any.
    TCoords — The name of the active texture coordinates array, if any.
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
connectivity. All the cells’ point lists are concatenated together. The second
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
    // coordinates. Read it in as a Vector
    Vector_<Vec3> coords = 
        points.getRequiredElementValueAs< Vector_<Vec3> >("DataArray");

    SimTK_ERRCHK2_ALWAYS(coords.size() == numPoints, method,
        "Expected coordinates for %d points but got %d.",
        numPoints, coords.size());

    for (int i=0; i < numPoints; ++i)
        addVertex(coords[i]);


    //const String& coordString = points.getRequiredElementValue("DataArray");
    //std::stringstream cstream(coordString);
    //for (int i=1; i <= numPoints; ++i) {
    //    Vec3 pt; cstream >> pt; 
    //    SimTK_ERRCHK_ALWAYS(!cstream.fail(), method,
    //        "Couldn't parse the <Points> coordinates DataArray.");
    //    // Eof is only allowed when reading the last point.
    //    if (i != numPoints)
    //        SimTK_ERRCHK2_ALWAYS(!cstream.eof(), method,
    //            "Expected coordinates for %d points but ran out after %d.",
    //            numPoints, i);
    //    addVertex(pt);
    //}

    // Polys are given by a connectivity array which lists the points forming
    // each polygon in a long unstructured list, then an offsets array, one per
    // polygon, which gives the index+1 of the *last* connectivity entry for
    // each polygon.
    Xml::Element polys = piece.getRequiredElement("Polys");

    // Find the connectivity and offset DataArrays.
    Xml::Element econnectivity, eoffsets;
    for (Xml::element_iterator p = polys.element_begin("DataArray");
         p != polys.element_end(); ++p) 
    {   const String& name = p->getRequiredAttributeValue("Name");
        if (name == "connectivity") econnectivity = *p;
        else if (name == "offsets") eoffsets = *p; }

    SimTK_ERRCHK_ALWAYS(econnectivity.isValid() && eoffsets.isValid(), method, 
        "Expected to find a DataArray with name='connectivity' and one with"
        " name='offsets' in the VTK PolyData file's <Polys> element but at"
        " least one of them was missing.");

    // Read in the arrays.
    Array_<int> offsets(numPolys);
    eoffsets.getValueAs<Array_<int> >(offsets);
    // Size may have changed if file is bad.
    SimTK_ERRCHK2_ALWAYS(offsets.size() == numPolys, method,
        "The number of offsets (%d) should have matched the stated "
        " NumberOfPolys value (%d).", offsets.size(), numPolys);

    // We expect that the last entry in the offsets array is one past the
    // end of the last polygon described in the connectivity array and hence
    // is the size of the connectivity array.
    const int expectedSize = numPolys ? offsets.back() : 0;
    Array_<int> connectivity(expectedSize);
    econnectivity.getValueAs<Array_<int> >(connectivity);

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

} // namespace SimTK
