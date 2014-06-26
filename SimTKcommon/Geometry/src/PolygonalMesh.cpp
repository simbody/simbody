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
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Pathname.h"

#include <cassert>
#include <sstream>
#include <string>
#include <set>
#include <map>
#include <fstream>

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

//------------------------------------------------------------------------------
//                                 LOAD FILE
//------------------------------------------------------------------------------
void PolygonalMesh::loadFile(const String& pathname) {
    std::string dir,fn,ext;
    bool isAbsolutePath;
    Pathname::deconstructPathname(pathname,isAbsolutePath,dir,fn,ext);
    String lext = String::toLower(ext);

    if (lext==".obj") loadObjFile(pathname);
    else if (lext==".vtp") loadVtpFile(pathname);
    else if (lext==".stl"||lext==".stla") loadStlFile(pathname);
    else {
        SimTK_ERRCHK1_ALWAYS(!"unrecognized extension",
            "PolygonalMesh::loadFile()",
            "Unrecognized file extension on mesh file '%s':\n"
            "  expected .obj, .stl, .stla, or .vtp.", pathname.c_str());
    }
}

//------------------------------------------------------------------------------
//                              LOAD OBJ FILE
//------------------------------------------------------------------------------

// For the pathname signature just open and punt to the istream signature.
void PolygonalMesh::loadObjFile(const String& pathname) {
    std::ifstream ifs(pathname);
    SimTK_ERRCHK1_ALWAYS(ifs.good(), "PolygonalMesh::loadObjFile()",
        "Failed to open file '%s'", pathname.c_str());
    loadObjFile(ifs);
    ifs.close();
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


//------------------------------------------------------------------------------
//                              LOAD VTP FILE
//------------------------------------------------------------------------------

/* Use our XML reader to parse VTK's PolyData file format and add the polygons
found there to whatever is currently in this PolygonalMesh object. OpenSim uses
this format for its geometric objects. 

Here is a somewhat stripped down and annotated version of Kitware's description
from vtk.org:

All the metadata is case sensitive.

PolyData -- Each PolyData piece specifies a set of points and cells 
independently from the other pieces. [Simbody Note: we will read in only the
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
    Scalars - The name of the active scalars array, if any.
    Vectors - The name of the active vectors array, if any.
    Normals - The name of the active normals array, if any.
    Tensors - The name of the active tensors array, if any.
    TCoords - The name of the active texture coordinates array, if any.
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
connectivity. All the cells' point lists are concatenated together. The second
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
one or more components per value. [Simbody Note: there are also "binary" and
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
        file. This is "ascii", "binary", or "appended". [Simbody only supports
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
//                              VERTEX MAP
//------------------------------------------------------------------------------
// This is a local utility class for use in weeding out duplicate vertices.
// Set an appropriate tolerance in the constructor, then vertices all of
// whose coordinates are within tol can be considered the same vertex.
// This uses a map to make the complexity n log n.
namespace {
struct VertKey {
    VertKey(const Vec3& v, Real tol=SignificantReal) : v(v), tol(tol) {}
    bool operator<(const VertKey& other) const {
        const Vec3 diff = v - other.v;
        if (diff[0] < -tol) return true;
        if (diff[0] >  tol) return false;
        if (diff[1] < -tol) return true;
        if (diff[1] >  tol) return false;
        if (diff[2] < -tol) return true;
        if (diff[2] >  tol) return false;
        return false; // they are numerically equal
    }
    Vec3 v;
    Real tol;
};
typedef std::map<VertKey,int> VertMap;
}

//------------------------------------------------------------------------------
//                              LOAD STL FILE
//------------------------------------------------------------------------------
namespace {

class STLFile {
public:
    STLFile(const String& pathname, const PolygonalMesh& mesh) 
    :   m_pathname(pathname), m_pathcstr(pathname.c_str()),
        m_vertexTol(NTraits<float>::getSignificant()),
        m_lineNo(0), m_sigLineNo(0) 
    {   preLoadVertMap(mesh); }

    // Examine file contents to determine whether this is an ascii-format 
    // STL; otherwise it is binary.
    bool isStlAsciiFormat();

    void loadStlAsciiFile(PolygonalMesh& mesh);
    void loadStlBinaryFile(PolygonalMesh& mesh);

private:
    bool getSignificantLine(bool eofOK);

    // If we're appending to an existing mesh we'll need to preload the
    // vertex map with the existing vertices.
    void preLoadVertMap(const PolygonalMesh& mesh) {
        for (int i=0; i < mesh.getNumVertices(); ++i) {
            const Vec3& v = mesh.getVertexPosition(i);
            m_vertMap.insert(std::make_pair(VertKey(v,m_vertexTol), i));
        }
    }

    // Look for a vertex close enough to this one and return its index if found,
    // otherwise add to the mesh.
    int getVertex(const Vec3& v, PolygonalMesh& mesh) {
        const VertKey key(v, m_vertexTol);
        VertMap::const_iterator p = m_vertMap.find(key);
        int ix;
        if (p != m_vertMap.end()) 
            ix = p->second;
        else {
            ix = mesh.addVertex(v);
            m_vertMap.insert(std::make_pair(key,ix));
        }
        return ix;
    }

    // The ascii/binary determination reads some lines; counts must restart.
    void resetLineCounts() {m_lineNo=m_sigLineNo=0;}

    const String&     m_pathname;
    const char* const m_pathcstr;
    const Real        m_vertexTol;
    VertMap           m_vertMap;

    std::ifstream     m_ifs;
    int               m_lineNo;         // current line in file
    int               m_sigLineNo;      // line # not counting blanks, comments
    String            m_keyword;        // first non-blank token on line
    std::stringstream m_restOfLine;     // full line except first token
};

}

void PolygonalMesh::loadStlFile(const String& pathname) {
    bool isAbsolutePath;
    std::string directory, fileName, extension;
    Pathname::deconstructPathname(pathname, isAbsolutePath,
                                  directory, fileName, extension);
    const bool hasAsciiExt = String::toLower(extension) == ".stla";

    STLFile stlfile(pathname, *this);

    if (hasAsciiExt || stlfile.isStlAsciiFormat()) {
        stlfile.loadStlAsciiFile(*this);
    } else {
        stlfile.loadStlBinaryFile(*this);
    }
}

// The standard format for an ASCII STL file is:
//
//   solid name
//   facet normal ni nj nk
//       outer loop
//           vertex v1x v1y v1z
//           vertex v2x v2y v2z
//           vertex v3x v3y v3z
//       endloop
//   endfacet
//   ...
//   endsolid name
//
// The "..." indicates that the facet/endfacet block repeats for each face. We
// will ignore the normal on the facet line and don't care if it is present.
// The 'name' is optional and we ignore it. Extra whitespace and case are
// ignored. We'll recognize "facetnormal" and "outerloop" if the spaces are
// missing; apparently that was once allowed.
//
// Extensions (mostly for compatibility with Gazebo's STL reader):
// - Allow comment lines that begin with #, !, or $; they are skipped.
// - Allow lines beginning with 'color'; they are ignored.
// - Allow negative numbers in vertices (stl standard says only +ve).
// - Allow more than three vertices per face.
// - Allow 'outer loop'/'endloop' to be left out.
// 
// If there are multiple solids in the STL file we'll just read the first one.

// We have to decide if this is really an ascii format stl; it might be
// binary. Unfortunately, some binary stl files also start with 'solid' so
// that isn't enough. We will simply try to parse the file as ascii and then
// if that leads to an inconsistency will try binary instead.
bool STLFile::isStlAsciiFormat() {
    m_ifs.open(m_pathname);
    SimTK_ERRCHK1_ALWAYS(m_ifs.good(), "PolygonalMesh::loadStlFile()",
        "Can't open file '%s'", m_pathcstr);

    bool isAscii = false;
    if (getSignificantLine(true) && m_keyword == "solid") {
        // Still might be binary. Look for a "facet" or "endsolid" line.
        while (getSignificantLine(true)) {
            if (m_keyword=="color") continue; // ignore
            isAscii = (   m_keyword=="facet" 
                       || m_keyword=="facetnormal"
                       || m_keyword=="endsolid");
            break;
        }
    }

    m_ifs.close();
    resetLineCounts();
    return isAscii;
}


void STLFile::loadStlAsciiFile(PolygonalMesh& mesh) {
    m_ifs.open(m_pathname);
    SimTK_ERRCHK1_ALWAYS(m_ifs.good(), "PolygonalMesh::loadStlFile()",
        "Can't open file '%s'", m_pathcstr);

    Array_<int> vertices;

    // Don't allow EOF until we've seen two significant lines.
    while (getSignificantLine(m_sigLineNo >= 2)) {
        if (m_sigLineNo==1 && m_keyword == "solid") continue;
        if (m_sigLineNo>1 && m_keyword == "endsolid")
            break;
        if (m_keyword == "color") continue;

        if (m_keyword == "facet" || m_keyword == "facetnormal") {
            // We're ignoring the normal on the facet line.
            getSignificantLine(false);

            bool outerLoopSeen=false;
            if (m_keyword=="outer" || m_keyword=="outerloop") {
                outerLoopSeen = true;
                getSignificantLine(false);
            }

            // Now process vertices.
            vertices.clear();
            while (m_keyword == "vertex") {
                Vec3 vertex;
                m_restOfLine >> vertex;
                SimTK_ERRCHK2_ALWAYS(m_restOfLine.eof(), 
                    "PolygonalMesh::loadStlFile()",
                    "Error at line %d in ASCII STL file '%s':\n"
                    "  badly formed vertex.", m_lineNo, m_pathcstr);
                vertices.push_back(getVertex(vertex, mesh));
                getSignificantLine(false);
            }

            // Next keyword is not "vertex".
            SimTK_ERRCHK3_ALWAYS(vertices.size() >= 3, 
                "PolygonalMesh::loadStlFile()",
                "Error at line %d in ASCII STL file '%s':\n"
                "  a facet had %d vertices; at least 3 required.", 
                m_lineNo, m_pathcstr, vertices.size());

            mesh.addFace(vertices);

            // Vertices must end with 'endloop' if started with 'outer loop'.
            if (outerLoopSeen) {
                SimTK_ERRCHK3_ALWAYS(m_keyword=="endloop", 
                    "PolygonalMesh::loadStlFile()",
                    "Error at line %d in ASCII STL file '%s':\n"
                    "  expected 'endloop' but got '%s'.",
                    m_lineNo, m_pathcstr, m_keyword.c_str());
                getSignificantLine(false);
            }

            // Now we expect 'endfacet'.
            SimTK_ERRCHK3_ALWAYS(m_keyword=="endfacet", 
                "PolygonalMesh::loadStlFile()",
                "Error at line %d in ASCII STL file '%s':\n"
                "  expected 'endfacet' but got '%s'.",
                m_lineNo, m_pathcstr, m_keyword.c_str());
        }
    }

    // We don't care if there is extra stuff in the file.
    m_ifs.close();
}

// This is the binary STL format:
//   uint8[80] - Header (ignored)
//   uint32    - Number of triangles
//   for each triangle
//      float[3]    - normal vector (we ignore this)
//      float[3]    - vertex 1  (counterclockwise order about the normal)
//      float[3]    - vertex 2
//      float[3]    - vertex 3
//      uint16      - "attribute byte count" (ignored)
//   end
//
// TODO: the STL binary format is always little-endian, like an Intel chip.
// The code here won't work properly on a big endian machine!
void STLFile::loadStlBinaryFile(PolygonalMesh& mesh) {
    // This should never fail since the above succeeded, but we'll check.
    m_ifs.open(m_pathname, std::ios_base::binary);
    SimTK_ERRCHK1_ALWAYS(m_ifs.good(), "PolygonalMesh::loadStlFile()",
        "Can't open file '%s'", m_pathcstr);

    unsigned char header[80];
    m_ifs.read((char*)header, 80);
    SimTK_ERRCHK1_ALWAYS(m_ifs.good() && m_ifs.gcount()==80, 
        "PolygonalMesh::loadStlFile()", "Bad binary STL file '%s':\n"
        "  couldn't read header.", m_pathcstr);

    unsigned nFaces;
    m_ifs.read((char*)&nFaces, sizeof(unsigned));
    SimTK_ERRCHK1_ALWAYS(m_ifs.good() && m_ifs.gcount()==sizeof(unsigned), 
        "PolygonalMesh::loadStlFile()", "Bad binary STL file '%s':\n"
        "  couldn't read triangle count.", m_pathcstr);

    Array_<int> vertices(3);
    float vbuf[3]; unsigned short sbuf;
    const unsigned vz = 3*sizeof(float);
    for (unsigned fx=0; fx < nFaces; ++fx) {
        m_ifs.read((char*)vbuf, vz); // normal ignored
        for (int vx=0; vx < 3; ++vx) {
            m_ifs.read((char*)vbuf, vz);
            SimTK_ERRCHK3_ALWAYS(m_ifs.good() && m_ifs.gcount()==vz, 
                "PolygonalMesh::loadStlFile()", "Bad binary STL file '%s':\n"
                "  couldn't read vertex %d for face %d.", m_pathcstr, vx, fx);
            const Vec3 vertex((Real)vbuf[0], (Real)vbuf[1], (Real)vbuf[2]);
            vertices[vx] = getVertex(vertex, mesh);
        }
        mesh.addFace(vertices);
        // Now read and toss the "attribute byte count".
        m_ifs.read((char*)&sbuf,sizeof(short));
        SimTK_ERRCHK2_ALWAYS(m_ifs.good() && m_ifs.gcount()==sizeof(short), 
            "PolygonalMesh::loadStlFile()", "Bad binary STL file '%s':\n"
            "  couldn't read attribute for face %d.", m_pathcstr, fx);
    }

    // We don't care if there is extra stuff in the file.
    m_ifs.close();
}

// Return the next line from the formatted input stream, ignoring blank
// lines and comment lines, and downshifting the returned keyword. Sets
// m_keyword and m_restOfLine and increments line counts. If eofOK==false,
// issues an error message if we hit EOF, otherwise it will quitely return
// false at EOF.
bool STLFile::getSignificantLine(bool eofOK) {
    std::string line;
    std::getline(m_ifs, line);
    while (m_ifs.good()) {
        ++m_lineNo;
        m_keyword = String::trimWhiteSpace(line); // using keyword as a temp
        if (   m_keyword.empty() 
            || m_keyword[0]=='#' || m_keyword[0]=='!' || m_keyword[0]=='$') 
        {
            std::getline(m_ifs, line);
            continue; // blank or comment
        }
        // Found a significant line.
        ++m_sigLineNo;
        m_keyword.toLower();
        m_restOfLine.clear();
        m_restOfLine.str(m_keyword);
        m_restOfLine >> m_keyword; // now it's the keyword at beginning of line
        return true;
    }

    SimTK_ERRCHK2_ALWAYS(!(m_ifs.fail()||m_ifs.bad()),
        "PolygonalMesh::loadStlFile()",
        "Error at line %d in ASCII STL file '%s':\n"
        "  error while reading file.", m_lineNo, m_pathcstr);

    // Must be EOF.
    SimTK_ERRCHK2_ALWAYS(eofOK, "PolygonalMesh::loadStlFile()",
        "Error at line %d in ASCII STL file '%s':\n"
        "  unexpected end of file.", m_lineNo, m_pathcstr);
    return false;
}

//------------------------------------------------------------------------------
//                            CREATE SPHERE MESH
//------------------------------------------------------------------------------

// Use unnamed namespace to keep VertKey class and VertMap type private to 
// this file, as well as a few helper functions.
namespace {


    /* Search a list of vertices for one close enough to this one and
    return its index if found, otherwise add to the end. */
    int getVertex(const Vec3& v, VertMap& vmap, Array_<Vec3>& verts) {
        VertMap::const_iterator p = vmap.find(VertKey(v));
        if (p != vmap.end()) return p->second;
        const int ix = (int)verts.size();
        verts.push_back(v);
        vmap.insert(std::make_pair(VertKey(v),ix));
        return ix;
    }

    /* Each face comes in as below, with vertices 0,1,2 on the surface
    of a sphere or radius r centered at the origin. We bisect the edges to get
    points a',b',c', then move out from the center to make points a,b,c
    on the sphere.
             1
            /\        
           /  \
        c /____\ b      Then construct new triangles
         /\    /\            [0,b,a]
        /  \  /  \           [a,b,c]
       /____\/____\          [c,2,a]
      2      a     0         [b,1,c]
    */
    void refineSphere(Real r, VertMap& vmap, 
                             Array_<Vec3>& verts, Array_<int>&  faces) {
        assert(faces.size() % 3 == 0);
        const int nVerts = faces.size(); // # face vertices on entry
        for (int i=0; i < nVerts; i+=3) {
            const int v0=faces[i], v1=faces[i+1], v2=faces[i+2];
            const Vec3 a = r*UnitVec3(verts[v0]+verts[v2]);
            const Vec3 b = r*UnitVec3(verts[v0]+verts[v1]);
            const Vec3 c = r*UnitVec3(verts[v1]+verts[v2]);
            const int va=getVertex(a,vmap,verts), 
                      vb=getVertex(b,vmap,verts), 
                      vc=getVertex(c,vmap,verts);
            // Replace the existing face with the 0ba triangle, then add the
            // rest. Refer to the above picture.
            faces[i+1] = vb; faces[i+2] = va;
            faces.push_back(va); faces.push_back(vb); faces.push_back(vc);//abc
            faces.push_back(vc); faces.push_back(v2); faces.push_back(va);//c2a
            faces.push_back(vb); faces.push_back(v1); faces.push_back(vc);//b1c
        }
    }


    void makeOctahedralMesh(const Vec3& r, Array_<Vec3>& vertices,
                                   Array_<int>&  faceIndices) {
        vertices.push_back(Vec3( r[0],  0,  0));   //0
        vertices.push_back(Vec3(-r[0],  0,  0));   //1
        vertices.push_back(Vec3( 0,  r[1],  0));   //2
        vertices.push_back(Vec3( 0, -r[1],  0));   //3
        vertices.push_back(Vec3( 0,  0,  r[2]));   //4
        vertices.push_back(Vec3( 0,  0, -r[2]));   //5
        int faces[8][3] = {{0, 2, 4}, {4, 2, 1}, {1, 2, 5}, {5, 2, 0}, 
                           {4, 3, 0}, {1, 3, 4}, {5, 3, 1}, {0, 3, 5}};
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 3; j++)
                faceIndices.push_back(faces[i][j]);
    }

} // end unnamed namespace

/*static*/ PolygonalMesh PolygonalMesh::
createSphereMesh(Real radius, int resolution) {
    SimTK_ERRCHK1_ALWAYS(radius > 0, "PolygonalMesh::createSphereMesh()",
        "Radius %g illegal.", radius);
    SimTK_ERRCHK1_ALWAYS(resolution >= 0, "PolygonalMesh::createSphereMesh()",
        "Resolution %d illegal.", resolution);

    Array_<Vec3> vertices;
    Array_<int> faceIndices;
    makeOctahedralMesh(Vec3(radius), vertices, faceIndices);

    VertMap vmap;
    for (unsigned i=0; i < vertices.size(); ++i)
        vmap[vertices[i]] = i;

    int level = resolution;
    while (level > 0) {
        refineSphere(radius, vmap, vertices, faceIndices);
        --level;
    }

    PolygonalMesh sphere;
    for (unsigned i=0; i < vertices.size(); ++i)
        sphere.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 3) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+3);
        sphere.addFace(verts);
    }

    return sphere; // just a shallow copy
}


//------------------------------------------------------------------------------
//                            CREATE BRICK MESH
//------------------------------------------------------------------------------
// Resolution 0 is just the four corners and six quads.
// Resolution n means divide longest edge with n extra vertices; then make the
// other faces about the same size.

// Extend the unnamed namespace with another local function.
namespace {

    // Given the number of vertices along the x,y,z edges and the three (i,j,k) 
    // coordinates of a particular vertex, return a packed representation of the
    // (i,j,k) location that we can use as an index. Using a long long here
    // just to avoid trouble, an int probably would have been fine.
    long long locCode(int nv[3], int i, int j, int k)
    {   return (long long)i*nv[1]*nv[2] + (long long)j*nv[2] + (long long)k; }

} // end of unnamed namespace

/*static*/ PolygonalMesh PolygonalMesh::
createBrickMesh(const Vec3& halfDims, int resolution) {
    SimTK_ERRCHK3_ALWAYS(halfDims > 0, "PolygonalMesh::createBrickMesh()",
        "Bad brick dimensions %g %g %g.", halfDims[0], halfDims[1], halfDims[2]);
    SimTK_ERRCHK1_ALWAYS(resolution >= 0, "PolygonalMesh::createBrickMesh()",
        "Resolution %d illegal.", resolution);

    const Vec3 dims(2*halfDims);

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

