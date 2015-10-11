#ifndef SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
#define SimTK_SimTKCOMMON_POLYGONAL_MESH_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
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

#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/PrivateImplementation.h"

namespace SimTK {

class PolygonalMesh;
class PolygonalMeshImpl;

// We only want the template instantiation to occur once. This symbol is defined
// in the SimTKcommon compilation unit that defines the PolygonalMesh class but
// should not be defined any other time.
#ifndef SimTK_SIMTKCOMMON_DEFINING_POLYGONALMESH
    extern template class PIMPLHandle<PolygonalMesh, PolygonalMeshImpl, true>;
#endif

/** This class provides a description of a mesh made of polygonal faces (not
limited to triangles). Its primary purpose is for loading geometry from files,
which can then be used for visualization or collision detection. For example,
the following lines load a mesh from a Wavefront OBJ file, then create a
DecorativeMesh from it.
@code
    PolygonalMesh mesh;
    mesh.loadObjFile("teapot.obj");
    DecorativeMesh decoration(mesh);
@endcode
You can also read a polygon mesh from a VTK PolyData (.vtp) file, or an STL
file (.stl) that is in ascii or binary format. You can also build meshes
programmatically, and some static methods are provided here for generating some
common shapes. If you don't know what kind of file you have, you can attempt to
read it with the loadFile() method which will examine the file extension to
determine the expected format.

The mesh has its own local frame and vertex locations are given in that
frame. You can scale and transform the vertices relative to that frame
(changing the values stored in the mesh) but more commonly the mesh will be
placed on a body relative to that body's frame, meaning you can re-use the
same mesh in various places.

We expect this to be a large object so give it shared (reference) semantics;
that is, the copy constructor and copy assignment default to shallow copies
(both handles will refer to the same data). If you want to make a deep
(non-shared) copy of a PolygonalMesh, use the copyAssign() method provided by
the PIMPLHandle base class. **/
class SimTK_SimTKCOMMON_EXPORT PolygonalMesh
:   public PIMPLHandle<PolygonalMesh, PolygonalMeshImpl, true> {
public:
    /** Create an empty %PolygonalMesh, with no vertices or faces. **/
    PolygonalMesh() {}

    /** Create a sphere-shaped mesh, with roughly uniform mesh elements.

    @param[in]  radius
        The radius of the underlying sphere. Vertices of the mesh will be on
        the sphere, with mesh elements somewhat inside.
    @param[in]  resolution
        Control for how dense a mesh to produce. Resolution 0 will produce an
        octahedron (8 triangular faces). Resolution 1 (the default) gives
        32 faces, resolution 2 gives 128. In general for resolution n there
        will be 2*4^(n+1) faces.
    @return A %PolygonalMesh representing a sphere of the specified radius. **/
    static PolygonalMesh createSphereMesh(Real  radius,
                                          int   resolution = 1);

    /** Create a brick-shaped mesh. A brick is a rectangular solid (a box)
    centered at and aligned with the mesh local frame. Note that its size is
    given with \e half dimensions. By default you will just get two mesh faces
    along the longest edge of the brick, with all other edges roughly the same
    size. You can control the mesh density with the \a resolution parameter.

    @param[in]  halfDims
        The half-dimensions of the brick. The extreme vertices are at
        -halfDims and +halfDims, so the brick is centered around the mesh
        local frame.
    @param[in]  resolution
        Control for how dense a mesh to produce. For this shape, \a resolution
        is interpreted as the number of extra vertices to insert in the
        \e longest edge of the brick. Extra vertices are inserted into the
        shorter edges if needed to keep the edge lengths approximately
        uniform for every mesh face. \a resolution=0 gives only
        vertices at the corners; the default is 1 meaning that the longest
        edge is split once.
    @return A %PolygonalMesh representing a brick of the requested size.

    <h3>Controlling the mesh density:</h3>
    If you want a brick mesh where all the edges in the mesh are roughly the
    same length, say \c wantEdgeLength, set \a resolution like this:
    @code
    Real wantEdgeLength = ...;
    Vec3 halfDims = ...;
    int resolution = (int)(max(halfDims)/wantEdgeLength + 0.5);
    @endcode

    If you want a brick mesh where all the edges are roughly the same length
    as the shortest edge of the brick, just set
    <code>wantEdgeLength=min(halfDims)</code> in the above calculation. **/
    static PolygonalMesh createBrickMesh(const Vec3& halfDims,
                                         int resolution = 1);

    /** Create a cylinder-shaped mesh, with the long axis in a given
    direction. By default you'll get a 12 sided polygon as the base and
    elements of roughly similar dimension along the edges. You can control the
    mesh density with the \a resolution parameter.

    @param[in]  axis
        The central axis direction of the cylinder, in the mesh local frame.
        This can be provided using the constants XAxis, YAxis, or ZAxis, or
        you can provide a unit vector in any direction.
    @param[in]  radius
        The cylinder radius.
    @param[in]  halfLength
        Half the length of the cylinder along its axis. The bases are at
        -halfLength and +halfLength along the \a axis, so the cylinder is
        centered around the mesh local frame origin.
    @param[in]  resolution
        Control for how dense a mesh to produce (see below for details).
    @return A %PolygonalMesh representing a cylinder of the requested dimensions
        and orientation.

    <h3>Controlling the mesh density:</h3>
    At resolution 0 the base is a hexagon with six triangular faces, and the
    tube is meshed with quad faces that are about as long
    as the diameter of the base. Resolution 1 (the default) makes the base
    a 12-sided polygon and introduces an intermediate 12-sided polygon of
    have the diameter. There will be triangles in the center still, but
    quad faces between the polygons. The length of the tube faces will be
    reduced to match. Higher resolutions refine the mesh similarly. **/
    static PolygonalMesh createCylinderMesh(const UnitVec3& axis,
                                            Real            radius,
                                            Real            halfLength,
                                            int             resolution=1);

    /** Restore this %PolygonalMesh to its default-constructed state, meaning
    that it will contain no vertices or faces after this call. **/
    void clear();

    /** Get the number of faces in the mesh. **/
    int getNumFaces() const;
    /** Get the number of vertices in the mesh. **/
    int getNumVertices() const;

    /** Get the position of a vertex in the mesh.
    @param[in]  vertex  The index of the vertex (as returned by addVertex()).
    @return The position of the specified vertex, measured and expressed in
    the mesh local frame. **/
    const Vec3& getVertexPosition(int vertex) const;
    /** Get the number of vertices that make up a particular face.
    @param[in]  face    The index of the face (as returned by addFace()). **/
    int getNumVerticesForFace(int face) const;
    /** Get the index of one of the vertices of a face.
    @param[in]  face    The index of the face (as returned by addFace()).
    @param[in]  vertex  The index of the vertex within the face (from 0, 1, or 2
                        for a triangular face, etc.) These are ordered the same
                        way as when the face was defined.
    @return The index of the specified vertex. **/
    int getFaceVertex(int face, int vertex) const;

    /** Add a vertex to the mesh.
    @param[in]  position   The position of the vertex to add, measured and
                           expressed in the mesh local frame.
    @return The index of the newly added vertex. **/
    int addVertex(const Vec3& position);

    /** Add a face to the mesh. Note that the ordering of the vertices defines
    the outward normal for the face; they must be counterclockwise around the
    desired normal.

    @param[in]  vertices    Indices of the vertices which make up the new face,
                            in counterclockwise order with respect to the face
                            normal.
    @return The index of the newly added face. **/
    int addFace(const Array_<int>& vertices);

    /** Scale a mesh by multiplying every vertex by a fixed value. Note that
    this permanently modifies the vertex locations within the mesh. Since the
    vertices are measured in the mesh local frame, scaling will appear to
    occur around the mesh origin (that is, the origin will remain where it
    was while everything else changes.
    @param[in]  scale   The scale factor. Can be any value except zero.
    @return A reference to this now-scaled mesh object. **/
    PolygonalMesh& scaleMesh(Real scale);

    /** %Transform a mesh by applying the given Transform to every vertex,
    leaving the mesh permanently changed. This has the effect of replacing the
    mesh local frame M with a new frame A.
    @param[in]  X_AM   The transform giving the pose of the mesh local frame in
                       the new frame A. Every vertex v_M becomes v_A=X_AM*v_M.
    @return A reference to this now-transformed mesh object. **/
    PolygonalMesh&  transformMesh(const Transform& X_AM);

    /** Attempt to interpret the given file as a mesh file, with the format
    determined from the file name extension. If we recognize the extension
    we'll call one of the specialized methods below; see the descriptions for
    more information. Ignoring case, we recognize:
        - <tt>.obj </tt>: Wavefront OBJ file
        - <tt>.stl </tt>: 3D Systems Stereolithography file (ascii or binary)
        - <tt>.stla</tt>: ascii-only stl extension
        - <tt>.vtp </tt>: VTK PolyData file (we can only read the ascii version)

    @param[in]  pathname    The name of a mesh file with a recognized extension.
    **/
    void loadFile(const String& pathname);

    /** Load a Wavefront OBJ (.obj) file, adding the vertices and faces it
    contains to this mesh, and ignoring anything else in the file. The suffix
    for these files is typically ".obj" but we don't check here.
    @param[in]  pathname    The name of a .obj file. **/
    void loadObjFile(const String& pathname);

    /** Alternate signature for Wavefront OBJ format that takes an already-open
    istream rather than a pathname. This is useful for testing since it
    can be supplied by a stringstream rather than a file.
    @param[in,out]  file    An input stream from which to load the file
                            contents. **/
    void loadObjFile(std::istream& file);

    /** Load a VTK PolyData (.vtp) file, adding the vertices and faces it
    contains to this mesh and ignoring anything else in the file. The suffix
    for these files is typically ".vtp" but we don't check here.
    @param[in]  pathname    The name of a .vtp file. **/
    void loadVtpFile(const String& pathname);

    /** Load an STL file, adding the vertices and faces it contains to this
    mesh and ignoring anything else in the file. The file may be in ascii or
    binary format. If the suffix is ".stla" then it can only be ascii.
    Otherwise, including ".stl" or anything else, we'll examine the contents to
    determine which format is used. STL files include many repeated vertices;
    we will collapse any that coincide to within a small tolerance so that there
    is some hope of getting a connected surface.
    @param[in]  pathname    The name of a .stl or .stla file. **/
    void loadStlFile(const String& pathname);

private:
    explicit PolygonalMesh(PolygonalMeshImpl* impl) : HandleBase(impl) {}
    void initializeHandleIfEmpty();
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
