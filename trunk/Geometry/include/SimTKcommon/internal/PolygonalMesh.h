#ifndef SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
#define SimTK_SimTKCOMMON_POLYGONAL_MESH_H_

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

#include "SimTKcommon/Simmatrix.h"

namespace SimTK {

class PolygonalMeshImpl;

/**
 * This class provides a description of a mesh made of polygonal faces.  Its primary purpose is for
 * loading geometry from files, which can then be used for visualization or collision detection.
 * For example, the following lines load a mesh from a Wavefront OBJ file, then create a DecorativeMesh
 * from it.
 *
 * <pre>
 * PolygonalMesh mesh;
 * std::ifstream file;
 * file.open("teapot.obj");
 * mesh.loadObjFile(file);
 * file.close();
 * DecorativeMesh decoration(mesh);
 * </pre>
 */
class SimTK_SimTKCOMMON_EXPORT PolygonalMesh {
public:
    /**
     * Create a PolygonalMesh, which initially contains no vertices or faces.
     */
    PolygonalMesh();
    PolygonalMesh(const PolygonalMesh& copy);
    PolygonalMesh& operator=(const PolygonalMesh& copy);
    ~PolygonalMesh();
    /**
     * Get the number of faces in the mesh.
     */
    int getNumFaces() const;
    /**
     * Get the number of vertices in the mesh.
     */
    int getNumVertices() const;
    /**
     * Get the position of a vertex in the mesh.
     *
     * @param vertex  the index of the vertex to get
     * @return the position of the specified vertex
     */
    const Vec3& getVertexPosition(int vertex) const;
    /**
     * Get the number of vertices that make up a face.
     *
     * @param face    the index of the face
     */
    int getNumVerticesForFace(int face) const;
    /**
     * Get the index of one of the vertices of a face.
     *
     * @param face    the index of the face
     * @param vertex  the index of the vertex within the face (from 0, 1, or 2 for a triangular face, etc.)
     * @return the index of the specified vertex
     */
    int getFaceVertex(int face, int vertex) const;
    /**
     * Add a vertex to the mesh.
     *
     * @param position    the position of the vertex to add
     * @return the index of the newly added vertex
     */
    int addVertex(const Vec3& position);
    /**
     * Add a face to the mesh.
     *
     * @param vertices    the indices of the vertices which make up the new face
     * @return the index of the newly added face
     */
    int addFace(const Array_<int>& vertices);
    /**
     * Scale a mesh by multiplying every vertex by a fix value.
     */
    void scaleMesh(Real scale);
    /**
     * Transform a mesh by applying a Transform to every vertex.
     */
    void transformMesh(const Transform& transform);
    /**
     * Load a Wavefront OBJ file, adding the vertices and faces it contains to this mesh.
     *
     * @param file    an input stream from which to load the file contents
     */
    void loadObjFile(std::istream& file);
    const PolygonalMeshImpl& getImpl() const;
    PolygonalMeshImpl& updImpl();
private:
    PolygonalMeshImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_POLYGONAL_MESH_H_
