#ifndef SimTK_SIMBODY_CONTACT_GEOMETRY_H_
#define SimTK_SIMBODY_CONTACT_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
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

#include "simbody/internal/common.h"
#include <vector>

namespace SimTK {

class ContactGeometryImpl;

/**
 * A ContactGeometry object describes the physical shape of a body.  It is used with GeneralContactSubsystem
 * for doing collision detection and contact modeling.
 */
class SimTK_SIMBODY_EXPORT ContactGeometry {
public:
    class HalfSpace;
    class Sphere;
    class TriangleMesh;
    class HalfSpaceImpl;
    class SphereImpl;
    class TriangleMeshImpl;
    ContactGeometry() : impl(0) {
    }
    ContactGeometry(const ContactGeometry& src);
    explicit ContactGeometry(ContactGeometryImpl* impl);
    virtual ~ContactGeometry();
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;
    ContactGeometry& ContactGeometry::operator=(const ContactGeometry& src);
    bool hasImpl() const {
        return impl != 0;
    }
    const ContactGeometryImpl& getImpl() const {
        assert(impl);
        return *impl;
    }
    ContactGeometryImpl& updImpl() {
        assert(impl);
        return *impl;
    }
    /**
     * Get a string which uniquely identifies the type of geometry this object represents.
     * Typically each subclass of ContactGeometry defines its own value.
     */
    const std::string& getType() const;
    /**
     * Get an integer which uniquely identifies the type of geometry this object represents.
     * A unique index is generated automatically for each unique type value as returned by getType().
     */
    int getTypeIndex() const;
protected:
    ContactGeometryImpl* impl;
};

/**
 * This ContactGeometry subclass represents an object that occupies the entire half-space x>0.
 * This is useful for representing walls and floors.
 */
class SimTK_SIMBODY_EXPORT ContactGeometry::HalfSpace : public ContactGeometry {
public:
    HalfSpace();
};

/**
 * This ContactGeometry subclass represents a sphere centered at the origin.
 */
class SimTK_SIMBODY_EXPORT ContactGeometry::Sphere : public ContactGeometry {
public:
    Sphere(Real radius);
    Real getRadius() const;
    void setRadius(Real radius);
    const SphereImpl& getImpl() const;
    SphereImpl& updImpl();
};

/**
 * This ContactGeometry subclass represents an arbitrary shape described by a mesh of triangular faces.
 * The mesh surface must satisfy the following requirements:
 *
 * <ol>
 * <li>It must be closed, so that any point can unambiguously be classified as either inside or outside.</li>
 * <li>It may not intersect itself anywhere, even at a single point.</li>
 * <li>It must be an oriented manifold.</li>
 * <li>The vertices for each face must be ordered counter-clockwise when viewed from the outside.  That is, if
 * v0, v1, and v2 are the locations of the three vertices for a face, the cross product (v1-v0)%(v2-v0) must
 * point outward.</li>
 * <li>The length of every edge must be non-zero.</li>
 * </ol>
 *
 * It is your responsibility to ensure that any mesh you create meets these requirements.  The constructor
 * will detect many incorrect meshes and signal them by throwing an exception, but it is not guaranteed to
 * detect all possible problems.  If a mesh fails to satisfy any of these requirements, the results of calculations
 * performed with it are undefined.  For example, collisions involving it might fail to be detected, or contact
 * forces on it might be calculated incorrectly.
 */
class SimTK_SIMBODY_EXPORT ContactGeometry::TriangleMesh : public ContactGeometry {
public:
    /**
     * Create a TriangleMesh.
     *
     * @param vertices     the positions of all vertices in the mesh
     * @param faceIndices  the indices of the vertices that make up each face.  The first three
     *                     elements are the vertices in the first face, the next three elements are
     *                     the vertices in the second face, etc.
     */
    TriangleMesh(const std::vector<Vec3>& vertices, const std::vector<int>& faceIndices);
    /**
     * Get the number of edges in the mesh.
     */
    int getNumEdges() const;
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
     * @param index  the index of the vertex to get
     * @return the position of the specified vertex
     */
    const Vec3& getVertexPosition(int index) const;
    /**
     * Get the index of one of the edges of a face.  Edge 0 connects vertices 0 and 1.
     * Edge 1 connects vertices 1 and 2.  Edge 2 connects vertices 0 and 2.
     *
     * @param face    the index of the face
     * @param edge    the index of the edge within the face (0, 1, or 2)
     * @return the index of the specified edge
     */
    int getFaceEdge(int face, int edge) const;
    /**
     * Get the index of one of the vertices of a face.
     *
     * @param face    the index of the face
     * @param vertex  the index of the vertex within the face (0, 1, or 2)
     * @return the index of the specified vertex
     */
    int getFaceVertex(int face, int vertex) const;
    /**
     * Get the index of one of the faces shared by an edge
     *
     * @param edge    the index of the edge
     * @param face    the index of the face within the edge (0 or 1)
     * @return the index of the specified face
     */
    int getEdgeFace(int edge, int face) const;
    /**
     * Get the index of one of the vertices shared by an edge
     *
     * @param edge    the index of the edge
     * @param vertex  the index of the vertex within the edge (0 or 1)
     * @return the index of the specified vertex
     */
    int getEdgeVertex(int edge, int vertex) const;
    /**
     * Find all edges that intersect a vertex.
     *
     * @param vertex  the index of the vertex
     * @param edges   the indices of all edges intersecting the vertex will be added to this
     */
    void findVertexEdges(int vertex, std::vector<int>& edges) const;
    const TriangleMeshImpl& getImpl() const;
    TriangleMeshImpl& updImpl();
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_GEOMETRY_H_
