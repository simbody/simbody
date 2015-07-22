#ifndef SimTK_SIMMATH_COLLISION_DETECTION_ALGORITHM_H_
#define SimTK_SIMMATH_COLLISION_DETECTION_ALGORITHM_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
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
#include "simmath/internal/common.h"
#include "simmath/internal/ContactGeometry.h"
#include "simmath/internal/Contact.h"

#include <map>

namespace SimTK {

/** A CollisionDetectionAlgorithm implements an algorithm for detecting overlaps
between pairs of ContactGeometry objects, and creating Contact objects based on
them. This class is used internally by GeneralContactSubsystem, and there
usually is no reason to access it directly. The exception is if you are
defining a new ContactGeometry subclass. In that case, you will also need to
define one or more CollisionDetectionAlgorithms to detect collisions with your
new geometry type, then register it calling registerAlgorithm(). **/
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm {
public:
    class HalfSpaceSphere;
    class SphereSphere;
    class HalfSpaceEllipsoid;
    class HalfSpaceTriangleMesh;
    class SphereTriangleMesh;
    class TriangleMeshTriangleMesh;
    class ConvexConvex;
    virtual ~CollisionDetectionAlgorithm() {}
    /**
     * Identify contacts between a pair of bodies.
     *
     * @param index1     the index of the first body within its contact set
     * @param object1    the ContactGeometry for the first body
     * @param transform1 the location and orientation of the first body in the
     *                   ground frame
     * @param index2     the index of the second body within its contact set
     * @param object2    the ContactGeometry for the second body
     * @param transform2 the location and orientation of the second body in the
     *                   ground frame
     * @param contacts   if the bodies overlap, a Contact should be added to
     *                   this for each distinct contact between them. (Multiple
     *                   contacts may exist if one of the bodies is concave.)
     */
    virtual void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const = 0;
    /**
     * Register a CollisionDetectionAlgorithm to be used for identifying
     * contacts between bodies of two specific types.
     *
     * @param type1      the type identifier for the ContactGeometry subclass
     *                   the algorithm expects as the first body
     * @param type2      the type identifier for the ContactGeometry subclass
     *                   the algorithm expects as the second body
     * @param algorithm  the algorithm to use for bodies of the specified types
     */
    static void registerAlgorithm(ContactGeometryTypeId type1,
                                  ContactGeometryTypeId type2,
                                  CollisionDetectionAlgorithm* algorithm);
    /**
     * Get the CollisionDetectionAlgorithm to use for identifying contacts
     * between bodies of two specific types.
     *
     * @param type1     the type id of the first body's ContactGeometry
     * @param type2     the type id of the second body's ContactGeometry
     * @return the CollisionDetectionAlgorithm to use, or NULL if no suitable
     *         algorithm has been registered
     */
    static CollisionDetectionAlgorithm*
        getAlgorithm(ContactGeometryTypeId type1, ContactGeometryTypeId type2);
private:
    struct AlgorithmMap
    :   public std::map<std::pair<ContactGeometryTypeId, ContactGeometryTypeId>,
                        CollisionDetectionAlgorithm*>
    {
        ~AlgorithmMap();
    };

    static AlgorithmMap algorithmMap;
};

/**
 * This algorithm detects contacts between a ContactGeometry::HalfSpace and a
 * ContactGeometry::Sphere.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::HalfSpaceSphere
:   public CollisionDetectionAlgorithm {
public:
    virtual ~HalfSpaceSphere() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
};

/**
 * This algorithm detects contacts between a ContactGeometry::HalfSpace and a
 * ContactGeometry::Ellipsoid.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::HalfSpaceEllipsoid
:   public CollisionDetectionAlgorithm {
public:
    virtual ~HalfSpaceEllipsoid() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
};

/**
 * This algorithm detects contacts between two ContactGeometry::Sphere objects.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::SphereSphere
:   public CollisionDetectionAlgorithm {
public:
    virtual ~SphereSphere() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
};

/**
 * This algorithm detects contacts between a ContactGeometry::HalfSpace and a
 * ContactGeometry::TriangleMesh.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::HalfSpaceTriangleMesh
:   public CollisionDetectionAlgorithm {
public:
    virtual ~HalfSpaceTriangleMesh() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
private:
    void processBox(const ContactGeometry::TriangleMesh& mesh,
                    const ContactGeometry::TriangleMesh::OBBTreeNode& node,
                    const Transform& transform, const Vec3& axisDir,
                    Real xoffset, std::set<int>& insideFaces) const;
    void addAllTriangles(const ContactGeometry::TriangleMesh::OBBTreeNode& node,
                         std::set<int>& insideFaces) const;
};

/**
 * This algorithm detects contacts between a ContactGeometry::Sphere and a
 * ContactGeometry::TriangleMesh.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::SphereTriangleMesh
:   public CollisionDetectionAlgorithm {
public:
    virtual ~SphereTriangleMesh() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
private:
    void processBox(const Vec3& center, Real radius2,
                    const ContactGeometry::TriangleMesh& mesh,
                    const ContactGeometry::TriangleMesh::OBBTreeNode& node,
                    std::set<int>& insideFaces) const;
};

/**
 * This algorithm detects contacts between two ContactGeometry::TriangleMesh
 * objects.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::TriangleMeshTriangleMesh
:   public CollisionDetectionAlgorithm {
public:
    virtual ~TriangleMeshTriangleMesh() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
private:
    void processNodes(const ContactGeometry::TriangleMesh& mesh1,
                      const ContactGeometry::TriangleMesh& mesh2,
                      const ContactGeometry::TriangleMesh::OBBTreeNode& node1,
                      const ContactGeometry::TriangleMesh::OBBTreeNode& node2,
                      const OrientedBoundingBox& node2Bounds,
                      const Transform& transform, std::set<int>& triangles1,
                      std::set<int>& triangles2) const;
    void findInsideTriangles(const ContactGeometry::TriangleMesh& mesh,
                             const ContactGeometry::TriangleMesh& otherMesh,
                             const Transform& transform,
                             std::set<int>& triangles) const;
    void tagFaces(const ContactGeometry::TriangleMesh& mesh,
                  Array_<int>& faceType, std::set<int>& triangles,
                  int index, int depth) const;
    static const int OUTSIDE = -1;
    static const int UNKNOWN = 0;
    static const int BOUNDARY = 1;
    static const int INSIDE = 2;
};

/**
 * This algorithm detects contacts between two ContactGeometry::Convex objects.
 */
class SimTK_SIMMATH_EXPORT CollisionDetectionAlgorithm::ConvexConvex
:   public CollisionDetectionAlgorithm {
public:
    virtual ~ConvexConvex() {}
    void processObjects
       (ContactSurfaceIndex index1, const ContactGeometry& object1,
        const Transform& transform1,
        ContactSurfaceIndex index2, const ContactGeometry& object2,
        const Transform& transform2,
        Array_<Contact>& contacts) const;
private:
    static Vec3 computeSupport(const ContactGeometry& object1,
                               const ContactGeometry& object2,
                               const Transform& transform, UnitVec3 direction);
    static void addContact
       (ContactSurfaceIndex index1, ContactSurfaceIndex index2,
        const ContactGeometry& object1,
        const ContactGeometry& object2,
        const Transform& transform1, const Transform& transform2,
        const Transform& transform12,
        Vec3 point1, Vec3 point2, Array_<Contact>& contacts);
    static Vec6 computeErrorVector(const ContactGeometry& object1,
                                   const ContactGeometry& object2,
                                   Vec3 pos1, Vec3 pos2,
                                   const Transform& transform12);
    static Mat66 computeJacobian(const ContactGeometry& object1,
                                 const ContactGeometry& object2,
                                 Vec3 pos1, Vec3 pos2,
                                 const Transform& transform12);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_COLLISION_DETECTION_ALGORITHM_H_
