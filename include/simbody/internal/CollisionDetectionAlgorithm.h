#ifndef SimTK_SIMBODY_COLLISION_DETECTION_ALGORITHM_H_
#define SimTK_SIMBODY_COLLISION_DETECTION_ALGORITHM_H_

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
#include <map>

namespace SimTK {

class ContactGeometry;
class Contact;

/**
 * A CollisionDetectionAlgorithm implements an algorithm for detecting overlaps between pairs of
 * ContactGeometry objects, and creating Contact objects based on them.  This class is used internally
 * by GeneralContactSubsystem, and there usually is no reason to access it directly.  The exception is
 * if you are defining a new ContactGeometry subclass.  In that case, you will also need to define one
 * or more CollisionDetectionAlgorithms to detect collisions with your new geometry type, then register it
 * calling registerAlgorithm().
 */
class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm {
public:
    class HalfSpaceSphere;
    class SphereSphere;
    /**
     * Identify contacts between a pair of bodies.
     *
     * @param index1     the index of the first body within its contact set
     * @param object1    the ContactGeometry for the first body
     * @param transform1 the location and orientation of the first body in the ground frame
     * @param index2     the index of the second body within its contact set
     * @param object2    the ContactGeometry for the second body
     * @param transform2 the location and orientation of the second body in the ground frame
     * @param contacts   if the bodies overlap, a Contact should be added to this for each distinct
     *                   contact between them.  (Multiple contacts may exist if one of the bodies
     *                   is concave.)
     */
    virtual void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const = 0;
    /**
     * Register a CollisionDetectionAlgorithm to be used for identifying contacts between bodies of two specific types.
     *
     * @param type1      the type identifier for the ContactGeometry subclass the algorithm expects as the first body
     * @param type2      the type identifier for the ContactGeometry subclass the algorithm expects as the second body
     * @param algorithm  the algorithm to use for bodies of the specified types
     */
    static void registerAlgorithm(const std::string& type1, const std::string& type2, CollisionDetectionAlgorithm* algorithm);
    /**
     * Get the CollisionDetectionAlgorithm to use for identifying contacts between bodies of two specific types.
     *
     * @param typeIndex1     the type index of the first body's ContactGeometry
     * @param typeIndex2     the type index of the second body's ContactGeometry
     * @return the CollisionDetectionAlgorithm to use, or NULL if no suitable algorithm has been registered
     */
    static CollisionDetectionAlgorithm* getAlgorithm(int typeIndex1, int typeIndex2);
private:
    static std::map<std::pair<int, int>, CollisionDetectionAlgorithm*>& getAlgorithmMap();
};

/**
 * This algorithm detects contacts between a ContactGeometry::HalfSpace and a ContactGeometry::Sphere.
 */
class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm::HalfSpaceSphere : public CollisionDetectionAlgorithm {
public:
    void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const;
};

/**
 * This algorithm detects contacts between two ContactGeometry::Sphere objects.
 */
class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm::SphereSphere : public CollisionDetectionAlgorithm {
public:
    void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_COLLISION_DETECTION_ALGORITHM_H_
