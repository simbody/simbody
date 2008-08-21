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

class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm {
public:
    class HalfSpaceSphere;
    class SphereSphere;
    virtual void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const = 0;
    static void registerAlgorithm(const std::string& type1, const std::string& type2, CollisionDetectionAlgorithm* algorithm);
    static CollisionDetectionAlgorithm* getAlgorithm(int typeIndex1, int typeIndex2);
private:
    static std::map<std::pair<int, int>, CollisionDetectionAlgorithm*>& getAlgorithmMap();
};

class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm::HalfSpaceSphere : public CollisionDetectionAlgorithm {
public:
    void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const;
};

class SimTK_SIMBODY_EXPORT CollisionDetectionAlgorithm::SphereSphere : public CollisionDetectionAlgorithm {
public:
    void processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
            int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_COLLISION_DETECTION_ALGORITHM_H_
