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


#include "simbody/internal/CollisionDetectionAlgorithm.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometryImpl.h"
#include <vector>

using std::map;
using std::pair;
using std::vector;

namespace SimTK {

static int registerStandardAlgorithms() {
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::HalfSpaceImpl::Type(), ContactGeometry::SphereImpl::Type(), new CollisionDetectionAlgorithm::HalfSpaceSphere());
    CollisionDetectionAlgorithm::registerAlgorithm(ContactGeometry::SphereImpl::Type(), ContactGeometry::SphereImpl::Type(), new CollisionDetectionAlgorithm::SphereSphere());
    return 1;
}

static int staticInitializer = registerStandardAlgorithms();

map<pair<int, int>, CollisionDetectionAlgorithm*>& CollisionDetectionAlgorithm::getAlgorithmMap() {
    static map<pair<int, int>, CollisionDetectionAlgorithm*> algorithmMap;
    return algorithmMap;
}

void CollisionDetectionAlgorithm::registerAlgorithm(const std::string& type1, const std::string& type2, CollisionDetectionAlgorithm* algorithm) {
    int typeIndex1 = ContactGeometryImpl::getIndexForType(type1);
    int typeIndex2 = ContactGeometryImpl::getIndexForType(type2);
    getAlgorithmMap()[pair<int, int>(typeIndex1, typeIndex2)] = algorithm;
}

CollisionDetectionAlgorithm* CollisionDetectionAlgorithm::getAlgorithm(int typeIndex1, int typeIndex2) {
    map<pair<int, int>, CollisionDetectionAlgorithm*> algorithmMap = getAlgorithmMap();
    map<pair<int, int>, CollisionDetectionAlgorithm*>::iterator iter = algorithmMap.find(pair<int, int>(typeIndex1, typeIndex2));
    if (iter == algorithmMap.end())
        return NULL;
    return iter->second;

}

void CollisionDetectionAlgorithm::HalfSpaceSphere::processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
        int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::SphereImpl& sphere = dynamic_cast<const ContactGeometry::SphereImpl&>(object2.getImpl());
    Vec3 location = (~transform1)*transform2.T(); // Location of the sphere in the half-space's coordinate frame
    Real r = sphere.getRadius();
    Real depth = r+location[0];
    if (depth < 0)
        return; // No intersection.
    Real contactRadius = std::sqrt(depth*r);
    Vec3 normal = transform1.R()*Vec3(-1, 0, 0);
    Vec3 contactLocation = transform1*Vec3(0.5*depth, location[1], location[2]);
    contacts.push_back(Contact(index1, index2, contactLocation, normal, contactRadius, depth));
}

void CollisionDetectionAlgorithm::SphereSphere::processObjects(int index1, const ContactGeometry object1, const Transform& transform1,
        int index2, const ContactGeometry object2, const Transform& transform2, std::vector<Contact>& contacts) const {
    const ContactGeometry::SphereImpl& sphere1 = dynamic_cast<const ContactGeometry::SphereImpl&>(object1.getImpl());
    const ContactGeometry::SphereImpl& sphere2 = dynamic_cast<const ContactGeometry::SphereImpl&>(object2.getImpl());
    Vec3 delta = transform2.T()-transform1.T();
    Real dist = delta.norm();
    if (dist == 0)
        return; // No sensible way to deal with this.
    Real r1 = sphere1.getRadius();
    Real r2 = sphere2.getRadius();
    Real depth = r1+r2-dist;
    if (depth > 0) {
        // They are overlapping.
        
        Real curvature = r1*r2/(r1+r2);
        Real contactRadius = std::sqrt(depth*curvature);
        Vec3 normal = delta/dist;
        Vec3 location = transform1.T()+(r1-0.5*depth)*normal;
        contacts.push_back(Contact(index1, index2, location, normal, contactRadius, depth));
    }
}

} // namespace SimTK

