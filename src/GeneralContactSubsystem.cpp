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


#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometryImpl.h"
#include "simbody/internal/CollisionDetectionAlgorithm.h"
#include <algorithm>
#include <map>
#include <vector>

using std::map;
using std::pair;
using std::sort;
using std::vector;

namespace SimTK {

// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, const vector<vector<Contact> >&) {
    assert(false);
    return o;
}

class ContactSet {
public:
    vector<MobilizedBody> bodies;
    vector<ContactGeometry> geometry;
    vector<Transform> transforms;
    mutable vector<Vec3> sphereCenters;
    mutable vector<Real> sphereRadii;
};

class ContactBodyExtent {
public:
    ContactBodyExtent(Real start, Real end, int index) : start(start), end(end), index(index) {
    }
    ContactBodyExtent() {
    }
    bool operator<(const ContactBodyExtent& e) const {
        return start < e.start;
    }
    Real start, end;
    int index;
};

class GeneralContactSubsystemImpl : public Subsystem::Guts {
public:
    GeneralContactSubsystemImpl() : contactsCacheIndex(-1), contactsValidCacheIndex(-1)
    {
    }

    GeneralContactSubsystemImpl* cloneImpl() const {
        return new GeneralContactSubsystemImpl(*this);
    }
    
    ContactSetIndex createContactSet() {
        invalidateSubsystemTopologyCache();
        int size = sets.size();
        sets.resize(size+1);
        return ContactSetIndex(size);
    }
    
    int getNumContactSets() const {
        return sets.size();
    }
    
    void addBody(ContactSetIndex setIndex, const MobilizedBody& body, const ContactGeometry& geom, Transform& transform) {
        assert(setIndex >= 0 && setIndex < sets.size());
        invalidateSubsystemTopologyCache();
        ContactSet& set = sets[setIndex];
        set.bodies.push_back(body);
        set.geometry.push_back(geom);
        set.transforms.push_back(transform);
    }

    int getNumBodies(ContactSetIndex set) const {
        return sets[set].bodies.size();
    }

    const MobilizedBody& getBody(ContactSetIndex set, int index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < sets[set].bodies.size());
        return sets[set].bodies[index];
    }

    const ContactGeometry& getBodyGeometry(ContactSetIndex set, int index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < sets[set].geometry.size());
        return sets[set].geometry[index];
    }

    ContactGeometry& updBodyGeometry(ContactSetIndex set, int index) {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < sets[set].geometry.size());
        invalidateSubsystemTopologyCache();
        return sets[set].geometry[index];
    }

    const Transform& getBodyTransform(ContactSetIndex set, int index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < sets[set].transforms.size());
        return sets[set].transforms[index];
    }

    Transform& updBodyTransform(ContactSetIndex set, int index) {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < sets[set].transforms.size());
        invalidateSubsystemTopologyCache();
        return sets[set].transforms[index];
    }

    const vector<Contact>& getContacts(const State& state, ContactSetIndex set) const {
        assert(set >= 0 && set < sets.size());
        SimTK_STAGECHECK_GE_ALWAYS(state.getSubsystemStage(getMySubsystemIndex()), Stage::Dynamics, "GeneralContactSubsystemImpl::getContacts()");
        vector<vector<Contact> >& contacts = Value<vector<vector<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
        return contacts[set];
    }
    
    int realizeSubsystemTopologyImpl(State& state) const {
        contactsCacheIndex = state.allocateCacheEntry(getMySubsystemIndex(), Stage::Dynamics, new Value<vector<vector<Contact> > >());
        contactsValidCacheIndex = state.allocateCacheEntry(getMySubsystemIndex(), Stage::Position, new Value<bool>());
        for (int i = 0; i < sets.size(); ++i) {
            const ContactSet& set = sets[i];
            int numBodies = set.bodies.size();
            set.sphereCenters.resize(numBodies);
            set.sphereRadii.resize(numBodies);
            for (int j = 0; j < numBodies; j++) {
                set.geometry[j].getBoundingSphere(set.sphereCenters[j], set.sphereRadii[j]);
                set.sphereCenters[j] = set.transforms[j]*set.sphereCenters[j];
            }
        }
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& state) const {
        Value<bool>::downcast(state.updCacheEntry(getMySubsystemIndex(), contactsValidCacheIndex)).upd() = false;
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& state) const {
        bool& contactsValid = Value<bool>::downcast(state.updCacheEntry(getMySubsystemIndex(), contactsValidCacheIndex)).upd();
        if (contactsValid)
            return 0;
        vector<vector<Contact> >& contacts = Value<vector<vector<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
        int numSets = getNumContactSets();
        contacts.resize(numSets);
        
        // Loop over all contact sets.
        
        for (int setIndex = 0; setIndex < numSets; setIndex++) {
            contacts[setIndex].clear();
            const ContactSet& set = sets[setIndex];
            int numBodies = set.bodies.size();
            
            // Perform a sweep-and-prune on a single axis to identify potential contacts.  First, find which
            // axis has the most variation in body locations.  That is the axis we will use.
            
            Vector_<Vec3> centers(numBodies);
            for (int i = 0; i < numBodies; i++)
                centers[i] = set.bodies[i].getBodyTransform(state)*set.sphereCenters[i];
            Vec3 average = mean(centers);
            Vec3 var(0);
            for (int i = 0; i < numBodies; i++)
                var += abs(centers[i]-average);
            int axis = (var[0] > var[1] ? 0 : 1);
            if (var[2] > var[axis])
                axis = 2;
            
            // Find the extent of each body along the axis and sort them by starting location.
            
            vector<ContactBodyExtent> extents(numBodies);
            for (int i = 0; i < numBodies; i++)
                extents[i] = ContactBodyExtent(centers[i][axis]-set.sphereRadii[i], centers[i][axis]+set.sphereRadii[i], i);
            sort(extents.begin(), extents.end());
            
            // Now sweep along the axis, finding potential contacts.
            
            for (int i = 0; i < numBodies; i++) {
                const int index1 = extents[i].index;
                const Transform transform1 = set.bodies[index1].getBodyTransform(state)*set.transforms[index1];
                const ContactGeometry& geom1 = set.geometry[index1];
                const int typeIndex1 = geom1.getTypeIndex();
                for (int j = i+1; extents[j].start <= extents[i].end && j < numBodies; j++) {
                    // They overlap along this axis.  See if the bounding spheres overlap.
                    
                    const int index2 = extents[j].index;
                    const Real sumRadius = set.sphereRadii[index1]+set.sphereRadii[index2];
                    if ((centers[index1]-centers[index2]).normSqr() <= sumRadius*sumRadius) {
                        // Do a full collision detection.

                        const Transform transform2 = set.bodies[index2].getBodyTransform(state)*set.transforms[index2];
                        const ContactGeometry& geom2 = set.geometry[index2];
                        const int typeIndex2 = geom2.getTypeIndex();
                        CollisionDetectionAlgorithm* algorithm = CollisionDetectionAlgorithm::getAlgorithm(typeIndex1, typeIndex2);
                        if (algorithm == NULL) {
                            algorithm = CollisionDetectionAlgorithm::getAlgorithm(typeIndex2, typeIndex1);
                            if (algorithm == NULL)
                                continue; // No algorithm available for detecting collisions between these two objects.
                            algorithm->processObjects(index2, geom2, transform2, index1, geom1, transform1, contacts[setIndex]);
                        }
                        else {
                            algorithm->processObjects(index1, geom1, transform1, index2, geom2, transform2, contacts[setIndex]);
                        }
                    }
                }
            }
        }
        contactsValid = true;
        return 0;
    }

    SimTK_DOWNCAST(GeneralContactSubsystemImpl, Subsystem::Guts);

private:
    mutable int contactsCacheIndex;
    mutable int contactsValidCacheIndex;
    vector<ContactSet> sets;
};

ContactSetIndex GeneralContactSubsystem::createContactSet() {
    return updImpl().createContactSet();
}

int GeneralContactSubsystem::getNumContactSets() const {
    return getImpl().getNumContactSets();
}

void GeneralContactSubsystem::addBody(ContactSetIndex index, const MobilizedBody& body, const ContactGeometry& geom, Transform transform) {
    return updImpl().addBody(index, body, geom, transform);
}

int GeneralContactSubsystem::getNumBodies(ContactSetIndex set) const {
    return getImpl().getNumBodies(set);
}

const MobilizedBody& GeneralContactSubsystem::getBody(ContactSetIndex set, int index) const {
    return getImpl().getBody(set, index);
}

const ContactGeometry& GeneralContactSubsystem::getBodyGeometry(ContactSetIndex set, int index) const {
    return getImpl().getBodyGeometry(set, index);
}

ContactGeometry& GeneralContactSubsystem::updBodyGeometry(ContactSetIndex set, int index) {
    return updImpl().updBodyGeometry(set, index);
}

const Transform& GeneralContactSubsystem::getBodyTransform(ContactSetIndex set, int index) const {
    return getImpl().getBodyTransform(set, index);
}

Transform& GeneralContactSubsystem::updBodyTransform(ContactSetIndex set, int index) {
    return updImpl().updBodyTransform(set, index);
}

const vector<Contact>& GeneralContactSubsystem::getContacts(const State& state, ContactSetIndex set) const {
    return getImpl().getContacts(state, set);
}

bool GeneralContactSubsystem::isInstanceOf(const Subsystem& s) {
    return GeneralContactSubsystemImpl::isA(s.getSubsystemGuts());
}
const GeneralContactSubsystem& GeneralContactSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const GeneralContactSubsystem&>(s);
}
GeneralContactSubsystem& GeneralContactSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<GeneralContactSubsystem&>(s);
}

const GeneralContactSubsystemImpl& GeneralContactSubsystem::getImpl() const {
    return dynamic_cast<const GeneralContactSubsystemImpl&>(getSubsystemGuts());
}
GeneralContactSubsystemImpl& GeneralContactSubsystem::updImpl() {
    return dynamic_cast<GeneralContactSubsystemImpl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
GeneralContactSubsystem::GeneralContactSubsystem() {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
}

GeneralContactSubsystem::GeneralContactSubsystem(MultibodySystem& mbs) {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
    mbs.setContactSubsystem(*this); // steal ownership
}

} // namespace SimTK

