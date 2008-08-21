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
#include "simbody/internal/ForceSubsystemGuts.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometryImpl.h"
#include "simbody/internal/CollisionDetectionAlgorithm.h"
#include <map>
#include <vector>

using std::map;
using std::pair;
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
};
    
class GeneralContactSubsystemImpl : public ForceSubsystemRep {
public:
    GeneralContactSubsystemImpl() : ForceSubsystemRep("GeneralContactSubsystem", "1.0"), contactsCacheIndex(-1), energyCacheIndex(-1)
    {
    }

    GeneralContactSubsystemImpl* cloneImpl() const {
        return new GeneralContactSubsystemImpl(*this);
    }
    
    ContactSetIndex createContactSet() {
        int size = sets.size();
        sets.resize(size+1);
        return ContactSetIndex(size);
    }
    
    int getNumContactSets() const {
        return sets.size();
    }
    
    void addBody(ContactSetIndex index, MobilizedBody& body, ContactGeometry geom, Transform& transform) {
        ContactSet& set = sets[index];
        set.bodies.push_back(body);
        set.geometry.push_back(geom);
        set.transforms.push_back(transform);
    }

    const MobilizedBody& getBody(ContactSetIndex set, int index) const {
        return sets[set].bodies[index];
    }

    const ContactGeometry& getBodyGeometry(ContactSetIndex set, int index) const {
        return sets[set].geometry[index];
    }

    const Transform& getBodyTransform(ContactSetIndex set, int index) const {
        return sets[set].transforms[index];
    }

    const vector<Contact>& getContacts(const State& state, ContactSetIndex set) const {
        vector<vector<Contact> >& contacts = Value<vector<vector<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
        return contacts[set];
    }
    
    int realizeSubsystemTopologyImpl(State& s) const {
        contactsCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), Stage::Dynamics, new Value<vector<vector<Contact> > >());
        energyCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(), Stage::Dynamics, new Value<Real>());
        return 0;
    }

    int realizeSubsystemModelImpl(State& s) const {
        // Sorry, no choices available at the moment.
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& state) const {
        vector<vector<Contact> >& contacts = Value<vector<vector<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
        int numSets = getNumContactSets();
        contacts.resize(numSets);
        
        // Loop over all contact sets.
        
        for (int setIndex = 0; setIndex < numSets; setIndex++) {
            contacts[setIndex].clear();
            const ContactSet& set = sets[setIndex];
            int numBodies = set.bodies.size();
            
            // Loop over all pairs of objects in this set.
            
            for (int i = 1; i < numBodies; i++) {
                Transform transform1 = set.bodies[i].getBodyTransform(state)*set.transforms[i];
                const ContactGeometry& geom1 = set.geometry[i];
                int typeIndex1 = geom1.getTypeIndex();
                for (int j = 0; j < i; j++) {
                    Transform transform2 = set.bodies[j].getBodyTransform(state)*set.transforms[j];
                    const ContactGeometry& geom2 = set.geometry[j];
                    int typeIndex2 = geom2.getTypeIndex();
                    CollisionDetectionAlgorithm* algorithm = CollisionDetectionAlgorithm::getAlgorithm(typeIndex1, typeIndex2);
                    if (algorithm == NULL) {
                        algorithm = CollisionDetectionAlgorithm::getAlgorithm(typeIndex2, typeIndex1);
                        if (algorithm == NULL)
                            continue; // No algorithm available for detecting collisions between these two objects.
                        algorithm->processObjects(j, geom2, transform2, i, geom1, transform1, contacts[setIndex]);
                    }
                    else {
                        algorithm->processObjects(i, geom1, transform1, j, geom2, transform2, contacts[setIndex]);
                    }
                }
            }
        }
        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    Real calcPotentialEnergy(const State& state) const {
        return Value<Real>::downcast(state.getCacheEntry(getMySubsystemIndex(), energyCacheIndex)).get();
    }

private:
    mutable int contactsCacheIndex;
    mutable int energyCacheIndex;
    vector<ContactSet> sets;
};

ContactSetIndex GeneralContactSubsystem::createContactSet() {
    return updImpl().createContactSet();
}

int GeneralContactSubsystem::getNumContactSets() const {
    return getImpl().getNumContactSets();
}

void GeneralContactSubsystem::addBody(ContactSetIndex index, MobilizedBody& body, ContactGeometry geom, Transform transform) {
    return updImpl().addBody(index, body, geom, transform);
}

const MobilizedBody& GeneralContactSubsystem::getBody(ContactSetIndex set, int index) const {
    return getImpl().getBody(set, index);
}

const ContactGeometry& GeneralContactSubsystem::getBodyGeometry(ContactSetIndex set, int index) const {
    return getImpl().getBodyGeometry(set, index);
}

const Transform& GeneralContactSubsystem::getBodyTransform(ContactSetIndex set, int index) const {
    return getImpl().getBodyTransform(set, index);
}

const vector<Contact>& GeneralContactSubsystem::getContacts(const State& state, ContactSetIndex set) const {
    return getImpl().getContacts(state, set);
}

bool GeneralContactSubsystem::isInstanceOf(const ForceSubsystem& s) {
    return GeneralContactSubsystemImpl::isA(s.getRep());
}
const GeneralContactSubsystem& GeneralContactSubsystem::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const GeneralContactSubsystem&>(s);
}
GeneralContactSubsystem& GeneralContactSubsystem::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<GeneralContactSubsystem&>(s);
}

const GeneralContactSubsystemImpl& GeneralContactSubsystem::getImpl() const {
    return dynamic_cast<const GeneralContactSubsystemImpl&>(ForceSubsystem::getRep());
}
GeneralContactSubsystemImpl& GeneralContactSubsystem::updImpl() {
    return dynamic_cast<GeneralContactSubsystemImpl&>(ForceSubsystem::updRep());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
GeneralContactSubsystem::GeneralContactSubsystem() : ForceSubsystem() {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
}

GeneralContactSubsystem::GeneralContactSubsystem(MultibodySystem& mbs) : ForceSubsystem() {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
    mbs.addForceSubsystem(*this); // steal ownership
}

} // namespace SimTK

