/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
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

#include "SimTKmath.h"

#include "simbody/internal/common.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <algorithm>

namespace SimTK {

// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, const Array_<Array_<Contact> >&) {
    assert(false);
    return o;
}

class ContactSet {
public:
    Array_<MobilizedBody,ContactSurfaceIndex>   bodies;
    Array_<ContactGeometry,ContactSurfaceIndex> geometry;
    Array_<Transform,ContactSurfaceIndex>       transforms;
    mutable Array_<Vec3,ContactSurfaceIndex>    sphereCenters;
    mutable Array_<Real,ContactSurfaceIndex>    sphereRadii;
};

class ContactBodyExtent {
public:
    ContactBodyExtent(Real start, Real end, ContactSurfaceIndex index)
    :   start(start), end(end), index(index) {}
    ContactBodyExtent() {}
    bool operator<(const ContactBodyExtent& e) const
    {   return start < e.start; }
    Real start, end;
    ContactSurfaceIndex index;
};


//==============================================================================
//                      GENERAL CONTACT SUBSYSTEM IMPL
//==============================================================================
class GeneralContactSubsystemImpl : public Subsystem::Guts {
public:
    GeneralContactSubsystemImpl() {}

    GeneralContactSubsystemImpl* cloneImpl() const override {
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

    const MobilizedBody& getBody(ContactSetIndex set, ContactSurfaceIndex index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < (int) sets[set].bodies.size());
        return sets[set].bodies[index];
    }

    const ContactGeometry& getBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < (int) sets[set].geometry.size());
        return sets[set].geometry[index];
    }

    ContactGeometry& updBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index) {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < (int) sets[set].geometry.size());
        invalidateSubsystemTopologyCache();
        return sets[set].geometry[index];
    }

    const Transform& getBodyTransform(ContactSetIndex set, ContactSurfaceIndex index) const {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < (int) sets[set].transforms.size());
        return sets[set].transforms[index];
    }

    Transform& updBodyTransform(ContactSetIndex set, ContactSurfaceIndex index) {
        assert(set >= 0 && set < sets.size());
        assert(index >= 0 && index < (int) sets[set].transforms.size());
        invalidateSubsystemTopologyCache();
        return sets[set].transforms[index];
    }

    const Array_<Contact>& getContacts(const State& state, ContactSetIndex set) const {
        assert(set >= 0 && set < sets.size());
        SimTK_STAGECHECK_GE_ALWAYS(state.getSubsystemStage(getMySubsystemIndex()), Stage::Dynamics, "GeneralContactSubsystemImpl::getContacts()");
        Array_<Array_<Contact> >& contacts = Value<Array_<Array_<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
        return contacts[set];
    }

    int realizeSubsystemTopologyImpl(State& state) const override {
        contactsCacheIndex = state.allocateCacheEntry(getMySubsystemIndex(), Stage::Dynamics, new Value<Array_<Array_<Contact> > >());
        contactsValidCacheIndex = state.allocateCacheEntry(getMySubsystemIndex(), Stage::Position, new Value<bool>());
        for (int i = 0; i < (int) sets.size(); ++i) {
            const ContactSet& set = sets[i];
            int numBodies = set.bodies.size();
            set.sphereCenters.resize(numBodies);
            set.sphereRadii.resize(numBodies);
            for (ContactSurfaceIndex j(0); j < numBodies; j++) {
                set.geometry[j].getBoundingSphere(set.sphereCenters[j], set.sphereRadii[j]);
                set.sphereCenters[j] = set.transforms[j]*set.sphereCenters[j];
            }
        }
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& state) const override {
        Value<bool>::downcast(state.updCacheEntry(getMySubsystemIndex(), contactsValidCacheIndex)).upd() = false;
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& state) const override {
        bool& contactsValid = Value<bool>::downcast(state.updCacheEntry(getMySubsystemIndex(), contactsValidCacheIndex)).upd();
        if (contactsValid)
            return 0;
        Array_<Array_<Contact> >& contacts = Value<Array_<Array_<Contact> > >::downcast(updCacheEntry(state, contactsCacheIndex)).upd();
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
            for (ContactSurfaceIndex i(0); i < numBodies; i++)
                centers[i] = set.bodies[i].getBodyTransform(state)*set.sphereCenters[i];
            Vec3 average = mean(centers);
            Vec3 var(0);
            for (int i = 0; i < numBodies; i++)
                var += abs(centers[i]-average);
            int axis = (var[0] > var[1] ? 0 : 1);
            if (var[2] > var[axis])
                axis = 2;

            // Find the extent of each body along the axis and sort them by starting location.

            Array_<ContactBodyExtent> extents(numBodies);
            for (ContactSurfaceIndex i(0); i < numBodies; i++)
                extents[i] = ContactBodyExtent(centers[i][axis]-set.sphereRadii[i], centers[i][axis]+set.sphereRadii[i], i);
            std::sort(extents.begin(), extents.end());

            // Now sweep along the axis, finding potential contacts.

            for (int i = 0; i < numBodies; i++) {
                const ContactSurfaceIndex index1 = extents[i].index;
                const Transform transform1 = set.bodies[index1].getBodyTransform(state)*set.transforms[index1];
                const ContactGeometry& geom1 = set.geometry[index1];
                const ContactGeometryTypeId typeId1 = geom1.getTypeId();
                for (int j = i+1; j < numBodies && extents[j].start <= extents[i].end; j++) {
                    // They overlap along this axis.  See if the bounding spheres overlap.

                    const ContactSurfaceIndex index2 = extents[j].index;
                    const Real sumRadius = set.sphereRadii[index1]+set.sphereRadii[index2];
                    if ((centers[index1]-centers[index2]).normSqr() <= sumRadius*sumRadius) {
                        // Do a full collision detection.

                        const Transform transform2 = set.bodies[index2].getBodyTransform(state)*set.transforms[index2];
                        const ContactGeometry& geom2 = set.geometry[index2];
                        const ContactGeometryTypeId typeId2 = geom2.getTypeId();
                        CollisionDetectionAlgorithm* algorithm =
                            CollisionDetectionAlgorithm::getAlgorithm
                                                            (typeId1, typeId2);
                        if (algorithm == NULL) {
                            algorithm = CollisionDetectionAlgorithm::
                                                getAlgorithm(typeId2, typeId1);
                            if (algorithm == NULL)
                                continue; // No algorithm available for detecting collisions between these two objects.
                            algorithm->processObjects(index2, geom2, transform2,
                                                      index1, geom1, transform1,
                                                      contacts[setIndex]);
                        }
                        else {
                            algorithm->processObjects(index1, geom1, transform1,
                                                      index2, geom2, transform2,
                                                      contacts[setIndex]);
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
    Array_<ContactSet>      sets;

    mutable CacheEntryIndex contactsCacheIndex;
    mutable CacheEntryIndex contactsValidCacheIndex;
};



//==============================================================================
//                         GENERAL CONTACT SUBSYSTEM
//==============================================================================
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

const MobilizedBody& GeneralContactSubsystem::getBody(ContactSetIndex set, ContactSurfaceIndex index) const {
    return getImpl().getBody(set, index);
}

const ContactGeometry& GeneralContactSubsystem::getBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index) const {
    return getImpl().getBodyGeometry(set, index);
}

ContactGeometry& GeneralContactSubsystem::updBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index) {
    return updImpl().updBodyGeometry(set, index);
}

const Transform& GeneralContactSubsystem::getBodyTransform(ContactSetIndex set, ContactSurfaceIndex index) const {
    return getImpl().getBodyTransform(set, index);
}

Transform& GeneralContactSubsystem::updBodyTransform(ContactSetIndex set, ContactSurfaceIndex index) {
    return updImpl().updBodyTransform(set, index);
}

const Array_<Contact>& GeneralContactSubsystem::getContacts(const State& state, ContactSetIndex set) const {
    return getImpl().getContacts(state, set);
}

bool GeneralContactSubsystem::isInstanceOf(const Subsystem& s) {
    return GeneralContactSubsystemImpl::isA(s.getSubsystemGuts());
}
const GeneralContactSubsystem& GeneralContactSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<const GeneralContactSubsystem&>(s);
}
GeneralContactSubsystem& GeneralContactSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<GeneralContactSubsystem&>(s);
}

const GeneralContactSubsystemImpl& GeneralContactSubsystem::getImpl() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const GeneralContactSubsystemImpl&>(getSubsystemGuts());
}
GeneralContactSubsystemImpl& GeneralContactSubsystem::updImpl() {
    return SimTK_DYNAMIC_CAST_DEBUG<GeneralContactSubsystemImpl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::Array_'s, which require a default constructor to be available.
GeneralContactSubsystem::GeneralContactSubsystem() {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
}

GeneralContactSubsystem::GeneralContactSubsystem(MultibodySystem& mbs) {
    adoptSubsystemGuts(new GeneralContactSubsystemImpl());
    mbs.setContactSubsystem(*this); // steal ownership
}

} // namespace SimTK

