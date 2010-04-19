/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometryImpl.h"

#include "simbody/internal/ContactTrackerSubsystem.h"

#include <algorithm>
using std::pair; using std::make_pair;
#include <iostream>
using std::cout; using std::endl;


namespace SimTK {

// We keep a list of all the ContactSurfaces being tracked by this
// subsystem, and a separate list of "bubble wrap" spheres that are 
// used in broad phase determination of which contact surfaces need
// to be examined more closely by ContactTracker objects.

SimTK_DEFINE_UNIQUE_INDEX_TYPE(BubbleIndex);

//TODO: should these be organized by body rather than surface?
//      surfaces on ground certainly should!

struct Surface {
    const MobilizedBody*    mobod;
    const ContactSurface*   surface; // in its own frame S
    Transform               X_BS;    // surface's pose on body
};

struct Bubble {
    const Vec3& getCenter() const {return sphere.getSubVec<3>(0);}
    const Real& getRadius() const {return sphere[3];}
    Vec3& updCenter() {return sphere.updSubVec<3>(0);}
    Real& updRadius() {return sphere[3];}

    Vec4                sphere;  // x,y,z,r in body frame
    ContactSurfaceIndex surface; // associated surface
};
static std::ostream& operator<<(std::ostream& o, const Bubble& bb) {
    o << bb.sphere << "->" << bb.surface;
    return o;
}

// For spatial sorting of bubbles.
class BubbleExtent {
public:
    BubbleExtent(Real start, Real end, BubbleIndex index) 
    :   start(start), end(end), index(index) {}
    BubbleExtent() {}
    bool operator<(const BubbleExtent& e) const 
    {   return start < e.start; }
    Real        start, end; // span along chosen sort axis
    BubbleIndex index;
};
static std::ostream& operator<<(std::ostream& o, const BubbleExtent& bx) {
    o << bx.index << ":(" << bx.start << "," << bx.end << ")";
    return o;
}

typedef std::map< pair<ContactGeometryTypeId,ContactGeometryTypeId>,
                  pair<ContactTracker*,bool> > TrackerMap;

// This type maps a contact surface onto a set of all the higher-numbered
// contact surfaces it might be touching, and for each of those we keep
// a pointer to that pair's Contact object if it is currently being tracked
// (null if it is new). You must use the *lower* numbered surface as the key 
// when inserting a pair so that any given pair of surfaces appears just once.
// However, the surface order in the Contact object will be determined by
// the order required by the corresponding tracker.
// TODO: replace with a fast data structure
typedef std::map<ContactSurfaceIndex,const Contact*> ContactSurfaceSet;
typedef std::map<ContactSurfaceIndex,ContactSurfaceSet> PairMap;

static std::ostream& operator<<(std::ostream& o, const ContactSurfaceSet& css) {
    ContactSurfaceSet::const_iterator p = css.begin();
    o << "{";
    for (; p != css.end(); ++p) {
        if (p!=css.begin()) o << ", ";
        o << p->first << ":0x" << p->second;
    }
    return o << "}";
}
static std::ostream& operator<<(std::ostream& o, const PairMap& pm) {
    PairMap::const_iterator p = pm.begin();
    for (; p != pm.end(); ++p) {
        if (p != pm.begin()) o << endl;
        o << p->first << "->" << p->second;
    }
    return o;
}



//==============================================================================
//                       CONTACT TRACKER SUBSYSTEM IMPL
//==============================================================================
class ContactTrackerSubsystemImpl : public Subsystem::Guts {
public:
// Constructor registers a default set of Trackers to use with geometry
// we know about. These can be overridden later.
ContactTrackerSubsystemImpl() : m_defaultTracker(0) {
    adoptContactTracker(new ContactTracker::HalfSpaceSphere());
    adoptContactTracker(new ContactTracker::SphereSphere());
}

ContactTrackerSubsystemImpl* cloneImpl() const 
{   return new ContactTrackerSubsystemImpl(*this); }

void adoptContactTracker(ContactTracker* tracker) {
    assert(tracker);
    // This is the order required by the tracker's trackContact() method.
    const pair<ContactGeometryTypeId,ContactGeometryTypeId>&
        types = tracker->getContactGeometryTypeIds();
    ContactGeometryTypeId low=types.first, high=types.second;
    const bool mustReverse = (low > high);
    if (mustReverse) std::swap(low,high);
    m_contactTrackers[make_pair(low,high)] = make_pair(tracker, mustReverse);
}

// Return the MultibodySystem which owns this ContactTrackerSubsystem.
const MultibodySystem& getMultibodySystem() const 
{   return MultibodySystem::downcast(getSystem()); }

// Return the SimbodyMatterSubsystem from which this ContactTrackerSubsystem
// gets the bodies to track.
const SimbodyMatterSubsystem& getMatterSubsystem() const 
{   return getMultibodySystem().getMatterSubsystem(); }

// Get access to state variables and cache entries.

const ContactSnapshot& getPrevActiveContacts(const State& state) const {
    const ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
            (getDiscreteVariable(state, m_activeContactsIx));
    return contacts;
}
const ContactSnapshot& getPrevPredictedContacts(const State& state) const {
    const ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
            (getDiscreteVariable(state, m_predictedContactsIx));
    return contacts;
}
const ContactSnapshot& getNextActiveContacts(const State& state) const {
    const ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
        (getDiscreteVarUpdateValue(state, m_activeContactsIx));
    return contacts;
}
const ContactSnapshot& getNextPredictedContacts(const State& state) const {
    const ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
        (getDiscreteVarUpdateValue(state, m_predictedContactsIx));
    return contacts;
}
ContactSnapshot& updNextActiveContacts(const State& state) const {
    ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
        (updDiscreteVarUpdateValue(state, m_activeContactsIx));
    return contacts;
}
ContactSnapshot& updNextPredictedContacts(const State& state) const {
    ContactSnapshot& contacts = Value<ContactSnapshot>::downcast
        (updDiscreteVarUpdateValue(state, m_predictedContactsIx));
    return contacts;
}

// Run through all the bodies to find the contact surfaces, assigning each
// a unique ContactSurfaceIndex. Then for each surface, get its geometry
// and create a Bubble from each of its bubble wrap spheres; each of those
// gets a unique BubbleIndex that maps back to the associated surface.
int realizeSubsystemTopologyImpl(State& state) const {
    // Briefly allow writing into the Topology cache; after this the
    // Topology cache is const.
    ContactTrackerSubsystemImpl* wThis = 
        const_cast<ContactTrackerSubsystemImpl*>(this);
    wThis->m_activeContactsIx = allocateAutoUpdateDiscreteVariable
        (state, Stage::Dynamics, new Value<ContactSnapshot>(), 
         Stage::Position);      // update depends on positions
    wThis->m_predictedContactsIx = allocateAutoUpdateDiscreteVariable
        (state, Stage::Dynamics, new Value<ContactSnapshot>(), 
         Stage::Acceleration);  // update depends on accelerations

    const SimbodyMatterSubsystem& matter = getMatterSubsystem();

    const int numBodies = matter.getNumBodies();
    wThis->m_surfaces.clear();
    wThis->m_bubbles.clear();

    for (MobilizedBodyIndex mbx(0); mbx < numBodies; ++mbx) {
        const MobilizedBody& mobod = matter.getMobilizedBody(mbx);
        const Body&          body  = mobod.getBody();
        for (int i=0; i < body.getNumContactSurfaces(); ++i) {
            const ContactSurfaceIndex surfx(m_surfaces.size());
            wThis->m_surfaces.push_back(); // default construct
            Surface& surf = wThis->m_surfaces.back();
            surf.mobod   = &mobod; 
            surf.surface = &body.getContactSurface(i);
            surf.X_BS    =  body.getContactSurfaceTransform(i);
            const ContactGeometry& geo  = surf.surface->getShape();
            // for each bubble:
            wThis->m_bubbles.push_back(); // default construct
            Bubble& bubble = wThis->m_bubbles.back();
            Vec3 center_S; 
            geo.getBoundingSphere(center_S, bubble.updRadius());
            bubble.updCenter() = surf.X_BS*center_S; // in B
            bubble.surface = surfx;
        }
    }
    //cout << "Bubbles:" << m_bubbles << "\n";
    return 0;
}

int realizeSubsystemPositionImpl(const State& state) const {
    return 0;
}

int realizeSubsystemVelocityImpl(const State& state) const {
    return 0;
}

// Adds new pairs to the existing set, if not already present.
void addInBroadPhasePairs(const State& state, PairMap& pairs) const {
    const int numBubbles = getNumBubbles();
    
    // Perform a sweep-and-prune on a single axis to identify potential 
    // contacts. First, find which axis has the most variation in body 
    // locations. That is the axis we will use.
    // TODO: this one-axis method is not good enough in general
    
    Vector_<Vec3> centers(numBubbles);
    for (BubbleIndex bbx(0); bbx < numBubbles; ++bbx) {
        const Bubble&  bubb = m_bubbles[bbx];
        const Surface& surf = m_surfaces[bubb.surface];
        centers[bbx] = surf.mobod->getBodyTransform(state) 
                        * bubb.getCenter();
    }
    Vec3 average = mean(centers);
    Vec3 var(0);
    for (BubbleIndex bbx(0); bbx < numBubbles; ++bbx)
        var += abs(centers[bbx]-average);
    int axis = (var[0] > var[1] ? 0 : 1);
    if (var[2] > var[axis])
        axis = 2;
    
    // Find the extent of each bubble along the axis and sort them by 
    // starting location.
    Array_<BubbleExtent,int> extents(numBubbles);
    for (BubbleIndex bbx(0); bbx < numBubbles; ++bbx) {
        const Real radius = m_bubbles[bbx].getRadius();
        const Real center = centers[bbx][axis];
        extents[bbx] = BubbleExtent(center-radius, center+radius, bbx);
    }
    //cout << "Bubble extents:" << extents << "\n";

    // Expensive: O(n log n)
    std::sort(extents.begin(), extents.end());
    
    // Now sweep along the axis, finding potential contacts.
    
    for (int ex1=0; ex1 < numBubbles; ++ex1) {
        const BubbleExtent& extent1 = extents[ex1];
        const Bubble&       bubb1   = m_bubbles[extent1.index];
        const Real          radius1 = bubb1.getRadius();
        const Vec3&         center1 = centers[extent1.index];

        // Loop over just the overlapping bubbles.
        for (int ex2(ex1+1); ex2 < numBubbles; ++ex2) {
            const BubbleExtent& extent2 = extents[ex2];
            if (extent2.start > extent1.end)
                break;  // no more bubbles can overlap with extent1

            // These bubbles do overlap along this axis. See if they are
            // actually touching.
            const Bubble& bubb2   = m_bubbles[extent2.index];
            const Real    radius2 = bubb2.getRadius();
            const Vec3&   center2 = centers[extent2.index];
            if ((center1-center2).normSqr() > square(radius1+radius2))
                continue; // nope

            // The bubbles are touching. We'll add the corresponding surfaces
            // to the narrow-phase list unless there are relevant exclusions.
            const Surface& surf1 = m_surfaces[bubb1.surface];
            const Surface& surf2 = m_surfaces[bubb2.surface];
            // Ignore if on the same body.
            if (surf1.mobod == surf2.mobod) continue;
            assert(bubb1.surface != bubb2.surface); // duh!
            // Ignore if surfaces are in a common clique.
            if (surf1.surface->isInSameClique(*surf2.surface)) continue;
            // We'll need to do a narrow phase investigation of these two
            // surfaces; use the lower-numbered one as the index to avoid
            // duplicates.
            ContactSurfaceIndex low=bubb1.surface, high=bubb2.surface;
            if (low > high) std::swap(low,high);
            ContactSurfaceSet& surfSet = pairs[low];
            // Insert this pair with null Contact if the pair isn't already
            // in the PairMap.
            surfSet.insert(make_pair(high,(Contact*)0));
        }
    }
}

// Call this any time after positions are known, to ensure that the active
// contact set has been updated for those positions. We can use three
// sources of information to compute the update:
//   - The previously-known set of active contacts
//   - The previously-known set of predicted contacts
//   - The current contact surface positions
// Note that we cannot update the set of predicted contacts here because
// that requires velocity and acceleration information.
//
// Algorithm:
//   for all "interesting" surface pairs (surf1,surf2):
//      prev = getPrevActiveContact(surf1,surf2) (might be empty)
//      trackContact(prev, boundsNow?, next)
//      updNextActiveContact(surf1,surf2) = next  (might be empty)
//   "interesting" means
//      - previously active, or
//      - previously predicted, or
//      - broad phase position bounds intersect
void ensureActiveContactsUpdated(const State& state) const {
    if (isDiscreteVarUpdateValueRealized(state, m_activeContactsIx))
        return; // already done

    const ContactSnapshot& active     = getPrevActiveContacts(state);
    const ContactSnapshot& predicted  = getPrevPredictedContacts(state);
    ContactSnapshot&       nextActive = updNextActiveContacts(state);

    // TODO: Can we reuse heap space in this cache entry?
    nextActive.clear();

    PairMap interesting;
    for (int i=0; i < active.getNumContacts(); ++i) {
        const Contact& contact = active.getContact(i);
        ContactSurfaceIndex low=contact.getSurface1(), 
                            high=contact.getSurface2();
        if (low > high) std::swap(low,high);
        ContactSurfaceSet& others = interesting[low];
        assert(others.find(high)==others.end());
        others[high] = &contact;
    }
    for (int i=0; i < predicted.getNumContacts(); ++i) {
        const Contact& contact = predicted.getContact(i);
        ContactSurfaceIndex low=contact.getSurface1(), 
                            high=contact.getSurface2();
        if (low > high) std::swap(low,high);
        ContactSurfaceSet& others = interesting[low];
        assert(others.find(high)==others.end());
        others[high] = &contact;
    }
    // This will ignore pairs that we already inserted above; new ones
    // will be inserted with null Contact object pointers.
    addInBroadPhasePairs(state, interesting);
    //cout << "Interesting pairs:\n" << interesting << "\n";

    PairMap::const_iterator p = interesting.begin();
    for (; p != interesting.end(); ++p) {
        const ContactSurfaceIndex index1 = p->first;
        const Transform transform1 = 
            m_surfaces[index1].mobod->getBodyTransform(state)
                * m_surfaces[index1].X_BS;
        const ContactGeometry& geom1 = m_surfaces[index1].surface->getShape();
        const ContactGeometryTypeId typeId1 = geom1.getTypeId();

        const ContactSurfaceSet& others = p->second;
        ContactSurfaceSet::const_iterator q = others.begin();
        for (; q != others.end(); ++q) {
            const ContactSurfaceIndex index2 = q->first;
            const Transform transform2 = 
                m_surfaces[index2].mobod->getBodyTransform(state)
                    * m_surfaces[index2].X_BS;
            const ContactGeometry& geom2 = 
                m_surfaces[index2].surface->getShape();
            const ContactGeometryTypeId typeId2 = geom2.getTypeId();
            if (!hasContactTracker(typeId1,typeId2))
                continue; // No algorithm available for detecting collisions between these two objects.
            bool mustReverse;
            const ContactTracker& tracker = 
                getContactTracker(typeId1, typeId2, mustReverse);

            // Put the surfaces in the order required by the tracker.
            const ContactSurfaceIndex trackSurf1 = (mustReverse? index2:index1);
            const ContactSurfaceIndex trackSurf2 = (mustReverse? index1:index2);

            UntrackedContact untracked; // empty handle in case we need it
            const Contact* prev = q->second;
            if (prev && prev->getCondition() == Contact::Broken)
                prev = 0; // that contact expired
            if (!prev) { 
                untracked = UntrackedContact(trackSurf1, trackSurf2);
                prev = &untracked;
            }
            Contact next; // empty handle
            if (mustReverse)
                tracker.trackContact
                   (*prev, transform2,geom2, transform1,geom1, 0/*TODO*/, next);
            else
                tracker.trackContact
                   (*prev, transform1,geom1, transform2,geom2, 0/*TODO*/, next);

            if (!next.isEmpty()) {
                next.setSurfaces(trackSurf1,trackSurf2);
                next.setContactId(prev->getCondition()==Contact::Untracked
                                    ? Contact::createNewContactId()
                                    : prev->getContactId()); // persistent
                if (   prev->getCondition()==Contact::Untracked
                    || prev->getCondition()==Contact::Anticipated)
                    next.setCondition(Contact::NewContact);
                else { // was NewContact or Ongoing; now Ongoing
                    assert(prev->getCondition()==Contact::NewContact
                           || prev->getCondition()==Contact::Ongoing);
                    next.setCondition(Contact::Ongoing);
                }
                nextActive.adoptContact(next);
            }
        }
    }

    markDiscreteVarUpdateValueRealized(state, m_activeContactsIx);
}

// Call this any time after accelerations are known, to ensure that the
// predicted contact set has been updated for new velocities and accelerations.
// We can use three sources of information to compute the update:
//   - The current (updated) set of active contacts
//   - The previously-known set of impending contacts
//   - The current contact surface positions, velocities, and accelerations
// The active contact update cannot be modified here although we can
// initiate its computation if it hasn't been done yet; after that it is 
// frozen.
//
// Algorithm:
//   for all "interesting" surface pairs (surf1,surf2):
//      prev = getPrevPredictedContact(surf1,surf2) (might be empty)
//      predictContact(prev, projBoundsNow?, next)
//      updNextImpendingContact(surf1,surf2) = next  (might be empty)
//   "interesting" means not currently active and:
//      - previously impending, or
//      - fast(surf1)||fast(surf2) and 
//           fast object broad phase projected bounds intersect
void ensurePredictedContactsUpdated(const State& state) const {
    if (isDiscreteVarUpdateValueRealized(state, m_predictedContactsIx))
        return; // already done

    const ContactSnapshot& nextActive    = getNextActiveContacts(state);
    const ContactSnapshot& prevPredicted = getPrevPredictedContacts(state);
    ContactSnapshot& nextPredicted = updNextPredictedContacts(state);
    // TODO
    markDiscreteVarUpdateValueRealized(state, m_predictedContactsIx);
}

int realizeSubsystemDynamicsImpl(const State& state) const {
    ensureActiveContactsUpdated(state);
    return 0;
}

int realizeSubsystemAccelerationImpl(const State& state) const {
    ensurePredictedContactsUpdated(state);
    return 0;
}

bool hasContactTracker(ContactGeometryTypeId id1, 
                       ContactGeometryTypeId id2) const 
{   if (id1 > id2) std::swap(id1,id2); // (low,high) order for lookup
    return m_contactTrackers.find(make_pair(id1,id2))
                != m_contactTrackers.end();
}

const ContactTracker& 
getContactTracker(ContactGeometryTypeId id1, ContactGeometryTypeId id2,
                  bool& mustReverse) const 
{   const bool inputSwapped = id1 > id2;
    if (inputSwapped) std::swap(id1,id2); // (low,high) order for lookup
    TrackerMap::const_iterator p = m_contactTrackers.find(make_pair(id1,id2));
    const bool trackerSwapped = p->second.second;
    if (p != m_contactTrackers.end()) {
        mustReverse = (inputSwapped != trackerSwapped); // xor
        assert(p->second.first);
        return *p->second.first;
    }
    // Couldn't find a Tracker.

    SimTK_ERRCHK2_ALWAYS(m_defaultTracker,
        "ContactTrackerSubsystem::getContactTracker()",
        "There was no registered ContactTracker for surface geometry"
        " type ids (%d,%d) and there was no default tracker.",
        (int)id1, (int)id2);

    return *m_defaultTracker;
}

int getNumSurfaces() const {return m_surfaces.size();}
int getNumBubbles()  const {return m_bubbles.size();}

SimTK_DOWNCAST(ContactTrackerSubsystemImpl, Subsystem::Guts);

private:
friend class ContactTrackerSubsystem;

    // TOPOLOGY STATE
// Always order the key with the lower numbered geometry type first but
// if that is the reverse from how the tracker is defined then the bool 
// should be set to true meaning "reverse arguments".
// This map owns the heap space for the ContactTrackers; be sure to 
// delete it when replacing or destructing.
TrackerMap          m_contactTrackers;
ContactTracker*     m_defaultTracker;

    // TOPOLOGY CACHE
Array_<Surface,ContactSurfaceIndex> m_surfaces;
Array_<Bubble,BubbleIndex>          m_bubbles;
DiscreteVariableIndex               m_activeContactsIx;
DiscreteVariableIndex               m_predictedContactsIx;
};



//==============================================================================
//                        CONTACT TRACKER SUBSYSTEM
//==============================================================================

bool ContactTrackerSubsystem::isInstanceOf(const Subsystem& s) {
    return ContactTrackerSubsystemImpl::isA(s.getSubsystemGuts());
}
const ContactTrackerSubsystem& ContactTrackerSubsystem::
downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const ContactTrackerSubsystem&>(s);
}
ContactTrackerSubsystem& ContactTrackerSubsystem::
updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<ContactTrackerSubsystem&>(s);
}

const ContactTrackerSubsystemImpl& ContactTrackerSubsystem::
getImpl() const {
    return dynamic_cast<const ContactTrackerSubsystemImpl&>(getSubsystemGuts());
}
ContactTrackerSubsystemImpl& ContactTrackerSubsystem::
updImpl() {
    return dynamic_cast<ContactTrackerSubsystemImpl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use
// except for making std::vectors, which require a default constructor to be 
// available.
ContactTrackerSubsystem::ContactTrackerSubsystem() 
{   adoptSubsystemGuts(new ContactTrackerSubsystemImpl()); }

ContactTrackerSubsystem::ContactTrackerSubsystem(MultibodySystem& mbs) 
{   adoptSubsystemGuts(new ContactTrackerSubsystemImpl());
    mbs.adoptSubsystem(*this); } // steal ownership

int ContactTrackerSubsystem::getNumSurfaces() const
{   return getImpl().getNumSurfaces(); }

const MobilizedBody& ContactTrackerSubsystem::
getMobilizedBody(ContactSurfaceIndex surfIx) const
{   return *getImpl().m_surfaces[surfIx].mobod; }

const ContactSurface& ContactTrackerSubsystem::
getContactSurface(ContactSurfaceIndex surfIx) const
{   return *getImpl().m_surfaces[surfIx].surface; }

const Transform& ContactTrackerSubsystem::
getContactSurfaceTransform(ContactSurfaceIndex surfIx) const
{   return getImpl().m_surfaces[surfIx].X_BS; }

void ContactTrackerSubsystem::
adoptContactTracker(ContactTracker* tracker)
{   updImpl().adoptContactTracker(tracker); }

bool ContactTrackerSubsystem::
hasContactTracker(ContactGeometryTypeId surface1, 
                  ContactGeometryTypeId surface2) const
{   return getImpl().hasContactTracker(surface1,surface2); }

const ContactTracker& ContactTrackerSubsystem::
getContactTracker(ContactGeometryTypeId surface1, 
                  ContactGeometryTypeId surface2,
                  bool& reverseOrder) const
{   return getImpl().getContactTracker(surface1,surface2,reverseOrder); }

const ContactSnapshot& ContactTrackerSubsystem::
getPreviousActiveContacts(const State& state) const
{   return getImpl().getPrevActiveContacts(state); }

const ContactSnapshot& ContactTrackerSubsystem::
getPreviousPredictedContacts(const State& state) const
{   return getImpl().getPrevPredictedContacts(state); }

const ContactSnapshot& ContactTrackerSubsystem::
getActiveContacts(const State& state) const {
    Real dummy;
    realizeActiveContacts(state,true,dummy);
    return getImpl().getNextActiveContacts(state);
}

const ContactSnapshot& ContactTrackerSubsystem::
getPredictedContacts(const State& state) const {
    Real dummy;
    realizePredictedContacts(state,true,dummy);
    return getImpl().getNextPredictedContacts(state);
}

bool ContactTrackerSubsystem::
realizeActiveContacts(const State& state, 
                      bool         lastTry,
                      Real&        stepAdvice) const
{   // TODO: errors and step advice
    getImpl().ensureActiveContactsUpdated(state);
    return true;
}

bool ContactTrackerSubsystem::
realizePredictedContacts(const State& state, 
                         bool         lastTry,
                         Real&        stepAdvice) const
{   // TODO: errors and step advice
    getImpl().ensurePredictedContactsUpdated(state);
    return true;
}



//==============================================================================
//                     HALFSPACE-SPHERE CONTACT TRACKER
//==============================================================================
// Cost is 21 flops if no contact, 43 with contact.
bool ContactTracker::HalfSpaceSphere::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GS, 
    const ContactGeometry& geoSphere,
    Real                   rangeOfInterest,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoSphere.getTypeId()==ContactGeometry::Sphere::classTypeId(),
       "ContactTracker::HalfSpaceSphere::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::SphereImpl& sphere = 
        reinterpret_cast<const ContactGeometry::SphereImpl&>
            (geoSphere.getImpl());
    
    // Location of the sphere center C in the half-space's frame "H", in
    // which the half space occupies all of x>0 space.
    const Vec3 p_HC=~X_GH*X_GS.p();   // vec from OH to C  (18 flops)
    const Real r = sphere.getRadius();
    const Real depth = p_HC[0] + r;   // 1 flop

    if (depth <= -rangeOfInterest) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    const UnitVec3 normal_G = -X_GH.x(); // 3 flops; no normalization required
    const Vec3 origin_H = Vec3(depth/2, p_HC[1], p_HC[2]); // 1 flop
    const Vec3 origin_G = X_GH*origin_H; // 18 flops
    // The surfaces are contacting (or close enough to be interesting).
    // The sphere's radius is also the effective radius.
    currentStatus = CircularPointContact(priorStatus.getSurface1(), Infinity,
                                         priorStatus.getSurface2(), r,
                                         r, depth, origin_G, normal_G);
    return true; // success
}

bool ContactTracker::HalfSpaceSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   rangeOfInterest,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceSphere::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceSphere::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   rangeOfInterest,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceSphere::initializeContact() not implemented yet."); 
    return false; }


//==============================================================================
//                       SPHERE-SPHERE CONTACT TRACKER
//==============================================================================
// Cost is 12 flops if no contact, 82 if contact
bool ContactTracker::SphereSphere::trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& geoSphere1,
    const Transform& X_GS2, 
    const ContactGeometry& geoSphere2,
    Real                   cutoff, // 0 for contact
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoSphere1.getTypeId()==ContactGeometry::Sphere::classTypeId()
        && geoSphere2.getTypeId()==ContactGeometry::Sphere::classTypeId(),
       "ContactTracker::SphereSphere::trackContact()");

    // No need for an expensive dynamic casts here; we know what we have.
    const ContactGeometry::SphereImpl& sphere1 = 
        reinterpret_cast<const ContactGeometry::SphereImpl&>
            (geoSphere1.getImpl());
    const ContactGeometry::SphereImpl& sphere2 = 
        reinterpret_cast<const ContactGeometry::SphereImpl&>
            (geoSphere2.getImpl());

    currentStatus.clear();

    // Find the vector from sphere center C1 to C2, expressed in G.
    const Vec3 v_C12_G = X_GS2.p() - X_GS1.p(); // 3 flops
    const Real d2 = v_C12_G.normSqr();          // 5 flops
    const Real r1 = sphere1.getRadius();
    const Real r2 = sphere2.getRadius();
    const Real rr = r1 + r2;                    // 1 flop

    // Quick check. If separated we can return nothing, unless we were
    // in contact last time in which case we have to return one last
    // Contact indicating that contact has been broken and by how much.
    if (d2 > square(rr+cutoff)) {       // 3 flops
        if (!priorStatus.getContactId().isValid())
            return true; // successful return: still separated
        const Real separation = std::sqrt(d2) - rr; // > cutoff
        currentStatus = BrokenContact(priorStatus.getSurface1(),
                                      priorStatus.getSurface2(),
                                      separation);
        return true;
    }

    const Real d = std::sqrt(d2); // ~20 flops
    if (d < SignificantReal) {    // 1 flop
        // TODO: If the centers are coincident we should use past information
        // to determine the most likely normal. For now just fail.
        return false;
    }

    const Real depth = rr - d; // >0 for penetration (1 flop)
    const Real r     = r1*r2/rr; // r=r1r2/(r1+r2)=1/(1/r1+1/r2) ~20 flops

    const UnitVec3 normal_G(v_C12_G/d, true); // 1/ + 3* = ~20 flops
    const Vec3     origin_G = X_GS1.p() + (r1 - depth/2)*normal_G; // 8 flops

    // The surfaces are contacting (or close enough to be interesting).
    currentStatus = CircularPointContact(priorStatus.getSurface1(), r1,
                                         priorStatus.getSurface2(), r2,
                                         r, depth, origin_G, normal_G);
    return true; // success
}

bool ContactTracker::SphereSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   rangeOfInterest,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereSphere::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::SphereSphere::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   rangeOfInterest,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereSphere::initializeContact() not implemented yet."); 
    return false; }



} // namespace SimTK

