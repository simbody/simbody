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

#include "SimTKmath.h"
#include "simbody/internal/common.h"

#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
//#include "SimTKcommon/internal/SubsystemGuts.h"
#include "simbody/internal/Contact.h"

#include "simbody/internal/ContactTrackerSubsystem.h"

#include <algorithm>
using std::pair; using std::make_pair;
#include <iostream>
using std::cout; using std::endl;
#include <set>


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
    adoptContactTracker(new ContactTracker::HalfSpaceEllipsoid());
    adoptContactTracker(new ContactTracker::HalfSpaceTriangleMesh());
    adoptContactTracker(new ContactTracker::SphereTriangleMesh());
    adoptContactTracker(new ContactTracker::TriangleMeshTriangleMesh());
}

~ContactTrackerSubsystemImpl() {
    delete m_defaultTracker;
    TrackerMap::iterator p = m_contactTrackers.begin();
    for (; p != m_contactTrackers.end(); ++p)
        delete p->second.first; // the tracker
    // The map itself gets deleted automatically.
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
                else { // was NewContact or Ongoing; now Ongoing or Broken
                    assert(prev->getCondition()==Contact::NewContact
                           || prev->getCondition()==Contact::Ongoing);
                    if (next.getTypeId() != BrokenContact::classTypeId())
                        next.setCondition(Contact::Ongoing);
                    // Condition will already by Broken for a BrokenContact
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
// Cost is 21 flops if no contact, 67 with contact.
bool ContactTracker::HalfSpaceSphere::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GS, 
    const ContactGeometry& geoSphere,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoSphere.getTypeId()==ContactGeometry::Sphere::classTypeId(),
       "ContactTracker::HalfSpaceSphere::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Sphere& sphere = 
        ContactGeometry::Sphere::getAs(geoSphere);

    const Rotation R_HG = ~X_GH.R(); // inverse rotation; no flops

    // p_HC is vector from H origin to S's center C
    const Vec3 p_HC = R_HG*(X_GS.p() - X_GH.p()); // 18 flops
    
    // Calculate depth of sphere center C given that the halfspace occupies 
    // all of x>0 space.
    const Real r = sphere.getRadius();
    const Real depth = p_HC[0] + r;   // 1 flop

    if (depth <= -cutoff) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // Calculate the rest of the X_HS transform as required by Contact.
    const Transform X_HS(R_HG*X_GS.R(), p_HC); // 45 flops
    const UnitVec3 normal_H(Vec3(-1,0,0), true); // 0 flops
    const Vec3     origin_H = Vec3(depth/2, p_HC[1], p_HC[2]); // 1 flop
    // The surfaces are contacting (or close enough to be interesting).
    // The sphere's radius is also the effective radius.
    currentStatus = CircularPointContact(priorStatus.getSurface1(), Infinity,
                                         priorStatus.getSurface2(), r,
                                         X_HS, r, depth, origin_H, normal_H);
    return true; // success
}

bool ContactTracker::HalfSpaceSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
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
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceSphere::initializeContact() not implemented yet."); 
    return false; }



//==============================================================================
//                     HALFSPACE-ELLIPSOID CONTACT TRACKER
//==============================================================================
// Cost is ~135 flops if no contact, ~425 with contact.
// The contact point on the ellipsoid must be the unique point that has its
// outward-facing normal in the opposite direction of the half space normal.
// We can find that point very fast and see how far it is from the half
// space surface. If it is close enough, we'll evaluate the curvatures at
// that point in preparation for generating forces with Hertz theory.
bool ContactTracker::HalfSpaceEllipsoid::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GE, 
    const ContactGeometry& geoEllipsoid,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoEllipsoid.getTypeId()==ContactGeometry::Ellipsoid::classTypeId(),
       "ContactTracker::HalfSpaceEllipsoid::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Ellipsoid& ellipsoid = 
        ContactGeometry::Ellipsoid::getAs(geoEllipsoid);

    // Our half space occupies the +x half so the normal is -x.
    const Transform X_HE = ~X_GH*X_GE; // 63 flops
    // Halfspace normal is -x, so the ellipsoid normal we're looking for is
    // in the half space's +x direction.
    const UnitVec3& n_E = (~X_HE.R()).x(); // halfspace normal in E
    const Vec3 Q_E = ellipsoid.findPointWithThisUnitNormal(n_E); // 40 flops
    const Vec3 Q_H = X_HE*Q_E; // Q measured from half space origin (18 flops)
    const Real depth = Q_H[0]; // x > 0 is penetrated

    if (depth <= -cutoff) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // The surfaces are contacting (or close enough to be interesting).
    // The ellipsoid's principal curvatures k at the contact point are also
    // the curvatures of the contact paraboloid since the half space doesn't
    // add anything interesting.
    Transform X_EQ; Vec2 k;
    ellipsoid.findParaboloidAtPointWithNormal(Q_E, n_E, X_EQ, k); // 220 flops

    // We have the contact paraboloid expressed in frame Q but Qz=n_E has the
    // wrong sign since we have to express it using the half space normal.
    // We have to end up with a right handed frame, so one of x or y has
    // to be negated too. (6 flops)
    Rotation& R_EQ = X_EQ.updR();
    R_EQ.setRotationColFromUnitVecTrustMe(ZAxis, -R_EQ.z()); // changing X_EQ
    R_EQ.setRotationColFromUnitVecTrustMe(XAxis, -R_EQ.x());

    // Now the frame is pointing in the right direction. Measure and express in 
    // half plane frame, then shift origin to half way between contact point 
    // Q on the undeformed ellipsoid and the corresponding contact point P 
    // on the undeformed half plane surface. It's easier to do this shift
    // in H since it is in the -Hx direction.
    Transform X_HC = X_HE*X_EQ; X_HC.updP()[0] -= depth/2; // 65 flops

    currentStatus = EllipticalPointContact(priorStatus.getSurface1(),
                                           priorStatus.getSurface2(),
                                           X_HE, X_HC, k, depth);
    return true; // success
}

bool ContactTracker::HalfSpaceEllipsoid::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceEllipsoid::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceEllipsoid::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceEllipsoid::initializeContact() not implemented yet."); 
    return false; }



//==============================================================================
//                       SPHERE-SPHERE CONTACT TRACKER
//==============================================================================
// Cost is 12 flops if no contact, 139 if contact
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
    const ContactGeometry::Sphere& sphere1 = 
        ContactGeometry::Sphere::getAs(geoSphere1);
    const ContactGeometry::Sphere& sphere2 = 
        ContactGeometry::Sphere::getAs(geoSphere2);

    currentStatus.clear();

    // Find the vector from sphere center C1 to C2, expressed in G.
    const Vec3 p_12_G = X_GS2.p() - X_GS1.p(); // 3 flops
    const Real d2 = p_12_G.normSqr();          // 5 flops
    const Real r1 = sphere1.getRadius();
    const Real r2 = sphere2.getRadius();
    const Real rr = r1 + r2;                    // 1 flop

    // Quick check. If separated we can return nothing, unless we were
    // in contact last time in which case we have to return one last
    // Contact indicating that contact has been broken and by how much.
    if (d2 > square(rr+cutoff)) {       // 3 flops
        if (!priorStatus.getContactId().isValid())
            return true; // successful return: still separated

        const Real separation = std::sqrt(d2) - rr;   // > cutoff, ~25 flops
        const Transform X_S1S2(~X_GS1.R()*X_GS2.R(), 
                               ~X_GS1.R()*p_12_G);    // 60 flops
        currentStatus = BrokenContact(priorStatus.getSurface1(),
                                      priorStatus.getSurface2(),
                                      X_S1S2, separation);
        return true;
    }

    const Real d = std::sqrt(d2); // ~20 flops
    if (d < SignificantReal) {    // 1 flop
        // TODO: If the centers are coincident we should use past information
        // to determine the most likely normal. For now just fail.
        return false;
    }

    const Transform X_S1S2(~X_GS1.R()*X_GS2.R(), 
                           ~X_GS1.R()*p_12_G);// 60 flops
    const Vec3& p_12 = X_S1S2.p(); // center-to-center vector in S1

    const Real depth = rr - d; // >0 for penetration (1 flop)
    const Real r     = r1*r2/rr; // r=r1r2/(r1+r2)=1/(1/r1+1/r2) ~20 flops

    const UnitVec3 normal_S1(p_12/d, true); // 1/ + 3* = ~20 flops
    const Vec3     origin_S1 = (r1 - depth/2)*normal_S1; // 5 flops

    // The surfaces are contacting (or close enough to be interesting).
    currentStatus = CircularPointContact(priorStatus.getSurface1(), r1,
                                         priorStatus.getSurface2(), r2,
                                         X_S1S2, r, depth, 
                                         origin_S1, normal_S1);
    return true; // success
}

bool ContactTracker::SphereSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
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
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereSphere::initializeContact() not implemented yet."); 
    return false; }




//==============================================================================
//                  HALFSPACE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::HalfSpaceTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GM, 
    const ContactGeometry& geoMesh,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::HalfSpace::isInstance(geoHalfSpace)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh),
       "ContactTracker::HalfSpaceTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::HalfSpaceTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::TriangleMesh& mesh = 
        ContactGeometry::TriangleMesh::getAs(geoMesh);

    // Transform giving mesh (S2) frame in the halfspace (S1) frame.
    const Transform X_HM = (~X_GH)*X_GM; 

    // Normal is halfspace -x direction; xdir is first column of R_MH.
    // That's a unit vector and -unitvec is also a unit vector so this
    // doesn't require normalization.
    const UnitVec3 hsNormal_M = -(~X_HM.R()).x();
    // Find the height of the halfspace face along the normal, measured
    // from the mesh origin.
    const Real hsFaceHeight_M = dot((~X_HM).p(), hsNormal_M);
    // Now collect all the faces that are all or partially below the 
    // halfspace surface.
    std::set<int> insideFaces;
    processBox(mesh, mesh.getOBBTreeNode(), X_HM, 
               hsNormal_M, hsFaceHeight_M, insideFaces);
    
    if (insideFaces.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_HM, 
                                        std::set<int>(), insideFaces);
    return true; // success
}

bool ContactTracker::HalfSpaceTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceTriangleMesh::initializeContact() not implemented yet."); 
    return false; }



// Check a single OBB and its contents (recursively) against the halfspace,
// appending any penetrating faces to the insideFaces list.
void ContactTracker::HalfSpaceTriangleMesh::processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Transform& X_HM, const UnitVec3& hsNormal_M, Real hsFaceHeight_M, 
    std::set<int>& insideFaces) const 
{   // First check against the node's bounding box.
    
    const OrientedBoundingBox& bounds = node.getBounds();
    const Transform& X_MB = bounds.getTransform(); // box frame in mesh
    const Vec3 p_BC = bounds.getSize()/2; // from box origin corner to center
    // Express the half space normal in the box frame, then reflect it into
    // the first (+,+,+) quadrant where it is the normal of a different 
    // but symmetric and more convenient half space.
    const UnitVec3 octant1hsNormal_B = (~X_MB.R()*hsNormal_M).abs();
    // Dot our octant1 radius p_BC with our octant1 normal to get
    // the extent of the box from its center in the direction of the octant1
    // reflection of the halfspace.
    const Real extent = dot(p_BC, octant1hsNormal_B);
    // Compute the height of the box center over the mesh origin,
    // measured along the real halfspace normal.
    const Vec3 boxCenter_M       = X_MB*p_BC;
    const Real boxCenterHeight_M = dot(boxCenter_M, hsNormal_M);
    // Subtract the halfspace surface position to get the height of the 
    // box center over the halfspace.
    const Real boxCenterHeight = boxCenterHeight_M - hsFaceHeight_M;
    if (boxCenterHeight >= extent)
        return;                             // no penetration
    if (boxCenterHeight <= -extent) {
        addAllTriangles(node, insideFaces); // box is entirely in halfspace
        return;
    }
    
    // Box is partially penetrated into halfspace. If it is not a leaf node, 
    // check its children.
    if (!node.isLeafNode()) {
        processBox(mesh, node.getFirstChildNode(), X_HM, hsNormal_M, 
                   hsFaceHeight_M, insideFaces);
        processBox(mesh, node.getSecondChildNode(), X_HM, hsNormal_M, 
                   hsFaceHeight_M, insideFaces);
        return;
    }
    
    // This is a leaf OBB node that is penetrating, so some of its triangles
    // may be penetrating.
    const Array_<int>& triangles = node.getTriangles();
    for (int i = 0; i < (int) triangles.size(); i++) {
        for (int vx=0; vx < 3; ++vx) {
            const int   vertex         = mesh.getFaceVertex(triangles[i], vx);
            const Vec3& vertexPos      = mesh.getVertexPosition(vertex);
            const Real  vertexHeight_M = dot(vertexPos, hsNormal_M);
            if (vertexHeight_M < hsFaceHeight_M) {
                insideFaces.insert(triangles[i]);
                break; // done with this face
            }
        }
    }
}

void ContactTracker::HalfSpaceTriangleMesh::addAllTriangles
   (const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    std::set<int>& insideFaces) const 
{
    if (node.isLeafNode()) {
        const Array_<int>& triangles = node.getTriangles();
        for (int i = 0; i < (int) triangles.size(); i++)
            insideFaces.insert(triangles[i]);
    }
    else {
        addAllTriangles(node.getFirstChildNode(), insideFaces);
        addAllTriangles(node.getSecondChildNode(), insideFaces);
    }
}





//==============================================================================
//                  SPHERE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::SphereTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GS, 
    const ContactGeometry& geoSphere,
    const Transform&       X_GM, 
    const ContactGeometry& geoMesh,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::Sphere::isInstance(geoSphere)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh),
       "ContactTracker::SphereTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::SphereTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Sphere&          sphere = 
        ContactGeometry::Sphere::getAs(geoSphere);
    const ContactGeometry::TriangleMesh&    mesh = 
        ContactGeometry::TriangleMesh::getAs(geoMesh);

    // Transform giving mesh (M) frame in the sphere (S) frame.
    const Transform X_SM = ~X_GS*X_GM; 

    // Want the sphere center measured and expressed in the mesh frame.
    const Vec3 p_MC = (~X_SM).p();
    std::set<int> insideFaces;
    processBox(mesh, mesh.getOBBTreeNode(), p_MC, square(sphere.getRadius()), 
               insideFaces);
    
    if (insideFaces.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_SM, 
                                        std::set<int>(), insideFaces);
    return true; // success
}

// Check a single OBB and its contents (recursively) against the sphere
// whose center location in M and radius squared is given, appending any 
// penetrating faces to the insideFaces list.
void ContactTracker::SphereTriangleMesh::processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Vec3& center_M, Real radius2, 
    std::set<int>& insideFaces) const 
{   // First check against the node's bounding box.

    const Vec3 nearest_M = node.getBounds().findNearestPoint(center_M);
    if ((nearest_M-center_M).normSqr() >= radius2)
        return; // no intersection possible
    
    // Bounding box is penetrating. If it's not a leaf node, check its children.
    if (!node.isLeafNode()) {
        processBox(mesh, node.getFirstChildNode(), center_M, radius2,
                   insideFaces);
        processBox(mesh, node.getSecondChildNode(), center_M, radius2,
                   insideFaces);
        return;
    }
    
    // This is a leaf node that may be penetrating; check the triangles.
    const Array_<int>& triangles = node.getTriangles();
    for (unsigned i = 0; i < triangles.size(); i++) {
        Vec2 uv;
        Vec3 nearest_M = mesh.findNearestPointToFace
                                    (center_M, triangles[i], uv);
        if ((nearest_M-center_M).normSqr() < radius2)
            insideFaces.insert(triangles[i]);
    }
}

bool ContactTracker::SphereTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::SphereTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereTriangleMesh::initializeContact() not implemented yet."); 
    return false; }






//==============================================================================
//               TRIANGLE MESH - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::TriangleMeshTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GM1, 
    const ContactGeometry& geoMesh1,
    const Transform&       X_GM2, 
    const ContactGeometry& geoMesh2,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::TriangleMesh::isInstance(geoMesh1)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh2),
       "ContactTracker::TriangleMeshTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::TriangleMeshTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::TriangleMesh& mesh1 = 
        ContactGeometry::TriangleMesh::getAs(geoMesh1);
    const ContactGeometry::TriangleMesh& mesh2 = 
        ContactGeometry::TriangleMesh::getAs(geoMesh2);

    // Transform giving mesh2 (M2) frame in the mesh1 (M1) frame.
    const Transform X_M1M2 = ~X_GM1*X_GM2; 
    std::set<int> insideFaces1, insideFaces2;

    // Get M2's bounding box in M1's frame.
    const OrientedBoundingBox 
        mesh2Bounds_M1 = X_M1M2*mesh2.getOBBTreeNode().getBounds();

    // Find the faces that are actually intersecting faces on the other
    // surface (this doesn't yet include faces that may be completely buried).
    findIntersectingFaces(mesh1, mesh2, 
                          mesh1.getOBBTreeNode(), mesh2.getOBBTreeNode(), 
                          mesh2Bounds_M1, X_M1M2, insideFaces1, insideFaces2);
    
    // It should never be the case that one set of faces is empty and the
    // other isn't, however it is conceivable that roundoff error could cause
    // this to happen so we'll check both lists.
    if (insideFaces1.empty() && insideFaces2.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    // There was an intersection. We now need to identify every triangle and 
    // vertex of each mesh that is inside the other mesh. We found the border
    // intersections above; now we have to fill in the buried faces.
    findBuriedFaces(mesh1, mesh2, ~X_M1M2, insideFaces1);
    findBuriedFaces(mesh2, mesh1,  X_M1M2, insideFaces2);

    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_M1M2, 
                                        insideFaces1, insideFaces2);
    return true; // success
}

void ContactTracker::TriangleMeshTriangleMesh::
findIntersectingFaces
   (const ContactGeometry::TriangleMesh&                mesh1, 
    const ContactGeometry::TriangleMesh&                mesh2,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node1, 
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node2, 
    const OrientedBoundingBox&                          node2Bounds_M1,
    const Transform&                                    X_M1M2, 
    std::set<int>&                                      triangles1, 
    std::set<int>&                                      triangles2) const 
{   // See if the bounding boxes intersect.
    
    if (!node1.getBounds().intersectsBox(node2Bounds_M1))
        return;
    
    // If either node is not a leaf node, process the children recursively.
    
    if (!node1.isLeafNode()) {
        if (!node2.isLeafNode()) {
            const OrientedBoundingBox firstChildBounds = 
                X_M1M2*node2.getFirstChildNode().getBounds();
            const OrientedBoundingBox secondChildBounds = 
                X_M1M2*node2.getSecondChildNode().getBounds();
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        }
        else {
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2, node2Bounds_M1, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2, node2Bounds_M1, X_M1M2, triangles1, triangles2);
        }
        return;
    }
    else if (!node2.isLeafNode()) {
        const OrientedBoundingBox firstChildBounds = 
            X_M1M2*node2.getFirstChildNode().getBounds();
        const OrientedBoundingBox secondChildBounds = 
            X_M1M2*node2.getSecondChildNode().getBounds();
        findIntersectingFaces(mesh1, mesh2, node1, node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
        findIntersectingFaces(mesh1, mesh2, node1, node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        return;
    }
    
    // These are both leaf nodes, so check triangles for intersections.
    
    const Array_<int>& node1triangles = node1.getTriangles();
    const Array_<int>& node2triangles = node2.getTriangles();
    for (unsigned i = 0; i < node2triangles.size(); i++) {
        const int face2 = node2triangles[i];
        Vec3 a1 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 0));
        Vec3 a2 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 1));
        Vec3 a3 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 2));
        const Geo::Triangle A(a1,a2,a3);
        for (unsigned j = 0; j < node1triangles.size(); j++) {
            const int face1 = node1triangles[j];
            const Vec3& b1 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 0));
            const Vec3& b2 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 1));
            const Vec3& b3 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 2));
            const Geo::Triangle B(b1,b2,b3);
            if (A.overlapsTriangle(B)) 
            {   // The triangles intersect.
                triangles1.insert(face1);
                triangles2.insert(face2);
            }
        }
    }
}

static const int Outside  = -1;
static const int Unknown  =  0;
static const int Boundary =  1;
static const int Inside   =  2;

void ContactTracker::TriangleMeshTriangleMesh::
findBuriedFaces(const ContactGeometry::TriangleMesh&    mesh,       // M 
                const ContactGeometry::TriangleMesh&    otherMesh,  // O
                const Transform&                        X_OM, 
                std::set<int>&                          insideFaces) const 
{  
    // Find which faces are inside.
    // We're passed in the list of Boundary faces, that is, those faces of
    // "mesh" that intersect faces of "otherMesh".
    Array_<int> faceType(mesh.getNumFaces(), Unknown);
    for (std::set<int>::iterator iter = insideFaces.begin(); 
                                 iter != insideFaces.end(); ++iter)
        faceType[*iter] = Boundary;

    for (int i = 0; i < (int) faceType.size(); i++) {
        if (faceType[i] == Unknown) {
            // Trace a ray from its center to determine whether it is inside.           
            const Vec3     origin_O    = X_OM    * mesh.findCentroid(i);
            const UnitVec3 direction_O = X_OM.R()* mesh.getFaceNormal(i);
            Real distance;
            int face;
            Vec2 uv;
            if (   otherMesh.intersectsRay(origin_O, direction_O, distance, 
                                           face, uv) 
                && ~direction_O*otherMesh.getFaceNormal(face) > 0) 
            {
                faceType[i] = Inside;
                insideFaces.insert(i);
            } else
                faceType[i] = Outside;
            
            // Recursively mark adjacent inside or outside Faces.           
            tagFaces(mesh, faceType, insideFaces, i, 0);
        }
    }
}

//TODO: the following method uses depth-first recursion to iterate through
//unmarked faces. For a large mesh this was observed to produce a stack
//overflow in OpenSim. Here we limit the recursion depth; after we get that
//deep we'll pop back out and do another expensive intersectsRay() test in
//the method above.
static const int MaxRecursionDepth = 500;

void ContactTracker::TriangleMeshTriangleMesh::
tagFaces(const ContactGeometry::TriangleMesh&   mesh, 
         Array_<int>&                           faceType,
         std::set<int>&                         triangles, 
         int                                    index,
         int                                    depth) const 
{
    for (int i = 0; i < 3; i++) {
        const int edge = mesh.getFaceEdge(index, i);
        const int face = (mesh.getEdgeFace(edge, 0) == index 
                            ? mesh.getEdgeFace(edge, 1) 
                            : mesh.getEdgeFace(edge, 0));
        if (faceType[face] == Unknown) {
            faceType[face] = faceType[index];
            if (faceType[index] > 0)
                triangles.insert(face);
            if (depth < MaxRecursionDepth)
                tagFaces(mesh, faceType, triangles, face, depth+1);
        }
    }
}


bool ContactTracker::TriangleMeshTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::TriangleMeshTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::TriangleMeshTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::TriangleMeshTriangleMesh::initializeContact() not implemented yet."); 
    return false; }



} // namespace SimTK

