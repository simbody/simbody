/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Implementations of non-inline methods of SubsystemRep.
 */

#include "simbody/SimbodyCommon.h"
#include "SubsystemRep.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::cout;
using std::endl;
using std::ostream;


// Returns -1, 0, 1 according to key {<,==,>} test ignoring case.
static int caseInsensitiveCompare(const std::string& key, const std::string& test);

namespace simtk {

    // PLACEMENT SLOT //



// We have just copied a Subsystem tree and this PlacementSlot is the new copy. If
// it had a valueSlot, that valueSlot is still pointing into the old Subsystem tree
// and needs to be repaired to point to the corresponding valueSlot in the new tree.
void PlacementSlot::repairValueReference(const Subsystem& oldRoot, const Subsystem& newRoot) {
    if (valueSlot) {
        valueSlot = const_cast<PlacementValueSlot*>
                        (FeatureRep::findCorrespondingPlacementValueSlot
                                    (oldRoot,*valueSlot,newRoot));
        if (valueSlot)
            valueSlot->setClientPlacementSlot(*this);
    }
}

String PlacementSlot::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "PlacementSlot ";
    if (hasOwner())
        s << getOwner().getFullName() << ":"
          << std::left << std::setw(2) << getIndexInOwner();
    else s << "NO OWNER";
    if (hasClientFeature())
        s << "[client:" << getClientFeature().getFullName() << "]";
    else s << "[NO CLIENT]";
    s << " " << getPlacement().toString(linePrefix);
    return s.str();
}

void PlacementSlot::checkPlacementConsistency(const Subsystem* expOwner, 
                                              int              expIndexInOwner,
                                              const Subsystem& expRoot) const
{
    cout << "CHECK PLACEMENT SLOT CONSISTENCY FOR PlacementSlot@" << this << endl;

    if (owner != expOwner)
        cout << "*** WRONG OWNER@" << owner << "; should have been " << expOwner << endl;
    if (indexInOwner != expIndexInOwner)
        cout << "*** WRONG INDEX " << indexInOwner << "; should have been " 
             << expIndexInOwner << endl;

    if (expOwner == 0) {
        if (client)
          cout << "*** UNOWNED PLACEMENT HAD CLIENT " << client->getFullName() << " ***" << endl;
    } else {
        if (!client) 
            cout << "*** NO CLIENT ***" << endl;
        else if (!client->getRep().hasPlacement())
            cout << "*** CLIENT " << client->getFullName() << " HAS NO PLACEMENT??? ***" << endl;
        else if (&(client->getRep().getPlacementSlot()) != this)
            cout << "*** CLIENT " << client->getFullName() << " HAS WRONG PLACEMENT SLOT@" 
                << &client->getRep().getPlacementSlot() << endl;
    }

    if (hasValueSlot()) {
        std::string nm = client ? client->getFullName() : "(NO CLIENT)";
        if (!getValueSlot().hasOwner())
            cout << "*** Referenced entry for " << nm << "'s value slot is unowned." << endl;
        else if (!getValueSlot().getOwner().getRep().findRootSubsystem().isSameSubsystem(expRoot))
            cout << "*** Referenced entry " << nm << "'s value slot is in wrong tree." << endl;
    }

    // Check the Placement
    if (!placement.hasRep()) 
        cout << "*** PLACEMENT HAS NO REP ***" << endl;
    else {
        if (&placement.getRep().getMyHandle() != &placement)
            cout << "*** placement.rep->handle=" << &placement.getRep().getMyHandle()
                << " which is *** WRONG ***" << endl;
        
        const Feature* offender;
        if (!placement.getRep().isLimitedToSubtree(expRoot, offender)) {
            cout << "*** Placement referenced Feature '" << offender->getFullName() 
                 << "' in wrong tree" << endl;
            cout << "*** Root should have been @" << &expRoot << ", not " << 
                &offender->getRep().findRootSubsystem() << endl;
        }
    }
}

    // SUBSYSTEM REP //

void SubsystemRep::realize(Stage g) const {
    // Always let the children go first.
    for (int i=0; i < getNSubsystems(); ++i)
        getSubsystem(i).realize(g);

    // For startup stage, we assume that we have recently modified the 
    // model, i.e., features and placements have changed. So we need to
    // (re)calculate the set of PlacementValues held by this Subsystem.
    // Note that we just did this for all the children, so all lower
    // PlacementValues should have been allocated. Nevertheless we have
    // to do this recursively to find all the unassigned Placements whose
    // highest Feature reference is one of the Placement expressions held here.
    if (g == Stage::Startup) {
        deleteAllPlacementValuesInThisSubsystem();
        const_cast<SubsystemRep*>(this)->finalizeStandardSubfeatures();
        allocatePlacementValueSlots(getMyHandle());
    }

    // TODO: for now we just assume all our values are invalid. Dependency
    // tracking would be better when costs are high.
    for (int i=0; i < getNPlacementValues(); ++i)
        updPlacementValueSlot(i).setValid(false);

    // TODO: current method will have side effect of realizing needed
    // PlacementValues as we go, so many of these will be evaluated
    // already by the time we get there.

    for (int i=0; i < getNPlacementValues(); ++i)
        if (!getPlacementValueSlot(i).isValid())
            getPlacementValueSlot(i).getClientPlacementSlot().realize(g);
}

const Subsystem*
SubsystemRep::findSubsystem(const std::string& nm) const {
    std::vector<std::string> segments;
    if (!isLegalSubsystemPathname(nm, &segments)) {
        if (segments.size())
            SIMTK_THROW2(Exception::IllegalSubsystemPathname,nm,segments.back());
        else
            SIMTK_THROW(Exception::EmptySubsystemPathname);
        //NOTREACHED
    }
    const Subsystem* found = &getMyHandle();
    for (size_t i=0; i<segments.size(); ++i) {
        size_t index;
        found = (found->getRep().findSubsystemIndex(segments[i],index) 
                    ? &found->getRep().childSubsystems[index] : 0);
        if (!found) break;
    }
    return found;
}

std::string 
SubsystemRep::getFullName() const { 
    std::string s;
    if (hasParentSubsystem())
        s = getParentSubsystem().getFullName() + "/";
    return s + getName(); 
}

const Subsystem& SubsystemRep::findRootSubsystem() const {
    return hasParentSubsystem() 
        ? getParentSubsystem().getRep().findRootSubsystem() : getMyHandle();
}
Subsystem& SubsystemRep::findUpdRootSubsystem() {
    return hasParentSubsystem() 
        ? updParentSubsystem().updRep().findUpdRootSubsystem() : updMyHandle();
}

// Wipe out all PlacementValues. Note that this routine is const because PlacementValues
// are mutable.
void SubsystemRep::deleteAllPlacementValuesInThisSubsystem() const {
    // Erase all references to the placement values (these may be in child subsystems).
    for (size_t i=0; i < placementValueSlots.size(); ++i) {
        PlacementValueSlot& pv = placementValueSlots[i];
        if (pv.hasClientPlacementSlot())
            pv.getClientPlacementSlot().clearValueSlot();
    }
    placementValueSlots.clear();
}

void SubsystemRep::deleteAllPlacementValuesFromHereToRoot() const {
    const SubsystemRep* nxt = this;
    do {
        nxt->deleteAllPlacementValuesInThisSubsystem();
        nxt = nxt->parent ? nxt->parent->rep : 0;
    } while (nxt);
}

void SubsystemRep::deletePlacementSlot(size_t i) {
    assert(i < placementSlots.size());
    deleteAllPlacementValuesFromHereToRoot(); // a little extreme, I know
    updPlacementSlot(i).updClientFeature().updRep().clearPlacementSlot();
    placementSlots.erase(i);
}

// Run around this subsystem and its children looking for Placements which can be
// evaluated without going any higher in the Subsystem tree than the target Subsysytem, and
// allocate PlacementValueSlots for them in the target.
void SubsystemRep::allocatePlacementValueSlots(const Subsystem& target) const {
    // Always let the children go first.
    for (int i=0; i < getNSubsystems(); ++i)
        getSubsystem(i).getRep().allocatePlacementValueSlots(target);

    for (int i=0; i < getNPlacements(); ++i) {
        const PlacementSlot& ps = getPlacementSlot(i);
        assert(ps.hasOwner());
        if (ps.hasValueSlot())
            continue;
        // can we get it a slot?
        const Subsystem* s = ps.getPlacement().getRep()
                                .findPlacementValueOwnerSubsystem(target);
        if (s && s->isSameSubsystem(target)) {
            PlacementValueSlot& pvs = 
                s->getRep().addPlacementValueLike(
                        ps.getPlacement().getRep().createEmptyPlacementValue());
            pvs.setClientPlacementSlot(ps);
            ps.assignValueSlot(pvs);
        }
    }
}
Subsystem& 
SubsystemRep::addSubsystemLike(const Subsystem& s, const std::string& nm) {
    assert(nm.size() > 0);
    const int index = (int)childSubsystems.size();
    childSubsystems.push_back(SubSubsystem()); // an empty handle
    Subsystem& newSubsystem = childSubsystems[index];
    s.getRep().cloneWithoutParentOrExternalPlacements(newSubsystem);
    newSubsystem.updRep().setParentSubsystem(updMyHandle(), index);
    newSubsystem.updRep().setName(nm);
    postProcessNewSubsystem(newSubsystem);
    return newSubsystem;
}

Feature& 
SubsystemRep::addFeatureLike(const Subsystem& f, const std::string& nm) {
    if (!Feature::isInstanceOf(f)) {
        SIMTK_THROW3(Exception::ExpectedFeaturePrototypeButGotSubsystem,
            getFullName(), nm, f.getName());
        //NOTREACHED
    }
    return reinterpret_cast<Feature&>(addSubsystemLike(f,nm));
}

void 
SubsystemRep::cloneWithoutParentOrExternalPlacements(Subsystem& newHandle) const
{
    SubsystemRep* copy = clone();
    copy->setMyHandle(newHandle);
    copy->parent = 0; copy->indexInParent = -1;
    newHandle.setRep(copy);

    // Re-parent all the copied child Subsystems to their new parent,
    // and fix the owned Placements to acknowledge their new owner.
    copy->reparentMyChildren();

    // Fix up all the internal placement references and delete the
    // external ones.
    copy->fixPlacements(this->getMyHandle(), copy->getMyHandle());
}


// Note that we can only allow placements involving internal features, e.g. children,
// grandchildren, etc. -- no external references. Otherwise someone further
// up the tree should own the new placement.
PlacementSlot& 
SubsystemRep::addPlacementLike(const Placement& p) {
    assert(p.hasRep());

    const Feature* offender;
    if (!p.getRep().isLimitedToSubtree(getMyHandle(),offender)) {
        SIMTK_THROW3(Exception::PlacementMustBeLocal,"SubsystemRep::addPlacementLike",
            this->getFullName(),offender->getFullName());
    }

    const int index = (int)placementSlots.size();
    placementSlots.push_back(PlacementSlot(p));
    PlacementSlot& newPlacementSlot = placementSlots[index];
    newPlacementSlot.setOwner(updMyHandle(), index);
    return newPlacementSlot;
}

PlacementValueSlot& 
SubsystemRep::addPlacementValueLike(const PlacementValue& v) const {
    assert(v.hasRep());

    // 'placementValueSlots' is mutable.
    const int index = (int)placementValueSlots.size();
    placementValueSlots.push_back(PlacementValueSlot(v));
    PlacementValueSlot& newPlacementValueSlot = placementValueSlots[index];
    newPlacementValueSlot.setOwner(getMyHandle(), index);
    return newPlacementValueSlot;
}

// Is Subsystem s in the tree rooted at oldRoot? If so, optionally return the 
// series of indices required to get to this Subsystem from the root.
// Complexity is O(log n) where n is tree depth.
/*static*/ bool 
SubsystemRep::isSubsystemInSubsystemTree(const Subsystem& oldRoot, const Subsystem& s,
                                         std::vector<int>* trace)
{
    if (trace) trace->clear();
    const Subsystem* const oldp = &oldRoot;
    const Subsystem*       sp   = &s;

    while (sp != oldp) {
        const Subsystem* const spParent = getParentPtr(*sp);
        if (!spParent) {
            if (trace) trace->clear(); // never mind ...
            return false;
        }
        if (trace) trace->push_back(sp->rep->getIndexInParent());
        sp = spParent;
    }

    return true;
}

// Is PlacementSlot p owned by a Feature in the tree rooted at oldRoot?
/*static*/ bool 
SubsystemRep::isPlacementInSubsystemTree(const Subsystem& oldRoot, 
                                         const PlacementSlot& p)
{
    if (!p.hasOwner())
        return false;   // a disembodied Placement
    return isSubsystemInSubsystemTree(oldRoot, p.getOwner());
}

// If Subsystem s is a member of the Subsystem tree rooted at oldRoot, find
// the corresponding Subsystem in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const Subsystem* 
SubsystemRep::findCorrespondingSubsystem
    (const Subsystem& oldRoot, const Subsystem& s, const Subsystem& newRoot)
{
    std::vector<int> trace;
    if (!isSubsystemInSubsystemTree(oldRoot,s,&trace))
        return 0;

    // Trace holds the indices needed to step from newRoot down to
    // the corresponding Feature (in reverse order).
    const Subsystem* newTreeRef = &newRoot;
    for (size_t i=trace.size(); i >=1; --i)
        newTreeRef = &newTreeRef->getRep().getSubsystem(trace[i-1]);
    return newTreeRef;
}

// Given two Subsystems, run up the tree towards the root to find
// their "least common denominator", i.e. the first shared node
// on the path back to the root. Return a pointer to that node
// if found, otherwise NULL meaning that the Subsystems aren't on
// the same tree. If the Subsystems are the same, then
// that Subsystem is the answer.
// Complexity is O(log n) (3 passes) where n is depth of Subsystem tree.

/*static*/ const Subsystem* 
SubsystemRep::findYoungestCommonAncestor(const Subsystem& s1, const Subsystem& s2)
{
    std::vector<const Subsystem*> s1path, s2path; // paths from nodes to their roots
    const Subsystem* s1p = &s1;
    const Subsystem* s2p = &s2;
    while (s1p) {s1path.push_back(s1p); s1p = getParentPtr(*s1p);}
    while (s2p) {s2path.push_back(s2p); s2p = getParentPtr(*s2p);}

    // If there is a common ancestor, we can find it by searching down from
    // the root (last element in each path). As soon as there is a difference,
    // the previous element is the common ancestor.
    const Subsystem* ancestor = 0;
    size_t i1 = s1path.size(), i2 = s2path.size();  // index of element just past end
    while (i1 && i2) {
        if (s1path[--i1] != s2path[--i2])
            break;
        ancestor = s1path[i1];
    }
    return ancestor;
}

/*static*/ Subsystem* 
SubsystemRep::findUpdYoungestCommonAncestor(Subsystem& s1, const Subsystem& s2) {
    return const_cast<Subsystem*>(findYoungestCommonAncestor(s1,s2));
}

// Debugging routine
void SubsystemRep::checkSubsystemConsistency(const Subsystem* expParent, 
                                         int expIndexInParent,
                                         const Subsystem& root) const {
    cout << "CHECK SUBSYSTEM CONSISTENCY FOR SubsystemRep@" << this 
         << "(" << getFullName() << ")" << endl;

    if (!myHandle) 
        cout << "*** NO HANDLE ***" << endl;
    else if (myHandle->rep != this)
        cout << "*** Handle->rep=" << myHandle->rep << " which is *** WRONG ***" << endl;

    if (parent != expParent)
        cout << " WRONG PARENT@" << parent << "; should have been " << expParent << endl;
    if (indexInParent != expIndexInParent)
        cout << "*** WRONG INDEX " << indexInParent << "; should have been " << expIndexInParent << endl;

    if (!findRootSubsystem().isSameSubsystem(root)) {
        cout << " WRONG ROOT@" << &findRootSubsystem() << "(" << findRootSubsystem().getFullName() << ")";
        cout << "; should have been " << &root << "(" << root.getFullName() << ")" << endl;
    }
    for (size_t i=0; i<(size_t)getNSubsystems(); ++i) 
        getSubsystem(i).checkSubsystemConsistency(&getMyHandle(), (int)i, root);
    for (size_t i=0; i<(size_t)getNPlacements(); ++i) 
        getPlacementSlot(i).checkPlacementConsistency(&getMyHandle(), (int)i, root);
    for (size_t i=0; i < (size_t)getNPlacementValues(); ++i)
        getPlacementValueSlot(i).checkPlacementValueConsistency(&getMyHandle(), (int)i, root);

    if (FeatureRep::isA(*this)) {
        const FeatureRep& fr = FeatureRep::downcast(*this);
        if (fr.hasPlacement()) {   // must be a Feature
            if (!fr.getPlacementSlot().hasOwner())
                cout << "*** Feature " << getFullName() << "'s placement is unowned." << endl;
            else if (!fr.getPlacementSlot().getOwner().getRep().findRootSubsystem().isSameSubsystem(root))
                cout << "*** Feature " << getFullName() << "'s placement is in wrong tree." << endl;
        }
    }
}

// Return true and ix==subsystem index if a subsystem of the given name is found.
// Otherwise return false and ix==childSubsystems.size().
bool 
SubsystemRep::findSubsystemIndex(const std::string& nm, size_t& ix) const {
    for (ix=0; ix < (size_t)getNSubsystems(); ++ix)
        if (caseInsensitiveCompare(nm, childSubsystems[ix].getName())==0)
            return true;
    return false;   // not found
}

// We have just copied a Subsystem subtree so all the parent pointers are
// wrong. Recursively repair them to point into the new tree.
void SubsystemRep::reparentMyChildren() {
    for (size_t i=0; i < (size_t)getNSubsystems(); ++i) {
        assert(childSubsystems[i].getRep().hasParentSubsystem());
        assert(childSubsystems[i].getRep().getIndexInParent() == i);    // shouldn't change
        childSubsystems[i].updRep().setParentSubsystem(updMyHandle(), i);
        childSubsystems[i].updRep().reparentMyChildren();               // recurse
    }
    for (size_t i=0; i < (size_t)getNPlacements(); ++i) {
        assert(placementSlots[i].hasOwner());
        assert(placementSlots[i].getIndexInOwner() == i);
        placementSlots[i].setOwner(updMyHandle(), i);
    }
    for (size_t i=0; i < (size_t)getNPlacementValues(); ++i) {
        assert(placementValueSlots[i].hasOwner());
        assert(placementValueSlots[i].getIndexInOwner() == i);
        placementValueSlots[i].setOwner(getMyHandle(), i);
    }
}

// We have just created at newRoot a copy of the tree rooted at oldRoot, and the
// current Subsystem (for which this is the Rep) is a node in the newRoot tree
// (with correct myHandle). However, the 'placement' pointers 
// still retain the values they had in the oldRoot tree; they must be 
// changed to point to the corresponding entities in the newRoot tree.
// If these pointers point outside the oldRoot tree, however, we'll just
// set them to 0 in the newRoot copy.
void SubsystemRep::fixPlacements(const Subsystem& oldRoot, const Subsystem& newRoot) {
    for (size_t i=0; i < (size_t)getNSubsystems(); ++i)
        childSubsystems[i].updRep().fixPlacements(oldRoot, newRoot);    // recurse

    for (size_t i=0; i < (size_t)getNPlacements(); ++i) {
        PlacementSlot& ps = placementSlots[i];
        ps.updPlacement().updRep().repairFeatureReferences(oldRoot,newRoot);
        ps.repairValueReference(oldRoot,newRoot);
    }

    if (FeatureRep::isA(*this))
        FeatureRep::downcast(*this).fixFeaturePlacement(oldRoot,newRoot);
}

// If Placement p's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
// find the corresponding Placement in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const PlacementSlot* 
SubsystemRep::findCorrespondingPlacementSlot
    (const Subsystem& oldRoot, const PlacementSlot& p, const Subsystem& newRoot)
{
    if (!p.hasOwner()) return 0;
    const Subsystem* corrOwner = findCorrespondingSubsystem(oldRoot,p.getOwner(),newRoot);
    if (!corrOwner) return 0;
    assert(corrOwner->hasRep());

    const PlacementSlot* newTreeRef = 
        &corrOwner->getRep().getPlacementSlot(p.getIndexInOwner());
    assert(newTreeRef);
    assert(&newTreeRef->getOwner() == corrOwner);
    assert(newTreeRef->getIndexInOwner() == p.getIndexInOwner());
    return newTreeRef;
}

// If PlacementValue v's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
// find the corresponding PlacementValue in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const PlacementValueSlot* 
SubsystemRep::findCorrespondingPlacementValueSlot
    (const Subsystem& oldRoot, const PlacementValueSlot& v, const Subsystem& newRoot)
{
    if (!v.hasOwner()) return 0;
    const Subsystem* corrOwner = findCorrespondingSubsystem(oldRoot,v.getOwner(),newRoot);
    if (!corrOwner) return 0;
    assert(corrOwner->hasRep());

    const PlacementValueSlot* newTreeRef = 
        &corrOwner->getRep().getPlacementValueSlot(v.getIndexInOwner());
    assert(newTreeRef);
    assert(&newTreeRef->getOwner() == corrOwner);
    assert(newTreeRef->getIndexInOwner() == v.getIndexInOwner());
    return newTreeRef;
}

// For now we'll allow only letters, digits, and underscore in names. Case is retained
// for display but otherwise insignificant.
bool SubsystemRep::isLegalSubsystemName(const std::string& n) {
    if (n.size()==0) return false;
    for (size_t i=0; i<n.size(); ++i)
        if (!(isalnum(n[i]) || n[i]=='_'))
            return false;
    return true;
}

// Take pathname of the form xxx/yyy/zzz, check its validity and optionally
// return as a list of separate subsystem names. We return true if we're successful,
// false if the pathname is malformed in some way. In that case the last segment
// returned will be the one that caused trouble.
bool SubsystemRep::isLegalSubsystemPathname(const std::string& pathname, 
                                            std::vector<std::string>* segments)
{
    std::string t;
    const size_t end = pathname.size();
    size_t nxt = 0;
    if (segments) segments->clear();
    bool foundAtLeastOne = false;
    // for each segment
    while (nxt < end) {
        // for each character of a segment
        while (nxt < end && pathname[nxt] != '/')
            t += pathname[nxt++];
        foundAtLeastOne = true;
        if (segments) segments->push_back(t);
        if (!isLegalSubsystemName(t))
            return false;
        t.clear();
        ++nxt; // skip '/' (or harmless extra increment at end)
    }
    return foundAtLeastOne;
}

} // namespace simtk



static int caseInsensitiveCompare(const std::string& key, const std::string& test) {
    const size_t minlen = std::min(key.size(), test.size());
    for (size_t i=0; i < minlen; ++i) {
        const int k = tolower(key[i]), t = tolower(test[i]);
        if (k < t) return -1;
        else if (k > t) return 1;
    }
    // caution -- size() is unsigned, don't get clever here
    if (key.size() > minlen) return 1;
    else if (test.size() > minlen) return -1;
    return 0;
}

