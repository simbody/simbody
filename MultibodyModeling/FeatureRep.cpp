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
 * Implementations of non-inline methods of FeatureRep.
 */

#include "SimbodyCommon.h"
#include "Feature.h"
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

    // SUBSYSTEM REP //

void SubsystemRep::realize(/*State,*/Stage g) const {
    // Always let the children go first.
    for (int i=0; i < getNSubsystems(); ++i)
        getSubsystem(i).realize(/*State,*/ g);

    // For startup stage, we assume that we have recently modified the 
    // model, i.e., features and placements have changed. So we need to
    // (re)calculate the set of PlacementValues held by this Subsystem.
    // Note that we just did this for all the children, so all lower
    // PlacementValues should have been allocated. Nevertheless we have
    // to do this recursively to find all the unassigned Placements whose
    // highest Feature reference is one of the Placement expressions held here.
    if (g == Stage::Startup) {
        deleteAllPlacementValues();
        allocatePlacementValueSlots(getMyHandle());
    }

    // TODO: for now we just assume all our values are invalid. Dependency
    // tracking would be better when costs are high.
    for (int i=0; i < getNPlacementValues(); ++i)
        updPlacementValue(i).updRep().setValid(false);

    // TODO: current method will have side effect of realizing needed
    // PlacementValues as we go, so many of these will be evaluated
    // already by the time we get there.
    for (int i=0; i < getNPlacementValues(); ++i)
        if (!getPlacementValue(i).isValid())
            getPlacementValue(i).getRep().getClientPlacement().realize(/*State,*/g);

    //if (FeatureRep::isA(*this)) {
    //    const FeatureRep& fr = FeatureRep::downcast(*this);
    //    if (fr.hasPlacement()) 
    //        fr.getPlacement().realize(/*State,*/ g);
    //}
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
Placement& 
SubsystemRep::addPlacementLike(const Placement& p) {
    assert(p.hasRep());

    const Feature* offender;
    if (!p.getRep().isLimitedToSubtree(getMyHandle(),offender)) {
        SIMTK_THROW3(Exception::PlacementMustBeLocal,"SubsystemRep::addPlacementLike",
            this->getFullName(),offender->getFullName());
    }

    const int index = (int)placementExpressions.size();
    placementExpressions.push_back(SubPlacement());
    Placement& newPlacement = placementExpressions[index];
    p.getRep().cloneUnownedWithNewHandle(newPlacement);
    newPlacement.updRep().setOwner(getMyHandle(), index);
    return newPlacement;
}

PlacementValue& 
SubsystemRep::addPlacementValueLike(const PlacementValue& v) const {
    assert(v.hasRep());

    // 'placementValues' is mutable.
    const int index = (int)placementValues.size();
    placementValues.push_back(PlacementValue());
    PlacementValue& newPlacementValue = placementValues[index];
    v.getRep().cloneUnownedWithNewHandle(newPlacementValue);
    newPlacementValue.updRep().setOwner(getMyHandle(), index);
    return newPlacementValue;
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

// Is Placement p owned by a Feature in the tree rooted at oldRoot?
/*static*/ bool 
SubsystemRep::isPlacementInSubsystemTree(const Subsystem& oldRoot, const Placement& p)
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
    for (size_t i=0; i<(size_t)getNPlacementExpressions(); ++i) 
        getPlacementExpression(i).checkPlacementConsistency(&getMyHandle(), (int)i, root);
    for (size_t i=0; i < (size_t)getNPlacementValues(); ++i)
        getPlacementValue(i).checkPlacementValueConsistency(&getMyHandle(), (int)i, root);
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
    for (size_t i=0; i < (size_t)getNPlacementExpressions(); ++i) {
        assert(placementExpressions[i].getRep().hasOwner());
        assert(placementExpressions[i].getRep().getIndexInOwner() == i);
        placementExpressions[i].updRep().setOwner(getMyHandle(), i);
    }
    for (size_t i=0; i < (size_t)getNPlacementValues(); ++i) {
        assert(placementValues[i].getRep().hasOwner());
        assert(placementValues[i].getRep().getIndexInOwner() == i);
        placementValues[i].updRep().setOwner(getMyHandle(), i);
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

    for (size_t i=0; i < (size_t)getNPlacementExpressions(); ++i) {
        PlacementRep& pr = placementExpressions[i].updRep();
        pr.repairFeatureReferences(oldRoot,newRoot);
        pr.repairValueReference(oldRoot,newRoot);
    }

    if (FeatureRep::isA(*this))
        FeatureRep::downcast(*this).fixFeaturePlacement(oldRoot,newRoot);
}

// If Placement p's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
// find the corresponding Placement in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const Placement* 
SubsystemRep::findCorrespondingPlacement
    (const Subsystem& oldRoot, const Placement& p, const Subsystem& newRoot)
{
    if (!p.hasOwner()) return 0;
    const Subsystem* corrOwner = findCorrespondingSubsystem(oldRoot,p.getOwner(),newRoot);
    if (!corrOwner) return 0;
    assert(corrOwner->hasRep());

    const Placement* newTreeRef = 
        &corrOwner->getRep().getPlacementExpression(p.getIndexInOwner());
    assert(newTreeRef);
    assert(&newTreeRef->getOwner() == corrOwner);
    assert(newTreeRef->getIndexInOwner() == p.getIndexInOwner());
    return newTreeRef;
}

// If PlacementValue v's owner Subsystem is a member of the Subsystem tree rooted at oldRoot,
// find the corresponding PlacementValue in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const PlacementValue* 
SubsystemRep::findCorrespondingPlacementValue
    (const Subsystem& oldRoot, const PlacementValue& v, const Subsystem& newRoot)
{
    if (!v.hasOwner()) return 0;
    const Subsystem* corrOwner = findCorrespondingSubsystem(oldRoot,v.getOwner(),newRoot);
    if (!corrOwner) return 0;
    assert(corrOwner->hasRep());

    const PlacementValue* newTreeRef = 
        &corrOwner->getRep().getPlacementValue(v.getIndexInOwner());
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

    // FEATURE REP //

// These are default implementations. Derived features which can actually
// be used as a placement of the given type should override.

/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsRealPlacement(RealPlacement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Real");
    //NOTREACHED
    return 0;
}
/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsVec3Placement(Vec3Placement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Vec3");
    //NOTREACHED
    return 0;
}
/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsStationPlacement(StationPlacement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Station");
    //NOTREACHED
    return 0;
}
/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsDirectionPlacement(DirectionPlacement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Direction");
    //NOTREACHED
    return 0;
}    
/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsOrientationPlacement(OrientationPlacement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Orientation");
    //NOTREACHED
    return 0;
} 
/*virtual*/ PlacementRep*
FeatureRep::useFeatureAsFramePlacement(FramePlacement&) const {
    SIMTK_THROW3(Exception::FeatureCantBeUsedAsPlacement,
                    getFullName(), getFeatureTypeName(), "Frame");
    //NOTREACHED
    return 0;
} 


// Use a Placement like p (possibly recast to something else) for this
// feature. Concrete FeatureRep's are responsible for interpreting the
// Placement and possibly converting it to something usable.
//
// We have to decide on an owner feature for the placement expression.
// That is the youngest common ancestor of this feature and all features
// mentioned explicitly in the placement expression.
//
// If this placement is currently evalutable (meaning it is a constant
// or references only features with evaluatable placements) then we
// can allocate a value slot for it in the owner feature. In addition,
// this may have enabled evaluation of any number of additional placements
// which were dependent (directly or indirectly) on the placement of 
// this feature. Value slots for a given placement expression x are always
// owned by the oldest owner of any of the placements on which x recursively
// depends. Note: we defer allocation of value slots until we're done 
// constructing the subsystem (i.e., not here).
void FeatureRep::place(const Placement& p) {
    if (hasPlacement()) {
        SIMTK_THROW1(Exception::FeatureHasAlreadyBeenPlaced,
            getFullName());
        //NOTREACHED  
    }

    assert(p.hasRep());

    // If possible, create a fixed-up copy of p which is suitable for
    // use as a Placement for this concrete FeatureRep.
    Placement pTweaked = 
        p.getRep().getPlacementType() == getRequiredPlacementType()
        ? p : convertToRequiredPlacementType(p);
    if (!pTweaked.hasRep()) {
        SIMTK_THROW3(Exception::PlacementCantBeUsedForThisFeature,
            PlacementRep::getPlacementTypeName(p.getRep().getPlacementType()),
            getFullName(), getFeatureTypeName());
        //NOTREACHED             
    }

    assert(pTweaked.getRep().getPlacementType() == getRequiredPlacementType());
    
    // If the Placement references any features, all its references
    // must be on the same feature tree as this feature (although not necessarily
    // *below* this feature). We will make the placement owner be the youngest
    // common ancestor of this feature and all the features referenced (directly)
    // by the placement. Note that this is not a recursive search through the
    // referenced features' placements -- we only care about direct feature 
    // references, not how they are placed (they may not even have placements
    // yet at all).

    const Feature* offender;
    if (!pTweaked.getRep().isLimitedToSubtree(findRootSubsystem(), offender)) {
        SIMTK_THROW2(Exception::FeatureAndPlacementOnDifferentTrees,
            getFullName(), offender->getFullName());
        //NOTREACHED
    }

    // If the Placement doesn't reference any features, it is a constant
    // value and can be owned by anyone. If the current Feature is a
    // prototype (has no parent) then we are "locking down" a value in
    // the prototype and the current Feature can own the placement itself.
    // If on the other hand the current Feature has a parent, then we
    // want the parent to own the placement (making it external). This is a 
    // significant difference because in the self-placement case the placement 
    // would remain in place after a copy, whereas external placements are
    // removed by copy (or assign). So either this Feature (if alone) or its
    // parent will be the youngest conceivable owner for the new Placement.
    const Subsystem& youngestAllowed = hasParentSubsystem() 
        ? getParentSubsystem() : SubsystemRep::getMyHandle();

    const Subsystem* commonAncestor = 
        pTweaked.getRep().findAncestorSubsystem(youngestAllowed);
    assert(commonAncestor); // there has to be one since they are on the same tree!

    // Please look the other way for a moment while we make a small change to
    // this const Feature ...
    Placement& good = 
        const_cast<Subsystem*>(commonAncestor)->updRep().addPlacementLike(pTweaked);

    // Some sanity (insanity?) checks.
    assert(good.hasOwner());
    assert(good.isConstant() || !good.getOwner().isSameSubsystem(getMyHandle()));
    assert(SubsystemRep::isSubsystemInSubsystemTree(good.getOwner(), getMyHandle()));
    assert(!good.dependsOn(getMyHandle())); // depends on *is* recursive
    placement = &good;
    good.updRep().setClientFeature(getMyHandle());
    postProcessNewPlacement();
}

// This is for use by SubsystemRep after a copy to fix the placement pointer and back pointer.
void FeatureRep::fixFeaturePlacement(const Subsystem& oldRoot, const Subsystem& newRoot) {
    if (placement) {
        placement = findCorrespondingPlacement(oldRoot,*placement,newRoot);
        if (placement) { 
            const_cast<Placement*>(placement)->updRep().setClientFeature(getMyHandle());
            if (placement->getRep().hasValueSlot())
                const_cast<Placement*>(placement)->updRep().updValueSlot().updRep()
                                                                .setClientPlacement(*placement);
        }
    }
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

