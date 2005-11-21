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
using std::endl;
using std::ostream;


// Returns -1, 0, 1 according to key {<,==,>} test ignoring case.
static int caseInsensitiveCompare(const std::string& key, const std::string& test);

namespace simtk {

    // FEATURE REP //

std::string 
FeatureRep::getFullName() const { 
    std::string s;
    if (hasParentFeature())
        s = getParentFeature().getFullName() + "/";
    return s + getName(); 
}

const Feature& FeatureRep::findRootFeature() const {
    return hasParentFeature() 
        ? getParentFeature().getRep().findRootFeature() : getMyHandle();
}
Feature& FeatureRep::findUpdRootFeature() {
    return hasParentFeature() 
        ? updParentFeature().updRep().findUpdRootFeature() : updMyHandle();
}

void 
FeatureRep::cloneWithoutParentOrExternalPlacements(Feature& newHandle) const
{
    FeatureRep* copy = clone();
    copy->setMyHandle(newHandle);
    copy->parent = 0; copy->indexInParent = -1;
    newHandle.setRep(copy);

    // Re-parent all the copied child Features to their new parent,
    // and fix the owned Placements to acknowledge their new owner.
    copy->reparentMyChildren();

    // Fix up all the internal placement references and delete the
    // external ones.
    copy->fixPlacements(this->getMyHandle(), copy->getMyHandle());
}

// Use a Placement like p (possibly recast to something else) for this
// feature. Concrete FeatureRep's are responsible for interpreting the
// Placement and possibly converting it to something usable.
void FeatureRep::place(const Placement& p) {
    assert(p.hasRep());

    // If possible, create a fixed-up copy of p which is suitable for
    // use as a Placement for this concrete FeatureRep.
    Placement pTweaked = 
        p.getRep().getPlacementType() == getRequiredPlacementType()
        ? p : recastPlacement(p);
    if (!pTweaked.hasRep()) {
        SIMTK_THROW3(Exception::PlacementCantBeUsedForThisFeature,
            PlacementRep::getPlacementTypeName(p.getRep().getPlacementType()),
            getFullName(), getFeatureTypeName());
        //NOTREACHED             
    }

    assert(pTweaked.getRep().getPlacementType() == getRequiredPlacementType());

    // If the Placement is a constant, it is a self-placement if this 
    // Feature has no parent; otherwise, we'll put the constant Placement
    // in the parent (making it external). This is a significant difference
    // because in the self-placement case the placement remains after a
    // copy, whereas external placements are removed by copy (or assign).
    if (pTweaked.isConstant()) {
        FeatureRep& ownerFeatureRep =
            hasParentFeature() ? updParentFeature().updRep() : *this;
        placement = &(ownerFeatureRep.addPlacementLike(pTweaked));
        return;
    } 
    
    // The placement is not a constant. In that case all its feature references
    // must be on the same feature tree as this feature (although not necessarily
    // *below* this feature). We will make the placement owner be the youngest
    // common ancestor of this feature and all the features referenced by the
    // placement.

    const Feature& myRoot = findRootFeature();
    const Feature* offender;
    if (!pTweaked.getRep().isLimitedToSubtree(myRoot, offender)) {
        SIMTK_THROW2(Exception::FeatureAndPlacementOnDifferentTrees,
            getFullName(), offender->getFullName());
        //NOTREACHED
    }

    const Feature& placementAncestor = 
        *pTweaked.getRep().findAncestorFeature(myRoot);

    Feature* commonAncestor = 
        findUpdYoungestCommonAncestor(updMyHandle(), placementAncestor);
    assert(commonAncestor); // there has to be one since they are on the same tree!
    Placement& good = commonAncestor->updRep().addPlacementLike(pTweaked);

    // Some sanity (insanity?) checks.
    assert(good.hasOwner());
    assert(!good.getOwner().isSameFeature(getMyHandle()));
    assert(FeatureRep::isFeatureInFeatureTree(good.getOwner(), getMyHandle()));
    assert(!good.dependsOn(getMyHandle()));
    placement = &good; 
}

Feature& 
FeatureRep::addSubfeatureLike(const Feature& f, const std::string& nm) {
    assert(nm.size() > 0);
    const int index = (int)subfeatures.size();
    subfeatures.push_back(SubFeature()); // an empty handle
    Feature& newFeature = subfeatures[index];
    f.getRep().cloneWithoutParentOrExternalPlacements(newFeature);
    newFeature.updRep().setParentFeature(updMyHandle(), index);
    newFeature.updRep().setName(nm);
    return newFeature;
}

// Note that we can only allow placements involving this feature, its children,
// grandchildren, etc. -- no external references. Otherwise someone further
// up the tree should own the new placement.
Placement& 
FeatureRep::addPlacementLike(const Placement& p) {
    assert(p.hasRep());

    const Feature* offender;
    if (!p.getRep().isLimitedToSubtree(getMyHandle(),offender)) {
        SIMTK_THROW3(Exception::PlacementMustBeLocal,"FeatureRep::addPlacementLike",
            this->getFullName(),offender->getFullName());
    }

    const int index = (int)placementExpressions.size();
    placementExpressions.push_back(SubPlacement());
    Placement& newPlacement = placementExpressions[index];
    p.getRep().cloneUnownedWithNewHandle(newPlacement);
    newPlacement.updRep().setOwner(getMyHandle(), index);
    return newPlacement;
}


// Is Feature f in the tree rooted at oldRoot? If so, optionally return the 
// series of indices required to get to this Feature from the root.
// Complexity is O(log n) where n is tree depth.
/*static*/ bool 
FeatureRep::isFeatureInFeatureTree(const Feature& oldRoot, const Feature& f,
                                std::vector<int>* trace)
{
    if (trace) trace->clear();
    const Feature* const oldp = &oldRoot;
    const Feature*       fp   = &f;

    while (fp != oldp) {
        const Feature* const fpParent = getParentPtr(*fp);
        if (!fpParent) {
            if (trace) trace->clear(); // never mind ...
            return false;
        }
        if (trace) trace->push_back(fp->rep->getIndexInParent());
        fp = fpParent;
    }

    return true;
}

// Is Placement p owned by a Feature in the tree rooted at oldRoot?
/*static*/ bool 
FeatureRep::isPlacementInFeatureTree(const Feature& oldRoot, const Placement& p)
{
    if (!p.hasOwner())
        return false;   // a disembodied Placement
    return isFeatureInFeatureTree(oldRoot, p.getOwner());
}

// If Feature f is a member of the Feature tree rooted at oldRoot, find
// the corresponding Feature in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const Feature* 
FeatureRep::findCorrespondingFeature
    (const Feature& oldRoot, const Feature& f, const Feature& newRoot)
{
    std::vector<int> trace;
    if (!isFeatureInFeatureTree(oldRoot,f,&trace))
        return 0;

    // Trace holds the indices needed to step from newRoot down to
    // the corresponding Feature (in reverse order).
    const Feature* newTreeRef = &newRoot;
    for (size_t i=trace.size(); i >=1; --i)
        newTreeRef = &newTreeRef->getRep().getSubfeature(trace[i-1]);
    return newTreeRef;
}

// Given two features, run up the tree towards the root to find
// their "least common denominator", i.e. the first shared node
// on the path back to the root. Return a pointer to that node
// if found, otherwise NULL meaning that the features aren't on
// the same tree. If the features are the same, then
// that feature is the answer.
// Complexity is O(log n) (3 passes) where n is depth of Feature tree.

/*static*/ const Feature* 
FeatureRep::findYoungestCommonAncestor(const Feature& f1, const Feature& f2)
{
    std::vector<const Feature*> f1path, f2path; // paths from nodes to their roots
    const Feature* f1p = &f1;
    const Feature* f2p = &f2;
    while (f1p) {f1path.push_back(f1p); f1p = getParentPtr(*f1p);}
    while (f2p) {f2path.push_back(f2p); f2p = getParentPtr(*f2p);}

    // If there is a common ancestor, we can find it by searching down from
    // the root (last element in each path). As soon as there is a difference,
    // the previous element is the common ancestor.
    const Feature* ancestor = 0;
    size_t i1 = f1path.size(), i2 = f2path.size();  // index of element just past end
    while (i1 && i2) {
        if (f1path[--i1] != f2path[--i2])
            break;
        ancestor = f1path[i1];
    }
    return ancestor;
}

/*static*/ Feature* 
FeatureRep::findUpdYoungestCommonAncestor(Feature& f1, const Feature& f2) {
    return const_cast<Feature*>(findYoungestCommonAncestor(f1,f2));
}

// Return true and ix==feature index if a feature of the given name is found.
// Otherwise return false and ix==childFeatures.size().
bool 
FeatureRep::findSubfeatureIndex(const std::string& nm, size_t& ix) const {
    for (ix=0; ix < (size_t)getNSubfeatures(); ++ix)
        if (caseInsensitiveCompare(nm, subfeatures[ix].getName())==0)
            return true;
    return false;   // not found
}

// We have just copied a Feature subtree so all the parent pointers are
// wrong. Recursively repair them to point into the new tree.
void FeatureRep::reparentMyChildren() {
    for (size_t i=0; i < (size_t)getNSubfeatures(); ++i) {
        assert(subfeatures[i].getRep().hasParentFeature());
        assert(subfeatures[i].getRep().getIndexInParent() == i);    // shouldn't change
        subfeatures[i].updRep().setParentFeature(updMyHandle(), i);
        subfeatures[i].updRep().reparentMyChildren();               // recurse
    }
    for (size_t i=0; i < (size_t)getNPlacementExpressions(); ++i) {
        assert(placementExpressions[i].getRep().hasOwner());
        assert(placementExpressions[i].getRep().getIndexInOwner() == i);
        placementExpressions[i].updRep().setOwner(getMyHandle(), i);
    }
}

// We have just created at newRoot a copy of the tree rooted at oldRoot, and the
// current Feature (for which this is the Rep) is a node in the newRoot tree
// (with correct myHandle). However, the 'placement' pointers 
// still retain the values they had in the oldRoot tree; they must be 
// changed to point to the corresponding entities in the newRoot tree.
// If these pointers point outside the oldRoot tree, however, we'll just
// set them to 0 in the newRoot copy.
void FeatureRep::fixPlacements(const Feature& oldRoot, const Feature& newRoot) {
    for (size_t i=0; i < (size_t)getNSubfeatures(); ++i)
        subfeatures[i].updRep().fixPlacements(oldRoot, newRoot);    // recurse

    for (size_t i=0; i < (size_t)getNPlacementExpressions(); ++i)
        placementExpressions[i].updRep().repairFeatureReferences(oldRoot,newRoot);

    if (placement)
        placement = findCorrespondingPlacement(oldRoot,*placement,newRoot);
}

// If Placement p's owner Feature is a member of the Feature tree rooted at oldRoot,
// find the corresponding Placement in the tree rooted at newRoot (which is expected
// to be a copy of oldRoot). Return NULL if not found for any reason.
/*static*/ const Placement* 
FeatureRep::findCorrespondingPlacement
    (const Feature& oldRoot, const Placement& p, const Feature& newRoot)
{
    if (!p.hasOwner()) return 0;
    const Feature* corrOwner = findCorrespondingFeature(oldRoot,p.getOwner(),newRoot);
    if (!corrOwner) return 0;
    assert(corrOwner->hasRep());

    const Placement* newTreeRef = 
        &corrOwner->getRep().getPlacementExpression(p.getIndexInOwner());
    assert(newTreeRef);
    assert(&newTreeRef->getOwner() == corrOwner);
    assert(newTreeRef->getIndexInOwner() == p.getIndexInOwner());
    return newTreeRef;
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
