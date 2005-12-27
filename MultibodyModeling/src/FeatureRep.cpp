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

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Feature.h"

#include "SubsystemRep.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::cout;
using std::endl;
using std::ostream;

namespace simtk {

    // FEATURE REP //

// Use a Placement like p (possibly recast to something else) for this
// feature. Concrete FeatureRep's are responsible for interpreting the
// Placement and possibly converting it to something usable.
//
// We have to decide on an owner Subsystem for the placement expression.
// That is the youngest common ancestor of this feature and all features
// mentioned explicitly in the placement expression.
//
// Note: we defer allocation of value slots until we're done 
// constructing the subsystem (i.e., not here).
void FeatureRep::place(const Placement& p) {
    if (hasPlacement()) {
        SIMTK_THROW1(Exception::FeatureHasAlreadyBeenPlaced,
            getFullName());
        //NOTREACHED  
    }
    replace(p);
}

// If this Feature already has a placement, we have to delete it and
// wipe out everything that might depend on it (for now, all placement
// values in the current feature, its parent and ancestors).
//
// If the new placement should be held by the same Subsystem as the
// current one, we'll re-use the same slot. Otherwise we'll mark the
// current slot unused and get a new one somewhere else.
void FeatureRep::replace(const Placement& p) {
    assert(p.hasRep());

    // If possible, create a fixed-up copy of p which is suitable for
    // use as a Placement for this concrete FeatureRep.

    Placement pTweaked;
    const Placement* sourcePlacement = 0;
    if (getSamplePlacement().getRep().isSamePlacementType(p))
        sourcePlacement = &p;
    else {
        // This will throw an exception if it doesn't work.
        pTweaked = convertToRequiredPlacementType(p);
        sourcePlacement = &pTweaked;
    }

    assert(isRequiredPlacementType(*sourcePlacement));

    // convenient abbreviations
    const Placement&    src    = *sourcePlacement;
    const PlacementRep& srcRep = src.getRep();
    
    // If the Placement references any features, all its references
    // must be on the same feature tree as this feature (although not necessarily
    // *below* this feature). We will make the placement owner be the youngest
    // common ancestor of this feature and all the features referenced (directly)
    // by the placement. Note that this is not a recursive search through the
    // referenced features' placements -- we only care about direct feature 
    // references, not how they are placed (they may not even have placements
    // yet at all).

    const Feature* offender;
    if (!srcRep.isLimitedToSubtree(findRootSubsystem(), offender)) {
        SIMTK_THROW2(Exception::FeatureAndPlacementOnDifferentTrees,
            getFullName(), offender->getFullName());
        //NOTREACHED
    }

    // The Placement cannot depend directly or indirectly on this 
    // Feature or we have a very awkward dependency of a Feature on
    // itself. This *might* be solvable by iteration, but we certainly
    // aren't going to attempt that! Much more likely, this is an
    // error on the part of the person putting the model together.
    if (srcRep.dependsOn(getMyHandle())) {
        SIMTK_THROW1(Exception::AFeatureCantBePlacedOnItself,
            getFullName());
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
        srcRep.findAncestorSubsystem(youngestAllowed);
    assert(commonAncestor); // there has to be one since they are on the same tree!

    // Now we have to decide what to do with the current placement's slot (if any).
    // If its owner subsystem is the same as commonAncestor we just calculated,
    // then we want to keep the same slot. Otherwise we'll delete the old one
    // and assign a new one.

    if (hasPlacement()) {
        if (commonAncestor->isSameSubsystem(getPlacementSlot().getOwner()))
            updPlacementSlot().updPlacement() = src;
        else
            removePlacement();
    }

    if (!hasPlacement()) {
        // Please look the other way for a moment while we make a small change to
        // this const Feature ...
        PlacementSlot& newSlot = 
            const_cast<Subsystem*>(commonAncestor)->updRep().addPlacementLike(src);
        placement = &newSlot;
        newSlot.setClientFeature(updMyHandle());
    }

    // Some sanity (insanity?) checks.
    assert(hasPlacement());
    const Subsystem& owner = getPlacementSlot().getOwner();
    assert(getPlacement().isConstant() || !owner.isSameSubsystem(getMyHandle()));
    assert(SubsystemRep::isSubsystemInSubsystemTree(owner, getMyHandle()));

    postProcessNewPlacement();
}

void FeatureRep::removePlacement() {
    if (!hasPlacement()) return;
    Subsystem& owner = updPlacementSlot().updOwner();
    owner.updRep().deletePlacementSlot(getPlacementSlot().getIndexInOwner());
    // that should have cleared our reference here
    assert(!hasPlacement());
}

// This is for use by SubsystemRep after a copy to fix the placement pointer and back pointer.
void FeatureRep::fixFeaturePlacement(const Subsystem& oldRoot, const Subsystem& newRoot) {
    if (placement) {
        placement = findCorrespondingPlacementSlot(oldRoot,*placement,newRoot);
        if (placement) 
            const_cast<PlacementSlot*>(placement)->setClientFeature(updMyHandle());
    }
}

} // namespace simtk

