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
#include "SubsystemRep.h"
#include "Feature.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::cout;
using std::endl;
using std::ostream;

namespace simtk {

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
// Note: we defer allocation of value slots until we're done 
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
    PlacementSlot& good = 
        const_cast<Subsystem*>(commonAncestor)->updRep().addPlacementLike(pTweaked);

    // Some sanity (insanity?) checks.
    assert(good.hasOwner());
    assert(good.getPlacement().isConstant() || !good.getOwner().isSameSubsystem(getMyHandle()));
    assert(SubsystemRep::isSubsystemInSubsystemTree(good.getOwner(), getMyHandle()));
    assert(!good.getPlacement().dependsOn(getMyHandle())); // depends on *is* recursive
    placement = &good;
    good.setClientFeature(getMyHandle());
    postProcessNewPlacement();
}

// This is for use by SubsystemRep after a copy to fix the placement pointer and back pointer.
void FeatureRep::fixFeaturePlacement(const Subsystem& oldRoot, const Subsystem& newRoot) {
    if (placement) {
        placement = findCorrespondingPlacementSlot(oldRoot,*placement,newRoot);
        if (placement) 
            const_cast<PlacementSlot*>(placement)->setClientFeature(getMyHandle());
    }
}

} // namespace simtk

