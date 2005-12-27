#ifndef SIMTK_FEATURE_REP_H_
#define SIMTK_FEATURE_REP_H_

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
 * Declarations for the *real* Multibody Modeling objects. These are opaque to
 * users.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Feature.h"
#include "simbody/internal/Placement.h"

#include "SubsystemRep.h"
#include "PlacementRep.h"

#include <string>
#include <cassert>
#include <sstream>
#include <cctype>

namespace simtk {

/**
 * FeatureRep is a still-abstract SubsystemRep which adds handling of the Feature's
 * placement to the basic SubsystemRep capabilities.
 */
class FeatureRep : public SubsystemRep {
public:
    FeatureRep(Feature& p, const std::string& nm, const Placement& placementSample)
        : SubsystemRep(p,nm), placement(0), sample(placementSample) { }
    virtual ~FeatureRep() { }

    // let's be more precise
    const Feature& getMyHandle() const
      { return reinterpret_cast<const Feature&>(SubsystemRep::getMyHandle()); }
    Feature&       updMyHandle()
      { return reinterpret_cast<Feature&>(SubsystemRep::updMyHandle()); }

    // This routine offers control after the feature has
    // been placed (that doesn't mean you can necessarily get a *value* for
    // that placement; just the expression defining that value).
    virtual void postProcessNewPlacement() { }

    // These allow the feature to weigh in on the suitability of a proposed
    // placement for the feature.

    // Generate a sample Placement of the correct type for this Feature.
    // Methods of this Placement can be used to probe the suitability of
    // proposed Placements.
    const Placement& getSamplePlacement() const {return sample;}
    bool isRequiredPlacementType(const Placement& p) const {
        return getSamplePlacement().getRep().isSamePlacementType(p);
    }

    // Given a proposed placement for this Feature which may not be of the 
    // required type, use it to generate a placement which is of the
    // right type. Throws an exception if it doesn't work.
    Placement convertToRequiredPlacementType(const Placement& p) const
      { return Placement(getSamplePlacement().getRep().createPlacementFrom(p)); }

    virtual std::string   getFeatureTypeName()            const = 0;

    // Create the appropriate concrete PlacementRep for a reference to the 
    // Placement of this kind of Feature, or to one of its Placement elements
    // if we're given an index (-1 means the whole Placement).
    virtual PlacementRep* createFeatureReference(Placement&, int i = -1) const = 0;

    bool             hasPlacement() const {return placement != 0;}

    const PlacementSlot& getPlacementSlot() const {
        if (!placement) 
            SIMTK_THROW1(Exception::RepLevelException, 
            "Feature has no placement");
        return *placement;
    }
    const Placement& getPlacement() const {
        return getPlacementSlot().getPlacement();
    }
    PlacementSlot& updPlacementSlot() {
        return const_cast<PlacementSlot&>(getPlacementSlot());
    }

    // Someone is deleting our placement. Erase the pointer.
    void clearPlacementSlot() {
        placement=0;
    }

    void place(const Placement& p);
    void replace(const Placement& p);
    void removePlacement();

    // This is for use by SubsystemRep after a copy to fix the placement pointer.
    void fixFeaturePlacement(const Subsystem& oldRoot, const Subsystem& newRoot);

    SIMTK_DOWNCAST(FeatureRep, SubsystemRep);
private:
    // If this Feature has been placed, this is the placement information.
    // If present, this PlacementSlot must be owned by this Feature, its parent
    // Subsystem or one of its ancestors.
    const PlacementSlot* placement;

    const Placement sample;
};

} // namespace simtk

#endif // SIMTK_FEATURE_REP_H_
