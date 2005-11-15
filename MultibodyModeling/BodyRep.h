#ifndef SIMTK_BODY_REP_H_
#define SIMTK_BODY_REP_H_

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
 * Declarations for the *real* Body objects. These are opaque to
 * users.
 */

#include "Body.h"
#include "FeatureRep.h"

namespace simtk {

/**
 * This is a still-abstract Feature representation, common to all
 * the Body features.
 */
class BodyRep : public FrameRep {
public:
    BodyRep(Body& m, const std::string& nm) : FrameRep(m,nm) { 
        initializeSubfeatures();
    }

    const RealMeasure& getMassMeasure() const {
        return RealMeasure::downcast(getSubfeature(massMeasureIndex));
    }
    const StationMeasure& getCentroidMeasure() const {
        return StationMeasure::downcast(getSubfeature(centroidMeasureIndex));
    }

    // virtuals getFeatureTypeName() && clone() still missing

    SIMTK_DOWNCAST2(BodyRep,FrameRep,FeatureRep);
private:
    // Every Body defines some mass-oriented measures.
    void initializeSubfeatures() {
        Feature& mm = addSubfeatureLike(RealMeasure("massMeasure"), "massMeasure");
        Feature& cm = addSubfeatureLike(StationMeasure("centroidMeasure"), "centroidMeasure");
        massMeasureIndex     = mm.getIndexInParent();
        centroidMeasureIndex = cm.getIndexInParent();
    }

    int massMeasureIndex, centroidMeasureIndex;
};

class RigidBodyRep : public BodyRep {
public:
    RigidBodyRep(RigidBody& pm, const std::string& nm)
      : BodyRep(pm,nm) { }

    std::string   getFeatureTypeName() const { return "RigidBody"; }
    PlacementType getRequiredPlacementType() const { return InvalidPlacementType; }
    FeatureRep*   clone() const { return new RigidBodyRep(*this); }

    SIMTK_DOWNCAST2(RigidBodyRep,BodyRep,FeatureRep);
private:
    void initializeFeatures() {
        addSubfeatureLike(RealParameter("mass"), "mass");
        addSubfeatureLike(Station("station"), "station");
/* TODO these measures need to have a FunctionPlacement
        updChildFeature("massMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("mass"))));
        updChildFeature("centroidMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("station"))));
            */
    }
};



} // namespace simtk


#endif // SIMTK_BODY_REP_H_
