#ifndef SIMTK_MASS_ELEMENT_REP_H_
#define SIMTK_MASS_ELEMENT_REP_H_

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
 * Declarations for the *real* MassElement objects. These are opaque to
 * users.
 */

#include "MassElement.h"
#include "FeatureRep.h"

namespace simtk {

/**
 * This is a still-abstract Feature representation, common to all
 * the MassElement features.
 */
class MassElementRep : public FeatureRep {
public:
    MassElementRep(MassElement& m, const std::string& nm) : FeatureRep(m,nm) { 
        initializeFeatures();
    }

    const RealMeasure& getMassMeasure() const {
        return RealMeasure::downcast(getChildFeature(0));
    }
    const StationMeasure& getCentroidMeasure() const {
        return StationMeasure::downcast(getChildFeature(1));
    }

    // virtuals getFeatureTypeName() && cloneWithoutPlacement() still missing

    SIMTK_DOWNCAST(MassElementRep,FeatureRep);
private:
    // Every MassElement defines some mass-oriented measures.
    void initializeFeatures() {
        addFeatureLike(RealMeasure("massMeasure"), "massMeasure");
        addFeatureLike(StationMeasure("centroidMeasure"), "centroidMeasure");
    }
};

class PointMassElementRep : public MassElementRep {
public:
    PointMassElementRep(PointMassElement& pm, const std::string& nm) 
      : MassElementRep(pm,nm) { 
        initializeFeatures();
    }

    // some self-placements
    void setMass(const Real& m) {
        const Placement& p = addPlacementLike(RealPlacement(m));
        updChildFeature("mass")->setPlacement(p);
    }

    void placePoint(const Vec3& v) {
        const Placement& p = addPlacementLike(StationPlacement(v));
        updChildFeature("station")->setPlacement(p);
    }

    std::string getFeatureTypeName() const { return "PointMassElement"; }
    // PointMass has subfeatures that need Placement, but does not
    // have its own Placement.
    PlacementType getRequiredPlacementType() const { return InvalidPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        PointMassElementRep* copy = new PointMassElementRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST2(PointMassElementRep,MassElementRep,FeatureRep);
private:
    void initializeFeatures() {
        addFeatureLike(RealParameter("mass"), "mass");
        addFeatureLike(Station("station"), "station");

        updChildFeature("massMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("mass"))));
        updChildFeature("centroidMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("station"))));
    }
};

class CylinderMassElementRep : public MassElementRep {
public:
    CylinderMassElementRep(CylinderMassElement& cm, const std::string& nm) 
      : MassElementRep(cm,nm) { 
        initializeFeatures();
    }

    // some self-placements
    void setMass(const Real& m) {
        const Placement& p = addPlacementLike(RealPlacement(m));
        updChildFeature("mass")->setPlacement(p);
    }
    void setRadius(const Real& r) {
        const Placement& p = addPlacementLike(RealPlacement(r));
        updChildFeature("radius")->setPlacement(p);
    }
    void setHalfLength(const Real& h) {
        const Placement& p = addPlacementLike(RealPlacement(h));
        updChildFeature("halfLength")->setPlacement(p);
    }
    void placeCenter(const Vec3& c) {
        const Placement& p = addPlacementLike(StationPlacement(c));
        updChildFeature("center")->setPlacement(p);
    }
    void placeAxis(const Vec3& a) {
        const Placement& p = addPlacementLike(DirectionPlacement(a));
        updChildFeature("axis")->setPlacement(p);
    }

    std::string getFeatureTypeName() const { return "CylinderMassElement"; }

    // Cylinder has subfeatures that need Placement, but does not
    // have its own Placement.
    PlacementType getRequiredPlacementType() const { return InvalidPlacementType; }
    FeatureRep* cloneWithoutPlacement(Feature& f) const {
        CylinderMassElementRep* copy = new CylinderMassElementRep(*this);
        copy->cleanUpAfterClone(f);
        return copy;
    }

    SIMTK_DOWNCAST2(CylinderMassElementRep,MassElementRep,FeatureRep);
private:
    void initializeFeatures() {
        addFeatureLike(RealParameter("mass"),       "mass");
        addFeatureLike(RealParameter("radius"),     "radius");
        addFeatureLike(RealParameter("halfLength"), "halfLength");
        addFeatureLike(Station      ("center"),     "center");
        addFeatureLike(Direction    ("axis"),       "axis");

        updChildFeature("massMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("mass"))));
        updChildFeature("centroidMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("center"))));
    }
};

} // namespace simtk


#endif // SIMTK_MASS_ELEMENT_REP_H_
