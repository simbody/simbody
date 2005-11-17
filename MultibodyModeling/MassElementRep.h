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
    MassElementRep(MassElement& m, const std::string& nm) : FeatureRep(m,nm) { }
    // must call initializeStandardSubfeatures() to complete construction.


    const RealMeasure& getMassMeasure() const {
        return RealMeasure::downcast(getSubfeature(massMeasureIndex));
    }
    const StationMeasure& getCentroidMeasure() const {
        return StationMeasure::downcast(getSubfeature(centroidMeasureIndex));
    }

    // virtuals getFeatureTypeName() && clone() still missing

    SIMTK_DOWNCAST(MassElementRep,FeatureRep);

protected:
    // Every MassElement defines some mass-oriented measures.
    virtual void initializeStandardSubfeatures() {
        Feature& mm = addSubfeatureLike(RealMeasure("massMeasure"), "massMeasure");
        Feature& cm = addSubfeatureLike(StationMeasure("centroidMeasure"), "centroidMeasure");
        massMeasureIndex     = mm.getIndexInParent();
        centroidMeasureIndex = cm.getIndexInParent();
    }

private:
    int massMeasureIndex, centroidMeasureIndex;
};

class PointMassElementRep : public MassElementRep {
public:
    PointMassElementRep(PointMassElement& pm, const std::string& nm) 
      : MassElementRep(pm,nm) { }
    // must call initializeStandardSubfeatures() to complete construction.

    Placement recastPlacement(const Placement& p) const {
        return p.getRep().castToStationPlacement();
    }

    // some self-placements
    void setMass(const Real& m) {
        const Placement& p = addPlacementLike(RealPlacement(m));
        findUpdSubfeature("mass")->place(p);
    }

    std::string getFeatureTypeName() const { return "PointMassElement"; }
    PlacementType getRequiredPlacementType() const { return StationPlacementType; }
    FeatureRep* clone() const { return new PointMassElementRep(*this); }

    SIMTK_DOWNCAST2(PointMassElementRep,MassElementRep,FeatureRep);

protected:
    // This will be called from initializeStandardSubfeatures() in MassElement.
    virtual void initializeStandardSubfeatures() {
        MassElementRep::initializeStandardSubfeatures();

        addSubfeatureLike(RealParameter("mass"), "mass");

        findUpdSubfeature("massMeasure")->place(
            addPlacementLike(FeaturePlacement(*findSubfeature("mass"))));
        findUpdSubfeature("centroidMeasure")->place(
            addPlacementLike(FeaturePlacement(getMyHandle())));
    }
};

class CylinderMassElementRep : public MassElementRep {
public:
    CylinderMassElementRep(CylinderMassElement& cm, const std::string& nm) 
      : MassElementRep(cm,nm) { }
    // must call initializeStandardSubfeatures() to complete construction.

    // no placement for the cylinder as a whole
    Placement recastPlacement(const Placement& p) const {
        return Placement();
    }

    // some self-placements
    void setMass(const Real& m) {
        const Placement& p = addPlacementLike(RealPlacement(m));
        findUpdSubfeature("mass")->place(p);
    }
    void setRadius(const Real& r) {
        const Placement& p = addPlacementLike(RealPlacement(r));
        findUpdSubfeature("radius")->place(p);
    }
    void setHalfLength(const Real& h) {
        const Placement& p = addPlacementLike(RealPlacement(h));
        findUpdSubfeature("halfLength")->place(p);
    }
    void placeCenter(const Vec3& c) {
        const Placement& p = addPlacementLike(StationPlacement(c));
        findUpdSubfeature("center")->place(p);
    }
    void placeAxis(const Vec3& a) {
        const Placement& p = addPlacementLike(DirectionPlacement(a));
        findUpdSubfeature("axis")->place(p);
    }

    std::string getFeatureTypeName() const { return "CylinderMassElement"; }

    // Cylinder has subfeatures that need Placement, but does not
    // have its own Placement.
    PlacementType getRequiredPlacementType() const { return InvalidPlacementType; }
    FeatureRep* clone() const { return new CylinderMassElementRep(*this); }

    SIMTK_DOWNCAST2(CylinderMassElementRep,MassElementRep,FeatureRep);
private:
    // This will be called from initializeStandardSubfeatures() in MassElement.
    virtual void initializeStandardSubfeatures() {
        MassElementRep::initializeStandardSubfeatures();

        addSubfeatureLike(RealParameter("mass"),       "mass");
        addSubfeatureLike(RealParameter("radius"),     "radius");
        addSubfeatureLike(RealParameter("halfLength"), "halfLength");
        addSubfeatureLike(Station      ("center"),     "center");
        addSubfeatureLike(Direction    ("axis"),       "axis");

        findUpdSubfeature("massMeasure")->place(
            addPlacementLike(FeaturePlacement(*findSubfeature("mass"))));
        findUpdSubfeature("centroidMeasure")->place(
            addPlacementLike(FeaturePlacement(*findSubfeature("center"))));
    }
};

} // namespace simtk


#endif // SIMTK_MASS_ELEMENT_REP_H_
