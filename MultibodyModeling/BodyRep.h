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
#include "BasicFeaturesRep.h"

namespace simtk {

/**
 * This is a still-abstract Feature representation, common to all
 * the Body subsystems.
 */
class BodyRep : public FrameRep {
public:
    BodyRep(Body& m, const std::string& nm) : FrameRep(m,nm) { }
    // must call initializeStandardSubfeatures to complete construction

    const RealMeasure& getMassMeasure() const {
        return RealMeasure::downcast(getFeature(massMeasureIndex));
    }
    RealMeasure& updMassMeasure() {
        return RealMeasure::downcast(updFeature(massMeasureIndex));
    }
    const StationMeasure& getCentroidMeasure() const {
        return StationMeasure::downcast(getFeature(centroidMeasureIndex));
    }
    StationMeasure& updCentroidMeasure() {
        return StationMeasure::downcast(updFeature(centroidMeasureIndex));
    }
    PlacementType getRequiredPlacementType() const { return FramePlacementType; }

    // virtuals getFeatureTypeName() && clone() still missing

    SIMTK_DOWNCAST(BodyRep,SubsystemRep);

protected:
    // Every Body defines some mass-oriented measures.
    virtual void initializeStandardSubfeatures() {
        FrameRep::initializeStandardSubfeatures();

        Feature& mm = addFeatureLike(RealMeasure("massMeasure"), "massMeasure");
        Feature& cm = addFeatureLike(StationMeasure("centroidMeasure"), "centroidMeasure");
        massMeasureIndex     = mm.getIndexInParent();
        centroidMeasureIndex = cm.getIndexInParent();
    }

private:
    int massMeasureIndex, centroidMeasureIndex;
};

class RigidBodyRep : public BodyRep {
public:
    RigidBodyRep(RigidBody& pm, const std::string& nm) : BodyRep(pm,nm) { }
    // must call initializeStandardSubfeatures to complete construction

    std::string   getFeatureTypeName() const { return "RigidBody"; }
    FeatureRep*   clone() const { return new RigidBodyRep(*this); }

    SIMTK_DOWNCAST(RigidBodyRep,SubsystemRep);

protected:
    virtual void initializeStandardSubfeatures() {
        BodyRep::initializeStandardSubfeatures();

    }

    virtual void finalizeStandardSubfeatures() {
        BodyRep::finalizeStandardSubfeatures();

        // Add up the masses and place the mass measure on the resulting Placement.
        Placement totalMass(0.);
        for (int i=0; i < getNSubsystems(); ++i) {
            if (MassElement::isInstanceOf(getSubsystem(i))) {
                const MassElement& me = MassElement::downcast(getSubsystem(i));
                totalMass = add(totalMass, me.getMassMeasure());
            }
        }

        updMassMeasure().place(totalMass);            
    }
};

class MultibodyRep : public SubsystemRep {
public:
    MultibodyRep(Multibody& handle, const std::string& nm)
      : SubsystemRep(handle,nm) { }
    // must call initializeStandardSubfeatures to complete construction

    std::string   getFeatureTypeName() const { return "Multibody"; }
    SubsystemRep* clone() const { return new MultibodyRep(*this); }

    SIMTK_DOWNCAST(MultibodyRep,SubsystemRep);

protected:
    virtual void initializeStandardSubfeatures() {
/* TODO these measures need to have a FunctionPlacement (with runtime states)
        updChildFeature("massMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("mass"))));
        updChildFeature("centroidMeasure")->setPlacement(
            addPlacementLike(FeaturePlacement(*getChildFeature("station"))));
            */
    }
};


class JointRep : public SubsystemRep {
public:
    JointRep(Joint& j, JointType jt, const std::string& nm) 
      : SubsystemRep(j,nm), refIndex(-1), movIndex(-1) { }
    // must call initializeStandardSubfeatures to complete construction

    // no placement for the joint as a whole
    Placement convertToRequiredPlacementType(const Placement& p) const {
        return Placement();
    }

    std::string   getFeatureTypeName() const { return "Joint"; }
    PlacementType getRequiredPlacementType() const { return VoidPlacementType; }
    SubsystemRep* clone() const { return new JointRep(*this); }

    const Frame& getReferenceFrame() const {return Frame::downcast(getFeature(refIndex));}
    const Frame& getMovingFrame()    const {return Frame::downcast(getFeature(movIndex)); }

    SIMTK_DOWNCAST(JointRep,SubsystemRep);
protected:
    virtual void initializeStandardSubfeatures() {
        Frame& R = Frame::downcast(addFeatureLike(Frame("R"), "reference"));
        Frame& M = Frame::downcast(addFeatureLike(Frame("M"), "moving"));

        refIndex = R.getIndexInParent();
        movIndex = M.getIndexInParent();
    }

    int refIndex, movIndex; // feature indices
};

} // namespace simtk


#endif // SIMTK_BODY_REP_H_
