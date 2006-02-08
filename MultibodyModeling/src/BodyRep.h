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

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Mechanics.h"
#include "simbody/internal/Placement.h"
#include "simbody/internal/BasicPlacements.h"
#include "simbody/internal/Body.h"

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

    const RealMeasure& getMass() const {
        return RealMeasure::downcast(getFeature(massMeasureIndex));
    }
    RealMeasure& updMass() {
        return RealMeasure::downcast(updFeature(massMeasureIndex));
    }
    const StationMeasure& getMassCenter() const {
        return StationMeasure::downcast(getFeature(COMMeasureIndex));
    }
    StationMeasure& updMassCenter() {
        return StationMeasure::downcast(updFeature(COMMeasureIndex));
    }
    const InertiaMeasure& getInertia() const {
        return InertiaMeasure::downcast(getFeature(inertiaMeasureIndex));
    }
    InertiaMeasure& updInertia() {
        return InertiaMeasure::downcast(updFeature(inertiaMeasureIndex));
    }
    const InertiaMeasure& getCentralInertia() const {
        return InertiaMeasure::downcast(getFeature(centralInertiaMeasureIndex));
    }
    InertiaMeasure& updCentralInertia() {
        return InertiaMeasure::downcast(updFeature(centralInertiaMeasureIndex));
    }

    // virtuals getFeatureTypeName() && clone() still missing

    SIMTK_DOWNCAST(BodyRep,SubsystemRep);

protected:
    // Every Body defines some mass-oriented measures.
    virtual void initializeStandardSubfeatures() {
        FrameRep::initializeStandardSubfeatures();

        // Place the body frame on itself.
        place(Frame());

        Feature& mm  = addFeatureLike(RealMeasure("mass"), "mass");
        Feature& cm  = addFeatureLike(StationMeasure("centroid"), "centroid");
        Feature& im  = addFeatureLike(InertiaMeasure("inertia"), "inertia");
        Feature& cim = addFeatureLike(InertiaMeasure("centralInertia"), "centralInertia");
        massMeasureIndex            = mm.getIndexInParent();
        COMMeasureIndex             = cm.getIndexInParent();
        inertiaMeasureIndex         = im.getIndexInParent();
        centralInertiaMeasureIndex  = cim.getIndexInParent();
    }

    virtual void finalizeStandardSubfeatures() {
        FrameRep::finalizeStandardSubfeatures();

        // Add up the masses and place the mass measure on the resulting Placement.
        RealPlacement    totalMass(0.);
        StationPlacement com(Vec3(0));
        InertiaPlacement inertia(InertiaMat(0.));
        for (int i=0; i < getNSubsystems(); ++i) {
            if (MassElement::isInstanceOf(getSubsystem(i))) {
                const MassElement& me = MassElement::downcast(getSubsystem(i));
                totalMass = RealPlacement::downcast   (add(totalMass, me.getMassMeasure()));
                com       = StationPlacement::downcast(add(com,  me.getMassMeasure()*me.getCentroidMeasure()));
                inertia   = InertiaPlacement::downcast(add(inertia,   me.getInertiaMeasure()));
            }
        }
        updMass().replace(totalMass);
        updMassCenter().replace(com / getMass());
        updInertia().replace(inertia);
        updCentralInertia().replace(InertiaPlacement(getInertia())
                                        .shiftToCOM(getMassCenter(),getMass()));
    }

private:
    int massMeasureIndex, COMMeasureIndex, inertiaMeasureIndex, centralInertiaMeasureIndex;
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
    JointRep(Joint& j, Joint::JointType jt, const std::string& nm) 
      : SubsystemRep(j,nm), jointType(jt), refIndex(-1), movIndex(-1) { }
    // must call initializeStandardSubfeatures to complete construction

    std::string   getFeatureTypeName() const { return "Joint"; }

    SubsystemRep* clone() const { return new JointRep(*this); }

    const FrameFeature& getReferenceFrame() const {return FrameFeature::downcast(getFeature(refIndex));}
    const FrameFeature& getMovingFrame()    const {return FrameFeature::downcast(getFeature(movIndex)); }
    Joint::JointType    getJointType()      const {return jointType;}

    SIMTK_DOWNCAST(JointRep,SubsystemRep);
protected:
    virtual void initializeStandardSubfeatures() {
        FrameFeature& R = FrameFeature::downcast(addFeatureLike(FrameFeature("R"), "reference"));
        FrameFeature& M = FrameFeature::downcast(addFeatureLike(FrameFeature("M"), "moving"));

        refIndex = R.getIndexInParent();
        movIndex = M.getIndexInParent();
    }

    Joint::JointType    jointType;
    int                 refIndex, movIndex; // feature indices
};

} // namespace simtk


#endif // SIMTK_BODY_REP_H_
