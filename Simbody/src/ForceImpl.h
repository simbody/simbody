#ifndef SimTK_SIMBODY_FORCE_IMPL_H_
#define SimTK_SIMBODY_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors: Nabeel Allana                                                *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/Force_BuiltIns.h"

namespace SimTK {

// This is what a Force handle points to.
class ForceImpl : public PIMPLImplementation<Force, ForceImpl> {
public:
    ForceImpl() : forces(0), defaultDisabled(false) {}
    ForceImpl(const ForceImpl& clone) {*this = clone;}

    void setDisabledByDefault(bool shouldBeDisabled) 
    {   invalidateTopologyCache();
        defaultDisabled = shouldBeDisabled; }

    bool isDisabledByDefault() const 
    {   return defaultDisabled; }

    virtual ~ForceImpl() {}
    virtual ForceImpl* clone() const = 0;
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    virtual bool shouldBeParallelized() const { 
        return false;
    }
    ForceIndex getForceIndex() const {return index;}
    const GeneralForceSubsystem& getForceSubsystem() const 
    {   assert(forces); return *forces; }
    void setForceSubsystem(GeneralForceSubsystem& frcsub, ForceIndex ix) {
        forces = &frcsub;
        index  = ix;
    }
    void invalidateTopologyCache() const {
        if (forces) forces->invalidateSubsystemTopologyCache();
    }

    // Every force element must provide the next two methods. Note that 
    // calcForce() must *add in* (+=) its forces to the given arrays.
    virtual void calcForce(const State&         state, 
                           Vector_<SpatialVec>& bodyForces, 
                           Vector_<Vec3>&       particleForces, 
                           Vector&              mobilityForces) const = 0;
    virtual Real calcPotentialEnergy(const State& state) const = 0;

    virtual void realizeTopology    (State& state) const {}
    virtual void realizeModel       (State& state) const {}
    virtual void realizeInstance    (const State& state) const {}
    virtual void realizeTime        (const State& state) const {}
    virtual void realizePosition    (const State& state) const {}
    virtual void realizeVelocity    (const State& state) const {}
    virtual void realizeDynamics    (const State& state) const {}
    virtual void realizeAcceleration(const State& state) const {}
    virtual void realizeReport      (const State& state) const {}

    virtual void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const {}

private:
        // CONSTRUCTION
    GeneralForceSubsystem* forces;  // just a reference; no delete on destruction
    ForceIndex             index;

        // TOPOLOGY "STATE"
    // Changing anything here invalidates the topology of the containing
    // force Subsystem and thus of the whole System. 

    // This says whether the Instance-stage "disabled" flag for this force 
    // element should be initially on or off. Most force elements are enabled 
    // by default.
    bool                   defaultDisabled;

        // TOPOLOGY "CACHE"
    // Nothing in the base Impl class.
};



//------------------------------------------------------------------------------
//                    TWO POINT LINEAR SPRING IMPL
//------------------------------------------------------------------------------
class Force::TwoPointLinearSpringImpl : public ForceImpl {
public:
    TwoPointLinearSpringImpl(const MobilizedBody& body1, const Vec3& station1, 
                             const MobilizedBody& body2, const Vec3& station2, 
                             Real k, Real x0);

    TwoPointLinearSpringImpl* clone() const override {
        return new TwoPointLinearSpringImpl(*this);
    }
    bool dependsOnlyOnPositions() const override {
        return true;
    }
    bool shouldBeParallelized() const override {
        return true;
    }
    void calcForce(const State&         state, 
                   Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>&       particleForces, 
                   Vector&              mobilityForces) const override;
    Real calcPotentialEnergy(const State& state) const override;

    void calcDecorativeGeometryAndAppend(const State& s, Stage stage, 
                                         Array_<DecorativeGeometry>& geom) 
                                         const override;

private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real k, x0;
};



//------------------------------------------------------------------------------
//                    TWO POINT LINEAR DAMPER IMPL
//------------------------------------------------------------------------------
class Force::TwoPointLinearDamperImpl : public ForceImpl {
public:
    TwoPointLinearDamperImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real damping);
    TwoPointLinearDamperImpl* clone() const {
        return new TwoPointLinearDamperImpl(*this);
    }
    bool shouldBeParallelized() const override {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real damping;
};



//------------------------------------------------------------------------------
//                    TWO POINT CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
class Force::TwoPointConstantForceImpl : public ForceImpl {
public:
    TwoPointConstantForceImpl(const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real force);
    TwoPointConstantForceImpl* clone() const {
        return new TwoPointConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    bool shouldBeParallelized() const override {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body1, body2;
    Vec3 station1, station2;
    Real force;
};



//------------------------------------------------------------------------------
//                    MOBILITY LINEAR SPRING IMPL
//------------------------------------------------------------------------------
class Force::MobilityLinearSpringImpl : public ForceImpl {
friend class MobilityLinearSpring;

    MobilityLinearSpringImpl(const MobilizedBody& mobod, 
                             MobilizerQIndex      whichQ, 
                             Real                 defaultStiffness, 
                             Real                 defaultQZero);

    const std::pair<Real,Real>& getParams(const State& state) const {
        return Value< std::pair<Real,Real> >::downcast
           (getForceSubsystem().getDiscreteVariable(state, m_paramsIx));
    }
    std::pair<Real,Real>& updParams(State& state) const {
        return Value< std::pair<Real,Real> >::updDowncast
           (getForceSubsystem().updDiscreteVariable(state, m_paramsIx));
    }

    MobilityLinearSpringImpl* clone() const override
    {   return new MobilityLinearSpringImpl(*this); }
    bool dependsOnlyOnPositions() const override {return true;}
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override;
    Real calcPotentialEnergy(const State& state) const override;

    // Allocate the discrete state variable for the parameters. 
    void realizeTopology(State& s) const override {
        m_paramsIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Dynamics, 
                 new Value< std::pair<Real,Real> >
                        (std::make_pair(m_defaultStiffness, m_defaultQZero)));
    }

//------------------------------------------------------------------------------
    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_mobodIx;
    const MobilizerQIndex           m_whichQ;
    Real                            m_defaultStiffness;
    Real                            m_defaultQZero;

    // TOPOLOGY CACHE (Set only once in realizeTopology(); const thereafter.)
    // The value is a pair<Real,Real> containing k,q0
    mutable DiscreteVariableIndex   m_paramsIx;
};



//------------------------------------------------------------------------------
//                       MOBILITY LINEAR DAMPER IMPL
//------------------------------------------------------------------------------
class Force::MobilityLinearDamperImpl : public ForceImpl {
friend class MobilityLinearDamper;

    MobilityLinearDamperImpl(const MobilizedBody&   mobod, 
                             MobilizerUIndex        whichU, 
                             Real                   defaultDamping);

    const Real& getDamping(const State& state) const {
        return Value<Real>::downcast
           (getForceSubsystem().getDiscreteVariable(state, m_dampingIx));
    }
    Real& updDamping(State& state) const {
        return Value<Real>::updDowncast
           (getForceSubsystem().updDiscreteVariable(state, m_dampingIx));
    }

    MobilityLinearDamperImpl* clone() const override 
    {   return new MobilityLinearDamperImpl(*this); }

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override;
    Real calcPotentialEnergy(const State& state) const override;

    // Allocate the discrete state variable for the parameters. 
    void realizeTopology(State& s) const override {
        m_dampingIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Dynamics, 
                 new Value<Real>(m_defaultDamping));
    }

//------------------------------------------------------------------------------
    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_mobodIx;
    MobilizerUIndex                 m_whichU;
    Real                            m_defaultDamping;

    // TOPOLOGY CACHE (Set only once in realizeTopology(); const thereafter.)
    // The value is a Real containing c.
    mutable DiscreteVariableIndex   m_dampingIx;

};



//------------------------------------------------------------------------------
//                       MOBILITY CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
// Hidden implementation class for a MobilityConstantForce force element.
class Force::MobilityConstantForceImpl : public ForceImpl {
friend class MobilityConstantForce;

    MobilityConstantForceImpl(const MobilizedBody&  mobod, 
                              MobilizerUIndex       whichU, 
                              Real                  defaultForce);

    Real getForce(const State& state) const {
        return Value<Real>::downcast
           (getForceSubsystem().getDiscreteVariable(state, m_forceIx));
    }
    Real& updForce(State& state) const {
        return Value<Real>::updDowncast
           (getForceSubsystem().updDiscreteVariable(state, m_forceIx));
    }

    // Implementation of virtual methods from ForceImpl:

    MobilityConstantForceImpl* clone() const override 
    {   return new MobilityConstantForceImpl(*this); }

    // Has to wait for Dynamics stage because that's all that gets invalidated
    // if the constant force is changed.
    bool dependsOnlyOnPositions() const override {return false;}

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override;

    Real calcPotentialEnergy(const State& state) const override {return 0;}

    // Allocate the discrete state variable for the force. 
    void realizeTopology(State& s) const override {
        m_forceIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Dynamics, 
                                      new Value<Real>(m_defaultForce));
    }

//------------------------------------------------------------------------------
    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_mobodIx;
    const MobilizerUIndex           m_whichU;
    Real                            m_defaultForce;

    // TOPOLOGY CACHE (Set only once in realizeTopology(); const thereafter.)
    mutable DiscreteVariableIndex   m_forceIx;
};

//------------------------------------------------------------------------------
//                       MOBILITY LINEAR STOP IMPL
//------------------------------------------------------------------------------
// Hidden implementation class for a MobilityLinearStop force element.
class Force::MobilityLinearStopImpl : public ForceImpl {
friend class MobilityLinearStop;

    // Type of the discrete state variable that holds values for this
    // stop's changeable parameters in a State. Since these affect only forces
    // the variable invalidates Dynamics stage and later when it changes.
    struct Parameters {
        Parameters(Real defStiffness, Real defDissipation,
                     Real defQLow, Real defQHigh)
        :   k(defStiffness), d(defDissipation), 
            qLow(defQLow), qHigh(defQHigh) {}

        Real    k, d, qLow, qHigh;
    };

    MobilityLinearStopImpl(const MobilizedBody&      mobod, 
                           MobilizerQIndex           whichQ, 
                           Real                      defaultStiffness,
                           Real                      defaultDissipation,
                           Real                      defaultQLow,
                           Real                      defaultQHigh);

    // Implementation of virtual methods from ForceImpl:
    MobilityLinearStopImpl* clone() const override 
    {   return new MobilityLinearStopImpl(*this); }
    bool dependsOnlyOnPositions() const override {return false;}

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override; 

    // We're not bothering to cache P.E. -- just recalculate it when asked.
    Real calcPotentialEnergy(const State& state) const override; 

    // Allocate the state variables and cache entry. 
    void realizeTopology(State& s) const override {
        // Allocate the discrete variable for dynamics parameters.
        const Parameters dv(m_defStiffness, m_defDissipation, 
                            m_defQLow, m_defQHigh);
        m_parametersIx = getForceSubsystem()
            .allocateDiscreteVariable(s, Stage::Dynamics, 
                                      new Value<Parameters>(dv));
    }

    const Parameters& getParameters(const State& s) const
    {   return Value<Parameters>::downcast
           (getForceSubsystem().getDiscreteVariable(s,m_parametersIx)); }
    Parameters& updParameters(State& s) const
    {   return Value<Parameters>::updDowncast
           (getForceSubsystem().updDiscreteVariable(s,m_parametersIx)); }

//------------------------------------------------------------------------------
    // TOPOLOGY STATE
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_mobodIx;
    const MobilizerQIndex           m_whichQ;

    Real m_defStiffness, m_defDissipation, m_defQLow, m_defQHigh;

    // TOPOLOGY CACHE (Set only once in realizeTopology(); const thereafter.)
    mutable DiscreteVariableIndex   m_parametersIx; // k, d, qLow, qHigh
};



//------------------------------------------------------------------------------
//                    MOBILITY DISCRETE FORCE IMPL
//------------------------------------------------------------------------------
// This Force element allocates a scalar discrete variable and applies its
// value as a generalized force at a particular mobility. The value can be
// set externally in an event handler or between steps.
class Force::MobilityDiscreteForceImpl : public ForceImpl {
public:
    MobilityDiscreteForceImpl(const MobilizedBody&  mobod, 
                              MobilizerUIndex       whichU, 
                              Real                  defaultForce);

    // Change the force value to be applied. The force will remain at this
    // value until changed again.
    void setMobilityForce(State& state, Real f) const;

    // Get the value of the generalized force to be applied.
    Real getMobilityForce(const State& state) const;

    // Override five virtuals from base class:

    // This is called at Simbody's realize(Dynamics) stage.
    void calcForce( const State&         state, 
                    Vector_<SpatialVec>& /*bodyForces*/, 
                    Vector_<Vec3>&       /*particleForces*/, 
                    Vector&              mobilityForces) const override;

    // This force element does not store potential energy.
    Real calcPotentialEnergy(const State& state) const override {return 0;}

    // Allocate the needed state variable and record its index.
    void realizeTopology(State& state) const override;

    MobilityDiscreteForceImpl* clone() const override 
    {   return new MobilityDiscreteForceImpl(*this); }

    // Force this to wait for Dynamics stage before calculating, because that's
    // all that gets invalidated when a new forces is applied.
    bool dependsOnlyOnPositions() const override {return false;}

private:
friend class Force::MobilityDiscreteForce;

    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_mobodIx;
    const MobilizerUIndex           m_whichU;
    Real                            m_defaultVal;

    mutable DiscreteVariableIndex   m_forceIx;
};



//------------------------------------------------------------------------------
//                         DISCRETE FORCES IMPL
//------------------------------------------------------------------------------
// This Force element allocates a Vector discrete variables and applies their
// values as generalized and body spatial forces. The values can be
// set externally in an event handler or between steps.
class Force::DiscreteForcesImpl : public ForceImpl {
public:
    DiscreteForcesImpl(const SimbodyMatterSubsystem& matter);

    const Vector& getAllMobilityForces(const State& state) const;
    Vector& updAllMobilityForces(State& state) const;

    const Vector_<SpatialVec>& getAllBodyForces(const State& state) const;
    Vector_<SpatialVec>& updAllBodyForces(State& state) const;

    // Override five virtuals from base class:

    // This is called at Simbody's realize(Dynamics) stage.
    void calcForce( const State&         state, 
                    Vector_<SpatialVec>& bodyForces, 
                    Vector_<Vec3>&       particleForces, 
                    Vector&              mobilityForces) const override;

    // This force element does not store potential energy.
    Real calcPotentialEnergy(const State& state) const override {return 0;}

    // Allocate the needed state variable and record its index.
    void realizeTopology(State& state) const override;

    DiscreteForcesImpl* clone() const override 
    {   return new DiscreteForcesImpl(*this); }

    // Force this to wait for Dynamics stage before calculating, because that's
    // all that gets invalidated when a new forces is applied.
    bool dependsOnlyOnPositions() const override {return false;}

private:
friend class Force::DiscreteForces;

    const SimbodyMatterSubsystem&   m_matter;

    mutable DiscreteVariableIndex   m_mobForcesIx;  // Vector(n)
    mutable DiscreteVariableIndex   m_bodyForcesIx; // Vector_<SpatialVec>(nb)
};



//------------------------------------------------------------------------------
//                           CONSTANT FORCE IMPL
//------------------------------------------------------------------------------
class Force::ConstantForceImpl : public ForceImpl {
public:
    ConstantForceImpl(const MobilizedBody& body, const Vec3& station, const Vec3& force);
    ConstantForceImpl* clone() const {
        return new ConstantForceImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 station, force;
};



//------------------------------------------------------------------------------
//                          CONSTANT TORQUE IMPL
//------------------------------------------------------------------------------
class Force::ConstantTorqueImpl : public ForceImpl {
public:
    ConstantTorqueImpl(const MobilizedBody& body, const Vec3& torque);
    ConstantTorqueImpl* clone() const {
        return new ConstantTorqueImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    const MobilizedBodyIndex body;
    Vec3 torque;
};



//------------------------------------------------------------------------------
//                          GLOBAL DAMPER IMPL
//------------------------------------------------------------------------------
class Force::GlobalDamperImpl : public ForceImpl {
public:
    GlobalDamperImpl(const SimbodyMatterSubsystem& matter, Real damping);
    GlobalDamperImpl* clone() const {
        return new GlobalDamperImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
private:
    const SimbodyMatterSubsystem& matter;
    Real damping;
};



//------------------------------------------------------------------------------
//                          UNIFORM GRAVITY IMPL
//------------------------------------------------------------------------------
class Force::UniformGravityImpl : public ForceImpl {
public:
    UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight);
    UniformGravityImpl* clone() const {
        return new UniformGravityImpl(*this);
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    Vec3 getGravity() const {
        return g;
    }
    void setGravity(const Vec3& gravity) {
        g = gravity;
        invalidateTopologyCache();
    }
    Real getZeroHeight() const {
        return zeroHeight;
    }
    void setZeroHeight(Real height) {
        zeroHeight = height;
        invalidateTopologyCache();
    }
private:
    const SimbodyMatterSubsystem& matter;
    Vec3 g;
    Real zeroHeight;
};



//------------------------------------------------------------------------------
//                              CUSTOM IMPL
//------------------------------------------------------------------------------
class Force::CustomImpl : public ForceImpl {
public:
    CustomImpl(Force::Custom::Implementation* implementation);
    CustomImpl* clone() const {
        return new CustomImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return implementation->dependsOnlyOnPositions();
    }
    bool shouldBeParallelized() const {
        return implementation->shouldBeParallelized();
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) 
                   const override;
    Real calcPotentialEnergy(const State& state) const override;
    ~CustomImpl() {
        delete implementation;
    }
    const Force::Custom::Implementation& getImplementation() const {
        return *implementation;
    }
    Force::Custom::Implementation& updImplementation() {
        return *implementation;
    }
protected:
    void realizeTopology(State& state) const override {
        implementation->realizeTopology(state);
    }
    void realizeModel(State& state) const override {
        implementation->realizeModel(state);
    }
    void realizeInstance(const State& state) const override {
        implementation->realizeInstance(state);
    }
    void realizeTime(const State& state) const override {
        implementation->realizeTime(state);
    }
    void realizePosition(const State& state) const override {
        implementation->realizePosition(state);
    }
    void realizeVelocity(const State& state) const override {
        implementation->realizeVelocity(state);
    }
    void realizeDynamics(const State& state) const override {
        implementation->realizeDynamics(state);
    }
    void realizeAcceleration(const State& state) const override {
        implementation->realizeAcceleration(state);
    }
    void realizeReport(const State& state) const override {
        implementation->realizeReport(state);
    }
    void calcDecorativeGeometryAndAppend(const State& state, Stage stage, 
        Array_<DecorativeGeometry>& geom) const override 
    {
        implementation->calcDecorativeGeometryAndAppend(state,stage,geom);
    }
private:
    Force::Custom::Implementation* implementation;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_IMPL_H_
