#ifndef SimTK_SIMBODY_MOTION_IMPL_H_
#define SimTK_SIMBODY_MOTION_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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
#include "simbody/internal/Motion.h"

#include "SimbodyTreeState.h"

namespace SimTK {

class MobilizedBodyImpl;


//==============================================================================
//                               MOTION IMPL
//==============================================================================
// This is the hidden implementation class for Motion objects. This is the
// abstract base class to which every Motion handle points.
class MotionImpl : public PIMPLImplementation<Motion, MotionImpl> {
public:
    MotionImpl()
    :   m_isDisabledByDefault(false) {}

    // Default copy constructor, copy assignment, and destructor; note that the
    // pointer to the mobilized body is not copied or deleted.

    bool hasMobilizedBody() const {return m_mobodImpl != nullptr;}
    const MobilizedBodyImpl& getMobilizedBodyImpl() const
    {   assert(m_mobodImpl); return *m_mobodImpl; }
    MobilizedBodyIndex getMobilizedBodyIndex() const;

    void setMobilizedBodyImpl(MobilizedBodyImpl* mbi)
    {   assert(!m_mobodImpl); m_mobodImpl = mbi; }
    void invalidateTopologyCache() const;

    const SimbodyMatterSubsystem& getMatterSubsystem() const;

    const AbstractValue&
    getDiscreteVariable(const State& s, DiscreteVariableIndex vx) const;

    AbstractValue&
    updDiscreteVariable(State& s, DiscreteVariableIndex vx) const;

    DiscreteVariableIndex
    allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const;

    template <class T> const T&
    getVar(const State& s, DiscreteVariableIndex vx) const {
        return Value<T>::downcast(getDiscreteVariable(s, vx));
    }

    template <class T> T&
    updVar(State& s, DiscreteVariableIndex vx) const {
        return Value<T>::updDowncast(updDiscreteVariable(s, vx));
    }

    template <class T> DiscreteVariableIndex
    allocVar(State& state, const T& initVal,
             const Stage& stage=Stage::Instance) const
    {
        return allocateDiscreteVariable(state, stage, new Value<T>(initVal));
    }


    virtual ~MotionImpl() {}
    virtual MotionImpl* clone() const = 0;

    // This reports whether this Motion is holonomic (Level::Position),
    // nonholonomic (Level::Velocity), or acceleration (Level::Acceleration).
    Motion::Level getLevel(const State& s) const {
        if (isDisabled(s)) return Motion::NoLevel;
        return getLevelVirtual(s); // ask concrete class
    }

    Motion::Method getLevelMethod(const State& s) const {
        if (isDisabled(s)) return Motion::NoMethod;
        return getLevelMethodVirtual(s);
    }


    void disable(State& s) const;
    void enable(State& s) const;
    bool isDisabled(const State&) const;

    void setDisabledByDefault(bool shouldBeDisabled)
    {   invalidateTopologyCache(); m_isDisabledByDefault=shouldBeDisabled; }
    bool isDisabledByDefault() const {return m_isDisabledByDefault;}

    // These operators calculate prescribed positions, velocities, or
    // accelerations given a State realized to the previous Stage.
    void calcPrescribedPosition      (const State& s, int nq, Real* q)      const
    {   calcPrescribedPositionVirtual(s,nq,q); }
    void calcPrescribedPositionDot   (const State& s, int nq, Real* qdot)   const
    {   calcPrescribedPositionDotVirtual(s,nq,qdot); }
    void calcPrescribedPositionDotDot(const State& s, int nq, Real* qdotdot)const
    {   calcPrescribedPositionDotDotVirtual(s,nq,qdotdot); }
    void calcPrescribedVelocity      (const State& s, int nu, Real* u)      const
    {   calcPrescribedVelocityVirtual(s,nu,u); }
    void calcPrescribedVelocityDot   (const State& s, int nu, Real* udot)   const
    {   calcPrescribedVelocityDotVirtual(s,nu,udot); }
    void calcPrescribedAcceleration  (const State& s, int nu, Real* udot)   const
    {   calcPrescribedAccelerationVirtual(s,nu,udot); }

    void realizeTopology(State& state)              const
    {   realizeTopologyVirtual(state); }
    void realizeModel(State& state)                 const
    {   realizeModelVirtual(state); }
    void realizeInstance(const State& state)        const
    {   realizeInstanceVirtual(state); }
    void realizeTime(const State& state)            const
    {   realizeTimeVirtual(state); }
    void realizePosition(const State& state)        const
    {   realizePositionVirtual(state); }
    void realizeVelocity(const State& state)        const
    {   realizeVelocityVirtual(state); }
    void realizeDynamics(const State& state)        const
    {   realizeDynamicsVirtual(state); }
    void realizeAcceleration(const State& state)    const
    {   realizeAccelerationVirtual(state); }
    void realizeReport(const State& state)          const
    {   realizeReportVirtual(state); }

    virtual Motion::Level getLevelVirtual(const State&) const = 0;
    virtual Motion::Method getLevelMethodVirtual(const State&) const
    {   return Motion::Prescribed; }

    virtual void calcPrescribedPositionVirtual
                   (const State&, int nq, Real* q)          const;
    virtual void calcPrescribedPositionDotVirtual
                   (const State&, int nq, Real* qdot)       const;
    virtual void calcPrescribedPositionDotDotVirtual
                   (const State&, int nq, Real* qdotdot)    const;
    virtual void calcPrescribedVelocityVirtual
                   (const State&, int nu, Real* u)          const;
    virtual void calcPrescribedVelocityDotVirtual
                   (const State&, int nu, Real* udot)       const;
    virtual void calcPrescribedAccelerationVirtual
                   (const State&, int nu, Real* udot)       const;

    virtual void realizeTopologyVirtual    (State&)         const {}
    virtual void realizeModelVirtual       (State&)         const {}
    virtual void realizeInstanceVirtual    (const State&)   const {}
    virtual void realizeTimeVirtual        (const State&)   const {}
    virtual void realizePositionVirtual    (const State&)   const {}
    virtual void realizeVelocityVirtual    (const State&)   const {}
    virtual void realizeDynamicsVirtual    (const State&)   const {}
    virtual void realizeAccelerationVirtual(const State&)   const {}
    virtual void realizeReportVirtual      (const State&)   const {}
private:
    ReferencePtr<MobilizedBodyImpl>     m_mobodImpl;
    bool                                m_isDisabledByDefault;
};


//------------------------------------------------------------------------------
//                               SINUSOID IMPL
//------------------------------------------------------------------------------
class Motion::SinusoidImpl : public MotionImpl {
public:
    // no default constructor
    SinusoidImpl(Motion::Level level,
                 Real amplitude, Real rate, Real phase)
    :   level(level), defAmplitude(amplitude), defRate(rate),
        defPhase(phase)
    {
    }

    SinusoidImpl* clone() const override {
        SinusoidImpl* copy = new SinusoidImpl(*this);
        return copy;
    }

    Motion::Level  getLevelVirtual (const State&) const override
    {   return level; }
    Motion::Method getLevelMethodVirtual(const State&) const override
    {   return Motion::Prescribed; }

    // Allocate variables if needed.
    void realizeTopologyVirtual(State& state) const override {
        // None yet.
    }


    void calcPrescribedPositionVirtual
       (const State& state, int nq, Real* q) const override {
        assert(level==Motion::Position); assert(nq==0 || q);
        const Real t = state.getTime();
        const Real out = defAmplitude*std::sin(defRate*t + defPhase);
        for (int i=0; i<nq; ++i)
            q[i] = out;
    }

    void calcPrescribedPositionDotVirtual
       (const State& state, int nq, Real* qdot) const override {
        assert(level==Motion::Position); assert(nq==0 || qdot);
        const Real t = state.getTime();
        const Real outd = defAmplitude*defRate*std::cos(defRate*t + defPhase);
        for (int i=0; i<nq; ++i)
            qdot[i] = outd;
    }

    void calcPrescribedPositionDotDotVirtual
       (const State& state, int nq, Real* qdotdot) const override {
        assert(level==Motion::Position); assert(nq==0 || qdotdot);
        const Real t = state.getTime();
        const Real outdd =
            -defAmplitude*defRate*defRate*std::sin(defRate*t + defPhase);
        for (int i=0; i<nq; ++i)
            qdotdot[i] = outdd;
    }

    void calcPrescribedVelocityVirtual
       (const State& state, int nu, Real* u) const override {
        assert(level==Motion::Velocity);
        assert(nu==0 || u);
        const Real t = state.getTime();
        const Real out = defAmplitude*std::sin(defRate*t + defPhase);
        for (int i=0; i<nu; ++i)
            u[i] = out;
    }

    void calcPrescribedVelocityDotVirtual
       (const State& state, int nu, Real* udot) const override {
        assert(level==Motion::Velocity);
        assert(nu==0 || udot);
        const Real t = state.getTime();
        const Real outd = defAmplitude*defRate*std::cos(defRate*t + defPhase);
        for (int i=0; i<nu; ++i)
            udot[i] = outd;
    }

    void calcPrescribedAccelerationVirtual
       (const State& state, int nu, Real* udot) const override {
        assert(level==Motion::Acceleration); assert(nu==0 || udot);
        const Real t = state.getTime();
        const Real out = defAmplitude*std::sin(defRate*t + defPhase);
        for (int i=0; i<nu; ++i)
            udot[i] = out;
    }
private:
        // TOPOLOGY "STATE"
    Motion::Level       level;
    Real                defAmplitude, defRate, defPhase;

        // TOPOLOGY "CACHE"
    // None yet.
};


//------------------------------------------------------------------------------
//                               STEADY IMPL
//------------------------------------------------------------------------------
class Motion::SteadyImpl : public MotionImpl {
public:
    // no default constructor
    explicit SteadyImpl(const Vec6& u) : defaultU(u) {}

    SteadyImpl* clone() const override {
        SteadyImpl* copy = new SteadyImpl(*this);
        copy->currentU.invalidate(); // no sharing state variables
        return copy;
    }

    void setDefaultRates(const Vec6& u) {
        invalidateTopologyCache();
        defaultU = u;
    }
    void setOneDefaultRate(MobilizerUIndex ux, Real u) {
        invalidateTopologyCache();
        defaultU[ux] = u;
    }
    Real getOneDefaultRate(MobilizerUIndex ux) const {
        return defaultU[ux];
    }

    const Vec6& getDefaultRates() const {return defaultU;}

    void setRates(State& s, const Vec6& u) const {
        updVar<Vec6>(s, currentU) = u;
    }
    void setOneRate(State& s, MobilizerUIndex ux, Real u) const {
        updVar<Vec6>(s, currentU)[ux] = u;
    }

    Real getOneRate(const State& s, MobilizerUIndex ux) const {
        return getVar<Vec6>(s, currentU)[ux];
    }

    Motion::Level  getLevelVirtual (const State&) const override
    {   return Motion::Velocity; }
    Motion::Method getLevelMethodVirtual(const State&) const override
    {   return Motion::Prescribed; }

    // Allocate a discrete variable to hold the constant rates.
    void realizeTopologyVirtual(State& state) const override {
        // This is in the Topology-stage "cache" so we can write to it,
        // but only here.
        const_cast<DiscreteVariableIndex&>(currentU) =
            allocVar(state, defaultU);
    }

    void calcPrescribedVelocityVirtual
       (const State& state, int nu, Real* u) const override
    {
        assert(0 <= nu && nu <= 6);
        assert(nu==0 || u);
        const Vec6& uval = getVar<Vec6>(state, currentU);
        for (int i=0; i<nu; ++i)
            u[i] = uval[i];
    }

    void calcPrescribedVelocityDotVirtual
       (const State& state, int nu, Real* udot) const override
    {
        assert(0 <= nu && nu <= 6);
        assert(nu==0 || udot);
        for (int i=0; i<nu; ++i)
            udot[i] = 0;
    }

private:
        // TOPOLOGY "STATE"
    Vec6                  defaultU;

        // TOPOLOGY "CACHE"
    DiscreteVariableIndex currentU;
};


//------------------------------------------------------------------------------
//                               CUSTOM IMPL
//------------------------------------------------------------------------------
class Motion::CustomImpl : public MotionImpl {
public:
    // Take over ownership of the supplied heap-allocated object.
    explicit CustomImpl(Motion::Custom::Implementation* implementation);

    CustomImpl(const CustomImpl& src) : implementation(0) {
        if (src.implementation)
            implementation = src.implementation->clone();
    }

    CustomImpl* clone() const override { return new CustomImpl(*this); }

    ~CustomImpl() {
        delete implementation;
    }

    Motion::Level getLevelVirtual(const State& s) const override {
        return getImplementation().getLevel(s);
    }
    Motion::Method getLevelMethodVirtual(const State& s) const override {
        return getImplementation().getLevelMethod(s);
    }

    const Motion::Custom::Implementation& getImplementation() const {
        assert(implementation); return *implementation;
    }
    Motion::Custom::Implementation& updImplementation() {
        assert(implementation); return *implementation;
    }

    void calcPrescribedPositionVirtual
       (const State& s, int nq, Real* q) const override
    {   getImplementation().calcPrescribedPosition(s,nq,q); }
    void calcPrescribedPositionDotVirtual
       (const State& s, int nq, Real* qdot) const override
    {   getImplementation().calcPrescribedPositionDot(s,nq,qdot); }
    void calcPrescribedPositionDotDotVirtual
       (const State& s, int nq, Real* qdotdot) const override
    {   getImplementation().calcPrescribedPositionDotDot(s,nq,qdotdot); }

    void calcPrescribedVelocityVirtual
       (const State& s, int nu, Real* u) const override
    {   getImplementation().calcPrescribedVelocity(s,nu,u); }
    void calcPrescribedVelocityDotVirtual
       (const State& s, int nu, Real* udot) const override
    {   getImplementation().calcPrescribedVelocityDot(s,nu,udot); }

    void calcPrescribedAccelerationVirtual
       (const State& s, int nu, Real* udot) const override
    {   getImplementation().calcPrescribedAcceleration(s,nu,udot); }

    void realizeTopologyVirtual(State& state) const override {
        getImplementation().realizeTopology(state);
    }
    void realizeModelVirtual(State& state) const override {
        getImplementation().realizeModel(state);
    }
    void realizeInstanceVirtual(const State& state) const override {
        getImplementation().realizeInstance(state);
    }
    void realizeTimeVirtual(const State& state) const override {
        getImplementation().realizeTime(state);
    }
    void realizePositionVirtual(const State& state) const override {
        getImplementation().realizePosition(state);
    }
    void realizeVelocityVirtual(const State& state) const override {
        getImplementation().realizeVelocity(state);
    }
    void realizeDynamicsVirtual(const State& state) const override {
        getImplementation().realizeDynamics(state);
    }
    void realizeAccelerationVirtual(const State& state) const override {
        getImplementation().realizeAcceleration(state);
    }
    void realizeReportVirtual(const State& state) const override {
        getImplementation().realizeReport(state);
    }
private:
    Motion::Custom::Implementation* implementation;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOTION_IMPL_H_
