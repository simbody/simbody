/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include "MotionImpl.h"
#include "MobilizedBodyImpl.h"

namespace SimTK {

const MobilizedBody& Motion::getMobilizedBody() const 
{   return getImpl().getMobilizedBodyImpl().getMyHandle(); }

Motion::Level Motion::getLevel(const State& s) const
{   return getImpl().getLevel(s); }

Motion::Method Motion::getLevelMethod(const State& s) const
{   return getImpl().getLevelMethod(s); }

/*static*/ const char*
Motion::nameOfLevel(Level l) {
    switch (l) {
      case Acceleration: return "Acceleration";
      case Velocity:     return "Velocity";
      case Position:     return "Position";
      default: 
        return "*** UNRECOGNIZED Motion::Level ***";
    }
}

/*static*/ const char*
Motion::nameOfMethod(Method m) {
    switch (m) {
      case Zero:       return "Zero";
      case Discrete:   return "Discrete";
      case Prescribed: return "Prescribed";
      case Free:       return "Free";
      case Fast:       return "Fast";
      default: 
        return "*** UNRECOGNIZED Motion::Method ***";
    }
}


const SimbodyMatterSubsystem&  
MotionImpl::getMatterSubsystem() const {
    return getMobilizedBodyImpl().getMySimbodyMatterSubsystem(); 
}

const AbstractValue&
MotionImpl::getDiscreteVariable(const State& s, DiscreteVariableIndex vx) const {
    return getMatterSubsystem().getDiscreteVariable(s, vx);
}

AbstractValue&
MotionImpl::updDiscreteVariable(State& s, DiscreteVariableIndex vx) const {
    return getMatterSubsystem().updDiscreteVariable(s, vx);
}

DiscreteVariableIndex
MotionImpl::allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const {
    return getMatterSubsystem().allocateDiscreteVariable(s, g, v);
}


void MotionImpl::invalidateTopologyCache() const {
    if (hasMobilizedBody()) 
        getMatterSubsystem().invalidateSubsystemTopologyCache();
}

// These are the default implementations for these virtual methods. They
// will throw exceptions if called.

void MotionImpl::
calcPrescribedPositionVirtual(const State& s, int nq, Real* q) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedPositionVirtual()",
        "A built-in Motion class did not supply a calcPrescribedPositionVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void MotionImpl::
calcPrescribedPositionDotVirtual(const State& s, int nq, Real* qdot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedPositionDotVirtual()",
        "A built-in Motion class did not supply a calcPrescribedPositionDotVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void MotionImpl::
calcPrescribedPositionDotDotVirtual(const State& s, int nq, Real* qdotdot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedPositionDotDotVirtual()",
        "A built-in Motion class did not supply a calcPrescribedPositionDotDotVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void MotionImpl::
calcPrescribedVelocityVirtual(const State& s, int nu, Real* u) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedVelocityVirtual()",
        "A built-in Motion class did not supply a calcPrescribedVelocityVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void MotionImpl::
calcPrescribedVelocityDotVirtual(const State& s, int nu, Real* udot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedVelocityDotVirtual()",
        "A built-in Motion class did not supply a calcPrescribedVelocityDotVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void MotionImpl::
calcPrescribedAccelerationVirtual(const State& s, int nu, Real* udot) const  {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "MotionImpl::calcPrescribedAccelerationVirtual()",
        "A built-in Motion class did not supply a calcPrescribedAccelerationVirtual() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}


//-------------------------------- Sinusoid ------------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS
   (Motion::Sinusoid, Motion::SinusoidImpl, Motion);

Motion::Sinusoid::Sinusoid(MobilizedBody& mobod, Motion::Level level,
                           Real amplitude, Real rate, Real phase)
:   Motion(new SinusoidImpl(level, amplitude, rate, phase)) {
    mobod.adoptMotion(*this);
}

//--------------------------------- Steady -------------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Motion::Steady, Motion::SteadyImpl, Motion);


Motion::Steady::Steady(MobilizedBody& mobod, Real u) 
  : Motion(new SteadyImpl(Vec<6>(u))) {
    mobod.adoptMotion(*this);
}

// These are specializations for each of the six possible instantiations
// of the templatized constructor. Anything else will be unresolved.

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<1>& u)
  : Motion(new SteadyImpl(Vec<6>(u[0],0,0,0,0,0))) {
    mobod.adoptMotion(*this);
}

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<2>& u)
  : Motion(new SteadyImpl(Vec<6>(u[0],u[1],0,0,0,0))) {
    mobod.adoptMotion(*this);
}

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<3>& u)
  : Motion(new SteadyImpl(Vec<6>(u[0],u[1],u[2],0,0,0))) {
    mobod.adoptMotion(*this);
}

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<4>& u)
  : Motion(new SteadyImpl(Vec<6>(u[0],u[1],u[2],u[3],0,0))) {
    mobod.adoptMotion(*this);
}

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<5>& u)
  : Motion(new SteadyImpl(Vec<6>(u[0],u[1],u[2],u[3],u[4],0))) {
    mobod.adoptMotion(*this);
}

template <> SimTK_SIMBODY_EXPORT
Motion::Steady::Steady(MobilizedBody& mobod, const Vec<6>& u)
  : Motion(new SteadyImpl(u)) {
    mobod.adoptMotion(*this);
}

Motion::Steady&
Motion::Steady::setDefaultRate(Real u)
{   updImpl().setDefaultRates(Vec6(u)); return *this; }
Motion::Steady& 
Motion::Steady::setOneDefaultRate(MobilizerUIndex ux, Real u)
{   updImpl().setOneDefaultRate(ux,u); return *this; }

void Motion::Steady::setRate(State& s, Real u) const
{   getImpl().setRates(s, Vec6(u)); }
void Motion::Steady::setOneRate(State& s, MobilizerUIndex ux, Real u) const
{   getImpl().setOneRate(s,ux,u); }

template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<1>& u)
{   updImpl().setDefaultRates(Vec6(u[0],0,0,0,0,0)); return *this; }
template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<2>& u)
{   updImpl().setDefaultRates(Vec6(u[0],u[1],0,0,0,0)); return *this; }
template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<3>& u)
{   updImpl().setDefaultRates(Vec6(u[0],u[1],u[2],0,0,0)); return *this; }
template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<4>& u)
{   updImpl().setDefaultRates(Vec6(u[0],u[1],u[2],u[3],0,0)); return *this; }
template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<5>& u)
{   updImpl().setDefaultRates(Vec6(u[0],u[1],u[2],u[3],u[4],0)); return *this; }
template <> SimTK_SIMBODY_EXPORT Motion::Steady&
Motion::Steady::setDefaultRates(const Vec<6>& u)
{   updImpl().setDefaultRates(u); return *this; }


//---------------------------------- Custom ------------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Motion::Custom, Motion::CustomImpl, Motion);

Motion::Custom::Custom(MobilizedBody& mobod, Implementation* implementation)
  : Motion(new CustomImpl(implementation)) {
    mobod.adoptMotion(*this);
}

const Motion::Custom::Implementation& Motion::Custom::getImplementation() const {
    return getImpl().getImplementation();
}

Motion::Custom::Implementation& Motion::Custom::updImplementation() {
    return updImpl().updImplementation();    
}


Motion::CustomImpl::CustomImpl(Motion::Custom::Implementation* implementation)
  : implementation(implementation) {}


// These are the default implementations for these virtual methods. They
// will throw exceptions if called.

void Motion::Custom::Implementation::
calcPrescribedPosition(const State& s, int nq, Real* q) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedPosition()",
        "Concrete Implementation did not supply a calcPrescribedPosition() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void Motion::Custom::Implementation::
calcPrescribedPositionDot(const State& s, int nq, Real* qdot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedPositionDot()",
        "Concrete Implementation did not supply a calcPrescribedPositionDot() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void Motion::Custom::Implementation::
calcPrescribedPositionDotDot(const State& s, int nq, Real* qdotdot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedPositionDotDot()",
        "Concrete Implementation did not supply a calcPrescribedPositionDotDot() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void Motion::Custom::Implementation::
calcPrescribedVelocity(const State& s, int nu, Real* u) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedVelocity()",
        "Concrete Implementation did not supply a calcPrescribedVelocity() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void Motion::Custom::Implementation::
calcPrescribedVelocityDot(const State& s, int nu, Real* udot) const {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedVelocityDot()",
        "Concrete Implementation did not supply a calcPrescribedVelocityDot() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

void Motion::Custom::Implementation::
calcPrescribedAcceleration(const State& s, int nu, Real* udot) const  {
    SimTK_ERRCHK_ALWAYS(!"unimplemented",
        "Motion::Custom::Implementation::calcPrescribedAcceleration()",
        "Concrete Implementation did not supply a calcPrescribedAcceleration() method "
        "but was apparently expected to do so.");
    /*NOTREACHED*/
}

} // namespace SimTK

