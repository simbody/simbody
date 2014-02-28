#ifndef SimTK_SIMBODY_SEMI_EXPLICIT_EULER_TIME_STEPPER_H_
#define SimTK_SIMBODY_SEMI_EXPLICIT_EULER_TIME_STEPPER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Thomas Uchida                                    *
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

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ImpulseSolver.h"
#include "simbody/internal/PGSImpulseSolver.h"
#include "simbody/internal/PLUSImpulseSolver.h"

namespace SimTK {

/** TODO: A low-accuracy, high performance, velocity-level time stepper for
models containing unilateral rigid contact or other conditional constraints.

At each step, this solver determines which subset of the multibody system's
conditional constraints may affect behavior. It then formulates the resulting
equations and invokes an ImpulseSolver to find constraint-space impulses that
solve the constraint equations and satisfy all the constraint inequalities.
Generally there are multiple solutions possible, and different ImpulseSolver
objects use different criteria.
**/

class SimTK_SIMBODY_EXPORT SemiExplicitEulerTimeStepper {
public:
    enum RestitutionModel {Poisson=0, PoissonOnce=1, Newton=2};
    enum Solver {PLUS=0, PGS=1};
    enum PositionProjectionMethod {Bilateral=0,Unilateral=1};

    explicit SemiExplicitEulerTimeStepper(const MultibodySystem& mbs)
    :   m_mbs(mbs), 
        m_accuracy(1e-2), m_consTol(1e-3), m_restitutionModel(Poisson),
        m_solver(PLUS), m_projectionMethod(Bilateral),
        m_defaultMinCORVelocity(1),
        m_defaultCaptureVelocity(m_consTol/10),
        m_defaultTransitionVelocity(m_consTol/10),
        m_defaultSignificantForce(SignificantReal)
    {  }

    const MultibodySystem& getMultibodySystem() const {return m_mbs;}

    void setRestitutionModel(RestitutionModel restModel)
    {   m_restitutionModel = restModel; }
    RestitutionModel getRestitutionModel() const {return m_restitutionModel;}

    /** Set integration accuracy; requires variable length steps. **/
    void setAccuracy(Real accuracy) {m_accuracy=accuracy;}
    /** Set the tolerance to which constraints must be satisfied. **/
    void setConstraintTol(Real consTol) {m_consTol=consTol;}

    /** Set the impact capture velocity to be used by default when a contact
    does not provide its own. This is the impact velocity below which the
    coefficient of restitution is to be treated as zero. This avoids a Zeno's
    paradox of ever-tinier time-wasting impacts. A capture velocity should
    be significantly greater than the velocity constraint tolerance
    since speeds below that are indistinguishable from zero anyway. In 
    practice we will use as the capture velocity the \e larger of this value 
    and the velocity constraint tolerance. **/
    void setDefaultImpactCaptureVelocity(Real vCapture) {
        SimTK_ERRCHK1_ALWAYS(vCapture>=0,
        "SemiExplicitEulerTimeStepper::setDefaultImpactCaptureVelocity()",
        "The impact capture velocity must be nonnegative but was %g.",
        vCapture);
        m_defaultCaptureVelocity = vCapture;
    }

    /** Set the minimum coefficient of restitution (COR) velocity to be used by
    default when a unilateral contact does not provide its own. This is
    the velocity at which the COR reaches its specified value; below this
    velocity the COR increases until the capture velocity is reached, at which
    point the COR is set to zero. This velocity cannot be less than the 
    capture velocity, but we don't check here. Instead we'll check at run time
    and use the larger of the minimum COR velocity and the capture velocity. **/
    void setDefaultImpactMinCORVelocity(Real vMinCOR) {
        SimTK_ERRCHK1_ALWAYS(vMinCOR>=0,
        "SemiExplicitEulerTimeStepper::setDefaultImpactMinCORVelocity()",
        "The velocity at which the minimum coefficient of restitution "
        " is reached must be nonnegative but was %g.", vMinCOR);
        m_defaultMinCORVelocity = vMinCOR;
    }

    /** Set the friction sliding-to-rolling transition velocity to be used 
    by default when a frictional contact does not provide its own. This is the
    slip velocity below which we are permitted to consider that the contact 
    may be rolling, that is, in stiction. The transition velocity should
    be significantly greater than the velocity constraint tolerance
    since speeds below that are indistinguishable from zero anyway. In 
    practice we will use as the transition velocity the \e larger of this value 
    and the velocity constraint tolerance. **/
    void setDefaultFrictionTransitionVelocity(Real vTransition) {
        SimTK_ERRCHK1_ALWAYS(vTransition>=0,
        "SemiExplicitEulerTimeStepper::setDefaultFrictionTransitionVelocity()",
        "The friction transition velocity must be nonnegative but was %g.",
        vTransition);
        m_defaultTransitionVelocity = vTransition;
    }

    Real getAccuracy() const {return m_accuracy;}
    Real getConstraintTol() const {return m_consTol;}

    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the velocity constraint 
    tolerance. **/
    Real getDefaultImpactCaptureVelocity() const 
    {   return m_defaultCaptureVelocity; }
    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the impact cature velocity. **/
    Real getDefaultImpactMinCORVelocity() const 
    {   return m_defaultMinCORVelocity; }
    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the velocity constraint 
    tolerance. **/
    Real getDefaultFrictionTransitionVelocity() const 
    {   return m_defaultTransitionVelocity; }

    /** Return the value actually being used as the default impact capture
    velocity. **/
    Real getDefaultImpactCaptureVelocityInUse() const 
    {   return std::max(m_defaultCaptureVelocity, m_consTol); }
    /** Return the value actually being used as the default impact minimum
    coefficient of restitution velocity. **/
    Real getDefaultImpactMinCORVelocityInUse() const 
    {   return std::max(m_defaultMinCORVelocity, 
                        getDefaultImpactCaptureVelocityInUse()); }
    /** Return the value actually being used as the default sliding-to-rolling
    friction transition velocity. **/
    Real getDefaultFrictionTransitionVelocityInUse() const 
    {   return std::max(m_defaultTransitionVelocity, m_consTol); }

    /** Initialize the TimeStepper's internally maintained state to a copy
    of the given state. **/
    void initialize(const State& initState);
    const State& getState() const {return m_state;}
    State& updState() {return m_state;}
    Real getTime() const {return m_state.getTime();}

    /** Advance to the indicated time in one or more steps, using repeated
    induced impacts. **/
    Integrator::SuccessfulStepStatus stepTo(Real time);

private:
    // Determine which constraints will be involved for this step.
    void findProximalConstraints(const State&);
    // Enable all proximal constraints, disable all distal constraints, 
    // reassigning multipliers if needed. Returns true if anything changed.
    bool enableProximalConstraints(State&);
    // After constraints are enabled, gather up useful info about them.
    void collectConstraintInfo(const State& s);
    // Calculate velocity-dependent coefficients of restitution and friction
    // and apply combining rules for dissimilar materials.
    void calcCoefficientsOfFriction(const State&, const Vector& verr);
    void calcCoefficientsOfRestitution(const State&, const Vector& verr);

    // Easy if there are no constraints active.
    void takeUnconstrainedStep(State& s, Real h);

    // Given a velocity constraint error, determine if any of its entries
    // indicate that an impact is occurring.
    bool isImpact(const State& s, const Vector& verr) const;

    // Adjust given verr to reflect Newton restitution. 
    bool applyNewtonRestitutionIfAny(const State&, Vector& verr) const;

    // Adjust given compression impulse to include Poisson restitution impulse.
    // Note which contacts are expanding.
    bool applyPoissonRestitutionIfAny(const State&, Vector& impulse,
                                      Array_<int>& expanders) const;

    bool calcExpansionImpulseIfAny(const State& s, const Array_<int>& impacters,
                                   const Vector& compressionImpulse,
                                   Vector& expansionImpulse,
                                   Array_<int>& expanders) const; 

    // This phase uses all the proximal constraints and should use a starting
    // guess for impulse saved from the last step if possible.
    bool doCompressionPhase(const State&, const Vector& eps,
                            Vector& compressionImpulse);
    // This phase uses all the proximal constraints, but we expect the result
    // to be zero unless expansion causes new violations.
    bool doExpansionPhase(const State&, const Vector& eps,
                          Vector& reactionImpulse);
    bool doInducedImpactRound(const State&, const Vector& eps,
                              Vector& impulse);
    bool anyPositionErrorsViolated(const State&, const Vector& perr) const;

    // This phase uses only holonomic constraints, and zero is a good initial
    // guess for the (hopefully small) position correction.
    bool doPositionCorrectionPhase(const State&, const Vector& eps,
                                   Vector& positionImpulse);


private:
    const MultibodySystem&      m_mbs;
    Real                        m_accuracy;
    Real                        m_consTol;
    RestitutionModel            m_restitutionModel;
    Solver                      m_solver;
    PositionProjectionMethod    m_projectionMethod;

    Real                        m_defaultCaptureVelocity;
    Real                        m_defaultMinCORVelocity;
    Real                        m_defaultTransitionVelocity;
    Real                        m_defaultSignificantForce;

    //TODO: make this selectable
    PGSImpulseSolver m_pgsSolver;

    // Persistent runtime data.
    State                       m_state;

    // Step temporaries.
    Matrix                      m_GMInvGt; // G M\ ~G
    Vector                      m_D; // soft diagonal

    Array_<UnilateralContactIndex>      m_proximalUniContacts, 
                                        m_distalUniContacts;
    Array_<StateLimitedFrictionIndex>   m_proximalStateLtdFriction,
                                        m_distalStateLtdFriction;

    // This is for use in the no-impact phase where all proximals participate.
    Array_<MultiplierIndex>                         m_allParticipating;

    // These lists are for use in impact phases.
    Array_<MultiplierIndex>                         m_participating;
    Array_<ImpulseSolver::UncondRT>                 m_unconditional;
    Array_<ImpulseSolver::UniContactRT>             m_uniContact; // with fric
    Array_<ImpulseSolver::UniSpeedRT>               m_uniSpeed;
    Array_<ImpulseSolver::BoundedRT>                m_bounded;
    Array_<ImpulseSolver::ConstraintLtdFrictionRT>  m_consLtdFriction;
    Array_<ImpulseSolver::StateLtdFrictionRT>       m_stateLtdFriction;

    // These lists are for use in position projection and include only
    // holonomic constraints (the unilateral contacts have friction off).
    Array_<MultiplierIndex>                         m_posParticipating;
    Array_<ImpulseSolver::UncondRT>                 m_posUnconditional;
    Array_<ImpulseSolver::UniContactRT>             m_posUniContact; // no fric
    Array_<ImpulseSolver::UniSpeedRT>               m_posNoUniSpeed;
    Array_<ImpulseSolver::BoundedRT>                m_posNoBounded;
    Array_<ImpulseSolver::ConstraintLtdFrictionRT>  m_posNoConsLtdFriction;
    Array_<ImpulseSolver::StateLtdFrictionRT>       m_posNoStateLtdFriction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_SEMI_EXPLICIT_EULER_TIME_STEPPER_H_

