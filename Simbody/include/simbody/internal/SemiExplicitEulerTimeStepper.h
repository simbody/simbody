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

/** A low-accuracy, high performance, velocity-level time stepper for
models containing unilateral rigid contacts or other conditional constraints.

At each step, this TimeStepper determines which subset of the multibody system's
conditional constraints may affect behavior; those are called "proximal"
constraints. It then formulates the resulting
equations and invokes an ImpulseSolver to find constraint-space impulses that
solve the constraint equations and satisfy all the constraint inequalities.
Generally there are multiple solutions possible, and different ImpulseSolver
objects use different criteria.

A variety of options are available for different methods of handling impacts,
to facilitate comparison of methods. For production, we recommend using the
default options.
**/

class SimTK_SIMBODY_EXPORT SemiExplicitEulerTimeStepper {
public:
    /** If an impact occurs at a contact where the coefficient of restitution
    (COR) is non zero, this option determines how we process the restitution
    impulse. The default \c Poisson uses Poisson's hypothesis, which treates the
    COR as a ratio of impulses, providing good physical behavior. The
    alternative \c Newton uses Newton's hypothesis which treats the COR as a 
    ratio of velocities; this is more commonly used in other packages and can be
    somewhat faster, but may exhibit blatantly nonphysical behavior such as 
    energy gain during an impact. Poisson and Newton are equivalent in most
    simple impact circumstances but can differ substantially in a multibody 
    impact involving multiple coupled simultaneous collisions. The final 
    alternative is \c NoRestitution, meaning that all CORs are 
    treated as zero so impacts are always maximally dissipative. **/
    enum RestitutionModel {Poisson=0, Newton=1, NoRestitution=2};
    enum InducedImpactModel {Simultaneous=0, Sequential=1, Mixed=2};
    enum PositionProjectionMethod {Bilateral=0,Unilateral=1,
                                   NoPositionProjection=2};
    enum ImpulseSolverType {PLUS=0, PGS=1};


    explicit SemiExplicitEulerTimeStepper(const MultibodySystem& mbs);

    /** The contained ImpulseSolver will be destructed here; don't reference 
    it afterwards! **/
    ~SemiExplicitEulerTimeStepper() {
        clearImpulseSolver();
    }

    /** Initialize the TimeStepper's internally maintained state to a copy
    of the given state; allocate and initialize the ImpulseSolver if there
    isn't one already. **/
    void initialize(const State& initState);
    /** Get access to the TimeStepper's internally maintained State. **/
    const State& getState() const {return m_state;}
    /** Get writable access to the TimeStepper's internally maintained
    State. Modifying this will change the course of the simulation. **/
    State& updState() {return m_state;}
    /** Shortcut to getting the current time from the TimeStepper's internally
    maintained State. **/
    Real getTime() const {return m_state.getTime();}

    /** synonym for getState. **/
    const State& getAdvancedState() const {return m_state;}

    /** synonym for updState. **/
    State& updAdvancedState() {return m_state;}

    /** synonym for getTime. **/
    Real getAdvancedTime() const {return m_state.getTime();}

    /** Advance to the indicated time in one or more steps, using repeated
    induced impacts. **/
    Integrator::SuccessfulStepStatus stepTo(Real time);

    /** Set integration accuracy; requires variable length steps. **/
    void setAccuracy(Real accuracy) {m_accuracy=accuracy;}
    /** Set the tolerance to which constraints must be satisfied. **/
    void setConstraintTolerance(Real consTol) {m_consTol=consTol;}
   
    void setRestitutionModel(RestitutionModel restModel)
    {   m_restitutionModel = restModel; }
    RestitutionModel getRestitutionModel() const 
    {   return m_restitutionModel; }
     
    void setInducedImpactModel(InducedImpactModel indModel)
    {   m_inducedImpactModel = indModel; }
    InducedImpactModel getInducedImpactModel() const 
    {   return m_inducedImpactModel; }
    
    /** Limit the number of induced impact rounds per time step. The 
    effect depends on the InducedImpactModel. If Simultaneous, this setting
    is ignored since there is always just a single round that fully resolves
    the impact. If Sequential, after this many rounds the remainder of the
    impact is left unresolved to be dealt with at the next time step (meaning
    penetration will be permitted). If Mixed, this many sequential impact rounds
    are executed, followed by a final round which resolves any remaining 
    penetration as for Simultaneous, meaning that all remaining induced impacts
    are treated as though they have zero coefficient of restitution. **/
    void setMaxInducedImpactsPerStep(int maxInduced) {
        SimTK_APIARGCHECK1_ALWAYS(maxInduced>=0, "SemiExplicitEulerTimeStepper",
            "setMaxInducedImpactsPerStep", "Illegal argument %d", maxInduced);
        m_maxInducedImpactsPerStep = maxInduced;
    }
    int getMaxInducedImpactsPerStep() const 
    {   return m_maxInducedImpactsPerStep; }
     
    void setPositionProjectionMethod(PositionProjectionMethod projMethod)
    {   m_projectionMethod = projMethod; }
    PositionProjectionMethod getPositionProjectionMethod() const 
    {   return m_projectionMethod; }
    
    void setImpulseSolverType(ImpulseSolverType solverType) {
        if (m_solverType != solverType) {
            // The new solver will get allocated in initialize().
            clearImpulseSolver();
            m_solverType = solverType; 
        }
    }
    ImpulseSolverType getImpulseSolverType() const 
    {   return m_solverType; }

    /** Set the impact capture velocity to be used by default when a contact
    does not provide its own. This is the impact velocity below which the
    coefficient of restitution is to be treated as zero. This avoids a Zeno's
    paradox of ever-tinier time-wasting impacts. A capture velocity should
    be significantly greater than the velocity constraint tolerance
    since speeds below that are indistinguishable from zero anyway. In 
    practice we will use as the capture velocity the \e larger of this value 
    and twice the velocity constraint tolerance currently in effect. **/
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
    and twice the velocity constraint tolerance currently in effect. **/
    void setDefaultFrictionTransitionVelocity(Real vTransition) {
        SimTK_ERRCHK1_ALWAYS(vTransition>=0,
        "SemiExplicitEulerTimeStepper::setDefaultFrictionTransitionVelocity()",
        "The friction transition velocity must be nonnegative but was %g.",
        vTransition);
        m_defaultTransitionVelocity = vTransition;
    }

    /** Set the threshold below which we can ignore forces. This gives the
    %TimeStepper to consider any force with smaller magnitude to be zero. 
    The default value for this is SimTK::SignificantReal, about 1e-14 in double
    precision. **/
    void setMinSignificantForce(Real minSignificantForce) {
        SimTK_ERRCHK1_ALWAYS(minSignificantForce>0,
        "SemiExplicitEulerTimeStepper::setMinSignificantForce()",
        "The minimum significant force magnitude must be greater than zero "
        "but was %g.", minSignificantForce);
        m_minSignificantForce = minSignificantForce;
    }
    Real getMinSignificantForce() const 
    {   return m_minSignificantForce; }

    /** Return the integration accuracy setting. This has no effect unless
    you are running in variable time step mode. **/
    Real getAccuracyInUse() const {return m_accuracy;}

    /** Return the tolerance to which we require constraints to be satisfied.
    This applies even if we are not controlling overall integration accuracy.
    It is also used as a "minimum meaningful velocity" value that overrides
    other velocity thresholds if they have been set to smaller values. **/
    Real getConstraintToleranceInUse() const {return m_consTol;}

    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the velocity constraint 
    tolerance. @see getDefaultImpactCaptureVelocityInUse() **/
    Real getDefaultImpactCaptureVelocity() const 
    {   return m_defaultCaptureVelocity; }
    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the impact cature velocity. 
    @see getDefaultImpactMinCORVelocityInUse() **/
    Real getDefaultImpactMinCORVelocity() const 
    {   return m_defaultMinCORVelocity; }
    /** Return the value set for this parameter, but the actual value used
    during execution will be no smaller than the velocity constraint 
    tolerance. @see getDefaultFrictionTransitionVelocityInUse() **/
    Real getDefaultFrictionTransitionVelocity() const 
    {   return m_defaultTransitionVelocity; }

    /** Return the value actually being used as the default impact capture
    velocity. **/
    Real getDefaultImpactCaptureVelocityInUse() const 
    {   return std::max(m_defaultCaptureVelocity, 2*m_consTol); }
    /** Return the value actually being used as the default impact minimum
    coefficient of restitution velocity. **/
    Real getDefaultImpactMinCORVelocityInUse() const 
    {   return std::max(m_defaultMinCORVelocity, 
                        getDefaultImpactCaptureVelocityInUse()); }
    /** Return the value actually being used as the default sliding-to-rolling
    friction transition velocity. **/
    Real getDefaultFrictionTransitionVelocityInUse() const 
    {   return std::max(m_defaultTransitionVelocity, 2*m_consTol); }

    /** Get access to the MultibodySystem for which this %TimeStepper was
    constructed. **/
    const MultibodySystem& getMultibodySystem() const {return m_mbs;}


    /** (Advanced) Get direct access to the ImpulseSolver. **/
    const ImpulseSolver& getImpulseSolver() const {
        SimTK_ERRCHK_ALWAYS(m_solver!=0, 
            "SemiExplicitEulerTimeStepper::getImpulseSolver()",
            "No solver is currently allocated.");
        return *m_solver;
    }
    /** (Advanced) Set your own ImpulseSolver; the %TimeStepper takes over
    ownership so don't delete afterwards! **/
    void setImpulseSolver(ImpulseSolver* impulseSolver) {
        clearImpulseSolver();
        m_solver = impulseSolver;
    }
    /** (Advanced) Delete the existing ImpulseSolver if any. **/
    void clearImpulseSolver() {
        delete m_solver; m_solver=0;
    }

    /** Get human-readable string representing the given enum value. **/
    static const char* getRestitutionModelName(RestitutionModel rm);
    /** Get human-readable string representing the given enum value. **/
    static const char* getInducedImpactModelName(InducedImpactModel iim);
    /** Get human-readable string representing the given enum value. **/
    static const char* getPositionProjectionMethodName
       (PositionProjectionMethod ppm);
    /** Get human-readable string representing the given enum value. **/
    static const char* getImpulseSolverTypeName(ImpulseSolverType ist);

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
    void calcCoefficientsOfRestitution(const State&, const Vector& verr,
                                       bool disableRestitution);

    // Easy if there are no constraints active.
    void takeUnconstrainedStep(State& s, Real h);

    // If we're in Newton restitution mode, calculating the verr change
    // that is needed to represent restitution. Output must already be
    // the same size as verr on entry if we're in Newton mode.
    bool calcNewtonRestitutionIfAny(const State&, const Vector& verr,
                                    Vector& newtonVerr) const;

    // Adjust given compression impulse to include Poisson restitution impulse.
    // Note which contacts are expanding.
    bool applyPoissonRestitutionIfAny(const State&, Vector& impulse,
                                      Array_<int>& expanders) const;

    bool calcExpansionImpulseIfAny(const State& s, const Array_<int>& impacters,
                                   const Vector& compressionImpulse,
                                   Vector& expansionImpulse,
                                   Array_<int>& expanders) const;

    // Perform a simultaneous impact if needed. All proximal constraints are 
    // dealt with so after this call there will be no more impacters, and no 
    // unapplied expansion impulses. For Poisson restitution this may be a 
    // compression/expansion impulse pair (and rarely a final compression
    // round to correct expanders that were forced back into impacting). 
    // For Newton restitution only a single impulse round is calculated.
    // Returns the number of impulse rounds actually taken, usually zero.
    int performSimultaneousImpact(const State& state, 
                                   Vector&      verr, // in/out
                                   Vector&      totalImpulse);

    // We identify impacters, observers, and expanders then perform a single
    // impulse calculation that ignores the observers. On return there may
    // be former observers and expanders that now have impacting approach
    // velocities so will be impacters on the next round. For Poisson
    // restitution, there may be expansion impulses that have not yet been 
    // applied; those contacts will be expanders on the next round.
    bool performInducedImpactRound(const Vector& verr, 
                                   const Vector& expansionImpulse);

    void classifyUnilateralContactsForSequentialImpact
       (const Vector&                           verr,
        const Vector&                           expansionImpulse,
        bool                                    includeAllProximals,
        bool                                    expansionOnly,
        Array_<ImpulseSolver::UniContactRT>&    uniContacts, 
        Array_<int>&                            impacters,
        Array_<int>&                            expanders,
        Array_<int>&                            observers,
        Array_<MultiplierIndex>&                participaters,
        Array_<MultiplierIndex>&                expanding) const;


    void classifyUnilateralContactsForSimultaneousImpact
       (const Vector&                           verr,
        const Vector&                           expansionImpulse,
        Array_<ImpulseSolver::UniContactRT>&    uniContacts, 
        Array_<int>&                            impacters,
        Array_<int>&                            expanders,
        Array_<int>&                            observers,
        Array_<MultiplierIndex>&                participaters,
        Array_<MultiplierIndex>&                expanding) const;

    // This phase uses all the proximal constraints and should use a starting
    // guess for impulse saved from the last step if possible.
    bool doCompressionPhase(const State&    state, 
                            Vector&         verrStart, // in/out
                            Vector&         verrApplied, // in/out
                            Vector&         compressionImpulse);
    // This phase uses all the proximal constraints, but we expect the result
    // to be zero unless expansion causes new violations.
    bool doExpansionPhase(const State&  state, 
                          const Array_<MultiplierIndex>& expanding,
                          Vector&       expansionImpulse,
                          Vector&       verrStart, // in/out
                          Vector&       reactionImpulse);
    bool doInducedImpactRound(const State&  state, 
                              const Array_<MultiplierIndex>& expanding,
                              Vector&       expansionImpulse,
                              Vector&       verrStart, // in/out
                              Vector&       impulse);
    bool anyPositionErrorsViolated(const State&, const Vector& perr) const;

    // This phase uses only holonomic constraints, and zero is a good initial
    // guess for the (hopefully small) position correction.
    bool doPositionCorrectionPhase(const State& state, 
                                   Vector&      pverr, // in/out
                                   Vector&      positionImpulse);


private:
    const MultibodySystem&      m_mbs;
    Real                        m_accuracy;
    Real                        m_consTol;
    RestitutionModel            m_restitutionModel;
    InducedImpactModel          m_inducedImpactModel;
    int                         m_maxInducedImpactsPerStep;
    PositionProjectionMethod    m_projectionMethod;
    ImpulseSolverType           m_solverType;


    Real                        m_defaultCaptureVelocity;
    Real                        m_defaultMinCORVelocity;
    Real                        m_defaultTransitionVelocity;
    Real                        m_minSignificantForce;

    ImpulseSolver*              m_solver;

    // Persistent runtime data.
    State                       m_state;
    Vector                      m_emptyVector; // don't change this!

    // Step temporaries.
    Matrix                      m_GMInvGt; // G M\ ~G
    Vector                      m_D; // soft diagonal
    Vector                      m_deltaU;
    Vector                      m_verr;
    Vector                      m_totalImpulse;
    Vector                      m_impulse;
    Vector                      m_genImpulse; // ~G*impulse

    Array_<UnilateralContactIndex>      m_proximalUniContacts, 
                                        m_distalUniContacts;
    Array_<StateLimitedFrictionIndex>   m_proximalStateLtdFriction,
                                        m_distalStateLtdFriction;

    // This is for use in the no-impact phase where all proximals participate.
    Array_<MultiplierIndex>                         m_allParticipating;

    // These are for use in impact phases.
    Array_<MultiplierIndex>                         m_participating;
    Array_<MultiplierIndex>                         m_expanding;
    Vector                                          m_expansionImpulse;
    Vector                                          m_newtonRestitutionVerr;
    Array_<int>                                     m_impacters;
    Array_<int>                                     m_expanders;
    Array_<int>                                     m_observers;

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

