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

#include "simbody/internal/common.h"
#include "simbody/internal/SemiExplicitEulerTimeStepper.h"
#include "simbody/internal/ConditionalConstraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <iostream>
using std::cout; using std::endl;

/* Implementation of SemiExplicitEulerTimeStepper. */

namespace SimTK {


//------------------------------------------------------------------------------
//                               STEP TO
//------------------------------------------------------------------------------
Integrator::SuccessfulStepStatus SemiExplicitEulerTimeStepper::
stepTo(Real time) {
    // Abbreviations.
    const MultibodySystem&              mbs    = m_mbs;
    const SimbodyMatterSubsystem&       matter = mbs.getMatterSubsystem();
    State&                              s      = m_state;

    const Real t0 = m_state.getTime();
    const Real h = time - t0;    // max timestep

    // Kinematics should already be realized so this won't do anything.
    mbs.realize(s, Stage::Position); 
    // Determine which constraints will be involved for this step.
    findProximalConstraints(s);
    // Enable all proximal constraints, reassigning multipliers if needed.
    enableProximalConstraints(s);
    collectConstraintInfo(s);

    mbs.realize(s, Stage::Velocity);
    const Vector u0 = s.getU(); // save

    const Vector& verr0 = s.getUErr();
    if (!verr0.size()) {
        takeUnconstrainedStep(s, h);
        return Integrator::ReachedScheduledEvent;
    }

    const int m = verr0.size();

    // Friction coefficient is fixed by initial slip velocity and doesn't
    // change during impact processing even the the slip velocity will change.
    // The logic is that it takes time for surface asperities to engage or
    // separate and no time is going by during an impact.
    calcCoefficientsOfFriction(s, verr0);
    calcCoefficientsOfRestitution(s, verr0);

    // If we're in Newton mode, or if a contact specifies Newton restitution,
    // then modify the appropriate verr's here.
    Vector verr = verr0;
    const bool anyNewton = applyNewtonRestitutionIfAny(s, verr);

    // Calculate the constraint compliance matrix GM\~G.
    matter.calcProjectedMInv(s, m_GMInvGt); // m X m

    Array_<int> expanders;
    Array_<int> impacters;
    Array_<int> whichFrictional;

    //cout << "Precomp verr=" << verr << endl;
    Vector totalImpulse(m, Real(0));
    Vector impulse(m), expansionImpulse(m, Real(0));

    bool anyImpact = false;
    while (true) {
        // Calculate participating constraints. All but:
        //   - ignore "observers" (nonnegative impact velocity)
        //   - ignore normal constraints for "expanders" (friction stays)
        // Then impacters=all contacts with negative impact velocity
        //      whichFrictional=all impacter & expander frictional elements
        m_participating.clear();;
        impacters.clear();
        whichFrictional.clear();
        for (unsigned i=0; i < m_uniContact.size(); ++i) {
            ImpulseSolver::UniContactRT& rt = m_uniContact[i];
            const MultiplierIndex mx = rt.m_Nk;
            if (expansionImpulse[mx] != 0) {
                rt.m_type = ImpulseSolver::Known;
                rt.m_knownPi = expansionImpulse[mx];
                if (rt.hasFriction()) { // friction still participates
                    for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                        m_participating.push_back(rt.m_Fk[j]);
                }
                continue; //expander
            }
            if (verr[mx] > -m_consTol) {
                rt.m_type = ImpulseSolver::Observe;
                continue; //observer
            }
            rt.m_type = ImpulseSolver::Participate;
            m_participating.push_back(mx);
            if (rt.hasFriction()) {
                for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                    m_participating.push_back(rt.m_Fk[j]);
            }
            impacters.push_back(i);
        }

        // TODO: unconditionals, unispeed, bounded, cons-ltd, state-ltd friction

        if (m_participating.empty())
            break;

        //--------------------------------------------------------------------------
        anyImpact = true;
        impulse = expansionImpulse;
        #ifndef NDEBUG
        printf("impacters: "); cout << impacters << endl;
        printf("whichFrictional: "); cout << whichFrictional << endl;
        printf("participating: "); cout << m_participating << endl;
        cout << "Preimpact impulse=" << impulse << endl;
        cout << "             verr=" << verr << endl;
        #endif
        doInducedImpactRound(s, verr, impulse);
        //--------------------------------------------------------------------------
        impulse -= expansionImpulse; // already applied it

        // Calculate the new expansion impulse.
        const bool anyExpansion = calcExpansionImpulseIfAny(s, impacters, impulse,
                                    expansionImpulse,expanders);
        #ifndef NDEBUG
        cout << "Postimpact impulse=" << impulse << endl;
        cout << "Expansion impulse=" << expansionImpulse << endl;
        cout << "Expanders: " << expanders << endl;
        #endif

        if (anyExpansion) impulse += expansionImpulse;
        totalImpulse += impulse;
        verr -= m_GMInvGt*impulse;
        #ifndef NDEBUG
        cout << "Postimpact verr=" << verr << endl;
        #endif

        if (!anyExpansion)
            break;
        calcCoefficientsOfRestitution(s, verr);
    }

    if (anyImpact) {
        Vector genImpulse, deltaU, expVerr; 
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, -totalImpulse, genImpulse);
        // gen impulse to delta-u
        matter.multiplyByMInv(s, genImpulse, deltaU);
        s.updU() += deltaU;
        mbs.realize(s, Stage::Velocity); // updates rotational forces
        #ifndef NDEBUG
        cout << "Impact: imp " << totalImpulse << "-> du " << deltaU << endl;
        #endif
    } else {
        #ifndef NDEBUG
        cout << "No impact.\n";
        #endif
    }

    // Evaluate applied forces and get reference to them. These include
    // gravity but not centrifugal forces.
    mbs.realize(s, Stage::Dynamics);
    const Vector&              f = mbs.getMobilityForces(s, Stage::Dynamics);
    const Vector_<SpatialVec>& F = mbs.getRigidBodyForces(s, Stage::Dynamics);

    // Calculate udotExt = M\(f + ~J*(F-C)) where C are centrifugal forces.
    // This is the unconstrained acceleration.
    Vector udotExt; Vector_<SpatialVec> A_GB;
    matter.calcAccelerationIgnoringConstraints(s,f,F,udotExt,A_GB);

    // Calculate verrExt = G*deltaU; the end-of-step constraint error due to 
    // external and centrifugal forces.
    Vector verrExt;
    matter.multiplyByG(s, h*udotExt, verrExt);
    verrExt += verr;

    // Make all unilateral contacts participate for this phase.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        rt.m_type = ImpulseSolver::Participate;
    }
    doCompressionPhase(s, verrExt, impulse);
    #ifndef NDEBUG
    cout << "Dynamics verr=" << verrExt << endl;
    cout << "      impulse=" << impulse << endl;
    #endif
    totalImpulse += impulse;
    s.updMultipliers() = totalImpulse/h;
    // Calculate constraint forces ~G*lambda (body frcs Fc, mobility frcs fc).
    Vector_<SpatialVec> Fc; Vector fc; 
    matter.calcConstraintForcesFromMultipliers(s,impulse/h, Fc, fc);

    // Now calculate udot = M\(f-fc + ~J*(F-Fc-C)).
    Vector udot;
    matter.calcAccelerationIgnoringConstraints(s,f-fc,F-Fc,udot,A_GB);

    // Update auxiliary states z, invalidating Stage::Dynamics.
    s.updZ() += h*s.getZDot();

    // Update u from deltaU, invalidating Stage::Velocity. 
    s.updU() += h*udot;

    s.updUDot() = (s.getU()-u0)/h;
    //TODO: need to calculate reaction forces from (udot,lambda), and
    // raise state's stage to Acceleration.

    // Done with velocity update. Now calculate qdot, possiblity including
    // an additional position error correction term.
    Vector qdot;
    const Vector& perr0 = s.getQErr();
    if (!anyPositionErrorsViolated(s, perr0)) {
        matter.multiplyByN(s,false,s.getU(),qdot);
    } else {
        // Don't include quaternions for position correction. 
        Vector posImpulse, genImpulse, posVerr, deltaU;
        posVerr.resize(s.getNUErr()); posVerr.setToZero();
        const int nQuat = matter.getNumQuaternionsInUse(s);
        posVerr(0, perr0.size()-nQuat) = perr0(0, perr0.size()-nQuat)/h;

        //----------------------------------------------------------------------
        // Calculate impulse and then deltaU=M\~G*impulse such that -h*deltaU 
        // will eliminate position errors, respecting only position constraints.
        doPositionCorrectionPhase(s, posVerr, posImpulse);
        //----------------------------------------------------------------------
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, posImpulse, genImpulse);
        // gen impulse to deltaU (watch sign)
        matter.multiplyByMInv(s, genImpulse, deltaU);

        // convert corrected u to qdot (note we're not changing u)
        matter.multiplyByN(s,false,s.getU()-deltaU, qdot);
    }

    // We have qdot, now update q, fix quaternions, update time.
    s.updQ() += h*qdot; // invalidates Stage::Position
    matter.normalizeQuaternions(s);
    s.updTime() += h;   // invalidates Stage::Time

    // Return from step with kinematics realized. Note that we may have
    // broken the velocity constraints by updating q, but we won't fix that
    // until the next step. Also position constraints are only imperfectly
    // satisfied by the correction above.
    mbs.realize(s, Stage::Velocity);
    #ifndef NDEBUG
    printf("end of step (%g,%g): verr=", t0,s.getTime());
    cout << s.getUErr() << endl;
    #endif

    return Integrator::ReachedScheduledEvent;
}

//------------------------------------------------------------------------------
//                             INITIALIZE
//------------------------------------------------------------------------------
void SemiExplicitEulerTimeStepper::initialize(const State& initState) {
    m_state = initState;
    m_mbs.realize(m_state, Stage::Acceleration);
}

//------------------------------------------------------------------------------
//                         FIND PROXIMAL CONSTRAINTS
//------------------------------------------------------------------------------
// The goal here is to examine the given State to classify as "proximal" or 
// "distal":
//      - each unilateral contact constraint and associated friction 
//        constraint, and
//      - each state-limited friction constraint.
// Only proximal constraints may participate in the generation of constraint 
// forces during a step; distal constraints have no effect whatsoever. In 
// general only a subset of the proximal unilateral constraints will actually 
// generate forces, but all must be considered.
//
// Velocity-level (nonholonomic) unilateral constraints are always considered
// proximal, as are unconditional constraints and their associated friction
// constraints.
//
// Simbody will be asked to generate constraint equations for all the proximal
// constraints, but not for the distal constraints.
//
// Algorithm:
// (1) Inspect each unilateral contact constraint in the MultibodySystem, 
// marking it proximal if its perr() value is within tolerance of violating 
// its defining unilateral inequality, distal otherwise. Typically, perr() is 
// a signed distance function and the defining inequality is perr() >= 0, so 
// the proximal condition is perr() <= tol. Associated friction constraints,
// if any, are marked proximal or distal to match the contact constraint.
//
// (2) Inspect each StateLimitedFriction constraint. Mark the friction element
// proximal if its known limiting force is non-zero (or above a small 
// threshold), otherwise it is distal.
//
// TODO: Constraint-limited friction should only be proximal if the limiting
// constraint is enabled.
void SemiExplicitEulerTimeStepper::
findProximalConstraints(const State& s) { //TODO: redo
    m_proximalUniContacts.clear();      m_distalUniContacts.clear();
    m_proximalStateLtdFriction.clear(); m_distalStateLtdFriction.clear();

    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    const int nUniContacts  = matter.getNumUnilateralContacts();
    const int nLtdFrictions = matter.getNumStateLimitedFrictions();

    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        if (contact.isProximal(s, m_consTol)) // may be scaled
            m_proximalUniContacts.push_back(ux);
        else m_distalUniContacts.push_back(ux);
    }

    for (StateLimitedFrictionIndex fx(0); fx < nLtdFrictions; ++fx) {
        const StateLimitedFriction& fric = matter.getStateLimitedFriction(fx);
        if (fric.getNormalForceMagnitude(s) >= m_defaultSignificantForce) 
            m_proximalStateLtdFriction.push_back(fx);
        else m_distalStateLtdFriction.push_back(fx);
    }

    #ifndef NDEBUG
    cout<<"Proximal unilateral contacts: "<< m_proximalUniContacts << endl;
    cout<<"Distal unilateral contacts: "  << m_distalUniContacts   << endl;
    cout<<"Proximal state-ltd friction: " << m_proximalStateLtdFriction << endl;
    cout<<"Distal state-ltd friction: "   << m_distalStateLtdFriction   << endl;
    #endif
}


// Enable all proximal constraints, disable all distal constraints, 
// reassigning multipliers if needed. Returns true if any change was made.
// If not, the multipliers retain their meaning from the last time this
// State was advanced.
bool SemiExplicitEulerTimeStepper::
enableProximalConstraints(State& s) {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    // Record friction application points. This has to be done while Position 
    // stage is still valid. TODO: This is a hack and shouldn't be necessary.
    Array_<Vec3> uniPosInfo(m_proximalUniContacts.size());
    Array_<Vec3> ltdPosInfo(m_proximalStateLtdFriction.size());
    for (unsigned i=0; i < m_proximalUniContacts.size(); ++i) {
        const UnilateralContactIndex cx = m_proximalUniContacts[i];
        uniPosInfo[i] = matter.getUnilateralContact(cx).getPositionInfo(s);
    }
    for (unsigned i=0; i < m_proximalStateLtdFriction.size(); ++i) {
        const StateLimitedFrictionIndex fx = m_proximalStateLtdFriction[i];
        ltdPosInfo[i] = matter.getStateLimitedFriction(fx).getPositionInfo(s);
    }

    bool changed = false;

    // Disable non-proximal constraints if they were previously enabled.
    for (unsigned i=0; i < m_distalUniContacts.size(); ++i) {
        const UnilateralContactIndex cx = m_distalUniContacts[i];
        const UnilateralContact& contact = matter.getUnilateralContact(cx);
        if (contact.disable(s)) changed = true;
    }
    for (unsigned i=0; i < m_distalStateLtdFriction.size(); ++i) {
        const StateLimitedFrictionIndex fx = m_distalStateLtdFriction[i];
        const StateLimitedFriction& fric = matter.getStateLimitedFriction(fx);
        if (fric.disable(s)) changed = true;
    }

    // Now enable the proximal constraints if they were disabled.
    for (unsigned i=0; i < m_proximalUniContacts.size(); ++i) {
        const UnilateralContactIndex cx = m_proximalUniContacts[i];
        const UnilateralContact& contact = matter.getUnilateralContact(cx);
        contact.setInstanceParameter(s, uniPosInfo[i]);
        if (contact.enable(s)) changed = true;
    }
    for (unsigned i=0; i < m_proximalStateLtdFriction.size(); ++i) {
        const StateLimitedFrictionIndex fx = m_proximalStateLtdFriction[i];
        const StateLimitedFriction& fric = matter.getStateLimitedFriction(fx);
        fric.setInstanceParameter(s, ltdPosInfo[i]);
        if (fric.enable(s)) changed = true;
    }

    // TODO: Note that we always have to move the friction application points
    // which is an Instance stage change; shouldn't be.
    m_mbs.realize(s, Stage::Instance); // assign multipliers

    return changed;
}

// Create lists of low-level constraint properties suited for analysis by
// an Impulse solver. Only proximal constraints are considered, and those
// must already have been enabled so that we can collect their Simbody-assigned
// multipliers here.
void SemiExplicitEulerTimeStepper::
collectConstraintInfo(const State& s) { //TODO: redo
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    // These are invalid now but we won't fill them out here.
    m_participating.clear();  m_posParticipating.clear();

    // This will be filled out here since it can't change later.
    m_allParticipating.clear();

    // These are the constraints for contact & impact.
    m_unconditional.clear();
    m_uniContact.clear();
    m_uniSpeed.clear();
    m_bounded.clear();
    m_consLtdFriction.clear();
    m_stateLtdFriction.clear();
    
    // These are just the holonomic constraints. 
    m_posUnconditional.clear();
    m_posUniContact.clear();

    const int nConstraints = matter.getNumConstraints();
    for (ConstraintIndex cx(0); cx < nConstraints; ++cx) {
        // TODO: sort out the enabled, unconditional ones and fill in
        // m_unconditional and m_posUnconditional.
        // m_allParticipating.push_back(...);
    }

    for (unsigned puc=0; puc < m_proximalUniContacts.size(); ++puc) {
        const UnilateralContactIndex cx = m_proximalUniContacts[puc];
        const UnilateralContact& contact = matter.getUnilateralContact(cx);
        m_uniContact.push_back(); // default constructed
        ImpulseSolver::UniContactRT& rt = m_uniContact.back();
        rt.m_ucx = cx;
        rt.m_Nk = contact.getContactMultiplierIndex(s);
        m_allParticipating.push_back(rt.m_Nk);
        if (contact.hasFriction(s)) {
            rt.m_Fk.resize(2);
            contact.getFrictionMultiplierIndices(s, rt.m_Fk[0], rt.m_Fk[1]);
            m_allParticipating.push_back(rt.m_Fk[0]);
            m_allParticipating.push_back(rt.m_Fk[1]);
        }

        // Only the normal constraint is included for position projection.
        m_posUniContact.push_back(); // default constructed
        ImpulseSolver::UniContactRT& posRt = m_posUniContact.back();
        posRt.m_ucx = cx;
        posRt.m_Nk = rt.m_Nk;
    }

    //TODO: speed, bounded, consLtd, proximal stateLtd friction
}

void SemiExplicitEulerTimeStepper::takeUnconstrainedStep(State& s, Real h) {
    const SimbodyMatterSubsystem&   matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Acceleration);
    const Vector& udot = s.getUDot(); // grab before invalidated
    s.updZ() += h*s.getZDot(); // invalidates Stage::Dynamics
    s.updU() += h*udot;        // invalidates Stage::Velocity
    Vector qdot;
    matter.multiplyByN(s,false,s.getU(),qdot);
    s.updQ() += h*qdot;         // invalidates Stage::Position
    matter.normalizeQuaternions(s);
    s.updTime() += h;           // invalidates Stage::Time
    m_mbs.realize(s, Stage::Velocity);

    //TODO: prescribed motion, auto update state variables, events(?)
}


bool SemiExplicitEulerTimeStepper::
isImpact(const State& s, const Vector& verr) const {
    // These are the proximal unilateral constraint runtimes.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        const ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        if (rt.m_sign*verr[rt.m_Nk] < -m_consTol) { // TODO: sign?
            printf("IMPACT cuz verr[%d]=%g\n", rt.m_Nk, verr[rt.m_Nk]);
            return true;
        }
    }
    return false;
}


// Calculate velocity-dependent coefficients of friction for all the
// proximal friction constraints.
void SemiExplicitEulerTimeStepper::
calcCoefficientsOfFriction(const State& s, const Vector& verr) {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    // These are the proximal unilateral constraint runtimes.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        if (!rt.hasFriction())
            continue;
        const UnilateralContact& uni = matter.getUnilateralContact(rt.m_ucx);
        rt.m_effMu = uni.calcEffectiveCOF(s, 
                        getDefaultFrictionTransitionVelocity());
        SimTK_DEBUG3("Uni fric %d speed=%g -> mu=%g\n", (int)rt.m_ucx, 
                     uni.getSlipVelocity(s).norm(), rt.m_effMu);
    }

    // These are the proximal state-limited friction elements.
    for (unsigned i=0; i < m_stateLtdFriction.size(); ++i) {
        ImpulseSolver::StateLtdFrictionRT& rt = m_stateLtdFriction[i];
        const StateLimitedFriction& fric = 
            matter.getStateLimitedFriction(rt.m_slfx);
        rt.m_effMu = fric.calcEffectiveCOF(s, 
                        getDefaultFrictionTransitionVelocity());
        SimTK_DEBUG3("State ltd fric %d speed=%g -> mu=%g\n", (int)rt.m_slfx, 
                     fric.getSlipSpeed(s), rt.m_effMu);
    }

    // These are the constraint-limited friction elements, which are 
    // always proximal.
    //for (unsigned i=0; i < m_consLtdFriction.size(); ++i) {
    //    ImpulseSolver::ConstraintLtdFrictionRT& rt = m_consLtdFriction[i];
    //    const ConstraintLimitedFriction& fric = 
    //        matter.getConstraintLimitedFriction(rt.m_clfx);
    //    rt.m_effMu = fric.calcEffectiveCOF(s, 
    //                    getDefaultFrictionTransitionVelocity());
    //    SimTK_DEBUG3("State ltd fric %d speed=%g -> mu=%g\n", (int)rt.m_slfx, 
    //                 fric.getSlipSpeed(s), rt.m_effMu);
    //}

}

// Calculate velocity-dependent coefficients of restitution.
// TODO: apply combining rules for dissimilar materials.
void SemiExplicitEulerTimeStepper::
calcCoefficientsOfRestitution(const State& s, const Vector& verr) {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    // These are the proximal unilateral constraint runtimes.
    // TODO: only need this for participating contacts.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        const UnilateralContact& uni = matter.getUnilateralContact(rt.m_ucx);
        rt.m_effCOR = uni.calcEffectiveCOR(s, 
                        getDefaultImpactCaptureVelocityInUse(),
                        getDefaultImpactMinCORVelocityInUse());
        SimTK_DEBUG3("Uni contact %d verr=%g -> cor=%g\n", (int)rt.m_ucx, 
                     uni.getVerr(s), rt.m_effCOR);
    }
}

bool SemiExplicitEulerTimeStepper::
applyNewtonRestitutionIfAny(const State& s, Vector& verr) const {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    if (getRestitutionModel() != Newton) 
        return false; //TODO: check individual contacts

    bool anyRestitution = false;
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        const ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        const MultiplierIndex              mx = rt.m_Nk;
        Real& v = verr[mx];
        if (rt.m_effCOR > 0) {
            v *= (1+rt.m_effCOR);
            anyRestitution = true;
        }
    }
    return anyRestitution;
}

bool SemiExplicitEulerTimeStepper::
applyPoissonRestitutionIfAny(const State& s, Vector& impulse,
                             Array_<int>& expanders) const {
    expanders.clear();
    if (getRestitutionModel() == Newton) 
        return false; //TODO: check individual contacts

    bool anyRestitution = false;
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        const ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        const MultiplierIndex              mx = rt.m_Nk;
        Real& pi = impulse[mx];
        if (rt.m_effCOR > 0 && std::abs(pi) >= SignificantReal) {
            pi *= (1+rt.m_effCOR);
            anyRestitution = true;
            expanders.push_back(i);
        }
    }
    return anyRestitution;
}


bool SemiExplicitEulerTimeStepper::
calcExpansionImpulseIfAny(const State& s, const Array_<int>& impacters,
                          const Vector& compressionImpulse,
                          Vector& expansionImpulse,
                          Array_<int>& expanders) const 
{
    expansionImpulse.resize(compressionImpulse.size());
    expansionImpulse.setToZero();
    expanders.clear();
    if (getRestitutionModel() == Newton) 
        return false; //TODO: check individual contacts

    bool anyRestitution = false;
    for (unsigned i=0; i < impacters.size(); ++i) {
        const int which = impacters[i];
        const ImpulseSolver::UniContactRT& rt = m_uniContact[which];
        const MultiplierIndex              mx = rt.m_Nk;
        const Real& pi = compressionImpulse[mx];
        if (rt.m_effCOR > 0 && std::abs(pi) >= SignificantReal) {
            expansionImpulse[mx] = pi*rt.m_effCOR;
            anyRestitution = true;
            expanders.push_back(impacters[i]);
        }
    }
    return anyRestitution;
}

// This phase uses all the proximal constraints and should use a starting
// guess for impulse saved from the last step if possible.
bool SemiExplicitEulerTimeStepper::
doCompressionPhase(const State& s, const Vector& eps, Vector& compImpulse) {
#ifndef NDEBUG
    printf("DYN t=%.15g verr=", s.getTime()); cout << eps << endl;
#endif
    // TODO: improve initial guess
    compImpulse.resize(m_GMInvGt.ncol()); compImpulse.setToZero();
    bool converged = m_pgsSolver.solve(0,
        m_allParticipating,m_GMInvGt,m_D,eps,compImpulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase uses all the proximal constraints, but we expect the result
// to be zero unless expansion causes new violations.
bool SemiExplicitEulerTimeStepper::
doExpansionPhase(const State&, const Vector& eps, Vector& reactionImpulse) {
    // TODO: improve initial guess
    reactionImpulse.resize(m_GMInvGt.ncol()); reactionImpulse.setToZero();
    bool converged = m_pgsSolver.solve(1,
        m_participating,m_GMInvGt,m_D,eps,reactionImpulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase includes only impacting contacts, plus the frictional constraints
// from expanders. It does not include any constraints from observers, nor the
// normal constraint from expanders.
bool SemiExplicitEulerTimeStepper::
doInducedImpactRound(const State& s, const Vector& eps, Vector& impulse)
{
#ifndef NDEBUG
    printf("IMP t=%.15g verr=", s.getTime()); cout << eps << endl;
#endif
    bool converged = m_pgsSolver.solve(0,
        m_participating,m_GMInvGt,m_D,eps,impulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase uses only holonomic constraints, and zero is a good initial
// guess for the (hopefully small) position correction.
bool SemiExplicitEulerTimeStepper::
doPositionCorrectionPhase(const State&, const Vector& eps,
                          Vector& positionImpulse) {
    positionImpulse.resize(m_GMInvGt.ncol()); positionImpulse.setToZero();
    bool converged = m_pgsSolver.solve(2,
        m_posParticipating,m_GMInvGt,m_D,eps,positionImpulse,
        m_posUnconditional,m_posUniContact,m_posNoUniSpeed,m_posNoBounded,
        m_posNoConsLtdFriction, m_posNoStateLtdFriction);
    return converged;
}

bool SemiExplicitEulerTimeStepper::
anyPositionErrorsViolated(const State&, const Vector& perr) const {
    // TODO: no need to fix if large perrs satisfy inequalities.
    bool anyViolated = perr.normInf() > m_consTol;
    SimTK_DEBUG2("maxAbs(perr)=%g -> %s\n", perr.normInf(),
                anyViolated ? "VIOLATED" : "OK");
    return anyViolated;
}



} // namespace SimTK