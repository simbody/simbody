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

using namespace SimTK;

/* Implementation of SemiExplicitEulerTimeStepper. */

namespace {
    using SXTS = SemiExplicitEulerTimeStepper;
    const Real                     DefAccuracy            = 1e-2;
    const Real                     DefConstraintTol       = 1e-3;
    const Real                     DefMinCORVelocity      = 1;
    const Real                     DefMinSignificantForce = SignificantReal;
    const int                      DefMaxInducedImpactsPerStep = 5;
    const SXTS::RestitutionModel   DefRestitutionModel    = SXTS::Poisson;
    const SXTS::InducedImpactModel DefInducedImpactModel  = SXTS::Simultaneous;
    const SXTS::ImpulseSolverType  DefImpulseSolverType   = SXTS::PLUS;
    const SXTS::PositionProjectionMethod DefPosProjMethod = SXTS::Bilateral;
}

namespace SimTK {
//------------------------------------------------------------------------------
//                              CONSTRUCTOR
//------------------------------------------------------------------------------

SemiExplicitEulerTimeStepper::
SemiExplicitEulerTimeStepper(const MultibodySystem& mbs)
:   m_mbs(mbs), 
    m_accuracy(DefAccuracy), m_consTol(DefConstraintTol), 
    m_restitutionModel(DefRestitutionModel),
    m_inducedImpactModel(DefInducedImpactModel),
    m_maxInducedImpactsPerStep(DefMaxInducedImpactsPerStep),
    m_projectionMethod(DefPosProjMethod),
    m_solverType(DefImpulseSolverType), 
    m_defaultMinCORVelocity(DefMinCORVelocity),
    m_defaultCaptureVelocity(m_consTol),
    m_defaultTransitionVelocity(m_consTol),
    m_minSignificantForce(DefMinSignificantForce),
    m_solver(0)
{}


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
    // m is the total number of proximal constraint equation. This won't 
    // change during the step.
    const int m = verr0.size();

    if (m==0) {
        takeUnconstrainedStep(s, h);
        return Integrator::ReachedScheduledEvent;
    }

    // Friction coefficient is fixed by initial slip velocity and doesn't change
    // during impact processing even though the slip velocity will change.
    // The logic is that it takes time for surface asperities to engage or
    // separate and no time is going by during an impact.
    calcCoefficientsOfFriction(s, verr0);

    // The CORs, on the other hand, start out with values based on the 
    // presenting impact velocities, but will change to reflect new velocities
    // produced during induced impacts.
    calcCoefficientsOfRestitution(s, verr0);

    // If we're in Newton mode, or if a contact specifies Newton restitution,
    // then modify the appropriate verr's here.
    Vector verr = verr0;
    const bool anyNewton = applyNewtonRestitutionIfAny(s, verr);

    // Calculate the constraint compliance matrix A=GM\~G.
    matter.calcProjectedMInv(s, m_GMInvGt); // m X m

    // TODO: this is for soft constraints. D >= 0.
    m_D.resize(m); m_D.setToZero();

    Vector totalImpulse(m, Real(0));
    Vector impulse(m);
    m_expansionImpulse.resize(m); 
    m_expansionImpulse.setToZero();

    int numImpactRounds = 0; // how many impact rounds?
    while (true) {
        classifyUnilateralContactsForSequentialImpact
           (verr, m_expansionImpulse, m_uniContact,
            m_impacters, m_expanders, m_observers, 
            m_participating, m_expanding);

        // TODO: unconditionals, unispeed, bounded, cons-ltd, state-ltd friction

        if (!(m_impacters.size() || m_expanders.size())) {
            SimTK_DEBUG1("No impacters or expanders after %d Poisson rounds.\n",
                         numImpactRounds);
            break;
        }

        //--------------------PERFORM ONE IMPACT ROUND--------------------------
        ++numImpactRounds;
        #ifndef NDEBUG
        printf("\nIMPACT ROUND %d:\n", numImpactRounds);
        printf("impacters: "); cout << m_impacters << endl;
        cout << "participating multx: " << m_participating << endl;
        cout << "expanding multx=" << m_expanding << endl;
        cout << "expansion impulse=" << m_expansionImpulse << endl;
        cout << "             verr=" << verr << endl;
        #endif
        doInducedImpactRound(s, m_expanding, m_expansionImpulse, verr, impulse);
        // verr has been updated here to verr-A*(expansion+impulse)
        //----------------------------------------------------------------------
        totalImpulse += impulse;
        if (!m_expanding.empty()) totalImpulse += m_expansionImpulse;

        // Calculate the new expansion impulse.
        const bool anyExpansion = 
            calcExpansionImpulseIfAny(s, m_impacters, impulse,
                                      m_expansionImpulse, m_expanders);
        #ifndef NDEBUG
        cout << "Postimpact verr=" << verr << endl;
        cout << "Postimpact impulse=" << impulse << endl;
        cout << "Next expansion impulse=" << m_expansionImpulse << endl;
        cout << "Expanders: " << m_expanders << endl;
        #endif

        if (numImpactRounds == m_maxInducedImpactsPerStep)
            break; // that's enough!

        calcCoefficientsOfRestitution(s, verr);
    }

    // If we processed an impact, update the state now to reflect the new
    // velocities and then recalculate the applied and rotational forces.

    if (numImpactRounds) {
        Vector genImpulse, deltaU, expVerr; 
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, -totalImpulse, genImpulse);
        // gen impulse to delta-u
        matter.multiplyByMInv(s, genImpulse, deltaU);
        s.updU() += deltaU;
        mbs.realize(s, Stage::Velocity); // updates rotational forces
        #ifndef NDEBUG
        cout << "Impact: imp " << totalImpulse << "-> du " << deltaU << endl;
        cout << "  Now verr=" << s.getUErr() << endl;
        #endif
    } else {
        SimTK_DEBUG("No impact.");
    }

    // Evaluate applied forces and get reference to them. These include
    // gravity but not centrifugal forces.
    mbs.realize(s, Stage::Dynamics);
    const Vector&              f = mbs.getMobilityForces(s, Stage::Dynamics);
    const Vector_<SpatialVec>& F = mbs.getRigidBodyForces(s, Stage::Dynamics);

    // Calculate udotExt = M\(f + ~J*(F-C)) where C are rotational forces.
    // This is the unconstrained acceleration.
    Vector udotExt; Vector_<SpatialVec> A_GB;
    matter.calcAccelerationIgnoringConstraints(s,f,F,udotExt,A_GB);

    // Calculate verrExt = G*deltaU; the end-of-step constraint error due to 
    // external and rotational forces.
    Vector verrExt;
    matter.multiplyByG(s, h*udotExt, verrExt);
    verrExt += verr;

    // Make all unilateral contacts participate for this phase.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        rt.m_type = ImpulseSolver::Participating;
    }

    #ifndef NDEBUG
    printf("\nDYNAMICS PHASE after %d impact rounds:\n", numImpactRounds);
    cout << "   verrExt=" << verrExt << endl;
    #endif
    doCompressionPhase(s, verrExt, impulse);
    #ifndef NDEBUG
    cout << "   dynamics impulse=" << impulse << endl;
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

    // Done with velocity update. Now calculate qdot, possibly including
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

        #ifndef NDEBUG
        printf("\nPOSITION CORRECTION PHASE:\n");
        cout << "   posVerr=" << posVerr << endl;
        #endif

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

    #ifndef NDEBUG
    mbs.realize(s, Stage::Velocity);
    cout << "Before position update, verr=" << s.getUErr() << endl;
    cout << "  perr=" << s.getQErr() << endl;
    #endif

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
    printf("END OF STEP (%g,%g): verr=", t0,s.getTime());
    cout << s.getUErr() << endl;
    cout << "  perr=" << s.getQErr() << endl;
    cout << "  constraint status:\n";
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        printf("  %d: cont %s fric %s\n", rt.m_ucx,
               ImpulseSolver::getUniCondName(rt.m_contactCond),
               ImpulseSolver::getFricCondName(rt.m_frictionCond));
    }
    #endif

    return Integrator::ReachedScheduledEvent;
}

//------------------------------------------------------------------------------
//                     PERFORM SIMULTANEOUS IMPACT
//------------------------------------------------------------------------------
bool SemiExplicitEulerTimeStepper::
performSimultaneousImpact(const Vector& verr, const Vector& expansionImpulse) {
    return false;
}

//------------------------------------------------------------------------------
//                    PERFORM INDUCED IMPACT ROUND
//------------------------------------------------------------------------------
bool SemiExplicitEulerTimeStepper::
performInducedImpactRound(const Vector& verr, const Vector& expansionImpulse) {
    classifyUnilateralContactsForSequentialImpact
       (verr, expansionImpulse, m_uniContact,
        m_impacters, m_expanders, m_observers, 
        m_participating, m_expanding);
    return false;
}

//------------------------------------------------------------------------------
//         CLASSIFY UNILATERAL CONTACTS FOR SEQUENTIAL IMPACT
//------------------------------------------------------------------------------
// Calculate participating constraints. Include all proximals except:
//   - ignore "observers" (unilateral contacts with nonnegative impact 
//     velocity)
//   - ignore normal constraints for "expanders" (friction stays)
// Then impacters=all contacts with negative impact velocity
void SemiExplicitEulerTimeStepper::
classifyUnilateralContactsForSequentialImpact
   (const Vector&                           verr,
    const Vector&                           expansionImpulse,
    Array_<ImpulseSolver::UniContactRT>&    uniContacts, 
    Array_<int>&                            impacters,
    Array_<int>&                            expanders,
    Array_<int>&                            observers,
    Array_<MultiplierIndex>&                participaters,
    Array_<MultiplierIndex>&                expanding) const
{
    assert(verr.size());
    assert(verr.size() == expansionImpulse.size());
    impacters.clear(); expanders.clear(); observers.clear(); 
    participaters.clear(); expanding.clear();

    for (unsigned i=0; i < uniContacts.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = uniContacts[i];
        const MultiplierIndex mz = rt.m_Nk;
        if (expansionImpulse[mz] != 0) {
            rt.m_type = ImpulseSolver::Known;
            if (rt.hasFriction()) // friction still participates
                for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                    participaters.push_back(rt.m_Fk[j]);
            expanders.push_back(i);  // uni contact index
            expanding.push_back(mz); // multiplier index
            continue; // Expander
        }

        // Observe even if we're slightly (but insignificantly) negative.
        // We're using 1/10 of constraint tol as the cutoff.
        // Don't use zero, or you'll loop forever due to roundoff!
        if (verr[mz] >= (Real(-0.1))*m_consTol) {
            rt.m_type = ImpulseSolver::Observing;
            observers.push_back(i);
            continue; // Observer
        }

        // Impacter (Participator)
        rt.m_type = ImpulseSolver::Participating;
        participaters.push_back(mz);
        if (rt.hasFriction())
            for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                participaters.push_back(rt.m_Fk[j]);
        impacters.push_back(i);
    }

    assert(   impacters.size()+expanders.size()+observers.size()
           == uniContacts.size());
}

//------------------------------------------------------------------------------
//                             INITIALIZE
//------------------------------------------------------------------------------
void SemiExplicitEulerTimeStepper::initialize(const State& initState) {
    m_state = initState;
    m_mbs.realize(m_state, Stage::Acceleration);

    if (!m_solver) {
        m_solver = m_solverType==PLUS 
            ? (ImpulseSolver*)new PLUSImpulseSolver(m_defaultTransitionVelocity)
            : (ImpulseSolver*)new PGSImpulseSolver(m_defaultTransitionVelocity);
    }

    SimTK_ERRCHK_ALWAYS(m_solver!=0,
                        "SemiExplicitEulerTimeStepper::initialize()",
                        "No ImpulseSolver available.");

    // Make sure the impulse solve knows our tolerance for slip velocity
    // during rolling.
    m_solver->setMaxRollingSpeed(getDefaultFrictionTransitionVelocityInUse());
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
        if (fric.getNormalForceMagnitude(s) >= m_minSignificantForce) 
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
    // This is invalid now but we can't fill it out here since the unilateral
    // contact part changes for each Poisson impact round.
    // TODO: should split up into fixed and changing parts.
    m_participating.clear();  

    // This will be filled out here since it can't change later.
    m_allParticipating.clear();

    // These are the proximal constraints for contact & impact.
    m_unconditional.clear();
    m_uniContact.clear();
    m_uniSpeed.clear();
    m_bounded.clear();
    m_consLtdFriction.clear();
    m_stateLtdFriction.clear();
    
    // These include only the proximal holonomic constraints. 
    m_posUnconditional.clear();
    m_posUniContact.clear();
    m_posParticipating.clear();

    const int nConstraints = matter.getNumConstraints();
    for (ConstraintIndex cx(0); cx < nConstraints; ++cx) {
        // TODO: sort out the enabled, unconditional ones and fill in
        // m_unconditional and m_posUnconditional.
        // m_allParticipating.push_back(...);
        // if holnomic: m_posParticipating.push_back(...);
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
        posRt.m_type = ImpulseSolver::Participating;
        posRt.m_ucx = cx;
        posRt.m_Nk = rt.m_Nk;
        m_posParticipating.push_back(rt.m_Nk); // normal is holonomic
    }

    //TODO: speed, bounded, consLtd, proximal stateLtd friction 
    // (all nonholonomic)
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
            printf("IMPACT cuz verr[%d]=%g\n", (int)rt.m_Nk, verr[rt.m_Nk]);
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
        Real slipSpeed = 0;
        for (unsigned j=0; j < rt.m_Fk.size(); ++j)
            slipSpeed += square(verr[rt.m_Fk[j]]);
        slipSpeed = std::sqrt(slipSpeed);
        const UnilateralContact& uni = matter.getUnilateralContact(rt.m_ucx);
        rt.m_effMu = uni.calcEffectiveCOF(s, 
                        getDefaultFrictionTransitionVelocity(),
                        slipSpeed);
        SimTK_DEBUG3("Uni fric %d verr speed=%g -> mu=%g\n", (int)rt.m_ucx, 
                     slipSpeed, rt.m_effMu);
    }

    // These are the proximal state-limited friction elements.
    for (unsigned i=0; i < m_stateLtdFriction.size(); ++i) {
        ImpulseSolver::StateLtdFrictionRT& rt = m_stateLtdFriction[i];
        Real slipSpeed = 0;
        for (unsigned j=0; j < rt.m_Fk.size(); ++j)
            slipSpeed += square(verr[rt.m_Fk[j]]);
        slipSpeed = std::sqrt(slipSpeed);
        const StateLimitedFriction& fric = 
            matter.getStateLimitedFriction(rt.m_slfx);
        rt.m_effMu = fric.calcEffectiveCOF(s, 
                        getDefaultFrictionTransitionVelocity(),
                        slipSpeed);
        SimTK_DEBUG3("State ltd fric %d speed=%g -> mu=%g\n", (int)rt.m_slfx, 
                     slipSpeed, rt.m_effMu);
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
        Real impactVel = verr[rt.m_Nk];
        if (impactVel > 0) { //TODO: sign
            rt.m_effCOR = 0; // separating
            continue;
        }
        if (impactVel >= -m_consTol) {
            rt.m_effCOR = 0;
            SimTK_DEBUG3(
                "Uni contact %d verr vel=%g at or below tol %g -> cor=0\n", 
                 (int)rt.m_ucx, impactVel, m_consTol);
            continue;
        }
        const UnilateralContact& uni = matter.getUnilateralContact(rt.m_ucx);
        rt.m_effCOR = uni.calcEffectiveCOR(s, 
                        getDefaultImpactCaptureVelocityInUse(),
                        getDefaultImpactMinCORVelocityInUse(),
                        -impactVel);
        SimTK_DEBUG3("Uni contact %d verr vel=%g -> cor=%g\n", (int)rt.m_ucx, 
                     impactVel, rt.m_effCOR);
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
        const MultiplierIndex              mz = rt.m_Nk;
        const Real& pi = compressionImpulse[mz];
        if (rt.m_effCOR > 0 && std::abs(pi) >= SignificantReal) {
            expansionImpulse[mz] = pi*rt.m_effCOR;
            anyRestitution = true;
            expanders.push_back(impacters[i]);
        }
    }
    return anyRestitution;
}

// This phase uses all the proximal constraints and should use a starting
// guess for impulse saved from the last step if possible.
bool SemiExplicitEulerTimeStepper::
doCompressionPhase(const State& s, Vector& verr, Vector& compImpulse) {
#ifndef NDEBUG
    printf("DYN t=%.15g verr=", s.getTime()); cout << verr << endl;
#endif
    // TODO: improve initial guess
    bool converged = m_solver->solve(0,
        m_allParticipating,m_GMInvGt,m_D,
        Array_<MultiplierIndex>(), verr,//dummy for piExpand(m); ignored
        verr,compImpulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase uses all the proximal constraints, but we expect the result
// to be zero unless expansion causes new violations.
bool SemiExplicitEulerTimeStepper::
doExpansionPhase(const State&   state,
                 const Array_<MultiplierIndex>& expanding,
                 Vector&        expansionImpulse, 
                 Vector&        verr, 
                 Vector&        reactionImpulse) {
    // TODO: improve initial guess
    bool converged = m_solver->solve(1,
        m_participating,m_GMInvGt,m_D,
        expanding,expansionImpulse,
        verr,reactionImpulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase includes only impacting contacts, plus the frictional constraints
// from expanders. It does not include any constraints from observers, nor the
// normal constraint from expanders.
bool SemiExplicitEulerTimeStepper::
doInducedImpactRound(const State& s, 
                     const Array_<MultiplierIndex>& expanding,
                     Vector& expansionImpulse, 
                     Vector& verr, Vector& impulse)
{
#ifndef NDEBUG
    printf("IMP t=%.15g verr=", s.getTime()); cout << verr << endl;
#endif
    bool converged = m_solver->solve(0,
        m_participating,m_GMInvGt,m_D,
        expanding,expansionImpulse,
        verr,impulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
    return converged;
}
// This phase uses only holonomic constraints, and zero is a good initial
// guess for the (hopefully small) position correction.
bool SemiExplicitEulerTimeStepper::
doPositionCorrectionPhase(const State& state, Vector& verr,
                          Vector& positionImpulse) {
    bool converged = m_solver->solve(2,
        m_posParticipating,m_GMInvGt,m_D,
        Array_<MultiplierIndex>(), verr,//dummy for piExpand; ignored
        verr, positionImpulse,
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