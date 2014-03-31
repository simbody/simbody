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

#include "SimbodyMatterSubsystemRep.h"

#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

/* Implementation of SemiExplicitEulerTimeStepper. */

namespace {
    const Real  DefAccuracy            = 1e-2;
    const Real  DefConstraintTol       = DefAccuracy/10;
    const Real  DefMinSignificantForce = SignificantReal;
    const int   DefMaxInducedImpactsPerStep = 5;
    const SemiExplicitEulerTimeStepper::RestitutionModel   
        DefRestitutionModel    = SemiExplicitEulerTimeStepper::Poisson;
    const SemiExplicitEulerTimeStepper::InducedImpactModel 
        DefInducedImpactModel  = SemiExplicitEulerTimeStepper::Simultaneous;
    const SemiExplicitEulerTimeStepper::ImpulseSolverType  
        DefImpulseSolverType   = SemiExplicitEulerTimeStepper::PLUS;
    const SemiExplicitEulerTimeStepper::PositionProjectionMethod 
        DefPosProjMethod = SemiExplicitEulerTimeStepper::Bilateral;
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
    m_defaultCaptureVelocity(0),    // means: use 2 x constraintTol
    m_defaultMinCORVelocity(0),     // means: use capture velocity
    m_defaultTransitionVelocity(0), // means: use 2 x constraintTol
    m_minSignificantForce(DefMinSignificantForce),
    m_solver(0)
{}


//------------------------------------------------------------------------------
//                                 STEP TO
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

    // Note that verr0 is a reference into State s so if we update the
    // velocities in s and realize them, verr0 will also be updated.
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

    // Calculate the constraint compliance matrix A=GM\~G.
    matter.calcProjectedMInv(s, m_GMInvGt); // m X m

    // TODO: this is for soft constraints. D >= 0.
    m_D.resize(m); m_D.setToZero();
    m_totalImpulse.resize(m); m_totalImpulse.setToZero();
    m_expansionImpulse.resize(m); m_expansionImpulse.setToZero();
    m_impulse.resize(m);

    m_verr = verr0;

    int numImpactRounds = 0;
    if (getInducedImpactModel() == Simultaneous) {
      numImpactRounds = performSimultaneousImpact(s, m_verr, m_totalImpulse);
    } else {
      bool disableRestitution = (getRestitutionModel()==NoRestitution);
      while (true) {
        // CORs start out with values based on the presenting impact velocities,
        // but change for each round to reflect new velocities produced during
        // induced impacts.
        calcCoefficientsOfRestitution(s, m_verr, disableRestitution);

        // If we're in Newton mode, calculate the restitution verr.
        const bool anyNewton = (getRestitutionModel()==Newton) &&
            calcNewtonRestitutionIfAny(s, m_verr, m_newtonRestitutionVerr);

        const bool isLastRound = 
            numImpactRounds >= m_maxInducedImpactsPerStep-1;

        const bool includeAllProximals =
            (getInducedImpactModel()!=Sequential) && isLastRound;

        classifyUnilateralContactsForSequentialImpact
           (m_verr, m_expansionImpulse, includeAllProximals, false,
            m_uniContact, m_impacters, m_expanders, m_observers, 
            m_participating, m_expanding);

        // TODO: unconditionals, unispeed, bounded, cons-ltd, state-ltd friction

        if (!(m_impacters.size() || m_expanders.size())) {
            SimTK_DEBUG1("No impacters or expanders after %d impact rounds.\n",
                         numImpactRounds);
            break;
        }

        //--------------------PERFORM ONE IMPACT ROUND--------------------------
        ++numImpactRounds;
        if (anyNewton)
            m_verr += m_newtonRestitutionVerr;
        #ifndef NDEBUG
        printf("\nIMPACT ROUND %d:\n", numImpactRounds);
        printf( "  impacters: "); cout << m_impacters << endl;
        cout << "  participating multx: " << m_participating << endl;
        cout << "  expanding multx=" << m_expanding << endl;
        cout << "  expansion impulse=" << m_expansionImpulse << endl;
        cout << "  verr=" << m_verr << endl;
        if (anyNewton)
            cout << "    incl. Newton verr=" << m_newtonRestitutionVerr << endl;
        #endif
        doInducedImpactRound(s,m_expanding,m_expansionImpulse,m_verr,m_impulse);
        // verr has been updated here to verr-A*(expansion+impulse)
        if (anyNewton) // remove any fake verrs
            m_verr -= m_newtonRestitutionVerr;
        m_totalImpulse += m_impulse;
        if (!m_expanding.empty()) m_totalImpulse += m_expansionImpulse;
        //----------------------------------------------------------------------

        // Calculate the next expansion impulse from the compression impulses
        // we just generated.
        const bool anyExpansion = 
            calcExpansionImpulseIfAny(s, m_impacters, m_impulse,
                                      m_expansionImpulse, m_expanders);
        #ifndef NDEBUG
        cout << "Postimpact verr=" << m_verr << endl;
        cout << "Postimpact impulse=" << m_impulse << endl;
        cout << "Next expansion impulse=" << m_expansionImpulse << endl;
        cout << "Next expanders: " << m_expanders << endl;
        #endif

        if (isLastRound) {
            SimTK_DEBUG1("\nTOO MANY IMPACT ROUNDS (%d); bailing.\n", 
                         numImpactRounds);
                         
            if (anyExpansion) {
                #ifndef NDEBUG
                cout << "  final round has leftover expansion impulse: " 
                     << m_expansionImpulse << endl;
                #endif
                if (getInducedImpactModel() == Sequential) {
                    SimTK_DEBUG("  Sequential mode; allow final penetration.\n");
                    // Do an expansion-only round that generates no reaction 
                    // except friction forces at the expanding contacts.
                    classifyUnilateralContactsForSequentialImpact
                        (m_verr, m_expansionImpulse, false, true, m_uniContact,
                        m_impacters, m_expanders, m_observers, 
                        m_participating, m_expanding);
                } else { // Mixed mode
                    // Do final expansion/compression round with everything 
                    // participating and all COR effectively 0.
                    SimTK_DEBUG1("  %s mode: do final COR=0 expansion round.\n",
                        getInducedImpactModelName(getInducedImpactModel()));
                    classifyUnilateralContactsForSequentialImpact
                       (m_verr, m_expansionImpulse, true, false,
                        m_uniContact, m_impacters, m_expanders, m_observers, 
                        m_participating, m_expanding);
                }
                //----------------PERFORM FINAL ROUND---------------
                ++numImpactRounds;
                #ifndef NDEBUG
                printf("FINAL EXPANSION ROUND:\n");
                cout << "   verr=" << m_verr << endl;
                #endif
                doInducedImpactRound(s, m_expanding, m_expansionImpulse, 
                                     m_verr, m_impulse);
                m_totalImpulse += m_impulse;
                m_totalImpulse += m_expansionImpulse;
                //--------------------------------------------------------------
                #ifndef NDEBUG
                cout << "  after final expansion, verr=" << m_verr << endl;
                #endif
            }
            break; // that's enough!
        }
      }
    }

    // If we processed an impact, update the state now to reflect the new
    // velocities and then recalculate the applied and rotational forces.

    if (numImpactRounds) {
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, -m_totalImpulse, m_genImpulse);
        // gen impulse to delta-u
        matter.multiplyByMInv(s, m_genImpulse, m_deltaU);
        s.updU() += m_deltaU;
        // This updates velocity-dependent forces including the rotational
        // (Coriolis & gyroscopic) forces. This also updates verr0 which is
        // a reference into the State s.
        mbs.realize(s, Stage::Velocity);
        #ifndef NDEBUG
        cout << "Impact: imp " << m_totalImpulse << "-> du " << m_deltaU << endl;
        cout << "  Now verr0=" << verr0 << endl;
        #endif
    } else {
        SimTK_DEBUG("No impact.\n");
    }

    // Evaluate applied forces and get reference to them. These include
    // gravity but not centrifugal forces.
    mbs.realize(s, Stage::Dynamics);
    const Vector&              f  = mbs.getMobilityForces(s, Stage::Dynamics);
    const Vector_<Vec3>&       Fp = mbs.getParticleForces(s, Stage::Dynamics);
    const Vector_<SpatialVec>& F  = mbs.getRigidBodyForces(s, Stage::Dynamics);

    // Get the acceleration-level cache entries and begin filling them in.
    // At end of step we'll make these available but they will represent
    // accelerations and reaction forces with updated u and old q.
    const SimbodyMatterSubsystemRep& matterRep = matter.getRep();
    const SBModelCache&              mc = matterRep.getModelCache(s);
    SBTreeAccelerationCache&         tac = 
        matterRep.updTreeAccelerationCache(s);
    SBConstrainedAccelerationCache&  cac = 
        matterRep.updConstrainedAccelerationCache(s);
    Vector&                 qdot    = s.updQDot();
    Vector&                 zdot    = s.updZDot();
    Vector&                 udot    = s.updUDot();
    Vector&                 qdotdot = s.updQDotDot();
    Vector&                 udotErr = s.updUDotErr();
    Vector&                 lambda  = s.updMultipliers();

    // Calculate udot = M\(f + ~J*(F-C)) where C are rotational forces.
    // This is the unconstrained acceleration; we'll be overwriting udot and
    // A_GB with their final values below.
    matterRep.calcTreeForwardDynamicsOperator
       (s, f, Fp, F, 0, 0, tac, udot, qdotdot, udotErr);
    m_deltaU = h*udot;

    // Calculate verr = G*deltaU; the end-of-step constraint error due to 
    // external and rotational forces.
    matter.multiplyByG(s, m_deltaU, m_verr);

    // Make all proximal unilateral contacts participate for this phase.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        rt.m_type = ImpulseSolver::Participating;
    }

    #ifndef NDEBUG
    printf("\nDYNAMICS PHASE after %d impact rounds:\n", numImpactRounds);
    cout << "   verr0  =" << verr0 << endl;
    cout << "   verrExt=" << m_verr << endl;
    cout << "   total verr=" << verr0+m_verr << endl;
    #endif

    // Add in verr0=Gu-b because we want G (u+du) - b = 0 when we're done
    // (except for sliding friction). The impulse solver needs to know the
    // total sliding velocity for proper friction classification, that
    // will be in verr=verr0+hGM\f.
    m_verr += verr0;

    // Use lambda as a temp here; we are really calculating lambda*h.
    doCompressionPhase(s, m_verr, lambda);
    #ifndef NDEBUG
    cout << "   dynamics impulse=" << lambda << endl;
    cout << "   updated verr=" << m_verr << endl;
    #endif
    m_totalImpulse += lambda;

    // Convert multipliers from impulses to forces. These are the multipliers
    // reported at end of step.
    lambda /= h;

    // Calculate constraint forces ~G*lambda (body frcs Fc, mobility frcs fc).
    Vector_<SpatialVec> Fc; Vector fc; 
    matterRep.calcConstraintForcesFromMultipliers(s,lambda,Fc,fc,
        cac.constrainedBodyForcesInG, cac.constraintMobilityForces);

    // Now calculate final udot = M\(f-fc + ~J*(F-Fc-C)) and corresponding
    // body accelerations A_GB.
    matterRep.calcTreeForwardDynamicsOperator
       (s, f, Fp, F, &fc, &Fc, tac, udot, qdotdot, udotErr);

    // Update auxiliary states z, invalidating Stage::Dynamics.
    s.updZ() += h*zdot;

    // Update u from deltaU, invalidating Stage::Velocity. 
    s.updU() += m_deltaU;

    // Done with velocity update. Now calculate qdot, possibly including
    // an additional position error correction term.
    // TODO: Do one linear position correction iteration unconditionally.
    // Over several steps this should perform the nonlinear correction needed
    // to get nice positions.
    const Vector& perr0 = s.getQErr();
    if (   m_projectionMethod==NoPositionProjection
        /*|| !anyPositionErrorsViolated(s, perr0)*/) //TODO: unconditional now
    {
        matter.multiplyByN(s,false,s.getU(),qdot);
    } else {
        // Perform position projection.
        // Don't include quaternions for position correction. 
        const int nQuat = matter.getNumQuaternionsInUse(s);
        m_verr.setToZero();
        m_verr(0, perr0.size()-nQuat) = perr0(0, perr0.size()-nQuat)/h;

        #ifndef NDEBUG
        printf("\nPOSITION CORRECTION PHASE:\n");
        cout << "   posVerr=" << m_verr << endl;
        #endif

        //----------------------------------------------------------------------
        // Calculate impulse and then deltaU=M\~G*impulse such that -h*deltaU 
        // will eliminate position errors, respecting only position constraints.
        doPositionCorrectionPhase(s, m_verr, m_impulse);
        //----------------------------------------------------------------------
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, m_impulse, m_genImpulse);
        // gen impulse to deltaU (watch sign)
        matter.multiplyByMInv(s, m_genImpulse, m_deltaU);

        // convert corrected u to qdot (note we're not changing u)
        matter.multiplyByN(s,false,s.getU()-m_deltaU, qdot);
    }

    #ifndef NDEBUG
    mbs.realize(s, Stage::Velocity);
    cout << "Before position update, verr=" << s.getUErr() << endl;
    cout << "  perr=" << s.getQErr() << endl;
    #endif

    // We have qdot, now update q, fix quaternions, update time.
    s.updQ() += h*qdot; // invalidates Stage::Position
    s.updTime() += h;   // invalidates Stage::Time
    mbs.realize(s, Stage::Time);
    matter.normalizeQuaternions(s);

    // Return from step with kinematics realized. Note that we may have
    // broken the velocity constraints by updating q, but we won't fix that
    // until the next step. Also position constraints are only imperfectly
    // satisfied by the correction above.
    // TODO: have to recalculate forces here just to get past Dynamics stage
    // so that reaction forces can be obtained.
    mbs.realize(s, Stage::Dynamics);

    // Mark the acceleration-level calculations valid now, despite the fact
    // that they don't reflect the latest time, q, u, and forces.
    matterRep.markCacheValueRealized(s, mc.constrainedAccelerationCacheIndex);
    matterRep.markCacheValueRealized(s, mc.treeAccelerationCacheIndex);


    #ifndef NDEBUG
    printf("END OF STEP (%g,%g): verr=", t0,s.getTime());
    cout << s.getUErr() << endl;
    cout << "  perr=" << s.getQErr() << endl;
    cout << "  constraint status:\n";
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        printf("  %d: cont %s fric %s\n", (int)rt.m_ucx,
               ImpulseSolver::getUniCondName(rt.m_contactCond),
               ImpulseSolver::getFricCondName(rt.m_frictionCond));
    }
    #endif

    return Integrator::ReachedScheduledEvent;
}

//------------------------------------------------------------------------------
//                     PERFORM SIMULTANEOUS IMPACT
//------------------------------------------------------------------------------
int SemiExplicitEulerTimeStepper::
performSimultaneousImpact(const State&  s, 
                          Vector&       verr, // in/out
                          Vector&       totalImpulse) 
{
    static const Vector noExpansion;
    const int m = verr.size();
    const bool disableRestitution = (getRestitutionModel()==NoRestitution);
    int numImpactRounds = 0; // how many impact rounds?

    classifyUnilateralContactsForSimultaneousImpact
       (verr, noExpansion, m_uniContact,
        m_impacters, m_expanders, m_observers, 
        m_participating, m_expanding);  
    assert(m_expanders.empty() && m_expanding.empty());

    // If we didn't see any impacters then we don't have to do anything.
    if (!m_impacters.size()) { // there can't be any expanders here
        SimTK_DEBUG(
            "performSimultaneousImpact(): No impacters; nothing to do.\n");
        return numImpactRounds; // zero
    }   

    // For simultaneous impact, CORs are determined only for the initial
    // impact velocities. In the reaction round they are treated as zero.
    calcCoefficientsOfRestitution(s, verr, disableRestitution);

    // If we're in Newton mode, calculate the restitution verr.
    const bool anyNewton = (getRestitutionModel()==Newton) &&
        calcNewtonRestitutionIfAny(s, verr, m_newtonRestitutionVerr);

    //----------------------PERFORM COMPRESSION ROUND---------------------------
    ++numImpactRounds;
    if (anyNewton)
        verr += m_newtonRestitutionVerr;
    #ifndef NDEBUG
    printf("\nperformSimultaneousImpact(): COMPRESSION ROUND:\n");
    printf( "  impacters: "); cout << m_impacters << endl;
    printf( "  observers: "); cout << m_observers << endl;
    cout << "  participating multx: " << m_participating << endl;
    cout << "  verr=" << verr << endl;
    if (anyNewton)
        cout << "    incl. Newton verr=" << m_newtonRestitutionVerr << endl;
    #endif
    // There are no expansion impulses yet; these are dummies.
    doInducedImpactRound(s, m_expanding, m_expansionImpulse, verr, 
                         totalImpulse); // the only impulse so far
    // verr has been updated here to verr-A*totalImpulse
    if (anyNewton) // remove any fake verrs
        verr -= m_newtonRestitutionVerr;
    //--------------------------------------------------------------------------

    // Calculate the next expansion impulse from the compression impulses
    // we just generated. This only occurs if we are doing Poisson restitution.
    const bool anyExpansion = 
        calcExpansionImpulseIfAny(s, m_impacters, totalImpulse,
                                    m_expansionImpulse, m_expanders);
    #ifndef NDEBUG
    cout << "  Postcompression verr=" << verr << endl;
    cout << "    impulse=" << totalImpulse << endl;
    cout << "  Next expansion impulse=" << m_expansionImpulse << endl;
    cout << "  Next expanders: " << m_expanders << endl;
    #endif  

    if (!anyExpansion) {
        SimTK_DEBUG("performSimultaneousImpact(): No expansion; done.\n");
        return numImpactRounds; // one
    }

    classifyUnilateralContactsForSimultaneousImpact
       (verr, m_expansionImpulse, m_uniContact,
        m_impacters, m_expanders, m_observers, 
        m_participating, m_expanding); 

    //-----------------------PERFORM EXPANSION ROUND----------------------------
    ++numImpactRounds;
    Vector reactionImpulse(m);
    #ifndef NDEBUG
    printf("\nperformSimultaneousImpact(): EXPANSION ROUND:\n");
    printf( "  impacters: "); cout << m_impacters << endl;
    printf( "  observers: "); cout << m_observers << endl;
    printf( "  expanders: "); cout << m_expanders << endl;
    cout << "  participating multx: " << m_participating << endl;
    cout << "  expanding multx: " << m_expanding << endl;
    cout << "  verr=" << verr << endl;
    #endif
    doInducedImpactRound(s, m_expanding, m_expansionImpulse, verr, 
                         reactionImpulse);
    // verr has been updated here to verr-A*(expansionImpulse+reactionImpulse)
    #ifndef NDEBUG
    printf( "  expansion impulse: "); cout << m_expansionImpulse << endl;
    printf( "  reaction impulse: "); cout << reactionImpulse << endl;
    #endif    
    totalImpulse += m_expansionImpulse;
    totalImpulse += reactionImpulse;
    //--------------------------------------------------------------------------

    // That should be the end but it is possible that one of the expanders
    // is now impacting.
    classifyUnilateralContactsForSimultaneousImpact
       (verr, noExpansion, m_uniContact,
        m_impacters, m_expanders, m_observers, 
        m_participating, m_expanding);  
    assert(m_expanders.empty() && m_expanding.empty());

    if (m_impacters.empty()) {
        SimTK_DEBUG("performSimultaneousImpact(): Expansion successful; done.\n");
        return numImpactRounds; // two
    }

    //-----------------------PERFORM CORRECTION ROUND----------------------------
    ++numImpactRounds;
    #ifndef NDEBUG
    printf("\nperformSimultaneousImpact(): CORRECTION ROUND:\n");
    printf( "  impacters: "); cout << m_impacters << endl;
    printf( "  observers: "); cout << m_observers << endl;
    cout << "  participating multx: " << m_participating << endl;
    cout << "  verr=" << verr << endl;
    #endif
    // No expansion here.
    doInducedImpactRound(s, m_expanding, m_expansionImpulse, verr, 
                         reactionImpulse);
    // verr has been updated here to verr-A*reactionImpulse
    #ifndef NDEBUG
    printf( "  correction impulse: "); cout << reactionImpulse << endl;
    printf( "  post-correction verr: "); cout << verr << endl;
    #endif
    totalImpulse += reactionImpulse;
    //--------------------------------------------------------------------------

    return numImpactRounds;
}

//------------------------------------------------------------------------------
//                    PERFORM INDUCED IMPACT ROUND
//------------------------------------------------------------------------------
bool SemiExplicitEulerTimeStepper::
performInducedImpactRound(const Vector& verr, const Vector& expansionImpulse) {
    classifyUnilateralContactsForSequentialImpact
       (verr, expansionImpulse, false, false, m_uniContact,
        m_impacters, m_expanders, m_observers, 
        m_participating, m_expanding);
    return false;
}

//------------------------------------------------------------------------------
//             CLASSIFY UNILATERAL CONTACTS FOR SEQUENTIAL IMPACT
//------------------------------------------------------------------------------
// Calculate participating constraints. Include all proximals except:
//   - ignore "observers" (unilateral contacts with nonnegative impact 
//     velocity)
//   - ignore normal constraints for "expanders" (friction stays)
// Then "impacters" are all contacts with negative impact velocity.
// Options: 
//   - includeAllProximals: observers also are participating; that way they
//       can never be violated when the impact or expansion impulses are
//       applied.
//   - expansion only: both observers and impacters are ignored; the only
//       participating constraints are the friction constraints at expanding
//       contacts.
//
// If includeAllProximals, then everyone participates meaning there won't be
// any penetrating contacts after the next set of impulses is generated.
// If expansionOnly, then everyone is an observer if they aren't expanding;
// that is used to apply a final expansion at the end of the step, with no
// reaction other than the expanding contacts' friction forces.
void SemiExplicitEulerTimeStepper::
classifyUnilateralContactsForSequentialImpact
   (const Vector&                           verr,
    const Vector&                           expansionImpulse,
    bool                                    includeAllProximals,
    bool                                    expansionOnly,
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
            rt.m_type = ImpulseSolver::Known; // normal does not participate
            if (rt.hasFriction()) // but friction still participates
                for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                    participaters.push_back(rt.m_Fk[j]);
            expanders.push_back(i);  // uni contact index
            expanding.push_back(mz); // multiplier index
            continue; // Expander
        }

        // Observe even if we're slightly (but insignificantly) negative.
        // We're using 1/10 of constraint tol as the cutoff.
        // Don't use zero, or you'll loop forever due to roundoff!
        // Also, if we're only doing expanders then everyone else is observing.
        if (expansionOnly || rt.m_sign*verr[mz] >= (Real(-0.1))*m_consTol) {
            if (includeAllProximals) {
                rt.m_type = ImpulseSolver::Participating;
                if (rt.hasFriction()) // but friction still participates
                    for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                        participaters.push_back(rt.m_Fk[j]);
                participaters.push_back(mz);
            } else 
                rt.m_type = ImpulseSolver::Observing;
            observers.push_back(i);
            continue; // Observer
        }

        // (We don't get here if expansionOnly flag is set.)
        // The normal velocity is substantially negative, so this is an
        // impacter and always participates.
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
//             CLASSIFY UNILATERAL CONTACTS FOR SEQUENTIAL IMPACT
//------------------------------------------------------------------------------
// Classify the unilateral contacts into impacters, observers, and expanders.
// and determine which constraint equations should participate assuming we
// are in Sequential ("one shot") impact mode. All impacter and observer 
// constraint equations should participate, as well as the friction constraints
// (but not the normal constraints) for expanders.
// If you don't want to process expanders, pass the expansionImpulse as a
// zero-length Vector, otherwise it must have length m.
void SemiExplicitEulerTimeStepper::
classifyUnilateralContactsForSimultaneousImpact
   (const Vector&                           verr,
    const Vector&                           expansionImpulse,
    Array_<ImpulseSolver::UniContactRT>&    uniContacts, 
    Array_<int>&                            impacters,
    Array_<int>&                            expanders,
    Array_<int>&                            observers,
    Array_<MultiplierIndex>&                participaters,
    Array_<MultiplierIndex>&                expanding) const
{
    const bool hasExpansion = expansionImpulse.size() > 0;
    const int m = verr.size(); assert(m); 
    assert(!hasExpansion || expansionImpulse.size() == m);

    impacters.clear(); expanders.clear(); observers.clear(); 
    participaters.clear(); expanding.clear();

    for (unsigned i=0; i < uniContacts.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = uniContacts[i];
        const MultiplierIndex mz = rt.m_Nk;
        if (hasExpansion && expansionImpulse[mz] != 0) {
            rt.m_type = ImpulseSolver::Known; // normal does not participate
            expanders.push_back(i);  // uni contact index
            expanding.push_back(mz); // multiplier index
        } else if (rt.m_sign*verr[mz] >= (Real(-0.1))*m_consTol) {
            // Observe even if we're slightly (but insignificantly) negative.
            // We're using 1/10 of constraint tol as the cutoff.
            // Don't use zero, or you'll loop forever due to roundoff!
            rt.m_type = ImpulseSolver::Participating;
            participaters.push_back(mz);
            observers.push_back(i);
        } else {
            // The normal velocity is substantially negative, so this is an
            // impacter.
            rt.m_type = ImpulseSolver::Participating;
            participaters.push_back(mz);
            impacters.push_back(i);
        }

        // Regardless of classification, friction constraints participate here.
        if (rt.hasFriction())
            for (unsigned j=0; j < rt.m_Fk.size(); ++j)
                participaters.push_back(rt.m_Fk[j]);
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
        const Real transVel = getDefaultFrictionTransitionVelocityInUse();
        m_solver = m_solverType==PLUS 
            ? (ImpulseSolver*)new PLUSImpulseSolver(transVel)
            : (ImpulseSolver*)new PGSImpulseSolver(transVel);
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
// forces during a step; distal constraints have no effect whatsoever. Simbody 
// will be asked to generate constraint equations for all the proximal
// constraints, but not for the distal constraints. In general only a subset of
// the proximal unilateral constraints will actually generate forces, but all 
// must be considered.
//
// Velocity-level (nonholonomic) unilateral constraints are always considered
// proximal, as are unconditional constraints and their associated friction
// constraints.
//
// Algorithm:
// (1) Inspect each unilateral contact constraint in the MultibodySystem, 
// marking it proximal if its perr() value is within tolerance of violating 
// its defining unilateral inequality, distal otherwise. Typically, perr() is 
// a signed distance function and the defining inequality is perr() >= 0, so 
// the proximal condition is perr() <= tol (that includes all negative values
// of perr). Associated friction constraints, if any, are marked proximal or 
// distal to match the contact constraint.
//
// (2) Inspect each StateLimitedFriction constraint. Mark the friction element
// proximal if its known limiting force is non-zero (or above a small 
// threshold), otherwise it is distal.
//
// TODO: ConstraintLimitedFriction constraints (i.e., friction that acts at
// bilateral constraints) should only be proximal if the limiting constraint is
// enabled.
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


//------------------------------------------------------------------------------
//                        ENABLE PROXIMAL CONSTRAINTS
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//                         COLLECT CONSTRAINT INFO
//------------------------------------------------------------------------------
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
        rt.m_sign = (Real)contact.getSignConvention();
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
        posRt.m_sign = (Real)contact.getSignConvention();
        posRt.m_ucx = cx;
        posRt.m_Nk = rt.m_Nk;
        m_posParticipating.push_back(rt.m_Nk); // normal is holonomic
    }

    //TODO: speed, bounded, consLtd, proximal stateLtd friction 
    // (all nonholonomic)
}

//------------------------------------------------------------------------------
//                        TAKE UNCONSTRAINED STEP
//------------------------------------------------------------------------------
void SemiExplicitEulerTimeStepper::
takeUnconstrainedStep(State& s, Real h) {
    const SimbodyMatterSubsystem&    matter = m_mbs.getMatterSubsystem();
    const SimbodyMatterSubsystemRep& matterRep = matter.getRep();
    const SBModelCache&              mc = matterRep.getModelCache(s);

    m_mbs.realize(s, Stage::Acceleration);
    const Vector& udot = s.getUDot(); // grab before invalidated
    s.updZ() += h*s.getZDot(); // invalidates Stage::Dynamics
    s.updU() += h*udot;        // invalidates Stage::Velocity
    Vector qdot;
    matter.multiplyByN(s,false,s.getU(),qdot);
    s.updQ() += h*qdot;         // invalidates Stage::Position
    matter.normalizeQuaternions(s);
    s.updTime() += h;           // invalidates Stage::Time
    m_mbs.realize(s, Stage::Dynamics);

    // Mark the acceleration-level calculations valid now, despite the fact
    // that they don't reflect the latest time, q, u, and forces.
    matterRep.markCacheValueRealized(s, mc.constrainedAccelerationCacheIndex);
    matterRep.markCacheValueRealized(s, mc.treeAccelerationCacheIndex);

    //TODO: prescribed motion, auto update state variables, events(?)
}


//------------------------------------------------------------------------------
//                      CALC COEFFICIENTS OF FRICTION
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//                     CALC COEFFICIENTS OF RESTITUTION
//------------------------------------------------------------------------------
// Calculate velocity-dependent coefficients of restitution.
void SemiExplicitEulerTimeStepper::
calcCoefficientsOfRestitution(const State& s, const Vector& verr,
                              bool disableRestitution) 
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    const Real defCaptureVel = getDefaultImpactCaptureVelocityInUse();
    const Real defMinCORVel  = getDefaultImpactMinCORVelocityInUse();

    if (disableRestitution) {
        SimTK_DEBUG("Restitution disabled; all COR=0\n");
        for (unsigned i=0; i < m_uniContact.size(); ++i) {
            ImpulseSolver::UniContactRT& rt = m_uniContact[i];
            rt.m_effCOR = 0;
        }
        return;
    }

    // These are the proximal unilateral constraint runtimes.
    // TODO: only need this for participating contacts.
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        Real impactVel = rt.m_sign*verr[rt.m_Nk];
        if (impactVel > 0) {
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
        rt.m_effCOR = uni.calcEffectiveCOR(s, defCaptureVel, defMinCORVel,
                                           -impactVel);
        SimTK_DEBUG3("Uni contact %d verr vel=%g -> cor=%g\n", (int)rt.m_ucx, 
                     impactVel, rt.m_effCOR);
    }
}

//------------------------------------------------------------------------------
//                      CALC NEWTON RESTITUTION IF ANY
//------------------------------------------------------------------------------
bool SemiExplicitEulerTimeStepper::
calcNewtonRestitutionIfAny(const State& s, const Vector& verr,
                           Vector& newtonVerr) const 
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();

    newtonVerr.resize(verr.size());
    newtonVerr.setToZero();
    bool anyRestitution = false;
    for (unsigned i=0; i < m_uniContact.size(); ++i) {
        const ImpulseSolver::UniContactRT& rt = m_uniContact[i];
        const MultiplierIndex              mx = rt.m_Nk;
        if (rt.m_effCOR > 0) {
            newtonVerr[mx] = rt.m_effCOR * verr[mx];
            anyRestitution = true;
        }
    }
    return anyRestitution;
}

//------------------------------------------------------------------------------
//                     APPLY POISSON RESTITUTION IF ANY
//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
//                      CALC EXPANSION IMPULSE IF ANY
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
//                         DO COMPRESSION PHASE
//------------------------------------------------------------------------------
// This phase uses all the proximal constraints and should use a starting
// guess for impulse saved from the last step if possible.
bool SemiExplicitEulerTimeStepper::
doCompressionPhase(const State& s, Vector& verr, Vector& compImpulse) {
#ifndef NDEBUG
    printf("DYN t=%.15g verr=", s.getTime()); cout << verr << endl;
#endif
    // TODO: improve initial guess
    m_expansionImpulse.setToZero(); //TODO: shouldn't need to zero this
    bool converged = m_solver->solve(0,
        m_allParticipating,m_GMInvGt,m_D,
        Array_<MultiplierIndex>(), m_expansionImpulse, 
        verr,compImpulse,
        m_unconditional,m_uniContact,m_uniSpeed,m_bounded,
        m_consLtdFriction, m_stateLtdFriction);
#ifndef NDEBUG
    m_solver->dumpUniContacts("Post-dynamics", m_uniContact);
#endif
    return converged;
}


//------------------------------------------------------------------------------
//                            DO EXPANSION PHASE
//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
//                        DO INDUCED IMPACT ROUND
//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
//                      DO POSITION CORRECTION PHASE
//------------------------------------------------------------------------------
// This phase uses only holonomic constraints.
// We require that m_uniContact has up-to-date information from the just-
// completed contact step so we can tell which contacts are active.
bool SemiExplicitEulerTimeStepper::
doPositionCorrectionPhase(const State& state, Vector& verr,
                          Vector& positionImpulse) {
    bool converged;
    if (m_projectionMethod == Unilateral) {
        SimTK_DEBUG1("UNILATERAL POSITION CORRECTION, %d participators\n",
                     (int)m_posParticipating.size());
        m_expansionImpulse.setToZero(); //TODO: shouldn't need to zero this
        converged = m_solver->solve(2,
            m_posParticipating,m_GMInvGt,m_D,
            Array_<MultiplierIndex>(), m_expansionImpulse,
            verr, positionImpulse,
            m_posUnconditional,m_posUniContact,m_posNoUniSpeed,m_posNoBounded,
            m_posNoConsLtdFriction, m_posNoStateLtdFriction);
    } else {
        // Bilateral. TODO: unconditionals
        m_participating.clear();
        for (unsigned k=0; k < m_uniContact.size(); ++k) {
            const ImpulseSolver::UniContactRT& rt = m_uniContact[k];
            const MultiplierIndex mx = rt.m_Nk;
            // Anything currently active, or inactive but penetrating, gets
            // a correction. TODO: is this the right strategy?
            if (rt.m_contactCond == ImpulseSolver::UniOff 
                && rt.m_sign*verr[mx] >= 0)
                continue;
            m_participating.push_back(rt.m_Nk);
        }
        SimTK_DEBUG1("BILATERAL POSITION CORRECTION, %d participators\n",
                    (int)m_participating.size());
        converged = m_solver->solveBilateral(m_participating,m_GMInvGt,m_D,
                                             verr, positionImpulse);
    }
    return converged;
}

//------------------------------------------------------------------------------
//                       ANY POSITION ERRORS VIOLATED
//------------------------------------------------------------------------------
bool SemiExplicitEulerTimeStepper::
anyPositionErrorsViolated(const State&, const Vector& perr) const {
    // TODO: no need to fix if large perrs satisfy inequalities.
    bool anyViolated = perr.normInf() > m_consTol;
    SimTK_DEBUG2("maxAbs(perr)=%g -> %s\n", perr.normInf(),
                anyViolated ? "VIOLATED" : "OK");
    return anyViolated;
}

//------------------------------------------------------------------------------
//                            DEBUGGING METHODS
//------------------------------------------------------------------------------
const char* SemiExplicitEulerTimeStepper::
getRestitutionModelName(RestitutionModel rm) {
    static const char* nm[]={"Poisson", "Newton", "NoRestitution"};
    return Poisson<=rm&&rm<=NoRestitution ? nm[rm] 
        : "UNKNOWNRestitutionModel";
}
const char* SemiExplicitEulerTimeStepper::
getInducedImpactModelName(InducedImpactModel iim) {
    static const char* nm[]={"Simultaneous", "Sequential", "Mixed"};
    return Simultaneous<=iim&&iim<=Mixed ? nm[iim] 
        : "UNKNOWNInducedImpactModel";
}
const char* SemiExplicitEulerTimeStepper::
getPositionProjectionMethodName(PositionProjectionMethod ppm) {
    static const char* nm[]={"Bilateral", "Unilateral", "NoPositionProjection"};
    return Bilateral<=ppm&&ppm<=NoPositionProjection ? nm[ppm] 
        : "UNKNOWNPositionProjectionMethod";
}
const char* SemiExplicitEulerTimeStepper::
getImpulseSolverTypeName(ImpulseSolverType ist) {
    static const char* nm[]={"PLUS", "PGS"};
    return PLUS<=ist&&ist<=PGS ? nm[ist] : "UNKNOWNImpulseSolverType";
}

} // namespace SimTK
