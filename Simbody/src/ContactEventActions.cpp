/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/ConditionalConstraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/ImpulseSolver.h"
#include "simbody/internal/PGSImpulseSolver.h"
#include "simbody/internal/PLUSImpulseSolver.h"

#include "ContactEventActions.h"
#include "SimbodyMatterSubsystemRep.h"

using namespace SimTK;

//==============================================================================
//                           IMPACT EVENT ACTION
//==============================================================================
//------------------------------------------------------------------------------
//                                CONSTRUCTOR
//------------------------------------------------------------------------------
ImpactEventAction::ImpactEventAction(const SimbodyMatterSubsystemRep& matter) 
:   EventAction(Change), m_matter(matter),
    m_solver(new PGSImpulseSolver(0.0)) {} // fill transition vel later
    //m_solver(new PLUSImpulseSolver(0.0)) {} // fill transition vel later

//------------------------------------------------------------------------------
//                         FIND PROXIMAL CONSTRAINTS
//------------------------------------------------------------------------------

/*static*/ void ImpactEventAction::
findProximalConstraints
   (const MultibodySystem&          mbs,
    const State&                    state,
    Real                            proximityTol,
    Array_<UnilateralContactIndex>& proximalUniContacts,
    Array_<UnilateralContactIndex>& needToActivate)
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    proximalUniContacts.clear();

    const int nUniContacts  = matter.getNumUnilateralContacts();

    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        if (contact.isActive(state)) {
            proximalUniContacts.push_back(ux);
            continue;
        }
        if (contact.isProximal(state, proximityTol)) { // may be scaled
            proximalUniContacts.push_back(ux);
            needToActivate.push_back(ux);
        }
    }
}

//------------------------------------------------------------------------------
//                         COLLECT CONSTRAINT INFO
//------------------------------------------------------------------------------
// Create lists of low-level constraint properties suited for analysis by
// an Impulse solver. Only proximal constraints are considered. We collect their
// Simbody-assigned multipliers here.
void ImpactEventAction::
collectConstraintInfo
   (const MultibodySystem&                  mbs,
    const State&                            state,
    const Array_<UnilateralContactIndex>&   proximalUniContacts,
    Array_<ImpulseSolver::UncondRT>&        unconditional,
    Array_<ImpulseSolver::UniContactRT>&    uniContact,
    Array_<MultiplierIndex>&                allParticipating) const
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // This will be filled out here since it can't change later.
    allParticipating.clear();

    // These are the proximal constraints for contact & impact.
    unconditional.clear();
    uniContact.clear();

    const int nConstraints = matter.getNumConstraints();
    for (ConstraintIndex cx(0); cx < nConstraints; ++cx) {
        // TODO: sort out the enabled, unconditional ones and fill in
        // m_unconditional and m_posUnconditional.
        // m_allParticipating.push_back(...);
        // if holnomic: m_posParticipating.push_back(...);
    }

    for (auto cx : proximalUniContacts) {
        const UnilateralContact& contact = matter.getUnilateralContact(cx);
        uniContact.emplace_back(); // default constructed
        ImpulseSolver::UniContactRT& rt = uniContact.back();
        rt.m_sign = (Real)contact.getSignConvention();
        rt.m_ucx = cx;
        rt.m_Nk = contact.getContactMultiplierIndex(state);
        allParticipating.push_back(rt.m_Nk);
        if (contact.hasFriction(state)) {
            rt.m_Fk.resize(2);
            contact.getFrictionMultiplierIndices(state, rt.m_Fk[0], rt.m_Fk[1]);
            allParticipating.push_back(rt.m_Fk[0]);
            allParticipating.push_back(rt.m_Fk[1]);
        }
    }

    //TODO: speed, bounded, consLtd, proximal stateLtd friction 
    // position projection stuff
    // (all nonholonomic)
}

//------------------------------------------------------------------------------
//            CLASSIFY UNILATERAL CONTACTS FOR SIMULTANEOUS IMPACT
//------------------------------------------------------------------------------
// Classify the unilateral contacts into impacters, observers, and expanders.
// and determine which constraint equations should participate assuming we
// are in Simultaneous ("one shot") impact mode. All impacter and observer 
// constraint equations should participate, as well as the friction constraints
// (but not the normal constraints) for expanders.
// If you don't want to process expanders, pass the expansionImpulse as a
// zero-length Vector, otherwise it must have length m.
void ImpactEventAction::
classifyUnilateralContactsForSimultaneousImpact
   (Real                                    posTol,
    Real                                    velTol,
    const Vector&                           verr,
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
        } else if (rt.m_sign*verr[mz] >= (Real(-0.1))*velTol) {
            //printf("rt.m_sign*verr[mz]=%.15g\n", rt.m_sign*verr[mz]);
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
//                     CALC COEFFICIENTS OF RESTITUTION
//------------------------------------------------------------------------------
// Calculate velocity-dependent coefficients of restitution.
void ImpactEventAction::
calcCoefficientsOfRestitution
   (const SimbodyMatterSubsystem& matter,
    const State& s, const Vector& verr,
    bool disableRestitution,
    Real velTol, Real defCaptureVel, Real defMinCORVel,
    Array_<ImpulseSolver::UniContactRT>& uniContact) const
{
    if (disableRestitution) {
        SimTK_DEBUG("Restitution disabled; all COR=0\n");
        for (unsigned i=0; i < uniContact.size(); ++i) {
            ImpulseSolver::UniContactRT& rt = uniContact[i];
            rt.m_effCOR = 0;
        }
        return;
    }

    // These are the proximal unilateral constraint runtimes.
    // TODO: only need this for participating contacts.
    for (unsigned i=0; i < uniContact.size(); ++i) {
        ImpulseSolver::UniContactRT& rt = uniContact[i];
        Real impactVel = rt.m_sign*verr[rt.m_Nk];
        if (impactVel > 0) {
            rt.m_effCOR = 0; // separating
            continue;
        }
        if (impactVel >= -velTol) {
            rt.m_effCOR = 0;
            SimTK_DEBUG3(
                "Uni contact %d verr vel=%g at or below tol %g -> cor=0\n", 
                 (int)rt.m_ucx, impactVel, velTol);
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
//                      CALC EXPANSION IMPULSE IF ANY
//------------------------------------------------------------------------------
/* We just generated a compression impulse, i.e., an impulse that brings all
the impact velocities to zero. Figure out if we need to apply an expansion
impulse and if so how much and to whom. */
bool ImpactEventAction::
calcExpansionImpulseIfAny
   (const Array_<ImpulseSolver::UniContactRT>&  uniContacts,
    const Array_<int>&                          impacters,
    const Vector&                               compressionImpulse,
    Vector&                                     expansionImpulse,
    Array_<int>&                                expanders) const 
{
    expansionImpulse.resize(compressionImpulse.size());
    expansionImpulse.setToZero();
    expanders.clear();

    bool anyRestitution = false;
    for (unsigned i=0; i < impacters.size(); ++i) {
        const int which = impacters[i];
        const ImpulseSolver::UniContactRT& rt = uniContacts[which];
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
//                              CHANGE VIRTUAL
//------------------------------------------------------------------------------
void ImpactEventAction::
changeVirtual(Study&                  study,
              const Event&            event,
              const EventTriggers&    triggers,
              EventChangeResult&      result) const 
{
    const MultibodySystem& mbs = MultibodySystem::downcast(study.getSystem());
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    const Real posTol = study.getConstraintToleranceInUse();
    const Real velTol = posTol / study.getTimeScaleInUse();

    const State& low = study.getCurrentState();
    State& high = study.updInternalState();

#ifndef NDEBUG
    printf("\n**********\nImpact Action (%.15g,%.15g], w=%.15g (p,vtol=%g,%g)\n", 
            low.getTime(), high.getTime(), high.getTime()-low.getTime(),
            posTol, velTol);
    for (auto& tp : triggers) {
        const EventWitness& w = 
            dynamic_cast<const EventWitness&>(*tp);
        const EventWitness::Value lv = w.calcWitnessValue(study, low);
        const EventWitness::Value hv = w.calcWitnessValue(study, high);
        printf("  Witness %d (trigger %d) '%s'=(%.15g(%d),%.15g(%d)]\n", 
                (int)w.getEventWitnessIndex(), (int)w.getEventTriggerId(),
                w.getTriggerDescription().c_str(),
                lv.getValue(), lv.getSign(), hv.getValue(), hv.getSign());
    }
#endif

    // TODO: Above this speed the solver should assume sliding.
    m_solver->setMaxRollingSpeed(2*velTol);

    Array_<UnilateralContactIndex> proximalUniContacts, needToActivate;
    findProximalConstraints(mbs, high, posTol, 
                            proximalUniContacts, needToActivate);
    #ifndef NDEBUG
    cout<<"Proximal unilateral contacts: "<< proximalUniContacts << endl;
    cout<<"  need to activate: "<< needToActivate << endl;
    #endif

    for (auto ucx : needToActivate) { 
        const auto& uni = m_matter.getUnilateralContact(ucx);
        uni.setCondition(high, CondConstraint::Active); // friction unknown
    }

    Array_<MultiplierIndex>                 allParticipating;
    Array_<ImpulseSolver::UncondRT>         unconditional;
    Array_<ImpulseSolver::UniContactRT>     uniContact;

    Array_<ImpulseSolver::UniSpeedRT>               uniSpeedDummy;
    Array_<ImpulseSolver::BoundedRT>                boundedDummy;
    Array_<ImpulseSolver::ConstraintLtdFrictionRT>  consLtdFrictionDummy;
    Array_<ImpulseSolver::StateLtdFrictionRT>       stateLtdFrictionDummy;

    collectConstraintInfo(mbs, high, proximalUniContacts,
                          unconditional, uniContact, allParticipating);

    mbs.realize(high, Stage::Velocity);

    // Note that verr0 is a reference into State s so if we update the
    // velocities in s and realize them, verr0 will also be updated.
    const Vector& verr0 = high.getUErr();

    // m is the total number of proximal constraint equations. This won't 
    // change during the step.
    const int m = verr0.size();

    // Calculate the constraint compliance matrix A=GM\~G and soft constraint
    // diagonal D.
    Matrix A(m,m); Vector D(m);
    matter.calcProjectedMInv(high, A); // m X m
    D.setToZero();

    // CLASSIFY FOR COMPRESSION
    const Vector noExpansion(m, 0.);
    Array_<int> impacters, expanders, observers;
    Array_<MultiplierIndex> participaters, expanding;
    classifyUnilateralContactsForSimultaneousImpact(posTol,velTol,
        verr0, noExpansion, uniContact,
        impacters, expanders, observers,
        participaters, expanding);

    #ifndef NDEBUG
    cout << "CLASSIFIED FOR COMPRESSION:\n";
    cout << "  impacters=" << impacters << endl;
    cout << "  expanders=" << expanders << endl;
    cout << "  observers=" << observers << endl;
    #endif

    // For simultaneous impact, CORs are determined only for the initial
    // impact velocities. In the reaction round they are treated as zero.
    calcCoefficientsOfRestitution(matter, high, verr0, 
        false, /*disableRestitution*/
        velTol, 2*velTol, 2*velTol, /*def capture/min COR vel*/
        uniContact);

    #ifndef NDEBUG
    ImpulseSolver::dumpUniContacts("UniContacts", uniContact);
    #endif


    //----------------------PERFORM COMPRESSION ROUND---------------------------
    Vector verr(m), verrApplied(m);
    Vector totalImpulse(m);
    verr = verr0; verrApplied.setToZero();
    bool converged = m_solver->solve(0,
        participaters, A, D, expanding, 
        noExpansion, verr, verrApplied, // in/out
        totalImpulse, // compression impulse is all so far
        unconditional, uniContact,
        uniSpeedDummy, boundedDummy, consLtdFrictionDummy, 
                                     stateLtdFrictionDummy);

    // Calculate the next expansion impulse from the compression impulses
    // we just generated. This only occurs if we are doing Poisson restitution.
    Vector expansionImpulse(m);
    const bool anyExpansion = 
        calcExpansionImpulseIfAny(uniContact, impacters, totalImpulse,
                                  expansionImpulse, expanders);
    #ifndef NDEBUG
    cout << "  Postcompression verr=" << verr << endl;
    cout << "    impulse=" << totalImpulse << endl;
    if (!anyExpansion) cout << "  No expansion needed." << endl;
    #endif  

    if(anyExpansion) {
        //---------------------PERFORM EXPANSION ROUND--------------------------
        #ifndef NDEBUG
        cout << "  Next expansion impulse=" << expansionImpulse << endl;
        cout << "  Next expanders: " << expanders << endl;
        #endif 
        classifyUnilateralContactsForSimultaneousImpact(posTol,velTol,
            verr, expansionImpulse, uniContact,
            impacters, expanders, observers,
            participaters, expanding);

        #ifndef NDEBUG
        cout << "CLASSIFIED FOR EXPANSION:\n";
        cout << "  impacters=" << impacters << endl;
        cout << "  expanders=" << expanders << endl;
        cout << "  observers=" << observers << endl;
        #endif

        Vector reactionImpulse(m);
        bool converged = m_solver->solve(0,
            participaters, A, D, expanding, 
            expansionImpulse, verr, verrApplied, // in/out
            reactionImpulse, // result
            unconditional, uniContact,
            uniSpeedDummy, boundedDummy, consLtdFrictionDummy, 
                                         stateLtdFrictionDummy);
        // verr has been updated here to verr-A*(expansionImpulse+reactionImpulse)
        #ifndef NDEBUG
        cout << "  expansion impulse: " << expansionImpulse << endl;
        cout << "  reaction impulse: " << reactionImpulse << endl;
        #endif    

        totalImpulse += expansionImpulse;
        totalImpulse += reactionImpulse;

        // That should be the end but it is possible that one of the expanders
        // is now impacting.
        classifyUnilateralContactsForSimultaneousImpact(posTol,velTol,
            verr, noExpansion, uniContact,
            impacters, expanders, observers, 
            participaters, expanding);  
        assert(expanders.empty() && expanding.empty());

        if (!impacters.empty()) {
            //-------------------PERFORM CORRECTION ROUND-----------------------
            #ifndef NDEBUG
            cout << "CLASSIFIED FOR CORRECTION:\n";
            cout << "  impacters=" << impacters << endl;
            cout << "  expanders=" << expanders << endl;
            cout << "  observers=" << observers << endl;
            #endif

            bool converged = m_solver->solve(0,
                participaters, A, D, expanding, 
                noExpansion, verr, verrApplied, // in/out
                reactionImpulse, // result
                unconditional, uniContact,
                uniSpeedDummy, boundedDummy, consLtdFrictionDummy, 
                                             stateLtdFrictionDummy);
            // verr has been updated here to verr-A*reactionImpulse
            #ifndef NDEBUG
            cout << "  correction impulse: " << reactionImpulse << endl;
            cout << "  post-correction verr: " << verr << endl;
            #endif

            totalImpulse += reactionImpulse;
        }
    }


    // constraint impulse -> gen impulse
    Vector genImpulse, deltaU;
    matter.multiplyByGTranspose(high, -totalImpulse, genImpulse);
    // gen impulse to delta-u
    matter.multiplyByMInv(high, genImpulse, deltaU);
    high.updU() += deltaU;
    mbs.realize(high, Stage::Velocity);
#ifndef NDEBUG
    cout << "deltaU=" << deltaU << " now verr=" << verr0 << endl;
#endif

    // We have potentially modified all system velocities during this
    // impact so it no longer matters which constraints were active on entry
    // to this handler. We want updateActiveSet() to consider only those
    // constraints that currently meet tolerances.
    // TODO 1: should still consider pre-active position constraints satisfied
    // if we did on entry; but only their position violations, not their 
    // velocity violations, can exceed tolerance.
    // TODO 2: constraints that were active on entry but not affected *at all*
    // by the impact (like they were fully decoupled from the impacters)
    // should probably NOT get turned off here.
    for (UnilateralContactIndex ucx(0); 
         ucx < m_matter.getNumUnilateralContacts(); ++ucx) 
    { 
        const auto& uni = m_matter.getUnilateralContact(ucx);
        uni.setCondition(high, CondConstraint::Off);
    }

    // Decide which constraints should be active. This must be the same 
    // method as is used in the Contact Action. It chooses an active set and
    // handles projection.
    ContactEventAction::updateActiveSet(mbs, high, posTol, velTol);

    SimTK_DEBUG("\n********** End of ImpactEventAction\n");
}

//==============================================================================
//                           CONTACT EVENT ACTION
//==============================================================================


//------------------------------------------------------------------------------
//                              CHANGE VIRTUAL
//------------------------------------------------------------------------------

void ContactEventAction::
changeVirtual(Study&                  study,
              const Event&            event,
              const EventTriggers&    triggers,
              EventChangeResult&      result) const 
{
    const MultibodySystem& mbs = m_matter.getMultibodySystem();
    const Real posTol = study.getConstraintToleranceInUse();
    const Real velTol = posTol / study.getTimeScaleInUse();

    const State& low = study.getCurrentState();
    State& high = study.updInternalState();

#ifndef NDEBUG
    printf("\n##########\n"
           "Contact Action (%.15g,%.15g], w=%.15g ptol=%g vtol=%g\n", 
            low.getTime(), high.getTime(), high.getTime()-low.getTime(),
           posTol, velTol);
    for (auto& tp : triggers) {
        const EventWitness& w = 
            dynamic_cast<const EventWitness&>(*tp);
        mbs.realize(low, w.getDependsOnStage());
        mbs.realize(high, w.getDependsOnStage());
        const EventWitness::Value lv = w.calcWitnessValue(study, low);
        const EventWitness::Value hv = w.calcWitnessValue(study, high);
        printf("  Trigger %d '%s'=(%.15g(%d),%.15g(%d)]\n", 
                (int)w.getEventTriggerId(),
                w.getTriggerDescription().c_str(),
                lv.getValue(), lv.getSign(),
                hv.getValue(), hv.getSign());
    }
#endif

    updateActiveSet(mbs, high, posTol, velTol);

#ifndef NDEBUG
    printf("\n########## end of ContactEventAction\n");
#endif
}

//------------------------------------------------------------------------------
//                            UPDATE ACTIVE SET
//------------------------------------------------------------------------------

/*static*/ void ContactEventAction::updateActiveSet
       (const MultibodySystem&          mbs,
        State&                          state, // in/out
        Real                            posTol,
        Real                            velTol) 
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Here we can only activate more contacts; lingering includes everything
    // currently active plus any inactive lingerers (perr=verr=0).
    Array_<UnilateralContactIndex> lingeringUniContacts;
    findLingeringConstraints(mbs, state, posTol, velTol, 
                             lingeringUniContacts);
    #ifndef NDEBUG
    cout <<"updateActiveSet() from " << findActiveSet(mbs,state) << endl;    
    cout <<"  All candidates: "<< lingeringUniContacts << endl;
    if (lingeringUniContacts.empty())
        cout <<"  NO LINGERING CONTACTS????\n";  
    #endif

    // If there aren't any candidates, that also means there were no conditional
    // constraints active, so there is nothing to do.
    if (lingeringUniContacts.empty())
        return;

    unsigned n = lingeringUniContacts.size();

    SimTK_ASSERT_ALWAYS(n<=64, "Max for combinatoric search is 64 constraints.");

    // Activate everything first.
    Real norm2Lambda;
    std::pair<UnilateralContactIndex,Real> worstForce, worstAcc;
    if (scoreActiveSet(mbs, state, lingeringUniContacts, lingeringUniContacts,
                       norm2Lambda, worstForce, worstAcc)) {
        SimTK_DEBUG1("  All active was a good solution (norm=%g)\n", 
                     norm2Lambda);
        mbs.realize(state, Stage::Position);
        matter.getSystem().project(state);
        return;
    }

    unsigned long long bestCombo=0;
    Real bestNorm2Lambda = Infinity;

    // Search all non-zero combinations for the best one.
    Array_<UnilateralContactIndex> actives;
    for (unsigned k=n-1; k > 0; --k) {
        SimTK_DEBUG2("  NOT GOOD @%d actives: Try choosing %d.\n", k+1, k);
        unsigned long long combo = 0;
        while (nextBitCombination(n, k, combo)) {
            actives.clear();
            for (unsigned long long bits=combo; bits; bits=clearLowBit(bits))
                actives.push_back(lingeringUniContacts[lowBitIndex(bits)]);
            if (scoreActiveSet(mbs, state, lingeringUniContacts,
                               actives, norm2Lambda, worstForce, worstAcc)) {
                if (norm2Lambda < bestNorm2Lambda) {
                    SimTK_DEBUG2("New winner with %d actives, norm=%g\n",
                                 k, norm2Lambda);
                    bestCombo=combo;
                    bestNorm2Lambda = norm2Lambda;
                }
            }
            if (bestCombo) {
                SimTK_DEBUG2("  Good solution w/%d actives, norm=%g\n",
                             k, bestNorm2Lambda);
                mbs.realize(state, Stage::Position);
                matter.getSystem().project(state);
                return;
            }
        }
    }

    // No good solution with any constraints on. Deactivate everything.
    if (scoreActiveSet(mbs, state, 
                       lingeringUniContacts, Array_<UnilateralContactIndex>(),
                       norm2Lambda, worstForce, worstAcc)) {
            SimTK_DEBUG("  No actives was a good solution\n");
    } else {
            SimTK_DEBUG("  COULDN'T FIND A GOOD SOLUTION; deactivated\n");
    }

    mbs.realize(state, Stage::Position);
    matter.getSystem().project(state);
}

//------------------------------------------------------------------------------
//                        FIND LINGERING CONSTRAINTS
//------------------------------------------------------------------------------

/*static*/ void ContactEventAction::
findLingeringConstraints
   (const MultibodySystem&          mbs,
    const State&                    state,
    Real                            proximityTol,
    Real                            velocityTol,
    Array_<UnilateralContactIndex>& lingeringUniContacts) 
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    lingeringUniContacts.clear();

    const int nUniContacts  = matter.getNumUnilateralContacts();

    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        if (contact.isActive(state)) {
            lingeringUniContacts.push_back(ux);
            continue;
        }

        if (    contact.isProximal(state, proximityTol)
            && !contact.isSeparating(state, velocityTol)) {
            lingeringUniContacts.push_back(ux);
        }
    }
}

//------------------------------------------------------------------------------
//                            FIND ACTIVE SET
//------------------------------------------------------------------------------
/*static*/ Array_<UnilateralContactIndex> ContactEventAction::
findActiveSet
   (const MultibodySystem&          mbs,
    const State&                    state) 
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    Array_<UnilateralContactIndex> activeContacts;

    const int nUniContacts  = matter.getNumUnilateralContacts();
    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        if (contact.isActive(state))
            activeContacts.push_back(ux);
    }
    return activeContacts;
}

//------------------------------------------------------------------------------
//                           ACTIVATE ACTIVE SET
//------------------------------------------------------------------------------
/*static*/ void ContactEventAction::
activateActiveSet
       (const MultibodySystem&                  mbs,
        State&                                  state,
        const Array_<UnilateralContactIndex>&   fullSet,
        const Array_<UnilateralContactIndex>&   activeSubset)
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Deactivate everything.
    for (auto ux : fullSet) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        contact.setCondition(state, CondConstraint::Off);
    }

    // Activate the active ones.
    for (auto ux : activeSubset) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        contact.setCondition(state, CondConstraint::Active);
    }
}

//------------------------------------------------------------------------------
//                           SCORE ACTIVE SUBSET
//------------------------------------------------------------------------------
// Given an active subset to try, activates them, calculates multipliers and
// acclerations, and looks for resulting violations of inequality conditions.
// The 2-norm of all the multipliers is used as a score since we want the least
// squares solution.
/*static*/ bool ContactEventAction::
scoreActiveSet
       (const MultibodySystem&                  mbs,
        State&                                  state,
        const Array_<UnilateralContactIndex>&   lingeringUniContacts,
        const Array_<UnilateralContactIndex>&   activeSubset,
        Real&                                   norm2Lambda,
        std::pair<UnilateralContactIndex,Real>& worstForce,
        std::pair<UnilateralContactIndex,Real>& worstAcc)
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();


    activateActiveSet(mbs,state,lingeringUniContacts,activeSubset);

    mbs.realize(state, Stage::Acceleration);
    const Vector& mults = state.getMultipliers();
    norm2Lambda = mults.norm();

    worstForce.first.invalidate(); worstForce.second=0;
    worstAcc.first.invalidate(); worstAcc.second=0;

    bool isGood = true;
    Array_<Real> forces, accs;
    for (auto ux : lingeringUniContacts) {
        const auto& uni = matter.getUnilateralContact(ux);
        const auto cond = uni.getCondition(state);
        const Real sign = uni.getSignConvention();
        if (cond == CondConstraint::Off) {
            const Real acc = sign*uni.getAerr(state);
            accs.push_back(acc);
            if (acc < worstAcc.second) {
                worstAcc.first = ux; worstAcc.second = acc;
                isGood = false;
            }
        } else {
            assert(cond == CondConstraint::Active);
            const Real frc = -sign*mults[uni.getContactMultiplierIndex(state)];
            forces.push_back(frc);
            if (frc < worstForce.second) {
                worstForce.first = ux; worstForce.second = frc;
                isGood = false;
            }
        }
    }

#ifndef NDEBUG
    cout << "scoreActiveSet: " << activeSubset << endl;
    cout << "  active forces=" << forces << endl;
    cout << "  inactive accs=" << accs << endl;

    cout << std::boolalpha << "  isGood=" << isGood 
         << " norm2=" << norm2Lambda << endl;
    if (worstForce.first.isValid())
        cout << "  worstForce uni " << worstForce.first 
             << ": " << worstForce.second << endl;
    if (worstAcc.first.isValid())
        cout << "  worstAcc uni " << worstAcc.first 
             << ": " << worstAcc.second << endl;
#endif

    return isGood;
}

