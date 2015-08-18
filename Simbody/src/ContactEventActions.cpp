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
    Array_<UnilateralContactIndex>& proximalUniContacts) 
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    proximalUniContacts.clear();

    const int nUniContacts  = matter.getNumUnilateralContacts();

    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = matter.getUnilateralContact(ux);
        if (   contact.isEnabled(state)
            || contact.isProximal(state, proximityTol)) // may be scaled
            proximalUniContacts.push_back(ux);
    }
}
//------------------------------------------------------------------------------
//                         COLLECT CONSTRAINT INFO
//------------------------------------------------------------------------------
// Create lists of low-level constraint properties suited for analysis by
// an Impulse solver. Only proximal constraints are considered, and those
// must already have been enabled so that we can collect their Simbody-assigned
// multipliers here.
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
        uniContact.push_back(); // default constructed
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
   (Real                                    consTol,
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
        } else if (rt.m_sign*verr[mz] >= (Real(-0.1))*consTol) {
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
    Real consTol, Real defCaptureVel, Real defMinCORVel,
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
        if (impactVel >= -consTol) {
            rt.m_effCOR = 0;
            SimTK_DEBUG3(
                "Uni contact %d verr vel=%g at or below tol %g -> cor=0\n", 
                 (int)rt.m_ucx, impactVel, consTol);
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

    const State& low = study.getCurrentState();
    State& high = study.updInternalState();

    printf("\n**********\nImpact Action (%.15g,%.15g], w=%.15g\n", 
            low.getTime(), high.getTime(), high.getTime()-low.getTime());
    for (auto& tp : triggers) {
        const EventTrigger::Witness& w = 
            dynamic_cast<const EventTrigger::Witness&>(*tp);
        printf("  Trigger %d '%s'=(%.15g,%.15g]\n", 
                (int)w.getEventTriggerId(),
                w.getTriggerDescription().c_str(),
                w.calcWitnessValue(study.getSystem(), low),
                w.calcWitnessValue(study.getSystem(), high));
    }

    const Real consTol = study.getConstraintToleranceInUse();

    // TODO: Above this speed the solver should assume sliding.
    m_solver->setMaxRollingSpeed(2*consTol);

    Array_<UnilateralContactIndex> proximalUniContacts;
    findProximalConstraints(mbs, high, consTol, proximalUniContacts);
    #ifndef NDEBUG
    cout<<"Proximal unilateral contacts: "<< proximalUniContacts << endl;
    #endif

    for (auto ucx : proximalUniContacts) { 
        const auto& uni = m_matter.getUnilateralContact(ucx);
        uni.enable(high);
    }
    mbs.realize(high, Stage::Instance); // assign multipliers

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
    const Vector noExpansion;
    Array_<int> impacters, expanders, observers;
    Array_<MultiplierIndex> participaters, expanding;
    classifyUnilateralContactsForSimultaneousImpact(consTol,
        verr0, noExpansion, uniContact,
        impacters, expanders, observers,
        participaters, expanding);

    cout << "CLASSIFIED:\n";
    cout << "  impacters=" << impacters << endl;
    cout << "  expanders=" << expanders << endl;
    cout << "  observers=" << observers << endl;

    // For simultaneous impact, CORs are determined only for the initial
    // impact velocities. In the reaction round they are treated as zero.
    calcCoefficientsOfRestitution(matter, high, verr0, 
        false, /*disableRestitution*/
        consTol, 2*consTol, 2*consTol, /*def capture/min COR vel*/
        uniContact);

    #ifndef NDEBUG
    ImpulseSolver::dumpUniContacts("UniContacts", uniContact);
    #endif


    // CALC COMPRESSION IMPULSE
    Vector verrStart(m), verrApplied(m);
    Vector compressionImpulse(m);
    Vector expansionImpulse(m, 0.);
    verrStart = verr0; verrApplied.setToZero();
    bool converged = m_solver->solve(0,
        participaters, A, D, expanding, 
        expansionImpulse, verrStart, verrApplied, // in/out
        compressionImpulse, // result
        unconditional, uniContact,
        uniSpeedDummy, boundedDummy, consLtdFrictionDummy, 
                                     stateLtdFrictionDummy);
    cout << "verr0=" << verr0 << " -> impulse " << compressionImpulse << endl;

    // constraint impulse -> gen impulse
    Vector genImpulse, deltaU;
    matter.multiplyByGTranspose(high, -compressionImpulse, genImpulse);
    // gen impulse to delta-u
    matter.multiplyByMInv(high, genImpulse, deltaU);
    high.updU() += deltaU;
    mbs.realize(high, Stage::Velocity);

    cout << "deltaU=" << deltaU << " now verr=" << verr0 << endl;
    mbs.realize(high, Stage::Acceleration);
    Vector mults = high.getMultipliers();
    cout << "mults=" << mults << endl;
    Array_<const UnilateralContact*> toDisable;
    for (auto ucx : proximalUniContacts) { 
        const auto& uni = m_matter.getUnilateralContact(ucx);
        const Real sign = uni.getSignConvention();
        const Real frc = -sign*mults[uni.getContactMultiplierIndex(high)];
        printf("uni contact %d: frc=%g (%s)\n", (int)ucx, frc,
               frc<0?"disable":"OK");
        if (frc < 0)
            toDisable.push_back(&uni);
    }
    for (auto ucp : toDisable)
        ucp->disable(high);
    mbs.realize(high, Stage::Position);

    m_matter.getSystem().project(high);
    printf("**********\n\n");
}

//==============================================================================
//                           CONTACT EVENT ACTION
//==============================================================================


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
        if (   contact.isEnabled(state)
            || (    contact.isProximal(state, proximityTol)
                && !contact.isSeparating(state, velocityTol)))
            lingeringUniContacts.push_back(ux);
    }
}

//------------------------------------------------------------------------------
//                              CHANGE VIRTUAL
//------------------------------------------------------------------------------
void ContactEventAction::
changeVirtual(Study&                  study,
              const Event&            event,
              const EventTriggers&    triggers,
              EventChangeResult&      result) const {
    //TODO: write this!
    const State& low = study.getCurrentState();
    State& high = study.updInternalState();
    const MultibodySystem& mbs = m_matter.getMultibodySystem();

    printf("\n##########\nContact Action (%.15g,%.15g], w=%.15g\n", 
            low.getTime(), high.getTime(), high.getTime()-low.getTime());
    for (auto& tp : triggers) {
        const EventTrigger::Witness& w = 
            dynamic_cast<const EventTrigger::Witness&>(*tp);
        mbs.realize(low, w.getDependsOnStage());
        mbs.realize(high, w.getDependsOnStage());
        printf("  Trigger %d '%s'=(%.15g,%.15g]\n", 
                (int)w.getEventTriggerId(),
                w.getTriggerDescription().c_str(),
                w.calcWitnessValue(study.getSystem(), low),
                w.calcWitnessValue(study.getSystem(), high));
    }

    const Real consTol = study.getConstraintToleranceInUse();
    //TODO: figure out pos & vel tol
    Array_<UnilateralContactIndex> lingeringUniContacts;
    findLingeringConstraints(mbs, high, consTol, consTol, lingeringUniContacts);
    #ifndef NDEBUG
    cout<<"Lingering unilateral contacts: "<< lingeringUniContacts << endl;
    #endif

    for (auto ucx : lingeringUniContacts) { 
        const auto& uni = m_matter.getUnilateralContact(ucx);
        uni.enable(high);
    }
    mbs.realize(high, Stage::Instance); // assign multipliers

    mbs.realize(high, Stage::Acceleration);
    const Vector& mults = high.getMultipliers();

    Array_<const UnilateralContact*> toDisable;
    for (auto ucx : lingeringUniContacts) {
        const auto& uni = m_matter.getUnilateralContact(ucx);
        const Real sign = uni.getSignConvention();
        const Real frc = -sign*mults[uni.getContactMultiplierIndex(high)];
        printf("uni contact %d: frc=%g (%s)\n", (int)ucx, frc,
               frc<0?"disable":"OK");
        if (frc < 0) // TODO: Not sufficient!!
            toDisable.push_back(&uni);
    }
    for (auto up : toDisable)
        up->disable(high);
    mbs.realize(high, Stage::Position);
    m_matter.getSystem().project(high);
    printf("##########\n");
}
