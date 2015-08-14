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
    const State&                            s,
    const Array_<UnilateralContactIndex>&   proximalUniContacts) 
{ //TODO: redo
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

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


    const int nConstraints = matter.getNumConstraints();
    for (ConstraintIndex cx(0); cx < nConstraints; ++cx) {
        // TODO: sort out the enabled, unconditional ones and fill in
        // m_unconditional and m_posUnconditional.
        // m_allParticipating.push_back(...);
        // if holnomic: m_posParticipating.push_back(...);
    }

    for (auto cx : proximalUniContacts) {
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
    }

    //TODO: speed, bounded, consLtd, proximal stateLtd friction 
    // position projection stuff
    // (all nonholonomic)
}

//------------------------------------------------------------------------------
//                              CHANGE VIRTUAL
//------------------------------------------------------------------------------
void ImpactEventAction::
changeVirtual(Study&                  study,
              const Event&            event,
              const EventTriggers&    triggers,
              EventChangeResult&      result) const {
    //TODO: write this!
    const MultibodySystem& mbs = MultibodySystem::downcast(study.getSystem());
    const State& low = study.getCurrentState();
    State& high = study.updInternalState();

    printf("Impact Action (%.15g,%.15g], w=%.15g\n", 
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



    m_matter.getSystem().project(high);
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
              EventChangeResult&      result) const {
    //TODO: write this!
    const State& low = study.getCurrentState();
    State& high = study.updInternalState();
    const MultibodySystem& mbs = m_matter.getMultibodySystem();

    printf("Contact Action (%.15g,%.15g], w=%.15g\n", 
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

    Array_<const UnilateralContact*> toDisable;
    for (UnilateralContactIndex ucx(0); 
        ucx < m_matter.getNumUnilateralContacts(); ++ucx) {
        const auto& uni = m_matter.getUnilateralContact(ucx);
        if (true)
            toDisable.push_back(&uni);
    }
    for (auto up : toDisable)
        up->disable(high);
    m_matter.getSystem().project(high);
}
