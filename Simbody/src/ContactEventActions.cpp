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
#include "simbody/internal/ConditionalConstraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include "ContactEventActions.h"
#include "SimbodyMatterSubsystemRep.h"

using namespace SimTK;

//==============================================================================
//                           IMPACT EVENT ACTION
//==============================================================================
void ImpactEventAction::
changeVirtual(Study&                  study,
              const Event&            event,
              const EventTriggers&    triggers,
              EventChangeResult&      result) const {
    //TODO: write this!
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

    Array_<const UnilateralContact*> toEnable;
    for (UnilateralContactIndex ucx(0); 
        ucx < m_matter.getNumUnilateralContacts(); ++ucx) {
        const auto& uni = m_matter.getUnilateralContact(ucx);
        if (uni.isProximal(high, study.getConstraintToleranceInUse()))
            toEnable.push_back(&uni);
    }
    for (auto up : toEnable)
        up->enable(high);
    m_matter.getSystem().project(high);
}

//==============================================================================
//                           CONTACT EVENT ACTION
//==============================================================================

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
