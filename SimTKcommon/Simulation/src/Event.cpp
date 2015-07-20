/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-14 Stanford University and the Authors.        *
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


/**@file
 *
 * Implementation of non-inline methods from the Event classes.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Event_Defs.h"
#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventAction.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventTrigger_Timer.h"
#include "SimTKcommon/internal/EventTrigger_Witness.h"

#include <cassert>
#include <string>

using namespace SimTK;

// Not inline because needs a definition for EventAction.
Event::Event(const std::string& eventDescription) 
:   m_eventDescription(eventDescription) {initialize();}
Event::Event(std::string&& eventDescription) 
:   m_eventDescription(std::move(eventDescription)) {initialize();}

EventActionIndex Event::adoptEventAction(EventAction* actionp) {
    SimTK_APIARGCHECK_ALWAYS(actionp != nullptr, "Event", "adoptEventAction",
                             "Action pointer can't be null.");
    m_actions.emplace_back(actionp);
    if (actionp->getActionMask() & EventAction::AnyChange)
        m_hasChangeAction = true;
    if (actionp->getActionMask() & EventAction::AnyReport)
        m_hasReportAction = true;
    return EventActionIndex(m_actions.size()-1);
}

void Event::performReportActions
   (const Study&            study,
    const EventTriggers&    triggers) const {
    for (auto&& action : m_actions) {
        if (action->getActionMask() & EventAction::AnyReport)
            action->report(study, *this, triggers);
    }
}

void Event::performChangeActions
   (Study&                  study,
    const EventTriggers&    triggers,
    EventChangeResult&      result) const {
    // Set up to provide a successful return if there are no change actions
    // to be performed.

    for (auto& action : m_actions) {
        if (!(action->getActionMask() & EventAction::AnyChange))
            continue;

        action->change(study, *this, triggers, result);
    }
}


const char* Event::getCauseName(EventCause cause) {
    switch(cause) {
    case EventCause::Initialization:    return "Initialization";
    case EventCause::Triggered:         return "Triggered";
    case EventCause::Scheduled:         return "Scheduled";
    case EventCause::TimeAdvanced:      return "TimeAdvanced";
    case EventCause::Signaled:          return "Signaled";
    case EventCause::Termination:       return "Termination";
    case EventCause::Invalid:           return "Invalid";
    }
    return "UNRECOGNIZED EVENT CAUSE";
}

std::string Event::eventTriggerDirectionString(TriggerDirection e) {
    // Catch special combos first
    if (e==NoEventTrigger)        return "NoEventTrigger";
    if (e==Falling)               return "Falling";
    if (e==Rising)                return "Rising";
    if (e==AnySignChange)         return "AnySignChange";

    // Not a special combo; unmask one at a time.
    const TriggerDirection triggerDirs[] =
     { PositiveToNegative,NegativeToPositive,NoEventTrigger };
    const char *triggerDirNames[] =
     { "PositiveToNegative","NegativeToPositive" };

    String s;
    for (int i=0; triggerDirs[i] != NoEventTrigger; ++i)
        if (e & triggerDirs[i]) {
            if (s.size()) s += "|";
            s += triggerDirNames[i];
            e = TriggerDirection((unsigned)e & ~((unsigned)triggerDirs[i])); 
        }

    // should have accounted for everything by now
    if (e != NoEventTrigger) {
        if (s.size()) s += " + ";
        s += "UNRECOGNIZED EVENT TRIGGER GARBAGE ";
        s += String((unsigned)e, "0x%x");
    }
    return s;
}


