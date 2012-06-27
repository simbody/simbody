/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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
#include "SimTKcommon/internal/Event.h"

#include <cassert>
#include <string>

namespace SimTK {


const char* Event::getCauseName(Cause cause) {
    switch(cause) {
    case Cause::Initialization:    return "Initialization";
    case Cause::Triggered:         return "Triggered";
    case Cause::Scheduled:         return "Scheduled";
    case Cause::TimeAdvanced:      return "TimeAdvanced";
    case Cause::Signaled:          return "Signaled";
    case Cause::Termination:       return "Termination";
    case Cause::Invalid:           return "Invalid";
    }
    return "UNRECOGNIZED EVENT CAUSE";
}


std::string Event::eventTriggerString(Trigger e) {
    // Catch special combos first
    if (e==NoEventTrigger)        return "NoEventTrigger";
    if (e==Falling)               return "Falling";
    if (e==Rising)                return "Rising";
    if (e==AnySignChange)         return "AnySignChange";

    // Not a special combo; unmask one at a time.
    const Trigger triggers[] =
     { PositiveToNegative,NegativeToPositive,NoEventTrigger };
    const char *triggerNames[] =
     { "PositiveToNegative","NegativeToPositive" };

    String s;
    for (int i=0; triggers[i] != NoEventTrigger; ++i)
        if (e & triggers[i]) {
            if (s.size()) s += "|";
            s += triggerNames[i];
            e = Trigger((unsigned)e & ~((unsigned)triggers[i])); 
        }

    // should have accounted for everything by now
    if (e != NoEventTrigger) {
        char buf[128];
        std::sprintf(buf, "0x%x", (unsigned)e);
        if (s.size()) s += " + ";
        s += "UNRECOGNIZED EVENT TRIGGER GARBAGE ";
        s += buf;
    }
    return s;
}



////////////////////////////
// EVENT TRIGGER INFO REP //
////////////////////////////

class EventTriggerInfo::EventTriggerInfoRep {
public:
    explicit EventTriggerInfoRep(EventTriggerInfo* h)
      : myHandle(h), eventId(EventId(InvalidIndex)), triggerOnRising(true), triggerOnFalling(true), localizationWindow(0.1)
    {
        assert(h);
    }

private:
    EventTriggerInfo* myHandle;
    friend class EventTriggerInfo;

    EventId  eventId;
    bool triggerOnRising;
    bool triggerOnFalling;
    Real localizationWindow;
};



    ////////////////////////
    // EVENT TRIGGER INFO //
    ////////////////////////

EventTriggerInfo::EventTriggerInfo() : rep(0) {
    rep = new EventTriggerInfoRep(this);
}
EventTriggerInfo::~EventTriggerInfo() {
    if (getRep().myHandle == this)
        delete rep;
    rep = 0;
}

EventTriggerInfo::EventTriggerInfo(EventId eventId) : rep(0) {
    rep = new EventTriggerInfoRep(this);
    rep->eventId = eventId;
}

EventTriggerInfo::EventTriggerInfo(const EventTriggerInfo& src) : rep(0) {
    rep = new EventTriggerInfoRep(src.getRep());
    rep->myHandle = this;
}

EventTriggerInfo& 
EventTriggerInfo::operator=(const EventTriggerInfo& src) {
    if (&src != this) {
        if (getRep().myHandle == this)
            delete rep;
        rep = new EventTriggerInfoRep(src.getRep());
        rep->myHandle = this;
    }
    return *this;
}

EventId EventTriggerInfo::getEventId() const {
    return getRep().eventId;
}
bool EventTriggerInfo::shouldTriggerOnRisingSignTransition() const {
    return getRep().triggerOnRising;
}
bool EventTriggerInfo::shouldTriggerOnFallingSignTransition() const {
    return getRep().triggerOnFalling;
}
Real EventTriggerInfo::getRequiredLocalizationTimeWindow()    const {
    return getRep().localizationWindow;
}

EventTriggerInfo& 
EventTriggerInfo::setEventId(EventId id) {
    updRep().eventId = id; 
    return *this;
}
EventTriggerInfo& 
EventTriggerInfo::setTriggerOnRisingSignTransition(bool shouldTrigger) {
    updRep().triggerOnRising = shouldTrigger; 
    return *this;
}
EventTriggerInfo& 
EventTriggerInfo::setTriggerOnFallingSignTransition(bool shouldTrigger) {
    updRep().triggerOnFalling = shouldTrigger; 
    return *this;
}
EventTriggerInfo& 
EventTriggerInfo::setRequiredLocalizationTimeWindow(Real w) {
    assert(w > 0);
    updRep().localizationWindow = w; 
    return *this;
}
    ////////////////////////////
    // EVENT TRIGGER INFO REP //
    ////////////////////////////

// All inline currently.

} // namespace SimTK

