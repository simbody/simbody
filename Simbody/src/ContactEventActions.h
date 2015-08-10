#ifndef SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_
#define SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_

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

// Collects declarations of event actions associated with event-based 
// conditional constraint handling. Currently:
//      ImpactEventAction
//      ContactEventAction
// 

#include "simbody/internal/common.h"

namespace SimTK {
class SimbodyMatterSubsystemRep;

//==============================================================================
//                           IMPACT EVENT ACTION
//==============================================================================
// This is the EventAction we will perform when an ImpactEvent triggers.
class ImpactEventAction : public EventAction {
public:
    explicit ImpactEventAction(const SimbodyMatterSubsystemRep& matter) 
    :   EventAction(Change), m_matter(matter) {}

    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override;

private:
    ImpactEventAction* cloneVirtual() const override 
    {   return new ImpactEventAction(*this); }

    const SimbodyMatterSubsystemRep&    m_matter;
};


//==============================================================================
//                           CONTACT EVENT ACTION
//==============================================================================
// This is the EventAction we will perform when an ContactEvent triggers.
class ContactEventAction : public EventAction {
public:
    explicit ContactEventAction(const SimbodyMatterSubsystemRep& matter) 
    :   EventAction(Change), m_matter(matter) {}

    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override;

private:
    ContactEventAction* cloneVirtual() const override 
    {   return new ContactEventAction(*this); }

    const SimbodyMatterSubsystemRep&    m_matter;
};

} // namespace


#endif // SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_


