#ifndef SimTK_SimTKCOMMON_EVENT_ACTION_H_
#define SimTK_SimTKCOMMON_EVENT_ACTION_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Event_Defs.h"

#include <string>

namespace SimTK {

class System;
class State;
class Study;

//==============================================================================
//                             EVENT ACTION
//==============================================================================
/** An %EventAction is an object that is invoked when a particular Event
occurs while a System is advanced through time. The action may consist of a 
read-only Report action, just reporting that the Event occurred, or a Change
action, in which case an event handler is invoked that may change the current 
state and thus affect the subsequent evolution of the System trajectory.

A Change action always causes the time stepper to stop the simulation at the
time of Event occurrence, and then restart afterwards with the possibly-modified
state. For a Report action, however, it is permissible for the time stepper
to provide an interpolated state that is an approximation of the actual 
trajectory. A Report action may optionally insist that the time stepper 
interrupt the simulation in the same manner as a Change action would do, so
that the report occurs exactly at one of the integration steps.

Every %EventAction is owned by a particular Event object; if you want the same
action performed for different Events you will need to make duplicate 
%EventAction objects. There is a clone() method available. **/
class SimTK_SimTKCOMMON_EXPORT EventAction {
public:
    /** Destructor is public and virtual for good cleanup behavior. **/
    virtual ~EventAction() = default;

    /** Change the action or actions to be performed; this is usually set on
    construction. The concrete %EventAction must have implemented the 
    appropriate method(s) for the actions you select; otherwise you will get a 
    runtime error when that method is needed. **/
    EventAction& setActionMask(unsigned actionMask) 
    {   m_actionMask = actionMask; return *this; }

    /** Set this to `false` to force the time stepper to complete a step prior 
    to invoking a Report action here; normally interpolation is permitted for 
    Reports. This flag has no effect on Change actions, which never allow
    interpolation because they must be able to cause an immediate change in
    the simulation trajectory. **/
    EventAction& setAllowInterpolationForReport(bool allowInterpolation) 
    {   m_allowInterpolationForReport = allowInterpolation; return *this; }

    /** Create a new instance of this concrete EventAction object. **/
    EventAction* clone() const {return cloneVirtual();}

    /** Return the set of actions to be performed. You can "or" this with the
    bitmask constants provided by the Reponse enum to see what is set. **/
    unsigned getActionMask() const 
    {   return m_actionMask; }

    /** Return whether we are allowing the time stepper to provide
    interpolated states when doing a Report action (that is the default). **/
    bool getAllowInterpolationForReport() const 
    {   return m_allowInterpolationForReport; }

    /** What methods does this EventAction want called when the event occurs? 
    These can be or-ed together. **/
    enum Response {
        Nothing     = 0x00, ///< Don't do anything.
        Report      = 0x01, ///< Report state that triggered.
        Change      = 0x02, ///< Given triggered state, fix it.
        PreReport   = 0x04, ///< Report state just prior to triggering.
        PostReport  = 0x08, ///< Report cleaned-up, post-handled state.
        AnyReport   = PreReport|Report|PostReport,
        AnyChange   = Change
    };

    /** Report the occurrence of an Event at the Study's "current state". **/
    void report(const Study&            study,
                const Event&            event,
                const EventTriggers&    triggers) const
    {
        ++m_numReports;
        reportVirtual(study, event, triggers);
    }
    
    /** Handle the occurrence of an Event by making a change to the 
    Study's "internal state". **/
    void change(Study&                  study,
                const Event&            event,
                const EventTriggers&    triggers,
                EventChangeResult&      result) const
    {
        ++m_numChanges;
        changeVirtual(study,event,triggers,result);
    }

protected:
    /** The constructor is for use by derived classes only. **/
    explicit EventAction(unsigned actionMask) 
    :   m_actionMask(actionMask), m_allowInterpolationForReport(true) 
    {   resetStats(); }

    /** Derived classes must implement this method, like this:
    @code
    class MyAction : public EventAction {
        MyAction* cloneVirtual() const override {
            return new MyAction(*this);
        }
    };
    @endcode
    **/
    virtual EventAction* cloneVirtual() const = 0;

    /** Derived classes must implement this if a Report action is selected. 
    Use `study.getCurrentState()` to obtain the State at which the report should
    occur, and `study.getSystem()` to access the System. **/
    virtual void reportVirtual
                   (const Study&            study,
                    const Event&            event,
                    const EventTriggers&    triggers) const
    {
        SimTK_THROW2(Exception::UnimplementedVirtualMethod, "EventAction", 
                     "reportVirtual");
    }

    /** Derived classes must implement this if a Change action is selected. 
    Use `study.updInternalState()` to obtain a writable reference to the 
    State that needs to be modified, and `study.getSystem()` to access 
    the System. **/
    virtual void changeVirtual
                   (Study&                  study,
                    const Event&            event,
                    const EventTriggers&    triggers,
                    EventChangeResult&      result) const
    {
        SimTK_THROW2(Exception::UnimplementedVirtualMethod, "EventAction", 
                     "changeVirtual");
    }

private:
    unsigned            m_actionMask;
    bool                m_allowInterpolationForReport;

    mutable long long   m_numReports;
    mutable long long   m_numChanges; 

    void resetStats() {m_numReports=m_numChanges=0;}
};


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_ACTION_H_
