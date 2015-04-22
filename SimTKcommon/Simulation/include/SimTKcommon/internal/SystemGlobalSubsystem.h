#ifndef SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_
#define SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"

#include <cassert>

namespace SimTK {
class ScheduledEventHandler;
class ScheduledEventReporter;
class TriggeredEventHandler;
class TriggeredEventReporter;

//==============================================================================
//                         SYSTEM GLOBAL SUBSYSTEM
//==============================================================================
/** This Subsystem manages resources that are global to the System rather than
associated with any particular Subsystem. For example, we maintain lists of
event handlers and reporters here. There is one of these automatically defined
in every System after it is constructed.

To obtain the global subsystem for a System, call getSystemGlobalSubsystem() or 
getSystemGlobalSubsystem() on it. Also, a System can be implicitly converted
to a Subsystem, in which case it actually returns a reference to
this subsystem. 

@par Realization order
For the early Stage realizations, the %SystemGlobalSubsystem is realized 
*first* so that resources it allocates can be used at the same Stage by other
Subsystems. Early means Stage::Topology, Stage::Model, and Stage::Instance.
For the later stages, this subsystem is realized *last*, after all
other subsystems. That way if Event::Signal functions defined here to trigger
events depend on computations done in other subsystems, those computations
will have been completed by the time they are needed. Later here means
Stage::Time and above, through Stage::Report. **/
class SimTK_SimTKCOMMON_EXPORT SystemGlobalSubsystem : public Subsystem {
public:
    explicit SystemGlobalSubsystem(System& sys);
    void addEventHandler(ScheduledEventHandler* handler);
    void addEventHandler(TriggeredEventHandler* handler);
    void addEventReporter(ScheduledEventReporter* handler) const;
    void addEventReporter(TriggeredEventReporter* handler) const;

    /** @cond **/  // don't let doxygen see this private class
    class Guts;
    /** @endcond **/
private:
    const Guts& getGuts() const;
    Guts& updGuts();
};


//==============================================================================
//                                 SYSTEM
//==============================================================================
// These methods can't be defined until SystemGlobalSubsystem has been
// declared.

inline void System::addEventHandler(ScheduledEventHandler* handler)
{   updSystemGlobalSubsystem().addEventHandler(handler); }
inline void System::addEventHandler(TriggeredEventHandler* handler)
{   updSystemGlobalSubsystem().addEventHandler(handler); }
inline void System::addEventReporter(ScheduledEventReporter* handler) const
{   getSystemGlobalSubsystem().addEventReporter(handler); }
inline void System::addEventReporter(TriggeredEventReporter* handler) const
{   getSystemGlobalSubsystem().addEventReporter(handler); }

inline System::operator const Subsystem&() const 
{   return getSystemGlobalSubsystem(); }
inline System::operator Subsystem&() 
{   return updSystemGlobalSubsystem(); }

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GLOBAL_SUBSYSTEM_H_
