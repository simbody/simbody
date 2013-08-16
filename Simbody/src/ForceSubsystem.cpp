/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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
 * Implementation of ForceSubsystem, a still-abstract Subsystem.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include "simbody/internal/ForceSubsystemGuts.h"

namespace SimTK {

    /////////////////////
    // FORCE SUBSYSTEM //
    /////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

/*static*/ bool 
ForceSubsystem::isInstanceOf(const Subsystem& s) {
    return ForceSubsystemRep::isA(s.getSubsystemGuts());
}
/*static*/ const ForceSubsystem&
ForceSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<const ForceSubsystem&>(s);
}
/*static*/ ForceSubsystem&
ForceSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<ForceSubsystem&>(s);
}


const ForceSubsystemRep& 
ForceSubsystem::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const ForceSubsystemRep&>(getSubsystemGuts());
}
ForceSubsystemRep&       
ForceSubsystem::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<ForceSubsystemRep&>(updSubsystemGuts());
}

} // namespace SimTK

