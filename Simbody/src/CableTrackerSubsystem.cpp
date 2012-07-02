/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Ian Stavness, Andreas Scholz                     *
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

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/CableTrackerSubsystem.h"

#include "CablePath_Impl.h"
#include "CableTrackerSubsystem_Impl.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//==============================================================================
//                        CABLE TRACKER SUBSYSTEM
//==============================================================================

bool CableTrackerSubsystem::isInstanceOf(const Subsystem& s) {
    return Impl::isA(s.getSubsystemGuts());
}

const CableTrackerSubsystem& CableTrackerSubsystem::
downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const CableTrackerSubsystem&>(s);
}
CableTrackerSubsystem& CableTrackerSubsystem::
updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<CableTrackerSubsystem&>(s);
}

const CableTrackerSubsystem::Impl& CableTrackerSubsystem::
getImpl() const {
    return dynamic_cast<const Impl&>(getSubsystemGuts());
}
CableTrackerSubsystem::Impl& CableTrackerSubsystem::
updImpl() {
    return dynamic_cast<Impl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use
// except for making std::vectors, which require a default constructor to be 
// available.
CableTrackerSubsystem::CableTrackerSubsystem() 
{   adoptSubsystemGuts(new Impl()); }

CableTrackerSubsystem::CableTrackerSubsystem(MultibodySystem& mbs) 
{   adoptSubsystemGuts(new Impl());
    mbs.adoptSubsystem(*this); } // steal ownership

int CableTrackerSubsystem::getNumCablePaths() const
{   return getImpl().getNumCablePaths(); }

const CablePath& CableTrackerSubsystem::
getCablePath(CablePathIndex cableIx) const
{   return getImpl().getCablePath(cableIx); }

CablePath& CableTrackerSubsystem::
updCablePath(CablePathIndex cableIx)
{   return updImpl().updCablePath(cableIx); }

