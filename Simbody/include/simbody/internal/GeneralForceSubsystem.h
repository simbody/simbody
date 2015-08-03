#ifndef SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
#define SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_

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

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MultibodySystem;
class SimbodyMatterSubsystem;
class Force;

/** This is a concrete subsystem which can apply arbitrary forces to a 
MultibodySystem. Each force element is represented by a Force object. For 
example, to add a spring between two bodies, you would write
@code
    GeneralForceSubsystem forces(system);
    // ...
    Force::TwoPointLinearSpring(forces, body1, station1, body2, station2, k, x0);
@endcode
**/
class SimTK_SIMBODY_EXPORT GeneralForceSubsystem : public ForceSubsystem {
public:
    GeneralForceSubsystem();
    explicit GeneralForceSubsystem(MultibodySystem&);

    /** Attach a new force to this subsystem. The subsystem takes over 
    ownership of the force, leaving the passed in handle as a reference to 
    it. This is normally called by the Force handle constructor.
    @param      force
        A writable reference to a Force handle whose referenced force element 
        is not yet owned by any subsystem. We will take over ownership of the 
        ForceImpl implementation objected referenced by the handle, bumping the
        reference count and leaving the reference in place so that the original
        handle can continue to be used to reference and modify the force 
        element.
    @return A unique integer index for the adopted force element that can be
    used to retrieve this force element from the subsystem later if needed. **/
    ForceIndex adoptForce(Force& force);
    
    /** Get the number of force elements which have been added to this 
    Subsystem. Legal ForceIndex values range from 0 to getNumForces()-1. **/
    int getNumForces() const;

    /** Get a const reference to a force element by index. **/
    const Force& getForce(ForceIndex index) const;

    /** Get a writable reference to a force element by index. **/
    Force& updForce(ForceIndex index);

    /** Find out whether a particular force element contained in this
    subsystem is currently disabled in the given state. **/
    bool isForceDisabled(const State& state, ForceIndex index) const;
    
    /** Disable or enable a particular force element contained in this 
    subsystem. This can usually be done more conveniently through the Force
    handle's disable() and enable() methods. Note that although force elements
    are normally enabled when created, it is possible that the force element
    will have been constructed to be disabled by default in which case it must
    be explicitly enabled. **/
    void setForceIsDisabled
       (State& state, ForceIndex index, bool shouldBeDisabled) const;

    /** Every Subsystem is owned by a System; a GeneralForceSubsystem expects
    to be owned by a MultibodySystem. This method returns a const reference
    to the containing MultibodySystem and will throw an exception if there is
    no containing System or it is not a MultibodySystem. **/
    const MultibodySystem& getMultibodySystem() const;

	/** Set the number of threads to distribute dynamics calculations among. 
	Defaults to number of cores on machine. **/
	void setNumberOfThreads(int nThreads) const;

    /** @cond **/   // don't show in Doxygen docs
    SimTK_PIMPL_DOWNCAST(GeneralForceSubsystem, ForceSubsystem);
    /** @endcond **/
private:
    // OBSOLETE; use getNumForces() instead.
    int getNForces() const {return getNumForces();}

    class GeneralForceSubsystemRep& updRep();
    const GeneralForceSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
