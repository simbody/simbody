#ifndef SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
#define SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MultibodySystem;
class SimbodyMatterSubsystem;
class Force;

/**
 * This is a concrete subsystem which can apply arbitrary forces to a MultibodySystem.
 * Each force is represented by a Force object.  For example, to add a spring between two
 * bodies, you would write
 * 
 * <pre>
 * GeneralForceSubsystem forces(system);
 * ...
 * Force::TwoPointLinearSpring(forces, body1, station1, body2, station2, k, x0);
 * </pre>
 */
class SimTK_SIMBODY_EXPORT GeneralForceSubsystem : public ForceSubsystem {
public:
    GeneralForceSubsystem();
    explicit GeneralForceSubsystem(MultibodySystem&);

    /// Attach a new force to this subsystem.  The subsystem takes over ownership of the force,
    /// leaving the passed in handle as a reference to it.
    ForceIndex adoptForce(Force& force);
    
    /// Get the number of Forces which have been added.
    int getNForces() const;
    
    /// Get a const reference to a Force by index.
    const Force& getForce(ForceIndex index) const;

    /// Get a modifiable reference to a Force by index.
    Force& updForce(ForceIndex index);
    
    SimTK_PIMPL_DOWNCAST(GeneralForceSubsystem, ForceSubsystem);
private:
    class GeneralForceSubsystemRep& updRep();
    const GeneralForceSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
