#ifndef SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H
#define SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-6 Stanford University and the Authors.         *
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

/** @file
 * Define the extendible library-side implementation of the ForceSubsystem.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SubsystemGuts.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/// Public declaration of internals for ForceSubsystem extension
class ForceSubsystem::Guts : public Subsystem::Guts {
public:
    Guts(const String& name, const String& version) 
      : Subsystem::Guts(name,version)
    {
    }

    // Make sure the virtual destructor in Subsystem::Guts remains
    // virtual in this intermediate class.
    virtual ~Guts() { }

    // All the other Subsystem::Guts virtuals remain unresolved.

    // Return the MultibodySystem which owns this ForceSubsystem.
    const MultibodySystem& getMultibodySystem() const {
        return MultibodySystem::downcast(getSystem());
    }

    SimTK_DOWNCAST(ForceSubsystem::Guts, Subsystem::Guts);
};

// typedef ForceSubsystem::Guts ForceSubsystemRep;

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_SUBSYSTEM_GUTS_H
