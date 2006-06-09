#ifndef SimTK_SIMBODY_FORCE_SUBSYSTEM_H_
#define SimTK_SIMBODY_FORCE_SUBSYSTEM_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Subsystem.h"

namespace SimTK {

class MatterSubsystem;

/**
 * This is logically an abstract class.
 */
class SimTK_SIMBODY_API ForceSubsystem : public Subsystem {
public:
    ForceSubsystem() { }

    /// This is a PIMPL virtual method. This is a Configured stage operator.
    Real calcPotentialEnergy(const State&) const;  // =0

    /// This is a PIMPL virtual method. This is a Dynamics stage operator.
    void addInForces(const State&, const MatterSubsystem&,
                     Vector_<SpatialVec>& rigidBodyForces,
                     Vector_<Vec3>&       particleForces,
                     Vector&              mobilityForces) const; // =0

    void setMatterSubsystemIndex(int subsys);
    int  getMatterSubsystemIndex() const;

    SimTK_PIMPL_DOWNCAST(ForceSubsystem, Subsystem);
    class ForceSubsystemRep& updRep();
    const ForceSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_SUBSYSTEM_H_
