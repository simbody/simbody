#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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


#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"

namespace SimTK {

class GeneralContactSubsystem;
class HuntCrossleyForceImpl;

class SimTK_SIMBODY_EXPORT HuntCrossleyForce : public Force {
public:
    /**
     * Create a Hunt-Crossley contact model.
     * 
     * @param contacts       the subsystem to which this contact model should be applied
     * @param contactSet     the index of the contact set to which this contact model will be applied
     */
    HuntCrossleyForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex contactSet);
    void setBodyParameters(int bodyIndex, Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction);
    Real getTransitionVelocity() const;
    void setTransitionVelocity(Real v);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(HuntCrossleyForce, HuntCrossleyForceImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_H_
