#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
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
#include "simbody/internal/HuntCrossleyForce.h"
#include "ForceImpl.h"

namespace SimTK {

class HuntCrossleyForceImpl : public ForceImpl {
public:
    class Parameters;
    HuntCrossleyForceImpl(GeneralContactSubsystem& subystem, ContactSetIndex set);
    HuntCrossleyForceImpl* clone() const {
        return new HuntCrossleyForceImpl(*this);
    }
    void setBodyParameters(int bodyIndex, Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction);
    const Parameters& getParameters(int bodyIndex) const;
    Parameters& updParameters(int bodyIndex);
    Real getTransitionVelocity() const;
    void setTransitionVelocity(Real v);
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    void realizeTopology(State& state) const;
private:
    const GeneralContactSubsystem& subsystem;
    const ContactSetIndex set;
    std::vector<Parameters> parameters;
    Real transitionVelocity;
    mutable int energyCacheIndex;
};

class HuntCrossleyForceImpl::Parameters {
public:
    Parameters() : stiffness(1), dissipation(0), staticFriction(0), dynamicFriction(0), viscousFriction(0) {
    }
    Parameters(Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction) :
            stiffness(std::pow(stiffness, 2.0/3.0)), dissipation(dissipation), staticFriction(staticFriction), dynamicFriction(dynamicFriction), viscousFriction(viscousFriction) {
    }
    Real stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
