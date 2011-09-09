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
#include "simbody/internal/ElasticFoundationForce.h"
#include "ForceImpl.h"

namespace SimTK {

class ElasticFoundationForceImpl : public ForceImpl {
public:
    class Parameters;
    ElasticFoundationForceImpl(GeneralContactSubsystem& subystem, 
                               ContactSetIndex set);
    ElasticFoundationForceImpl* clone() const {
        return new ElasticFoundationForceImpl(*this);
    }
    void setBodyParameters
       (ContactSurfaceIndex bodyIndex, Real stiffness, Real dissipation, 
        Real staticFriction, Real dynamicFriction, Real viscousFriction);
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    void realizeTopology(State& state) const;
    void processContact(const State& state, ContactSurfaceIndex meshIndex, 
                        ContactSurfaceIndex otherBodyIndex, 
                        const Parameters& param, 
                        const std::set<int>& insideFaces,
                        Real areaScale,
                        Vector_<SpatialVec>& bodyForces, Real& pe) const;
private:
    friend class ElasticFoundationForce;
    const GeneralContactSubsystem& subsystem;
    const ContactSetIndex set;
    std::map<ContactSurfaceIndex, Parameters> parameters;
    Real transitionVelocity;
    mutable CacheEntryIndex energyCacheIndex;
};

class ElasticFoundationForceImpl::Parameters {
public:
    Parameters() : stiffness(1), dissipation(0), staticFriction(0), dynamicFriction(0), viscousFriction(0) {
    }
    Parameters(Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction) :
            stiffness(stiffness), dissipation(dissipation), staticFriction(staticFriction), dynamicFriction(dynamicFriction), viscousFriction(viscousFriction) {
    }
    Real stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction;
    Array_<Vec3> springPosition;
    Array_<UnitVec3> springNormal;
    Array_<Real> springArea;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
