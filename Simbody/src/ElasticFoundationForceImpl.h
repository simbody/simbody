#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_FORCE_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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
