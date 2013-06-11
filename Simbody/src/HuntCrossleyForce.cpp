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

#include "SimTKmath.h"

#include "simbody/internal/common.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MobilizedBody.h"

#include "HuntCrossleyForceImpl.h"

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(HuntCrossleyForce, HuntCrossleyForceImpl, Force);

HuntCrossleyForce::HuntCrossleyForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex set) :
        Force(new HuntCrossleyForceImpl(contacts, set)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

void HuntCrossleyForce::setBodyParameters
   (ContactSurfaceIndex surfIndex, Real stiffness, Real dissipation, 
    Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    updImpl().setBodyParameters(surfIndex, stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction);
}

Real HuntCrossleyForce::getTransitionVelocity() const {
    return getImpl().getTransitionVelocity();
}

void HuntCrossleyForce::setTransitionVelocity(Real v) {
    updImpl().setTransitionVelocity(v);
}

ContactSetIndex HuntCrossleyForce::getContactSetIndex() const {
    return getImpl().getContactSetIndex();
}


HuntCrossleyForceImpl::HuntCrossleyForceImpl(GeneralContactSubsystem& subsystem, ContactSetIndex set) : 
        subsystem(subsystem), set(set), transitionVelocity(Real(0.01)) {
}

void HuntCrossleyForceImpl::setBodyParameters
   (ContactSurfaceIndex bodyIndex, Real stiffness, Real dissipation, 
    Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    updParameters(bodyIndex) = Parameters(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction);
    subsystem.invalidateSubsystemTopologyCache();
}

const HuntCrossleyForceImpl::Parameters& HuntCrossleyForceImpl::
getParameters(ContactSurfaceIndex bodyIndex) const {
    assert(bodyIndex >= 0 && bodyIndex < subsystem.getNumBodies(set));
    // This fills in the default values which the missing entries implicitly 
    // had already.
    if (bodyIndex >= parameters.size())
        const_cast<Array_<Parameters,ContactSurfaceIndex>&>(parameters)
            .resize(bodyIndex+1);
    return parameters[bodyIndex];
}

HuntCrossleyForceImpl::Parameters& HuntCrossleyForceImpl::
updParameters(ContactSurfaceIndex bodyIndex) {
    assert(bodyIndex >= 0 && bodyIndex < subsystem.getNumBodies(set));
    subsystem.invalidateSubsystemTopologyCache();
    if (bodyIndex >= (int) parameters.size())
        parameters.resize(bodyIndex+1);
    return parameters[bodyIndex];
}


Real HuntCrossleyForceImpl::getTransitionVelocity() const {
    return transitionVelocity;
}

void HuntCrossleyForceImpl::setTransitionVelocity(Real v) {
    transitionVelocity = v;
    subsystem.invalidateSubsystemTopologyCache();
}

void HuntCrossleyForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                                      Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Array_<Contact>& contacts = subsystem.getContacts(state, set);
    Real& pe = Value<Real>::downcast(state.updCacheEntry(subsystem.getMySubsystemIndex(), energyCacheIndex)).upd();
    pe = 0.0;
    for (int i = 0; i < (int) contacts.size(); i++) {
        if (!PointContact::isInstance(contacts[i]))
            continue;
        const PointContact& contact = static_cast<const PointContact&>(contacts[i]);
        const Parameters& param1 = getParameters(contact.getSurface1());
        const Parameters& param2 = getParameters(contact.getSurface2());
        
        // Adjust the contact location based on the relative stiffness of the two materials.
        
        const Real s1 = param2.stiffness/(param1.stiffness+param2.stiffness);
        const Real s2 = 1-s1;
        const Real depth = contact.getDepth();
        const Vec3& normal = contact.getNormal();
        const Vec3 location = contact.getLocation()+(depth*(Real(0.5)-s1))*normal;
        
        // Calculate the Hertz force.

        const Real k = param1.stiffness*s1;
        const Real c = param1.dissipation*s1 + param2.dissipation*s2;
        const Real radius = contact.getEffectiveRadiusOfCurvature();
        const Real fH = Real(4./3.)*k*depth*std::sqrt(radius*k*depth);
        pe += Real(2./5.)*fH*depth;
        
        // Calculate the relative velocity of the two bodies at the contact point.
        
        const MobilizedBody& body1 = subsystem.getBody(set, contact.getSurface1());
        const MobilizedBody& body2 = subsystem.getBody(set, contact.getSurface2());
        const Vec3 station1 = body1.findStationAtGroundPoint(state, location);
        const Vec3 station2 = body2.findStationAtGroundPoint(state, location);
        const Vec3 v1 = body1.findStationVelocityInGround(state, station1);
        const Vec3 v2 = body2.findStationVelocityInGround(state, station2);
        const Vec3 v = v1-v2;
        const Real vnormal = dot(v, normal);
        const Vec3 vtangent = v-vnormal*normal;
        
        // Calculate the Hunt-Crossley force.
        
        const Real f = fH*(1+Real(1.5)*c*vnormal);
        Vec3 force = (f > 0 ? f*normal : Vec3(0));
        
        // Calculate the friction force.
        
        const Real vslip = vtangent.norm();
        if (vslip != 0) {
            const bool hasStatic = (param1.staticFriction != 0 || param2.staticFriction != 0);
            const bool hasDynamic= (param1.dynamicFriction != 0 || param2.dynamicFriction != 0);
            const bool hasViscous = (param1.viscousFriction != 0 || param2.viscousFriction != 0);
            const Real us = hasStatic ? 2*param1.staticFriction*param2.staticFriction/(param1.staticFriction+param2.staticFriction) : 0;
            const Real ud = hasDynamic ? 2*param1.dynamicFriction*param2.dynamicFriction/(param1.dynamicFriction+param2.dynamicFriction) : 0;
            const Real uv = hasViscous ? 2*param1.viscousFriction*param2.viscousFriction/(param1.viscousFriction+param2.viscousFriction) : 0;
            const Real vrel = vslip/getTransitionVelocity();
            const Real ffriction = f*(std::min(vrel, Real(1))*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            force += ffriction*vtangent/vslip;
        }
        
        // Apply the force to the bodies.
        
        body1.applyForceToBodyPoint(state, station1, -force, bodyForces);
        body2.applyForceToBodyPoint(state, station2, force, bodyForces);
    }
}

Real HuntCrossleyForceImpl::calcPotentialEnergy(const State& state) const {
    return Value<Real>::downcast(state.getCacheEntry(subsystem.getMySubsystemIndex(), energyCacheIndex)).get();
}

void HuntCrossleyForceImpl::realizeTopology(State& state) const {
        energyCacheIndex = state.allocateCacheEntry(subsystem.getMySubsystemIndex(), Stage::Dynamics, new Value<Real>());
}

} // namespace SimTK

