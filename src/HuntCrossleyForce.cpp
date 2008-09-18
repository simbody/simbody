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
#include "simbody/internal/Contact.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/MobilizedBody.h"

#include "HuntCrossleyForceImpl.h"

using std::vector;

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(HuntCrossleyForce, HuntCrossleyForceImpl, Force);

HuntCrossleyForce::HuntCrossleyForce(GeneralForceSubsystem& forces, GeneralContactSubsystem& contacts, ContactSetIndex set) :
        Force(new HuntCrossleyForceImpl(contacts, set)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

void HuntCrossleyForce::setBodyParameters(int bodyIndex, Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    updImpl().setBodyParameters(bodyIndex, stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction);
}

Real HuntCrossleyForce::getTransitionVelocity() const {
    return getImpl().getTransitionVelocity();
}

void HuntCrossleyForce::setTransitionVelocity(Real v) {
    updImpl().setTransitionVelocity(v);
}

HuntCrossleyForceImpl::HuntCrossleyForceImpl(GeneralContactSubsystem& subsystem, ContactSetIndex set) : 
        subsystem(subsystem), set(set), transitionVelocity(0.001) {
}

void HuntCrossleyForceImpl::setBodyParameters(int bodyIndex, Real stiffness, Real dissipation, Real staticFriction, Real dynamicFriction, Real viscousFriction) {
    updParameters(bodyIndex) = Parameters(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction);
    subsystem.invalidateSubsystemTopologyCache();
}

const HuntCrossleyForceImpl::Parameters& HuntCrossleyForceImpl::getParameters(int bodyIndex) const {
    assert(bodyIndex >= 0 && bodyIndex < subsystem.getNumBodies(set));
    if (bodyIndex >= parameters.size())
        const_cast<vector<Parameters>&>(parameters).resize(bodyIndex+1); // This fills in the default values which the missing entries implicitly had already.
    return parameters[bodyIndex];
}

HuntCrossleyForceImpl::Parameters& HuntCrossleyForceImpl::updParameters(int bodyIndex) {
    assert(bodyIndex >= 0 && bodyIndex < subsystem.getNumBodies(set));
    subsystem.invalidateSubsystemTopologyCache();
    if (bodyIndex >= parameters.size())
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

void HuntCrossleyForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const vector<Contact>& contacts = subsystem.getContacts(state, set);
    for (int i = 0; i < contacts.size(); i++) {
        if (!PointContact::isInstance(contacts[i]))
            continue;
        const PointContact& contact = static_cast<const PointContact&>(contacts[i]);
        const Parameters& param1 = getParameters(contact.getFirstBody());
        const Parameters& param2 = getParameters(contact.getSecondBody());
        
        // Adjust the contact location based on the relative stiffness of the two materials.
        
        const Real s1 = param2.stiffness/(param1.stiffness+param2.stiffness);
        const Real s2 = 1-s1;
        const Real depth = contact.getDepth();
        const Vec3& normal = contact.getNormal();
        const Vec3 location = contact.getLocation()+(depth*(0.5-s1))*normal;
        
        // Calculate the Hertz force.

        const Real k = param1.stiffness*s1;
        const Real c = param1.dissipation*s1 + param2.dissipation*s2;
        const Real radius = contact.getRadius();
        const Real curvature = radius*radius/depth;
        const Real fH = (4.0/3.0)*k*depth*std::sqrt(curvature*k*depth);
        
        // Calculate the relative velocity of the two bodies at the contact point.
        
        const MobilizedBody& body1 = subsystem.getBody(set, contact.getFirstBody());
        const MobilizedBody& body2 = subsystem.getBody(set, contact.getSecondBody());
        const Vec3 station1 = body1.findStationAtGroundPoint(state, location);
        const Vec3 station2 = body2.findStationAtGroundPoint(state, location);
        const Vec3 v1 = body1.findStationVelocityInGround(state, station1);
        const Vec3 v2 = body2.findStationVelocityInGround(state, station2);
        const Vec3 v = v1-v2;
        const Real vnormal = dot(v, normal);
        const Vec3 vtangent = v-vnormal*normal;
        
        // Calculate the Hunt-Crossley force.
        
        const Real f = fH*(1+1.5*c*vnormal);
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
            const Real ffriction = f*(std::min(vrel, 1.0)*(ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
            force += ffriction*vtangent/vslip;
        }
        
        // Apply the force to the bodies.
        
        body1.applyForceToBodyPoint(state, station1, -force, bodyForces);
        body2.applyForceToBodyPoint(state, station2, force, bodyForces);
    }
}

Real HuntCrossleyForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

} // namespace SimTK

