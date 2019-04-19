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
 * Contributors: Antoine Falisse, Gil Serrancoli                              *
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

#include "HuntCrossleyForceImpl_smooth.h"

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(HuntCrossleyForce_smooth,
    HuntCrossleyForceImpl_smooth, Force);

HuntCrossleyForce_smooth::HuntCrossleyForce_smooth
    (GeneralForceSubsystem& forces) :
    Force(new HuntCrossleyForceImpl_smooth(forces)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

HuntCrossleyForceImpl_smooth::HuntCrossleyForceImpl_smooth
    (GeneralForceSubsystem& subsystem) :
    subsystem(subsystem){
}

void HuntCrossleyForceImpl_smooth::realizeTopology(State& state) const {
    energyCacheIndex=state.allocateCacheEntry(subsystem.getMySubsystemIndex(),
        Stage::Dynamics, new Value<Real>());
}

void HuntCrossleyForce_smooth::setParameters
    (Real stiffness, Real dissipation, Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity) {
    updImpl().setParameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity);
}

void HuntCrossleyForceImpl_smooth::setParameters
    (Real stiffness, Real dissipation,  Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity) {
    updParameters() = Parameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity);
}

void HuntCrossleyForce_smooth::setStiffness(Real stiffness) {
    updImpl().parameters.stiffness = stiffness;
}

void HuntCrossleyForceImpl_smooth::setStiffness(Real stiffness) {
    parameters.stiffness = stiffness;
}

void HuntCrossleyForce_smooth::setDissipation(Real dissipation) {
    updImpl().parameters.dissipation = dissipation;
}

void HuntCrossleyForceImpl_smooth::setDissipation(Real dissipation) {
    parameters.dissipation = dissipation;
}

void HuntCrossleyForce_smooth::setStaticFriction(Real staticFriction) {
    updImpl().parameters.staticFriction = staticFriction;
}

void HuntCrossleyForceImpl_smooth::setStaticFriction(Real staticFriction) {
    parameters.staticFriction = staticFriction;
}

void HuntCrossleyForce_smooth::setDynamicFriction(Real dynamicFriction) {
    updImpl().parameters.dynamicFriction = dynamicFriction;
}

void HuntCrossleyForceImpl_smooth::setDynamicFriction(Real dynamicFriction) {
    parameters.dynamicFriction = dynamicFriction;
}

void HuntCrossleyForce_smooth::setViscousFriction(Real viscousFriction) {
    updImpl().parameters.viscousFriction = viscousFriction;
}

void HuntCrossleyForceImpl_smooth::setViscousFriction(Real viscousFriction) {
    parameters.viscousFriction = viscousFriction;
}

void HuntCrossleyForce_smooth::setTransitionVelocity(Real transitionVelocity) {
    updImpl().parameters.transitionVelocity = transitionVelocity;
}

void HuntCrossleyForceImpl_smooth::setTransitionVelocity
    (Real transitionVelocity) {
    parameters.transitionVelocity = transitionVelocity;
}

void HuntCrossleyForce_smooth::setContactPlane(Vec3 normal, Real offset) {
    updImpl().setContactPlane(normal, offset);
}

void HuntCrossleyForceImpl_smooth::setContactPlane(Vec3 normal, Real offset) {
    ContactPlane = Plane(normal, offset);
}

void HuntCrossleyForce_smooth::setContactSphere(MobilizedBody bodyInput) {
    updImpl().BodySphere = bodyInput;
}

void HuntCrossleyForceImpl_smooth::setContactSphere(MobilizedBody bodyInput) {
    BodySphere = bodyInput;
}

void HuntCrossleyForce_smooth::setLocContactSphere(Vec3 LocContactSphere) {
    updImpl().LocContactSphere = LocContactSphere;
}

void HuntCrossleyForceImpl_smooth::setLocContactSphere(Vec3 LocContactSphere) {
    LocContactSphere = LocContactSphere;
}

void HuntCrossleyForce_smooth::setRadiusContactSphere(Real radius) {
    updImpl().RadiusContactSphere = radius;
}

void HuntCrossleyForceImpl_smooth::setRadiusContactSphere(Real radius) {
    RadiusContactSphere = radius;
}

MobilizedBody HuntCrossleyForce_smooth::getBodySphere() {
    return updImpl().BodySphere;
}

MobilizedBody HuntCrossleyForceImpl_smooth::getBodySphere() {
    return BodySphere;
}

Vec3 HuntCrossleyForce_smooth::getLocContactSphere() {
    return updImpl().LocContactSphere;
}

Vec3 HuntCrossleyForceImpl_smooth::getLocContactSphere() {
    return LocContactSphere;
}

Real HuntCrossleyForce_smooth::setRadiusContactSphere() {
    return updImpl().RadiusContactSphere;
}

Real HuntCrossleyForceImpl_smooth::getRadiusContactSphere() {
    return RadiusContactSphere;
}

const HuntCrossleyForceImpl_smooth::Parameters& HuntCrossleyForceImpl_smooth::
    getParameters() const {
    return parameters;
}

HuntCrossleyForceImpl_smooth::Parameters& HuntCrossleyForceImpl_smooth::
    updParameters() {
    return parameters;
}

void HuntCrossleyForceImpl_smooth::getContactPointSphere(const State& state,
    Vec3& contactPointPos) const {
    Vec3 posSphereInGround =
        BodySphere.findStationLocationInGround(state, LocContactSphere);
    contactPointPos = posSphereInGround -
        RadiusContactSphere*ContactPlane.getNormal();
}

Vec3 HuntCrossleyForce_smooth::getContactPointInBody(const State& state) {
    return updImpl().getContactPointInBody(state);
}

Vec3 HuntCrossleyForceImpl_smooth::getContactPointInBody(const State& state) {
    Vec3 posSphereInGround =
        BodySphere.findStationLocationInGround(state, LocContactSphere);
    Vec3 contactPointPos = posSphereInGround - RadiusContactSphere*ContactPlane.getNormal();
    return BodySphere.findStationAtGroundPoint(state, contactPointPos);
}

void HuntCrossleyForceImpl_smooth::calcForce(const State& state,
    Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces,
    Vector& mobilityForces) const {
    // Calculate the indentation based on the contact point location.
    Vec3 contactPointPos;
    getContactPointSphere(state, contactPointPos);
    const Real Indentation = - ContactPlane.getDistance(contactPointPos);
    // Adjust the contact location based on the relative stiffness of the two
    // materials. Here we assume, as in the original Simbody Hunt-Crossley
    // contact model, that both materials have the same relative stiffness.
    // As described in Sherman(2011), the point of contact will then be
    // located midway between the two surfaces. We therefore need to add half
    // the indentation to the contact location that was determined as the
    // location of the contact sphere center minus its radius.
    const Vec3 normal = ContactPlane.getNormal();
    const Vec3 contactPointPosAdj =
        contactPointPos+Real(1./2.)*Indentation*normal;
    const Vec3 contactPointPosAdjInB =
        BodySphere.findStationAtGroundPoint(state, contactPointPosAdj);
    // Calculate the contact point velocity.
    const Vec3 contactPointVel =
        BodySphere.findStationVelocityInGround(state, contactPointPosAdjInB);
    // Calculate the tangential and indentation velocities.
    const Vec3 v = contactPointVel;
    const Real vnormal = dot(v, normal);
    const Vec3 vtangent = v - vnormal*normal;
    const Real IndentationVel = -vnormal;
    // Get the contact model parameters.
    const Parameters parameters = getParameters();
    const Real stiffness = parameters.stiffness;
    const Real dissipation = parameters.dissipation;
    const Real vt = parameters.transitionVelocity;
    const Real us = parameters.staticFriction;
    const Real ud = parameters.dynamicFriction;
    const Real uv = parameters.viscousFriction;
    // Set the parameters for the smooth approximations.
    double eps = 1e-5;
    double bv = 50;
    double bd = 300;
    // Calculate the Hertz force.
    const Real k = (1./2.)*std::pow(stiffness, (2./3.));
    const Real fH = (4./3.)*k*std::sqrt(RadiusContactSphere*k)*
        std::pow(std::sqrt(Indentation*Indentation+eps),(3./2.));
    // Calculate the Hunt-Crossley force.
    const Real c = dissipation;
    const Real fHd = fH*(1.+(3./2.)*c*IndentationVel);
    const Real fn = fHd*(1./2.+(1./2.)*std::tanh(bd*Indentation))*
        (1./2.+(1./2.)*std::tanh(bv*(IndentationVel+(2./(3.*c)))));
    Vec3 force = fn*normal;
    // Calculate the friction force.
    const Real aux = pow(vtangent[0],2) +
        pow(vtangent[1],2)+pow(vtangent[2],2)+eps;
    const Real vslip = pow(aux,1./2.);
    const Real vrel = vslip / vt;
    const Real ffriction = fn*(std::min(vrel,Real(1))*
        (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
    force += ffriction*(-vtangent) / vslip;
    // Apply the force to the bodies.
    BodySphere.applyForceToBodyPoint(state, contactPointPosAdjInB,
        force, bodyForces);
}

// TODO
Real HuntCrossleyForceImpl_smooth::calcPotentialEnergy(const State& state)
    const {
    Real PotentialEnergy = 0;
    return PotentialEnergy;
}

} // namespace SimTK

