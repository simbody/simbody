/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-19 Stanford University and the Authors.        *
 * Authors: Antoine Falisse, Gil Serrancoli                                   *
 * Contributors: Peter Eastman                                                *
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
#include "simbody/internal/MobilizedBody.h"

#include "SmoothSphereHalfplaneForceImpl.h"

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(SmoothSphereHalfplaneForce,
    SmoothSphereHalfplaneForceImpl, Force);

SmoothSphereHalfplaneForce::SmoothSphereHalfplaneForce
    (GeneralForceSubsystem& forces) :
    Force(new SmoothSphereHalfplaneForceImpl(forces)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

SmoothSphereHalfplaneForceImpl::SmoothSphereHalfplaneForceImpl
    (GeneralForceSubsystem& subsystem) :
    subsystem(subsystem){
}

void SmoothSphereHalfplaneForceImpl::realizeTopology(State& state) const {
    energyCacheIndex=state.allocateCacheEntry(subsystem.getMySubsystemIndex(),
        Stage::Dynamics, new Value<Real>());
}

void SmoothSphereHalfplaneForce::setParameters
    (Real stiffness, Real dissipation, Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updImpl().setParameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfplaneForceImpl::setParameters
    (Real stiffness, Real dissipation,  Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updParameters() = Parameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfplaneForce::setStiffness(Real stiffness) {
    updImpl().parameters.stiffness = stiffness;
}

void SmoothSphereHalfplaneForceImpl::setStiffness(Real stiffness) {
    parameters.stiffness = stiffness;
}

void SmoothSphereHalfplaneForce::setDissipation(Real dissipation) {
    updImpl().parameters.dissipation = dissipation;
}

void SmoothSphereHalfplaneForceImpl::setDissipation(Real dissipation) {
    parameters.dissipation = dissipation;
}

void SmoothSphereHalfplaneForce::setStaticFriction(Real staticFriction) {
    updImpl().parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfplaneForceImpl::setStaticFriction(Real staticFriction) {
    parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfplaneForce::setDynamicFriction(Real dynamicFriction) {
    updImpl().parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfplaneForceImpl::setDynamicFriction(Real dynamicFriction) {
    parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfplaneForce::setViscousFriction(Real viscousFriction) {
    updImpl().parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfplaneForceImpl::setViscousFriction(Real viscousFriction) {
    parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfplaneForce::setTransitionVelocity(
    Real transitionVelocity) {
    updImpl().parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfplaneForceImpl::setTransitionVelocity
    (Real transitionVelocity) {
    parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfplaneForce::setConstantContactForce(Real cf) {
    updImpl().parameters.cf = cf;
}

void SmoothSphereHalfplaneForceImpl::setConstantContactForce(Real cf) {
    parameters.cf = cf;
}

void SmoothSphereHalfplaneForce::setParameterTanhHertzForce(Real bd) {
    updImpl().parameters.bd = bd;
}

void SmoothSphereHalfplaneForceImpl::setParameterTanhHertzForce(Real bd) {
    parameters.bd = bd;
}

void SmoothSphereHalfplaneForce::setParameterTanhHuntCrossleyForce(Real bv) {
    updImpl().parameters.bv = bv;
}

void SmoothSphereHalfplaneForceImpl::setParameterTanhHuntCrossleyForce(
    Real bv) {
    parameters.bv = bv;
}

void SmoothSphereHalfplaneForce::setContactSphereInBody(
    MobilizedBody bodyInput1) {
    updImpl().bodySphere = bodyInput1;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereInBody(
    MobilizedBody bodyInput1){
    bodySphere = bodyInput1;
}

void SmoothSphereHalfplaneForce::setContactPlaneInBody(
    MobilizedBody bodyInput2) {
    updImpl().bodyPlane = bodyInput2;
}

void SmoothSphereHalfplaneForceImpl::setContactPlaneInBody(
    MobilizedBody bodyInput2){
    bodyPlane = bodyInput2;
}

void SmoothSphereHalfplaneForce::setContactSphereLocationInBody(
    Vec3 contactSphereLocation) {
    updImpl().contactSphereLocation = contactSphereLocation;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereLocationInBody(
    Vec3 contactSphereLocation) {
    contactSphereLocation = contactSphereLocation;
}

void SmoothSphereHalfplaneForce::setContactPlaneFrame(
    Transform planeFrame) {
    updImpl().contactPlaneFrame = planeFrame;
}

void SmoothSphereHalfplaneForceImpl::setContactPlaneFrame(
    Transform planeFrame) {
    contactPlaneFrame = planeFrame;
}

void SmoothSphereHalfplaneForce::setContactSphereRadius(Real radius) {
    updImpl().contactSphereRadius = radius;
}

void SmoothSphereHalfplaneForceImpl::setContactSphereRadius(Real radius) {
    contactSphereRadius = radius;
}

MobilizedBody SmoothSphereHalfplaneForce::getBodySphere() {
    return updImpl().bodySphere;
}

MobilizedBody SmoothSphereHalfplaneForceImpl::getBodySphere() {
    return bodySphere;
}

MobilizedBody SmoothSphereHalfplaneForce::getBodyPlane() {
    return updImpl().bodyPlane;
}

MobilizedBody SmoothSphereHalfplaneForceImpl::getBodyPlane() {
    return bodyPlane;
}

Vec3 SmoothSphereHalfplaneForce::getContactSphereLocationInBody() {
    return updImpl().contactSphereLocation;
}

Vec3 SmoothSphereHalfplaneForceImpl::getContactSphereLocationInBody() {
    return contactSphereLocation;
}

Real SmoothSphereHalfplaneForce::getContactSphereRadius() {
    return updImpl().contactSphereRadius;
}

Real SmoothSphereHalfplaneForceImpl::getContactSphereRadius() {
    return contactSphereRadius;
}

Transform SmoothSphereHalfplaneForce::getContactPlaneTransform() {
    return updImpl().contactPlaneFrame;
}

Transform SmoothSphereHalfplaneForceImpl::getContactPlaneTransform() {
    return contactPlaneFrame;
}

const SmoothSphereHalfplaneForceImpl::Parameters&
    SmoothSphereHalfplaneForceImpl::getParameters() const {
    return parameters;
}

SmoothSphereHalfplaneForceImpl::Parameters& SmoothSphereHalfplaneForceImpl::
    updParameters() {
    return parameters;
}

void SmoothSphereHalfplaneForceImpl::getNormalContactPlane(const State& state,
    UnitVec3& normalContactPlane) const {
    // TODO check for minus and for component. Seems weird wrt existing code
    normalContactPlane =
        -(bodyPlane.getBodyRotation(state)*contactPlaneFrame.y());
}

void SmoothSphereHalfplaneForceImpl::getContactSphereOrigin(const State& state,
    Vec3& contactSphereOrigin) const {
    contactSphereOrigin =
        bodySphere.getBodyTransform(state) * contactSphereLocation;
}

void SmoothSphereHalfplaneForceImpl::calcForce(const State& state,
    Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces,
    Vector& mobilityForces) const {
    // Calculate the indentation based on the contact point location.
    Vec3 contactSphereOrigin;
    getContactSphereOrigin(state, contactSphereOrigin);
    UnitVec3 normal;
    getNormalContactPlane(state, normal);
    const Vec3 contactPointPosition = contactSphereOrigin -
        contactSphereRadius*normal;
    // TODO not sure it is correct see SphereOnPlaneContact::findSeparation
    const Vec3 contactSphereLocationInPlane =
        bodySphere.findStationLocationInAnotherBody(
            state,contactSphereLocation,bodyPlane);
    const Vec3 p_PO_F = contactSphereLocationInPlane - contactPlaneFrame.p();
    const Real indentation =
        -(dot(p_PO_F, contactPlaneFrame.y()) - contactSphereRadius);
    // Initialize the potential energy.
    Real& pe = Value<Real>::updDowncast(state.updCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).upd();
    pe = 0.0;
    // Adjust the contact location based on the relative stiffness of the two
    // materials. Here we assume, as in the original Simbody Hunt-Crossley
    // contact model, that both materials have the same relative stiffness.
    // As described in Sherman(2011), the point of contact will then be
    // located midway between the two surfaces. We therefore need to add half
    // the indentation to the contact location that was determined as the
    // location of the contact sphere center minus its radius.
    ////////////////const Vec3 normal = contactPlane.getNormal(); // TODO
    const Vec3 contactPointPositionSphereAdjustedInGround =
        contactPointPosition+Real(1./2.)*indentation*normal;
    // Calculate the contact point velocity.
    const Vec3 station1 = bodySphere.findStationAtGroundPoint(state,
        contactPointPositionSphereAdjustedInGround);
    const Vec3 station2 = bodyPlane.findStationAtGroundPoint(state,
        contactPointPositionSphereAdjustedInGround);
    const Vec3 v1 = bodySphere.findStationVelocityInGround(state, station1);
    const Vec3 v2 = bodyPlane.findStationVelocityInGround(state, station2);
    const Vec3 v = v1-v2;
    // Calculate the normal and tangential velocities.
    // TODO is the sign correct?
    const Real vnormal = dot(v, normal);
    const Vec3 vtangent = v - vnormal*normal;
    // Get the contact model parameters.
    const Parameters parameters = getParameters();
    const Real stiffness = parameters.stiffness;
    const Real dissipation = parameters.dissipation;
    const Real vt = parameters.transitionVelocity;
    const Real us = parameters.staticFriction;
    const Real ud = parameters.dynamicFriction;
    const Real uv = parameters.viscousFriction;
    const Real cf = parameters.cf;
    const Real bd = parameters.bd;
    const Real bv = parameters.bv;
    // Calculate the Hertz force.
    const Real k = (1./2.)*std::pow(stiffness, (2./3.));
    const Real fH = (4./3.)*k*std::sqrt(contactSphereRadius*k)*
        std::pow(std::sqrt(indentation*indentation+cf),(3./2.));
    pe += Real(2./5.)*fH*indentation;
    // Calculate the Hunt-Crossley force.
    const Real c = dissipation;
    const Real fHd = fH*(1.+(3./2.)*c*vnormal);
    const Real fn = fHd*(1./2.+(1./2.)*std::tanh(bd*indentation))*
        (1./2.+(1./2.)*std::tanh(bv*(vnormal+(2./(3.*c)))));
    Vec3 force = fn*normal;
    // Calculate the friction force.
    const Real aux = vtangent.normSqr()+cf;
    const Real vslip = pow(aux,1./2.);
    const Real vrel = vslip / vt;
    const Real ffriction = fn*(std::min(vrel,Real(1))*
        (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
    force += ffriction*(vtangent) / vslip;
    // Apply the force to the bodies.
    bodySphere.applyForceToBodyPoint(state, station1, -force, bodyForces);
    // TODO is this necessary?
    bodyPlane.applyForceToBodyPoint(state, station2, force, bodyForces);
}

Real SmoothSphereHalfplaneForceImpl::calcPotentialEnergy(const State& state)
    const { return Value<Real>::downcast(state.getCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).get();
}

} // namespace SimTK

