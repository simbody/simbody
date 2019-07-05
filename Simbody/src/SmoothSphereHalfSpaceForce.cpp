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

#include "SmoothSphereHalfSpaceForceImpl.h"

namespace SimTK {

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(SmoothSphereHalfSpaceForce,
    SmoothSphereHalfSpaceForceImpl, Force);

SmoothSphereHalfSpaceForce::SmoothSphereHalfSpaceForce
    (GeneralForceSubsystem& forces) :
    Force(new SmoothSphereHalfSpaceForceImpl(forces)) {
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

SmoothSphereHalfSpaceForceImpl::SmoothSphereHalfSpaceForceImpl
    (GeneralForceSubsystem& subsystem) :
    subsystem(subsystem){
}

void SmoothSphereHalfSpaceForceImpl::realizeTopology(State& state) const {
    energyCacheIndex=state.allocateCacheEntry(subsystem.getMySubsystemIndex(),
        Stage::Dynamics, new Value<Real>());
}

void SmoothSphereHalfSpaceForce::setParameters
    (Real stiffness, Real dissipation, Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updImpl().setParameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfSpaceForceImpl::setParameters
    (Real stiffness, Real dissipation,  Real staticFriction,
    Real dynamicFriction, Real viscousFriction, Real transitionVelocity, Real
    cf, Real bd, Real bv) {
    updParameters() = Parameters(stiffness, dissipation, staticFriction,
        dynamicFriction, viscousFriction, transitionVelocity, cf, bd, bv);
}

void SmoothSphereHalfSpaceForce::setStiffness(Real stiffness) {
    updImpl().parameters.stiffness = stiffness;
}

void SmoothSphereHalfSpaceForceImpl::setStiffness(Real stiffness) {
    parameters.stiffness = stiffness;
}

void SmoothSphereHalfSpaceForce::setDissipation(Real dissipation) {
    updImpl().parameters.dissipation = dissipation;
}

void SmoothSphereHalfSpaceForceImpl::setDissipation(Real dissipation) {
    parameters.dissipation = dissipation;
}

void SmoothSphereHalfSpaceForce::setStaticFriction(Real staticFriction) {
    updImpl().parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfSpaceForceImpl::setStaticFriction(Real staticFriction) {
    parameters.staticFriction = staticFriction;
}

void SmoothSphereHalfSpaceForce::setDynamicFriction(Real dynamicFriction) {
    updImpl().parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfSpaceForceImpl::setDynamicFriction(Real dynamicFriction) {
    parameters.dynamicFriction = dynamicFriction;
}

void SmoothSphereHalfSpaceForce::setViscousFriction(Real viscousFriction) {
    updImpl().parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfSpaceForceImpl::setViscousFriction(Real viscousFriction) {
    parameters.viscousFriction = viscousFriction;
}

void SmoothSphereHalfSpaceForce::setTransitionVelocity(
    Real transitionVelocity) {
    updImpl().parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfSpaceForceImpl::setTransitionVelocity
    (Real transitionVelocity) {
    parameters.transitionVelocity = transitionVelocity;
}

void SmoothSphereHalfSpaceForce::setConstantContactForce(Real cf) {
    updImpl().parameters.cf = cf;
}

void SmoothSphereHalfSpaceForceImpl::setConstantContactForce(Real cf) {
    parameters.cf = cf;
}

void SmoothSphereHalfSpaceForce::setHertzSmoothing(Real bd) {
    updImpl().parameters.bd = bd;
}

void SmoothSphereHalfSpaceForceImpl::setHertzSmoothing(Real bd) {
    parameters.bd = bd;
}

void SmoothSphereHalfSpaceForce::setHuntCrossleySmoothing(Real bv) {
    updImpl().parameters.bv = bv;
}

void SmoothSphereHalfSpaceForceImpl::setHuntCrossleySmoothing(Real bv) {
    parameters.bv = bv;
}

void SmoothSphereHalfSpaceForce::setContactSphereBody(
    MobilizedBody bodyInput1) {
    updImpl().bodySphere = bodyInput1;
}

void SmoothSphereHalfSpaceForceImpl::setContactSphereBody(
    MobilizedBody bodyInput1){
    bodySphere = bodyInput1;
}

void SmoothSphereHalfSpaceForce::setContactHalfSpaceBody(
    MobilizedBody bodyInput2) {
    updImpl().bodyHalfSpace = bodyInput2;
}

void SmoothSphereHalfSpaceForceImpl::setContactHalfSpaceBody(
    MobilizedBody bodyInput2){
    bodyHalfSpace = bodyInput2;
}

void SmoothSphereHalfSpaceForce::setContactSphereLocationInBody(
    Vec3 contactSphereLocation) {
    updImpl().contactSphereLocation = contactSphereLocation;
}

void SmoothSphereHalfSpaceForceImpl::setContactSphereLocationInBody(
    Vec3 contactSphereLocation) {
    contactSphereLocation = contactSphereLocation;
}

void SmoothSphereHalfSpaceForce::setContactHalfSpaceFrame(
    Transform halfSpaceFrame) {
    updImpl().contactHalfSpaceFrame = halfSpaceFrame;
}

void SmoothSphereHalfSpaceForceImpl::setContactHalfSpaceFrame(
    Transform halfSpaceFrame) {
    contactHalfSpaceFrame = halfSpaceFrame;
}

void SmoothSphereHalfSpaceForce::setContactSphereRadius(Real radius) {
    updImpl().contactSphereRadius = radius;
}

void SmoothSphereHalfSpaceForceImpl::setContactSphereRadius(Real radius) {
    contactSphereRadius = radius;
}

MobilizedBody SmoothSphereHalfSpaceForce::getBodySphere() {
    return updImpl().bodySphere;
}

MobilizedBody SmoothSphereHalfSpaceForceImpl::getBodySphere() {
    return bodySphere;
}

MobilizedBody SmoothSphereHalfSpaceForce::getBodyHalfSpace() {
    return updImpl().bodyHalfSpace;
}

MobilizedBody SmoothSphereHalfSpaceForceImpl::getBodyHalfSpace() {
    return bodyHalfSpace;
}

Vec3 SmoothSphereHalfSpaceForce::getContactSphereLocationInBody() {
    return updImpl().contactSphereLocation;
}

Vec3 SmoothSphereHalfSpaceForceImpl::getContactSphereLocationInBody() {
    return contactSphereLocation;
}

Real SmoothSphereHalfSpaceForce::getContactSphereRadius() {
    return updImpl().contactSphereRadius;
}

Real SmoothSphereHalfSpaceForceImpl::getContactSphereRadius() {
    return contactSphereRadius;
}

Transform SmoothSphereHalfSpaceForce::getContactHalfSpaceTransform() {
    return updImpl().contactHalfSpaceFrame;
}

Transform SmoothSphereHalfSpaceForceImpl::getContactHalfSpaceTransform() {
    return contactHalfSpaceFrame;
}

const SmoothSphereHalfSpaceForceImpl::Parameters&
    SmoothSphereHalfSpaceForceImpl::getParameters() const {
    return parameters;
}

SmoothSphereHalfSpaceForceImpl::Parameters& SmoothSphereHalfSpaceForceImpl::
    updParameters() {
    return parameters;
}

void SmoothSphereHalfSpaceForceImpl::getNormalContactHalfSpace(
    const State& state, UnitVec3& normalContactHalfSpace) const {
    normalContactHalfSpace =
        (bodyHalfSpace.getBodyRotation(state)*contactHalfSpaceFrame.x());
}

void SmoothSphereHalfSpaceForceImpl::getContactSphereOrigin(const State& state,
    Vec3& contactSphereOrigin) const {
    contactSphereOrigin =
        bodySphere.getBodyTransform(state) * contactSphereLocation;
}

void SmoothSphereHalfSpaceForceImpl::calcForce(const State& state,
    Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces,
    Vector& mobilityForces) const {
    // Calculate the indentation.
    const Vec3 contactSphereLocationInBodyHalfSpace =
        bodySphere.findStationLocationInAnotherBody(
            state,contactSphereLocation,bodyHalfSpace);
    const Vec3 distanceSphereHalfSpaceInBodyHalfSpace =
        contactSphereLocationInBodyHalfSpace - contactHalfSpaceFrame.p();
    const Real indentation = -(dot(distanceSphereHalfSpaceInBodyHalfSpace,
        -contactHalfSpaceFrame.x()) - contactSphereRadius);

    // Initialize the potential energy.
    Real& pe = Value<Real>::updDowncast(state.updCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).upd();
    pe = 0.0;

    // Calculate the contact point location (in ground).
    Vec3 contactSphereOriginInGround;
    getContactSphereOrigin(state, contactSphereOriginInGround);
    UnitVec3 normal; // the normal is pointing in the direction of contact
    getNormalContactHalfSpace(state, normal);
    const Vec3 contactPointPositionInGround = contactSphereOriginInGround +
        contactSphereRadius*normal;
    // Adjust the contact location based on the relative stiffness of the two
    // materials. Here we assume, as in the original Simbody Hunt-Crossley
    // contact model, that both materials have the same relative stiffness.
    // As described in Sherman et al. (2011), the point of contact will then be
    // located midway between the two surfaces. We therefore need to subtract
    // (subtraction due to normal direction) half the indentation to the
    // contact location that was determined as the location of the contact
    // sphere center plus its radius (plus due to normal direction).
    const Vec3 contactPointPositionAdjustedInGround =
        contactPointPositionInGround-Real(1./2.)*indentation*normal;
    // Calculate the contact point velocity.
    const Vec3 station1 = bodySphere.findStationAtGroundPoint(state,
        contactPointPositionAdjustedInGround);
    const Vec3 station2 = bodyHalfSpace.findStationAtGroundPoint(state,
        contactPointPositionAdjustedInGround);
    const Vec3 v1 = bodySphere.findStationVelocityInGround(state, station1);
    const Vec3 v2 = bodyHalfSpace.findStationVelocityInGround(state, station2);
    const Vec3 v = v1-v2;
    // Calculate the normal and tangential velocities.
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
    const Real k = (1./2.)*std::pow(stiffness,2./3.);
    const Real fh_pos = (4./3.)*k*std::sqrt(contactSphereRadius*k)*
        std::pow(std::sqrt(indentation*indentation+cf),3./2.);
    const Real fh_smooth = fh_pos*(1./2.+(1./2.)*std::tanh(bd*indentation));
    // Calculate the potential energy.
    // The potential energy is the integral of the Hertz force. Due to the
    // smooth approximation, there is no exact expression for the potential
    // energy. Here we provide an approximation based on the original
    // expression (i.e., pe = Real(2./5.)*fHertz*indentation) where we replace
    // fHertz by the smooth approximation (i.e., fh_smooth).
    pe += Real(2./5.)*fh_smooth*indentation;
    // Calculate the Hunt-Crossley force.
    const Real c = dissipation;
    const Real fhc_pos = fh_smooth*(1.+(3./2.)*c*vnormal);
    const Real fhc_smooth = fhc_pos*(1./2.+(1./2.)*std::tanh(
        bv*(vnormal+(2./(3.*c)))));
    Vec3 force = fhc_smooth*normal;
    // Calculate the friction force.
    const Real aux = vtangent.normSqr()+cf;
    const Real vslip = pow(aux,1./2.);
    const Real vrel = vslip / vt;
    const Real ff = fhc_smooth*(std::min(vrel,Real(1))*
        (ud+2*(us-ud)/(1+vrel*vrel))+uv*vslip);
    force += ff*(vtangent) / vslip;
    // Apply the force to the bodies.
    bodySphere.applyForceToBodyPoint(state, station1, -force, bodyForces);
    bodyHalfSpace.applyForceToBodyPoint(state, station2, force, bodyForces);
}

Real SmoothSphereHalfSpaceForceImpl::calcPotentialEnergy(const State& state)
    const { return Value<Real>::downcast(state.getCacheEntry(
        subsystem.getMySubsystemIndex(), energyCacheIndex)).get();
}

} // namespace SimTK

