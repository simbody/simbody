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
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Force.h"

#include "ForceImpl.h"

namespace SimTK {

ForceIndex Force::getForceIndex() const {
    return getImpl().getForceIndex();
}

// TwoPointLinearSpring

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointLinearSpring, Force::TwoPointLinearSpringImpl, Force);

Force::TwoPointLinearSpring::TwoPointLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real k, Real x0) : Force(new TwoPointLinearSpringImpl(
        body1, station1, body2, station2, k, x0)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::TwoPointLinearSpringImpl::TwoPointLinearSpringImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real k, Real x0) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), k(k), x0(x0) {
}

void Force::TwoPointLinearSpringImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.T() + s2_G;

    const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
    const Real d   = r_G.norm();  // distance between the points
    const Real stretch   = d - x0; // + -> tension, - -> compression
    const Real frcScalar = k*stretch; // k(x-x0)

    const Vec3 f1_G = (frcScalar/d) * r_G;
    bodyForces[body1] +=  SpatialVec(s1_G % f1_G, f1_G);
    bodyForces[body2] -=  SpatialVec(s2_G % f1_G, f1_G);
}

Real Force::TwoPointLinearSpringImpl::calcPotentialEnergy(const State& state) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.T() + s2_G;

    const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
    const Real d   = r_G.norm();  // distance between the points
    const Real stretch   = d - x0; // + -> tension, - -> compression

    return 0.5*k*stretch*stretch; // 1/2 k (x-x0)^2
}

// TwoPointLinearDamper

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointLinearDamper, Force::TwoPointLinearDamperImpl, Force);

Force::TwoPointLinearDamper::TwoPointLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real damping) : Force(new TwoPointLinearDamperImpl(
        body1, station1, body2, station2, damping)) {
    assert(damping >= 0);
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::TwoPointLinearDamperImpl::TwoPointLinearDamperImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real damping) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), damping(damping) {
}

void Force::TwoPointLinearDamperImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.T() + s2_G;

    const Vec3 v1_G = matter.getMobilizedBody(body1).findStationVelocityInGround(state, station1);
    const Vec3 v2_G = matter.getMobilizedBody(body2).findStationVelocityInGround(state, station2);
    const Vec3 vRel = v2_G - v1_G; // relative velocity

    const UnitVec3 d(p2_G - p1_G); // direction from point1 to point2
    const Real frc = damping*dot(vRel, d); // c*v

    const Vec3 f1_G = frc*d;
    bodyForces[body1] +=  SpatialVec(s1_G % f1_G, f1_G);
    bodyForces[body2] -=  SpatialVec(s2_G % f1_G, f1_G);
}

Real Force::TwoPointLinearDamperImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// TwoPointConstantForce

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::TwoPointConstantForce, Force::TwoPointConstantForceImpl, Force);

Force::TwoPointConstantForce::TwoPointConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real force) : Force(new TwoPointConstantForceImpl(
        body1, station1, body2, station2, force)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::TwoPointConstantForceImpl::TwoPointConstantForceImpl(const MobilizedBody& body1, const Vec3& station1,
        const MobilizedBody& body2, const Vec3& station2, Real force) : matter(body1.getMatterSubsystem()),
        body1(body1.getMobilizedBodyIndex()), station1(station1),
        body2(body2.getMobilizedBodyIndex()), station2(station2), force(force) {
}

void Force::TwoPointConstantForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB1 = matter.getMobilizedBody(body1).getBodyTransform(state);
    const Transform& X_GB2 = matter.getMobilizedBody(body2).getBodyTransform(state);

    const Vec3 s1_G = X_GB1.R() * station1;
    const Vec3 s2_G = X_GB2.R() * station2;

    const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
    const Vec3 p2_G = X_GB2.T() + s2_G;

    const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
    const Real x   = r_G.norm();  // distance between the points
    const UnitVec3 d(r_G/x, true);

    const Vec3 f2_G = force * d;
    bodyForces[body1] -=  SpatialVec(s1_G % f2_G, f2_G);
    bodyForces[body2] +=  SpatialVec(s2_G % f2_G, f2_G);
}

Real Force::TwoPointConstantForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// MobilityLinearSpring

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityLinearSpring, Force::MobilityLinearSpringImpl, Force);

Force::MobilityLinearSpring::MobilityLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate,
        Real k, Real x0) : Force(new MobilityLinearSpringImpl(body, coordinate, k, x0)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::MobilityLinearSpringImpl::MobilityLinearSpringImpl(const MobilizedBody& body, int coordinate,
        Real k, Real x0) : matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), coordinate(coordinate), k(k), x0(x0) {
}

void Force::MobilityLinearSpringImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const MobilizedBody& mb = matter.getMobilizedBody(body);
    const Real q = mb.getOneQ(state, coordinate);
    const Real frc = -k*(q-x0);
    mb.applyOneMobilityForce(state, coordinate, frc, mobilityForces);
}

Real Force::MobilityLinearSpringImpl::calcPotentialEnergy(const State& state) const {
    const MobilizedBody& mb = matter.getMobilizedBody(body);
    const Real q = mb.getOneQ(state, coordinate);
    const Real frc = -k*(q-x0);
    return 0.5*k*(q-x0)*(q-x0);
}

// MobilityLinearDamper

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityLinearDamper, Force::MobilityLinearDamperImpl, Force);

Force::MobilityLinearDamper::MobilityLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate,
        Real damping) : Force(new MobilityLinearDamperImpl(body, coordinate, damping)) {
    assert(damping >= 0);
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::MobilityLinearDamperImpl::MobilityLinearDamperImpl(const MobilizedBody& body, int coordinate,
        Real damping) : matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), coordinate(coordinate), damping(damping) {
}

void Force::MobilityLinearDamperImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const MobilizedBody& mb = matter.getMobilizedBody(body);
    const Real u = mb.getOneU(state, coordinate);
    const Real frc = -damping*u;
    mb.applyOneMobilityForce(state, coordinate, frc, mobilityForces);
}

Real Force::MobilityLinearDamperImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// MobilityConstantForce

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::MobilityConstantForce, Force::MobilityConstantForceImpl, Force);

Force::MobilityConstantForce::MobilityConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate,
        Real force) : Force(new MobilityConstantForceImpl(body, coordinate, force)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::MobilityConstantForceImpl::MobilityConstantForceImpl(const MobilizedBody& body, int coordinate,
        Real force) : matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), coordinate(coordinate), force(force) {
}

void Force::MobilityConstantForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const MobilizedBody& mb = matter.getMobilizedBody(body);
    const Real q = mb.getOneQ(state, coordinate);
    mb.applyOneMobilityForce(state, coordinate, force, mobilityForces);
}

Real Force::MobilityConstantForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// ConstantForce

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::ConstantForce, Force::ConstantForceImpl, Force);

Force::ConstantForce::ConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& station, const Vec3& force) :
        Force(new ConstantForceImpl(body, station, force)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::ConstantForceImpl::ConstantForceImpl(const MobilizedBody& body, const Vec3& station, const Vec3& force) :
        matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), station(station), force(force) {
}

void Force::ConstantForceImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const Transform& X_GB = matter.getMobilizedBody(body).getBodyTransform(state);
    const Vec3 station_G = X_GB.R() * station;
    bodyForces[body] += SpatialVec(station_G % force, force);
}

Real Force::ConstantForceImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// ConstantTorque

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::ConstantTorque, Force::ConstantTorqueImpl, Force);

Force::ConstantTorque::ConstantTorque(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& torque) :
        Force(new ConstantTorqueImpl(body, torque)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::ConstantTorqueImpl::ConstantTorqueImpl(const MobilizedBody& body, const Vec3& torque) :
        matter(body.getMatterSubsystem()), body(body.getMobilizedBodyIndex()), torque(torque) {
}

void Force::ConstantTorqueImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    bodyForces[body][0] += torque;
}

Real Force::ConstantTorqueImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// GlobalDamper

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::GlobalDamper, Force::GlobalDamperImpl, Force);

Force::GlobalDamper::GlobalDamper(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter,
        Real damping) : Force(new GlobalDamperImpl(matter, damping)) {
    assert(damping >= 0);
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Force::GlobalDamperImpl::GlobalDamperImpl(const SimbodyMatterSubsystem& matter, Real damping) : matter(matter), damping(damping) {
}

void Force::GlobalDamperImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    mobilityForces -= damping*matter.getU(state);
}

Real Force::GlobalDamperImpl::calcPotentialEnergy(const State& state) const {
    return 0;
}

// UniformGravity

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::UniformGravity, Force::UniformGravityImpl, Force);

Force::UniformGravity::UniformGravity(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter,
        const Vec3& g, Real zeroHeight) : Force(new UniformGravityImpl(matter, g, zeroHeight)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

Vec3 Force::UniformGravity::getGravity() const {
    return getImpl().getGravity();
}

void Force::UniformGravity::setGravity(const Vec3& g) {
    updImpl().setGravity(g);
}

Real Force::UniformGravity::getZeroHeight() const {
    return getImpl().getZeroHeight();
}

void Force::UniformGravity::setZeroHeight(Real height) {
    updImpl().setZeroHeight(height);
}

Force::UniformGravityImpl::UniformGravityImpl(const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight) : matter(matter), g(g), zeroHeight(zeroHeight) {
}

void Force::UniformGravityImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    const int nBodies    = matter.getNBodies();
    const int nParticles = matter.getNParticles();

    if (nParticles) {
        const Vector& m = matter.getAllParticleMasses(state);
        const Vector_<Vec3>& loc_G = matter.getAllParticleLocations(state);
        for (int i=0; i < nParticles; ++i) {
            particleForces[i] += g * m[i];
        }
    }

    // no need to apply gravity to Ground!
    for (MobilizedBodyIndex i(1); i < nBodies; ++i) {
        const MassProperties& mprops = matter.getMobilizedBody(i).getBodyMassProperties(state);
        const Real&      m       = mprops.getMass();
        const Vec3&      com_B   = mprops.getMassCenter();
        const Transform& X_GB    = matter.getMobilizedBody(i).getBodyTransform(state);
        const Vec3       com_B_G = X_GB.R()*com_B;
        const Vec3       frc_G   = m*g;

        bodyForces[i] += SpatialVec(com_B_G % frc_G, frc_G); 
    }
}

Real Force::UniformGravityImpl::calcPotentialEnergy(const State& state) const {
    const int nBodies    = matter.getNBodies();
    const int nParticles = matter.getNParticles();
    Real pe = 0.0;

    if (nParticles) {
        const Vector& m = matter.getAllParticleMasses(state);
        const Vector_<Vec3>& loc_G = matter.getAllParticleLocations(state);
        for (int i=0; i < nParticles; ++i) {
            pe -= m[i]*(~g*loc_G[i] + zeroHeight); // odd signs because height is in -g direction
        }
    }

    // no need to apply gravity to Ground!
    for (MobilizedBodyIndex i(1); i < nBodies; ++i) {
        const MassProperties& mprops = matter.getMobilizedBody(i).getBodyMassProperties(state);
        const Real&      m       = mprops.getMass();
        const Vec3&      com_B   = mprops.getMassCenter();
        const Transform& X_GB    = matter.getMobilizedBody(i).getBodyTransform(state);
        const Vec3       com_B_G = X_GB.R()*com_B;
        const Vec3       com_G   = X_GB.T() + com_B_G;

        pe -= m*(~g*com_G + zeroHeight); // odd signs because height is in -g direction
    }
    return pe;
}

void Force::UniformGravityImpl::invalidateTopologyCache() {
    matter.invalidateSubsystemTopologyCache();
}

// Custom

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::Custom, Force::CustomImpl, Force);

Force::Custom::Custom(GeneralForceSubsystem& forces, Implementation* implementation) : 
        Force(new CustomImpl(implementation)) {
    updImpl().setForceIndex(forces.adoptForce(*this));
}

const Force::Custom::Implementation& Force::Custom::getImplementation() const {
    return getImpl().getImplementation();
}

Force::Custom::Implementation& Force::Custom::updImplementation() {
    return updImpl().updImplementation();    
}


Force::CustomImpl::CustomImpl(Force::Custom::Implementation* implementation) : implementation(implementation) {
}

void Force::CustomImpl::calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
    implementation->calcForce(state, bodyForces, particleForces, mobilityForces);
}

Real Force::CustomImpl::calcPotentialEnergy(const State& state) const {
    return implementation->calcPotentialEnergy(state);
}

} // namespace SimTK

