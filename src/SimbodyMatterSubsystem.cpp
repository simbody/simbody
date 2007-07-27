/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/**@file
 *
 * Implementation of SimbodyMatterSubsystem, a concrete 
 * MatterSubsystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/MobilizedBody.h"

#include "MobilizedBodyRep.h"
#include "SimbodyMatterSubsystemRep.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {


/*static*/ bool 
SimbodyMatterSubsystem::isInstanceOf(const Subsystem& s) {
    return SimbodyMatterSubsystemRep::isA(s.getSubsystemGuts());
}
/*static*/ const SimbodyMatterSubsystem&
SimbodyMatterSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const SimbodyMatterSubsystem&>(s);
}
/*static*/ SimbodyMatterSubsystem&
SimbodyMatterSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<SimbodyMatterSubsystem&>(s);
}

const SimbodyMatterSubsystemRep& 
SimbodyMatterSubsystem::getRep() const {
    return dynamic_cast<const SimbodyMatterSubsystemRep&>(getSubsystemGuts());
}
SimbodyMatterSubsystemRep&       
SimbodyMatterSubsystem::updRep() {
    return dynamic_cast<SimbodyMatterSubsystemRep&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
SimbodyMatterSubsystem::SimbodyMatterSubsystem() 
  : Subsystem()
{
    adoptSubsystemGuts(new SimbodyMatterSubsystemRep());
    updRep().createGroundBody(); //TODO: handle this differently
}

SimbodyMatterSubsystem::SimbodyMatterSubsystem(MultibodySystem& mbs) 
  : Subsystem()
{
    adoptSubsystemGuts(new SimbodyMatterSubsystemRep());
    updRep().createGroundBody(); //TODO: handle this differently
    mbs.setMatterSubsystem(*this);
}

MobilizedBodyId SimbodyMatterSubsystem::adoptMobilizedBody(MobilizedBodyId parent, MobilizedBody& child) {
    return updRep().adoptMobilizedBody(parent,child);
}
const MobilizedBody& SimbodyMatterSubsystem::getMobilizedBody(MobilizedBodyId id) const {
    return getRep().getMobilizedBody(id);
}
MobilizedBody& SimbodyMatterSubsystem::updMobilizedBody(MobilizedBodyId id) {
    return updRep().updMobilizedBody(id);
}
const MobilizedBody::Ground& SimbodyMatterSubsystem::getGround() const {
    return getRep().getGround();
}
MobilizedBody::Ground& SimbodyMatterSubsystem::updGround() {
    return updRep().updGround();
}

ConstraintId SimbodyMatterSubsystem::adoptConstraint(Constraint& child) {
    return updRep().adoptConstraint(child);
}
const Constraint& SimbodyMatterSubsystem::getConstraint(ConstraintId id) const {
    return getRep().getConstraint(id);
}
Constraint& SimbodyMatterSubsystem::updConstraint(ConstraintId id) {
    return updRep().updConstraint(id);
}

void SimbodyMatterSubsystem::calcAcceleration(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot,
    Vector_<SpatialVec>&       A_GB) const
{
    Vector_<Vec3> particleForces; // TODO
    SBAccelerationCache ac;
    ac.allocate(getRep().topologyCache);

    Vector udotErr(getNUDotErr(s)); // unwanted return value
    Vector multipliers(getNMultipliers(s)); // unwanted return value

    getRep().calcLoopForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                             ac, udot, multipliers, udotErr);

    A_GB = ac.bodyAccelerationInGround;
}


void SimbodyMatterSubsystem::calcAccelerationIgnoringConstraints(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot,
    Vector_<SpatialVec>&       A_GB) const
{
    Vector              netHingeForces; // unwanted side effect

    getRep().calcTreeAccelerations(s,mobilityForces,bodyForces,
        netHingeForces, A_GB, udot);
}


void SimbodyMatterSubsystem::calcMInverseV(const State& s,
    const Vector&        v,
    Vector&              MinvV,
    Vector_<SpatialVec>& A_GB) const
{
    getRep().calcMInverseF(s,v, A_GB, MinvV);
}

// Convert spatial forces to internal equivalent, ignoring velocity and
// constraints.
void SimbodyMatterSubsystem::calcInternalGradientFromSpatial(const State& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    getRep().calcInternalGradientFromSpatial(s,dEdR,dEdQ);
}

// Convert spatial forces and centrifugal forces to an equivalent set
// of mobilizer forces, ignoring constraints.
void SimbodyMatterSubsystem::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobForces) const
{
    getRep().calcTreeEquivalentMobilityForces(s,bodyForces,mobForces);
}

Real SimbodyMatterSubsystem::calcKineticEnergy(const State& s) const {
    return getRep().calcKineticEnergy(s);
}


void SimbodyMatterSubsystem::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    getRep().calcQDot(s, u, qdot);
}

void SimbodyMatterSubsystem::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    getRep().calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodyMatterSubsystem::getNBodies()        const {return getRep().getNBodies();}
int SimbodyMatterSubsystem::getNMobilities()    const {return getRep().getNMobilities();}
int SimbodyMatterSubsystem::getNConstraints()   const {return getRep().getNConstraints();}
int SimbodyMatterSubsystem::getNParticles()     const {return getRep().getNParticles();}

int SimbodyMatterSubsystem::getTotalQAlloc()    const {return getRep().getTotalQAlloc();}

// Modeling info.
void SimbodyMatterSubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { getRep().setUseEulerAngles(s,useAngles); }
void SimbodyMatterSubsystem::setMobilizerIsPrescribed(State& s, MobilizedBodyId body, bool prescribed) const
  { getRep().setMobilizerIsPrescribed(s,body,prescribed); }
void SimbodyMatterSubsystem::setConstraintIsEnabled(State& s, ConstraintId constraint, bool enabled) const
  { getRep().setConstraintIsEnabled(s,constraint,enabled); }
bool SimbodyMatterSubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodyMatterSubsystem::isMobilizerPrescribed(const State& s, MobilizedBodyId body) const
  { return getRep().isMobilizerPrescribed(s,body); }
bool SimbodyMatterSubsystem::isConstraintEnabled(const State& s, ConstraintId constraint) const
  { return getRep().isConstraintEnabled(s,constraint); }

int SimbodyMatterSubsystem::getNQuaternionsInUse(const State& s) const {
    return getRep().getNQuaternionsInUse(s);
}
bool SimbodyMatterSubsystem::isUsingQuaternion(const State& s, MobilizedBodyId body) const {
    return getRep().isUsingQuaternion(s, body);
}
int SimbodyMatterSubsystem::getQuaternionIndex(const State& s, MobilizedBodyId body) const {
    return getRep().getQuaternionIndex(s, body);
}

void SimbodyMatterSubsystem::enforcePositionConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
  { getRep().enforcePositionConstraints(s, requiredTol, desiredTol); }
void SimbodyMatterSubsystem::enforceVelocityConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
  { getRep().enforceVelocityConstraints(s, requiredTol, desiredTol); }

const SpatialVec&
SimbodyMatterSubsystem::getCoriolisAcceleration(const State& s, MobilizedBodyId body) const {
    return getRep().getCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCoriolisAcceleration(const State& s, MobilizedBodyId body) const {
    return getRep().getTotalCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getGyroscopicForce(const State& s, MobilizedBodyId body) const {
    return getRep().getGyroscopicForce(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getCentrifugalForces(const State& s, MobilizedBodyId body) const {
    return getRep().getCentrifugalForces(s,body);
}

const SpatialMat& 
SimbodyMatterSubsystem::getArticulatedBodyInertia(const State& s, MobilizedBodyId body) const {
    return getRep().getArticulatedBodyInertia(s,body);
}

const Vector& 
SimbodyMatterSubsystem::getAllParticleMasses(const State& s) const { 
    return getRep().getAllParticleMasses(s); 
}
Vector& SimbodyMatterSubsystem::updAllParticleMasses(State& s) const {
    return getRep().updAllParticleMasses(s); 
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleLocations(const State& s) const { 
    return getRep().getAllParticleLocations(s); 
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleVelocities(const State& s) const {
    return getRep().getAllParticleVelocities(s);
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleAccelerations(const State& s) const {
    return getRep().getAllParticleAccelerations(s);
}

void SimbodyMatterSubsystem::addInStationForce(const State& s, MobilizedBodyId body, const Vec3& stationInB, 
                                        const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNBodies());
    const Rotation& R_GB = getRep().getBodyTransform(s,body).R();
    bodyForces[body] += SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}
void SimbodyMatterSubsystem::addInBodyTorque(const State& s, MobilizedBodyId body, const Vec3& torqueInG,
                                      Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNBodies());
    bodyForces[body][0] += torqueInG; // no force
}
void SimbodyMatterSubsystem::addInMobilityForce(const State& s, MobilizedBodyId body, int index, Real d,
                                         Vector& mobilityForces) const 
{ 
    assert(mobilityForces.size() == getRep().getNMobilities());
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= index && index < nu);
    mobilityForces[uStart+index] += d;
}

Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleLocations(State& s) const {
    return getRep().updAllParticleLocations(s);
}
Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleVelocities(State& s) const {
    return getRep().updAllParticleVelocities(s);
}

bool SimbodyMatterSubsystem::projectQConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return getRep().projectQConstraints(s,y_err,tol,targetTol); 
}
bool SimbodyMatterSubsystem::projectUConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return getRep().projectUConstraints(s,y_err,tol,targetTol); 
}

/// Calculate the total system mass.
/// TODO: this should be precalculated.
///
/// @par Required stage
///   \c Stage::Instance
Real SimbodyMatterSubsystem::calcSystemMass(const State& s) const {
    Real mass = 0;
    for (MobilizedBodyId b(1); b < getNBodies(); ++b)
        mass += getMobilizedBody(b).getBodyMassProperties(s).getMass();
    return mass;
}


/// Return the location r_OG_C of the system mass center C, measured from the ground
/// origin OG, and expressed in Ground. 
///
/// @par Required stage
///   \c Stage::Position
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterLocationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);

    for (MobilizedBodyId b(1); b < getNBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Transform&      X_GB    = getMobilizedBody(b).getBodyTransform(s);
        const Real            mb      = MB_OB_B.getMass();
        const Vec3            r_OG_CB = X_GB * MB_OB_B.getMassCenter();
        mass += mb;
        com  += mb * r_OG_CB; // weighted by mass
    }

    if (mass != 0) 
        com /= mass;

    return com;
}


/// Return total system mass, mass center location measured from the Ground origin,
/// and system inertia taken about the Ground origin, expressed in Ground.
///
/// @par Required stage
///   \c Stage::Position
MassProperties SimbodyMatterSubsystem::calcSystemMassPropertiesInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);
    Inertia I    = Inertia(0);

    for (MobilizedBodyId b(1); b < getNBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Transform&      X_GB    = getMobilizedBody(b).getBodyTransform(s);
        const MassProperties  MB_OG_G = MB_OB_B.calcTransformedMassProps(X_GB);
        const Real            mb      = MB_OG_G.getMass();
        mass += mb;
        com  += mb * MB_OG_G.getMassCenter();
        I    += MB_OG_G.getInertia();   // already has mass built in
    }

    if (mass != 0) {
        com /= mass;
        I   /= mass;
    }

    return MassProperties(mass, com, I);
}

/// Return the system inertia matrix taken about the system center of mass,
/// expressed in Ground.
///
/// @par Required stage
///   \c Stage::Position
Inertia SimbodyMatterSubsystem::calcSystemCentralInertiaInGround(const State& s) const {
    const MassProperties M_OG_G = calcSystemMassPropertiesInGround(s);
    return M_OG_G.calcCentralInertia();
}


/// Return the velocity V_G_C = d/dt r_OG_C of the system mass center C in the Ground frame G,
/// expressed in G.
///
/// @par Required stage
///   \c Stage::Velocity
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterVelocityInGround(const State& s) const {
    Real    mass = 0;
    Vec3    comv = Vec3(0);

    for (MobilizedBodyId b(1); b < getNBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 v_G_CB = getMobilizedBody(b).calcBodyFixedPointVelocityInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        comv += mb * v_G_CB; // weighted by mass
    }

    if (mass != 0) 
        comv /= mass;

    return comv;
}

/// Return the acceleration A_G_C = d^2/dt^2 r_OG_C of the system mass center C in
/// the Ground frame G, expressed in G.
///
/// @par Required stage
///   \c Stage::Acceleration
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterAccelerationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    coma = Vec3(0);

    for (MobilizedBodyId b(1); b < getNBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 a_G_CB = getMobilizedBody(b).calcBodyFixedPointAccelerationInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        coma += mb * a_G_CB; // weighted by mass
    }

    if (mass != 0) 
        coma /= mass;

    return coma;
}

/// Return the momentum of the system as a whole (angular, linear) measured
/// in the ground frame, taken about the ground origin and expressed in ground.
/// (The linear component is independent of the "about" point.)
///
/// @par Required stage
///   \c Stage::Velocity
SpatialVec SimbodyMatterSubsystem::calcSystemMomentumAboutGroundOrigin(const State& s) const {
    SpatialVec mom(Vec3(0), Vec3(0));
    for (MobilizedBodyId b(1); b < getNBodies(); ++b) {
        const SpatialVec mom_CB_G = getMobilizedBody(b).calcBodyMomentumAboutBodyMassCenterInGround(s);
        const Vec3&      Iw = mom_CB_G[0];
        const Vec3&      mv = mom_CB_G[1];
        const Vec3       r = getMobilizedBody(b).locateBodyMassCenterOnGround(s);
        mom[0] += (Iw + r % mv); // add central angular momentum plus contribution from mass center location
        mom[1] += mv;            // just add up central linear momenta
    }
    return mom;
}


} // namespace SimTK

