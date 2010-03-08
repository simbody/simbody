/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-9 Stanford University and the Authors.         *
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
 * Implementation of SimbodyMatterSubsystem, a concrete Subsystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/MobilizedBody.h"

#include "MobilizedBodyImpl.h"
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

// Create Subsystem but don't associate it with any System. This isn't much 
// use except for making std::vector's, which require a default constructor 
// to be available.
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

MobilizedBodyIndex SimbodyMatterSubsystem::adoptMobilizedBody(MobilizedBodyIndex parent, MobilizedBody& child) {
    return updRep().adoptMobilizedBody(parent,child);
}
const MobilizedBody& SimbodyMatterSubsystem::getMobilizedBody(MobilizedBodyIndex id) const {
    return getRep().getMobilizedBody(id);
}
MobilizedBody& SimbodyMatterSubsystem::updMobilizedBody(MobilizedBodyIndex id) {
    return updRep().updMobilizedBody(id);
}
const MobilizedBody::Ground& SimbodyMatterSubsystem::getGround() const {
    return getRep().getGround();
}
MobilizedBody::Ground& SimbodyMatterSubsystem::updGround() {
    return updRep().updGround();
}

bool SimbodyMatterSubsystem::getShowDefaultGeometry() const {
    return getRep().getShowDefaultGeometry();
}

void SimbodyMatterSubsystem::setShowDefaultGeometry(bool show) {
    updRep().setShowDefaultGeometry(show);
}


ConstraintIndex SimbodyMatterSubsystem::adoptConstraint(Constraint& child) {
    return updRep().adoptConstraint(child);
}
const Constraint& SimbodyMatterSubsystem::getConstraint(ConstraintIndex id) const {
    return getRep().getConstraint(id);
}
Constraint& SimbodyMatterSubsystem::updConstraint(ConstraintIndex id) {
    return updRep().updConstraint(id);
}

//TODO: should allow zero-length force arrays to stand for zeroes.
void SimbodyMatterSubsystem::calcAcceleration
   (const State&                state,
    const Vector&               appliedMobilityForces,
    const Vector_<SpatialVec>&  appliedBodyForces,
    Vector&                     udot,
    Vector_<SpatialVec>&        A_GB) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcAcceleration",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcAcceleration",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());

    Vector_<Vec3> appliedParticleForces; // TODO

    // Create a dummy acceleration cache to hold the result.
    const SBModelCache&    mc = getRep().getModelCache(state);
    const SBInstanceCache& ic = getRep().getInstanceCache(state);
    SBTreeAccelerationCache tac;
    tac.allocate(getRep().topologyCache, mc, ic);

    Vector udotErr(getNUDotErr(state)); // unwanted return value
    Vector multipliers(getNMultipliers(state)); // unwanted return value

    getRep().calcLoopForwardDynamicsOperator(state, 
        appliedMobilityForces, appliedParticleForces, appliedBodyForces,
        tac, udot, multipliers, udotErr);

    A_GB = tac.bodyAccelerationInGround;
}

//TODO: should allow zero-length force arrays to stand for zeroes.
void SimbodyMatterSubsystem::calcAccelerationIgnoringConstraints
   (const State&                state,
    const Vector&               appliedMobilityForces,
    const Vector_<SpatialVec>&  appliedBodyForces,
    Vector&                     udot, // output only; returns pres. accels
    Vector_<SpatialVec>&        A_GB) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcAccelerationIgnoringConstraints",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcAccelerationIgnoringConstraints",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());

    Vector netHingeForces(getNumMobilities()); // unwanted side effects
    Vector tau;

    getRep().calcTreeAccelerations(state,
        appliedMobilityForces, appliedBodyForces,
        netHingeForces, A_GB, udot, tau);
}


void SimbodyMatterSubsystem::calcMInverseV(const State& s,
    const Vector&        v,
    Vector&              MinvV) const
{
	Vector_<SpatialVec> A_GB;
    getRep().calcMInverseF(s,v, A_GB, MinvV);
}

void SimbodyMatterSubsystem::calcResidualForceIgnoringConstraints
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    Vector&                    residualMobilityForces) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==0 || appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==0 || appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());
    SimTK_APIARGCHECK2_ALWAYS(
        knownUdot.size()==0 || knownUdot.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d knownUdots but there are %d mobilities.",
        knownUdot.size(), getNumMobilities());

    residualMobilityForces.resize(getNumMobilities());

	Vector_<SpatialVec> A_GB(getNumBodies());
    getRep().calcTreeResidualForces(state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        A_GB, residualMobilityForces);
}

void SimbodyMatterSubsystem::calcMV(const State& s, 
	const Vector& v, 
	Vector& MV) const
{
	Vector_<SpatialVec> A_GB;
    getRep().calcMV(s,v, A_GB, MV);
}

void SimbodyMatterSubsystem::
calcPNInv(const State& s, Matrix& PNInv) const {
    return getRep().calcHolonomicConstraintMatrixPNInv(s,PNInv);
}

void SimbodyMatterSubsystem::
calcP(const State& s, Matrix& P) const {
    return getRep().calcHolonomicVelocityConstraintMatrixP(s,P);
}

void SimbodyMatterSubsystem::
calcPt(const State& s, Matrix& Pt) const {
    return getRep().calcHolonomicVelocityConstraintMatrixPt(s,Pt);
}

void SimbodyMatterSubsystem::calcCompositeBodyInertias(const State& s,
    Vector_<SpatialMat>& R) const
{
    getRep().calcCompositeBodyInertias(s,R);
}

// Convert internal kinematics to spatial equivalent, ignoring velocity and constraints.
void SimbodyMatterSubsystem::calcSpatialKinematicsFromInternal(const State& s,
    const Vector&        v,
    Vector_<SpatialVec>& Jv) const
{
    getRep().calcSpatialKinematicsFromInternal(s,v,Jv);
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

void SimbodyMatterSubsystem::calcMobilizerReactionForces(const State& s, Vector_<SpatialVec>& forces) const {
    getRep().calcMobilizerReactionForces(s, forces);
}

void SimbodyMatterSubsystem::calcConstraintForcesFromMultipliers
   (const State& s, const Vector& lambda,
    Vector_<SpatialVec>& bodyForcesInG,
    Vector&              mobilityForces) const
{
    getRep().calcConstraintForcesFromMultipliers(s,lambda,bodyForcesInG,mobilityForces);
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

void SimbodyMatterSubsystem::multiplyByN
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   getRep().multiplyByN(s,matrixOnRight,in,out); }
void SimbodyMatterSubsystem::multiplyByNInv
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   getRep().multiplyByNInv(s,matrixOnRight,in,out); }
void SimbodyMatterSubsystem::multiplyByNDot
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   getRep().multiplyByNDot(s,matrixOnRight,in,out); }


// Topological info. Note the lack of a State argument.
int SimbodyMatterSubsystem::getNumBodies()        const {return getRep().getNumBodies();}
int SimbodyMatterSubsystem::getNumMobilities()    const {return getRep().getNumMobilities();}
int SimbodyMatterSubsystem::getNumConstraints()   const {return getRep().getNumConstraints();}
int SimbodyMatterSubsystem::getNumParticles()     const {return getRep().getNumParticles();}

int SimbodyMatterSubsystem::getTotalQAlloc()    const {return getRep().getTotalQAlloc();}

// Modeling info.
void SimbodyMatterSubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { getRep().setUseEulerAngles(s,useAngles); }
void SimbodyMatterSubsystem::setMobilizerIsPrescribed(State& s, MobilizedBodyIndex body, bool prescribed) const
  { getRep().setMobilizerIsPrescribed(s,body,prescribed); }
void SimbodyMatterSubsystem::setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disabled) const
  { getRep().setConstraintIsDisabled(s,constraint,disabled); }
bool SimbodyMatterSubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodyMatterSubsystem::isMobilizerPrescribed(const State& s, MobilizedBodyIndex body) const
  { return getRep().isMobilizerPrescribed(s,body); }
bool SimbodyMatterSubsystem::isConstraintDisabled(const State& s, ConstraintIndex constraint) const
  { return getRep().isConstraintDisabled(s,constraint); }
void SimbodyMatterSubsystem::convertToEulerAngles(const State& inputState, State& outputState) const
  { return getRep().convertToEulerAngles(inputState, outputState); }
void SimbodyMatterSubsystem::convertToQuaternions(const State& inputState, State& outputState) const
  { return getRep().convertToQuaternions(inputState, outputState); }

int SimbodyMatterSubsystem::getNumQuaternionsInUse(const State& s) const {
    return getRep().getNumQuaternionsInUse(s);
}
bool SimbodyMatterSubsystem::isUsingQuaternion(const State& s, MobilizedBodyIndex body) const {
    return getRep().isUsingQuaternion(s, body);
}
QuaternionPoolIndex SimbodyMatterSubsystem::getQuaternionPoolIndex(const State& s, MobilizedBodyIndex body) const {
    return getRep().getQuaternionPoolIndex(s, body);
}
AnglePoolIndex SimbodyMatterSubsystem::getAnglePoolIndex(const State& s, MobilizedBodyIndex body) const {
    return getRep().getAnglePoolIndex(s, body);
}
const SpatialVec&
SimbodyMatterSubsystem::getCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
    return getRep().getCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
    return getRep().getTotalCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getGyroscopicForce(const State& s, MobilizedBodyIndex body) const {
    return getRep().getGyroscopicForce(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getCentrifugalForces(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getTotalCentrifugalForces(s,body);
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

void SimbodyMatterSubsystem::addInStationForce(const State& s, MobilizedBodyIndex body, const Vec3& stationInB, 
                                               const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNumBodies());
    const Rotation& R_GB = getRep().getBodyTransform(s,body).R();
    bodyForces[body] += SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void SimbodyMatterSubsystem::realizeCompositeBodyInertias(const State& s) const {
    getRep().realizeCompositeBodyInertias(s);
}

void SimbodyMatterSubsystem::realizeArticulatedBodyInertias(const State& s) const {
    getRep().realizeArticulatedBodyInertias(s);
}

const SpatialMat& 
SimbodyMatterSubsystem::getCompositeBodyInertia(const State& s, MobilizedBodyIndex mbx) const {
    return getRep().getCompositeBodyInertias(s)[mbx]; // will lazy-evaluate if necessary
}

const SpatialMat& 
SimbodyMatterSubsystem::getArticulatedBodyInertia(const State& s, MobilizedBodyIndex mbx) const {
    return getRep().getArticulatedBodyInertias(s)[mbx]; // will lazy-evaluate if necessary
}

void SimbodyMatterSubsystem::addInBodyTorque(const State& s, MobilizedBodyIndex body, const Vec3& torqueInG,
                                             Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNumBodies());
    bodyForces[body][0] += torqueInG; // no force
}
void SimbodyMatterSubsystem::addInMobilityForce(const State& s, MobilizedBodyIndex body, MobilizerUIndex which, Real d,
                                                Vector& mobilityForces) const 
{ 
    assert(mobilityForces.size() == getRep().getNumMobilities());
    UIndex uStart; int nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= which && which < nu);
    mobilityForces[uStart+which] += d;
}

Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleLocations(State& s) const {
    return getRep().updAllParticleLocations(s);
}
Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleVelocities(State& s) const {
    return getRep().updAllParticleVelocities(s);
}

bool SimbodyMatterSubsystem::prescribe(State& s, Stage g) const {
    return getRep().prescribe(s,g);
}

bool SimbodyMatterSubsystem::projectQConstraints(State& s, Real consAccuracy, const Vector& yWeights,
                                                 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{ 
    return getRep().projectQConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts); 
}
bool SimbodyMatterSubsystem::projectUConstraints(State& s, Real consAccuracy, const Vector& yWeights,
												 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{ 
    return getRep().projectUConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts); 
}

/// Calculate the total system mass.
/// TODO: this should be precalculated.
Real SimbodyMatterSubsystem::calcSystemMass(const State& s) const {
    Real mass = 0;
    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b)
        mass += getMobilizedBody(b).getBodyMassProperties(s).getMass();
    return mass;
}


// Return the location r_OG_C of the system mass center C, measured from the ground
// origin OG, and expressed in Ground. 
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterLocationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
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


// Return total system mass, mass center location measured from the Ground origin,
// and system inertia taken about the Ground origin, expressed in Ground.
MassProperties SimbodyMatterSubsystem::calcSystemMassPropertiesInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);
    Inertia I    = Inertia(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Transform&      X_GB    = getMobilizedBody(b).getBodyTransform(s);
        const MassProperties  MB_OG_G = MB_OB_B.calcTransformedMassProps(~X_GB);
        const Real            mb      = MB_OG_G.getMass();
        mass += mb;
        com  += mb * MB_OG_G.getMassCenter();
        I    += MB_OG_G.getInertia();   // already has mass built in
    }

    if (mass != 0)
        com /= mass;

    return MassProperties(mass, com, I);
}

// Return the system inertia matrix taken about the system center of mass,
// expressed in Ground.
Inertia SimbodyMatterSubsystem::calcSystemCentralInertiaInGround(const State& s) const {
    const MassProperties M_OG_G = calcSystemMassPropertiesInGround(s);
    return M_OG_G.calcCentralInertia();
}


// Return the velocity V_G_C = d/dt r_OG_C of the system mass center C in the Ground frame G,
// expressed in G.
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterVelocityInGround(const State& s) const {
    Real    mass = 0;
    Vec3    comv = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 v_G_CB = getMobilizedBody(b).findStationVelocityInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        comv += mb * v_G_CB; // weighted by mass
    }

    if (mass != 0) 
        comv /= mass;

    return comv;
}

// Return the acceleration A_G_C = d^2/dt^2 r_OG_C of the system mass center C in
// the Ground frame G, expressed in G.
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterAccelerationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    coma = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 a_G_CB = getMobilizedBody(b).findStationAccelerationInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        coma += mb * a_G_CB; // weighted by mass
    }

    if (mass != 0) 
        coma /= mass;

    return coma;
}

// Return the momentum of the system as a whole (angular, linear) measured
// in the ground frame, taken about the ground origin and expressed in ground.
// (The linear component is independent of the "about" point.)
SpatialVec SimbodyMatterSubsystem::calcSystemMomentumAboutGroundOrigin(const State& s) const {
    SpatialVec mom(Vec3(0), Vec3(0));
    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const SpatialVec mom_CB_G = getMobilizedBody(b).calcBodyMomentumAboutBodyMassCenterInGround(s);
        const Vec3&      Iw = mom_CB_G[0];
        const Vec3&      mv = mom_CB_G[1];
        const Vec3       r = getMobilizedBody(b).findMassCenterLocationInGround(s);
        mom[0] += (Iw + r % mv); // add central angular momentum plus contribution from mass center location
        mom[1] += mv;            // just add up central linear momenta
    }
    return mom;
}

// Return the momentum of the system as a whole (angular, linear) measured
// in the ground frame, taken about the current system center of mass
// location and expressed in ground.
// (The linear component is independent of the "about" point.)
SpatialVec SimbodyMatterSubsystem::calcSystemCentralMomentum(const State& s) const {
    SpatialVec mom(Vec3(0), Vec3(0));
    Real mtot = 0;  // system mass
    Vec3 com(0);    // system mass center

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MobilizedBody& mobod = getMobilizedBody(b);
        const Real m    = mobod.getBodyMass(s);
        const Vec3 CB_G = mobod.findMassCenterLocationInGround(s);
        mtot += m;
        com  += m * CB_G;
        const SpatialVec mom_CB_G = mobod.calcBodyMomentumAboutBodyMassCenterInGround(s);
        const Vec3&      Iw = mom_CB_G[0];
        const Vec3&      mv = mom_CB_G[1];
        mom[0] += (Iw + CB_G % mv); // add central angular momentum plus contribution from mass center location
        mom[1] += mv;               // just add up central linear momenta
    }
    if (mtot != 0)
        com /= mtot;

    // Shift momentum from ground origin to system COM (only angular affected).
    mom[0] -= com % mom[1];
    return mom;
}

} // namespace SimTK

