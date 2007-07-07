/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of MatterSubsystem, a still-abstract Subsystem handle,
 * and MatterSubsystemRep the actual C++ abstract class representing
 * the implementation of a MatterSubsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MatterSubsystem.h"

#include "MatterSubsystemRep.h"
#include "MobilizedBodyRep.h"

namespace SimTK {

    /////////////////////
    // MATTER SUBSYTEM //
    /////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

/*static*/ bool 
MatterSubsystem::isInstanceOf(const Subsystem& s) {
    return MatterSubsystemRep::isA(s.getRep());
}
/*static*/ const MatterSubsystem&
MatterSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MatterSubsystem&>(s);
}
/*static*/ MatterSubsystem&
MatterSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MatterSubsystem&>(s);
}

const MatterSubsystemRep& 
MatterSubsystem::getRep() const {
    return dynamic_cast<const MatterSubsystemRep&>(*rep);
}
MatterSubsystemRep&       
MatterSubsystem::updRep() {
    return dynamic_cast<MatterSubsystemRep&>(*rep);
}


int MatterSubsystem::getNBodies() const {
    return getRep().getNBodies();
}
int MatterSubsystem::getNParticles() const {
    return getRep().getNParticles();
}
int MatterSubsystem::getNMobilities() const {
    return getRep().getNMobilities();
}
int MatterSubsystem::getNConstraints() const {
    return getRep().getNConstraints();
}
MobilizedBodyId MatterSubsystem::getParent(MobilizedBodyId bodyNum) const { 
    return getRep().getParent(bodyNum); 
}
Array<MobilizedBodyId> 
MatterSubsystem::getChildren(MobilizedBodyId bodyNum) const { 
    return getRep().getChildren(bodyNum); 
}

const MassProperties& MatterSubsystem::getDefaultBodyMassProperties(MobilizedBodyId body) const {
    return getRep().getDefaultBodyMassProperties(body);
}
const Transform&  
MatterSubsystem::getDefaultMobilizerFrame(MobilizedBodyId bodyNum) const { 
    return getRep().getDefaultMobilizerFrame(bodyNum); 
}
const Transform& 
MatterSubsystem::getDefaultMobilizerFrameOnParent(MobilizedBodyId bodyNum) const { 
    return getRep().getDefaultMobilizerFrameOnParent(bodyNum); 
}

const MassProperties& MatterSubsystem::getBodyMassProperties(const State& s, MobilizedBodyId body) const {
    return getRep().getBodyMassProperties(s,body);
}

const Vector& 
MatterSubsystem::getAllParticleMasses(const State& s) const { 
    return getRep().getAllParticleMasses(s); 
}

const Transform&  
MatterSubsystem::getMobilizerFrame(const State& s, MobilizedBodyId bodyNum) const { 
    return getRep().getMobilizerFrame(s, bodyNum); 
}
const Transform& 
MatterSubsystem::getMobilizerFrameOnParent(const State& s, MobilizedBodyId bodyNum) const { 
    return getRep().getMobilizerFrameOnParent(s, bodyNum); 
}

MassProperties& MatterSubsystem::updBodyMassProperties(State& s, MobilizedBodyId b) const {
    return getRep().updBodyMassProperties(s,b); 
}
Transform& MatterSubsystem::updMobilizerFrame(State& s, MobilizedBodyId b) const {
    return getRep().updMobilizerFrame(s,b); 
}
Transform& MatterSubsystem::updMobilizerFrameOnParent(State& s, MobilizedBodyId b) const {
    return getRep().updMobilizerFrameOnParent(s,b); 
}
Vector& MatterSubsystem::updAllParticleMasses(State& s) const {
    return getRep().updAllParticleMasses(s); 
}

const Transform& 
MatterSubsystem::getBodyTransform(const State& s, MobilizedBodyId bodyNum) const { 
    return getRep().getBodyTransform(s,bodyNum); 
}
const Vector_<Vec3>& 
MatterSubsystem::getAllParticleLocations(const State& s) const { 
    return getRep().getAllParticleLocations(s); 
}

const SpatialVec& 
MatterSubsystem::getBodyVelocity(const State& s, MobilizedBodyId bodyNum) const { 
    return getRep().getBodyVelocity(s,bodyNum); 
}
const Vector_<Vec3>& 
MatterSubsystem::getAllParticleVelocities(const State& s) const {
    return getRep().getAllParticleVelocities(s);
}

const SpatialVec& 
MatterSubsystem::getBodyAcceleration(const State& s, MobilizedBodyId bodyNum) const { 
    return getRep().getBodyAcceleration(s,bodyNum); 
}
const Vector_<Vec3>& 
MatterSubsystem::getAllParticleAccelerations(const State& s) const {
    return getRep().getAllParticleAccelerations(s);
}

void MatterSubsystem::addInStationForce(const State& s, MobilizedBodyId body, const Vec3& stationInB, 
                                        const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNBodies());
    const Rotation& R_GB = getRep().getBodyTransform(s,body).R();
    bodyForces[body] += SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}
void MatterSubsystem::addInBodyTorque(const State& s, MobilizedBodyId body, const Vec3& torqueInG,
                                      Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNBodies());
    bodyForces[body][0] += torqueInG; // no force
}
void MatterSubsystem::addInMobilityForce(const State& s, MobilizedBodyId body, int index, Real d,
                                         Vector& mobilityForces) const 
{ 
    assert(mobilityForces.size() == getRep().getNMobilities());
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= index && index < nu);
    mobilityForces[uStart+index] += d;
}

int MatterSubsystem::getNMobilizerCoords(const State& s, MobilizedBodyId body) const {
    int qStart, nq;
    getRep().findMobilizerQs(s,body,qStart,nq);
    return nq;
}
int MatterSubsystem::getNMobilizerSpeeds(const State& s, MobilizedBodyId body) const {
    int uStart, nu;
    getRep().findMobilizerUs(s,body,uStart,nu);
    return nu;
}

Real MatterSubsystem::getMobilizerQ(const State& s, MobilizedBodyId body, int index) const {
    int qStart, nq; getRep().findMobilizerQs(s,body,qStart,nq);
    assert(0 <= index && index < nq);
    return getRep().getQ(s)[qStart+index];
}
Real MatterSubsystem::getMobilizerU(const State& s, MobilizedBodyId body, int index) const { 
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= index && index < nu);
    return getRep().getU(s)[uStart+index];
}

void MatterSubsystem::setMobilizerQ(State& s, MobilizedBodyId body, int index, Real q) const { 
    int qStart, nq; getRep().findMobilizerQs(s,body,qStart,nq);
    assert(0 <= index && index < nq);
    getRep().updQ(s)[qStart+index] = q;
}
void MatterSubsystem::setMobilizerU(State& s, MobilizedBodyId body, int index, Real u) const { 
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= index && index < nu);
    getRep().updU(s)[uStart+index] = u;
}

const Vector& MatterSubsystem::getAllMobilizerCoords(const State& s) const {
    return getRep().getAllMobilizerCoords(s);
}

const Vector& MatterSubsystem::getAllMobilizerSpeeds(const State& s) const {
    return getRep().getAllMobilizerSpeeds(s);
}

const Vector& MatterSubsystem::getAllMobilizerAppliedForces(const State& s) const {
    return getRep().getAllMobilizerAppliedForces(s);
}

const Vector_<SpatialVec>& MatterSubsystem::getAllBodyAppliedForces(const State& s) const {
    return getRep().getAllBodyAppliedForces(s);
}


const Vector_<Vec3>& MatterSubsystem::getAllParticleAppliedForces(const State& s) const {
    return getRep().getAllParticleAppliedForces(s);
}

Real MatterSubsystem::getMobilizerCoord(const State& s, MobilizedBodyId body) const {
    int qStart, nq; getRep().findMobilizerQs(s,body,qStart,nq);
    assert(nq == 1);
    return getRep().getQ(s)[qStart];
}
Vector MatterSubsystem::getMobilizerCoords(const State& s, MobilizedBodyId body) const {
    int qStart, nq; getRep().findMobilizerQs(s,body,qStart,nq);
    return getRep().getQ(s)(qStart,nq);
}
Real MatterSubsystem::getMobilizerSpeed(const State& s, MobilizedBodyId body) const {
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(nu == 1);
    return getRep().getU(s)[uStart];
}
Vector MatterSubsystem::getMobilizerSpeeds(const State& s, MobilizedBodyId body) const {
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    return getRep().getU(s)(uStart,nu);
}
Real MatterSubsystem::getMobilizerAppliedForce(const State& s, MobilizedBodyId body) const {
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(nu == 1);
    return getRep().getAllMobilizerAppliedForces(s)[uStart];
}
Vector MatterSubsystem::getMobilizerAppliedForces(const State& s, MobilizedBodyId body) const {
    int uStart, nu; getRep().findMobilizerUs(s,body,uStart,nu);
    return getRep().getAllMobilizerAppliedForces(s)(uStart,nu);
}

const Vec<2>& MatterSubsystem::getMobilizerCoordsAsVec2(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 2);
    return Vec<2>::getAs(&r.getQ(s)[qStart]);
}
const Vec<3>& MatterSubsystem::getMobilizerCoordsAsVec3(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 3);
    return Vec<3>::getAs(&r.getQ(s)[qStart]);
}
const Vec<4>& MatterSubsystem::getMobilizerCoordsAsVec4(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 4);
    return Vec<4>::getAs(&r.getQ(s)[qStart]);
}
const Vec<5>& MatterSubsystem::getMobilizerCoordsAsVec5(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 5);
    return Vec<5>::getAs(&r.getQ(s)[qStart]);
}
const Vec<6>& MatterSubsystem::getMobilizerCoordsAsVec6(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 6);
    return Vec<6>::getAs(&r.getQ(s)[qStart]);
}
const Vec<7>& MatterSubsystem::getMobilizerCoordsAsVec7(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 7);
    return Vec<7>::getAs(&r.getQ(s)[qStart]);
}


const Vec<2>& MatterSubsystem::getMobilizerSpeedsAsVec2(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 2);
    return Vec<2>::getAs(&r.getU(s)[uStart]);
}
const Vec<3>& MatterSubsystem::getMobilizerSpeedsAsVec3(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 3);
    return Vec<3>::getAs(&r.getU(s)[uStart]);
}
const Vec<4>& MatterSubsystem::getMobilizerSpeedsAsVec4(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 4);
    return Vec<4>::getAs(&r.getU(s)[uStart]);
}
const Vec<5>& MatterSubsystem::getMobilizerSpeedsAsVec5(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 5);
    return Vec<5>::getAs(&r.getU(s)[uStart]);
}
const Vec<6>& MatterSubsystem::getMobilizerSpeedsAsVec6(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 6);
    return Vec<6>::getAs(&r.getU(s)[uStart]);
}



const Vec<2>& MatterSubsystem::getMobilizerAppliedForceAsVec2(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 2);
    return Vec<2>::getAs(&r.getAllMobilizerAppliedForces(s)[uStart]);
}
const Vec<3>& MatterSubsystem::getMobilizerAppliedForceAsVec3(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 3);
    return Vec<3>::getAs(&r.getAllMobilizerAppliedForces(s)[uStart]);
}
const Vec<4>& MatterSubsystem::getMobilizerAppliedForceAsVec4(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 4);
    return Vec<4>::getAs(&r.getAllMobilizerAppliedForces(s)[uStart]);
}
const Vec<5>& MatterSubsystem::getMobilizerAppliedForceAsVec5(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 5);
    return Vec<5>::getAs(&r.getAllMobilizerAppliedForces(s)[uStart]);
}
const Vec<6>& MatterSubsystem::getMobilizerAppliedForceAsVec6(const State& s, MobilizedBodyId body) const {
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 6);
    return Vec<6>::getAs(&r.getAllMobilizerAppliedForces(s)[uStart]);
}

Vector_<SpatialVec>& MatterSubsystem::updAllBodyAppliedForces(State& s) const {
    return getRep().updAllBodyAppliedForces(s);
}
Vector_<Vec3>& MatterSubsystem::updAllParticleLocations(State& s) const {
    return getRep().updAllParticleLocations(s);
}
Vector_<Vec3>& MatterSubsystem::updAllParticleVelocities(State& s) const {
    return getRep().updAllParticleVelocities(s);
}
Vector_<Vec3>& MatterSubsystem::updAllParticleAppliedForces(State& s) const {
    return getRep().updAllParticleAppliedForces(s);
}

void MatterSubsystem::setMobilizerCoord(State& s, MobilizedBodyId body, Real q) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 1);
    r.updQ(s)[qStart] = q;
}
void MatterSubsystem::setMobilizerCoordsAsVec2(State& s, MobilizedBodyId body, const Vec<2>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 2);
    Vec<2>::updAs(&r.updQ(s)[qStart]) = qs;
}
void MatterSubsystem::setMobilizerCoordsAsVec3(State& s, MobilizedBodyId body, const Vec<3>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 3);
    Vec<3>::updAs(&r.updQ(s)[qStart]) = qs;
}
void MatterSubsystem::setMobilizerCoordsAsVec4(State& s, MobilizedBodyId body, const Vec<4>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 4);
    Vec<4>::updAs(&r.updQ(s)[qStart]) = qs;
}
void MatterSubsystem::setMobilizerCoordsAsVec5(State& s, MobilizedBodyId body, const Vec<5>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 5);
    Vec<5>::updAs(&r.updQ(s)[qStart]) = qs;
}
void MatterSubsystem::setMobilizerCoordsAsVec6(State& s, MobilizedBodyId body, const Vec<6>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 6);
    Vec<6>::updAs(&r.updQ(s)[qStart]) = qs;
}
void MatterSubsystem::setMobilizerCoordsAsVec7(State& s, MobilizedBodyId body, const Vec<7>& qs) const {
    const MatterSubsystemRep& r = getRep();
    int qStart, nq; r.findMobilizerQs(s,body,qStart,nq);
    assert(nq == 7);
    Vec<7>::updAs(&r.updQ(s)[qStart]) = qs;
}

void MatterSubsystem::setMobilizerSpeed(State& s, MobilizedBodyId body, Real u) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 1);
    r.updU(s)[uStart] = u;
}

void MatterSubsystem::setMobilizerSpeedsAsVec2(State& s, MobilizedBodyId body, const Vec<2>& us) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 2);
    Vec<2>::updAs(&getRep().updU(s)[uStart]) = us;
}
void MatterSubsystem::setMobilizerSpeedsAsVec3(State& s, MobilizedBodyId body, const Vec<3>& us) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 3);
    Vec<3>::updAs(&getRep().updU(s)[uStart]) = us;
}
void MatterSubsystem::setMobilizerSpeedsAsVec4(State& s, MobilizedBodyId body, const Vec<4>& us) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 4);
    Vec<4>::updAs(&getRep().updU(s)[uStart]) = us;
}
void MatterSubsystem::setMobilizerSpeedsAsVec5(State& s, MobilizedBodyId body, const Vec<5>& us) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 5);
    Vec<5>::updAs(&getRep().updU(s)[uStart]) = us;
}
void MatterSubsystem::setMobilizerSpeedsAsVec6(State& s, MobilizedBodyId body, const Vec<6>& us) const{
    const MatterSubsystemRep& r = getRep();
    int uStart, nu; r.findMobilizerUs(s,body,uStart,nu);
    assert(nu == 6);
    Vec<6>::updAs(&getRep().updU(s)[uStart]) = us;
}

Real MatterSubsystem::calcQConstraintNorm(const State& s) const { 
    return getRep().calcQConstraintNorm(s); 
}
Real MatterSubsystem::calcUConstraintNorm(const State& s) const { 
    return getRep().calcUConstraintNorm(s); 
}
Real MatterSubsystem::calcUDotConstraintNorm(const State& s) const { 
    return getRep().calcUDotConstraintNorm(s); 
}

bool MatterSubsystem::projectQConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return getRep().projectQConstraints(s,y_err,tol,targetTol); 
}
bool MatterSubsystem::projectUConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return getRep().projectUConstraints(s,y_err,tol,targetTol); 
}



} // namespace SimTK

