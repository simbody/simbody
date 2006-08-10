/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of MatterSubsystem, a still-abstract Subsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MatterSubsystem.h"

#include "MatterSubsystemRep.h"

namespace SimTK {

    /////////////////////
    // MatterSubsystem //
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
int MatterSubsystem::getParent(int bodyNum) const { 
    return getRep().getParent(bodyNum); 
}
Array<int> 
MatterSubsystem::getChildren(int bodyNum) const { 
    return getRep().getChildren(bodyNum); 
}
const Transform&  
MatterSubsystem::getMobilizerFrame(const State& s, int bodyNum) const { 
    return getRep().getMobilizerFrame(s, bodyNum); 
}
const Transform& 
MatterSubsystem::getMobilizerFrameOnParent(const State& s, int bodyNum) const { 
    return getRep().getMobilizerFrameOnParent(s, bodyNum); 
}
const Real&  
MatterSubsystem::getBodyMass(const State& s, int bodyNum) const { 
    return getRep().getBodyMass(s,bodyNum); 
}
const Vec3&  
MatterSubsystem::getBodyCenterOfMassStation(const State& s, int bodyNum) const { 
    return getRep().getBodyCenterOfMassStation(s,bodyNum); 
}

const Vector& 
MatterSubsystem::getParticleMasses(const State& s) const { 
    return getRep().getParticleMasses(s); 
}

const Vector_<Vec3>& 
MatterSubsystem::getParticleLocations(const State& s) const { 
    return getRep().getParticleLocations(s); 
}

const Transform& 
MatterSubsystem::getBodyConfiguration(const State& s, int bodyNum) const { 
    return getRep().getBodyConfiguration(s,bodyNum); 
}
const SpatialVec& 
MatterSubsystem::getBodyVelocity(const State& s, int bodyNum) const { 
    return getRep().getBodyVelocity(s,bodyNum); 
}

void MatterSubsystem::addInStationForce(const State& s, int body, const Vec3& stationInB, 
                                        const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const {
    getRep().addInStationForce(s,body,stationInB,forceInG,bodyForces); 
}
void MatterSubsystem::addInBodyTorque(const State& s, int body, const Vec3& torqueInG,
                                      Vector_<SpatialVec>& bodyForces) const {
    getRep().addInBodyTorque(s,body,torqueInG,bodyForces); 
}
void MatterSubsystem::addInMobilityForce(const State& s, int body, int axis, const Real& d,
                                         Vector& mobilityForces) const { 
    getRep().addInMobilityForce(s,body,axis,d,mobilityForces); 
}

const Real&
MatterSubsystem::getMobilizerQ(const State& s, int body, int axis) const { 
    return getRep().getMobilizerQ(s,body,axis); 
}
const Real&
MatterSubsystem::getMobilizerU(const State& s, int body, int axis) const { 
    return getRep().getMobilizerU(s,body,axis); 
}

void MatterSubsystem::setMobilizerQ(State& s, int body, int axis, const Real& q) const { 
    getRep().setMobilizerQ(s,body,axis,q); 
}
void MatterSubsystem::setMobilizerU(State& s, int body, int axis, const Real& u) const { 
    getRep().setMobilizerU(s,body,axis,u); 
}


const Transform& MatterSubsystem::getMobilizerConfiguration(const State& s, int body) const { 
    return getRep().getMobilizerConfiguration(s,body); 
}
const SpatialVec& MatterSubsystem::getMobilizerVelocity(const State& s, int body) const { 
    return getRep().getMobilizerVelocity(s,body); 
}
void MatterSubsystem::setMobilizerConfiguration(State& s, int body, const Transform& X_JbJ) const { 
    getRep().setMobilizerConfiguration(s,body,X_JbJ); 
}
void MatterSubsystem::setMobilizerVelocity(State& s, int body, const SpatialVec& V_JbJ) const { 
    getRep().setMobilizerVelocity(s,body,V_JbJ); 
}


const Vector& MatterSubsystem::getQConstraintErrors(const State& s) const { 
    return getRep().getQConstraintErrors(s); 
}
Real MatterSubsystem::calcQConstraintNorm(const State& s) const { 
    return getRep().calcQConstraintNorm(s); 
}
const Vector& MatterSubsystem::getUConstraintErrors(const State& s) const { 
    return getRep().getUConstraintErrors(s); 
}
Real MatterSubsystem::calcUConstraintNorm(const State& s) const { 
    return getRep().calcUConstraintNorm(s); 
}
const Vector& MatterSubsystem::getUDotConstraintErrors(const State& s) const { 
    return getRep().getUDotConstraintErrors(s); 
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

