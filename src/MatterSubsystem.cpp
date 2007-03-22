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
BodyId MatterSubsystem::getParent(BodyId bodyNum) const { 
    return getRep().getParent(bodyNum); 
}
Array<BodyId> 
MatterSubsystem::getChildren(BodyId bodyNum) const { 
    return getRep().getChildren(bodyNum); 
}

const MassProperties& MatterSubsystem::getBodyMassProperties(const State& s, BodyId body) const {
    return getRep().getBodyMassProperties(s,body);
}

const Vector& 
MatterSubsystem::getParticleMasses(const State& s) const { 
    return getRep().getParticleMasses(s); 
}

const Transform&  
MatterSubsystem::getMobilizerFrame(const State& s, BodyId bodyNum) const { 
    return getRep().getMobilizerFrame(s, bodyNum); 
}
const Transform& 
MatterSubsystem::getMobilizerFrameOnParent(const State& s, BodyId bodyNum) const { 
    return getRep().getMobilizerFrameOnParent(s, bodyNum); 
}


const Transform& 
MatterSubsystem::getBodyPosition(const State& s, BodyId bodyNum) const { 
    return getRep().getBodyPosition(s,bodyNum); 
}
const Vector_<Vec3>& 
MatterSubsystem::getParticleLocations(const State& s) const { 
    return getRep().getParticleLocations(s); 
}

const SpatialVec& 
MatterSubsystem::getBodyVelocity(const State& s, BodyId bodyNum) const { 
    return getRep().getBodyVelocity(s,bodyNum); 
}
const SpatialVec& 
MatterSubsystem::getBodyAcceleration(const State& s, BodyId bodyNum) const { 
    return getRep().getBodyAcceleration(s,bodyNum); 
}

void MatterSubsystem::addInStationForce(const State& s, BodyId body, const Vec3& stationInB, 
                                        const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const {
    getRep().addInStationForce(s,body,stationInB,forceInG,bodyForces); 
}
void MatterSubsystem::addInBodyTorque(const State& s, BodyId body, const Vec3& torqueInG,
                                      Vector_<SpatialVec>& bodyForces) const {
    getRep().addInBodyTorque(s,body,torqueInG,bodyForces); 
}
void MatterSubsystem::addInMobilityForce(const State& s, BodyId body, int axis, const Real& d,
                                         Vector& mobilityForces) const { 
    getRep().addInMobilityForce(s,body,axis,d,mobilityForces); 
}

const Real&
MatterSubsystem::getMobilizerQ(const State& s, BodyId body, int axis) const { 
    return getRep().getMobilizerQ(s,body,axis); 
}
const Real&
MatterSubsystem::getMobilizerU(const State& s, BodyId body, int axis) const { 
    return getRep().getMobilizerU(s,body,axis); 
}

void MatterSubsystem::setMobilizerQ(State& s, BodyId body, int axis, const Real& q) const { 
    getRep().setMobilizerQ(s,body,axis,q); 
}
void MatterSubsystem::setMobilizerU(State& s, BodyId body, int axis, const Real& u) const { 
    getRep().setMobilizerU(s,body,axis,u); 
}


const Transform& MatterSubsystem::getMobilizerPosition(const State& s, BodyId body) const { 
    return getRep().getMobilizerPosition(s,body); 
}
const SpatialVec& MatterSubsystem::getMobilizerVelocity(const State& s, BodyId body) const { 
    return getRep().getMobilizerVelocity(s,body); 
}
void MatterSubsystem::setMobilizerPosition(State& s, BodyId body, const Transform& X_JbJ) const { 
    getRep().setMobilizerPosition(s,body,X_JbJ); 
}
void MatterSubsystem::setMobilizerVelocity(State& s, BodyId body, const SpatialVec& V_JbJ) const { 
    getRep().setMobilizerVelocity(s,body,V_JbJ); 
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

