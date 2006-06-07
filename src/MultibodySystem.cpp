/* Copyright (c) 2006 Stanford University and Michael Sherman.
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
 * Implementation of MultibodySystem, a concrete System.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

#include "MultibodySystemRep.h"

namespace SimTK {

    /////////////////////////
    // MechanicalSubsystem //
    /////////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

/*static*/ bool 
MechanicalSubsystem::isInstanceOf(const Subsystem& s) {
    return MechanicalSubsystemRep::isA(s.getRep());
}
/*static*/ const MechanicalSubsystem&
MechanicalSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MechanicalSubsystem&>(s);
}
/*static*/ MechanicalSubsystem&
MechanicalSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MechanicalSubsystem&>(s);
}

const MechanicalSubsystemRep& 
MechanicalSubsystem::getRep() const {
    return dynamic_cast<const MechanicalSubsystemRep&>(*rep);
}
MechanicalSubsystemRep&       
MechanicalSubsystem::updRep() {
    return dynamic_cast<MechanicalSubsystemRep&>(*rep);
}

int MechanicalSubsystem::getNBodies() const {
    return getRep().getNBodies();
}
int MechanicalSubsystem::getNConstraints() const {
    return MechanicalSubsystemRep::downcast(*rep).getNConstraints();
}
int MechanicalSubsystem::getParent(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getParent(bodyNum); 
}
Array<int> 
MechanicalSubsystem::getChildren(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getChildren(bodyNum); 
}
const Transform&  
MechanicalSubsystem::getJointFrame(const State& s, int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointFrame(s, bodyNum); 
}
const Transform& 
MechanicalSubsystem::getJointFrameOnParent(const State& s, int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointFrameOnParent(s, bodyNum); 
}
const Vec3&  
MechanicalSubsystem::getBodyCenterOfMass(const State& s, int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getBodyCenterOfMass(s,bodyNum); 
}
const Transform& 
MechanicalSubsystem::getBodyConfiguration(const State& s, int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getBodyConfiguration(s,bodyNum); 
}
const SpatialVec& 
MechanicalSubsystem::getBodyVelocity(const State& s, int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getBodyVelocity(s,bodyNum); 
}

const Real&
MechanicalSubsystem::getJointQ(const State& s, int body, int axis) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointQ(s,body,axis); 
}
const Real&
MechanicalSubsystem::getJointU(const State& s, int body, int axis) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointU(s,body,axis); 
}

void MechanicalSubsystem::setJointQ(State& s, int body, int axis, const Real& q) const { 
    MechanicalSubsystemRep::downcast(*rep).setJointQ(s,body,axis,q); 
}
void MechanicalSubsystem::setJointU(State& s, int body, int axis, const Real& u) const { 
    MechanicalSubsystemRep::downcast(*rep).setJointU(s,body,axis,u); 
}

const Vector& MechanicalSubsystem::getQConstraintErrors(const State& s) const { 
    return MechanicalSubsystemRep::downcast(*rep).getQConstraintErrors(s); 
}
const Real& MechanicalSubsystem::getQConstraintNorm(const State& s) const { 
    return MechanicalSubsystemRep::downcast(*rep).getQConstraintNorm(s); 
}
const Vector& MechanicalSubsystem::getUConstraintErrors(const State& s) const { 
    return MechanicalSubsystemRep::downcast(*rep).getUConstraintErrors(s); 
}
const Real& MechanicalSubsystem::getUConstraintNorm(const State& s) const { 
    return MechanicalSubsystemRep::downcast(*rep).getUConstraintNorm(s); 
}
bool MechanicalSubsystem::projectQConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return MechanicalSubsystemRep::downcast(*rep).projectQConstraints(
        s,y_err,tol,targetTol); 
}
bool MechanicalSubsystem::projectUConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return MechanicalSubsystemRep::downcast(*rep).projectUConstraints(
        s,y_err,tol,targetTol); 
}


    ///////////////////////////////
    // MechanicalForcesSubsystem //
    ///////////////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.


/*static*/ bool 
MechanicalForcesSubsystem::isInstanceOf(const Subsystem& s) {
    return MechanicalForcesSubsystemRep::isA(s.getRep());
}
/*static*/ const MechanicalForcesSubsystem&
MechanicalForcesSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MechanicalForcesSubsystem&>(s);
}
/*static*/ MechanicalForcesSubsystem&
MechanicalForcesSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MechanicalForcesSubsystem&>(s);
}

    /////////////////////
    // MultibodySystem //
    /////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

MultibodySystem::MultibodySystem() {
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);
}

MultibodySystem::MultibodySystem(MechanicalSubsystem& m, 
                                 MechanicalForcesSubsystem& f)
{
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);

    setMechanicalSubsystem(m);
    setMechanicalForcesSubsystem(f);
}

bool MultibodySystem::project(State& s, Vector& y_err, 
             const Real& tol,
             const Real& dontProjectFac,
             const Real& targetTol
             ) const
{
    return MultibodySystemRep::downcast(*rep).project(
                s,y_err,tol,dontProjectFac,targetTol);
}

MechanicalSubsystem&       
MultibodySystem::setMechanicalSubsystem(MechanicalSubsystem& m) {
    return MultibodySystemRep::downcast(*rep).setMechanicalSubsystem(m);
}
MechanicalForcesSubsystem& 
MultibodySystem::setMechanicalForcesSubsystem(MechanicalForcesSubsystem& f) {
    return MultibodySystemRep::downcast(*rep).setMechanicalForcesSubsystem(f);
}

const MechanicalSubsystem&       
MultibodySystem::getMechanicalSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMechanicalSubsystem();
}

const MechanicalForcesSubsystem& 
MultibodySystem::getMechanicalForcesSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMechanicalForcesSubsystem();
}

MechanicalSubsystem&       
MultibodySystem::updMechanicalSubsystem() {
    return MultibodySystemRep::downcast(*rep).updMechanicalSubsystem();
}

MechanicalForcesSubsystem& 
MultibodySystem::updMechanicalForcesSubsystem() {
    return MultibodySystemRep::downcast(*rep).updMechanicalForcesSubsystem();
}

// TODO: camera facing, screen fixed, calculated geometry (e.g. line between stations
// on two different bodies, marker at system COM)
void MultibodySystem::addAnalyticGeometry(int body, const Transform& X_BG, const AnalyticGeometry& g) {
    MultibodySystemRep::downcast(*rep).addAnalyticGeometry(body,X_BG,g);
}

void MultibodySystem::addDecorativeGeometry(int body, const Transform& X_BG, const DecorativeGeometry& g) {
    MultibodySystemRep::downcast(*rep).addDecorativeGeometry(body,X_BG,g);
}

const Array<AnalyticGeometry>&   
MultibodySystem::getBodyAnalyticGeometry(int body) {
    return MultibodySystemRep::downcast(*rep).getBodyAnalyticGeometry(body);
}

const Array<DecorativeGeometry>& 
MultibodySystem::getBodyDecorativeGeometry(int body) {
    return MultibodySystemRep::downcast(*rep).getBodyDecorativeGeometry(body);
}

} // namespace SimTK

