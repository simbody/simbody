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
 * Implementation of MatterSubsystem, a still-abstract Subsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MatterSubsystem.h"

#include "MatterSubsystemRep.h"

namespace SimTK {

    /////////////////////////
    // MatterSubsystem //
    /////////////////////////

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
int MatterSubsystem::getNConstraints() const {
    return MatterSubsystemRep::downcast(*rep).getNConstraints();
}
int MatterSubsystem::getParent(int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getParent(bodyNum); 
}
Array<int> 
MatterSubsystem::getChildren(int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getChildren(bodyNum); 
}
const Transform&  
MatterSubsystem::getJointFrame(const State& s, int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getJointFrame(s, bodyNum); 
}
const Transform& 
MatterSubsystem::getJointFrameOnParent(const State& s, int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getJointFrameOnParent(s, bodyNum); 
}
const Vec3&  
MatterSubsystem::getBodyCenterOfMass(const State& s, int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getBodyCenterOfMass(s,bodyNum); 
}
const Transform& 
MatterSubsystem::getBodyConfiguration(const State& s, int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getBodyConfiguration(s,bodyNum); 
}
const SpatialVec& 
MatterSubsystem::getBodyVelocity(const State& s, int bodyNum) const { 
    return MatterSubsystemRep::downcast(*rep).getBodyVelocity(s,bodyNum); 
}

const Real&
MatterSubsystem::getJointQ(const State& s, int body, int axis) const { 
    return MatterSubsystemRep::downcast(*rep).getJointQ(s,body,axis); 
}
const Real&
MatterSubsystem::getJointU(const State& s, int body, int axis) const { 
    return MatterSubsystemRep::downcast(*rep).getJointU(s,body,axis); 
}

void MatterSubsystem::setJointQ(State& s, int body, int axis, const Real& q) const { 
    MatterSubsystemRep::downcast(*rep).setJointQ(s,body,axis,q); 
}
void MatterSubsystem::setJointU(State& s, int body, int axis, const Real& u) const { 
    MatterSubsystemRep::downcast(*rep).setJointU(s,body,axis,u); 
}


const Vector& MatterSubsystem::getQConstraintErrors(const State& s) const { 
    return MatterSubsystemRep::downcast(*rep).getQConstraintErrors(s); 
}
const Real& MatterSubsystem::getQConstraintNorm(const State& s) const { 
    return MatterSubsystemRep::downcast(*rep).getQConstraintNorm(s); 
}
const Vector& MatterSubsystem::getUConstraintErrors(const State& s) const { 
    return MatterSubsystemRep::downcast(*rep).getUConstraintErrors(s); 
}
const Real& MatterSubsystem::getUConstraintNorm(const State& s) const { 
    return MatterSubsystemRep::downcast(*rep).getUConstraintNorm(s); 
}
bool MatterSubsystem::projectQConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return MatterSubsystemRep::downcast(*rep).projectQConstraints(
        s,y_err,tol,targetTol); 
}
bool MatterSubsystem::projectUConstraints(State& s, Vector& y_err, Real tol, Real targetTol) const { 
    return MatterSubsystemRep::downcast(*rep).projectUConstraints(
        s,y_err,tol,targetTol); 
}

} // namespace SimTK

