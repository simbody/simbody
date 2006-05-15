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

int MechanicalSubsystem::getNBodies() const {
    return MechanicalSubsystemRep::downcast(*rep).getNBodies();
}
int MechanicalSubsystem::getNConstraints() const {
    return MechanicalSubsystemRep::downcast(*rep).getNConstraints();
}
int MechanicalSubsystem::getParent(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getParent(bodyNum); 
}
const Array<int>& 
MechanicalSubsystem::getChildren(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getChildren(bodyNum); 
}
const Transform&  
MechanicalSubsystem::getJointFrame(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointFrame(bodyNum); 
}
const Transform& 
MechanicalSubsystem::getJointFrameOnParent(int bodyNum) const { 
    return MechanicalSubsystemRep::downcast(*rep).getJointFrameOnParent(bodyNum); 
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

void MechanicalSubsystem::realizeParameters(const State& s) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeParameters(s); 
}
void MechanicalSubsystem::realizeTime(const State& s) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeParameters(s); 
}
void MechanicalSubsystem::realizeConfiguration(const State& s) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeTime(s); 
}
void MechanicalSubsystem::realizeMotion(const State& s) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeMotion(s); 
}
void MechanicalSubsystem::realizeDynamics(const State& s, const MechanicalForcesSubsystem& f) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeDynamics(s,f); 
}
void MechanicalSubsystem::realizeReaction(const State& s, const MechanicalForcesSubsystem& f) const { 
    MechanicalSubsystemRep::downcast(*rep).realizeReaction(s,f); 
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

    ///////////////////////////////
    // MechanicalForcesSubsystem //
    ///////////////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

void MechanicalForcesSubsystem::realizeParameters(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeParameters(s,m); 
}
void MechanicalForcesSubsystem::realizeTime(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeParameters(s,m); 
}
void MechanicalForcesSubsystem::realizeConfiguration(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeTime(s,m); 
}
void MechanicalForcesSubsystem::realizeMotion(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeMotion(s,m); 
}
void MechanicalForcesSubsystem::realizeDynamics(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeDynamics(s,m); 
}
void MechanicalForcesSubsystem::realizeReaction(const State& s, const MechanicalSubsystem& m) const { 
    MechanicalForcesSubsystemRep::downcast(*rep).realizeReaction(s,m); 
}

    /////////////////////
    // MultibodySystem //
    /////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

MultibodySystem::MultibodySystem(const MechanicalSubsystem& m, 
                                 const MechanicalForcesSubsystem& f)
{
    rep = new MultibodySystemRep(m,f);
    rep->setMyHandle(*this);
}

const MechanicalSubsystem&       
MultibodySystem::getMechanicalSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMechanicalSubsystem();
}

const MechanicalForcesSubsystem& 
MultibodySystem::getMechanicalForcesSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMechanicalForcesSubsystem();
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

void MultibodySystem::realizeConstruction (State& s)       const {MultibodySystemRep::downcast(*rep).realizeConstruction(s); }
void MultibodySystem::realizeModeling     (State& s)       const {MultibodySystemRep::downcast(*rep).realizeModeling(s);     }
void MultibodySystem::realizeParameters   (const State& s) const {MultibodySystemRep::downcast(*rep).realizeParameters(s);   }
void MultibodySystem::realizeTime         (const State& s) const {MultibodySystemRep::downcast(*rep).realizeTime(s);         }
void MultibodySystem::realizeConfiguration(const State& s) const {MultibodySystemRep::downcast(*rep).realizeConfiguration(s);}
void MultibodySystem::realizeMotion       (const State& s) const {MultibodySystemRep::downcast(*rep).realizeMotion(s);       }
void MultibodySystem::realizeDynamics     (const State& s) const {MultibodySystemRep::downcast(*rep).realizeDynamics(s);     }
void MultibodySystem::realizeReaction     (const State& s) const {MultibodySystemRep::downcast(*rep).realizeReaction(s);     }

} // namespace SimTK

