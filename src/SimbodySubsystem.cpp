/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * Implementation of SimbodySubsystem, a concrete MechanicalSubsystem.
 */

#include "Simbody.h"
#include "RigidBodyTree.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {

SimbodySubsystem::SimbodySubsystem() : MechanicalSubsystem() {
    rep = new RigidBodyTree();
    rep->setMyHandle(*this);
}


const RigidBodyTree& 
SimbodySubsystem::getRep() const {
    return dynamic_cast<const RigidBodyTree&>(*rep);
}
RigidBodyTree&       
SimbodySubsystem::updRep() {
    return dynamic_cast<RigidBodyTree&>(*rep);
}

int SimbodySubsystem::addRigidBody(
    const MassProperties&     mp,
    const Transform&          bodyJointFrameInB,    // X_BJ
    int                       parent,
    const Transform&          parentJointFrameInP,  // X_PJb
    const JointSpecification& joint)
{
    const int save = getRep().nextUSlot;

    RigidBodyNode& pn = updRep().updRigidBodyNode(parent);
    const int rbIndex = updRep().addRigidBodyNode(pn,
        mp, parentJointFrameInP, bodyJointFrameInB, joint.getJointType(), joint.isReversed(),
        updRep().nextUSlot, updRep().nextUSqSlot, updRep().nextQSlot);

    //cout << "CREATED BODY " << rbIndex << ": U states " << save << "-" << getRep().nextUSlot-1 << endl;
    return rbIndex;
}

int SimbodySubsystem::addConstantDistanceConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC,
     const Real& distance)
{
    return updRep().addConstantDistanceConstraint(
       getRep().getRigidBodyNode(parent), stationInP,
       getRep().getRigidBodyNode(child),  stationInC,
       distance);
}

int SimbodySubsystem::addCoincidentStationsConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC)
{
    return updRep().addCoincidentStationsConstraint(
       getRep().getRigidBodyNode(parent), stationInP,
       getRep().getRigidBodyNode(child),  stationInC);
}

int SimbodySubsystem::addWeldConstraint
    (int parent, const Transform& frameInP,
     int child,  const Transform& frameInC)
{
    return updRep().addWeldConstraint(
       getRep().getRigidBodyNode(parent), frameInP,
       getRep().getRigidBodyNode(child),  frameInC);
}

// Note the lack of a State argument when completing construction.
void SimbodySubsystem::endConstruction() {updRep().endConstruction();}

void SimbodySubsystem::calcInternalGradientFromSpatial(const State& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    getRep().calcInternalGradientFromSpatial(s,dEdR,dEdQ);
}

void SimbodySubsystem::calcTreeEquivalentJointForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    jointForces) const
{
    getRep().calcTreeEquivalentJointForces(s,bodyForces,jointForces);
}

Real SimbodySubsystem::calcKineticEnergy(const State& s) const {
    return getRep().calcKineticEnergy(s);
}

void SimbodySubsystem::calcTreeUDot(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot) const
{
    Vector              netHingeForces; // unwanted side effects
    Vector_<SpatialVec> A_GB;

    getRep().calcTreeAccelerations(s,jointForces,bodyForces,
        netHingeForces, A_GB, udot);
}

void SimbodySubsystem::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    getRep().calcQDot(s, u, qdot);
}

void SimbodySubsystem::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    getRep().calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodySubsystem::getNBodies()        const {return getRep().getNBodies();}
int SimbodySubsystem::getTotalDOF()       const {return getRep().getTotalDOF();}
int SimbodySubsystem::getTotalQAlloc()    const {return getRep().getTotalQAlloc();}
int SimbodySubsystem::getNConstraints()   const {return getRep().getNConstraints();}
int SimbodySubsystem::getTotalMultAlloc() const {return getRep().getTotalMultAlloc();}

int SimbodySubsystem::getQIndex(int body) const {return getRep().getQIndex(body);}
int SimbodySubsystem::getQAlloc(int body) const {return getRep().getQAlloc(body);}
int SimbodySubsystem::getUIndex(int body) const {return getRep().getUIndex(body);}
int SimbodySubsystem::getDOF   (int body) const {return getRep().getDOF(body);}

int SimbodySubsystem::getMultIndex(int constraint) const {return getRep().getMultIndex(constraint);}
int SimbodySubsystem::getMaxNMult (int constraint) const {return getRep().getMaxNMult(constraint);}

// Modeling info.
void SimbodySubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { getRep().setUseEulerAngles(s,useAngles); }
void SimbodySubsystem::setJointIsPrescribed(State& s, int joint, bool prescribed) const
  { getRep().setJointIsPrescribed(s,joint,prescribed); }
void SimbodySubsystem::setConstraintIsEnabled(State& s, int constraint, bool enabled) const
  { getRep().setConstraintIsEnabled(s,constraint,enabled); }
bool SimbodySubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodySubsystem::isJointPrescribed(const State& s, int joint) const
  { return getRep().isJointPrescribed(s,joint); }
bool SimbodySubsystem::isConstraintEnabled(const State& s, int constraint) const
  { return getRep().isConstraintEnabled(s,constraint); }


const Vector& SimbodySubsystem::getQ(const State& s) const {return getRep().getQ(s);}
const Vector& SimbodySubsystem::getU(const State& s) const {return getRep().getU(s);}

const Vector&
SimbodySubsystem::getAppliedJointForces(const State& s) const {
    return getRep().getAppliedJointForces(s);
}
const Vector_<SpatialVec>&
SimbodySubsystem::getAppliedBodyForces(const State& s) const {
    return getRep().getAppliedBodyForces(s);
}

void SimbodySubsystem::setQ(State& s, const Vector& q) const {getRep().setQ(s,q);}
void SimbodySubsystem::setU(State& s, const Vector& u) const {getRep().setU(s,u);}
Vector& SimbodySubsystem::updQ(State& s) const {return getRep().updQ(s);}
Vector& SimbodySubsystem::updU(State& s) const {return getRep().updU(s);}

void SimbodySubsystem::setJointQ(State& s, int body, int axis, const Real& r) const
  { return getRep().setJointQ(s,body,axis,r); }
void SimbodySubsystem::setJointU(State& s, int body, int axis, const Real& r) const
  { return getRep().setJointU(s,body,axis,r); }

const Real& SimbodySubsystem::getJointQ(const State& s, int body, int axis) const
  { return getRep().getJointQ(s,body,axis); }
const Real& SimbodySubsystem::getJointU(const State& s, int body, int axis) const
  { return getRep().getJointU(s,body,axis); }

void SimbodySubsystem::setPrescribedUdot(State& s, int body, int axis, const Real& r) const
  { return getRep().setPrescribedUdot(s,body,axis,r); }

void SimbodySubsystem::clearAppliedForces(State& s) const {getRep().clearAppliedForces(s);}
void SimbodySubsystem::applyGravity(State& s, const Vec3& g) const {getRep().applyGravity(s,g);}
void SimbodySubsystem::applyPointForce(State& s, int body, const Vec3& stationInB, 
                                                        const Vec3& forceInG) const 
  { getRep().applyPointForce(s,body,stationInB,forceInG); }

void SimbodySubsystem::applyBodyTorque(State& s, int body, const Vec3& torqueInG) const 
  { getRep().applyBodyTorque(s,body,torqueInG); }
void SimbodySubsystem::applyJointForce(State& s, int body, int axis, const Real& d) const
  { getRep().applyJointForce(s,body,axis,d); }

void SimbodySubsystem::enforceConfigurationConstraints(State& s) const
  { getRep().enforceConfigurationConstraints(s); }
void SimbodySubsystem::enforceMotionConstraints(State& s) const
  { getRep().enforceMotionConstraints(s); }

const Transform&
SimbodySubsystem::getBodyConfiguration(const State& s, int body) const
  { return getRep().getBodyConfiguration(s,body); }

const SpatialVec&
SimbodySubsystem::getBodyVelocity(const State& s, int body) const
  { return getRep().getBodyVelocity(s,body); }

const SpatialVec&
SimbodySubsystem::getBodyAcceleration(const State& s, int body) const
  { return getRep().getBodyAcceleration(s,body); }

const Vector& SimbodySubsystem::getQDot   (const State& s) const {return getRep().getQDot(s);}
const Vector& SimbodySubsystem::getUDot   (const State& s) const {return getRep().getUDot(s);}
const Vector& SimbodySubsystem::getQDotDot(const State& s) const {return getRep().getQDotDot(s);}

} // namespace SimTK

