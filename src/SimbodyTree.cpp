/**@file
 *
 * Implementation of SimbodyTree.
 */

#include "Simbody.h"
#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;


SimbodyTree::SimbodyTree() {
    rep = new RigidBodyTree();
}

SimbodyTree::~SimbodyTree() {
    delete rep;
}

int SimbodyTree::addRigidBody(
    int                       parent,
    const Transform&       parentJointFrameInP,  // X_PJb
    const JointSpecification& joint,
    const Transform&       bodyJointFrameInB,    // X_BJ
    const MassProperties&     mp)
{
    const int save = rep->nextUSlot;

    RigidBodyNode& pn = rep->updRigidBodyNode(parent);
    const int rbIndex = rep->addRigidBodyNode(pn,
        mp, parentJointFrameInP, bodyJointFrameInB, joint.getJointType(), joint.isReversed(),
        rep->nextUSlot, rep->nextUSqSlot, rep->nextQSlot);

    cout << "CREATED BODY " << rbIndex << ": U states " << save << "-" << rep->nextUSlot-1 << endl;
    return rbIndex;
}

int SimbodyTree::addConstantDistanceConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC,
     const Real& distance)
{
    return rep->addConstantDistanceConstraint(
       rep->getRigidBodyNode(parent), stationInP,
       rep->getRigidBodyNode(child),  stationInC,
       distance);
}

int SimbodyTree::addCoincidentStationsConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC)
{
    return rep->addCoincidentStationsConstraint(
       rep->getRigidBodyNode(parent), stationInP,
       rep->getRigidBodyNode(child),  stationInC);
}

int SimbodyTree::addWeldConstraint
    (int parent, const Transform& frameInP,
     int child,  const Transform& frameInC)
{
    return rep->addWeldConstraint(
       rep->getRigidBodyNode(parent), frameInP,
       rep->getRigidBodyNode(child),  frameInC);
}

// Note the lack of a State argument when completing construction.
void SimbodyTree::realizeConstruction (State& s)             {rep->realizeConstruction(s);}
void SimbodyTree::realizeModeling     (State& s)       const {rep->realizeModeling(s);}

void SimbodyTree::realizeParameters   (const State& s) const {rep->realizeParameters(s);}
void SimbodyTree::realizeTime         (const State& s) const {rep->realizeTime(s);}
void SimbodyTree::realizeConfiguration(const State& s) const {rep->realizeConfiguration(s);}
void SimbodyTree::realizeMotion       (const State& s) const {rep->realizeMotion(s);}
void SimbodyTree::realizeDynamics     (const State& s) const {rep->realizeDynamics(s);}
void SimbodyTree::realizeReaction     (const State& s) const {rep->realizeReaction(s);}


void SimbodyTree::calcInternalGradientFromSpatial(const State& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    rep->calcInternalGradientFromSpatial(s,dEdR,dEdQ);
}

void SimbodyTree::calcTreeEquivalentJointForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    jointForces) const
{
    rep->calcTreeEquivalentJointForces(s,bodyForces,jointForces);
}

Real SimbodyTree::calcKineticEnergy(const State& s) const {
    return rep->calcKineticEnergy(s);
}

void SimbodyTree::calcTreeUDot(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot) const
{
    Vector              netHingeForces; // unwanted side effects
    Vector_<SpatialVec> A_GB;

    rep->calcTreeAccelerations(s,jointForces,bodyForces,
        netHingeForces, A_GB, udot);
}

void SimbodyTree::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    rep->calcQDot(s, u, qdot);
}

void SimbodyTree::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    rep->calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodyTree::getNBodies()        const {return rep->getNBodies();}
int SimbodyTree::getTotalDOF()       const {return rep->getTotalDOF();}
int SimbodyTree::getTotalQAlloc()    const {return rep->getTotalQAlloc();}
int SimbodyTree::getNConstraints()   const {return rep->getNConstraints();}
int SimbodyTree::getTotalMultAlloc() const {return rep->getTotalMultAlloc();}

int SimbodyTree::getQIndex(int body) const {return rep->getQIndex(body);}
int SimbodyTree::getQAlloc(int body) const {return rep->getQAlloc(body);}
int SimbodyTree::getUIndex(int body) const {return rep->getUIndex(body);}
int SimbodyTree::getDOF   (int body) const {return rep->getDOF(body);}

int SimbodyTree::getMultIndex(int constraint) const {return rep->getMultIndex(constraint);}
int SimbodyTree::getMaxNMult (int constraint) const {return rep->getMaxNMult(constraint);}

// Modeling info.
void SimbodyTree::setUseEulerAngles(State& s, bool useAngles) const
  { rep->setUseEulerAngles(s,useAngles); }
void SimbodyTree::setJointIsPrescribed(State& s, int joint, bool prescribed) const
  { rep->setJointIsPrescribed(s,joint,prescribed); }
void SimbodyTree::setConstraintIsEnabled(State& s, int constraint, bool enabled) const
  { rep->setConstraintIsEnabled(s,constraint,enabled); }
bool SimbodyTree::getUseEulerAngles(const State& s) const
  { return rep->getUseEulerAngles(s); }
bool SimbodyTree::isJointPrescribed(const State& s, int joint) const
  { return rep->isJointPrescribed(s,joint); }
bool SimbodyTree::isConstraintEnabled(const State& s, int constraint) const
  { return rep->isConstraintEnabled(s,constraint); }


const Vector& SimbodyTree::getQ(const State& s) const {return rep->getQ(s);}
const Vector& SimbodyTree::getU(const State& s) const {return rep->getU(s);}

const Vector&
SimbodyTree::getAppliedJointForces(const State& s) const {
    return rep->getAppliedJointForces(s);
}
const Vector_<SpatialVec>&
SimbodyTree::getAppliedBodyForces(const State& s) const {
    return rep->getAppliedBodyForces(s);
}

void SimbodyTree::setQ(State& s, const Vector& q) const {rep->setQ(s,q);}
void SimbodyTree::setU(State& s, const Vector& u) const {rep->setU(s,u);}
Vector& SimbodyTree::updQ(State& s) const {return rep->updQ(s);}
Vector& SimbodyTree::updU(State& s) const {return rep->updU(s);}

void SimbodyTree::setJointQ(State& s, int body, int axis, const Real& r) const
  { return rep->setJointQ(s,body,axis,r); }
void SimbodyTree::setJointU(State& s, int body, int axis, const Real& r) const
  { return rep->setJointU(s,body,axis,r); }

void SimbodyTree::setPrescribedUdot(State& s, int body, int axis, const Real& r) const
  { return rep->setPrescribedUdot(s,body,axis,r); }

void SimbodyTree::clearAppliedForces(State& s) const {rep->clearAppliedForces(s);}
void SimbodyTree::applyGravity(State& s, const Vec3& g) const {rep->applyGravity(s,g);}
void SimbodyTree::applyPointForce(State& s, int body, const Vec3& stationInB, 
                                                        const Vec3& forceInG) const 
  { rep->applyPointForce(s,body,stationInB,forceInG); }

void SimbodyTree::applyBodyTorque(State& s, int body, const Vec3& torqueInG) const 
  { rep->applyBodyTorque(s,body,torqueInG); }
void SimbodyTree::applyJointForce(State& s, int body, int axis, const Real& d) const
  { rep->applyJointForce(s,body,axis,d); }

void SimbodyTree::enforceConfigurationConstraints(State& s) const
  { rep->enforceConfigurationConstraints(s); }
void SimbodyTree::enforceMotionConstraints(State& s) const
  { rep->enforceMotionConstraints(s); }

const Transform&
SimbodyTree::getBodyConfiguration(const State& s, int body) const
  { return rep->getRigidBodyNode(body).getX_GB(s); }

const SpatialVec&
SimbodyTree::getBodyVelocity(const State& s, int body) const
  { return rep->getRigidBodyNode(body).getV_GB(s); }

const SpatialVec&
SimbodyTree::getBodyAcceleration(const State& s, int body) const
  { return rep->getRigidBodyNode(body).getA_GB(s); }

const Vector& SimbodyTree::getQDot   (const State& s) const {return rep->getQDot(s);}
const Vector& SimbodyTree::getUDot   (const State& s) const {return rep->getUDot(s);}
const Vector& SimbodyTree::getQDotDot(const State& s) const {return rep->getQDotDot(s);}
