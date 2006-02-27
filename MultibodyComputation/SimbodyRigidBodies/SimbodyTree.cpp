/**@file
 *
 * Implementation of SimbodyTree.
 */

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;


SimbodyTree::SimbodyTree() {
    rep = new RigidBodyTree();
}

int SimbodyTree::addRigidBody(
    int                       parent,
    const TransformMat&       parentJointFrameInP,  // X_PJb
    const JointSpecification& joint,
    const TransformMat&       bodyJointFrameInB,    // X_BJ
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

// Note the lack of a State argument when completing construction.
void SimbodyTree::realizeConstruction() {
    rep->realizeConstruction();
}
void SimbodyTree::realizeModeling     (const SBState& s) const {rep->realizeModeling(s.getRep());}
void SimbodyTree::realizeParameters   (const SBState& s) const {rep->realizeParameters(s.getRep());}
void SimbodyTree::realizeTime         (const SBState& s) const {rep->realizeTime(s.getRep());}
void SimbodyTree::realizeConfiguration(const SBState& s) const {rep->realizeConfiguration(s.getRep());}
void SimbodyTree::realizeMotion       (const SBState& s) const {rep->realizeMotion(s.getRep());}
void SimbodyTree::realizeReaction     (const SBState& s) const {rep->realizeReaction(s.getRep());}
void SimbodyTree::realize(const SBState& s, SBStage g) const {
    rep->realize(s.getRep(), g);
}

void SimbodyTree::calcInternalGradientFromSpatial(const SBState& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    rep->calcInternalGradientFromSpatial(s.getRep(),dEdR,dEdQ);
}

Real SimbodyTree::calcKineticEnergy(const SBState& s) const {
    return rep->calcKineticEnergy(s.getRep());
}

void SimbodyTree::calcTreeUDot(const SBState& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot) const
{
    rep->calcTreeUDot(s.getRep(),jointForces,bodyForces,udot);
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

const SBState& SimbodyTree::getInitialState() const {return rep->getInitialState();}

// Modeling info.
void SimbodyTree::setUseEulerAngles(SBState& s, bool useAngles) const
  { rep->setUseEulerAngles(s.updRep(),useAngles); }
void SimbodyTree::setJointIsPrescribed(SBState& s, int joint, bool prescribed) const
  { rep->setJointIsPrescribed(s.updRep(),joint,prescribed); }
void SimbodyTree::setConstraintIsEnabled(SBState& s, int constraint, bool enabled) const
  { rep->setConstraintIsEnabled(s.updRep(),constraint,enabled); }
bool SimbodyTree::getUseEulerAngles(const SBState& s) const
  { return rep->getUseEulerAngles(s.getRep()); }
bool SimbodyTree::isJointPrescribed(const SBState& s, int joint) const
  { return rep->isJointPrescribed(s.getRep(),joint); }
bool SimbodyTree::isConstraintEnabled(const SBState& s, int constraint) const
  { return rep->isConstraintEnabled(s.getRep(),constraint); }


const Vector& SimbodyTree::getQ(const SBState& s) const {return rep->getQ(s.getRep());}
const Vector& SimbodyTree::getU(const SBState& s) const {return rep->getU(s.getRep());}

const Vector&
SimbodyTree::getAppliedJointForces(const SBState& s) const {
    return rep->getAppliedJointForces(s.getRep());
}
const Vector_<SpatialVec>&
SimbodyTree::getAppliedBodyForces(const SBState& s) const {
    return rep->getAppliedBodyForces(s.getRep());
}

void SimbodyTree::setQ(SBState& s, const Vector& q) const {rep->setQ(s.updRep(),q);}
void SimbodyTree::setU(SBState& s, const Vector& u) const {rep->setU(s.updRep(),u);}
VectorView& SimbodyTree::updQ(SBState& s) const {return rep->updQ(s.updRep());}
VectorView& SimbodyTree::updU(SBState& s) const {return rep->updU(s.updRep());}

void SimbodyTree::setJointQ(SBState& s, int body, int axis, const Real& r) const
  { return rep->setJointQ(s.updRep(),body,axis,r); }
void SimbodyTree::setJointU(SBState& s, int body, int axis, const Real& r) const
  { return rep->setJointU(s.updRep(),body,axis,r); }

void SimbodyTree::setPrescribedUdot(SBState& s, int body, int axis, const Real& r) const
  { return rep->setPrescribedUdot(s.updRep(),body,axis,r); }

void SimbodyTree::clearAppliedForces(SBState& s) const {rep->clearAppliedForces(s.updRep());}
void SimbodyTree::applyGravity(SBState& s, const Vec3& g) const {rep->applyGravity(s.updRep(),g);}
void SimbodyTree::applyPointForce(SBState& s, int body, const Vec3& stationInB, 
                                                        const Vec3& forceInG) const 
  { rep->applyPointForce(s.updRep(),body,stationInB,forceInG); }

void SimbodyTree::applyBodyTorque(SBState& s, int body, const Vec3& torqueInG) const 
  { rep->applyBodyTorque(s.updRep(),body,torqueInG); }
void SimbodyTree::applyJointForce(SBState& s, int body, int axis, const Real& d) const
  { rep->applyJointForce(s.updRep(),body,axis,d); }

void SimbodyTree::enforceConfigurationConstraints(SBState& s) const
  { rep->enforceConfigurationConstraints(s.updRep()); }
void SimbodyTree::enforceMotionConstraints(SBState& s) const
  { rep->enforceMotionConstraints(s.updRep()); }

const TransformMat&
SimbodyTree::getBodyConfiguration(const SBState& s, int body) const
  { return rep->getRigidBodyNode(body).getX_GB(s.getRep()); }

const SpatialVec&
SimbodyTree::getBodyVelocity(const SBState& s, int body) const
  { return rep->getRigidBodyNode(body).getV_GB(s.getRep()); }

const SpatialVec&
SimbodyTree::getBodyAcceleration(const SBState& s, int body) const
  { return rep->getRigidBodyNode(body).getA_GB(s.getRep()); }

const Vector& SimbodyTree::getQDot   (const SBState& s) const {return rep->getQDot(s.getRep());}
const Vector& SimbodyTree::getUDot   (const SBState& s) const {return rep->getUDot(s.getRep());}
const Vector& SimbodyTree::getQDotDot(const SBState& s) const {return rep->getQDotDot(s.getRep());}
