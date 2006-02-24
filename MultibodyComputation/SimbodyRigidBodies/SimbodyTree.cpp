/**@file
 *
 * Implementation of SimbodyTree.
 */

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
#include "RigidBodyTree.h"

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
    rep->realizeConstruction(1e-6,0);   // TODO don't *even* ask me about those arguments
}
void SimbodyTree::realizeModeling(const SBState& s) const {
    rep->realizeModeling(s);
}
void SimbodyTree::realizeParameters(const SBState& s) const {
    rep->realizeParameters(s);
}
void SimbodyTree::realizeConfiguration(const SBState& s) const {
    rep->realizeConfiguration(s);
}
void SimbodyTree::realizeMotion(const SBState& s) const {
    rep->realizeMotion(s);
}
void SimbodyTree::realizeReaction(const SBState& s) const {
    rep->realizeReaction(s);
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
void SimbodyTree::setUseEulerAngles(SBState& s, bool useAngles) const {
    rep->setUseEulerAngles(s,useAngles);
}
void SimbodyTree::setJointIsPrescribed(SBState& s, int joint, bool prescribed) const {
    rep->setJointIsPrescribed(s,joint,prescribed);
}
void SimbodyTree::setConstraintIsEnabled(SBState& s, int constraint, bool enabled) const {
    rep->setConstraintIsEnabled(s,constraint,enabled);
}

bool SimbodyTree::getUseEulerAngles(const SBState& s) const {
    return rep->getUseEulerAngles(s);
}
bool SimbodyTree::isJointPrescribed(const SBState& s, int joint) const {
    return rep->isJointPrescribed(s,joint);
}
bool SimbodyTree::isConstraintEnabled(const SBState& s, int constraint) const {
    return rep->isConstraintEnabled(s,constraint);
}
const SBState&
SimbodyTree::getDefaultState() const {
    return rep->getDefaultState();
}

const Vector&
SimbodyTree::getQ(const SBState& s) const {
    return rep->getQ(s);
}
const Vector&
SimbodyTree::getU(const SBState& s) const {
    return rep->getU(s);
}
const Vector&
SimbodyTree::getAppliedJointForces(const SBState& s) const {
    return rep->getAppliedJointForces(s);
}
const Vector_<SpatialVec>&
SimbodyTree::getAppliedBodyForces(const SBState& s) const {
    return rep->getAppliedBodyForces(s);
}

void SimbodyTree::setQ(SBState& s, const Vector& q) const {
    rep->setQ(s,q);
}

void SimbodyTree::setU(SBState& s, const Vector& u) const {
    rep->setU(s,u);
}

VectorView& SimbodyTree::updQ(SBState& s) const {
    return rep->updQ(s);
}

VectorView& SimbodyTree::updU(SBState& s) const {
    return rep->updU(s);
}

void SimbodyTree::setJointQ(SBState& s, int body, int axis, const Real& r) const {
    return rep->setJointQ(s,body,axis,r);
}
void SimbodyTree::setJointU(SBState& s, int body, int axis, const Real& r) const {
    return rep->setJointU(s,body,axis,r);
}


void SimbodyTree::setPrescribedUdot(SBState& s, int body, int axis, const Real& r) const {
    return rep->setPrescribedUdot(s,body,axis,r);
}

void SimbodyTree::clearAppliedForces(SBState& s) const {
    rep->clearAppliedForces(s);
}
void SimbodyTree::applyGravity(SBState& s, const Vec3& g) const {
    rep->applyGravity(s,g);
}
void SimbodyTree::applyPointForce(SBState& s, int body, const Vec3& stationInB, 
                                  const Vec3& forceInG) const 
{
    rep->applyPointForce(s,body,stationInB,forceInG);
}
void SimbodyTree::applyBodyTorque(SBState& s, int body, 
                                  const Vec3& torqueInG) const 
{
    rep->applyBodyTorque(s,body,torqueInG);
}
void SimbodyTree::applyJointForce(SBState& s, int body, int axis, const Real& d) const {
    rep->applyJointForce(s,body,axis,d);
}

const TransformMat&
SimbodyTree::getBodyConfiguration(const SBState& s, int body) const {
    return rep->getRigidBodyNode(body).getX_GB(s);
}

const SpatialVec&
SimbodyTree::getBodyVelocity(const SBState& s, int body) const {
    return rep->getRigidBodyNode(body).getV_GB(s);
}

const SpatialVec&
SimbodyTree::getBodyAcceleration(const SBState& s, int body) const {
    return rep->getRigidBodyNode(body).getA_GB(s);
}
