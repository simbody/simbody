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
    const TransformMat&       parentJointFrameInP,
    const JointSpecification& joint,
    const TransformMat&       bodyJointFrameInB,
    const MassProperties&     mp)
{
    const int save = rep->nextUSlot;
    RigidBodyNode* rb = RigidBodyNode::create(
        mp,
        bodyJointFrameInB,
        joint.getJointType(),
        joint.isReversed(),
        false,
        rep->nextUSlot, rep->nextUSqSlot, rep->nextQSlot);
    cout << "CREATED: states " << save << "-" << rep->nextUSlot-1 << endl;

    RigidBodyNode& pn = rep->updRigidBodyNode(parent);
    const int rbIndex = rep->addRigidBodyNode(pn,TransformMat()/*XXX*/,rb);
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
