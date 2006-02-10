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
        rep->nextUSlot, rep->nextQSlot);
    cout << "CREATED: states " << save << "-" << rep->nextUSlot-1 << endl;

    RigidBodyNode& pn = rep->updRigidBodyNode(parent);
    const int rbIndex = rep->addRigidBodyNode(pn,TransformMat()/*XXX*/,rb);
    return rbIndex;
}
