/**@file
 *
 * Implementation of RigidBodyTree.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include "dint-loop.h"

#include "vec3.h"

#include <sthead.h>
#include <cdsMath.h>
#include <cdsString.h>
#include <cdsSStream.h>
#include <cdsVector.h>
#include <fixedVector.h>
#include <subVector.h>
#include <fixedMatrix.h>
#include <matrixTools.h>
#include <cdsIomanip.h>
#include <cdsFstream.h>

using namespace InternalDynamics;
using MatrixTools::inverse;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat66;

RigidBodyTree::~RigidBodyTree() {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) {
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) 
            delete rbNodeLevels[i][j];
        rbNodeLevels[i].resize(0);
    }
    rbNodeLevels.resize(0);
}

//
// destructNode - should also work on partially constructed objects
//
void
RigidBodyTree::destructNode(RigidBodyNode* n) {
    for (int i=0; i < n->getNChildren(); i++)
        destructNode( n->updChild(i) );
    rbNodeLevels[n->getLevel()].remove( rbNodeLevels[n->getLevel()].getIndex(n) );
    delete n;
}


//
// Add node to tree. Takes over ownership of the node's heap space.
//
void
RigidBodyTree::addRigidBodyNode(RigidBodyNode* node)
{
    int level = node->level;
    if ( rbNodeLevels.size()<=level )
        rbNodeLevels.resize(level+1);
    rbNodeLevels[level].append(node);
}


// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
RigidBodyTree::calcP() {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcP();
}


// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
RigidBodyTree::calcZ() {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcZ();
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void
RigidBodyTree::calcPandZ() {
    // level 0 for atoms whose position is fixed
    // InternalDynamics::RecurseTipToBase<HingeNode::calcPandZ> d1;
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcPandZ();
}

// Calc acceleration: sweep from base to tip.
RVec
RigidBodyTree::calcGetAccel() {
    RVec acc( ivm->dim() );
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel();
    for (l_int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getAccel(acc);
    return acc;
}

// Calc P quantities: sweep from tip to base.
RVec
RigidBodyTree::getAccel() {
    calcPandZ();
    RVec acc = calcGetAccel();

    if ( ivm->lConstraints->fixAccel() ) 
        for (int i=0 ; i<rbNodeLevels.size() ; i++)
            for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
                rbNodeLevels[i][j]->getAccel(acc);

    return acc;
}

void
RigidBodyTree::updateAccel() {
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcZ();
    for (l_int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel();
}

// Calc internal force: sweep from tip to base.
RVec
RigidBodyTree::getInternalForce() {
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcInternalForce();

    RVec T( ivm->dim() );
    for (l_int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getInternalForce(T);

    ivm->lConstraints->fixGradient(T);

    return T;
}

void 
RigidBodyTree::setPosVel(const RVec& pos, const RVec& vel) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setPosVel(pos,vel); 
}

void 
RigidBodyTree::setVel(const RVec& vel)  {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setVel(vel); 
}

void 
RigidBodyTree::enforceConstraints(RVec& pos, RVec& vel) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->enforceConstraints(pos,vel);

    ivm->lConstraints->enforce(pos,vel); //FIX: previous constraints still obeyed?
}

RVec
RigidBodyTree::getPos() const {
    RVec pos(ivm->dim());
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getPos(pos); 
    return pos;
}

RVec
RigidBodyTree::getVel() const {
    RVec vel(ivm->dim());
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getVel(vel); 
    return vel;
}
