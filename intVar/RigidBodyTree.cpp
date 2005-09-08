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
typedef CDSList<Vec6>         VecVec6;

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
void RigidBodyTree::destructNode(RigidBodyNode* n) {
    for (int i=0; i < n->getNChildren(); i++)
        destructNode( n->updChild(i) );
    rbNodeLevels[n->getLevel()].remove( rbNodeLevels[n->getLevel()].getIndex(n) );
    delete n;
}

void RigidBodyTree::setPos(const RVec& pos)  {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setPos(pos); 
}

// setPos() & calcConfigurationKinematics() must have been called
void RigidBodyTree::setVel(const RVec& vel)  {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setVel(vel); 
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcP() {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcP();
}


// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcZ(const VecVec6& spatialForces) {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(spatialForces[node.getNodeNum()]);
        }
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcPandZ(const VecVec6& spatialForces) {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcPandZ(spatialForces[node.getNodeNum()]);
        }
}

//
// Y is used for length constraints.
// Recursion is base to tip.
//
void RigidBodyTree::calcY() {
    for (int i=0; i < rbNodeLevels.size(); i++)
        for (int j=0; j < rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcY();
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcGetAccel(RVec& acc) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel();
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getAccel(acc);
}

// Calc P quantities: sweep from tip to base.
void RigidBodyTree::getAccel(const VecVec6& spatialForces, RVec& acc) {
    calcPandZ(spatialForces);
    calcGetAccel(acc);

    if ( ivm->lConstraints->fixAccel() ) 
        for (int i=0 ; i<rbNodeLevels.size() ; i++)
            for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
                rbNodeLevels[i][j]->getAccel(acc);
}

void
RigidBodyTree::updateAccel(const VecVec6& spatialForces) {
    calcZ(spatialForces);

    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel();
}

// Calc internal force: sweep from tip to base.
void RigidBodyTree::getInternalForce(const VecVec6& spatialForces, RVec& T) {
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalForce(spatialForces[node.getNodeNum()]);
        }

    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getInternalForce(T);

    ivm->lConstraints->fixGradient(T);
}

void 
RigidBodyTree::enforceConstraints(RVec& pos, RVec& vel) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->enforceConstraints(pos,vel);

    ivm->lConstraints->enforce(pos,vel); //FIX: previous constraints still obeyed?
}

void RigidBodyTree::getPos(RVec& pos) const {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getPos(pos); 
}

void RigidBodyTree::getVel(RVec& vel) const {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getVel(vel); 
}
