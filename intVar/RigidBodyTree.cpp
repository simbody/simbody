/**@file
 *
 * Implementation of RigidBodyTree.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include "LengthConstraints.h"

#include "vec3.h"

#include "sthead.h"
#include "cdsMath.h"
#include "cdsString.h"
#include "cdsSStream.h"
#include "cdsVector.h"
#include "fixedVector.h"
#include "subVector.h"
#include "fixedMatrix.h"
#include "matrixTools.h"
#include "cdsIomanip.h"
#include "cdsFstream.h"

using namespace InternalDynamics;
using MatrixTools::inverse;

typedef FixedVector<double,6> Vec6;
typedef FixedMatrix<double,6> Mat66;
typedef CDSList<Vec6>         VecVec6;

RigidBodyTree::~RigidBodyTree() {
    delete lConstraints;

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
        destructNode( n->getChild(i) );
    rbNodeLevels[n->getLevel()].remove( rbNodeLevels[n->getLevel()].getIndex(n) );
    delete n;
}

// Add a new node, taking over the heap space.
int RigidBodyTree::addRigidBodyNode(RigidBodyNode&  parent,
                                    const Frame&    referenceConfig,    // body frame in parent
                                    RigidBodyNode*& nodep)
{
    RigidBodyNode* n = nodep; nodep=0;  // take ownership
    const int level = parent.getLevel() + 1;
    n->setLevel(level);

    // Put node in tree at the right level
    if (rbNodeLevels.size()<=level) rbNodeLevels.resize(level+1);
    const int nxt = rbNodeLevels[level].size();
    rbNodeLevels[level].append(n);

    // Assign a unique reference integer to this node, for use by caller
    const int nodeNum = nodeNum2NodeMap.size();
    nodeNum2NodeMap.append(RigidBodyNodeIndex(level,nxt));
    n->setNodeNum(nodeNum);

    // Link in to the tree topology (bidirectional).
    parent.addChild(n);

    return nodeNum;
}

// Add a distance constraint and allocate slots to hold the runtime information for
// its stations. Return the assigned distance constraint index for caller's use.
int RigidBodyTree::addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d)
{
    const int indx = distanceConstraints.size();
    distanceConstraints.append(RBDistanceConstraint(s1,s2,d));
    RBDistanceConstraint& dc = distanceConstraints[indx];
    int rtnxt = stationRuntimeInfo.size();
    distanceConstraints[indx].grabStationRuntimes(rtnxt);
    stationRuntimeInfo.resize(rtnxt);
    return indx;
}

void RigidBodyTree::finishConstruction(IVM* ivm) {
    lConstraints = new LengthConstraints(ivm);
    lConstraints->construct(distanceConstraints);
}

// Set generalized coordinates: sweep from base to tips.
void RigidBodyTree::setPos(const RVec& pos)  {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setPos(pos); 
}

// Set generalized speeds: sweep from base to tip.
// setPos() must have been called already.
void RigidBodyTree::setVel(const RVec& vel)  {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->setVel(vel); 
}

// Enforce coordinate constraints -- order doesn't matter.
void RigidBodyTree::enforceTreeConstraints(RVec& pos, RVec& vel) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->enforceConstraints(pos,vel);
}

// Enforce loop constraints.
void RigidBodyTree::enforceConstraints(RVec& pos, RVec& vel) {
    lConstraints->enforce(pos,vel); //FIX: previous constraints still obeyed? (CDS)
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

// Y is used for length constraints: sweep from base to tip.
void RigidBodyTree::calcY() {
    for (int i=0; i < rbNodeLevels.size(); i++)
        for (int j=0; j < rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcY();
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcTreeAccel() {
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel();
}

// Acceleration fixup to account for constraints.
void RigidBodyTree::fixAccelForConstraints() {
    lConstraints->fixAccel();
}

void RigidBodyTree::fixVel0(RVec& vel) {
    lConstraints->fixVel0(vel);
}

// Calc unconstrained internal forces from spatial forces: sweep from tip to base.
void RigidBodyTree::calcTreeInternalForces(const VecVec6& spatialForces) {
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalForce(spatialForces[node.getNodeNum()]);
        }
}

// Retrieve already-computed internal forces (order doesn't matter).
void RigidBodyTree::getInternalForces(RVec& T) {
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getInternalForce(T);
}

void RigidBodyTree::getConstraintCorrectedInternalForces(RVec& T) {
    lConstraints->fixGradient(T);
}

// Get current generalized coordinates (order doesn't matter).
void RigidBodyTree::getPos(RVec& pos) const {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getPos(pos); 
}

// Get current generalized speeds (order doesn't matter).
void RigidBodyTree::getVel(RVec& vel) const {
    for (int i=0 ; i<rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getVel(vel); 
}

// Retrieve already-calculated accelerations (order doesn't matter)
void RigidBodyTree::getAcc(RVec& acc) const {
    for (int i=0 ; i<rbNodeLevels.size() ; i++)
        for (int j=0 ; j<rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getAccel(acc);
}
