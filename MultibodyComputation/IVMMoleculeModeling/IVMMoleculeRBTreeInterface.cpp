/**@file
 * Sherm 060126: This is support code for the abstract IVMMoleculeRBTreeInterface class.
 */

#include "IVMMoleculeRBTreeInterface.h"


    // INTERFACE TO ORIGINAL IVM RIGID BODY CODE

#include "IVMRigidBodyTree.h"

class IVMRigidBodyTreeMolIfc : public  IVMMoleculeRBTreeInterface {
public:
    IVMRigidBodyTreeMolIfc() { }
    ~IVMRigidBodyTreeMolIfc() { }
    IVMMoleculeRBTreeInterface* clone() const 
      { return new IVMRigidBodyTreeMolIfc(*this); }

    int addRigidBodyNode(int parentNodeNum,
                         const IVMFrame&          refConfig,    // body frame in parent
                         const IVMMassProperties& massProps,    // mass properties in body frame
                         const IVMFrame&          jointFrame,   // inboard joint frame J in body frame
                         const IVMJointType&      jtype,
                         bool                     isReversed,
                         bool                     useEuler,
                         int&                     nxtStateOffset) 
    {
        IVMRigidBodyNode* nodep = IVMRigidBodyNode::create(
                                        massProps, jointFrame, jtype, 
                                        isReversed, useEuler, nxtStateOffset); 
        return ivmTree.addRigidBodyNode(ivmTree.updRigidBodyNode(parentNodeNum),
                                        refConfig, nodep);
    }

    int addGroundNode() {
        int dummy;
        IVMRigidBodyNode* gnodep = IVMRigidBodyNode::create(
                                    IVMMassProperties(),
                                    IVMFrame(),
                                    IVMThisIsGround,
                                    false, false, dummy);
        return ivmTree.addGroundNode(gnodep);
    }

    void RBNodeSetVelFromSVel(int nodeNum, const CDSVec6& sVel) {
        IVMRigidBodyNode& rb = ivmTree.updRigidBodyNode(nodeNum);
        rb.setVelFromSVel(sVel);
    }
    void RBNodeGetConfig(int nodeNum, CDSMat33& R_GB, CDSVec3& OB_G) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        R_GB = rb.getR_GB(); OB_G = rb.getOB_G();
    }

    CDSVec6 RBNodeGetSpatialVel(int nodeNum) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        return rb.getSpatialVel();
    }

    double RBNodeGetMass(int nodeNum) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        return rb.getMass();
    }

    CDSVec3 RBNodeGetCOM_G(int nodeNum) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        return rb.getCOM_G();
    }

    double RBNodeCalcKineticEnergy(int nodeNum) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        return rb.calcKineticEnergy();
    }

    void RBNodeDump(int nodeNum, std::ostream& s) const {
        const IVMRigidBodyNode& rb = ivmTree.getRigidBodyNode(nodeNum);
        return rb.nodeDump(s);
    }

    int addDistanceConstraint(int nodeNum1, const CDSVec3& station1,
                              int nodeNum2, const CDSVec3& station2,
                              double distance)
    {
        IVMStation s1(ivmTree.updRigidBodyNode(nodeNum1), station1);
        IVMStation s2(ivmTree.updRigidBodyNode(nodeNum2), station2);
        return ivmTree.addDistanceConstraint(s1,s2,distance);
    }

    void finishConstruction(const double& ctol, int verbose) {
        ivmTree.finishConstruction(ctol,verbose);
    }

    int getNBodies() const {return ivmTree.getNBodies();}

    int getDOF() const {return ivmTree.getDOF();} 
    int getDim() const {return ivmTree.getDim();} 

    void setPos(const RVec& pos) {ivmTree.setPos(pos);}
    void setVel(const RVec& vel) {ivmTree.setVel(vel);}

    void getPos(RVec& pos) const {ivmTree.getPos(pos);}
    void getVel(RVec& vel) const {ivmTree.getVel(vel);}
    void getAcc(RVec& acc) const {ivmTree.getAcc(acc);}
    
    void velFromCartesian(const RVec& pos, RVec& vel) {
        ivmTree.velFromCartesian(pos,vel);
    }
    void enforceTreeConstraints(RVec& pos, RVec& vel) {
        ivmTree.enforceTreeConstraints(pos,vel);
    }
    void enforceConstraints(RVec& pos, RVec& vel) {
        ivmTree.enforceConstraints(pos,vel);
    }
    void prepareForDynamics() {
        ivmTree.prepareForDynamics();
    }
    void calcTreeForwardDynamics(const CDSVecVec6& spatialForces) {
        ivmTree.calcTreeForwardDynamics(spatialForces);
    }
    void calcLoopForwardDynamics(const CDSVecVec6& spatialForces) {
        ivmTree.calcLoopForwardDynamics(spatialForces);
    }

    void calcP() {ivmTree.calcP();}
    void calcZ(const CDSVecVec6& spatialForces) {
        ivmTree.calcZ(spatialForces);
    }
    void calcTreeAccel()    {ivmTree.calcTreeAccel();}
    void fixVel0(RVec& vel) {ivmTree.fixVel0(vel);}
    void calcY()            {ivmTree.calcY();}
    void calcTreeInternalForces(const CDSVecVec6& spatialForces) {
        ivmTree.calcTreeInternalForces(spatialForces);
    }
    void getInternalForces(RVec& T) {
        ivmTree.getInternalForces(T);
    }

    void getConstraintCorrectedInternalForces(RVec& T) {
        ivmTree.getConstraintCorrectedInternalForces(T);
    }
private:
    IVMRigidBodyTree ivmTree;
};

/*static*/ IVMMoleculeRBTreeInterface*
IVMMoleculeRBTreeInterface::create(bool newStyle) {
    assert(!newStyle);
    return new IVMRigidBodyTreeMolIfc();
}


