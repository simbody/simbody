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

      // TODO: THESE NEED WORK TO GET RID OF NODE
    int addRigidBodyNode(IVMRigidBodyNode&  parent,
                         const IVMFrame&    referenceConfig,    // body frame in parent
                         IVMRigidBodyNode*& nodep)
    {
        return ivmTree.addRigidBodyNode(parent,referenceConfig,nodep);
    }
    int addGroundNode(IVMRigidBodyNode*& gnodep) {
        return ivmTree.addGroundNode(gnodep);
    }
    const IVMRigidBodyNode& getRigidBodyNode(int nodeNum) const {
        return ivmTree.getRigidBodyNode(nodeNum);
    }
    IVMRigidBodyNode& updRigidBodyNode(int nodeNum) {
        return ivmTree.updRigidBodyNode(nodeNum);
    }

    int addDistanceConstraint(const IVMStation& s1, const IVMStation& s2,
                              const double& d) {
        return ivmTree.addDistanceConstraint(s1,s2,d);
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


