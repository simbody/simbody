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

    int getInitialStateOffset() const {return 1;}

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


    // INTERFACE TO NEW SIMBODY RIGID BODY CODE

#include "RigidBodyTree.h"


static CDSVec3 toCDSVec3(const Vec3& v) {return CDSVec3(v[0],v[1],v[2]);}
static Vec3    toVec3(const CDSVec3& v) {return Vec3(v[0],v[1],v[2]);}
static CDSVec6 toCDSVec6(const SpatialVec& v) { 
    CDSVec6 cv;
    for (int i=0; i<3; ++i) cv[i] = v[0][i], cv[i+3] = v[1][i];
    return cv;
}
static SpatialVec toSpatialVec(const CDSVec6& v) 
  { return SpatialVec(Vec3(v[0],v[1],v[2]),Vec3(v[3],v[4],v[5])); }

static CDSVecVec6 toCDSVecVec6(const Vector_<SpatialVec>& a) {
    CDSVecVec6 vv(a.size());
    for (int i=0; i < (int)a.size(); ++i)
        vv(i) = toCDSVec6(a[i]);
    return vv;
}

static Vector_<SpatialVec> toSpatialVecList(const CDSVecVec6& vv) {
    Vector_<SpatialVec> svl(vv.size());
    for (int i=0; i < (int)vv.size(); ++i)
        svl[i] = toSpatialVec(vv(i));
    return svl;
}

static CDSMat33 toCDSMat33(const Mat33& m) {
    return CDSMat33(m(0,0), m(0,1), m(0,2),
                    m(1,0), m(1,1), m(1,2),
                    m(2,0), m(2,1), m(2,2));
}
static RVec toRVec(const Vector& v) {
    RVec r(v.size());   // a 1-based vector
    for (int i=0; i < v.size(); ++i)
        r[i+1] = v[i];
    return r;
}
static Vector toVector(const RVec& r) {
    return Vector(r.size(), r.pointer());
}
static Mat33 toMat33(const CDSMat33& m) {
    return Mat33(Row3(m(0,0), m(0,1), m(0,2)),
                 Row3(m(1,0), m(1,1), m(1,2)),
                 Row3(m(2,0), m(2,1), m(2,2)));
}
static IVMInertia toIVMInertia(const InertiaMat& i) {
    return IVMInertia(toCDSMat33(i.toMat33()));
}

static InertiaMat toMatInertia(const IVMInertia& i) {
    return InertiaMat(toMat33(i));
}
static RotationMat toMatRotation(const CDSMat33& m) {
    const Mat33 m33 = toMat33(m);
    return reinterpret_cast<const RotationMat&>(m33);
}

static IVMMassProperties toIVMMassProperties(const Real& m, const Vec3& c, const InertiaMat& i) {
    return IVMMassProperties(m, toCDSVec3(c), toIVMInertia(i));
}

static MassProperties toMassProperties(const IVMMassProperties& mp) {
    return MassProperties(mp.getMass(), toVec3(mp.getCOM()), 
                          toMatInertia(mp.getInertia()));
}

static IVMFrame toIVMFrame(const TransformMat& f) {
    return IVMFrame(toCDSMat33(f.R().asMat33()), toCDSVec3(f.T()));
}

static TransformMat toFrame(const IVMFrame& f) {
    return TransformMat(toMatRotation(f.getRot_RF()), toVec3(f.getLoc_RF()));
}

static Joint::JointType
toJointType(IVMJointType jt) {
    switch (jt) {
    case IVMUnknownJointType:        return Joint::UnknownJointType;
    case IVMThisIsGround:            return Joint::ThisIsGround;
    case IVMWeldJoint:               return Joint::Weld;
    case IVMTorsionJoint:            return Joint::Torsion;  // aka PinJoint
    case IVMSlidingJoint:            return Joint::Sliding;
    case IVMUJoint:                  return Joint::Universal;
    case IVMCylinderJoint:           return Joint::Cylinder;
    case IVMPlanarJoint:             return Joint::Planar;
    case IVMGimbalJoint:             return Joint::Gimbal;
    case IVMOrientationJoint:        return Joint::Orientation; // aka BallJoint
    case IVMCartesianJoint:          return Joint::Cartesian;
    case IVMFreeLineJoint:           return Joint::FreeLine;
    case IVMFreeJoint:               return Joint::Free;
    default: assert(false);
    }
    //NOTREACHED
    return Joint::UnknownJointType;
}

class SimbodyRigidBodyTreeMolIfc : public  IVMMoleculeRBTreeInterface {
public:
    SimbodyRigidBodyTreeMolIfc() { }
    ~SimbodyRigidBodyTreeMolIfc() { }

    IVMMoleculeRBTreeInterface* clone() const 
      { return new SimbodyRigidBodyTreeMolIfc(*this); }

    int getInitialStateOffset() const {return 0;}

    int addRigidBodyNode(int parentNodeNum,
                         const IVMFrame&          IVMrefConfig,    // body frame in parent
                         const IVMMassProperties& IVMmassProps,    // mass properties in body frame
                         const IVMFrame&          IVMjointFrame,   // inboard joint frame J in body frame
                         const IVMJointType&      IVMjtype,
                         bool                     isReversed,
                         bool                     useEuler,
                         int&                     nxtStateOffset) 
    {
        const TransformMat     refConfig  = toFrame(IVMrefConfig);
        const MassProperties   massProps  = toMassProperties(IVMmassProps);
        const TransformMat     jointFrame = toFrame(IVMjointFrame);
        const Joint::JointType jtype = toJointType(IVMjtype);

        RigidBodyNode* nodep = RigidBodyNode::create(
                                        massProps, jointFrame, jtype, 
                                        isReversed, useEuler, nxtStateOffset); 
        return simTree.addRigidBodyNode(simTree.updRigidBodyNode(parentNodeNum),
                                        refConfig, nodep);
    }

    int addGroundNode() {
        int dummy;
        RigidBodyNode* gnodep = RigidBodyNode::create(
                                    MassProperties(),
                                    TransformMat(),
                                    Joint::ThisIsGround,
                                    false, false, dummy);
        return simTree.addGroundNode(gnodep);
    }

    void RBNodeSetVelFromSVel(int nodeNum, const CDSVec6& sVel) {
        RigidBodyNode& rb = simTree.updRigidBodyNode(nodeNum);
        rb.setVelFromSVel(toSpatialVec(sVel));
    }
    void RBNodeGetConfig(int nodeNum, CDSMat33& R_GB, CDSVec3& OB_G) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        R_GB = toCDSMat33(rb.getR_GB().asMat33()); OB_G = toCDSVec3(rb.getOB_G());
    }

    CDSVec6 RBNodeGetSpatialVel(int nodeNum) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        return toCDSVec6(rb.getSpatialVel());
    }

    double RBNodeGetMass(int nodeNum) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        return rb.getMass();
    }

    CDSVec3 RBNodeGetCOM_G(int nodeNum) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        return toCDSVec3(rb.getCOM_G());
    }

    double RBNodeCalcKineticEnergy(int nodeNum) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        return rb.calcKineticEnergy();
    }

    void RBNodeDump(int nodeNum, std::ostream& s) const {
        const RigidBodyNode& rb = simTree.getRigidBodyNode(nodeNum);
        return rb.nodeDump(s);
    }

    int addDistanceConstraint(int nodeNum1, const CDSVec3& station1,
                              int nodeNum2, const CDSVec3& station2,
                              double distance)
    {
        RBStation s1(simTree.updRigidBodyNode(nodeNum1), toVec3(station1));
        RBStation s2(simTree.updRigidBodyNode(nodeNum2), toVec3(station2));
        return simTree.addDistanceConstraint(s1,s2,distance);
    }

    void finishConstruction(const double& ctol, int verbose) {
        simTree.realizeConstruction(ctol,verbose);
    }

    int getNBodies() const {return simTree.getNBodies();}

    int getDOF() const {return simTree.getDOF();} 
    int getDim() const {return simTree.getDim();} 

    void setPos(const RVec& pos) {simTree.setPos(toVector(pos));}
    void setVel(const RVec& vel) {simTree.setVel(toVector(vel));}

    void getPos(RVec& pos) const {
        Vector p(pos.size());
        simTree.getPos(p);
        pos = toRVec(p);
    }
    void getVel(RVec& vel) const {
        Vector v(vel.size());
        simTree.getVel(v);
        vel = toRVec(v);
    }
    void getAcc(RVec& acc) const {
        Vector a(acc.size());
        simTree.getAcc(a);
        acc = toRVec(a);
    }
    
    void velFromCartesian(const RVec& pos, RVec& vel) {
        Vector v(vel.size());
        simTree.velFromCartesian(toVector(pos),v);
        vel = toRVec(v);
    }
    void enforceTreeConstraints(RVec& pos, RVec& vel) {
        Vector p(pos.size()), v(vel.size());
        p = toVector(pos); v = toVector(vel);
        simTree.enforceTreeConstraints(p,v);
        pos = toRVec(p); vel = toRVec(v);
    }
    void enforceConstraints(RVec& pos, RVec& vel) {
        Vector p(pos.size()), v(vel.size());
        p = toVector(pos); v = toVector(vel);
        simTree.enforceConstraints(p,v);
        pos = toRVec(p); vel = toRVec(v);
    }
    void prepareForDynamics() {
        simTree.prepareForDynamics();
    }
    void calcTreeForwardDynamics(const CDSVecVec6& spatialForces) {
        const Vector_<SpatialVec> svl = toSpatialVecList(spatialForces);
        simTree.calcTreeForwardDynamics(svl);
    }
    void calcLoopForwardDynamics(const CDSVecVec6& spatialForces) {
        const Vector_<SpatialVec> svl = toSpatialVecList(spatialForces);
        simTree.calcLoopForwardDynamics(svl);
    }

    void calcP() {simTree.calcP();}
    void calcZ(const CDSVecVec6& spatialForces) {
        const Vector_<SpatialVec> svl = toSpatialVecList(spatialForces);
        simTree.calcZ(svl);
    }
    void calcTreeAccel()    {simTree.calcTreeAccel();}
    void fixVel0(RVec& vel) {
        Vector v(vel.size());
        v = toVector(vel);
        simTree.fixVel0(v);
        vel = toRVec(v);
    }
    void calcY()            {simTree.calcY();}
    void calcTreeInternalForces(const CDSVecVec6& spatialForces) {
        const Vector_<SpatialVec> svl = toSpatialVecList(spatialForces);
        simTree.calcTreeInternalForces(svl);
    }
    void getInternalForces(RVec& T) {
        Vector t(T.size());
        simTree.getInternalForces(t);
        T = toRVec(t);
    }

    void getConstraintCorrectedInternalForces(RVec& T) {
        Vector t(T.size());
        simTree.getConstraintCorrectedInternalForces(t);
        T = toRVec(t);
    }
private:
    RigidBodyTree simTree;
};


/*static*/ IVMMoleculeRBTreeInterface*
IVMMoleculeRBTreeInterface::create(bool newStyle) {
    return newStyle ? (IVMMoleculeRBTreeInterface*)new SimbodyRigidBodyTreeMolIfc()
                    : (IVMMoleculeRBTreeInterface*)new IVMRigidBodyTreeMolIfc();
}


