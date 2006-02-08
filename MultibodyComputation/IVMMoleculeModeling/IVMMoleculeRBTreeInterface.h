#ifndef IVM_MOLECULE_RB_TREE_INTERFACE_H_
#define IVM_MOLECULE_RB_TREE_INTERFACE_H_

/**@file
 * Sherm 060126: this is an abstraction of the IVM rigid body tree interface.
 * The idea is to have the original code implement the interface in as
 * straightforward a way possible, and make sure the code still works. Then,
 * implement the same interface using the new Simbody version of this code 
 * (which will require much translation between vector formats and so on)
 * to show that it still yields the original results.
 */


#include "cdsList.h"
#include "cdsVector.h"
#include "cdsVec3.h"
#include "IVMMassProperties.h"

typedef FixedVector<double,6>        CDSVec6;
typedef CDSVector<double,1>          RVec;   // first element has index 1
typedef CDSList<CDSVec6>             CDSVecVec6;

/**
 * This class provides as virtual methods all the operations of the IVMRigidBodyTree,
 * plus those needed for manipulation the IVMRigidBodyNodes as well, without exposing
 * any internals.
 */
class IVMMoleculeRBTreeInterface {
public:
    IVMMoleculeRBTreeInterface() { }
    virtual ~IVMMoleculeRBTreeInterface() { }

    static IVMMoleculeRBTreeInterface* create(bool newStyle);
    virtual IVMMoleculeRBTreeInterface* clone() const=0;
    virtual int getInitialStateOffset() const=0; // either 0 or 1

    /// Create a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.

    virtual int addRigidBodyNode(int parentNodeNum,
                         const IVMFrame&          refConfig,    // body frame in parent
                         const IVMMassProperties& massProps,    // mass properties in body frame
                         const IVMFrame&          jointFrame,   // inboard joint frame J in body frame
                         const IVMJointType&      jtype,
                         bool                     isReversed,
                         bool                     useEuler,
                         int&                     nxtStateOffset)=0;

    /// Same as addRigidBodyNode but special-cased for ground.
    virtual int addGroundNode()=0;

    virtual void    RBNodeSetVelFromSVel(int nodeNum, const CDSVec6& sVel)=0;
    virtual void    RBNodeGetConfig(int nodeNum, CDSMat33& R_GB, CDSVec3& OB_G) const=0;
    virtual CDSVec6 RBNodeGetSpatialVel(int nodeNum) const=0;
    virtual double  RBNodeGetMass(int nodeNum) const=0;
    virtual CDSVec3 RBNodeGetCOM_G(int nodeNum) const=0;
    virtual double  RBNodeCalcKineticEnergy(int nodeNum) const=0;
    virtual void    RBNodeDump(int nodeNum, std::ostream&) const=0;

    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    virtual int addDistanceConstraint(int nodeNum1, const CDSVec3& station1,
                                      int nodeNum2, const CDSVec3& station2,
                                      double distance) = 0;

    /// Call this after all bodies & constraints have been added.
    virtual void finishConstruction(const double& ctol, int verbose)=0;

    // includes ground
    virtual int getNBodies() const=0;

    virtual int getDOF() const=0;
    virtual int getDim() const=0;

    // Kinematics -- calculate spatial quantities from internal states.
    virtual void setPos(const RVec& pos)=0;
    virtual void setVel(const RVec& vel)=0;

    virtual void getPos(RVec& pos) const=0;
    virtual void getVel(RVec& vel) const=0;
    virtual void getAcc(RVec& acc) const=0;
    
    /// This is a solver which generates internal velocities from spatial ones.
    virtual void velFromCartesian(const RVec& pos, RVec& vel)=0;

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints (just quaternions constraints; ignores loops).
    virtual void enforceTreeConstraints(RVec& pos, RVec& vel)=0;

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaternion constraints).
    virtual void enforceConstraints(RVec& pos, RVec& vel)=0;

    /// Prepare for dynamics by calculating position-dependent quantities
    /// like the articulated body inertias P.
    virtual void prepareForDynamics()=0;

    /// Given a set of spatial forces, calculate accelerations ignoring
    /// constraints. Must have already called prepareForDynamics().
    /// TODO: also applies stored internal forces (hinge torques) which
    /// will cause surprises if non-zero.
    virtual void calcTreeForwardDynamics(const CDSVecVec6& spatialForces)=0;

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    virtual void calcLoopForwardDynamics(const CDSVecVec6& spatialForces)=0;


    /// Unconstrained (tree) dynamics 
    virtual void calcP()=0;                             // articulated body inertias
    virtual void calcZ(const CDSVecVec6& spatialForces)=0; // articulated body remainder forces
    virtual void calcTreeAccel()=0;                     // accels with forces from last calcZ

    virtual void fixVel0(RVec& vel)=0; // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    virtual void calcY()=0;

    /// Convert spatial forces to internal (joint) forces, ignoring constraints.
    virtual void calcTreeInternalForces(const CDSVecVec6& spatialForces)=0;

    /// Retrieve last-computed internal (joint) forces.
    virtual void getInternalForces(RVec& T)=0;

    virtual void getConstraintCorrectedInternalForces(RVec& T)=0; // TODO has to move elsewhere

    friend ostream& operator<<(ostream&, const IVMMoleculeRBTreeInterface&);
};

ostream& operator<<(ostream&, const IVMMoleculeRBTreeInterface&);

#endif // IVM_MOLECULE_RB_TREE_INTERFACE_H_
