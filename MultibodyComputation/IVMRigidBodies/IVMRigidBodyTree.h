#ifndef IVM_RIGID_BODY_TREE_H_
#define IVM_RIGID_BODY_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"
#include "cdsVec3.h"

#include "IVMRigidBodyNode.h"

#include <cassert>

typedef CDSList<IVMRigidBodyNode*>   IVMNodePtrList;
typedef CDSVector<double,1>          RVec;   // first element has index 1
typedef CDSList<CDSVec6>             CDSVecVec6;

class IVM;
class IVMLengthConstraints;
class IVMStationRuntime;
class IVMDistanceConstraintRuntime;

/**
 * A station is a point located on a particular rigid body. A station is
 * measured from the body frame origin and expressed in the body frame.
 */
class IVMStation {
public:
    IVMStation() : rbNode(0) { } // so we can have arrays of these
    IVMStation(IVMRigidBodyNode& n, const CDSVec3& pos) : rbNode(&n), station_B(pos) { }
    // default copy, assignment, destructor

    void calcPosInfo(IVMStationRuntime&) const;
    void calcVelInfo(IVMStationRuntime&) const;
    void calcAccInfo(IVMStationRuntime&) const;

    IVMRigidBodyNode&    getNode()    const { assert(isValid()); return *rbNode; }
    const CDSVec3&       getStation() const { assert(isValid()); return station_B; }
    bool                 isValid()    const { return rbNode != 0; }
private:
    IVMRigidBodyNode*    rbNode;
    CDSVec3              station_B;
};
ostream& operator<<(ostream&, const IVMStation&);

class IVMStationRuntime {
public:
    CDSVec3 station_G;    // vector from body origin OB to station, reexpressed in G
    CDSVec3 stationVel_G; // velocity of station relative to velocity of OB, expr. in G
    CDSVec3 pos_G;        // spatial quantities
    CDSVec3 vel_G;
    CDSVec3 acc_G;

    CDSVec3 force_G;      // the constraint force (calculated)
};

/**
 * This class requests that two stations, one on each of two rigid bodies,
 * be maintained at a certain separation distance at all times.
 */
class IVMDistanceConstraint {
public:
    IVMDistanceConstraint() : distance(-1.), runtimeIndex(-1) {}
    IVMDistanceConstraint(const IVMStation& s1, const IVMStation& s2, const double& d) {
        assert(s1.isValid() && s2.isValid() && d >= 0.);
        stations[0] = s1; stations[1] = s2; distance = d;
        runtimeIndex = -1;
    }

    void calcPosInfo(IVMDistanceConstraintRuntime&) const;
    void calcVelInfo(IVMDistanceConstraintRuntime&) const;
    void calcAccInfo(IVMDistanceConstraintRuntime&) const;

    void setRuntimeIndex(int ix) {assert(ix>=0); runtimeIndex=ix;}
    int  getRuntimeIndex() const {assert(isValid()&&runtimeIndex>=0); return runtimeIndex;}

    const double&     getDistance()     const { return distance; }
    const IVMStation& getStation(int i) const { assert(isValid() && (i==1||i==2)); return stations[i-1]; }
    bool              isValid()         const { return distance >= 0.; }

protected:
    double       distance;
    IVMStation   stations[2];
    int          runtimeIndex;
};

class IVMDistanceConstraintRuntime {
public:
    IVMDistanceConstraintRuntime() { }

    IVMStationRuntime stationRuntimes[2];

    CDSVec3 fromTip1ToTip2_G;    // tip2.pos - tip1.pos
    CDSVec3 unitDirection_G;     // fromTip1ToTip2/|fromTip1ToTip2|

    CDSVec3 relVel_G;            // spatial relative velocity tip2.vel-tip1.vel

    CDSReal posErr;
    CDSReal velErr;
    CDSReal accErr;
};

/**
 * The IVMRigidBodyTree class owns the tree of joint-connected rigid bodies, called
 * IVMRigidBodyNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * IVMRigidBodyTree is the owner of the IVMRigidBodyNode objects (which are abstract),
 * pointers to which are stored in the tree.
 */
class IVMRigidBodyTree {
public:
    IVMRigidBodyTree() : lConstraints(0), DOFTotal(-1), dimTotal(-1) { }
    ~IVMRigidBodyTree();

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addRigidBodyNode(IVMRigidBodyNode&  parent,
                         const IVMFrame&    referenceConfig,    // body frame in parent
                         IVMRigidBodyNode*& nodep);

    /// Same as addRigidBodyNode but special-cased for ground.
    int addGroundNode(IVMRigidBodyNode*& gnodep);

    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    int addDistanceConstraint(const IVMStation& s1, const IVMStation& s2, const double& d);

    /// Call this after all bodies & constraints have been added.
    void finishConstruction(const double& ctol, int verbose);

    // includes ground
    int getNBodies() const { return nodeNum2NodeMap.size(); }

    int getDOF() const { return DOFTotal; } 
    int getDim() const { return dimTotal; } 

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos);
    void setVel(const RVec& vel);

    void getPos(RVec& pos) const;
    void getVel(RVec& vel) const;
    void getAcc(RVec& acc) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const RVec& pos, RVec& vel) {assert(false);/*TODO*/}

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints (just quaternions constraints; ignores loops).
    void enforceTreeConstraints(RVec& pos, RVec& vel);

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaternion constraints).
    void enforceConstraints(RVec& pos, RVec& vel);

    /// Prepare for dynamics by calculating position-dependent quantities
    /// like the articulated body inertias P.
    void prepareForDynamics();

    /// Given a set of spatial forces, calculate accelerations ignoring
    /// constraints. Must have already called prepareForDynamics().
    /// TODO: also applies stored internal forces (hinge torques) which
    /// will cause surprises if non-zero.
    void calcTreeForwardDynamics(const CDSVecVec6& spatialForces);

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    void calcLoopForwardDynamics(const CDSVecVec6& spatialForces);


    /// Unconstrained (tree) dynamics 
    void calcP();                             // articulated body inertias
    void calcZ(const CDSVecVec6& spatialForces); // articulated body remainder forces
    void calcTreeAccel();                     // accels with forces from last calcZ

    void fixVel0(RVec& vel); // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY();

    /// Convert spatial forces to internal (joint) forces, ignoring constraints.
    void calcTreeInternalForces(const CDSVecVec6& spatialForces);

    /// Retrieve last-computed internal (joint) forces.
    void getInternalForces(RVec& T);

    void getConstraintCorrectedInternalForces(RVec& T); // TODO has to move elsewhere

    const IVMRigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }
    IVMRigidBodyNode& updRigidBodyNode(int nodeNum) {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

private:
    struct RigidBodyNodeIndex {
        RigidBodyNodeIndex(int l, int o) : level(l), offset(o) { }
        int level, offset;
    };

    // set by finishConstruction
    int DOFTotal;   // summed over all nodes
    int dimTotal;

    // This holds pointers to nodes and serves to map (level,offset) to nodeSeqNo.
    CDSList<IVMNodePtrList> rbNodeLevels;
    // Map nodeNum to (level,offset).
    CDSList<RigidBodyNodeIndex> nodeNum2NodeMap;

    CDSList<IVMDistanceConstraint>        distanceConstraints;
    // TODO: later this moves to state cache (sherm)
    CDSList<IVMDistanceConstraintRuntime> dcRuntimeInfo;
    
    IVMLengthConstraints* lConstraints;
    friend ostream& operator<<(ostream&, const IVMRigidBodyTree&);
};

ostream& operator<<(ostream&, const IVMRigidBodyTree&);

#endif // IVM_RIGID_BODY_TREE_H_
