#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "simbody/internal/SimbodyCommon.h"
using namespace simtk;

#include "cdsList.h"
#include "cdsVec3.h"

#include "RigidBodyNode.h"

#include <cassert>

typedef CDSList<RigidBodyNode*>   RBNodePtrList;
typedef CDSList<CDSVec6>          CDSVecVec6;

class IVM;
class LengthConstraints;
class RBStationRuntime;
class RBDistanceConstraintRuntime;

/**
 * A station is a point located on a particular rigid body. A station is
 * measured from the body frame origin and expressed in the body frame.
 */
class RBStation {
public:
    RBStation() : rbNode(0) { } // so we can have arrays of these
    RBStation(RigidBodyNode& n, const CDSVec3& pos) : rbNode(&n), station_B(pos) { }
    // default copy, assignment, destructor

    void calcPosInfo(RBStationRuntime&) const;
    void calcVelInfo(RBStationRuntime&) const;
    void calcAccInfo(RBStationRuntime&) const;

    RigidBodyNode&       getNode()    const { assert(isValid()); return *rbNode; }
    const CDSVec3&       getStation() const { assert(isValid()); return station_B; }
    bool                 isValid()    const { return rbNode != 0; }
private:
    RigidBodyNode*       rbNode;
    CDSVec3              station_B;
};
ostream& operator<<(ostream&, const RBStation&);

class RBStationRuntime {
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
class RBDistanceConstraint {
public:
    RBDistanceConstraint() : distance(-1.), runtimeIndex(-1) {}
    RBDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d) {
        assert(s1.isValid() && s2.isValid() && d >= 0.);
        stations[0] = s1; stations[1] = s2; distance = d;
        runtimeIndex = -1;
    }

    void calcPosInfo(RBDistanceConstraintRuntime&) const;
    void calcVelInfo(RBDistanceConstraintRuntime&) const;
    void calcAccInfo(RBDistanceConstraintRuntime&) const;

    void setRuntimeIndex(int ix) {assert(ix>=0); runtimeIndex=ix;}
    int  getRuntimeIndex() const {assert(isValid()&&runtimeIndex>=0); return runtimeIndex;}

    const double&    getDistance()     const { return distance; }
    const RBStation& getStation(int i) const { assert(isValid() && (i==1||i==2)); return stations[i-1]; }
    bool             isValid()         const { return distance >= 0.; }

protected:
    double       distance;
    RBStation    stations[2];
    int          runtimeIndex;
};

class RBDistanceConstraintRuntime {
public:
    RBDistanceConstraintRuntime() { }

    RBStationRuntime stationRuntimes[2];

    CDSVec3 fromTip1ToTip2_G;    // tip2.pos - tip1.pos
    CDSVec3 unitDirection_G;     // fromTip1ToTip2/|fromTip1ToTip2|

    CDSVec3 relVel_G;            // spatial relative velocity tip2.vel-tip1.vel

    CDSReal posErr;
    CDSReal velErr;
    CDSReal accErr;
};

/**
 * The RigidBodyTree class owns the tree of joint-connected rigid bodies, called
 * RigidBodyNodes. The tree is stored by levels, with level 0 being ground, level 1
 * being bodies which are connected to ground (base bodies), level 2 connected to
 * level 1 and so on. Nodes at the same level are stored together in an array,
 * but the order does not reflect the logical tree structure; that is maintained
 * via parent & children pointers kept in the nodes.
 * 
 * RigidBodyTree is the owner of the RigidBodyNode objects (which are abstract), pointers to
 * which are stored in the tree.
 */
class RigidBodyTree {
public:
    RigidBodyTree() : lConstraints(0), DOFTotal(-1), dimTotal(-1) { }
    ~RigidBodyTree();

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addRigidBodyNode(RigidBodyNode&  parent,
                         const RBFrame&  referenceConfig,    // body frame in parent
                         RigidBodyNode*& nodep);

    /// Same as addRigidBodyNode but special-cased for ground.
    int addGroundNode(RigidBodyNode*& gnodep);

    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    int addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d);

    /// Call this after all bodies & constraints have been added.
    void finishConstruction(const double& ctol, int verbose);

    // includes ground
    int getNBodies() const { return nodeNum2NodeMap.size(); }

    int getDOF() const { return DOFTotal; } 
    int getDim() const { return dimTotal; } 

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const Vector& pos);
    void setVel(const Vector& vel);

    void getPos(Vector& pos) const;
    void getVel(Vector& vel) const;
    void getAcc(Vector& acc) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel);

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints (just quaternions constraints; ignores loops).
    void enforceTreeConstraints(Vector& pos, Vector& vel);

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaternion constraints).
    void enforceConstraints(Vector& pos, Vector& vel);

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

    void fixVel0(Vector& vel); // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY();

    /// Convert spatial forces to internal (joint) forces, ignoring constraints.
    void calcTreeInternalForces(const CDSVecVec6& spatialForces);

    /// Retrieve last-computed internal (joint) forces.
    void getInternalForces(Vector& T);

    void getConstraintCorrectedInternalForces(Vector& T); // TODO has to move elsewhere

    const RigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }
    RigidBodyNode& updRigidBodyNode(int nodeNum) {
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
    CDSList<RBNodePtrList> rbNodeLevels;
    // Map nodeNum to (level,offset).
    CDSList<RigidBodyNodeIndex> nodeNum2NodeMap;

    CDSList<RBDistanceConstraint>        distanceConstraints;
    // TODO: later this moves to state cache (sherm)
    CDSList<RBDistanceConstraintRuntime> dcRuntimeInfo;
    
    LengthConstraints* lConstraints;
    friend ostream& operator<<(ostream&, const RigidBodyTree&);
};

ostream& operator<<(ostream&, const RigidBodyTree&);

#endif /* RIGID_BODY_TREE_H_ */
