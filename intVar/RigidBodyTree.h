#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "cdsList.h"
#include "cdsVector.h"
#include "vec3.h"

#include "RigidBodyNode.h"

#include <cassert>

typedef CDSList<RigidBodyNode*>   RBNodePtrList;
typedef CDSVector<double,1>       RVec;   // first element has index 1
typedef CDSList<Vec6>             VecVec6;

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
    RBStation(RigidBodyNode& n, const Vec3& pos) : rbNode(&n), station_B(pos) { }
    // default copy, assignment, destructor

    void calcPosInfo(RBStationRuntime&) const;
    void calcVelInfo(RBStationRuntime&) const;
    void calcAccInfo(RBStationRuntime&) const;

    RigidBodyNode&       getNode()    const { assert(isValid()); return *rbNode; }
    const Vec3&          getStation() const { assert(isValid()); return station_B; }
    bool                 isValid()    const { return rbNode != 0; }
private:
    RigidBodyNode*       rbNode;
    Vec3                 station_B;
};
ostream& operator<<(ostream&, const RBStation&);

class RBStationRuntime {
public:
    Vec3 station_G;    // vector from body origin OB to station, reexpressed in G
    Vec3 stationVel_G; // velocity of station relative to velocity of OB, expr. in G
    Vec3 pos_G;        // spatial quantities
    Vec3 vel_G;
    Vec3 acc_G;

    Vec3 force_G;      // the constraint force (calculated)
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

    Vec3 fromTip1ToTip2_G;    // tip2.pos - tip1.pos
    Vec3 unitDirection_G;     // fromTip1ToTip2/|fromTip1ToTip2|

    Vec3 relVel_G;            // spatial relative velocity tip2.vel-tip1.vel

    Real posErr;
    Real velErr;
    Real accErr;
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
    RigidBodyTree() : lConstraints(0) { }
    ~RigidBodyTree();

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addRigidBodyNode(RigidBodyNode&  parent,
                         const Frame&    referenceConfig,    // body frame in parent
                         RigidBodyNode*& nodep);

    /// Same as addRigidBodyNode but special-cased for ground.
    int addGroundNode(RigidBodyNode*& gnodep);

    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    int addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d);

    /// Call this after all bodies & constraints have been added.
    /// TODO: this "ivm" has to go.
    void finishConstruction(IVM* ivm);

    // deallocate subtree rooted at the indicated node
    void destructNode(RigidBodyNode*); 

    // includes ground
    int getNBodies() const { return nodeNum2NodeMap.size(); }

//    int getLevel(int nodeNum) const { return getRigidBodyNode(nodeNum)->getLevel(); }

    int getDOF() const; 
    int getDim() const; 

    // Kinematics -- calculate spatial quantities from internal states.
    void setPos(const RVec& pos);
    void setVel(const RVec& vel);

    void getPos(RVec& pos) const;
    void getVel(RVec& vel) const;
    void getAcc(RVec& acc) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const RVec& pos, RVec& vel);

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints (just quaternions constraints; ignores loops).
    void enforceTreeConstraints(RVec& pos, RVec& vel);

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaterion constraints).
    void enforceConstraints(RVec& pos, RVec& vel);

    /// Prepare for dynamics by calculating position-dependent quantities
    /// like the articulated body inertias P.
    void prepareForDynamics();

    /// Given a set of spatial forces, calculate accelerations ignoring
    /// constraints. Must have already called prepareForDynamics().
    /// TODO: also applies stored internal forces (hinge torques) which
    /// will cause surprises if non-zero.
    void calcTreeForwardDynamics(const VecVec6& spatialForces);

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    void calcLoopForwardDynamics(const VecVec6& spatialForces);


    /// Unconstrained (tree) dynamics 
    void calcP();                             // articulated body inertias
    void calcZ(const VecVec6& spatialForces); // articulated body remainder forces
    void calcTreeAccel();                     // accels with forces from last calcZ

    void fixAccelForConstraints();            // call after calcTreeAccel

    void fixVel0(RVec& vel); // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY();

    /// Convert spatial forces to internal (joint) forces, ignoring constraints.
    void calcTreeInternalForces(const VecVec6& spatialForces);

    /// Retrieve last-computed internal (joint) forces.
    void getInternalForces(RVec& T);

    void getConstraintCorrectedInternalForces(RVec& T); // TODO has to move elsewhere

    void propagateSVel();

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
