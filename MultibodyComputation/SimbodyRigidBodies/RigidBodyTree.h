#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
using namespace simtk;

#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"

#include <cassert>
#include <vector>
#include <iostream>

typedef std::vector<RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>           SpatialVecList;

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

    RigidBodyNode&    getNode()    const { assert(isValid()); return *rbNode; }
    const Vec3&       getStation() const { assert(isValid()); return station_B; }
    bool              isValid()    const { return rbNode != 0; }
private:
    RigidBodyNode*    rbNode;
    Vec3              station_B;
};
std::ostream& operator<<(std::ostream&, const RBStation&);

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
    RigidBodyTree() 
      : nextUSlot(0), nextQSlot(0), DOFTotal(-1), SqDOFTotal(-1), maxNQTotal(-1), 
        built(false), lConstraints(0) 
      { addGroundNode(); }

    RigidBodyTree(const RigidBodyTree&);
    RigidBodyTree& operator=(const RigidBodyTree&);
    ~RigidBodyTree();

    /// Take ownership of a new node, add it to the tree, and assign it
    /// a node number. This is NOT a regular labeling; it is just
    /// for reference. You can depend on nodeNum being (a) unique, and (b) a
    /// small enough integer to make it a reasonable index, but don't depend
    /// on it having any particular value or being sequential or even
    /// monotonically increasing.
    int addRigidBodyNode(RigidBodyNode&      parent,
                         const TransformMat& referenceConfig, // body frame in parent
                         RigidBodyNode*&     nodep);


    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    int addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d);

    /// Call this after all bodies & constraints have been added.
    void realizeConstruction(const double& ctol, int verbose);

        // CALLABLE AFTER realizeConstruction()

    const SBState& getDefaultState() const {return defaultState;}

    // includes ground
    int getNBodies()      const {assert(built); return nodeNum2NodeMap.size();}
    int getNConstraints() const {assert(built); return distanceConstraints.size();}

    int getTotalDOF()    const {assert(built); return DOFTotal;}
    int getTotalSqDOF()  const {assert(built); return SqDOFTotal;}
    int getTotalQAlloc() const {assert(built); return maxNQTotal;}

    int getTotalMultAlloc() const {assert(built); assert(false); return -1;} // TODO

    int getQIndex(int body) const {assert(built);return getRigidBodyNode(body).getQIndex();}
    int getQAlloc(int body) const {assert(built);return getRigidBodyNode(body).getMaxNQ();}
    int getUIndex(int body) const {assert(built);return getRigidBodyNode(body).getUIndex();}
    int getDOF   (int body) const {assert(built);return getRigidBodyNode(body).getDOF();}

    int getMultIndex(int constraint) const {assert(built);assert(false); return -1;} //TODO
    int getMaxNMult (int constraint) const {assert(built);assert(false); return -1;} //TODO


    // Modeling info.
    void setDefaultModelingValues(SBState& s) {
        assert(s.getStage() >= BuiltStage);
        s.cache->stage = BuiltStage; // trim back if necessary
        s.vars->useEulerAngles = false;
        s.vars->prescribed.assign(getNBodies(), false);
        s.vars->enabled.assign(getNConstraints(), false);
    }
    void setUseEulerAngles(SBState& s, bool useAngles) const {
        assert(s.getStage() >= BuiltStage);
        if (s.getStage() >= ModeledStage && (s.vars->useEulerAngles == useAngles))
            return; // no change
        s.cache->stage = BuiltStage; // back up if necessary
        s.vars->useEulerAngles = useAngles;
    }
    void setJointIsPrescribed(SBState& s, int joint, bool prescribe) const {
        assert(s.getStage() >= BuiltStage);
        if (s.getStage() >= ModeledStage && (s.vars->prescribed[joint] == prescribe))
            return; // no change
        s.cache->stage = BuiltStage; // back up if necessary
        s.vars->prescribed[joint] = prescribe;
    }
    void setConstraintIsEnabled(SBState& s, int constraint, bool enable) const {
        assert(s.getStage() >= BuiltStage);
        if (s.getStage() >= ModeledStage && (s.vars->enabled[constraint] == enable))
            return; // no change
        s.cache->stage = BuiltStage; // back up if necessary
        s.vars->enabled[constraint] = enable;   
    }

    bool getUseEulerAngles(const SBState& s) const {
        assert(s.getStage() >= ModeledStage);
        return s.vars->useEulerAngles;
    }
    bool isJointPrescribed(const SBState& s, int joint) const {
        assert(s.getStage() >= ModeledStage);
        return s.vars->prescribed[joint];
    }
    bool isConstraintEnabled(const SBState& s, int constraint) const {
        assert(s.getStage() >= ModeledStage);
        return s.vars->enabled[constraint];
    }

    /// Call this after all modeling choices have been made, such as whether
    /// to use quaternions, what joints to prescribe, etc.
    void realizeModeling(const SBState&) const;

        // CALLABLE AFTER realizeModeling()

    /// Call this after all available parameters have been given values,
    /// e.g. body masses.
    void realizeParameters(const SBState&) const;


    // Kinematics -- calculate spatial quantities from internal states.
    void realizeConfiguration(const SBState&);
    void realizeMotion(const SBState&);

    // Dynamics -- calculate accelerations and internal forces from 
    // forces and prescribed accelerations supplied in the State.
    void realizeReaction(const SBState&) const;


    void getDefaultParameters   (SBState&) const;
    void getDefaultConfiguration(SBState&) const;
    void getDefaultVelocity     (SBState&) const;


    void getAcc(Vector& acc) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

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
    void calcTreeForwardDynamics(const SpatialVecList& spatialForces);

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    void calcLoopForwardDynamics(const SpatialVecList& spatialForces);


    /// Unconstrained (tree) dynamics 
    void calcP();                             // articulated body inertias
    void calcZ(const SpatialVecList& spatialForces); // articulated body remainder forces
    void calcTreeAccel();                     // accels with forces from last calcZ

    void fixVel0(Vector& vel); // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY();

    /// Convert spatial forces to internal (joint) forces, ignoring constraints.
    void calcTreeInternalForces(const SpatialVecList& spatialForces);

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

    // Initialize to 0 at beginning of construction. These are for doling
    // out Q & U state variables to the nodes.
    int nextUSlot;
    int nextQSlot;

    // set by realizeConstruction
    int DOFTotal;   // summed over all nodes
    int SqDOFTotal; // sum of squares of ndofs per node
    int maxNQTotal; // sum of dofs with room for quaternions

    bool built;

    // This is a complete State available immediately after realizeConstruction().
    // It contains default Modeling values, and everything else is allocated in
    // accordance with those. It has been realized through ModelingStage.
    SBState defaultState;

    // This holds pointers to nodes and serves to map (level,offset) to nodeSeqNo.
    std::vector<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

    std::vector<RBDistanceConstraint>        distanceConstraints;
    // TODO: later this moves to state cache (sherm)
    std::vector<RBDistanceConstraintRuntime> dcRuntimeInfo;
    
    LengthConstraints* lConstraints;

    void addGroundNode();
    friend std::ostream& operator<<(std::ostream&, const RigidBodyTree&);
    friend class SimbodyTree;
};

std::ostream& operator<<(std::ostream&, const RigidBodyTree&);

#endif // RIGID_BODY_TREE_H_
