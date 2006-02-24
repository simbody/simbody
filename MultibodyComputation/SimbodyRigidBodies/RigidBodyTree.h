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

    void calcPosInfo(const SBState&, RBStationRuntime&) const;
    void calcVelInfo(const SBState&, RBStationRuntime&) const;
    void calcAccInfo(const SBState&, RBStationRuntime&) const;

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

    void calcPosInfo(const SBState&, RBDistanceConstraintRuntime&) const;
    void calcVelInfo(const SBState&, RBDistanceConstraintRuntime&) const;
    void calcAccInfo(const SBState&, RBDistanceConstraintRuntime&) const;

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
      : nextUSlot(0), nextUSqSlot(0), nextQSlot(0), DOFTotal(-1), SqDOFTotal(-1), maxNQTotal(-1), 
        built(false), lConstraints(0) 
      { addGroundNode(); }

    RigidBodyTree(const RigidBodyTree&);
    RigidBodyTree& operator=(const RigidBodyTree&);
    ~RigidBodyTree();

    /// Create a new node, add it to the tree, and assign it
    /// a node number, which is a regular labeling starting with node 0 which is ground.
    int addRigidBodyNode
        (RigidBodyNode&          parent,
         const MassProperties&   m,            // mass properties in body frame
         const TransformMat&     X_PJb,        // parent's frame for attaching this joint
         const TransformMat&     X_BJ,         // inboard joint frame J in body frame
         JointSpecification::JointType        
                                 type,
         bool                    isReversed,   // child-to-parent orientation?
         int&                    nxtU,
         int&                    nxtUSq,
         int&                    nxtQ); 


    /// Add a distance constraint and allocate slots to hold the runtime information for
    /// its stations. Return the assigned distance constraint index for caller's use.
    int addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d);

    /// Call this after all bodies & constraints have been added.
    void realizeConstruction(const double& ctol, int verbose);

    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // includes ground
    int getNBodies()      const {return nodeNum2NodeMap.size();}
    int getNConstraints() const {return distanceConstraints.size();}

        // CALLABLE AFTER realizeConstruction()

    const SBState& getDefaultState() const {return defaultState;}


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
        s.vars->prescribed[0] = true; // ground
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

    void setQ(SBState&, const Vector& q) const;
    void setU(SBState&, const Vector& u) const;

    void setJointQ(SBState& s, int body, int axis, const Real& r) const {
        assert(s.cache->stage >= ModeledStage);
        if (s.cache->stage > TimedStage)
            s.cache->stage = TimedStage; // back up if necessary

        const RigidBodyNode& n = getRigidBodyNode(body);
        assert(0 <= axis && axis < n.getNQ(s));
        s.vars->q[n.getQIndex()+axis] = r;
    }

    void setJointU(SBState& s, int body, int axis, const Real& r) const {
        assert(s.cache->stage >= ModeledStage);
        if (s.cache->stage > ConfiguredStage)
            s.cache->stage = ConfiguredStage; // back up if necessary

        const RigidBodyNode& n = getRigidBodyNode(body);
        assert(0 <= axis && axis < n.getDOF());
        s.vars->u[n.getUIndex()+axis] = r;
    }

    void setPrescribedUdot(SBState& s, int body, int axis, const Real& r) const {
        assert(s.cache->stage >= ModeledStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary

        const RigidBodyNode& n = getRigidBodyNode(body);
        assert(0 <= axis && axis < n.getDOF());
        s.vars->prescribedUdot[n.getUIndex()+axis] = r;
    }

    const Vector& getQ(const SBState&) const;
    VectorView&   updQ(SBState&)       const;

    const Vector& getU(const SBState&) const;
    VectorView&   updU(SBState&)       const;

    const Vector& getAppliedJointForces(const SBState&) const;
    const Vector_<SpatialVec>& getAppliedBodyForces(const SBState&) const;

    const Vector& getQdot(const SBState&) const;
    const Vector& getUdot(const SBState&) const;
    const Vector& getQdotDot(const SBState&) const;


    // Kinematics -- calculate spatial quantities from internal states.
    void realizeConfiguration(const SBState&) const;
    void realizeMotion(const SBState&) const;

    // Dynamics -- calculate accelerations and internal forces from 
    // forces and prescribed accelerations supplied in the State.
    void realizeReaction(const SBState&) const;

    void clearAppliedForces(SBState& s) const {
        assert(s.cache->stage >= ModeledStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary
        s.vars->appliedJointForces.setToZero();
        s.vars->appliedBodyForces.setToZero();
    }

    void applyGravity(SBState& s, const Vec3& g) const {
        assert(s.cache->stage >= ConfiguredStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary
        for (int body=1; body<getNBodies(); ++body) {
            const RigidBodyNode& n = getRigidBodyNode(body);
            applyPointForce(s, body, n.getCOM_B(s), n.getMass(s)*g);
        }
    }

    void applyPointForce(SBState& s, int body, const Vec3& stationInB, 
                         const Vec3& forceInG) const
    {
        assert(s.cache->stage >= ConfiguredStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary

        const RotationMat& R_GB = getRigidBodyNode(body).getX_GB(s).R();
        s.vars->appliedBodyForces[body] += 
            SpatialVec((R_GB*stationInB) % forceInG, forceInG);
    }

    void applyBodyTorque(SBState& s, int body, const Vec3& torqueInG) const {
        assert(s.cache->stage >= ConfiguredStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary

        s.vars->appliedBodyForces[body] += SpatialVec(torqueInG, Vec3(0));
    }

    void applyJointForce(SBState& s, int body, int axis, const Real& r) const {
        assert(s.cache->stage >= ConfiguredStage);
        if (s.cache->stage > MovingStage)
            s.cache->stage = MovingStage; // back up if necessary

        const RigidBodyNode& n = getRigidBodyNode(body);
        assert(0 <= axis && axis < n.getDOF());
        s.vars->appliedJointForces[n.getUIndex()+axis] = r;
    }

    void getDefaultParameters   (SBState&) const;
    void getDefaultConfiguration(SBState&) const;
    void getDefaultVelocity     (SBState&) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    /// This is a solver which tweaks the state to make it satisfy position
    /// and velocity constraints (just quaternions constraints; ignores loops).
    void enforceQuaternionConstraints(SBState&) const;

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaternion constraints).
    void enforceLengthConstraints(SBState&) const;

    /// Prepare for dynamics by calculating position-dependent quantities
    /// like the articulated body inertias P.
    void prepareForDynamics(const SBState&) const;

    /// Given a set of spatial forces, calculate accelerations ignoring
    /// constraints. Must have already called prepareForDynamics().
    /// TODO: also applies stored internal forces (hinge torques) which
    /// will cause surprises if non-zero.
    void calcTreeForwardDynamics(const SBState&, 
                                 const SpatialVecList& spatialForces) const;

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    void calcLoopForwardDynamics(const SBState&, 
                                 const SpatialVecList& spatialForces) const;


    /// Unconstrained (tree) dynamics 
    void calcP(const SBState&) const;                        // articulated body inertias
    void calcZ(const SBState&, const SpatialVecList& spatialForces) const; // articulated body remainder forces
    void calcTreeAccel(const SBState&) const;                // accels with forces from last calcZ

    void fixVel0(SBState&, Vector& vel) const; // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY(const SBState&) const;

    /// Calculate the product J*X where J is the partial velocity Jacobian dV/du
    /// and X is a vector of SpatialVec's, one per body. See Eq. 76&77 in
    /// Schwieters' paper, and see 81a & b for a use of this routine to compute
    /// energy gradient in internal coordinates. In that case X=dE/dR, that is
    /// the gradient of the energy w.r.t. atomic positions, summed and shifted
    /// to body origins. There we are pretending dR/dq is the same as dV/du, which
    /// will be true if dq/dt = u, which works for all cases except quaternions.
    /// Schwieters handles that by using Euler angles for orientation coordinates
    /// when doing minimizations. But note that the routine works in terms of u,
    /// not q, so it produces a meaningful result in all cases, just not one that
    /// can be mapped directly back to quaternion coordinates. This is an O(n)
    /// operator which can be called after realizeConfiguration().
    /// It has no effect on the cache.
    void calcInternalGradientFromSpatial(const SBState&, const SpatialVecList& X, 
                                         Vector& JX);

    /// Pass in internal forces in T; they will be adjusted by this routine.
    void calcConstraintCorrectedInternalForces(const SBState&, Vector& T); 

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
    int nextUSqSlot;
    int nextQSlot;

    // set by realizeConstruction
    int DOFTotal;   // summed over all nodes
    int SqDOFTotal; // sum of squares of ndofs per node
    int maxNQTotal; // sum of dofs with room for quaternions

    bool built;

    // This is a complete State available immediately after realizeConstruction().
    // It contains default Modeling values, and everything else is allocated in
    // accordance with those.
    SBState defaultState;

    // This holds pointers to nodes and serves to map (level,offset) to nodeSeqNo.
    std::vector<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

    std::vector<RBDistanceConstraint>        distanceConstraints;
    // TODO: later this moves to state cache (sherm)
    mutable std::vector<RBDistanceConstraintRuntime> dcRuntimeInfo;
    
    LengthConstraints* lConstraints;

    void addGroundNode();
    friend std::ostream& operator<<(std::ostream&, const RigidBodyTree&);
    friend class SimbodyTree;
};

std::ostream& operator<<(std::ostream&, const RigidBodyTree&);

#endif // RIGID_BODY_TREE_H_
