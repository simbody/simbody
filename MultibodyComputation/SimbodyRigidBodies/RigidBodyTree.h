#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "simbody/Simbody.h"
#include "SimbodyTree.h"
using namespace simtk;

class RigidBodyNode;

namespace simtk {
class SBStateRep;
class SBModelingVars;
class SBParameterVars;
class SBTimeVars;
class SBConfigurationVars;
class SBMotionVars;
class SBDynamicVars;
}

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

    void calcPosInfo(const SBStateRep&, RBStationRuntime&) const;
    void calcVelInfo(const SBStateRep&, RBStationRuntime&) const;
    void calcAccInfo(const SBStateRep&, RBStationRuntime&) const;

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

    void calcPosInfo(const SBStateRep&, RBDistanceConstraintRuntime&) const;
    void calcVelInfo(const SBStateRep&, RBDistanceConstraintRuntime&) const;
    void calcAccInfo(const SBStateRep&, RBDistanceConstraintRuntime&) const;

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

    /// This is available any time.
    SBStage getStage(const SBStateRep&) const;

    /// This will realize the state up to the indicated stage, possibly advancing
    /// multiple stages internally by calling each "realizeWhatever()" routine in order.
    void realize(const SBStateRep& s, SBStage stage) const;

    /// Call this after all bodies & constraints have been added.
    void realizeConstruction (); // will set built==true
    void realizeModeling     (const SBStateRep&) const;
    void realizeParameters   (const SBStateRep&) const;
    void realizeTime         (const SBStateRep&) const;
    void realizeConfiguration(const SBStateRep&) const;
    void realizeMotion       (const SBStateRep&) const;
    void realizeReaction     (const SBStateRep&) const;

    const SBState& getInitialState() const {
        assert(built); return initialState;
    }
    void setDefaultModelingValues     (const SBStateRep&, SBModelingVars&)      const;
    void setDefaultParameterValues    (const SBStateRep&, SBParameterVars&)     const;
    void setDefaultTimeValues         (const SBStateRep&, SBTimeVars&)          const;
    void setDefaultConfigurationValues(const SBStateRep&, SBConfigurationVars&) const;
    void setDefaultMotionValues       (const SBStateRep&, SBMotionVars&)        const;
    void setDefaultDynamicValues      (const SBStateRep&, SBDynamicVars&)       const;



    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // includes ground
    int getNBodies()      const {return nodeNum2NodeMap.size();}
    int getNConstraints() const {return distanceConstraints.size();}

        // CALLABLE AFTER realizeConstruction()

    int getTotalDOF()    const {assert(built); return DOFTotal;}
    int getTotalSqDOF()  const {assert(built); return SqDOFTotal;}
    int getTotalQAlloc() const {assert(built); return maxNQTotal;}

    int getTotalMultAlloc() const {assert(built); assert(false); return -1;} // TODO

    int getQIndex(int body) const;
    int getQAlloc(int body) const;
    int getUIndex(int body) const;
    int getDOF   (int body) const;

    int getMultIndex(int constraint) const {assert(built);assert(false); return -1;} //TODO
    int getMaxNMult (int constraint) const {assert(built);assert(false); return -1;} //TODO


    // Modeling info.

    void setUseEulerAngles(SBStateRep& s, bool useAngles) const;
    void setJointIsPrescribed(SBStateRep& s, int joint, bool prescribe) const;
    void setConstraintIsEnabled(SBStateRep& s, int constraint, bool enable) const;
    bool getUseEulerAngles(const SBStateRep& s) const;
    bool isJointPrescribed(const SBStateRep& s, int joint) const;
    bool isConstraintEnabled(const SBStateRep& s, int constraint) const;

        // CALLABLE AFTER realizeModeling()

    void setQ(SBStateRep&, const Vector& q) const;
    void setU(SBStateRep&, const Vector& u) const;

    void setJointQ(SBStateRep& s, int body, int axis, const Real& r) const;
    void setJointU(SBStateRep& s, int body, int axis, const Real& r) const;
    void setPrescribedUdot(SBStateRep& s, int body, int axis, const Real& r) const;

    const Vector& getQ(const SBStateRep&) const;
    Vector&       updQ(SBStateRep&)       const;

    const Vector& getU(const SBStateRep&) const;
    Vector&       updU(SBStateRep&)       const;

    const Vector& getAppliedJointForces(const SBStateRep&) const;
    const Vector_<SpatialVec>& getAppliedBodyForces(const SBStateRep&) const;

    const Vector& getQDot(const SBStateRep&) const;
    const Vector& getUDot(const SBStateRep&) const;
    const Vector& getQDotDot(const SBStateRep&) const;




    // Dynamics -- calculate accelerations and internal forces from 
    // forces and prescribed accelerations supplied in the State.


    void clearAppliedForces(SBStateRep& s) const;
    void applyGravity(SBStateRep& s, const Vec3& g) const;
    void applyPointForce(SBStateRep& s, int body, const Vec3& stationInB, 
                         const Vec3& forceInG) const;
    void applyBodyTorque(SBStateRep& s, int body, const Vec3& torqueInG) const;
    void applyJointForce(SBStateRep& s, int body, int axis, const Real& r) const;
    
    /// This is a solver which generates internal velocities from spatial ones.
    void velFromCartesian(const Vector& pos, Vector& vel) {assert(false);/*TODO*/}

    /// This is a solver which tweaks the state to make it satisfy
    /// quaternions constraints; ignores loops. If any change is made
    /// the stage will be backed up to just before ConfiguredStage.
    /// You should always call realize(ConfiguredStage) after this.
    void enforceQuaternionConstraints(SBStateRep&) const;

    /// This is a solver which tweaks the state to make it satisfy general
    /// constraints (other than quaternion constraints).
    void enforceLengthConstraints(SBStateRep&) const;

    /// Prepare for dynamics by calculating position-dependent quantities
    /// like the articulated body inertias P.
    void prepareForDynamics(const SBStateRep&) const;

    /// Given a set of spatial forces, calculate accelerations ignoring
    /// constraints. Must have already called prepareForDynamics().
    /// TODO: also applies stored internal forces (hinge torques) which
    /// will cause surprises if non-zero.
    void calcTreeForwardDynamics(const SBStateRep&, 
                                 const SpatialVecList& spatialForces) const;

    /// Given a set of spatial forces, calculate acclerations resulting from
    /// those forces and enforcement of acceleration constraints.
    void calcLoopForwardDynamics(const SBStateRep&, 
                                 const SpatialVecList& spatialForces) const;


    /// Unconstrained (tree) dynamics 
    void calcP(const SBStateRep&) const;                        // articulated body inertias
    void calcZ(const SBStateRep&, const SpatialVecList& spatialForces) const; // articulated body remainder forces
    void calcTreeAccel(const SBStateRep&) const;                // accels with forces from last calcZ

    void fixVel0(SBStateRep&, Vector& vel) const; // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY(const SBStateRep&) const;

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
    void calcInternalGradientFromSpatial(const SBStateRep&, const SpatialVecList& X, 
                                         Vector& JX);

    /// Pass in internal forces in T; they will be adjusted by this routine.
    void calcConstraintCorrectedInternalForces(const SBStateRep&, Vector& T); 

    const RigidBodyNode& getRigidBodyNode(int nodeNum) const {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }
    RigidBodyNode& updRigidBodyNode(int nodeNum) {
        const RigidBodyNodeIndex& ix = nodeNum2NodeMap[nodeNum];
        return *rbNodeLevels[ix.level][ix.offset];
    }

    bool isBuilt() const {return built;}

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
    SBState initialState;

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
