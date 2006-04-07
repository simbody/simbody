#ifndef RIGID_BODY_TREE_H_
#define RIGID_BODY_TREE_H_

#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyTree.h"
using namespace SimTK;

class RigidBodyNode;
class ConstraintNode;
class RBDistanceConstraint;
class RBStation;

namespace SimTK {
class SBStateRep;
class SBModelingVars;
class SBParameterVars;
class SBTimeVars;
class SBConfigurationVars;
class SBMotionVars;
class SBDynamicsVars;
class SBReactionVars;
}

#include <cassert>
#include <vector>
#include <iostream>

typedef std::vector<RigidBodyNode*>   RBNodePtrList;
typedef Vector_<SpatialVec>           SpatialVecList;

class IVM;
class LengthConstraints;

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

    // Create a new node, add it to the tree, and assign it
    // a node number, which is a regular labeling starting with node 0 which is ground.
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


    // Constrain stations on each of two distinct bodies to remain a
    // particular distance apart at all times. Distance must be
    // significantly greater than 0 so that this can be implemented as a
    // single constraint force acting along the instantaneous line between
    // the stations. Parent and child distinction here is meaningless.
    int addConstantDistanceConstraint(const RigidBodyNode& parent, const Vec3& stationInP,
                                      const RigidBodyNode& child,  const Vec3& stationInC,
                                      const Real& distance);

    // Constrain stations on each of two distinct bodies to remain superimposed.
    // This restricts all translation but no rotation so adds three constraint
    // equations. Parent and child distinction here is meaningless.
    int addCoincidentStationsConstraint(const RigidBodyNode& parent, const Vec3& stationInP,
                                        const RigidBodyNode& child,  const Vec3& stationInC);

    // Constrain frames fixed to each of two distinct bodies to remain
    // superimposed. Parent and child here mean nothing! This adds six
    // constraint equations.
    int addWeldConstraint(const RigidBodyNode& parent, const TransformMat& frameInP,
                          const RigidBodyNode& child,  const TransformMat& frameInC);


    /// This is available any time.
    Stage getStage(const SBStateRep&) const;

    /// This will realize the state up to the indicated stage, possibly advancing
    /// multiple stages internally by calling each "realizeWhatever()" routine in order.
    void realize(const SBStateRep& s, Stage stage) const;

    /// Call this after all bodies & constraints have been added.
    void realizeConstruction (); // will set built==true
    void realizeModeling     (const SBStateRep&) const;
    void realizeParameters   (const SBStateRep&) const;
    void realizeTime         (const SBStateRep&) const;
    void realizeConfiguration(const SBStateRep&) const;
    void realizeMotion       (const SBStateRep&) const;
    void realizeDynamics     (const SBStateRep&) const;
    void realizeReaction     (const SBStateRep&) const;

    Real calcKineticEnergy(const SBStateRep&) const;


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
    void calcInternalGradientFromSpatial(const SBStateRep&, 
        const SpatialVecList& X, 
        Vector&               JX);

    // Given a set of body forces, return the equivalent set of joint torques 
    // IGNORING CONSTRAINTS.
    // Must be in DynamicsStage so that articulated body inertias are available,
    // however, velocities are ignored. This operator has NO effect on the state
    // cache. It makes a single O(N) pass.
    void calcTreeEquivalentJointForces(const SBStateRep&, 
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    jointForces);

    void calcTreeAccelerations(const SBStateRep& s,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    netHingeForces,
        Vector_<SpatialVec>&       A_GB,
        Vector&                    udot) const; 

    // Must be in Stage::Configured to calculate qdot = Q*u.
    void calcQDot(const SBStateRep& s,
        const Vector& u,
        Vector&       qdot) const;

    // Must be in MovingStage to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const SBStateRep& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    const SBState& getInitialState() const {
        assert(built); return initialState;
    }
    void setDefaultModelingValues     (const SBStateRep&, SBModelingVars&)      const;
    void setDefaultParameterValues    (const SBStateRep&, SBParameterVars&)     const;
    void setDefaultTimeValues         (const SBStateRep&, SBTimeVars&)          const;
    void setDefaultConfigurationValues(const SBStateRep&, SBConfigurationVars&) const;
    void setDefaultMotionValues       (const SBStateRep&, SBMotionVars&)        const;
    void setDefaultDynamicsValues     (const SBStateRep&, SBDynamicsVars&)      const;
    void setDefaultReactionValues     (const SBStateRep&, SBReactionVars&)      const;

    // These counts can be obtained even during construction, where they
    // just return the current counts.
    // includes ground
    int getNBodies()      const {return nodeNum2NodeMap.size();}
    int getNConstraints() const {return constraintNodes.size();}

    // A single constraint may generate multiple of these.
    int getNDistanceConstraints() const {return distanceConstraints.size();}

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

    void enforceConfigurationConstraints(SBStateRep&) const;
    void enforceMotionConstraints(SBStateRep&) const;

    /// Unconstrained (tree) dynamics 
    void calcArticulatedBodyInertias(const SBStateRep&) const;                        // articulated body inertias
    void calcZ(const SBStateRep&, const SpatialVecList& spatialForces) const; // articulated body remainder forces
    void calcTreeAccel(const SBStateRep&) const;                // accels with forces from last calcZ

    void fixVel0(SBStateRep&, Vector& vel) const; // TODO -- yuck

    /// Part of constrained dynamics (TODO -- more to move here)
    void calcY(const SBStateRep&) const;


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

    // Add a distance constraint equation and assign it a particular multiplier
    // slot to use. Return the assigned distance constraint index for caller's use.
    int addOneDistanceConstraintEquation(
        const RBStation& s1, const RBStation& s2, const Real& d,
        int multIndex);

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

    // This holds pointers to nodes and serves to map (level,offset) to nodeNum.
    Array<RBNodePtrList>      rbNodeLevels;
    // Map nodeNum to (level,offset).
    std::vector<RigidBodyNodeIndex> nodeNum2NodeMap;

    // This holds pointers to the abstract constraint nodes which correspond
    // the the user's idea of constraints in a manner analogous to the
    // linked bodies represented by RigidBodyNodes. Each of these may generate
    // several constraint equations.
    Array<ConstraintNode*>       constraintNodes;
    Array<RBDistanceConstraint*> distanceConstraints;
    
    LengthConstraints* lConstraints;

    void addGroundNode();
    int addConstraintNode(ConstraintNode*&);

    // Given a forces in the state, calculate accelerations ignoring
    // constraints, and leave the results in the state. 
    // Must have already called realizeDynamics().
    // We also allow some extra forces to be supplied, with the intent
    // that these will be used to deal with internal forces generated
    // by constraints. 
    void calcTreeForwardDynamics (const SBStateRep& s,
        const Vector*              extraJointForces,
        const Vector_<SpatialVec>* extraBodyForces) const;

    // Given a set of forces in the state, calculate acclerations resulting from
    // those forces and enforcement of acceleration constraints, and update the state.
    void calcLoopForwardDynamics(const SBStateRep&) const;

    friend std::ostream& operator<<(std::ostream&, const RigidBodyTree&);
    friend class SimbodyTree;
};

std::ostream& operator<<(std::ostream&, const RigidBodyTree&);

#endif // RIGID_BODY_TREE_H_
