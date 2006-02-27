#ifndef SIMTK_SIMBODY_TREE_H_
#define SIMTK_SIMBODY_TREE_H_

#include "simbody/Simbody.h"

#include <cassert>
#include <vector>
#include <iostream>

class RigidBodyTree;

namespace simtk {

class TransformMat;
class JointSpecification;
class InertiaMat;
class SimbodyTreeRep;
class MassProperties;

/** 
 * Coordinate allocation:
 * Body 0 is ground and has no coordinates. Welded bodies have a higher body number
 * but no coordinates. We still set qoff and uoff as we would if these had dofs,
 * but nu==nu==0 for these so the next body has the same qoff & uoff.
 *
 *
 *   Body 0     1       2    3      4      ...
 *         ---------- ----- --- ---------- ---
 *      q |          |     |   |          |
 *         ---------- ----- --- ---------- ---
 *                    ^
 *                    qoff[2], nq[2], nqmax[2]
 *
 *   Body 0     1     2    3    4      ...
 *         -------- ----- --- -------- ---
 *      u |        |     |   |        |
 *         -------- ----- --- -------- ---
 *                  ^
 *                  uoff[2], nu[2] (==ndof[2])
 *
 * qdot, qdotdot are allocated exactly as for q. Derivatives of unused coordinates
 * are set to 0.
 * udot, prescribedUdot are allocated exactly as for u. All udots are set.
 * Entries of prescribedUdot which do not correspond to prescribed joints will
 * be ignored and need not be set.
 */

enum SBStage {
    UninitializedStage  = 0, // these are ordered
    BuiltStage          = 1,
    ModeledStage        = 2,
    ParametrizedStage   = 3,
    TimedStage          = 4,
    ConfiguredStage     = 5,
    MovingStage         = 6,
    DynamicsStage       = 7, // dynamic properties & operators available
    ReactingStage       = 8  // accelerations & reaction forces available
};

// This is the handle class for the hidden SBState implementation.
class SBStateRep;
class SBState {
public:
    SBState() : rep(0) { }
    ~SBState();
    SBState(const SBState&);
    SBState& operator=(const SBState&);

    const SBStateRep& getRep() const {assert(rep); return *rep;}
    SBStateRep&       updRep()       {assert(rep); return *rep;}
    SBStateRep* rep;
};


/**
 * The Simbody low-level multibody tree interface.
 * Equations represented:
 *
 *               qdot = Q u
 *   M udot + ~A mult = f
 *             A udot = b(t,q,u)
 *
 *           v(t,q,u) = 0
 *             g(t,q) = 0
 * 
 * where            f = T - R*(F + C(q,u))
 * where T is user-applied joint forces, F is user-applied body
 * forces and torques and gravity, and C includes coriolis and
 * gyroscopic forces. R* is the operator that maps spatial forces
 * to joint forces.
 *
 * We calculate the constraint multipliers like this:
 *           AM~A mult = A udot0-b, udot0=inv(M)f
 * using the pseudo inverse to give a least squares solution for
 * mult: mult = pinv(AM~A)(A inv(M)f - b). Then the real udot is
 * udot = udot0 - udotC, with udotC = inv(M)(~A mult).
 */
class SimbodyTree {
public:
    /// Create a tree containing only the ground body (body 0).
    SimbodyTree();

    /// Add a general rigid body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addRigidBody(int                       parent,
                     const TransformMat&       parentJointFrameInP,
                     const JointSpecification& joint,
                     const TransformMat&       bodyJointFrameInB,
                     const MassProperties&);

    /// Add a massless body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addMasslessBody(int                       parent,
                        const TransformMat&       parentJointFrameInP,
                        const JointSpecification& joint,
                        const TransformMat&       bodyJointFrameInB);

    /// Special case for convenience: attach a general rigid body to
    /// a body (ground by default) using a free joint and only the
    /// body frames.
    int addFreeRigidBody(const MassProperties&, int parent=0);

    /// Special case: add a free particle (point mass) to the tree
    /// by connecting it to a body (ground by default) using a
    /// Cartesian joint (3d translation) with fixed frame the parent's
    /// body frame and the point location as the moving "frame".
    int addFreeParticle (const Real& mass,      int parent=0);

    /// Constrain stations on each of two distinct bodies to remain
    /// a particular distance apart at all times. Distance must be
    /// significantly greater than 0 so that this can be implemented
    /// as a single constraint force acting along the instantaneous
    /// line between the stations.
    int addConstantDistanceConstraint(int parent, const Vec3& stationInP,
                                      int child,  const Vec3& stationInC,
                                      const Real& distance);

    /// Constrain stations on each of two distinct bodies to remain
    /// superimposed. This restricts all translation but no rotation
    /// so adds three constraint equations.
    int addCoincidentStationsConstraint(int parent, const Vec3& stationInP,
                                        int child,  const Vec3& stationInC);

    /// Constrain frames fixed to each of two distinct bodies to
    /// remain superimposed. Parent and child here mean nothing!
    /// This adds six constraint equations.
    int addWeldConstraint(int parent, const TransformMat& frameInP,
                          int child,  const TransformMat& frameInC);

    /// Topology and default values are frozen after this call.
    // TODO: "formulation" instead of "modeling" stage?
    void realizeConstruction();
    const SBState& getInitialState() const;

    void realizeModeling     (const SBState&) const;
    void realizeParameters   (const SBState&) const;
    void realizeTime         (const SBState&) const;
    void realizeConfiguration(const SBState&) const;
    void realizeMotion       (const SBState&) const;
    void realizeDynamics    (const SBState&) const;
    void realizeReaction     (const SBState&) const;

    void realize(const SBState&, SBStage) const;

    // Operators

    /// Requires realization through ConfiguredStage.
    void calcInternalGradientFromSpatial(const SBState&,
        const Vector_<SpatialVec>& dEdR,
        Vector&                    dEdQ) const; // really Qbar

    /// Requires realization through MovingStage.
    Real calcKineticEnergy(const SBState&) const;

    /// Requires realization through DynamicsStage.
    void calcTreeUDot(const SBState&,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot) const;

    // Constraint projections.

    void enforceConfigurationConstraints(SBState&) const;
    void enforceMotionConstraints(SBState&) const;

    // These are available after realizeConstruction().

    /// The number of bodies includes all rigid bodies, particles, massless
    /// bodies and ground. Bodies and their inboard joints have the same 
    /// number, starting with ground at 0 with a regular labeling such
    /// that children have higher body numbers than their parents. Joint 0
    /// is meaningless, but otherwise joint n is the inboard joint of body n.
    int getNBodies() const;

    /// This is the total number of defined constraints, each of which may
    /// generate more than one constraint equation.
    int getNConstraints() const;

    /// The sum of all the joint degrees of freedom. This is also the length
    /// of state variable vector u.
    int getTotalDOF() const; 

    /// The sum of all the q vector allocations for each joint. These may not
    /// all be in use.
    int getTotalQAlloc() const;

    /// This is the sum of all the allocations for constraint multipliers.
    int getTotalMultAlloc() const;

    // Per-body info.
    int getQIndex(int body) const;
    int getQAlloc(int body) const; // must wait for modeling for actual NQ
    int getUIndex(int body) const;
    int getDOF   (int body) const; // always same as # u's

    // Per-constraint info;
    int getMultIndex(int constraint) const;
    int getMaxNMult (int constraint) const;  // wait for modeling to get actual NMult


    /// For all ball and free joints, decide what method we should use
    /// to model their orientations. Choices are: quaternions (best
    /// for dynamics), or rotation angles (3-2-1 Euler sequence, good for
    /// optimization). TODO: allow settable zero rotation for Euler sequence,
    /// with convenient way to say "this is zero".
    void setUseEulerAngles(SBState&, bool) const;
    void setJointIsPrescribed(SBState&, int joint, bool) const;
    void setConstraintIsEnabled(SBState&, int constraint, bool) const;

    // Return modeling information from the SBState.
    bool getUseEulerAngles  (const SBState&) const;
    bool isJointPrescribed  (const SBState&, int joint)      const;
    bool isConstraintEnabled(const SBState&, int constraint) const;


    // Parameter setting & getting would go here 

    void setJointQ(SBState&, int body, int axis, const Real&) const;
    void setJointU(SBState&, int body, int axis, const Real&) const;
    void setPrescribedUdot(SBState&, int body, int axis, const Real&) const;
    void clearAppliedForces(SBState&) const;

    /// Configuration Stage. 

    void applyGravity    (SBState&, const Vec3& g) const;
    void applyPointForce (SBState&, int body, const Vec3& stationInB, 
                          const Vec3& forceInG) const;
    void applyBodyTorque (SBState&, int body, 
                          const Vec3& torqueInG) const;
    void applyJointForce(SBState&, int body, int axis, const Real&) const;


    const TransformMat& getBodyConfiguration(const SBState&, int body) const;
    const SpatialVec&   getBodyVelocity     (const SBState&, int body) const;
    const SpatialVec&   getBodyAcceleration (const SBState&, int body) const;

    const Vector& getQ(const SBState&) const;
    const Vector& getU(const SBState&) const;
    const Vector& getAppliedJointForces(const SBState&) const;
    const Vector_<SpatialVec>& getAppliedBodyForces(const SBState&) const;

    void setQ(SBState&, const Vector& q) const;
    void setU(SBState&, const Vector& u) const;
    VectorView& updQ(SBState&) const;
    VectorView& updU(SBState&) const;

    const Vector& getQDot   (const SBState&) const;
    const Vector& getUDot   (const SBState&) const;
    const Vector& getQDotDot(const SBState&) const;

private:
    RigidBodyTree* rep;
    //SimbodyTreeRep* rep;
};

std::ostream& operator<<(std::ostream&, const SimbodyTree&);

};

#endif // SIMTK_SIMBODY_TREE_H_
