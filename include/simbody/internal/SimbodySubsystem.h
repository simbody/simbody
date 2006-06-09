#ifndef SimTK_SIMBODY_SUBSYSTEM_H_
#define SimTK_SIMBODY_SUBSYSTEM_H_

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/MatterSubsystem.h"

#include <cassert>
#include <vector>
#include <iostream>

class RigidBodyTree;

namespace SimTK {

class Transform;
class JointSpecification;
class InertiaMat;
class SimbodyTreeRep;
class MassProperties;



/**
 * The Simbody low-level multibody tree interface.
 * Equations represented:
 *
 *                     qdot = Q u
 *   M udot + ~A mult + R*C = T + R*F
 *                   A udot = b(t,q,u)
 *
 *                 v(t,q,u) = 0
 *                   g(t,q) = 0
 * 
 * where M(q) is the mass matrix, A(q) the constraint matrix, C(q,u)
 * the coriolis and gyroscopic forces, T is user-applied joint forces,
 * F is user-applied body forces and torques and gravity. 
 * R* is the operator that maps spatial forces to joint forces (partial
 * velocity matrix in Kane's terminology).
 *
 * We calculate the constraint multipliers like this:
 *           AM~A mult = A udot0-b, udot0=inv(M)f
 * using the pseudo inverse of AM~A to give a least squares solution for
 * mult: mult = pinv(AM~A)(A inv(M)f - b). Then the real udot is
 * udot = udot0 - udotC, with udotC = inv(M)(~A mult).
 *
 * NOTE: none of the above matrices are actually formed or factored!
 */
class SimTK_SIMBODY_API SimbodySubsystem : public MatterSubsystem {
public:
    /// Create a tree containing only the ground body (body 0).
    SimbodySubsystem();

    // These are the same as the compiler defaults but are handy to
    // have around explicitly for debugging.
    ~SimbodySubsystem() {
    }
    SimbodySubsystem(const SimbodySubsystem& ss) : MatterSubsystem(ss) {
    }
    SimbodySubsystem& operator=(const SimbodySubsystem& ss) {
        MatterSubsystem::operator=(ss);
        return *this;
    }

    /// Add a general rigid body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addRigidBody(const MassProperties&     massProps,
                     const Transform&          bodyJointFrameInB,
                     int                       parent,
                     const Transform&          parentJointFrameInP,
                     const JointSpecification& joint);

    /// Add a massless body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addMasslessBody(int                       parent,
                        const Transform&          parentJointFrameInP,
                        const JointSpecification& joint,
                        const Transform&          bodyJointFrameInB);

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
    int addWeldConstraint(int parent, const Transform& frameInP,
                          int child,  const Transform& frameInC);

    /// Topology and default values are frozen after this call. If you don't
    /// call it then it will be called automatically by realizeConstruction().
    void endConstruction();

    // Operators

    /// Requires realization through Stage::Configured.
    void calcInternalGradientFromSpatial(const State&,
        const Vector_<SpatialVec>& dEdR,
        Vector&                    dEdQ) const; // really Qbar

    /// Requires realization through MovingStage.
    Real calcKineticEnergy(const State&) const;

    /// Requires realization through DynamicsStage although
    /// velocities are irrelevant.
    void calcTreeEquivalentJointForces(const State&, 
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    jointForces) const;

    /// Requires realization through DynamicsStage.
    void calcTreeUDot(const State&,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot) const;

    /// Must be in Stage::Configured to calculate qdot = Q*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    /// Must be in Stage::Moving to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    // Constraint projections.

    /// Project position coordinates (q's) so that they satisfy their 
    /// constraints.
    void enforceConfigurationConstraints(State&) const;

    /// Project velocity coordinates (u's) so that they satisfy their
    /// constraints.
    void enforceMotionConstraints(State&) const;

    // These are available after realizeConstruction().

    /// The number of bodies includes all rigid bodies, particles, massless
    /// bodies and ground. Bodies and their inboard joints have the same 
    /// number, starting with ground at 0 with a regular labeling such
    /// that children have higher body numbers than their parents. Joint 0
    /// is meaningless (or I suppose you could think of it as the weld
    /// joint that attaches ground to the universe), but otherwise 
    /// joint n is the inboard joint of body n.
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
    void setUseEulerAngles(State&, bool) const;
    void setJointIsPrescribed(State&, int joint, bool) const;
    void setConstraintIsEnabled(State&, int constraint, bool) const;

    // Return modeling information from the State.
    bool getUseEulerAngles  (const State&) const;
    bool isJointPrescribed  (const State&, int joint)      const;
    bool isConstraintEnabled(const State&, int constraint) const;


    // Parameter setting & getting TODO
    void setGravity(State&, const Vec3& g) const; // any time after modeling
    const Vec3& getGravity(const State&) const;

    void setJointQ(State&, int body, int axis, const Real&) const;
    void setJointU(State&, int body, int axis, const Real&) const;

    // Configuration Stage. 


    /// Obtain the current orientation and position of the body frame of
    /// the indicated body. Must be in Configured stage. The configuration
    /// is provided as the Transform X_GB from the ground frame to the
    /// body frame.
    const Transform&  getBodyConfiguration(const State&, int body) const;

    /// Obtain the current spatial angular and linear velocity of the body frame of
    /// the indicated body. Must be in Moving stage. This is the velocity 
    /// V_GB of the body frame measured and expressed in the ground frame.
    const SpatialVec& getBodyVelocity     (const State&, int body) const;
    const SpatialVec& getBodyAcceleration (const State&, int body) const;

    const Real& getJointQ(const State&, int body, int axis) const;
    const Real& getJointU(const State&, int body, int axis) const;
    const Real& getJointQDot(const State&, int body, int axis) const;
    const Real& getJointUDot(const State&, int body, int axis) const;
    const Real& getJointQDotDot(const State&, int body, int axis) const;

    /// Get the location in space of a station (point) fixed on a body. This
    /// just makes use of the transform associated with the body's current
    /// configuration. Naturally the station is provided in the body frame, and
    /// the result is the vector from the ground origin to the station, expressed
    /// in the ground frame.
    const Vec3 getStationLocation(const State& s, int body, const Vec3& station_B) const {
        const Transform& X_GB = getBodyConfiguration(s, body);
        return X_GB.T() + X_GB.R() * station_B;
    }

    const Vector& getQ(const State&) const;
    const Vector& getU(const State&) const;
    const Vector& getAppliedMobilityForces(const State&) const;
    const Vector_<SpatialVec>& getAppliedBodyForces(const State&) const;

    void setQ(State&, const Vector& q) const;
    void setU(State&, const Vector& u) const;
    Vector& updQ(State&) const;
    Vector& updU(State&) const;

    const Vector& getQDot   (const State&) const;
    const Vector& getUDot   (const State&) const;
    const Vector& getQDotDot(const State&) const;

    SimTK_PIMPL_DOWNCAST(SimbodySubsystem, MatterSubsystem);
private:
    const RigidBodyTree& getRep() const;
    RigidBodyTree&       updRep();
};

SimTK_SIMBODY_API std::ostream& 
operator<<(std::ostream&, const SimbodySubsystem&);

};

#endif // SimTK_SIMBODY_SUBSYSTEM_H_
