#ifndef SimTK_SIMBODY_TREE_H_
#define SimTK_SIMBODY_TREE_H_

#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyState.h"

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
 * This class holds a full set of applied forces, including body
 * forces (by which we mean forces and torques), and joint forces..
 *
 * Body forces are expressed in the ground frame. Joints constitute
 * their own frames, and joint forces are applied as per-joint-DOF
 * scalars whose meanings depend on the quirks of the individual
 * joints.
 *
 * Note that this is a low-level container class; it does not know
 * the joint types or current configuration of the bodies. It just
 * knows how many items of each type there are. It must be given to a
 * SimbodyTree, along with the current state, in order to be
 * interpreted correctly.
 */
class SimbodyForceSet {
public:
    SimbodyForceSet() { }
    SimbodyForceSet(int nBody, int nDOFs) : bodyForces(nBody), jointForces(nDOFs) { }
    // default copy, assignment, destructor

    void resize(int nBody, int nDOFs) 
       {bodyForces.resize(nBody); jointForces.resize(nDOFs);}

    /// Set all forces to zero.
    void clear() {clearBodyForces(); clearJointForces();}

    /// Clear only the body forces, leaving joint forces alone.
    void clearBodyForces()  {bodyForces.setToZero();}

    /// Clear only the joint forces, leaving the body forces alone.
    void clearJointForces() {jointForces.setToZero();}

    const Vector_<SpatialVec>& getBodyForces() const {return bodyForces;}
    Vector_<SpatialVec>&       updBodyForces()       {return bodyForces;}

    const Vector& getJointForces() const {return jointForces;}
    Vector&       updJointForces()       {return jointForces;}

    // Named assignment operators; prefer the actual operators in C++
    SimbodyForceSet& assign  (const SimbodyForceSet& f) {return (*this = f);}
    SimbodyForceSet& add     (const SimbodyForceSet& f) 
       {bodyForces+=f.bodyForces; jointForces+=f.jointForces; return *this;}
    SimbodyForceSet& subtract(const SimbodyForceSet& f) 
       {bodyForces-=f.bodyForces; jointForces-=f.jointForces; return *this;}
    SimbodyForceSet& scale   (const Real& s) {bodyForces*=s; jointForces*=s; return *this;}

    SimbodyForceSet& operator+=(const SimbodyForceSet& f) {return add(f);}
    SimbodyForceSet& operator-=(const SimbodyForceSet& f) {return subtract(f);}
    SimbodyForceSet& operator*=(const Real& s)            {return scale(s);}
private:
    Vector_<SpatialVec> bodyForces;
    Vector_<Vec3>       particleForces;
    Vector_<Real>       jointForces;
};

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
 * using the pseudo inverse to give a least squares solution for
 * mult: mult = pinv(AM~A)(A inv(M)f - b). Then the real udot is
 * udot = udot0 - udotC, with udotC = inv(M)(~A mult).
 */
class SimTK_SIMBODY_API SimbodyTree {
public:
    /// Create a tree containing only the ground body (body 0).
    SimbodyTree();
    ~SimbodyTree();
    SimbodyTree(const SimbodyTree&);
    SimbodyTree& operator=(const SimbodyTree&);

    /// Add a general rigid body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addRigidBody(int                       parent,
                     const Transform&          parentJointFrameInP,
                     const JointSpecification& joint,
                     const Transform&          bodyJointFrameInB,
                     const MassProperties&);

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

    /// Topology and default values are frozen after this call, and all
    /// Modeling variables have been allocated slots in the State and
    /// given appropriate default values. The State's stage cannot be
    /// any higher than Allocated.
    void realizeConstruction(State&);

    /// All Modeling choices are frozen after this call, and all remaining
    /// State variables have been allocated and given appropriate initial
    /// values. Modeling variables are NOT changed by this call. The
    /// State's stage cannot be any higher than Built.
    void realizeModeling    (State&) const;

    void realizeParameters   (const State&) const;
    void realizeTime         (const State&) const;
    void realizeConfiguration(const State&) const;
    void realizeMotion       (const State&) const;
    void realizeDynamics     (const State&) const;

    /// To generate a reaction (i.e., accelerations) from this State we
    /// will need the help of a ForceSubsystem to extract its forces
    /// from the State.
    /// TODO: an alternative design would be not to allow this here
    /// but to have the common parent System generate the reactions.
    void realizeReaction     (const State&/*, const ForceSubsystem&*/) const;

    void realize(const State&, Stage) const;

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

    // Must be in Stage::Configured to calculate qdot = Q*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    // Must be in Stage::Moving to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    // Constraint projections.

    void enforceConfigurationConstraints(State&) const;
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

    // Return modeling information from the SBState.
    bool getUseEulerAngles  (const State&) const;
    bool isJointPrescribed  (const State&, int joint)      const;
    bool isConstraintEnabled(const State&, int constraint) const;


    // Parameter setting & getting
    void setGravity(State&, const Vec3& g) const; // any time after modeling
    const Vec3& getGravity(const State&) const;

    void setJointQ(State&, int body, int axis, const Real&) const;
    void setJointU(State&, int body, int axis, const Real&) const;
    void setPrescribedUdot(State&, int body, int axis, const Real&) const;
    void clearAppliedForces(State&) const;

    // Configuration Stage. 

    // OBSOLETE
    void applyGravity    (State&, const Vec3& g) const;


    void applyPointForce (State&, int body, const Vec3& stationInB, 
                          const Vec3& forceInG) const;
    void applyBodyTorque (State&, int body, 
                          const Vec3& torqueInG) const;
    void applyJointForce(State&, int body, int axis, const Real&) const;


    const Transform&  getBodyConfiguration(const State&, int body) const;
    const SpatialVec& getBodyVelocity     (const State&, int body) const;
    const SpatialVec& getBodyAcceleration (const State&, int body) const;

    const Vec3 getStationLocation(const State& s, int body, const Vec3& station_B) const {
        const Transform& X_GB = getBodyConfiguration(s, body);
        return X_GB.T() + X_GB.R() * station_B;
    }

    const Vector& getQ(const State&) const;
    const Vector& getU(const State&) const;
    const Vector& getAppliedJointForces(const State&) const;
    const Vector_<SpatialVec>& getAppliedBodyForces(const State&) const;

    void setQ(State&, const Vector& q) const;
    void setU(State&, const Vector& u) const;
    Vector& updQ(State&) const;
    Vector& updU(State&) const;

    const Vector& getQDot   (const State&) const;
    const Vector& getUDot   (const State&) const;
    const Vector& getQDotDot(const State&) const;

private:
    RigidBodyTree* rep;
};

SimTK_SIMBODY_API std::ostream& 
operator<<(std::ostream&, const SimbodyTree&);

};

#endif // SimTK_SIMBODY_TREE_H_
