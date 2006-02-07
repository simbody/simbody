#ifndef SIMTK_SIMBODY_TREE_H_
#define SIMTK_SIMBODY_TREE_H_

#include "simbody/Simbody.h"

#include <cassert>
#include <vector>
#include <iostream>

namespace simtk {

class TransformMat;
class JointSpecification;
class InertiaMat;
class ForceSystem;
class SimbodyTreeRep;
class MassProperties;

struct SimbodyTreeResults {
    Stage realizationLevel;     // must be kept up to date by State changes

    // TODO: constraint runtimes

    // TODO: Modeling
    //   counts of various things resulting from modeling choices,
    //   constraint enabling, prescribed motion

    // TODO: Parameters
    //   body mass props; particle masses
    //   X_BJ, X_PJi transforms
    //   distance constraint distances & station positions

    // Configuration
    std::vector<TransformMat> bodyConfigInParent; // nb (joint config)
    std::vector<TransformMat> bodyConfigInGround; // nb
    Vector_<SpatialMat>       bodySpatialInertia; // nb

    Vector positionConstraintErrors;              // npc
    Matrix_<Vec3> storageForHt;                   // 2 x ndof

    // Motion
    Vector_<SpatialVec> bodyVelocityInParent;     // nb (joint velocity)
    Vector_<SpatialVec> bodyVelocityInGround;     // nb

    Vector velocityConstraintErrors;              // nvc
    Vector qdot;                                  // nq

    // Dynamics
    Vector_<SpatialMat> articulatedBodyInertia;   // nb (P)
    Vector_<SpatialVec> bodyAccelerationInGround; // nb
    Vector_<SpatialVec> coriolisForces;           // nb (& gyroscopic, Pa+b)

    Vector accelerationConstraintErrors;          // nac
    Vector udotAndLambda;                         // nu+nac
    Vector netHingeForces;                        // nu (T-(~Am+R(F+C))
    Vector qdotdot;                               // nq

    // dynamic temporaries
    Vector_<Real>       storageForDI;   // sum(nu[j]^2)
    Matrix_<Vec3>       storageForG;    // 2 X ndof
    Vector              nu;
    Vector              epsilon;

};

struct SimbodyTreeState {
    // Modeling
    bool              useEulerAngles;
    std::vector<bool> prescribed;           // nb
    std::vector<bool> enabled;              // nac

    // Parametrization
    // TODO: body masses, etc.

    // Configuration
    Vector q;                               // nq

    // Motion
    Vector u;                               // nu

    // Dynamics
    Vector_<SpatialVec> appliedBodyForces;  // nb
    Vector              appliedJointForces; // nu
    Vector              prescribedUdot;     // nu

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
    int addRigidBody(int parent,
                     const TransformMat&       parentJointFrameInP,
                     const JointSpecification& joint,
                     const TransformMat&       bodyJointFrameInB,
                     const MassProperties&);

    /// Add a massless body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    int addMasslessBody(int parent,
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
    void         finishConstruction();
    const State& getDefaultState() const;
    void         realize(const State&, Stage) const;

    // These require realize(Modeling) afterwards.

    /// For all ball and free joints, decide what method we should use
    /// to model their orientations. Choices are: quaternions (best
    /// for dynamics), or rotation angles (1-2-3 Euler sequence, good for
    /// optimization). TODO: allow settable zero rotation for Euler sequence,
    /// with convenient way to say "this is zero".
    void setUseRotationAngles(State&, bool) const;
    void setJointIsPrescribed(State&, int joint, bool) const;
    void setConstraintIsEnabled(State&, int constraint, bool) const;


    void setJointQ(State&, int joint, const Real*) const;
    void setJointU(State&, int joint, const Real*) const;
    void setPrescribedUdot(State&, int joint, const Real*) const;

    // Return modeling information from the State.
    bool isJointPrescribed  (const State&, int joint) const;
    bool isConstraintEnabled(const State&, int constraint) const;


    void clearAppliedForces(State&) const;

    void applyGravity    (State&) const;
    void applyPointForce (State&, int body, const Vec3& stationInB, 
                          const Vec3& forceInG) const;
    void applyBodyTorque (State&, int body, 
                          const Vec3& torqueInG) const;
    void applyJointForce(State&, int joint, const Real*) const;


    const TransformMat& getBodyConfiguration(const State&, int body) const;
    const SpatialVec&   getBodyVelocity     (const State&, int body) const;
    const SpatialVec&   getBodyAcceleration (const State&, int body) const;

    Real&        updTime(State&) const;
    Vector&      updQ   (State&) const;
    Vector&      updU   (State&) const;
    ForceSystem& updForceSystem(State&) const;

    const Vector& getQDot   (const State&) const;
    const Vector& getUDot   (const State&) const;
    const Vector& getQDotDot(const State&) const;

private:
    SimbodyTreeRep* rep;
};

std::ostream& operator<<(std::ostream&, const SimbodyTree&);

};

#endif // SIMTK_SIMBODY_TREE_H_
