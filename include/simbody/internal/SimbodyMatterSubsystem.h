#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_H_

/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/MatterSubsystem.h"

#include <cassert>
#include <vector>
#include <iostream>

class SimbodyMatterSubsystemRep;

namespace SimTK {

class Transform;
class Mobilizer;
class Inertia;
class SimbodyTreeRep;
class MassProperties;



/**
 * The Simbody low-level multibody tree interface.
 * Equations represented:
 *
 *                               qdot = Q u
 *                               n(q) = 0
 *
 *                   M udot + ~G mult = f
 *                  G udot + b(t,q,u) = 0
 *
 *   [A]    [ba]
 * G=[V]  b=[bv]  f=T+R*(F-C)
 *   [P]    [bp]
 *
 * a(t,q,u,udot) = A udot + ba(t,q,u) = 0
 *          vdot = V udot + bv(t,q,u) = 0
 *       pdotdot = P udot + bp(t,q,u) = 0
 *           
 *                           v(t,q,u) = 0
 *                pdot = P u + c(t,q) = 0
 *
 *                             p(t,q) = 0
 * 
 * where M(q) is the mass matrix, G(q) the acceleration constraint matrix, C(q,u)
 * the coriolis and gyroscopic forces, T is user-applied joint forces,
 * F is user-applied body forces and torques and gravity. 
 * R* is the operator that maps spatial forces to joint forces. p() are the
 * holonomic (position) constraints, v() the non-holonomic (velocity) constraints,
 * and a() the reaction (acceleration) constraints, which must be linear with A
 * the coefficient matrix for a(). pdot, pdotdot are obtained
 * by differentiation of p(), vdot by differentiation of v().
 * P=partial(pdot)/partial(u), V=partial(v)/partial(u).
 * n() is the set of quaternion normalization constraints.
 *
 * We calculate the constraint multipliers like this:
 *           G M^-1 ~G mult = G udot0 - b, udot0=M^-1 f
 * using the pseudo inverse of G M^-1 ~G to give a least squares solution for
 * mult: mult = pinv(G M^-1 ~G)(G M^-1 f - b). Then the real udot is
 * udot = udot0 - udotC, with udotC = M^-1 ~G mult.
 *
 * NOTE: only the constraint matrices have to be formed and factored:
 *     G M^-1 ~G    to calculate multipliers (square, symmetric: LDL' if
 *                  well conditioned, else pseudoinverse)
 *
 *     P            for projection onto position manifold (pseudoinverse)
 *
 *    [V]           for projection onto velocity manifold (pseudoinverse)
 *    [P]
 *
 * In many cases these matrices consisted of decoupled blocks which can
 * be solved independently; we try to take advantage of that whenever possible
 * to solve a set of smaller systems rather than one large one. Also, in the
 * majority of biosimulation applications we are likely to have only holonomic
 * (position) constraints, so there is no V or A so P is the whole story.
 */
class SimTK_SIMBODY_EXPORT SimbodyMatterSubsystem : public MatterSubsystem {
public:
    /// Create a tree containing only the ground body (body 0).
    SimbodyMatterSubsystem();

    // These are the same as the compiler defaults but are handy to
    // have around explicitly for debugging.
    ~SimbodyMatterSubsystem() {
    }
    SimbodyMatterSubsystem(const SimbodyMatterSubsystem& ss) : MatterSubsystem(ss) {
    }
    SimbodyMatterSubsystem& operator=(const SimbodyMatterSubsystem& ss) {
        MatterSubsystem::operator=(ss);
        return *this;
    }

    /// Add a general rigid body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    BodyId addRigidBody(const MassProperties&  massProps,
                     const Transform&          bodyJointFrameInB,
                     BodyId                    parent,
                     const Transform&          parentJointFrameInP,
                     const Mobilizer&          mobilizer);

    /// Add a massless body to the growing tree by connecting it
    /// to one of the bodies already in the tree.
    BodyId addMasslessBody(const Transform&       bodyJointFrameInB,
                        BodyId                    parent,
                        const Transform&          parentJointFrameInP,
                        const Mobilizer&          mobilizer);

    /// Special case for convenience: attach a general rigid body to
    /// a body (ground by default) using a free joint and only the
    /// body frames.
    /// NOTE: it is *NOT* allowed to add bodies outboard to this one.
    BodyId addFreeRigidBody(const MassProperties&, BodyId parent=GroundId);

    /// Special case: add a free particle (point mass) to the tree
    /// by connecting it to a body (ground by default) using a
    /// Cartesian joint (3d translation) with fixed frame the parent's
    /// body frame and the point location as the moving "frame".
    /// NOTE: it is *NOT* allowed to add bodies outboard to this one.
    BodyId addFreeParticle (const Real& mass,      BodyId parent=GroundId);

    /// Constrain stations on each of two distinct bodies to remain
    /// a particular distance apart at all times. Distance must be
    /// significantly greater than 0 so that this can be implemented
    /// as a single constraint force acting along the instantaneous
    /// line between the stations.
    ConstraintId addConstantDistanceConstraint(BodyId parent, const Vec3& stationInP,
                                               BodyId child,  const Vec3& stationInC,
                                               const Real& distance);

    /// Constrain stations on each of two distinct bodies to remain
    /// superimposed. This restricts all translation but no rotation
    /// so adds three constraint equations.
    ConstraintId addCoincidentStationsConstraint(BodyId parent, const Vec3& stationInP,
                                                 BodyId child,  const Vec3& stationInC);

    /// Constrain frames fixed to each of two distinct bodies to
    /// remain superimposed. Parent and child here mean nothing!
    /// This adds six constraint equations.
    ConstraintId addWeldConstraint(BodyId parent, const Transform& frameInP,
                                   BodyId child,  const Transform& frameInC);

    /// Topology and default values are frozen after this call. If you don't
    /// call it then it will be called automatically by realizeTopology().
    void endConstruction();

    // Operators

    /// Requires realization through Stage::Position.
    void calcInternalGradientFromSpatial(const State&,
        const Vector_<SpatialVec>& dEdR,
        Vector&                    dEdQ) const; // really Qbar

    /// Requires realization through MovingStage.
    Real calcKineticEnergy(const State&) const;

    /// Requires realization through DynamicsStage although
    /// velocities are irrelevant.
    void calcTreeEquivalentMobilityForces(const State&, 
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    mobilityForces) const;

    /// Requires realization through DynamicsStage.
    void calcTreeUDot(const State&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot,
        Vector_<SpatialVec>&       A_GB) const;

    /// Requires realization through DynamicsStage.
    void calcMInverseF(const State&,
        const Vector&        f,
        Vector&              udot,
        Vector_<SpatialVec>& A_GB) const;

    /// Must be in Stage::Position to calculate qdot = Q*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    /// Must be in Stage::Velocity to calculate qdotdot = Qdot*u + Q*udot.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    // Constraint projections.

    /// Project position coordinates (q's) so that they satisfy their 
    /// constraints to at least tol.
    void enforcePositionConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    /// Project velocity coordinates (u's) so that they satisfy their
    /// constraints to at least tol.
    void enforceVelocityConstraints(State&, const Real& requiredTol, const Real& desiredTol) const;

    // These are available after realizeTopology().

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
    int getQIndex(BodyId) const;
    int getQAlloc(BodyId) const; // must wait for modeling for actual NQ
    int getUIndex(BodyId) const;
    int getDOF   (BodyId) const; // always same as # u's


    // Per-constraint info;
    int getMultIndex(int constraint) const;
    int getMaxNMult (int constraint) const;  // wait for modeling to get actual NMult


    /// For all ball and free joints, decide what method we should use
    /// to model their orientations. Choices are: quaternions (best
    /// for dynamics), or rotation angles (3-2-1 Euler sequence, good for
    /// optimization). TODO: allow settable zero rotation for Euler sequence,
    /// with convenient way to say "this is zero".
    void setUseEulerAngles(State&, bool) const;
    void setMobilizerIsPrescribed(State&, BodyId, bool) const;
    void setConstraintIsEnabled(State&, int constraint, bool) const;

    // Return modeling information from the State.
    bool getUseEulerAngles  (const State&) const;
    bool isMobilizerPrescribed  (const State&, BodyId) const;
    int  getNQuaternionsInUse(const State&) const;
    bool isUsingQuaternion(const State&, BodyId) const;
    int  getQuaternionIndex(const State&, BodyId) const;

    bool isConstraintEnabled(const State&, int constraint) const;

    // Position Stage. 

    // Dynamics stage responses.

    // Cross joint
    const SpatialVec& getCoriolisAcceleration(const State&, BodyId) const;

    // Including parent
    const SpatialVec& getTotalCoriolisAcceleration(const State&, BodyId) const;

    const SpatialVec& getGyroscopicForce(const State&, BodyId) const;
    const SpatialVec& getCentrifugalForces(const State&, BodyId) const;
    const SpatialMat& getArticulatedBodyInertia(const State& s, BodyId) const;

    const Real& getMobilizerQDot(const State&, BodyId, int axis) const;
    const Real& getMobilizerUDot(const State&, BodyId, int axis) const;
    const Real& getMobilizerQDotDot(const State&, BodyId, int axis) const;

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

    SimTK_PIMPL_DOWNCAST(SimbodyMatterSubsystem, Subsystem);
private:
    const SimbodyMatterSubsystemRep& getRep() const;
    SimbodyMatterSubsystemRep&       updRep();
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem&);

};

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
