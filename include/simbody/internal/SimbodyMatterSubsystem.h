#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"

#include <cassert>
#include <vector>
#include <iostream>

class SimbodyMatterSubsystemRep;

namespace SimTK {

class MultibodySystem;
class Constraint;

/**
 * The Simbody low-level multibody tree interface.
 * Equations represented:                                  @verbatim
 *
 *                  qdot = Q u
 *                  zdot = zdot(t,q,u,z)
 *
 *      M udot + ~G mult = f(t,q,u,z)
 *      G udot           = b(t,q,u) 
 *
 *              where
 *
 *       [P]    [bp]
 *     G=[V]  b=[bv]  f=T+J*(F-C)
 *       [A]    [ba]
 *
 * a(t,q,u,udot) = A udot - ba(t,q,u) = 0
 *          vdot = V udot - bv(t,q,u) = 0
 *       pdotdot = P udot - bp(t,q,u) = 0
 *           
 *                           v(t,q,u) = 0
 *                pdot = P u - c(t,q) = 0
 *
 *                             p(t,q) = 0
 *                               n(q) = 0                  @endverbatim
 *
 * 
 * where M(q) is the mass matrix, G(q) the acceleration constraint matrix, C(q,u)
 * the coriolis and gyroscopic forces, T is user-applied joint mobility forces,
 * F is user-applied body forces and torques and gravity. 
 * J* is the operator that maps spatial forces to joint mobility forces. p() are the
 * holonomic (position) constraints, v() the non-holonomic (velocity) constraints,
 * and a() the acceleration-only constraints, which must be linear, with A
 * the coefficient matrix for a(). pdot, pdotdot are obtained
 * by differentiation of p(), vdot by differentiation of v().
 * P=partial(pdot)/partial(u) (yes, that's u, not q), V=partial(v)/partial(u).
 * n(q) is the set of quaternion normalization constraints, which exist only at the
 * position level and are uncoupled from everything else.
 *
 * We calculate the constraint multipliers like this:
 *           G M^-1 ~G mult = G udot0 - b, udot0=M^-1 f
 * using the pseudo inverse of G M^-1 ~G to give a least squares solution for
 * mult: mult = pinv(G M^-1 ~G)(G M^-1 f - b). Then the real udot is
 * udot = udot0 - udotC, with udotC = M^-1 ~G mult. Note: M^-1* is an
 * O(N) operator that provides the desired result; it *does not* require
 * forming or factoring M.
 *
 * NOTE: only the following constraint matrices have to be formed and factored:
 * @verbatim
 *    [G M^-1 ~G]   to calculate multipliers (square, symmetric: LDL' if
 *                  well conditioned, else pseudoinverse)
 *
 *    [P]           for projection onto position manifold (pseudoinverse)
 *
 *    [P]           for projection onto velocity manifold (pseudoinverse)
 *    [V]
 * @endverbatim
 *
 * In many cases these matrices consist of decoupled blocks which can
 * be solved independently; we try to take advantage of that whenever possible
 * to solve a set of smaller systems rather than one large one. Also, in the
 * majority of biosimulation applications we are likely to have only holonomic
 * (position) constraints, so there is no V or A and G=P is the whole story.
 */
class SimTK_SIMBODY_EXPORT SimbodyMatterSubsystem : public Subsystem {
public:
    /// Create a tree containing only the ground body (body 0).
    SimbodyMatterSubsystem();
    explicit SimbodyMatterSubsystem(MultibodySystem&);

    class Subtree; // used for working with a connected subgraph of the MobilizedBody tree
    class SubtreeResults;

    // These are the same as the compiler defaults but are handy to
    // have around explicitly for debugging.
    ~SimbodyMatterSubsystem() {
    }
    SimbodyMatterSubsystem(const SimbodyMatterSubsystem& ss) : Subsystem(ss) {
    }
    SimbodyMatterSubsystem& operator=(const SimbodyMatterSubsystem& ss) {
        Subsystem::operator=(ss);
        return *this;
    }


        ///////////////////////////////
        // PAUL'S FRIENDLY INTERFACE //
        ///////////////////////////////

    /// Calculate the total system mass.
    ///
    /// @par Required stage
    ///   \c Stage::Instance
    Real calcSystemMass(const State& s) const;


    /// Return the location r_OG_C of the system mass center C, measured from the ground
    /// origin OG, and expressed in Ground. 
    ///
    /// @par Required stage
    ///   \c Stage::Position
    Vec3 calcSystemMassCenterLocationInGround(const State& s) const;


    /// Return total system mass, mass center location measured from the Ground origin,
    /// and system inertia taken about the Ground origin, expressed in Ground.
    ///
    /// @par Required stage
    ///   \c Stage::Position
    MassProperties calcSystemMassPropertiesInGround(const State& s) const;

    /// Return the system inertia matrix taken about the system center of mass,
    /// expressed in Ground.
    ///
    /// @par Required stage
    ///   \c Stage::Position
    Inertia calcSystemCentralInertiaInGround(const State& s) const;

    /// Return the velocity V_G_C = d/dt r_OG_C of the system mass center C in the Ground frame G,
    /// expressed in G.
    ///
    /// @par Required stage
    ///   \c Stage::Velocity
    Vec3 calcSystemMassCenterVelocityInGround(const State& s) const;

    /// Return the acceleration A_G_C = d^2/dt^2 r_OG_C of the system mass center C in
    /// the Ground frame G, expressed in G.
    ///
    /// @par Required stage
    ///   \c Stage::Acceleration
    Vec3 calcSystemMassCenterAccelerationInGround(const State& s) const;

    /// Return the momentum of the system as a whole (angular, linear) measured
    /// in the ground frame, taken about the ground origin and expressed in ground.
    /// (The linear component is independent of the "about" point.)
    ///
    /// @par Required stage
    ///   \c Stage::Velocity
    SpatialVec calcSystemMomentumAboutGroundOrigin(const State& s) const;

        //////////////////
        // CONSTRUCTION //
        //////////////////

    // Attach new matter using the indicated parent body as the reference
    // frame, with the mobilizer and mass properties provided by 'child'.
    // We take over ownership of child's representation from the given
    // handle, leaving that handle as a reference to our new matter object.
    // It is an error if the given handle wasn't the owner of the
    // matter representation.
    MobilizedBodyId      adoptMobilizedBody(MobilizedBodyId parent, MobilizedBody& child);
    const MobilizedBody& getMobilizedBody(MobilizedBodyId) const;
    MobilizedBody&       updMobilizedBody(MobilizedBodyId);


    // Note: topology is not marked invalid upon returning a writable reference
    // here; that will be done only if a non-const method of the returned MobilizedBody
    // is called. That means it is OK to use Ground() to satisfy a const argument;
    // it won't have an "invalidate topology" side effect.
    const MobilizedBody::Ground& getGround() const;
    MobilizedBody::Ground&       updGround();
    MobilizedBody::Ground&       Ground() {return updGround();}

    ConstraintId      adoptConstraint(Constraint&);
    const Constraint& getConstraint(ConstraintId) const;
    Constraint&       updConstraint(ConstraintId);

        ///////////////
        // OPERATORS //
        ///////////////

        // Operators make use of the State but do not write their results back
        // into the State, not even into the State cache.

    /// This is the primary dynamics operator. It takes a state which
    /// has been realized to the Dynamics stage, a complete set of forces to apply,
    /// and returns the accelerations that result. Only the forces supplied here,
    /// and those resulting from centrifugal effects, affect the results. Everything
    /// in the matter subsystem is accounted for including velocities and 
    /// acceleration constraints, which will always be satisified as long as the
    /// constraints are consistent. If the position and velocity constraints
    /// aren't already satisified in the State, these accelerations
    /// are harder to interpret physically, but they will still be calculated and
    /// the acceleration constraints will still be satisified. No attempt
    /// will be made to satisfy position and velocity constraints, or even to check
    /// whether they are statisfied.
    /// This is an O(n*nc^2) operator worst case where all nc constraint equations
    /// are coupled.
    /// Requires prior realization through Stage::Dynamics.
    void calcAcceleration(const State&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot,
        Vector_<SpatialVec>&       A_GB) const;


    /// This operator is similar to calcAcceleration but ignores the effects of
    /// acceleration constraints. The supplied forces and velocity-induced centrifugal
    /// effects are properly accounted for, but any forces that would have resulted
    /// from enforcing the contraints are not present.
    /// This is an O(N) operator.
    /// Requires prior realization through Stage::Dynamics.
    void calcAccelerationIgnoringConstraints(const State&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot,
        Vector_<SpatialVec>&       A_GB) const;

    /// This operator calculates M^-1 v where M is the system mass matrix and v
    /// is a supplied vector with one entry per mobility. If v is a set of 
    /// mobility forces f, the result is a generalized acceleration (udot=M^-1 f). Only 
    /// the supplied vector is used, and M depends only on position states,
    /// so the result here is not affected by velocities in the State.
    /// However, this fast O(N) operator requires that the Dynamics stage operators
    /// are already available, so the State must be realized to Stage::Dynamics
    /// even though velocities are ignored.
    /// Requires prior realization through Stage::Dynamics.
    void calcMInverseV(const State&,
        const Vector&        v,
        Vector&              MinvV,
        Vector_<SpatialVec>& A_GB) const; // <-- TODO: get rid of this A

    /// Requires realization through Stage::Position.
    void calcSpatialKinematicsFromInternal(const State&,
        const Vector&        v,
        Vector_<SpatialVec>& Jv) const;

    /// Requires realization through Stage::Position.
    void calcInternalGradientFromSpatial(const State&,
        const Vector_<SpatialVec>& dEdR,
        Vector&                    dEdQ) const; // really Qbar

    /// Requires realization through Stage::Velocity.
    Real calcKineticEnergy(const State&) const;

    /// Accounts for applied forces
    /// and centrifugal forces produced by non-zero velocities in the State. Returns
    /// a set of mobility forces which replace both the applied bodyForces and the
    /// centrifugal forces.
    /// Requires prior realization through Stage::Dynamics. 
    void calcTreeEquivalentMobilityForces(const State&, 
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    mobilityForces) const;

    /// Must be in Stage::Position to calculate qdot = Q(q)*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    /// Must be in Stage::Velocity to calculate qdotdot = Q(q)*udot + Qdot(q,u)*u.
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

    /// The number of bodies includes all rigid bodies, massless
    /// bodies and ground but not particles. Bodies and their inboard mobilizers
    /// have the same number since they are grouped together as a MobilizedBody
    /// MobilizedBody numbering starts with ground at 0 with a regular labeling such
    /// that children have higher body numbers than their parents. Mobilizer 0
    /// is meaningless (or I suppose you could think of it as the weld
    /// joint that attaches ground to the universe), but otherwise 
    /// mobilizer n is the inboard mobilizer of body n.
    int getNBodies() const;

    /// This is the total number of defined constraints, each of which may
    /// generate more than one constraint equation.
    int getNConstraints() const;

    /// TODO: total number of particles.
    int getNParticles() const;

    /// The sum of all the joint degrees of freedom. This is also the length
    /// of state variable vector u.
    int getNMobilities() const; 

    /// The sum of all the q vector allocations for each joint. There may be
    /// some that are not in use for particular modeling options.
    int getTotalQAlloc() const;

    /// This is the sum of all the allocations for constraint multipliers,
    /// one per acceleration constraint equation.
    int getTotalMultAlloc() const;

    /// For all mobilizers offering unrestricted orientation, decide what
    /// method we should use to model their orientations. Choices are: quaternions (best
    /// for dynamics), or rotation angles (1-2-3 Euler sequence, good for
    /// optimization). TODO: (1) other Euler sequences, (2) allow settable zero
    /// rotation for Euler sequence, with convenient way to say "this is zero".
    void setUseEulerAngles(State&, bool) const;
    bool getUseEulerAngles  (const State&) const;
    int  getNQuaternionsInUse(const State&) const;

    // TODO: these are obsolete. Their functions should be
    // taken over by methods in the MobilizedBody and Constraint classes.
    void setMobilizerIsPrescribed(State&, MobilizedBodyId, bool) const;
    bool isMobilizerPrescribed  (const State&, MobilizedBodyId) const;
    bool isUsingQuaternion(const State&, MobilizedBodyId) const;
    int  getQuaternionIndex(const State&, MobilizedBodyId) const;
    void setConstraintIsEnabled(State&, ConstraintId constraint, bool) const;
    bool isConstraintEnabled(const State&, ConstraintId constraint) const;

    // Position Stage. 

    // Dynamics stage responses.

    // Cross joint
    const SpatialVec& getCoriolisAcceleration(const State&, MobilizedBodyId) const;

    // Including parent
    const SpatialVec& getTotalCoriolisAcceleration(const State&, MobilizedBodyId) const;

    const SpatialVec& getGyroscopicForce(const State&, MobilizedBodyId) const;
    const SpatialVec& getCentrifugalForces(const State&, MobilizedBodyId) const;
    const SpatialMat& getArticulatedBodyInertia(const State& s, MobilizedBodyId) const;

        // PARTICLES

    // The generalized coordinates for a particle are always the three measure numbers
    // (x,y,z) of the particle's Ground-relative Cartesian location vector. The generalized
    // speeds are always the three corresponding measure numbers of the particle's
    // Ground-relative Cartesian velocity. The generalized applied forces are
    // always the three measure numbers of a Ground-relative force vector.
    const Vector_<Vec3>& getAllParticleLocations    (const State&) const;
    const Vector_<Vec3>& getAllParticleVelocities   (const State&) const;

    const Vec3& getParticleLocation(const State& s, ParticleId p) const {
        return getAllParticleLocations(s)[p];
    }
    const Vec3& getParticleVelocity(const State& s, ParticleId p) const {
        return getAllParticleVelocities(s)[p];
    }

    Vector& updAllParticleMasses(State& s) const;

    void setAllParticleMasses(State& s, const Vector& masses) const {
        updAllParticleMasses(s) = masses;
    }

    // Note that particle generalized coordinates, speeds, and applied forces
    // are defined to be the particle Cartesian locations, velocities, and
    // applied force vectors, so can be set directly at Stage::Model or higher.

    // These are the only routines that must be provided by the concrete MatterSubsystem.
    Vector_<Vec3>& updAllParticleLocations(State&)     const;
    Vector_<Vec3>& updAllParticleVelocities(State&)    const;

    // The following inline routines are provided by the generic MatterSubsystem class
    // for convenience.

    Vec3& updParticleLocation(State& s, ParticleId p) const {
        return updAllParticleLocations(s)[p];
    }
    Vec3& updParticleVelocity(State& s, ParticleId p) const {
        return updAllParticleVelocities(s)[p];
    }

    void setParticleLocation(State& s, ParticleId p, const Vec3& r) const {
        updAllParticleLocations(s)[p] = r;
    }
    void setParticleVelocity(State& s, ParticleId p, const Vec3& v) const {
        updAllParticleVelocities(s)[p] = v;
    }

    void setAllParticleLocations(State& s, const Vector_<Vec3>& r) const {
        updAllParticleLocations(s) = r;
    }
    void setAllParticleVelocities(State& s, const Vector_<Vec3>& v) const {
        updAllParticleVelocities(s) = v;
    }

    /// TODO: not implemented yet; particles must be treated as rigid bodies for now.
    const Vector& getAllParticleMasses(const State&) const;

    const Vector_<Vec3>& getAllParticleAccelerations(const State&) const;

    const Vec3& getParticleAcceleration(const State& s, ParticleId p) const {
        return getAllParticleAccelerations(s)[p];
    }

        // POSITION STAGE responses //

        // POSITION STAGE operators //

    /// Apply a force to a point on a body (a station). Provide the
    /// station in the body frame, force in the ground frame. Must
    /// be realized to Position stage prior to call.
    void addInStationForce(const State&, MobilizedBodyId bodyB, const Vec3& stationOnB, 
                           const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const;

    /// Apply a torque to a body. Provide the torque vector in the
    /// ground frame.
    void addInBodyTorque(const State&, MobilizedBodyId, const Vec3& torqueInG, 
                         Vector_<SpatialVec>& bodyForces) const;

    /// Apply a scalar joint force or torque to an axis of the
    /// indicated body's mobilizer.
    void addInMobilityForce(const State&, MobilizedBodyId, int axis, Real f,
                            Vector& mobilityForces) const;

        // POSITION STAGE solvers //

    /// This is a solver you can call after the State has been realized
    /// to stage Position. It will project the Q constraints
    /// along the error norm so that getQConstraintNorm() <= tol, and will
    /// project out the corresponding component of y_err so that y_err's Q norm
    /// is reduced. Returns true if it does anything at all to State or y_err.
    bool projectQConstraints(State&, Vector& y_err, Real tol, Real targetTol) const;

        // VELOCITY STAGE responses //

        // VELOCITY STAGE operators //

        // VELOCITY STAGE solvers //

    /// This is a solver you can call after the State has been realized
    /// to stage Velocity. It will project the U constraints
    /// along the error norm so that getUConstraintNorm() <= tol, and will
    /// project out the corresponding component of y_err so that y_err's U norm
    /// is reduced.
    bool projectUConstraints(State&, Vector& y_err, Real tol, Real targetTol) const;

        // ACCELERATION STAGE reponses


        // Bookkeeping
    SimTK_PIMPL_DOWNCAST(SimbodyMatterSubsystem, Subsystem);
    const SimbodyMatterSubsystemRep& getRep() const;
    SimbodyMatterSubsystemRep&       updRep();
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem&);


} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
