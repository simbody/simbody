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
 * Portions copyright (c) 2006-9 Stanford University and the Authors.         *
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

class MobilizedBody;
class MultibodySystem;
class Constraint;

/**
 * The Simbody low-level multibody tree interface.
 * Equations represented:
 * <pre>
 *                  qdot = N u
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
 *       pdotdot = P udot - bp(t,q,u) = 0
 *          vdot = V udot - bv(t,q,u) = 0
 * a(t,q,u,udot) = A udot - ba(t,q,u) = 0
 *           
 *                pdot = P u - c(t,q) = 0
 *                           v(t,q,u) = 0
 *
 *                             p(t,q) = 0
 *                               n(q) = 0
 * </pre>
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
 * (We can get partial(p)/partial(q) when we need it as P*N^-1.)
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
 *    [P N^-1]      for projection onto position manifold (pseudoinverse)
 *
 *    [P;V]         for projection onto velocity manifold (pseudoinverse)
 *                  (using Matlab notation meaning rows of P over rows of V)
 * @endverbatim
 *
 * When working in a weighted norm with weights W on the state variables and
 * weights T (1/tolerance) on the constraint errors, the matrices we need are
 * actually [Tp PN^-1 Wq^-1], [Tpv [P;V] Wu^-1], etc. with T and W diagonal
 * weighting matrices. These can then be used to find least squares solutions
 * in the weighted norms.
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
    
    /// Get whether default decorative geometry is displayed for bodies in this system.
    bool getShowDefaultGeometry() const;

    /// Set whether default decorative geometry is displayed for bodies in this system.
    void setShowDefaultGeometry(bool show);


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
    MobilizedBodyIndex   adoptMobilizedBody(MobilizedBodyIndex parent, MobilizedBody& child);
    const MobilizedBody& getMobilizedBody(MobilizedBodyIndex) const;
    MobilizedBody&       updMobilizedBody(MobilizedBodyIndex);


    // Note: topology is not marked invalid upon returning a writable reference
    // here; that will be done only if a non-const method of the returned MobilizedBody
    // is called. That means it is OK to use Ground() to satisfy a const argument;
    // it won't have an "invalidate topology" side effect.
    const MobilizedBody::Ground& getGround() const;
    MobilizedBody::Ground&       updGround();
    MobilizedBody::Ground&       Ground() {return updGround();}

    ConstraintIndex   adoptConstraint(Constraint&);
    const Constraint& getConstraint(ConstraintIndex) const;
    Constraint&       updConstraint(ConstraintIndex);

        ///////////////
        // OPERATORS //
        ///////////////

        // Operators make use of the State but do not write their results back
        // into the State, not even into the State cache.

    /**
     * This is the primary forward dynamics operator. It takes a state which
     * has been realized to the Dynamics stage, a complete set of forces to apply,
     * and returns the accelerations that result. Only the forces supplied here,
     * and those resulting from centrifugal effects, affect the results. Everything
     * in the matter subsystem is accounted for including velocities and 
     * acceleration constraints, which will always be satisified as long as the
     * constraints are consistent. If the position and velocity constraints
     * aren't already satisified in the State, these accelerations
     * are harder to interpret physically, but they will still be calculated and
     * the acceleration constraints will still be satisfied. No attempt
     * will be made to satisfy position and velocity constraints, or even to check
     * whether they are statisfied.
     * This operator solves the two equations
     * <pre>
     *      M udot + G^T lambda = F
     *      G udot              = b
     * </pre>
     * for udot and lambda, but does not return lambda. F includes both the 
     * applied forces and the "bias" forces due to rigid body rotations.
     * This is an O(n*nc^2) operator worst case where all nc constraint equations
     * are coupled.
     * Requires prior realization through Stage::Dynamics.
     */
    void calcAcceleration(const State&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot,
        Vector_<SpatialVec>&       A_GB) const;

    /**
     * This operator is similar to calcAcceleration but ignores the effects of
     * acceleration constraints. The supplied forces and velocity-induced centrifugal
     * effects are properly accounted for, but any forces that would have resulted
     * from enforcing the contraints are not present.
     * This operator solves the equation
     * <pre>
     *     M udot = F
     * </pre>
     * for udot. F includes both the applied forces and the "bias" forces due
     * to rigid body rotations, but does not include any constraint forces.
     * This is an O(N) operator.
     * Requires prior realization through Stage::Dynamics.
     */
    void calcAccelerationIgnoringConstraints(const State&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector&                    udot,
        Vector_<SpatialVec>&       A_GB) const;

    /**
     * This operator calculates M^-1 v where M is the system mass matrix and v
     * is a supplied vector with one entry per mobility. If v is a set of 
     * mobility forces f, the result is a generalized acceleration (udot=M^-1 f). Only 
     * the supplied vector is used, and M depends only on position states,
     * so the result here is not affected by velocities in the State.
     * However, this fast O(N) operator requires that the Dynamics stage operators
     * are already available, so the State must be realized to Stage::Dynamics
     * even though velocities are ignored.
     * Requires prior realization through Stage::Dynamics.
     */
    void calcMInverseV(const State&,
        const Vector&        v,
        Vector&              MinvV) const;

    /**
     * NOT IMPLEMENTED YET --
     * This is the primary inverse dynamics operator. Using position and velocity
     * from the given state, a set of applied forces, and a known set of mobility
     * accelerations and constraint multipliers, it calculates the additional
     * mobility forces that would be required to satisfy Newton's 2nd law.
     * That is, this operator returns
     * <pre>
     *     f_residual = M udot + G^T lambda + f_inertial - f_applied
     * </pre>
     * where f_applied is the mobility-space equivalent to all the
     * applied forces (including mobility and body forces), f_inertial 
     * is the mobility-space equivalent of the velocity-dependent
     * inertial forces due to rigid body rotations (coriolis  
     * and gyroscopic forces), and the udots and lambdas are given values of the 
     * generalized accelerations and constraint multipliers, resp.
     * TODO
     * @see calcResidualForceIgnoringConstraints()
     */
    void calcResidualForce
       (const State&               state,
        const Vector&              appliedMobilityForces,
        const Vector_<SpatialVec>& appliedBodyForces,
        const Vector&              knownUdot,
        const Vector&              knownMultipliers,
        Vector&                    residualMobilityForces) const;

    /**
     * This is the inverse dynamics operator for the tree system; if there are
     * any constraints they are ignored. This method solves
     * <pre>
     *      f_residual = M udot + f_inertial - f_applied
     * </pre> 
     * in O(n) time, meaning that the mass matrix M is never formed. Inverse
     * dynamics is considerably faster than forward dynamics (even though that
     * is also O(n) in Simbody).
     *
     * In the above equation we solve for the residual forces \c f_residual given
     * desired accelerations and (optionally) a set of applied forces. Here 
     * \c f_applied is the mobility-space equivalent of all the applied forces
     * (including mobility and body forces), \c f_inertial is the mobility-space
     * equivalent of the velocity-dependent inertial forces due to rigid 
     * body rotations (coriolis and gyroscopic forces), and \c udot is the 
     * given set of values for the desired generalized accelerations. The returned 
     * \c f_residual is the additional generalized force (that is, mobilizer 
     * force) that would have to be applied at each mobility to give the desired
     * \c udot. The inertial forces depend on the velocities \c u already realized 
     * in the State. Otherwise, only the explicitly-supplied forces affect the 
     * results of this operator; any forces that may be present elsewhere in 
     * the System are ignored.
     *
     * @param[in] state
     *      A State valid for the containing System, already realized to
     *      Stage::Velocity.
     * @param[in] appliedMobilityForces
     *      One scalar generalized force applied per mobility. Can be zero
     *      length if there are no mobility forces; otherwise must have exactly 
     *      one entry per mobility in the matter subsystem.
     * @param[in] appliedBodyForces
     *      One spatial force for each body. A spatial force is a force applied
     *      to the body origin and a torque on the body, each expressed in the 
     *      Ground frame. The supplied Vector must be either zero length or have 
     *      exactly one entry per body in the matter subsystem.
     * @param[in] knownUdot
     *      These are the desired generalized accelerations, one per mobility. 
     *      If this is zero length it will be treated as all-zero; otherwise 
     *      it must have exactly one entry per mobility in the matter subsystem.
     * @param[out] residualMobilityForces
     *      These are the residual generalized forces which, if applied, would 
     *      produce the \p knownUdot. This will be resized if necessary to have 
     *      one scalar entry per mobility. 
     *
     * @see calcResidualForce(), calcMV()
     * @see calcAcceleration(), calcAccelerationIgnoringConstraints()
     */
    void calcResidualForceIgnoringConstraints
       (const State&               state,
        const Vector&              appliedMobilityForces,
        const Vector_<SpatialVec>& appliedBodyForces,
        const Vector&              knownUdot,
        Vector&                    residualMobilityForces) const;

    /// This operator calculate M*v where M is the system mass matrix and
    /// v is a supplied vector with one entry per mobility. If v is a set
    /// of mobility accelerations (generalized accelerations), then the
    /// result is a generalized force (f=M*a). Only the supplied vector
    /// is used, and M depends only on position states, so the result
    /// here is not affected by velocities in the State.
    /// Requires prior realization through Stage::Position.
    void calcMV(const State&, const Vector& v, Vector& MV) const;

    /// Returns
    ///     constraintErr = G udot - b
    /// the residual error in the acceleration constraint equation given
    /// a generalized acceleration udot.
    /// Requires velocites to have been realized so that b(t,q,u) is 
    /// available.
    void calcAccConstraintErr(const State&,
        const Vector&   knownUdot,
        Vector&         constraintErr) const;

    /// Returns G*v, the product of the mXn acceleration constraint Jacobian
    /// and a vector of length n. m is the number of active acceleration constraint equations,
    /// n is the number of mobilities.
    /// This is an O(n) operation.
    void calcGV(const State&,
        const Vector&   v,
        Vector&         Gv) const;

    /// Returns G^T*v, the product of the nXm transpose of the acceleration constraint Jacobian
    /// G and a vector v of length m. m is the number of active acceleration constraint equations,
    /// n is the number of mobilities. If v is a set of constraint multipliers, then
    /// f=G^T*v is the set of equivalent generalized forces they generate.
    /// This is an O(n) operation.
    void calcGtV(const State&,
        const Vector&   v,
        Vector&         GtV) const;


    /// This operator explicitly calcuates the n X n mass matrix M. Note that
    /// this is inherently an O(n^2) operation since the mass matrix has
    /// n^2 elements (although only n(n+1)/2 are unique due to symmetry).
    /// <em>DO NOT USE THIS CALL DURING NORMAL DYNAMICS</em>. To
    /// do so would change an O(n) operation into an O(n^2) one. Instead,
    /// see if you can accomplish what you need with O(n) operators like
    /// calcMV() which calculates the matrix-vector product M*v in O(n)
    /// without explicitly forming M.
    /// Also, don't invert this matrix numerically to get M^-1. Instead, call
    /// the method calcMInv() which can produce M^-1 directly.
    /// @see calcMV()
    /// @see calcMInv()
    void calcM(const State&, Matrix& M) const;

    /// This operator explicitly calculates the n X n mass matrix inverse
    /// M^-1. This is an O(n^2) operation, which is of course within a 
    /// constant factor of optimal for returning a matrix with n^2 
    /// elements explicitly. (There are actually only n(n+1)/2 unique
    /// elements since the matrix is symmetric.)
    /// <em>DO NOT USE THIS CALL DURING NORMAL DYNAMICS</em>. To
    /// do so would change an O(n) operation into an O(n^2) one. Instead,
    /// see if you can accomplish what you need with O(n) operators like
    /// calcMInvV() which calculates the matrix-vector product M^-1*v in O(n)
    /// without explicitly forming M or M^-1.
    /// If you need M explicitly, you can get it with the calcM() method.
    /// @see calcMInvV()
    /// @see calcM()
    void calcMInv(const State&, Matrix& MInv) const;

    /// This O(nm) operator explicitly calculates the m X n acceleration-level constraint
    /// Jacobian G = [P;V;A] which appears in the system equations of motion.
    /// This method generates G columnwise use the acceleration-level
    /// constraint error equations.
    /// To within numerical error, this should be identical to the transpose of
    /// the matrix returned by calcGt() which uses a different method.
    /// Consider using the calcGV() method instead of this one, which forms the
    /// matrix-vector product G*v in O(n) time without explicitly forming G.
    /// @see calcGt()
    /// @see calcGV()
    void calcG(const State&, Matrix& G) const;

    /// This O(nm) operator explicitly calculates the n X m transpose of the acceleration-
    /// level constraint Jacobian G = [P;V;A] which appears in the system equations
    /// of motion. This method generates G^T columnwise use the constraint force
    /// generating methods which map constraint multipliers to constraint forces.
    /// To within numerical error, this should be identical to the transpose of
    /// the matrix returned by calcG() which uses a different method.
    /// Consider using the calcGtV() method instead of this one, which forms the
    /// matrix-vector product G^T*v in O(n) time without explicitly forming G^T.
    /// @see calcG()
    /// @see calcGtV()
    void calcGt(const State&, Matrix& Gt) const;



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

    /// Must be in Stage::Position to calculate qdot = N(q)*u.
    void calcQDot(const State& s,
        const Vector& u,
        Vector&       qdot) const;

    /// Must be in Stage::Velocity to calculate qdotdot = N(q)*udot + Ndot(q,u)*u.
    void calcQDotDot(const State& s,
        const Vector& udot,
        Vector&       qdotdot) const;

    /// Must be in Stage::Position to calculate out_q = N(q)*in_u (e.g., qdot=N*u)
    /// or out_u = ~N*in_q. Note that one of "in" and "out" is always "q-like" while
    /// the other is "u-like", but which is which changes if the matrix is transposed.
    /// Note that the transposed operation here is the same as multiplying by N on
    /// the right, with the Vectors viewed as RowVectors instead.
    /// This is an O(n) operator since N is block diagonal.
    void multiplyByN(const State& s, bool transpose, const Vector& in, Vector& out) const;

    /// Must be in Stage::Position to calculate out_u = NInv(q)*in_q (e.g., u=NInv*qdot)
    /// or out_q = ~NInv*in_u. Note that one of "in" and "out" is always "q-like" while
    /// the other is "u-like", but which is which changes if the matrix is transposed.
    /// Note that the transposed operation here is the same as multiplying by NInv on
    /// the right, with the Vectors viewed as RowVectors instead.
    /// This is an O(N) operator since NInv is block diagonal.
    void multiplyByNInv(const State& s, bool transpose, const Vector& in, Vector& out) const;

    /// Calculate the mobilizer reaction force generated by each MobilizedBody.  This is the
    /// constraint force that would be required to make the system move in the same way if that
    /// MobilizedBody were converted to a Free body.  A mobilizer exerts equal and opposite reaction
    /// forces on the parent and child bodies.  This method reports the force on the child body.
    /// The force is applied at the origin of the outboard frame M (fixed to the child), and expressed
    /// in the ground frame.
    ///
    /// The State must have been realized to Stage::Acceleration to use this method.
    void calcMobilizerReactionForces(const State& s, Vector_<SpatialVec>& forces) const;

    // These are available after realizeTopology().

    /// The number of bodies includes all rigid bodies, massless
    /// bodies and ground but not particles. Bodies and their inboard mobilizers
    /// have the same number since they are grouped together as a MobilizedBody
    /// MobilizedBody numbering starts with ground at 0 with a regular labeling such
    /// that children have higher body numbers than their parents. Mobilizer 0
    /// is meaningless (or I suppose you could think of it as the weld
    /// joint that attaches ground to the universe), but otherwise 
    /// mobilizer n is the inboard mobilizer of body n.
    int getNumBodies() const;

    /// This is the total number of defined constraints, each of which may
    /// generate more than one constraint equation.
    int getNumConstraints() const;

    /// TODO: total number of particles.
    int getNumParticles() const;

    /// The sum of all the mobilizer degrees of freedom. This is also the length
    /// of the state variable vector u and the mobility forces array.
    int getNumMobilities() const;

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

    int  getNumQuaternionsInUse(const State&) const;

    void setMobilizerIsPrescribed(State&, MobilizedBodyIndex, bool) const;
    bool isMobilizerPrescribed  (const State&, MobilizedBodyIndex) const;
    bool isUsingQuaternion(const State&, MobilizedBodyIndex) const;
    QuaternionPoolIndex getQuaternionPoolIndex(const State&, MobilizedBodyIndex) const;
    AnglePoolIndex      getAnglePoolIndex(const State&, MobilizedBodyIndex) const;
    void setConstraintIsDisabled(State&, ConstraintIndex constraint, bool) const;
    bool isConstraintDisabled(const State&, ConstraintIndex constraint) const;
    
    /// Given a State which is modeled using quaternions, convert it to a
    /// representation based on Euler angles and store the result in another state.
    void convertToEulerAngles(const State& inputState, State& outputState) const;
    
    /// Given a State which is modeled using Euler angles, convert it to a
    /// representation based on quaternions and store the result in another state.
    void convertToQuaternions(const State& inputState, State& outputState) const; 

    // Dynamics stage responses.

    // Cross joint
    const SpatialVec& getCoriolisAcceleration(const State&, MobilizedBodyIndex) const;

    // Including parent
    const SpatialVec& getTotalCoriolisAcceleration(const State&, MobilizedBodyIndex) const;

    const SpatialVec& getGyroscopicForce(const State&, MobilizedBodyIndex) const;
    const SpatialVec& getCentrifugalForces(const State&, MobilizedBodyIndex) const;
    const SpatialMat& getArticulatedBodyInertia(const State& s, MobilizedBodyIndex) const;

        // PARTICLES
        // TODO: not currently implemented. Use a point mass with a Cartesian (translation)
        // mobilizer to Ground instead. The idea here would be to special-case particles
        // to make them faster; there would be no additional functionality.

    // The generalized coordinates for a particle are always the three measure numbers
    // (x,y,z) of the particle's Ground-relative Cartesian location vector. The generalized
    // speeds are always the three corresponding measure numbers of the particle's
    // Ground-relative Cartesian velocity. The generalized applied forces are
    // always the three measure numbers of a Ground-relative force vector.
    const Vector_<Vec3>& getAllParticleLocations    (const State&) const;
    const Vector_<Vec3>& getAllParticleVelocities   (const State&) const;

    const Vec3& getParticleLocation(const State& s, ParticleIndex p) const {
        return getAllParticleLocations(s)[p];
    }
    const Vec3& getParticleVelocity(const State& s, ParticleIndex p) const {
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

    Vec3& updParticleLocation(State& s, ParticleIndex p) const {
        return updAllParticleLocations(s)[p];
    }
    Vec3& updParticleVelocity(State& s, ParticleIndex p) const {
        return updAllParticleVelocities(s)[p];
    }

    void setParticleLocation(State& s, ParticleIndex p, const Vec3& r) const {
        updAllParticleLocations(s)[p] = r;
    }
    void setParticleVelocity(State& s, ParticleIndex p, const Vec3& v) const {
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

    const Vec3& getParticleAcceleration(const State& s, ParticleIndex p) const {
        return getAllParticleAccelerations(s)[p];
    }

        // POSITION STAGE responses //

        // POSITION STAGE operators //

    /// Apply a force to a point on a body (a station). Provide the
    /// station in the body frame, force in the ground frame. Must
    /// be realized to Position stage prior to call.
    void addInStationForce(const State&, MobilizedBodyIndex bodyB, const Vec3& stationOnB, 
                           const Vec3& forceInG, Vector_<SpatialVec>& bodyForcesInG) const;

    /// Apply a torque to a body. Provide the torque vector in the
    /// ground frame.
    void addInBodyTorque(const State&, MobilizedBodyIndex, const Vec3& torqueInG, 
                         Vector_<SpatialVec>& bodyForcesInG) const;

    /// Apply a scalar joint force or torque to an axis of the
    /// indicated body's mobilizer.
    void addInMobilityForce(const State&, MobilizedBodyIndex, MobilizerUIndex which, Real f,
                            Vector& mobilityForces) const;

        // POSITION STAGE solvers //

    /// This is a solver you can call after the State has been realized
    /// to stage Position. It will project the q constraints
    /// along the error norm so that getQConstraintNorm() <= consAccuracy, and will
    /// project out the corresponding component of yErrest so that yErrest's q norm
    /// is reduced. Returns true if it does anything at all to State or yErrest.
    bool projectQConstraints(State& s, Real consAccuracy, const Vector& yWeights,
							 const Vector& ooTols, Vector& yErrest, System::ProjectOptions) const;

        // VELOCITY STAGE responses //

        // VELOCITY STAGE operators //

        // VELOCITY STAGE solvers //

    /// This is a solver you can call after the State has been realized
    /// to stage Velocity. It will project the u constraints
    /// along the error norm so that getUConstraintNorm() <= consAccuracy, and will
    /// project out the corresponding component of yErrest so that yErrest's u norm
    /// is reduced. Returns true if it does anything at all to State or yErrest.
    bool projectUConstraints(State& s, Real consAccuracy, const Vector& yWeights,
							 const Vector& ooTols, Vector& yErrest, System::ProjectOptions) const;

        // ACCELERATION STAGE reponses


        // Bookkeeping
    SimTK_PIMPL_DOWNCAST(SimbodyMatterSubsystem, Subsystem);
    const SimbodyMatterSubsystemRep& getRep() const;
    SimbodyMatterSubsystemRep&       updRep();

    // OBSOLETE; TODO: remove in SimTK 2.0
    void multiplyByQMatrix(const State& s, bool transposeMatrix, const Vector& in, Vector& out) const
    {   multiplyByN(s,transposeMatrix,in,out);}
    // OBSOLETE; TODO: remove in SimTK 2.0
    void multiplyByQMatrixInverse(const State& s, bool transposeMatrix, const Vector& in, Vector& out) const
    {   multiplyByNInv(s,transposeMatrix,in,out);}
    // OBSOLETE; TODO: remove in SimTK 2.0
    int getNBodies() const {return getNumBodies();}
    // OBSOLETE; TODO: remove in SimTK 2.0
    int getNConstraints() const {return getNumConstraints();}
    // OBSOLETE; TODO: remove in SimTK 2.0
    int getNMobilities() const {return getNumMobilities();}
    // OBSOLETE; TODO: remove in SimTK 2.0
    int getNParticles() const {return getNumParticles();}
private:
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem&);


} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
