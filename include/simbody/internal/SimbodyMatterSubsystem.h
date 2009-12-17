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
 * where M(q) is the mass matrix, G(q) the acceleration constraint matrix, 
 * C(q,u) the coriolis and gyroscopic forces, T is user-applied joint mobility
 * forces, F is user-applied body forces and torques and gravity. J* is the 
 * operator that maps spatial forces to joint mobility forces. p() are the
 * holonomic (position) constraints, v() the non-holonomic (velocity) 
 * constraints, and a() the acceleration-only constraints, which must be 
 * linear, with A the coefficient matrix for a(). pdot, pdotdot are obtained
 * by differentiation of p(), vdot by differentiation of v().
 * P=partial(pdot)/partial(u) (yes, that's u, not q), V=partial(v)/partial(u).
 * (We can get partial(p)/partial(q) when we need it as P*N^-1.) n(q) is the 
 * set of quaternion normalization constraints, which exist only at the
 * position level and are uncoupled from everything else.
 *
 * We calculate the constraint multipliers like this:
 * <pre>
 *           G M^-1 ~G mult = G udot0 - b, udot0=M^-1 f
 * </pre>
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

/// Get whether default decorative geometry is displayed for bodies in 
/// this system.
bool getShowDefaultGeometry() const;

/// Set whether default decorative geometry is displayed for bodies in 
/// this system.
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

/// Return the momentum of the system as a whole (angular, linear) measured
/// in the ground frame, taken about the current system center of mass
/// location and expressed in ground.
/// (The linear component is independent of the "about" point.)
///
/// @par Required stage
///   \c Stage::Velocity
SpatialVec calcSystemCentralMomentum(const State& s) const;

    //////////////////
    // CONSTRUCTION //
    //////////////////

/// Attach new matter by attaching it to the indicated parent body. The 
/// mobilizer and mass properties provided by \a child. A new 
/// MobilizedBodyIndex is assigned for the child; it is guaranteed to be 
/// numerically larger than the MobilizedBodyIndex of the parent. We take 
/// over ownership of child's implementation object from the given 
/// MobilizedBody handle, leaving that handle as a reference to the 
/// implementation object now owned by the matter subsystem. It is an 
/// error if the given MobilizedBody handle wasn't the owner of the 
/// implementation object to which it refers.
/// @note
/// This method is usually called by concrete MobilizedBody constructors;
/// it does not normally need to be called by end users.
MobilizedBodyIndex   adoptMobilizedBody(MobilizedBodyIndex  parent, 
                                        MobilizedBody&      child);

/// Given a MobilizedBodyIndex, return a read-only (const) reference to the 
/// corresponding MobilizedBody within this matter subsystem. This method 
/// will fail if no such MobilizedBody is present.
const MobilizedBody& getMobilizedBody(MobilizedBodyIndex) const;

/// Given a MobilizedBodyIndex, return a writable reference to the 
/// corresponding MobilizedBody within this matter subsystem. This method 
/// will fail if no such MobilizedBody is present.
MobilizedBody&       updMobilizedBody(MobilizedBodyIndex);


// Note: topology is not marked invalid upon returning a writable reference
// here; that will be done only if a non-const method of the returned 
// MobilizedBody is called. That means it is OK to use Ground() to satisfy 
// a const argument; it won't have an "invalidate topology" side effect.

/// Return a read-only (const) reference to the Ground MobilizedBody
/// within this matter subsystem.
const MobilizedBody::Ground& getGround() const;
/// Return a writable reference to the Ground MobilizedBody within this
/// matter subsystem.
MobilizedBody::Ground&       updGround();
/// This is a synonym for updGround() that makes for nicer-looking examples.
/// @see updGround()
MobilizedBody::Ground&       Ground() {return updGround();}

/// Add a new Constraint object to the matter subsystem. The details of
/// the Constraint are opaque here. A new ConstraintIndex is assigned.
/// We take  over ownership of the implementation object from the given 
/// Constraint handle, leaving that handle as a reference to the 
/// implementation object now owned by the matter subsystem. It is an 
/// error if the given Constraint handle wasn't the owner of the 
/// implementation object to which it refers.
/// @note
/// This method is usually called by concrete Constraint constructors;
/// it does not normally need to be called by end users.
ConstraintIndex   adoptConstraint(Constraint&);

/// Given a ConstraintIndex, return a read-only (const) reference to the 
/// corresponding Constraint within this matter subsystem. This method 
/// will fail if no such Constraint is present.
const Constraint& getConstraint(ConstraintIndex) const;

/// Given a ConstraintIndex, return a writable reference to the 
/// corresponding Constraint within this matter subsystem. This method 
/// will fail if no such Constraint is present.
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
 * and those calculated internally from prescribed motion, constraints, and
 * centrifugal effects, affect the results. Acceleration constraints are 
 * always satisfied on return as long as the constraints are consistent. 
 * If the position and velocity constraints aren't already satisified in the 
 * State, results are harder to interpret physically, but they will still be 
 * calculated and the acceleration constraints will still be satisfied. No 
 * attempt will be made to satisfy position and velocity constraints, or to 
 * set prescribed positions and velocities, nor even to check whether these 
 * are satisfied; position and velocity constraint and prescribed positions 
 * and velocities are simply irrelevant here.
 *
 * Given applied forces f_applied, this operator solves this set of equations:
 * <pre>
 *      M udot + G^T lambda + f_bias = f_applied    (1)
 *      G udot              - b      = 0            (2)
 * </pre>
 * for udot and lambda (although it does not return lambda). 
 * f_applied is the set of generalized (mobility) forces equivalent to the 
 * \a mobilityForces and \a bodyForces arguments supplied here (in particular, 
 * it does \e not include forces due to prescribed motion).
 * M, G, and b are defined by the mobilized bodies, constraints, and prescribed 
 * motions present in the System. f_bias includes 
 * the velocity-dependent gyroscopic and coriolis forces due to rigid body 
 * rotations and is extracted internally from the already-realized state. 
 *
 * Prescribed accelerations are treated logically as constraints included in 
 * equation (2), although the corresponding part of G is an identity matrix. 
 * Thus the generalized forces used to enforce prescribed motion are a subset 
 * of the lambdas. Note that this method does not allow you to specify your 
 * own prescribed udots; those are calculated from the mobilizers' 
 * state-dependent Motion specifications that are already part of the System.
 *
 * This is an O(n*m^2) operator worst case where all m constraint equations
 * are coupled (prescribed motions are counted in n, not m). Requires prior 
 * realization through Stage::Dynamics.
 */
void calcAcceleration(const State&,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot,    // output only; no prescribed motions
    Vector_<SpatialVec>&       A_GB) const;

/**
 * This operator solves this set of equations:
 * <pre>
 *      M udot + G^T lambda = f_applied - f_bias    (1)
 *      G udot              = b                     (2)
 * </pre>
 * for udot and lambda given udot_p and f_applied, where udot={udot_p, udot_r} 
 * and lambda={lambda_p, lambda_r}, with suffix "_p" meaning "prescribed" and 
 * "_r" meaning "regular".
 *
 * We're ignoring any forces and motions that are part of the System, except
 * to determine which of the mobilizers are prescribed, in which case their
 * accelerations must be supplied here.
 *
 * Notes:
 *  - the complete set of udots is returned (udot_p will
 *    have been copied into the appropriate slots).
 *  - lambda is returned in separate lambda_p and lambda_r Vectors.
 *  - lambda_p's are mobility forces but have to be mapped to
 *    the mobilities
 *  - lambda_p's are \e residual forces; they have the opposite sign from
 *    \e applied mobility forces
 *  - lambda_r's must be mapped to residual mobility forces via G^T
 */
void calcDynamicsFromForcesAndMotions(const State&,
    const Vector&              appliedMobilityForces,  // [nu] or [0]
    const Vector_<SpatialVec>& appliedBodyForces,      // [nb] or [0]
    const Vector&              udot_p,                 // [npres] or [0]
    Vector&                    udot,                   // [nu]    (output only)
    Vector&                    lambda_p,               // [npres] (output only)
    Vector&                    lambda_r) const;        // [m]     (output only)

/**
 * This is the same as calcDynamicsFromForcesAndMotions() except that the
 * Motions are taken from those in the System rather than from arguments.
 * All applied forces in the System are ignored; they are overridden by
 * the arguments here. Otherwise everything in the discussion of 
 * calcDynamicsFromForcesAndMotions() applies here too.
 * @see calcDynamicsFromForcesAndMotions()
 */
void calcDynamicsFromForces(const State&,
    const Vector&              appliedMobilityForces,  // [nu] or [0]
    const Vector_<SpatialVec>& appliedBodyForces,      // [nb] or [0]
    Vector&                    udot,                   // [nu]    (output only)
    Vector&                    lambda_p,               // [npres] (output only)
    Vector&                    lambda_r) const;        // [m]     (output only) 

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
 * This operator calculates in O(N) time the product M^-1*v where M is the 
 * system mass matrix and v is a supplied vector with one entry per 
 * mobility. If v is a set of generalized forces f, the result is a 
 * generalized acceleration (udot=M^-1*f). Only the supplied vector is 
 * used, M depends only on position states, so the result here is not 
 * affected by velocities in the State. In particular, you'll have to 
 * calculate your own inertial forces and put them in f if you want 
 * them included.
 *
 * If the supplied State does not already contain realized values for the
 * articulated body inertias, then they will be realized when this operator
 * is first called for a new set of positions. When there are prescribed 
 * accelerations, articulated body inertias are created using mixed 
 * articulated and composite inertia calculations depending on the 
 * prescribed status of the corresponding mobilizers. In that case this
 * method works only with the "free" (non-prescribed) mobilities. Only
 * the entries in v corresponding to free mobilities are examined, and
 * only the entries in MinvV corresponding to free mobilities are written.
 *
 * Once the appropriate articulated body inertias are available, repeated
 * calls to this operator are very fast, with worst case around 80*n flops
 * when all mobilizers have 1 dof.
 *
 * Requires prior realization through Stage::Position.
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

/**
 * This operator calculates in O(N) time the product M*v where M is the 
 * system mass matrix and v is a supplied vector with one entry per 
 * mobility. If v is a set of mobility accelerations (generalized 
 * accelerations), then the result is a generalized force (f=M*a). 
 * Only the supplied vector is used, and M depends only on position 
 * states, so the result here is not affected by velocities in the State.
 * Constraints and prescribed motions are ignored.
 *
 * The current implementation requires about 120*n flops and does
 * not require realization of composite-body or articulated-body
 * inertias.
 *
 * Requires prior realization through Stage::Position.
 */
void calcMV(const State&, const Vector& v, Vector& MV) const;

/// This operator calculates the composite body inertias R given
/// a State realized to Position stage. Composite body inertias
/// are the spatial mass properties of the rigid body formed by
/// a particular body and all bodies outboard of that body if
/// all the outboard mobilizers were welded in their current
/// orientations.
void calcCompositeBodyInertias(const State&,
    Vector_<SpatialMat>& R) const;

/**
Given a complete set of generalized accelerations, this kinematic operator
calculates the resulting body accelerations, including velocity-dependent 
terms taken from the supplied State.
@pre \a state must already be realized to Velocity stage
@param[in] state
    The State from which position- and velocity- related terms are taken; 
    must already have been realized to Velocity stage.
@param[in] knownUDot
    A complete set of generalized accelerations. Must have the same length 
    as the number of mobilities, or if length zero the udots will be taken 
    as all zero in which case only velocity-dependent accelerations will be
    returned in \a A_GB.
@param[out] A_GB
    Spatial accelerations of all the body frames measured and expressed in
    the Ground frame, resulting from supplied generalized accelerations 
    \a knownUDot and velocity-dependent acceleration terms taken from 
    \a state. This will be resized if necessary to the number of bodies 
    <em>including</em> Ground so that the returned array may be indexed by 
    MobilizedBodyIndex with A_GB[0]==0 always. The angular acceleration
    vector for MobilizedBody i is A_GB[i][0]; linear acceleration of the
    body's origin is A_GB[i][1].
*/
void calcAccelerationFromUDot(const State&         state,
                              const Vector&        knownUDot,
                              Vector_<SpatialVec>& A_GB) const;

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
/// and a vector of length n. m is the number of active acceleration 
/// constraint equations, n is the number of mobilities.
/// This is an O(m+n) operation.
void calcGV(const State&,
    const Vector&   v,
    Vector&         Gv) const;

/// Returns G^T*v, the product of the nXm transpose of the acceleration 
/// constraint Jacobian G and a vector v of length m. m is the number of 
/// active acceleration constraint equations, n is the number of 
/// mobilities. If v is a set of constraint multipliers, then
/// f=G^T*v is the set of equivalent generalized forces they generate.
/// This is an O(m+n) operation.
void calcGtV(const State&,
    const Vector&   v,
    Vector&         GtV) const;


/// This operator explicitly calculates the n X n mass matrix M. Note that
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

/// Treating all constraints together, given a comprehensive set of 
/// multipliers lambda, generate the complete set of body and mobility 
/// forces applied by all the constraints; watch the sign -- normally 
/// constraint forces have opposite sign from applied forces. If you want 
/// to take Simbody-calculated multipliers and use them to generate forces
/// that look like applied forces, negate the multipliers before making 
/// this call.
/// 
/// State must be realized to Stage::Position to call this operator 
/// (although typically the multipliers are obtained by realizing to 
/// Stage::Acceleration).
void calcConstraintForcesFromMultipliers
  (const State& s, const Vector& multipliers,
   Vector_<SpatialVec>& bodyForcesInG,
   Vector&              mobilityForces) const;

/// Calculate the mobilizer reaction force generated by each MobilizedBody.  
/// This is the constraint force that would be required to make the system 
/// move in the same way if that MobilizedBody were converted to a Free body.
///
/// A simple way to think of the reaction force is to think of cutting the 
/// mobilizer, then imagine the force required to make the system move in 
/// the same manner as when the mobilizer was present. This is what the 
/// reaction forces accomplish. With that definition, mobility forces (as 
/// opposed to body forces) are \e included in the reactions. Some 
/// conventions do not include the mobility forces in the definition of a 
/// reaction force. We chose to include them since this preserves Newton’s 
/// 3rd law of equal and opposite reactions between bodies. Ours is the same 
/// convention as used in SD/FAST.
/// 
/// A mobilizer exerts equal and opposite reaction forces on the parent 
/// (inboard) and child (outboard) bodies. This method reports the force on 
/// the child body, as though it were applied at the origin of the outboard 
/// mobilizer frame M (fixed to the child), and expressed in the Ground frame.
///
/// @param[in] state        
///     A State compatible with this System that has already been realized 
///     to Stage::Acceleration.
/// @param[out] forcesAtMInG    
///     A Vector of spatial force vectors, indexed by MobilizedBodyIndex 
///     (beginning with 0 for Ground), giving the reaction moment and force
///     applied by each body's unique inboard mobilizer to that body. The
///     force is returned as though it were applied at the origin of the 
///     body's mobilizer frame M. The returned force is expressed in the
///     Ground frame. Note that applied mobility (generalized) forces are 
///     \e included in the returned reaction forces (see above).
///                           
void calcMobilizerReactionForces
   (const State&         state, 
    Vector_<SpatialVec>& forcesAtMInG) const;

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

/// Accounts for applied forces and inertial forces produced by non-
/// zero velocities in the State. Returns a set of mobility forces which
/// replace both the applied bodyForces and the inertial forces.
/// Requires prior realization through Stage::Dynamics. 
void calcTreeEquivalentMobilityForces(const State&, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const;



/// Must be in Stage::Position to calculate qdot = N(q)*u.
void calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const;

/// Must be in Stage::Velocity to calculate 
/// qdotdot = N(q)*udot + Ndot(q,u)*u.
void calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const;

/// Must be in Stage::Position to calculate out_q = N(q)*in_u (e.g., 
/// qdot=N*u) or out_u = ~N*in_q. Note that one of "in" and "out" is always
/// "q-like" while the other is "u-like", but which is which changes if the 
/// matrix is transposed. Note that the transposed operation here is the 
/// same as multiplying by N on the right, with the Vectors viewed as 
/// RowVectors instead. This is an O(n) operator since N is block diagonal.
void multiplyByN(const State& s, bool transpose, 
                 const Vector& in, Vector& out) const;

/// Must be in Stage::Position to calculate out_u = NInv(q)*in_q (e.g., 
/// u=NInv*qdot) or out_q = ~NInv*in_u. Note that one of "in" and "out" is
/// always "q-like" while the other is "u-like", but which is which changes
/// if the matrix is transposed. Note that the transposed operation here is
/// the same as multiplying by NInv on the right, with the Vectors viewed 
/// as RowVectors instead. This is an O(N) operator since NInv is block 
/// diagonal.
void multiplyByNInv(const State& s, bool transpose, 
                    const Vector& in, Vector& out) const;

/// Must be in Stage::Velocity to calculate out_q = NDot(q,u)*in_u
/// or out_u = ~NDot(q,u)*in_q. This is used, for example, as part of the 
/// conversion between udot and qdotdot. Note that one of "in" and "out" is
/// always "q-like" while the other is "u-like", but which is which changes
/// if the matrix is transposed. Note that the transposed operation here is
/// the same as multiplying by NDot on the right, with the Vectors viewed 
/// as RowVectors instead. This is an O(N) operator since NDot is block 
/// diagonal.
void multiplyByNDot(const State& s, bool transpose, 
                    const Vector& in, Vector& out) const;

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
/// method we should use to model their orientations. Choices are: 
/// quaternions (best for dynamics), or rotation angles (1-2-3 Euler 
/// sequence, good for optimization). TODO: (1) other Euler sequences, 
/// (2) allow settable zero rotation for Euler sequence, with convenient
/// way to say "this is zero".
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

    // POSITION STAGE realizations //

/// This method checks whether composite body inertias have already
/// been computed since the last change to a Position stage state 
/// variable (q) and if so returns immediately at little cost; otherwise,
/// it initiates computation of composite body inertias for all of
/// the mobilized bodies. These are not otherwise
/// computed unless specifically requested. You cannot call
/// this method unless the State has already been realized through
/// Position stage.
void realizeCompositeBodyInertias(const State&) const;

/// This method checks whether articulated body inertias have already
/// been computed since the last change to a Position stage state 
/// variable (q) and if so returns immediately at little cost; otherwise,
/// it initiates the relatively expensive computation of articulated 
/// body inertias for all of the mobilized bodies. These are not otherwise
/// computed until they are needed at Dynamics stage. You cannot call
/// this method unless the State has already been realized through
/// Position stage.
void realizeArticulatedBodyInertias(const State&) const;


    // POSITION STAGE responses //

/// Return the composite body inertia for a particular mobilized body. You
/// can call this any time after the State has been realized to Position
/// stage, however it will first trigger realization of all the composite body
/// inertias if they have not already been calculated. Ground is mobilized body 
/// zero; its composite body inertia has infinite mass and principle moments of
/// inertia, and zero center of mass.
/// @see realizeCompositeBodyInertias()
const SpatialMat& getCompositeBodyInertia(const State&, MobilizedBodyIndex) const;

/// Return the articulated body inertia for a particular mobilized body. You
/// can call this any time after the State has been realized to Position
/// stage, however it will first trigger expensive realization of all the articulated body
/// inertias if they have not already been calculated. Ground is mobilized body 
/// zero; its articulated body inertia is the same as its composite body inertia --
/// an ordinary Spatial Inertia but with infinite mass and principle moments of
/// inertia, and zero center of mass.
/// @see realizeArticulatedBodyInertias()
const SpatialMat& getArticulatedBodyInertia(const State&, MobilizedBodyIndex) const;

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

/// This solver transfers prescribed positions presQ or prescribed velocities
/// presU into the supplied State. When called with Stage::Position, the State
/// must have already been realized to Time stage; when called with 
/// Stage::Velocity, the State must already have been realized to Position
/// Stage. Returns true if it makes any State changes.
bool prescribe(State&, Stage) const;

/// This is a solver you can call after the State has been realized
/// to stage Position. It will project the q constraints
/// along the error norm so that getQConstraintNorm() <= consAccuracy, and will
/// project out the corresponding component of yErrest so that yErrest's q norm
/// is reduced. Returns true if it does anything at all to State or yErrest.
bool projectQConstraints(State& s, Real consAccuracy, const Vector& yWeights,
						 const Vector& ooTols, Vector& yErrest, System::ProjectOptions) const;

    // VELOCITY STAGE responses //

// Cross joint
const SpatialVec& getCoriolisAcceleration(const State&, MobilizedBodyIndex) const;

// Including parent
const SpatialVec& getTotalCoriolisAcceleration(const State&, MobilizedBodyIndex) const;

const SpatialVec& getGyroscopicForce(const State&, MobilizedBodyIndex) const;
const SpatialVec& getCentrifugalForces(const State&, MobilizedBodyIndex) const;

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
