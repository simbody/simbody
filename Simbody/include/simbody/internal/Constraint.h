#ifndef SimTK_SIMBODY_CONSTRAINT_H_
#define SimTK_SIMBODY_CONSTRAINT_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/** @file
This defines the base Constraint class and related classes, which are used to 
specify limitations on the mobility of the mobilized bodies in a 
SimbodyMatterSubsystem. **/


#include "SimTKmath.h"
#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;
class Constraint;
class ConstraintImpl;

// We only want the template instantiation to occur once. This symbol is 
// defined in the SimTK core compilation unit that defines the Constraint 
// class but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_CONSTRAINT
    extern template class PIMPLHandle<Constraint, ConstraintImpl, true>;
#endif

    ///////////////////////////
    // CONSTRAINT BASE CLASS //
    ///////////////////////////

/** This is the base class for all %Constraint classes, which is just a handle 
for the underlying hidden implementation. There is a set of built-in 
constraints and a generic "Custom" constraint (an abstract base class) from 
which advanced users may derive their own constraints. Each built-in constraint
type is a local subclass within %Constraint, and is also derived from 
%Constraint.

%Constraint is a PIMPL-style abstract base class, with concrete classes defined 
for each kind of constraint. **/
class SimTK_SIMBODY_EXPORT Constraint 
:   public PIMPLHandle<Constraint, ConstraintImpl, true> {
public:
/** Default constructor creates an empty %Constraint handle that can be used
to reference any %Constraint. **/
Constraint() { }
/** For internal use: construct a new %Constraint handle referencing a 
particular implementation object. **/
explicit Constraint(ConstraintImpl* r) : HandleBase(r) { }

/** Disable this %Constraint, effectively removing it from the system. This
is an Instance-stage change and affects the allocation of %Constraint-
related resources in the supplied State. **/
void disable(State&) const;

/** Enable this %Constraint, without necessarily satisfying it. This is an 
Instance-stage change and affects the allocation of %Constraint-related  
resources in the supplied State. Note that merely enabling a constraint does 
not ensure that the State's positions and velocities satisfy that constraint; 
initial satisfaction requires use of an appropriate project() solver.
@see SimTK::System::project() **/
void enable(State&) const;
/** Test whether this constraint is currently disabled in the supplied 
State. **/
bool isDisabled(const State&) const;
/** Test whether this %Constraint is disabled by default in which case it must 
be explicitly enabled before it will take effect.
@see setDisabledByDefault(), enable() **/
bool isDisabledByDefault() const;

/** Normally Constraints are enabled when defined and can be disabled later. If
you want to define this constraint but have it be off by default, use
this method.
@see isDisabledByDefault(), enable(), disable(), isDisabled() **/
void setDisabledByDefault(bool shouldBeDisabled);

/** This is an implicit conversion from Constraint to ConstraintIndex when 
needed. This will fail if the %Constraint is not contained in a subsystem. **/
operator ConstraintIndex() const {return getConstraintIndex();}

/** Get a const reference to the matter subsystem that contains this 
%Constraint. This will throw an exception if the %Constraint has not yet been
added to any subsystem; if you aren't sure use isInSubsystem() first to
check.  
@see updMatterSubsystem(), isInSubsystem(), getConstraintIndex() **/
const SimbodyMatterSubsystem& getMatterSubsystem()      const;

/** Assuming you have writable access to this %Constraint, get a writable 
reference to the containing matter subsystem. This will throw an exception if 
the %Constraint has not yet been added to any subsystem; if you aren't sure 
use isInSubsystem() first to check.  
@see getMatterSubsystem(), isInSubsystem(), getConstraintIndex() **/
SimbodyMatterSubsystem& updMatterSubsystem();

/** Get the ConstraintIndex that was assigned to this %Constraint when it was
added to the matter subsystem. This will throw an exception if the %Constraint 
has not yet been added to any subsystem; if you aren't sure use isInSubsystem() 
first to check. There is also an implicit conversion from %Constraint to
ConstraintIndex, so you don't normally need to call this directly.
@see getMatterSubsystem(), isInSubsystem() **/
ConstraintIndex getConstraintIndex() const;

/** Test whether this %Constraint is contained within a matter subsystem. **/
bool isInSubsystem() const;
/** Test whether the supplied MobilizedBody is in the same matter subsystem
as this %Constraint. Also returns false if either the %Constraint or the 
%MobilizedBody is not in any subsystem, or if neither is. **/
bool isInSameSubsystem(const MobilizedBody& mobod) const;

    // TOPOLOGY STAGE (i.e., post-construction) //

/** Return the number of unique bodies \e directly restricted by this 
constraint. Included are any bodies to which this %Constraint may apply a body 
force (i.e., torque or point force). The Ancestor body is not included unless
it was specified as a Constrained Body. This is the length of the bodyForces 
array for this %Constraint. **/
int getNumConstrainedBodies() const;

/** Return a const reference to the actual MobilizedBody corresponding to one
of the Constrained Bodies included in the count returned by 
getNumConstrainedBodies(). The index must be in the range 
0 <= \a consBodyIx < getNumConstrainedBodies(). **/
const MobilizedBody& getMobilizedBodyFromConstrainedBody
   (ConstrainedBodyIndex consBodyIx) const;

/** Return a const reference to the actual MobilizedBody which is serving as
the Ancestor body for the constrained bodies in this Constraint. This
will fail if there are no constrained bodies (i.e., if 
getNumConstrainedBodies()==0). **/
const MobilizedBody& getAncestorMobilizedBody() const;

/** Return the number of unique mobilizers \e directly restricted by this
%Constraint. Included are any mobilizers to which the %Constraint may
apply any mobility force. Like bodies, mobilizers are referenced using the 
MobilizedBody containing them. Note that all the mobilities of a Constrained
Mobilizer are included in the set of constrainable Qs or constrainable Us for 
this %Constraint even if not all of them are constrained. **/
int getNumConstrainedMobilizers() const;

/** Return a const reference to the actual MobilizedBody corresponding to one
of the Constrained Mobilizers included in the count returned by 
getNumConstrainedMobilizers(). The index must be in the range 
0 <= \a consMobilizerIx < getNumConstrainedMobilizers(). **/
const MobilizedBody& getMobilizedBodyFromConstrainedMobilizer
   (ConstrainedMobilizerIndex consMobilizerIx) const;

/** Return a subtree object indicating which parts of the multibody tree
are potentially affected by this %Constraint. **/
const SimbodyMatterSubtree& getSubtree() const;

    // MODEL STAGE //
// nothing in base class currently

    // INSTANCE STAGE //

/** Return the number of constrainable generalized coordinates q associated 
with a particular constrained mobilizer. This is just the number of generalized
coordinates for that mobilizer; any or all of them may actually be 
unconstrained. **/
int getNumConstrainedQ(const State&, ConstrainedMobilizerIndex) const;
/** Return the number of constrainable mobilities u associated with a 
particular constrained mobilizer. This is just the number of generalized speeds
for that mobilizer; any or all of them may actually be unconstrained. The 
number of constrainable udots is the same. **/
int getNumConstrainedU(const State&, ConstrainedMobilizerIndex) const;

/** Return the index into the constrained mobilities u array corresponding to a
particular mobility of the indicated ConstrainedMobilizer. Don't confuse this 
with the set of \e participating mobilities which also includes all mobilities
on each branch between the ancestor and a constrained body. The \e constrained
mobilities are just those belonging to the mobilizers which are directly 
constrained. **/
ConstrainedUIndex getConstrainedUIndex
    (const State&, ConstrainedMobilizerIndex, MobilizerUIndex which) const;
/** Return the index into the constrained coordinates q array corresponding to 
a particular coordinate of the indicated ConstrainedMobilizer. Don't confuse 
this with the set of \e participating coordinates which also includes all 
coordinates on each branch between the ancestor and a constrained body. The 
\e constrained coordinates are just those belonging to the mobilizers which are
directly constrained. **/
ConstrainedQIndex getConstrainedQIndex
    (const State&, ConstrainedMobilizerIndex, MobilizerQIndex which) const;

/** Return the sum of the number of coordinates q associated with each of
the constrained mobilizers. **/
int getNumConstrainedQ(const State&) const;

/** Return the sum of the number of mobilities u associated with each of the 
constrained mobilizers. These are the only mobilities to which the constraint 
may directly apply a force, so this is also the dimension of the mobilityForces 
array. **/
int getNumConstrainedU(const State&) const;

/** Map one of this %Constraint's constrained q's to the corresponding index 
within the matter subsystem's whole q vector. **/
QIndex getQIndexOfConstrainedQ(const State&      state,
                               ConstrainedQIndex consQIndex) const;
/** Map one of this %Constraint's constrained U's (or mobilities) to the 
corresponding index within the matter subsystem's whole u vector. **/
UIndex getUIndexOfConstrainedU(const State&      state,
                               ConstrainedUIndex consUIndex) const;

/** Find out how many holonomic (position), nonholonomic (velocity), and
acceleration-only constraint equations are currently being generated 
by this Constraint. 

@param[in]      state
    The State from which the current status of this %Constraint is obtained.
    Must have been realized through Instance stage.
@param[out]     mp  The number of holonomic constraint equations.
@param[out]     mv  The number of nonholonomic constraint equations.
@param[out]     ma  The number of acceleration-only constraint equations.

Note that the counts here do not include the derivatives of the higher-order
constraint equations, just the number at the level they were defined. **/
void getNumConstraintEquationsInUse(const State& state, 
                                    int& mp, int& mv, int& ma) const;

/** Return the start of the blocks of multipliers (or acceleration errors)
assigned to this %Constraint. Separate blocks are allocated for holonomic
(position), nonholonomic (velocity), and acceleration-only constraint
equations. The size of each block is given by getNumConstraintEquationsInUse();
this %Constraint's multipliers are assigned contiguously within each block.
If any size is zero, the corresponding index is returned invalid. 
    
@param[in]      state
    The State from which the current status of this %Constraint is obtained.
    Must have been realized through Instance stage.
@param[out]     px0  
    The index of the first slot for the second time derivatives of position 
    (holonomic) constraint equations.
@param[out]     vx0 
    The index of the first slot for the time derivatives of velocity 
    (nonholonomic) constraint equations.
@param[out]     ax0
    The index of the first slot for the acceleration-only constraint equations. 

For position and velocity constraints, the multiplier slots correspond to
the time derivatives of these constraints that are used to create acceleration
constraints from them. **/
void getIndexOfMultipliersInUse(const State& state,
                                MultiplierIndex& px0, 
                                MultiplierIndex& vx0, 
                                MultiplierIndex& ax0) const;

/** Set the part of a complete constraint-space vector that belongs to this
constraint. The full vector has dimension m=mp+mv+ma, that is, one entry
per acceleration-level constraint equation.
@param[in]      state
    The State from which the current status of this %Constraint is obtained.
    Must have been realized through Instance stage.
@param[in]      myPart
    The constraint-space scalars for this %Constraint in the order
    position, velocity, acceleration if this %Constraint produces constraint
    equations of different types. The number of entries must match the number
    of constraint equations generated by this %Constraint.
@param[in,out]  constraintSpace
    An array of full constraint space dimension m. If it has length zero on
    entry we'll resize it to m and initialize it to zero, otherwise the size 
    must be exactly m and we'll only modify the slots that belong to this
    %Constraint.

Note that we're writing only to the output argument; this method does not
calculate or modify anything else.

@see getMyPartFromConstraintSpaceVector(), getIndexOfMultipliersInUse() **/
void setMyPartInConstraintSpaceVector(const State& state,
                                      const Vector& myPart,
                                      Vector& constraintSpace) const;

/** Get the part of a complete constraint-space vector that belongs to this
constraint. The full vector has dimension m=mp+mv+ma, that is, one entry
per acceleration-level constraint equation.
@param[in]      state
    The State from which the current status of this %Constraint is obtained.
    Must have been realized through Instance stage.
@param[in]      constraintSpace
    An array of full constraint space dimension m. We will only examine the
    entries belonging to this %Constraint.
@param[out]     myPart
    The constraint-space scalars for this %Constraint in the order
    position, velocity, acceleration if this %Constraint produces constraint
    equations of different types. The number of entries will match the number
    of constraint equations generated by this %Constraint and the argument
    will be resized if necessary.

@see setMyPartInConstraintSpaceVector(), getIndexOfMultipliersInUse() **/
void getMyPartFromConstraintSpaceVector(const State& state,
                                        const Vector& constraintSpace,
                                        Vector& myPart) const;

    // POSITION STAGE //
/** Get a Vector containing the position errors. Many subclasses provide 
their own methods for getting this information in a more specific form. **/
Vector getPositionErrorsAsVector(const State&) const;   // mp of these
Vector calcPositionErrorFromQ(const State&, const Vector& q) const;

// Matrix P = partial(perr_dot)/partial(u). (just the holonomic constraints)
Matrix calcPositionConstraintMatrixP(const State&) const; // mp X nu
Matrix calcPositionConstraintMatrixPt(const State&) const; // nu X mp

// Matrix PNInv = partial(perr)/partial(q) = P*N^-1
Matrix calcPositionConstraintMatrixPNInv(const State&) const; // mp X nq

/** This operator calculates this constraint's body and mobility forces given 
the complete set of multipliers lambda for this Constraint. We expect that 
lambda has been packed to include multipliers associated with the second 
time derivatives of the position (holonomic) constraints, the first time
derivatives of the velocity (nonholonomic) constraints, and the 
acceleration-only constraints, in that order.

The state must be realized already to Stage::Velocity. Returned body forces 
correspond only to the <em>constrained bodies</em> and the mobility forces 
correspond only to the <em>constrained mobilities</em>; they must be unpacked 
by the caller into the actual system mobilized bodies and actual system 
mobilities. Note that the body forces are in the ancestor body frame A, not 
necessarily the Ground frame G, and that they are opposite in sign from
applied forces. If you want to calculate forces you can treat as applied
forces, negate \a lambda before the call. **/
void calcConstraintForcesFromMultipliers(const State&,
    const Vector&        lambda,                // mp+mv+ma of these
    Vector_<SpatialVec>& bodyForcesInA,         // numConstrainedBodies
    Vector&              mobilityForces) const; // numConstrainedU

    // VELOCITY STAGE //
/** Get a Vector containing the velocity errors. Many subclasses provide 
their own methods for getting this information in a more specific form. **/
Vector getVelocityErrorsAsVector(const State&) const;   // mp+mv of these
Vector calcVelocityErrorFromU(const State&,     // mp+mv of these
                              const Vector& u) const;   // numParticipatingU u's

// Matrix V = partial(verr)/partial(u) for just the non-holonomic 
// constraints.
Matrix calcVelocityConstraintMatrixV(const State&) const;  // mv X nu
Matrix calcVelocityConstraintMatrixVt(const State&) const; // nu X mv

    // DYNAMICS STAGE //
// nothing in base class currently

    // ACCELERATION STAGE //
/** Get a Vector containing the acceleration errors. Many subclasses 
provide their own methods for getting this information in a more 
specific form. **/
Vector getAccelerationErrorsAsVector(const State&) const;   // mp+mv+ma of these
Vector calcAccelerationErrorFromUDot(const State&,  // mp+mv+ma of these
                                     const Vector& udot) const; // numParticipatingU udot's

/** Get a Vector containing the Lagrange multipliers. Many subclasses 
provide their own methods for getting this information in a more 
specific form. **/
Vector getMultipliersAsVector(const State&) const;  // mp+mv+ma of these   

/** Given a State realized through Acceleration stage, return the forces
that were applied to the system by this %Constraint, with body forces
expressed in Ground. Note that the sign convention for constraint forces
is opposite that of applied forces, because constraints appear on the left
hand side in Simbody's equations of motion, while applied forces are on 
the right hand side.

These forces are the same as what you would get if you get the multipliers 
from this \a state using getMultipliersAsVector(), call 
calcConstraintForcesFromMultipliers(), and re-express the constrained body 
forces in the Ground frame. However, the ones returned here are already 
calculated so require only copying out of the \a state cache. **/
void getConstraintForcesAsVectors
   (const State&         state,
    Vector_<SpatialVec>& bodyForcesInG, // numConstrainedBodies
    Vector&              mobilityForces) const; // numConstrainedU

/** For convenience, returns constrained body forces as the function return. 
@see getConstraintForcesAsVectors() **/
Vector_<SpatialVec> getConstrainedBodyForcesAsVector(const State& state) const {
    Vector_<SpatialVec> bodyForcesInG;
    Vector              mobilityForces;
    getConstraintForcesAsVectors(state,bodyForcesInG,mobilityForces);
    return bodyForcesInG;
}
/** For convenience, returns constrained mobility forces as the function
return. 
@see getConstraintForcesAsVectors() **/
Vector getConstrainedMobilityForcesAsVector(const State& state) const {
    Vector_<SpatialVec> bodyForcesInG;
    Vector              mobilityForces;
    getConstraintForcesAsVectors(state,bodyForcesInG,mobilityForces);
    return mobilityForces;
}

/** Calculate the power being applied by this %Constraint to the system.
The \a state must be realized through Acceleration stage so that the 
applied constraint forces are known. Then power is calculated as the
dot product of the \e applied body spatial forces and body spatial velocities, 
plus the dot product of the \e applied mobility forces and corresponding 
mobilities (generalized speeds) u. I emphasized \e applied here because the
sign convention is opposite for constraint forces, so the power calculation
requires negating the constraint forces.

For any non-working %Constraint, power should always be within machine
precision of zero. This is a very useful test when debugging new Constraints.
For working Constraints, you can calculate work done as the time integral of 
the power. Then if you embed the %Constraint in an otherwise conservative
system, the sum of system potential and kinetic energy, minus the work done
by this constraint, should be constant to within integration accuracy.
Power and work here are signed quantities with positive sign meaning that
the %Constraint is adding energy to the system and negative meaning it is 
removing energy from the system. 

Computational cost here is low because the forces and velocities are already
known. Only the dot product need be computed, at a cost of about 
11 ncb + 2 ncu flops, where ncb is the number of constrained bodies and ncu
is the number of constrained mobilities for this %Constraint. **/
Real calcPower(const State& state) const;

// Matrix A = partial(aerr)/partial(udot) for just the acceleration-only 
// constraints.
Matrix calcAccelerationConstraintMatrixA(const State&) const;  // ma X nu
Matrix calcAccelerationConstraintMatrixAt(const State&) const; // nu X ma

/** (Advanced) Mark this constraint as one that is only conditionally active.
The conditions under which it is active must be evaluated elsewhere. This should
be set immediately after construction of the Constraint and invalidates
Stage::Topology. **/
void setIsConditional(bool isConditional);
/** (Advanced) Get the value of the isConditional flag. **/
bool isConditional() const;
                  
//------------------------------------------------------------------------------
// These are the built-in Constraint types, and some synonyms. Each built in 
// Constraint type is declared in its own header file using naming convention 
// Constraint_Rod.h, for example. All the built-in headers are collected in 
// Constraint_BuiltIns.h; you should include new ones there also.
class Rod;  
typedef Rod  ConstantDistance;  ///< Synonym for Rod constraint.

class Ball; 
typedef Ball CoincidentPoints;  ///< Synonym for Ball constraint.
typedef Ball Spherical;         ///< Synonym for Ball constraint.

class Weld; 
typedef Weld CoincidentFrames;

class PointInPlane;  // translations perpendicular to plane normal only
class PointOnLine;   // translations along a line only
class ConstantAngle; // prevent rotation about common normal of two vectors
class ConstantOrientation; // allows any translation but no rotation
class NoSlip1D; // same velocity at a point along a direction
class ConstantCoordinate; // prescribe generalized coordinate value
class ConstantSpeed; // prescribe generalized speed value
class ConstantAcceleration; // prescribe generalized acceleration value
class Custom;
class CoordinateCoupler;
class SpeedCoupler;
class PrescribedMotion;
class PointOnPlaneContact; 
class SphereOnPlaneContact; // ball in contact with plane (sliding or rolling)
class SphereOnSphereContact; // ball in contact with ball (sliding or rolling)
class LineOnLineContact;    // edge/edge contact

// Internal use only.
class RodImpl;
class BallImpl;
class WeldImpl;
class PointInPlaneImpl;
class PointOnLineImpl;
class ConstantAngleImpl;
class ConstantOrientationImpl;
class NoSlip1DImpl;
class ConstantCoordinateImpl;
class ConstantSpeedImpl;
class ConstantAccelerationImpl;
class CustomImpl;
class CoordinateCouplerImpl;
class SpeedCouplerImpl;
class PrescribedMotionImpl;
class PointOnPlaneContactImpl; 
class SphereOnPlaneContactImpl;
class SphereOnSphereContactImpl;
class LineOnLineContactImpl;

};

    //////////////////////////////
    // POINT ON LINE CONSTRAINT //
    //////////////////////////////

/**
 *  Two constraint equations. This constraint enforces that a point fixed to
 *  one body (the "follower body") must travel along a line fixed on another body (the
 *  "line body"). The constraint is enforced by an internal (non-working)
 *  scalar force acting at the spatial location of the follower point, directed in the
 *  plane for which the line is a normal, and equal and opposite on the two bodies.
 * 
 *  The assembly condition is the same as the run-time constraint: the point
 *  has to be moved onto the line.
 */
class SimTK_SIMBODY_EXPORT Constraint::PointOnLine : public Constraint  {
public:
    // no default constructor
    PointOnLine(MobilizedBody& lineBody_B, const UnitVec3& defaultLineDirection_B, const Vec3& defaultPointOnLine_B,
                MobilizedBody& followerBody_F, const Vec3& defaultFollowerPoint_F);
    
    /** Default constructor creates an empty handle. **/
    PointOnLine() {}

    // These affect only generated decorative geometry for visualization;
    // the line is really infinite in extent and the
    // point is really of zero radius.
    PointOnLine& setLineDisplayHalfLength(Real);
    PointOnLine& setPointDisplayRadius(Real);
    Real getLineDisplayHalfLength() const;
    Real getPointDisplayRadius() const;

    // Defaults for Instance variables.
    PointOnLine& setDefaultLineDirection(const UnitVec3&);
    PointOnLine& setDefaultPointOnLine(const Vec3&);
    PointOnLine& setDefaultFollowerPoint(const Vec3&);

    // Stage::Topology
    MobilizedBodyIndex getLineMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const UnitVec3& getDefaultLineDirection() const;
    const Vec3&     getDefaultPointOnLine() const;
    const Vec3&     getDefaultFollowerPoint() const;

    // Stage::Instance
    const UnitVec3& getLineDirection(const State&) const;
    const Vec3&     getPointOnLine(const State&) const;
    const Vec3&     getFollowerPoint(const State&) const;

    // Stage::Position, Velocity
    Vec2 getPositionErrors(const State&) const;
    Vec2 getVelocityErrors(const State&) const;

    // Stage::Acceleration
    Vec2 getAccelerationErrors(const State&) const;
    Vec2 getMultipliers(const State&) const;
    const Vec2& getForceOnFollowerPoint(const State&) const; // in normal direction
    
    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (PointOnLine, PointOnLineImpl, Constraint);
    /** @endcond **/
};

    ///////////////////////////////
    // CONSTANT ANGLE CONSTRAINT //
    ///////////////////////////////

/**
 * This constraint consists of a single constraint equation that enforces that
 * a unit vector v1 fixed to one body (the "base body") must maintain a fixed 
 * angle theta with respect to a unit vector v2 fixed on the other body (the 
 * "follower body"). This can be done with a single constraint equation as long 
 * as theta is sufficiently far away from 0 and +/-Pi (180 degrees), with the 
 * numerically best performance at theta=Pi/2 (90 degrees).
 *
 * @warning
 * Do not use this constraint to \e align the vectors, that is for angles near 
 * 0 or +/- Pi; performance will noticeably degrade within a few degrees of 
 * these limits and numerical integration will eventually fail at the limits.
 * 
 * If you want to enforce that two axes are aligned with one another (that 
 * is, the angle between them is 0 or +/-Pi), that takes \e two constraint 
 * equations since the only remaining rotation is about the common axis. (That 
 * is, two rotational degrees of freedom are removed; that can't be done with 
 * one constraint equation -- the situation is analogous to the inability of
 * a Rod (distance) constraint to keep two points at 0 distance.) Instead,
 * you can use two ConstantAngle constraints on pairs of vectors perpendicular 
 * to the aligned ones, so that each ConstantAngle is set to the optimal 90 degrees.
 * 
 * This constraint is enforced by an internal scalar torque applied equal and
 * opposite on each body, about the mutual perpendicular to the two vectors.
 * 
 * The assembly condition is the same as the run-time constraint: the 
 * bodies must be rotated until the vectors have the right angle between them.
 */
class SimTK_SIMBODY_EXPORT Constraint::ConstantAngle : public Constraint {
public:
    // no default constructor
    ConstantAngle(MobilizedBody& baseBody_B,     const UnitVec3& defaultAxis_B,
                  MobilizedBody& followerBody_F, const UnitVec3& defaultAxis_F, 
                  Real angle = Pi/2);
    
    /** Default constructor creates an empty handle. **/
    ConstantAngle() {}

    // These affect only generated decorative geometry for visualization.
    ConstantAngle& setAxisDisplayLength(Real);
    ConstantAngle& setAxisDisplayWidth(Real);
    Real getAxisDisplayLength() const;
    Real getAxisDisplayWidth() const;

    // Defaults for Instance variables.
    ConstantAngle& setDefaultBaseAxis(const UnitVec3&);
    ConstantAngle& setDefaultFollowerAxis(const UnitVec3&);
    ConstantAngle& setDefaultAngle(Real);

    // Stage::Topology
    MobilizedBodyIndex getBaseMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const UnitVec3& getDefaultBaseAxis() const;
    const UnitVec3& getDefaultFollowerAxis() const;
    Real getDefaultAngle() const;

    // Stage::Instance
    const UnitVec3& getBaseAxis(const State&) const;
    const UnitVec3& getFollowerAxis(const State&) const;
    Real getAngle(const State&) const;

    // Stage::Position, Velocity
    Real getPositionError(const State&) const;
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getTorqueOnFollowerBody(const State&) const; // about f X b
    
    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ConstantAngle, ConstantAngleImpl, Constraint);
    /** @endcond **/
};


    /////////////////////////////////////
    // CONSTANT ORIENTATION CONSTRAINT //
    /////////////////////////////////////

/**
 *  Three constraint equations. This constraint enforces that a reference frame
 *  fixed to one body (the "follower body") must have the same orientation as another
 *  reference frame fixed on another body (the "base body"). That is, we have three
 *  constraint equations that collectively prohibit any relative rotation between
 *  the base and follower. The run time equations we use are just three "constant angle"
 *  constraints enforcing perpendicularity between follower's x,y,z axes with the base
 *  y,z,x axes respectively.
 * 
 *  This constraint is enforced by an internal (non-working) torque vector applied equal and
 *  opposite on each body.
 * 
 *  TODO: The assembly condition is not the same as the run-time constraint, because the
 *  perpendicularity conditions can be satisfied with antiparallel axes. For assembly
 *  we must have additional (redundant) constraints requiring parallel axes.
 */
class SimTK_SIMBODY_EXPORT Constraint::ConstantOrientation : public Constraint
{
public:
    // no default constructor
    ConstantOrientation(MobilizedBody& baseBody_B,     const Rotation& defaultRB,
                        MobilizedBody& followerBody_F, const Rotation& defaultRF); 
    
    /** Default constructor creates an empty handle. **/
    ConstantOrientation() {}

    //TODO: default visualization geometry?

    // Defaults for Instance variables.
    ConstantOrientation& setDefaultBaseRotation(const Rotation&);
    ConstantOrientation& setDefaultFollowerRotation(const Rotation&);

    // Stage::Topology
    MobilizedBodyIndex getBaseMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const Rotation& getDefaultBaseRotation() const;
    const Rotation& getDefaultFollowerRotation() const;

    // Stage::Instance
    const Rotation& getBaseRotation(const State&) const;
    const Rotation& getFollowerRotation(const State&) const;

    // Stage::Position, Velocity
    Vec3 getPositionErrors(const State&) const;
    Vec3 getVelocityErrors(const State&) const;

    // Stage::Acceleration
    Vec3 getAccelerationErrors(const State&) const;
    Vec3 getMultipliers(const State&) const;
    Vec3 getTorqueOnFollowerBody(const State&) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (ConstantOrientation, ConstantOrientationImpl, Constraint);
    /** @endcond **/
};

    ///////////////////////////
    // NO SLIP 1D CONSTRAINT //
    ///////////////////////////

/**
 * One non-holonomic constraint equation. There is a contact point P and a no-slip 
 * direction n fixed in a case body C. There are two moving bodies B0 and B1. The 
 * material point of B0 and the material point of B1 which are each coincident 
 * with the contact point P must have identical velocities in C, along the direction n.
 * This can be used to implement simple rolling contact between disks, such as occurs
 * in gear trains.
 * 
 * The assembly condition is the same as the run-time constraint: the velocities must
 * be made to match.
 */
class SimTK_SIMBODY_EXPORT Constraint::NoSlip1D : public Constraint {
public:
    // no default constructor

    /** Define the up to three bodies involved in this constraint: the two
    "moving" bodies and a Case body, and a default contact point and no-slip
    direction in the Case body frame C. (If you are modeling gears then the
    Case is the gearbox.) The case serves to define the 
    contact geometry but no forces are applied to it. It is OK for the Case
    body to be the same body as one of the moving bodies. **/
    NoSlip1D(MobilizedBody& caseBodyC, const Vec3& P_C, const UnitVec3& n_C,
             MobilizedBody& movingBody0, MobilizedBody& movingBody1);
    
    /** Default constructor creates an empty handle. **/
    NoSlip1D() {}

    /** Change the contact point at which this %Constraint acts.
    Provide the station location in the Case body local frame.
    This overrides the default point that was supplied on construction. This
    is an Instance-stage change. **/
    void setContactPoint(State& state, const Vec3& point_C) const;
    /** Change the no-slip direction along which this %Constraint acts.
    Provide the direction unit vector in the Case body local frame.
    This overrides the default direction that was supplied on construction. This
    is an Instance-stage change. **/
    void setDirection(State& state, const UnitVec3& direction_C) const;

    /** Return from the given \a state the contact point, in the Case body 
    frame. **/
    const Vec3& getContactPoint(const State& state) const;
    /** Return from the given \a state the no-slip direction, in  the Case 
    body frame. **/
    const UnitVec3& getDirection(const State& state) const;

    // These affect only generated decorative geometry for visualization;
    // the plane is really infinite in extent with zero depth and the
    // point is really of zero radius.

    /** For visualization only, set the length of the line used to show the
    no-slip direction. **/
    NoSlip1D& setDirectionDisplayLength(Real);
    /** For visualization only, set the radius of the sphere used to show
    the contact point location. **/
    NoSlip1D& setPointDisplayRadius(Real);
    /** Return the current value of the visualization line length for the 
    no-slip direction. **/
    Real getDirectionDisplayLength() const;
    /** Return the current value of the radius for visualization of the
    contact point. **/
    Real getPointDisplayRadius() const;

    // Defaults for Instance variables.

    /** Change the default contact point; this is the initial value for
    for the actual contact point and is a topological change. **/
    NoSlip1D& setDefaultContactPoint(const Vec3&); 
    /** Change the default no-slip direction; this is the initial value for
    for the actual direction and is a topological change. **/
    NoSlip1D& setDefaultDirection(const UnitVec3&);


    // Stage::Topology

    /** Get the mobilized body index of the Case body that was set during
    construction. **/
    MobilizedBodyIndex getCaseMobilizedBodyIndex() const;
    /** Get the mobilized body index of moving body 0 or moving body 1 that 
    was set during construction. Set \a which to 0 or 1 accordingly. **/
    MobilizedBodyIndex getMovingBodyMobilizedBodyIndex(int which) const;

    /** Obtain the default value for the no-slip direction, expressed in the
    Case body frame. **/
    const UnitVec3& getDefaultDirection() const;
    /** Obtain the default value for the contact point, in the Case body
    frame. **/
    const Vec3&     getDefaultContactPoint() const;


    // Stage::Position, Velocity
        // no position error

    /** Get the velocity error for this constraint equation, using configuration
    and velocity information from the given \a state, which must already have
    been realized through Velocity stage. **/
    Real getVelocityError(const State& state) const;

    // Stage::Acceleration

    /** Get the acceleration error for this constraint equation, using 
    configuration, velocity, and acceleration information from the given 
    \a state, which must already have been realized through Acceleration 
    stage. **/
    Real getAccelerationError(const State&) const;

    /** Get the Lagrange multiplier for this constraint equation, using 
    configuration, velocity, and acceleration information from the given 
    \a state, which must already have been realized through Acceleration 
    stage. While this is linearly related to the constraint force it may have 
    arbitrary sign and scaling; if you want an actual force use 
    getForceAtContactPoint() instead. **/
    Real getMultiplier(const State&) const;

    /** Determine the constraint force currently being generated by this 
    constraint. The force is as applied to the second moving body, that is,
    moving body 1, and is applied along the no-slip direction vector. **/
    Real getForceAtContactPoint(const State&) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(NoSlip1D, NoSlip1DImpl, Constraint);
    /** @endcond **/
};

    /////////////////////////
    // CONSTANT COORDINATE //
    /////////////////////////

/** Constrain a single mobilizer coordinate q to have a particular value.

Generates one position-level (holonomic) constraint equation. Some generalized 
coordinate q is required to remain at a particular position value p. This may
be an angle, a length, or some other unit depending on how the mobilizer is
defined.

Consider using the lock() or lockAt() feature of mobilizers (see MobilizedBody 
description) instead of this constraint; if applicable, locking is more 
efficient since it does not require adding a constraint equation to the system.
 
The assembly condition is the same as the run-time constraint: q must be set
to position p. 
@see MobilizedBody::lock(), MobilizedBody::lockAt() **/
class SimTK_SIMBODY_EXPORT Constraint::ConstantCoordinate : public Constraint {
public:
    /** Construct a constant coordinate constraint on a particular generalized
    coordinate q of the given mobilizer. Provide a default position value to 
    which the q should be locked; you can change it later via setPosition(). **/
    ConstantCoordinate(MobilizedBody& mobilizer, MobilizerQIndex whichQ, 
                       Real defaultPosition);

    /** Construct a constant coordinate constraint on the generalized
    coordinate q of the given mobilizer, assuming there is only one 
    coordinate. (Constrains the first coordinate if there are several.) Provide 
    a default position value to which the q should be locked; you can change it 
    later via setPosition(). **/
    ConstantCoordinate(MobilizedBody& mobilizer, Real defaultPosition); 
    
    /** Default constructor creates an empty handle you can use to reference
    any existing %ConstantCoordinate Constraint. **/
    ConstantCoordinate() {}

    /** Return the index of the mobilized body to which this constant coordinate
    constraint is being applied (to \e one of its coordinates). This is set on
    construction of the %ConstantCoordinate constraint. **/
    MobilizedBodyIndex getMobilizedBodyIndex() const;

    /** Return the particular coordinate whose position is controlled by
    this %ConstantCoordinate constraint. This is set on construction. **/
    MobilizerQIndex    getWhichQ() const;

    /** Return the default value for the position to be enforced. This is set on
    construction or via setDefaultPosition(). This is used to initialize the 
    position when a default State is created, but it can be overriden by 
    changing the value in the State using setPosition(). **/
    Real getDefaultPosition() const;

    /** Change the default value for the position to be enforced by this 
    constraint. This is a topological change, meaning you'll have to call
    realizeTopology() on the containing System and obtain a new State before
    you can use it. If you just want to make a runtime change in the State,
    see setPosition(). **/
    ConstantCoordinate& setDefaultPosition(Real position);

    /** Override the default position with this one whose value is stored in the
    given State. This invalidates the Position stage in the state. Don't 
    confuse this with setDefaultPosition() -- the value set here overrides that
    one. **/
    void setPosition(State& state, Real position) const;

    /** Get the current value of the position set point from the indicated 
    State. This is the value currently in effect, either from the default or 
    from a previous call to setPosition(). **/
    Real getPosition(const State& state) const;

    /** Return the amount by which the given State fails to satisfy this
    %ConstantCoordinate constraint. This is a signed value, q-p where p is
    the currently effective desired position as returned by getPosition()
    on this same \a state. The \a state must already be realized 
    through Stage::Position. **/
    Real getPositionError(const State& state) const;

    /** Return the amount by which the given State fails to satisfy the time
    derivative of this %ConstantCoordinate constraint, which should be zero. 
    This is a signed value equal to the current value of qdot (=d/dt q). The 
    \a state must already be realized through Stage::Velocity. **/
    Real getVelocityError(const State& state) const;

    /** Return the amount by which the accelerations in the given State fail
    to satify the second time derivative of this constraint, which should be 
    zero. This is a signed value equal to the current value of qdotdot
    (=d^2/dt^2 q). The \a state must already be realized through 
    Stage::Acceleration. **/
    Real getAccelerationError(const State& state) const;

    /** Get the value of the Lagrange multiplier generated to satisfy this
    constraint. For a %ConstantCoordinate constraint, the multiplier has the
    same magnitude as the generalized force although by convention constraint 
    multipliers have the opposite sign from applied forces. The \a state must 
    already be realized through Stage::Acceleration.**/
    Real getMultiplier(const State& state) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (ConstantCoordinate, ConstantCoordinateImpl, Constraint);
    /** @endcond **/
};

    ////////////////////
    // CONSTANT SPEED //
    ////////////////////

/** Constrain a single mobility to have a particular speed.

One non-holonomic constraint equation. Some mobility u is required to be at a
particular value s.

Consider using the lock() or lockAt() feature of mobilizers (see MobilizedBody 
description) instead of this constraint; if applicable, locking is more 
efficient since it does not require adding a constraint equation to the system.
 
The assembly condition is the same as the run-time constraint: u must be set
to s. 
@see MobilizedBody::lock(), MobilizedBody::lockAt() **/
class SimTK_SIMBODY_EXPORT Constraint::ConstantSpeed : public Constraint {
public:
    /** Construct a constant speed constraint on a particular mobility
    of the given mobilizer. **/
    ConstantSpeed(MobilizedBody& mobilizer, MobilizerUIndex whichU, 
                  Real defaultSpeed);
    /** Construct a constant speed constraint on the mobility
    of the given mobilizer, assuming there is only one mobility. **/
    ConstantSpeed(MobilizedBody& mobilizer, Real defaultSpeed); 
    
    /** Default constructor creates an empty handle you can use to reference
    any existing %ConstantSpeed Constraint. **/
    ConstantSpeed() {}

    /** Return the index of the mobilized body to which this constant speed
    constraint is being applied (to \e one of its mobilities). This is set on
    construction of the %ConstantSpeed constraint. **/
    MobilizedBodyIndex getMobilizedBodyIndex() const;

    /** Return the particular mobility whose generalized speed is controlled by
    this %ConstantSpeed constraint. This is set on construction. **/
    MobilizerUIndex    getWhichU() const;

    /** Return the default value for the speed to be enforced. This is set on
    construction or via setDefaultSpeed(). This is used to initialize the speed
    when a default State is created, but it can be overriden by changing the
    value in the State using setSpeed(). **/
    Real               getDefaultSpeed() const;

    /** Change the default value for the speed to be enforced by this 
    constraint. This is a topological change, meaning you'll have to call
    realizeTopology() on the containing System and obtain a new State before
    you can use it. If you just want to make a runtime change in the State,
    see setSpeed(). **/
    ConstantSpeed&     setDefaultSpeed(Real speed);

    /** Override the default speed with this one whose value is stored in the
    given State. This invalidates the Velocity stage in the state. Don't 
    confuse this with setDefaultSpeed() -- the value set here overrides that
    one. **/
    void setSpeed(State& state, Real speed) const;

    /** Get the current value of the speed set point from the indicated State.
    This is the value currently in effect, either from the default or from a
    previous call to setSpeed(). **/
    Real getSpeed(const State& state) const;

    // no position error

    /** Return the amount by which the given State fails to satisfy this
    %ConstantSpeed constraint. The \a state must already be realized through 
    Stage::Velocity. **/
    Real getVelocityError(const State& state) const;

    /** Return the amount by which the accelerations in the given State fail
    to satify the time derivative of this constraint (which must be zero). 
    The \a state must already be realized through Stage::Acceleration. **/
    Real getAccelerationError(const State& state) const;

    /** Get the value of the Lagrange multiplier generated to satisfy this
    constraint. For a %ConstantSpeed constraint, that is the same as the
    generalized force although by convention constraint multipliers have the
    opposite sign from applied forces. The \a state must already be realized 
    through Stage::Acceleration.**/
    Real getMultiplier(const State& state) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (ConstantSpeed, ConstantSpeedImpl, Constraint);
    /** @endcond **/
};

    ///////////////////////////
    // CONSTANT ACCELERATION //
    ///////////////////////////

/** Constrain a single mobility to have a particular acceleration.

One acceleration-only constraint equation. Some generalized acceleration
udot is required to be at a particular value a.

Consider using the lock() feature of mobilizers (see MobilizedBody description)
instead of this constraint; if applicable, locking is more efficient since it 
does not require adding a constraint equation to the system.

There is no assembly condition because this does not involve state
variables q or u, just u's time derivative udot. 
@see MobilizedBody::lock() **/
class SimTK_SIMBODY_EXPORT Constraint::ConstantAcceleration : public Constraint
{
public:
    /** Construct a constant acceleration constraint on a particular mobility
    of the given mobilizer. **/
    ConstantAcceleration(MobilizedBody& mobilizer, MobilizerUIndex whichU, 
                         Real defaultAcceleration);

    /** Construct a constant acceleration constraint on the mobility
    of the given mobilizer, assuming there is only one mobility. **/
    ConstantAcceleration(MobilizedBody& mobilizer, 
                         Real defaultAcceleration);
    
    /** Default constructor creates an empty handle you can use to reference
    any existing %ConstantAcceleration Constraint. **/
    ConstantAcceleration() {}

    /** Return the index of the mobilized body to which this constant 
    acceleration constraint is being applied (to \e one of its mobilities). 
    This is set on construction of the %ConstantAcceleration constraint. **/
    MobilizedBodyIndex getMobilizedBodyIndex() const;

    /** Return the particular mobility whose generalized acceleration is 
    controlled by this %ConstantAcceleration constraint. This is set on 
    construction. **/
    MobilizerUIndex    getWhichU() const;

    /** Return the default value for the acceleration to be enforced. This is 
    set on construction or via setDefaultAcceleration(). This is used to 
    initialize the acceleration when a default State is created, but it can be 
    overriden by changing the value in the State using setAcceleration(). **/
    Real               getDefaultAcceleration() const;

    /** Change the default value for the acceleration to be enforced by this 
    constraint. This is a topological change, meaning you'll have to call
    realizeTopology() on the containing System and obtain a new State before
    you can use it. If you just want to make a runtime change in the State,
    see setAcceleration(). **/
    ConstantAcceleration& setDefaultAcceleration(Real accel);

    /** Override the default acceleration with this one whose value is stored 
    in the given State. This invalidates the Acceleration stage in the state. 
    Don't confuse this with setDefaultAcceleration() -- the value set here 
    overrides that one. **/
    void setAcceleration(State& state, Real accel) const;

    /** Get the current value of the acceleration set point from the indicated 
    State. This is the value currently in effect, either from the default or 
    from a previous call to setAcceleration(). **/
    Real getAcceleration(const State& state) const;

    // no position or velocity error

    /** Return the amount by which the accelerations in the given State fail
    to satify this constraint. The \a state must already be realized through 
    Stage::Acceleration. **/
    Real getAccelerationError(const State&) const;

    /** Get the value of the Lagrange multipler generated to satisfy this
    constraint. For a %ConstantAcceleration constraint, that is the same as the
    generalized force although by convention constraint multipliers have the
    opposite sign from applied forces. The \a state must already be realized 
    through Stage::Acceleration.**/
    Real getMultiplier(const State&) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (ConstantAcceleration, ConstantAccelerationImpl, Constraint);
    /** @endcond **/
};

//==============================================================================
//                                 CUSTOM
//==============================================================================
/** The handle class Constraint::Custom (dataless) and its companion class 
Constraint::Custom::Implementation can be used together to define new 
Constraint types with arbitrary properties. To use it, create a class that 
extends Constraint::Custom::Implementation. You can then create an instance 
of it and pass it to the Constraint::Custom constructor:

<pre>
Constraint::Custom myConstraint(new MyConstraintImplementation( args ));
</pre>

Alternatively, you can also create a new Handle class which is a subclass of 
Constraint::Custom and which creates the Implementation itself in its 
constructors.

<pre>
class MyConstraint : public Constraint::Custom {
public:
  MyConstraint( args ) : Constraint::Custom(new MyForceImplementation( args )) {
  }
}
</pre>

This allows an end user to simply write

<pre>
MyConstraint( args );
</pre>

and not worry about implementation classes or creating objects on the heap.  
If you do this, your Constraint::Custom subclass must not have any data 
members or virtual methods.  If it does, it will not work correctly. Instead, 
store all data in the Implementation subclass. **/
class SimTK_SIMBODY_EXPORT Constraint::Custom : public Constraint {
public:
    class Implementation;
    class ImplementationImpl;

    /** Create a Custom Constraint.
     * 
     * @param implementation
     *      The object which implements the custom constraint. The 
     *      Constraint::Custom takes over ownership of the implementation 
     *      object, and deletes it when the Constraint itself is deleted.
     */
    explicit Custom(Implementation* implementation);

    
    /** Default constructor creates an empty handle. **/
    Custom() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, Constraint);
protected:
    const Implementation& getImplementation() const;
    Implementation&       updImplementation();
};

//==============================================================================
//                           CUSTOM::IMPLEMENTATION
//==============================================================================

// We only want the template instantiation to occur once. This symbol is 
// defined in the SimTK core compilation unit that defines the Constraint 
// class but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_CONSTRAINT
    extern template class PIMPLHandle<Constraint::Custom::Implementation, 
                                      Constraint::Custom::ImplementationImpl>;
#endif

/** This is the abstract base class for the implementation of custom 
constraints.\ See Constraint::Custom for more information. **/
class SimTK_SIMBODY_EXPORT Constraint::Custom::Implementation 
:   public PIMPLHandle<Implementation,ImplementationImpl> {
public:
// No default constructor because you have to supply at least the 
// SimbodyMatterSubsystem to which this Constraint belongs.

/** Destructor is virtual so derived classes get a chance to clean up if 
necessary. **/
virtual ~Implementation() { }

/** This method should produce a deep copy identical to the concrete derived 
Implementation object underlying this Implementation base class object. Note 
that the result is new heap space; the caller must be sure to take ownership of
the returned pointer and call delete on it when done. **/
virtual Implementation* clone() const = 0;

/** This Implementation base class constructor sets the topological defaults 
for the number of position level (holonomic), velocity level (nonholonomic), 
and acceleration-only constraint equations to be generated. **/
Implementation(SimbodyMatterSubsystem&, int mp, int mv, int ma);

/** The default constructor for the Implementation base class sets the number
of generated equations to zero for this constraint, meaning the Constraint 
won't do anything by default. The actual number can be changed using 
setNumConstraintEquationsInUse() prior to realizeModel(). **/ 
explicit Implementation(SimbodyMatterSubsystem&);

/** Return a reference to the matter subsystem containing this constraint. **/
const SimbodyMatterSubsystem& getMatterSubsystem() const;

    // Topological information//

/** Call this if you want to make sure that the next realizeTopology() call 
does something. This is done automatically when you modify the constraint in 
ways understood by Simbody, such as adding a ConstrainedBody. But if you are 
just changing some of your own topology and want to make sure you get a chance 
to recompute something in realizeTopology(), make this call at the time of 
modification. **/
void invalidateTopologyCache() const;

/** This is an alternate way to set the default number of equations to be 
generated if you didn't specify them in the base class constructor. A reference
to this Implementation is returned so that this can be used in a sequence like
an assignment operator. **/
Implementation& setDefaultNumConstraintEquations(int mp, int mv, int ma);

/** Normally Constraints are enabled when defined and can be disabled later. If
you want to define this constraint but have it be off by default, use this
method. A reference to this Implementation is returned so that this can be
used in a sequence like an assignment operator. **/
Implementation& setDisabledByDefault(bool shouldBeDisabled);

/** Call this during construction phase to add a body to the topological 
structure of this Constraint. This body's mobilizer's mobilities are \e not 
part of the constraint; mobilizers must be added separately. Numbering starts 
from 0 for each Constraint. The supplied MobilizedBody must be in the Matter 
Subsystem of which this Constraint is a part. **/
ConstrainedBodyIndex addConstrainedBody(const MobilizedBody&);

/** Call this during construction phase to add a mobilizer to the topological 
structure of this Constraint. All the coordinates q and mobilities u for this 
mobilizer are added also, but we don't know how many of those there will be 
until Stage::Model. Numbering starts from 0 for each Constraint. The supplied 
MobilizedBody must be in the Matter Subsystem of which this Constraint is 
a part. **/
ConstrainedMobilizerIndex addConstrainedMobilizer(const MobilizedBody&);

/** Map a constrained body for this constraint to the mobilized body to which
it corresponds in the matter subsystem. You should not use this to extract
any information in the constraint error or forces methods; always work with
the constrained bodies and constrained mobilities instead. **/
MobilizedBodyIndex 
getMobilizedBodyIndexOfConstrainedBody(ConstrainedBodyIndex) const;

/** Map a constrained mobilizer for this constraint to the mobilized body to 
which it corresponds in the matter subsystem. You should not use this to 
extract any information in the constraint error or forces methods; always work
with the constrained bodies and constrained mobilities instead. **/
MobilizedBodyIndex 
getMobilizedBodyIndexOfConstrainedMobilizer(ConstrainedMobilizerIndex) const;


/** @name       Methods for use with ConstrainedMobilizers
When a constraint acts directly on generalized coordinates q or generalized
speeds u (or their time derivatives), use methods in this section to access
those values in your constraint error and force methods. The "from state"
methods should only be used to pull information from the state that is at a
higher level than the method being written. For example, if you are calculating
velocity errors you can get positions from the state, but not velocities.
Instead, the velocities will be passed as an argument. **/
/**@{**/


/** Use this method in your calcPositionErrors() implementation to extract the
value of a particular generalized coordinate q selected by (mobilizer,whichQ),
from the "constrained q" argument that is passed to the method from Simbody. 
@param[in]      state
    Supplied state which is used only for modeling information; generalized
    coordinates q within \a state are ignored.
@param[in]      constrainedQ
    This is the argument that is supplied to calcPositionErrors() from which
    we will extract the particular q value selected by the next two arguments.
@param[in]      mobilizer
    The constrained mobilizer one of whose generalized coordinates is of
    interest.
@param[in]      whichQ
    The particular generalized coordinate of \a mobilizer whose value we
    want. The actual value will be selected from \a constrainedQ.
@return
    The value of the generalized coordinate q of interest, extracted from
    the \a constrainedQ argument. **/
Real getOneQ(const State&                           state,
             const Array_<Real,ConstrainedQIndex>&  constrainedQ,
             ConstrainedMobilizerIndex              mobilizer, 
             MobilizerQIndex                        whichQ) const;

/** Same as the getOneQ() method but for use in methods to which no explicit
"constrained q" argument is supplied. The desired q value is obtained
from \a state. You can call this from any constraint implementation method
\e except calcPositionError(). **/
Real getOneQFromState(const State&              state, 
                      ConstrainedMobilizerIndex mobilizer, 
                      MobilizerQIndex           whichQ) const;

/** Use this method in your calcPositionDotErrors() implementation to extract
the value of a particular generalized coordinate derivative qdot selected by 
(mobilizer,whichQ), from the "constrained qdot" argument that is passed to the 
method from Simbody. 
@param[in]      state
    Supplied state which is used only for modeling information; qdots within 
    \a state are ignored.
@param[in]      constrainedQDot
    This is the argument that is supplied to calcPositionDotErrors() from which
    we will extract the particular qdot value selected by the next two 
    arguments.
@param[in]      mobilizer
    The constrained mobilizer one of whose generalized coordinates is of
    interest.
@param[in]      whichQ
    The particular generalized coordinate of \a mobilizer whose qdot value we
    want. The actual value will be selected from \a constrainedQDot.
@return
    The value of the generalized coordinate derivative qdot of interest, 
    extracted from the \a constrainedQDot argument. **/
Real getOneQDot(const State&                           state,
                const Array_<Real,ConstrainedQIndex>&  constrainedQDot,
                ConstrainedMobilizerIndex              mobilizer, 
                MobilizerQIndex                        whichQ) const;

/** Same as the getOneQDot() method above but for use in velocity- or 
acceleration-level methods to which no explicit "constrained qdot" argument 
is supplied. The desired qdot value is obtained from \a state. You can call 
this from calcPositionDotDotError(). State must already be realized to the 
Velocity stage. **/
Real getOneQDotFromState(const State&             state, 
                        ConstrainedMobilizerIndex mobilizer, 
                        MobilizerQIndex           whichQ) const;


/** Use this method in your calcPositionDotDotErrors() implementation to extract
the value of a particular generalized coordinate second derivative qdotdot 
selected by (mobilizer,whichQ), from the "constrained qdotdot" argument that 
is passed to the method from Simbody. 
@param[in]      state
    Supplied state which is used only for modeling information; qdotdots within 
    \a state are ignored.
@param[in]      constrainedQDotDot
    This is the argument that is supplied to calcPositionDotDotErrors() from 
    which we will extract the particular qdotdot value selected by the next two 
    arguments.
@param[in]      mobilizer
    The constrained mobilizer one of whose generalized coordinates is of
    interest.
@param[in]      whichQ
    The particular generalized coordinate of \a mobilizer whose qdotdot value
    we want. The actual value will be selected from \a constrainedQDotDot.
@return
    The value of the generalized coordinate second derivative qdotdot of 
    interest, extracted from the \a constrainedQDotDot argument. 
    
There is no getOneQDotDotFromState() method because all the acceleration-
level methods are passed qdotdot or udot as an explicit argument. **/
Real getOneQDotDot(const State&                           state,
                   const Array_<Real,ConstrainedQIndex>&  constrainedQDotDot,
                   ConstrainedMobilizerIndex              mobilizer, 
                   MobilizerQIndex                        whichQ) const;

/** Use this method in your calcVelocityErrors() implementation to extract the
value of a particular generalized speed u selected by (mobilizer,whichU),
from the "constrained u" argument that is passed to the method from Simbody. 
@param[in]      state
    Supplied state which is used only for modeling information; generalized
    speeds u within \a state are ignored.
@param[in]      constrainedU
    This is the argument that is supplied to calcVelocityErrors() from which
    we will extract the particular u value selected by the next two arguments.
@param[in]      mobilizer
    The constrained mobilizer one of whose generalized speeds is of interest.
@param[in]      whichU
    The particular generalized speed of \a mobilizer whose value we
    want. The actual value will be selected from \a constrainedU.
@return
    The value of the generalized speed u of interest, extracted from
    the \a constrainedU argument. **/
Real getOneU(const State&                           state,
             const Array_<Real,ConstrainedUIndex>&  constrainedU,
             ConstrainedMobilizerIndex              mobilizer, 
             MobilizerUIndex                        whichU) const;

/** Same as the getOneU() method but for use in velocity- or acceleration-
level methods to which no explicit "constrained u" argument is supplied. The 
desired u value is obtained from \a state. You can call this only from
calcVelocityDotError(), calcAccelerationError(), and any constraint force
method. The State needs to be realized only as high as Model stage, but don't 
use this value in calcPositionError() or addInPositionConstraintForces(). 
Those must be limited to dependencies on time and configuration only. **/
Real getOneUFromState(const State&              state,
                      ConstrainedMobilizerIndex mobilizer, 
                      MobilizerUIndex           whichU) const;

/** Use this method in your calcVelocityDotErrors() and calcAccelerationErrors()
implementations to extract the value of a particular generalized speed 
derivative udot selected by (mobilizer,whichU), from the "constrained udot" 
argument that is passed to these two methods from Simbody. 
@param[in]      state
    Supplied state which is used only for modeling information; udots within 
    \a state are ignored.
@param[in]      constrainedUDot
    This is the argument that is supplied to calcVelocityDotErrors() and
    calcAccelerationErrros() from which we will extract the particular udot 
    value selected by the next two arguments.
@param[in]      mobilizer
    The constrained mobilizer one of whose generalized speeds is of interest.
@param[in]      whichU
    The particular generalized speed of \a mobilizer whose udot value we
    want. The actual value will be selected from \a constrainedUDot.
@return
    The value of the generalized speed derivative udot of interest, 
    extracted from the \a constrainedUDot argument.    

There is no getOneUDotFromState() method because all the acceleration-
level methods are passed qdotdot or udot as an explicit argument. **/
Real getOneUDot(const State&                           state,
                const Array_<Real,ConstrainedUIndex>&  constrainedUDot,
                ConstrainedMobilizerIndex              mobilizer, 
                MobilizerUIndex                        whichU) const;

/** Apply a scalar generalized (mobility-space) force \a fu to a particular 
mobility of one of this %Constraint's Constrained Mobilizers, \e adding it in to
the appropriate slot of the mobilityForces vector, which is of length 
getNumConstrainedU() for this %Constraint. State need only have been realized 
to Model stage, but this is intended for use in Velocity-stage calls to
addInXXXConstraintForce() methods for nonholonomic (velocity) or acceleration-only
constraint equations. 
@see addInOneQForce() for use in position (holonomic) constraints **/
void addInOneMobilityForce
   (const State&                    state, 
    ConstrainedMobilizerIndex       mobilizer, 
    MobilizerUIndex                 whichU,
    Real                            fu, 
    Array_<Real,ConstrainedUIndex>& mobilityForces) const;

/** For use with holonomic (position) constraints, this method allows 
generalized forces to be applied in "q-space" rather than "u-space". 
A scalar q-space generalized force \a fq is applied to a particular 
generalized coordinate (q) of one of this position (holonomic) %Constraint's 
Constrained Mobilizers, \e adding it in to the appropriate slot of the qForces 
vector, which must be of length getNumConstrainedQ() for this %Constraint. 
State need only have been realized to Model stage, but this is intended for 
Position-stage use in the addInPositionConstraintForce() method for position
constraint equations. 

Simbody will convert these automatically to mobility (u) space as needed
via fu = ~N * fq, where N is block-diagonal kinematic coupling matrix that 
appears in the equation qdot = N*u.
@see addInOneMobilityForces() for velocity and acceleration-only constraint
equations **/
void addInOneQForce
   (const State&                    state, 
    ConstrainedMobilizerIndex       mobilizer, 
    MobilizerQIndex                 whichQ,
    Real                            fq, 
    Array_<Real,ConstrainedQIndex>& qForces) const;
/**@}**/


/** @name          Methods for use with Constrained Bodies
When a constraint is enforced (at least in part) by applying forces to bodies, 
use the methods in this section to access position, velocity, and acceleration
information about those constrained bodies. Note that you can pull higher-level
information from the state, but information at the current level for a method 
must be taken from the supplied arguments instead. For example, if you are 
writing an acceleration error routine, you can get time, position, and velocity
information from the state but must get acceleration information from the
body accelerations that are supplied by Simbody as arguments. **/
/**@{**/


/** Extract from the \a allX_AB argument the spatial transform X_AB giving the 
pose (orientation and location) of a Constrained Body B's body frame B in this 
constraint's Ancestor frame A. **/
const Transform& getBodyTransform
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            bodyB) const;
/** Convenient inline interface to getBodyTransform() that returns just the
orientation as the Rotation matrix R_AB. **/
const Rotation& getBodyRotation
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyTransform(allX_AB,bodyB).R(); }
/** Convenient inline interface to getBodyTransform() that returns just the
location part of B's pose in A, that is the vector p_AB from A's origin Ao
to B's origin Bo, expressed in A. **/
const Vec3& getBodyOriginLocation
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyTransform(allX_AB,bodyB).p(); }

/** Extract from the State cache the spatial transform X_AB giving the 
pose (orientation and location) of a Constrained Body B's body frame B in this 
constraint's Ancestor frame A. Do not use this method in a routine that has
an explicit argument providing the transforms X_AB; use the above 
getBodyTransform() method instead. **/
const Transform& getBodyTransformFromState
   (const State& state, ConstrainedBodyIndex B)    const; // X_AB
/** Convenient inline interface to getBodyTransformFromState() that returns 
just the orientation as the Rotation matrix R_AB. **/
const Rotation& getBodyRotationFromState
   (const State& state, ConstrainedBodyIndex bodyB) const
{   return getBodyTransformFromState(state,bodyB).R(); }
/** Convenient inline interface to getBodyTransformFromState() that returns 
just the location part of B's pose in A, that is the vector p_AB from A's 
origin Ao to B's origin Bo, expressed in A. **/
const Vec3& getBodyOriginLocationFromState
   (const State& state, ConstrainedBodyIndex bodyB) const
{   return getBodyTransformFromState(state,bodyB).p(); }

/** Extract from the \a allV_AB argument the spatial velocity V_AB giving the 
angular and linear velocity of a Constrained Body B's body frame B measured and 
expressed in this constraint's Ancestor frame A. **/
const SpatialVec& getBodyVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            bodyB) const;
/** Convenient inline interface to getBodyVelocity() that returns just the
angular velocity vector w_AB. **/
const Vec3& getBodyAngularVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyVelocity(allV_AB,bodyB)[0]; }
/** Convenient inline interface to getBodyVelocity() that returns just the
linear velocity vector v_AB. **/
const Vec3& getBodyOriginVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyVelocity(allV_AB,bodyB)[1]; }

/** Extract from the State cache the spatial velocity V_AB giving the angular
and linear velocity of a Constrained Body B's body frame B measured and 
expressed in this Constraint's Ancestor frame A. Do not use this method in a 
routine that has an explicit argument providing the spatial velocities V_AB; 
use the above getBodyVelocity() method instead. **/
const SpatialVec& getBodyVelocityFromState
   (const State& state, ConstrainedBodyIndex bodyB)     const; // V_AB
/** Convenient inline interface to getBodyVelocityFromState() that returns just
the angular velocity vector w_AB. **/
const Vec3& getBodyAngularVelocityFromState
   (const State& state, ConstrainedBodyIndex bodyB) const
{   return getBodyVelocityFromState(state,bodyB)[0]; }
/** Convenient inline interface to getBodyVelocityFromState() that returns just
the linear velocity vector v_AB. **/
const Vec3& getBodyOriginVelocityFromState
   (const State& state, ConstrainedBodyIndex bodyB) const 
{   return getBodyVelocityFromState(state,bodyB)[1]; }

/** Extract from the \a allA_AB argument the spatial acceleration A_AB giving 
the angular and linear acceleration of a Constrained Body B's body frame B 
measured and expressed in this constraint's Ancestor frame A. Note that there
is no getBodyAccelerationFromState() method because all acceleration-level
methods will be passed body accelerations explicitly. **/
const SpatialVec& getBodyAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            bodyB) const;
/** Convenient inline interface to getBodyAcceleration() that returns just the
angular acceleration vector b_AB. **/
const Vec3& getBodyAngularAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyAcceleration(allA_AB,bodyB)[0]; }
/** Convenient inline interface to getBodyAcceleration() that returns just the
linear acceleration vector a_AB. **/
const Vec3& getBodyOriginAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            bodyB) const
{   return getBodyAcceleration(allA_AB,bodyB)[1]; }

    // Calculate location, velocity, and acceleration for a given station.

/** Calculate the position p_AS in the Ancestor frame of a station S of a 
Constrained Body B, specified with the position vector p_BS (or more 
explicitly, p_BoS) from the B frame origin Bo to the point S, expressed
in the B frame. The return value is a position vector from the Ancestor 
frame's origin Ao to the location of the point S, expressed in the 
Ancestor frame. Cost is 18 flops. **/
Vec3 findStationLocation
   (const Array_<Transform, ConstrainedBodyIndex>&  allX_AB, 
    ConstrainedBodyIndex                            bodyB, 
    const Vec3&                                     p_BS) const 
{
    const Transform& X_AB = allX_AB[bodyB];
    return X_AB * p_BS; // re-measure and re-express
}
/** Same as findStationLocation() but for when you have to get the position
information from the \a state rather than from an explicit argument. 
Cost is 18 flops. **/
Vec3 findStationLocationFromState
   (const State&            state, 
    ConstrainedBodyIndex    bodyB, 
    const Vec3&             p_BS) const 
{   
    const Transform& X_AB = getBodyTransformFromState(state,bodyB);
    return X_AB * p_BS; // re-measure and re-express
}

/** Calculate the velocity v_AS in the Ancestor frame of a station S of a 
Constrained Body B, specified with the position vector p_BS (or more 
explicitly, p_BoS) from the B frame origin Bo to the point S, expressed
in the B frame. The return value v_AS is a vector expressed in the Ancestor 
frame, and is the time derivative taken in A of the position vector p_AS. 
Cost is 27 flops. **/
Vec3 findStationVelocity
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            bodyB, 
    const Vec3&                                     p_BS) const 
{   // p_BS_A is p_BS rexpressed in A but not shifted to Ao
    const Rotation&   R_AB   = getBodyRotationFromState(state,bodyB);
    const Vec3        p_BS_A = R_AB * p_BS; 
    const SpatialVec& V_AB   = allV_AB[bodyB];
    return V_AB[1] + (V_AB[0] % p_BS_A);    // v + w X r
}
/** Same as findStationVelocity() but for when you have to get the velocity
information from the \a state rather than from an explicit argument. 
Cost is 27 flops. **/
Vec3 findStationVelocityFromState
   (const State&            state, 
    ConstrainedBodyIndex    bodyB, 
    const Vec3&             p_BS) const 
{   // p_BS_A is p_BS rexpressed in A but not shifted to Ao
    const Rotation&   R_AB   = getBodyRotationFromState(state,bodyB);
    const Vec3        p_BS_A = R_AB * p_BS;
    const SpatialVec& V_AB   = getBodyVelocityFromState(state,bodyB);
    return V_AB[1] + (V_AB[0] % p_BS_A);    // v + w X r
}

/** Calculate the acceleration a_AS in the Ancestor frame of a station S of
a Constrained Body B, specified with the position vector p_BS (or more 
explicitly, p_BoS) from the B frame origin Bo to the point S, expressed in the 
B frame. The return value a_AS is a vector expressed in the Ancestor frame,
and is the time derivative taken in A of the velocity vector v_AS and hence the
second derivative taken in A of the position vectory p_AS. Note that there
is no findStationAccelerationFromState() method because all acceleration-level
routines here are provided acceleration information in explicit arguments. 
Cost is 48 flops. **/
Vec3 findStationAcceleration
   (const State&                                    state, 
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            bodyB, 
    const Vec3&                                     p_BS) const 
{   // p_BS_A is p_BS rexpressed in A but not shifted to Ao
    const Rotation&   R_AB   = getBodyRotationFromState(state,bodyB);
    const Vec3        p_BS_A = R_AB * p_BS; 
    const Vec3&       w_AB   = getBodyAngularVelocityFromState(state,bodyB);
    const SpatialVec& A_AB   = allA_AB[bodyB];

    // Result is a + b X r + w X (w X r).
    // ("b" is angular acceleration; w is angular velocity).
    const Vec3 a_AS = A_AB[1] + (A_AB[0] % p_BS_A) 
                              + w_AB % (w_AB % p_BS_A); // % is not associative
    return a_AS;
}

    // Utilities for applying constraint forces to ConstrainedBodies.

/** Apply an Ancestor-frame force to a B-frame station S given by the position 
vector p_BS (or more explicitly, p_BoS) from the B frame origin Bo to the point
S, expressed in the B frame, <em>adding to</em> the appropriate 
\p bodyForcesInA entry for this ConstrainedBody B. **/
void addInStationForce
   (const State&                                state,  
    ConstrainedBodyIndex                        bodyB,
    const Vec3&                                 p_BS, 
    const Vec3&                                 forceInA, 
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA) const;

/** Apply an Ancestor-frame torque to body B, <em>adding to</em> the 
appropriate \p bodyForcesInA entry for this ConstrainedBody B. **/
void addInBodyTorque
   (const State&                                state, 
    ConstrainedBodyIndex                        bodyB,
    const Vec3&                                 torqueInA, 
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA) const;
/**@}**/

/** @name             Utility methods
These provide access to quantities associated with this constraint,
suitable for use in the optional realize() virtual methods. **/
/**@{**/
/** Given a \a state as passed to your realizeAcceleration() implementation,
obtain the multipliers that Simbody just calculated for this %Constraint. **/
void getMultipliers(const State&  state, 
                    Array_<Real>& multipliers) const;
/**@}**/

protected:
/** @name             Optional realize() virtual methods
Provide implementations of these methods if you want to allocate State 
variables (such as modeling options or parameters) or want to pre-calculate 
some expensive quantities and store them in the State cache for your future 
use. Note that the Position, Velocity, and Acceleration-stage realize methods 
will be called <em>after</em> the constraint error calculating methods 
associated with this Constraint's constraint equations have been used by
Simbody to perform any constraint calculations. That means, for example, you 
can access calculated multipliers from your realizeAcceleration() method. **/

/**@{**/
/** The Matter Subsystem's realizeTopology() method will call this method after
all MobilizedBody topology has been processed. This gives the Constraint a 
chance to 
  - calculate Topology stage "cache" values (mutable values which are 
    stored in the derived Implementation class directly), and
  - allocate Model-stage state variables for later use, and
  - allocate Model-stage cache entries in the State.
The indices to the Model-stage state & cache entries must be stored locally as 
part of the Topology-stage cache. **/
virtual void realizeTopology(State&) const { }

/** The Matter Subsystem's realizeModel() method will call this method after 
all MobilizedBody Model-stage processing has been done. This gives the 
Constraint a chance to 
  - calculate Model stage cache values according to the settings of the 
    Model variables,
  - allocate any later-Stage variables that may be needed (typically these 
    will be Instance stage variables containing geometric information or 
    constraint parameters like lengths or velocities.
The indices to any of the State entries allocated here must be stored in the 
State as part of the Model-stage cache. **/
virtual void realizeModel(State&) const { }

/** The Matter Subsystem's realizeInstance() method will call this method after
all MobilizedBody Instance-stage processing has been done. This gives the 
Constraint a chance to 
  - calculate Instance stage cache values according to the settings of the 
    Instance variables. **/
virtual void realizeInstance(const State&) const { }

/** The Matter Subsystem's realizeTime() method will call this method after any
MobilizedBody Time-stage processing has been done. This gives the Constraint
a chance to 
  - calculate Time stage cache values according to the current value of 
    time found in the State. **/
virtual void realizeTime(const State&) const { }

/** The Matter Subsystem's realizePosition() method will call this method after
any MobilizedBody Position-stage processing has been done, and \e after the
call has been made to your calcPositionErrors() operator. This gives the 
Constraint a chance to 
  - calculate Position stage cache values according to the current values of 
    positions and position errors found in the State. **/
virtual void realizePosition(const State&) const { }

/** The Matter Subsystem's realizeVelocity() method will call this method after
any MobilizedBody Velocity-stage processing has been done, and \e after your
calcVelocityErrors() and calcPositionDotErrors() operators have been called. 
This gives the Constraint a chance to 
  - calculate Velocity stage cache values according to the current values of 
    velocities and velocity errors found in the State. **/
virtual void realizeVelocity(const State&) const { }

/** The Matter Subsystem's realizeDynamics() method will call this method after
any MobilizedBody Dynamics-stage processing has been done. This gives the 
Constraint a chance to 
  - calculate Dynamics stage cache values according to the current values found
    in the State. **/
virtual void realizeDynamics(const State&) const { }

/** The Matter Subsystem's realizeAcceleration() method will call this method
after any MobilizedBody Acceleration-stage processing has been done, and 
\e after your calcAccelerationErrors(), calcVelocityDotErrors(), and
calcPositionDotDotErrors() operators have been called. This gives the 
Constraint a chance to 
  - calculate Acceleration stage cache values according to the current values 
    of body and mobility accelerations, acceleration errors, and multiplier
    values found in the state. **/
virtual void realizeAcceleration(const State&) const { }

/** The Matter Subsystem's realizeReport() method will call this method after 
any MobilizedBody Report-stage processing has been done. This gives the 
Constraint a chance to 
  - calculate Report stage cache values according to the current values found
    in the State. **/
virtual void realizeReport(const State&) const { }
/**@}**/

/** @name           Position (Holonomic) Constraint Virtuals
These must be defined if there are any position (holonomic) constraint
equations generated by this %Constraint. **/
/**@{**/

/** Calculate the \e mp position-constraint errors due to the position-level 
specification of a holonomic constraint and write them to \a perr, which will
have been allocated to length \e mp; do not reallocate it. When this is called,
\a state will already have been realized to Stage::Time; all position 
information used in your implementation must be taken from the passed-in 
arguments \a X_AB and \a constrainedQ, not from \a state. **/
virtual void calcPositionErrors     
   (const State&                                    state,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)       // mp of these
    const;

/** Calculate the \e mp velocity errors arising from the first time derivative
of the position-level holonomic constraint function calcPositionErrors(), and 
write them to \a pverr, which will have been allocated to length \e mp; do not 
reallocate it. When this is called, \a state will have already been realized to
Stage::Position; all velocity information used in your implementation must be 
taken from the passed-in arguments \a V_AB and \a constrainedQDot, not from 
\a state. However, you can obtain position information for the constrained 
bodies and constrained mobilizers from \a state using getOneQFromState(), 
getBodyTransformFromState(), and related methods. The implementation of this 
method must produce \e exactly the time derivative of the implementation of 
calcPositionErrors(). **/
virtual void calcPositionDotErrors      
   (const State&                                    state, // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) // mp of these
    const;

/** Calculate the \e mp errors arising from the second time derivative of the 
position-level holonomic constraint function calcPositionErrors(), and write 
them to \a paerr, which will have been allocated to length \e mp; do not 
reallocate it. When this is called, \a state will already have been realized to
Stage::Velocity; all acceleration-level information used in your implementation
must be taken from the passed-in arguments \a A_AB and \a constrainedQDotDot,
\e not from \a state. However, you can obtain position and velocity information
for the constrained bodies and constrained mobilizers from \a state using 
getOneQFromState(), getOneQDotFromState(), getBodyTransformFromState(), 
getBodyVelocityFromState(), and related methods. The implementation of this 
method must produce \e exactly the time derivative of the implementation of 
calcPositionDotErrors(). **/
virtual void calcPositionDotDotErrors     
   (const State&                                    state, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) // mp of these
    const;

/** From the \e mp supplied Lagrange multipliers provided in \a multipliers,
calculate the forces produced by this Constraint on its Constrained Bodies and
Constrained Qs. Body spatial forces are applied at the body origin and 
expressed in the Ancestor frame and written to an array \a bodyForcesInA of 
length getNumConstrainedBodies(). Q forces are written to an array \a qForces 
of length getNumConstrainedQ(), that is, the number of constrained 
<em>generalized coordinates</em> q, not the number of constrained \e mobilizers
or constrained \e mobilities u. When this is called, \a state will already have
been realized to Stage::Position and all position-stage cache information is 
available including any that may have been calculated during the prior call to 
this Constraint's calcPositionErrors() method and realizePosition() method. 
Simbody will already have ensured that the force-return arrays have been 
allocated to the right size and properly initialized; you need update only 
those to which you are applying forces.

@note Don't forget that you must <em>add in</em> your force contributions;
don't just write them or you'll wipe out all preceding constraints'
contributions! **/
virtual void addInPositionConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&       qForces) const;
/**@}**/

/** @name         Velocity (Nonholonomic) Constraint Virtuals
These must be defined if there are any velocity (nonholonomic) constraint 
equations generated by this Constraint. **/
/**@{**/

/** Calculate the \e mv velocity-constraint errors due to the velocity-level 
specification of a nonholonomic constraint and write them to \a verr, which will
already have been allocated to length \e mv; do not reallocate it. When this
is called, \a state will have been realized to Stage::Position; all 
velocity-level information used in your implementation must be taken from the 
passed-in arguments \a V_AB and \a constrainedU, not from \a state. However,
you may obtain time or any position-related information from \a state. 
A nonholonomic constraint may depend on \e any position information; you do
not have to limit that to constrained bodies and mobilizers as you do for
velocity-level information. **/
virtual void calcVelocityErrors     
   (const State&                                    state,  // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const;

/** Calculate the \e mv errors arising from the first time derivative
of the velocity-level specification of a nonholonomic constraint and write them
to \a vaerr, which will already have been allocated to length \e mv; do not 
reallocate it. When this is called, \a state will have been realized to 
Stage::Velocity; all acceleration-level information used in your implementation
must be taken from the passed-in arguments \a A_AB and \a constrainedUDot, 
\e not from \a state. However, you can obtain from \a state time, and any 
needed position and velocity information. The implementation of this method 
must produce \e exactly the time derivative of the implementation of 
calcVelocityErrors(). **/
virtual void calcVelocityDotErrors     
   (const State&                                    state,  // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const;

/** From the \p mv supplied Lagrange multipliers provided in \a multipliers,
calculate the forces produced by this Constraint on its Constrained Bodies and
Constrained Mobilities due to its velocity-level (nonholonomic) constraints. 
Body spatial forces are applied at the body origin and expressed in the 
Ancestor frame and written to an array \a bodyForcesInA of length 
getNumConstrainedBodies(). Mobility (generalized) forces are written to an 
array \a mobilityForces of length getNumConstrainedU(), that is, the number of 
constrained \e mobilities, not the number of constrained \e mobilizers. The 
supplied \a state will have been realized to Stage::Velocity and all 
position- and velocity-stage cache information is available including any that 
may have been calculated during the prior call to this constraint's 
realizePosition() and realizeVelocity() methods. Simbody will already have 
ensured that the force-return arrays have been allocated to the right 
size and initialized properly; you need only update the non-zero ones.

@note Don't forget that you must <em>add in</em> your force contributions;
don't just write them or you'll wipe out all preceding constraints'
contributions! **/
virtual void addInVelocityConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&       mobilityForces) const;
/**@}**/

/** @name           Acceleration-Only Constraint Virtuals
These must be defined if there are any acceleration-only constraint equations
generated by this Constraint. **/
/**@{**/

/** Calculate the \e ma acceleration-constraint errors due to the 
specification of an acceleration-only constraint and write them to \a aerr, 
which will already have been allocated to length \e ma; do not reallocate it. 
When this is called, \a state will have been realized to Stage::Velocity; all 
acceleration-level information used in your implementation must be taken from 
the passed-in arguments \a A_AB and \a constrainedUDot, \e not from \a state. 
However, an acceleration-only constraint may depend arbitrarily on time,
position, and velocity information which you may obtain freely from \a state;
you do not have to limit that to constrained bodies and mobilizers as you do 
for acceleration-level information. 

@note This method \e must be linear in the accelerations; Simbody has no
way to enforce that so it is up to you to do this correctly. **/
virtual void calcAccelerationErrors      
   (const State&                                    state, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr)  // ma of these
    const;

/** From the \e ma supplied Lagrange multipliers provided in \a multipliers,
calculate the forces produced by this Constraint on its Constrained Bodies and
Constrained Mobilities due to its acceleration-only constraints. Body spatial 
forces are applied at the body origin and expressed in the Ancestor frame and 
written to an array \a bodyForcesInA of length getNumConstrainedBodies(). 
Mobility forces are written to an array \a mobilityForces of length 
getNumConstrainedU(), that is, the number of constrained \e mobilities, not the
number of constrained \e mobilizers. The \a state will have been realized to 
Stage::Velocity and all position- and velocity-stage cache information is 
available including any that may have been calculated during the prior call to
this constraint's realizePosition() and realizeVelocity() methods. Simbody will
already have ensured that the force-return arrays have been allocated to the 
right size and initialized properly; you need only update the non-zero ones. 

@note Don't forget that you must <em>add in</em> your force contributions;
don't just write them or you'll wipe out all preceding constraints'
contributions! **/
virtual void addInAccelerationConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&       mobilityForces) const;
/**@}**/

/** Implement this optional method if you would like your constraint to 
generate any suggestions for geometry that could be used as default 
visualization as an aid to understanding a system containing this constraint. 
For example, if your constraint connects two points, you might want to draw a 
line between those points. You can also generate text labels, and you can
provide methods for controlling the presence or appearance of your generated 
geometry. If you don't implement this routine no geometry will be generated. **/
virtual void calcDecorativeGeometryAndAppend
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
}

friend class Constraint::CustomImpl;
};



//==============================================================================
//                           COORDINATE COUPLER
//==============================================================================
/** This is a %Constraint that uses a Function object to define a single 
holonomic (position) constraint equation acting to relate a set of generalized 
coordinates q. You provide a Function which takes some subset of the system's 
generalized coordinates as arguments, and returns a single value. It also must 
support partial derivatives up to second order. The constraint enforces that 
the value of the function should equal 0 at all times. For example, if you
wanted q1 and q2 to be constrained to have the same value you could define your 
function f as f=q1-q2. **/
class SimTK_SIMBODY_EXPORT Constraint::CoordinateCoupler 
:   public Constraint::Custom {
public:
    /** Create a CoordinateCoupler. You specify a Function and a list of 
    generalized coordinates to pass to it as arguments. Each generalized 
    coordinate is specified by a MobilizedBody and the index of the coordinate
    within its mobilizer. For example <pre>
        matter.getMobilizedBody(coordMobod[2]).getOneQ(state, coordQIndex[2])
    </pre> will be passed to the Function as the value of the second argument.
    
    @param matter      
        The matter subsystem to which this constraint will be added.
    @param function    
        The Function whose value should be maintained at zero by this 
        constraint at all times. The constraint takes over ownership of this 
        object, and automatically deletes in when the constraint is deleted.
    @param coordMobod   
        The MobilizedBody corresponding to each generalized coordinate that 
        should be passed as a function argument.
    @param coordQIndex  
        The index corresponding to each generalized coordinate that should be 
        passed as a function argument.
    **/
    CoordinateCoupler(SimbodyMatterSubsystem&           matter, 
                      const Function*                   function, 
                      const Array_<MobilizedBodyIndex>& coordMobod, 
                      const Array_<MobilizerQIndex>&    coordQIndex);

    /** For compatibility with std::vector; no copying is done. **/
    CoordinateCoupler(SimbodyMatterSubsystem& matter, const Function* function, 
                      const std::vector<MobilizedBodyIndex>& coordBody, 
                      const std::vector<MobilizerQIndex>& coordIndex) 
    {   // Invoke the above constructor with converted arguments.
        new (this) CoordinateCoupler(matter,function,
            ArrayViewConst_<MobilizedBodyIndex>(coordBody),
            ArrayViewConst_<MobilizerQIndex>(coordIndex));
    }
    
    /** Default constructor creates an empty handle. **/
    CoordinateCoupler() {}
};



//==============================================================================
//                              SPEED COUPLER
//==============================================================================
/** This is a %Constraint that uses a Function object to 
define a nonholonomic (velocity) constraint. You provide a Function which takes
some subset of the system's generalized speeds as arguments, and returns a 
single value. It also must support partial derivatives up to second order. The
constraint enforces that the value of the function should equal 0 at all times.
 
The Function may optionally depend on coordinates (q) as well as speeds (u), 
but it only acts as a constraint on the speeds. The constraint takes the 
current values of the coordinates as constants, then tries to modify only the 
speeds so as to satisfy the constraint. **/
class SimTK_SIMBODY_EXPORT Constraint::SpeedCoupler : public Constraint::Custom {
public:
    /** Create a SpeedCoupler.  You specify a Function and a list of generalized 
    speeds to pass to it as arguments. Each generalized speed is specified by a 
    MobilizedBody and the index of the speeds within that body.  For example
    matter.getMobilizedBody(bodies[2]).getOneU(state, speeds[2]) will be passed 
    to the function as the value of the second argument.
     
    @param      matter      
        The matter subsystem to which this constraint will be added.
    @param      function    
        The Function whose value should equal 0 at all times. The constraint 
        takes over ownership of this object, and automatically deletes it when 
        the constraint is deleted.
    @param      speedBody   
        The MobilizedBody corresponding to each generalized speed that should be
        passed as a function argument.
    @param      speedIndex  
        The index corresponding to each generalized speed that should be passed 
        as a function argument. **/
    SpeedCoupler(SimbodyMatterSubsystem&            matter, 
                 const Function*                    function, 
                 const Array_<MobilizedBodyIndex>&  speedBody, 
                 const Array_<MobilizerUIndex>&     speedIndex);

    /** For compatibility with std::vector; no copying is done. **/
    SpeedCoupler(SimbodyMatterSubsystem&                matter, 
                 const Function*                        function, 
                 const std::vector<MobilizedBodyIndex>& speedBody, 
                 const std::vector<MobilizerUIndex>&    speedIndex) 
    {   // Invoke above constructor with converted arguments.
        new (this) SpeedCoupler(matter, function,
                                ArrayViewConst_<MobilizedBodyIndex>(speedBody),
                                ArrayViewConst_<MobilizerUIndex>(speedIndex));
    }

    /** Create a SpeedCoupler. You specify a Function and a list of generalized 
    coordinates and speeds to pass to it as arguments. Each generalized speed is
    specified by a MobilizedBody and the index of the speeds within that body.
    For example matter.getMobilizedBody(bodies[2]).getOneU(state, speeds[2]) will
    be passed to the function as the value of the second argument. Generalized 
    coordinates come after generalized speeds in the argument list. For example, 
    if you specify three generalized speeds and two generalized coordinates, the
    Function must take a total of five arguments. The first three are the speeds,
    and the last two are the coordinates.
     
    @param matter      
        The matter subsystem to which this constraint will be added.
    @param function    
        The Function whose value should equal 0 at all times. The constraint 
        takes over ownership of this object, and automatically deletes it when 
        the constraint is deleted.
    @param speedBody   
        The MobilizedBody corresponding to each generalized speed that should be 
        passed as a function argument.
    @param speedIndex  
        The index corresponding to each generalized speed that should be passed
        as a function argument.
    @param coordBody   
        The MobilizedBody corresponding to each generalized coordinate that 
        should be passed as a function argument.
    @param coordIndex  
        The index corresponding to each generalized coordinate that should be 
        passed as a function argument. **/
    SpeedCoupler(SimbodyMatterSubsystem&            matter, 
                 const Function*                    function, 
                 const Array_<MobilizedBodyIndex>&  speedBody, 
                 const Array_<MobilizerUIndex>&     speedIndex,
                 const Array_<MobilizedBodyIndex>&  coordBody, 
                 const Array_<MobilizerQIndex>&     coordIndex);

    /** For compatibility with std::vector; no copying is done. **/
    SpeedCoupler(SimbodyMatterSubsystem&                matter, 
                 const Function*                        function, 
                 const std::vector<MobilizedBodyIndex>& speedBody, 
                 const std::vector<MobilizerUIndex>&    speedIndex,
                 const std::vector<MobilizedBodyIndex>& coordBody, 
                 const std::vector<MobilizerQIndex>&    coordIndex)
    {   // Invoke above constructor with converted arguments.
        new (this) SpeedCoupler(matter, function,
                                ArrayViewConst_<MobilizedBodyIndex>(speedBody),
                                ArrayViewConst_<MobilizerUIndex>(speedIndex),
                                ArrayViewConst_<MobilizedBodyIndex>(coordBody),
                                ArrayViewConst_<MobilizerQIndex>(coordIndex));
    }

    /** Default constructor creates an empty handle. **/
    SpeedCoupler() {}
};



//==============================================================================
//                             PRESCRIBED MOTION
//==============================================================================
/** This is a %Constraint that uses a Function to prescribe
the behavior of a single generalized coordinate as a function of time. You 
provide a Function which takes the current time as its argument and returns the
required value of the generalized coordinate. It also must support derivatives 
up to second order. **/
class SimTK_SIMBODY_EXPORT Constraint::PrescribedMotion 
:   public Constraint::Custom {
public:
    /** Create a PrescribedMotion constraint. You specify a Function that takes 
    time as its single argument, and returns the required value for the 
    constrained coordinate.
     
    @param      matter      
        The matter subsystem to which this constraint will be added.
    @param      function    
        The Function which specifies the value of the constrained coordinate.  
        The constraint takes over ownership of this object, and automatically 
        deletes it when the constraint is deleted.
    @param      coordBody   
        The MobilizedBody corresponding to the generalized coordinate which will
        be constrained.
    @param      coordIndex  
        The index of the generalized coordinate which will be constrained. **/
    PrescribedMotion(SimbodyMatterSubsystem&    matter, 
                     const Function*            function, 
                     MobilizedBodyIndex         coordBody, 
                     MobilizerQIndex            coordIndex);

    
    /** Default constructor creates an empty handle. **/
    PrescribedMotion() {}
};



} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_H_



