#ifndef SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_H_
#define SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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
Declares the Constraint::SphereOnPlaneContact class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                        SPHERE ON PLANE CONTACT
//==============================================================================
/** This constraint represents a \e bilateral connection between a sphere on one
body and a plane on another.

On construction you may choose whether the connection enforces rolling;
otherwise the sphere will slip along the plane. There is always one position
(holonomic) constraint equation: the sphere must touch the plane. That leaves
the sphere free to take on any orientation and to be touching at any point in
the plane. Optionally, there are also two velocity (nonholonomic) constraint
equations that prevent relative slip between the sphere and the plane and thus
enforce rolling.

Note that this is a bilateral, unconditional connection and will push or pull
as necessary to keep the ball in contact with the plane. If rolling is being
enforced then whatever tangential forces are necessary to keep the ball rolling
will be generated, regardless of the normal force. These constraints can form
the basis for unilateral contact, but additional conditions must be added to
determine when they are active.

There are two mobilized bodies involved, we'll call them F ("floor") and B
("ball"). (These names are just mnemonics; both bodies can be free to move
or either one could be Ground.) F is a body to which a plane P has been fixed,
and B is a body to which a sphere S is fixed. The plane is defined by a frame P
given relative to body F by the transform X_FP. The coordinate axes of P are
used for parameterization of the %Constraint. The z direction Pz of that frame
is the plane normal; the Px,Py directions are used to express the tangential
slip velocity and tangential forces in case of rolling; the P frame origin Po
(given by the vector p_FP) provides the height h=p_FP.Pz of the plane over the
floor body's origin (Fo) in the direction of the plane normal Pz.

MobilizedBody B has a sphere fixed to it with center So and radius r, with So
given by the vector p_BS in the B frame. Call the point at the bottom (w.r.t.
plane normal) of the sphere C, with C=So-r*Pz. The position constraint equation
for contact then enforces that the location of the bottom of the sphere (at C)
is in the plane P, and all contact forces are applied at the location of point
C. Because position constraints cannot be enforced perfectly, contact will occur
slightly above or slightly below the plane surface, depending where C is.
However, contact will always occur at a radius of exactly r from the sphere
center S.

The contact constraints here are enforced by a normal multiplier acting along
Pz, and optionally two tangential multipliers acting along Px and Py
respectively. Together these can be interpreted as a force acting in a frame C
that is always aligned with P, but whose origin is at the contact point Co at
the bottom of the ball.

The assembly condition is the same as the position constraint: some point on
the surface of the sphere must be touching the plane. There is no assembly
condition for the tangential constraints since they do not restrict the
allowable pose during assembly. **/
class SimTK_SIMBODY_EXPORT Constraint::SphereOnPlaneContact
:   public Constraint  {
public:
    // no default constructor

//------------------------------------------------------------------------------
/** @name                       Construction
Methods in this section refer both to constructors, and to methods that can
be used to set or change contruction (Topology-stage) parameters; these
specify the values assigned by default to the corresponding state variables.
Note:
  - Changing one of these default parameters invalidates the containing
    System's topology, meaning that realizeTopology() will have to be called
    and a new State obtained before subsequent use.
  - The set...() methods return a reference to "this" %SphereOnPlaneContact
    element (in the manner of an assignment operator) so they can be chained in
    a single expression. **/
/*@{*/

/** Construct a sphere-on-plane constraint as described in the
Constraint::SphereOnPlaneContact class documentation.

@param      planeMobod
    The "floor" MobilizedBody F to which the plane P is attached.
@param      defaultPlaneFrame
    The Transform X_FP that defines the plane relative to the body
    frame of \a planeMobod F. This is the value that will be present in
    a default State; you can modify it later.
@param      sphereMobod
    The "ball" MobilizedBody B to which the sphere S is attached.
@param      defaultSphereCenter
    The center of sphere S defined relative to the body frame of the
    \a sphereMobod B. This is a vector p_BS from body B's origin Bo to
    the sphere center So, expressed in the B frame. This is the value that
    will be present in a default State; you can modify it later.
@param      defaultSphereRadius
    The radius r of the sphere S. This is the value that will be present in
    a default State; you can modify it later.
@param      enforceRolling
    Whether to generate tangential forces to make the sphere roll on the
    plane. Otherwise only a normal force is generated and the sphere is
    free to slip.
**/
SphereOnPlaneContact(   MobilizedBody&      planeMobod,
                        const Transform&    defaultPlaneFrame,
                        MobilizedBody&      sphereMobod,
                        const Vec3&         defaultSphereCenter,
                        Real                defaultSphereRadius,
                        bool                enforceRolling);

/** Default constructor creates an empty handle that can be used to
reference any %SphereOnPlaneContact %Constraint.**/
SphereOnPlaneContact() {}

/** Return a reference to the MobilizedBody to which the plane is
attached. This refers to the \a planeMobod that was given in the
constructor and cannot be changed after construction. **/
const MobilizedBody& getPlaneMobilizedBody() const;

/** Return a reference to the MobilizedBody to which the sphere is
attached. This refers to the \a sphereMobod that was given in the
constructor and cannot be changed after construction. **/
const MobilizedBody& getSphereMobilizedBody() const;

/** Report whether this %Constraint was constructed to generate rolling
constraints (otherwise it is frictionless). This cannot be changed after
construction. **/
bool isEnforcingRolling() const;


/** Replace the default plane frame that was supplied on construction.
This is a topological change; you'll have to call realizeTopology() again
if you call this method. **/
SphereOnPlaneContact&
setDefaultPlaneFrame(const Transform& defaultPlaneFrame);
/** Replace the default center point that was supplied on construction.
This is a topological change; you'll have to call realizeTopology() again
if you call this method. **/
SphereOnPlaneContact&
setDefaultSphereCenter(const Vec3& defaultSphereCenter);
/** Replace the default sphere radius that was supplied on construction.
This is a topological change; you'll have to call realizeTopology() again
if you call this method. **/
SphereOnPlaneContact&
setDefaultSphereRadius(Real defaultSphereRadius);


/** Return the default plane frame as set during construction or by the
most recent call to setDefaultPlaneFrame(). **/
const Transform& getDefaultPlaneFrame() const;
/** Return the default center point as set during construction or by the
most recent call to setDefaultSphereCenter(). **/
const Vec3& getDefaultSphereCenter() const;
/** Return the default sphere radius as set during construction or by the
most recent call to setDefaultSphereRadius(). **/
Real getDefaultSphereRadius() const;
/*@}............................ Construction ................................*/

//------------------------------------------------------------------------------
/** @name                       Visualization
If you are allowing the SimbodyMatterSubsystem to generate default geometry,
this constraint element will attempt to represent the plane and sphere that
it connects. Methods here give you some limited control over the generated
default geometry; if you need more control you should disable the generation
of default geometry and make your own. **/
/*@{*/
/** This affects only generated decorative geometry for default
visualization; the plane is really infinite in extent. If you don't set
this the default half width is 1 length unit. Set this to zero to
disable any attempt to generate default visualization for the plane. **/
SphereOnPlaneContact& setPlaneDisplayHalfWidth(Real halfWidth);
/** Return the plane half-width that will be used if we're asked to generate
default visualization geometry.  If this is zero we won't generate any plane
visualization geometry. **/
Real getPlaneDisplayHalfWidth() const;
/*@}........................... Visualization ................................*/


//------------------------------------------------------------------------------
/** @name                      Runtime Changes
These refer to Position-stage discrete state variables that determine the sphere
and plane parameters to be used to calculate constraint forces from a given
State object. If these are not set explicitly, the parameters are set to those
provided in the constructor or via the correponding setDefault...() methods.
Note:
  - Changing one of these parameters invalidates the given State's
    Stage::Position, meaning that the State's stage will be no higher than
    Stage::Time after the parameter change. That ensures that position
    constraint errors will be recalculated before they are used.
  - The set...() methods here return a const reference to "this"
    %SphereOnPlaneContact element (in the manner of an assignment operator,
    except read-only) so they can be chained in a single expression. **/
/*@{*/

/** Modify the location of the plane in this \a state by providing a new
transform X_FP giving the plane frame relative to the plane body F's body
frame. This overrides the defaultPlaneFrame in the given \a state, whose
Stage::Position is invalidated. **/
const SphereOnPlaneContact&
setPlaneFrame(State& state, const Transform& planeFrame) const;

/** Modify the location of the sphere in this \a state by providing a new
vector p_BS giving the sphere center location relative to the sphere body B's
body frame origin, and expressed in the B frame. This overrides the
defaultSphereCenter in the given \a state, whose Stage::Position is
invalidated. **/
const SphereOnPlaneContact&
setSphereCenter(State& state, const Vec3& sphereCenter) const;

/** Modify the radius of the sphere in this \a state. This overrides the
defaultSphereRadius in the given \a state, whose Stage::Position is
invalidated. **/
const SphereOnPlaneContact&
setSphereRadius(State& state, Real sphereRadius) const;

/** Return the plane frame X_FP that is currently in effect for this
%Constraint. Note that the origin of the returned frame will be exactly on
the plane surface; that is not necessarily where contact occurs since the
ball may be above or below the surface by position constraint
tolerance. **/
const Transform& getPlaneFrame(const State& state) const;
/** Return the sphere's center point location p_BO that is current in effect
for this %Constraint. The location is measured and expressed in the ball
body's frame. **/
const Vec3& getSphereCenter(const State& state) const;
/** Return the sphere radius that is currently in effect for this
%Constraint. **/
Real getSphereRadius(const State& state) const;
/*@}.......................... Runtime Changes ...............................*/

//------------------------------------------------------------------------------
/** @name                      Computations
Methods here provide access to values already calculated by this constraint
element, and provide operators you can call to calculate related values. **/
/*@{*/
/** The returned position error can be viewed as the signed distance from
the lowest point of the sphere to the plane surface. It is positive when the
sphere is above the plane (along the plane normal Pz) and negative when it is
penetrating the plane. The given \a state must have already been realized
through Stage::Position. **/
Real getPositionError(const State& state) const;

/** The returned velocity error vector has the time derivative of the quantity
returned by getPositionError() in its z coordinate, and violation of the
rolling constraints in its x and y coordinates. If rolling is not being
enforced then the x and y components are returned zero; they will not contain
the slip velocity in that case since any slip velocity is acceptable. Note
that the returned vector is expressed in the plane frame P, that is, this is
the velocity of the contact point of the sphere body, measured with respect
to the plane, and expressed along Px,Py, and Pz. The given \a state must
have already been realized through Stage::Velocity. **/
Vec3 getVelocityErrors(const State& state) const;

/** This vector is the time derivative of the value returned by
getVelocityError(). Note that this is different than the acceleration of
the point of the sphere at the contact point beceause the contact point moves
with respect to the sphere. The given \a state must have already been realized
through Stage::Acceleration. **/
Vec3 getAccelerationErrors(const State& state) const;

/** This are the Lagrange multipliers required to enforce the constraint
equations generated here. For this %Constraint it has units of force, but
recall that the sign convention for multipliers is the opposite of that for
applied forces. Thus the returned value is the negative of the force being
applied to the sphere at the contact point, expressed in the plane frame P.
The x,y coordinates are the forces in the plane used to enforce rolling (or
zero if rolling is not being enforced), and the z coordinate is the force
needed to enforce contact. Since this is an unconditional, bilateral
constraint the multipliers may have any sign and magnitude. The given \a state
must already be realized to Stage::Acceleration. **/
Vec3 getMultipliers(const State& state) const;

/** Return the force vector currently being applied by this constraint to
the contact point on the sphere body, expressed in the Ground frame. An equal
and opposite force is applied to the plane body, to its material point that
is at that same location in space. This is zero if the constraint is not
currently enabled. Tangential forces are generated only when rolling is being
enforced, but since this result is in the Ground frame all the vector measure
numbers may be non-zero regardless.  The given \a state must already be realized
to Stage::Acceleration. **/
Vec3 findForceOnSphereInG(const State& state) const;

/** Return the contact point location in the Ground frame. We define this
to be the Ground location coincident with the bottom of the sphere with
respect to the plane normal Pz. Note that this does not imply that the
constraint is satisifed; that point could be far from the plane or deeply
below the plane. This calculates a valid value even if this constraint is
currently disabled. The given \a state must already be realized to
Stage::Position. **/
Vec3 findContactPointInG(const State& state) const;

/** Calculate the separation distance or penetration depth of the sphere
and the plane. This is positive if the lowest point of the sphere is above
the plane (in the sense of the plane normal Pz), and negative if the
lowest point is below the plane, in which case it measures the maximum
penetration depth. This calculates a valid value even if this constraint is
currently disabled. The given \a state must be realized to
Stage::Position. **/
Real findSeparation(const State& state) const;
/*@}............................ Computations ................................*/

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (SphereOnPlaneContact, SphereOnPlaneContactImpl, Constraint);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_PLANE_CONTACT_H_



