#ifndef SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_H_
#define SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_H_

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
Declares the Constraint::SphereOnSphereContact class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                        SPHERE ON SPHERE CONTACT
//==============================================================================
/** This constraint represents a \e bilateral connection between a sphere on one
body and a sphere on another.

On construction you may choose whether the connection enforces rolling;
otherwise the spheres will slip against each other. There is always one position 
(holonomic) constraint equation: the sphere surfaces must touch at some point. 
That leaves the spheres with five unconstrained degrees of freedom (dofs) to 
take on any relative orientation and to be touching at any point on their 
surfaces. Optionally, there are also two velocity (nonholonomic) constraint 
equations that prevent relative slip between the spheres and thus enforce 
rolling. In that case there can be five position dofs but only three 
unconstrained velocity dofs.

Note that this is a bilateral, unconditional connection and will push or pull 
as necessary to keep the spheres in contact. If rolling is being enforced then 
whatever tangential forces are necessary to keep the spheres rolling 
will be generated, regardless of the normal force. These constraints can form 
the basis for unilateral contact with Coulomb friction, but additional 
conditions must be added to determine when they are active.

There are two mobilized bodies involved, we'll call them F and B. Either body 
can move or either can be Ground; however, we will orient the signs in our 
calculations so that we treat F as though it were "fixed"; for example, the 
contact normal will be considered to point from F towards B, as though F were 
Ground. F has a sphere attached to it with center point Sf and radius rf. B has
a sphere attached to it with center Sb and radius rb. The line between the 
centers is p_SfSb, oriented from body F's sphere center to body B's sphere 
center. The contact normal will be a unit vector aligned with p_SfSb.

Let d=||p_SfSb|| be the separation between the sphere centers. We define the 
contact point Co to be located on the center-to-center line so that it divides
the line into two segments whose lengths are proportional to the size of the
sphere nearest each segment: <pre>
    Co = Sf + (rf/(rf+rb)) * p_SfSb
</pre>
If the contact constraint is satisfied (so that the spheres are touching), then
d=rf+rb and Co will be located at the point of contact, at a distance rf from 
Sf and rb from Sb. Otherwise Co will be located along the center line with the
error distributed in proportion to the radii. Since position constraints
cannot be satisifed perfectly, this ensures that the fractional torque error 
produced by applying tangential forces at a point away from the sphere surfaces
is the same for each sphere.

For any given pose, we will define a contact frame C with origin Co whose z
axis Cz is aligned with the center-to-center line, pointing from Co towards Sb.
The x and y axes Cx and Cy define the contact plane. They are chosen arbitrarily
but consistently for any relative pose of B in F, and held constant during 
subsequent velocity and acceleration calculations using the same pose. Cx and Cy
are used to parameterize the slip velocity and the tangential forces; the two 
tangential multipliers are measure numbers in the Cx and Cy directions. The
instantaneous contact frame in Ground is available as the transform X_GC.

The contact constraints here are enforced by a normal multiplier acting along
Cz, and optionally two tangential multipliers acting along Cx and Cy
respectively. Together these three multipliers can be interpreted as a force 
expressed in frame C, and applied equal and opposite to bodies F and B at
their material points coincident with contact point Co.

The assembly condition is the same as the position constraint: the surfaces
of the spheres must meet at some point. There is no assembly
condition for the tangential constraints since they do not restrict the
allowable pose during assembly. **/
class SimTK_SIMBODY_EXPORT Constraint::SphereOnSphereContact
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
  - The set...() methods return a reference to "this" %SphereOnSphereContact
    element (in the manner of an assignment operator) so they can be chained in
    a single expression. 
 
The default parameters can be overridden in any given State, and modified
without affecting Topology; see the "Runtime Changes" section below. **/
/*@{*/

/** Construct a sphere-on-sphere constraint as described in the 
Constraint::SphereOnSphereContact class documentation.

@param      mobod_F   
    The first MobilizedBody object to which a contacting sphere is attached.
    We'll call it F, the "fixed" body just to orient the contact; actually
    either or both bodies can be moving.
@param      defaultCenter_F
    This is the location of the sphere center on \a mobod_F, given as a vector 
    from the F frame origin Fo to the sphere center Sf, expressed in F. This is
    the center location that will be present in a default State; you can modify
    it later.
@param      defaultRadius_F
    The radius rf of the sphere attached to \a mobod_F. This is the value that
    will be present in a default State; you can modify it later.
@param      mobod_B
    The second MobilizedBody object to which a second contacting sphere is 
    attached. We'll call this mobilized body B. It can be Ground or a moving
    body.
@param      defaultCenter_B
    This is the location of the sphere center on \a mobod_B, given as a vector 
    from the B frame origin Bo to the sphere center Sb, expressed in B. This is
    the center location that will be present in a default State; you can modify
    it later.
@param      defaultRadius_B
    The radius rb of the sphere attached to \a mobod_B. This is the value that
    will be present in a default State; you can modify it later.
@param      enforceRolling
    Whether to generate tangential forces to make the sphere roll on the
    plane. Otherwise only a normal force is generated and the sphere is
    free to slip.
**/
SphereOnSphereContact(  MobilizedBody&      mobod_F, 
                        const Vec3&         defaultCenterOnF, 
                        Real                defaultRadiusOnF, 
                        MobilizedBody&      mobod_B, 
                        const Vec3&         defaultCenterOnB,
                        Real                defaultRadiusOnB,
                        bool                enforceRolling);
    
/** Default constructor creates an empty handle that can be used to
reference any %SphereOnSphereContact %Constraint.**/
SphereOnSphereContact() {}

/** Return a reference to the first MobilizedBody to which a sphere is
attached. This refers to the \a mobod_F that was given in the 
constructor and cannot be changed after construction. **/
const MobilizedBody& getMobilizedBodyF() const;

/** Return a reference to the second MobilizedBody to which a sphere is
attached. This refers to the \a mobod_B that was given in the 
constructor and cannot be changed after construction. **/
const MobilizedBody& getMobilizedBodyB() const;

/** Report whether this %Constraint was constructed to generate rolling
constraints (otherwise it is frictionless). This cannot be changed after
construction. **/
bool isEnforcingRolling() const;

/** Replace the default center point for the sphere attached to the first body,
\a mobod_F, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setCenterOnF() **/
SphereOnSphereContact& 
setDefaultCenterOnF(const Vec3& defaultCenter);
/** Replace the default radius for the sphere attached to the first body,
\a mobod_F, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setRadiusOnF() **/
SphereOnSphereContact& 
setDefaultRadiusOnF(Real defaultRadius);
/** Replace the default center point for the sphere attached to the second
body, \a mobod_B, that was supplied on construction. This is a topological 
change; you'll have to call realizeTopology() again if you call this method. 
@see setCenterOnB() **/
SphereOnSphereContact& 
setDefaultCenterOnB(const Vec3& defaultCenter);
/** Replace the default radius for the sphere attached to the second body,
\a mobod_B, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setRadiusOnB() **/
SphereOnSphereContact& 
setDefaultRadiusOnB(Real defaultRadius);

/** Return the default center point for the sphere attached to the first
body, \a mobod_F, as set during construction or by the most recent call to 
setDefaultCenterOnF(). 
@see getCenterOnF() **/
const Vec3& getDefaultCenterOnF() const;
/** Return the default radius for the sphere attached to the first body, 
\a mobod_F, as set during construction or by the most recent call to 
setDefaultRadiusF(). 
@see getRadiusOnF() **/
Real getDefaultRadiusOnF() const;
/** Return the default center point for the sphere attached to the second
body, \a mobod_B, as set during construction or by the most recent call to 
setDefaultCenterOnB(). 
@see getCenterOnB() **/
const Vec3& getDefaultCenterOnB() const;
/** Return the default radius for the sphere attached to the second body, 
\a mobod_B, as set during construction or by the most recent call to 
setDefaultRadiusB(). 
@see getRadiusOnB() **/
Real getDefaultRadiusOnB() const;
/*@}............................ Construction ................................*/

//------------------------------------------------------------------------------
/** @name                       Visualization
If you are allowing the SimbodyMatterSubsystem to generate default geometry,
this constraint element will attempt to represent the spheres that
it connects. Methods here give you some limited control over the generated
default geometry; if you need more control you should disable the generation
of default geometry and make your own. **/
/*@{*/
    // nothing yet
/*@}........................... Visualization ................................*/


//------------------------------------------------------------------------------
/** @name                      Runtime Changes
These refer to Position-stage discrete state variables that determine the sphere
parameters to be used to calculate constraint forces from a given State object.
If these are not set explicitly, the parameters are set to those provided in the
constructor or via the correponding setDefault...() methods. 
Note:
  - Changing one of these parameters invalidates the given State's 
    Stage::Position, meaning that the State's stage will be no higher than
    Stage::Time after the parameter change. That ensures that position 
    constraint errors will be recalculated before they are used.
  - The set...() methods here return a const reference to "this" 
    %SphereOnSphereContact element (in the manner of an assignment operator, 
    except read-only) so they can be chained in a single expression. 
    
You can also modify and examine the default parameters; see the "Construction"
section above. **/
/*@{*/

/** Modify the location of the sphere on the first body, \a mobod_F, in this 
\a state by providing a new vector p_FSf giving the sphere center location Sf
relative to body F's body frame origin Fo, and expressed in the F frame. This 
overrides the \a defaultCenterOnF in the given \a state, whose Stage::Position 
is invalidated. **/
const SphereOnSphereContact& 
setCenterOnF(State& state, const Vec3& sphereCenter) const;

/** Modify the radius of the sphere on the first body, \a mobod_F, in this 
\a state. This overrides the \a defaultRadiusOnF in the given \a state, whose 
Stage::Position is invalidated. **/
const SphereOnSphereContact& 
setRadiusOnF(State& state, Real sphereRadius) const;

/** Modify the location of the sphere on the second body, \a mobod_B, in this 
\a state by providing a new vector p_BSb giving the sphere center location Sb
relative to body B's body frame origin Bo, and expressed in the B frame. This 
overrides the \a defaultCenterOnB in the given \a state, whose Stage::Position 
is invalidated. **/
const SphereOnSphereContact& 
setCenterOnB(State& state, const Vec3& sphereCenter) const;

/** Modify the radius of the sphere on the second body, \a mobod_B, in this 
\a state. This overrides the \a defaultRadiusOnB in the given \a state, whose 
Stage::Position is invalidated. **/
const SphereOnSphereContact& 
setRadiusOnB(State& state, Real sphereRadius) const;

/** Return the position of the center point of the sphere attached to the
first body, \a mobod_F, as currently set in the given \a state. The value 
is returned as a vector p_FSf from body F's origin Fo to its sphere center Sf,
expressed in the F frame. **/
const Vec3& getCenterOnF(const State& state) const;
/** Return the radius of the sphere attached to the first body, \a mobod_F, as 
currently set in the given \a state. **/
Real getRadiusOnF(const State& state) const;

/** Return the position of the center point of the sphere attached to the
second body, \a mobod_B, as currently set in the given \a state. The value 
is returned as a vector p_BSb from body B's origin Bo to its sphere center Sb,
expressed in the B frame. **/
const Vec3& getCenterOnB(const State& state) const;
/** Return the radius of the sphere attached to the first body, \a mobod_F, as 
currently set in the given \a state. **/
Real getRadiusOnB(const State& state) const;
/*@}.......................... Runtime Changes ...............................*/

//------------------------------------------------------------------------------
/** @name                      Computations
Methods here provide access to values already calculated by this constraint
element, and provide operators you can call to calculate related values. **/
/*@{*/
/** The returned position error can be viewed as the signed distance between
the spheres. It is positive when the spheres are separated and measures their
closest approach distance. It is negative when the spheres are interpenetrating
and measures the penetration depth. The given \a state must have already been 
realized through Stage::Position. **/
Real getPositionError(const State& state) const;

/** The returned velocity error vector has the time derivative of the quantity
returned by getPositionError() in its z coordinate, and violation of the
rolling constraints in its x and y coordinates. If rolling is not being 
enforced then the x and y components are returned zero; they will not contain
the slip velocity in that case since any slip velocity is acceptable. Note
that the returned vector is expressed in the instantaneous contact frame C,
considered as fixed on body F. That is, this is the velocity of the contact
point on body B's sphere in the F frame, expressed in C. The given \a state must
have already been realized through Stage::Velocity. **/
Vec3 getVelocityErrors(const State& state) const;

/** This vector is the time derivative of the value returned by
getVelocityError(). Note that this is different than the acceleration of
the point of sphere B at the contact point because the contact point moves
with respect to that sphere. The given \a state must have already been realized
through Stage::Acceleration. **/
Vec3 getAccelerationErrors(const State& state) const;

/** This are the Lagrange multipliers required to enforce the constraint
equations generated here. For this %Constraint it has units of force, but 
recall that the sign convention for multipliers is the opposite of that for 
applied forces. Thus the returned value is the negative of the force being
applied to sphere B at the contact point, expressed in the contact frame C.
The x,y coordinates are the tangential force used to enforce rolling (or
zero if rolling is not being enforced), and the z coordinate is the force 
needed to enforce contact. Since this is an unconditional, bilateral 
constraint the multipliers may have any sign and magnitude. The given \a state 
must already be realized to Stage::Acceleration. **/
Vec3 getMultipliers(const State& state) const;

/** Return the force vector currently being applied by this constraint to
the contact point Co on the sphere attached to body B (the second body in the
constructor), expressed in the Ground frame. An equal and opposite force is 
applied to the sphere attached to body F, to its material point that
is at that same location in space. This is zero if the constraint is not
currently enabled. Tangential forces are generated only when rolling is being
enforced, but since this result is in the Ground frame all the vector measure
numbers may be non-zero regardless. The given \a state must already be realized
to Stage::Acceleration. **/
Vec3 findForceOnSphereBInG(const State& state) const;

/** Return the instantaneous contact frame C in the Ground frame. The actual
contact point is the origin Co of the returned frame, which is placed at a 
point along the center-to-center line segment as described in the class
documentation for this SphereOnSphereContact class. Note that this does not 
imply that the normal constraint is satisifed; point Co could be far from the 
sphere surfaces or embedded beneath them. The x-y directions Cx and Cy are
somewhat arbitrary, but they are stable for a given relative pose between 
bodies F and B. This calculates a valid value even 
if this constraint is currently disabled. The given \a state must already be 
realized to Stage::Position. **/
Transform findContactFrameInG(const State& state) const;

/** Calculate the separation distance or penetration depth of the two spheres.
It is positive when the spheres are separated and measures their closest 
approach distance. It is negative when the spheres are interpenetrating
and measures the penetration depth. This calculates a valid value even if this 
constraint is currently disabled. The given \a state must be realized to 
Stage::Position. **/
Real findSeparation(const State& state) const;
/*@}............................ Computations ................................*/

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (SphereOnSphereContact, SphereOnSphereContactImpl, Constraint);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_H_



