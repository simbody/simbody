#ifndef SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_H_
#define SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_H_

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
Declares the Constraint::PointOnPlaneContact class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                         POINT ON PLANE CONTACT
//==============================================================================

/** (Advanced) This is the underlying constraint for unilateral contact with
friction but must be combined with contact and friction conditions. This 
enforces the same normal condition as Constraint::PointInPlane but also adds
two velocity-level no-slip constraints. Thus there are three
constraints, one position (holonomic) constraint and two velocity 
(nonholonomic) constraints. Note that this is a bilateral
constraint and will push or pull as necessary to keep the point in contact
with the plane, and that sticking is enforced regardless of the amount of
normal force being generated. If you want to make this unilateral, you must
handle switching it on and off separately; when this constraint is enabled it
always enforces the contact and no-slip conditions.

There are two mobilized bodies involved. MobiliedBody S, the plane "surface" 
body, has a plane P fixed to it, with the plane defined by a frame P given 
relative to body S by the transform X_SP. MobiliedBody B, the "follower" body 
has a station point F (a vertex) fixed to it, given with respect to body B at 
location p_BF. The coordinate axes of the plane frame P (fixed in S) are used 
for parameterization of the %Constraint. The z direction Pz of that frame is the 
plane normal; the Px,Py directions are used to express the tangential velocity; 
the P frame origin Po provides the height h=dot(Po,Pz) of the plane
over the surface body's origin in the direction of the plane normal Pz.

The position constraint equation for contact enforces that the follower point F 
must always be in the plane P, that is, pz_PF=p_PF[2]=0, or equivalently <pre>
    (1) perr = pz_PF = dot(p_SF,Pz)-h 
</pre> 
That constraint equation is enforced by an internal 
(non-working) scalar force acting at the spatial location of the follower point, 
directed along the plane normal, and equal and opposite on the two bodies at F 
on B and at the instantaneously coincident point p_SF on S. Because position 
constraints are not enforced perfectly, contact will occur slightly above or 
slightly below the plane surface, wherever the follower point is.

The two velocity constraint equations enforce that the follower point has no
relative tangential velocity when measured in the plane frame P, that is we
want to enforce vx_PF=vy_PF=0. The time derivative of the normal constraint is 
also enforced at the velocity level so we have vz_PF=0. Taken together the
velocity-level constraints are just <pre>
    (2) verr = v_PF = ~R_SP * v_SF
</pre>
Equation (2) is the same as the velocity-level constraints for a ball joint 
between frame P on S and point F on B. The acceleration-level constraints are 
just the time derivative of (2), i.e. <pre>
    (3) aerr = a_PF = ~R_SP * a_SF
</pre>
The assembly condition is that the follower point must be in the plane; there
is no assembly condition in the tangential direction since those constraints
are at the velocity level only. **/

class SimTK_SIMBODY_EXPORT Constraint::PointOnPlaneContact 
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
  - The set...() methods return a reference to "this" %PointOnPlaneContact
    element (in the manner of an assignment operator) so they can be chained in
    a single expression. **/
/*@{*/
/** Construct a point-in-plane normal constraint and two no-slip friction
constraints as described in the Constraint::PointOnPlaneContact
class documentation. **/
PointOnPlaneContact(MobilizedBody&     planeBody, 
                    const Transform&   defaultPlaneFrame, 
                    MobilizedBody&     followerBody, 
                    const Vec3&        defaultFollowerPoint);
    
/** Default constructor creates an empty handle that can be used to
reference any %PointOnPlaneContact %Constraint.**/
PointOnPlaneContact() {}

/** Replace the default plane frame that was supplied on construction.
This is a topological change; you'll have to call realizeTopology() again
if you call this method. **/
PointOnPlaneContact& 
setDefaultPlaneFrame(const Transform& defaultPlaneFrame);
/** Replace the default follower point that was supplied on construction.
This is a topological change; you'll have to call realizeTopology() again
if you call this method. **/
PointOnPlaneContact& 
setDefaultFollowerPoint(const Vec3& defaultFollowerPoint);

/** Return the MobilizedBodyIndex of the plane MobilizedBody. **/
MobilizedBodyIndex getPlaneMobilizedBodyIndex() const;
/** Return the MobilizedBodyIndex of the follower MobilizedBody. **/
MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

/** Return the default plane frame as set during construction or by the
most recent call to setDefaultPlaneFrame(). **/
const Transform& getDefaultPlaneFrame() const;
/** Return the default follower point as set during construction or by the
most recent call to setDefaultFollowerPoint(). **/
const Vec3&      getDefaultFollowerPoint() const;
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
PointOnPlaneContact& setPlaneDisplayHalfWidth(Real halfWidth);
/** This affects only generated decorative geometry for default 
visualization; the point is really zero radius. If you don't set
this the default radius is .05 length unit. Set this to zero to 
disable any attempt to generate default visualization for the point. **/
PointOnPlaneContact& setPointDisplayRadius(Real radius);
/** Return the plane half-width that will be used if we're asked to generate
default visualization geometry. If this is zero we won't generate any plane
visualization geometry. **/
Real getPlaneDisplayHalfWidth() const;
/** Return the sphere radius that will be used to visualize the follower
point if we're asked to generate default visualization geometry. If this is 
zero we won't generate any point visualization geometry. **/
Real getPointDisplayRadius() const;
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

/** Return the plane frame X_SP that is currently in effect for this
%Constraint. Note that the origin of the returned frame will be exactly on
the plane surface; that is not necessarily where contact occurs since the
follower point may be above or below the surface by position constraint
tolerance. **/
const Transform& getPlaneFrame(const State& state) const;
const Vec3&      getFollowerPoint(const State& state) const;
/*@}.......................... Runtime Changes ...............................*/


//------------------------------------------------------------------------------
/** @name                      Computations
Methods here provide access to values already calculated by this constraint
element, and provide operators you can call to calculate related values. **/
/*@{*/
/** The returned position error can be viewed as a signed distance. It is
positive when the follower point is above the plane and negative when
it is below the plane. The given \a state must have already been 
realized through Stage::Position. **/
Real getPositionError(const State& state) const;

/** The velocity error vector is the velocity of the follower point in
the contact frame. The contact frame is parallel to the plane frame but
with its origin shifted to the spatial location as the follower point. 
Note that the coordinates are ordered x-y-z so the first two numbers are
the tangential slip velocity and the third is the velocity in the normal
direction, with positive indicating separation and negative indicating
approach. The z value here is the time derivative of the quantity returned
by getPositionError(). The given \a state must have already been 
realized through Stage::Velocity. **/
Vec3 getVelocityErrors(const State& state) const;

/** This vector is the time derivative of the value returned by
getVelocityError(). The given \a state must have already been realized
through Stage::Acceleration. **/
Vec3 getAccelerationErrors(const State& state) const;

/** These are the Lagrange multipliers required to enforce the three
constraint equations generated here. For this %Constraint they have units
of force, but the sign convention for multipliers is the opposite of that
for applied forces. Thus the returned vector may be considered the force
applied by the follower point to the coincident point of the plane body,
expressed in the contact frame. **/
Vec3 getMultipliers(const State& state) const;

/** This is the force applied by the plane body to the follower point,
expressed in the contact frame. For this %Constraint, the value returned
here is identical to the vector returned by getMultipliers(), but with the 
opposite sign. **/
Vec3 getForceOnFollowerPoint(const State& state) const;
/*@}............................ Computations ................................*/

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (PointOnPlaneContact, PointOnPlaneContactImpl, Constraint);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_POINT_ON_PLANE_CONTACT_H_



