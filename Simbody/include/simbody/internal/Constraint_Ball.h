#ifndef SimTK_SIMBODY_CONSTRAINT_BALL_H_
#define SimTK_SIMBODY_CONSTRAINT_BALL_H_

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
Declares the Constraint::Ball class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                  CONSTRAINT::  BALL (COINCIDENT POINTS)
//==============================================================================

/** Enforce that a fixed station on one body remains coincident with a fixed
station on a second body, as though there were a ball joint connecting them at
those points. Uses three position-level (holonomic) constraint equations to
prevent relative translation in three orthogonal directions.

At construction you specify the two bodies to be connected by the %Constraint,
and give default values for a station on each body. The State is initialized
with those default stations, but you can change them later.

The constraint is enforced by an internal (non-working) force applied at the
spatial location of the point on body 2, on material points of each body that
are coincident with that spatial location. Note that this is somewhat asymmetric
when the ball is not properly assembled -- it acts as though the contact occurs
at the point on body 2, \e not at the point on body 1. That is critical to
ensure that Newton's 3rd law is satisified -- the action and reaction must
occur at the same point.

The assembly condition is the same as the runtime constraint -- the two points
can be brought together by driving the perr to zero. **/
class SimTK_SIMBODY_EXPORT Constraint::Ball : public Constraint {
public:
/** Connect the origin of \a body1 to the origin of \a body2. That is,
the default stations will both be (0,0,0). You can change those later
in the State using setPointOnBody1() and setPointOnBody2(). **/
Ball(MobilizedBody& body1, MobilizedBody& body2);
/** Connect \a body1 and \a body2 at given station points, given in the
body frame of the corresponding body. You can change
those later in the State using setPointOnBody1() and setPointOnBody2(). **/
Ball(MobilizedBody& body1, const Vec3& defaultPoint1,
        MobilizedBody& body2, const Vec3& defaultPoint2);

/** Default constructor creates an empty handle. **/
Ball() {}

/** Change the station point on body 1 at which this %Constraint acts.
Provide the station location in the body 1 local frame.
This overrides the default point that was supplied on construction. This
is an Instance-stage change. **/
void setPointOnBody1(State& state, const Vec3& point_B1) const;
/** Change the station point on body 2 at which this %Constraint acts.
Provide the station location in the body 2 local frame.
This overrides the default point that was supplied on construction. This
is an Instance-stage change. **/
void setPointOnBody2(State& state, const Vec3& point_B2) const;

/** Return from the given \a state the constrained station on body 1, in
the body 1 frame. **/
const Vec3& getPointOnBody1(const State& state) const;
/** Return from the given \a state the constrained station on body 2, in
the body 2 frame. **/
const Vec3& getPointOnBody2(const State& state) const;

/** Change the default station location on body 1. This is a topological
change meaning you'll have to call realizeTopology() again after changing
the default point. If you want to change the constrained station during
a simulation, use setPointOnBody1() instead to override it in a State. **/
Ball& setDefaultPointOnBody1(const Vec3& defaultPoint_B1);
/** Change the default station location on body 2. This is a topological
change meaning you'll have to call realizeTopology() again after changing
the default point. If you want to change the constrained station during
a simulation, use setPointOnBody2() instead to override it in a State. **/
Ball& setDefaultPointOnBody2(const Vec3& defaultPoint_B2);

/** Return the default location for the station point on body 1, as a
vector in the body 1 frame. Note that
this is not necessarily the station point being used for any given State;
use getPointOnBody1() for that. **/
const Vec3& getDefaultPointOnBody1() const;
/** Return the default location for the station point on body 2, as a
vector in the body 2 frame. Note that
this is not necessarily the station point being used for any given State;
use getPointOnBody2() for that. **/
const Vec3& getDefaultPointOnBody2() const;


/** For visualization only, you can override the default radius used by
this %Constraint to draw itself. **/
Ball& setDefaultRadius(Real r);
/** Retrieve the radius being used for visualization of
this %Constraint. **/
Real getDefaultRadius() const;

/** Return the MobilizedBodyIndex corresponding to body 1. **/
MobilizedBodyIndex getBody1MobilizedBodyIndex() const;
/** Return the MobilizedBodyIndex corresponding to body 2. **/
MobilizedBodyIndex getBody2MobilizedBodyIndex() const;


/** Return the current position-level constraint error for this %Constraint.
This is the vector between the constrained stations on body 1 and body 2,
which would be zero if this constraint were perfectly satisfied. The
returned vector is measured in the Ancestor body frame. The given
\a state must be realized through Position stage. **/
Vec3 getPositionErrors(const State& state) const;

/** Return the current velocity-level constraint error for this %Constraint.
This is the relative velocity between the material points of body 1 and
body 2 that are coincident with the constrained station point on body 2;
note that this is subtly different from the time derivative of the
position error vector. The returned vector is measured in the Ancestor
body frame. The given \a state must be realized through Velocity stage. **/
Vec3 getVelocityErrors(const State& state) const;

/** Return the current acceleration-level constraint error for this
%Constraint. This is the relative acceleration between the material points
of body 1 and body 2 that are coincident with the constrained station point
on body 2; this is precisely the time derivative of the
velocity error vector (but not exactly the second time derivative of the
position error). The returned vector is measured in the Ancestor
body frame. The given \a state must be realized through Acceleration
stage. **/
Vec3 getAccelerationErrors(const State&) const;

/** Return the force currently being applied by this %Constraint to the
point of body 1 that is coincident in space with the constrained point on
body 2. The force vector is expressed in body 1's local frame. **/
Vec3 getBallReactionForceOnBody1(const State&) const;
/** Return the force currently being applied by this %Constraint to body 2,
at its constrained station point. The force vector is expressed in body 2's
local frame. **/
Vec3 getBallReactionForceOnBody2(const State&) const;

/** Return the three Lagrange multipliers associated with the three
accleration-level constraint equations generated by this %Constraint.
Although these are related to reaction forces, if that's what you're
interested in you should use getBallReactionForcesOnBody1() or
getBallReactionForceOnBody2() instead; the definition of the multipliers
is somewhat arbitrary and will not always be easy to interpret as forces.
The given \a state must be realized through Acceleration stage. **/
Vec3 getMultipliers(const State& state) const;

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ball, BallImpl, Constraint);
/** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_BALL_H_



