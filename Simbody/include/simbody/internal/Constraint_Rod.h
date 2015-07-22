#ifndef SimTK_SIMBODY_CONSTRAINT_ROD_H_
#define SimTK_SIMBODY_CONSTRAINT_ROD_H_

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
Declares the Constraint::Rod class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                  ROD (CONSTANT DISTANCE) CONSTRAINT
//==============================================================================

/** This constraint consists of one constraint equation that enforces a constant
distance between a point on one body and a point on another body. This is
like connecting them by a rigid, massless rod with ball joints at either end.
The constraint is enforced by a force acting along the rod with opposite
signs at either end. When positive, this represents tension in the rod
pulling the points together; when negative it represents compression keeping
the points separated.

@warning
You can't use this to enforce a distance of zero between two points.
That takes three constraints because there is no restriction on the force
direction. For a distance of zero (i.e., you want the points to be
coincident) use a Ball constraint, a.k.a. CoincidentPoints constraint.
**/
class SimTK_SIMBODY_EXPORT Constraint::Rod : public Constraint {
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
  - The set...() methods return a reference to "this" %Rod constraint
    element (in the manner of an assignment operator) so they can be chained in
    a single expression.

The default parameters can be overridden in any given State, and modified
without affecting Topology; see the "Runtime Changes" section below. **/
/*@{*/
/** Construct a %Rod (constant distance) constraint as described in the
Constraint::Rod class documentation.

@param      mobod1
    The first MobilizedBody object to which one of the rod endpoints is fixed.
@param      defaultPointOnBody1
    This is the location of a point P on \a mobod1, given as a vector p_B1P
    from the B1 frame origin to the point P, expressed in the \a mobod1 body
    frame B1. This is the point location that will be present in a default
    State; you can modify it later.
@param      mobod2
    The second MobilizedBody object to which the other rod endpoint is fixed.
@param      defaultPointOnBody2
    This is the location of a point Q on \a mobod2, given as a vector p_B2Q
    from the B2 frame origin to the point Q, expressed in the \a mobod2 body
    frame B2. This is the point location that will be present in a default
    State; you can modify it later.
@param      defaultRodLength
    The rod length (required distance between points P and Q). This is the value
    that will be present in a default State; you can modify it later.
**/
Rod(MobilizedBody& mobod1, const Vec3& defaultPointOnBody1,
    MobilizedBody& mobod2, const Vec3& defaultPointOnBody2,
    Real defaultRodLength);

/** Construct a %Rod (constant distance) constraint as described in the
Constraint::Rod class documentation, but using the body origins as the %Rod
end points.

@param      mobod1
    The first MobilizedBody object whose origin is used as one of the end
    points. This is the same as specifying %Vec3(0) for \a defaultPointOnBody1
    in the preceding constructor. This is the value that will be present in a
    default State; you can modify it later.
@param      mobod2
    The second MobilizedBody object whose origin is used as one of the end
    points. This is the same as specifying %Vec3(0) for \a defaultPointOnBody2
    in the preceding constructor. This is the value that will be present in a
    default State; you can modify it later.
@param      defaultRodLength
    The rod length (required distance between points P and Q). This is the value
    that will be present in a default State; you can modify it later.
**/
Rod(MobilizedBody& mobod1, MobilizedBody& mobod2,
    Real defaultRodLength);

/** Default constructor creates an empty handle that can be used to
reference any %Rod %Constraint.**/
Rod() {}

/** Return a reference to the first MobilizedBody given in the
constructor. This cannot be changed after construction. **/
const MobilizedBody& getMobilizedBody1() const;

/** Return a reference to the second MobilizedBody given in the
constructor. This cannot be changed after construction. **/
const MobilizedBody& getMobilizedBody2() const;

/** Replace the default location for the point attached to the first body,
\a mobod1, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method.
@see setPointOnBody1() **/
Rod& setDefaultPointOnBody1(const Vec3& defaultPoint);
/** Replace the default location for the point attached to the second body,
\a mobod2, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method.
@see setPointOnBody2() **/
Rod& setDefaultPointOnBody2(const Vec3& defaultPoint);
/** Replace the default rod length (distance), replacing whatever value was
supplied on construction. This is a topological change; you'll have to call
realizeTopology() again if you call this method.
@see setRodLength() **/
Rod& setDefaultRodLength(Real defaultRodLength);

/** Return the default location for the point attached to the first mobilized
body, \a mobod1, as set during construction or by the most recent call to
setDefaultPointOnBody1().
@see getPointOnBody1() **/
const Vec3& getDefaultPointOnBody1() const;
/** Return the default location for the point attached to the second mobilized
body, \a mobod2, as set during construction or by the most recent call to
setDefaultPointOnBody2().
@see getPointOnBody2() **/
const Vec3& getDefaultPointOnBody2() const;
/** Return the default rod length (distance) as set during construction or by
the most recent call to setDefaultRodLength().
@see getRodLength() **/
Real getDefaultRodLength() const;
/*@}............................ Construction ................................*/

//------------------------------------------------------------------------------
/** @name                       Visualization
If you are allowing the SimbodyMatterSubsystem to generate default geometry,
this constraint element will attempt to represent the distance constraint by
drawing a Green point representation at the point on body 1, an Orange point
representation at the point on body 2, and a Gray line along the point-to-point
line segment, beginning at point 1 but with length restricted to the rod length
(distance). Methods here give you some limited control over the generated
default geometry; if you need more control you should disable the generation of
default geometry and make your own. Methods are available to ask the constraint
element for its geometry for that purpose. **/
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
/** Modify the location of the point P on the first body, \a mobod1, in this
\a state by providing a new vector p_B1P giving the location of P measured from
body 1's origin, and expressed in the B1 frame. This overrides the
\a defaultPointOnBody1 in the given \a state, whose Stage::Position
is invalidated.
@see getPointOnBody1(), setDefaultPointOnBody1() **/
const Rod& setPointOnBody1(State& state, const Vec3& p_B1P) const;
/** Modify the location of the point Q on the second body, \a mobod2, in this
\a state by providing a new vector p_B2Q giving the location of Q measured from
body 2's origin, and expressed in the B2 frame. This overrides the
\a defaultPointOnBody2 in the given \a state, whose Stage::Position
is invalidated.
@see getPointOnBody2(), setDefaultPointOnBody2() **/
const Rod& setPointOnBody2(State& state, const Vec3& p_B2Q) const;
/** Modify the rod length (distance) in this \a state. This overrides the
\a defaultRodLength in the given \a state, whose Stage::Position is
invalidated.
@see getRodLength(), setDefaultRodLength() **/
const Rod& setRodLength(State& state, Real rodLength) const;
/** Return the position on body 1 of the rod end point P on the
first body, \a mobod1, as currently set in the given \a state. The value
is returned as a vector p_B1P from body 1's origin to the point P,
expressed in the B1 frame.
@see setPointOnBody1(), getDefaultPointOnBody1() **/
const Vec3& getPointOnBody1(const State& state) const;
/** Return the position on body 2 of the rod end point Q attached to the
second body, \a mobod2, as currently set in the given \a state. The value
is returned as a vector p_B2Q from body 2's origin to the point Q,
expressed in the B2 frame.
@see setPointOnBody2(), getDefaultPointOnBody2() **/
const Vec3& getPointOnBody2(const State& state) const;
/** Return the rod length (distance) as currently set in the given \a state.
@see setRodLength(), getDefaultRodLength() **/
Real getRodLength(const State& state) const;
/*@}.......................... Runtime Changes ...............................*/

//------------------------------------------------------------------------------
/** @name                      Computations
Methods here provide access to values already calculated by this constraint
element, and provide operators you can call to calculate related values. **/
/*@{*/
/** This is the signed violation of the position constraint, in length units.
It is positive when the two points are separated by more than the required
distance; negative when closer than the required distance. The given \a state
must have already been realized through Stage::Position. **/
Real getPositionError(const State&) const;

/** This is the time derivative of the value returned by getPositionError(); in
this case it is the relative velocity of the two points projected onto the
direction of the line between them. **/
Real getVelocityError(const State&) const;

/** This is the time derivative of the value returned by getVelocityError(). **/
Real getAccelerationError(const State&) const;

/** This is the Lagrange multiplier required to enforce the constraint
equation generated here. For this %Constraint it has units of force, but recall
that the sign convention for multipliers is opposite that of forces. Use the
other provided methods if you want to get meaningful forces. The given \a
state must already be realized to Stage::Acceleration. **/
Real getMultiplier(const State&) const;

/** This returns the tension in the %Rod being used to enforce the constraint.
It is positive if the %Rod is in tension, negative if in compression. The
result is zero if the constraint is currently disabled.  The given \a state
must already be realized to Stage::Acceleration. **/
Real getRodTension(const State&) const;

/** Return the instantaneous orientation of the %Rod in the Ground frame. This
is a unit vector along the direction from point1 (on body 1) to point2 (on
body 2), unless the points overlap in which case it is an arbitrary direction.
This method calculates a valid value even if this constraint is currently
disabled. The given \a state must already be realized to Stage::Position. **/
UnitVec3 findRodOrientationInG(const State& state) const;

/** Calculate the amount by which this constraint is violated.
It is positive when the distance between the points is greater than the
specified %Rod length, negative when less. If the constraint is currently
enabled then this returns the same value as getPositionError(); however, this
method calculates a valid value even if the constraint is currently disabled.
The given \a state must be realized to Stage::Position. **/
Real findLengthViolation(const State& state) const;

/*@}............................ Computations ................................*/

//------------------------------------------------------------------------------
/** @name               Advanced/Obscure/Obsolete/Debugging
You probably don't want to use the methods in this section. **/
/*@{*/
/** (Obscure) Use getMobilizedBody1() instead. **/
MobilizedBodyIndex getBody1MobilizedBodyIndex() const;
/** (Obscure) Use getMobilizedBody2() instead. **/
MobilizedBodyIndex getBody2MobilizedBodyIndex() const;
/*@}................... Advanced/Obscure/Obsolete/Debugging ..................*/

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Rod, RodImpl, Constraint);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_ROD_H_



