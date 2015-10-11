#ifndef SimTK_SIMBODY_CABLE_SPRING_H_
#define SimTK_SIMBODY_CABLE_SPRING_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-13 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/CablePath.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force
subclass CableSpring, providing an elastic band force element that follows a
given CablePath. **/

namespace SimTK {

/** This force element implements a passive elastic element (like a
rubber band) that follows a frictionless CablePath across a set of "obstacles".
The element calculates a uniform nonnegative tension that is used to apply
forces at the end points of the CablePath and to each intermediate obstacle.
The model provides stiffness and dissipation, and has a slack length below
which the tension is zero and no forces are generated. Dissipated power is
calculated and integrated so that work lost to dissipation is available for
conservation of energy calculations.

\par Theory:

Given current CablePath length L with time derivative Ldot, and slack length
L0 for this force element, define stretch x=max(0,L-L0) and xdot=Ldot. Then
calculate the nonnegative tension f(x,xdot), the potential energy pe(x) stored
in this force element, and dissipating power powerLoss(x,xdot) due to
rate-dependent resistance to length change, as follows:
<pre>
    f_stretch   = k*x
    f_rate      = max(-f_stretch, f_stretch*c*xdot)
    f           = f_stretch + f_rate
    pe          = k*x^2/2
    powerLoss   = f_rate * xdot
    dissipation = integ(powerLoss, dt)
</pre>
where \c k >= 0 is the stiffness coefficient (force per unit length) and
\c c >= 0 is the dissipation coefficient (1/velocity) for this force element.

Note that in this model the tension component \c f_rate, due to rate-dependent
dissipation in the elastic element, is calculated as a fraction of the
stretch-dependent tension component \c f_stretch, similar to a Hunt and Crossley
contact model. That is why the dissipation coefficient \c c has units of
1/velocity, rather than force/velocity. This makes the generated tension grow
smoothly from zero as the slack length is exceeded, but even so this is not
necessarily a good model for any particular physically-realizable elastic
element. A ligament, for example, while perhaps qualitatively similar to this
element, is likely to require a more complicated relationship between
\c x, \c xdot, and \c f to yield quantitative agreement with experimental data.

While the total tension \c f must be nonnegative, the tension \c f_rate due to
stretch rate can be positive (resists stretching) or negative (resists
shortening) but can never be less than \c -f_stretch since \c f can't be
negative. When a stretched cable spring is shortening so rapidly that force due
to dissipation cancels the force due to stiffness, there will be no tension
generated but we will still dissipate power at a rate that exactly accounts for
the continuing loss of potential energy in the spring (<code>k*x*xdot</code>)
as it shortens to its slack length. In that way total energy
<code>E=pe+ke+dissipation</code> is conserved, where \c dissipation is the
time integral of \c powerLoss calculated above (note that \c powerLoss >= 0
as defined).

@see CablePath, CableTrackerSubsystem **/
class SimTK_SIMBODY_EXPORT CableSpring : public Force {
public:

/** Create an elastic force element that follows a given CablePath and add it
to a GeneralForceSubsystem. Default values for the cable spring properties
are given here; you can override them in the State. See the %CableSpring class
description for a detailed explanation of these parameters.

@param[in,out]      forces
    The subsystem to which this force element should be added.
@param[in]          path
    The CablePath that defines the routing of this elastic element over
    geometric obstacles.
@param[in]          defaultStiffness
    A nonnegative spring constant representing the stiffness of this element,
    in units of force/length, where the force represents a uniform tension along
    the element that results from stretching it beyond its slack length.
@param[in]          defaultSlackLength
    The maximum length this elastic element can have before it begins to
    generate force. At or below this length the element is slack and has zero
    tension and zero power dissipation.
@param[in]          defaultDissipationCoef
    A nonnegative dissipation coefficient for this elastic element in units of
    1/velocity.
**/
CableSpring(GeneralForceSubsystem&  forces,
            const CablePath&        path,
            Real                    defaultStiffness,
            Real                    defaultSlackLength,
            Real                    defaultDissipationCoef);

/** Default constructor creates an empty handle. **/
CableSpring() {}


//------------------------------------------------------------------------------
/** @name                      Default Parameters

These refer to Topology-stage parameters normally set in the constructor; these
determine the initial values assigned to the corresponding Instance-stage state
variables.

\par Notes
- Changing one of these parameters invalidates the containing System's
  topology, meaning that realizeTopology() will have to be called
  before subsequent use.
- The set...() methods return a reference to "this" CableSpring (in
  the manner of an assignment operator) so they can be chained in a
  single expression. **/
/**@{**/

/** Set the stiffness (spring constant) k that will be used by default
for this cable spring; must be nonnegative.
@param[in]  stiffness       The spring constant k to be used by default.
@return A writable reference to "this" %CableSpring which will now have the new
default stiffness. **/
CableSpring& setDefaultStiffness(Real stiffness);
/** Set the slack length L0 that will be used by default for this cable spring;
must be nonnegative.
@param[in]  slackLength     The slack length L0 to be used by default.
@return A writable reference to "this" %CableSpring which will now have the new
default slack length. **/
CableSpring& setDefaultSlackLength(Real slackLength);
/** Set the dissipation coefficient c that will be used by default for this
cable spring; must be nonnegative.
@param[in]  dissipation     The dissipation coefficient c to be used by default.
@return A writable reference to "this" %CableSpring which will now have the new
default dissipation coefficient. **/
CableSpring& setDefaultDissipationCoef(Real dissipation);

/** Return the stiffnesses k that will be used by default for this cable spring.
@return The default spring constant. **/
Real getDefaultStiffness() const;
/** Return the slack length L0 that will be used by default for this cable spring.
@return The default slack length. **/
Real getDefaultSlackLength() const;
/** Return the dissipation coefficient c that will be used by default for this
cable spring.
@return The default dissipation coefficient. **/
Real getDefaultDissipationCoef() const;
/**@}**/



//------------------------------------------------------------------------------
/**@name                     Instance Parameters

These refer to the Instance-stage state variables that determine the geometry
and material properties that will be used in computations involving this cable
spring when performed with the given State. If these are not set explicitly,
the default values are set to those provided in the constructor or via the
corresponding setDefault...() methods.

\par Notes
- Changing one of these parameters invalidates the given State's Instance stage,
  meaning that realize(Instance) will have to be called (explicitly or
  implicitly by realizing a higher stage) before subsequent use.
- The set...() methods here return a const reference to "this" %CableSpring (in
  the manner of an assignment operator, except read-only) so they can be chained
  in a single expression.
**/
/**@{**/
/** Set the stiffness (spring constant) k that will be used for this cable
spring when evaluated using this State.
@pre \a state realized to Stage::Topology
@param[in,out]  state       The State object to be modified by this method.
@param[in]      stiffness   The spring constant k.
@return A const reference to "this" %CableSpring for convenient chaining of
set...() methods in a single expression. **/
const CableSpring& setStiffness(State&      state,
                                Real        stiffness) const;
/** Set the slack length L0 that will be used for this cable spring when
evaluated using this State.
@pre \a state realized to Stage::Topology
@param[in,out]  state       The State object to be modified by this method.
@param[in]      slackLength The slack length L0.
@return A const reference to "this" %CableSpring for convenient chaining of
set...() methods in a single expression. **/
const CableSpring& setSlackLength(State&    state,
                                  Real      slackLength) const;
/** Set the dissipation coefficient c that will be used for this cable spring
when evaluated using this State.
@pre \a state realized to Stage::Topology
@param[in,out]  state       The State object that is modified by this method.
@param[in]      dissipationCoef The dissipation coefficient c.
@return A const reference to "this" %CableSpring for convenient chaining of
set...() methods in a single expression. **/
const CableSpring& setDissipationCoef(State&    state,
                                      Real      dissipationCoef) const;

/** Return the stiffness (spring constant) k currently being used for this cable
spring by this State.
@pre \a state realized to Stage::Topology
@param[in]  state   The State object from which to obtain the stiffness.
@return The current value in the given state of the spring constant k. **/
Real getStiffness(const State& state) const;
/** Return the slack length L0 currently being used for this cable spring by
this State.
@pre \a state realized to Stage::Topology
@param[in]  state   The State object from which to obtain the slack length.
@return The current value in the given state of the slack length L0. **/
Real getSlackLength(const State& state) const;
/** Return the dissipation coefficient c currently being used for this cable
spring by this State.
@pre \a state realized to Stage::Topology
@param[in]  state   The State object from which to obtain the dissipation
                    coefficient.
@return The current value in the given state of the dissipation
coefficient c. **/
Real getDissipationCoef(const State& state)   const;
/**@}**/



//------------------------------------------------------------------------------
/** @name               Position-related Quantities

These methods return position-dependent quantities that are calculated
by the cable spring as part of its force calculations and stored in the
State cache. These can be obtained at no cost, although the first call
after a position change may initiate computation if forces haven't
already been computed.

@pre These methods may be called only after the supplied state has
     been realized to Stage::Position. **/
/**@{**/
/** Return the current length of the CablePath that underlies this cable spring
element. Note that this may return a value less than the spring's slack length.
@pre \a state realized to Stage::Position
@param  state   The State from which the path length is retrieved.
@return The current length of the underlying CablePath (nonnegative).
@see getLengthDot() **/
Real getLength(const State& state) const;
/**@}**/



//------------------------------------------------------------------------------
/** @name                 Velocity-related Quantities

These methods return velocity-dependent quantities that are calculated
by the cable spring as part of its force calculations and stored in the
State cache. These can be obtained at no cost, although the first call
after a velocity change may initiate computation if forces haven't
already been computed.

@pre These methods may be called only after the supplied state has
     been realized to Stage::Velocity. **/
/**@{**/
/** Return the current rate of length change (time derivative of length) of the
CablePath that underlies this cable spring element.
@pre \a state realized to Stage::Velocity
@param  state   The State from which the path length rate of change is
                retrieved.
@return The current length change rate of the underlying CablePath.
@see getLength() **/
Real getLengthDot(const State& state) const;

/**@}**/



//------------------------------------------------------------------------------
/** @name                          Forces

These methods return the forces being applied by this cable spring in the
configuration and velocities contained in the supplied State. These are
evaluated during force calculation and available at no computational cost
afterwards, although the first call after a velocity change may initiate
computation if forces haven't already been computed.

@pre These methods may be called only after the supplied state has been
     realized to Stage::Velocity. **/
/**@{**/
/** Return the current level of tension in the cable spring.
@pre \a state realized to Stage::Velocity
@param[in]          state
    The State from which the current tension is retrieved.
@return
    The current tension in the cable spring in the configuration and velocity
    contained in \a state (a nonnegative scalar). **/
Real getTension(const State& state) const;
/**@}**/



//------------------------------------------------------------------------------
/** @name                   Energy, Power, and Work

These methods return the energy, power, and work-related quantities associated
with this %CableSpring element for the values in the supplied State. **/
/**@{**/
/** Obtain the potential energy stored in this cable spring in the current
configuration.
@pre \a state realized to Stage::Position
@param[in]          state
    The State from whose cache the potential energy is retrieved.
@return
    The potential energy currently contained in the cable spring in the
    configuration contained in \a state (a nonnegative scalar). **/
Real getPotentialEnergy(const State& state) const;

/** Obtain the rate at which energy is being dissipated by this cable spring,
that is, the power being lost (presumably due to heat). This is in units of
energy/time which is watts in MKS.
@pre \a state realized to Stage::Velocity
@param[in]          state
    The State from which to obtain the current value of the power
    dissipation.
@return
    The instantaneous power dissipation (a nonnegative scalar).
@see getDissipatedEnergy() for the time-integrated power loss **/
Real getPowerDissipation(const State& state) const;

/** Obtain the total amount of energy dissipated by this cable spring since some
arbitrary starting point. This is the time integral of the power dissipation.
For a system whose only non-conservative forces are cable springs, the sum of
potential, kinetic, and dissipated energies should be conserved. This is a
State variable so you can obtain its value any time after it is allocated.
@pre \a state realized to Stage::Model
@param[in]          state
    The State from which to obtain the current value of the dissipated
    energy.
@return
    The total dissipated energy (a nonnegative scalar).
@see getPowerDissipation() for the instantaneous power loss **/
Real getDissipatedEnergy(const State& state) const;

/** Set the accumulated dissipated energy to an arbitrary value. Typically
this is used only to reset the dissipated energy to zero, but non-zero values
can be useful if you are trying to match some existing data or continuing a
simulation. This is a State variable so you can set its value any time after it
is allocated.
@pre \a state realized to Stage::Model
@param[in,out]      state
    The State whose dissipated energy variable for this cable spring is to
    be modified.
@param[in]          energy
    The new value for the accumulated dissipated energy (must be a nonnegative
    scalar). **/
void setDissipatedEnergy(State& state, Real energy) const;
/**@}**/

/**@name                Advanced/obscure/debugging
Miscellaneous methods that you probably aren't going to need. **/
/**@{**/
/** Get a reference to the CablePath that is used to determine the routing of
this cable spring force element. This is set during %CableSpring construction
and cannot be changed subsequently. **/
const CablePath& getCablePath() const;
/**@}**/

/** @cond **/ // Don't show this in Doxygen.
class Impl;
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (CableSpring, CableSpring::Impl, Force);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CABLE_SPRING_H_
