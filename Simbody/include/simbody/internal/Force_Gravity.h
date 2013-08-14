#ifndef SimTK_SIMBODY_FORCE_GRAVITY_H_
#define SimTK_SIMBODY_FORCE_GRAVITY_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-13 Stanford University and the Authors.        *
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
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force 
subclass Force::Gravity and is logically part of Force.h. The file assumes that
Force.h will have included all necessary declarations. **/

namespace SimTK {

/** This force element represents a uniform gravitational field applied to a
set of bodies. This can be used instead of the more limited 
Force::UniformGravity class if you want more control.

Force::Gravity is given a magnitude g in units of acceleration, and a "down" 
direction (unit vector) d, expressed in the Ground frame. For example, in SI 
units "standard gravity" g at the Earth's surface is defined to be exactly 
g=9.80665 m/s^2. If gravity acts in Ground's -Y direction, then its direction 
vector is d=[0,-1,0]; together this yields a gravity vector 
v=g*d=[0,-9.80665,0]. By default gravity will be applied to each MobilizedBody
in the given SimbodyMatterSubsystem (except for Ground). You can specify that
some bodies are excluded; those will not have a gravity force applied. You can
include multiple Force::Gravity elements with different exclusion lists in the
same System if you want different bodies to feel different effects.

\par Force
Each body B that has not been explicitly excluded will experience a force 
fb = mb*g*d, applied to its center of mass, where mb is the mass of body B. 
Although this is a pure force, note that when it is measured in the body frame 
B there will also be a moment unless the body frame origin Bo is located at the 
body's mass center. You can obtain the applied forces if you need them, for
example for gravity compensation; see getBodyForce().

\par Potential Energy
Gravitational potential energy for a body B is mb*g*hb where 
hb is the height of body B's mass center over an arbitrary "zero"
height hz (default is hz=0), measured along the "up" direction -d. If pb is 
the Ground frame vector giving the position of body B's mass center, its
height over or under hz is hb=pb*(-d) - hz. Note that this is a signed
quantity so the potential energy is also signed. **/
class SimTK_SIMBODY_EXPORT Force::Gravity : public Force {
public:


//------------------------------------------------------------------------------
/** @name                       Construction
Methods in this section refer both to constructors, and to methods that can
be used to set or change contruction (Topology-stage) parameters; these
specify the values assigned by default to the corresponding state variables. 
Note:
  - Changing one of these default parameters invalidates the containing 
    System's topology, meaning that realizeTopology() will have to be called 
    and a new State obtained before subsequent use. 
  - The set...() methods return a reference to "this" Gravity element (in
    the manner of an assignment operator) so they can be chained in a single 
    expression. **/
/*@{*/


/** This is the most general constructor for creating a Gravity force element 
within a particular force subsystem and affecting all the bodies of a particular
matter subsystem, given gravity magnitude and direction separately. 

If you want to exclude by default some bodies from the effects of this force 
element, call setDefaultBodyIsExcluded(). You can also make changes on the 
fly with setBodyIsExcluded().

@param[in,out]      forces       
    The subsystem to which this force should be added.
@param[in]          matter        
    The subsystem containing the bodies that will be affected.
@param[in]          down
    The default gravity direction vector d (i.e., the "down" direction),
    expressed in the Ground frame. Gravity will be directed along this direction
    unless explicitly changed within a particular State via the setDirection() 
    or setGravityVector() methods. This is a unit vector; if you provide an 
    ordinary Vec3 it will be normalized before use and its magnitude ignored.
    It is an error to provide a zero vector here; provide a zero magnitude
    instead or use the alternate constructor.
@param[in]          g
    The default magnitude g to be used for gravity, given as a nonnegative
    scalar with units of acceleration. The gravity vector is calculated as 
    v=g*d where d is the current direction vector, typically the default 
    direction as given by the \a down argument to this constructor. 
    Gravity will have the magnitude given here unless explicitly changed within
    a particular State via the setMagnitude() or setGravityVector() methods.
@param[in]          zeroHeight
    This is an optional specification of the default value for the height
    up the gravity vector that is considered to be "zero" for purposes of
    calculating the gravitational potential energy. The default is 
    \a zeroHeight == 0, i.e., a body's potential energy is defined to be zero
    when the height of its mass center is the same as the height of the Ground
    origin. The zero height will have the value specified here unless
    explicitly changed within a particular State use the setZeroHeight() 
    method. **/
Gravity(GeneralForceSubsystem&          forces, 
        const SimbodyMatterSubsystem&   matter,
        const UnitVec3&                 down,
        Real                            g,
        Real                            zeroHeight = 0);

/** Convenience constructor to create a %Gravity force element by specifying 
only a gravity vector, that is, the product of the gravity magnitude and the 
"down" direction vector. We have to extract both a magnitude and direction from 
this since they are maintained separately by this force element. If you provide 
a zero vector here (exactly zero, that is), then the magnitude is zero but there
is no direction. In that case we default to the negative of the containing 
System's default "up" direction. 

If you want to exclude by default some bodies from the effects of this force 
element, call setDefaultBodyIsExcluded(). You can also make changes on the 
fly with setBodyIsExcluded().

@param[in,out]      forces       
    The subsystem to which this force should be added.
@param[in]          matter        
    The subsystem containing the bodies that will be affected.
@param[in]          gravity
    The default gravity vector v, interpreted as v=g*d where g=|\a gravity| is 
    a positive scalar and d is the "down" direction unit vector d=\a gravity/g.
    If the magnitude is exactly zero we'll set the "down" direction to the 
    opposite of the containing System's "up" direction, otherwise the direction
    is extracted from the given vector.

This constructor sets the default zero height hz to zero (for use in calculating 
gravitational potential energy). If you want a different default value, use
the more general constructor or the setDefaultZeroHeight() method. 
@see SimTK::System::setUpDirection() **/
Gravity(GeneralForceSubsystem&          forces, 
        const SimbodyMatterSubsystem&   matter,
        const Vec3&                     gravity);

/** Convenience constructor to create a %Gravity force element by specifying 
only gravity's magnitude, with the direction chosen to oppose the containing 
System's "up" direction. Note that the direction is extracted from the System 
and recorded when this %Gravity element is contructed; it will not be affected
by subsequent changes to the System's up direction.

If you want to exclude by default some bodies from the effects of this force 
element, call setDefaultBodyIsExcluded(). You can also make changes on the 
fly with setBodyIsExcluded().

@param[in,out]      forces       
    The subsystem to which this force should be added.
@param[in]          matter        
    The subsystem containing the bodies that will be affected.
@param[in]          g
    The default magnitude g to be used for gravity, given as a nonnegative
    scalar with units of acceleration. %Gravity will have the magnitude given 
    here unless explicitly changed within a particular State via the 
    setMagnitude() or setGravityVector() methods.

This constructor sets the default zero height hz to zero (for use in calculating 
gravitational potential energy). If you want a different default value, use
the more general constructor or the setDefaultZeroHeight() method.  
@see SimTK::System::setUpDirection() **/
Gravity(GeneralForceSubsystem&          forces, 
        const SimbodyMatterSubsystem&   matter,
        Real                            g);
   
/** Default constructor creates an empty handle. **/
Gravity() {}

/** Set how the indicated body is to be treated by default. Normally all bodies
except Ground are affected by this force element so you only need to call this
if you want to exclude a body. This is a topological change since the flag will 
be stored with the System; for on-the-fly changes you should call
setBodyIsExcluded() instead. Ground is always excluded; you can call this 
method on Ground but the call will be ignored.  

@param[in]      mobod
    This is a MobilizedBody or MobilizedBodyIndex for the body whose 
    default exclusion status you want to change.
@param[in]      isExcluded
    Set to true if the default should be that this body is \e not to be affected
    by this %Gravity force element.
@return A writable reference to "this" Gravity element which will now have
    the indicated body excluded.

@see setBodyIsExcluded() **/
Gravity& setDefaultBodyIsExcluded(MobilizedBodyIndex mobod, bool isExcluded);

/** Set the default value for the gravity vector, that is, the direction and
magnitude with which gravity will act. From the given vector we extract 
separately the magnitude g and "down" direction d (a unit vector), unless the 
magnitude is exactly zero in which case we'll leave the direction unchanged.
@param[in]      gravity    
    A vector giving the magnitude and direction with which gravity will act
    on the bodies. If exactly zero only the magnitude will be changed.
@return A writable reference to "this" Gravity element which will now have
    the new default magnitude and direction. **/
Gravity& setDefaultGravityVector(const Vec3& gravity);

/** Set the default "down" direction d, that is, the direction along which 
gravity will act.
@param[in] down    
    A unit vector giving the direction along which gravity will act ("down").
    If you pass an ordinary Vec3, it will be normalized before use and its
    length will be ignored. Do not pass a zero-length Vec3, however.
@return A writable reference to "this" Gravity element which will now have
    the new default \a direction. **/
Gravity& setDefaultDownDirection(const UnitVec3& down);
/** Convenience overload that takes the down direction as a Vec3 and 
normalizes it (throwing away the magnitude) to create the required unit 
vector for the down direction. It is an error if the supplied Vec3 has zero 
length. **/
Gravity& setDefaultDownDirection(const Vec3& down)
{   return setDefaultDownDirection(UnitVec3(down)); }

/** Set the default magnitude of gravity (a nonegative scalar). This will be 
combined with the default \a direction unit vector d to calculate the gravity 
vector v = g*d that is used to generate body forces mb*v where mb is the
mass of body b and the force is applied to b's mass center.
@param[in]      g
    The magnitude of the uniform gravitational field, given as a nonnegative
    scalar in units of acceleration. If this is zero no forces will be applied
    by this element.
@return A writable reference to "this" Gravity element which will now have
    the new default \a strength. **/
Gravity& setDefaultMagnitude(Real g);

/** Set the default zero height (a signed scalar), for use in potential energy
calculation. See the Gravity class comments for how this is used.
@param[in]      zeroHeight
    The height that is considered to be "zero" for purposes of potential
    energy calculation. An affected body whose mass center is at this height
    will have zero gravitational potential energy.
@return A writable reference to "this" Gravity element which will now have
    the new default \a zeroHeight. **/
Gravity& setDefaultZeroHeight(Real zeroHeight);

/** Return the current setting of the "is excluded by default" flag for the 
given body. This is the status that the flag will have in the default State
returned by System::realizeTopology().
@param[in]      mobod
    This is a MobilizedBody or MobilizedBodyIndex for the body whose 
    exclusion status you want to obtain.
@return The current value of the "is excluded by default" flag for this body; 
    this will always be \c true for Ground.    
@see getBodyIsExcluded() **/
bool getDefaultBodyIsExcluded(MobilizedBodyIndex mobod) const;
/** Return the default gravity vector being used for this Gravity force 
element, calculated from the default magnitude and direction. **/
Vec3 getDefaultGravityVector() const;
/** Return the default down direction (a unit vector) for this Gravity force 
element. **/
const UnitVec3& getDefaultDownDirection() const;
/** Return the default gravity magnitude g (a nonnegative scalar). **/
Real getDefaultMagnitude() const;
/** Return the default zero height used in the calculation of graviational
potential energy. **/
Real getDefaultZeroHeight() const;

/*@}............................ Construction ................................*/



//------------------------------------------------------------------------------
/** @name                      Instance Parameters

These refer to the Instance-stage state variables that determine the actual
values to be used to calculate gravitational forces and energy from a given
State object, and which bodies are to be excluded from that calculation. If 
these are not set explicitly, the values are set to those provided in the 
constructor or via the correponding setDefault...() methods. Note:
  - Changing one of these parameters invalidates the given State's 
    Instance stage, meaning that realize(Instance) will have to be called 
    (explicitly or implicitly by realizing a higher stage) before 
    subsequent use.
  - The set...() methods here return a const reference to "this" Gravity force
    element (in the manner of an assignment operator, except read-only) so they
    can be chained in a single expression. **/
/*@{*/

/** Within a given State, selectively exclude (or include) a body from the 
effects of this Gravity force. The default State will inherit the "is excluded
by default" values set during construction, but then you can change them. Ground 
is always excluded; you can call this method on Ground but the call will
be ignored. 
@param[in,out]      state
    The State object that is modified by this method.
@param[in]          mobod
    This is a MobilizedBody or MobilizedBodyIndex for the body whose 
    exclusion status you want to change.
@param[in]          isExcluded
    Set to true if this body should \e not be affected by this force element.
@return A const reference to "this" Gravity force element for convenient 
    chaining of  set...() methods in a single expression.
@pre \a state must already be realized to Stage::Topology. 
@see setDefaultBodyIsExcluded(), getBodyIsExcluded() **/
const Gravity& setBodyIsExcluded(State& state, MobilizedBodyIndex mobod, 
                                bool isExcluded) const;

/** Set the gravity vector v, that is, the magnitude and direction with which 
gravitational forces will act, overriding the default magnitude and direction 
that are stored in this %Gravity force element. If the given vector is exactly
zero, then only the magnitude will be changed here.
@param[in,out]      state
    The State object that is modified by this method.
@param[in]          gravity    
    The gravity vector including both the magnitude and direction. If this is
    exactly zero only the magnitude will be changed.
@return A const reference to "this" Gravity element for convenient chaining of 
    set...() methods in a single expression.
@pre \a state must already be realized to Stage::Topology. 
@see setDownDirection(), setMagnitude() **/
const Gravity& setGravityVector(State& state, const Vec3& gravity) const;

/** Set the "down" direction d (a unit vector), that is, the direction along 
which gravitational forces will act. Magnitude can be changed independently
with setMagnitude().
@param[in,out]      state
    The State object that is modified by this method.
@param[in]          down    
    A unit vector giving the "down" direction along which gravity will act.
@return A const reference to "this" Gravity element for convenient chaining of 
    set...() methods in a single expression.
@pre \a state must already be realized to Stage::Topology. 
@see setMagnitude() **/
const Gravity& setDownDirection(State&           state,
                                const UnitVec3&  down) const;
/** Convenience overload that takes the down direction as a Vec3 and 
normalizes it (throwing away the magnitude) to create the required unit 
vector. It is an error if the supplied Vec3 has zero length. **/
const Gravity& setDownDirection(State&          state,
                                const Vec3&     down) const
{   return setDownDirection(state, UnitVec3(down)); }

/** Set the gravity magnitude g (a nonnegative scalar) in this \a state,
overriding the default gravity magnitude that is stored in this Gravity force
element.
@param[in,out] state    The State object that is modified by this method.
@param[in]     g        The magnitude of the uniform gravitational field.
@return A const reference to "this" Gravity element for convenient chaining of 
    set...() methods in a single expression.
@pre \a state must already be realized to Stage::Topology. 
@see setDownDirection() **/
const Gravity& setMagnitude(State& state, Real g) const;
/** Set the potential energy zero height hz (a scalar) in this \a state,
overriding the default zero height that is stored in this Gravity force
element.
@param[in,out]      state
    The State object that is modified by this method.
@param[in]          hz    
    The height at which potential energy is considered to be zero.
@return A const reference to "this" Gravity element for convenient chaining of 
    set...() methods in a single expression.
@pre \a state must already be realized to Stage::Topology. **/
const Gravity& setZeroHeight(State& state, Real hz) const;

/** Return the current setting of the "is excluded" flag for a given body in a 
given State. This is not necessarily the same as the default value for this 
body's exclusion status.
@param[in]      state
    The state containing the "is excluded" flag to be obtained.
@param[in]      mobod
    A MobilizedBody or MobilizedBodyIndex for the body whose exclusion status 
    you want to obtain.
@return The current value of the "is excluded" flag for this body in this 
    \a state; this will always be \c true for Ground. 
@see getDefaultBodyIsExcluded(), setBodyIsExcluded() **/
bool getBodyIsExcluded(const State& state, MobilizedBodyIndex mobod) const;
/** Get the gravity vector v that will be used for computations done with this
\a state. This is calculated as v=g*d from the current values of the magnitude
g and direction d Instance variables in the \a state.
@param[in]      state
    The state containing the gravity Instance variables used to calculate the
    gravity vector.
@return The current value of the gravity vector.
@pre \a state must have been realized to Topology stage. **/
Vec3 getGravityVector(const State& state) const;
/** Get the gravity "down" direction d (a unit vector) that is currently set 
in this \a state.
@param[in]      state
    The state containing the gravity direction Instance variable whose value is
    returned.
@return The current value of the gravity direction Instance variable.
@pre \a state must have been realized to Topology stage. **/
const UnitVec3& getDownDirection(const State& state) const;
/** Get the gravity magnitude g (a nonnegative scalar) that is currently set 
in this \a state.
@param[in]      state
    The state containing the gravity magnitude Instance variable whose value is
    returned.
@return The current value of the gravity magnitude Instance variable.
@pre \a state must have been realized to Topology stage. **/
Real getMagnitude(const State& state) const;
/** Get the zero height hz that is currently set in this \a state for use in
calculating gravitational potential energy. See the class documentation for 
information about hz is used.
@param[in]      state
    The state containing the zero height Instance variable whose value is
    returned.
@return The current value of the zero height Instance variable.
@pre \a state must have been realized to Topology stage. **/
Real getZeroHeight(const State& state) const;

/*@}........................ Instance Parameters .............................*/



//------------------------------------------------------------------------------
/** @name                       Energy and Forces
These methods return the gravitational potential energy and the forces being 
applied by this Gravity element in the configuration contained in the supplied 
State. These are evaluated during force calculation and available at no 
computational cost afterwards, although the first call after a position or 
instance-variable change will initiate computation if forces haven't already 
been computed. **/
/*@{*/

/** Obtain the gravitational potential energy contained in this Gravity force
element in the given configuration. Note that gravitational potential energy
is m*g*h for each body where the height h is measured from an arbitrary zero
height. By default the zero height is defined to be height of the Ground 
origin but that can be changed arbitrarily.
@param[in]          state    
    The State from whose cache the potential energy is retrieved.
@return The gravitational potential energy of this Gravity force element in the 
    configuration contained in \a state (a signed scalar).
@pre \a state must be realized to Stage::Position **/
Real getPotentialEnergy(const State& state) const;

/** Obtain a reference to the set of gravitational forces currently being 
applied by this %Gravity force element, as a Vector of spatial forces indexed
by MobilizedBodyIndex. The force on Ground or on any explicitly excluded body 
is zero, but there is an entry for every mobilized body (starting with Ground) 
in the returned result. The return vector of spatial vectors is in a form suitable for direct use with
the SimbodyMatterSubsystem Jacobian operators, which can be used to transform
these into generalized forces. This can be useful for gravity compensation.

@param[in]  state   
    The State from whose cache the force is retrieved.
@return A reference to the Vector of spatial forces (meaning the gravitational 
    force and moment about the body origin) currently being applied to each 
    of the mobilized bodies, expressed in the Ground frame.
@pre \a state must be realized to Stage::Position 

Forces are returned as though applied at the body origin, which is
\e not necessarily the same place as the body center of mass. That means that
in general there will be both a force and a moment returned for each body.

If gravity forces have not yet been calculated for the configuration given
in \a state, the computation will be initiated here and cached for use later.
**/
const Vector_<SpatialVec>& getBodyForces(const State& state) const;

/** Convenience method to extract the gravitational force on just one body;
see getBodyForces() to get the whole set at once, and for more information.

@param[in]  state   
    The State from whose cache the force is retrieved.
@param[in]  mobod   
    The MobilizedBody or MobilizedBodyIndex whose force is wanted.
@return The spatial force (meaning the gravitational force and moment about
    the body origin) currently being applied to the given mobilized body, 
    expressed in the Ground frame.
@pre \a state must be realized to Stage::Position **/
const SpatialVec& 
getBodyForce(const State& state, MobilizedBodyIndex mobod) const
{   return getBodyForces(state)[mobod]; }

// Particles aren't supported yet so don't show this in Doxygen.
/** @cond **/
/** Obtain the gravitational forces currently being applied by this Gravity
force element to the particles.
@param[in]  state       The State from whose cache the force is retrieved.
@return The current value of the forces, expressed in the Ground frame.
@pre \a state must be realized to Stage::Position **/
const Vector_<Vec3>& getParticleForces(const State& state) const;
/** @endcond **/

/*@}............................ Energy and Forces ...........................*/

//------------------------------------------------------------------------------
/** @name                     Debugging/Testing
This is information for use by developers of this class for debugging and
testing -- please ignore. **/
/*@{*/
/** Return a count of the number of times the set of gravitational forces
or potential energy was calculated since this force element was constructed. 
Note that this is not a per-body count. This count is only incremented when an 
actual computation is performed, not when forces are requested but are already 
valid. Also, we don't consider it a computation if the gravity magnitude is 
zero. This is intended for use in testing that caching and invalidation is being
done properly. **/
long long getNumEvaluations() const;

/** Return true if the gravitational forces for this configuration have already
been calculated and are up to date. That means they can be retrieved with no
further evaluation. **/
bool isForceCacheValid(const State& state) const;

/** Invalidate the stored gravitational forces if they have already been 
calculated at this configuration. That will force a new evaluation the next
time they are requested (unless the gravity magnitude is currently zero). **/
void invalidateForceCache(const State& state) const;
/*@}*/

// Don't show this in Doxygen.
/** @cond **/
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Gravity, GravityImpl, Force);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_GRAVITY_H_
