#ifndef SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_
#define SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021-22 Authors.
 Authors: Frank C. Anderson
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ----------------------------------------------------------------------------*/

#include "SimTKmath.h"
#include "simbody/internal/ExponentialSpringParameters.h"
#include "simbody/internal/Force.h"

namespace SimTK {

class MultibodySystem;
class MultibodySubsystem;
class MobilizedBody;
class ExponentialSpringForceImpl;

//=============================================================================
/** Class ExponentialSpringForce uses an "exponential spring" as a means
of modeling contact of a specified point on a MobilizedBody with a contact
plane that is fixed to Ground. In this documentation, this specified point
is referred to as the "body station". Each ExponentialSpringForce instance
acts at only one body station. In practice, you should choose a number of
body stations strategically located across the surface of a MobilizedBody,
and construct an ExponentialSpringForce for each of those body stations. For
example, if the body were a cube, you would likely choose the body stations
to be the corners of the cube and construct an ExponentialSpringForce
instance for each corner of the cube (so a total of 8 instances). The contact
plane is typically used to model interactions with a floor, but need not be
limited to this use case. The contact plane can be rotated and displaced
relative to the Ground frame and so can be used to model a wall, ramp, or
some other planar structure.

A distinguishing feature of the exponential spring, relative to other
contact models, is that it ALWAYS applies a force to the body; there is
never a time in a simulation when the spring force is not applied. This
seemingly non-physical feature works because (assuming default parameters are
in use) the force becomes small (less than 0.01 N) as soon as the body station
is more than 1.0 cm above the contact plane and extremely small (less than
0.000001 N) as soon as the body station is more than 2.0 cm above the contact
plane. The sizes of these forces are likely small compared to the errors made
in neglecting air resistance, assuming a uniform gravitational field, or
estimating inertial parameters.

As a side note, electrostatic forces, which are the fundamental force from
which all conventional contact forces arise, are never actually zero either.
They obey Coulomb's Law (Fₑ = kₑq₁q₂/r²) which states that an
electrostatic force only approaches zero as the distance between two
charges (r) approaches infinity. It is certainly not the case that exponential
springs accurately model contact force at the level of electrostatic forces;
we only suggest that the use of a contact force field that acts over a great
distance is reasonable as long as the force gets sufficiently small
sufficiently quickly.

Introducing any additional modeling approximation into a simulation is hard to
justify, however, unless there are significant benefits. In the case of
exponential springs there are two: simplicity and speed. With exponential
springs, there is no need to search for body intersections (the number of
active springs does not change during a simulation) and there is no need to
find the precise time of contact. In addition, because the function
describing the contact force is smooth, the integration step sizes taken
by variable-step explicit integrators are often well-behaved.

One of the challenges in modeling friction is handling the static
case. When the net force applied to a body in a direction tangent to a contact
plane is less than the limiting value of the static friction force (μₛ Fₙ),
the body should not drift but remain fixed in place. One approach to handling
this case is to apply conditional constraints that enforce zero acceleration
of the body in the contact plane. Unfortunately, implementation of this
approach is complex, particularly if there are multiple contact points on the
body. Another approach is to apply a velocity-dependent frictional force
that opposes any slide velocity (i.e., f = -cv). This latter approach is
relatively simple to implement, but will not entirely eliminate drift. By
making c very large, the drift velocity can be made small but it cannot be
brought to zero. And, unfortunately, as c is made large, the system equations
can become stiff, requiring integrator step size to be reduced, sometimes
considerably (at least in explicit integrators).

Part of the speed-up offered by ExponentialSpringForce derives, not from the
use of an exponential for the normal force, but from the friction model
employed. ExponentialSpringForce uses a spring-based frictional model that
includes a velocity-dependent (damping) term AND a position-dependent
(elastic) term. By including the elastic term, drift is entirely eliminated
and relatively large integration step sizes are maintained.

In initial comparisons, using class ExponentialSpringForce to model
contact resulted in simulation cpu times that were typically 4 times faster
(and sometimes 100+ times faster) than when using class
CompliantContactSubsystem, and yet the simulated motions were similar.
These comparisons can be reproduced by building and running the Test_Adhoc -
ExponentialSpringsComparison project that is included with Simbody.

Aspects of the exponential spring contact model are described in the following
publication:

        Anderson F.C. and Pandy M.G. (1999). A dynamics optimization
        solution for vertical jumping in three dimensions. Computer Methods
        in Biomechanics and Biomedical Engineering 2(3):201-231.

The current class makes several improvements to that contact model, most
notably including 1) the ability to rotate and translate the contact plane and
2) the ability to specify both a static and a kinetic coefficient of friction. 
The computational details of the contact model implemented by this class
follow below.

----------------------------------
Computations and Coordinate Frames
----------------------------------
The positive z-axis of the contact plane defines its normal. The positive
z-axis is the axis along which the repelling normal force (modeled using an
exponential) is applied. The x-axis and y-axis of the contact plane together
define the tangent plane. Friction forces will always be tangent to the x-y
plane.

In the equations below, a variable with a "z" suffix (e.g., pz, vz, or cz)
refers to a quantity that is normal to the contact plane or that pertains
to calculation of the normal force. A variable with an "xy" suffix
(e.g., pxy, vxy, or cxy) refers to a quantity that lies in or is tangent to
the contact plane or that pertains to calculation of the friction force.

### Normal Force (positive z-axis)

The elastic part of the normal force is computed using an exponential
whose shape is a function of three parameters (d₀, d₁, and d₂):

        fzElastic = d₁exp(−d₂(pz−d₀))

Note that pz is the displacement of the body station above (pz > 0.0)
or below (pz < 0.0) the contact plane. The default values of the shape
parameters were chosen to maximize integration step size while maintaining a
number of constraints (e.g., the normal force should fall below 0.01 Newtons
when pz > 1.0 cm). The damping part of the normal force is linear in velocity
and scaled by the elastic part:

        fzDamping = −cz vz fzElastic,

where vz is the normal component of the velocity of the body station and
cz is the damping coefficient for the normal direction. All together, the
spring force in the normal direction is given by

        fz  = fzElastic + fzDamping
            = d₁exp(d₂(py−d₀)) − cz vz d₁exp(d₂(pz−d₀)))
            = d₁exp(d₂(pz−d₀)) (1 − cz vz)

which has the form of the Hunt & Crossley damping model:

        K. H. Hunt and F. R. E. Crossley (1975). Coefficient of Restitution
        Interpreted as Damping in Vibroimpact. ASME Journal of Applied
        Mechanics, pp. 440-445.

### Friction Force (x-y plane)

The friction force is computed by blending two different friction models.
The blending is performed based on the 'Sliding' State of the
ExponentialSpringForce class. 'Sliding' is a continuous state variable (a Z in
Simbody vocabulary) that characterizes the extent to which either static or
kinetic conditions are present.

        Sliding = 0.0     static- fully fixed in place (lower bound)
        Sliding = 1.0     kinetic- fully sliding (upper bound)

More details about the Sliding State are given in sections below.

#### Friction Model 1 - Pure Damping (Sliding = 1.0)
When the body station is moving with respect to the contact plane, the
friction force is computed using a simple damping term:

        fricDamp = −cxy vxy

where cxy is the damping coefficient in the contact plane and vxy is the
velocity of the body station in the contact plane expressed in the contact
plane. However, the magnitude of the total frictional force is not allowed to
exceed the frictional limit:

        fricLimit = μ fz
        if (|fricDamp| > fricLimit)
            fricDamp = −fricLimit vxy / |vxy| = −μ fz vxy / |vxy|

where μ is the instantaneous coefficient of friction (more below). Note that
fz is always positive and so fricLimit is a positive scalar. Thus, for
velocities in the contact plane above some threshold velocity, which is
typically small (i.e., less than 0.1 m/s), this model is consistent with a
standard Coulomb Friction model.

#### Friction Model 2 - Damped Linear Spring (Sliding = 0.0)
When the body station is anchored with respect to the contact plane, the
friction force is represented by a damped linear spring. The viscous term is
given by the same damping expression as above:

        fricDampSpr = −cxy vxy

and the elastic term is given by

        fricElasSpr = −kxy (pxy−p₀)

where kxy is the friction spring elasticity, pxy is the position of the body
station projected onto the contact plane, and p₀ is the current spring zero
(i.e., the elastic anchor point of the friction spring). Note that p₀ always
resides in the contact plane.

The total friction spring force is then given by the sum of the elastic and
viscous terms:

        fricSpr = fricElasSpr + fricDampSpr

If the magnitude of the fricSpr exceeds the magnitude of the friction limit,
the terms are scaled down:

        if(|fricSpr| > fricLimit)
            scaleFactor = fricLimit / |fricSpr|
            fricDampSpr = scaleFactor * fricDampSpr
            fricElasSpr = scaleFactor * fricElasSpr
            fricSpr = fricElasSpr + fricDampSpr

Note that scaling down the friction spring force does not alter its direction.

#### Blending the Friction Models
Blending Model 1 and Model 2 is accomplished using linear expressions of the
Sliding State:

        fricElasBlend = fricElasSpr * (1.0 − Sliding)
        fricDampBlend = fricDampSpr + (fricDamp − fricDampSpr)*Sliding
        fricBlend = fricElasBlend + fricDampBlend

Regarding the elastic term, when Sliding = 0.0, fricElasBlend is given
entirely by fricElasSpr, and, as Sliding → 1.0, fricElasBlend → 0.0.
Regarding the damping term, when Sliding = 0.0, fricDampBlend is
given entirely by fricDampSpr, and, as Sliding → 1.0, fricDampBlend →
fricDamp. Any difference between fricDampSpr and fricDamp is due to the
difference in the ways the friction limit is enforced.

Thus, Model 1 (Pure Damping) dominates as Sliding → 1.0, and
Model 2 (Damped Linear Spring) dominates as Sliding → 0.0. The blending is
well behaved and smooth for all values of Sliding between 0.0 and 1.0.

#### Moving the Friction Spring Zero
The friction spring zero (p₀) (the elastic anchor point) is always made
to be consistent with the final value of the blended elastic force
(fricElasBlend):

        p₀ = pxy + fricElasBlend / kxy;
        p₀[2] = 0.0;  // ensures that p₀ lies in the contact plane

When fricElasBlend = 0.0, notice that p₀ = pxy, p₀[2] = 0.0. In other words,
when Sliding = 1.0, p₀ is simply the position of the body station projected
onto the contact plane.

In Simbody, p₀ is handled as an Auto Update Discrete %State. See
State::allocateAutoUpdateDiscreteVariable() for a detailed description. Any
change to p₀ is made to the Update Cache (not to the State directly), and the
integrator moves this cache value to the actual p₀ State after a successful
integration step is obtained.

#### Coefficients of Friction and SlidingDot
Coefficients of kinetic (sliding) and static (fixed) friction can be specified
separately for the spring, subject to the following constraints:

        0.0 ≤ μₖ ≤ μₛ

Note that there is no upper bound on μₛ. The instantaneous coefficient of
friction (μ) is calculated based on the value of the Sliding State:

        μ = μₛ − Sliding*(μₛ − μₖ)

The time derivative of Sliding, SlidingDot, is used to drive Sliding toward
the extremes of 0.0 or 1.0, depending on the following criteria:

- If the frictional spring force (fricSpr) exceeded the frictional limit at
any point during its calculation, Sliding is driven toward 1.0 (rise):

        SlidingDot = (1.0 − Sliding)/tau

- If the frictional spring force (fricSpr) does not exceed the frictional
limit at any point during its calculation and if the kinematics of the body
station are near static equilibrium, Sliding is driven toward 0.0 (decay):

        SlidingDot = −Sliding/tau

The threshold for being "near" static equilibrium is established by two
parameters, vSettle and aSettle. When the velocity and acceleration of the
body station relative to the contact plane are below vSettle and aSettle,
respectively, static equilibrium is considered effectively reached.

During a simulation, once a rise or decay is triggered, the rise or decay
continues uninterrupted until Sliding crosses 0.95 or 0.05, respectively. Once
these thresholds have been crossed, the criteria for rise and decay are again
monitored.

In the above equations for SlidingDot, tau is the characteristic time it takes
for the Sliding State to rise or decay. The default value of tau is 0.01 sec.
The motivation for using a continuous state variable is that, although the
transition between static and kinetic may happen quickly, it does not happen
instantaneously. Evolving Sliding based on a differential equation ensures
that μ is continuous and that the blending of friction models is well behaved.

### Future Enhancements

Enhancements might include the capacity to
- Fix the contact plane to a body other than Ground,
- Use a polygonal mesh as the contact surface, and
- Move the body station in manner that adapts to contact circumstances. For
example, to model a coin rolling on a table a small number of body stations
could be continually moved to the portion of the coin that is closest to the
table.

------
STATES
------
Each instance of ExponentialSpringForce possesses 5 states, which are listed
below in the appropriate category:

### DISCRETE STATES (parameters)
- μₛ (Real) = Static coefficient of friction.  0.0 ≤ μₛ
- μₖ (Real) = Kinetic coefficient of friction.  0.0 ≤ μₖ ≤ μₛ

As discrete states, μₛ and μₖ can be changed during the course of a simulation
without invalidating the System topology. This feature allows μₛ and μₖ to be
altered during a simulation to model, for example, a slippery spot on the
floor.

### AUTO UPDATE DISCRETE STATES
- p₀ (Vec3) = Zero point  of the frictional spring. p₀ always lies in the
contact plane.
- SlidingAction (int) = Action to take when setting SlidingDot. Possible
actions are Check, Rise, or Decay.

### CONTINUOUS STATE
- Sliding (Real) = Indicator of whether the spring zero (p₀) is moving or
fixed in the contact plane.  0.0 ≤ Sliding ≤ 1.0.
The "Sliding" state is used to transition the instantaneous coefficient of
friction (μ) between μₖ and μₛ. A value of 0.0 indicates that p₀ is fixed in
place, in which case μ = μₛ. A value of 1.0 indicates that p₀ is sliding, in
which case μ = μₖ. A value between 0.0 and 1.0 indicates that a transition
from sliding to fixed or from fixed to sliding is underway, in which case
μₖ ≤ μ ≤ μₛ. Sliding is also used to blend between friction Model 1 and Model 2
(see above for details).

----------
PARAMETERS
----------
Customizable Topology-stage parameters specifying the characteristics of the
exponential spring are managed using ExponentialSpringParameters.

----
DATA
----
Calculated quantities are cached when the System is realized through
Stage::Dynamics. The cached data is made accessible via conventional "get"
methods (see below). */
class SimTK_SIMBODY_EXPORT ExponentialSpringForce : public Force {
public:
    /** Construct an exponential spring force object with customized
    parameters.
    @param forces The general force subsystem that encapsulates many (and
    possibly all) of the body forces, torques, and generalized forces applied
    to the system during a simulation.
    @param contactPlane Transform specifying the location and orientation of
    the contact plane with respect to the Ground frame. The positive z-axis
    of the contact plane defines the normal direction; the x- and y-axes
    define the tangent (or friction) plane.
    @param body MobilizedBody that will interact / collide with the contact
    plane.
    @param station Point on the specified body at which the contact force
    will be applied. The position and velocity of this point relative to
    the contact plane determine the magnitude and direction of the contact
    force.
    @param params Customized Topology-stage parameters. If params is omitted
    from the argument list, the default set of parameters is used. */
    ExponentialSpringForce(GeneralForceSubsystem& forces,
        const Transform& contactPlane,
        const MobilizedBody& body, const Vec3& station,
        ExponentialSpringParameters params = ExponentialSpringParameters());

    /** Get the Transform specifying the location and orientation of the
     Contact Plane with respect to the Ground frame. This transform can be
     used to transform quantities expressed in the Contact Plane to the Ground
     frame and vice versa. */
    const Transform& getContactPlaneTransform() const;

    /** Get the body (i.e., the MobilizedBody) to which the resultant contact
    force is applied and for which this exponential spring was instantiated. */
    const MobilizedBody& getBody() const;

    /** Get the point (body station) that interacts with the contact plane
    and at which the resulting contact force is applied. The body station is
    expressed in the frame of the body for which this exponential spring was
    instantiated. */
    const Vec3& getStation() const;

    /** Set the customizable Topology-stage parameters on this exponential
    spring instance. To set the customizable Topology-stage parameters, create
    an ExponentialSpringParameters object, set the desired parameters on that
    object, and then call this method to modify the parameter values owned by
    this ExponentialSpringForce instance. Calling this method will invalidate
    the System at Stage::Topology; therefore, following a call to this method,
    System::realizeTopology() must be called before simulation can proceed.
    @param params Parameters object. */
    void setParameters(const ExponentialSpringParameters& params);

    /** Get a const reference to the parameters object owned by this
    exponential spring. The returned reference can be used to inspect the
    current values of the parameters or to create a copy of the parameters
    that is not owned by this exponential spring. The returned reference can
    also be used as the argument to setParameters() on a different spring
    object in order to assign its parameters to the values of this one. A
    non-const reference is not made available to the user in order to hinder
    casual changes made directly to the parameters owned by this subsystem.
    All changes made to the parameters owned by this subsystem must go through
    setParameters(). In this way, managing when system needs to be re-realized
    at Stage::Topology is streamlined. */
    const ExponentialSpringParameters& getParameters() const;

    /** Set the static coefficient of friction (μₛ) for this exponential
    spring. The value of μₛ is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, μₛ can be set at any
    time during a simulation. A change to μₛ will invalidate the System at
    Stage::Dynamics.
    @param state State object that will be modified.
    @param mus %Value of the static coefficient of friction. No upper bound.
    0.0 ≤ μₛ. If μₛ < μₖ, μₖ is set equal to μₛ. */
    void setMuStatic(State& state, Real mus);

    /** Get the static coefficient of friction (μₛ) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve μₛ. */
    Real getMuStatic(const State& state) const;

    /** Set the kinetic coefficient of friction (μₖ) for this exponential
    spring. The value of μₖ is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, μₖ can be set at any
    time during a simulation. A change to μₖ will invalidate the System at
    Stage::Dynamics.
    @param state State object that will be modified.
    @param muk %Value of the kinetic coefficient of friction. No upper bound.
    0.0 ≤ μₖ. If μₖ > μₛ, μₛ is set equal to μₖ. */
    void setMuKinetic(State& state, Real muk);

    /** Get the kinetic coefficient of friction (μₖ) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve μₖ. */
    Real getMuKinetic(const State& state) const;

    /** Set the Sliding state of the spring.
    @param state State object on which the new value will be set.
    @param sliding New value of Sliding.  0.0 ≤ sliding ≤ 1.0 */
    void setSliding(State& state, Real sliding);

    /** Get the Sliding state of the spring.
    @param state State object from which to retrieve Sliding. */
    Real getSliding(const State& state) const;

    /** Reset the friction spring zero (elastic anchor point) so that it
    coincides with the projection of the spring station onto the contact
    plane. This step is often needed at the beginning of a simulation to
    ensure that a simulation does not begin with large friction forces.
    After this call, the elastic portion of the friction force (fxyElas)
    should be 0.0. Calling this method will invalidate the System at
    Stage::Dynamics.
    @param state State object on which to base the reset. */
    void resetSpringZero(State& state) const;

    /** Get the elastic part of the normal force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getNormalForceElasticPart(const State& state,
        bool inGround = true) const;

    /** Get the damping part of the normal force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getNormalForceDampingPart(const State& state,
        bool inGround = true) const;

    /** Get the normal force. The system must be realized to Stage::Dynamics
    to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getNormalForce(const State& state,
        bool inGround = true) const;

    /** Get the instantaneous coefficient of friction (μ). The system must be
    realized to Stage::Dynamics to access this data. μ is obtained by using
    the Sliding state to transition between μₖ and μₛ:

            μ = μₛ - Sliding*(μₛ - μₖ)

    Because 0.0 ≤ Sliding ≤ 1.0, μₖ ≤ μ ≤ μₛ.
    @param state State object on which to base the calculations. */
    Real getMu(const State& state) const;

    /** Get the friction limit. The system must be realized to Stage::Dynamics
    to access this data.
    @param state State object on which to base the calculations. */
    Real getFrictionForceLimit(const State& state) const;

    /** Get the elastic part of the friction force. The system must be
    realized to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getFrictionForceElasticPart(
        const State& state, bool inGround = true) const;

    /** Get the damping part of the friction force. The system must be
    realized to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getFrictionForceDampingPart(
        const State& state, bool inGround = true) const;

    /** Get the total friction force. The total frictional force is always
    just the sum of the elastic part and the damping part of the frictional
    force, which may be obtained separately by calling
    getFrictionalForceElasticPart() and getFrictionalForceDampingPart().
    The system must be realized to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getFrictionForce(const State& state, bool inGround = true) const;

    /** Get the total force applied to the body by this
    ExponentialSpringForce instance. The total force is the vector sum of the
    friction force, which may be obtained by a call to getFrictionForce(), and
    the normal force, which may be obtained by a call to getNormalForce().
    The system must be realized to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getForce(const State& state, bool inGround = true) const;

    /** Get the position of the body station (i.e., the point on the body at
    which the force generated by this ExponentialSpringForce is applied). This
    method differs from getStation() in terms of the frame in which the
    station is expressed. getStation() expresses the point in the frame of
    the MobilizedBody. getStationPosition() expresses the point either in the
    Ground frame or in the frame of the Contact Plane. The system must be
    realized to Stage::Position to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getStationPosition(const State& state, bool inGround = true) const;

    /** Get the velocity of the body station (i.e., the point on the body at
    which the force generated by this ExponentialSpringForce is applied).
    The system must be realized to Stage::Velocity to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getStationVelocity(const State& state, bool inGround = true) const;

    /** Get the position of the friction spring zero (p₀), which will always
    lie in the Contact Plane. p₀ is the elastic anchor point of the damped
    linear spring used in Friction Model 2. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to base the calculations.
    @param inGround Flag for choosing the frame in which the returned quantity
    will be expressed. If true (the default), the quantity will be expressed
    in the Ground frame. If false, the quantity will be expressed in the frame
    of the contact plane. */
    Vec3 getFrictionSpringZeroPosition(const State& state,
        bool inGround=true) const;

private:
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(
        ExponentialSpringForce, ExponentialSpringForceImpl, Force);
    Stage getStage(const State& state) const {
        return getForceSubsystem().getStage(state);
    }
};


} // namespace SimTK
#endif // SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

