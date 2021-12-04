#ifndef SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_
#define SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021 Authors.
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
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ExponentialSpringParameters.h"

namespace SimTK {

class MultibodySystem;
class MultibodySubsystem;
class MobilizedBody;
class ExponentialSpringForceImpl;

//=============================================================================
/** Class ExponentialSpringForce uses an 'exponential spring' as a means
of modeling contact of a point on a MobilizedBody with a contact plane that is
fixed to Ground. In practice, you should choose a number of body-fixed points
('Stations' in Simbody vocabulary), strategically located across the
surface of a MobilizedBody, and create an ExponentialSpringForce that acts at
each of those points. For example, if the body were a cube, you would likely
choose to place an exponential spring at each corner of the cube. The contact
plane is typically used to model interactions with a floor, but need not be
limited to this use case. The contact plane can be rotated and displaced
relative to the ground frame and so can be used to model a wall or ramp, for
example.

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
springs accurately model contact force at the level of fundamental forces; we
only suggest that the use of a force field that acts over a great distance is
reasonable as long as the force gets sufficiently small sufficiently quickly.

Introducing any additional modeling approximation into a simulation is hard to
justify, however, unless there are significant benefits. In the case of
exponential springs there are two: simplicity and speed. With exponential
springs, there is no need to search for body intersections (the number active
of springs does not change during a simulation) and there is no need to
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
that opposes any slide velocity (i.e., f = -kᵥv). This latter approach is
relatively simple to implement, but will not entirely eliminate drift. By
making kᵥ very large, the drift velocity can be made small but it cannot be
brought to zero. And, unfortunately, as kᵥ is made large, the system equations
can become stiff, requiring integrator step size to be reduced, sometimes
considerably (at least in explicit integrators).

Part of the speed-up offered by ExponentialSpringForce derives, not from the
use of an exponential for the normal force, but from the friction model
employed. ExponentialSpringForce uses a spring-based frictional model that
includes a velocity-depending (damping) term AND a position-dependent
(elastic) term. By including the elastic term, drift is entirely eliminated
and relatively large integration step sizes are maintained.

In initial comparisons, using class ExponentialSpringForce to model
contact resulted in simulation cpu times that were typically 4 times less
(and sometimes 50 times less) than when using class CompliantContactSubsystem,
and yet the simulated motions were similar. These comparisons can be
reproduced by building and running the Test_Adhoc -
ExponentialSpringsComparison project that is included with Simbody.

Details of the exponential spring contact model are available in the
following publication:

        Anderson F.C. and Pandy M.G. (1999). A dynamics optimization
        solution for vertical jumping in three dimensions. Computer Methods
        in Biomechanics and Biomedical Engineering 2(3):201-231.

The current class makes several improvements to that contact model, most
notably including 1) the ability to rotate and translate the contact plane
and 2) the ability to specify both a static and a kinetic coefficient of
friction. The computational details of the contact model implemented by
this class follow below.

----------------------------------
Computations and Coordinate Frames
----------------------------------
The positive z-axis of the contact plane defines its normal. The positive
z-axis is the axis along which the repelling normal force (modeled using an
exponential) is applied. The x-axis and y-axis of the contact plane together
define the tangent plane. Friction forces will always be tangent to the x-y
plane.

In the equations below, a variable with a "z" suffix (e.g., pz or vz) refers
to a quantity that is normal to the contact plane. A variable with an "xy"
suffix refers to a quantity that lies in the contact plane.

### Normal Force (positive z-axis)

The elastic part of the normal force is computed using an exponential
whose shape is a function of three parameters (d₀, d₁, and d₂):

        fzElastic = d₁exp(-d₂(pz-d₀))

-Note that pz is the displacement of the body spring station above (pz > 0.0)
or below (pz < 0.0) the contact plane. The default values of the shape
parameters were chosen to maximize integration step size while maintaining a
number of constraints (e.g., the normal force must fall below 0.01 Newtons
when pz > 1.0 cm). The damping part of the normal force is linear in velocity
and scaled by the elastic part:

        fzDamping = - kzᵥ vz fzElastic,

where vz is the normal component of the velocity of the specified Station and
kzᵥ is the damping coefficient for the normal direction. All together, the
spring force in the normal direction is given by

        fz  = fzElastic + fzDamping
            = d₁exp(d₂(py-d₀)) - kzᵥ vz d₁exp(d₂(pz-d₀)))
            = d₁exp(d₂(pz-d₀)) (1 - kzᵥ vz)

which has the form of the Hunt & Crossley damping model:

        K. H. Hunt and F. R. E. Crossley (1975). Coefficient of Restitution
        Interpreted as Damping in Vibroimpact. ASME Journal of Applied
        Mechanics, pp. 440-445.

### Friction Force (x-y plane)

Friction force is implemented using a damped linear spring. Note that the
friction force cannot be represented by a scalar (like the normal force)
because it is free to point any direction in the xy plane. The elastic and
viscous terms are given by

        fricElas = -kₚ (pxy-p₀)

        fricDamp = -kᵥ vxy

and the total friction force is given by

        friction = fricElas + fricDamp

where kₚ is the spring elasticity, kᵥ is the spring viscosity, pxy is the
position of the body station projected onto the contact plane, p₀ is the
current spring zero, which always resides in the contact plane, and vxy is
the velocity of the body station in the contact plane expressed in the
contact plane.

By default, given a value for kₚ, the value of kᵥ is computed so as to
result in critical damping for a specified mass:

        kᵥ = 2.0 * sqrt(kₚ*mass)

Valid values of kₚ can range widely (e.g., kₚ = 1,000 to kₚ = 1,000,000)
depending on the material properties of the MobilizedBody and the contact
plane. In addition to being set by the above equation for critical damping,
kᵥ can be set to any positive value, independent of kₚ. In general, the higher
kₚ and kᵥ, the smaller the integration step size will need to be in order to
produce an accurate integration.

### A Moving Spring Zero

When the computed frictional force exceeds the allowed limit, it is capped at
the limit (fxzLimit):

        fxyLimit = μ fz
        if(friction.norm() > fxyLimit) {
        friction = fxyLimit * friction.normalize()

where μ is the instantaneous coefficient of friction (more below).

If the magnitude of the elastic part of the friction force by itself exceeds
fxyLimit, a new spring zero is found such that the magnitude of the
elastic part would be equal to fxyLimit:

        if(fricElas.norm() > fxyLimit)  p₀ = pxy + (pxy-p₀)/kₚ

In Simbody, p₀ is handled as an Auto Update Discrete State. See
State::allocateAutoUpdateDiscreteVariable() for a detailed description. Any
change to p₀ is made to the Update Cache (not to the State directly), and the
integrator copies this cache value to the actual p₀ State after a successful
integration step is obtained.

### Coefficients of Friction

Coefficients of kinetic (sliding) and static (fixed) friction can be
specified for the spring subject to the following constraints:

       0.0 ≤ μₖ ≤ μₛ

Note that there is no upper bound on μₛ.

A continuous state variable (Z = Sliding) is used to characterize the
sliding state of the spring.

        Sliding = 0     means fully fixed in place (lower bound)
        Sliding = 1     means fully sliding (upper bound)

The instantaneous coefficient of friction (mu) is calculated based on the
value of Sliding:

        μ = μₛ - Sliding*(μₛ - μₖ)

The time derivative of Sliding is used to drive Sliding toward the
extreme of 0.0 or 1.0, depending on the following criteria:

When the elastic part of the frictional force exceeds its limit,
Sliding is driven toward 1:

        if (fricElas.norm() > fxyLimit)  SlidingDot = (1.0-Sliding)/tau

When vxy falls below some specified "settle" velocity (e.g., 0.01 m/s) AND
the frictional force is less than its limit, Sliding is driven toward 0:

        else if (vxy.norm < 0.01)  SlidingDot = -Sliding/tau

Otherwise, no change to the Sliding state is made:

        else  SlidingDot = 0.0

In the above equations, tau is the characteristic time it takes for the Sliding
state to rise or decay. The motivation for using a continuous state variable
for Sliding is that, although the transition between fixed and sliding may
happen quickly, it does not happen instantaneously.  And, modeling Sliding
based on a differential equation ensures that μ is continuous.

### Future Enhancements

Future enhancements could include the capacity to
1) fix the contact plane to a body other than ground
2) use a polygonal mesh as the contact surface, and
3) move the specified MobilizedBody Station in manner that adapts to the
contact circumstances. For example, to model a coin rolling on a table
a small number of stations could be continually moved to the portion of the
coin that is closest to the table.

------
STATES
------
Each instance of ExponentialSpringForce posseses 4 states, which are listed
below in the appropriate category:

### DISCRETE STATES (parameters)
- μₛ = Static coefficient of friction.  0.0 ≤ μₛ
- μₖ = Kinetic coefficient of friction.  0.0 ≤ μₖ ≤ μₛ

### AUTO UPDATE DISCRETE STATE
- p₀ = Zero point (Vec3) of the frictional spring. p₀ always lies in the
contact plane.

### CONTINUOUS STATE
- Sliding = Indicator of whether the spring zero (p₀) has been moving
(sliding) or fixed in the contact plane.  0.0 ≤ Sliding ≤ 1.0.
The "Sliding" state is used to transition the instantaneous coefficient of
frction (μ) between μₖ and μₛ. A value of 0.0 indicates that p₀ is fixed in
place, in which case μ = μₛ. A value of 1.0 indicates that p₀ is sliding, in
which case μ = μₖ. A value between 0.0 and 1.0 indicates that a transition
from sliding to fixed or from fixed to sliding is undwerway, in which case
μₖ ≤ μ ≤ μₛ.

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
class SimTK_SIMBODY_EXPORT ExponentialSpringForce : public ForceSubsystem {
public:
    /** Construct an exponential spring force object with customized
    parameters.
    @param system The system being simulated or studied.
    @param contactPlane Transform specifying the location and orientation of
    the contact plane with respect to the Ground frame. The positive z-axis
    defines the normal of the contact plane; friction forces lie in (are
    tangent to) the x-y plane.
    @param body MobilizedBody that will interact / collide with the contact
    plane.
    @param station Point on the specified body at which the contact force
    will be applied. The position and velocity of this point relative to
    the contact plane determine the magnitude and direction of the contact
    force.
    @param mus Static coefficient of friction. No upper bound. 0.0 ≤ μₛ
    @param muk Kinetic coefficient of friction. 0.0 ≤ μₖ ≤ μₛ. During
    construction, if μₖ > μₛ, μₖ is set equal to μₛ.
    @param params Customized Topology-stage parameters. If params is omitted
    from the argument list, the default set of parameters is used. */
    ExponentialSpringForce(MultibodySystem& system,
        const Transform& contactPlane,
        const MobilizedBody& body, const Vec3& station,
        Real mus, Real muk,
        ExponentialSpringParameters params = ExponentialSpringParameters());

    /** Get the Transform specifying the location and orientation of the
     Contact Plane with respect to the Ground frame. This transform can be
     used to transform quantities expressed in the Contact Plane to the Ground
     frame and vice versa. The z-axis of the contact plane frame specifies the
     normal. */
    const Transform& getContactPlane() const;

    /** Get the body (i.e., the MobilizedBody) for which this exponential
    spring was instantiated. */
    const MobilizedBody& getBody() const;

    /** Get the point (station) that interacts with the contact plane and at
    which the resulting contact force is applied. The station is expressed in
    the frame of the body for which this exponential spring was
    instantiated. */
    const Vec3& getStation() const;

    /** Set the customizable Topology-stage parameters on this exponential
    spring instance. To do this, create an ExponentialSpringParameters object,
    set the desired parameters on that object, and then call this method to
    modify the parameter values owned by this ExponentialSpringForce instance.
    Calling this method will invalidate the System at Stage::Topology;
    therefore, following a call to this method, System::realizeTopology() must
    be called before simulation can proceed.
    @param params Parameters object. */
    void setParameters(const ExponentialSpringParameters& params);

    /** Get a const reference to the parameters object owned by this
    exponential spring. The returned reference can be used to inspect the
    current values of the parameters or to create a copy of the paramters that
    is not owned by this exponential spring. The returned reference can also be
    used as the argument to setParameters() on a different spring object in
    order to assign its parameters to the values of this one. A non-const
    reference is not made available to the user in order to hinder casual
    changes made directly to the parameters owned this subsystem. All changes
    made to the parameters owned by this subsystem must go through
    setParameters(). In this way, managing when system needs to be re-realized
    at Stage::Topology is streamlined. */
    const ExponentialSpringParameters& getParameters() const;

    /** Set the static coefficient of friction (μₛ) for this exponential
    spring.
    The value of μₛ is held in the System's State object. Unlike the
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
    spring.
    The value of μₖ is held in the System's State object. Unlike the
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

    /** Get the Sliding state of the spring.
    @param state State object from which to retrieve Sliding. */
    Real getSliding(const State& state) const;

    /** Reset the spring zero so that it coincides with the projection of the
    spring station onto the contact plane. This step is often needed at the
    beginning of a simulation to ensure that a simulation does not begin with
    large friction forces. After this call, the elastic portion of the
    friction force (fxzElas) should be 0.0. Calling this method will
    invalidate the System at Stage::Dynamics.
    @param state State object on which to base the reset. */
    void resetSpringZero(State& state) const;

    /** Get the elastic part of the normal force. The system must be realized to
    Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getNormalForceElasticPart(const State& state,
        bool inGround = true) const;

    /** Get the damping part of the normal force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getNormalForceDampingPart(const State& state,
        bool inGround = true) const;

    /** Get the normal force. The system must be realized to Stage::Dynamics
    to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getNormalForce(const State& state,
        bool inGround = true) const;

    /** Get the instantaneous coefficient of friction (μ). The system must be
    realized to Stage::Dynamics to access this data. μ is obtained by using
    the Sliding state to transition between μₖ and μₛ:
            μ = μₛ - Sliding*(μₛ - μₖ)
    Because 0.0 ≤ Sliding ≤ 1.0, μₖ ≤ μ ≤ μₛ.
    @param state State object on which to based the calculations. */
    Real getMu(const State& state) const;

    /** Get the friction limit. The system must be realized to Stage::Dynamics
    to access this data.
    @param state State object on which to based the calculations. */
    Real getFrictionForceLimit(const State& state) const;

    /** Get the elastic part of the friction force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getFrictionForceElasticPart(
        const State& state, bool inGround = true) const;

    /** Get the damping part of the friction force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getFrictionForceDampingPart(
        const State& state, bool inGround = true) const;

    /** Get the total friction force. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getFrictionForce(const State& state, bool inGround = true) const;

    /** Get the total spring force applied to the MobilizedBody. The system
    must be realized to Stage::Dynamics to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getForce(const State& state, bool inGround = true) const;

    /** Get the position of the spring station (i.e., the point at which the
    spring force is applied to the body). The system must be realized to
    Stage::Position to access this data.
    This method differs from getStation() in terms of the frame in which the
    station is expressed. getStation() expresses the point in the frame of
    the MobilizedBody. getStationPosition() expresses the point either in the
    Ground frame or in the frame of the Contact Plane.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getStationPosition(const State& state, bool inGround = true) const;

    /** Get the velocity of the spring station. The system must be realized to
    Stage::Velocity to access this data.
    @param state State object on which to based the calculations.
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getStationVelocity(const State& state, bool inGround = true) const;

    /** Get the position of the spring zero. The system must be realized
    to Stage::Dynamics to access this data.
    @param state State object from which to retrieve the data. 
    @param inGround Flag for choosing the frame in which the returned
    quantity will be expressed. If true, the quantity will be expressed in the
    Ground frame. If false, the quantity will be expressed in the frame of
    the contact plane. */
    Vec3 getSpringZeroPosition(const State& state, bool inGround = true) const;

private:
    /** Retrieve a writable reference to the underlying implementation. */
    ExponentialSpringForceImpl& updImpl();

    /** Retrieve a read-only reference to the underlying implementation. */
    const ExponentialSpringForceImpl& getImpl() const;
};


} // namespace SimTK
#endif // SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

