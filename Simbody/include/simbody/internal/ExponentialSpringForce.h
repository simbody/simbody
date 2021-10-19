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
/** ExponentialSpringData is a helper class that is used to store key
quantities associated with the ExponentialSpringForce Subsystem during a
simulation. An instance of this class serves as the data Cache Entry for the
ExponentialSpringForceImpl Subsystem. All of its member variables are
guaranteed to be calculated and set once the System has been realized
to Stage::Dynamics.

To understand what the quantities organized by this class represent, a basic
description of the contact problem that is solved, along with a description
of coordinate frame conventions, will be helpful.

Class ExponentialSpringForce computes and applies a contact force at a
specified point on a MobilizedBody (i.e., a Station) due to interaction of
that point with a specified contact plane. That plane is typically used to
model interactions with a floor, but need not be limited to this use case.
The contact plane can be rotated and displaced relative to the ground frame
and so can be used to model a wall or ramp, for example.

%Contact force computations are carried out in the frame of the contact plane.
The positive y-axis of the contact frame defines the normal of the
contact plane. The positive y-axis is the axis along which a repelling
normal force (modeled using an exponential) is applied. The x-axis and
z-axis of the contact frame are tangent to the contact plane. The friction
force will always lie in x-z plane.

Member variables with a "y" suffix (e.g., py, vy, or fy) indicate that
these quantities are directed normal to the contact plane.  Member varaibles
with an "xz" suffix (e.g., pxz, vxz, or fxz) indicate that these quantities
lie in the contact plane (or tangent to it) and are associated with the
friction force.

Member variables with a "_G" suffix are expressed in the ground frame. Member
variables without a "_G" suffix are expressed in frame of the contact plane.

This class has public variable members for easy access, which would normally
cause trouble. The typical user, however, will only have access to a const
reference to the ExponentialSpringData object owned by the State and so will
only be able to read its member variables not change them. A non-const refernce
can be obtained only through the underlying implementation
(ExponentialSpringForceImpl); only ExponentialSpringForceImpl has reason
and knows how to assign values to the member variables.

A copy constructor and an assignment operator are provided in the event a user
wants to store the calculated data independent of the State.

@see ExponentialSpringForce */
class SimTK_SIMBODY_EXPORT ExponentialSpringData {
public:
    /** Default Constructor. */
    ExponentialSpringData();

    /** Copy Constructor.
    @param source Const reference to the source object that should be copied.*/
    ExponentialSpringData(const ExponentialSpringData& source);

    /** Copy Assignment Operator.
    @param source Const reference to the source object to which this instance
    should be assigned. */
    ExponentialSpringData& operator=(const ExponentialSpringData& source);

    /** Position of the body spring station in the ground frame. */
    Vec3 p_G;

    /** Velocity of the body spring station in the ground frame. */
    Vec3 v_G;

    /** Position of the body spring station in the frame of the contact
    plane. */
    Vec3 p;

    /** Velocity of the body spring station in the frame of the contact
    plane. */
    Vec3 v;

    /** Displacement of the body spring station normal to the floor expressed
    in the frame of the contact plane. */
    Real py;

    /** Velocity of the body spring station normal to the floor expressed
    in the frame of the contact plane. */
    Real vy;

    /** Position of the body spring station projected onto the floor expressed
    in the frame of the contact plane. */
    Vec3 pxz;

    /** Velocity of the body spring station in the plane of the floor
    expressed in the frame of the contact plane. */
    Vec3 vxz;

    /** Position of the body spring station projected onto the floor expressed
    in the ground frame.  This quantity is not needed to calculate the spring
    force, but is likely desired for visualization. */
    Vec3 pxz_G;

    /** Elastic force in the normal direction expressed in the frame of the
    contact plane. */
    Real fyElas;

    /** Damping force in the normal direction expressed in the frame of the
    contact plane. */
    Real fyDamp;

    /** Total normal force expressed in the frame of the contact plane. */
    Real fy;

    /** Instantaneous coefficient of friction. */
    Real mu;

    /** Limit of the frictional force. */
    Real fxyLimit;

    /** Elastic frictional force expressed in the frame of the contact plane.*/
    Vec3 fricElas;

    /** Damping frictional force expressed in the frame of the contact plane.*/
    Vec3 fricDamp;

    /** Total frictional force (elastic + damping) expressed in the frame of
    the contact plane. */
    Vec3 fric;

    /** Magnitude of the frictional force. */
    Real fxy;

    /** Resultant spring force (normal + friction) expressed in the floor
    frame. */
    Vec3 f;

    /** Resultant spring force (normal + friction) expressed in the ground
    frame. */
    Vec3 f_G;
};


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
is more than 1.0 cm above the ground and extremely small (less than 0.000001 N) as soon as the body
station is more than 2.0 cm above the ground. The sizes of these forces are
likely small compared to the errors made in neglecting air resistance,
assuming a uniform gravitational field, or estimating inertial parameters.

As a side note, electrostatic forces, which are the fundamental force from
which all conventional contact forces arise, are never actually zero either.
They obey Coulomb's Law (Fe = ke*q1*q2/r^2) which states that an
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
plane is less than the limiting value of the static friction force (μₛ fₙ),
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
(and sometimes 50 times less) than when using class CompliantContactSubsystem.
The simulated motions were similar. These comparisons can be reproduced by
building and running the Test_Adhoc - ExponentialSpringsComparison project
that is included with Simbody.

Details of the the exponential spring contact model are available in the
following publication:

        Anderson F.C. and Pandy M.G. (1999). A dynamics optimization
        solution for vertical jumping in three dimensions. Computer Methods
        in Biomechanics and Biomedical Engineering 2(3):201-231.

The current class makes several improvements to that contact model, most
notably including 1) the ability to rotate and translate the contact plane
and 2) the ability to specify static and kinetic coefficients of friction.
The computational details of the contact model implemented by this class
follow below.

----------------------------------
Computations and Coordinate Frames
----------------------------------
%Contact force computations are carried out in the frame of the contact plane.
The positive y-axis of the contact plane defines its normal. The positive
y-axis is the axis along which the repelling normal force (modeled using an
exponential) is applied. The x-axis and z-axis of the contact plane together
define the tangent plane. Friction forces will always lie in x-z plane.

In the equations below, a variable with a "y" suffix (e.g., py or vy) refers
to a quantity that is normal to the contact plane. A variables with an "xz"
suffix refers to a quantity that lies in the contact plane.

### Normal Force (positive y-axis)

The elastic part of the normal force is computed using an exponential
whose shape is a function of three parameters (d₀, d₁, and d₂):

        fyElastic = d₁exp(d₂(py-d₀)),

where py is the height above the contact plane of the specified Station
on the MobilizedBody. The damping part is linear in velocity and scaled by the
elastic part:

        fyDamping = - kyᵥ vy fyElastic,

where vy is the normal component of the velocity of the specified Station and
kyᵥ is the damping coefficient for the normal direction. All together, the
spring force in the normal direction is given by

        fy  = fyElastic + fyDamping
            = d₁exp(d₂(py-d₀)) - kyᵥ vy d₁exp(d₂(py-d₀)))
            = d₁exp(d₂(py-d₀)) (1 - kyᵥ vy)

which has the form of the Hunt & Crossley damping model. See K. H. Hunt and
F. R. E. Crossley (1975). Coefficient of Restitution Interpreted as Damping in
Vibroimpact. ASME Journal of Applied Mechanics, pp. 440-445.

### Friction Force (x-z plane)

Friction force (a vector because its direction in the x-z plane can change)
is implemented using a damped linear spring. The elastic and viscous terms
are given by

        fricElas = -kₚ (pxy-p₀)

        fricDamp = -kᵥ vxy

and the total friction force is given by

        friction = fricElas + fricDamp

where kₚ is the spring elasticity, kᵥ is the spring viscosity, pxy is the
position of the body station projected onto the contact plane, p₀ is the
current spring zero, which always resides in the contact plane, and vxy is
the velocity of the body station in the contact plane.

By default, the frictional parameters kₚ and kᵥ are chosen to result in
critical damping for a specified mass:

        kᵥ = 2.0 * sqrt(kₚ*mass)

Valid values of kₚ can range widely (e.g., kₚ = 1,000 to kₚ = 1,000,000)
depending on the material properties of the MobilizedBody and the contact
plane. In addition to being set by the above equation for critical damping,
kᵥ can be set to 0.0 or to some other independent quantity. In general, the
higher kₚ and kᵥ, the smaller the integration step size will need to be in
order to produce an accurate integration.

### A Moving Spring Zero

When the computed frictional force exceeds the allowed limit, it is capped at
the limit (fxyLimit):

       fxyLimit = μ fy
       if(friction.norm() > fxyLimit) {
           friction = fxyLimit * friction.normalize()
       }

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

       0.0 ≤ μₖ ≤ μₛ ≤ 1.0

A continuous state variable (Z = Sliding) is used to characterize the
sliding state of the spring.

        Sliding = 0     means fully fixed in place (lower bound)
        Sliding = 1     means fully sliding (upper bound)

The instantaneous coefficient of friction (mu) is calculated based on the
value of Sliding:

        μ = μₛ - Sliding (μₛ - μₖ)

The time derivative of Sliding is used to drive Sliding toward the
extreme of 0.0 or 1.0, depending on the following criteria:

When the frictional force exceeds its limit, Sliding is driven toward 1:

        if (friction.norm() >= fxyLimit)  SlidingDot = (1.0-Sliding)/tau

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

Future enhancements might include the capacity to
1) use a polygonal mesh as the contact surface, and
2) move the specified MobilizedBody Station in manner that adapts to the
contact circumstances. For example, to model a coin rolling on a table
a small number of stations could be continually moved to the portion of the
coin that is closest to the table.

------
STATES
------
Each instance of ExponentialSpringForce posseses 6 states. These are listed
below in the appropriate category:

### DISCRETE STATES (parameters)
μₛ = Static coefficient of friction.  0.0 ≤ μₛ ≤ 1.0  
μₖ = Kinetic coefficient of friction.  0.0 ≤ μₖ ≤ μₛ  
station = Point (Vec3) expressed in the body frame at which the force is
exerted on the MobilizedBody.

### AUTO UPDATE DISCRETE STATE
p₀ = Zero point (Vec3) of the frictional spring in the contact plane. p₀
always lies in the contact plane and is expressed in the frame of the contact
plane.

### CONTINUOUS STATE
Sliding = Indicator of whether the MobilizedBody is sliding relative to the
contact plane.

----------
PARAMETERS
----------
Customizable Topology-stage parameters specifying the characteristics of the
exponential spring are managed using ExponentialSpringParameters.

----
DATA
----
Calculated quantities are cached when the System is realized through
Stage::Dynamics. The cached data is made accessible via the helper class
ExponentialSpringData. */
class SimTK_SIMBODY_EXPORT ExponentialSpringForce : public ForceSubsystem {
public:
    /** Construct an exponential spring force object with customized
    parameters.
    @param system The system being simulated or studied.
    @param contactPlane Transform specifying the location and orientation of
    the contact plane. The positive y-axis defines the normal of the contact
    plane; friction forces occur in the x-z plane.
    @param body %Body that will interact / collide with the contact plane.
    @param station Point on the specified body at which the contact force
    will be applied. The position and velocity of this point relative to
    the contact plane determine the magnitude and direction of the contact
    force.
    @param mus Static coefficient of friction. No upper bound. 0.0 ≤ μₛ
    @param muk Kinetic coefficient of friction. 0.0 ≤ μₖ ≤ μₛ. If μₖ > μₛ,
    μₖ is set equal to μₛ.
    @param params Customized parameters. */
    ExponentialSpringForce(MultibodySystem& system,
        const Transform& contactPlane,
        const MobilizedBody& body, const Vec3& station,
        Real mus, Real muk,
        ExponentialSpringParameters params = ExponentialSpringParameters());

    /** Set the customizable parameters on this exponential spring instance.
    To do this, create an ExponentialSpringParameters object, set the desired
    parameters on that object, and then call this method to modify the
    parameter values owned by this ExponentialSpringForce instance. Calling
    this method will invalidate the System at Stage::Topology; therefore,
    following a call to this method, System::realizeTopology() must be called
    before simulation can proceed.
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
    is streamlined. */
    const ExponentialSpringParameters& getParameters() const;

    /** Set the static coefficient of friction (μₛ) for this exponential
    spring.
    The value of μₛ is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, μₛ can be set at any
    time during a simulation. A change to μₛ will invalidate the System at
    Stage::Dynamics.
    @param state State object that will be modified.
    @param mus Value of the static coefficient of friction. No upper bound.
    0.0 ≤ μₛ. If μₛ < μₖ, μₖ is set equal to μₛ. */
    void setMuStatic(State& state, Real mus);

    /** Get the static coefficient of friction (μₛ) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve μₛ. */
    const Real& getMuStatic(const State& state) const;
 
    /** Set the kinetic coefficient of friction (μₖ) for this exponential
    spring.
    The value of μₖ is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, μₖ can be set at any
    time during a simulation. A change to μₖ will invalidate the System at
    Stage::Dynamics.
    @param state State object that will be modified.
    @param muk Value of the kinetic coefficient of friction. No upper bound.
    0.0 ≤ μₖ. If μₖ > μₛ, μₛ is set equal to μₖ. */
    void setMuKinetic(State& state, Real muk);

    /** Get the kinetic coefficient of friction (μₖ) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve μₖ. */
    const Real& getMuKinetic(const State& state) const;

    /** Reset the spring zero so that it coincides with the projection of the
    spring station onto the contact plane. This step is often needed at the
    beginning of a simulation to ensure that a simulation does not begin with
    large friction forces. After this call, the elastic portion of the
    friction force (fxyElas) should be 0.0. Calling this method will
    invalidate the System at Stage::Dynamics.
    @param state State object on which to base the reset. */
    void resetSpringZero(State& state) const;

    /** Retrieve a const reference to this exponential spring's data cache.
    Before a call to this method, the System should be realized through
    Stage::Dynamics; an exception will be thrown otherwise.  The returned
    reference is const to prevent changes to the underlying data cache.
    Only the underlying implementation (ExponentialSpringForceImpl) has
    write access to the data cache.
    @param state State object from which to retrieve the data. */
    const ExponentialSpringData& getData(const State& state) const;

    /** Retrieve a writable reference to the underlying implementation. */
    ExponentialSpringForceImpl& updImpl();

    /** Retrieve a read-only reference to the underlying implementation. */
    const ExponentialSpringForceImpl& getImpl() const;
};


} // namespace SimTK
#endif // SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

