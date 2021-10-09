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

#include "Simbody.h"

namespace SimTK {


//=============================================================================
/** ExponentialSpringParameters is a helper class used to customize the
characteristics of an ExponentialSpringForce instance. It manages all
parameters that are implemented as Stage::Topology variables. To customize any
of the Topology-stage parameters on an ExponentialSpringForce instance, the
user should

1) Create an ExponentialSpringParameters object. For example,

        ExponentialSpringParameters myParams;

2) Use the available 'set' methods in this class to change the parameters
of that object. For example,

        myParams.setNormalViscosity(0.25);

3) Use ExponentialSpringForce::setParameters() to alter the parameters of one
(or many) ExponentialSpringForce instances. For example,

        ExponentialSpringForce spr1, spr2;
        spr1.setParameters(myParams);
        spr2.setParameters(myParams);

4) Realize the system to Stage::Topology.  When a new set of parameters is
set on an ExponentialSpringForce instance, as above in step 3, the System will
be invalidated at Stage::Topology.  The System must therefore be realized at
Stage::Topology before a simulation can proceed.

        system.realizeTopology();

Note that each ExponentialSpringForce instance owns its own private
ExponentialSpringParameters object. The myParams object is just used to set
the desired parameter values of the privately owned parameters object. It is
fine for objects like myParams to go out out of scope or for myParams objects
allocated from the heap to be deleted.

Therefore, also note that the parameter values possessed by an
ExponentialSpringForce instance do not necessarily correspond to the values
held by an instance of this class until a call to
ExponentialSpringForce::setParameters() is made.

The default values of the parameters held by ExponentialSpringParameters
work well for typical contact interactions, but clearly may not be
appropriate for simulating many contact interactions. For example, one might
want to simulate an interaction in which very little energy is dissipated
during a bounce.

The default values of the parameters are expressed in units of Newtons,
meters, seconds, and kilograms; however, you may use an alternate set of
self-consistent units by re-specifying all parameters.

Finally, there are 2 quantities not managed by this class that can be used
to further customize the behavior of an ExponentialSpringForce instance: the
static coefficient of friction (mus) and the kinetic coefficient of friction
(muk). mus and muk are parameters, but they are not implemented as
Topology-stage variables; they are implemented as Dynamics-stage discrete
state variables. Thus, they can be changed during the course of a simulation
without invalidating the System at Stage::Topology. This functionality allows
a contact plane to posses non-uniform frictional characteristics across its
surface. For example, a patch of ice on a sidewalk could be modeled. mus and
muk can be set during a simulation using ExponentialSpringForce::setMuStatic()
and ExponentialSpringForce::setMuKinetic().

For an explanation of the equations that underlie exponential spring forces,
refer to class ExponentialSpringForce.*/
class SimTK_SIMBODY_EXPORT ExponentialSpringParameters {
public:
    /** Default Constructor. Member variables are initialized to default
    values. */
    ExponentialSpringParameters();

    /** Copy constructor. Member variables are initialized to the values of
    the specified source object.
    @param source Const reference to the source object to be copied. */
    ExponentialSpringParameters(const ExponentialSpringParameters& source);

    /** Copy assignment operator.
    @param source Const reference to the source object to which this instance
    should be assigned. */
    ExponentialSpringParameters&
        operator=(const ExponentialSpringParameters& source);

    /** Set the parameters that determine the shape of the exponential.
    @param d0 shifts the exponential function up and down with respect to
    the contact plane. Its default value is 0.0065905 m (~7 mm above the
    contact plane). This slight upward shift eliminates significant
    penetration into the contact plane. d0 can positive, have a value of 0.0,
    or be negative.
    @param d1 linearly scales the applied force up or down. Its default
    value is 0.5336 (unitless). d1 should be positive to generate a repulsive
    force.
    @param d2 linearly scales the exponent. Its default value is 1150.0 / m.
    Larger values of d2 make the exponential curve rise more rapidly as py
    gets less than d0. d1 should be positive. */
    void setShapeParameters(Real d0, Real d1 = 0.5336, Real d2 = 1150.0);

    /** Get the parameters that determine the shape of the exponential.
    @see setShapeParameters() */
    void getShapeParameters(Real& d0, Real& d1, Real& d2) const;

    /** Set the viscosity of the exponential part of an ExponentialSpringForce.
    This viscosity only multiplies the component of the velocity of the
    spring station that is normal to the contact plane. To eliminate energy
    dissipation due to motion in the normal direction, set kvNorm equal to 0.0.
    @param kvNorm Viscosity of a spring in the normal direction.
    Its default value is 0.5 N*s/m. kvNorm should be positive or zero. */
    void setNormalViscosity(Real& kvNorm);

    /** Get the viscosity of the exponential part of an ExponentialSpringForce.
    The exponential part of the contact force is always normal to the contact
    plane.
    @returns Viscosity in the normal direction */
    Real getNormalViscosity() const;

    /** Set both the elasticity and viscosity of the friction spring
    associated with an exponential spring. The value that is set for the
    viscosity is computed so that critical damping would result for a
    specified mass (i.e., kv = 2*sqrt(kp*mass)). A call to this method
    overrides any values set previously by calls to setElasticity() or
    setViscosity().
    @param kp Elasticity of the friction spring. kp should be positive.
    @param mass Mass of the body for which critical damping would be achieved.
    Articulated bodies generally don't have a constant mass as it relates to
    acceleration in a particular direction, so think of this mass as a kind
    of average mass. A default mass of 1.0 kg is used if a mass is not
    specified. */
    void setElasticityAndComputeViscosity(Real kp, Real mass = 1.0);

    /** Set the elasticity of the friction spring associated with an
    exponential spring. A call to this method overrides any value of
    elasticity previously set by a call to setElasticityAndComputeViscosity().
    @param kp Elasticity of the friction spring. Its default value is
    2000.0 N/m. kp should be positive. */
    void setElasticity(Real kp);

    /** Get the elasticity of the friction spring. The value of the elasticity
    is the default value (2000.0 N/m) or the value set by a call to either
    setElasticity() or setElasticityAndComputeVicosity(), whichever was called
    most recently.
    @returns Elasticity of the friction spring */
    Real getElasticity() const;

    /** Set the viscosity of the friction spring. A call to this method
    overrides any value of viscosity previously set by a call to
    setElasticityAndComputeViscosity(). Setting the viscosity equal
    to 0.0 is fine, but may not have the expected result.  If a body is not
    sliding, setting kv = 0.0 will simply allow the body to vibrate in place
    indefinitely. If a body is sliding, even if kv = 0.0, the kinetic energy
    of the body will still be dissipated because the frictional force will
    not be zero. (The elastic part of the friction spring is stretched and so
    still applies a force, but the potential energy stored in the spring is
    not increased because the spring zero is continually released.) The only
    way to eliminate energy dissipation entirely is to set the coefficients
    of friction equal to 0.0, which can be done by a call to
    ExponentialSpringForce::setMuStatic(0.0).
    @param kv Viscosity of the friction spring. Its default value is
    2.0*sqrt(kp*mass) = 2.0*sqrt(2000*1) ~= 89.4427 N*s/m. kv should be 0.0
    or positive. */
    void setViscosity(Real kv);

    /** Get the viscosity of the friction spring. The value of the viscosity
    is the default value (~89.4427 N*s/m) or the value set by a call to either
    setViscosity() or setElasticityAndComputeViscosity(), whichever was called
    most recently.
    @returns Viscosity of the friction spring */
    Real getViscosity() const;

    /** Set the time constant for transitioning back and forth between the
    static and kinetic coefficients of friction. The transition is mediated
    by a rising or falling exponential that is asymptotic to mu static or
    or mu kinetic respectively.
    @param tau Time constant for sliding transitions. The default value of tau
    is 0.01 s. tau must be positive. */
    void setSlidingTimeConstant(Real tau);

    /** Get the time constant for transitioning back and forth between the
    static and kinetic coefficients of friction.
    @returns time constant for sliding transitions. */
    Real getSlidingTimeConstant() const;

    /** Set the velocity below which the coefficient of friction transitions
    to the static coefficient of friction.
    @param vSettle Settle velocity. It's default value is 0.01 m/s. vSettle
    must be positive.*/
    void setSettleVelocity(Real vSettle);

    /** Get the velocity below which the coefficient of friction transitions
    to the static coefficient of friction.
    @returns Settle velocity */
    Real getSettleVelocity() const;

protected:
    Real d0;
    Real d1;
    Real d2;
    Real kvNorm;
    Real kpFric;
    Real kvFric;
    Real kTau;
    Real vSettle;

    /** The member variables managed by this class are used by the underlying
    implementation of the ExponentialSpringForce Subsystem.  Direct access
    is given to ExponentialSpringForceImpl for reasons of convenience. Direct
    access is not given to other classes to ensure that parameters are set
    only to valid values. */
    friend class ExponentialSpringForceImpl;
};


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
making kᵥ very large, the drift velocity can be made negligibly small but it
cannot be brought exactly to zero. And, unfortunately, as kᵥ is made very
large, the system equations become stiff, requiring integrator step size to
be reduced, sometimes considerably (at least in explicit integrators).

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

Normal %Force (positive y-axis) \n
The elastic part of the normal force is computed using an exponential
whose shape is a function of three parameters (d0, d1, and d2):

        fyElastic = d1 * exp(d2 * (py - d0)),

where py is the height above the contact plane of the specified Station
on the MobilizedBody. The damping part is linear in velocity and scaled by the
elastic part:

        fyDamping = - kv*vy * fyElastic,

where vy is the normal component of the velocity of the specified Station.
All together, the spring force in the normal direction is given by

        fy  = fyElastic + fyDamping
            = d1*exp(d2*(py-d0)) - kv*vy * d1*exp(d2*(py-d0))

Friction %Force (x-z plane) \n
Friction force (a vector because its direction in the x-z plane can change)
is implemented using a damped linear spring. The elastic and viscous terms
are given by

        fricElas = -kp*(pxy-p0)

        fricDamp = -kv*(vxy)

and the total friction force is given by

        friction = fricElas + fricDamp

where kp is the spring elasticity, kv is the spring viscosity, pxy is the
position of the body station projected onto the contact plane, p0 is the
current spring zero, which always resides in the contact plane, and vxy is
the velocity of the body station in the contact plane.

By default, the frictional parameters kp and kv are chosen to result in
critical damping for a specified mass:

        kv = 2.0 * sqrt(kp*mass)

Valid values of kp can range widely (e.g., kp = 1,000 to kp = 1,000,000)
depending on the material properties of the MobilizedBody and the contact
plane. In addition to being set by the above equation for critical damping,
kv can be set to 0.0 or to some other independent quantity. In general, the
higher kp and kv, the smaller the integration step size will need to be in
order to produce an accurate integration.

A Moving Spring Zero\n
When the computed frictional force exceeds the allowed limit, it is capped at
the limit (fxyLimit):

       fxyLimit = mu*fy
       if(friction.norm() > fxyLimit) {
           friction = fxyLimit * friction.normalize()
       }

where mu is the instantaneous coefficient of friction (more below).

If the magnitude of the elastic part of the friction force by itself exceeds
fxyLimit, a new spring zero (p0) is found such that the magnitude of the
elastic part would be equal to fxyLimit:

       if(fricElas.norm() > fxyLimit)  p0New = pxy + (pxy-p0)/kp

In Simbody, p0 is handled as an Auto Update Discrete State.  Any change
to p0 is made to the Update Cache (not to the State directly), and the
integrator swaps this cache value with the actual p0 State after a
successful integration step is obtained.

Coefficients of Friction \n
Coefficients of kinetic (sliding) and static (fixed) friction can be
specified for the spring subject to the following constraints:

       0.0 ≤ muk ≤ mus ≤ 1.0

A continuous state variable (Z = Sliding) is used to characterize the
sliding state of the spring.

        Sliding = 0     means fully fixed in place (lower bound)
        Sliding = 1     means fully sliding (upper bound)

The instantaneous coefficient of friction (mu) is calculated based on the
value of Sliding:

        mu = mus - Sliding*(mus - muk)

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
based on a differential equation ensures that mu is continuous.

Future Enhancements \n
Future enhancements might include the capacity to
1) use a polygonal mesh as the contact surface, and
2) move the specified MobilizedBody Station in manner that adapts to the
contact circumstances. For example, to model a coin rolling on a table
a small number of stations could be continually moved to the portion of the
coin that is closest to the table.

------
STATES
------
PARAMETERS (implemented as discrete states)\n
    mus = Static coefficient of friction.  0.0  mus ≤ 1.0 \n
    muk = Kinetic coefficient of friction.  0.0 ≤ muk ≤ mus \n
    station = Point on the MobilizedBody expressed in the body frame.


AUTO UPDATE DISCRETE STATE\n
    p0 = Zero of the frictional spring.


CONTINUOUS STATE\n
    Sliding = Indicator of whether the MobilizedBody is sliding
    relative to the contact plane.


----------
PARAMETERS
----------
Customizable parameters specifying the characteristics of the exponential
spring are handled using ExponentialSpringParameters.

----
DATA
----
Calculated quantities are cached when the System is realized through
Stage::Dynamics. The cached data is made accessible via the helper class
ExponentialSpringData. */
class SimTK_SIMBODY_EXPORT ExponentialSpringForce : public ForceSubsystem {
public:
    /** Construct an exponential spring force object with default values for
    parameters like elasticity, viscosity, etc.
    @param system The system being simulated or studied.
    @param contactPlane Transform specifying the location and orientation of
    the contact plane. The positive y-axis defines the normal of the contact
    plane; friction forces occur in the x-z plane.
    @param body Body that will interact / collide with the contact plane.
    @param station Point on the specified body at which the contact force
    will be applied. The position and velocity of this point relative to
    the contact plane determine the magnitude and direction of the contact
    force.
    @param mus Initial value of the static coefficient of friction.
    0.0 ≤ mus ≤ 1.0
    @param muk Initial value of the kinetic coefficient of friction.
    0.0 ≤ muk ≤ mus */
    ExponentialSpringForce(MultibodySystem& system,
        const Transform& contactPlane,
        const MobilizedBody& body,const Vec3& station,
        Real mus, Real muk);

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
    @param mus Static coefficient of friction.   0.0 ≤ mus ≤ 1.0
    @param muk Kinetic coefficient of friction.  0.0 ≤ muk ≤ mus
    @param params Customized parameters. */
    ExponentialSpringForce(MultibodySystem& system,
        const Transform& contactPlane,
        const MobilizedBody& body, const Vec3& station,
        Real mus, Real muk, const ExponentialSpringParameters& params);

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

    /** Set the static coefficient of friction (mus) for this exponential
    spring.
    The value of mus is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, mus can be set at any
    time during a simulation. A change to mus will invalidate the System at
    Stage::Dynamics.
    @param state State object that will be modified.
    @param mus New value of the static coefficient of friction.
    0.0 ≤ mus ≤ 1.0 */
    void setMuStatic(State& state, const Real& mus);

    /** Get the static coefficient of friction (mus) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve mus. */
    const Real& getMuStatic(const State& state) const;
 
    /** Set the kinetic coefficient of friction (muk) for this exponential
    spring.
    The value of muk is held in the System's State object. Unlike the
    parameters managed by ExponentialSpringParameters, muk can be set at any
    time during a simulation. A change to muk will invalidate the System at
    Stage::Dynamics. @param state State object that will be modified.
    @param muk New value of the kinetic coefficient of friction.
    0.0 ≤ muk ≤ mus */
    void setMuKinetic(State& state, const Real& muk);

    /** Get the kinetic coefficient of friction (muk) held by the specified
    state for this exponential spring.
    @param state State object from which to retrieve muk. */
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
#endif SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

