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
 ------------------------------------------------------------------------------*/

#include "Simbody.h"

namespace SimTK {

//=============================================================================
/** Class ExponentialSpringParameters manages the parameters used to specify
the characteristics of an ExponentialSpringForce instance.  Changing a
parameter will invalidate the system at Stage::Topology, although the
invalidation will not occur until the parameters are set on an actual
instance of ExponentialSpringForce.

The force in the normal direction is computed using an exponential:

       fNormal = d1 * exp(-d2*(py-d0)) * (1.0 - kvNorm)

where py is the distance of the spring point on the Body above the plane
of the floor expressied in the frame of the contact plane.

The friction force is modeled using a linear spring with damping:

       fFric = -kpFric*(pxz-p0) - kvFric*vxz,

where pxz is the position of the spring on the Body projected onto the plane
of the floor and expressed in the floor frame, and p0 is the spring zero
in the plane of the floor.

d0 is the only parameter that can be positive, negative, or zero.  All other
parameters should be positive (and sometimes 0.0).

Note that there are 4 additional quantities that determine the
behavior of an exponential spring that are not handled by this class.
They are p0, mus, muk, and sliding.  These quantities are handled as
as states in order to allow them to change during the course of a simulation.
See ExponentialSpringForce for details on these states.

@see ExponentialSpringForce::setParameters()
@see ExponentialSpringForce */
class SimTK_SIMBODY_EXPORT ExponentialSpringParameters {
protected:
    /** d0 shifts the exponential function up and down with respect to
    the floor plane.Its default value is 0.0065905 (~7 mm above the floor
    plane). This slight upward shift prevents significant penetration
    into the floor plane.  d0 can have a value of 0.0. */
    Real d0;

    /** d1 linearly scales the applied force up or down. Its default value
    is 0.5336. */
    Real d1;

    /** d2 multiplies the exponent.  Its default value is 1150.0. Larger values
    of d2 makes the exponential curve rise more rapidly. */
    Real d2;

    /** kvNorm is the viscosity of the spring in the normal direction. Its
    default value is 0.5. kvNorm can have a value of 0.0. */
    Real kvNorm;

    /** kpFric is the elasticity of the friction spring. Its default value
    is 2000. */
    Real kpFric;

    /** kvFric is the viscosity of the friction spring. By default, its value
    is set relative to kpFric so that critical damping would occur for a
    specified mass (i.e., kvFric = 2 * sqrt(kpFric * mass)).
    kvFric can be 0.0. */
    Real kvFric;

    /** kTau is equal to 1.0/tau. tau is the characteristic time to transition
    between the static and kinetic coefficents of friction. The user sets tau.
    kTau is maintained as the member variable to avoid the cost of
    frequently dividing by tau.  The default value of tau is 0.01 seconds.
    @see setSlidingTimeConstant() */
    Real kTau;

    /** vSettle is the velocity of pxz below which the friction spring
    transitions to using the static coefficient of friction. Its default value
    is 0.01 m/s. */
    Real vSettle;

public:
    /** Default Constructor */
    ExponentialSpringParameters();

    /** Copy constructor.
    @param params Const reference to the object to be copied. */
    ExponentialSpringParameters(const ExponentialSpringParameters& params);

    /** Assignment operator.
    @param params Const reference to the object to which this instance should
    be assigned. */
    ExponentialSpringParameters& operator =
        (const ExponentialSpringParameters& params);

    /** Set the parameters conrolling the shape of the exponential.
    @param d0 shifts the exponential function up and down with respect to
    the floor plane. Its default value is 0.0065905 (~7 mm above the floor). 
    This slight upward shift eliminates significant penetration into the floor.
    @param d1 linearly scales the applied force up or down. Its default
    value is 0.5336.
    @param d2 multiplies the exponent. Its default value is 1150.0. Larger
    values of d2 make the exponential curve rise more rapidly as py gets less
    than d0. */
    void setShapeParameters(Real d0, Real d1 = 0.5336, Real d2 = 1150.0);

    /** Get the paramters that control the shape of the exponential.
    @see setShapeParameters() */
    void getShapeParameters(Real& d0, Real& d1, Real& d2) const;

    /** Set the viscosity of the exponential part of the spring. This viscosity
    only applies to the component of the velocity of the body spring point
    that is normal to the floor plane.
    @param kvNorm Viscosity of the spring in the normal direction. Its default
    value is 0.5. */
    void setNormalViscosity(Real& kvNorm);

    /** Get the viscosity of the exponential part of the spring.
    @returns Spring viscosity in the normal direction */
    Real getNormalViscosity() const;

    /** Set the elasticity of the friction spring and compute and set a
    viscosity that will result in critical damping for the specified mass
    (i.e., kv = 2*sqrt(kp*mass)).
    @param kp Elasticity of the friction spring.
    @param mass Mass of the body for which critical damping would be achieved.
    Articulated bodies effectively don't have a constant mass as it relates to
    acceleration in a particular direction, so think of this mass as a kind
    of average effective mass. A default mass of 1.0 kg is used if a mass
    is not specified. */
    void setElasticityAndComputeViscosity(Real kp, Real mass = 1.0);

    /** Set the elasticity of the friction spring. It's default value is 2000.
    @param kp Elastcity of the friction spring */
    void setElasticity(Real kp);

    /** Set the vicosity of the friction spring.  Its value can be set to 0.0.
    @param kv Viscosity of the friction spring */
    void setViscosity(Real kv);

    /** Get the elasticity of the friction spring.
    @returns Elasticity of the friction spring */
    Real getElasticity() const;

    /** Get the viscosity of the friction spring.
    @returns Viscosity of the friction spring */
    Real getViscosity() const;

    /** Set the time constant for transitioning back and forth between the
    static and kinetic coefficents of friction.  The transition is mediated
    by a rising or falling exponential that is asymptotic to mu static or
    or mu kinetic respectively. The default value of tau is 0.01 s.
    @param tau Time constant for sliding transitions */
    void setSlidingTimeConstant(Real tau);

    /** Get the time constant for transitioning back and forth between the
    static and kinetic coefficents of friction.
    @returns time constant for sliding transitions
    @see setSlidingTimeConstant() */
    Real getSlidingTimeConstant() const;

    /** Set the velocity below which the coefficient of friction transitions
    to the static coefficient of friction. It's default value is 0.01 m/s.
    @param vSettle Settle velocity */
    void setSettleVelocity(Real vSettle);

    /** Get the velocity below which the coefficient of friction transitions
    to the static coefficient of friction.
    @returns Settle velocity */
    Real getSettleVelocity() const;

    /** The member variables managed by this class are used by the underlying
    implementation of the %ExponentialSpringForce Subsystem.  Direct access
    is given to %ExponentialSpringForceImpl for reasons of speed and
    convenience.  Direct access is not given to other classes to ensure that
    these paramters are set only to valid values. */
    friend class ExponentialSpringForceImpl;
};


//=============================================================================
/** ExponentialSpringData is a helper class that is used to store key
quantities associated with the ExponentialSpringForce Subsystem during a
simulation.  An instance of this class serves as the data Cache Entry for the
ExponentialSpringForceImpl Subsystem. All of its member variables are
guarranteed to be calculated and set once the System has been realized
to the Dynamics Stage.

To understand what the quantities organized by this class represent, a basic
description of the contact problem that is solved, along with a description
of coordinate frame conventions, will be helpful.

Class ExponentialSpringForce computes and applies a contact force at a
specified point on a MoblizedBody (i.e., a Station) due to interaction of
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
these quantities pertain to the normal of the contact plane.  Member
variables with a "xz" suffix (e.g., pxz, vxz, or fxz) indicate that
these quanties lie in the contact plane (tangent to it) and are
associated with the friction force.

Member variables with a "_G" suffix are expressed in the ground frame. Member
variables without a "_G" suffix are expressed in frame of the contact plane. 

Note- It is recognized that a class whose member variables are public does
not generally conform to good coding practices. In the case of this class,
however, the typical user will only have access to a const reference to
the instance owned by the State and, thus, the user will not be able to
alter the member variables (unless extraordinary measures are taken). Only
the underlying implementation (ExponentialSpringForceImpl) has writable
access to the instance of this class that is owned by the State.

A copy constructor and assignment operator are provided in the event a user
wants to create a copy of the data that is not owned by the State.

@see ExponentialSpringForce */
class SimTK_SIMBODY_EXPORT ExponentialSpringData {
public:
    /** Default Constructor. */
    ExponentialSpringData();

    /** Copy Constructor.
    @param data Const refernce to the object that should be copied. */
    ExponentialSpringData(const ExponentialSpringData& data);

    /** Assignment Operator.
    @param data Const reference to the object to which this instance should
    be assigned. */
    ExponentialSpringData& operator = (const ExponentialSpringData& data);

    /** Position of the body spring station in the ground frame. */
    Vec3 p_G;

    /** Velocity of the body spring station in the ground frame. */
    Vec3 v_G;

    /** Position of the body spring station expressed in the frame of the
    contact plane. */
    Vec3 p;

    /** Velocity of the body spring station expressed in the frame of the
    contact plane. */
    Vec3 v;

    /** Displacment of the body spring station normal to the floor expressed
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
    force, but is likely desired for vizualization. */
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

    /** Elastic frictional force expressed in the frame of the contact plane. */
    Vec3 fricElas;

    /** Damping frictional force expressed in the frame of the contact plane. */
    Vec3 fricDamp;

    /** Total frictional force (elastic + damping) expressed in the frame of
    the contact plane. */
    Vec3 fric;

    /** Magnitude of the frictional force. */
    Real fxy;

    /** Resultant spring force (normal + friction) expressed in the floor
    frame. */
    Vec3 f;

    /** Resultant spring force (normal + frcition) expressed in the ground
    frame. */
    Vec3 f_G;
};


//=============================================================================
/** Class ExponentialSpringForce implements an 'exponential spring' as a means
of modeling contact of a point on a MobilizedBody with a contact plane. In
practice, a number of body-fixed points ('Stations' in Simbody vocabulary) are
strategically chosen across the surface of the MobilizedBody so as to represent
the geometry of that body reasonably well, and an ExponentialSpringForce is
created for each of those points.  For example, if the body were a cube, one
would likely choose to place an exponential spring at each corner of the cube.
The contact plane is typically used to model interactions with a floor, but
need not be limited to this use case. The contact plane can be rotated and
displaced relative to the ground frame and so can be used to model a wall or
ramp, for example.

A distinguishing feature of the exponential spring, relative to other
contact models, is that it ALWAYS applies a force to the body; there is
never a time in a simulation when the spring force is not applied.  This
seemingly non-physical feature works because the force becomes small
(less than 0.01 N) as soon as the body station is more than 1.0 cm above
the ground and extremely small (less than 0.000001 N) as soon as the body
station is more than 2.0 cm above the ground. The sizes of these forces is
below the precision with which the acceleration due to gravity
(e.g., g = (0.,-9.8,0.)) is typically specified in a simulation. The sizes
of these is also small compared to the errors made in neglecting air
resistance or to errors made in the estimation of inertial properties.

An advantage of this distinguishing feature is that there is no need to
search for intersections (the number active of springs does not change
during a simulation) and there is no need to find the precise time of
contact.  In addition, because the function describing the contact force
is smooth, adjustments in integration step size made by variable-step
explicit integrators are often well-behaved.

As a side note, electrostatic forces, which are the fundamental force from
which all conventional contact forces arise, are never actually zero either.
They obey Coulomb's Law (Fe = ke*q1*q2/r^2) which states that an
electrostatic force only approaches zero as the distance between two
charges (r) approaches infinity.  It is certainly not the claim here that
Exponential Springs accurately models contact force at the level of
fundamental forces; it is only suggested that the use of a force field
that acts over a great distance is not that far fetched as long as the
force gets sufficiently small sufficiently quickly.

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

        fyDamping = - kv*vy * fElastic,

where vy is the normal component of the velocity of the specified Station.
All together, the spring force in the normal direction is given by

        fy  = fElastic + fDamping
            = d1*exp(d2*(py-d0)) - kv*vy * d1*exp(d2*(py-d0))

Friction %Force (x-z plane) \n
Friction force (type Vec3) is implemented using a damped linear spring.
The elastic and viscous terms are given by

        fxyElas = -kp*(pxy-p0)

        fxyDamp = -kv*(vxy)

and the total friction force is given by

        friction = fxyElas + fxyDamp

where kp is the spring elasticity, kv is the spring viscosity, pxy is the
position of the spring Station projected onto the contact plane, p0 is the
current spring zero, which always resides in the contact plane, and vxy is
the velocity of the Station in the contact plane.

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

       if(fxyElas.norm() > fxyLimit)  p0New = pxy + (pxy-p0)/kp

In Simbody, p0 is handled as an Auto Update Discrete State.  Any change
to p0 is made to the Update Cache (not to the State directly), and the
integrator swaps this cache value with the actual p0 State after a
successful integration step is obtained.

Coefficients of Friction \n
Coefficients of kinetic (sliding) and static (fixed) friction can be
specified for the spring subject to the following constraints:

       0.0 <= muk <= mus <= 1.0

A continuous state variable (Z = Sliding) is used to characterize the
sliding state of the spring.

        Sliding = 0     means fully fixed in place (lower bound)
        Sliding = 1     means fully sliding (upper bound)

The instantaneous coefficient of friction (mu) is calculated based on the
value of Sliding:

        mu = mus - Sliding*(mus - muk)

The time derivative of Sliding is used to drive Sliding toward the
extreme of 0.0 or 1.0, depending on the following criteria...

When the frictional force exceeds its limit, Sliding is driven toward 1:

        if (friction.norm() >= fxyLimit)  SlidingDot = (1.0 - Sliding)/tau

When vxy falls below some specified magnitude (e.g., 0.01 m/s) AND the
frictional force is less than its limit, Sliding is driven toward 0:

        else if (vxy.norm < 0.01)  SlidingDot = -Sliding/tau

Otherwise, no change to the Sliding state is made:

        else  SlidingDot = 0.0

In the above equations, tau is the characteristic rate at which the Sliding
rises or decays. The motivation for using a continuous state variable for
Sliding is that, although the transition between fixed and sliding may happen
quickly, it does not happen instantaneously.  And, modeling Sliding based on a
differential equation ensures that mu is continuous.

Future Enhancements \n
Future enhancements might include the capacity to
1) use a polygonal mesh as the contact plane, and
2) move the specified MobilizedBody Station in manner that adapts to the 
contact circumstances. For example, to model a coin rolling on a table
a small number of stations could be continually moved to the portion of the
coin that is closest to the table.

------
STATES
------
DISCRETE STATES\n
@param mus static coefficient of friction.  0.0 <= mus <= 1.0
@param muk kinetic coefficient of friction.  0.0 <= muk <= mus
@param station Point on the MobilizedBody expressed in the body frame.

AUTO UPDATE DISCRETE STATE\n
@param p0 Zero of the frictional spring.

CONTINUOUS STATE\n
@param Sliding Characterizes whether the spring zero is sliding (1.0) or fixed
in place (0.0). This state is used for transitioning back and forth between
mus and muk.

----------
PARAMETERS
----------
Parameters specifying the characteristics of the exponential spring are
handled using ExponentialSpringParameters.

----
DATA
----
Calculated quantities are cached during the Dynamics Stage.  The cached
data is made accessible by the helper class ExponentialSpringData. */
class SimTK_SIMBODY_EXPORT ExponentialSpringForce : public ForceSubsystem {
public:
    ExponentialSpringForce(MultibodySystem& system,
        const Transform& floor,const MobilizedBody& body,const Vec3& station);
    ExponentialSpringForce(MultibodySystem& system, const ExponentialSpringParameters& params,
        const Transform& floor, const MobilizedBody& body, const Vec3& station);
    // Parameters (set these before realizing at the Topology Stage)
    void setParameters(const ExponentialSpringParameters& params);
    const ExponentialSpringParameters& getParameters() const;
    // Coefficients of Friction (discrete states - can change during a simulation)
    // Static
    void setMuStatic(State& state, const Real& mus);
    const Real& getMuStatic(const State& state) const;
    // Kinetic
    void setMuKinetic(State& state, const Real& muk);
    const Real& getMuKinetic(const State& state) const;
    /// Reset the spring zero (auto update discrete state) so that it
    /// coincides in the floor plane with the specified station on the body.
    /// This is often done at the beginning of a simulation to ensure that
    /// a simulation does not begin with large friction forces.
    void resetSprZero(State& state) const;
    /// Retrieve access to a const reference to the data cache.  The System
    /// should be realized through Stage::Dynamics before a call to this method;
    /// an exception will be thrown otherwise.
    const ExponentialSpringData& getData(const State& state) const;
    // Shortcuts for accessing the underlying implementation
    ExponentialSpringForceImpl& updImpl();
    const ExponentialSpringForceImpl& getImpl() const;
};


} // namespace SimTK
#endif SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

