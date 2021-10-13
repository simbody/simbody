#ifndef SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_
#define SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_

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
static coefficient of friction (μₛ) and the kinetic coefficient of friction
(μₖ). μₛ and μₖ are parameters, but they are not implemented as
Topology-stage variables; they are implemented as Dynamics-stage discrete
state variables. Thus, they can be changed during the course of a simulation
without invalidating the System at Stage::Topology. This functionality allows
a contact plane to posses non-uniform frictional characteristics across its
surface. For example, a patch of ice on a sidewalk could be modeled. μₛ and
μₖ can be set during a simulation using ExponentialSpringForce::setMuStatic()
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
    specified mass (i.e., kᵥ = 2*sqrt(kₚ*mass)). A call to this method
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
    sliding, setting kᵥ = 0.0 will simply allow the body to vibrate in place
    indefinitely. If a body is sliding, even if kᵥ = 0.0, the kinetic energy
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
    @param τ Time constant for sliding transitions. The default value of τ
    is 0.01 s. τ must be positive. */
    void setSlidingTimeConstant(Real τ);

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

private:
    Real d0;
    Real d1;
    Real d2;
    Real kvNorm;
    Real kpFric;
    Real kvFric;
    Real tau;
    Real vSettle;

};


} // namespace SimTK
#endif SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_

