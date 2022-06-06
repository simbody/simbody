#ifndef SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_
#define SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_

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

namespace SimTK {

//=============================================================================
/** ExponentialSpringParameters is a helper class used to customize the
force-producing characteristics of an ExponentialSpringForce instance.
ExponentialSpringParameters manages all parameters that are implemented as
Stage::Topology variables. To customize any of the Topology-stage parameters
on an ExponentialSpringForce instance, you should

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

4) Realize the system to Stage::Topology. When a new set of parameters is
set on an ExponentialSpringForce instance, as above in step 3, the System will
be invalidated at Stage::Topology.  The System must therefore be realized at
Stage::Topology (and hence Stage::Model) before a simulation can proceed.

        system.realizeTopology();

Note that each ExponentialSpringForce instance owns its own private
ExponentialSpringParameters object. The myParams object is just used to set
the desired parameter values of the privately owned parameters object. It is
fine for objects like myParams to go out of scope or for myParams objects
allocated from the heap to be deleted.

Therefore, also note that the parameter values possessed by an
ExponentialSpringForce instance do not necessarily correspond to the values
held by an instance of this class until a call to
ExponentialSpringForce::setParameters() is made.

The default values of the parameters held by ExponentialSpringParameters
work well for typical contact interactions, but clearly may not be
appropriate for simulating many contact interactions. For example, you might
want to simulate an interaction in which very little energy is dissipated
during a contact event, in which case you'd reduce the normal viscosity and
friction spring viscosity, as well as decrease the coefficients of
friction. Valid values of the friction spring elasticity and viscosity
(kxy and cxy) can range widely (e.g., kxy = 1,000 to kxy = 1,000,000)
depending on the material properties of the objects represented by a
particular MobilizedBody and contact plane (e.g., bare feet on a yoga mat vs.
a steel bearing on a marble floor). In general, higher values of kxy and cxy
will result in smaller integration step sizes.

The default values of the parameters are expressed in units of Newtons,
meters, seconds, and kilograms; however, you may use an alternate set of
self-consistent units by re-specifying all parameters.

Finally, there are 2 parameters in this class that possess a different
status than the others. They are 1) the initial value for the static
coefficient of friction (μₛ) and 2) the initial value for the kinetic
coefficient of friction (μₖ).

These two parameters, like the others, are used to set the force-producing
characteristics of an ExponentialSpringForce instance at the beginning of a
simulation when a model is put together. And, if the initial μₛ or the
initial μₖ are changed, those new values will not be realized for an
ExponentialSpringForce instance until ExponentialSpringForce::setParameters()
on the instance is called and, further, the System is re-realized to
Stage::Topology.

Unlike the other parameters, however, there are Dynamics-stage discrete state
variables that correspond to the initial μₛ and the default μₖ. These
discrete-state versions (think of them as the instantaneous μₛ and the
instantaneous μₖ) can be changed during the course of a simulation without
having to call setParameters() and re-realize the System at Stage::Toplogy.
This functionality allows a contact plane to posses non-uniform frictional
characteristics across its surface. For example, a patch of ice on a sidewalk
could be modeled. The discrete-state versions of μₛ and μₖ (the instantaneous
versions) should NOT be changed via the ExponentialSpringParameters class, but
by setting them directly on an ExponentialSpringForce instance by calling
ExponentialSpringForce::setMuStatic() and
ExponentialSpringForce::setMuKinetic().

In summary, regarding coefficients of friction, use
ExponentialSpringParameters to establish the initial values for μₛ and μₖ;
use ExponentialSpringForce::setMuStatic() and
ExponentialSpringForce::setMuKinetic() to change μₛ and μₖ during the course
of a simulation.

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
    ExponentialSpringParameters(const ExponentialSpringParameters& source)
        = default;

    /** Copy assignment operator.
    @param source Const reference to the source object to which this instance
    should be assigned. */
    ExponentialSpringParameters&
        operator=(const ExponentialSpringParameters& source) = default;

    /** Equality operator. All member variables of this object must be
    equal to the corresponding members of the other object, in order for the
    equality operator to return true.
    @param other The other object with which to test equality. */
    bool operator==(const ExponentialSpringParameters& other) const;

    /** Inquality operator. If any member of this object is not equal
    to its corresponding member in the other, true is returned. */
    bool operator!=(const ExponentialSpringParameters& other) const;

    /** Set the parameters that determine the shape of the exponential. The
    exponential, which models the elastic response of the spring in the normal
    direction, is a function of 3 parameters:

            fzElastic = d₁exp(−d₂(pz−d₀))

    Note that pz is the displacement of the body spring station above
    (pz > 0.0) or below (pz < 0.0) the contact plane. The default values of
    the shape parameters were chosen to maximize integration step size while
    maintaining a number of constraints (e.g., the normal force should fall
    below 0.01 Newtons when pz > 1.0 cm).
    @param d0 shifts the exponential function up and down with respect to the
    contact plane. Its default value is 0.0065905 m (~7 mm above the contact
    plane). This slight upward shift reduces penetration of the body spring
    station below the contact plane. That is, unless there is an extreme
    impact event, the repulsive force applied to the body will generally be
    large enough to keep pz from going negative. There is no issue with pz
    going negative (nothing special happens); the shift just facilitates an
    interpretation of the contact interaction that is conceptually appealing.
    d0 can be positive, have a value of 0.0, or be negative.
    @param d1 linearly scales the applied force up or down. Its default
    value is 0.5336 Newtons. d1 should be positive to generate a repulsive
    force directed along the positive z-axis of the contact plane.
    @param d2 linearly scales the exponent. Its default value is 1150.0/m.
    Larger values of d2 make the exponential curve rise more rapidly as pz
    gets small or becomes negative. d1 should be positive. */
    void setShapeParameters(Real d0, Real d1 = 0.5336, Real d2 = 1150.0);

    /** Get the parameters that determine the shape of the exponential. See
    setShapeParameters() for a detailed description of these parameters. */
    void getShapeParameters(Real& d0, Real& d1, Real& d2) const;

    /** Set the viscosity of the exponential part of an ExponentialSpringForce.
    This viscosity only multiplies the component of the velocity of the
    spring station that is normal to the contact plane. To eliminate energy
    dissipation due to motion in the normal direction, set cz equal to 0.0.
    @param cz Viscosity of a spring in the normal direction.
    Its default value is 0.5 N*s/m. cz should be positive or zero. */
    void setNormalViscosity(Real cz);

    /** Get the viscosity of the exponential part of an ExponentialSpringForce.
    The exponential part of the contact force is always normal to the contact
    plane.
    @returns Viscosity in the normal direction (cz) */
    Real getNormalViscosity() const;

    /** Set the maximum normal force allowed for an ExponentialSpringForce. If
    the normal force exceeds the maximum normal force, it is clamped at the
    maximum normal force before being applied to the MobilizedBody. Setting the
    maximum force can be important for avoiding numerical isses and integration
    failures for motions that involve extreme impact events (e.g., a mass
    striking the floor after falling from a high height).
    @param maxFz Maximum normal force. Its default value is 100,000.0 N.
    maxFz should be positive. */
    void setMaxNormalForce(Real maxFz);

    /** Get the maximum normal force allowed for an ExponentialSpringForce. */
    Real getMaxNormalForce() const;

    /** Set both the elasticity and viscosity of the damped linear spring
    used to model friction in Friction Model 2. The value that is set for the
    viscosity is computed so that critical damping would result for a
    specified mass (i.e., cxy = 2*sqrt(kxy*mass)). A call to this method
    overrides any values set previously by calls to setElasticity() or
    setViscosity().
    @param kxy Elasticity of the friction spring. kxy should be positive.
    @param mass Mass of the body for which critical damping would be achieved.
    Articulated bodies generally don't have an effective mass that is constant
    as it relates to acceleration in a particular direction, so think of this
    mass as a kind of average mass. A default mass of 1.0 kg is used if a mass
    is not specified. */
    void setElasticityAndViscosityForCriticalDamping(Real kxy, Real mass=1.0);

    /** Set the elasticity of the friction spring (kxy). A call to this method
    overrides any value of elasticity previously set by a call to
    setElasticityAndViscosityForCriticalDamping().
    @param kxy Elasticity of the friction spring. Its default value is
    20000.0 N/m. kxy should be positive. */
    void setFrictionElasticity(Real kxy);

    /** Get the elasticity of the friction spring. The value of the elasticity
    is the default value (20000.0 N/m) or the value set by a call to either
    setElasticity() or setElasticityAndComputeVicosity(), whichever was called
    most recently.
    @returns Elasticity of the friction spring. */
    Real getFrictionElasticity() const;

    /** Set the viscosity of the friction spring (cxy). A call to this method
    overrides any value of viscosity previously set by a call to
    setElasticityAndViscosityForCriticalDamping(). Setting the viscosity equal
    to 0.0 is fine. Be aware, however, that if a body is not sliding,
    setting cxy = 0.0 will simply allow the body to vibrate in place
    indefinitely. If a body is sliding, even if cxy = 0.0, the kinetic energy
    of the body will still be dissipated because the frictional force is
    directed opposite the sliding velocity, and the elastic part of the
    friction spring will not store additional potential energy because the
    friction spring zero (p₀) is continually released. The only way to
    eliminate energy dissipation entirely is to also set the coefficients of
    friction equal to 0.0.
    @param cxy Viscosity of the friction spring. Its default value is
    2.0*sqrt(kxy*mass) = 2.0*sqrt(20000*1) ~= 282.8427 N*s/m. cxy should be
    0.0 or positive. */
    void setFrictionViscosity(Real cxy);

    /** Get the viscosity of the friction spring. The value of the viscosity
    is the default value (~282.8427 N*s/m) or the value set by a call to either
    setViscosity() or setElasticityAndViscosityForCriticalDamping(), whichever
    was called most recently.
    @returns Viscosity of the friction spring. */
    Real getFrictionViscosity() const;

    /** Set the time constant for transitioning back and forth between the
    static (μₛ) and kinetic (μₖ) coefficients of friction. The transition is
    mediated by a rising or decaying exponential that is asymptotic to μₛ or
    μₖ, respectively.
    @param tau Time constant (τ) for sliding transitions. The default value of
    τ is 0.01 s. τ must be positive. */
    void setSlidingTimeConstant(Real tau);

    /** Get the time constant for transitioning back and forth between the
    static (μₛ) and kinetic (μₖ) coefficients of friction.
    @returns Time constant for sliding transitions. */
    Real getSlidingTimeConstant() const;

    /** Set the velocity below which the coefficient of friction transitions
    to the static coefficient of friction (μₛ).
    @param vSettle Settle velocity. It's default value is 0.01 m/s. vSettle
    must be positive. */
    void setSettleVelocity(Real vSettle);

    /** Get the velocity magnitude below which the coefficient of friction
    transitions to the static coefficient of friction.
    @returns Settle velocity */
    Real getSettleVelocity() const;

    /** Set the acceleration magnitude below which the coefficient of friction
    transitions to the static coefficient of friction (μₛ).
    @param aSettle Settle acceleration. It's default value is 1.0 m/s².
    aSettle must be positive.*/
    void setSettleAcceleration(Real aSettle);

    /** Get the acceleration magnitude below which the coefficient of friction
    transitions to the static coefficient of friction.
    @returns Settle acceleration. */
    Real getSettleAcceleration() const;

    /** Set the initial value of the static coefficient of friction (μₛ).
    Note- use ExponentialSpringForce::setMuStatic(), not this method, to set
    the instantaneous value of μₛ during the course of a simulation.
    @param mus Initial static coefficient of friction (μₛ). It's default
    value is 0.7 (unitless). 0.0 ≤ μₖ ≤ μₛ. If μₛ < μₖ, μₖ is set to μₛ. */
    void setInitialMuStatic(Real mus);

    /** Get the initial value of the static coefficient of friction (μₛ).
    @returns Initial μₛ. */
    Real getInitialMuStatic() const;

    /** Set the initial value of the kinetic coefficient of friction (μₖ).
    Note- use ExponentialSpringForce::setMuKinetic(), not this method, to set
    the instantaneous value of μₖ during the course of a simulation.
    @param muk Initial kinetic coefficient of friction (μₖ). It's default
    value is 0.5 (unitless). 0.0 ≤ μₖ ≤ μₛ. If μₖ > μₛ, μₛ is set to μₖ.*/
    void setInitialMuKinetic(Real mus);

    /** Get the initial value of the kinetic coefficient of friction (μₖ).
    @returns Initial μₖ. */
    Real getInitialMuKinetic() const;

private:
    Real d0;
    Real d1;
    Real d2;
    Real cz;
    Real maxFz;
    Real kxy;
    Real cxy;
    Real tau;
    Real vSettle;
    Real aSettle;
    Real initMus;
    Real initMuk;
};

}       // namespace SimTK
#endif  // SimTK_SIMBODY_EXPONENTIAL_SPRING_PARAMETERS_H_

