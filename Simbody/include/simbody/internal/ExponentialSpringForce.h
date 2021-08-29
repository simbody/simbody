#ifndef SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_
#define SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * Copyright (c) 2021 Authors.                                                *
 * Authors: Frank C. Anderson                                                 *
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

#include "Simbody.h"
#include "simbody/internal/ForceSubsystemGuts.h"

namespace SimTK {

//=============================================================================
/** Class ExponentialSpringParameters manages the parameters used to specify the
* characteristics of an ExponentialSpringForce instance.  Changing a parameter
* will invalidate the system at %Stage::Topology, although the invalidation
* will not occur until the parameters are set on an actual instance of
* %ExponentialSpringForce.
*
* The force in the normal direction is computed using an exponential:
*
*		fNormal = d1 * exp(-d2*(py-d0)) * (1.0 - kvNorm)
*
* where py is the distance of the spring point on the Body above the plane
* of the floor expressied in the floor frame.
*
* The friction force is modeled using a linear spring with damping:
*
*		fFric = -kpFric*(pxz-p0) - kvFric*v,
*
* where pxz is the position of the spring on the Body projected onto the plane
* of the floor and expressed in the floor frame, and p0 is the spring zero
* in the plane of the floor.
*
* d0 is the only parameter that can be positive, negative, or zero.  All other
* parameters should be positive (and sometimes 0.0).
*
* Note that there are 4 additional quantities that determine the
* behavior of an exponential spring that are not handled by this class.
* They are p0, mus, muk, and sliding.  These quantities are handled as
* as states in order to allow them to change during the course of a simulation.
* See ExponentialSpringForce for details on these states.
*
* @see ExponentialSpringForce::setParameters()
* @see ExponentialSpringForce
*/
class SimTK_SIMBODY_EXPORT ExponentialSpringParameters
{
protected:
	/// <summary>
	/// d0 shifts the exponential function up and down with respect to
	/// the floor plane.Its default value is 0.0065905 (~7 mm above the floor
	/// plane). This slight upward shift prevents significant penetration
	/// into the floor plane.  d0 can have a value of 0.0.
	/// </summary>
	Real d0;

	/// <summary>
	/// d1 linearly scales the applied force up or down. Its default value is 0.5336.
	/// </summary>
	Real d1;

	/// <summary>
	/// d2 multiplies the exponent.  Its default value is 1150.0. Larger values
	/// of d2 makes the exponential curve rise more rapidly.
	/// </summary>
	Real d2;

	/// <summary>
	/// kvNorm is the viscosity of the spring in the normal direction. Its default
	/// value is 0.5. kvNorm can have a value of 0.0.
	/// </summary>
	Real kvNorm;

	/// <summary>
	/// kpFric is the elasticity of the friction spring. Its default value is 2000.
	/// </summary>
	Real kpFric;

	/// <summary>
	/// kvFric is the viscosity of the friction spring. By default, its value is
	/// set relative to kpFric so that critical damping would occur for a specified
	/// mass (i.e., kvFric = 2 * sqrt(kpFric * mass)).  kvFric can be 0.0.
	/// </summary>
	Real kvFric;

	/// <summary>
	/// kTau is equal to 1.0/tau. tau is the characteristic time to transition
	/// between the static and kinetic coefficents of friction.  The user sets tau.
	/// kTau is maintained as the member variable to avoid the cost of
	/// frequently dividing by tau.  The default value of tau is 0.01 seconds.
	/// @see setSlidingTimeConstant()
	/// </summary>
	Real kTau;

	/// <summary>
	/// vSettle is the velocity of pxz below which the friction spring transitions
	/// to using the static coefficient of friction. Its default value is 0.01 m/s.
	/// </summary>
	Real vSettle;

public:
	// Constructors and operators
	/// Default Constructor
	ExponentialSpringParameters();

	/// <summary>
	///  Copy constructor.
	/// </summary>
	/// <param name="params"> Reference to the object to be copied. </param>
	ExponentialSpringParameters(const ExponentialSpringParameters& params);

	/// <summary>
	/// Assignment operator.
	/// </summary>
	void operator = (const ExponentialSpringParameters& params);

	/// <summary>
	/// Set the parameters conrolling the shape of the exponential.
	/// </summary>
	/// <param name="d0">
	/// d0 shifts the exponential function up and down with respect to the floor
	/// plane. Its default value is 0.0065905 (~7 mm above the floor plane). 
	/// This slight upward shift prevents significant penetration into the floor.
	/// </param>
	/// <param name="d1">
	/// d1 linearly scales the applied force up or down. Its default value is 0.5336.
	/// </param>
	/// d2 multiplies the exponent. Its default value is 1150.0. Larger values of
	/// d2 make the exponential curve rise more rapidly as py gets less than d0.
	/// <param name="d2"></param>
	void setShapeParameters(Real d0, Real d1 = 0.5336, Real d2 = 1150.0);

	/// <summary>
	/// Get the paramters that control the shape of the exponential.
	/// </summary>
	void getShapeParameters(Real& d0, Real& d1, Real& d2) const;

	/// <summary>
	/// Set the viscosity of the exponential part of the spring.  This viscosity
	/// only applies to the component of the velocity of the body spring point
	/// that is normal to the floor plane.
	/// </summary>
	/// <param name="kvNorm">
	/// Viscosity of the spring in the normal direction. Its default value is 0.5.
	/// </param>

	void setNormalViscosity(Real& kvNorm);
	/// <summary>
	/// Get the viscosity of the exponential part of the spring.
	/// </summary>
	/// <returns>normal viscosity</returns>
	Real getNormalViscosity() const;

	/// <summary>
	/// Set the elasticity of the friction spring and compute and set a viscosity
	/// that will result in critical damping for the specified mass
	/// (i.e., kv = 2*sqrt(kp*mass)), which is done by default.
	/// Friction forces lie in the plane tangent to the floor.
	/// </summary>
	/// <param name="kp">Elasticity of the friction spring.</param>
	/// <param name="mass">
	/// Mass of the body for which critical damping would be achieved.  Articulated
	/// bodies effectively don't have a constant mass as it relates to acceleration in
	/// a particular direction, so think of this mass as a kind of average
	/// effective mass.  A default mass of 1.0 kg is used if a mass is not specified.
	/// </param>
	void setElasticityAndComputeViscosity(Real kp, Real mass = 1.0);

	/// <summary>
	/// Set the elasticity of the friction spring.  It's default value is 2000.0.
	/// </summary>
	/// <param name="kp">Elastcity of the frictgion spring.</param>
	void setElasticity(Real kp);

	/// <summary>
	/// Set the vicosity of the friction spring.  Its value can be set to 0.0.
	/// </summary>
	/// <param name="kv">Viscosity of the friction spring.</param>
	void setViscosity(Real kv);

	/// <summary>
	/// Get the elasticity of the friction spring.
	/// </summary>
	/// <returns>Elasticity of the friction spring.</returns>
	Real getElasticity() const;

	/// <summary>
	/// Get the viscosity of the friction spring.
	/// </summary>
	/// <returns>Viscosity of the friction spring.</returns>
	Real getViscosity() const;

	/// <summary>
	/// Set the time constant for transitioning back and forth between
	/// the static and kinetic coefficents of friction.  The transition
	/// is mediated by a rising or falling exponential that is asymptotic to
	/// mu static or mu kinetic respectively.
	/// </summary>
	/// <param name="tau">time constant</param>
	void setSlidingTimeConstant(Real tau);

	/// <summary>
	/// Get the time constant for transitioning back and forth between
	/// the static and kinetic coefficents of friction.  The default value is
	/// 0.01 seconds.
	/// </summary>
	/// <returns>time constant</returns>
	Real getSlidingTimeConstant() const;

	/// <summary>
	/// Set the velocity below which the coefficient of friction transitions
	/// to the static coefficient of friction.  It's default value is 0.01 m/s.
	/// </summary>
	/// <param name="vSettle">Settle velocity</param>
	void setSettleVelocity(Real vSettle);

	/// <summary>
	/// Get the velocity below which the coefficient of friction transitions
	/// to the static coefficient of friction.
	/// </summary>
	/// <returns>Settle velocity</returns>
	Real getSettleVelocity() const;

	/// <summary>
	/// The member variables managed by this class are used by the underlying
	/// implementation of the %ExponentialSpringForce Subsystem.  Direct access
	/// is given to %ExponentialSpringForceImpl for reasons of speed.  Direct access
	/// is not given to other classes to ensure that these paramters are set only
	/// to valid values.
	/// </summary>
	friend class ExponentialSpringForceImpl;
};


//=============================================================================
/** Class ExponentialSpringData is use to store key quantities associated with
* the %ExponentialSpringForce Subsystem during a simulation.  An instance
* of this class serves as the Cache Entry for the
* %ExponentialSpringForceImpl Subsystem.
*
* All of its member variables are guarranteed to be calculated and set once
* the System has been realized to the Dynamics Stage.
*
* At this time, contact for points specified on a body occurs only with a
* plane representing a "floor".  The floor plane can be rotated and
* displaced relative to the ground frame.
*
* The positive y axis of the floor frame defines the normal of the floor
* plane.  It is the axis along which the spring force is modeled as an
* exponential.
*
* The x and z axes of the floor frame lie in the plane of the floor.
* The friction force lies in the floor plane.
*
* A copy constructor and assignment operator are provided for this class
* in case a user wants to create a copy of the data that is not owned
* by the State.
*
* Member variables with an _G suffix are expressed in the ground frame.
* Member veriables without a suffix are expressed in the floor frame.
*/
class SimTK_SIMBODY_EXPORT ExponentialSpringData
{
public:
	// Default Constructor
	ExponentialSpringData();

	// Copy Constructor
	ExponentialSpringData(const ExponentialSpringData& data);

	// Assignment Operator
	void operator = (const ExponentialSpringData& data);

	// Members ----
	// Position of the body spring station in the ground frame.
	Vec3 p_G;

	// Velocity of the body spring station in the ground frame.
	Vec3 v_G;

	// Position of the body spring station expressed in the floor frame.
	Vec3 p;

	// Velocity of the body spring station expressed in the floor frame.
	Vec3 v;

	// Displacment of the body spring station normal to the floor
	// expressed in the floor frame.
	Real py;

	// Velocity of the body spring station normal to the floor
	// expressed in the floor frame.
	Real vy;

	// Position of the body spring station projected onto the floor
	// expressed in the floor frame.
	Vec3 pxz;

	// Velocity of the body spring station in the plane of the floor
	// expressed in the floor frame.
	Vec3 vxz;

	// Position of the body spring station projected onto the floor
	// expressed in the ground frame.  This quantity is not needed to
	// calculate the spring force, but is likely desired for vizualization.
	Vec3 pxz_G;

	// Elastic force in the normal direction.
	Real fyElas;

	// Damping force in the normal direction.
	Real fyDamp;

	// Total normal force.
	Real fy;

	// Instantaneous coefficient of friction.
	Real mu;

	// Limit of the frictional force.
	Real fxyLimit;

	// Elastic frictional force.
	Vec3 fricElas;

	// Damping frictional force.
	Vec3 fricDamp;

	// Resultant frictional force
	Vec3 fric;

	// Magnitude of the frictional force
	Real fxy;

	// Total spring force
	Vec3 f;

	// Total spring force expressed in the ground frame.
	Vec3 f_G;
};


//=============================================================================
/** Class ExponentialSpringForce implements an exponential spring as a means of
* modeling contact of a body with a floor.  At present, the floor is modeled
* simply as a horizontal plane.  In practice, a number of body-fixed points
* ('Stations' in Simbody vocabulary) are strategically chosen across the
* surface of the body so as to represent the geometry of that body reasonably
* well, and an exponential spring force is applied at each of those points.
* For example, if the body were a cube, one would likely choose to place a
* spring at each corner of the cube.
*
* A distinguishing feature of the exponential spring, relative to other
* contact models, is that it ALWAYS applies a force to the body; there is
* never a time in a simulation when the spring force is not applied.  This
* seemingly non-physical feature works because the force becomes small
* (less than 0.01 N) as soon as the body station is more than 1.0 cm above
* the ground and extremely small (less than 0.000001 N) as soon as the body
* station is more than 2.0 cm above the ground.  The sizes of these forces
* is well below the accuracy with which the acceleration due to gravity
* (e.g., g = (0.,-9.8,0.)) is typically specified in a simulation.
*
* An advantage of this distinguishing feature is that there is no need to
* search for intersections (the number active of springs does not change
* during a simulation) and there is no need to find the precise time of
* contact.  In addition, because the function describing the contact force
* is smooth, adjustments in integration step size made by variable-step
* explicit integrators are often well-behaved.
*
* As a side note, electrostatic forces, which are the fundamental force from
* which all conventional contact forces arise, are never actually zero either.
* They obey Coulomb's Law (Fe = ke*q1*q2/r^2) which states that an
* electrostatic force only approaches zero as the distance between two
* charges (r) approaches infinity.  It is certainly not the claim here that
* Exponential Springs accurately models contact force at the level of
* fundamental forces; it is only suggested that the use of a force field
* that acts over a great distance is not that far fetched as long as the
* force gets sufficiently small sufficiently quickly.
*
* Normal Force \n
* The elastic part of the normal force is computed using an exponential
* whose shape is a function of three parameters (d0, d1, and d2):
*
*	fElastic = d1 * exp(d2 * (py - d0))
*
* where py is the height above the floor plane of the specified Station
* on the Body. The default values for the parameters are
*
*		d0 = 0.0065905
*		d1 = 0.5336
*		d2 = -1150.0
*
* The damping part is linear in velocity and scaled by the elastic part:
*
*		fDamping = - kv*vy * fElastic,
*
* where kv = 0.5 by default and vy is the component of the velocity of the
* Body Station in the normal direction relative to the floor plane.
*
* So, all together, the spring force in the normal direction is given by
*
*		fNormal	= fElastic + fDamping
*				= d1*exp(d2*(py-d0)) - kv*vy*d1*exp(d2*(py-d0))
*
* Friction \n
* Friction is implemented using a simple linear, damped spring.
*
*		friction = -kp*(p-p0) - kv*(v)
*
* where kp is the elasticity, kv is the viscosity, p is the position of the
* spring Station on the Body projected onto the plane of the floor, p0 is the
* spring zero in the plane of the floor, and v is the velocity of the spring
* Station projected onto the floor.  p, p0, and v are all expressed in the
* frame of the floor.
*
* By default, the frictional parameters kp and kv are chosen to result in critical
* damping for a specified mass:
*
*		kv = 2.0 * sqrt(kp*mass)
*
* Valid values of kp can range widely (e.g., kp = 1,000 to kp = 1,000,000)
* depending on the material properties of the Body and Floor.  In general, the
* higher kp and kv, the smaller the integration step size will need to be
* in order to produce an accurate integration.
*
* When the frictional force is greater than its allowed limit, the frictional
* force is capped at this limit:
*
*		fLimit = mu*fNormal
*		if(friction.norm() > fLimit) {
*			friction = fLimit * friction.normalize()
*		}	
* 
* where mu is the instantaneous coefficient of friction.
* 
* In addition, if the magnitude of the elastic part of the frictional
* force by itself exceeds fLimit, a new spring zero (p0) is found such
* that the magnitude of the elastic part is equal to fLimit:
*
*		if( (kp*(p-p0)).norm() > fLimit)  p0New = p + (p-p0)/kp
*
* In Simbody, p0 is handled as an Auto Update Discrete State.  Any change
* to p0 is made to the Update Cache (not to the State directly), and the
* integrator swaps this cache value with the actual p0 State after a
* successful integration step is obtained.
*
* Coefficients of kinetic (sliding) and static (fixed) friction can be
* specified for the spring subject to the following constraints:
*
*		0.0 <= muk <= mus <= 1.0
*
* A continuous state variable (Z = Sliding) is used to characterize the
* sliding state of the spring.
*
* 		Sliding = 0		means fully fixed in place (lower bound)
* 		Sliding = 1		means fully sliding (upper bound)
* 
* The instantaneous coefficient of friction (mu) is calculated based on the
* value of Sliding:
*
*		mu = mus - Sliding*(mus - muk)
* 
* The time derivative of Sliding is used to drive Sliding toward the
* extreme of 0.0 or 1.0, depending on the following criteria...
* 
* When the frictional force exceeds its limit, Sliding is driven toward 1:
* 
* 		if (friction.norm() >= fLimit)  SlidingDot = (1.0 - Sliding)/T
* 
* When the velocity of p falls below some specified velocity (e.g., 0.01 m/s)
* AND the frictional force is less than its limit, Sliding is driven toward 0:
*
* 		else if (v.norm < 0.01)  SlidingDot = -Sliding/T
* 
* Otherwise, no change to the Sliding state is made at all:
*
* 		else  SlidingDot = 0.0
*
* In the above equations, T is the characteristic rate at which the Sliding
* rises or decays (e.g., T = 0.01 sec).
* 
* The motivation for using a continuous state variable for Sliding is that,
* although the transition between fixed and sliding may happen quickly,
* it does not happen instantaneously.  And, modeling Sliding based on a
* differential equation ensures that mu is continuous.
*
* A future enhancement might include the capacity to model the ground
* or floor with a polygonal mesh, as opposed to an infinite plane.
*
* ------------------------------------------------------------------
* PARAMETERS
* ------------------------------------------------------------------
* Parameters specifying the characteristics of the exponential
* spring are handled using %ExponentialSpringParameters.
*
* @see ExponentialSpringParameters.
*
* Whenever the paramters are set, the system is invalated at
* the %Stage::Topology.
*
* ------------------------------------------------------------------
* STATES
* ------------------------------------------------------------------
* DISCRETE STATES (Changeable during simulation)\n
* mus:	static coefficient of friction.
*		(0.0 <= mus <= 1.0; default = 0.7)
* muk:	kinetic coefficient of friction.
*		(0.0 <= muk <= mus; default = 0.5)
* Station:	%Station on the Body expressed in the Body frame at which the spring
*			force is applied.  (NOT A STATE YET. @todo make station a state)
*
* AUTO UPDATE DISCRETE STATE\n
* (Computed and put in cache during %System::realize(). Swapped to the actual state
* after a successful integration step.)
* p0:	zero of the frictional spring.
*		Calling resetSprZero() will set p0 equal to Station and then
*		project p0 on to the floor. To project on to the fllor, the
*		y component is simpmly set to zero (p0[1] = 0.0).
*
* CONTINUOUS STATE (Governed by differential equation. Evolved by integrator.)\n
* Sliding:	characterizes whether the spring zero is sliding (1.0) or fixed in place (0.0).
*			This state is used for transitioning back and forth between mus and muk.
*/
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
	// Reset the spring zero (auto update discrete state) so that it
	// coincides in the floor plane with the specified station on the body.
	// This is often done at the beginning of a simulation to ensure that
	// a simulation does not begin with large friction forces.
	void resetSprZero(State& state) const;
	// Data Cache
	const ExponentialSpringData& calcData(const State& state) const;
	const ExponentialSpringData& getData(const State& state) const;
	// Shortcuts for accessing the underlying implementation
	ExponentialSpringForceImpl& updImpl();
	const ExponentialSpringForceImpl& getImpl() const;
};


} // namespace SimTK
#endif SimTK_SIMBODY_EXPONENTIAL_SPRING_FORCE_H_

