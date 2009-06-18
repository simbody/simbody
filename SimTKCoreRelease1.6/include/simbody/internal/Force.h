#ifndef SimTK_SIMBODY_FORCE_H_
#define SimTK_SIMBODY_FORCE_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/GeneralForceSubsystem.h"

namespace SimTK {

class SimbodyMatterSubsystem;
class GeneralForceSubsystem;
class MobilizedBody;
class Force;
class ForceImpl;

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the Force class but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_FORCE
    extern template class PIMPLHandle<Force, ForceImpl, true>;
#endif

/**
 * A Force object applies forces to the bodies in a system.  There are subclasses for various standard types
 * of forces, or you can create your own forces by using Custom.
 */
    
class SimTK_SIMBODY_EXPORT Force : public PIMPLHandle<Force, ForceImpl, true> {
public:
    Force() { }
    explicit Force(ForceImpl* r) : HandleBase(r) { }

	/// Get the GeneralForceSubsystem of which this Force is an element.
	const GeneralForceSubsystem& getForceSubsystem() const;

    /// Get the index of this Force in its GeneralForceSubsystem.
    ForceIndex getForceIndex() const;
    
    class TwoPointLinearSpring;
    class TwoPointLinearDamper;
    class TwoPointConstantForce;
    class MobilityLinearSpring;
    class MobilityLinearDamper;
    class MobilityConstantForce;
    class ConstantForce;
    class ConstantTorque;
    class GlobalDamper;
    class Thermostat;
    class UniformGravity;
    class Custom;
    
    class TwoPointLinearSpringImpl;
    class TwoPointLinearDamperImpl;
    class TwoPointConstantForceImpl;
    class MobilityLinearSpringImpl;
    class MobilityLinearDamperImpl;
    class MobilityConstantForceImpl;
    class ConstantForceImpl;
    class ConstantTorqueImpl;
    class GlobalDamperImpl;
    class ThermostatImpl;
    class UniformGravityImpl;
    class CustomImpl;
};

/**
 * A linear spring between two points, specified as a station on
 * each of two bodies. The stiffness k and 
 * unstretched length x0 are provided. Then if d is the unit vector
 * from point1 to point2, and x the current separation, we have
 * f = k(x-x0) and we apply a force f*d to point1 and -f*d to point2.
 * This contributes to potential energy: pe = 1/2 k (x-x0)^2.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearSpring : public Force {
public:
    /**
     * Create a TwoPointLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param k          the spring constant
     * @param x0         the distance at which the force is 0
     */
    TwoPointLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real k, Real x0);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointLinearSpring, TwoPointLinearSpringImpl, Force);
};

/**
 * A force which resists changes in the distance between
 * two points, acting along the line between those points. The points
 * are specified as a station on each of two bodies. A damping constant
 * c >= 0 is given. If the relative (scalar) velocity between the points
 * is v, then we apply a force of magnitude f=c*|v| to each point in
 * a direction which opposes their separation. This is not a potential
 * force and thus does not contribute to the potential energy calculation.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearDamper: public Force {
public:
    /**
     * Create a TwoPointLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param damping    the damping constant
     */
    TwoPointLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real damping);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointLinearDamper, TwoPointLinearDamperImpl, Force);
};

/**
 * A constant force f (a signed scalar) which acts along the line between
 * two points, specified as a station on each of two bodies. A positive
 * force acts to separate the points; negative pulls them together. The
 * force magnitude is independent of the separation between the points.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointConstantForce: public Force {
public:
    /**
     * Create a TwoPointConstantForce.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param force      the magnitude of the force to apply
     */
    TwoPointConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real force);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointConstantForce, TwoPointConstantForceImpl, Force);
};

/**
 * A linear spring along or around a mobility coordinate. The
 * stiffness k is provided, along with an arbitrary "zero" 
 * coordinate value q0 at which the spring generates no force.
 * The generated force is k*(q-q0), and potential energy is 
 * pe = 1/2 k (q-q0)^2.
 * This is not meaningful unless the mobility coordinate is such that qdot=u 
 * for that coordinate.  In particular, do not use this on a coordinate
 * which is part of a quaternion.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearSpring : public Force {
public:
    /**
     * Create a MobilityLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param k          the spring constant
     * @param q0         the value of the coordinate at which the force is 0
     */
    MobilityLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real k, Real q0);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearSpring, MobilityLinearSpringImpl, Force);
};

/**
 * A linear damper on a mobility coordinate. The
 * damping constant c is provided, with the generated force
 * being -c*u where u is the mobility's generalize speed.
 * This is meaningful on any mobility, since all our
 * generalized speeds have physical meaning. This is not
 * a potential force and hence does not contribute to
 * potential energy.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearDamper : public Force {
public:
    /**
     * Create a MobilityLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param damping    the damping constant
     */
    MobilityLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real damping);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearDamper, MobilityLinearDamperImpl, Force);
};

/**
 * A constant (scalar) "force" f applied to a mobility. The mobility here selects a
 * generalized speed (u), not a generalized coordinate (q), and the
 * meaning depends on the definition of the generalized speed. If that
 * speed is a translation then this is a force; if a rotation then
 * this is a torque.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityConstantForce : public Force {
public:
    /**
     * Create a MobilityConstantForce.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param force      the force to apply
     */
    MobilityConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real force);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityConstantForce, MobilityConstantForceImpl, Force);
};

/**
 * A constant force applied to a body station. The force is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantForce: public Force {
public:
    ConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& station, const Vec3& force);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ConstantForce, ConstantForceImpl, Force);
};

/**
 * A constant torque to a body. The torque is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantTorque: public Force {
public:
    ConstantTorque(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& torque);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ConstantTorque, ConstantTorqueImpl, Force);
};

/**
 * A general energy "drain" on the system. This is 
 * done by effectively adding a damper to every generalized
 * speed (mobility) in the system. Each generalized speed 
 * u_i feels a force -dampingFactor*u_i.
 * This usually is not physically meaningful, but it can be useful in some
 * circumstances just to drain energy out of the model when
 * the specific energy-draining mechanism is not important.
 * You can have more than one of these in which case the
 * dampingFactors are added. No individual dampingFactor is
 * allowed to be negative. This is not a potential force and
 * hence does not contribute to potential energy.
 */

class SimTK_SIMBODY_EXPORT Force::GlobalDamper : public Force {
public:
    GlobalDamper(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, Real damping);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(GlobalDamper, GlobalDamperImpl, Force);
};

/**
 * This is a feedback-controlled force that maintains a particular temperature Tb, as though
 * the system were immersed in an infinite heat bath at that temperature. There are two
 * parameters, the temperature Tb and a "relaxation time" t which controls how tightly
 * the temperature is maintained. This thermostat is particularly useful in molecular 
 * simulations but can be applied to any mechanical system also.
 *
 * Temperature is defined here as T = (2*KE) / (N*Kb) where KE is the system
 * kinetic energy, N is the number of coupled degrees of freedom (mobilities minus 
 * active, nonredundant constraints), and Kb is Boltzmann's constant in appropriate
 * units. 
 *
 * We use a Nose'-Hoover chain to achieve excellent statistical mechanics properties
 * with a continuous force. At equilibrium the temperature will have a Boltzmann
 * distribution; the relaxation time controls how long it takes the system to reach
 * equilibrium with the bath. Smaller values of relaxation time produce faster
 * response but can make the system stiff and will normally require smaller step
 * sizes; larger values will take longer to equilibrate but will run faster.
 *
 * This Force does not produce any potential energy. However, there is a "bath energy"
 * available through a separate call which can be used in combination with the
 * system energy to construct a conserved quantity; this is described further below.
 *
 * \par Theory:
 *
 * For an m-chain Nose'-Hoover chain, we will define m auxiliary 
 * "thermostat" state variables ci, 0<=i<=m-1, with units of 1/time.
 * The 0'th thermostat variable c0 is used to generate a 
 * force f applied to the system mobilities u:
 * <pre>
 *		f = -c0 * M * u
 * </pre>
 * where M is the system mass matrix and u is the vector of generalized
 * speeds. The c variables should be initialized to zero at the start of a 
 * simulation. Ideally, you should initialize the u's so that they
 * are already at the right temperature, but if not you should still
 * make them non-zero -- you can see above that if you have no
 * velocities you will get no Nose'-Hoover forces.
 *
 * If m==1, we have the standard Nose'-Hoover method except with
 * a relaxation time specified instead of the thermal mass parameter, as
 * in reference [2]:
 * <pre>
 *      cdot[0]   =   (T/Tb - 1)/t^2
 * </pre>
 * Otherwise, for m > 1 we have:
 * <pre>
 *      cdot[0]   =   (T/Tb - 1)/t^2 - c0*c1
 *      cdot[1]   =   N*c0^2 - 1/t^2 - c1*c2
 *      cdot[i]   = c[i-1]^2 - 1/t^2 - c[i]*c[i+1]   (2<=i<m-1)
 *      cdot[m-1] = c[m-2]^2 - 1/t^2 
 * </pre>
 * For comparision with the literature where thermal mass
 * parameters Qi are used, we use Q0 = N Kb Tb t^2 and
 * Qi = Kb Tb t^2, i > 0. That is, the first thermostat
 * that controls the N degrees of freedom is N times "heavier"
 * than the subsequent ones, each of which controls only
 * the one dof of its next-lower thermostat. See refs [1]
 * and [2].
 * 
 * In addition there is a set of state variables si given by
 * sdot[i] = c[i]. Together these permit us to define a "bath energy"
 * which can be combined with system energy to produce a conserved
 * quantity. Bath energy is KEb + PEb where
 * <pre>
 *		KEb = 1/2 Kb Tb t^2 (N c0^2 + sum(ci^2))
 *		PEb = Kb Tb (N s0 + sum(si))
 * </pre>
 * where Kb is Boltzmann's constant, Tb the bath temperature, N the
 * number of degrees of freedom in the temperature definition, and
 * the sums run from 1 to m-1.
 * Note that you must request the bath energy separately; we do
 * not return any potential energy for this force otherwise.
 *
 * \par References:
 *
 * [1] Martyna, GJ; Klien, ML; Tuckerman, M. Nose'-Hoover chains:
 * The canonical ensemble via continuous dynamics. J. Chem. Phys.
 * 97(4):2635-2643 (1992).
 *
 * [2] Vaidehi, N; Jain, A; Goddard, WA. Constant Temperature
 * Constrained Molecular Dynamics: The Newton-Euler Inverse Mass
 * Operator Method. J. Phys. Chem. 100:10508-10517 (1996).
 */

class SimTK_SIMBODY_EXPORT Force::Thermostat : public Force {
public:
	/// Define a global thermostat (one that affects all degrees of freedom) at
	/// a given default temperature and relaxation time. The number of Nose'-Hoover
	/// chains is given a default value.
    Thermostat(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, 
			   Real boltzmannsConstant, Real bathTemperature, Real relaxationTime);

	/// TODO: not implemented yet. Remove a body from consideration in
	/// the thermostat. Typically this would be the system base body so
	/// that overall rigid body translation and orientation is not counted
	/// as part of the temperature.
	Thermostat& excludeMobilizedBody(MobilizedBodyIndex);

	/// Set the default (state independent) number of Nose'-Hoover chains.
	/// This is a Topology-stage change.
	Thermostat& setDefaultNumChains(int numChains);

	/// Set the default (state independent) bath temperature. This will be 
	/// interpreted using the value of Boltzmann's constant Kb provided on
	/// construction. The units will be Kb/energy, typically Kelvins.
	Thermostat& setDefaultBathTemperature(Real bathTemperature);

	/// Set the default (state independent) relaxation time.
	Thermostat& setDefaultRelaxationTime(Real relaxationTime);

	/// Get the initial value for the number of chains that will be used for
	/// the "number of chains" State variable. A new value may be set in 
	/// a particular State.
	int getDefaultNumChains() const;
	/// Get the initial value for the bath temperature that will be use for
	/// the "bath temperature" State variable. A new value may be set in 
	/// a particular State.
	Real getDefaultBathTemperature() const;
	/// Get the initial value for the bath temperature that will be use for
	/// the "bath temperature" State variable.
	Real getDefaultRelaxationTime() const;
	/// Can't change the value of Boltzmann's constant after construction.
	Real getBoltzmannsConstant() const;

	/// Set the actual number of Nose'-Hoover chains to be used. This variable
	/// controls the number of auxiliary state variables allocated by the Thermostat
	/// so invalidates Model stage (TODO: should be Instance).
	void setNumChains(State&, int numChains) const;
	/// Set the bath temperature which serves as the target temperature for
	/// the thermostat. This is given in units defined by the value of 
	/// Boltzmann's constant (which has units of energy/temperature) that was 
	/// set on construction. This sets an Instance-stage state variable so
	/// invalidates Instance and higher stages in the given State.
	void setBathTemperature(State&, Real Tb) const;
	/// Set the relaxation time which determines how long the system will
	/// take to equilibrate to the bath temperature. This sets an 
	/// Instance-stage state variable so invalidates Instance and higher 
	/// stages in the given State.
	void setRelaxationTime(State&, Real t) const;

	/// Obtain the current number of Nose'-Hoover chains in use. This is a
	/// state variable so can be obtained any time after realization of
	/// the Model stage.
	int getNumChains(const State&) const;
	/// Obtain the current bath temperature, in units which are determined
	/// by the value of Boltzmann's constant that was supplied on construction
	/// of this Thermostat force element. This is a state variable so can 
	/// be obtained any time after realization of the Model stage.
	Real getBathTemperature(const State&) const;
	/// Obtain the current relaxation time. This is a state variable so can 
	/// be obtained any time after realization of the Model stage.
	Real getRelaxationTime(const State&) const;

	/// Return the number of degrees of freedom being used in the definition
	/// of temperature for this thermostat. This is the net of the total 
	/// number of mobilities selected minus nonredundant constraints.
	int getNumDegreesOfFreedom(const State&) const;

	/// Return the temperature of the controlled degrees of freedom
	/// via the definition T = 2*ke / (N*Kb) where N is the number
	/// of degrees of freedom. You can call this after Stage::Velocity
	/// has been realized.
	Real getCurrentTemperature(const State&) const;

	/// This is a solver that initializes the thermostat state variables to zero.
	void initializeChainState(State&) const;
	/// Set the thermostat state variables to particular values. The Vector's
	/// length must be the same as twice the current number of chains called for by
	/// the State.
	void setChainState(State&, const Vector&) const;

	/// Return the current values of the thermostat chain variables. The 
	/// returned vector will have twice the length that getNumChains(s) would return
	/// if called on this same State.
	Vector getChainState(const State&) const;

	/// Calculate the total "bath energy" which, when added to the system
	/// energy, should yield a conserved quantity (assuming all other forces
	/// are conservative).
    Real calcBathEnergy(const State& state) const;

	/// Set the controlled system to a set of randomized velocities which
	/// yields the bath temperature. This ignores the current system velocities.
	/// TODO: not implemented yet.
	void initializeSystemToBathTemperature(State&) const;

	/// Set the controlled system to a set of randomized velocities which
	/// yields a particular temperature. This ignores the current system velocities.
	/// The temperature is interpreted using the value of Boltzmann's constant
	/// that was provided on construction of this Thermostat.
	/// TODO: not implemented yet.
	void setSystemToTemperature(State&, Real T) const;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Thermostat, ThermostatImpl, Force);
};

/**
 * A uniform gravitational force applied to every body in the system.  The force
 * is specified by a vector in the Ground frame.  You can optionally specify
 * a height at which the gravitational potential energy is zero.
 */

class SimTK_SIMBODY_EXPORT Force::UniformGravity : public Force {
public:
    UniformGravity(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight=0);
    Vec3 getGravity() const;
    void setGravity(const Vec3& g);
    Real getZeroHeight() const;
    void setZeroHeight(Real height);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(UniformGravity, UniformGravityImpl, Force);
};

/**
 * This class can be used to define new forces.  To use it, create a class that extends Force::Custom::Implementation.
 * You can then create an instance of it and pass it to the Force::Custom constructor:
 * 
 * <pre>
 * Force::Custom myForce(forces, new MyForceImplementation());
 * </pre>
 * 
 * Alternatively, you can create a subclass of Force::Custom which creates the Implementation itself:
 * 
 * <pre>
 * class MyForce : public Force::Custom {
 * public:
 *   MyForce(GeneralForceSubsystem& forces) : Force::Custom(forces, new MyForceImplementation()) {
 *   }
 * }
 * </pre>
 * 
 * This allows a user to simply write
 * 
 * <pre>
 * MyForce(forces);
 * </pre>
 * 
 * and not worry about implementation classes or creating objects on the heap.  If you do this, your Force::Custom
 * subclass must not have any data members or virtual methods.  If it does, it will not work correctly.  Instead,
 * store all data in the Implementation subclass.
 */

class SimTK_SIMBODY_EXPORT Force::Custom : public Force {
public:
    class Implementation;
    /**
     * Create a Custom force.
     * 
     * @param forces         the subsystem to which this force should be added
     * @param implementation the object which implements the custom force.  The Force::Custom takes over
     *                       ownership of the implementation object, and deletes it when the Force itself
     *                       is deleted.
     */
    Custom(GeneralForceSubsystem& forces, Implementation* implementation);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, Force);
protected:
    const Implementation& getImplementation() const;
    Implementation& updImplementation();
};

class SimTK_SIMBODY_EXPORT Force::Custom::Implementation {
public:
    virtual ~Implementation() { }
    /**
     * Calculate the force for a given state.
     * 
     * @param state          the State for which to calculate the force
     * @param bodyForces     spatial forces on MobilizedBodies are accumulated in this.  To apply a force to a body,
     *                       add it to the appropriate element of this vector.
     * @param particleForces forces on particles are accumulated in this.  Since particles are not yet implemented,
     *                       this is ignored.
     * @param mobilityForces forces on individual mobilities (elements of the state's u vector) are accumulated in this.
     *                       To apply a force to a mobility, add it to the appropriate element of this vector.
     */
    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const = 0;
    /**
     * Calculate this force's contribution to the potential energy of the System.
     * 
     * @param state          the State for which to calculate the potential energy
     */
    virtual Real calcPotentialEnergy(const State& state) const = 0;
    /**
     * Get whether this force depends only on the position variables (q), not on the velocies (u) or auxiliary variables (z).
     * The default implementation returns false.  If the force depends only on positions, you should override this to return
     * true.  This allows force calculations to be optimized in some cases.
     */
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    /**
     * The following methods may optionally be overridden to do specialized realization for a Force.
     */
    //@{
    virtual void realizeTopology(State& state) const {
    }
    virtual void realizeModel(State& state) const {
    }
    virtual void realizeInstance(const State& state) const {
    }
    virtual void realizeTime(const State& state) const {
    }
    virtual void realizePosition(const State& state) const {
    }
    virtual void realizeVelocity(const State& state) const {
    }
    virtual void realizeDynamics(const State& state) const {
    }
    virtual void realizeAcceleration(const State& state) const {
    }
    virtual void realizeReport(const State& state) const {
    }
    //@}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_H_
