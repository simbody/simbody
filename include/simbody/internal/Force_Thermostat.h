#ifndef SimTK_SIMBODY_FORCE_THERMOSTAT_H_
#define SimTK_SIMBODY_FORCE_THERMOSTAT_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Chris Bruns                                                  *
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

/** @file
 * This contains the user-visible API ("handle" class) for the SimTK::Force 
 * subclass Force::Thermostat and is logically part of Force.h. The file 
 * assumes that Force.h will have included all necessary declarations.
 */

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

namespace SimTK {

/**
 * This is a feedback-controlled force that uses Nose'-Hoover chains to 
 * maintain a particular temperature Tb, as though the system were immersed in 
 * an infinite heat bath at that temperature. There are two parameters, the 
 * temperature Tb and a "relaxation time" t which controls how tightly the 
 * temperature is maintained. This thermostat is particularly useful in 
 * molecular simulations but can be applied to any mechanical system also.
 *
 * Temperature is defined here as T = (2*KE) / (N*kB) where KE is the system
 * kinetic energy, N is the number of coupled degrees of freedom (mobilities 
 * minus active, nonredundant constraints, minus up to 6 rigid body dofs for 
 * the system as a whole), and kB is Boltzmann's constant in appropriate units. 
 *
 * We use a Nose'-Hoover chain to achieve excellent statistical mechanics 
 * properties with a continuous force. At equilibrium the temperature will have
 * a Boltzmann distribution; the relaxation time controls how long it takes the
 * system to reach equilibrium with the bath. Smaller values of relaxation time
 * produce faster response but can make the system stiff and will normally 
 * require smaller step sizes; larger values will take longer to equilibrate 
 * but will run faster.
 *
 * This Force does not produce any potential energy. However, there is a "bath 
 * energy" available through a separate call which can be used in combination 
 * with the system energy to construct a conserved quantity; this is described 
 * further below.
 *
 * \par Theory:
 *
 * The current system temperature is defined 
 * <pre>
 *      T = (2*KE) / (N*kB)
 * </pre>
 * where KE is the kinetic energy of the moving bodies whose N degrees of
 * freedom are being controlled (not necessarily all the bodies in the system),
 * and kB is Boltzmann's constant. Our goal here is to control T so that it
 * follows a Boltzmann distribution around the specified bath temperature Tb.
 *
 * For an m-chain Nose'-Hoover chain, we will define m auxiliary "thermostat" 
 * state variables c[i], 0<=i<m, with units of 1/time. The 0'th thermostat 
 * variable c[0] is used to generate a generalized force f applied to the 
 * system mobilities u:
 * <pre>
 *		f = -c[0] * M * u
 * </pre>
 * where M is the system mass matrix and u is the vector of generalized speeds.
 * (Note that in Simbody the M*u product is formed in O(n) time; M itself
 * is never formed.) The c variables should be initialized to zero at the 
 * start of a simulation. Ideally, you should initialize the u's so that they 
 * are already at the right temperature, but if not you should still make them
 * non-zero -- you can see above that if you have no velocities you will get 
 * no Nose'-Hoover forces.
 *
 * If m==1, we have the standard Nose'-Hoover method except with a relaxation 
 * time specified instead of the thermal mass parameter, as in reference [2]:
 * <pre>
 *      cdot[0]   =   (T/Tb - 1)/t^2
 * </pre>
 * Otherwise, for m > 1 we have:
 * <pre>
 *      cdot[0]   =   (T/Tb - 1)/t^2 - c[0]*c[1]
 *      cdot[1]   = N*c[0]^2 - 1/t^2 - c[1]*c[2]
 *      cdot[i]   = c[i-1]^2 - 1/t^2 - c[i]*c[i+1]   (2<=i<m-1)
 *      cdot[m-1] = c[m-2]^2 - 1/t^2 
 * </pre>
 * For comparison with the literature where thermal mass parameters Qi are 
 * used, we use Q0 = N kB Tb t^2 and Qi = kB Tb t^2, i > 0. That is, the first 
 * thermostat that controls the N thermal degrees of freedom is N times 
 * "heavier" than the subsequent ones, each of which controls only the one 
 * dof of its next-lower thermostat. See refs [1] and [2].
 * 
 * In addition there is a set of state variables si given by sdot[i]=c[i]. 
 * Together these permit us to define a "bath energy" which can be combined 
 * with system energy to produce a conserved quantity. Bath energy is KEb+PEb
 * where
 * <pre>
 *		KEb = 1/2 kB Tb t^2 (N c[0]^2 + sum(c[i]^2))
 *		PEb = kB Tb (N s[0] + sum(s[i]))
 * </pre>
 * where kB is Boltzmann's constant, Tb the bath temperature, N the number of 
 * thermal degrees of freedom in the temperature definition, and the sums run 
 * from 1 to m-1. Note that you must request the bath energy separately; we do
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
	/// a given default temperature and relaxation time. The number of 
    /// Nose'-Hoover chains is given a default value.
    Thermostat(GeneralForceSubsystem&        forces, 
               const SimbodyMatterSubsystem& matter, 
			   Real                          boltzmannsConstant, 
               Real                          bathTemperature, 
               Real                          relaxationTime,
               int                           numExcludedDofs = 6);

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

    /// Set the default number of system rigid body degrees of freedom (0-6)
    /// to be excluded from the calculation of the number of thermal degrees 
    /// of freedom N; if you don't call this it is assumed that 6 dofs should
    /// be excluded.
    Thermostat& setDefaultNumExcludedDofs(int numExcludedDofs);

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
	/// Get the initial value for the number of system rigid body degrees
    /// of freedom (0-6) to be excluded from the calculation of the number of
    /// thermal degrees of freedom N. A new value may be set in a particular 
    /// State.
	int getDefaultNumExcludedDofs() const;

	/// Can't change the value of Boltzmann's constant after construction;
    /// this is the value being used.
	Real getBoltzmannsConstant() const;

	/// Set the actual number of Nose'-Hoover chains to be used. This variable
	/// controls the number of auxiliary state variables allocated by the 
    /// Thermostat so invalidates Model stage.
	const Thermostat& setNumChains(State&, int numChains) const;
	/// Set the bath temperature which serves as the target temperature for
	/// the thermostat. This is given in units defined by the value of 
	/// Boltzmann's constant (which has units of energy/temperature) that was 
	/// set on construction. This sets an Instance-stage state variable so
	/// invalidates Instance and higher stages in the given State.
	const Thermostat& setBathTemperature(State&, Real Tb) const;
	/// Set the relaxation time which determines how long the system will
	/// take to equilibrate to the bath temperature. This sets an 
	/// Instance-stage state variable so invalidates Instance and higher 
	/// stages in the given State.
	const Thermostat& setRelaxationTime(State&, Real t) const;
    /// Set the actual number of system rigid body degrees of freedom (0-6)
    /// to be excluded from the calculation of the number of thermal degrees 
    /// of freedom N. This sets an Instance-stage state variable so invalidates
    /// Instance and higher stages in the given State.
    const Thermostat& setNumExcludedDofs(State&, int numExcludedDofs) const;

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
	/// Get the current value for the number of system rigid body degrees
    /// of freedom (0-6) to be excluded from the calculation of the number of
    /// thermal degrees of freedom N. This is a state variable so can 
	/// be obtained any time after realization of the Model stage.
	int getNumExcludedDofs(const State&) const;

	/// Return the number of thermal degrees of freedom being used in the 
    /// definition of temperature for this thermostat. This is the net of the
    /// total number of mobilities minus nonredundant constraints minus
    /// the number of excluded system rigid body degrees of freedom (0-6).
	int getNumThermalDofs(const State&) const;

	/// Return the temperature of the controlled degrees of freedom via the 
    /// definition T = 2*ke / (N*Kb) where N is the number of thermal degrees 
    /// of freedom. You can call this after Stage::Velocity has been realized.
	Real getCurrentTemperature(const State&) const;

	/// This is a solver that initializes thermostat state variables to zero.
	void initializeChainState(State&) const;
	/// Set the thermostat state variables to particular values. The Vector's
	/// length must be the same as twice the current number of chains called 
    /// for by the State.
	void setChainState(State&, const Vector&) const;

	/// Return the current values of the thermostat chain variables. The 
	/// returned vector will have twice the length that getNumChains() would 
    /// return if called on this same State.
	Vector getChainState(const State&) const;

	/// Calculate the total "bath energy" which, when added to the system
	/// energy, should yield a conserved quantity (assuming all other forces
	/// are conservative).
    Real calcBathEnergy(const State& state) const;

    /// Get the amount of power being applied by the thermostat to the 
    /// system; sign is positive when energy is coming from the bath.
    Real getExternalPower(const State& state) const;

    /// Get the amount of work that has been done by the bath on the
    /// system since an arbitrary start time.
    Real getExternalWork(const State& state) const;

    /// Set the current value of the work done by the bath to an
    /// arbitrary value; normally zero for initialization.
    void setExternalWork(State& state, Real work) const;

	/// Set the controlled system to a set of randomized velocities which
	/// yields the bath temperature. This ignores the current system velocities.
	/// TODO: not implemented yet.
	void initializeSystemToBathTemperature(State&) const;

	/// Set the controlled system to a set of randomized velocities which
	/// yields a particular temperature. This ignores the current system 
    /// velocities. The temperature is interpreted using the value of 
    /// Boltzmann's constant that was provided on construction of this 
    /// Thermostat. TODO: not implemented yet.
	void setSystemToTemperature(State&, Real T) const;

    // Don't show this in Doxygen.
    /// @cond
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Thermostat, ThermostatImpl, Force);
    /// @endcond
};

} // namespace SimTK


#endif // SimTK_SIMBODY_FORCE_THERMOSTAT_H_
