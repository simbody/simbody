#ifndef SimTK_SIMBODY_VELOCITY_RESCALING_THERMOSTAT_H_
#define SimTK_SIMBODY_VELOCITY_RESCALING_THERMOSTAT_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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
#include "simbody/internal/MultibodySystem.h"

namespace SimTK {

/**
 * This is an event handler that acts as a thermostat for controlling the temperature of a simulation.
 * It does this through velocity rescaling.  At regular intervals, it calculates the total kinetic
 * energy of the system, then rescales all the velocities so the kinetic energy will exactly equal
 * kT/2 per degree of freedom.
 */

class SimTK_SIMBODY_EXPORT VelocityRescalingThermostat : public ScheduledEventHandler {
public:
    /**
     * Create a VelocityRescalingThermostat.
     * @param system the MultibodySystem to be simulated
     * @param boltzmannsConstant the value of Boltzmann's constant in whatever units are used by the system
     * @param temperature the temperature to maintain  the system at.  The default value is 293.15 Kelvin (20 C).
     * @param rescalingInterval the time interval at which to rescale velocities.  The default value is 1.0.
     */
    VelocityRescalingThermostat(const MultibodySystem& system, Real boltzmannsConstant, Real temperature=293.15, Real rescalingInterval=1.0);
    /**
     * Get the temperature this thermostat is set to maintain.
     */
    Real getTemperature();
    /**
     * Set the temperature this thermostat is set to maintain.
     */
    void setTemperature(Real temp);
    /**
     * Get the time interval at which the velocities will be rescaled.
     */
    Real getRescalingInterval();
    /**
     * Set the time interval at which the velocities will be rescaled.
     */
    void setRescalingInterval(Real interval);
    void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate);
    Real getNextEventTime(const State& state) const;
    ~VelocityRescalingThermostat();
private:
    class VelocityRescalingThermostatImpl;
    VelocityRescalingThermostatImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_VELOCITY_RESCALING_THERMOSTAT_H_
