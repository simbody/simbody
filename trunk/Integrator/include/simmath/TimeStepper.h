#ifndef SimTK_SIMMATH_TIMESTEPPER_H_
#define SimTK_SIMMATH_TIMESTEPPER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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

#include "simmath/internal/common.h"
#include "Integrator.h"

namespace SimTK {
class Integrator;

/**
 * This class uses an Integrator to advance a System through time.  For example:
 * 
 * <pre>
 * TimeStepper stepper(system, integrator);
 * stepper.initialize(initialState);
 * stepper.stepTo(finalTime);
 * </pre>
 * 
 * stepTo() invokes the Integrator to advance time.  It detects events which may occur, calls event handlers as
 * appropriate, then invokes the Integrator again to continue advancing time until the target time is reached.
 */
class SimTK_SIMMATH_EXPORT TimeStepper {
public:
    /**
     * Create a TimeStepper to advance a System.
     * 
     * This constructor leaves the Integrator unspecified.  You therefore must call setIntegrator()
     * before calling initialize().
     */
    explicit TimeStepper(const System& system);
    /**
     * Create a TimeStepper to advance a System using an Integrator.
     */
    TimeStepper(const System& system, Integrator& integrator);
    ~TimeStepper();
    /**
     * Set the Integrator this TimeStepper will use to advance the System.
     */
    void setIntegrator(Integrator& integrator);
    /**
     * Get the Integrator being used to advance the System.
     */
    const Integrator& getIntegrator() const;
    /**
     * Get a non-const reference to the Integrator being used to advance the System.
     */
    Integrator& updIntegrator();
    /**
     * Get whether the TimeStepper should report every significant state returned by the Integrator.
     * If this is true, stepTo() will return whenever the Integrator reports a significant state,
     * such as when an event occurs or the start of a new continuous interval.  If this is false,
     * stepTo() will only return when the specified time has been reached or when the simulation is
     * terminated.
     */
    bool getReportAllSignificantStates();
    /**
     * Set whether the TimeStepper should report every significant state returned by the Integrator.
     * If this is true, stepTo() will return whenever the Integrator reports a significant state,
     * such as when an event occurs or the start of a new continuous interval.  If this is false,
     * stepTo() will only return when the specified time has been reached or when the simulation is
     * terminated.
     */
    void setReportAllSignificantStates(bool b);
    /**
     * Supply the time stepper with a starting state.  This must be called after the Integrator has
     * been set, and before the first call to stepTo().
     * 
     * The specified state is copied into the Integrator's internally maintained "advanced" state;
     * subsequent changes to the State object passed in here will not affect the simulation.
     */
    void initialize(const State&);
    /**
     * Get the current State of the System being integrated.  Usually this will correspond to the
     * time specified in the most recent call to stepTo().  It may be an interpolated state which
     * is earlier than the Integrator's "advanced" state.
     */
    const State& getState() const;
    /**
     * Get the current time of the System being integrated.  This is identical to calling getState().getTime().
     */
    Real getTime() const {return getState().getTime();}
    /**
     * Use the Integrator to advance the System up to the specified time.  This method will repeatedly
     * invoke the Integrator as necessary, handling any events which occur.  The TimeStepper must be
     * initialized by calling initialize() before this method may be called.
     * 
     * When this method returns, the System will usually have been advanced all the way to the specified
     * final time.  There are situations where it may return sooner, however: if setReportAllSignificantStates()
     * was set to true, and a significant state occurred; if {@link Integrator#setFinalTime setFinalTime()} was invoked on the Integrator,
     * and the final time was reached; or if an event handler requested that the simulation terminate
     * immediately.
     */
    Integrator::SuccessfulStepStatus stepTo(Real time);
private:
    class TimeStepperRep* rep;
    friend class TimeStepperRep;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_TIMESTEPPER_H_
