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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath "time stepper" tools.
 */

#include "SimTKcommon.h"

#include "simmath/internal/common.h"

namespace SimTK {
class Integrator;

/**
 * Given a system of equations, an object of this class will advance that
 * system through time. The system is expected to be a mixture of continuous
 * and discrete equations, involving continuous state variables y and 
 * discrete state variables d. The continuous part is an ODE-on-a-manifold
 * system suitable for solution via coordinate projection[1], structured like
 * this:
 *         (1)  y' = f(d;t,y)         differential equations
 *         (2)  c  = c(d;t,y)         algebraic equations (manifold is c=0)
 *         (3)  e  = e(d;t,y)         event triggers (watch for zero crossings)
 * with initial conditions t0,y0,d0 such that c=0. By "ODE on a manifold" we
 * mean that the ODE (1) automatically satisfies the condition that IF c==0,
 * THEN c'=0, where c'=partial(c)/partial(t) + [partial(c)/partial(y)]*y'. This is
 * a less stringent condition than an ODE with invariant, in which c'=0 regardless
 * of whether c=0.
 *
 * [1] Hairer, Lubich, Wanner, "Geometric Numerical Integration: Structure-Preserving
 * Algorithms for Ordinary Differential Equations", 2nd ed., section IV.4, pg 109ff,
 * Springer, 2006.
 *
 * The discrete variables d are updated upon occurence of specific events, which are
 * detected using a set of scalar-valued event trigger functions (3).
 * An event trigger function for a particular event should be designed so that it has a
 * zero crossing when the event occurs. The integrator can thus watch for sign changes
 * in event triggers and terminate the current step when a zero crossing occurs,
 * notifying the system and giving it a chance to handle the event; that is, 
 * update its state variables discontinuously.
 *
 * The zero crossings of continuous event trigger functions will
 * be isolated quickly; discontinuous ones have to be "binary chopped" which
 * is more expensive. There are several special case events:
 *          - events which are simply a known function of t
 *          - "end of step" events
 *          - external events (e.g., clock time, user interrupt)
 *
 * We are given a set of weights W for the y's, and a set of tolerances T
 * for the constraint errors. Given an accuracy specification (like 0.1%),
 * the integrators here are expected to solve for y(t) such that the
 * local error |W*y|_RMS <= accuracy, and |T*c(t,y)|_RMS <= accuracy at
 * all times.
 *
 * TODO: isolation tolerances for witnesses; dealing with simultaneity.
 *
 */
class SimTK_SIMMATH_EXPORT TimeStepper {
public:
    explicit TimeStepper(const System&);
    TimeStepper(const System&, Integrator&);
    ~TimeStepper();

    void setIntegrator(Integrator&);
    const Integrator& getIntegrator() const;
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
    
    // Supply the time stepper with a starting state. This is *copied*
    // into the time stepper's internally maintained "advanced" state; 
    // subsequent changes to the State object passed in here will not
    // affect the simulation.
    void initialize(const State&);

    // Return a state which satisfies the caller's step request. This 
    // may be an interpolated value earlier than getAdvancedState().
    const State& getState() const;
    Real         getTime() const {return getState().getTime();}

    void stepTo(Real time);
    // opaque implementation for binary compatibility
    class TimeStepperRep* rep;
    friend class TimeStepperRep;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_TIMESTEPPER_H_
