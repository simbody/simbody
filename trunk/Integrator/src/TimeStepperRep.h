#ifndef SimTK_SIMMATH_TIMESTEPPER_REP_H_
#define SimTK_SIMMATH_TIMESTEPPER_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

/** @file
 * This is the declaration of the TimeStepperRep class which
 * represents the implementation of the TimeStepper class.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include "simmath/Integrator.h"
#include "simmath/TimeStepper.h"

#include <exception>

namespace SimTK {

    ////////////////////////////
    // CLASS TIME STEPPER REP //
    ////////////////////////////

// This is a concrete class -- there is only one kind of SimTK TimeStepper.
// The IntegratorRep on the other hand is abstract.
class TimeStepperRep {
public:
    TimeStepperRep(TimeStepper* handle, const System& system);
    // default destructor, no default constructor, no copy or copy assign

    void stepTo(Real time);

    const State& getState() const {return integ->getState();}

    const System getSystem() const {return system;}
    void setIntegrator(Integrator& integrator) {
        integ = &integrator;
    }
    const Integrator& getIntegrator() const {
        assert(integ);
        return *integ;
    }
    Integrator& updIntegrator() {
        assert(integ);
        return *integ;
    }
    bool getReportAllSignificantStates() {
        return reportAllSignificantStates;
    }
    void setReportAllSignificantStates(bool b) {
        reportAllSignificantStates = b;
    }

private:
    TimeStepper* myHandle;
    friend class TimeStepper;

    // This is the System to be advanced through time.
    const System& system;

    // This is the Integrator used to advance through continuous intervals. The Integrator
    // maintains the current State for the System.
    Integrator* integ;
    
    // Whether to report every significant state.
    bool reportAllSignificantStates;

    // suppress
    TimeStepperRep(const TimeStepperRep&);
    TimeStepperRep& operator=(const TimeStepperRep&);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_TIMESTEPPER_REP_H_


