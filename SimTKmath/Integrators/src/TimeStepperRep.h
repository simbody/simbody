#ifndef SimTK_SIMMATH_TIMESTEPPER_REP_H_
#define SimTK_SIMMATH_TIMESTEPPER_REP_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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

    Integrator::SuccessfulStepStatus stepTo(Real time);

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
    bool getReportAllSignificantStates() const {
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
    
    // The last time at events and reports were processed.
    Real lastEventTime, lastReportTime;

    // suppress
    TimeStepperRep(const TimeStepperRep&);
    TimeStepperRep& operator=(const TimeStepperRep&);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_TIMESTEPPER_REP_H_


