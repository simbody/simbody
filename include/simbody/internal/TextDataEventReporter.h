#ifndef SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_
#define SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_

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


/**
 * This is an event reporter which prints out numeric data at regular intervals in tabular form.
 * You provide it with a UserFunction, which calculates the values to be reported.  At every
 * reporting interval, it invokes the UserFunction, then prints out the current time along
 * with the value or values returned by the function.
 * 
 * After creating a TextDataEventReporter, invoke addEventReporter() on the System's default subsystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

class SimTK_SIMBODY_EXPORT TextDataEventReporter : public ScheduledEventReporter {
public:
    TextDataEventReporter(const System& system, UserFunction<Real>* function, Real reportInterval);
    TextDataEventReporter(const System& system, UserFunction<Vector>* function, Real reportInterval);
    ~TextDataEventReporter();
    Real getReportInterval() const;
    void setReportInterval(Real interval);
    Real getNextEventTime(const State&) const;
    void handleEvent(const State& state);
    class TextDataEventReporterRep;
protected:
    TextDataEventReporterRep* rep;
    const TextDataEventReporterRep& getRep() const {assert(rep); return *rep;}
    TextDataEventReporterRep&       updRep() const {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_
