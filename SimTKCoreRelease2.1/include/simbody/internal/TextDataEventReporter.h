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


#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

/**
 * This template class defines a standard interface for objects that calculate 
 * a function based on a System and State for use in a TextDataEventReporter.
 */
template <class T>
class UserFunction {
public:
    virtual T evaluate(const System& system, const State& state) = 0;
};

/**
 * This is an EventReporter which prints out numeric data at regular intervals 
 * in tabular form. You provide it with a UserFunction, which calculates the 
 * values to be reported.  At every reporting interval, it invokes the 
 * UserFunction, then prints out the current time along with the value or 
 * values returned by the function.
 * 
 * After creating a TextDataEventReporter, invoke addEventReporter() on the 
 * System's default subsystem.
 */
class SimTK_SIMBODY_EXPORT TextDataEventReporter : public PeriodicEventReporter {
public:
    /**
     * Create a TextDataEventReporter which calculates a single number at each 
     * reporting interval, and displays it along with the time.
     */
    TextDataEventReporter(const System& system, UserFunction<Real>* function, Real reportInterval);
    /**
     * Create a TextDataEventReporter which calculates a vector of numbers at 
     * each reporting interval, and displays them along with the time.
     */
    TextDataEventReporter(const System& system, UserFunction<Vector>* function, Real reportInterval);
    ~TextDataEventReporter();
    void handleEvent(const State& state) const;
    class TextDataEventReporterRep;
protected:
    TextDataEventReporterRep* rep;
    const TextDataEventReporterRep& getRep() const {assert(rep); return *rep;}
    TextDataEventReporterRep&       updRep() const {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_
