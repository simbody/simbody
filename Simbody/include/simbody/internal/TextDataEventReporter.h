#ifndef SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_
#define SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {


/** This is an EventReporter which prints out numeric data at regular intervals
in tabular form. You provide it with a UserFunction, which calculates the
values to be reported.  At every reporting interval, it invokes the
UserFunction, then prints out the current time along with the value or
values returned by the function.

After creating a TextDataEventReporter, add it to the System by calling the
addEventReporter() method. **/
class SimTK_SIMBODY_EXPORT TextDataEventReporter
:   public PeriodicEventReporter {
public:

    /** This template class defines a standard interface for objects that
    calculate a function based on a System and State for use in a
    TextDataEventReporter. **/
    template <class T> class UserFunction {
    public:
        virtual ~UserFunction() {}
        virtual T evaluate(const System& system, const State& state) = 0;
    };

    /** Create a TextDataEventReporter which calculates a single number at each
    reporting interval, and displays it along with the time. Takes ownership
    of the UserFunction object. **/
    TextDataEventReporter(const System&         system,
                          UserFunction<Real>*   function,
                          Real                  reportInterval);

    /** Create a TextDataEventReporter which calculates a vector of numbers at
    each reporting interval, and displays them along with the time. Takes
    ownership of the UserFunction object.  **/
    TextDataEventReporter(const System&         system,
                          UserFunction<Vector>* function,
                          Real                  reportInterval);
    /** The destructor will take care of deleting the UserFunction object;
    don't delete it yourself. **/
    ~TextDataEventReporter();

    /** This is the implementation of the EventReporter virtual. **/
    void handleEvent(const State& state) const override;

    class TextDataEventReporterRep;
protected:
    TextDataEventReporterRep* rep;
    const TextDataEventReporterRep& getRep() const {assert(rep); return *rep;}
    TextDataEventReporterRep&       updRep() const {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_TEXT_DATA_EVENT_REPORTER_H_
