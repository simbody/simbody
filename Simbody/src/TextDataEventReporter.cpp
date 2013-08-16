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

#include "simbody/internal/TextDataEventReporter.h"

using std::cout;
using std::endl;
using namespace SimTK;

/**
 * This class is the internal implementation for TextDataEventReporter.  It has two subclasses,
 * one for UserFunctions that return a Real and one for UserFunctions that return a Vector.
 */

class TextDataEventReporter::TextDataEventReporterRep {
public:
    TextDataEventReporterRep(const System& system) : system(system) {
    }
    virtual ~TextDataEventReporterRep() { }
    virtual void printValues(const State& state) const = 0;
    void handleEvent(const State& state) const {
        cout << state.getTime();
        printValues(state);
        cout << endl;
    }
    TextDataEventReporter* handle;
    const System& system;
    Real reportInterval;
    class RealFunction;
    class VectorFunction;
};

class TextDataEventReporter::TextDataEventReporterRep::RealFunction : public TextDataEventReporter::TextDataEventReporterRep {
public:
    RealFunction(const System& system, UserFunction<Real>* function) : TextDataEventReporterRep(system), function(function) {
    }
    ~RealFunction() {
        delete function;
    }
    void printValues(const State& state) const {
        Real value = function->evaluate(system, state);
        cout << "\t" << value;
    }
    UserFunction<Real>* function;
};

class TextDataEventReporter::TextDataEventReporterRep::VectorFunction : public TextDataEventReporter::TextDataEventReporterRep {
public:
    VectorFunction(const System& system, UserFunction<Vector>* function) : TextDataEventReporterRep(system), function(function) {
    }
    ~VectorFunction() {
        delete function;
    }
    void printValues(const State& state) const {
        Vector values = function->evaluate(system, state);
        for (int i = 0; i < values.size(); ++i)
            cout << "\t" << values[i];
    }
    UserFunction<Vector>* function;
};

TextDataEventReporter::TextDataEventReporter(const System& system, UserFunction<Real>* function, Real reportInterval) : PeriodicEventReporter(reportInterval) {
    rep = new TextDataEventReporterRep::RealFunction(system, function);
    updRep().handle = this;
}

TextDataEventReporter::TextDataEventReporter(const System& system, UserFunction<Vector>* function, Real reportInterval) : PeriodicEventReporter(reportInterval) {
    rep = new TextDataEventReporterRep::VectorFunction(system, function);
    updRep().handle = this;
}

TextDataEventReporter::~TextDataEventReporter() {
    if (rep->handle == this)
        delete rep;
}

void TextDataEventReporter::handleEvent(const State& state) const {
    updRep().handleEvent(state);
}
