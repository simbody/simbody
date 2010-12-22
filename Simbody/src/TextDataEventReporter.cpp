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
