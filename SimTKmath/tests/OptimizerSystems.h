#ifndef SimTK_SIMMATH_OPTIMIZER_SYSTEMS_H_
#define SimTK_SIMMATH_OPTIMIZER_SYSTEMS_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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

// These websites list test functions for optimization.
// https://en.wikipedia.org/wiki/Test_functions_for_optimization
// http://www.sfu.ca/~ssurjano/optimization.html

#include "SimTKmath.h"

using SimTK::Vector;
using SimTK::Real;
using SimTK::OptimizerSystem;
using SimTK::Pi;
using SimTK::square;
using SimTK::sqrt;

class TestOptimizerSystem : public OptimizerSystem {
public:
    TestOptimizerSystem(int nParameters)
        : OptimizerSystem(nParameters) {}
    virtual Real optimalValue() const {
        Real f;
        objectiveFunc(optimalParameters(), true, f);
        return f;
    }
    virtual Vector optimalParameters() const = 0;
};

// This function come from Nikolaus Hansen's source code
// (https://github.com/cma-es). I think it is supposed to be cigar-shaped.
class Cigtab : public TestOptimizerSystem {
public:
    Cigtab(int nParameters) : TestOptimizerSystem(nParameters) {}
    int objectiveFunc(const Vector& x, bool new_parameters,
            Real& f) const {
        f = 1e4 * x[0] * x[0] + 1e-4 * x[1] * x[1];
        for (int i = 0; i < getNumParameters(); ++i) {
            f += x[i] * x[i];
        }
        return 0;
    }
    Vector optimalParameters() const override {
        Vector x(getNumParameters());
        x.setToZero();
        return x;
    }
};

// A function with many local minima. http://www.sfu.ca/~ssurjano/ackley.html
class Ackley : public TestOptimizerSystem {
public:
    Ackley(int nParameters) : TestOptimizerSystem(nParameters),
            a(20), b(0.2), c(2 * Pi) {
        // The website above says this function usually has the following
        // bounds:
        Vector limits(nParameters);
        limits.setTo(32.768);
        setParameterLimits(-limits, limits);
    }
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        const Real n = getNumParameters();
        Real sumcos = 0;
        for (int i = 0; i < n; ++i) {
            sumcos += cos(c * x[i]);
        }
        f = -a * exp(-b * x.normRMS()) - exp(sumcos / n) + a + SimTK::E;
        return 0;
    }
    Vector optimalParameters() const {
        Vector x(getNumParameters());
        x.setToZero();
        return x;
    }
private:
    const Real a;
    const Real b;
    const Real c;
};

// A very complex 2D function http://www.sfu.ca/~ssurjano/drop.html
class DropWave : public TestOptimizerSystem {
public:
    DropWave() : TestOptimizerSystem(2) {
        // The website above says this function usually has the following
        // bounds:
        Vector limits(getNumParameters());
        limits.setTo(5.12);
        setParameterLimits(-limits, limits);
    }
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        const Real dotprod = x[0] * x[0] + x[1] * x[1];
        f = -(1 + cos(12 * sqrt(dotprod))) / (0.5 * dotprod + 2);
        return 0;
    }
    Vector optimalParameters() const {
        Vector x(getNumParameters());
        x.setToZero();
        return x;
    }
};

// Looks like a curved valley.
// https://en.wikipedia.org/wiki/Test_functions_for_optimization
class Rosenbrock : public TestOptimizerSystem {
public:
    Rosenbrock(int nParameters) : TestOptimizerSystem(nParameters) {}
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        f = 0;
        for (int i = 0; i < getNumParameters() - 1; ++i) {
            f += 100.0 * square(x[i+1] - square(x[i])) + square(x[i] - 1);
        }
        return 0;
    }
    Vector optimalParameters() const {
        Vector x(getNumParameters());
        x.setTo(1);
        return x;
    }
};

// http://www.sfu.ca/~ssurjano/schwef.html
class Schwefel : public TestOptimizerSystem {
public:
    Schwefel(int nParameters) : TestOptimizerSystem(nParameters) {
        // The website above says this function usually has the following
        // bounds:
        Vector limits(nParameters);
        limits.setTo(500);
        setParameterLimits(-limits, limits);
    }
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        Real sum = 0;
        for (int i = 0; i < getNumParameters(); ++i) {
            sum += x[i] * sin(sqrt(std::abs(x[i])));
        }
        f = 418.9829 * getNumParameters() -sum;
        return 0;
    }
    Vector optimalParameters() const {
        Vector x(getNumParameters());
        x.setTo(420.9687);
        return x;
    }
};

class Easom : public TestOptimizerSystem {
public:
    Easom() : TestOptimizerSystem(2) {
        // The website above says this function usually has the following
        // bounds:
        Vector limits(getNumParameters());
        limits.setTo(100);
        setParameterLimits(-limits, limits);
    }
    int objectiveFunc(const Vector& x, bool new_parameters, Real& f) const {
        f = -cos(x[0]) * cos(x[1]) * 
            exp(-pow(x[0] - Pi, 2) - pow(x[1] - Pi, 2));
        return 0;
    }
    Vector optimalParameters() const {
        Vector x(getNumParameters());
        x.setTo(Pi);
        return x;
    }
};


#endif // SimTK_SIMMATH_OPTIMIZER_SYSTEMS_H_
