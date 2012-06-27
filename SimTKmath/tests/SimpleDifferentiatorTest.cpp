/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * This is a test program which uses the Differentiator class in various ways.
 */

#include "SimTKmath.h"

// Just so we can get the version number:
#include "SimTKlapack.h"

#include <cstdio>
#include <iostream>

using SimTK::Real;
using SimTK::Vector;
using SimTK::Matrix;
using SimTK::Differentiator;
using std::printf;
using std::cout;
using std::endl;

// This is a system of functions of a particular set of parameters (state).
// The underlying function wants time also, so we provide that as data in
// the concrete class. Time should be set prior to calculation of the Jacobian.
class MyVectorFunc : public Differentiator::JacobianFunction {
public:
    MyVectorFunc(int nf, int ny) 
        : Differentiator::JacobianFunction(nf,ny), time(0) { }

    void setTime(Real t) {time=t;}
    Real getTime() const {return time;}

    // Must provide this pure virtual function.
    int f(const Vector& y, Vector& fy) const;
private:
    Real time;
};

// This is a single scalar function of a vector of parameters.
class MyObjectiveFunc : public Differentiator::GradientFunction {
public:
    MyObjectiveFunc(int ny) 
        : Differentiator::GradientFunction(ny), time(0) { }

    void setTime(Real t) {time=t;}
    Real getTime() const {return time;}

    // Must provide this pure virtual function.
    int f(const Vector& y, Real& fy) const;
private:
    Real time;
};

// This represents a generic scalar function of a scalar parameter,
// where the actual function has a simple C signature.
class GenericScalarFunc : public Differentiator::ScalarFunction {
    typedef Real (*CFunc)(Real);
public:
    GenericScalarFunc(CFunc cf) 
        : Differentiator::ScalarFunction(), cp(cf) { }

    // Must provide this pure virtual function.
    int f(Real x, Real& fx) const {
        fx = cp(x);
        return 0;
    }
    
    CFunc cp;
};

class SinOmegaX : public Differentiator::ScalarFunction {
public:
    SinOmegaX(Real omega) : w(omega) { }

    // Must provide this virtual function.
    int f(Real x, Real& fx) const {
        fx = std::sin(w*x);
        return 0; // success
    }
private:
    const Real w;
};


int main () {
    try {
        const Real w=3;
        SinOmegaX      sinwx(w);  // user-written class
        Differentiator dsinwx(sinwx);

        const Real x = 1.234;
        Real exact  = w*std::cos(w*x);
        Real approx = dsinwx.calcDerivative(x);

        std::printf("exact =%16.12f\n", exact);
        std::printf("approx=%16.12f err=%.3e\n", 
            approx, std::abs((approx-exact)/exact));

        return 0;
    } 
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}
