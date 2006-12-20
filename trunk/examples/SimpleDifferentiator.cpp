/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDERS, OR CONTRIBUTORS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * This is a test program which uses the Differentiator class in various ways.
 */

//#define SimTK_USE_STATIC_LIBRARIES

#include "Simmath.h"
#include "Differentiator.h"

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
