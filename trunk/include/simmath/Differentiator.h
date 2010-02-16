#ifndef SimTK_DIFFERENTIATOR_H_
#define SimTK_DIFFERENTIATOR_H_

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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath numerical differentiation tools.
 */

#include "SimTKcommon.h"
#include "internal/common.h"
#include "SimTKcommon/internal/BigMatrix.h"

namespace SimTK {


/**
 * Given a function f(y), where f, y or both can be vectors, 
 * calculate the derivative (gradient, Jacobian) df/dy.
 * 
 * Calculation is done using numerical differencing, which should be
 * considered a last resort for cases in which the analytic
 * derivative is unavailable. (Note that you can obtain an
 * analytic gradient automatically from the source code
 * for f using complex step derivatives, ADIFOR, etc.)
 */
class SimTK_SIMMATH_EXPORT Differentiator {
public:
    // This are local classes within Differentiator; defined below.
    class ScalarFunction;   // ordinary scalar function of a scalar
    class GradientFunction; // scalar function of vector
    class JacobianFunction; // vector function of vector
    class Function;         // abstraction of the above

    // These are the exceptions that can be thrown by this class.
    class OpNotAllowedForFunctionOfThisShape;
    class UserFunctionThrewAnException;
    class UserFunctionReturnedNonzeroStatus;
    class UnknownMethodSpecified;


    enum Method {
        UnspecifiedMethod=0,
        ForwardDifference=1,
        CentralDifference=2
    };
    static bool        isValidMethod(Method);
    static const char* getMethodName(Method);
    static int         getMethodOrder(Method);

    virtual ~Differentiator();
    explicit Differentiator(const Function& f, 
                            Method          defaultMethod=UnspecifiedMethod);

    // You can change the default method; normally it is ForwardDifference.
    // If you set it to 'UnspecifiedMethod' it goes back to the original default.
    Differentiator& setDefaultMethod(Method);
    Method          getDefaultMethod() const;

    // These are the real routines, which are efficient and flexible
    // but somewhat messy to use.
    void calcDerivative(Real y0, Real fy0, Real& dfdy, 
                        Method=UnspecifiedMethod) const;
    void calcGradient  (const Vector& y0, Real fy0, Vector& gf,
                        Method=UnspecifiedMethod) const;
    void calcJacobian  (const Vector& y0, const Vector& fy0, Matrix& dfdy,
                        Method=UnspecifiedMethod) const;

    // These provide a simpler though less efficient interface. They will
    // do some heap allocation, and will make an initial unperturbed call
    // to the user function.
    Real   calcDerivative(Real          y0, Method=UnspecifiedMethod) const;
    Vector calcGradient  (const Vector& y0, Method=UnspecifiedMethod) const;
    Matrix calcJacobian  (const Vector& y0, Method=UnspecifiedMethod) const;

    // Statistics (mutable)
    void resetAllStatistics();                 // reset all stats to zero
    int getNumDifferentiations() const;        // total # calls of calcWhatever
    int getNumDifferentiationFailures() const; // # of those that failed
    int getNumCallsToUserFunction() const;     // total # calls to user function

    // This is a local class.
    class DifferentiatorRep;
private:
    // opaque implementation for binary compatibility
    DifferentiatorRep* rep;

private:
    //OBSOLETE NAMES
    int getNDifferentiations() const {return getNumDifferentiations();}
    int getNDifferentiationFailures() const {return getNumDifferentiationFailures();}
    int getNCallsToUserFunction() const {return getNumCallsToUserFunction();}
};

/**
 * This abstract class defines a function to be differentiated (repeatedly)
 * by a Differentiator object. Users should not access this class directly;
 * instead, use one of the specialized function classes ScalarFunction,
 * GradientFunction, or JacobianFunction depending on the type of function
 * you want to differentiate.
 *
 * The Differentiator class will assume the function is calculated to
 * about machine accuracy unless told otherwise. But if f is the result of some
 * approximate calculation (for example, it came from another Differentiator
 * approximation, or from numerical integration),  we will need to know that in
 * order to have a reasonable crack at calculating df.
 */
class SimTK_SIMMATH_EXPORT Differentiator::Function {
public:
    Function& setNumFunctions(int);
    Function& setNumParameters(int);
    Function& setEstimatedAccuracy(Real);

    // These values are fixed after construction.
    int  getNumFunctions()  const;
    int  getNumParameters() const;
    Real getEstimatedAccuracy() const; // approx. "roundoff" in f calculation

    // Statistics (mutable)
    void resetAllStatistics();
    int getNumCalls()    const; // # evaluations of this function since reset
    int getNumFailures() const; // # of calls which failed

    // This is the declaration of a local class name.
    class FunctionRep;
protected:
    Function();
    ~Function();

    // opaque implementation for binary compatibility
    FunctionRep* rep;

private:
    // suppress copy constructor and copy assignment
    Function(const Function&);
    Function& operator=(const Function&);

private:
    //OBSOLETE NAMES
    Function& setNFunctions(int n) {return setNumFunctions(n);}
    Function& setNParameters(int n) {return setNumParameters(n);}
    int  getNFunctions()  const {return getNumFunctions();}
    int  getNParameters() const {return getNumParameters();}
    int getNCalls()    const {return getNumCalls();}
    int getNFailures() const {return getNumFailures();}



friend class Differentiator;
};

/**
 * Derive a concrete class from this one if you have a scalar function
 * of a single scalar variable that you want to differentiate.
 */
class SimTK_SIMMATH_EXPORT Differentiator::ScalarFunction : public Differentiator::Function {
public:
    virtual int f(Real x, Real& fx) const=0;

protected:
    explicit ScalarFunction(Real acc=-1);
    virtual ~ScalarFunction() { }

private:
    // suppress copy constructor and copy assignment
    ScalarFunction(const Function&);
    ScalarFunction& operator=(const Function&);
};

/**
 * Derive a concrete class from this one if you have a scalar function
 * of multiple variables that you want to differentiate. This is the typical
 * form for an optimization objective function, for example.
 */
class SimTK_SIMMATH_EXPORT Differentiator::GradientFunction : public Differentiator::Function {
public:
    virtual int f(const Vector& y, Real& fy) const=0;

protected:
    explicit GradientFunction(int ny=-1, Real acc=-1);
    virtual ~GradientFunction() { }

private:
    // suppress copy constructor and copy assignment
    GradientFunction(const GradientFunction&);
    GradientFunction& operator=(const GradientFunction&);
};

/**
 * Derive a concrete class from this one if you have a set of functions
 * (i.e., a vector-valued function) of multiple variables that you want
 * to differentiate. This is the typical form for a multibody system, for example.
 */
class SimTK_SIMMATH_EXPORT Differentiator::JacobianFunction : public Differentiator::Function {
public:
    virtual int f(const Vector& y, Vector& fy) const=0;

protected:
    explicit JacobianFunction(int nf=-1, int ny=-1, Real acc=-1); 
    virtual ~JacobianFunction() { }

private:
    // suppress copy constructor and copy assignment
    JacobianFunction(const JacobianFunction&);
    JacobianFunction& operator=(const JacobianFunction&);
};

} // namespace SimTK

#endif // SimTK_DIFFERENTIATOR_H_
