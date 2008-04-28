#ifndef SimTK_SIMMATH_FUNCTION_H_
#define SimTK_SIMMATH_FUNCTION_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

namespace SimTK {

/**
 * This abstract class represents a mathematical function that calculates an N component vector
 * based on M real arguments.  The number of output components is set as a template argument, while
 * the number of input components may be determined at runtime.
 * 
 * Subclasses define particular mathematical functions.  Predefined subclasses are provided for
 * several common function types: Function::Constant, Function::Linear, and Function::Polynomial.
 * You can define your own subclasses for other function types.  The Spline class also provides
 * a convenient way to create various types of Functions.
 */

template <int N>
class Function {
public:
    class Constant;
    class Linear;
    class Polynomial;
    virtual ~Function() {
    }
    /**
     * Calculate the value of this function at a particular point.
     * 
     * @param x     the Vector of input arguments.  Its size must equal the value returned by getArgumentSize().
     */
    virtual Vec<N> calcValue(const Vector& x) const = 0;
    /**
     * Calculate a partial derivative of this function at a particular point.  Which derivative to take is specified
     * by listing the input components with which to take it.  For example, if derivComponents=={0}, that indicates
     * a first derivative with respective to component 0.  If derivComponents=={0, 0, 0}, that indicates a third
     * derivative with respective to component 0.  If derivComponents=={4, 7}, that indicates a partial second derivative with
     * respect to components 4 and 7.
     * 
     * @param derivComponents  the input components with respect to which the derivative should be taken.  Its size must be
     *                         less than or equal to the value returned by getMaxDerivativeOrder().
     * @param x                the Vector of input arguments.  Its size must equal the value returned by getArgumentSize().
     */
    virtual Vec<N> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const = 0;
    /**
     * Get the number of components expected in the input vector.
     */
    virtual int getArgumentSize() const = 0;
    /**
     * Get the maximum derivative order this Function object can calculate.
     */
    virtual int getMaxDerivativeOrder() const = 0;
};

/**
 * This is a Function subclass which simply returns a fixed value, independent of its arguments.
 */

template <int N>
class Function<N>::Constant : public Function<N> {
public:
    /**
     * Create a Function::Constant object.
     * 
     * @param value        the value which should be returned by calcValue();
     * @param argumentSize the value which should be returned by getArgumentSize();
     */
    Constant(Vec<N> value, int argumentSize) : value(value), argumentSize(argumentSize) {
    }
    Vec<N> calcValue(const Vector& x) const {
        assert(x.size() == argumentSize);
        return value;
    }
    Vec<N> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        return Vec<N>(0.0);
    }
    virtual int getArgumentSize() const {
        return argumentSize;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }
private:
    const int argumentSize;
    const Vec<N> value;
};

/**
 * This is a Function subclass whose output value is a linear function of its arguments:
 * f(x, y, ...) = ax+by+...+c.
 */

template <int N>
class Function<N>::Linear : public Function<N> {
public:
    /**
     * Create a Function::Linear object.
     * 
     * @param coefficients  the coefficients of the linear function.  The number of arguments expected by the
     *                      function is equal to coefficients.size()-1.  coefficients[0] is the coefficient for
     *                      the first argument, coefficients[1] is the coefficient for the second argument, etc.
     *                      The final element of coefficients contains the constant term.
     */
    Linear(const Vector_<Vec<N> >& coefficients) : coefficients(coefficients) {
    }
    Vec<N> calcValue(const Vector& x) const {
        assert(x.size() == coefficients.size()-1);
        Vec<N> value(0);
        for (int i = 0; i < x.size(); ++i)
            value += x[i]*coefficients[i];
        value += coefficients[x.size()];
        return value;
    }
    Vec<N> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        assert(x.size() == coefficients.size()-1);
        assert(derivComponents.size() > 0);
        if (derivComponents.size() == 1)
            return coefficients(derivComponents[0]);
        return Vec<N>(0.0);
    }
    virtual int getArgumentSize() const {
        return coefficients.size()-1;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }
private:
    const Vector_<Vec<N> > coefficients;
};


/**
 * This is a Function subclass whose output value is a polynomial of its argument:
 * f(x) = ax^n+bx^(n-1)+...+c.
 */

template <int N>
class Function<N>::Polynomial : public Function<N> {
public:
    /**
     * Create a Function::Polynomial object.
     * 
     * @param coefficients the polynomial coefficients in order of decreasing powers
     */
    Polynomial(const Vector_<Vec<N> >& coefficients) : coefficients(coefficients) {
    }
    Vec<N> calcValue(const Vector& x) const {
        assert(x.size() == 1);
        Real arg = x[0];
        Vec<N> value(0);
        for (int i = 0; i < coefficients.size(); ++i)
            value = value*arg + coefficients[i];
        return value;
    }
    Vec<N> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        assert(x.size() == 1);
        assert(derivComponents.size() > 0);
        Real arg = x[0];
        Vec<N> value(0);
        int derivOrder = derivComponents.size();
        int polyOrder = coefficients.size()-1;
        for (int i = 0; i <= polyOrder-derivOrder; ++i) {
            Vec<N> coeff = coefficients[i];
            for (int j = 0; j < derivOrder; ++j)
                coeff *= polyOrder-i-j;
            value = value*arg + coeff;
        }
        return value;
    }
    virtual int getArgumentSize() const {
        return 1;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }
private:
    const Vector_<Vec<N> > coefficients;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_FUNCTION_H_


