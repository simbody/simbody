#ifndef SimTK_SimTKCOMMON_FUNCTION_H_
#define SimTK_SimTKCOMMON_FUNCTION_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

// Note: this file was moved from Simmath to SimTKcommon 20100601; see the
// Simmath repository for earlier history.

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include <cassert>

namespace SimTK {

/**
 * This abstract class represents a mathematical function that calculates a 
 * value of arbitrary type based on M real arguments.  The output type is set 
 * as a template argument, while the number of input components may be 
 * determined at runtime.  The name "Function" (with no trailing _) may be
 * used as a synonym for Function_<Real>.
 * 
 * Subclasses define particular mathematical functions.  Predefined subclasses 
 * are provided for several common function types: Function_<T>::Constant, 
 * Function_<T>::Linear, Function_<T>::Polynomial, and Function_<T>::Step.
 * You can define your own subclasses for other function types.  The 
 * Spline_ class also provides a convenient way to create various types of 
 * Functions.
 */
template <class T>
class Function_ {
public:
    class Constant;
    class Linear;
    class Polynomial;
    class Sinusoid;
    class Step;
    virtual ~Function_() {
    }
    /**
     * Calculate the value of this function at a particular point.
     * 
     * @param x     the Vector of input arguments. Its size must equal the 
     *              value returned by getArgumentSize().
     */
    virtual T calcValue(const Vector& x) const = 0;
    /**
     * Calculate a partial derivative of this function at a particular point.  
     * Which derivative to take is specified by listing the input components 
     * with which to take it. For example, if derivComponents=={0}, that 
     * indicates a first derivative with respective to component 0. If 
     * derivComponents=={0, 0, 0}, that indicates a third derivative with 
     * respective to component 0.  If derivComponents=={4, 7}, that indicates a 
     * partial second derivative with respect to components 4 and 7.
     * 
     * @param       derivComponents  
     *      The input components with respect to which the derivative should be
     *      taken.  Its size must be less than or equal to the value returned 
     *      by getMaxDerivativeOrder().
     * @param       x                
     *      The Vector of input arguments. Its size must equal the value 
     *      returned by getArgumentSize().
     * @return
     *      The value of the selected derivative, which is of type T.
     */
    virtual T calcDerivative(const Array_<int>& derivComponents, 
                             const Vector&      x) const = 0;

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }

    /**
     * Get the number of components expected in the input vector.
     */
    virtual int getArgumentSize() const = 0;
    /**
     * Get the maximum derivative order this Function_ object can calculate.
     */
    virtual int getMaxDerivativeOrder() const = 0;
};

/** This typedef is used for the very common case that the return type of
the Function object is Real. **/
typedef Function_<Real> Function;



/**
 * This is a Function_ subclass which simply returns a fixed value, independent
 * of its arguments.
 */
template <class T>
class Function_<T>::Constant : public Function_<T> {
public:
    /**
     * Create a Function_::Constant object.
     * 
     * @param value        the value which should be returned by calcValue();
     * @param argumentSize the value which should be returned by 
     *                     getArgumentSize(), with a default of 1.
     */
    explicit Constant(T value, int argumentSize=1) 
    :   argumentSize(argumentSize), value(value) {
    }
    T calcValue(const Vector& x) const {
        assert(x.size() == argumentSize);
        return value;
    }
    T calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
        return static_cast<T>(0);
    }
    virtual int getArgumentSize() const {
        return argumentSize;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }

private:
    const int argumentSize;
    const T value;
};

/**
 * This is a Function_ subclass whose output value is a linear function of its 
 * arguments: f(x, y, ...) = ax+by+...+c.
 */
template <class T>
class Function_<T>::Linear : public Function_<T> {
public:
    /**
     * Create a Function_::Linear object.
     * 
     * @param coefficients  
     *      The coefficients of the linear function. The number of arguments 
     *      expected by the function is equal to coefficients.size()-1.  
     *      coefficients[0] is the coefficient for the first argument, 
     *      coefficients[1] is the coefficient for the second argument, etc.
     *      The final element of coefficients contains the constant term.
     */
    explicit Linear(const Vector_<T>& coefficients) : coefficients(coefficients) {
    }
    T calcValue(const Vector& x) const {
        assert(x.size() == coefficients.size()-1);
        T value = static_cast<T>(0);
        for (int i = 0; i < x.size(); ++i)
            value += x[i]*coefficients[i];
        value += coefficients[x.size()];
        return value;
    }
    T calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
        assert(x.size() == coefficients.size()-1);
        assert(derivComponents.size() > 0);
        if (derivComponents.size() == 1)
            return coefficients(derivComponents[0]);
        return static_cast<T>(0);
    }
    virtual int getArgumentSize() const {
        return coefficients.size()-1;
    }
    int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }
private:
    const Vector_<T> coefficients;
};


/**
 * This is a Function_ subclass whose output value is a polynomial of its 
 * argument: f(x) = ax^n+bx^(n-1)+...+c.
 */
template <class T>
class Function_<T>::Polynomial : public Function_<T> {
public:
    /**
     * Create a Function_::Polynomial object.
     * 
     * @param coefficients the polynomial coefficients in order of decreasing 
     *        powers
     */
    Polynomial(const Vector_<T>& coefficients) : coefficients(coefficients) {
    }
    T calcValue(const Vector& x) const {
        assert(x.size() == 1);
        Real arg = x[0];
        T value = static_cast<T>(0);
        for (int i = 0; i < coefficients.size(); ++i)
            value = value*arg + coefficients[i];
        return value;
    }
    T calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
        assert(x.size() == 1);
        assert(derivComponents.size() > 0);
        Real arg = x[0];
        T value = static_cast<T>(0);
        const int derivOrder = (int)derivComponents.size();
        const int polyOrder = coefficients.size()-1;
        for (int i = 0; i <= polyOrder-derivOrder; ++i) {
            T coeff = coefficients[i];
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

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }
private:
    const Vector_<T> coefficients;
};


/**
 * This is a Function_ subclass whose output value is a sinusoid of its 
 * argument: f(x) = a*sin(w*x + p) where a is amplitude, w is frequency
 * in radians per unit of x, p is phase in radians.
 *
 * This is only defined for a scalar (Real) return value.
 */
template <>
class Function_<Real>::Sinusoid : public Function_<Real> {
public:
    /**
     * Create a Function::Sinusoid object, returning a*sin(w*x+p).
     * 
     * @param[in] amplitude 'a' in the above formula
     * @param[in] frequency 'w' in the above formula
     * @param[in] phase     'p' in the above formula
     */
    Sinusoid(Real amplitude, Real frequency, Real phase=0) 
    :   a(amplitude), w(frequency), p(phase) {}

    void setAmplitude(Real amplitude) {a=amplitude;}
    void setFrequency(Real frequency) {w=frequency;}
    void setPhase    (Real phase)     {p=phase;}

    Real getAmplitude() const {return a;}
    Real getFrequency() const {return w;}
    Real getPhase    () const {return p;}

    // Implementation of Function_<T> virtuals.

    virtual Real calcValue(const Vector& x) const {
        const Real t = x[0]; // we expect just one argument
        return a*std::sin(w*t + p);
    }

    virtual Real calcDerivative(const Array_<int>& derivComponents,
                                const Vector&      x) const {
        const Real t = x[0]; // time is the only argument
        const int order = derivComponents.size();
        // The n'th derivative is
        //    sign * a * w^n * sc
        // where sign is -1 if floor(order/2) is odd, else 1
        // and   sc is cos(w*t+p) if order is odd, else sin(w*t+p)
        switch(order) {
        case 0: return  a*      std::sin(w*t + p);
        case 1: return  a*w*    std::cos(w*t + p);
        case 2: return -a*w*w*  std::sin(w*t + p);
        case 3: return -a*w*w*w*std::cos(w*t + p);
        default:
            const Real sign = Real(((order/2) & 0x1) ? -1 : 1);
            const Real sc   = (order & 0x1) ? std::cos(w*t+p) : std::sin(w*t+p);
            const Real wn   = std::pow(w, order);
            return sign*a*wn*sc;
        }
    }

    virtual int getArgumentSize() const {return 1;} // just time
    virtual int getMaxDerivativeOrder() const {
        return std::numeric_limits<int>::max();
    }

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    Real calcDerivative(const std::vector<int>& derivComponents, 
                        const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }
private:
    Real a, w, p;
};

/**
 * This is a Function_ subclass whose output value y=f(x) is smoothly stepped
 * from y=y0 to y1 as its input argument goes from x=x0 to x1. This is 
 * an S-shaped function with first and second derivatives y'(x0)=y'(x1)=0
 * and y''(x0)=y''(x1)==0. The third derivative y''' exists and is continuous
 * but we cannot guarantee anything about it at the end points.
 */
template <class T>
class Function_<T>::Step : public Function_<T> {
public:
    /**
     * Create a Function_::Step object that smoothly interpolates its output
     * through a given range as its input moves through its range.
     * 
     * @param y0    Output value when (x-x0)*sign(x1-x0) <= 0.
     * @param y1    Output value when (x-x1)*sign(x1-x0) >= 0.
     * @param x0    Start of switching interval.
     * @param x1    End of switching interval.
     *
     * @tparam T    The template type is the type of y0 and y1. This must
     *              be a type that supports subtraction and scalar
     *              multiplication by a Real so that we can compute
     *              an expression like y=y0 + f*(y1-y0) for some Real scalar f.
     *
     * Note that the numerical values of x0 and x1 can be in either order
     * x0 < x1 or x0 > x1.
     */
    Step(const T& y0, const T& y1, Real x0, Real x1) 
    :   m_y0(y0), m_y1(y1), m_yr(y1-y0), m_zero(Real(0)*y0),
        m_x0(x0), m_x1(x1), m_ooxr(1/(x1-x0)), m_sign(sign(m_ooxr)) 
    {   SimTK_ERRCHK1_ALWAYS(x0 != x1, "Function_<T>::Step::ctor()",
        "A zero-length switching interval is illegal; both ends were %g.", x0);
    }

    T calcValue(const Vector& xin) const {
        SimTK_ERRCHK1_ALWAYS(xin.size() == 1,
            "Function_<T>::Step::calcValue()", 
            "Expected just one input argument but got %d.", xin.size());

        const Real x = xin[0];
        if ((x-m_x0)*m_sign <= 0) return m_y0;
        if ((x-m_x1)*m_sign >= 0) return m_y1;
        // f goes from 0 to 1 as x goes from x0 to x1.
        const Real f = stepAny(0,1,m_x0,m_ooxr, x);
        return m_y0 + f*m_yr;
    }

    T calcDerivative(const Array_<int>& derivComponents, const Vector& xin) const {
        SimTK_ERRCHK1_ALWAYS(xin.size() == 1,
            "Function_<T>::Step::calcDerivative()", 
            "Expected just one input argument but got %d.", xin.size());

        const int derivOrder = (int)derivComponents.size();
        SimTK_ERRCHK1_ALWAYS(1 <= derivOrder && derivOrder <= 3,
            "Function_<T>::Step::calcDerivative()",
            "Only 1st, 2nd, and 3rd derivatives of the step are available,"
            " but derivative %d was requested.", derivOrder);
        const Real x = xin[0];
        if ((x-m_x0)*m_sign <= 0) return m_zero;
        if ((x-m_x1)*m_sign >= 0) return m_zero;
        switch(derivOrder) {
          case 1: return dstepAny (1,m_x0,m_ooxr, x) * m_yr;
          case 2: return d2stepAny(1,m_x0,m_ooxr, x) * m_yr;
          case 3: return d3stepAny(1,m_x0,m_ooxr, x) * m_yr;
          default: assert(!"impossible derivOrder");
        }
        return NaN*m_yr; /*NOTREACHED*/
    }

    virtual int getArgumentSize() const {return 1;}
    int getMaxDerivativeOrder() const {return 3;}

    /** This provides compatibility with std::vector without requiring any 
    copying. **/
    T calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const 
    {   return calcDerivative(ArrayViewConst_<int>(derivComponents),x); }
private:
    const T    m_y0, m_y1, m_yr;   // precalculate yr=(y1-y0)
    const T    m_zero;             // precalculate T(0)
    const Real m_x0, m_x1, m_ooxr; // precalculate ooxr=1/(x1-x0)
    const Real m_sign;             // sign(x1-x0) is 1 or -1
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_FUNCTION_H_


