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
 * the number of input components may be determined at runtime.  Subclasses define particular
 * mathematical functions.
 */

template <int N>
class Function {
public:
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

} // namespace SimTK

#endif // SimTK_SIMMATH_FUNCTION_H_


