#ifndef SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_
#define SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-8 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

#include "AbstractIntegratorRep.h"

namespace SimTK {

/**
 * This is the private (library side) implementation of the 
 * RungeKuttaMersonIntegratorRep class which is a concrete class
 * implementing the abstract IntegratorRep.
 */

class RungeKuttaMersonIntegratorRep : public AbstractIntegratorRep {
public:
    RungeKuttaMersonIntegratorRep(Integrator* handle, const System& sys);
protected:
    bool attemptODEStep
       (Real t1, Vector& yErrEst, int& errOrder, int& numIterations);
private:    
    static const int NTemps = 3;
    Vector ytmp[NTemps];
};

} // namespace SimTK

#endif // SimTK_SIMMATH_RUNGE_KUTTA_MERSON_INTEGRATOR_REP_H_


