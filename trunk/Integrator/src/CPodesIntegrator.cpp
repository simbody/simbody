/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
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

#include "simmath/CPodesIntegrator.h"

using namespace SimTK;

/**
 * Create a CPodesIntegrator for integrating a System.  The nonlinear system iteration type is chosen automatically
 * based on the linear multistep method: Newton iteration for BDF (the default), and functional iteration for Adams.
 */

CPodesIntegrator::CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method) : rep(this, sys, method), Integrator(rep) {
}

/**
 * Create a CPodesIntegrator for integrating a System.
 */

CPodesIntegrator::CPodesIntegrator(const System& sys, CPodes::LinearMultistepMethod method, CPodes::NonlinearSystemIterationType iterationType) : rep(this, sys, method, iterationType), Integrator(rep) {
}

/**
 * CPODES provides its own mechanism for projecting the system onto the constraint manifold.  By default,
 * CPodesIntegrator uses the System's project() method for doing projection, which is usually more
 * efficient.  Invoking this method tells it to use the CPODES mechanism instead.
 * 
 * This method must be invoked before the integrator is initialized.  Invoking it after initialization
 * will produce an exception.
 */

void CPodesIntegrator::setUseCPodesProjection() {
    rep.setUseCPodesProjection();
}

