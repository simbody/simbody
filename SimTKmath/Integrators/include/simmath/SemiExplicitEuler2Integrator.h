#ifndef SimTK_SIMMATH_SEMI_EXPLICIT_EULER_2_INTEGRATOR_H_
#define SimTK_SIMMATH_SEMI_EXPLICIT_EULER_2_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/Integrator.h"

namespace SimTK {

class SemiExplicitEuler2IntegratorRep;

/** This is an implementation of a variable-step, first-order semi-explicit 
Euler method, also known as semi-implicit Euler or symplectic Euler.

The fixed-step, first-order semi-explicit Euler method is very popular for game 
engines such as ODE, and Simbody offers that method as 
SemiExplicitEulerIntegrator. See the documentation there for a discussion and 
theory. This integrator is an attempt to extend that method to give some
error control. For equivalent performance to the fixed-step method, this one
needs to run at twice the step size (because it takes two substeps). Be sure
to set its maximum allowable step large enough, and use a very loose accuracy
setting.

<h3>Theory</h3>

See SemiExplicitEulerIntegrator for the theory behind the underlying first
order method. Here we use Richardson Extrapolation (step doubling) to get an
error estimate. See Solving Ordinary Differential Equations I: Nonstiff 
Problems, 2nd rev. ed., Hairer, Norsett, & Wanner, section II.4, pp 164-165.
Note that we are not using the local extrapolation trick to get a higher-order
result; see below for why.

The idea is to take two steps of length h, and a single big step of length 2h, 
and combine the results. Starting at t=t0 we can get to t2=t0+2h two different
ways: <pre>
    u1 = u0 + h udot(q0,u0)              ub = u0 + 2h udot(q0,u0)
    q1 = q0 + h N(q0)u1         and      qb = q0 + 2h N(q0)ub
    u2 = u1 + h udot(q1,u1)                  (b = big step)
    q2 = q1 + h N(q1)u2
</pre> Note that we can use udot(q0,u0) in both methods so it doesn't cost much
to calculate both solutions y2(t2)=(q2,u2) and yb(t2)=(qb,ub).

Now assume the true solution at t2 is y(t2)=(q,u). Then <pre>
    y = y2 + 2C h^2  + O(h^3)
    y = yb + C(2h)^2 + O(h^3)
</pre> where C is some unknown constant of the underlying method. Ignoring the
3rd order terms, subtract these two equations to get: <pre>
    y2-yb = 2Ch^2
</pre> which is the 2nd order error term for y2. We can use that for error 
control. 

Hairer suggests using that error term to improve on y2 to get a better estimate 
for y (this is called "local extrapolation") : <pre>
    y2' = y2 + (y2-yb) = 2 y2 - yb
    y   = y2' + O(h^3)                  <-- we don't actually do this
</pre> y2' is a 2nd order accurate estimate for y(t2). But ... we don't have an
error estimate for the revised value and we don't actually know that it is more
accurate. In fact, although this sounds lovely in theory, in practice we will 
often be running at close to a stability boundary and using the "big step" value
to update the result causes stability problems. So we'll skip the local 
extrapolation and stick with the first order result from the two half-steps.

@author Michael Sherman
**/

class SimTK_SIMMATH_EXPORT SemiExplicitEuler2Integrator : public Integrator {
public:
    /** Create a SemiExplicitEuler2Integrator for integrating a System with 
    variable size steps. **/
    SemiExplicitEuler2Integrator(const System& sys);
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SEMI_EXPLICIT_EULER_2_INTEGRATOR_H_


