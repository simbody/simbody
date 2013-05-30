#ifndef SimTK_SIMMATH_SEMI_EXPLICIT_EULER_INTEGRATOR_H_
#define SimTK_SIMMATH_SEMI_EXPLICIT_EULER_INTEGRATOR_H_

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

class SemiExplicitEulerIntegratorRep;

/** This is an implementation of the fixed-step Semi-Explicit Euler method, also
known as Semi-Implicit Euler or Symplectic Euler.

This method is very popular for game engines such as ODE, because it is
a first order method that is more stable than Explicit Euler. However, the
implementation used in gaming is fully explicit, which is only correct if 
there are no velocity-dependent force elements. Coriolis forces, damping,
etc. would require an implicit solution which is generally too slow for
real time simulation.

@note There is no error estimator for this integrator so it cannot adjust the
step size. It is given a fixed step size on construction and will use that
unless you change it. The step size may be reduced to isolate events, and the
method is capable of interpolation (linear) if you want reports at shorter 
intervals than the step size.

<h3>Theory</h3>

The correct implementation of this method is described in Geometric Numerical
Integration, Hairer, Lubich & Wanner 2006, page 3. The method applies to
a system of differential equations that can be partitioned like
this: <pre>
    udot = udot(q,u)
    qdot = qdot(q,u)
         = N(q)u      in Simbody
</pre> Semi-Implicit Euler treats one variable implicitly, and the other
explicitly, producing two different possible forms: <pre>
Form 1:
    q1 = q0 + h qdot(q1,u0)   <-- implicit in q
    u1 = u0 + h udot(q1,u0)
Form 2:
    u1 = u0 + h udot(q0,u1)   <-- implicit in u
    q1 = q0 + h qdot(q0,u1)
</pre> If udot were udot(q) only, Form 2 would be <pre>
    u1 = u0 + h udot(q0)      <-- if no velocity dependence
    q1 = q0 + h N(q0)u1
</pre> That would be a correct, yet fully explicit implementation of Symplectic
Euler. Instead, the standard gaming implementation is <pre>
    u1 = u0 + h udot(q0,u0)   <-- wrong; should be implicit in u
    q1 = q0 + h N(q0)u1
</pre> This differs from a true Symplectic Euler by however much the velocity-
dependent forces change from u0 to u1. If they are very small, no great loss.
That may be less likely in internal coordinates though. Making this implicit
would be quite expensive since we don't have analytical partial derivatives
of udot available.

Form 1 on the other hand could be implemented efficiently as a true 
semi-implicit method: <pre>
    q1 = q0 + h N(q1)u0       <-- implicit in q
    u1 = u0 + h udot(q1,u0)
</pre> The directional partial derivatives D N(q)*u0 / D q are easily calculated
for the block diagonal matrix N which in many cases is just an identity matrix.
A disadvantage is that accelerations are unknown at (q0,u0) and (q1,u1); if 
those are needed (very common!) a second evaluation of the accelerations 
would be required. For that expense, a second order integration method could
have been used instead. I also don't know whether this form would have the same
stability properties in practice as the other one; it would be interesting to
experiment. However, for now we use Form 2 in the gamer's fully-explicit style.
This will not work well in the presence of huge damping forces and probably not
for very high rotation rates either!
**/

class SimTK_SIMMATH_EXPORT SemiExplicitEulerIntegrator : public Integrator {
public:
    /** Create a SemiExplicitEulerIntegrator for integrating a System with 
    fixed size steps. **/
    SemiExplicitEulerIntegrator(const System& sys, Real stepSize);
    ~SemiExplicitEulerIntegrator();
};

} // namespace SimTK

#endif // SimTK_SIMMATH_SEMI_EXPLICIT_EULER_INTEGRATOR_H_


