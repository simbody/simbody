#ifndef SimTK_SIMBODY_PGS_IMPULSE_SOLVER_H_
#define SimTK_SIMBODY_PGS_IMPULSE_SOLVER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "simbody/internal/ImpulseSolver.h"

namespace SimTK {

/** Projected Gauss Seidel impulse solver.
Finds a solution to
<pre>
    [A+D] (piExpand+piUnknown) = v
    subject to
    piExpand <= 0 (given)

    piUnknown_z <= 0         unilateral contact must push, not pull
    ||(piUnknown_x,piUnknown+y)|| <= -mu*pi_z   unilateral friction cone
    where pi=piExpand+piUnknown.

    piUnknown_speed <= 0           one-way ratchets
    lb <= piUnknown_bounded <= ub  torque-limited motor

    ||piUnknown_F|| <= mu*||piUnknown_N||  friction in bilateral constraint
    ||piUnknown_F|| <= mu*N                friction with known normal force N
</pre>
When piUnknown_z[k] hits its upper limit of 0, we must have v_z[k] >= 0
(contact surfaces separating). We don't explicitly enforce that here; it 
depends on all diag(A)[z[k]] > 0. That means that if v_z[k]<0 we could improve
the solution by making piUnknown_z[k] negative, so it wouldn't have hit the
limit.
**/

class SimTK_SIMBODY_EXPORT PGSImpulseSolver : public ImpulseSolver {
public:
    explicit PGSImpulseSolver(Real roll2slipTransitionSpeed) 
    :   ImpulseSolver(roll2slipTransitionSpeed,
                      1e-6, // default PGS convergence tolerance
                      100), // default PGS max number iterations
        m_SOR(1.2) {}

    /** Solve with conditional constraints. In the common underdetermined
    case (redundant contact) we will return the first solution encountered but
    it is unlikely to be the best possible solution. **/
    bool solve
       (int                                 phase,
        const Array_<MultiplierIndex>&      participating,
        const Matrix&                       A,
        const Vector&                       D, 
        const Array_<MultiplierIndex>&      expanding, // nx<=m of these 
        Vector&                             piExpand,
        Vector&                             verrStart, // in/out
        Vector&                             verrApplied, // in/out
        Vector&                             pi, 
        Array_<UncondRT>&                   unconditional,
        Array_<UniContactRT>&               uniContact,
        Array_<UniSpeedRT>&                 uniSpeed,
        Array_<BoundedRT>&                  bounded,
        Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
        Array_<StateLtdFrictionRT>&         stateLtdFriction
        ) const OVERRIDE_11;

    /** Solve with only unconditional constraints. In the underdetermined
    case we return one of the possible solutions but it will not in general
    be the least squares one. In the overdetermined, inconsistent case we
    will iterate for a long time and may converge on the least-error solution
    but cannot guarantee that. **/
    bool solveBilateral
       (const Array_<MultiplierIndex>&      participating, // p<=m of these 
        const Matrix&                       A,     // m X m, symmetric
        const Vector&                       D,     // m, diag>=0 added to A
        const Vector&                       rhs,   // m, RHS
        Vector&                             pi     // m, unknown result
        ) const OVERRIDE_11;

private:
    Real m_SOR; 
};

} // namespace SimTK

#endif // SimTK_SIMBODY_PGS_IMPULSE_SOLVER_H_
