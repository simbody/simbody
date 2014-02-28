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

/** TODO: Projected Gauss Seidel impulse solver.
**/

class SimTK_SIMBODY_EXPORT PGSImpulseSolver : public ImpulseSolver {
public:
    PGSImpulseSolver() 
    :   m_convergenceTol(1e-6), m_maxIters(100), m_SOR(1.2) {}

    /** Solve. **/
    bool solve
       (int                                 phase,
        const Array_<MultiplierIndex>&      participating, // p<=m of these 
        const Matrix&                       A,     // m X m, symmetric
        const Vector&                       D,     // m, diag >= 0 added to A
        const Vector&                       rhs,   // m, RHS
        Vector&                             pi,    // m, initial guess & result
        Array_<UncondRT>&                   unconditional,
        Array_<UniContactRT>&               uniContact, // with friction
        Array_<UniSpeedRT>&                 uniSpeed,
        Array_<BoundedRT>&                  bounded,
        Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
        Array_<StateLtdFrictionRT>&         stateLtdFriction
        ) const OVERRIDE_11;

private:
    Real m_convergenceTol;
    int  m_maxIters;
    Real m_SOR; 

};

} // namespace SimTK

#endif // SimTK_SIMBODY_PGS_IMPULSE_SOLVER_H_
