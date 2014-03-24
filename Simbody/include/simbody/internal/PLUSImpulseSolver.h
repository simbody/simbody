#ifndef SimTK_SIMBODY_PLUS_IMPULSE_SOLVER_H_
#define SimTK_SIMBODY_PLUS_IMPULSE_SOLVER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Thomas Uchida, Michael Sherman                                    *
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

/** TODO: PLUS (Poisson-Lankarani-Uchida-Sherman) impulse solver.
**/

class SimTK_SIMBODY_EXPORT PLUSImpulseSolver : public ImpulseSolver {
public:
    explicit PLUSImpulseSolver(Real roll2slipTransitionSpeed) 
    :   ImpulseSolver(roll2slipTransitionSpeed,
                      1e-10, // default PLUS convergence tol
                      20),   // default Newton iteration limit
        m_minSmoothness(SqrtEps), // sharpness of smoothed discontinuities
        m_cosMaxSlidingDirChange(std::cos(Pi/6)) // 30 degrees
    {}

    /** Solve with conditional constraints. **/
    bool solve
       (int                                 phase,
        const Array_<MultiplierIndex>&      participating,
        const Matrix&                       A,
        const Vector&                       D,
        const Array_<MultiplierIndex>&      expanding,
        Vector&                             piExpand, // in/out
        Vector&                             verr, // in/out
        Vector&                             pi,  
        Array_<UncondRT>&                   unconditional,
        Array_<UniContactRT>&               uniContact,
        Array_<UniSpeedRT>&                 uniSpeed,
        Array_<BoundedRT>&                  bounded,
        Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
        Array_<StateLtdFrictionRT>&         stateLtdFriction
        ) const OVERRIDE_11;

    /** Solve with only unconditional constraints. **/
    bool solveBilateral
       (const Array_<MultiplierIndex>&      participating, // p<=m of these 
        const Matrix&                       A,     // m X m, symmetric
        const Vector&                       D,     // m, diag>=0 added to A
        const Vector&                       rhs,   // m, RHS
        Vector&                             pi     // m, unknown result
        ) const OVERRIDE_11;

    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(PLUSImpulseSolver, ActiveIndex);

private:

    // Given point P and line segment AB, find the point closest to P that lies
    // on AB, which we call Q. Returns stepLength, the ratio AQ:AB. In our case,
    // P is the origin and AB is the line segment connecting the initial and
    // final tangential velocity vectors.
    // @author Thomas Uchida
    Real calcSlidingStepLengthToOrigin(const Vec2& A, const Vec2& B, Vec2& Q)
        const;
    Real calcSlidingStepLengthToOrigin(const Vec3& A, const Vec3& B, Vec3& Q)
        const;

    // Given vectors A and B, find step length alpha such that the angle between
    // A and A+alpha*(B-A) is MaxSlidingDirChange. The solutions were generated
    // in Maple using the law of cosines, then exported as optimized code.
    // @author Thomas Uchida
    Real calcSlidingStepLengthToMaxChange(const Vec2& A, const Vec2& B) const;
    Real calcSlidingStepLengthToMaxChange(const Vec3& A, const Vec3& B) const;

    void classifyFrictionals(Array_<UniContactRT>& uniContacts) const;

    // Go through the given set of active constraints and build a reverse map
    // from the multipliers to the active index.
    void fillMult2Active(const Array_<MultiplierIndex,ActiveIndex>& active,
                         Array_<ActiveIndex,MultiplierIndex>& mult2active) const;

    // Copy the active rows and columns of A into the Jacobian. These will
    // be the right values for the linear equations, but rows for nonlinear
    // equations (sliding, impending) will get overwritten. Initialize piActive 
    // from pi.
    void initializeNewton(const Matrix&          A,
                          const Vector&          piGuess,
                          const Array_<UniContactRT>& bounded) const;

    // Given a new piActive, update the impending slip directions and calculate
    // the new err(piActive).
    void updateDirectionsAndCalcCurrentError
       (const Matrix& A, Array_<UniContactRT>& uniContact,
        const Vector& piELeft, const Vector& piActive, 
        Vector& errActive) const;

    // Replace rows of Jacobian for constraints corresponding to sliding or
    // impending slip frictional elements. This is the partial derivative of the
    // constraint error w.r.t. pi. Also set rhs m_verrActive.
    void updateJacobianForSliding(const Matrix&             A,
                                  const Array_<UniContactRT>& uniContact,
                                  const Vector& piELeft) const;

    // These are set on construction.
    Real m_minSmoothness;
    Real m_cosMaxSlidingDirChange;

    // This starts out as verr and is then reduced during each interval.
    mutable Vector m_verrLeft; // m of these
    mutable Vector m_verrExpand; // -A*piExpand for not-yet-applied piE

    // This is a subset of the given participating constraints that are
    // presently active. Only the rows and columns of A that are listed here
    // can be used (and we'll replace some of those rows). Note that a 
    // "known" unilateral contact (typically one undergoing Poisson expansion)
    // is *not* active, although its friction constraints are.
    mutable Array_<MultiplierIndex,ActiveIndex> m_active;   // na of these

    // This is the inverse mapping from m_active. Given an index into the full
    // A matrix (all proximal constraint equations, each with a Simbody-assigned
    // multiplier), this returns either the corresponding index into the 
    // m_active array, or an invalid index if this proximal constraint is not 
    // active.
    mutable Array_<ActiveIndex,MultiplierIndex> m_mult2active; // m of these

    // Each of these is indexed by ActiveIndex; they have dimension na.
    mutable Matrix m_JacActive;  // Jacobian for Newton iteration
    mutable Vector m_rhsActive;  // per-interval RHS for Newton iteration
    mutable Vector m_piActive;   // Current impulse during Newton.
    mutable Vector m_errActive;  // Error(piActive)

    mutable Matrix m_bilateralActive;  // temp for use by solveBilateral()
};

} // namespace SimTK

#endif // SimTK_SIMBODY_PLUS_IMPULSE_SOLVER_H_