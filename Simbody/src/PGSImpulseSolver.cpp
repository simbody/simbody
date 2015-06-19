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
 * Contributors: Thomas Uchida, John Hsu                                      *
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

#include "simbody/internal/common.h"
#include "simbody/internal/ImpulseSolver.h"
#include "simbody/internal/PGSImpulseSolver.h"

#include <algorithm>

#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// Local utilities.
namespace {

// Calculate (A+D)[row]*pi, but only looking at the given columns.
// Vectors must be contiguous, Matrix must be packed and in column order (i.e.
// columns are contiguous.) So A(r,c) = A[r + c*m].
Real doRowSum(const Array_<MultiplierIndex>& columns,
              const MultiplierIndex&         row,
              const Matrix&                  A,
              const Vector&                  D,
              const Vector&                  pi)
{
    assert(pi.hasContiguousData());
    const Real* pip = &pi[0];
    const int m = A.nrow();

    assert(A.hasContiguousData()); // packed
    assert(A(0).hasContiguousData()); // in column order
    const Real* Ap = &A(0,0);

    Real rowSum = 0;
    for (unsigned c=0; c < columns.size(); ++c) {
        const MultiplierIndex cx = columns[c];
        const Real* cp = Ap + cx*m; // point to start of column
        rowSum += cp[row]*pip[cx];
    }
    if (D.size()) {
        assert(D.hasContiguousData());
        const Real* Dp = &D[0];
        rowSum += Dp[row]*pip[row];
    }
    return rowSum;
}

// Calculate sums=(A+D)[rows]*pi, but only looking at the given columns.
// We expect that A is in column order so we'll work down the rows before
// we switch columns.
// Vectors must be contiguous, Matrix must be packed and in column order (i.e.
// columns are contiguous.) So A(r,c) = A[r + c*m].
void doRowSums(const Array_<int>& columns, // these are MultiplierIndex ints
               const Array_<int>& rows,
               const Matrix&      A, 
               const Vector&      D,
               const Vector&      pi,
               Array_<Real>&      sums)
{
    assert(pi.hasContiguousData());
    const Real* pip = &pi[0];
    const int m = A.nrow(), n = A.ncol();

    assert(A.hasContiguousData()); // packed
    assert(A(0).hasContiguousData()); // in column order
    const Real* Ap = &A(0,0);

    sums.resize(rows.size()); sums.fill(Real(0));
    for (unsigned c=0; c < columns.size(); ++c) {
        const int cx = columns[c];
        const Real* cp = Ap + cx*m; // point to start of column
        for (unsigned i=0; i<rows.size(); ++i)
            sums[i] += cp[rows[i]]*pip[cx];
    }
    if (D.size()) {
        assert(D.hasContiguousData());
        const Real* Dp = &D[0];
        for (unsigned i=0; i<rows.size(); ++i)
            sums[i] += Dp[rows[i]]*pip[rows[i]];
    }
}

// Given a rowSum, update one element of pi and return the squared error.
// If the corresponding diagonal of A is nonpositive, we will quietly skip
// the update.
inline Real doUpdate(const MultiplierIndex& row,
                     const Matrix&          A,
                     const Vector&          D,
                     const Vector&          rhs,
                     const Real&            SOR, // successive over relaxation
                     const Real&            rowSum,
                     Vector&                pi)
{
    Real Arr = A(row,row);
    if (D.size()) Arr += D[row];
    const Real er = rhs[row]-rowSum;
    if (Arr > Real(0))
        pi[row] += SOR * er/Arr;
    return square(er);
}

// Same but now we're doing multiple row updates and return the sum of the
// squared errors for those rows.
Real doUpdates(const Array_<int>& rows, // These are MultiplierIndex ints
               const Matrix&                  A,
               const Vector&                  D,
               const Vector&                  rhs,
               const Real&                    SOR,
               const Array_<Real>&            rowSums,
               Vector&                        pi)
{
    const bool hasDiag = (D.size() > 0);
    Real er2 = 0;
    for (unsigned i=0; i<rows.size(); ++i) {
        const MultiplierIndex row(rows[i]);
        Real Arr = A(row,row);
        if (hasDiag) Arr += D[row];
        const Real er = rhs[row]-rowSums[i];
        if (Arr > Real(0))
            pi[row] += SOR * er/Arr;
        er2 += square(er);
    }
    return er2;
}

// Multiply the active entries of a row of the full matrix A (mXm) by a sparse,
// full-length (m) column containing only the indicated non-zero entries. 
// Useful for A[r]*piExpand.
static Real multRowTimesSparseCol(const Matrix& A, MultiplierIndex row, 
           const Array_<MultiplierIndex>& nonZero,
           const Vector& sparseCol) 
{
    const RowVectorView Ar = A[row];
    Real result = 0;
    for (unsigned nz(0); nz < nonZero.size(); ++nz) {
        const MultiplierIndex mx = nonZero[nz];
        result += Ar[mx] * sparseCol[mx];
    }
    return result;
}

/** Given a unilateral multiplier pi and its sign convention, ensure that
sign*pi<=0. Return true if any change is made. **/
inline ImpulseSolver::UniCond boundUnilateral(Real sign, Real& pi) {
    assert(sign==1 || sign==-1);
    if (sign*pi > 0) {pi=0; return ImpulseSolver::UniOff;}
    return ImpulseSolver::UniActive;
}

/** Given a scalar pi, ensure that lb <= pi <= ub by moving pi to the nearest
bound if necessary. Return true if any change is made. **/
inline ImpulseSolver::BndCond boundScalar(Real lb, Real& pi, Real ub) {
    assert(lb <= ub);
    if      (pi > ub) {pi=ub; return ImpulseSolver::SlipHigh;}
    else if (pi < lb) {pi=lb; return ImpulseSolver::SlipLow;}
    return ImpulseSolver::Engaged;
}

/** Given an index set IV, ensure that ||pi[IV]|| <= maxLen by scaling the
vector to that length if necessary. Return true if any change is made. **/
ImpulseSolver::FricCond 
boundVector(Real maxLen, const Array_<MultiplierIndex>& IV, Vector& pi) {
    assert(maxLen >= 0);
    const Real maxLen2 = square(maxLen);
    Real piNorm2 = 0;
    for (unsigned i=0; i<IV.size(); ++i) piNorm2 += square(pi[IV[i]]);
    if (piNorm2 <= maxLen2) 
        return ImpulseSolver::Rolling;
    const Real scale = std::sqrt(maxLen2/piNorm2); // 0 <= scale < 1
    for (unsigned i=0; i<IV.size(); ++i) pi[IV[i]] *= scale;
    return ImpulseSolver::Sliding;
}

/** Given index set IN identifying the components of the normal force vector,
and index set IF identifying the components of the friction vector, ensure
that ||pi[IF]|| <= mu*||pi[IN]|| by scaling the friction vector if necessary.
Return true if any change is made. **/
ImpulseSolver::FricCond 
boundFriction(Real mu, 
              const Array_<int>& IN, // these are MultiplierIndex ints 
              const Array_<int>& IF, 
              Vector& pi) {
    assert(mu >= 0);
    Real N2=0, F2=0; // squares of normal and friction force magnitudes
    for (unsigned i=0; i<IN.size(); ++i) N2 += square(pi[IN[i]]);
    for (unsigned i=0; i<IF.size(); ++i) F2 += square(pi[IF[i]]);
    const Real mu2N2 = mu*mu*N2;
    if (F2 <= mu2N2) 
        return ImpulseSolver::Rolling;
    const Real scale = std::sqrt(mu2N2/F2); // 0 <= scale < 1
    for (unsigned i=0; i<IF.size(); ++i) pi[IF[i]] *= scale;
    return ImpulseSolver::Sliding;
}
}

namespace SimTK {



//==============================================================================
//                   PROJECTED GAUSS SEIDEL IMPULSE SOLVER
//==============================================================================
/* 
We are given 
    - A, square matrix of dimension m 
    - v, constraint velocity vector (length m)
    - pi, solution vector with initial value pi=pi0 (length m)
representing m scalar constraint equations A[i]*pi=rhs[i].

A smaller square "participating" subset may be selected via
    - I, selection index set, a p-element subset of IA={1,...,m}

The selected subset I is partitioned into four disjoint index sets
    - IU: Unconditional
    - IC: Unilateral contact, optionally with planar friction
    - IS: Unilateral speed constraint
    - IB: Bounded scalar constraint
    - IS: State-limited friction
    - IF: Constraint-limited friction

Each unconditional constraint k provides
    - a unique index set of 1-6 multipliers IU_k from IU

Each unilateral contact k provides
    - a unique normal multiplier index
    - whether the normal force is known (expander) or unknown (participater)
    - if known, then the value of the normal force
    - optionally two friction multipliers
    - the effective coefficient of friction mu

Each unilateral speed constraint k provides
    - a single constraint index

Each bounded scalar constraint k provides
    - a single constraint index iB_k from IB, and 
    - effective lower and upper bounds lb_k, ub_k.

Each state-limited friction constraint k specifies 
    - a unique index set of 1-3 distinct constraints IS_k from IS,
    - a nonnegative scalar N_k specifing the limiting normal force, as
      determined from the state and passed in to this method
    - the effective coefficient of friction mu.

Each constraint-limited friction constraint k specifies 
    - a unique index set of 1-3 distinct friction constraints IF_k from IF, 
    - an index set of 1-3 distinct normal constraints IN_k from IU,
    - the effective coefficient of friction mu.

Given those inputs, we attempt to solve: 
    A[I,I] w[I] = b[I]
    subject to lb_k <= w[iB_k] <= ub_k       for bounded constraints k in IB
    and        ||w[IV_k]|| <= L_k            for vector constraints k in IV
    and        ||w[IF_k]|| <= mu*||w[IN_k]|| for friction constraints k in IF

Implicitly, complementarity conditions must hold:
    w_i in interior of constraint -> A[i]*w == b[i]
    w_i on boundary of constraint -> A[i]*w != b[i]

*/


//------------------------------------------------------------------------------
//                                 SOLVE
//------------------------------------------------------------------------------
bool PGSImpulseSolver::
solve(int                                 phase,
      const Array_<MultiplierIndex>&      participating, // p<=m of these 
      const Matrix&                       A,     // m X m, symmetric
      const Vector&                       D,     // m, diag >= 0 added to A
      const Array_<MultiplierIndex>&      expanding,
      Vector&                             piExpand,   // m
      Vector&                             verrStart, // m, in/out
      Vector&                             verrApplied, // m
      Vector&                             pi,         // m, piUnknown
      Array_<UncondRT>&                   unconditional,
      Array_<UniContactRT>&               uniContact, // with friction
      Array_<UniSpeedRT>&                 uniSpeed,
      Array_<BoundedRT>&                  bounded,
      Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
      Array_<StateLtdFrictionRT>&         stateLtdFriction
      ) const 
{
    SimTK_DEBUG("\n-----------------\n");
    SimTK_DEBUG(  "START PGS SOLVER:\n");
    ++m_nSolves[phase];

#ifndef NDEBUG
   {FactorQTZ fac(A);
    cout << "A=" << A; cout << "D=" << D; 
    cout << "verrStart=" << verrStart << endl;
    cout << "verrApplied=" << verrApplied << endl;
    cout << "expanding mx=" << expanding << endl;
    cout << "piExpand=" << piExpand << endl;
    Vector verrDbg, x;
    verrDbg = verrStart;
    if (verrApplied.size()) verrDbg += verrApplied;
    fac.solve(verrDbg, x); 
    cout << "x=" << x << endl;
    cout << "resid=" << A*x-verrDbg << endl;}
#endif

    const int m=A.nrow();
    assert(A.ncol()==m); assert(D.size()==m);
    assert(verrStart.size()==m); 
    assert(verrApplied.size()==0 || verrApplied.size()==m);
    assert(piExpand.size()==m); 

    const int p = (int)participating.size();
    const int nx = (int)expanding.size();
    assert(p<=m); assert(nx<=m);
    
    pi.resize(m);
    pi.setToZero(); // Use this for piUnknown

    // If there are applied forces, add them to the rhs.
    if (verrApplied.size()) 
        verrStart += verrApplied;

    // Move expansion impulse to RHS. We will always apply the full expansion
    // impulse in one interval in this solver.
    if (nx) 
        for (MultiplierIndex mx(0); mx < m; ++mx) {
            verrStart[mx] -= multRowTimesSparseCol(A,mx,expanding,piExpand)
                             + D[mx]*piExpand[mx];
        }

    // Now rhs = verrStart + verrApplied - [A+D]*piExpand.
    #ifndef NDEBUG
    printf("PGS::solve(): using verr="); cout << verrStart << endl;
    #endif


    // Partitions of selected subset.
    const int mUncond   = (int)unconditional.size();
    const int mUniSpeed = (int)uniSpeed.size();
    const int mBounded  = (int)bounded.size();
    // State limited friction has no dependence on unknown multipliers.
    const int mStateLtd = (int)stateLtdFriction.size();
    // Must do unilateral friction and constraint-limited friction last because
    // they depend on normal multipliers.
    const int mUniCont  = (int)uniContact.size();
    const int mConsLtd  = (int)consLtdFriction.size();

    // If debugging, check for consistent constraint equation count.
    #ifndef NDEBUG
    {int mCount = mUniSpeed + mBounded; // 1 each
    for (int k=0; k<mUncond; ++k)
        mCount += unconditional[k].m_mults.size();
    for (int k=0; k<mUniCont; ++k) {
        if (uniContact[k].m_type==Observing)
            continue; // neither normal nor friction participate
        if (uniContact[k].m_type==Participating)
            ++mCount; // normal participates
        if (uniContact[k].hasFriction())
            mCount += 2; // friction participates even if normal is Known
    }
    for (int k=0; k<mStateLtd; ++k)
        mCount += stateLtdFriction[k].m_Fk.size();
    for (int k=0; k<mConsLtd; ++k)
        mCount += consLtdFriction[k].m_Fk.size();
    assert(mCount == p);}
    #endif

    if (p == 0) {
        SimTK_DEBUG1("PGS %d: nothing to do; converged in 0 iters.\n", phase);
        // Returning pi=0; can still have piExpand!=0 so verr is updated.
        return true;
    }

    // Track total error for all included equations, and the error for just
    // those equations that are being enforced.
    bool converged = false;
    Real normRMSall = Infinity, normRMSenf = Infinity, sor = m_SOR;
    Real prevNormRMSenf = NaN;
    int its = 1;
    Array_<Real> rowSums; // handy temp
    for (; its <= m_maxIters; ++its) {
        ++m_nIters[phase];
        Real sum2all = 0, sum2enf = 0; // track solution errors
        prevNormRMSenf = normRMSenf;

        // UNCONDITIONAL: these are always on.
        for (int k=0; k < mUncond; ++k) {
            const UncondRT& rt = unconditional[k];
            doRowSums(participating,rt.m_mults,A,D,pi,rowSums);
            const Real er2=doUpdates(rt.m_mults,A,D,verrStart,sor,rowSums,pi);
            sum2all += er2; sum2enf += er2;
        }

        // UNILATERAL CONTACT NORMALS. Do all of these before any friction.
        for (int k=0; k < mUniCont; ++k) {
            UniContactRT& rt = uniContact[k];
            if (rt.m_type != Participating)
                continue;
            const MultiplierIndex Nk = rt.m_Nk;
            const Real rowSum=doRowSum(participating,Nk,A,D,pi);
            const Real er2=doUpdate(Nk,A,D,verrStart,sor,rowSum,pi);
            sum2all += er2;
            rt.m_contactCond = boundUnilateral(rt.m_sign, pi[Nk]);
            if (rt.m_contactCond == UniActive)
                sum2enf += er2;
        }

        // UNILATERAL CONTACT FRICTION. These are limited by the normal
        // multiplier or by a known normal force during Poisson expansion.
        for (int k=0; k < mUniCont; ++k) {
            UniContactRT& rt = uniContact[k];
            if (rt.m_type == Observing || !rt.hasFriction())
                continue;
            const MultiplierIndex Nk = rt.m_Nk;
            const Array_<MultiplierIndex>& Fk = rt.m_Fk;
            doRowSums(participating,Fk,A,D,pi,rowSums);
            const Real er2=doUpdates(Fk,A,D,verrStart,sor,rowSums,pi);
            sum2all += er2;
            Real N = std::abs(pi[Nk] + piExpand[Nk]);
            rt.m_frictionCond=boundVector(rt.m_effMu*N, Fk, pi);
            if (rt.m_frictionCond==Rolling)
                sum2enf += er2;
        }

        // BOUNDED: conditional scalar constraints with constant bounds
        // on resulting pi.
        for (int k=0; k < mBounded; ++k) {
            BoundedRT& rt = bounded[k];
            const MultiplierIndex rx = rt.m_ix;
            const Real rowSum=doRowSum(participating,rx,A,D,pi);
            const Real er2=doUpdate(rx,A,D,verrStart,sor,rowSum,pi);
            sum2all += er2;
            rt.m_boundedCond=boundScalar(rt.m_lb, pi[rx], rt.m_ub);
            if (rt.m_boundedCond == Engaged)
                sum2enf += er2;
        }

        // STATE LIMITED FRICTION: a set of constraint equations forming a 
        // vector whose maximum length is limited.
        for (int k=0; k < mStateLtd; ++k) {
            StateLtdFrictionRT& rt = stateLtdFriction[k];
            const Array_<MultiplierIndex>& Fk = rt.m_Fk;
            doRowSums(participating,Fk,A,D,pi,rowSums);
            const Real localEr2=doUpdates(Fk,A,D,verrStart,sor,rowSums,pi);
            sum2all += localEr2;
            rt.m_frictionCond=boundVector(rt.m_effMu*rt.m_knownN, Fk, pi);
            if (rt.m_frictionCond==Rolling)
                sum2enf += localEr2;
        }

        // CONSTRAINT LIMITED FRICTION: a set of constraint equations forming 
        // a vector whose maximum length is limited by the norm of other 
        // multipliers pi.
        for (int k=0; k < mConsLtd; ++k) {
            ConstraintLtdFrictionRT& rt = consLtdFriction[k];
            const Array_<int>& Fk = rt.m_Fk; // friction components
            const Array_<int>& Nk = rt.m_Nk; // normal components
            doRowSums(participating,Fk,A,D,pi,rowSums);
            const Real localEr2=doUpdates(Fk,A,D,verrStart,sor,rowSums,pi);
            sum2all += localEr2;
            rt.m_frictionCond=boundFriction(rt.m_effMu,Nk,Fk,pi);
            if (rt.m_frictionCond==Rolling)
                sum2enf += localEr2;
        }
        normRMSall = std::sqrt(sum2all/p);
        normRMSenf = std::sqrt(sum2enf/p);

        const Real rate = normRMSenf/prevNormRMSenf;

        if (rate > 1) {
            SimTK_DEBUG3("GOT WORSE@%d: sor=%g rate=%g\n", its, sor, rate);
            if (sor > .1)
                sor = std::max(.8*sor, .1);
        } 

        #ifndef NDEBUG
        printf("%d/%d: EST rmsAll=%g rmsEnf=%g rate=%g\n", phase, its,
                     normRMSall, normRMSenf, 
                     normRMSenf/prevNormRMSenf);
        #endif
        //#ifdef NDEBUG // i.e., NOT debugging (TODO)
        //if (its > 90)
        //    printf("%d/%d: EST rmsAll=%g rmsEnf=%g rate=%g\n", phase, its,
        //             normRMSall, normRMSenf, 
        //             normRMSenf/prevNormRMSenf);
        //#endif
        if (normRMSenf < m_convergenceTol) //TODO: add failure-to-improve check
        {
            SimTK_DEBUG3("PGS %d converged to %g in %d iters\n", 
                         phase, normRMSenf, its);
            converged = true;
            break;
        }
        #ifndef NDEBUG
        cout << "pi=" << pi << " err=" << normRMSenf << " rate=" << rate << endl;
        #endif
    }

    if (!converged) {
        SimTK_DEBUG3("PGS %d CONVERGENCE FAILURE: %d iters -> norm=%g\n",
               phase, its, normRMSenf);
        ++m_nFail[phase];
    }

    verrStart -= A*pi;
    verrStart -= D.elementwiseMultiply(pi);
    #ifndef NDEBUG
    cout << "FINAL@" << its << " pi=" << pi << " verr=" << verrStart
         <<  " resid=" << normRMSenf << endl;
    #endif
    return converged;
}


//------------------------------------------------------------------------------
//                           SOLVE BILATERAL
//------------------------------------------------------------------------------
bool PGSImpulseSolver::
solveBilateral
   (const Array_<MultiplierIndex>&  participating, // p<=m of these 
    const Matrix&                   A,     // m X m, symmetric
    const Vector&                   D,     // m, diag>=0 added to A
    const Vector&                   rhs,   // m, RHS
    Vector&                         pi     // m, unknown result
    ) const
{
    SimTK_DEBUG("--------------------------------\n");
    SimTK_DEBUG(  "PGS BILATERAL SOLVER:\n");
    ++m_nBilateralSolves;

    const int m=A.nrow(); 
    const int p = (int)participating.size();

    assert(A.ncol()==m); 
    assert(D.size()==0 || D.size()==m);
    assert(rhs.size()==m);
    assert(p<=m);
 
    pi.resize(m);
    pi.setToZero(); // That takes care of all non-participators.

    if (p == 0) {
        SimTK_DEBUG("  no bilateral participators. Nothing to do.\n");
        SimTK_DEBUG("--------------------------------\n");
        return true;
    }


    // Track total error for all included equations, and the error for just
    // those equations that are being enforced.
    bool converged = false;
    Real normRMSenf = Infinity, sor = m_SOR;
    Real prevNormRMSenf = NaN;
    int its = 1;
    Array_<Real> rowSums; // handy temp
    for (; its <= m_maxIters; ++its) {
        ++m_nBilateralIters;
        Real sum2enf = 0; // track solution errors
        prevNormRMSenf = normRMSenf;

        // All the participating constraints are unconditionally active.
        Array_<MultiplierIndex> mults(1);
        for (int k=0; k < p; ++k) {
            mults[0] = participating[k];
            doRowSums(participating,mults,A,D,pi,rowSums);
            const Real localEr2=doUpdates(mults,A,D,rhs,sor,rowSums,pi);
            sum2enf += localEr2;
        }

        normRMSenf = std::sqrt(sum2enf/p);
        const Real rate = normRMSenf/prevNormRMSenf;

        if (rate > 1) {
            SimTK_DEBUG3("GOT WORSE@%d: sor=%g rate=%g\n", its, sor, rate);
            if (sor > .1)
                sor = std::max(.8*sor, .1);
        } 

        #ifndef NDEBUG
        printf("iter %d: EST rmsEnf=%g rate=%g\n", its,
                     normRMSenf, normRMSenf/prevNormRMSenf);
        #endif

        if (normRMSenf < m_convergenceTol) //TODO: add failure-to-improve check
        {
            SimTK_DEBUG2("BILATERAL PGS converged to %g in %d iters\n", 
                         normRMSenf, its);
            converged = true;
            break;
        }
        #ifndef NDEBUG
        cout << "pi=" << pi << " err=" << normRMSenf << " rate=" << rate << endl;
        #endif
    }

    if (!converged) {
        SimTK_DEBUG2("BILATERAL PGS CONVERGENCE FAILURE: %d iters -> norm=%g\n",
              its, normRMSenf);
        ++m_nBilateralFail;
    }

    #ifndef NDEBUG
    cout << "A=" << A;
    cout << "D=" << D << endl;
    cout << "rhs=" << rhs << endl;
    cout << "active=" << participating << endl;
    cout << "-> pi=" << pi << endl;
    if (D.size()) cout << "resid=" << A*pi+D.elementwiseMultiply(pi)-rhs << endl;
    else cout << "resid=" << A*pi-rhs << endl;
    #endif
    SimTK_DEBUG("--------------------------------\n");
    return converged;

}


} // namespace SimTK