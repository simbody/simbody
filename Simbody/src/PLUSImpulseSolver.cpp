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

#include "simbody/internal/common.h"
#include "simbody/internal/ImpulseSolver.h"
#include "simbody/internal/PLUSImpulseSolver.h"

#include <algorithm>
#include <cassert>

#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;


// Local utilities.
namespace {

// Multiply the active entries of a row of the full matrix A by a packed
// column containing only active entries. Useful for A[r]*piActive.
static Real multRowTimesActiveCol(const Matrix& A, MultiplierIndex row, 
           const Array_<MultiplierIndex,PLUSImpulseSolver::ActiveIndex>& active,
           const Vector& colActive) 
{
    const RowVectorView Ar = A[row];
    Real result = 0;
    for (PLUSImpulseSolver::ActiveIndex ax(0); ax < active.size(); ++ax)
        result += Ar[active[ax]] * colActive[ax];
    return result;
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

// Unpack an active column vector and add its values into a full column.
static void addInActiveCol
   (const Array_<MultiplierIndex,PLUSImpulseSolver::ActiveIndex>& active,
    const Vector& colActive,
    Vector& colFull) 
{
    for (PLUSImpulseSolver::ActiveIndex ax(0); ax < active.size(); ++ax) 
        colFull[active[ax]] += colActive[ax];
}

// On return a<=b.
inline void sort2(int& a, int& b) {
    if (a>b) std::swap(a,b);
}
// On return a<=b<=c.
inline void sort3(int& a, int& b, int& c) {
    sort2(a,b); // a<=b
    sort2(b,c); // a<=c, b<=c
    sort2(a,b); // a<=b<=c
}

// Smooth, convex approximation to max(z,0); small eps is smoother.
inline Real softmax0(Real z, Real eps) {
    assert(eps>0);
    return (z+std::sqrt(z*z+eps))/2;
}
// Partial derivative of softmax0 with respect to z.
inline Real dsoftmax0(Real z, Real eps) {
    assert(eps>0);
    return (1+z/std::sqrt(z*z+eps))/2;
}

// Smooth, concave approximation to min(z,0); small eps is smoother.
inline Real softmin0(Real z, Real eps) {
    assert(eps>0);
    return (z-std::sqrt(z*z+eps))/2;
}
// Partial derivative of softmin0 with respect to z.
inline Real dsoftmin0(Real z, Real eps) {
    assert(eps>0);
    return (1-z/std::sqrt(z*z+eps))/2;
}

// Smooth, convex approximation to abs(z); small eps is smoother.
inline Real softabs(Real z, Real eps) {
    assert(eps>0);
    return std::sqrt(z*z+eps);
}
// Partial derivative of softabs with respect to z.
inline Real dsoftabs(Real z, Real eps) {
    assert(eps>0);
    return z/std::sqrt(z*z+eps);
}

}

namespace SimTK {

//==============================================================================
//                   PLUS SUCCESSIVE PRUNING IMPULSE SOLVER
//==============================================================================


//------------------------------------------------------------------------------
//                                SOLVE
//------------------------------------------------------------------------------
bool PLUSImpulseSolver::
solve(int                                 phase,
      const Array_<MultiplierIndex>&      participating, 
      const Matrix&                       A,
      const Vector&                       D,
      const Array_<MultiplierIndex>&      expanding,
      Vector&                             piExpand, // in/out
      Vector&                             verr,     // in/out
      Vector&                             pi, 
      Array_<UncondRT>&                   unconditional,
      Array_<UniContactRT>&               uniContact,
      Array_<UniSpeedRT>&                 uniSpeed,
      Array_<BoundedRT>&                  bounded,
      Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
      Array_<StateLtdFrictionRT>&         stateLtdFriction
      ) const 
{
    SimTK_DEBUG("\n--------------------------------\n");
    SimTK_DEBUG(  "START SUCCESSIVE PRUNING SOLVER:\n");
    ++m_nSolves[phase];

    const int m=A.nrow(); assert(A.ncol()==m); 
    assert(D.size()==m);
    assert(verr.size()==m);
    assert(piExpand.size()==m);

    // These are not mutually exclusive; a contact can be in both lists.
    const int p = (int)participating.size();
    const int nx = (int)expanding.size();
    assert(p<=m); assert(nx<=m);
 
    pi.resize(m);
    pi.setToZero(); // Use this for piUnknown


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


    // This is reduced with each completed sliding interval. We will eventually
    // eliminate all of it except for entries corresponding to friction that
    // remains Sliding throughout the impulse solution.
    m_verrLeft = verr; // what's left to solve TODO: get rid of this
    Vector piELeft = piExpand; // TODO: and this
    m_verrExpand.resize(m); m_verrExpand.setToZero();

    Vector piTotal(m, Real(0)), piGuess(m);
    Vector piSave, dpi; // temps

    // Track total error for all included equations, and the error for just
    // those equations that are being enforced.
    bool converged = false;
    Real normRMSall = Infinity, normRMSenf = Infinity;
    Real prevNormRMSenf = NaN;

    // Each sliding interval requires a complete restart, except that we 
    // continue to accumulate piTotal. We're done when we take an interval of 
    // length frac==1.
    int interval = 0;
    Real frac = 0;
    while (frac < 1) {
        ++interval; 
        m_active = participating; m_mult2active.resize(m);
        fillMult2Active(m_active, m_mult2active);

        // Calculate remaining expansion impulse part of RHS verrE=A*piE.
        // This is how much we'll change verr if we get to apply the full
        // expansion impulse in this sliding interval.
        if (nx)
            for (MultiplierIndex mx(0); mx < m; ++mx) {
                m_verrExpand[mx] = multRowTimesSparseCol(A,mx,expanding,piELeft)
                                   + D[mx]*piELeft[mx];
            }

        if (p == 0) {
            SimTK_DEBUG1("PLUS %d: nothing to do; converged in 0 iters.\n", phase);
            // Returning pi=0; can still have piExpand!=0 so verr is updated.
            if (nx) 
                verr -= m_verrExpand;
            return true;
        }

        #ifndef NDEBUG
        printf("\n***** Sliding interval %d start\n", interval);
        cout << "  active=" << m_active << endl;
        cout << "  mult2active=" << m_mult2active << endl;
        cout << "  piTotal=" << piTotal << endl;
        cout << "  verrLeft=" << m_verrLeft << endl;
        cout << "  expanding=" << expanding << endl;
        cout << "  piELeft=" << piELeft << endl;
        cout << "  verrExpand=" << m_verrExpand << endl;
        #endif

        piGuess = 0; // Hold the best-guess impulse for this interval.

        // Determine step begin Rolling vs. Sliding and get slip directions.
        // Sets all non-Observer uni contacts to active or known.
        classifyFrictionals(uniContact); // no Impendings at interval start

        int its = 1;
        for (; ; ++its) {

            #ifndef NDEBUG
            printf("\n....... Active set iter %d start\n", its); 
            cout << ": active=" << m_active << endl;
            for (unsigned uc=0; uc<uniContact.size(); ++uc) {
                const UniContactRT& rt = uniContact[uc];
                printf("%s UniCont %d (ix=%d): cond=%s/%s, vel=%g,%g, mag=%g\n",
                       getContactTypeName(rt.m_type),(int)uc,(int)rt.m_ucx,
                       getUniCondName(rt.m_contactCond),
                       getFricCondName(rt.m_frictionCond),
                       rt.m_slipVel[0],rt.m_slipVel[1],rt.m_slipMag);
            }
            #endif

            // piGuess has the best guess impulse from the previous active set,
            // unpacked into the associated multiplier slots. This will be
            // the actual piActive values projected to be in-bounds.

            m_mult2active.resize(m);
            fillMult2Active(m_active, m_mult2active);
            initializeNewton(A, piGuess, uniContact);
            updateDirectionsAndCalcCurrentError(A, uniContact, piELeft,
                                                m_piActive,m_errActive);
            
            if (m_active.empty())
                break;
          
            updateJacobianForSliding(A, uniContact, piELeft);
            Real prevNorm = NaN;
            Real errNorm = m_errActive.norm();
            int newtIter = 0;
            SimTK_DEBUG1(">>>> Start NEWTON solve with errNorm=%g...\n", errNorm);
            while (errNorm > m_convergenceTol) {
                ++newtIter;
                // Solve for deltaPi.
                FactorQTZ fac(m_JacActive);
                fac.solve(m_errActive, dpi);
                const Real deltaNorm = dpi.norm();

                #ifndef NDEBUG
                printf("> NEWTON iter %d: errNorm=%g(v) -> deltaNorm=%g(pi)\n", 
                       newtIter, errNorm, dpi.norm());
                cout << "> JacActive=" << m_JacActive;
                cout << "> piActive=" << m_piActive << endl;
                cout << "> errActive=" << m_errActive << endl;
                cout << "> deltaPi=" << dpi << endl;
                #endif

                // Backtracking line search.
                const Real MinFrac = 0.01; // take at least this much
                const Real SearchReduceFac = 0.5;
                
                Real frac = 1;
                int nsearch = 0;
                piSave = m_piActive;
                prevNorm = errNorm;
                while (true) {
                    ++nsearch;
                    SimTK_DEBUG3("Line search iter %d: frac=%g, prevNorm=%g.\n", 
                                 nsearch, frac, prevNorm);
                    m_piActive = piSave - frac*dpi;

                    updateDirectionsAndCalcCurrentError(A,uniContact,piELeft,
                                                        m_piActive,
                                                        m_errActive);
                    errNorm = m_errActive.norm();
                    #ifndef NDEBUG
                    cout << "> piNow=" << m_piActive << endl;
                    cout << "> errNow=" << m_errActive
                         << " normNow=" << errNorm << endl;
                    #endif
                    if (errNorm < prevNorm)
                        break;

                    frac *= SearchReduceFac;
                    if (frac*SearchReduceFac < MinFrac) {
                        SimTK_DEBUG3("LINE SEARCH STUCK at iter %d: accepting "
                            "small norm increase %g at frac=%g\n", nsearch,
                            errNorm - prevNorm, frac);
                        break;
                    }
                    SimTK_DEBUG2("GOT WORSE @iter %d: backtrack to frac=%g\n", 
                           nsearch, frac);
                }

                SimTK_DEBUG2("Improvement rate now/prev at iter %d is %g\n",
                                newtIter, errNorm / prevNorm);

                if (errNorm < m_convergenceTol)
                    break; // we have a winner

                // If we're not making sufficient progress after 3 steps,
                // we have likely converged to a local minimum and the problem 
                // is infeasible. Here we're arbitrarily saying that if we
                // can't even improve 5%, we'll just give up.
                if (newtIter >= 3 && errNorm > 0.95*prevNorm) {
                    SimTK_DEBUG2("PLUSImpulseSolver Newton: poor progress "
                    "after %d iters; errNorm=%g (infeasible?).\n", 
                    newtIter, errNorm);
                    break; // converged, but not to zero
                }

                if (newtIter >= m_maxIters) {
                    SimTK_DEBUG2("PLUSImpulseSolver Newton failed to converge "
                        "after %d iters; errNorm=%g.\n", m_maxIters, errNorm);
                    break; // we have a loser
                }

                updateJacobianForSliding(A, uniContact, piELeft);
                prevNorm = errNorm;
            }

            SimTK_DEBUG2("<<<< NEWTON done in %d iters; norm=%g.\n",
                         newtIter,errNorm);

            // UNCONDITIONAL: these are always on.
            for (int fx=0; fx < mUncond; ++fx) {
                const UncondRT& rt = unconditional[fx];
                for (unsigned i=0; i<rt.m_mults.size(); ++i) {
                    const MultiplierIndex mx = rt.m_mults[i];
                    piGuess[mx] = m_piActive[m_mult2active[mx]]; // unpack
                }
            }

            // BOUNDED: conditional scalar constraints with constant bounds
            // on resulting pi.
            int worstBounded=0; Real worstBoundedValue=0;
            for (int k=0; k < mBounded; ++k) {
                const BoundedRT&      rt = bounded[k];
                const MultiplierIndex mx = rt.m_ix;
                const ActiveIndex     ax = m_mult2active[mx];
                if (!ax.isValid())
                    continue; // not active
                // Only the in-bounds value gets saved in piGuess in case we
                // need to use it for an initial guess on the next iteration.
                piGuess[mx] = clamp(rt.m_lb, m_piActive[ax], rt.m_ub);
                const Real err=std::abs(m_piActive[ax] - piGuess[mx]);
                if (err>worstBoundedValue) 
                    worstBounded=k, worstBoundedValue=err;
            }

            // UNI CONTACT NORMAL: conditional scalar constraints with 
            // with restriction pi <= 0.
            int worstUniNormal=0; Real worstUniNormalValue=0;
            for (int k=0; k < mUniCont; ++k) {
                const UniContactRT&   rt = uniContact[k];
                const MultiplierIndex mx = rt.m_Nk;
                if (rt.m_contactCond==UniOff || rt.m_contactCond==UniKnown) {
                    piGuess[mx] = 0;
                    continue;
                }

                // Participating and active.
                assert(rt.m_contactCond == UniActive);
                const ActiveIndex ax = m_mult2active[mx];
                assert(ax.isValid());
                // Only the in-bounds value gets saved in piGuess in case we
                // need to use it for an initial guess on the next iteration.
                const Real piAdj = rt.m_sign*m_piActive[ax] < 0 ? m_piActive[ax] 
                                                                : Real(0);
                piGuess[mx] = piAdj; 
                const Real err=std::abs(m_piActive[ax] - piAdj);
                if (err>worstUniNormalValue) 
                    worstUniNormal=k, worstUniNormalValue=err;
            }


            // UNI CONTACT FRICTION: a set of constraint equations forming a 
            // vector whose maximum length is limited by the associated 
            // unilateral contact normal force.
            int worstFric=0; Real worstFricValue=0;
            for (int k=0; k < mUniCont; ++k) {
                const UniContactRT& rt = uniContact[k];
                if (rt.m_contactCond==UniOff || !rt.hasFriction())
                    continue;
                // Known, or Participating and active, and has friction.
                const Array_<MultiplierIndex>& Fk = rt.m_Fk; // friction components
                const MultiplierIndex Nk = rt.m_Nk; // normal component
                assert(m_mult2active[Fk[0]].isValid());
                const Real mu = rt.m_effMu;
                Real scale = 1; // might change if we're rolling

                // Only if rolling is there an inequality constraint that
                // must be satisfied; calculate its violation here.
                if (rt.m_frictionCond == Rolling) {
                    Real tmag=0, nmag=0;
                    for (unsigned i=0; i<Fk.size(); ++i) {
                        const MultiplierIndex mx = Fk[i];
                        const ActiveIndex ax = m_mult2active[mx];
                        tmag += square(m_piActive[ax]);
                    }
                    tmag = std::sqrt(tmag);

                    // "Sucking" normal forces are zero already in piGuess,
                    // and known normal force has been inserted if needed.
                    nmag = std::abs(piGuess[Nk] + piELeft[Nk]); 
                    if (tmag > mu*nmag) {
                        scale = mu*nmag/tmag;
                        const Real err = tmag - mu*nmag;
                        if (err > worstFricValue)
                            worstFric=k, worstFricValue=err;
                    }
                }

                // Copy the possibly-reduced value into piGuess.
                for (unsigned i=0; i<Fk.size(); ++i) {
                    const MultiplierIndex mx = Fk[i];
                    const ActiveIndex ax = m_mult2active[mx];
                    piGuess[mx] = scale*m_piActive[ax];
                }
            }

            // TODO: uni speed, constraint- and state-limited friction.

            if (   worstBoundedValue<=SignificantReal
                && worstUniNormalValue<=SignificantReal
                && worstFricValue<=SignificantReal) 
            {
                SimTK_DEBUG3("Bounded/Contact/Rolling OK: worst=%g/%g/%g. "
                    "Check sliding next.\n",
                    worstBoundedValue, worstUniNormalValue, worstFricValue);
                break;
            }

            //TODO: bounded

            bool mustReleaseFriction = true; // if we don't release a normal.
            if (worstUniNormalValue > worstFricValue) {
                SimTK_DEBUG2("Worst offender is normal contact %d err=%g ...\n", 
                    worstUniNormal, worstUniNormalValue);
                // A contact normal is the worst offender. However, if it has a
                // rolling friction constraint active we should release that 
                // first because doing so might fix the contact normal.
                UniContactRT& rt = uniContact[worstUniNormal];
                if (!rt.hasFriction() || rt.m_frictionCond != Rolling) {
                    const MultiplierIndex rx = rt.m_Nk;
                    rt.m_contactCond = UniOff;
                    // Update active set; must work from highest numbered to 
                    // lowest to avoid having to move a lot of entries.
                    if (!rt.hasFriction()) {
                        m_active.eraseFast(m_active.begin()+m_mult2active[rx]);
                    } else {
                        const Array_<MultiplierIndex>& Fk = rt.m_Fk;
                        int a=m_mult2active[rx],b=m_mult2active[Fk[0]],
                            c=m_mult2active[Fk[1]];
                        sort3(a,b,c);
                        m_active.eraseFast(m_active.begin()+c);
                        m_active.eraseFast(m_active.begin()+b);
                        m_active.eraseFast(m_active.begin()+a);
                    }
                    // mult2active is invalid now.
                    mustReleaseFriction = false;
                    SimTK_DEBUG1("... normal contact %d released.\n", 
                                 worstUniNormal);
                } else {
                    SimTK_DEBUG("... but it's Rolling, so that must go first.\n");
                    worstFric = worstUniNormal;
                    worstFricValue = NaN;
                    mustReleaseFriction = true;
                }
            }

            if (mustReleaseFriction) {
                UniContactRT& rt = uniContact[worstFric];
                const Array_<MultiplierIndex>& Fk = rt.m_Fk;
                const MultiplierIndex          Nk = rt.m_Nk;
                const ActiveIndex ax=m_mult2active[Fk[0]], ay=m_mult2active[Fk[1]], 
                                  az=m_mult2active[Nk];

                SimTK_DEBUG2("switch worst fric %d from roll->impend err=%g\n", 
                       worstFric, worstFricValue);
                rt.m_frictionCond = Impending;

                #ifndef NDEBUG
                // Oppose the last rolling force as a guess at the slip velocity.
                // Sign convention for multiplier is opposite velocity, so no 
                // explicit negation here.
                const Vec2 ft(piGuess[Fk[0]], piGuess[Fk[1]]);
                cout << "  rolling impulse was " << ft << endl;
                #endif
            }
        } 

        // Need to check what fraction of this interval we can accept. We are
        // only limited by frictional contacts that are currently Sliding;
        // Rolling and Impending-slip contacts don't restrict the interval.
        frac = 1;
        for (int k=0; k < mUniCont; ++k) {
            const UniContactRT& rt = uniContact[k];
            if (rt.m_contactCond==UniOff || rt.m_frictionCond != Sliding)
                continue;
            const Array_<MultiplierIndex>& Fk = rt.m_Fk;
            const MultiplierIndex          Nk = rt.m_Nk;
            assert(Fk.size()==2); //TODO: generalize
            // New velocity db=[Ax Ay]*(pi+piE). TODO: D?
            Vec2 db(  multRowTimesActiveCol(A,Fk[0],m_active,m_piActive)
                    + m_verrExpand[Fk[0]],
                      multRowTimesActiveCol(A,Fk[1],m_active,m_piActive)
                    + m_verrExpand[Fk[1]]);
            Vec2 bend = rt.m_slipVel - db;
            #ifndef NDEBUG
            cout << "slipVel " << k << " from " << rt.m_slipVel 
                 << " to " << bend << endl;
            #endif
            const Real bendMag = bend.norm();
            SimTK_ASSERT2_ALWAYS(rt.m_slipMag > m_maxRollingTangVel,
                "PLUSImpulseSolver::solve(): contact %d misclassified as "
                "Sliding; slip speed %g too small (Rolling at %g).",
                rt.m_slipMag, m_maxRollingTangVel);

            if (bendMag <= m_maxRollingTangVel) {
                SimTK_DEBUG2("Friction %d slowed to a halt, v=%g\n", k, bendMag);
                continue;
            }
            const Real cosTheta = 
                clamp(-1, dot(rt.m_slipVel,bend)/(rt.m_slipMag*bendMag), 1);
            if (cosTheta >= m_cosMaxSlidingDirChange) {
                SimTK_DEBUG3("Friction %d rotated %g deg, less than max %g\n", 
                       k, std::acos(cosTheta)*180/Pi,
                       std::acos(m_cosMaxSlidingDirChange)*180/Pi);
                continue;
            }
            SimTK_DEBUG4("TOO BIG: Sliding fric %d; endmag=%g, rot=%g deg > %g\n", 
                   k, bendMag, std::acos(cosTheta)*180/Pi,
                   std::acos(m_cosMaxSlidingDirChange)*180/Pi);

            Vec2 endPt;
            Real frac1 = calcSlidingStepLengthToOrigin(rt.m_slipVel,bend,endPt);
            const Real endPtMag = endPt.norm();
            if (endPtMag <= m_maxRollingTangVel) {
                SimTK_DEBUG2("  Frac=%g halts it, v=%g\n", frac1, endPtMag);
                frac = std::min(frac, frac1);
                continue;
            }
            Real frac2 = calcSlidingStepLengthToMaxChange(rt.m_slipVel,bend);
            SimTK_DEBUG2("  Frac=%g reduces angle to %g degrees.\n", 
                   frac2, std::acos(m_cosMaxSlidingDirChange)*180/Pi);
            frac = std::min(frac, frac2);
        }

        for (unsigned i=0; i < expanding.size(); ++i) {
            const MultiplierIndex mx = expanding[i];
            const Real alphaPiE = frac*piELeft[mx];
            piELeft[mx] -= alphaPiE; // How much piE left to do
        }
        m_piActive *= frac;
        addInActiveCol(m_active, m_piActive, piTotal); // accumulate in piTotal

        // Update rhs. TODO: D*piActive
        for (MultiplierIndex mx(0); mx < m; ++mx) {
            m_verrLeft[mx] -= multRowTimesActiveCol(A,mx,m_active,m_piActive)
                              + frac*m_verrExpand[mx];
        }

        #ifndef NDEBUG
        printf("SP interval %d end: frac=%g\n", interval, frac);
        cout << ": m_piActive=" << m_piActive << endl;
        cout << ": m_verrLeft=" << m_verrLeft << endl;
        cout << ": piELeft=" << piELeft << endl;
        #endif
    }

    // Return the result. TODO: don't copy 
    pi = piTotal; // doesn't include piE
    verr = m_verrLeft;

    // Check how we did on the original problem.
    SimTK_DEBUG("SP DONE. Check normal complementarity ...\n");
    for (unsigned k=0; k < uniContact.size(); ++k) {
        const UniContactRT& rt = uniContact[k];
        const MultiplierIndex mx = rt.m_Nk;
        SimTK_DEBUG4("%d: pi=%g verr=%g pi*v=%g\n", k, 
                     pi[mx], verr[mx], pi[mx]*verr[mx]);
    } 
    //TODO: printf("SP DONE. Check friction cones ...\n");

    #ifndef NDEBUG
    cout << "SP FINAL " << interval << " intervals, piTotal=" << piTotal 
         <<  " errNorm=" << m_errActive.norm() << endl;
    #endif
    return converged;
}

//------------------------------------------------------------------------------
//                           SOLVE BILATERAL
//------------------------------------------------------------------------------
bool PLUSImpulseSolver::
solveBilateral
   (const Array_<MultiplierIndex>&  participating, // p<=m of these 
    const Matrix&                   A,     // m X m, symmetric
    const Vector&                   D,     // m, diag>=0 added to A
    const Vector&                   rhs,   // m, RHS
    Vector&                         pi     // m, unknown result
    ) const
{
    SimTK_DEBUG("--------------------------------\n");
    SimTK_DEBUG(  "PLUS BILATERAL SOLVER:\n");
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

    m_nBilateralIters++; // Just one "iteration".

    // Set up the smaller problem containing only the participating constraints:
    //   bilateralActive = P*A*~P
    //   rhsActive       = P*rhs
    //   piActive = a place to put the result for just the active part
    m_bilateralActive.resize(p,p);
    m_rhsActive.resize(p); m_piActive.resize(p);
    const bool hasD = (D.size() > 0);
    for (ActiveIndex aj(0); aj < p; ++aj) {
        const MultiplierIndex mj = participating[aj];
        for (ActiveIndex ai(0); ai < p; ++ai) {
            const MultiplierIndex mi = participating[ai];
            m_bilateralActive(ai,aj) = A(mi,mj);
        }
        if (hasD)
            m_bilateralActive(aj,aj) += D[mj];
        m_rhsActive[aj] = rhs[mj];
    }

    // Calculate the pseudoinverse of P*A*~P, and then solve to get
    //     piActive = pinv(P*A*~P) * rhsActive
    FactorQTZ pinv(m_bilateralActive);
    pinv.solve(m_rhsActive, m_piActive);
    // Distribute the active result into the full impulse vector.
    for (ActiveIndex ai(0); ai < p; ++ai) {
        const MultiplierIndex mi = participating[ai];
        pi[mi] = m_piActive[ai];
    }

    #ifndef NDEBUG
    cout << "A=" << A;
    cout << "D=" << D << endl;
    cout << "rcond(A+D)=" << pinv.getRCondEstimate() 
         << " rank=" << pinv.getRank() << endl;
    cout << "rhs=" << rhs << endl;
    cout << "active=" << participating << endl;
    cout << "-> piActive=" << m_piActive << endl;
    cout << "-> pi=" << pi << endl;
    cout << "resid active=" << m_bilateralActive*m_piActive-m_rhsActive << endl;
    if (D.size()) cout << "resid=" << A*pi+D.elementwiseMultiply(pi)-rhs << endl;
    else cout << "resid=" << A*pi-rhs << endl;
    #endif
    SimTK_DEBUG("--------------------------------\n");
    return true;
}

//------------------------------------------------------------------------------
//              CALC SLIDING STEP LENGTH TO ORIGIN (Vec2 and Vec3)
//------------------------------------------------------------------------------
Real PLUSImpulseSolver::
calcSlidingStepLengthToOrigin(const Vec2& A, const Vec2& B, Vec2& Q) const
{
    // Check whether initial tangential velocity is small (impending slip).
    if (A.normSqr() < square(m_maxRollingTangVel)) {
        SimTK_DEBUG2("--> initial slip velocity small (%g<%g); stepLen=1\n",
                     A.norm(), m_maxRollingTangVel);
        Q = B;
        return 1;
    }

    const Vec2 P     = Vec2(0);
    const Vec2 AtoP  = P-A, AtoB  = B-A;
    const Real ABsqr = AtoB.normSqr();

    // Ensure line segment is of meaningful length.
    if (ABsqr < SimTK::SignificantReal) {
        SimTK_DEBUG1("-->ABsqr=%g short; returning stepLength=1\n", ABsqr);
        Q = B;
        return 1;
    }

    // Normalized distance from A to Q.
    const Real stepLength = clamp(0.0, dot(AtoP,AtoB)/ABsqr, 1.0);
    Q = A + stepLength*AtoB;

    SimTK_DEBUG2("--> returning stepLength=%g (dist to origin=%g)\n",
                    stepLength, Q.norm());

    return stepLength;
}

Real PLUSImpulseSolver::
calcSlidingStepLengthToOrigin(const Vec3& A, const Vec3& B, Vec3& Q) const
{
    // Check whether initial tangential velocity is small (impending slip).
    if (A.normSqr() < square(m_maxRollingTangVel)) {
        SimTK_DEBUG2("--> initial slip velocity small (%g<%g); stepLen=1\n",
                     A.norm(), m_maxRollingTangVel);
        Q = B;
        return 1;
    }

    const Vec3 P     = Vec3(0);
    const Vec3 AtoP  = P-A, AtoB  = B-A;
    const Real ABsqr = AtoB.normSqr();

    // Ensure line segment is of meaningful length.
    if (ABsqr < SimTK::SignificantReal) {
        SimTK_DEBUG1("-->ABsqr=%g short; returning stepLength=1\n", ABsqr);
        Q = B;
        return 1;
    }

    // Normalized distance from A to Q.
    const Real stepLength = clamp(0.0, dot(AtoP,AtoB)/ABsqr, 1.0);
    Q = A + stepLength*AtoB;

    SimTK_DEBUG2("--> returning stepLength=%g (dist to origin=%g)\n",
                    stepLength, Q.norm());

    return stepLength;
}


//------------------------------------------------------------------------------
//            CALC SLIDING STEP LENGTH TO MAX CHANGE (Vec2 and Vec3)
//------------------------------------------------------------------------------
Real PLUSImpulseSolver::
calcSlidingStepLengthToMaxChange(const Vec2& A, const Vec2& B) const
{
    // Temporary variables created by dsolve/numeric/optimize.
    Real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, sol1, sol2;
    const Vec2 v = B-A;

    // Optimized computation sequence generated in Maple.
    t1 = m_cosMaxSlidingDirChange;
    t1 *= t1;
    t2 = t1 - 1;
    t3 = A[0]*v[1] - A[1]*v[0];
    t3 = std::sqrt(-t1*t2*t3*t3);
    t4 = t2*v[0]*A[0];
    t5 = A[1]*v[1];
    t2 *= t5;
    t6 = v[1]*v[1];
    t7 = v[0]*v[0];
    t8 = t6 + t7;
    t9 = A[1]*A[1];
    t10 = A[0]*A[0];
    t1 = t1*(t10*t8 + t8*t9) - t10*t7 - t6*t9 - 2*t5*A[0]*v[0];
    t5 = t10 + t9;
    t1 = 1 / t1;

    sol1 = -t1*t5*(t2 + t4 + t3);
    sol2 = -t1*t5*(t2 + t4 - t3);
    assert(sol1>=0 || sol2>=0); //TODO: is this guaranteed?
    Real sol;
    if (sol1 < 0) sol=sol2;
    else if (sol2 < 0) sol=sol1;
    else sol = std::min(sol1, sol2);

    SimTK_DEBUG3("-->max change solutions: %g and %g; returning %g\n",
                 sol1,sol2,sol);

    return sol;
}

Real PLUSImpulseSolver::
calcSlidingStepLengthToMaxChange(const Vec3& A, const Vec3& B) const
{
    // Temporary variables created by dsolve/numeric/optimize.
    Real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    Real sol1, sol2;
    const Vec3 v = B-A;

    // Optimized computation sequence generated in Maple.
    t1 = m_cosMaxSlidingDirChange;
    t1 *= t1;
    t2 = t1 - 1;
    t3 = A[0] * A[0];
    t4 = v[0] * v[0];
    t5 = A[2] * A[2];
    t6 = v[1] * v[1];
    t7 = A[1] * A[1];
    t8 = A[1] * v[1];
    t9 = A[0] * v[0];
    t10 = std::sqrt(-(t1 * t2 * (t3 * t6 + t4 * t7 + t5 * (t6 + t4) 
          + (-2 * A[2] * (t9 + t8) + (t7 + t3) * v[2]) * v[2] - 2 * t8 * t9)));
    t11 = t9 * t2;
    t12 = t8 * t2;
    t13 = A[2] * v[2];
    t2 = t13 * t2;
    t14 = v[2] * v[2];
    t15 = t6 + t14 + t4;
    t1 = t1 * (t15 * t3 + t15 * t5 + t15 * t7) - t14 * t5 - t3 * t4 - t6 * t7
         + t9 * (-2 * t8 - 2 * t13) - 2 * t13 * t8;
    t3 = t7 + t3 + t5;
    t1 = 1 / t1;

    sol1 = -(t12 + t2 + t11 + t10) * t1 * t3;
    sol2 = -(t12 + t2 + t11 - t10) * t1 * t3;

    Real sol;
    if (sol1 < 0) sol=sol2;
    else if (sol2 < 0) sol=sol1;
    else sol = std::min(sol1, sol2);

    SimTK_DEBUG3("-->max change solutions: %g and %g; returning %g\n",
                 sol1,sol2,sol);

    return sol;
}

//------------------------------------------------------------------------------
//                           CLASSIFY FRICTIONALS
//------------------------------------------------------------------------------
// At the start of each sliding interval, classify all frictional contacts.
// For unilateral contact friction, if the unilateral normal contact is 
// Observing (passive) then its friction constraints are off also. Otherwise
// (normal is Participating or Known), every frictional contact is classified 
// as Rolling or Sliding depending on the current slip velocity as present
// in the remaining right hand side of the rolling equations in A. No frictional
// contact is marked Impending at the start of a sliding interval; that only
// occurs as a result of a transition from Rolling.
void PLUSImpulseSolver::
classifyFrictionals(Array_<UniContactRT>& uniContact) const {
    SimTK_DEBUG1("classifyFrictionals(): %d uni contacts\n", 
                 (int)uniContact.size());
    for (unsigned k=0; k < uniContact.size(); ++k) {
        UniContactRT& rt = uniContact[k];

        // Set contact condition.
        if (rt.m_type==Participating) rt.m_contactCond = UniActive;
        else if (rt.m_type==Known)  rt.m_contactCond = UniKnown;
        else {assert(rt.m_type==Observing); rt.m_contactCond = UniOff;}

        // Set friction condition and slip velocity.
        if (rt.m_type == Observing || !rt.hasFriction()) {
            rt.m_frictionCond = FricOff;
            rt.m_slipVel = Vec2(NaN); // for bug catching
            rt.m_slipMag = NaN;
        } else { // normal is Participating or Known and has friction.
            const Array_<MultiplierIndex>& Fk = rt.m_Fk; // friction components
            assert(Fk.size()==2); //TODO: generalize
            Real tmag=0;
            for (unsigned i=0; i<Fk.size(); ++i) {
                const MultiplierIndex mx = Fk[i];
                rt.m_slipVel[i] = m_verrLeft[mx];
                tmag += square(m_verrLeft[mx]);
            }
            tmag = std::sqrt(tmag);
            rt.m_slipMag = tmag;
            rt.m_frictionCond = tmag > m_maxRollingTangVel ? Sliding : Rolling;
        }

        #ifndef NDEBUG
        printf("  %s contact %d is %s; vel=%g,%g, mag=%g\n",
               getContactTypeName(rt.m_type), (int)k,
               getFricCondName(rt.m_frictionCond), rt.m_slipVel[0],
               rt.m_slipVel[1], rt.m_slipMag);
        #endif
    }
}

//------------------------------------------------------------------------------
//                            FILL MULT 2 ACTIVE
//------------------------------------------------------------------------------
// mult2active must already have been resized to size of A
void PLUSImpulseSolver::
fillMult2Active(const Array_<MultiplierIndex,ActiveIndex>& active,
                Array_<ActiveIndex,MultiplierIndex>& mult2active) const
{
    const int p = active.size();
    mult2active.fill(ActiveIndex()); // invalid
    for (ActiveIndex aj(0); aj < p; ++aj) {
        const MultiplierIndex mj = active[aj];
        mult2active[mj] = aj;
    }
    #ifndef NDEBUG
    printf("fillMult2Active:\n");
    cout << ": active=" << active << endl;
    cout << ": mult2active=" << mult2active << endl;
    #endif
}

//------------------------------------------------------------------------------
//                             INITIALIZE NEWTON
//------------------------------------------------------------------------------
// Initialize for a Newton iteration. Fill in the part of the Jacobian
// corresponding to linear equations since those won't change. Transfer
// previous impulses pi to new piActive. Assumes m_active and m_mult2active
// have been filled in.
void PLUSImpulseSolver::
initializeNewton(const Matrix& A, const Vector& pi, // m of these 
                 const Array_<UniContactRT>& uniContact) const { 
    const int na = m_active.size();
    m_JacActive.resize(na,na); m_rhsActive.resize(na); m_piActive.resize(na);
    m_errActive.resize(na);
    for (ActiveIndex aj(0); aj < na; ++aj) {
        const MultiplierIndex mj = m_active[aj];
        for (ActiveIndex ai(0); ai < na; ++ai) {
            const MultiplierIndex mi = m_active[ai];
            m_JacActive(ai,aj) = A(mi,mj);
        }
        m_rhsActive[aj] = m_verrLeft[mj] - m_verrExpand[mj];
        m_piActive[aj]  = pi[mj];
    }
    // For impacters, guess a small separating impulse. This improves
    // convergence because it puts the max() terms in the Jacobian on
    // the right branch.
    for (unsigned k=0; k < uniContact.size(); ++k) {
        const UniContactRT& rt = uniContact[k];
        if (rt.m_contactCond != UniActive)
            continue;

        const MultiplierIndex mz = rt.m_Nk;
        const ActiveIndex az = m_mult2active[mz];
        assert(az.isValid());
        // Don't use rt.m_sign here because it affects both pi & rhs; we just
        // want the signs to match.
        m_piActive[az] = .01*sign(m_rhsActive[az]); //-1,0,1
        SimTK_DEBUG3("  active normal %d has v=%g; guess pi=%g\n",
                    (int)az,m_rhsActive[az],m_piActive[az]);
    }
    #ifndef NDEBUG
    printf("initializeNewton:\n");
    cout << ": verrLeft was=" << m_verrLeft << endl;
    cout << ": verrExpand was=" << m_verrExpand << endl;
    cout << ": rhsActive=" << m_rhsActive << endl;
    cout << ": pi was=" << pi << endl;
    cout << ": piActive=" << m_piActive << endl;
    #endif
}

//------------------------------------------------------------------------------
//                  UPDATE DIRECTIONS AND CALC CURRENT ERROR
//------------------------------------------------------------------------------
// Calculate err(pi). For Impending slip frictional contacts we also revise
// the slip direction based on the current values of pi and piExpand.
void PLUSImpulseSolver::
updateDirectionsAndCalcCurrentError
   (const Matrix& A,  Array_<UniContactRT>& uniContact, 
    const Vector& piELeft, const Vector& piActive,
    Vector& errActive) const 
{
    const int na = m_active.size();
    assert(piActive.size() == na);
    errActive.resize(na);
    // Initialize as though all rolling.
    for (ActiveIndex ai(0); ai < na; ++ai) {
        const MultiplierIndex mi = m_active[ai];
        // err = A pi - rhs (piExpand included in rhs)
        errActive[ai] = multRowTimesActiveCol(A,mi,m_active,piActive)
                        - m_rhsActive[ai];
    }

    // Replace error equations for sliding and impending slip. For impending
    // slip we'll also update slipVel and slipMag since we'll need them again
    // when we calculate the Jacobian.
    for (unsigned k=0; k < uniContact.size(); ++k) {
        UniContactRT& rt = uniContact[k];
        if (rt.m_contactCond == UniOff || !rt.hasFriction()) 
            continue; // inactive, or no friction
        if (!(rt.m_frictionCond==Sliding || rt.m_frictionCond==Impending))
            continue; // no need to modify the equations

        const Array_<MultiplierIndex>& Fk = rt.m_Fk;
        const MultiplierIndex          Nk = rt.m_Nk;
        assert(Fk.size()==2); //TODO: generalize
        const MultiplierIndex mx=Fk[0], my=Fk[1], mz=Nk;

        if (rt.m_frictionCond==Impending) {
            // Update slip direction to [Ax Ay]*(pi+piE).
            Vec2 d(multRowTimesActiveCol(A,mx,m_active,piActive)
                   + m_verrExpand[mx],
                   multRowTimesActiveCol(A,my,m_active,piActive)
                   + m_verrExpand[my]);
            const Real dnorm = d.norm();
            rt.m_slipVel = d; rt.m_slipMag = dnorm;
            SimTK_DEBUG3("Updated impending slipVel %d to %g,%g\n",k,d[0],d[1]);
        }

        const Real mu = rt.m_effMu;
        const ActiveIndex ax=m_mult2active[mx], ay=m_mult2active[my], 
                          az=m_mult2active[mz];
        const Real pix = piActive[ax], piy=piActive[ay];

        // Applying the sign here makes sure pizE is negative.
        const Real pizE = rt.m_sign*piELeft[mz];

        // The slipping or impending slip equation is
        //     f = - mu N v/|v|
        // where v is the slip velocity or impending slip velocity (above), and
        // f is the tangential friction force vector, and N>=0 is the magnitude
        // of the normal force. Our convention for multipliers is opposite
        // applied forces, so we have
        //    -pi_xy = - mu (-piz) v/|v|
        // ==> pi_xy = - mu piz v/|v|
        // ==> |v| pi_xy + mu piz v = 0
        // We write the two error functions like this:
        //     errx=|v|pi_x + mu*vx*[piE+min(pi_z,0)] 
        //     erry=|v|pi_y + mu*vy*[piE+min(pi_z,0)] 

        // Calculate the terms common to Active and Known contacts.
        errActive[ax] = rt.m_slipMag*pix + mu*rt.m_slipVel[0]*pizE;
        errActive[ay] = rt.m_slipMag*piy + mu*rt.m_slipVel[1]*pizE;

        // Add in the additional term for Active contacts.
        if (rt.m_contactCond == UniActive) {
            const ActiveIndex az=m_mult2active[mz];
            assert(az.isValid());
            // Applying the sign here makes piz negative if it represents a 
            // valid "pushing" force.
            const Real piz=rt.m_sign*piActive[az];
            const Real minz = std::min(piz, Real(0));
            errActive[ax] += mu*rt.m_slipVel[0]*minz;
            errActive[ay] += mu*rt.m_slipVel[1]*minz;
        }
    }
}


//------------------------------------------------------------------------------
//                       UDPATE JACOBIAN FOR SLIDING
//------------------------------------------------------------------------------
// Calculate Jacobian 
//      J= D err(pi) / D pi
// See updateDirectionsAndCalcCurrentError() for err(pi). All rows of J 
// corresponding to linear equations, including rolling constraints, have 
// already been filled in since they can't change during the iteration. Only 
// sliding and impending friction rows are potentially nonlinear and thus
// subject to change during the Newton iterations.
void PLUSImpulseSolver::
updateJacobianForSliding(const Matrix& A,
                         const Array_<UniContactRT>& uniContact,
                         const Vector& piELeft) const {
    int nPairsChanged = 0;
    for (unsigned k=0; k < uniContact.size(); ++k) {
        const UniContactRT& rt = uniContact[k];
        if (!(rt.m_contactCond==UniActive||rt.m_contactCond==UniKnown)
            || !rt.hasFriction())
            continue;
        // Known, or Participating and active, and has friction.
        if (!(rt.m_frictionCond==Sliding || rt.m_frictionCond==Impending))
            continue;

        const Array_<MultiplierIndex>& Fk = rt.m_Fk;
        assert(Fk.size()==2); //TODO: generalize
        const MultiplierIndex mx=Fk[0], my=Fk[1];
        assert(m_mult2active[mx].isValid()); 
        assert(m_mult2active[my].isValid()); 

         // Handy abbreviations to better match equations.
        const Real mu = rt.m_effMu;
        const ActiveIndex ax=m_mult2active[mx], ay=m_mult2active[my];
        const Real pix = m_piActive[ax], piy=m_piActive[ay];
        const Vec2 d     = rt.m_slipVel;
        const Real dnorm = rt.m_slipMag;
        const Vec2 dhat = dnorm > TinyReal ? d/dnorm : Vec2(0);

        m_JacActive[ax] = m_JacActive[ay] = 0; // zero the rows
        if (rt.m_frictionCond==Impending) {
            // Calculate terms for derivative of norm(d) w.r.t. pi.
            const RowVectorView Ax = A[mx], Ay = A[my];
            const MultiplierIndex mz = rt.m_Nk;
            const Real pizE = rt.m_sign*piELeft[mz];

            if (rt.m_contactCond==UniActive) { // Impending normal is active
                const ActiveIndex az=m_mult2active[mz];
                assert(az.isValid());
                const Real piz=rt.m_sign*m_piActive[az];
                const Real Axz=Ax(mz), Ayz=Ay(mz);
                const Real minz  = softmin0(piz, m_minSmoothness);
                const Real dminz = dsoftmin0(piz, m_minSmoothness);
                // errx=|d|pix + dx*mu*(pizE+softmin0(piz))   [erry similar]
                // d/dpii errx = s*pix + mu*Axi*(pizE+softmin0(piz)), i!=x,z
                // d/dpix errx = s*pix + mu*Axx*(pizE+softmin0(piz)) + |d|
                // d/dpiz errx = s*pix + mu*Axz*(pizE+softmin0(piz))
                //                                       + mu*dx*dsoftmin0(piz)
                // where s = dot(d/|d|, [Axi Ayi]).

                // Fill in generic terms for unrelated constraints (not x,y,z)
                for (ActiveIndex ai(0); ai<m_active.size(); ++ai) {
                    const MultiplierIndex mi = m_active[ai];
                    const Real pii=m_piActive[ai];
                    const Real Axi=Ax(mi), Ayi=Ay(mi);
                    const Real s = ~dhat*Vec2(Axi,Ayi);
                    m_JacActive(ax,ai) = s*pix + mu*Axi*(pizE+minz);
                    m_JacActive(ay,ai) = s*piy + mu*Ayi*(pizE+minz);
                }
                // Add additional terms for related rows.
                m_JacActive(ax,ax) += dnorm;            // d errx / dx
                m_JacActive(ay,ay) += dnorm;            // d erry / dy
                m_JacActive(ax,az) += mu*d[0]*dminz;    // d errx / dz
                m_JacActive(ay,az) += mu*d[1]*dminz;    // d erry / dz

            } else { // Impending normal is an expander; piz is not variable
                assert(rt.m_contactCond==UniKnown);
                // errx=|d|pix + dx*mu*pizE   [erry similar]
                // d/dpii errx = s*pix + mu*Axi*pizE, for i != x
                // d/dpix errx = s*pix + mu*Axx*pizE + |d|
                // where s = dot(d/|d|, [Axi Ayi]).

                // Fill in generic terms for unrelated constraints (not x,y)
                for (ActiveIndex ai(0); ai<m_active.size(); ++ai) {
                    const MultiplierIndex mi = m_active[ai];
                    const Real pii=m_piActive[ai];
                    const Real Axi=Ax(mi), Ayi=Ay(mi);
                    const Real s = ~dhat*Vec2(Axi,Ayi);
                    m_JacActive(ax,ai) = s*pix + mu*Axi*pizE;
                    m_JacActive(ay,ai) = s*piy + mu*Ayi*pizE;
                }
                m_JacActive(ax,ax) += dnorm;
                m_JacActive(ay,ay) += dnorm;
            }
        } else { // Slipping
            m_JacActive(ax,ax) = m_JacActive(ay,ay) = dnorm;
            // That's all for an expander; active also has z derivs.
            if (rt.m_contactCond==UniActive) { // normal is active
                const ActiveIndex az=m_mult2active[rt.m_Nk];
                assert(az.isValid());
                const Real piz=m_piActive[az];
                // errx=|d|pi_x + mu*dx*softmin0(piz)   [erry similar]
                // d/dpi_x errx = |d|
                // d/dpi_z errx = mu*dx*dsoftmin0(piz)
                const Real dminz = dsoftmin0(piz, m_minSmoothness);
                m_JacActive(ax,az) = mu*d[0]*dminz;
                m_JacActive(ay,az) = mu*d[1]*dminz;
            } 
        }
        ++nPairsChanged;
    }
    #ifndef NDEBUG
    if (nPairsChanged) {
        printf("Updated %d pairs of rows in Jacobian:", nPairsChanged);
        //cout << m_JacActive;
    }
    // Calculate Jacobian numerically.
    Array_<UniContactRT> uniContactTmp = uniContact;
    Vector piActive = m_piActive;
    Vector errActive0, errActive1;
    Matrix numJac(piActive.size(), piActive.size());
    for (int i=0; i < piActive.size(); ++i) {
        const Real save = piActive[i];
        piActive[i] = save - 1e-6;
        updateDirectionsAndCalcCurrentError(A,uniContactTmp,piELeft,
                                            piActive,errActive0);
        piActive[i] = save + 1e-6;
        updateDirectionsAndCalcCurrentError(A,uniContactTmp,piELeft,
                                            piActive,errActive1);
        numJac(i) = (errActive1-errActive0)/2e-6;
        piActive[i] = save;
    }
    //cout << "JacErr=" << m_JacActive-numJac;
    cout << "Jacobian num vs. analytic norm=" << (m_JacActive-numJac).norm() << endl;
    //cout << "USING NUMERICAL JAC\n"; m_JacActive = numJac;
    #endif
}


} // namespace SimTK