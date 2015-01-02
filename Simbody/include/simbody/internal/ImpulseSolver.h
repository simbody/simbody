#ifndef SimTK_SIMBODY_IMPULSE_SOLVER_H_
#define SimTK_SIMBODY_IMPULSE_SOLVER_H_

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

#include "SimTKmath.h"
#include "simbody/internal/common.h"

namespace SimTK {

/** This is the abstract base class for impulse solvers, which solve an
important subproblem of the contact and impact equations.

Impact problem: <pre>
    M  du    + ~G (pi+piE) = 0
    G (u+du) -  D (pi+piE) = b - verrNewton
</pre> where verrNewton is constraint space velocity error due to Newton
restitution (not used for Poisson restitution). Moving knowns to the right:
<pre>
    M  du    + ~G pi = -~G piE
    G  du    -  D pi = b - G u - verrNewton + D piE
</pre> Substituting 2nd eqn into first gives this impact subproblem: 
<pre>
    [A+D] pi = verr0 + verrNewton + verrExpand
    where verr0 = Gu-b, verrExpand = -[A+D]piE, and A=G M\ ~G.
</pre> 
%Contact problem: <pre>
    M du + ~G pi  = h f
    G du -  D pi  = b - G u
</pre> Substituting gives this contact subproblem:
<pre>
    [A+D] pi = verr0 + verrApplied
    where verrApplied = h G M\ f.
</pre>

The form of the problems is the same, with different RHS and no expansion
impulse given for contact. The impulse solver thus solves this system of 
equations and inequalities:
<pre>
    [A+D] (piExpand + piUnknown) = verrStart + verrApplied
    where verrStart = verr0 + verrNewton
</pre>
subject to several inequalities and replacement of friction rows by sliding or 
impending slip equations
for the unknown impulse piUnknown. piExpand is the given Poisson expansion
impulse, non-zero only for expanding unilateral normal contacts, with 
piExpand_z[k]<=0 for each unilateral contact k. A is an mXm symmetric positive
semidefinite matrix, D a diagonal matrix with nonnegative elements, 
pi=piExpand+piUnknown and rhs are m-vectors. We write pi this way because
friction limits depend on pi, not just piUnknown. We require piExpand_z<=0
for unilateral normal contacts z.

Similarly, we must separate verrStart and verrApplied, because verrStart
contains the actual constraint-space velocities, especially sliding, which
we need in order to determine sliding direction and contact status (rolling
or sliding). verrApplied just represents what the applied forces would do if
there were no constraints; the constraints will react to that so that the 
actual velocity changes will be much different.

For each "sliding step", we classify frictional contacts based on the current
contents of verrStart, then solve:
<pre>
    [A+D] piUnknown = verrStart + verrApplied - [A+D]*piExpand 
    with inequalities
    piUnknown_z <= 0
    sqrt(piUnknown_x^2+piUnknown_y^2) <= -mu*pi_z
    where pi=piUnknown+piExpand
</pre>
We then choose a fraction s of this sliding step to accept, with 0<s<=1, then
update <pre>
    verrStart   += s*(verrApplied - [A+D]*pi)
    piExpand    -= s*piExpand
    verrApplied -= s*verrApplied
</pre>
If s < 1 then we are not done. In that case we have removed some of the verr
and used up some of the expansion impulse. Return to do 
another sliding step until we take one where s==1.

We return piUnknown and the updated verrStart which would be zero if all 
contacts were active and rolling.

There are often multiple solutions for piUnknown; consult the documentation for 
particular ImpulseSolver implementations to determine which solution is 
returned. Possibilities include: any solution (PGS), and the least squares 
solution (PLUS).
**/

class SimTK_SIMBODY_EXPORT ImpulseSolver {
public:
    struct UncondRT;
    struct UniContactRT;
    struct UniSpeedRT;
    struct BoundedRT;
    struct ConstraintLtdFrictionRT;
    struct StateLtdFrictionRT;

    // How to treat a unilateral contact (input to solver).
    enum ContactType {TypeNA=-1, Observing=0, Known=1, Participating=2};

    // These are the results after the solve is complete. Don't mess with the
    // enum numbering here; we're counting on it.
    enum UniCond  {UniNA=-1, UniOff=0, UniActive=1, UniKnown=2};
    enum FricCond {FricNA=-1, FricOff=0, Sliding=1, Impending=2, Rolling=3};
    enum BndCond  {BndNA=-1, SlipLow=0, ImpendLow=1, Engaged=2, 
                   ImpendHigh=3, SlipHigh=4};


    ImpulseSolver(Real roll2slipTransitionSpeed,
                  Real convergenceTol,
                  int maxIters) 
    :   m_maxRollingTangVel(roll2slipTransitionSpeed),
        m_convergenceTol(convergenceTol),
        m_maxIters(maxIters)
    {
        clearStats();
    }

    virtual ~ImpulseSolver() {}

    void setMaxRollingSpeed(Real roll2slipTransitionSpeed) {
        assert(roll2slipTransitionSpeed >= 0);
        m_maxRollingTangVel = roll2slipTransitionSpeed; 
    }
    Real getMaxRollingSpeed() const {return m_maxRollingTangVel;}

    void setConvergenceTol(Real tol) {
        assert(tol >= 0);
        m_convergenceTol = tol;
    }
    Real getConvergenceTol() const {return m_convergenceTol;}

    void setMaxIterations(int maxIts) {
        assert(maxIts > 0);
        m_maxIters = maxIts;
    }
    int getMaxIterations() const {return m_maxIters;}

    // We'll keep stats separately for different "phases". The meaning of a
    // phase is up to the caller.
    static const int MaxNumPhases = 3;

    void clearStats() const {
        for (int i=0; i < MaxNumPhases; ++i)
            clearStats(i);
        m_nBilateralSolves = m_nBilateralIters = m_nBilateralFail = 0;
    }

    void clearStats(int phase) const {
        SimTK_ERRCHK2(0<=phase&&phase<MaxNumPhases,
            "ImpulseSolver::clearStats(phase)",
            "Phase must be 0..%d but was %d\n", MaxNumPhases-1, phase);
        m_nSolves[phase] = m_nIters[phase] = m_nFail[phase] = 0;
    }

    /** Solve. **/
    virtual bool solve
       (int                                 phase,
        const Array_<MultiplierIndex>&      participating, // p<=m of these 
        const Matrix&                       A,     // m X m, symmetric
        const Vector&                       D,     // m, diag>=0 added to A
        const Array_<MultiplierIndex>&      expanding, // nx<=m of these 
        Vector&                             piExpand, // m
        Vector&                             verrStart,   // m, RHS (in/out)
        Vector&                             verrApplied, // m
        Vector&                             pi,       // m, known+unknown
        Array_<UncondRT>&                   unconditional,
        Array_<UniContactRT>&               uniContact, // with friction
        Array_<UniSpeedRT>&                 uniSpeed,
        Array_<BoundedRT>&                  bounded,
        Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
        Array_<StateLtdFrictionRT>&         stateLtdFriction
        ) const = 0;


    /** Solve a set of bilateral (unconditional) constraints for the impulse
    necessary to enforce them. This can be used for projecting a set of
    violated active constraints onto their manifold. This just solves the 
    linear system <pre>
        P*(A+D)*~P P*pi = P*rhs 
                Pbar*pi = 0
    </pre>
    where P is a pXm "participation" matrix such that P(i,j)=1 if constraint
    j is the i'th active constraint, zero otherwise, and Pbar is the
    "nonparticipation" matrix such that Pbar(i,j)=1 if constraint j is the i'th
    inactive constraint, zero otherwise. A is mXm symmetric, positive 
    semidefinite, but may be rank deficient. D is an mXm diagonal matrix with
    only nonnegative elements. The returned solution pi (mX1) should ideally be 
    the solution of minimum 2-norm ||pi|| if the system is underdetermined, or 
    the solution that minimizes the 2-norm of the error ||(A+D)pi-rhs|| 
    (participating part only) if the solution is overdetermined and 
    inconsistent. However, concrete ImpulseSolvers are free to return a 
    different solution provide their behavior is well documented. The method 
    used should be qualitatively similar to that used by the solve() method for
    the same concrete %ImpulseSolver. For example, if solve() uses an
    iterative method then this should also do so. **/
    virtual bool solveBilateral
       (const Array_<MultiplierIndex>&      participating, // p<=m of these 
        const Matrix&                       A,     // m X m, symmetric
        const Vector&                       D,     // m, diag>=0 added to A
        const Vector&                       rhs,   // m, RHS
        Vector&                             pi     // m, unknown result
        ) const = 0;

    // Printable names for the enum values for debugging.
    static const char* getContactTypeName(ContactType ct);
    static const char* getUniCondName(UniCond uc);
    static const char* getFricCondName(FricCond fc);
    static const char* getBndCondName(BndCond bc);

    // Show details for each uni contact in the array.
    static void dumpUniContacts(const String& msg,
                                const Array_<UniContactRT>& uniContacts);

protected:
    Real m_maxRollingTangVel; // Sliding above this speed if solver cares.
    Real m_convergenceTol;    // Meaning depends on concrete solver.
    int  m_maxIters;          // Meaning depends on concrete solver.

    mutable long long m_nSolves[MaxNumPhases];
    mutable long long m_nIters[MaxNumPhases];
    mutable long long m_nFail[MaxNumPhases];
    mutable long long m_nBilateralSolves;
    mutable long long m_nBilateralIters;
    mutable long long m_nBilateralFail;
};

struct ImpulseSolver::UncondRT {
    UncondRT() {}

    // Input to solver.
    ConstraintIndex m_constraint;    // Back pointer to Simbody element.
    Array_<MultiplierIndex> m_mults; // Which constraint multipliers?

    // Set by solver on return.
    Array_<Real>            m_impulse; // Same size as m_mults.
};

// A unilateral contact (possibly with friction), joint stop, rope.
// These are the only constraints that can undergo impacts. Note that the COR
// is here for the convenience of the time stepper; it doesn't affect the
// impulse solvers. "Known" here means the normal constraint does not 
// participate (that is, the constraint equation cannot be active), but an 
// expansion impulse has been supplied for it.
struct ImpulseSolver::UniContactRT {
    UniContactRT() 
    :   m_sign(1), m_type(TypeNA), m_effCOR(NaN), m_effMu(NaN),
        m_contactCond(UniNA), m_frictionCond(FricNA), 
        m_impulse(NaN)
    {}

    bool hasFriction() const {return !m_Fk.empty();}

    // Input to solver.
    UnilateralContactIndex m_ucx;   // Back pointer to Simbody element.
    MultiplierIndex        m_Nk;    // multiplier for the normal constraint
    Real                   m_sign;  // sign convention for normal multiplier

    Array_<MultiplierIndex> m_Fk;   // optional friction multipliers

    // These solver inputs can change during a step.
    ContactType     m_type;       // Observing, Known, Participating
    Real            m_effCOR;     // velocity-dependent COR
    Real            m_effMu;      // if there is friction, else NaN

    // Working values for use by the solver, with final values returned.
    UniCond         m_contactCond;
    FricCond        m_frictionCond;
    Vec2            m_slipVel;
    Real            m_slipMag;
    Vec3            m_impulse;
};

// Ratchet.
struct ImpulseSolver::UniSpeedRT {
    UniSpeedRT(MultiplierIndex ix, int sign) 
    :   m_ix(ix), m_sign(sign), m_speedCond(UniNA), m_impulse(NaN)
    {   assert(sign==-1 || sign==1); }

    // Input to solver.
    MultiplierIndex m_ix;         // which constraint multiplier
    Real            m_sign;       // allowable sign for non-zero multiplier

    // Set by solver on return.
    UniCond         m_speedCond;
    Vec3            m_impulse;
};

// Torque-limited motor.
struct ImpulseSolver::BoundedRT {
    BoundedRT(MultiplierIndex ix, Real lb, Real ub) 
    :   m_ix(ix), m_lb(lb), m_ub(ub), m_boundedCond(BndNA), m_impulse(NaN) 
    {   assert(m_lb<=m_ub); }

    // Input to solver.
    MultiplierIndex m_ix;         // which constraint multiplier
    Real            m_lb, m_ub;   // effective lower, upper bounds; lb <= ub

    // Set by solver on return.
    BndCond         m_boundedCond;
    Real            m_impulse;
};

// Friction acting at a joint-like constraint, bead-on-a-wire.
struct ImpulseSolver::ConstraintLtdFrictionRT {
    ConstraintLtdFrictionRT
       (const Array_<MultiplierIndex>& frictionComponents, 
        const Array_<MultiplierIndex>& normalComponents,
        Real                           effMu)
    :   m_Fk(frictionComponents), m_Nk(normalComponents), 
        m_effMu(effMu), m_frictionCond(FricNA), 
        m_Fimpulse(frictionComponents.size(), NaN) 
    {   assert(m_Fk.size()<=3 && m_Nk.size()<=3); 
        assert(isNaN(m_effMu) || m_effMu>=0); }

    // Inputs to solver.
    ConstraintLimitedFrictionIndex m_clfx;   // Back pointer to Simbody element.
    Array_<MultiplierIndex> m_Fk, m_Nk;
    Real                    m_effMu;

    // Set by solver on return.
    FricCond                m_frictionCond;
    Array_<Real>            m_Fimpulse; // same size as m_Fk
};

// Friction acting at a compliant contact.
struct ImpulseSolver::StateLtdFrictionRT {
    StateLtdFrictionRT(const Array_<MultiplierIndex>& frictionComponents, 
                       Real knownN, Real muEff)
    :   m_Fk(frictionComponents), m_knownN(knownN), m_effMu(muEff), 
        m_frictionCond(FricNA), m_Fimpulse(frictionComponents.size(), NaN) 
    {   assert(m_Fk.size()<=3); assert(m_knownN >= 0);
        assert(isNaN(m_effMu) || m_effMu>=0); }

    // Inputs to solver.
    StateLimitedFrictionIndex m_slfx;   // Back pointer to Simbody element.
    Array_<MultiplierIndex> m_Fk;
    Real                    m_knownN;
    Real                    m_effMu;

    // Set by solver on return.
    FricCond                m_frictionCond;
    Array_<Real>            m_Fimpulse; // same size as m_Fk
};

} // namespace SimTK

#endif // SimTK_SIMBODY_IMPULSE_SOLVER_H_
