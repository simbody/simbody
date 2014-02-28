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

/** This is the abstract base class for impulse solvers.
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
    enum ContactType {TypeNA=-1, Observe=0, Known=1, Participate=2};

    // These are the results after the solve is complete. Don't mess with the
    // enum numbering here; we're counting on it.
    enum UniCond  {UniNA=-1, UniOff=0, UniActive=1};
    enum FricCond {FricNA=-1, FricOff=0, Sliding=1, Impending=2, Rolling=3};
    enum BndCond  {BndNA=-1, SlipLow=0, ImpendLow=1, Engaged=2, 
                   ImpendHigh=3, SlipHigh=4};


    ImpulseSolver() {
        clearStats();
    }

    // We'll keep stats separately for different "phases". The meaning of a
    // phase is up to the caller.
    static const int MaxNumPhases = 3;

    void clearStats() const {
        for (int i=0; i < MaxNumPhases; ++i)
            clearStats(i);
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
        const Vector&                       D,     // m, diag >= 0 added to A
        const Vector&                       verr,  // m, RHS
        Vector&                             pi,    // m, initial guess & result
        Array_<UncondRT>&                   unconditional,
        Array_<UniContactRT>&               uniContact, // with friction
        Array_<UniSpeedRT>&                 uniSpeed,
        Array_<BoundedRT>&                  bounded,
        Array_<ConstraintLtdFrictionRT>&    consLtdFriction,
        Array_<StateLtdFrictionRT>&         stateLtdFriction
        ) const = 0;

    // Printable names for the enum values for debugging.
    static const char* getContactTypeName(ContactType ct);
    static const char* getUniCondName(UniCond uc);
    static const char* getFricCondName(FricCond fc);
    static const char* getBndCondName(BndCond bc);

protected:
    mutable long long m_nSolves[MaxNumPhases];
    mutable long long m_nIters[MaxNumPhases];
    mutable long long m_nFail[MaxNumPhases];
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
// impulse solvers.
struct ImpulseSolver::UniContactRT {
    UniContactRT() 
    :   m_sign(1), m_type(TypeNA), m_effCOR(NaN), m_effMu(NaN),
        m_knownPi(NaN), m_contactCond(UniNA), m_frictionCond(FricNA), 
        m_impulse(NaN)
    {}

    bool hasFriction() const {return !m_Fk.empty();}

    // Input to solver.
    UnilateralContactIndex m_ucx;   // Back pointer to Simbody element.
    MultiplierIndex        m_Nk;    // multiplier for the normal constraint
    Real                   m_sign;  // sign convention for normal multiplier

    Array_<MultiplierIndex> m_Fk;   // optional friction multipliers

    // These solver inputs can change during a step.
    ContactType     m_type;       // Observe, Known, Participate
    Real            m_effCOR;     // velocity-dependent COR
    Real            m_effMu;      // if there is friction, else NaN
    Real            m_knownPi;    // known impulse if status==Known

    // Set by solver on return.
    UniCond         m_contactCond;
    FricCond        m_frictionCond;
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
