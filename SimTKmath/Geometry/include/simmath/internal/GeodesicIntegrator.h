#ifndef SimTK_SIMMATH_GEODESIC_INTEGRATOR_H_
#define SimTK_SIMMATH_GEODESIC_INTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

/** @file
 * Contains a stripped-down low overhead numerical integrator for use with
 * small, fixed-sized sets of Differential-Algebraic equations such as arise
 * when computing geodesics over smooth surfaces.
 */

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

namespace SimTK {

/** This is a stripped-down numerical integrator for small ODE or DAE problems
whose size is known at compile time, with no provision for discrete variables,
event detection, or interpolation. You cannot use this integrator to advance
a Simbody System; see Integrator instead. Everything is defined in this header
file so that the integration can proceed with virtually no overhead. Templates
are used rather than run-time polymorphism, so there are no virtual function
calls. The system of equations is given as a template object that must
implement particular methods which the compiler may inline if they are simple
enough.

<h3>Class Equations</h3>
This integrator is instantiated with a class that encapsulates the system of
equations to be solved, and must provide compile time constants and methods 
with the following signatures:
@code
class MyEquations {
public:
    enum { NQ=5,       // number of 2nd order equations
           NC=3,       // total number of position & velocity constraints
           N = 2*NQ }; // total number of states y

    // Calculate state derivatives ydot given time and state y.
    void calcDerivs(Real t, const Vec<N>& y, Vec<N>& ydot) const;

    // Calculate amount by which the given time and state y violate the
    // constraints and return the constraint errors in cerr.
    void calcConstraintErrors(Real t, const Vec<N>& y, Vec<NC>& cerr) const;

    // Given a time and state y, ensure that the state satisfies the constraints
    // to within the indicated absolute tolerance, by performing the shortest
    // (i.e. least squares) projection of the state back to the constraint
    // manifold. Return false if the desired tolerance cannot be achieved. 
    // Otherwise (true return), a subsequent call to calcConstraintErrors() 
    // would return each |cerr[i]|<=consTol.
    bool projectIfNeeded(Real consTol, Real t, Vec<N>& y) const;
};

@endcode

<h3>Usage</h3>

@code
// Create an object of a type that satisfies the Equations description above.
MyEquations eqns(...); // ... is whatever arguments you require.
// Instantiate an integrator for this system of equations and specify the
// integration accuracy and constraint tolerance.
GeodesicIntegrator<MyEquations> integ(eqns, accuracy, constraintTol);
// Initialize; will project constraints if necessary.
integ.initialize(t0, y0);
// Integrate to time finalTime, getting output every completed step.
while (true) {
    std::cout << "t=" << integ.getTime() << " y=" << integ.getY() << "\n";
    if (integ.getTime() == finalTime)
        break;
    integ.takeOneStep(finalTime);
}
@endcode

<h3>Mathematical Overview</h3>

This is an explicit, variable-step integrator solving a 2nd-order DAE 
structured as an ODE-on-a-manifold system[1] like this:
<pre>
        (1)  udot = f(t,q,u)      NQ dynamic differential equations
        (2)  qdot = u             NQ kinematic differential equations
        (3)  0    = c(t,q,u)      NC constraints
</pre>
Here the "dot" suffix indicates differentiation with respect to the independent
variable t which we'll refer to as time here although it can be anything (for
geodesic calculations it is arc length). We'll  
call the second order variables q the "position variables", and their time 
derivatives u the "velocity variables". Collected together we call the state 
y={q,u}. At the beginning of a step, we expect to have been given initial 
conditions t0,q0,u0 such that |c(t0,q0,u0)|<=tol. The user provides the accuracy 
requirement and constraint tolerance. We solve the system to that accuracy while
keeping the constraints within tolerance. The integrator returns after taking
a successful step which may involve trial evaluations that are retracted.

By "ODE on a manifold" we mean that the ODE (1,2) automatically satisfies the 
condition that IF c==0, THEN cdot=0, where
<pre>
    cdot=Dc/Dt + Dc/Dq*qdot + Dc/Du*udot
</pre>
This means that satisfaction of the acceleration-level constraints is built
into the dynamic differential equations (1) so that we need only deal with
relatively slow drift of the solution away from the position and velocity
constraint manifolds.

To handle the constraint drift we use the method of coordinate projection and
expect the supplied Equations object to be able to perform a least-squares
projection of a state (q,u) to move it onto the constraint manifolds.

[1] Hairer, Lubich, Wanner, "Geometric Numerical Integration: 
Structure-Preserving Algorithms for Ordinary Differential Equations", 2nd ed., 
section IV.4, pg 109ff, Springer, 2006. **/
template <class Eqn>
class GeodesicIntegrator {
public:
    enum { NQ = Eqn::NQ, NC = Eqn::NC, N = 2*NQ };

    /** Construct an integrator for the given set of equations \a eqn, which are
    to be solved to the given \a accuracy, with constraints maintained to 
    within the given \a constraintTol. **/
    GeodesicIntegrator(const Eqn& eqn, Real accuracy, Real constraintTol) 
    :   m_eqn(eqn), m_accuracy(accuracy), m_consTol(constraintTol),
        m_hInit(NaN), m_hLast(NaN), m_hNext(Real(0.1)), m_nInitialize(0) 
    {   reinitializeCounters(); }

    /** Call this once before taking a series of steps. This sets the initial
    conditions, and calculates the starting derivatives and constraint errors.
    The constraints must be satisfied already by the given state; an error
    is thrown if not. **/
    void initialize(Real t, const Vec<N>& y) {
        ++m_nInitialize;
        reinitializeCounters();
        m_hInit = m_hLast = NaN;
        m_hNext = Real(0.1); // override if you have a better idea
        m_t = t; m_y = y;
        if (!m_eqn.projectIfNeeded(m_consTol, m_t, m_y)) {
            Vec<NC> cerr;
            m_eqn.calcConstraintErrors(m_t, m_y, cerr);
            const Real consErr = calcNormInf(cerr);
            SimTK_ERRCHK2_ALWAYS(!"projection failed",
                "GeodesicIntegrator::initialize()",
                "Couldn't project constraints to tol=%g;"
                " largest error was %g.", m_consTol, consErr);
        }
        m_eqn.calcDerivs(m_t, m_y, m_ydot);
    }

    /** Set initial time and state prior to integrating. State derivatives
    and constraint errors are calculated and an error is thrown if the 
    constraints are not already satisifed to the required tolerance. **/
    void setTimeAndState(Real t, const Vec<N>& y) {
        m_t = t; m_y = y;
        m_eqn.calcDerivs(m_t, m_y, m_ydot);
        Vec<NC> cerr;
        m_eqn.calcConstraintErrors(m_t, m_y, cerr);
        const Real consErr = calcNormInf(cerr);
        if (consErr > m_consTol)
            SimTK_ERRCHK2_ALWAYS(!"constraints not satisfied",
                "GeodesicIntegrator::setTimeAndState()",
                "Supplied state failed to satisfy constraints to tol=%g;"
                " largest error was %g.", m_consTol, consErr);
    }

    /** Use this if you think you know a better initial step size to try than
    the default. **/
    void setNextStepSizeToTry(Real h) {m_hNext=h;}
    /** Return the size of the next time step the integrator will attempt on
    the next call to takeOneStep(). **/
    Real getNextStepSizeToTry() const {return m_hNext;}

    /** Return the accuracy requirement as set in the constructor. **/
    Real getRequiredAccuracy() const {return m_accuracy;}
    /** Return the constraint tolerance as set in the constructor. **/
    Real getConstraintTolerance() const {return m_consTol;}

    /** Return the size of the first accepted step to be taken after the most
    recent initialize() call. **/
    Real getActualInitialStepSizeTaken() const {return m_hInit;}

    /** Return the number of successful time steps taken since the most recent
    initialize() call. **/
    int getNumStepsTaken() const {return m_nStepsTaken;}
    /** Return the total number of steps that were attempted since the most
    recent initialize() call. In general this will be more than the number
    of steps taken since some will be rejected. **/
    int getNumStepsAttempted() const {return m_nStepsAttempted;}
    /** How many steps were rejected because they did not satisfy the 
    accuracy requirement, since the most recent initalize() call. This is
    common but for non-stiff systems should be only a modest fraction of the
    number of steps taken. **/
    int getNumErrorTestFailures() const {return m_nErrtestFailures;}
    /** How many steps were rejected because the projectIfNeeded() method was
    unable to satisfy the constraint tolerance (since the most recent
    initialize() call). This should be very rare. **/
    int getNumProjectionFailures() const {return m_nProjectionFailures;}

    /** Return the number of calls to initialize() since construction of this
    integrator object. **/
    int getNumInitializations() const {return m_nInitialize;}

    /** Advance time and state by one error-controlled step and return, but 
    in no case advance past t=tStop. The integrator's internal time, 
    state, and state derivatives are advanced to the end of the step. If this
    step reaches \a tStop, the returned time will be \e exactly tStop. **/
    void takeOneStep(Real tStop);

    /** Return the current time. **/
    const Real& getTime() const {return m_t;}
    /** Return the complete current state as a Vec<N>. **/
    const Vec<N>&  getY() const {return m_y;}
    /** Return just the "position" variables q from the current state. **/
    const Vec<NQ>& getQ() const {return Vec<NQ>::getAs(&m_y[0]);}
    /** Return just the "velocity" variables u from the current state. **/
    const Vec<NQ>& getU() const {return Vec<NQ>::getAs(&m_y[NQ]);}
    /** Return the complete set of time derivatives of the current state. **/
    const Vec<N>&  getYDot() const {return m_ydot;}
    /** Return just the derivatives qdot of the "position" variables q. **/
    const Vec<NQ>& getQDot() const {return Vec<NQ>::getAs(&m_ydot[0]);}
    /** Return just the derivatives udot of the "velocity" variables u. **/
    const Vec<NQ>& getUDot() const {return Vec<NQ>::getAs(&m_ydot[NQ]);}

    /** This is a utility routine that returns the infinity norm (maximum
    absolute value) contained in a fixed-size, scalar Vec. **/
    template <int Z> static Real calcNormInf(const Vec<Z>& v) {
        Real norm = 0;
        for (int i=0; i < Z; ++i) {
            Real aval = std::abs(v[i]);
            if (aval > norm) norm = aval;
        }
        return norm;
    }

    /** This is a utility routine that returns the RMS norm of a fixed-size, 
    scalar Vec. **/
    template <int Z> static Real calcNormRMS(const Vec<Z>& v) {
        Real norm = 0;
        for (int i=0; i< Z; ++i) 
            norm += square(v[i]);
        return std::sqrt(norm/Z);
    }

private:
    void takeRKMStep(Real h, Vec<N>& y1, Vec<N>& y1err) const;

    const Eqn&      m_eqn;      // The DAE system to be solved.
    Real            m_accuracy; // Absolute accuracy requirement for y.
    Real            m_consTol;  // Absolute tolerance for constraint errors.

    Real            m_t;        // Current value of the independent variable.
    Vec<N>          m_y;        // Current q,u in that order.

    Real            m_hInit;    // Actual initial step taken.
    Real            m_hLast;    // Last step taken.
    Real            m_hNext;    // max step size to try next

    Vec<N>          m_ydot;     // ydot(t,y)

    // Counters.
    int m_nInitialize; // zeroed on construction only
    void reinitializeCounters() {
        m_nStepsTaken=m_nStepsAttempted=0;
        m_nErrtestFailures=m_nProjectionFailures=0;
    }
    int m_nStepsTaken;
    int m_nStepsAttempted;
    int m_nErrtestFailures;
    int m_nProjectionFailures;
};

template <class Eqn> void
GeodesicIntegrator<Eqn>::takeOneStep(Real tStop) {
    const Real Safety = Real(0.9), MinShrink = Real(0.1), MaxGrow = Real(5);
    const Real HysteresisLow =  Real(0.9), HysteresisHigh = Real(1.2);
    const Real MaxStretch = Real(0.1);
    const Real hMin = m_t <= 1 ? SignificantReal : SignificantReal*m_t;
    const Real hStretch = MaxStretch*m_hNext;

    // Figure out the target ending time for the next step. Choosing time
    // rather than step size lets us end at exactly tStop.
    Real t1 = m_t + m_hNext; // this is the usual case
    // If a small stretching of the next step would get us to tStop, try 
    // to make it all the way in one step.
    if (t1 + hStretch > tStop) 
        t1 = tStop;

    Real h, errNorm;
    Vec<N> y1, y1err;

    // Try smaller and smaller step sizes if necessary until we get one
    // that satisfies the error requirement, and for which projection
    // succeeds.
    while (true) {
        ++m_nStepsAttempted;
        h = t1 - m_t; assert(h>0);
        takeRKMStep(h, y1, y1err);
        errNorm = calcNormInf(y1err);
        if (errNorm > m_accuracy) {
            ++m_nErrtestFailures;
            // Failed to achieve required accuracy at this step size h.
            SimTK_ERRCHK4_ALWAYS(h > hMin,
                "GeodesicIntegrator::takeOneStep()", 
                "Accuracy %g worse than required %g at t=%g with step size"
                " h=%g; can't go any smaller.", errNorm, m_accuracy, m_t, h);

            // Shrink step by (acc/err)^(1/4) for 4th order.
            Real hNew = Safety * h * std::sqrt(std::sqrt(m_accuracy/errNorm));
            hNew = clamp(MinShrink*h, hNew, HysteresisLow*h);
            t1 = m_t + hNew;
            continue;
        }
        // Accuracy achieved. Can we satisfy the constraints?
        if (m_eqn.projectIfNeeded(m_consTol, t1, y1))
            break; // all good

        // Constraint projection failed. Hopefully that's because the
        // step was too big.
        ++m_nProjectionFailures;

        SimTK_ERRCHK3_ALWAYS(h > hMin,
            "GeodesicIntegrator::takeOneStep()", 
            "Projection failed to reach constraint tolerance %g at t=%g "
            "with step size h=%g; can't shrink step further.", 
            m_consTol, m_t, h);

        const Real hNew = MinShrink*h;
        t1 = m_t + hNew;
    }

    // We achieved desired accuracy at step size h, and satisfied the
    // constraints. (t1,y1) is now the final time and state; calculate 
    // state derivatives which will be used at the start of the next step.
    ++m_nStepsTaken;
    if (m_nStepsTaken==1) m_hInit = h; // that was the initial step
    m_t = t1; m_y = y1; m_hLast = h;
    m_eqn.calcDerivs(m_t, m_y, m_ydot);

    // If the step we just took ended right at tStop, don't use it to
    // predict a new step size; instead we'll just use the same hNext we
    // would have used here if it weren't for tStop.
    if (t1 < tStop) {
        // Possibly grow step for next time.
        Real hNew = errNorm == 0 ? MaxGrow*h
            :  Safety * h * std::sqrt(std::sqrt(m_accuracy/errNorm));
        if (hNew < HysteresisHigh*h) hNew = h; // don't bother
        hNew = std::min(hNew, MaxGrow*h);
        m_hNext = hNew;
    }
}

template <class Eqn> void
GeodesicIntegrator<Eqn>::takeRKMStep(Real h, Vec<N>& y1, Vec<N>& y1err) const {
    const Real h2=h/2, h3=h/3, h6=h/6, h8=h/8;
    const Real t0=m_t, t1=m_t+h;
    const Vec<N>& y0 = m_y;
    const Vec<N>& f0 = m_ydot;
    Vec<N> f1, f2, f3, f4, ysave;
    m_eqn.calcDerivs(t0+h3, y0 + h3* f0,         f1);
    m_eqn.calcDerivs(t0+h3, y0 + h6*(f0 + f1),   f2);
    m_eqn.calcDerivs(t0+h2, y0 + h8*(f0 + 3*f2), f3);

    // We'll need this for error estimation.
    ysave = y0 + h2*(f0 - 3*f2 + 4*f3);
    m_eqn.calcDerivs(t1, ysave, f4);

    // Final value. This is the 4th order accurate estimate for 
    // y1=y(t0+h)+O(h^5): y1 = y0 + (h/6)*(f0 + 4 f3 + f4). 
    // Don't evaluate here yet because we might reject this step or we
    // might need to do a projection.
    y1 = y0 + h6*(f0 + 4*f3 + f4);

    // This is an embedded 3rd-order estimate y1hat=y(t0+h)+O(h^4).
    //     y1hat = y0 + (h/10)*(f0 + 3 f2 + 4 f3 + 2 f4)
    // We don't actually have any need for y1hat, just its 4th-order
    // error estimate y1hat-y1=(1/5)(y1-ysave) (easily verified from the above).

    for (int i=0; i<N; ++i)
        y1err[i] = std::abs(y1[i]-ysave[i]) / 5;
}

} // namespace SimTK

#endif // SimTK_SIMMATH_GEODESIC_INTEGRATOR_H_
