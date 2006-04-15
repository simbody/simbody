#ifndef SimTK_SIMBODY_NUMERICAL_METHODS_H_
#define SimTK_SIMBODY_NUMERICAL_METHODS_H_

#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyState.h"

#include <cassert>

namespace SimTK {

/**
 * This interface is to be implemented whenever you want to treat
 * your mechanical system as a DAE to be solved with coordinate projection.
 * In this formulation, we expect to compute state derivatives
 * directly by embedding the acceleration-level constraints, forming
 * an ODE. We integrate this system using an ordinary ODE method,
 * but then the solution will violate position and velocity constraints.
 * The coordinate projection scheme moves the solution back onto
 * the constraint manifold. Typically there is a lot of exploitable
 * structure to the constraints so we expect the user function here
 * to be capable of performing the projection itself, resulting
 * in a *small* change to the state. If the calling ODE solver
 * is a multistep method which must retain past solution values,
 * it needs to update its internal solution memory using the 
 * *projected* result.
 *
 * The state y is internally partitioned into three segments:
 *          y = {q,u,z)
 * although the caller does not need to know about this partioning. The
 * q's are configuration coordinates and participate in all the constraints.
 * u's are velocity coordinates (the qdots are some linear function of the u's)
 * and participate in velocity constraints. z's are "dynamic" state variables
 * typically associated with forces and are not involved in any of the
 * constraints (TODO: we could add dynamic constraints as well).
 *
 * The system of equations represented here is:
 *          ydot=f(t,y)       derivative function
 *          perr=p(t,q)       position constraint errors
 *          verr=v(t,q,u)     velocity constraint errors
 *
 * Scaling
 * -------
 * In internal coordinates, different state variable can be expected
 * to have very different scales. We expect the system here to understand
 * the scaling of each of its state variables. That is, it should be
 * able to provide "units" ui for each yi in y. Unit ui (>0) is the amount
 * we would consider a "unit" change in yi, measuring about its current
 * value. NOTE: "unit" here does not mean "small"! These units are
 * intended as a statement about the model, not the numerical solution,
 * so they are independent of the requested integration accuracy. The 
 * default unit is 1 which would treat all variables as equally important.
 * Also we expect to know for each variable whether it is appropriate
 * to use relative tolerance on that variable. Angles, for example, are
 * periodic so the accuracy requirement is independent of the current
 * value which may have "wrapped" many times. We'll define ri=1 if
 * relative error is OK, ri=0 if not.
 *
 * Given this knowledge, we expect the system to implement a method which,
 * given a scalar relative tolerance rtol and absolute tolerance atol,
 * can return a multiplicative weighting to be applied to estimated 
 * state errors, such that a resulting WRMS norm <= 1 would indicate that
 * the requested tolerances had been achieved.
 *
 * If an integrator solution is being performed with relative tolerance rtol
 * and absolute tolerance atol, then the weight for an error in yi is
 *          Wi = 1/(ri*rtol*|yi| + ui*atol)
 *
 * Constraint tolerances
 * ---------------------
 * We do not expect to project constraints every step, since the ODE 
 * automatically satisfies them at the acceleration level. However, we
 * do need to project them if they have drifted too far, and the 
 * acceptable drift depends on user-requested accuracy. As for scaling
 * of the state variables, we expect the model to understand the units for
 * each constraint violation, so that a "unit" change in any constraint
 * would be considered equally as egregious as a "unit" change in any other.
 * Again note that "unit" does not mean small. Typical choices for 
 * a translational constraint might be meters for large systems and nanometers
 * for molecules and in both cases one might consider 0.001 unit to be a
 * "small" violation (a millimeter or hundredth of an Angstrom). Velocity
 * constriants do not need separate units, but they require a notion of
 * timescale. That is, what do you consider a significant amount of time?
 * Call that time ts; then a "unit" velocity constraint error is that velocity
 * which would result in a unit position error after a time ts.
 *
 * With the constraint units and timescale known, we are given a single
 * scalar constraint tolerance which is interpreted as the maximum number
 * of units of error we can tolerate in any constraint. If the current
 * constraint violations are each smaller than that (in their own units)
 * then we won't project for that step. When we do project, the projection
 * should be done so as to reduce the violations substantially.
 *
 */
class MechanicalDAESystem {
public:
    virtual long size()                const = 0;  // problem size (|y|)

    // The solver should communicate this information at the beginning and then
    // leave it. If you only have one number, e.g. 0.001, it is reasonable to 
    // use it for both values here because they are intended to have roughly
    // the same scale, at least qualitatively.
    virtual void setAccuracy(const Real& solution, const Real& constraint) = 0;

    virtual void setState(const Real& t, const Vector& y) = 0;
    virtual bool realize()                          const = 0; // ODE only
    virtual bool realizeAndProject(bool& anyChange) = 0; // DAE using coordinate projection

    virtual const Real&   getTimescale() const = 0;

    virtual const Real&   getT()       const = 0;   // These are available after setState(),
    virtual const Vector& getY()       const = 0;   // though y is modified by realizeAndProject().
    virtual const Vector& getWeights() const = 0;

    virtual const Vector& getYDot()    const = 0;   // Available after either realize() call.
    virtual const Vector& getPositionError() const=0;
    virtual const Vector& getVelocityError() const=0;
    virtual const Vector& getAccelerationError() const=0;
};

/**
 * Abstract ODEIntegrator. We expect this to be able to
 * take a single, error-controlled step advancing the solution
 * y(t) to an ODE initial value problem of the form
 *           y' = f(t,y)
 *        y(t0) = y0
 * In general an integrator will maintain a complex internal
 * state which must be allocated when the problem size is
 * known.
 *
 * The function
 */
class MechanicalDAEIntegrator {
public:
    explicit MechanicalDAEIntegrator(MechanicalDAESystem& s) : mech(s) {
        initializeUserStuff();
    }
    virtual ~MechanicalDAEIntegrator() { }

    virtual MechanicalDAEIntegrator* clone()       const = 0;
    virtual const Real&    getT() const = 0;
    virtual const Vector&  getY() const = 0;

    virtual bool setInitialConditions(const Real& t0, const Vector& y0) = 0;

    virtual bool step(const Real& tout) = 0;

    virtual Real getConstraintTolerance() const = 0;

    void setStopTime(const Real& tstop) {
        assert(tstop == -1. || (0. <= tstop));
        userStopTime = tstop;
    }

    void setInitialStepSize(const Real& z) {
        assert(z == -1. || z > 0.);
        assert(userMinStepSize==-1. || z >= userMinStepSize);
        assert(userMaxStepSize==-1. || z <= userMaxStepSize);
        userInitStepSize = z;
    }
    void setMinimumStepSize(const Real& z) { 
        assert(z == -1. || z > 0.);
        assert(userInitStepSize==-1. || z <= userInitStepSize);
        assert(userMaxStepSize ==-1. || z <= userMaxStepSize);
        userMinStepSize = z;
    }
    void setMaximumStepSize(const Real& z) {
        assert(z == -1. || z > 0.);
        assert(userInitStepSize==-1. || z >= userInitStepSize);
        assert(userMinStepSize ==-1. || z >= userMinStepSize);
        userMinStepSize = z;
    }
    void setAccuracy(const Real& accuracy) {
        assert(accuracy == -1. || (0. < accuracy && accuracy < 1.));
        userAccuracy = accuracy;
    }
    void setTolerances(const Real& relTol, const Real& absTol, 
                       const Real& consTol) {
        assert(relTol  == -1. || (0. < relTol  && relTol  <= 1.));
        assert(absTol  == -1. || (0. < absTol  && absTol  <= 1.));
        assert(consTol == -1. || (0. < consTol && consTol <= 1.));

        userRelTol=relTol; userAbsTol=absTol; userConsTol=consTol;
    }
    void setVelocityConstraintRescale(const Real& rescale) {
        assert(rescale == -1 || (0. <= rescale));
        userVConsRescale = rescale;
    }

protected:
    MechanicalDAESystem& mech;

    // collect user requests
    Real userInitStepSize, userMinStepSize, userMaxStepSize;
    Real userAccuracy; // use for relTol, absTol, constraintTol
    Real userRelTol, userAbsTol, userConsTol, userVConsRescale; // for fussy people
    Real userStopTime; // never go past this
    long userMaxNumSteps; // that is, in a single call to step(); 0=no limit

    // Mark all user-supplied options "not supplied by user".
    void initializeUserStuff() {
        userInitStepSize = userMinStepSize = userMaxStepSize = -1.;
        userAccuracy = userRelTol = userAbsTol = userConsTol = userVConsRescale = -1.;
        userStopTime = -1.;
        userMaxNumSteps = -1;
    }

};

class ExplicitEuler : public MechanicalDAEIntegrator {
public:
    ExplicitEuler(MechanicalDAESystem& s) : MechanicalDAEIntegrator(s) {
        construct();
    }

    MechanicalDAEIntegrator* clone() const {return new ExplicitEuler(*this);}

    const Real&   getT() const {return t;}
    const Vector& getY() const {return y;}

    Real getConstraintTolerance() const {
        return consTol;
    }

    bool setInitialConditions(const Real& t0, const Vector& y0) {
        initializeIntegrationParameters();
        mech.setAccuracy(absTol, consTol);

        mech.setState(t0,y0);
        bool stateWasChanged;
        if (!mech.realizeAndProject(stateWasChanged)) 
            return false;
        
        t = t0; y = mech.getY(); // pick up projection if any
        initialized = true;
        return true;
    }

    bool step(const Real& tOut) {
        assert(initialized && tOut >= t);
        const Real tMax = std::min(tOut, stopTime);

        bool wasLastStep;
        long stepsTaken = 0;
        do {
            if (maxNumSteps && stepsTaken >= maxNumSteps)
                return false; // too much work

            // We expect to be realized at (t,y) already.
            ydot0 = mech.getYDot(); // save so we can restart

            // Take the biggest step we can get up to maxStep
            Real hTry = std::min(maxStepSize, tMax-t);
            for (;;) {
                wasLastStep = false;
                Real nextT = t+hTry;
                if (tMax-nextT < minStepSize)
                    nextT = tMax, hTry = tMax-t, wasLastStep=true;
                mech.setState(nextT, y + hTry*ydot0);
                bool stateWasChanged;
                if (mech.realizeAndProject(stateWasChanged)) {
                    // Took a step of size hTry. Update solution 
                    // and try another step. Note that we are picking up
                    // any projection changes by extracting y.
                    t = nextT; y = mech.getY();
                    break;
                }
                // Failed

                if (hTry <= minStepSize)
                    return false; // can't proceed
                // Reduce step size and try again
                hTry = std::max(minStepSize, hTry/2);
            }

            // Took a successful step of size hTry.
            ++stepsTaken;
        } while (!wasLastStep);

        // We reached tMax successfully.
        return true;
    }
private:

    
    void initializeIntegrationParameters() {
        initializeStepSizes();
        initializeTolerances();
        if (userStopTime != -1.) stopTime = userStopTime;
        else stopTime = CNT<Real>::getInfinity();
        if (userMaxNumSteps != -1) 
            maxNumSteps = userMaxNumSteps; // 0 means infinity
        else maxNumSteps = 0;
    }

    void initializeTolerances() {
        accuracy     = (userAccuracy     != -1. ? userAccuracy     : 1e-3);
        relTol       = (userRelTol       != -1. ? userRelTol       : accuracy); 
        absTol       = (userAbsTol       != -1. ? userAbsTol       : accuracy); 
        consTol      = (userConsTol      != -1. ? userConsTol      : accuracy); 
        vconsRescale = (userVConsRescale != -1. ? userVConsRescale : 1.); 
    }

    void initializeStepSizes() {
        if (userMaxStepSize != -1.) { // got max
            maxStepSize = userMaxStepSize;
            if (userInitStepSize != -1.) { // got max & init
                initStepSize = userInitStepSize;
                if (userMinStepSize != -1.)
                    minStepSize = userMinStepSize; // got min,max,init
                else minStepSize = std::min(maxStepSize/3, initStepSize);
            } else { // got max, not init
                initStepSize = maxStepSize;
                if (userMinStepSize != -1.)
                    minStepSize = userMinStepSize; // max & min, not init
                else minStepSize = maxStepSize/3;
            }
        } else { // didn't get max
            if (userInitStepSize != -1.) { // got init, not max
                initStepSize = userInitStepSize;
                maxStepSize = initStepSize;
                minStepSize = (userMinStepSize != -1. ? userMinStepSize : maxStepSize/3);
            } else { // didn't get init or max
                if (userMinStepSize != -1.) { // got only min
                    minStepSize = userMinStepSize;
                    maxStepSize = std::max(mech.getTimescale()/10., minStepSize);
                    initStepSize = maxStepSize;
                } else { // didn't get anything
                    maxStepSize = initStepSize = mech.getTimescale()/10.;
                    minStepSize = maxStepSize/3;
                }
            }
        }
    }

    void construct() {
        initStepSize=minStepSize=maxStepSize
            =accuracy=relTol=absTol=consTol=vconsRescale
            =stopTime=t=CNT<Real>::getNaN();
        maxNumSteps = -1;
        y.resize(mech.size());
        ydot0.resize(mech.size());
        initialized = false;
    }

    Real   initStepSize, minStepSize, maxStepSize;
    Real   accuracy, relTol, absTol, consTol, vconsRescale;
    Real   stopTime;
    long   maxNumSteps;
    Real   t;
    Vector y;
    Vector ydot0;

    bool initialized;
};

/*
class RungeKuttaMerson : public MechanicalDAEIntegrator {
public:
    minStep
    oneStepOnly
    nextStepSize
    OKtoPassTout
    relTol (scalar)
    absTol (n-vector)

    bool step(const Real& tStop) {
        // We expect to be realized at (t,y) already.
        bool success = false;

        assert(tStop >= t);
        maxStep = (tStop - t);

        if (maxStep < minStep) {
            ode.setT(tStop); ode.setY(y);
            if (!ode.realizeAndProject()) return false;
            if (ode.stateWasChanged())
                y = ode.getY();
            t = tStop;
            return true;
        }

        acceptableError = y*relTol + absTol;

        Real hNext = predictedStepSize;
        while (t < tStop) {
            bool hWasShrunk    = false;

            // Near tStop we'll shrink or stretch as appropriate. If we shrink
            // substantially we'll make note of that and not attempt to 
            // grow the step size if we succeed with this tiny thing.
            if (t + 1.1*hNext > tStop) {
                const Real hAdj = tStop - t;
                hWasShrunk = hAdj < 0.9*hNext;
                hNext = hAdj;
            }


            Real trialH = hNext;
            // Now we're going to attempt to take as big a step as we can
            // take (up to hNext). We'll keep shrinking trialH until something
            // works or we die.
            for (;;) {
                takeAnRK4MStep(t, y, ode.getYDot(), trialH,
                              ytmp[0], ytmp[1], ytmp[2], 
                              ynew, errEst);

                // We'll use infinity norm just to be conservative (all entries
                // are >= 0).
                scalarError = errEst.elementwiseDivide(acceptableError).max();

                if (scalarError <= 1) {
                    // Passed error test; see if we can evaluate at ynew.
                    ode.setT(t+hTrial); ode.setY(ynew);
                    if (ode.realizeAndProject()) {
                        // Accept this step (we'll pick up projection changes here).
                        t = ode.getT(); y = ode.getY();
                        // Can we double the step? TODO: should be more subtle
                        hNext = trialH;
                        if (!hWasShrunk) {
                            predictedNextH = trialH;
                            if (scalarError*64. <= 1.) 
                                predictedNextH *= 2.;
                        }
                        break; // done with this step
                    }
                }

                // Step did not succeed.
                if (trialH <= minStep ) {
                    // TODO: is it just a discontinuity?
                    return false;
                }
                else 
                    trialH = std::max(trialH/2, minStep);
            }
            // Completed a step at h=trialH, hNext is set properly
        }



        if (success) {
            t = ode.getT(); y = ode.getY();
            return 0;
        }

        actualDt = 0.;
        return 1;
    }
private:
    bool f(const Real& t, const Vector& y, Vector& yd) const {
        ode.setState(t,y);
        if (!ode.realize()) 
            return false;
        yd = ode.getYDot();
        return true;
    }

    bool takeAnRK4MStep(const Real& t0, const Vector& y0, const Vector& yd0,
                        const Real& h, 
                        Vector& ytmp, Vector& yd1, Vector& yd2,
                        Vector& y, Vector& errEst)
    {
        ytmp = y0 + (h/3)*yd0;
        if (!f(t0+h/3, ytmp, yd1)) return false;

        ytmp = y0 + (h/6)*(yd0+yd1);
        if (!f(t0+h/3, ytmp, yd1)) return false;

        ytmp = y0 + (h/8)*(yd0 + 3.*yd1);
        if (!f(t0+h/2, ytmp, yd2)) return false;

        ytmp = y0 + (h/2)*(yd0 - 3.*yd1 + 4.*yd2);
        if (!f(t0+h,   ytmp, yd1)) return false;

        // final
        y = y0 + (h/6)*(yd0 + 4.*yd2 + yd1);
        errEst = 0.2*(y-ytmp).abs();
    }

private:
    Real   t;   // current solution
    Vector y;

    Vector ytmp[4];
};


class MultibodyTimeStepper {
public:
    MultibodyTimeStepper(const MultibodySystem&, State& asModeled);

    void setInitialConditions(const State&);

    void step(const Real& dt, const Real& 


private:
    class MultibodyTimeStepperRep* rep;
};
*/

} // namespace SimTK

#endif // SimTK_SIMBODY_NUMERICAL_METHODS_H_
