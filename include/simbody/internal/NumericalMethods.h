#ifndef SimTK_SIMBODY_NUMERICAL_METHODS_H_
#define SimTK_SIMBODY_NUMERICAL_METHODS_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "SimTKcpodes.h"

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
 * constraints do not need separate units, but they require a notion of
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


/**
 * Abstract Index3 DAE Integrator using the Coordinate Projection
 * method. We expect this to be able to
 * take a single, error-controlled step advancing the solution
 * y(t) to a DAE initial value problem of the form
 *           y' = f(t,y)
 *           0  = c(t,y)  constraints (position & velocity)
 *
 *        y(t0) = y0      initial conditions
 *     c(t0,y0) = 0
 *
 * We expect to be given a DAE System to integrate which can
 * calculate the derivatives from the state and perform the 
 * necessary projections to solve c(t,y)=0.
 *
 * In general an integrator will maintain a complex internal
 * state which must be allocated when the problem size is
 * known.
 */
class MechanicalDAEIntegrator {
public:
    explicit MechanicalDAEIntegrator(const MultibodySystem& mb, State& s) 
        : mbs(mb), state(s) {
        initializeUserStuff();
        zeroStats();
    }
    virtual ~MechanicalDAEIntegrator() { }

    const MultibodySystem& getMultibodySystem() const {return mbs;}
    const State&           getState()           const {return state;}
    State&                 updState()                 {return state;}

    virtual MechanicalDAEIntegrator* clone()       const = 0;

    virtual bool initialize() = 0;
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
        userMaxStepSize = z;
    }
    void setAccuracy(const Real& accuracy) {
        assert(accuracy == -1. || (0. < accuracy && accuracy < 1.));
        userAccuracy = accuracy;
    }
    void setRelativeTolerance(const Real& relTol) {
        assert(relTol  == -1. || (0. < relTol  && relTol  <= 1.));
        userRelTol=relTol;
    }
    void setAbsoluteTolerance(const Real& absTol) {
        assert(absTol  == -1. || (0. < absTol  && absTol  <= 1.));
        userAbsTol=absTol;
    }
    void setConstraintTolerance(const Real& consTol) {
        assert(consTol == -1. || (0. < consTol && consTol <= 1.));
        userConsTol=consTol;
    }
    void setVelocityConstraintRescale(const Real& rescale) {
        assert(rescale == -1 || (0. <= rescale));
        userVConsRescale = rescale;
    }
    void setProjectEveryStep(bool forceProject) {
        userProjectEveryStep = forceProject ? 1 : 0;
    }

    long getStepsAttempted()     const {return statsStepsAttempted;}
    long getStepsTaken()         const {return statsStepsTaken;}
    long getErrorTestFailures()  const {return statsErrorTestFailures;}
    long getRealizeFailures()    const {return statsRealizationFailures;}
    long getProjectionFailures() const {return statsProjectionFailures;}
    long getStepSizeChanges()    const {return statsStepSizeChanges;}

protected:
    const MultibodySystem& mbs;
    State&                 state;

    // collect user requests
    Real userInitStepSize, userMinStepSize, userMaxStepSize;
    Real userAccuracy; // use for relTol, absTol, constraintTol
    Real userRelTol, userAbsTol, userConsTol, userVConsRescale; // for fussy people
    Real userStopTime; // never go past this
    long userMaxNumSteps; // that is, in a single call to step(); 0=no limit
    int  userProjectEveryStep; // -1 (not supplied), 0(false), 1(true)

    // Mark all user-supplied options "not supplied by user".
    void initializeUserStuff() {
        userInitStepSize = userMinStepSize = userMaxStepSize = -1.;
        userAccuracy = userRelTol = userAbsTol = userConsTol = userVConsRescale = -1.;
        userStopTime = -1.;
        userMaxNumSteps = -1;
        userProjectEveryStep= -1;
    }

    // Statistics
    long statsStepsAttempted;
    long statsStepsTaken;
    long statsErrorTestFailures;
    long statsRealizationFailures;
    long statsProjectionFailures;
    long statsStepSizeChanges;

    void zeroStats() {
        statsStepsAttempted = statsStepsTaken = statsErrorTestFailures
            = statsRealizationFailures = statsProjectionFailures = statsStepSizeChanges = 0;
    }

};

// This class implements the abstract CPodesSystem interface understood by
// our C++ interface to CPodes.
class CPodesMultibodySystem : public CPodesSystem {
public:
    CPodesMultibodySystem(MechanicalDAEIntegrator* p) : mdae(p) 
    {
    }

    // Override default implementations of these virtual functions.
    //   explicitODE()
    //   project()


    // Calculate ydot = f(t,y).
    int explicitODE(Real t, const Vector& y, Vector& ydot) const {
        mdae->updState().updY() = y;
        mdae->updState().updTime() = t;

        try { 
            mdae->getMultibodySystem().realize(mdae->getState(), Stage::Acceleration); 
        }
        catch(...) { return CPodes::RecoverableError; } // assume recoverable
        ydot = mdae->getState().getYDot();
        return CPodes::Success;
    }

    // Calculate yerr = c(t,y).
    int constraint(Real t, const Vector& y, Vector& yerr) const {
        mdae->updState().updY() = y;
        mdae->updState().updTime() = t;

        try { 
            mdae->getMultibodySystem().realize(mdae->getState(), Stage::Velocity); 
        }
        catch(...) { return CPodes::RecoverableError; } // assume recoverable
        yerr = mdae->getState().getYErr();
        return CPodes::Success;
    }

    // Given a state (t,y) not on the constraint manifold, return ycorr
    // such that (t,y+ycorr+eps) is on the manifold, with ||eps||_wrms <= epsProj. 
    // 'err' passed in as the integrator's current error estimate for state y;
    // optionally project it to eliminate the portion normal to the manifold.
    int project(Real t, const Vector& y, Vector& ycorr, Real epsProj, 
                Vector& err) const 
    {
        State& s = mdae->updState();
        s.updY() = y;
        s.updTime() = t;

        const MultibodySystem& mbs = mdae->getMultibodySystem();

        try {
            const Real tol = mdae->getConstraintTolerance();
            Vector yUnitWeights, unitTolerances;
            mbs.realize(s, Stage::Position);
            mbs.calcYUnitWeights(s, yUnitWeights);
            mbs.calcYErrUnitTolerances(s, unitTolerances);

            mbs.project(s, tol, yUnitWeights, unitTolerances, err);
        }
        catch (...) { return CPodes::RecoverableError; } // assume recoverable

        ycorr = s.getY()-y;
        return CPodes::Success;
    }

    /*
    virtual int  quadrature(Real t, const Vector& y, 
                            Vector& qout) const;
    virtual int  root(Real t, const Vector& y, const Vector& yp,
                      Vector& gout) const;
    virtual int  weight(const Vector& y, Vector& weights) const;
    virtual void errorHandler(int error_code, const char* module,
                              const char* function, char* msg) const;
    */
    
    MechanicalDAEIntegrator *mdae;
};

class CPodesIntegrator : public MechanicalDAEIntegrator {
public:
    CPodesIntegrator(const MultibodySystem& mb, State& s) 
      : MechanicalDAEIntegrator(mb,s), cpodes(0), sys(0)
    {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "CPodesIntegrator::CPodesIntegrator()");
        reconstructForNewModel();
    }

    // Virtual function implementations
    ~CPodesIntegrator() {
        delete cpodes;
        delete sys;
    }
    CPodesIntegrator* clone() const {assert(false); return 0;}

    const CPodes& getCPodes() const {assert(cpodes); return *cpodes;}

    bool initialize() {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "CPodesIntegrator::initialize()");
        if (state.getSystemStage() < Stage::Model)
            reconstructForNewModel();
        initializeIntegrationParameters();

        initialized = true;

        mbs.realize(state, Stage::Velocity);

        const int ny = state.getY().size();
        const int nc = state.getYErr().size();

        Vector ycorr(ny);
        Vector err(ny, Real(0));
        if (sys->project(state.getTime(), state.getY(), ycorr, 0./*ignored*/, err)
                != CPodes::Success)
            return false;
        state.updY() = state.getY() + ycorr;

        mbs.realize(state, Stage::Velocity);

        Vector ydot(ny);
        if (sys->explicitODE(state.getTime(), state.getY(), ydot) 
               != CPodes::Success)
            return false;

        int retval;
        if ((retval=cpodes->init(*sys, state.getTime(), state.getY(), ydot,
            CPodes::ScalarScalar, relTol, &absTol)) != CPodes::Success) 
        {
            printf("init returned %d\n", retval);
            return false;
        }
        cpodes->lapackDense(ny);
        cpodes->setNonlinConvCoef(0.01); // TODO (default is 0.1)
        cpodes->setMaxNumSteps(50000);
        cpodes->projDefine();
        /*
        if (nc) {
            cpodes->projInit(CPodes::L2Norm, CPodes::Nonlinear,
                    Vector(nc, getConstraintTolerance()));
            //cpodes->setProjUpdateErrEst(false);
            //cpodes->setProjFrequency(long proj_freq);
            //cpodes->setProjTestCnstr(true);
            //cpodes->setProjLsetupFreq(long proj_lset_freq);
            //cpodes->setProjNonlinConvCoef(.1);
            cpodes->lapackDenseProj(nc, ny, 
                //CPodes::ProjectWithSchurComplement
                CPodes::ProjectWithQRPivot
                //CPodes::ProjectWithLU
                );
        }
        /**/
        return true;
    }

    bool step(const Real& tout) {
        // Re-parametrizing or remodeling requires a new call to initialize().
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Instance,
            "CPodesIntegrator::step()");
        assert(initialized);

        Real   tret; // ignored
        Vector ypout(state.getY().size()); // ignored

        // TODO: deal with stop time
        int res = cpodes->step(tout, &tret, state.updY(), ypout, CPodes::Normal);
        return res >= 0;
    }

    Real getConstraintTolerance() const {
        assert(initialized);
        return consTol;
    }

    Real getPredictedNextStep() const {
        assert(initialized);
        Real hnext;
        (void)cpodes->getCurrentStep(&hnext);
        return hnext;
    }
private:
    void initializeIntegrationParameters() {
        mbs.realize(state, Stage::Instance);

        initializeStepSizes();
        initializeTolerances();  

        if (userStopTime != -1.) 
            cpodes->setStopTime(userStopTime);
        if (userMaxNumSteps != -1) 
            cpodes->setMaxNumSteps(userMaxNumSteps);
        if (userProjectEveryStep != -1)
            if (userProjectEveryStep==1)
                cpodes->setProjFrequency(1); // every step
    }

    void initializeStepSizes() {
        if (userInitStepSize != -1)
            cpodes->setInitStep(userInitStepSize);
        if (userMinStepSize != -1)
            cpodes->setMinStep(userMinStepSize);
        if (userMaxStepSize != -1)
            cpodes->setMaxStep(userMaxStepSize);
    }

    void initializeTolerances() {
        accuracy     = (userAccuracy     != -1. ? userAccuracy     : 1e-3);
        relTol       = (userRelTol       != -1. ? userRelTol       : accuracy); 
        absTol       = (userAbsTol       != -1. ? userAbsTol       : accuracy); 
        consTol      = (userConsTol      != -1. ? userConsTol      : accuracy); 
        vconsRescale = (userVConsRescale != -1. ? userVConsRescale : 1.); 
    }
private:
    void reconstructForNewModel() {
        initialized = false;
        delete cpodes;
        delete sys;
        mbs.realize(state, Stage::Model);
        accuracy=relTol=absTol=consTol=vconsRescale
            = CNT<Real>::getNaN();

        cpodes = new CPodes(CPodes::ExplicitODE, CPodes::BDF, CPodes::Newton);
        //cpodes = new CPodes(CPodes::ExplicitODE, CPodes::BDF, CPodes::Functional);
        sys = new CPodesMultibodySystem(this);
    }

    CPodes*                cpodes;
    CPodesMultibodySystem* sys;
    Real accuracy, relTol, absTol, consTol, vconsRescale;
    bool initialized;
};

class ExplicitEuler : public MechanicalDAEIntegrator {
public:
    ExplicitEuler(const MultibodySystem& mb, State& s) 
      : MechanicalDAEIntegrator(mb,s) 
    {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "ExplicitEuler::ExplicitEuler()");
        reconstructForNewModel();
    }

    ExplicitEuler* clone() const {return new ExplicitEuler(*this);}

    Real getConstraintTolerance() const {
        assert(initialized);
        return consTol;
    }

    Real getPredictedNextStep() const {
        assert(initialized);
        return maxStepSize;
    }

    bool initialize() {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "ExplicitEuler::initialize()");
        if (state.getSystemStage() < Stage::Model)
            reconstructForNewModel();
        initializeIntegrationParameters();

        try { 
            Vector yUnitWeights, unitTolerances;
            mbs.realize(state, Stage::Position);
            mbs.calcYUnitWeights(state, yUnitWeights);
            mbs.calcYErrUnitTolerances(state, unitTolerances);
            Vector dummyErrest; //TODO
            mbs.project(state, consTol, yUnitWeights, unitTolerances,
                        dummyErrest);
        }
        catch (...) { ++statsProjectionFailures; return false; }

        try { mbs.realize(state, Stage::Acceleration); }
        catch (...) { ++statsRealizationFailures; return false; }

        initialized = true;
        return true;
    }

    bool step(const Real& tOut) {
        // Re-parametrizing or remodeling requires a new call to initialize().
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Instance,
            "ExplicitEuler::step()");

        assert(initialized && tOut >= state.getTime());
        const Real tMax = std::min(tOut, stopTime);

        bool wasLastStep;
        long stepsTaken = 0;
        do {
            if (maxNumSteps && stepsTaken >= maxNumSteps)
                return false; // too much work

            mbs.realize(state, Stage::Acceleration);
            t0    = state.getTime();
            ydot0 = state.getYDot(); // save so we can restart

            // Take the biggest step we can get up to maxStep
            Real hTry = std::min(maxStepSize, tMax-t0);
            for (;;) {
                ++statsStepsAttempted;
                wasLastStep = false;
                Real nextT = t0+hTry;
                if (tMax-nextT < minStepSize)
                    nextT = tMax, hTry = tMax-t0, wasLastStep=true;

                state.updTime() = nextT;
                state.updY()    = state.getY() + hTry*ydot0;

                bool projectOK = false, realizeOK = false;
                try {
                    projectOK = true;
                    Vector yUnitWeights, unitTolerances;
                    mbs.realize(state, Stage::Position);
                    mbs.calcYUnitWeights(state, yUnitWeights);
                    mbs.calcYErrUnitTolerances(state, unitTolerances);
                    bool needProjection = true;
                    if (!projectEveryStep) {
                        mbs.realize(state, Stage::Velocity); // so we can get UErrs
                        const Real err = 
                            mbs.calcWeightedRMSNorm(state.getYErr(), unitTolerances);
                        if (err < 0.9*consTol)
                            needProjection = false; // far from being violated
                    }
                    if (needProjection) {
                        Vector dummyErrest; //TODO
                        mbs.project(state, consTol, yUnitWeights, unitTolerances,
                                    dummyErrest);
                    }
                }
                catch (...) { projectOK=false; }

                if (projectOK) {
                    try {
                        realizeOK = true;
                        mbs.realize(state, Stage::Acceleration);
                    } catch (...) { realizeOK = false; }
                    if (!realizeOK)  ++statsRealizationFailures;
                } else ++statsProjectionFailures;

                if (projectOK && realizeOK)
                    break; // took a step of size hTry.

                // Failed

                if (hTry <= minStepSize)
                    return false; // can't proceed
                // Reduce step size and try again
                hTry = std::max(minStepSize, hTry/2);
            }

            // Took a successful step of size hTry.
            ++stepsTaken;
            ++statsStepsTaken;
        } while (!wasLastStep);

        // We reached tMax successfully.
        return true;
    }
private:

    
    void initializeIntegrationParameters() {
        mbs.realize(state, Stage::Instance);

        initializeStepSizes();
        initializeTolerances();
        if (userStopTime != -1.) stopTime = userStopTime;
        else stopTime = CNT<Real>::getInfinity();
        if (userMaxNumSteps != -1) 
            maxNumSteps = userMaxNumSteps; // 0 means infinity
        else maxNumSteps = 0;
        if (userProjectEveryStep != -1)
            projectEveryStep = (userProjectEveryStep==1);
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
                    maxStepSize = std::max(mbs.calcTimescale(state)/10., minStepSize);
                    initStepSize = maxStepSize;
                } else { // didn't get anything
                    maxStepSize = initStepSize = mbs.calcTimescale(state)/10.;
                    minStepSize = maxStepSize/3;
                }
            }
        }
    }

    void reconstructForNewModel() {
        mbs.realize(state, Stage::Model);
        initStepSize=minStepSize=maxStepSize
            =accuracy=relTol=absTol=consTol=vconsRescale
            =stopTime=CNT<Real>::getNaN();
        maxNumSteps = -1;
        projectEveryStep = false;
        ydot0.resize(state.getY().size());
        initialized = false;
    }

    Real   initStepSize, minStepSize, maxStepSize;
    Real   accuracy, relTol, absTol, consTol, vconsRescale;
    Real   stopTime;
    long   maxNumSteps;
    bool   projectEveryStep;
    Real   t0;
    Vector ydot0;

    bool initialized;
};


class RungeKuttaMerson : public MechanicalDAEIntegrator {
public:
    RungeKuttaMerson(const MultibodySystem& mb, State& s, bool noProject=false) 
      : MechanicalDAEIntegrator(mb,s) 
    {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "RungeKuttaMerson::RungeKuttaMerson()");
        reconstructForNewModel();
        if (noProject) suppressProject = true;
    }

    RungeKuttaMerson* clone() const {return new RungeKuttaMerson(*this);}

    Real getConstraintTolerance() const {
        assert(initialized);
        return consTol;
    }
    Real getPredictedNextStep() const {
        assert(initialized);
        return predictedNextStep;
    }

    bool initialize() {
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Topology,
            "RungeKuttaMerson::initialize()");
        if (state.getSystemStage() < Stage::Model)
            reconstructForNewModel();

        initializeIntegrationParameters();

        if (!evaluateAndProject(true))
            return false;
        
        initialized = true;
        return true;
    }

    ////////
    // TODO: this is not preserving the old state properly. Can't just save continuous states!
    ////////
    bool step(const Real& tOut) {
        // Re-parametrizing or remodeling requires a new call to initialize().
        SimTK_STAGECHECK_GE_ALWAYS(state.getSystemStage(), Stage::Instance,
            "RungeKuttaMerson::step()");

        assert(initialized && tOut >= state.getTime());
        const Real tMax = std::min(tOut, stopTime);

        // Duck out now if there isn't enough time left to take a step.
        if (tMax-state.getTime() < minStepSize) {
            ++statsStepsAttempted;
            state.updTime() = tMax;

            if (!evaluateAndProject())
                return false;

            ++statsStepsTaken;
            return true;
        }

        bool wasLastStep;
        long stepsTaken = 0;

        // hWanted tracks the step size we would like to have taken if the
        // end of the step interval had not intervened.
        Real hWanted = std::min(predictedNextStep, maxStepSize);
        do {
            if (maxNumSteps && stepsTaken >= maxNumSteps)
                return false; // too much work

            // Keep track of whether we had to use an excessively small step size
            // in order to stop at tOut. In that case we don't want to use the
            // shrunken step to predict the next one.
            bool trialHWasLimitedByTOut = false;

            wasLastStep = false;

            t0 = state.getTime(); // save these so we can restart
            y0 = state.getY();

            // Near tStop we'll shrink or stretch as appropriate. If we shrink
            // substantially we'll make note of that and not attempt to 
            // grow the step size if we succeed with this tiny thing.
            Real trialH = hWanted;
            if (t0 + 1.05*trialH > tMax) {
                const Real hAdj = tMax - t0;
                trialHWasLimitedByTOut = hAdj < 0.9*trialH;
                trialH = hAdj;
                wasLastStep = true;
            }

            mbs.realize(state, Stage::Acceleration);
            ydot0 = state.getYDot();    // save for faster restart

            // Now we're going to attempt to take as big a step as we can
            // take (up to hNext). We'll keep shrinking trialH until something
            // works or we die.
            Real scalarError;
            bool RKstepOK;
            for (;;) {
                ++statsStepsAttempted;

                RKstepOK = 
                    takeAnRK4MStep(t0, y0, ydot0, trialH,
                                   ytmp[0], ytmp[1], ytmp[2], 
                                   ynew, errEst);

                if (!RKstepOK) ++statsRealizationFailures;

                if (RKstepOK) {
                    Vector wts;
                    mbs.calcYUnitWeights(state, wts);
                    scalarError = mbs.calcWeightedRMSNorm(errEst, wts) / accuracy;
                    //scalarError = mbs.calcYErrorNorm(state, errEst) / accuracy;

                    if (scalarError <= 1.) {
                        // Passed error test; see if we can evaluate at ynew.
                        state.updTime() = t0+trialH;
                        state.updY()    = ynew;

                        if (evaluateAndProject()) {
                            // Accept this step (we'll pick up projection changes here).
                            t0 = state.getTime(); y0 = state.getY();

                            // Adjust step size. Restrict growth to 5x, shrinking to 0.1x.
                            // We won't grow at all if the step we just executed was
                            // artificially shrunk due, e.g., to us being near tMax.
                            // We'll still allow shrinkage in that case, although that
                            // is unlikely to be indicated.
                            // This is a 4th order method so 1/4 is the correct exponent.
                            // Note the 0.9 safety factor applied here -- we always
                            // underpredict the step size by 10% to avoid thrashing.

                            Real optimalStepScale = 
                                scalarError > 0 ? 0.9/pow(scalarError, 0.25) : 5.;
                            optimalStepScale = 
                                std::min(std::max(optimalStepScale, 0.1), 5.);

                            // Inject a little stability here -- don't grow the step
                            // by less than 20% or shrink by less than 10%.
                            if (0.9 < optimalStepScale && optimalStepScale < 1.2)
                                optimalStepScale = 1.;

                            // Change the value of hWanted unless we just succeeded with
                            // an unnaturally small end-of-interval step.
                            if (!(trialHWasLimitedByTOut && optimalStepScale >= 1)) {
                                hWanted = optimalStepScale * trialH;
                                hWanted = std::min(std::max(hWanted, minStepSize), maxStepSize);
                                if (hWanted != trialH)
                                    ++statsStepSizeChanges;
                            }
                            predictedNextStep = hWanted;
                            break; // done with this step
                        }
                        ++statsRealizationFailures;
                        ++statsProjectionFailures;
                    } else 
                        ++statsErrorTestFailures;
                }

                // Step did not succeed.
                if (trialH <= minStepSize )
                    return false;

                Real optimalStepScale;
                if (RKstepOK) { 
                    optimalStepScale = scalarError > 1. ? 0.9/pow(scalarError, 0.25) : 1.;
                    optimalStepScale = 
                        std::max(optimalStepScale, 0.1);
                } else
                    optimalStepScale = 0.5;

                hWanted = optimalStepScale * trialH;
                hWanted = std::min(std::max(hWanted, minStepSize), maxStepSize);
                if (hWanted != trialH)
                    ++statsStepSizeChanges;

                trialH = hWanted;
                trialHWasLimitedByTOut = false;
                wasLastStep = false;
            }
            // Completed a step at h=trialH, hNext is set properly
            ++stepsTaken;
            ++statsStepsTaken;
        } while (!wasLastStep);

        return true;
    }

private:
    bool f(const Real& t, const Vector& y, Vector& yd) const {
        state.updY() = y;
        state.updTime() = t;

        try { mbs.realize(state, Stage::Acceleration); }
        catch(...) { return false; }
        yd = state.getYDot();
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

        for (int i=0; i<y.size(); ++i)
            errEst[i] = 0.2*std::fabs(y[i]-ytmp[i]);

        return true;
    }

    // Realize through Position stage, project, then realize
    // through Acceleration stage.
    bool evaluateAndProject(bool forceProject=false) {
        bool projectOK = false, realizeOK = false;
        try {
            projectOK = true;
            Vector yUnitWeights, unitTolerances;
            mbs.realize(state, Stage::Position);
            mbs.calcYUnitWeights(state, yUnitWeights);
            mbs.calcYErrUnitTolerances(state, unitTolerances);
            bool needProjection = true;
            if (!projectEveryStep) {
                mbs.realize(state, Stage::Velocity); // so we can get UErrs
                const Real err = 
                    mbs.calcWeightedRMSNorm(state.getYErr(), unitTolerances);
                if (err < 0.9*consTol)
                    needProjection = false; // far from being violated
            }
            if ((needProjection && !suppressProject) || forceProject) {
                Vector dummyErrest; //TODO
                mbs.project(state, consTol, yUnitWeights, unitTolerances,
                            dummyErrest);
            }
        }
        catch (...) { projectOK=false; }

        if (projectOK) {
            try {
                realizeOK = true;
                mbs.realize(state, Stage::Acceleration);
            } catch (...) { realizeOK = false; }
            if (!realizeOK)  ++statsRealizationFailures;
        } else ++statsProjectionFailures;

        return projectOK && realizeOK;
    }

private:    
    void initializeIntegrationParameters() {
        mbs.realize(state, Stage::Instance);

        initializeStepSizes();
        predictedNextStep = initStepSize;
        initializeTolerances();        
        if (userStopTime != -1.) stopTime = userStopTime;
        else stopTime = CNT<Real>::getInfinity();
        if (userMaxNumSteps != -1) 
            maxNumSteps = userMaxNumSteps; // 0 means infinity
        else maxNumSteps = 0;
        if (userProjectEveryStep != -1)
            projectEveryStep = (userProjectEveryStep==1);
    }

    void initializeTolerances() {
        accuracy     = (userAccuracy     != -1. ? userAccuracy     : 1e-3);
        relTol       = (userRelTol       != -1. ? userRelTol       : accuracy); 
        absTol       = (userAbsTol       != -1. ? userAbsTol       : accuracy); 
        consTol      = (userConsTol      != -1. ? userConsTol      : accuracy); 
        vconsRescale = (userVConsRescale != -1. ? userVConsRescale : 1.); 
    }

    void initializeStepSizes() {
        const Real MinStep = NTraits<Real>::Eps_34; // e.g., 1e-12 in double
        if (userMaxStepSize != -1.) { // got max
            maxStepSize = userMaxStepSize;
            if (userInitStepSize != -1.) { // got max & init
                initStepSize = userInitStepSize;
                if (userMinStepSize != -1.)
                    minStepSize = userMinStepSize; // got min,max,init
                else minStepSize = std::min(MinStep, initStepSize);
            } else { // got max, not init
                if (userMinStepSize != -1.) {
                    minStepSize = userMinStepSize; // max & min, not init
                    initStepSize = std::max(minStepSize, mbs.calcTimescale(state)/10.);
                } else {
                    // got max only
                    initStepSize = std::min(mbs.calcTimescale(state)/10., maxStepSize);
                    minStepSize  = std::min(MinStep, initStepSize);
                }
            }
        } else { // didn't get max
            maxStepSize = CNT<Real>::getInfinity();
            if (userInitStepSize != -1.) { // got init, not max
                initStepSize = userInitStepSize;
                if (userMinStepSize != -1.) minStepSize = userMinStepSize;
                else minStepSize = std::min(MinStep, initStepSize);
            } else { // didn't get init or max
                if (userMinStepSize != -1.) { // got only min
                    minStepSize = userMinStepSize;
                    initStepSize = std::max(mbs.calcTimescale(state)/10., minStepSize);
                } else { // didn't get anything
                    initStepSize = mbs.calcTimescale(state)/10.;
                    minStepSize  = std::min(MinStep, initStepSize);
                }
            }
        }
    }

    static const int NTemps = 4;

    void reconstructForNewModel() {
        mbs.realize(state, Stage::Model);
        initStepSize=minStepSize=maxStepSize
            =accuracy=relTol=absTol=consTol=vconsRescale
            =stopTime=predictedNextStep=CNT<Real>::getNaN();
        maxNumSteps = -1;
        projectEveryStep = false;
        suppressProject = false;
        y0.resize(state.getY().size());
        ydot0.resize(state.getY().size());
        errEst.resize(state.getY().size());
        ynew.resize(state.getY().size());
        for (int i=0; i<NTemps; ++i)
            ytmp[i].resize(state.getY().size());
        initialized = false;
    }

    Real   initStepSize, minStepSize, maxStepSize;
    Real   accuracy, relTol, absTol, consTol, vconsRescale;
    Real   stopTime;
    long   maxNumSteps;
    bool   projectEveryStep;
    bool   suppressProject;
    Real   t0;
    Vector y0, ydot0, errEst, ynew;
    Vector ytmp[NTemps];
    Real predictedNextStep;

    bool initialized;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_NUMERICAL_METHODS_H_
