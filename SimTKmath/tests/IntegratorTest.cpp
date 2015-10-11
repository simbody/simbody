/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
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

/**@file
 * This is a test program which uses the Integrator class in various ways.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include "SimTKmath.h"

#include <cstdio>
#include <cassert>
#include <iostream>

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

static int qProj, qProjFail;
static int uProj, uProjFail;

// User-defined system to be integrated. This is a kind of SimTK::System.
class MyPendulum;
class MyPendulumGuts: public System::Guts {
    friend class MyPendulum;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex       massIndex, lengthIndex, gravityIndex;
    mutable QIndex                      q0;
    mutable UIndex                      u0;
    mutable QErrIndex                   qerr0;
    mutable UErrIndex                   uerr0;
    mutable UDotErrIndex                udoterr0;
    mutable EventTriggerByStageIndex    trigger0;
    mutable CacheEntryIndex             mgForceIndex; // a cache entry m*g calculated at Dynamics stage
    mutable EventId                     eventId0, eventId1, eventId2;
public:
    MyPendulumGuts() : Guts() {
        // Index types set themselves invalid on construction.
    }

    inline const MyPendulum& getMyPendulum() const;

    /*virtual*/MyPendulumGuts* cloneImpl() const override {return new MyPendulumGuts(*this);}

        /////////////////////////////////////////////////////////
        // Implementation of continuous DynamicSystem virtuals //
        /////////////////////////////////////////////////////////

    int realizeTopologyImpl(State&) const override;
    int realizeModelImpl(State&) const override;
    int realizeInstanceImpl(const State&) const override;
    int realizePositionImpl(const State&) const override;
    int realizeVelocityImpl(const State&) const override;
    int realizeDynamicsImpl(const State&) const override;
    int realizeAccelerationImpl(const State&) const override;

    // qdot==u here so these are just copies
    void multiplyByNImpl(const State& state, const Vector& u,
                         Vector& dq) const override {dq=u;}
    void multiplyByNTransposeImpl(const State& state, const Vector& fq,
                                  Vector& fu) const override {fu=fq;}
    void multiplyByNPInvImpl(const State& state, const Vector& dq,
                             Vector& u) const override {u=dq;}
    void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu,
                                      Vector& fq) const override {fq=fu;}

    // No prescribed motion.
    bool prescribeQImpl(State&) const override {return false;}
    bool prescribeUImpl(State&) const override {return false;}

    void projectQImpl(State&, Vector& qErrEst,
             const ProjectOptions& options, ProjectResults& results) const override;
    void projectUImpl(State&, Vector& uErrEst,
             const ProjectOptions& options, ProjectResults& results) const override;


        ////////////////////////////////////////////////
        // Implementation of discrete System virtuals //
        ////////////////////////////////////////////////

    int calcEventTriggerInfoImpl
       (const State& s, Array_<EventTriggerInfo>& eti) const override
    {
        eti.clear();
        eti.push_back(EventTriggerInfo(eventId0)
                      .setRequiredLocalizationTimeWindow(1)
                      .setTriggerOnRisingSignTransition(false));
        eti.push_back(EventTriggerInfo(eventId1));
        eti.push_back(EventTriggerInfo(eventId2)
                      .setTriggerOnFallingSignTransition(false));
        return 0;
    }

    int calcTimeOfNextScheduledEventImpl
       (const State& s, Real& tNextEvent,
        Array_<EventId>& eventIds, bool includeCurrentTime) const override
    {
        // Generate an event every 5.123 seconds.
        int nFives = (int)(s.getTime() / 5.123); // rounded down
        if (s.getTime()==0) nFives=1; // don't start with the event

        tNextEvent = nFives * Real(5.123);
        // Careful ...
        if (   tNextEvent < s.getTime()
            || (tNextEvent == s.getTime() && !includeCurrentTime))
            tNextEvent += Real(5.123);
        eventIds.push_back(eventId1); // event Id for scheduled pulse

        return 0;
    }


    // This should be called when the integrator returns indicating that
    // a discontinuity (event trigger) has been detected. The current
    // state is inconsistent in some way and we expect the event handlers
    // to correct that. Time will be the same before and after, but the
    // state may have changed discontinuously.
    void handleEventsImpl
       (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
        const HandleEventsOptions& options, HandleEventsResults& results) const
        override
    {
        cout << "===> t=" << s.getTime() << ": HANDLING "
             << Event::getCauseName(cause) << " EVENT!!!" << endl;
        if (eventIds.size())
            cout << "  EVENT IDS: " << eventIds << endl;
        if (cause != Event::Cause::TimeAdvanced) {
            std::swap(s.updQ()[0],s.updQ()[1]); // invalidates Position stage
            s.updU()=0;
        }
        results.setExitStatus(HandleEventsResults::Succeeded);
    }

};

// This is the handle class for a MyPendulum System.
// It must not have any data members. Data, if needed, is
// in the corresponding "Guts" class.

class MyPendulum: public System {
public:
    MyPendulum() : System()
    {
        adoptSystemGuts(new MyPendulumGuts());
        DefaultSystemSubsystem defsub(*this);
        updGuts().subsysIndex = defsub.getMySubsystemIndex();

        setHasTimeAdvancedEvents(false);
        (void)realizeTopology();
    }

    const MyPendulumGuts& getGuts() const {
        return dynamic_cast<const MyPendulumGuts&>(getSystemGuts());
    }

    MyPendulumGuts& updGuts() {
        return dynamic_cast<MyPendulumGuts&>(updSystemGuts());
    }

    // Instance variables are written to our defaultState.
    void setDefaultMass(Real mass) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.massIndex) = Value<Real>(mass);
    }

    void setDefaultLength(Real length) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.lengthIndex) = Value<Real>(length);
    }

    void setDefaultGravity(Real gravity) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.gravityIndex) = Value<Real>(gravity);
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        const MyPendulumGuts& guts = getGuts();
        updDefaultState().updU(guts.subsysIndex) = u;
        updDefaultState().updQ(guts.subsysIndex) = q;
        updDefaultState().updTime() = t;
    }

    Real getMass(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& m = s.getDiscreteVariable(guts.subsysIndex, guts.massIndex);
        return Value<Real>::downcast(m).get();
    }
    Real getDefaultMass() const {return getMass(getDefaultState());}

    Real getLength(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& d = s.getDiscreteVariable(guts.subsysIndex, guts.lengthIndex);
        return Value<Real>::downcast(d).get();
    }
    Real getDefaultLength() const {return getLength(getDefaultState());}

    Real getGravity(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& g = s.getDiscreteVariable(guts.subsysIndex, guts.gravityIndex);
        return Value<Real>::downcast(g).get();
    }
    Real getDefaultGravity() const {return getGravity(getDefaultState());}


    void dump(const char* msg) const {
        const MyPendulumGuts& guts = getGuts();
        cout << std::string(msg) << ": MyPendulum default state dump:" << endl;
        cout << "  mass   =" << getDefaultMass() << endl;
        cout << "  length =" << getDefaultLength() << endl;
        cout << "  gravity=" << getDefaultGravity() << endl;
        cout << "  time=" << getDefaultState().getTime() << endl;
        cout << "  q=" << getDefaultState().getQ(guts.subsysIndex) << endl;
        cout << "  u=" << getDefaultState().getU(guts.subsysIndex) << endl;
    }

};

inline const MyPendulum& MyPendulumGuts::getMyPendulum() const {
    return static_cast<const MyPendulum&>(getSystem());
}

static void printFinalStats(const Integrator& integ);

static void reportState(const char* msg,
                        const System& sys, const Integrator& integ) {
    if (*msg) printf("%s\n", msg);
    const State& s = integ.getState();

    sys.realize(s);
    printf(
        " -%6s- %9.6lf(%9.6lf) %14.10lf  %14.10lf  %14.10lf  %14.10lf | %14.10lf %14.10lf %14.10lf %14.10lf %14.10lf\n",
        integ.isStateInterpolated() ? "INTERP" : "------",
        s.getTime(), integ.getAdvancedTime(),
        s.getY()[0], s.getY()[1], s.getY()[2], s.getY()[3],
        s.getYErr()[0], s.getYErr()[1],
        s.getEventTriggersByStage(SubsystemIndex(0),Stage::Position)[0],
        s.getEventTriggersByStage(SubsystemIndex(0),Stage::Position)[1],
        s.getEventTriggersByStage(SubsystemIndex(0),Stage::Position)[2]);

    cout << "YDot:        " << s.getYDot() << endl;
    cout << "Multipliers: " << s.getMultipliers() << endl;
    cout << "UDotErrs:    " << s.getUDotErr() << endl;
}

int main () {
  try
  { MyPendulum sys;
    RungeKuttaMersonIntegrator integ(sys);
    //RungeKuttaFeldbergIntegrator integ(sys);
    //RungeKutta3Integrator integ(sys);
    //CPodesIntegrator integ(sys);
    //VerletIntegrator integ(sys);
    //ExplicitEulerIntegrator integ(sys);

    const Real t0=0;
    const Real qi[] = {1,0}; // (x,y)=(1,0)
    const Real ui[] = {0,0}; // v=0
    const Vector q0(2, qi);
    const Vector u0(2, ui);

    sys.setDefaultMass(10);
    sys.setDefaultTimeAndState(t0, q0, u0);
    sys.dump("Initial");

    integ.setAccuracy(1e-2);
    integ.setConstraintTolerance(1e-4);

    //integ.setAllowInterpolation(false);
    //integ.setProjectEveryStep(true);
    //integ.setProjectInterpolatedStates(false);
    //integ.setInitialStepSize(0.1);
    //integ.setUseInfinityNorm(true);
    //integ.setReturnEveryInternalStep(true);

    const Real tFinal = 30.003;
    const Real hReport = 1.;

    integ.setFinalTime(tFinal);
    integ.initialize(sys.getDefaultState());

    cout << "ACCURACY IN USE = " << integ.getAccuracyInUse() << endl;
    cout << "Initial y=" << integ.getState().getY() << endl;
    cout << "Initial ydot=" << integ.getState().getYDot() << endl;

    Real prevScheduledEventTime = -Infinity;
    for (int reportNo=0; !integ.isSimulationOver();
         reportNo += (integ.getTime() >= reportNo*hReport))
    {
        Array_<EventId> scheduledEventIds;
        Real nextScheduledEvent = NTraits<Real>::getInfinity();
        sys.calcTimeOfNextScheduledEvent(integ.getAdvancedState(),
            nextScheduledEvent, scheduledEventIds,
            integ.getAdvancedTime() > prevScheduledEventTime);

        HandleEventsOptions handleOpts(integ.getAccuracyInUse());
        HandleEventsResults handleResults;

        printf("----------------------------------------------------------\n");
        printf("stepTo(%g,%g)\n", reportNo*hReport, nextScheduledEvent);
        switch(integ.stepTo(reportNo*hReport, nextScheduledEvent)) {
            case Integrator::ReachedStepLimit: printf("STEP LIMIT\n"); break;
            case Integrator::ReachedReportTime: printf("REPORT TIME AT t=%.17g\n", integ.getTime()); break;
            case Integrator::StartOfContinuousInterval: printf("START OF CONTINUOUS INTERVAL\n"); break;

            case Integrator::ReachedScheduledEvent:  {
                Stage lowestModified = Stage::Empty;
                bool shouldTerminate;
                printf("SCHEDULED EVENT\n");
                reportState("BEFORE SCHEDULED EVENT:", sys, integ);
                sys.handleEvents(integ.updAdvancedState(),
                    Event::Cause::Scheduled,
                    scheduledEventIds,handleOpts,handleResults);
                shouldTerminate =
                    handleResults.getExitStatus()==HandleEventsResults::ShouldTerminate;
                lowestModified = handleResults.getLowestModifiedStage();
                integ.reinitialize(lowestModified, shouldTerminate);
                prevScheduledEventTime = integ.getAdvancedTime();
                break;
            }

            case Integrator::TimeHasAdvanced: {
                Stage lowestModified = Stage::Empty;
                bool shouldTerminate;
                printf("TIME HAS ADVANCED TO %g\n", integ.getTime());
                sys.handleEvents(integ.updAdvancedState(),
                    Event::Cause::TimeAdvanced,
                    Array_<EventId>(),handleOpts,handleResults);
                shouldTerminate =
                    handleResults.getExitStatus()==HandleEventsResults::ShouldTerminate;
                lowestModified = handleResults.getLowestModifiedStage();
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
                                              }

            case Integrator::ReachedEventTrigger: {
                Stage lowestModified = Stage::Empty;
                bool shouldTerminate;
                printf("EVENT TRIGGERED AT tLow=%.17g tHigh=%.17g!!\n",
                    integ.getTime(), integ.getAdvancedTime());
                cout << std::setprecision(17);
                cout << "Event window:     " << integ.getEventWindow() << endl;
                cout << "Triggered events: " << integ.getTriggeredEvents()<<"\n";
                cout << "Transitions seen:";
                for (int i=0; i<(int)integ.getEventTransitionsSeen().size(); ++i)
                    cout << " " << Event::eventTriggerString(integ.getEventTransitionsSeen()[i]);
                cout << endl;
                cout << "Est event times:  " << integ.getEstimatedEventTimes() << "\n";
                reportState("BEFORE TRIGGERED EVENT:", sys, integ);
                // state(t-) => state(t+)
                sys.handleEvents(integ.updAdvancedState(),
                    Event::Cause::Triggered,
                    integ.getTriggeredEvents(),handleOpts,handleResults);
                shouldTerminate =
                    handleResults.getExitStatus()==HandleEventsResults::ShouldTerminate;
                lowestModified = handleResults.getLowestModifiedStage();
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
            }

            case Integrator::EndOfSimulation: {
                Stage lowestModified = Stage::Empty;
                bool shouldTerminate;
                String reason = integ.getTerminationReasonString
                    (integ.getTerminationReason());

                printf("SIMULATION IS OVER. TERMINATION REASON=%s\n",
                       reason.c_str());
                sys.handleEvents(integ.updAdvancedState(),
                    Event::Cause::Termination,
                    Array_<EventId>(),handleOpts,handleResults);
                shouldTerminate =
                    handleResults.getExitStatus()==HandleEventsResults::ShouldTerminate;
                lowestModified = handleResults.getLowestModifiedStage();
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
            }

            default: assert(!"Unrecognized return from stepTo()");
        }

        // fall through to here to report
        reportState("", sys, integ);
    }

    printFinalStats(integ);

    return 0;
  }
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}

static void printFinalStats(const Integrator& integ)
{
  Real h0u;
  int nst, nattempt, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int nproj, nprojq, nproju, nce, nsetupsP, nprf, nprqf, npruf;

  h0u=NaN;
  nst=nattempt=nfe=nsetups=nje=nfeLS=nni=ncfn=netf=nge=-1;
  nproj=nprojq=nproju=nce=nsetupsP=nprf=nprqf=npruf=-1;

  /*
  flag = cpode.getActualInitStep(&h0u);
  flag = cpode.getNumSteps(&nst);
  flag = cpode.getNumFctEvals(&nfe);
  flag = cpode.getNumLinSolvSetups(&nsetups);
  flag = cpode.getNumErrTestFails(&netf);
  flag = cpode.getNumNonlinSolvIters(&nni);
  flag = cpode.getNumNonlinSolvConvFails(&ncfn);
  flag = cpode.dlsGetNumJacEvals(&nje);
  flag = cpode.dlsGetNumFctEvals(&nfeLS);
  flag = cpode.getProjStats(&nproj, &nce, &nsetupsP, &nprf);
  flag = cpode.getNumGEvals(&nge);
  */

  h0u   = integ.getActualInitialStepSizeTaken();
  nst   = integ.getNumStepsTaken();
  nattempt = integ.getNumStepsAttempted();
  nfe   = integ.getNumRealizations();
  netf  = integ.getNumErrorTestFailures();
  nproj = integ.getNumProjections();
  nprojq = integ.getNumQProjections();
  nproju = integ.getNumUProjections();
  nprf  = integ.getNumProjectionFailures();
  nprqf  = integ.getNumQProjectionFailures();
  npruf  = integ.getNumUProjectionFailures();

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6d nattempt = %-6d nfe  = %-6d nsetups = %-6d\n",
     nst, nattempt, nfe, nsetups);
  printf("nfeLS = %-6d nje = %d\n",
     nfeLS, nje);
  printf("nni = %-6d ncfn = %-6d netf = %-6d \n",
     nni, ncfn, netf);
  printf("nproj = %-6d nprojq = %-6d nproju = %-6d\n",
         nproj, nprojq, nproju);
  printf("nprf = %-6d nprqf = %-6d npruf = %-6d\n",
         nprf, nprqf, npruf);
  printf("nge = %d\n", nge);

  printf("qProj=%d qProjFail=%d\n", qProj, qProjFail);
  printf("uProj=%d uProjFail=%d\n", uProj, uProjFail);

}

/*
 * This system is a 2d pendulum swinging in gravity. It is modeled as
 * a point mass free in the plane, plus a distance constraint to model
 * the rod.
 *
 *    y       | g               O
 *    ^       v                  \  d
 *    |                           \
 *    |                            * m
 *     ------> x
 *
 * Gravity acts in the y direction, the rod is length d, mass m, pivot
 * location is the ground origin (0,0).
 *
 * The DAE for a generic multibody system is:
 *       qdot = Qu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r^2 = x^2  + y^2
 *       v^2 = x'^2 + y'^2
 * We will express the "rod length=d" constraint as
 *       (r^2 - d^2)/2 = 0    (perr)
 *           xx' + yy' = 0    (verr)
 *         xx'' + yy'' = -v^2 (aerr)
 *
 * So the matrix A = d perr/dq = [x y] and b = -v^2, and the
 * equations of motion are:
 *     [ m 0 x ] [ x'' ]   [  0  ]
 *     [ 0 m y ] [ y'' ] = [ -mg ]
 *     [ x y 0 ] [ L   ]   [-v^2 ]
 * where L (the Lagrange multiplier) is proportional to
 * the rod tension. You can solve this to get
 *     L   = (m*v^2 - mg*y)/(r^2)
 *     x'' = - x*L/m
 *     y'' = - y*L/m - g
 *
 */
int MyPendulumGuts::realizeTopologyImpl(State& s) const {
    // Instance variables mass, length, gravity
    massIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(1));
    lengthIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(1));
    gravityIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(13.7503716373294544));
    const Vector init(2, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    mgForceIndex = s.allocateCacheEntry(subsysIndex, Stage::Dynamics,
                                             new Value<Real>());
    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int MyPendulumGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int MyPendulumGuts::realizeInstanceImpl(const State& s) const {
    qerr0 = s.allocateQErr(subsysIndex, 1);
    uerr0 = s.allocateUErr(subsysIndex, 1);
    udoterr0 = s.allocateUDotErr(subsysIndex, 1); // and multiplier
    trigger0 = s.allocateEventTrigger(subsysIndex, Stage::Position, 3);
    eventId0 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    eventId1 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    eventId2 = getSystem().getDefaultSubsystem().createEventId(subsysIndex, s);
    System::Guts::realizeInstanceImpl(s);
    return 0;
}
int MyPendulumGuts::realizePositionImpl(const State& s) const {
    const Real    d = getMyPendulum().getLength(s);
    const Vector& q = s.getQ(subsysIndex);
    // This is the perr() equation.
    s.updQErr(subsysIndex)[0] = (q[0]*q[0] + q[1]*q[1] - d*d)/2;

    s.updEventTriggersByStage(subsysIndex, Stage::Position)[0] = 100*q[0]-q[1];

    // Make sure this boolean trigger *crosses* zero; it won't work right
    // if one end is actually zero. We'll use -.5 for false, .5 for true.
    s.updEventTriggersByStage(subsysIndex, Stage::Position)[1] =
        (s.getTime() > /*1.49552*/1.49545 && s.getTime() < 12.28937)-0.5;

    s.updEventTriggersByStage(subsysIndex, Stage::Position)[2] =
        s.getTime()-1.495508;
    System::Guts::realizePositionImpl(s);
    return 0;
}

int MyPendulumGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];

    // This is the verr() equation.
    s.updUErr(subsysIndex)[0]  = q[0]*u[0] + q[1]*u[1];
    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int MyPendulumGuts::realizeDynamicsImpl(const State& s) const {
    const Real m  = getMyPendulum().getMass(s);
    const Real g  = getMyPendulum().getGravity(s);

    Real& mg = Value<Real>::downcast(s.updCacheEntry(subsysIndex, mgForceIndex)).upd();
    // Calculate the force due to gravity.
    mg = m*g;
    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int MyPendulumGuts::realizeAccelerationImpl(const State& s) const {
    const Real m  = getMyPendulum().getMass(s);
    const Real g  = getMyPendulum().getGravity(s);
    // we're pretending we couldn't calculate this here!
    const Real mg = Value<Real>::downcast
                       (s.updCacheEntry(subsysIndex, mgForceIndex)).get();

    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       udot = s.updUDot(subsysIndex);
    Vector&       qdotdot = s.updQDotDot(subsysIndex);

    const Real r2 = q[0]*q[0] + q[1]*q[1];
    const Real v2 = u[0]*u[0] + u[1]*u[1];
    const Real L  = (m*v2 - mg*q[1])/r2;
    udot[0] = - q[0]*L/m;
    udot[1] = - q[1]*L/m - g;
    qdotdot = udot; // N=identity for this problem
    s.updMultipliers(subsysIndex)[0] = L;
    s.updUDotErr(subsysIndex)[0] = q[0]*udot[0] + q[1]*udot[1] + v2;
    System::Guts::realizeAccelerationImpl(s);
    return 0;
}

/*
 * Here we want to remove any constraint errors from the current state,
 * and project out any component of the integrator's error estimate
 * perpendicular to the constraint manifold. We will do this sequentially
 * rather than handling position and velocity simultaneously.
 *
 * For this system we have P = d perr/dq = V = d verr/du = [x y].
 * Weighted, we have PW=tp*[x/wx y/wy] VW=tv*[x/wxd y/wyd].
 * With pinv(A)=~A*(A*~A)^-1, we have:
 *
 *    pinv(P)  = ~[            x             y] /  (    x ^2+     y ^2)
 *    pinv(PW) = ~(1/tp)*[(wx *wy ^2)*x (wx ^2*wy) *y] / ((wy *x)^2+(wx *y)^2)
 *    pinv(VW) = ~(1/tv)*[(wxd*wyd^2)*x (wxd^2*wyd)*y] / ((wyd*x)^2+(wxd*y)^2)
 *      (the latter assuming x,y already projected on position manifold)
 *
 * We want to solve
 *    |perr(q0 - dq)|_TRMS <= accuracy, such that dq=min_WLS(dq)
 *    PW(q0) dq = Tp * perr(q0); q = q0-dq
 * Then
 *    |verr(q,u0 - du)|_TRMS <= accuracy, du=min_WLS(du)
 *    VW(q) du = Tv * verr(q,u0); u = u0-du
 *
 *
 * To remove the corresponding error estimates:
 *    PW(q) qperp = PW(q) qerrest; qerrest -= qperp
 *    VW(q) uperp = VW(q) uerrest; uerrest -= uperp
 *
 *
 */
static Real wrms(const Vector& y, const Vector& w) {
    Real sumsq = 0;
    for (int i=0; i<y.size(); ++i)
        sumsq += square(y[i]*w[i]);
    return std::sqrt(sumsq/y.size());
}

// qerrest is in/out
void MyPendulumGuts::projectQImpl(State& s, Vector& qerrest,
                                const ProjectOptions& opts,
                                ProjectResults& results) const

{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Real projLimit = opts.getProjectionLimit();
    const bool forceProj = opts.isOptionSet(ProjectOptions::ForceProjection);
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getQErrWeights(subsysIndex);
    // Since qdot=u here we can use uweights directly as qweights.
    const Vec2& wq = Vec2::getAs(&uweights[0]);
    const Real& tp = ctols[0]; // inverse tolerances 1/ti

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    Real& ep = s.updQErr(subsysIndex)[0]; // ep changes as we go

    results.setAnyChangeMade(false);

    //cout << "BEFORE wperr=" << tp*ep << endl;
    if (!forceProj && std::abs(tp*ep) <= consAccuracy) {
        results.setExitStatus(ProjectResults::Succeeded);
        return;
    }
    if (std::abs(tp*ep) > projLimit) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        ++qProjFail;
        return;
    }
    ++qProj;
    results.setAnyChangeMade(true);
    Real wqchg;
    do {
        // Position projection
        Real r2 = ~q*q; // x^2+y^2
        Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
        Row2 P(~q), PW(tp*q[0]/wq[0], tp*q[1]/wq[1]);
        Vec2 Pinv(q/r2);
        Vec2 PWinv = Vec2(square(wq[1])*wq[0]*q[0],
                            square(wq[0])*wq[1]*q[1]) / (tp*wqr2);
        Vec2 dq  = Pinv*(ep);      //cout << "dq=" << dq << endl;
        Vec2 wdq = PWinv*(tp*ep);  //cout << "wdq=" << wdq << endl;

        wqchg = std::sqrt(wdq.normSqr()/q.size()); // wrms norm

        s.updQ(subsysIndex)[0] -= wdq[0]/wq[0];
        s.updQ(subsysIndex)[1] -= wdq[1]/wq[1];
        realize(s, Stage::Position); // recalc QErr (ep)

        //cout << "AFTER q-=wdq/W wperr=" << tp*ep << " wqchg=" << wqchg << endl;
    } while (std::abs(tp*ep) > consAccuracy && wqchg >= 0.01*consAccuracy);

    //cout << "...AFTER wperr=" << tp*ep << endl;

    // Now do error estimates.

    if (qerrest.size()) {
        Vec2& eq = Vec2::updAs(&qerrest[0]);

        // Recalc PW, PWInv:
        const Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
        const Row2 PW = Row2(tp*q[0]/wq[0], tp*q[1]/wq[1]);
        const Vec2 PWinv = Vec2(wq[0]*square(wq[1])*q[0],
                                square(wq[0])*wq[1]*q[1]) / (tp*wqr2);

        Vec2 qperp = PWinv*(PW*eq);

        //cout << "ERREST before=" << yerrest
        //     << " wrms=" << wrms(qerrest,qweights) << endl;
        //cout << "PW*eq=" << PW*eq << endl;
        eq -= qperp;

        //cout << "ERREST after=" << yerrest
        //     << " wrms=" << wrms(qerrest,qweights) << endl;
        //cout << "PW*eq=" << PW*eq << endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

void MyPendulumGuts::projectUImpl(State& s, Vector& uerrest,
             const ProjectOptions& opts, ProjectResults& results) const
{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Real projLimit = opts.getProjectionLimit();
    const bool forceProj = opts.isOptionSet(ProjectOptions::ForceProjection);
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getUErrWeights(subsysIndex);

    const Vec2& wu = Vec2::getAs(&uweights[0]);
    const Real& tv = ctols[0];

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    const Vec2& u = Vec2::getAs(&s.getU(subsysIndex)[0]);
    Real& ev = s.updUErr(subsysIndex)[0]; // ev changes as we go

    results.setAnyChangeMade(false);

    //cout << "BEFORE wverr=" << tv*ev << endl;
    if (!forceProj && std::abs(tv*ev) <= consAccuracy) {
        results.setExitStatus(ProjectResults::Succeeded);
        return;
    }
    if (std::abs(tv*ev) > projLimit) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        ++uProjFail;
        return;
    }

    ++uProj;
    results.setAnyChangeMade(true);

    // Do velocity projection at current values of q, which should have
    // been projected already.
    Real r2 = ~q*q; // x^2+y^2
    Real wur2 = square(wu[1]*q[0]) + square(wu[0]*q[1]);
    Row2 V(~q), VW(tv*q[0]/wu[0], tv*q[1]/wu[1]);
    Vec2 Vinv(q/r2);
    Vec2 VWinv = Vec2(square(wu[1])*wu[0]*q[0],
                      square(wu[0])*wu[1]*q[1]) / (tv*wur2);
    realize(s, Stage::Velocity); // calculate UErr (ev)

    //cout << "BEFORE wverr=" << tv*ev << endl;
    Vec2 du  = Vinv*(ev);      //cout << "du=" << du << endl;
    Vec2 wdu = VWinv*(tv*ev);  //cout << "wdu=" << wdu << endl;

    s.updU(subsysIndex)[0] -= wdu[0]/wu[0];
    s.updU(subsysIndex)[1] -= wdu[1]/wu[1];

    realize(s, Stage::Velocity); // recalc UErr
    //cout << "AFTER u-=wdu wverr=" << tv*ev << endl;

    //cout << "...AFTER wverr=" << tv*ev << endl;

    // Now do error estimates.


    if (uerrest.size()) {
        Vec2& eu = Vec2::updAs(&uerrest[0]);
        Vec2 uperp = VWinv*(VW*eu);

        //cout << "ERREST before=" << uerrest
        //     << " wrms=" << wrms(uerrest,uweights) << endl;
        //cout << " VW*eu=" << VW*eu << endl;
        eu -= uperp;

        //cout << "ERREST after=" << yerrest
        //     << " wrms=" << wrms(uerrest,uweights) << endl;
        //cout << " VW*eu=" << VW*eu << endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}



