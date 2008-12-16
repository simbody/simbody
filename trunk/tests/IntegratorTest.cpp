/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
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

// User-defined system to be integrated. This is a kind of SimTK::System.
class MyPendulum;
class MyPendulumGuts: public System::Guts {
    friend class MyPendulum;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
    mutable int massIndex, lengthIndex, gravityIndex;
    mutable int q0, u0, qerr0, uerr0, udoterr0, event0;
    mutable int mgForceIndex; // a cache entry m*g calculated at Dynamics stage
    mutable EventId eventId0, eventId1, eventId2;
public:
    MyPendulumGuts() : Guts() {
        subsysIndex = InvalidSubsystemIndex;
        massIndex = lengthIndex = gravityIndex = -1;
        q0 = u0 = qerr0 = uerr0 = udoterr0 = event0 = -1;
        mgForceIndex = -1;
    }

    const MyPendulum& getMyPendulum() const {
        return reinterpret_cast<const MyPendulum&>(getSystem());
    }

    /*virtual*/MyPendulumGuts* cloneImpl() const {return new MyPendulumGuts(*this);}

        /////////////////////////////////////////////////////////
        // Implementation of continuous DynamicSystem virtuals //
        /////////////////////////////////////////////////////////

    /*virtual*/int realizeTopologyImpl(State&) const;
    /*virtual*/int realizeModelImpl(State&) const;
    /*virtual*/int realizeInstanceImpl(const State&) const;
    /*virtual*/int realizePositionImpl(const State&) const;
    /*virtual*/int realizeVelocityImpl(const State&) const;
    /*virtual*/int realizeDynamicsImpl(const State&) const;
    /*virtual*/int realizeAccelerationImpl(const State&) const;

    /*virtual*/Real calcTimescaleImpl(const State& s) const {
        assert(s.getSystemStage() >= Stage::Instance);
        return 1;
    }


    /*virtual*/int calcYUnitWeightsImpl(const State& s, Vector& wts) const {
        assert(s.getSystemStage() >= Stage::Position);
        wts.resize(s.getNY());
        wts = 1; wts[1]=1;
        return 0;
    }


    // Returns *inverse* tols 1/ti.
    /*virtual*/int calcYErrUnitTolerancesImpl(const State& s, Vector& ootols) const {
        assert(s.getSystemStage() >= Stage::Instance);
        ootols.resize(s.getNYErr());
        ootols=1; ootols[0]=1; ootols[1]=1;
        return 0;
    }


    /*virtual*/int projectImpl(State&, Real consAccuracy, const Vector& yweights,
                           const Vector& ctols, Vector& yerrest, System::ProjectOptions) const;


        ////////////////////////////////////////////////
        // Implementation of discrete System virtuals //
        ////////////////////////////////////////////////

    /*virtual*/int calcEventTriggerInfoImpl(const State& s, std::vector<System::EventTriggerInfo>& eti) const {
        eti.clear();
        eti.push_back(System::EventTriggerInfo(eventId0)
                      .setRequiredLocalizationTimeWindow(1)
                      .setTriggerOnRisingSignTransition(false));
        eti.push_back(System::EventTriggerInfo(eventId1));
        eti.push_back(System::EventTriggerInfo(eventId2) 
                      .setTriggerOnFallingSignTransition(false));
        return 0;
    }

    /*virtual*/int calcTimeOfNextScheduledEventImpl(const State& s, Real& tNextEvent, 
                                                    std::vector<int>& eventIds, bool includeCurrentTime) const
    {
        // Generate an event every 5.123 seconds.
        int nFives = (int)(s.getTime() / 5.123); // rounded down
        tNextEvent = (nFives+1) * Real(5.123);
        // Careful ...
        if (tNextEvent <= s.getTime())
            tNextEvent += Real(5.123);
        eventIds.push_back(eventId1); // event Id for scheduled pulse

        return 0;
    }


    // This should be called when the integrator returns indicating that
    // a discontinuity (event trigger) has been detected. The current
    // state is inconsistent in some way and we expect the event handlers
    // to correct that. Time will be the same before and after, but the
    // state may have changed discontinuously. At least the discrete state
    // variables which indicate unhandled triggers must be cleared.
    /*virtual*/int handleEventsImpl
       (State& s, System::EventCause cause, const std::vector<int>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const
    {
        cout << "===> t=" << s.getTime() << ": HANDLING " 
             << System::getEventCauseName(cause) << " EVENT!!!" << endl;
//        if (eventIds.size())
//            cout << "  EVENT IDS: " << eventIds << endl;

        return 0;
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

static void printFinalStats(const Integrator& integ);
int main () {
  try 
  { MyPendulum sys;
    RungeKuttaMersonIntegrator integ(sys);

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

    const Real tFinal = 30.003;
    const Real hReport = 1.;

    integ.setFinalTime(tFinal);
    integ.initialize(sys.getDefaultState());

    cout << "ACCURACY IN USE = " << integ.getAccuracyInUse() << endl;


    for (int reportNo=0; !integ.isSimulationOver(); 
         reportNo += (integ.getTime() >= reportNo*hReport))
    {
        std::vector<EventId> scheduledEventIds;
        Real nextScheduledEvent = NTraits<Real>::getInfinity();
        sys.calcTimeOfNextScheduledEvent(integ.getAdvancedState(), 
            nextScheduledEvent, scheduledEventIds, true);

        switch(integ.stepTo(reportNo*hReport, nextScheduledEvent)) {
            case Integrator::ReachedStepLimit: printf("STEP LIMIT\n"); break;
            case Integrator::ReachedReportTime: printf("REPORT TIME AT t=%.17g\n", integ.getTime()); break;
            case Integrator::StartOfContinuousInterval: printf("START OF CONTINUOUS INTERVAL"); break;

            case Integrator::ReachedScheduledEvent:  {
                Stage lowestModified = Stage::Report;
                bool shouldTerminate;
                printf("SCHEDULED EVENT\n");
                sys.handleEvents(integ.updAdvancedState(),
                    System::ScheduledEvents,
                    scheduledEventIds,
                    integ.getAccuracyInUse(),
                    integ.getStateWeightsInUse(),
                    integ.getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
            }

            case Integrator::TimeHasAdvanced: {
                Stage lowestModified = Stage::Report;
                bool shouldTerminate;
                printf("TIME HAS ADVANCED TO %g\n", integ.getTime()); 
                sys.handleEvents(integ.updAdvancedState(),
                    System::TimeAdvancedEvent,
                    std::vector<EventId>(),
                    integ.getAccuracyInUse(),
                    integ.getStateWeightsInUse(),
                    integ.getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
                                              }

            case Integrator::ReachedEventTrigger: {
                Stage lowestModified = Stage::Report;
                bool shouldTerminate;
                printf("EVENT TRIGGERED AT tLow=%.17g tHigh=%.17g!!\n", 
                    integ.getTime(), integ.getAdvancedTime());
                cout << std::setprecision(17);
                cout << "Event window:     " << integ.getEventWindow() << endl;
//                cout << "Triggered events: " << integ.getTriggeredEvents();
                cout << "Transitions seen:";
                for (int i=0; i<(int)integ.getEventTransitionsSeen().size(); ++i)
                    cout << " " << EventStatus::eventTriggerString(integ.getEventTransitionsSeen()[i]);
                cout << endl;
//                cout << "Est event times:  " << integ.getEstimatedEventTimes();

                // state(t-) => state(t+)
                sys.handleEvents(integ.updAdvancedState(),
                    System::TriggeredEvents,
                    integ.getTriggeredEvents(),
                    integ.getAccuracyInUse(),
                    integ.getStateWeightsInUse(),
                    integ.getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
            }
                                                  
            case Integrator::EndOfSimulation: {
                Stage lowestModified = Stage::Report;
                bool shouldTerminate;

                printf("SIMULATION IS OVER. TERMINATION REASON=<TODO>\n");
                sys.handleEvents(integ.updAdvancedState(),
                    System::TerminationEvent,
                    std::vector<EventId>(),
                    integ.getAccuracyInUse(),
                    integ.getStateWeightsInUse(),
                    integ.getConstraintWeightsInUse(),
                    lowestModified, shouldTerminate);
                integ.reinitialize(lowestModified, shouldTerminate);
                break;
            }

            default: assert(!"Unrecognized return from stepTo()");
        }

        // fall through to here to report
        const State& s = integ.getState();

        sys.realize(s);
        printf(
            " -%6s- %9.6lf(%9.6lf) %14.10lf  %14.10lf  %14.10lf  %14.10lf | %14.10lf %14.10lf %14.10lf %14.10lf %14.10lf\n",
            integ.isStateInterpolated() ? "INTERP" : "------",
            s.getTime(), integ.getAdvancedTime(),
            s.getY()[0], s.getY()[1], s.getY()[2], s.getY()[3],
            s.getYErr()[0], s.getYErr()[1], 
            s.getEventsByStage(SubsystemIndex(0),Stage::Position)[0], s.getEventsByStage(SubsystemIndex(0),Stage::Position)[1], 
            s.getEventsByStage(SubsystemIndex(0),Stage::Position)[2]);

        cout << "YDot:        " << s.getYDot() << endl;
        cout << "Multipliers: " << s.getMultipliers() << endl;
        cout << "UDotErrs:    " << s.getUDotErr() << endl;
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
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  long int nproj, nce, nsetupsP, nprf;

  h0u=NaN;
  nst=nfe=nsetups=nje=nfeLS=nni=ncfn=netf=nge=-1;
  nproj=nce=nsetupsP=nprf=-1;

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
  nst   = integ.getNStepsTaken();
  nfe   = integ.getNRealizations();
  netf  = integ.getNErrorTestFailures();
  nproj = integ.getNProjections();
  nprf  = integ.getNProjectionFailures();

  printf("\nFinal Statistics:\n");
  printf("h0u = %g\n",h0u);
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld\n",
	 nst, nfe, nsetups);
  printf("nfeLS = %-6ld nje = %ld\n",
	 nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld \n",
	 nni, ncfn, netf);
  printf("nproj = %-6ld nce = %-6ld nsetupsP = %-6ld nprf = %-6ld\n",
         nproj, nce, nsetupsP, nprf);
  printf("nge = %ld\n", nge);

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
    event0 = s.allocateEvent(subsysIndex, Stage::Position, 3);
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
    
    s.updEventsByStage(subsysIndex, Stage::Position)[0] = 100*q[0]-q[1];

    s.updEventsByStage(subsysIndex, Stage::Position)[1] = 
        s.getTime() > 1.49552 && s.getTime() < 12.28937;

    s.updEventsByStage(subsysIndex, Stage::Position)[2] = s.getTime()-1.495508;
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

    const Real r2 = q[0]*q[0] + q[1]*q[1];
    const Real v2 = u[0]*u[0] + u[1]*u[1];
    const Real L  = (m*v2 - mg*q[1])/r2;
    udot[0] = - q[0]*L/m;
    udot[1] = - q[1]*L/m - g;
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

int MyPendulumGuts::projectImpl(State& s, Real consAccuracy,
                                const Vector& yweights, const Vector& ctols,
                                Vector& yerrest, System::ProjectOptions opts) const // yerrest is in/out
{
    const Vec2& wq = Vec2::getAs(&yweights[0]);
    const Vec2& wu = Vec2::getAs(&yweights[2]);
    const Real& tp = ctols[0]; // inverse tolerances 1/ti
    const Real& tv = ctols[1];

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    const Vec2& u = Vec2::getAs(&s.getU(subsysIndex)[0]);
    Real& ep = s.updQErr(subsysIndex)[0];
    Real& ev = s.updUErr(subsysIndex)[0];

    //cout << "BEFORE wperr=" << tp*ep << endl;

    Real wqchg;
    if (opts.hasAnyPositionOptions()) {
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
    }

    // Do velocity projection at new values of q
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

    // Now do error estimates.


    if (yerrest.size()) {
        Vec2& eq = Vec2::updAs(&yerrest[0]);
        Vec2& eu = Vec2::updAs(&yerrest[2]);

        // Recalc PW, PWInv:
        const Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
        const Row2 PW = Row2(tp*q[0]/wq[0], tp*q[1]/wq[1]);
        const Vec2 PWinv = Vec2(wq[0]*square(wq[1])*q[0], 
                                square(wq[0])*wq[1]*q[1]) / (tp*wqr2);

        Vec2 qperp = PWinv*(PW*eq);
        Vec2 uperp = VWinv*(VW*eu);

        //cout << "ERREST before=" << yerrest 
        //     << " wrms=" << wrms(yerrest,yweights) << endl;
        //cout << "PW*eq=" << PW*eq << " VW*eu=" << VW*eu << endl;
        eq -= qperp; eu -= uperp;

        //cout << "ERREST after=" << yerrest 
        //     << " wrms=" << wrms(yerrest,yweights) << endl;
        //cout << "PW*eq=" << PW*eq << " VW*eu=" << VW*eu << endl;
    }

    return 0;
}



