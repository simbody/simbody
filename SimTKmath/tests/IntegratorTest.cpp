/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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

/* This is a test program which uses the Integrator class to advance a 
hybrid continuous-discrete system. The system is a planar pendulum modeled
as a DAE, but without using any of Simbody's multibody formulation.
The discrete event handling here is just a random assortment of triggers
and actions to test the event detection/localization capabilities of 
Integrator. We do event handling explicitly here without using a TimeStepper,
based on return values from the Integrator. 

This system is a 2d pendulum swinging in gravity. It is modeled as a point 
mass free in the plane, plus a distance constraint to model the rod.

   y       | g               O
   ^       v                  \  d
   |                           \
   |                            * m
    ------> x

Gravity acts in the -y direction, the rod is length d, mass m, pivot location is
the ground origin (0,0). */

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

static void reportState(const char* msg, const Integrator& integ);
static void printFinalStats(const Integrator& integ);

//==============================================================================
//                                  EVENTS
//==============================================================================
// There are several built-in events in addition to those defined here.

// This event is caused by a variety of triggers. Its actions change the 
// subsequent trajectory.
class TimeToMakeAChangeEvent : public Event {
public:
    TimeToMakeAChangeEvent() : Event("Time to make a change") {}
private:
    TimeToMakeAChangeEvent* cloneVirtual() const override
    {   return new TimeToMakeAChangeEvent(*this); }
};

// This event occurs when it is time to make a generic report. Associated 
// actions won't change the simulation trajectory.
class ReportStateEvent : public Event {
public:
    ReportStateEvent() : Event("Reporting time reached") {}
private:
    ReportStateEvent* cloneVirtual() const override 
    {   return new ReportStateEvent(*this); }
};


//==============================================================================
//                                  TRIGGERS
//==============================================================================
// There are several built-in triggers in addition to those defined here, and
// timer triggers such as the ones we'll use in the example don't require a
// separate user-written class.

//--------------------------------- Q RATIO ------------------------------------
// Trigger when q1 = ratio*q0. Witness is 
//         w(q)=ratio*q0 - q1
//       dw(qd)=ratio*qdot0 - qdot1
//     ddw(qdd)=ratio*qdotdot0 - qdotdot1
class QRatio : public EventTrigger::Witness {
public:
    QRatio(SystemQIndex q0, SystemQIndex q1, Real ratio)
    :   EventTrigger::Witness("q ratio achieved"), 
        m_q0(q0), m_q1(q1), m_ratio(ratio) {
        setAccuracyRelativeTimeLocalizationWindow(1);
        setTriggerOnRisingSignTransition(false); // falling only
    }

    SystemQIndex getQ0() const {return m_q0;}
    SystemQIndex getQ1() const {return m_q1;}
    Real getRatio() const {return m_ratio;}

private:
    QRatio* cloneVirtual() const override {return new QRatio(*this);}

    Real calcWitnessValueVirtual(const System& system,
                                 const State&  state,
                                 int           derivOrder) const override
    {
        Real val = NaN;
        switch(derivOrder) {
        case 0:
            val = m_ratio*state.getQ()[m_q0] - state.getQ()[m_q1];
            break;
        case 1:
            val = m_ratio*state.getQDot()[m_q0] - state.getQDot()[m_q1];
            break;
        case 2:
            val = m_ratio*state.getQDotDot()[m_q0] - state.getQDotDot()[m_q1];
            break;
        default:
            SimTK_ASSERT1_ALWAYS(!"illegal derivOrder",
                "QRatio:: calcValueVirtual(): "
                "Deriv order %d out of range 0-2.", derivOrder);
        }
        return val;
    }

    Stage getDependsOnStageVirtual(int derivOrder) const override {
        Stage stage = Stage::Empty;
        switch(derivOrder) {
        case 0:
            stage = Stage::Time; // q's are state variables
            break;
        case 1:
            stage = Stage::Velocity; // qdots are calculated from states u
            break;
        case 2:
            stage = Stage::Acceleration;
            break;
        default:
            SimTK_ASSERT1_ALWAYS(!"illegal derivOrder",
                "QRatio::getDependsOnStageVirtual(): "
                "Deriv order %d out of range 0-2.", derivOrder);
        }
        return stage;
    }

    int getNumTimeDerivativesVirtual() const override 
    {   return 2; }


private:
    SystemQIndex    m_q0, m_q1;
    Real            m_ratio;
};

//---------------------------- SQUARE WAVE WITNESS -----------------------------
// This nonsmooth witness is a square wave that is -.5 until a start time,
// then +.5 until an end time, then back to -.5 forever. No derivatives are
// provided.
class SquareWaveWitness : public EventTrigger::Witness {
public:
    SquareWaveWitness(double tstart, double tend)
    :   EventTrigger::Witness("square wave transition"), 
        m_startTime(tstart), m_endTime(tend) {
        assert(tend > tstart);
        setWitnessIsDiscontinuous(true);
        // take defaults: trigger on both transitions; default window
    }
private:
    SquareWaveWitness* cloneVirtual() const override
    {   return new SquareWaveWitness(*this); }

    Real calcWitnessValueVirtual(const System& system,
                                 const State&  state,
                                 int           derivOrder) const override {
        assert(derivOrder==0); // value only
        const double now = state.getTime();
        return Real((m_startTime < now && now < m_endTime)-0.5);
    }

    Stage getDependsOnStageVirtual(int derivOrder) const override
    {   return Stage::Time; }

    const double m_startTime, m_endTime;
};

//--------------------------- TIME CROSSING WITNESS ----------------------------
// This will trigger when time reaches a set value. In real life this should be
// done using a Timer rather than a Witness.
class TimeCrossingWitness : public EventTrigger::Witness {
public:
    explicit TimeCrossingWitness(double t)
    :   EventTrigger::Witness("time crossing reached"), m_triggerTime(t) {
        setTriggerOnFallingSignTransition(false); // rising only
    }

private:
    TimeCrossingWitness* cloneVirtual() const override
    {   return new TimeCrossingWitness(*this); }

    Real calcWitnessValueVirtual(const System& system,
                                 const State&  state,
                                 int           derivOrder) const override {
        assert(derivOrder==0); // value only
        return state.getTime() - m_triggerTime;
    }

    Stage getDependsOnStageVirtual(int derivOrder) const override
    {   return Stage::Time; }

    const double m_triggerTime;
};

//==============================================================================
//                                  ACTIONS
//==============================================================================
//-------------------------- REPORT INITIALIZATION -----------------------------
class ReportInitializationAction : public EventAction {
public:
    ReportInitializationAction() : EventAction(Report) {}

private:
    ReportInitializationAction* cloneVirtual() const override
    {   return new ReportInitializationAction(*this); }

    void reportVirtual(const Study&            study,
                       const Event&            event,
                       const EventTriggers&    triggers) const override {
        printf("ReportInitializationAction::reportVirtual()\n");
        printf("***** INITIALIZATION REPORT *****\n");
        const System& system = study.getSystem();
        const State& state = study.getCurrentState();
        cout << "===> t=" << state.getTime() << ": REPORTING " 
             << " EVENT id=" << event.getEventId() << endl;
        cout << "  descr: '" << event.getEventDescription() << "'" << endl;

        for (auto tp : triggers) {
            printf("  trigger id=%d desc='%s' type=%s\n",
                (int)tp->getEventTriggerId(), tp->getTriggerDescription().c_str(),
                typeid(*tp).name());
        }
    }

};

//------------------------- REPORT EVENT OCCURRENCE ----------------------------
class ReportEventOccurrenceAction : public EventAction {
public:
    ReportEventOccurrenceAction() : EventAction(Report) {}

private:
    ReportEventOccurrenceAction* cloneVirtual() const override
    {   return new ReportEventOccurrenceAction(*this); }

    void reportVirtual(const Study&            study,
                       const Event&            event,
                       const EventTriggers&    triggers) const override
    {
        printf("ReportEventOccurrenceAction::reportVirtual()\n");
        const System& system = study.getSystem();
        const State& state = study.getCurrentState();
        cout << "===> t=" << state.getTime() << ": REPORTING " 
             << " EVENT id=" << event.getEventId() << endl;
        cout << "  descr: '" << event.getEventDescription() << "'" << endl;

        for (auto tp : triggers) {
            printf("  trigger id=%d desc='%s' type=%s\n",
                (int)tp->getEventTriggerId(), tp->getTriggerDescription().c_str(),
                typeid(*tp).name());
        }
    }

};


//------------------------------ REPORT STATE ----------------------------------
class ReportStateAction : public EventAction {
public:
    explicit ReportStateAction(const std::string& msg) 
    :   EventAction(Report), m_message(msg) {}

private:
    ReportStateAction* cloneVirtual() const override
    {   return new ReportStateAction(*this); }

    void reportVirtual(const Study&            study,
                       const Event&            event,
                       const EventTriggers&    triggers) const override
    {
        printf("ReportStateAction(%s)::reportVirtual()\n", m_message.c_str());
        reportState("Reporting from ReportStateAction", 
                    Integrator::downcast(study));
    }

    std::string m_message;
};


//--------------------------- TIME ADVANCED ACTION -----------------------------
class TimeAdvancedAction : public EventAction {
public:
    TimeAdvancedAction() : EventAction(Report) {}

private:
    TimeAdvancedAction* cloneVirtual() const override
    {   return new TimeAdvancedAction(*this); }

    void reportVirtual(const Study&            study,
                       const Event&            event,
                       const EventTriggers&    triggers) const override
    {
        printf("TimeAdvancedAction::reportVirtual()\n");
        printf("    TIME ADVANCED TO %g\n", study.getCurrentState().getTime());
        printf("    %d trigger id=%d\n", (int)triggers.size(),
               (int)triggers[0]->getEventTriggerId());
    }
};


//------------------------------ SWAP Qs ACTION --------------------------------
// Unconditionally swap two pre-specified q's when this action is performed.
class SwapQsAction : public EventAction {
public:
    SwapQsAction(SystemQIndex qix, SystemQIndex qjx) 
    :   EventAction(Change), m_qix(qix), m_qjx(qjx) {}

private:
    SwapQsAction* cloneVirtual() const override
    {   return new SwapQsAction(*this); }

    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override
    {
        printf("SwapQsAction::changeVirtual()\n");
        const System& system = study.getSystem();
        State& state = study.updInternalState();

        cout << "SwapQsAction CHANGE: Event '" 
             << event.getEventDescription()
             << "' id=" << event.getEventId() << endl;

        for (auto tp : triggers) {
            printf("  trigger id=%d desc='%s' type=%s\n",
                (int)tp->getEventTriggerId(), tp->getTriggerDescription().c_str(),
                typeid(*tp).name());
        }

        printf("Swapping q%d with q%d; zeroing u's\n", (int)m_qix, (int)m_qjx);
        std::swap(state.updQ()[m_qix],state.updQ()[m_qjx]); //invalidates Position stage
        state.updU() = 0;
        //result.reportExitStatus(EventChangeResult::Succeeded);
        //result.reportExitStatus(EventChangeResult::ShouldTerminate);
        //result.reportExitStatus(EventChangeResult::Failed,
        //                        "OMG I fucked up!!!");
    }

    SystemQIndex m_qix, m_qjx;
};


//==============================================================================
//                           MY PENDULUM GUTS
//==============================================================================
// User-defined system to be integrated. This is a kind of SimTK::System.
class MyPendulum;
class MyPendulumGuts: public System::Guts {
public:
    MyPendulumGuts() : Guts() {
        // Index types set themselves invalid on construction.
    }

    inline const MyPendulum& getMyPendulum() const;

    MyPendulumGuts* cloneImpl() const override 
    {   return new MyPendulumGuts(*this); }

        //////////////////////////////////////////////////
        // Implementation of continuous System virtuals //
        //////////////////////////////////////////////////

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
                      const ProjectOptions& options, ProjectResults& results) 
                      const override;
    void projectUImpl(State&, Vector& uErrEst, 
                      const ProjectOptions& options, ProjectResults& results) 
                      const override;

private:
friend class MyPendulum;

    // Called from realizeTopology() with cache temporarily writable.
    void realizeToTopologyCache(State& state);

    // TOPOLOGY STATE
    SubsystemIndex              m_subsysIndex;

    // TOPOLOGY CACHE
    DiscreteVariableIndex       m_massIndex, m_lengthIndex, m_gravityIndex;
    QIndex                      m_q0;
    UIndex                      m_u0;
    QErrIndex                   m_qerr0;
    UErrIndex                   m_uerr0;
    UDotErrIndex                m_udoterr0;
    CacheEntryIndex             m_mgForceIndex; // m*g, calc. at Dynamics stage
    EventId                     m_timeToMakeAChangeEvent, 
                                m_reportStateEvent;
};


//==============================================================================
//                                MY PENDULUM
//==============================================================================
// This is the handle class for a MyPendulum System. It must not have any data 
// members. Data, if needed, is in the corresponding "Guts" class.
class MyPendulum: public System {
public:

    // Constructor allocates needed resources from System base class.
    MyPendulum() : System(new MyPendulumGuts())
    { 
        MyPendulumGuts& guts = updGuts();
        guts.m_subsysIndex = getSystemGlobalSubsystem().getMySubsystemIndex();

        // Set up TimeToMakeAChangeEvent and its Triggers.
       {auto ttmc = new TimeToMakeAChangeEvent();
        ttmc->adoptEventAction(new ReportEventOccurrenceAction());
        ttmc->adoptEventAction(new SwapQsAction(SystemQIndex(0),
                                                SystemQIndex(1)));
        guts.m_timeToMakeAChangeEvent = adoptEvent(ttmc);}

       {auto e5123 = new EventTrigger::Timer::Periodic(5.123);
        e5123->setShouldTriggerAtZero(false);
        e5123->addEvent(guts.m_timeToMakeAChangeEvent); 
        adoptEventTrigger(e5123);}

       // These triggers are designed to occur simultaneously.
       {auto sqw = new SquareWaveWitness(/*1.49552*/1.49545, 12.28937);
        sqw->addEvent(guts.m_timeToMakeAChangeEvent);
        adoptEventTrigger(sqw);}

       {// Trigger when q1=100*q0.
        auto qr = new QRatio(SystemQIndex(0),SystemQIndex(1),100.);
        qr->addEvent(guts.m_timeToMakeAChangeEvent);
        adoptEventTrigger(qr);}

       {auto tcw = new TimeCrossingWitness(1.495508);
        tcw->addEvent(guts.m_timeToMakeAChangeEvent);
        adoptEventTrigger(tcw);}

        // Set up ReportStateEvent and its Triggers.
       {auto rs = new ReportStateEvent();
        rs->adoptEventAction(new ReportEventOccurrenceAction());
        rs->adoptEventAction(new ReportStateAction("for ReportStateEvent"));
        guts.m_reportStateEvent = adoptEvent(rs);}

       {auto report1 = new EventTrigger::Timer::Periodic(1.);
        report1->addEvent(guts.m_reportStateEvent); 
        adoptEventTrigger(report1);}

       {auto report2 = new EventTrigger::Timer::Designated
           ("some designated times", {.1, 3.14159, 1.5, 9.99, 23.456});
        report2->addEvent(guts.m_reportStateEvent);
        adoptEventTrigger(report2);}


        // Add initialization action.
        updInitializationEvent()
            .adoptEventAction(new ReportInitializationAction());
        updInitializationEvent()
            .adoptEventAction(new SwapQsAction(SystemQIndex(0),
                                               SystemQIndex(1)));

        // Add a time advanced action.
        updTimeAdvancedEvent().adoptEventAction(new TimeAdvancedAction());
        //setHasTimeAdvancedEvents(true); // TODO: shouldn't be needed now

        // Add termination action.
        updTerminationEvent()
            .adoptEventAction(new ReportStateAction("for TerminationEvent"));


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
        MyPendulumGuts& guts = updGuts();
        guts.updDefaultState()
            .updDiscreteVariable(guts.m_subsysIndex, 
                                 guts.m_massIndex) = Value<Real>(mass);
    }

    void setDefaultLength(Real length) {
        MyPendulumGuts& guts = updGuts();
        guts.updDefaultState()
            .updDiscreteVariable(guts.m_subsysIndex, 
                                 guts.m_lengthIndex) = Value<Real>(length);
    }

    void setDefaultGravity(Real gravity) {
        MyPendulumGuts& guts = updGuts();
        guts.updDefaultState()
            .updDiscreteVariable(guts.m_subsysIndex, 
                                 guts.m_gravityIndex) = Value<Real>(gravity);
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        MyPendulumGuts& guts = updGuts();
        guts.updDefaultState().updU(guts.m_subsysIndex) = u;
        guts.updDefaultState().updQ(guts.m_subsysIndex) = q;
        guts.updDefaultState().updTime() = t;
    }

    Real getMass(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& m = s.getDiscreteVariable(guts.m_subsysIndex, 
                                                       guts.m_massIndex);
        return Value<Real>::downcast(m).get();
    }
    Real getDefaultMass() const {return getMass(getDefaultState());}

    Real getLength(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& d = s.getDiscreteVariable(guts.m_subsysIndex, 
                                                       guts.m_lengthIndex);
        return Value<Real>::downcast(d).get();
    }
    Real getDefaultLength() const {return getLength(getDefaultState());}

    Real getGravity(const State& s) const {
        const MyPendulumGuts& guts = getGuts();
        const AbstractValue& g = s.getDiscreteVariable(guts.m_subsysIndex, 
                                                       guts.m_gravityIndex);
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
        cout << "  q=" << getDefaultState().getQ(guts.m_subsysIndex) << endl;
        cout << "  u=" << getDefaultState().getU(guts.m_subsysIndex) << endl;
    }

};

// Had to wait until static_cast can see that MyPendulum derives from System.
inline const MyPendulum& MyPendulumGuts::getMyPendulum() const {
    return static_cast<const MyPendulum&>(getSystem());
}


//==============================================================================
//                                    MAIN
//==============================================================================

// Perform event change actions, handle errors, reinitialize integrator.
static void performEventChangeActions
   (Integrator&         integ,
    EventsAndCauses&    triggeredEvents) 
{
    const System& sys = integ.getSystem();

    printf("----> PERFORM EVENT CHANGE ACTIONS @t=%.15g "
           "(%d triggers)\n", 
           integ.getAdvancedTime(), (int)triggeredEvents.size());

    EventChangeResult result;
    sys.performEventChangeActions(integ, triggeredEvents, result);
    const bool shouldTerminate = 
        result.getExitStatus() == EventChangeResult::ShouldTerminate;
    const Stage lowestModified = 
        result.getLowestModifiedStage();

    printf("----> CHANGE RESULT: status=%s, lowestModified=%s, "
           "shouldTerminate=%s\n",
           EventChangeResult::getExitStatusName(result.getExitStatus()),
           lowestModified.getName().c_str(),
           shouldTerminate?"true":"false");

    if (result.getExitStatus() == EventChangeResult::Failed) {
        printf("----> HANDLER FAILED, message='%s'\n",
               result.getMessage().c_str());
        integ.terminate(Integrator::AnUnrecoverableErrorOccurred);
        throw std::runtime_error(result.getMessage());
    }

    integ.reinitialize(lowestModified, shouldTerminate);
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

    Real timeOfLastScheduledReport = -Infinity, 
         timeOfLastScheduledChange = -Infinity;

    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;

    while (!integ.isSimulationOver()) {
        EventTriggers scheduledReportTimers, scheduledChangeTimers;
        Real timeOfNextReport, timeOfNextChange;
        
        // Clear the lists of triggered events (otherwise they will get
        // appended to by noteEventOccurrences()).
        triggeredEvents.clear(); ignoredEventIds.clear();

        // The same Timer may appear both on the report and change lists since
        // an Event may have both report and change actions. In that case, the
        // integrator will return ReachedReportTime first, at which point 
        // report actions should be performed, then on reentry will immediately
        // return ReachedScheduledEvent at the same time.
        // TimeAdvanced events behave similarly, but witness-triggered events
        // get reported *before* the event, then the state that actually caused
        // the event must not appear as part of the trajectory, so its report
        // actions must be invoked after its change actions have fixed whatever
        // is wrong.
        sys.findNextScheduledEventTimes
           (integ, timeOfLastScheduledReport, timeOfLastScheduledChange,
            timeOfNextReport, scheduledReportTimers,
            timeOfNextChange, scheduledChangeTimers);

        printf("----------------------------------------------------------\n");
        printf("stepTo(%g,%g)\n", timeOfNextReport, timeOfNextChange);
        switch(integ.stepTo(timeOfNextReport, timeOfNextChange)) {
        case Integrator::ReachedStepLimit: 
            printf("STEP LIMIT\n"); // not an event
            break;
        case Integrator::StartOfContinuousInterval: 
            printf("START OF CONTINUOUS INTERVAL\n"); // not an event
            break;

        case Integrator::ReachedReportTime: {
            printf("ReachedReportTime t=%.17g\n", integ.getTime()); 
            timeOfLastScheduledReport = integ.getTime();
            sys.noteEventOccurrence(scheduledReportTimers,
                                    triggeredEvents, ignoredEventIds);
            sys.performEventReportActions(integ, triggeredEvents);
            break;
        }

        case Integrator::ReachedScheduledEvent:  {
            printf("ReachedScheduledEvent t=%.17g\n", integ.getTime());
            reportState("BEFORE SCHEDULED EVENT:", integ);
            timeOfLastScheduledChange = integ.getTime();
            sys.noteEventOccurrence(scheduledChangeTimers,
                                    triggeredEvents, ignoredEventIds);
            performEventChangeActions(integ, triggeredEvents);
            break;
        }

        // See above for discussion: report, then change.
        case Integrator::TimeHasAdvanced: {
            printf("TIME HAS ADVANCED TO %g\n", integ.getTime());
            sys.noteEventOccurrence({&sys.getTimeAdvancedTrigger()},
                                    triggeredEvents, ignoredEventIds);
            sys.performEventReportActions(integ, triggeredEvents);

            performEventChangeActions(integ, triggeredEvents);
            break;
        }

        // Here we are looking at an invalid state that needs to be changed
        // before it can be part of the trajectory. Don't report it; just
        // perform change actions.
        case Integrator::ReachedEventTrigger: {
            printf("EVENT TRIGGERED AT tLow=%.17g tHigh=%.17g!!\n", 
                integ.getTime(), integ.getAdvancedTime());
            cout << std::setprecision(17);
            cout << "Event window:     " << integ.getEventWindow() << endl;
            cout << "Triggered witnesses: " 
                    << integ.getTriggeredWitnesses()<< "\n";
            cout << "Transitions seen:";
            for (auto tr : integ.getWitnessTransitionsSeen())
                cout << " " << Event::eventTriggerDirectionString(tr);
            cout << "\nEst event times:  " 
                    << integ.getEstimatedTriggerTimes() << "\n";
            reportState("BEFORE TRIGGERED EVENT:", integ);

            sys.noteEventOccurrence(integ.getTriggeredWitnesses(),
                                    triggeredEvents, ignoredEventIds);
            performEventChangeActions     // state(t-) => state(t+)
                (integ, triggeredEvents);
            break;
        }
            
        // Like Initialization events, Termination events are issued by
        // the Integrator rather than the TimeStepper.
        case Integrator::EndOfSimulation: {
            String reason = integ.getTerminationReasonString
                (integ.getTerminationReason());

            printf("SIMULATION IS OVER. TERMINATION REASON=%s\n",
                    reason.c_str());
            break;
        }

        default: assert(!"Unrecognized return from stepTo()");
        }

        // fall through to here to report
        reportState("fall through", integ);
    }

    printFinalStats(integ);

    return 0;
  } 
  catch (std::exception& e) {
    std::printf("FAILED: %s\n", e.what());
    return 1;
  }
}


//------------------------------------------------------------------------------
//                             REPORT STATE
//------------------------------------------------------------------------------
static void reportState(const char* msg, const Integrator& integ) {
    if (*msg) printf("%s\n", msg);

    const System& sys = integ.getSystem();
    const State&  s   = integ.getState();

    sys.realize(s);
    printf(
        " -%6s- %9.6lf(%9.6lf) %14.10lf  %14.10lf  %14.10lf  %14.10lf | %14.10lf %14.10lf\n",
        integ.isStateInterpolated() ? "INTERP" : "------",
        s.getTime(), integ.getAdvancedTime(),
        s.getY()[0], s.getY()[1], s.getY()[2], s.getY()[3],
        s.getYErr()[0], s.getYErr()[1]);

    Array_<const EventTrigger::Witness*,ActiveWitnessIndex> witnesses;
    sys.findActiveEventWitnesses(integ, witnesses);
    printf("%d witness values:\n", witnesses.size());
    for (auto w : witnesses) {
        printf("  %d: %g\n", (int)w->getEventTriggerId(), 
               w->calcWitnessValue(sys,s,0));
    }

    cout << "YDot:        " << s.getYDot() << endl;
    cout << "Multipliers: " << s.getMultipliers() << endl;
    cout << "UDotErrs:    " << s.getUDotErr() << endl;
}

//------------------------------------------------------------------------------
//                           PRINT FINAL STATS
//------------------------------------------------------------------------------
static void printFinalStats(const Integrator& integ)
{
  Real h0u;
  int nst, nattempt, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int nproj, nprojq, nproju, nce, nsetupsP, nprf, nprqf, npruf;

  h0u=NaN;
  nst=nattempt=nfe=nsetups=nje=nfeLS=nni=ncfn=netf=nge=-1;
  nproj=nprojq=nproju=nce=nsetupsP=nprf=nprqf=npruf=-1;

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


//------------------------------------------------------------------------------
//                        REALIZE TOPOLOGY IMPL
//------------------------------------------------------------------------------
// Finalize System construction; allocate State resources; last time to write
// on System topology cache.

void MyPendulumGuts::realizeToTopologyCache(State& s) {
    const Subsystem& sub = getSubsystem(m_subsysIndex);
    // Instance variables mass, length, gravity
    m_massIndex = sub.allocateDiscreteVariable(s, Stage::Instance,
                                                      new Value<Real>(1));
    m_lengthIndex = sub.allocateDiscreteVariable(s, Stage::Instance,
                                                      new Value<Real>(1));
    m_gravityIndex = sub.allocateDiscreteVariable(s, Stage::Instance,
                                    new Value<Real>(13.7503716373294544));
    const Vector init(2, Real(0));
    m_q0 = sub.allocateQ(s, init);
    m_u0 = sub.allocateU(s, init);

    m_mgForceIndex = sub.allocateCacheEntry(s, Stage::Dynamics,
                                                            new Value<Real>());

    m_qerr0 = sub.allocateQErr(s, 1);
    m_uerr0 = sub.allocateUErr(s, 1);
    m_udoterr0 = sub.allocateUDotErr(s, 1); // and multiplier
}

int MyPendulumGuts::realizeTopologyImpl(State& s) const {
    const_cast<MyPendulumGuts*>(this)->realizeToTopologyCache(s);
    return 0;
}


int MyPendulumGuts::realizeModelImpl(State& s) const {
    return 0;
}
int MyPendulumGuts::realizeInstanceImpl(const State& s) const {

    return 0;
}

//------------------------------------------------------------------------------
//                        REALIZE POSITION IMPL
//------------------------------------------------------------------------------
// Calculate constraint position errors perr().
int MyPendulumGuts::realizePositionImpl(const State& s) const {
    const Subsystem& sub = getSubsystem(m_subsysIndex);
    const Real    d = getMyPendulum().getLength(s);
    const Vector& q = s.getQ(m_subsysIndex);

    // This is the perr() equation.
    sub.updQErr(s)[0] = (q[0]*q[0] + q[1]*q[1] - d*d)/2;

    return 0;
}

//------------------------------------------------------------------------------
//                         REALIZE VELOCITY IMPL
//------------------------------------------------------------------------------
// Calculate qdot=N*u and constraint velocity errors verr().
int MyPendulumGuts::realizeVelocityImpl(const State& s) const {
    const Subsystem& sub = getSubsystem(m_subsysIndex);
    const Vector& q    = s.getQ(m_subsysIndex);
    const Vector& u    = s.getU(m_subsysIndex);
    Vector&       qdot = s.updQDot(m_subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];

    // This is the verr() equation.
    sub.updUErr(s)[0] = q[0]*u[0] + q[1]*u[1];
    return 0;
}

//------------------------------------------------------------------------------
//                         REALIZE DYNAMICS IMPL
//------------------------------------------------------------------------------
// Calculate forces, in this case just gravity. The result goes into the
// State cache.
int MyPendulumGuts::realizeDynamicsImpl(const State& s) const {
    const Subsystem& sub = getSubsystem(m_subsysIndex);
    const Real m  = getMyPendulum().getMass(s);
    const Real g  = getMyPendulum().getGravity(s);

    Real& mg = Value<Real>::updDowncast
                                (sub.updCacheEntry(s, m_mgForceIndex)).upd();

    // Calculate the force due to gravity and store it in the cache.
    mg = m*g;
    return 0;
}

//------------------------------------------------------------------------------
//                        REALIZE ACCELERATION IMPL
//------------------------------------------------------------------------------
/* The DAE for a generic multibody system is:
        qdot = N u
      M udot = f - ~G lambda
      G udot = b
      perr(t,q) = 0
      verr(t,q,u) = 0

Let   r^2 = x^2  + y^2
      v^2 = x'^2 + y'^2
We will express the "rod length=d" constraint as 
      (r^2 - d^2)/2 = 0    (perr)
          xx' + yy' = 0    (verr)
        xx'' + yy'' = -v^2 (aerr)

So the matrix G = d perr/dq = [x y] and b = -v^2, and the equations of 
motion are:
    [ m 0 x ] [ x'' ]   [  0  ]
    [ 0 m y ] [ y'' ] = [ -mg ]
    [ x y 0 ] [ L   ]   [-v^2 ]
where L (the Lagrange multiplier) is proportional to the rod tension. You can 
solve this to get
    L   = (m*v^2 - mg*y)/(r^2)
    x'' = - x*L/m
    y'' = - y*L/m - g
*/ 
int MyPendulumGuts::realizeAccelerationImpl(const State& s) const {
    const Subsystem& sub = getSubsystem(m_subsysIndex);
    const Real m  = getMyPendulum().getMass(s);
    const Real g  = getMyPendulum().getGravity(s);

    // We're pretending we couldn't calculate this here!
    const Real mg = Value<Real>::downcast
                                (sub.updCacheEntry(s, m_mgForceIndex)).get();

    const Vector& q    = sub.getQ(s);
    const Vector& u    = sub.getU(s);
    Vector&       udot = sub.updUDot(s);
    Vector&       qdotdot = sub.updQDotDot(s);

    const Real r2 = q[0]*q[0] + q[1]*q[1];
    const Real v2 = u[0]*u[0] + u[1]*u[1];
    const Real L  = (m*v2 - mg*q[1])/r2;
    udot[0] = - q[0]*L/m;
    udot[1] = - q[1]*L/m - g;
    qdotdot = udot; // N=identity for this problem
    sub.updMultipliers(s)[0] = L;

    // This is the aerr() equation.
    sub.updUDotErr(s)[0] = q[0]*udot[0] + q[1]*udot[1] + v2;

    return 0;
}

//------------------------------------------------------------------------------
//                            PROJECT Q IMPL
//------------------------------------------------------------------------------
/* Here we want to remove any constraint errors from the current state,
and project out any component of the integrator's error estimate
perpendicular to the constraint manifold. We will do this sequentially
rather than handling position and velocity simultaneously.

For this system we have P = d perr/dq = V = d verr/du = [x y].
Weighted, we have PW=tp*[x/wx y/wy] VW=tv*[x/wxd y/wyd]. 
With pinv(A)=~A*(A*~A)^-1, we have:

   pinv(P)  = ~[            x             y] /  (    x ^2+     y ^2)
   pinv(PW) = ~(1/tp)*[(wx *wy ^2)*x (wx ^2*wy) *y] / ((wy *x)^2+(wx *y)^2)
   pinv(VW) = ~(1/tv)*[(wxd*wyd^2)*x (wxd^2*wyd)*y] / ((wyd*x)^2+(wxd*y)^2)
     (the latter assuming x,y already projected on position manifold)

We want to solve
   |perr(q0 - dq)|_TRMS <= accuracy, such that dq=min_WLS(dq)
   PW(q0) dq = Tp * perr(q0); q = q0-dq
Then
   |verr(q,u0 - du)|_TRMS <= accuracy, du=min_WLS(du)
   VW(q) du = Tv * verr(q,u0); u = u0-du


To remove the corresponding error estimates:
   PW(q) qperp = PW(q) qerrest; qerrest -= qperp
   VW(q) uperp = VW(q) uerrest; uerrest -= uperp
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
    const Vector& uweights = s.getUWeights(m_subsysIndex);
    const Vector& ctols = s.getQErrWeights(m_subsysIndex);
    // Since qdot=u here we can use uweights directly as qweights.
    const Vec2& wq = Vec2::getAs(&uweights[0]);
    const Real& tp = ctols[0]; // inverse tolerances 1/ti

    const Vec2& q = Vec2::getAs(&s.getQ(m_subsysIndex)[0]); // set up aliases
    Real& ep = s.updQErr(m_subsysIndex)[0]; // ep changes as we go

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
    
        s.updQ(m_subsysIndex)[0] -= wdq[0]/wq[0]; 
        s.updQ(m_subsysIndex)[1] -= wdq[1]/wq[1]; 
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

//------------------------------------------------------------------------------
//                            PROJECT U IMPL
//------------------------------------------------------------------------------
/* Velocity projection -- see comment above projectQImpl(). */
void MyPendulumGuts::projectUImpl(State& s, Vector& uerrest, 
             const ProjectOptions& opts, ProjectResults& results) const
{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Real projLimit = opts.getProjectionLimit();
    const bool forceProj = opts.isOptionSet(ProjectOptions::ForceProjection);
    const Vector& uweights = s.getUWeights(m_subsysIndex);
    const Vector& ctols = s.getUErrWeights(m_subsysIndex);

    const Vec2& wu = Vec2::getAs(&uweights[0]);
    const Real& tv = ctols[0];

    const Vec2& q = Vec2::getAs(&s.getQ(m_subsysIndex)[0]); // set up aliases
    const Vec2& u = Vec2::getAs(&s.getU(m_subsysIndex)[0]);
    Real& ev = s.updUErr(m_subsysIndex)[0]; // ev changes as we go

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

    s.updU(m_subsysIndex)[0] -= wdu[0]/wu[0]; 
    s.updU(m_subsysIndex)[1] -= wdu[1]/wu[1];

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



