/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-15 Stanford University and the Authors.        *
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

/* Here we'll build a very simple System containing a simple
Subsystem, and integrate with a simple integrator and manual event handling. 
This avoids all the usual Simmath and Simbody trappings and lets us just check 
out the underlying simulation architecture. */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include "SimTKcommon/internal/Event.h"
#include "SimTKcommon/internal/EventTrigger.h"
#include "SimTKcommon/internal/EventTimer.h"
#include "SimTKcommon/internal/EventWitness.h"

#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

//==============================================================================
//                             TEST SYSTEM
//==============================================================================
class TestSystemGuts : public System::Guts {
public:
    // implementations of System::Guts virtuals
    TestSystemGuts* cloneImpl() const override
    {   return new TestSystemGuts(*this); }

    bool prescribeQImpl(State& state) const override
    {
        return false;
    }
    bool prescribeUImpl(State& state) const override
    {
        return false;
    }
    void projectQImpl(State& state, Vector& yerrest, const ProjectOptions& opts, 
                     ProjectResults& result) const override
    {   
        result.setExitStatus(ProjectResults::Succeeded);
    }
    void projectUImpl(State& state, Vector& yerrest, const ProjectOptions& opts, 
                     ProjectResults& result) const override
    {
        result.setExitStatus(ProjectResults::Succeeded);
    }
};

class TestSystem : public System {
public:
    TestSystem() : System(new TestSystemGuts()) {}

    static const TestSystem& getAs(const System& sys)
    {   assert(dynamic_cast<const Guts*>(&sys.getSystemGuts()));
        return static_cast<const TestSystem&>(sys); }
    static TestSystem& updAs(System& sys)
    {   assert(dynamic_cast<const Guts*>(&sys.getSystemGuts()));
        return static_cast<TestSystem&>(sys); }
private:
    const TestSystemGuts& getGuts() const
    {   return dynamic_cast<const TestSystemGuts&>(getSystemGuts());}
    TestSystemGuts& updGuts()
    {   return dynamic_cast<TestSystemGuts&>(updSystemGuts());}
};


//==============================================================================
//                             TEST SUBSYSTEM
//==============================================================================

class TestSubsystemGuts : public Subsystem::Guts {
    struct StateVars {
        QIndex myQs;
        UIndex myUs;
    };
    struct CacheEntries {
        CacheEntryIndex qSumCacheIx, uSumCacheIx;
    };
    friend std::ostream& operator<<(std::ostream& o, const CacheEntries&);
    friend std::ostream& operator<<(std::ostream& o, const StateVars&);
public:
    TestSubsystemGuts() {}

    const Vec3& getQ3(const State& s) const 
    {   return Vec3::getAs(&getQ(s)[getStateVars(s).myQs]); }
    const Vec3& getU3(const State& s) const 
    {   return Vec3::getAs(&getU(s)[getStateVars(s).myQs]); }
    const Vec3& getQDot3(const State& s) const 
    {   return Vec3::getAs(&getQDot(s)[getStateVars(s).myQs]); }
    const Vec3& getUDot3(const State& s) const 
    {   return Vec3::getAs(&getUDot(s)[getStateVars(s).myUs]); }
    const Vec3& getQDotDot3(const State& s) const 
    {   return Vec3::getAs(&getQDotDot(s)[getStateVars(s).myQs]); }
    Real getQSum(const State& s) const 
    {   return Value<Real>::downcast
                            (getCacheEntry(s,getCacheEntries(s).qSumCacheIx)); }
    Real getUSum(const State& s) const 
    {   return Value<Real>::downcast(
                            getCacheEntry(s,getCacheEntries(s).uSumCacheIx)); }

    Vec3& updQ3(State& s) const 
    {return Vec3::updAs(&updQ(s)[getStateVars(s).myQs]);}
    Vec3& updU3(State& s) const 
    {return Vec3::updAs(&updU(s)[getStateVars(s).myUs]);}

private:
friend class TestSubsystem;
    // implementations of Subsystem::Guts virtuals

    TestSubsystemGuts* cloneImpl() const override
    {   return new TestSubsystemGuts(*this); }


    int realizeSubsystemTopologyImpl(State& s) const override {
        auto mThis = const_cast<TestSubsystemGuts*>(this);
        mThis->myStateVars = 
            allocateCacheEntry(s, Stage::Model, new Value<StateVars>());
        mThis->myCacheEntries = 
            allocateCacheEntry(s, Stage::Instance, new Value<CacheEntries>());
        return 0;
    }


    int realizeSubsystemModelImpl(State& s) const override {
        StateVars& vars = updStateVars(s);
        vars.myQs = allocateQ(s, Vector(Vec3(0)));
        vars.myUs = allocateU(s, Vector(Vec3(0)));
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const override {
        CacheEntries& cache = updCacheEntries(s);
        cache.qSumCacheIx = allocateCacheEntry(s, Stage::Position, 
                                               new Value<Real>(0));
        cache.uSumCacheIx = allocateCacheEntry(s, Stage::Velocity, 
                                               new Value<Real>(0));
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const override {
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const override {
        updQSum(s) = sum(getQ3(s));
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const override {
        updQDot3(s) = getU3(s);
        updUSum(s) = sum(getU3(s));
        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const override {
        updQDotDot3(s) = updUDot3(s) = Vec3(1,2,3);
        return 0;
    }


    Vec3& updQDot3(const State& s) const 
    {   return Vec3::updAs(&updQDot(s)[getStateVars(s).myQs]); }
    Vec3& updUDot3(const State& s) const 
    {   return Vec3::updAs(&updUDot(s)[getStateVars(s).myUs]); }
    Vec3& updQDotDot3(const State& s) const 
    {   return Vec3::updAs(&updQDotDot(s)[getStateVars(s).myQs]); }
    Real& updQSum(const State& s) const 
    {   return Value<Real>::updDowncast
                            (updCacheEntry(s,getCacheEntries(s).qSumCacheIx)); }
    Real& updUSum(const State& s) const 
    {   return Value<Real>::updDowncast
                            (updCacheEntry(s,getCacheEntries(s).uSumCacheIx)); }
    
    const StateVars& getStateVars(const State& s) const
    {   assert(myStateVars.isValid());
        return Value<StateVars>::downcast(getCacheEntry(s,myStateVars)); }
    const CacheEntries& getCacheEntries(const State& s) const
    {   assert(myCacheEntries.isValid());
        return Value<CacheEntries>::downcast(getCacheEntry(s,myCacheEntries)); }
    StateVars& updStateVars   (const State& s) const 
    {   assert(myStateVars.isValid());
        return Value<StateVars>::updDowncast(updCacheEntry(s,myStateVars)); }
    CacheEntries& updCacheEntries(const State& s) const
    {   assert(myCacheEntries.isValid());
        return Value<CacheEntries>::updDowncast(updCacheEntry(s,myCacheEntries)); }

    const TestSystem& getTestSystem() const {return TestSystem::getAs(getSystem());}
    TestSystem& updTestSystem() {return TestSystem::updAs(updSystem());}

        // TOPOLOGY STATE VARIABLES //
    EventId     m_somethingNeedsToHappenId;

        // TOPOLOGY CACHE //
    CacheEntryIndex myStateVars;
    CacheEntryIndex myCacheEntries;

};

// This Subsystem has 3 q's and 3 u's of its own, as well as whatever State
// variables its Measures require.
class TestSubsystem : public Subsystem {
public:
    class SomethingNeedsToHappen;
    class HandleEvent;
    class USumWitness;
    class TimeReachedWitness;

    TestSubsystem(System& sys);

    const SomethingNeedsToHappen& getMyEvent() const;

    Real getQSum(const State& s) const {return getGuts().getQSum(s);}
    Real getUSum(const State& s) const {return getGuts().getUSum(s);}

    static const TestSubsystem& downcast(const Subsystem& subsys) {
        assert(dynamic_cast<const TestSubsystemGuts*>
                                                (&subsys.getSubsystemGuts()));
        return static_cast<const TestSubsystem&>(subsys);
    }
private:
    const TestSubsystemGuts& getGuts() const
    {   return dynamic_cast<const TestSubsystemGuts&>(getSubsystemGuts());}
    TestSubsystemGuts& updGuts()
    {   return dynamic_cast<TestSubsystemGuts&>(updSubsystemGuts());}
};


// This event belongs to TestSubsystem. It is caused by several different 
// triggers. 
class TestSubsystem::SomethingNeedsToHappen : public Event {
public:
    SomethingNeedsToHappen() : Event("Something needs to happen") {}
private:
    SomethingNeedsToHappen* cloneVirtual() const override
    {   return new SomethingNeedsToHappen(*this); }
};


auto TestSubsystem::getMyEvent() const -> const SomethingNeedsToHappen& {
    return dynamic_cast<const SomethingNeedsToHappen&>
        (getSystem().getEvent(getGuts().m_somethingNeedsToHappenId));
}

// This is the Action we'll take whenever the SomethingNeedsToHappen event
// occurs. We just report the arguments we got, and then trash positions.
class TestSubsystem::HandleEvent : public EventAction {
public:
    HandleEvent(SubsystemIndex ix) 
    :   EventAction(Change), m_subsysIx(ix) {}

private:
    HandleEvent* cloneVirtual() const override
    {   return new HandleEvent(*this); }

    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override
    {
        auto& sys = study.getSystem();
        auto& s = study.updInternalState();

        // Make sure we are owned by a TestSubsystem.
        auto& testSub = TestSubsystem::downcast(sys.getSubsystem(m_subsysIx));

        cout << "**** TestSubsystem::HandleEvent t=" << s.getTime() 
             << " acc=" << study.getAccuracyInUse()
             << " tol=" << study.getConstraintToleranceInUse();
        cout << " Event id=" << event.getEventId()
             << " desc='" << event.getEventDescription() << "'\n";

        for (auto tp : triggers) {
            printf("  trigger id=%d desc='%s' type=%s\n",
                (int)tp->getEventTriggerId(), tp->getTriggerDescription().c_str(),
                typeid(*tp).name());
        }
        cout << "****" << endl;

        // Pretend we changed a position to test lowestModifiedStage
        // calculation. Try to hide our duplicity by realizing it again.
        s.invalidateAllCacheAtOrAbove(Stage::Position); 
        sys.realize(s, Stage::Velocity);
        result.reportExitStatus(EventChangeResult::Succeeded);
    }

    const SubsystemIndex m_subsysIx;
};

//------------------------------ U SUM WITNESS ---------------------------------
// Trigger when the sum of certain generalized speeds reaches a particular
// value. (This is treated as a Velocity-stage witness.)
class TestSubsystem::USumWitness : public EventWitness {
public:
    USumWitness(SubsystemIndex ix, Real targetSum)
    :   EventWitness("u sum witness", Bilateral, RisingAndFalling, Continuous),
        m_subsysIx(ix), m_targetSum(targetSum) {
        // default allows triggering both falling and rising transition
    }

private:
    USumWitness* cloneVirtual() const override {return new USumWitness(*this);}

    Value calcWitnessValueVirtual(const Study&  study,
                                  const State&  state,
                                  int           derivOrder) const override
    {
        assert(derivOrder==0);
        const System& system = study.getSystem();
        auto& testsub = TestSubsystem::downcast(system.getSubsystem(m_subsysIx));
        return Value(testsub.getUSum(state) - m_targetSum,
                     SignificantReal);
    }

    Stage getDependsOnStageVirtual(int derivOrder) const override {
        assert(derivOrder==0);
        return Stage::Velocity;
    }

    int getNumTimeDerivativesVirtual() const override 
    {   return 0; }


private:
    const SubsystemIndex    m_subsysIx;
    Real                    m_targetSum;
};

//--------------------------- TIME CROSSING WITNESS ----------------------------
// This will trigger when time reaches a set value. In real life this should be
// done using a Timer rather than a Witness.
class TestSubsystem::TimeReachedWitness : public EventWitness {
public:
    explicit TimeReachedWitness(double t)
    :   EventWitness("designated time reached", Bilateral, Rising, Continuous), 
        m_triggerTime(t) {}

private:
    TimeReachedWitness* cloneVirtual() const override
    {   return new TimeReachedWitness(*this); }

    Value calcWitnessValueVirtual(const Study&  study,
                                  const State&  state,
                                  int           derivOrder) const override {
        assert(derivOrder==0); // value only
        return Value(state.getTime() - m_triggerTime, SignificantReal);
    }

    Stage getDependsOnStageVirtual(int derivOrder) const override
    {   return Stage::Time; }

    const double m_triggerTime;
};

TestSubsystem::TestSubsystem(System& sys) {
    const Real TriggerTime1 = .6789, TriggerTime2 = 1.234;
    const Real TriggerUSum = 5;

    adoptSubsystemGuts(new TestSubsystemGuts());
    sys.adoptSubsystem(*this);

    // Create event and register with the System.
    {auto eventp = new SomethingNeedsToHappen();
        eventp->adoptEventAction(new HandleEvent(this->getMySubsystemIndex()));
        updGuts().m_somethingNeedsToHappenId = sys.adoptEvent(eventp);}

    // Register three witnesses that cause the above event.
    {auto time1p = new TimeReachedWitness(TriggerTime1);
        time1p->addEvent(getGuts().m_somethingNeedsToHappenId);
        sys.adoptEventTrigger(time1p);}

    {auto time2p = new TimeReachedWitness(TriggerTime2);
        time2p->addEvent(getGuts().m_somethingNeedsToHappenId);
        sys.adoptEventTrigger(time2p);}

    {auto sump = new USumWitness(this->getMySubsystemIndex(), TriggerUSum);
        sump->addEvent(getGuts().m_somethingNeedsToHappenId);
        sys.adoptEventTrigger(sump);}
}

//==============================================================================
//                         MY SIN COS MEASURE
//==============================================================================
// Computes sin(t) and cos(t) as a Vec2.
template <class T>
class MySinCos : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(MySinCos, Measure_<T>);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(MySinCos, Measure_<T>);
};

template <class T>
class MySinCos<T>::Implementation : public Measure_<T>::Implementation {
public:
    Implementation() 
    :   Measure_<T>::Implementation(T(Vec2(0)), 1) {}

    // Default copy constructor, destructor, copy assignment are fine.

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const override
    {   return new Implementation(*this); }
    int getNumTimeDerivativesVirtual() const override {return 0;}
    Stage getDependsOnStageVirtual(int order) const override
    {   return Stage::Time; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
        override
    {
        SimTK_ASSERT1_ALWAYS(derivOrder==0,
            "MySinCos::Implementation::calcCachedValueVirtual():"
            " derivOrder %d seen but only 0 allowed.", derivOrder);

        value[0] = std::sin(s.getTime());
        value[1] = std::cos(s.getTime());
    }
};


//==============================================================================
//                           FIXED TIME STEPPER
//==============================================================================
// Our time stepper must derive from Study so that we can satisfy the 
// signature for performing event actions.
class FixedTimeStepper : public Study {
public:
    explicit FixedTimeStepper(System& system)
    :   m_system(system), m_consTol(NaN) {}

    // Copy in the given state to this study's internal state and set the
    // constraint tolerance to be used when projecting.
    void initialize(const State& initState, Real consTol);

    // Take a fixed-size, explicit midpoint step of length h, advancing the 
    // internal state. Return true if we should terminate.
    bool step(double h);

    // Find all the triggering event witnesses that depend only on stage.
    void findEvents(Stage stage, EventTriggers& triggered) const;

    // Handle all the events that were triggered by witnesses at this stage.
    // This may modify the internal state. Return true if an event handler
    // says we should terminate the simulation.
    bool handleEvents(Stage stage, const EventTriggers& triggered);

    // Calcuate the values of all witnesses currently in m_witnesses. This 
    // state must already be realized to Stage::Acceleration.
    void calcWitnessValues(const State&                 state, 
                           Array_<EventWitness::Value>& values) {
        values.resize(m_witnesses.size());

        for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
            const EventWitness& w = *m_witnesses[awx];
            values[awx] = w.calcWitnessValue(*this, state);
        }
    }
private:
    System&     m_system;
    State       m_state;
    Real        m_consTol;

    Array_<const EventWitness*, 
           ActiveWitnessIndex>      m_witnesses;
    Array_<EventWitness::Value>     m_witnessValuesPrev; // same length

    // Implement the Study interface.
    const System& getSystemVirtual() const override
    {   return m_system; }
    const State& getCurrentStateVirtual() const override
    {   return m_state; } 
    const State& getInternalStateVirtual() const override
    {   return m_state; }  
    State& updInternalStateVirtual() override
    {   return m_state; } 
    Real getAccuracyInUseVirtual() const override
    {   return NaN; }
    Real getConstraintToleranceInUseVirtual() const override
    {   return m_consTol; }
};

//==============================================================================
//                                TEST ONE
//==============================================================================
void testOne() {
    TestSystem sys;
    TestSubsystem subsys(sys);

    // Add a Result measure to the system subsystem. This depends on 
    // Position stage and invalidates Dynamics and later stages.
    Measure_<Vector>::Result vectorResult(sys.updSystemGlobalSubsystem(), 
        Stage::Position, Stage::Dynamics);

    MeasureIndex vectorResultIx = vectorResult.getSubsystemMeasureIndex();
    cout << "vectorResult index=" 
         << vectorResultIx << endl;

    Measure_<Vector>::Result myVecRes =  Measure_<Vector>::Result::getAs(
        sys.updSystemGlobalSubsystem().getMeasure(vectorResultIx));

    Measure::Result result(sys.updSystemGlobalSubsystem(), 
        Stage::Time, Stage::Position);
    Measure::Result autoResult(sys.updSystemGlobalSubsystem(), 
        Stage::Time, Stage::Position);
    autoResult.setIsPresumedValidAtDependsOnStage(true);

    Measure::Zero zero(subsys);
    Measure::Constant three(subsys, 3);
    Measure_<Vec3>::Constant v3const(subsys, Vec3(1,2,3));
    Measure::Sinusoid cos2pit(subsys, 1, 2*Pi, Pi/2);

    // Integrate the cos(2pi*t) measure with IC=0; should give sin(2pi*t)/2pi.
    Measure::Integrate sin2pitOver2pi(subsys, cos2pit, zero);

    // These two compute -cos(t), sin(t) by integrating sin(t), cos(t) with
    // initial conditions -1,0.
    Measure_<Vec2>::Constant cossinInit(subsys, Vec2(-1,0));
    MySinCos<Vec2> mysincos(subsys);
    Measure_<Vec2>::Integrate cossin(subsys, mysincos, cossinInit);

    Measure_<Vector>::Constant vcossinInit(subsys, Vector(Vec2(-1,0)));
    MySinCos<Vector> vmysincos(subsys);
    Measure_<Vector>::Integrate vcossin(subsys, vmysincos, vcossinInit,
                                        Vector(2,Zero));
    Measure_<Vector>::Delay vcossin_delaypt1(subsys, vcossin, .1);

    Measure_<Real>::Minimum minCos2pit(subsys, cos2pit);
    Measure_<Real>::Maximum maxCos2pit(subsys, cos2pit);
    Measure_<Real>::MinAbs minAbsCos2pit(subsys, cos2pit);
    Measure_<Real>::MaxAbs maxAbsCos2pit(subsys, cos2pit);

    Measure::Differentiate dInteg(subsys, sin2pitOver2pi);
    dInteg.setForceUseApproximation(true);

    Measure::Time tMeasure;
    Measure::Time tSubMeas(subsys);

    Measure::Delay tDelayed(subsys, tSubMeas, 0.01);

    Measure::Scale t1000(subsys, 1000, tSubMeas);

    Measure::Variable mv(subsys, Stage::Position, 29);
    cout << "mv def value=" << mv.getDefaultValue() << endl;
    mv.setDefaultValue(-19);
    cout << "mv def value now=" << mv.getDefaultValue() << endl;

    cout << "Measure::Zero=" << Measure::Zero().getValue(State()) << endl;
    cout << "Measure::One=" << Measure::One().getValue(State()) << endl;

    Measure::Plus vplus(subsys, mv, cos2pit);
    Measure::Minus vminus(subsys, mv, cos2pit);

    Measure::Plus vplus2;
    vplus2.deepAssign(vplus);

    Measure m; // just a handle for referencing other measures
    m = cos2pit;



    cout << "vplus ref count=" << vplus.getRefCount() << endl;
    cout << "vplus2 ref count=" << vplus2.getRefCount() << endl;


    State state = sys.realizeTopology();
    cout << "sys topo version=" << sys.getSystemTopologyCacheVersion() << "\n";
    cout << "state topo version=" << state.getSystemTopologyStageVersion() << "\n";
    sys.invalidateSystemTopologyCache();
    sys.realizeTopology();
    cout << "sys topo version=" << sys.getSystemTopologyCacheVersion() << "\n";
    // Use sneaky loophole since we know state is still good.
    state.setSystemTopologyStageVersion(sys.getSystemTopologyCacheVersion());
    sys.realizeModel(state);

    m = zero;
    cout << "m=" << m.getValue(state) << endl;


    cout << "uWeights=" << state.getUWeights() << "\n";
    cout << "zWeights=" << state.getZWeights() << "\n";

    state.updUWeights()[1] = 9;
    state.updZWeights() = 21;

    cout << "uWeights=" << state.getUWeights() << "\n";
    cout << "zWeights=" << state.getZWeights() << "\n";

    sys.realize(state,Stage::Instance);

    cout << "qerrWeights=" << state.getQErrWeights() << "\n";
    cout << "uerrWeights=" << state.getUErrWeights() << "\n";

    State dupState = state;
    cout << "dup uWeights=" << state.getUWeights() << "\n";
    cout << "dup zWeights=" << state.getZWeights() << "\n";
    cout << "dup qerrWeights=" << state.getQErrWeights() << "\n";
    cout << "dup uerrWeights=" << state.getUErrWeights() << "\n";

    // Allocate vectorResult and initialize it. (Can't mark it valid yet.)
    vectorResult.updValue(state).resize(3);
    vectorResult.updValue(state) = Vector(Vec3(1,2,3));
    
    result.updValue(state) = 1.234;
    autoResult.updValue(state) = 4.321;


    state.setTime(1.234);
    sys.realize(state, Stage::Time);
    cout << "Initially, tMeasure=" << tMeasure.getValue(state)
         << " tSubMeas=" << tSubMeas.getValue(state) 
         << " 1000*tMeasure=" << t1000.getValue(state)
         << endl;

    Measure_<Mat22>::One m22Ident;
    cout << "Measure_<Mat22>::One=" << m22Ident.getValue(state) << endl;

    cout << "mv after realizeTopo=" << mv.getValue(state) << endl;

    cout << "cossinInit=" << cossinInit.getValue(state) << endl;


    State s2,s3;
    s2 = state; // new copies of variables
    s2 = state; // should do only assignments w/o heap allocation

    const Real h = .001;
    const int nSteps = 2000;
    const int outputInterval = 100;

    state.setTime(0);
    SimTK_TEST(state.getTime()==0);

    sys.realize(state, Stage::Position);
    cout << "mv+cos2pit=" << vplus.getValue(state) << endl;
    cout << "mv-cos2pit=" << vminus.getValue(state) << endl;

    cout << "Sys stage after realize(Pos):" 
         << state.getSystemStage().getName() << endl;

    mv.setValue(state, 1.234);

    cout << "Sys stage after mv=1.234:" 
         << state.getSystemStage().getName() << endl;

    cout << "mv is now=" << mv.getValue(state) << endl;

    sys.realize(state, Stage::Position);

    cout << "Realized Position:\n";
    vectorResult.markAsValid(state);
    cout << "vectorResult=" << vectorResult.getValue(state) << endl;
    cout << "myVecRes=" << myVecRes.getValue(state) << endl;
    result.markAsValid(state);
    cout << "result=" << result.getValue(state) << endl;
    // Shouldn't need to mark this one.
    cout << "autoResult=" << autoResult.getValue(state) << endl;

    sys.realize(state, Stage::Acceleration);

    cout << "Now stage=" << state.getSystemStage() << endl;
    vectorResult.setValue(state, Vector(5,9));
    cout << "After vectorResult.setValue(), vectorResult="
         << vectorResult.getValue(state) << endl;
    cout << "... but stage=" << state.getSystemStage() << endl;

    cossin.setValue(state, cossinInit.getValue(state));
    vcossin.setValue(state, vcossinInit.getValue(state));

    FixedTimeStepper ts(sys);
    ts.initialize(state, 1e-6); // initial state, constraint tolerance

    for (int i=0; i <= nSteps; ++i) {
        const State& state = ts.getInternalState();
        const double t = state.getTime();

        if (i % outputInterval == 0) {
            sys.realize(state, Stage::Report);
            cout << "\ntMeasure=" << tMeasure.getValue(state)
                 << " d/dt tMeasure=" << tMeasure.getValue(state,1)
                 << " d3/dt3 tMeasure=" << tMeasure.getValue(state,3)
                 << " 1000*tSubMeas=" << t1000.getValue(state)
                 << " t=" << t << endl;

            SimTK_TEST(tMeasure.getValue(state) == t);
            SimTK_TEST(tMeasure.getValue(state,1) == 1); // derivs
            SimTK_TEST(tMeasure.getValue(state,2) == 0);
            SimTK_TEST(tMeasure.getValue(state,3) == 0);

            SimTK_TEST(tSubMeas.getValue(state) == t);
            SimTK_TEST_EQ(t1000.getValue(state), 1000*t);

            cout << " tDelayed=" << tDelayed.getValue(state) << endl;
            if (t >= 0.01)
                SimTK_TEST_EQ(tDelayed.getValue(state), t-.01);

            cout << "q=" << state.getQ() << " u=" << state.getU() << endl;
            cout << "qSum=" << subsys.getQSum(state) << " uSum=" << subsys.getUSum(state) << endl;
            cout << "three=" << three.getValue(state) << " v3const=" << v3const.getValue(state) << endl;

            SimTK_TEST(zero.getValue(state)==0);
            SimTK_TEST(three.getValue(state)==3);
            SimTK_TEST(v3const.getValue(state)==Vec3(1,2,3));

            // This is calculated exactly.
            cout << "cos2pit=" << cos2pit.getValue(state) 
                 << " cos(2pi*t)=" << std::cos(2*Pi*t) << endl;
            SimTK_TEST_EQ(cos2pit.getValue(state), std::cos(2*Pi*t));

            cout << "Min(cos2pit)=" << minCos2pit.getValue(state) 
                 << " @t=" << minCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "Max(cos2pit)=" << maxCos2pit.getValue(state) 
                 << " @t=" << maxCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "MinAbs(cos2pit)=" << minAbsCos2pit.getValue(state) 
                 << " @t=" << minAbsCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "MaxAbs(cos2pit)=" << maxAbsCos2pit.getValue(state) 
                 << " @t=" << maxAbsCos2pit.getTimeOfExtremeValue(state) << endl;

            // This is calcuated by numerically integrating cos2pit.
            cout << "sin2pitOver2pi=" << sin2pitOver2pi.getValue(state) 
                 << " sin(2pi*t)/2pi=" << std::sin(2*Pi*t)/(2*Pi) << endl;
            SimTK_TEST_EQ_TOL(sin2pitOver2pi.getValue(state),
                              std::sin(2*Pi*t)/(2*Pi), 1e-3);

            // This should be an exact derivative, that is cos2pit.
            cout << "d/dt sin2pitOver2pi=" 
                 << sin2pitOver2pi.getValue(state,1) << endl;
            SimTK_TEST_EQ(sin2pitOver2pi.getValue(state,1), std::cos(2*Pi*t));

            // But this one is calculated by numerically differencing the
            // already-approximate integral (not good at t=0).
            cout << "dInteg=" 
                 << dInteg.getValue(state) << endl;
            if (i > 0) {
                SimTK_TEST_EQ_TOL(dInteg.getValue(state),
                                  std::cos(2*Pi*t), 1e-2);
            }

            // These should be -cos(t), sin(t), computed by integration.
            cout << "cossin=" << cossin.getValue(state) << "\n";
            cout << "vcossin=" << vcossin.getValue(state) << "\n";
            SimTK_TEST_EQ_TOL(cossin.getValue(state), 
                              Vec2(-std::cos(t),std::sin(t)), 1e-3);
            SimTK_TEST_EQ_TOL(vcossin.getValue(state), 
                              Vector(Vec2(-std::cos(t),std::sin(t))), 1e-3);

            // Delayed, so should be -cos(t-.1), sin(t-.1).
            cout << "vcossin delay .1=" << vcossin_delaypt1.getValue(state) << "\n";
            if (t >= 0.1)
                SimTK_TEST_EQ_TOL(vcossin_delaypt1.getValue(state), 
                                  Vector(Vec2(-std::cos(t-.1),std::sin(t-.1))), 
                                  1e-3);
        }

        if (i == nSteps)
            break;

        // Adjust step size to account for roundoff so that 2000 1ms steps
        // adds up to 2 seconds.
        const double hActual = (i+1)*h - state.getTime();
        if (ts.step(hActual))
            break;
    }

    auto& event = subsys.getMyEvent();

    printf("DONE. Time=%.15g, #events=%lld\n", ts.getCurrentState().getTime(),
           event.getNumOccurrences());

    // With the step size adjustments above we should end almost exactly at 2.
    SimTK_TEST_EQ(ts.getCurrentState().getTime(), 2);

    // The SomethingNeedsToHappen event should have occurred once for each of
    // the two designated-time witnesses, and once for the u-sum witness.
    SimTK_TEST(event.getNumOccurrences() == 3);

}

//==============================================================================
//                         TEST DESIGNATED TIMER
//==============================================================================
void testDesignatedTimer() {
    EventTrigger::Timer::Designated d("my timer", {.1,.3,4.2,.1999,1e-5});
    d.insertDesignatedTime(.707);
    d.insertDesignatedTimes({1,2,.3,3,.01,9,1,.9,2,1,2});
    d.insertDesignatedTime(.707);
    d.insertDesignatedTimes({1,2,3,4,5,6});

    // float times picked for exact representations so we can check them
    // to double precision below.
    std::set<float> goodTimes {99.875f,99.125f,25.75f};

    d.insertDesignatedTimes(goodTimes.begin(), goodTimes.end());

    // This is the right answer for the sorted/uniqued times.
    const double right[]={1e-005,0.01,0.1,0.1999,0.3,0.707,0.9,
                          1,2,3,4,4.2,5,6,9,25.75,99.125,99.875};

    for (unsigned i=0; i < d.getNumDesignatedTimes(); ++i)
        SimTK_TEST_EQ(d.getDesignatedTime(i), right[i]); // default tol

    auto e(d); // test copying of designated timer
    for (unsigned i=0; i < e.getNumDesignatedTimes(); ++i)
        SimTK_TEST_EQ(e.getDesignatedTime(i), right[i]); // default tol

    TestSystem sys; State s = sys.realizeTopology(); 

    s.setTime(.8);
    SimTK_TEST_EQ(d.calcTimeOfNextTrigger(sys,s,0), 0.9);

    s.setTime(2);
    SimTK_TEST(d.calcTimeOfNextTrigger(sys,s,0)==2);
    SimTK_TEST(d.calcTimeOfNextTrigger(sys,s,2)==3);

    s.setTime(0);
    SimTK_TEST_EQ(d.calcTimeOfNextTrigger(sys,s,-Infinity), 1e-005);

    s.setTime(100);
    SimTK_TEST(d.calcTimeOfNextTrigger(sys,s,-Infinity) == Infinity);
}


//==============================================================================
//                                   MAIN
//==============================================================================
int main() {              
    SimTK_START_TEST("TestSimulation");
        SimTK_SUBTEST(testOne);
        SimTK_SUBTEST(testDesignatedTimer);
    SimTK_END_TEST();
}

void FixedTimeStepper::initialize(const State& initState, Real consTol) {
    m_state = initState;
    m_consTol = consTol;

    // Handle the initialization event.
    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;
    getSystem().noteEventOccurrence({&getSystem().getInitializationTrigger()},
                                    triggeredEvents, ignoredEventIds);

    EventChangeResult result;
    getSystem().performEventChangeActions(*this, triggeredEvents, result);
    SimTK_ERRCHK_ALWAYS(
        result.getExitStatus()==EventChangeResult::Succeeded,
        "FixedTimeStepper::initialize()", 
        "An initialization event handler failed or requested termination.");

    // Realize the starting state to Acceleration stage; then take care of
    // auto-update discrete variables (which leave stage unchanged).
    m_system.realize(m_state, Stage::Acceleration);
    m_state.autoUpdateDiscreteVariables();

    // Issue initialization event report actions if any.
    getSystem().performEventReportActions(*this, triggeredEvents);

    // Remember start-of-step values for witness functions.
    getSystem().findActiveEventWitnesses(*this, m_witnesses);
    calcWitnessValues(m_state, m_witnessValuesPrev);
}

// On entry, m_state will typically already be realized to Acceleration stage. 
bool FixedTimeStepper::step(double h) {
    const Real h2 = h/2;
    auto& sys   = m_system;
    auto& state = m_state;

    sys.realize(state, Stage::Acceleration); // make sure

    // Commit the values for the discrete variable updates calculated
    // at the end of the previous step (no change to stage).
    state.autoUpdateDiscreteVariables();

    // This is an alias for ydot in the state, so will get updated as we go.
    const Vector& ydot = state.getYDot();

    // Record the witness values at start of step.
    calcWitnessValues(state, m_witnessValuesPrev);

    // First integrator stage: unconstrained continuous system only.
    state.updY()    += h2*ydot;
    state.updTime() += h2;
    sys.realize(state, Stage::Time);
    sys.prescribeQ(state);
    sys.realize(state, Stage::Position);
    sys.prescribeU(state);
    sys.realize(state, Stage::Acceleration); // updates ydot

    // Second (final) integrator stage.
    // 1. Unconstrained continuous system.
    state.updY() += h2*ydot; // that is, y = y0 + h*(ydot0+ydot)/2
    state.updTime() += h2;
    sys.realize(state, Stage::Time);
    sys.prescribeQ(state);

    EventTriggers   triggered;

    // 2. Deal with time-dependent events.
    findEvents(Stage::Time, triggered);
    if (handleEvents(Stage::Time, triggered))
        return true;

    sys.realize(state, Stage::Position);
    sys.prescribeU(state);

    // 3a. Project position-dependent constraints.
    Vector temp;
    ProjectOptions opts(m_consTol);
    ProjectResults results;
    sys.projectQ(state, temp, opts, results);

    // 3b. Handle position-dependent events.
    findEvents(Stage::Position, triggered);
    if (handleEvents(Stage::Position, triggered))
        return true;

    // 4a. Project velocity-dependent constraints.
    sys.realize(state, Stage::Velocity);
    sys.projectU(state, temp, opts, results);

    // 4b. Handle velocity-dependent events.
    findEvents(Stage::Velocity, triggered);
    if (handleEvents(Stage::Velocity, triggered))
        return true;

    // 5. Handle dynamics- and acceleration-dependent events.
    sys.realize(state, Stage::Dynamics);
    findEvents(Stage::Dynamics, triggered);
    if (handleEvents(Stage::Dynamics, triggered))
        return true;

    sys.realize(state, Stage::Acceleration);
    findEvents(Stage::Acceleration, triggered);
    if (handleEvents(Stage::Acceleration, triggered))
        return true;

    // Ensure State is realized through Acceleration Stage.
    sys.realize(state, Stage::Acceleration);
    return false;
}


// Find the event triggers at a particular stage that changed sign since
// they were last recorded in m_witnessValuesPrev.
void FixedTimeStepper::
findEvents(Stage stage, EventTriggers& triggered) const {
    triggered.clear();

    for (ActiveWitnessIndex awx(0); awx < m_witnesses.size(); ++awx) {
        const EventWitness& w = *m_witnesses[awx];
        if (w.getDependsOnStage() != stage)
            continue;
        const EventWitness::Value oldVal = m_witnessValuesPrev[awx];
        const EventWitness::Value newVal = 
            w.calcWitnessValue(*this, m_state);

        // This sign() function return -1,0,1. Transitions *to* zero are events;
        // transitions *from* zero are not. That avoids double-counting.
        if (oldVal.getSign() && (newVal.getSign() != oldVal.getSign())) {
            triggered.push_back(&w);
        }
    }
}

bool FixedTimeStepper::
handleEvents(Stage stage, const EventTriggers& triggered) {
    auto& sys = m_system;
    auto& state = m_state;

    if (triggered.empty())
        return false;

    EventsAndCauses triggeredEvents;
    Array_<EventId> ignoredEventIds;
    getSystem().noteEventOccurrence(triggered,
                                    triggeredEvents, ignoredEventIds);

    cout << "==> Handling " << triggered.size() 
         << " events at Stage " << stage << ":";
    for (unsigned i=0; i < triggered.size(); ++i)
        cout << " " << triggered[i];
    cout << endl;

    EventChangeResult result;

    Array_<StageVersion> stageVersions;
    state.getSystemStageVersions(stageVersions);
    cout << "BEFORE handling stage versions=\n";
    cout << stageVersions << "\n";

    sys.performEventChangeActions(*this, triggeredEvents, result);

    state.getSystemStageVersions(stageVersions);
    cout << "AFTER handling stage versions=\n";
    cout << stageVersions << "\n";
    cout << "Result lowestStage=" << result.getLowestModifiedStage() <<"\n";

    if (result.getExitStatus()==EventChangeResult::ShouldTerminate) {
        cout << "==> Event at Stage " << stage 
             << " requested termination at t=" << state.getTime() << endl;
        return true;
    }

    return false;
}

