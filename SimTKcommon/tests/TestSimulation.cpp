/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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

/*
 * Here we'll build a very simple System containing a simple
 * Subsystem, and integrate with a simple integrator. This avoids all the usual
 * Simmath and Simbody trappings and lets us just check out the underlying
 * simulation architecture.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

#include <map>
#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS((cond), "Assertion failed.");}
#define ASSERT_EQ(v1,v2)    \
    {SimTK_ASSERT_ALWAYS(numericallyEqual((v1),(v2)),   \
     "Values should have been numerically equivalent.");}

// Scale by the magnitude of the quantities being compared, so that we don't
// ask for unreasonable precision. For magnitudes near zero, we'll be satisfied
// if both are very small without demanding that they must also be relatively
// close. That is, we use a relative tolerance for big numbers and an absolute
// tolerance for small ones.
bool numericallyEqual(float v1, float v2) {
    const float scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1f);
    return std::abs(v1-v2) < scale*NTraits<float>::getSignificant();
}
bool numericallyEqual(double v1, double v2) {
    const double scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1);
    return std::abs(v1-v2) < scale*NTraits<double>::getSignificant();
}
template <class P>
bool numericallyEqual(const std::complex<P>& v1, const std::complex<P>& v2) {
    return numericallyEqual(v1.real(), v2.real())
            && numericallyEqual(v1.imag(), v2.imag());
}
template <class P>
bool numericallyEqual(const conjugate<P>& v1, const conjugate<P>& v2) {
    return numericallyEqual(v1.real(), v2.real())
            && numericallyEqual(v1.imag(), v2.imag());
}
template <class P>
bool numericallyEqual(const std::complex<P>& v1, const conjugate<P>& v2) {
    return numericallyEqual(v1.real(), v2.real())
            && numericallyEqual(v1.imag(), v2.imag());
}
template <class P>
bool numericallyEqual(const conjugate<P>& v1, const std::complex<P>& v2) {
    return numericallyEqual(v1.real(), v2.real())
            && numericallyEqual(v1.imag(), v2.imag());
}
template <class P>
bool numericallyEqual(const negator<P>& v1, const negator<P>& v2) {
    return numericallyEqual(-v1, -v2);  // P, P
}
template <class P>
bool numericallyEqual(const P& v1, const negator<P>& v2) {
    return numericallyEqual(-v1, -v2);  // P, P
}
template <class P>
bool numericallyEqual(const negator<P>& v1, const P& v2) {
    return numericallyEqual(-v1, -v2);  // P, P
}
template <class P>
bool numericallyEqual(const negator<std::complex<P> >& v1, const conjugate<P>& v2) {
    return numericallyEqual(-v1, -v2);  // complex, conjugate
}
template <class P>
bool numericallyEqual(const negator<conjugate<P> >& v1, const std::complex<P>& v2) {
    return numericallyEqual(-v1, -v2);  // conjugate, complex
}
template <class P>
bool numericallyEqual(const std::complex<P>& v1, const negator<conjugate<P> >& v2) {
    return numericallyEqual(-v1, -v2); // complex, conjugate
}
template <class P>
bool numericallyEqual(const conjugate<P>& v1, const negator<std::complex<P> >& v2) {
    return numericallyEqual(-v1, -v2); // conjugate, complex
}

namespace SimTK {
typedef std::map<EventId, SubsystemIndex> EventRegistry;
std::ostream& operator<<(std::ostream& o, const EventRegistry&) {return o;}
}

class SystemSubsystemGuts : public Subsystem::Guts {
public:
    SystemSubsystemGuts() {}

    const EventRegistry& getEventRegistry(const State& s) const
    {   assert(eventRegistry >= 0);
        return Value<EventRegistry>::downcast(getCacheEntry(s, eventRegistry)); }
    EventRegistry& updEventRegistry(const State& s) const
    {   assert(eventRegistry >= 0);
        return Value<EventRegistry>::updDowncast(updCacheEntry(s, eventRegistry)); }

    // implementations of Subsystem::Guts virtuals
    SystemSubsystemGuts* cloneImpl() const override
    {   return new SystemSubsystemGuts(*this); }
    int realizeSubsystemTopologyImpl(State& s) const override {
        eventRegistry = allocateCacheEntry(s, Stage::Instance, 
                                           new Value<EventRegistry>());
        return 0;
    }

private:
    mutable CacheEntryIndex eventRegistry;
};


class SystemSubsystem : public Subsystem {
public:
    SystemSubsystem() {
        adoptSubsystemGuts(new SystemSubsystemGuts());
    }

    void registerEventsToSubsystem(const State& s, const Subsystem::Guts& sub, 
                                   EventId start, int nEvents) const
    {   EventRegistry& er = getGuts().updEventRegistry(s);
        SubsystemIndex sx = sub.getMySubsystemIndex();
        for (int i=start; i < start+nEvents; ++i)
            er[EventId(i)] = sx;
    }

    const EventRegistry& getEventRegistry(const State& s) const 
    {   return getGuts().getEventRegistry(s); }

private:
    const SystemSubsystemGuts& getGuts() const
    {   return dynamic_cast<const SystemSubsystemGuts&>(getSubsystemGuts());}
    SystemSubsystemGuts& updGuts()
    {   return dynamic_cast<SystemSubsystemGuts&>(updSubsystemGuts());}
};


class TestSystemGuts : public System::Guts {
public:
    const SystemSubsystem& getSystemSubsystem() const {return syssub;}
    SystemSubsystem& updSystemSubsystem() {return syssub;}

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
    void handleEventsImpl
       (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
        const HandleEventsOptions& options, HandleEventsResults& results) const
        override
    {
        cout << "handleEventsImpl t=" << s.getTime() 
             << " cause=" << Event::getCauseName(cause) << endl;

        if (eventIds.empty()) {
            for (SubsystemIndex sx(0); sx < getNumSubsystems(); ++sx) {
                const Subsystem& sub = getSubsystem(sx);
                sub.getSubsystemGuts().handleEvents(s, cause, eventIds, 
                    options, results);
                if (results.getExitStatus()==HandleEventsResults::Failed)
                    break;
            }
            return;
        }

        // If there are EventIds, dole them out to the owning subsystem.

        const EventRegistry& registry = getSystemSubsystem().getEventRegistry(s);

        std::map<SubsystemIndex, Array_<EventId> > eventsPerSub;
        for (EventId eid(0); eid < eventIds.size(); ++eid)
            eventsPerSub[ registry.find(eid)->second ].push_back(eid);

        std::map<SubsystemIndex, Array_<EventId> >::const_iterator 
            i = eventsPerSub.begin();
        for (; i != eventsPerSub.end(); ++i) {
            const Subsystem& sub = getSubsystem(i->first);
            sub.getSubsystemGuts().handleEvents(s, cause, i->second, 
                options, results);
            if (results.getExitStatus()==HandleEventsResults::Failed)
                break;
        }
    }

    int reportEventsImpl(const State& s, Event::Cause cause, 
                         const Array_<EventId>& eventIds) const override
    {
        cout << "reportEventsImpl t=" << s.getTime() << " cause=" 
             << Event::getCauseName(cause) << endl;
        return 0;
    }

private:
    SystemSubsystem syssub;
};

class TestSystem : public System {
public:
    TestSystem() {
        adoptSystemGuts(new TestSystemGuts());
        adoptSubsystem(updGuts().updSystemSubsystem());
    }

    const SystemSubsystem& getSystemSubsystem() const 
    {   return getGuts().getSystemSubsystem(); }
    SystemSubsystem& updSystemSubsystem()
    {   return updGuts().updSystemSubsystem(); }

    void registerEventsToSubsystem(const State& s, const Subsystem::Guts& sub, 
                                   EventId start, int nEvents) const 
    {
        const SystemSubsystem& syssub = getGuts().getSystemSubsystem();
        syssub.registerEventsToSubsystem(s,sub,start,nEvents);
    }

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


class TestSubsystemGuts : public Subsystem::Guts {
    struct StateVars {
        QIndex myQs;
        UIndex myUs;
    };
    struct CacheEntries {
        CacheEntryIndex qSumCacheIx, uSumCacheIx;
        EventTriggerByStageIndex timeTriggerIx, velTriggerIx;
    };
public:
    TestSubsystemGuts() {}

    const Vec3& getQ3(const State& s) const {return Vec3::getAs(&getQ(s)[getStateVars(s).myQs]);}
    const Vec3& getU3(const State& s) const {return Vec3::getAs(&getU(s)[getStateVars(s).myQs]);}
    const Vec3& getQDot3(const State& s) const {return Vec3::getAs(&getQDot(s)[getStateVars(s).myQs]);}
    const Vec3& getUDot3(const State& s) const {return Vec3::getAs(&getUDot(s)[getStateVars(s).myUs]);}
    const Vec3& getQDotDot3(const State& s) const {return Vec3::getAs(&getQDotDot(s)[getStateVars(s).myQs]);}
    Real getQSum(const State& s) const {return Value<Real>::downcast(getCacheEntry(s,getCacheEntries(s).qSumCacheIx));}
    Real getUSum(const State& s) const {return Value<Real>::downcast(getCacheEntry(s,getCacheEntries(s).uSumCacheIx));}

    Real getTimeTrigger1(const State& s) const {return getEventTriggersByStage(s, Stage::Time)[getCacheEntries(s).timeTriggerIx];}
    Real getTimeTrigger2(const State& s) const {return getEventTriggersByStage(s, Stage::Time)[getCacheEntries(s).timeTriggerIx+1];}
    Real getVelTrigger(const State& s) const {return getEventTriggersByStage(s, Stage::Velocity)[getCacheEntries(s).velTriggerIx];}

    Vec3& updQ3(State& s) const {return Vec3::updAs(&updQ(s)[getStateVars(s).myQs]);}
    Vec3& updU3(State& s) const {return Vec3::updAs(&updU(s)[getStateVars(s).myUs]);}

    // implementations of Subsystem::Guts virtuals

    TestSubsystemGuts* cloneImpl() const override
    {   return new TestSubsystemGuts(*this); }


    int realizeSubsystemTopologyImpl(State& s) const override {
        myStateVars = allocateCacheEntry(s, Stage::Model, new Value<StateVars>());
        myCacheEntries = allocateCacheEntry(s, Stage::Instance, new Value<CacheEntries>());
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
        cache.qSumCacheIx = allocateCacheEntry(s, Stage::Position, new Value<Real>(0));
        cache.uSumCacheIx = allocateCacheEntry(s, Stage::Velocity, new Value<Real>(0));
        cache.timeTriggerIx = allocateEventTriggersByStage(s, Stage::Time, 2);
        cache.velTriggerIx  = allocateEventTriggersByStage(s, Stage::Velocity, 1);

        getTestSystem().registerEventsToSubsystem(s, *this, EventId(cache.timeTriggerIx), 2);
        getTestSystem().registerEventsToSubsystem(s, *this, EventId(cache.velTriggerIx), 1);
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const override {
        const Real TriggerTime1 = .6789, TriggerTime2 = 1.234;
        updTimeTrigger1(s) = s.getTime() - TriggerTime1;
        updTimeTrigger2(s) = s.getTime() - TriggerTime2;
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const override {
        updQSum(s) = sum(getQ3(s));
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const override {
        const Real TriggerUSum = 5;
        updQDot3(s) = getU3(s);
        const Real usum = updUSum(s) = sum(getU3(s));
        updVelTrigger(s) = usum - TriggerUSum; 
        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const override {
        updQDotDot3(s) = updUDot3(s) = Vec3(1,2,3);
        return 0;
    }

    void handleEventsImpl(State& s, Event::Cause cause, 
                          const Array_<EventId>& eventIds,
                          const HandleEventsOptions& options,
                          HandleEventsResults& results) const override
    {
        cout << "**** TestSubsystem::handleEventsImpl t=" << s.getTime() 
             << " acc=" << options.getAccuracy()
             << " eventIds=";
        for (unsigned i=0; i < eventIds.size(); ++i)
           cout << " " << eventIds[i];
        cout << " ****" << endl;

        // Pretend we changed a position to test lowestModifiedStage
        // calculation. Try to hide our duplicity by realizing it again.
        s.invalidateAllCacheAtOrAbove(Stage::Position); 
        getSystem().realize(s, Stage::Velocity);
        results.setExitStatus(HandleEventsResults::Succeeded);
    }

    void reportEventsImpl(const State&, Event::Cause, 
                          const Array_<EventId>& eventIds) const override
    {
    }


private:
    Vec3& updQDot3(const State& s) const {return Vec3::updAs(&updQDot(s)[getStateVars(s).myQs]);}
    Vec3& updUDot3(const State& s) const {return Vec3::updAs(&updUDot(s)[getStateVars(s).myUs]);}
    Vec3& updQDotDot3(const State& s) const {return Vec3::updAs(&updQDotDot(s)[getStateVars(s).myQs]);}
    Real& updQSum(const State& s) const {return Value<Real>::updDowncast(updCacheEntry(s,getCacheEntries(s).qSumCacheIx));}
    Real& updUSum(const State& s) const {return Value<Real>::updDowncast(updCacheEntry(s,getCacheEntries(s).uSumCacheIx));}
    Real& updTimeTrigger1(const State& s) const {return updEventTriggersByStage(s, Stage::Time)[getCacheEntries(s).timeTriggerIx];}
    Real& updTimeTrigger2(const State& s) const {return updEventTriggersByStage(s, Stage::Time)[getCacheEntries(s).timeTriggerIx+1];}
    Real& updVelTrigger(const State& s) const {return updEventTriggersByStage(s, Stage::Velocity)[getCacheEntries(s).velTriggerIx];}

    const StateVars& getStateVars(const State& s) const
    {   assert(myStateVars >= 0);
        return Value<StateVars>::downcast(getCacheEntry(s,myStateVars)); }
    const CacheEntries& getCacheEntries(const State& s) const
    {   assert(myCacheEntries >= 0);
        return Value<CacheEntries>::downcast(getCacheEntry(s,myCacheEntries)); }
    StateVars& updStateVars   (const State& s) const 
    {   assert(myStateVars >= 0);
        return Value<StateVars>::updDowncast(updCacheEntry(s,myStateVars)); }
    CacheEntries& updCacheEntries(const State& s) const
    {   assert(myCacheEntries >= 0);
        return Value<CacheEntries>::updDowncast(updCacheEntry(s,myCacheEntries)); }

    const TestSystem& getTestSystem() const {return TestSystem::getAs(getSystem());}
    TestSystem& updTestSystem() {return TestSystem::updAs(updSystem());}

        // TOPOLOGY STATE VARIABLES //

    Array_<EventHandler*>          eventHandlers;
    mutable Array_<EventReporter*> eventReporters;

        // TOPOLOGY CACHE //
    mutable CacheEntryIndex myStateVars;
    mutable CacheEntryIndex myCacheEntries;

};

// This Subsystem has 3 q's and 3 u's of its own, as well as whatever State
// variables its Measures require.
class TestSubsystem : public Subsystem {
public:
    TestSubsystem(System& sys) {
        adoptSubsystemGuts(new TestSubsystemGuts());
        sys.adoptSubsystem(*this);
    }

    Real getQSum(const State& s) {return getGuts().getQSum(s);}
    Real getUSum(const State& s) {return getGuts().getUSum(s);}
private:
    const TestSubsystemGuts& getGuts() const
    {   return dynamic_cast<const TestSubsystemGuts&>(getSubsystemGuts());}
    TestSubsystemGuts& updGuts()
    {   return dynamic_cast<TestSubsystemGuts&>(updSubsystemGuts());}
};

// Find the event triggers at a particular stage that changed sign since
// they were last recorded in events0.
static void findEvents(const State& state, Stage g, const Vector& triggers0,
                       Array_<EventId>& triggered)
{
    const int n     = state.getNEventTriggersByStage(g);
    const int start = state.getEventTriggerStartByStage(g); // location within triggers0 Vector

    const Vector& stageTriggers = state.getEventTriggersByStage(g);

    triggered.clear();
    for (int i=0; i < n; ++i) {
        const EventId allStageId = EventId(start + i);
        if (sign(triggers0[allStageId]) != sign(stageTriggers[i]))
            triggered.push_back(allStageId);
    }
}

static Real   accuracy = 1e-6;
static Real   timescale;

static bool handleEvents(const System& sys, State& state, Stage g,
                         const Array_<EventId>& triggered) 
{
    if (triggered.empty())
        return false;

    cout << "==> Handling " << triggered.size() << " events at Stage " << g << ":";
    for (unsigned i=0; i < triggered.size(); ++i)
        cout << " " << triggered[i];
    cout << endl;

    bool shouldTerminate = false; 
    HandleEventsOptions options(accuracy);
    HandleEventsResults results;

    Array_<StageVersion> stageVersions;
    state.getSystemStageVersions(stageVersions);
    cout << "BEFORE handling stage versions=\n";
    cout << stageVersions << "\n";
    sys.handleEvents(state, Event::Cause::Triggered, triggered,
                     options, results);
    state.getSystemStageVersions(stageVersions);
    cout << "AFTER handling stage versions=\n";
    cout << stageVersions << "\n";
    cout << "Results lowestStage=" << results.getLowestModifiedStage() <<"\n";

    if (results.getExitStatus()==HandleEventsResults::ShouldTerminate) {
        cout << "==> Event at Stage " << g 
             << " requested termination at t=" << state.getTime() << endl;
        shouldTerminate = true;
    }

    return shouldTerminate;
}

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
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    int getNumTimeDerivativesVirtual() const {return 0;}
    Stage getDependsOnStageVirtual(int order) const 
    {   return Stage::Time; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
    {
        SimTK_ASSERT1_ALWAYS(derivOrder==0,
            "MySinCos::Implementation::calcCachedValueVirtual():"
            " derivOrder %d seen but only 0 allowed.", derivOrder);

        value[0] = std::sin(s.getTime());
        value[1] = std::cos(s.getTime());
    }
};

template <class T>
class MyRealMeasure : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(MyRealMeasure, Measure);

    SimTK_MEASURE_HANDLE_POSTSCRIPT(MyRealMeasure, Measure);
};

template <class T>
class MyRealMeasure<T>::Implementation : public Measure_<T>::Implementation {
public:
    Implementation* cloneVirtual() const override
    {   return new Implementation(*this); }
    int getNumTimeDerivativesVirtual() const override 
    {   return 0; }
    Stage getDependsOnStageVirtual(int order) const override 
    {   return Stage::Time; }

};


void testOne() {
    TestSystem sys;
    TestSubsystem subsys(sys);

    // Add a Result measure to the system subsystem. This depends on 
    // Position stage and invalidates Dynamics and later stages.
    Measure_<Vector>::Result vectorResult(sys.updSystemSubsystem(), 
        Stage::Position, Stage::Dynamics);

    MeasureIndex vectorResultIx = vectorResult.getSubsystemMeasureIndex();
    cout << "vectorResult index=" 
         << vectorResultIx << endl;

    Measure_<Vector>::Result myVecRes =  Measure_<Vector>::Result::getAs(
        sys.updSystemSubsystem().getMeasure(vectorResultIx));

    Measure::Result result(sys.updSystemSubsystem(), 
        Stage::Time, Stage::Position);
    Measure::Result autoResult(sys.updSystemSubsystem(), 
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

    Measure m;
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

    // Explicit midpoint steps.
    const Real h = .001;
    const int nSteps = 2000;
    const int outputInterval = 100;
    state.setTime(0);

    ASSERT(state.getTime()==0);


    //initialize()
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

    // Fill in statics above.
    timescale = sys.getDefaultTimeScale();

    sys.realize(state, Stage::Acceleration);

    cout << "Now stage=" << state.getSystemStage() << endl;
    vectorResult.setValue(state, Vector(5,9));
    cout << "After vectorResult.setValue(), vectorResult="
         << vectorResult.getValue(state) << endl;
    cout << "... but stage=" << state.getSystemStage() << endl;

    cossin.setValue(state, cossinInit.getValue(state));
    vcossin.setValue(state, vcossinInit.getValue(state));

    // Handler is allowed to throw an exception if it fails since we don't
    // have a way to recover.
    HandleEventsOptions handleOpts;
    HandleEventsResults results;
    sys.handleEvents(state, Event::Cause::Initialization,
                     Array_<EventId>(), handleOpts, results);
    SimTK_ERRCHK_ALWAYS(
        results.getExitStatus()!=HandleEventsResults::ShouldTerminate,
        "Integrator::initialize()", 
        "An initialization event handler requested termination.");

    sys.realize(state, Stage::Acceleration);
    state.autoUpdateDiscreteVariables(); // ??

    for (int i=0; i <= nSteps; ++i) {

        if (i % outputInterval == 0) {
            sys.realize(state, Stage::Report);
            cout << "\ntMeasure=" << tMeasure.getValue(state)
                 << " d/dt tMeasure=" << tMeasure.getValue(state,1)
                 << " d3/dt3 tMeasure=" << tMeasure.getValue(state,3)
                 << " 1000*tSubMeas=" << t1000.getValue(state)
                 << " t=" << state.getTime() << endl;
            cout << " tDelayed=" << tDelayed.getValue(state) << endl;
            cout << "q=" << state.getQ() << " u=" << state.getU() << endl;
            cout << "qSum=" << subsys.getQSum(state) << " uSum=" << subsys.getUSum(state) << endl;
            cout << "three=" << three.getValue(state) << " v3const=" << v3const.getValue(state) << endl;
            cout << "cos2pit=" << cos2pit.getValue(state) 
                 << " cos(2pi*t)=" << std::cos(2*Pi*state.getTime()) << endl;
            cout << "Min(cos2pit)=" << minCos2pit.getValue(state) 
                 << " @t=" << minCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "Max(cos2pit)=" << maxCos2pit.getValue(state) 
                 << " @t=" << maxCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "MinAbs(cos2pit)=" << minAbsCos2pit.getValue(state) 
                 << " @t=" << minAbsCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "MaxAbs(cos2pit)=" << maxAbsCos2pit.getValue(state) 
                 << " @t=" << maxAbsCos2pit.getTimeOfExtremeValue(state) << endl;
            cout << "sin2pitOver2pi=" << sin2pitOver2pi.getValue(state) 
                 << " sin(2pi*t)/2pi=" << std::sin(2*Pi*state.getTime())/(2*Pi) << endl;
            cout << "d/dt sin2pitOver2pi=" 
                 << sin2pitOver2pi.getValue(state,1) << endl;
            cout << "dInteg=" 
                 << dInteg.getValue(state) << endl;
            cout << "cossin=" << cossin.getValue(state) << "\n";
            cout << "vcossin=" << vcossin.getValue(state) << "\n";
            cout << "vcossin delay .1=" << vcossin_delaypt1.getValue(state) << "\n";
        }

        if (i == nSteps)
            break;

        const Real h2 = h/2;
        const Vector ydot0 = state.getYDot();
        const Vector triggers0 = state.getEventTriggers();
        Array_<EventId> triggered;

        // Commit the values for the discrete variable updates calculated
        // at the end of the previous step. This includes both explicitly
        // discrete variables and continuous variables which are defined
        // by algebraic rather than differential equations, such as
        // prescribed motions. This requires that all calculations have
        // been performed already using the *updated* values, *not* the
        // state values; that permits us to perform this update without
        // invalidating any cache entries.
        state.autoUpdateDiscreteVariables();

        // First integrator stage: unconstrained continuous system only.
        state.updY()    += h2*ydot0;
        state.updTime() += h2;
        sys.realize(state, Stage::Time);
        sys.prescribeQ(state);
        sys.realize(state, Stage::Position);
        sys.prescribeU(state);
        sys.realize(state);

        // Second (final) integrator stage.
        // 1. Unconstrained continuous system.
        const Vector& ydot = state.getYDot();
        state.updY() += h2*ydot; // that is, y = y0 + h*(ydot0+ydot)/2
        state.updTime() += h2;
        sys.realize(state, Stage::Time);
        sys.prescribeQ(state);

        // 2. Deal with time-dependent events.
        findEvents(state, Stage::Time, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Time, triggered))
            break;

        sys.realize(state, Stage::Position);
        sys.prescribeU(state);

        // 3a. Project position-dependent constraints.
        Vector temp;
        ProjectOptions opts(accuracy);
        ProjectResults results;
        sys.projectQ(state, temp, opts, results);

        // 3b. Handle position-dependent events.
        findEvents(state, Stage::Position, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Position, triggered))
            break;

        // 4a. Project velocity-dependent constraints.
        sys.realize(state, Stage::Velocity);
        sys.projectU(state, temp, opts, results);

        // 4b. Handle velocity-dependent events.
        findEvents(state, Stage::Velocity, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Velocity, triggered))
            break;

        // 5. Handle dynamics- and acceleration-dependent events.
        sys.realize(state, Stage::Dynamics);
        findEvents(state, Stage::Dynamics, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Dynamics, triggered))
            break;
        sys.realize(state, Stage::Acceleration);
        findEvents(state, Stage::Acceleration, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Acceleration, triggered))
            break;

        // Ensure State is realized through Acceleration Stage.
        sys.realize(state, Stage::Acceleration);
    }
}

int main() {
    try {
        testOne();
    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

