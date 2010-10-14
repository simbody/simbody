/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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


/*
 * Here we'll build a very simple System containing a simple
 * Subsystem, and integrate with a simple integrator. This avoids all the usual
 * Simmath and Simbody trappings and lets us just check out the underlying
 * simulation architecture.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

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
        return Value<EventRegistry>::downcast(updCacheEntry(s, eventRegistry)); }

    // implementations of Subsystem::Guts virtuals
    SystemSubsystemGuts* cloneImpl() const {return new SystemSubsystemGuts(*this);}
    int realizeSubsystemTopologyImpl(State& s) const {
        eventRegistry = allocateCacheEntry(s, Stage::Instance, new Value<EventRegistry>());
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

    void registerEventsToSubsystem(const State& s, const Subsystem::Guts& sub, EventId start, int nEvents) const
    {   EventRegistry& er = getGuts().updEventRegistry(s);
        SubsystemIndex sx = sub.getMySubsystemIndex();
        for (int i=start; i < start+nEvents; ++i)
            er[EventId(i)] = sx;
    }

    const EventRegistry& getEventRegistry(const State& s) const {return getGuts().getEventRegistry(s);}

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
    TestSystemGuts* cloneImpl() const {return new TestSystemGuts(*this);}

    int projectImpl(State&, Real consAccuracy, const Vector& yweights,
                    const Vector& ootols, Vector& yerrest, System::ProjectOptions) const
    {
        return 0;
    }

    int handleEventsImpl
       (State& s, Event::Cause cause, const Array_<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const
    {
        cout << "handleEventsImpl t=" << s.getTime() << " cause=" << Event::getCauseName(cause) << endl;
        const EventRegistry& registry = getSystemSubsystem().getEventRegistry(s);

        std::map<SubsystemIndex, Array_<EventId> > eventsPerSub;
        for (EventId eid(0); eid < eventIds.size(); ++eid)
            eventsPerSub[ registry.find(eid)->second ].push_back(eid);

        std::map<SubsystemIndex, Array_<EventId> >::const_iterator i = eventsPerSub.begin();
        for (; i != eventsPerSub.end(); ++i) {
            Stage lowest = Stage::Report;
            bool  terminate = false;
            const Subsystem& sub = getSubsystem(i->first);
            sub.getSubsystemGuts().handleEvents(s, cause, i->second, accuracy, yWeights, ooConstraintTols, lowest, terminate);
            if (lowest < lowestModified)
                lowestModified = lowest;
            if (terminate) {
                shouldTerminate = true;
                break;
            }
        }
            
        return 0;
    }

    int reportEventsImpl(const State& s, Event::Cause cause, const Array_<EventId>& eventIds) const
    {
        cout << "reportEventsImpl t=" << s.getTime() << " cause=" << Event::getCauseName(cause) << endl;
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

    void registerEventsToSubsystem(const State& s, const Subsystem::Guts& sub, EventId start, int nEvents) const {
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
    friend std::ostream& operator<<(std::ostream& o, const CacheEntries&);
    friend std::ostream& operator<<(std::ostream& o, const StateVars&);
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

    TestSubsystemGuts* cloneImpl() const {return new TestSubsystemGuts(*this);}


    int realizeSubsystemTopologyImpl(State& s) const {
        myStateVars = allocateCacheEntry(s, Stage::Model, new Value<StateVars>());
        myCacheEntries = allocateCacheEntry(s, Stage::Instance, new Value<CacheEntries>());
        return 0;
    }


    int realizeSubsystemModelImpl(State& s) const {
        StateVars& vars = updStateVars(s);
        vars.myQs = allocateQ(s, Vector(Vec3(0)));
        vars.myUs = allocateU(s, Vector(Vec3(0)));
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        CacheEntries& cache = updCacheEntries(s);
        cache.qSumCacheIx = allocateCacheEntry(s, Stage::Position, new Value<Real>(0));
        cache.uSumCacheIx = allocateCacheEntry(s, Stage::Velocity, new Value<Real>(0));
        cache.timeTriggerIx = allocateEventTriggersByStage(s, Stage::Time, 2);
        cache.velTriggerIx  = allocateEventTriggersByStage(s, Stage::Velocity, 1);

        // TODO: this doesn't work
        getTestSystem().registerEventsToSubsystem(s, *this, EventId(cache.timeTriggerIx), 2);
        getTestSystem().registerEventsToSubsystem(s, *this, EventId(cache.velTriggerIx), 1);
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        const Real TriggerTime1 = .6789, TriggerTime2 = 1.234;
        updTimeTrigger1(s) = s.getTime() - TriggerTime1;
        updTimeTrigger2(s) = s.getTime() - TriggerTime2;
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        updQSum(s) = sum(getQ3(s));
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        const Real TriggerUSum = 5;
        updQDot3(s) = getU3(s);
        const Real usum = updUSum(s) = sum(getU3(s));
        updVelTrigger(s) = usum - TriggerUSum; 
        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const {
        updQDotDot3(s) = updUDot3(s) = Vec3(1,2,3);
        return 0;
    }

    void handleEvents(State& s, Event::Cause cause, const Array_<EventId>& eventIds,
        Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols,
        Stage& lowestModified, bool& shouldTerminate) const
    {
        cout << "**** TestSubsystem::handleEvents t=" << s.getTime() << " eventIds=";
        for (unsigned i=0; i < eventIds.size(); ++i)
           cout << " " << eventIds[i];
        cout << " ****" << endl;
    }

    void reportEvents(const State&, Event::Cause, const Array_<EventId>& eventIds) const
    {
    }


private:
    Vec3& updQDot3(const State& s) const {return Vec3::updAs(&updQDot(s)[getStateVars(s).myQs]);}
    Vec3& updUDot3(const State& s) const {return Vec3::updAs(&updUDot(s)[getStateVars(s).myUs]);}
    Vec3& updQDotDot3(const State& s) const {return Vec3::updAs(&updQDotDot(s)[getStateVars(s).myQs]);}
    Real& updQSum(const State& s) const {return Value<Real>::downcast(updCacheEntry(s,getCacheEntries(s).qSumCacheIx));}
    Real& updUSum(const State& s) const {return Value<Real>::downcast(updCacheEntry(s,getCacheEntries(s).uSumCacheIx));}
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
        return Value<StateVars>::downcast(updCacheEntry(s,myStateVars)); }
    CacheEntries& updCacheEntries(const State& s) const
    {   assert(myCacheEntries >= 0);
        return Value<CacheEntries>::downcast(updCacheEntry(s,myCacheEntries)); }

    const TestSystem& getTestSystem() const {return TestSystem::getAs(getSystem());}
    TestSystem& updTestSystem() {return TestSystem::updAs(updSystem());}

        // TOPOLOGY STATE VARIABLES //

    Array_<EventHandler*>          eventHandlers;
    mutable Array_<EventReporter*> eventReporters;

        // TOPOLOGY CACHE //
    mutable CacheEntryIndex myStateVars;
    mutable CacheEntryIndex myCacheEntries;

};

std::ostream& operator<<(std::ostream& o, const TestSubsystemGuts::StateVars&) {return o;}
std::ostream& operator<<(std::ostream& o, const TestSubsystemGuts::CacheEntries&) {return o;}

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
static Vector weights;
static Vector ooTols;

static bool handleEvents(const System& sys, State& state, Stage g,
                         const Array_<EventId>& triggered) 
{
    if (triggered.empty())
        return false;

    cout << "==> Handling " << triggered.size() << " events at Stage " << g << ":";
    for (unsigned i=0; i < triggered.size(); ++i)
        cout << " " << triggered[i];
    cout << endl;

    Stage lowestModified = Stage::Report;
    bool shouldTerminate = false; 
    sys.handleEvents(state, Event::Cause::Triggered, triggered, accuracy, weights, ooTols, 
                     lowestModified, shouldTerminate);

    if (shouldTerminate) {
        cout << "==> Event at Stage " << g << " requested termination at t=" << state.getTime() << endl;
    }

    return shouldTerminate;
}

void testOne() {
    TestSystem sys;
    TestSubsystem subsys(sys);

    Measure::Zero zero(subsys);
    Measure::Constant three(subsys, 3);
    Measure_<Vec3>::Constant v3const(subsys, Vec3(1,2,3));
    Measure::Sinusoid cos2pit(subsys, 1, 2*Pi, Pi/2);

    // Integrate the cos(2pi*t) measure with IC=0; should give sin(2pi*t)/2pi.
    Measure::Integrate sin2pitOver2pi(subsys, cos2pit, zero);

    Measure::Differentiate dInteg(subsys, sin2pitOver2pi);
    dInteg.setForceUseApproximation(true);

    Measure::Time tMeasure;
    Measure::Time tSubMeas(subsys);

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

    m = zero;
    cout << "m=" << m.getValue(State()) << endl;


    cout << "vplus ref count=" << vplus.getRefCount() << endl;
    cout << "vplus2 ref count=" << vplus2.getRefCount() << endl;


    State state = sys.realizeTopology();
    state.setTime(1.234);
    sys.realize(state, Stage::Time);
    cout << "Initially, tMeasure=" << tMeasure.getValue(state)
         << " tSubMeas=" << tSubMeas.getValue(state) 
         << " 1000*tMeasure=" << t1000.getValue(state)
         << endl;

    Measure_<Mat22>::One m22Ident;
    cout << "Measure_<Mat22>::One=" << m22Ident.getValue(state) << endl;

    cout << "mv after realizeTopo=" << mv.getValue(state) << endl;

    State s2,s3;
    s2 = state; // new copies of variables
    s2 = state; // should do only assignments w/o heap allocation

    // Explicit midpoint steps.
    const Real h = .0001;
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

    // Fill in statics above.
    timescale = sys.calcTimescale(state);
    sys.calcYUnitWeights(state, weights);
    sys.calcYErrUnitTolerances(state, ooTols);

    sys.realize(state, Stage::Acceleration);
    state.autoUpdateDiscreteVariables(); // ??y
    for (int i=0; i <= nSteps; ++i) {

        if (i % outputInterval == 0) {
            sys.realize(state, Stage::Report);
            cout << "\ntMeasure=" << tMeasure.getValue(state)
                 << " d/dt tMeasure=" << tMeasure.getValue(state,1)
                 << " d3/dt3 tMeasure=" << tMeasure.getValue(state,3)
                 << " 1000*tSubMeas=" << t1000.getValue(state)
                 << " t=" << state.getTime() << endl;
            cout << "q=" << state.getQ() << " u=" << state.getU() << endl;
            cout << "qSum=" << subsys.getQSum(state) << " uSum=" << subsys.getUSum(state) << endl;
            cout << "three=" << three.getValue(state) << " v3const=" << v3const.getValue(state) << endl;
            cout << "cos2pit=" << cos2pit.getValue(state) 
                 << " cos(2pi*t)=" << std::cos(2*Pi*state.getTime()) << endl;
            cout << "sin2pitOver2pi=" << sin2pitOver2pi.getValue(state) 
                 << " sin(2pi*t)/2pi=" << std::sin(2*Pi*state.getTime())/(2*Pi) << endl;
            cout << "d/dt sin2pitOver2pi=" 
                 << sin2pitOver2pi.getValue(state,1) << endl;
            cout << "dInteg=" 
                 << dInteg.getValue(state) << endl;
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
        sys.realize(state);

        // Second (final) integrator stage.
        // 1. Unconstrained continuous system.
        const Vector& ydot = state.getYDot();
        state.updY() += h2*ydot; // that is, y = y0 + h*(ydot0+ydot)/2
        state.updTime() += h2;

        // 2. Deal with time-dependent events.
        sys.realize(state, Stage::Time);
        findEvents(state, Stage::Time, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Time, triggered))
            break;

        // 3a. Project position-dependent constraints.
        sys.realize(state, Stage::Position);
        Vector temp;
        sys.project(state, accuracy, weights, ooTols, temp, System::ProjectOptions::PositionOnly);

        // 3b. Handle position-dependent events.
        findEvents(state, Stage::Position, triggers0, triggered);
        if (handleEvents(sys, state, Stage::Position, triggered))
            break;

        // 4a. Project velocity-dependent constraints.
        sys.realize(state, Stage::Velocity);
        sys.project(state, accuracy, weights, ooTols, temp, System::ProjectOptions::VelocityOnly);

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

