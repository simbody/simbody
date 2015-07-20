#ifndef SimTK_SimTKCOMMON_EVENT_TRIGGER_TIMER_H_
#define SimTK_SimTKCOMMON_EVENT_TRIGGER_TIMER_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/EventTrigger.h"

#include <initializer_list>

namespace SimTK {
    
SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventTimerIndex);

//==============================================================================
//                              EVENT TIMER
//==============================================================================
/** An %EventTrigger::Timer is an EventTrigger that causes and Event to occur at
a particular time or times that can be determined before they happen. This 
allows the time stepper to adjust the step size so that it completes a step 
exactly at the time of next event occurrence.

A %Timer provides a method that a time stepper can call with the current
time and state to determine when the next timed event should occur. This class
is still abstract; derive from it to define particular timers. Several 
predefined %Timer classes are provided. **/
class SimTK_SimTKCOMMON_EXPORT EventTrigger::Timer : public EventTrigger {
    using Super = EventTrigger;
public:
    class Periodic;
    class Designated;

    /** Determine the next time at which this timer will trigger. This
    method is intended for use by time steppers in deciding how big a step to
    take next.
    @param      system                
        The System to be used for computations.
    @param      state                
        The current time and state of `system`. This is the earliest possible
        trigger time.
    @param      timeOfLastTrigger    
        We will only return a time strictly later than this one. Set to 
        state.getTime() to prevent re-triggering at the current time. Set
        to -Infinity if you don't want a restriction.
    @returns The next time at which this timer should be trigger. Returns
        Infinity if the timer should not trigger any more. **/
    double calcTimeOfNextTrigger(const System& system, const State& state, 
                                 double timeOfLastTrigger) const {
        double t =
            calcTimeOfNextTriggerVirtual(system, state, timeOfLastTrigger);

        // Treat a returned time earlier than allowed as "don't trigger". This
        // allows for dumber concrete classes. For example, a function that
        // returns just a constant time like "5" will act like a trigger-once
        // timer.
        if (t < state.getTime() || t <= timeOfLastTrigger)
            t = Infinity;
        return t;
    }

    /** (Internal use only) The EventTimerIndex is set after all 
    EventTriggers are known and examined for timers. **/
    EventTimerIndex getEventTimerIndex() const {return m_timerIndex;}
protected:
    /** Create a %Timer with no associated Event objects. A description
    can be useful for reporting and debugging. **/
    explicit Timer(const std::string& triggerDescription) 
    :   Super(triggerDescription) {}

    /** This method must be implemented by any concrete %EventTrigger::Timer. 
    It is always called with `state` realized through Stage::Acceleration. **/
    virtual double calcTimeOfNextTriggerVirtual
       (const System& system, const State& state, 
        double timeOfLastTrigger) const = 0;

private:
friend class SystemGlobalSubsystem;

    // Set when the timers are being collected from within the triggers.
    EventTimerIndex     m_timerIndex;
};



//==============================================================================
//                           DESIGNATED EVENT TIMER
//==============================================================================
/** This is a predefined concrete EventTrigger::Timer that triggers at times 
which are given explicitly. You can initialize the list of times with code like 
this: @code
    EventTrigger::Timer::Designated 
        eventTimes("description", {.1, 1.2, 2, 3.12, 3.4};
    //   OR
    // myTimes is a container of numbers assignment compatible to double,
    // such as std::vector<double> or SimTK::Array_<int>.
    EventTrigger::Timer::Designated 
        eventTimes("description", myTimes.begin(), myTimes.end());
@endcode
You can also insert times after construction using provided methods. Any 
designated times will be sorted in increasing order, and duplicates will be
removed. **/
class SimTK_SimTKCOMMON_EXPORT EventTrigger::Timer::Designated 
:   public EventTrigger::Timer {
    using Super = EventTrigger::Timer;
public:
    /** Create a designated event timer with no times set and no associated
    events yet. A description can be useful for reporting and debugging. **/
    explicit Designated(const std::string& description)
    :   Super(description) {}

    /** Create a designated event timer with a single time set. The time must
    be nonnegative. **/
    Designated(const std::string& description, double t);

    /** Create a designated event timer from a list of specific event times,
    given by begin and one-past-the-end iterators for some container. The 
    container must contain values that are assignment compatible with `double`, 
    such as integers or floats. The times will be sorted and made unique prior 
    to use, at a cost of O(n log n) in the length of the list. An exception will
    be thrown if any of the times are negative. **/
    template <class InputIter>
    Designated(const std::string& description, 
               const InputIter& first, const InputIter& last1)
    :   Super(description), m_triggerTimes(first, last1) {constructHelper();}

    /** Same as previous constructor but takes an `std::initializer_list` as
    input rather than a pair of iterators. **/
    Designated(const std::string&              description, 
               std::initializer_list<double>   times)
    :   Designated(description, times.begin(), times.end()) {} // delegate

    /** Return the number of unique designated signal times contained here. **/
    unsigned getNumDesignatedTimes() const {return m_triggerTimes.size();}
    /** Return one of the designated signal times in sorted order. **/
    double getDesignatedTime(unsigned i) const {return m_triggerTimes[i];}

    /** Designate a time at which this timer should trigger. This
    will be inserted into its proper place in the sorted list of times, and 
    ignored if already present in the list. The time must be nonnegative. Cost
    is constant time (very fast) if `t` is later than any other time in the 
    list, otherwise O(n) in the current length of the list. **/
    void insertDesignatedTime(double t);

    /** Add a list of designated times at which this timer should trigger. The 
    list is provided with a pair of iterators in the usual fashion.
    The times will be inserted into their proper places in the sorted list 
    of times, and ignored if already present in the list. All times must be 
    nonnegative. Cost is O(n log n) to sort the input list, where n is the input
    length, and O(m) to merge it with the already-sorted internal list, where m 
    is the total number of designated times in the two lists. **/
    template <class InputIter>
    void insertDesignatedTimes(const InputIter& first, const InputIter& last1) 
    {   insertHelper(Array_<double>(first, last1)); }

    /** Same as insertDesignatedTimes() above, but takes an 
    `std::initializer_list` rather than a pair of iterators. **/
    void insertDesignatedTimes(std::initializer_list<double> times) 
    {   insertDesignatedTimes(times.begin(), times.end()); }

    /** Clear the current list of designated times for this timer. **/
    void clearDesignatedTimes() {m_triggerTimes.clear();}

private:
    // Check that t is nonnegative and throw an exception if not.
    void checkTime(const char* methodName, double t) const;

    // Assume we have stuffed an unwashed list of times in; we now need to 
    // check that the times are valid, sort them, and throw out duplicates.
    void constructHelper();

    // Given an unwashed list of times in a temporary variable that we can
    // write on, clean these up and insert them properly into the list.
    void insertHelper(Array_<double>&& times);

    Designated* cloneVirtual() const override {return new Designated(*this);}

    // The current implementation takes O(log n) time.
    double calcTimeOfNextTriggerVirtual
       (const System& system, const State& state, 
        double timeOfLastTrigger) const override;

private:
    Array_<double>  m_triggerTimes; // unique, sorted list
};


//==============================================================================
//                            PERIODIC EVENT TIMER
//==============================================================================
/** This is a predefined concrete EventTrigger::Timer that triggers 
periodically, at times n*period. Normally n starts from zero but you may
optionally specify that it should start from one. **/
class SimTK_SimTKCOMMON_EXPORT EventTrigger::Timer::Periodic
:   public EventTrigger::Timer {
    using Super = EventTrigger::Timer;
public:
    /** Create a periodic event trigger and specify the time interval between
    event occurrences.
    @param  description Text describing this particular timer (can be helpful
                        for reporting and debugging).
    @param  period      The time interval between triggerings of this timer. 
                        Must be strictly greater than zero. **/
    Periodic(const std::string& description, double period)
    :   Super(description), m_triggerAtZero(true) {
        setEventPeriod(period);
    }

    /** Create a periodic event trigger and specify the time interval between
    event occurrences. The description will be supplied as "Periodic 
    EventTrigger with period t", where t is replaced with the supplied
    `interval`.
    @param  period  The time interval between triggerings of this timer.  
                    Must be strictly greater than zero. **/
    explicit Periodic(double period)
    :   Super("Periodic EventTrigger with period " 
              + String(period, "%.15g")), m_triggerAtZero(true) {
        setEventPeriod(period);
    }

    /** Prevent this timer from triggering at the start. By default it will 
    trigger at times `n*period` for n starting with n=0; set this flag to 
    make it start at n=1 instead. **/
    Periodic& setShouldTriggerAtZero(bool triggerAtZero)
    {   m_triggerAtZero = triggerAtZero; return *this; }

    /** Return whether this periodic event trigger is set to trigger at
    zero; the default is that it will do so. **/
    bool getShouldTriggerAtZero() const {return m_triggerAtZero;}
    
    /** Get the time interval between triggerings of this timer. This was 
    specified in the constructor or setEventPeriod(). **/   
    double getEventPeriod() const {return m_period;}

    /** Change the event interval for this periodic event timer. **/
    Periodic& setEventPeriod(double period) {
        SimTK_APIARGCHECK1_ALWAYS(period > 0, 
            "EventTrigger::Timer::Periodic", "setEventPeriod", 
            "The period was %g but must be greater than zero.", period);
        m_period = period;
        return *this;
    }

private:
    Periodic* cloneVirtual() const override {return new Periodic(*this);}

    double calcTimeOfNextTriggerVirtual
       (const System& system, const State& state, 
        double timeOfLastTrigger) const override;
private:
    double m_period;
    bool   m_triggerAtZero;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_TRIGGER_TIMER_H_
