#ifndef SimTK_SimTKCOMMON_EVENT_TRIGGER_WITNESS_H_
#define SimTK_SimTKCOMMON_EVENT_TRIGGER_WITNESS_H_

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
#include "SimTKcommon/internal/EventTrigger.h"

namespace SimTK {

SimTK_DEFINE_UNIQUE_INDEX_TYPE(EventWitnessIndex);

//==============================================================================
//                               EVENT WITNESS
//==============================================================================
/** An %EventTrigger::Witness is a kind of EventTrigger that provides a scalar 
function of time and state that can be efficiently monitored to detect the 
occurence of a particular condition that causes an Event to occur. Zero 
crossings of this witness function cause the associated Event(s) to occur. 

A %Witness may be defined so that either or both a rising (negative to positive)
or a falling (positive to negative) zero crossing triggers the Event. 
Transitions *to* zero are consider zero "crossings"; transitions *from* zero
are not -- that prevents double-triggering.

%Witness functions are assumed to be continuous unless specifically designated
discontinuous. For continuous witnesses, time derivative information may be 
provided to improve prediction of zero crossings and to speed up localization.
Continuous witnesses may also enforce a localization requirement on the 
witness *value* as well as localization in time. 

This is an abstract class; concrete %Witness classes are provided with 
particular witness functions defined. **/
class SimTK_SimTKCOMMON_EXPORT EventTrigger::Witness : public EventTrigger {
    using Super = EventTrigger;
public:
    /** These are the available strategies for localizing the zero-crossings
    of an event witness. **/
    enum LocalizationStrategy {
        None,               ///< Don't localize; just process when observed.
        Time,               ///< Narrow occurrence to a given time window.
        Value,              ///< Require witness value very near zero.
        BothTimeAndValue,   ///< Require time *and* value localization (strict).
        EitherTimeOrValue   ///< OK if either time *or* value localized (loose).
    };

    /** Designate this witness function as discontinuous. The default is that
    we consider the witness continuous; if that assumption is wrong we will 
    waste time trying to predict and localize its occurrence. Setting this flag
    prevents attempts to predict, changes the localization algorithm to 
    binary-chop in time, and prevents any attempt to localize by value. The
    only allowable values for LocalizationStrategy for a discontinuous witness
    are `None` and `Time`. **/
    void setWitnessIsDiscontinuous(bool isDiscontinuous)
    {   m_witnessIsDiscontinuous = isDiscontinuous; }

    /** Return whether this witness function is known to be discontinuous. 
    Otherwise we are assuming it is continuous. **/
    bool getWitnessIsDiscontinuous() const
    {   return m_witnessIsDiscontinuous; }

    /** Calculate the current value of this %Witness or one of its time 
    derivatives.
    @param          system                
        The System to be used for computations.
    @param[in]      state
        The State of `system` whose time and contents are to be used. Must
        already be realized to at least Stage `getDependsOnStage(derivOrder)`.
    @param[in]      derivOrder  
        Which derivative value is to be obtained: 0 -> the witness function
        value, 1 -> its 1st time derivative, etc. Must not be higher than the 
        value returned by getNumTimeDerivatives(); defaults to 0.
    @return The function or derivative value (a scalar). **/                  
    Real calcWitnessValue(const System&     system,
                          const State&      state, 
                          int               derivOrder=0) const {
        SimTK_ASSERT2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "EventTrigger::Witness::calcWitnessValue(): "
            "The derivative order %d is out of the expected range 0..%d.",
            derivOrder, getNumTimeDerivatives());
        return calcWitnessValueVirtual(system, state, derivOrder);
    }

    /** Every %Witness can produce a value, and some can provide one or more
    total derivatives with respect to time of that value. This method reports 
    how many are available: 1 -> first derivative d/dt is available, 2 -> first
    and second derivative d^2/dt^2 are available, etc. We interpret the "0th 
    derivative" to be the witness function's value.
    @return The maximum available derivative order. **/
    int getNumTimeDerivatives() const 
    {   return getNumTimeDerivativesVirtual(); }

    /** At what Stage can we expect the value of this witness function or one
    of its time derivatives to be available?
    @param[in]      derivOrder
        Which derivative level is to be checked: 0 -> the value, 1 -> the 1st 
        time derivative, etc. Must not be higher than the value returned by 
        getNumTimeDerivatives() and defaults to 0.
    @return The Stage after which this value or derivative is available. **/                  
    Stage getDependsOnStage(int derivOrder=0) const {
        SimTK_ASSERT2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "EventTrigger::Witness::calcWitnessValue(): "
            "The derivative order %d is out of the expected range 0..%d.",
            derivOrder, getNumTimeDerivatives());
        return getDependsOnStageVirtual(derivOrder);
    }

    /** Specify whether a rising (negative to positive) transition of the 
    witness function should trigger the associated event(s). **/
    Witness& setTriggerOnRisingSignTransition(bool triggerOnRising) {
        m_triggerOnRising = triggerOnRising;
        return *this;
    }

    /** Return whether this witness will trigger on a transition from negative
    to non-negative. **/
    bool getTriggerOnRisingSignTransition() const
    {   return m_triggerOnRising; }

    /** Specify whether a falling (positive to negative) transition of the 
    witness function should trigger the associated event(s). **/
    Witness& setTriggerOnFallingSignTransition(bool triggerOnFalling) {
        m_triggerOnFalling = triggerOnFalling;
        return *this;
    }

    /** Return whether this witness will trigger on a transition from positive
    to non-positive. **/
    bool getTriggerOnFallingSignTransition() const
    {   return m_triggerOnFalling; }

    /** Specify the strategy to be used when localizing the zero crossing of
    this witness. The choices are: don't localize, localize by time window, 
    localize by value, require both, or either one. For discontinous witness
    functions only no localization or localization by time are allowed. **/
    Witness& setLocalizationStrategy(LocalizationStrategy strategy) 
    {   m_strategy = strategy; return *this; }

    /** Return the LocalizationStrategy to be used for this witness function.
    For discontinuous witness functions only `None` and `Time` are effective;
    we'll treat anything but `None` as though it were `Time`. **/
    LocalizationStrategy getLocalizationStrategy() const 
    {   return m_strategy; }

    /** Set the localization time window as a fraction of Study accuracy.
    During time stepping, this will be multiplied by the accuracy currently in
    effect to produce the absolute window width in time units. The default 
    value is `scale`=0.1 (10%*Accuracy*TimeScale). Set to Infinity to disable 
    accuracy-relative time localization. 
    @see setAbsoluteLocalizationTimeWindow()
    @see System::setDefaultTimeScale() **/
    Witness& setAccuracyRelativeTimeLocalizationWindow(Real scale) {
        SimTK_APIARGCHECK1_ALWAYS(scale > 0, "EventTrigger::Witness",
            "setAccuracyRelativeTimeLocalizationWindow",
            "Supplied scale was %g but must be positive.", scale);
        m_accuracyRelativeTimeWindow = scale;
        return *this;
    }

    /** Return the current setting of the accuracy-relative time localization
    window, as a fraction of Study accuracy. **/
    Real getAccuracyRelativeTimeLocalizationWindow() const 
    {   return m_accuracyRelativeTimeWindow; }

    /** Set the maximum tolerance `dt` for time localization of a reported zero
    crossing, in units of time. This applies regardless of the accuracy 
    currently in use by the time stepper. If time localization is enabled,
    events will be localized such that the actual event time is within `dt` of 
    the reported time of occurrence. The default value is `dt=Infinity`, that 
    is, we do not use an absolute tolerance but rather use a tolerance that
    is scaled by the requested integration accuracy. 
    @see setAccuracyRelativeLocalizationTimeWindow() **/    
    Witness& setAbsoluteTimeLocalizationWindow(double dt) {
        SimTK_APIARGCHECK1_ALWAYS(dt > 0, "EventTrigger::Witness",
            "setAbsoluteTimeLocalizationWindow",
            "Supplied dt was %g but must be positive.", dt);
        m_absoluteTimeWindow = dt;
        return *this;
    }

    /** Return the current setting of the absolute time localization window, in 
    units of time. **/
    double getAbsoluteTimeLocalizationWindow() const
    {   return m_absoluteTimeWindow; }

    /** Set the scale factor for determining the required tolerance when
    isolating a zero crossing of the witness function value. During time
    stepping, this scale factor will be multiplied by the accuracy currently in
    effect to produce the absolute tolerance in units of the witness 
    function's value. Events will be localized such that the witness function
    value is no more than +/- tolerance at the time the event is reported as
    having occurred. The default value is `scale`=0.1 (10%*Accuracy). Set to 
    Infinity to disable accuracy-relative time localization. **/
    Witness& setAccuracyRelativeValueLocalizationTolerance(Real scale) {
        SimTK_APIARGCHECK1_ALWAYS(scale > 0, "EventTrigger::Witness",
            "setAccuracyRelativeValueLocalizationTolerance",
            "Supplied scale was %g but must be positive.", scale);
        m_accuracyRelativeValueTolerance = scale;
        return *this;
    }

    /** Return the current setting of the accuracy-relative value localization
    tolerance, in units of the witness value but scaled by Study accuracy. **/
    Real getAccuracyRelativeValueLocalizationTolerance() const 
    {   return m_accuracyRelativeValueTolerance; }

    /** Set the absolute tolerance for error in the value of the witness 
    function at the reported zero crossing, in units of the witness function
    value. Events will be localized such that the witness function value is no
    more than +/- `tolerance` at the time the event is reported as having 
    occurred. The default value is `tolerance=Infinity`, that is, we do not use 
    an absolute tolerance. **/
    Witness& setAbsoluteValueLocalizationTolerance(Real tolerance) {
        SimTK_APIARGCHECK1_ALWAYS(tolerance > 0, "EventTrigger::Witness",
            "setAbsoluteValueLocalizationTolerance",
            "Supplied tolerance was %g but must be positive.", tolerance);
        m_absoluteValueTolerance = tolerance;
        return *this;
    }

    /** Return the current setting of the absolute tolerance for error in the
    value of the witness function at a reported zero crossing, in units of the
    witness function value. **/
    double getAbsoluteValueLocalizationTolerance() const
    {   return m_absoluteValueTolerance; }


    /** Return a bitmask for selecting sign transitions that are significant
    for this witness function. This will return Falling, Rising, or 
    AnySignChange=Falling|Rising. **/
    Event::TriggerDirection calcTransitionMask() const {
        unsigned mask = 0;
        if (m_triggerOnRising) {
            mask |= Event::NegativeToPositive;
        }
        if (m_triggerOnFalling) {
            mask |= Event::PositiveToNegative;
        }
        return Event::TriggerDirection(mask);
    }

    /** (Internal use only) The EventWitnessIndex is set after all 
    EventTriggers are known and examined for witnesses. **/
    EventWitnessIndex getEventWitnessIndex() const {return m_witnessIndex;}

    /** (Internal use only) This returns the assigned slot in the value array
    for the stage associated with the value, first derivative, or second
    derivative. **/
    int getValueIndex(int derivOrder) const {
        assert(derivOrder <= getNumTimeDerivatives());
        return m_valueIndex[derivOrder];
    }

    /** At most we are interested in the witness value, first, and second
    derivatives. We can't make use of any more than that. **/
    enum {MaxDeriv = 2};

protected:
    /** Create a %Witness with no associated Event objects. A description
    can be useful for reporting and debugging. **/
    explicit Witness(const std::string& description) 
    :   Super(description) {setDefaults();} 

    /** You must override this method to define a witness function. You do not
    need to error check `derivOrder`; the base class will ensure that it is in 
    range 0 to getNumTimeDerivatives(). **/
    virtual Real calcWitnessValueVirtual(const System& system,
                                         const State&  state, 
                                         int           derivOrder) const = 0; 

    /** You must override this method to provide the Stage dependency of you
    witness function so that the solver can avoid calling it until the State
    has been sufficiently realized. The dependency level will generally be
    different for the value and each of its time derivatives; `derivOrder`
    will be passed in as zero for the value's Stage, 1 for the 1st derivative's
    Stage, and so on. You do not need to error check `derivOrder`; the base 
    class will ensure that it is in range 0 to getNumTimeDerivatives(). **/
    virtual Stage getDependsOnStageVirtual(int derivOrder) const = 0;

    /** If you can provide derivatives, override this method to say how many.
    One is very helpful, two can help a little more, but any higher derivatives
    will be ignored. The default implementation returns 0, meaning that only 
    the witness function value is available but no derivatives. **/
    virtual int getNumTimeDerivativesVirtual() const {return 0;}

private:
friend class SystemGlobalSubsystem;
    void setDefaults() {
        m_triggerOnRising = m_triggerOnFalling = true;
        m_witnessIsDiscontinuous = false;
        m_strategy = BothTimeAndValue;
        m_accuracyRelativeTimeWindow     = 0.1; // 10% of System timescale
        m_accuracyRelativeValueTolerance = 0.1; // 10% of value
        m_absoluteTimeWindow     = Infinity;
        m_absoluteValueTolerance = Infinity;

        m_witnessIndex.invalidate();
        for (int i=0; i <= MaxDeriv; ++i) m_valueIndex[i] = -1; 
    }

    bool                    m_triggerOnRising;
    bool                    m_triggerOnFalling;
    bool                    m_witnessIsDiscontinuous;
    LocalizationStrategy    m_strategy;

    // Localization information
    Real                    m_accuracyRelativeTimeWindow;
    double                  m_absoluteTimeWindow;

    Real                    m_accuracyRelativeValueTolerance;
    Real                    m_absoluteValueTolerance;

    // Set when the witnesses are being collected from within the triggers.
    EventWitnessIndex       m_witnessIndex;

    // These are the assigned slots in the per-stage value arrays, for the
    // slot that goes with the value, first derivative, and second derivative.
    int                     m_valueIndex[1+MaxDeriv];
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_WITNESS_H_
