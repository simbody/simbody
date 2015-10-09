#ifndef SimTK_SimTKCOMMON_EVENT_WITNESS_H_
#define SimTK_SimTKCOMMON_EVENT_WITNESS_H_

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
/** An %EventWitness is a kind of EventTrigger that provides a scalar 
function of time and state that can be efficiently monitored to detect the 
occurence of an Event. Zero crossings of the witness function cause the 
associated Event(s) to occur. 

A witness may be unidirectional, triggering either on Rising (negative to 
positive) or Falling (positive to negative) zero crossings, or it may be
bidirectional and trigger on any zero crossings.

A witness may also be designated Unilateral, in which case one sign is 
"forbidden" and should always trigger an Event when seen. In that case we expect
the Event Action to cause the witness to reverse direction. The canonical 
example is an impact witness that is a signed distance with negative distance 
indicating forbidden interpenetration; the impact handler is expected to cause a
reversal and prevent the witness from proceeding into the forbidden 
negative-distance zone. Note that a Unilateral witness must also be
unidirectional.

Alternatively, a witness may be Bilateral, meaning that only zero *crossings* 
indicate events; after that the witness value may continue to grow in the same 
direction. To avoid double-triggering, Bilateral witness functions are 
considered to trigger when a transition *to* zero is observed, but not when a 
transition *from* zero is observed. Note that a Bilateral witness can still be 
unidirectional. For examples, a witness designed to isolate the occurrence
of maximum values would trigger on a positive-to-negative change in slope, but
could then continue to become more negative after that since there is nothing
forbidden about negative slope values.

witness functions can be handled most efficiently if they are continuous but
will still work if they are discontinuous. If you know your witness is
discontinuous you can mark it that way which will help the localization 
algorithm somewhat. For continuous witnesses, time derivative information may be 
provided to improve prediction of zero crossings and to speed up localization.

Deadband
--------
"Zero" for a witness means "within the deadband", meaning that the value is
numerically indistinguishable from zero. The width of the deadband
is accuracy-dependent, with the mapping from accuracy to deadband width
determined by the witness and possibly including state dependence. The study
provides the System and State, accuracy, constraint tolerance, and 
acceleration precision. The witness calculates its value and sign -1,0, or 1 
with 0 meaning "the value is in the deadband". If the witness has derivatives,
they have their own deadbands.

This is an abstract class; concrete witness classes are provided with 
particular witness functions defined. **/
class SimTK_SimTKCOMMON_EXPORT EventWitness : public EventTrigger {
    using Super = EventTrigger;
public:
    /** Specify whether witness is expected to take on positive and negative
    values, or just one sign. Set at construction. **/
    enum Range {Bilateral,Unilateral};
    /** Specify the transition direction that causes this witness to trigger
    an event. A unilateral witness must pick one direction. Set at 
    construction. **/
    enum Direction {Rising,Falling,RisingAndFalling};
    /** This is a hint to the localization algorithm about how best to 
    localize the triggering of this witness. Continuous witnesses can be
    localized much more quickly than discontinuous ones. Set at 
    construction. **/
    enum Continuity {Continuous,Discontinuous};

    /** Witnesses trigger events in response to sign transitions of their 
    associated witness function. This enum defines constants for use in 
    specifying which kind of transition has been seen, or which kinds are 
    considered interesting. For the latter purpose, these can be or'ed 
    together to make a mask. **/
    enum TransitionMask : unsigned {
        NoTransition            =0x00,    // must be 0

        // These constants are chosen to be easy to generate from (before,after)
        // pairs; don't change them! (See classifyTransition().)

        // Rising transitions
        NegativeToZero          =0x01, // base 3: 01= bit 1
        NegativeToPositive      =0x02, // 02=bit 2
        ZeroToPositive          =0x10, // 12=bit 5

        // Falling transitions.
        PositiveToNegative      =0x20, // base 3: 20=bit 6
        PositiveToZero          =0x40, // 21=bit 7
        ZeroToNegative          =0x04  // 10=bit 3

        // "or"ed combinations of the above are also considered TransitionMask
        // values
    };

    /** This is the value of a witness function or one of its derivatives.
    Both the actual value and the sign incorporating a finite-width deadband
    are included. **/
    class Value {
    public:
        Value() = default;
        Value(Real value, Real halfDeadband) 
            : m_value(value), m_sign(sign(deadband(value,halfDeadband))) {}
        Real getValue() const {return m_value;}
        int  getSign() const {return m_sign;}
    private:
        Real m_value{NaN};  ///< The actual value calculated.
        int  m_sign{0};     ///< -1,0,1; 0 means value is in the deadband.
    };

    /** Calculate the current value of this %Witness or one of its time 
    derivatives, and determine if it is in the deadband around zero.
    @param          study                
        The Study that is requesting this witness value, containing a reference
        to the System under study. This also supplies accuracy requirements
        that can be used to size the deadband appropriately.
    @param[in]      state
        The State of the Study's System whose time and contents are to be used.
        This will likely be one of the State objects maintained by `study`.
        This State must already be realized to at least Stage 
        `getDependsOnStage(derivOrder)`.
    @param[in]      derivOrder  
        Which derivative value is to be obtained: 0 -> the witness function
        value, 1 -> its 1st time derivative, etc. Must not be higher than the 
        value returned by getNumTimeDerivatives(); defaults to 0.
    @return The function or derivative value (a scalar), and its sign
        characterization -1,0,1 where 0 means "in the deadband." **/                  
    Value calcWitnessValue(const Study&      study,
                           const State&      state, 
                           int               derivOrder=0) const {
        SimTK_ASSERT2(0 <= derivOrder && derivOrder <= getNumTimeDerivatives(),
            "EventWitness::calcWitnessValue(): "
            "The derivative order %d is out of the expected range 0..%d.",
            derivOrder, getNumTimeDerivatives());
        return calcWitnessValueVirtual(study, state, derivOrder);
    }

    /** Given a sign as returned by calcWitnessValue(), determine if it is
    in the forbidden region of a unilateral witness. If the sign is 0 or the
    witness is bilateral, we return false. **/
    bool isSignForbidden(int sign) const {
        assert(std::abs(sign)<=1);
        return sign == m_forbiddenSign;
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
            "EventWitness::getDependsOnStage(): "
            "The derivative order %d is out of the expected range 0..%d.",
            derivOrder, getNumTimeDerivatives());
        return getDependsOnStageVirtual(derivOrder);
    }

    /** Return whether this witness will trigger on a transition from negative
    to non-negative. **/
    bool getTriggerOnRisingSignTransition() const
    {   return m_direction==Rising || m_direction==RisingAndFalling; }

    /** Return whether this witness will trigger on a transition from positive
    to non-positive. **/
    bool getTriggerOnFallingSignTransition() const
    {   return m_direction==Falling || m_direction==RisingAndFalling; }

    /** Set the localization time window as a fraction of Study accuracy.
    During time stepping, this will be multiplied by the accuracy currently in
    effect to produce the absolute window width in time units. The default 
    value is `scale`=0.1 (10%*Accuracy*TimeScale). Set to Infinity to disable 
    accuracy-relative time localization. 
    @see setAbsoluteLocalizationTimeWindow()
    @see System::setDefaultTimeScale() **/
    EventWitness& setAccuracyRelativeTimeLocalizationWindow(Real scale) {
        SimTK_APIARGCHECK1_ALWAYS(scale > 0, "EventWitness",
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
    EventWitness& setAbsoluteTimeLocalizationWindow(double dt) {
        SimTK_APIARGCHECK1_ALWAYS(dt > 0, "EventWitness",
            "setAbsoluteTimeLocalizationWindow",
            "Supplied dt was %g but must be positive.", dt);
        m_absoluteTimeWindow = dt;
        return *this;
    }

    /** Return the current setting of the absolute time localization window, in 
    units of time. **/
    double getAbsoluteTimeLocalizationWindow() const
    {   return m_absoluteTimeWindow; }

    /** (Internal use only) The EventWitnessIndex is set after all 
    EventTriggers are known and examined for witnesses. **/
    EventWitnessIndex getEventWitnessIndex() const {return m_witnessIndex;}

    /** Classify a before/after sign transition. Before and after must both be 
    -1,0, or 1 as returned by the SimTK::sign() function applied to the 
    (possibly deadbanded) witness function value at the beginning and end of a 
    step. **/
    static TransitionMask classifyTransition(int before, int after) {
        assert(std::abs(before)<=1 && std::abs(after)<=1);
        if (after==before) return NoTransition;
        // Treat as 2-digit base 3 number ba to pick out a bit.
        const int whichBit = 3*(before+1) + (after+1); // optimizer cleans up
        return TransitionMask(1 << (whichBit-1));
    }

    /** Given an observed transition, weed out ignorable ones using the supplied
    mask. That is, the return will indicate NoTransition unless the original
    transition was present in the mask. **/
    static TransitionMask maskTransition(TransitionMask transition, 
                                         TransitionMask mask) 
    {
        assert(atMostOneBitIsSet((unsigned)transition));
        // we're depending on NoTransition==0
        return TransitionMask(transition & mask); 
    }

    Range getRange() const {return m_range;}
    Direction getDirection() const {return m_direction;}
    Continuity getContinuity() const {return m_continuity;}
    TransitionMask getTransitionMask() const {return m_transitionMask;}

    /** Translate a Range to a human-readable string. **/
    static std::string toString(Range range);

    /** Translate a Direction to a human-readable string. **/
    static std::string toString(Direction direction);

    /** Translate a Continuity to a human-readable string. **/
    static std::string toString(Continuity continuity);

    /** Translate a TransitionMask into a human-readable string; if there
    are multiple transitions in the mask they will appear or-ed together. **/
    static std::string toString(TransitionMask transitions);


    /** At most we are interested in the witness value, first, and second
    derivatives. We can't make use of any more than that. **/
    enum {MaxDeriv = 2};

protected:
    /** Create a %Witness with no associated Event objects. A description
    can be useful for reporting and debugging. **/
    EventWitness(const std::string& description,
            Range range, Direction direction, Continuity continuity) 
    :   Super(description), 
        m_range(range), m_direction(direction), m_continuity(continuity), 
        m_transitionMask(calcTransitionMask()),
        m_forbiddenSign(calcForbiddenSign()),
        m_accuracyRelativeTimeWindow(0.1), // 10% of System timescale
        m_absoluteTimeWindow(Infinity)
    {
        SimTK_ERRCHK_ALWAYS
           (!(range==Unilateral && direction==RisingAndFalling),
            "EventWitness::EventWitness()",
            "A unilateral witness cannot trigger on both rising and falling "
            "transitions.");
    }

    /** You must override this method to define a witness function and to
    classify its value as in or out of the deadband. You do not need to error 
    check `derivOrder`; the base class will ensure that it is in range 0 to 
    getNumTimeDerivatives(). **/
    virtual Value calcWitnessValueVirtual(const Study&  study,
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

    // Calculate the set of transitions for this witness that should cause an
    // event to trigger. Bilaterals trigger going *to* zero but not *from* zero,
    // in order to avoid double triggering. Unilaterals also trigger going from
    // zero into their forbidden sign.
    TransitionMask calcTransitionMask() const {
        unsigned mask = 0;
        if (m_direction==Rising || m_direction==RisingAndFalling) {
            mask |= (NegativeToZero|NegativeToPositive);
            if (m_range==Unilateral) mask |= ZeroToPositive;
        }
        if (m_direction==Falling || m_direction==RisingAndFalling) {
            mask |= (PositiveToZero|PositiveToNegative);
            if (m_range==Unilateral) mask |= ZeroToNegative;
        }
        return TransitionMask(mask);
    }

    // Return -1 or 1 for unilateral witnesses; some non-sign value for
    // bilaterals so that no sign is forbidden.
    int calcForbiddenSign() const {
        if (m_range==Bilateral) return 2; // not a sign
        assert(m_direction != RisingAndFalling); // not allowed for unilaterals
        return m_direction==Rising ? 1 : -1;
    }

    // These properties are set at construction and can't be changed.
    const Range             m_range;
    const Direction         m_direction;
    const Continuity        m_continuity;

    // These are calculated in the constructor.
    const TransitionMask    m_transitionMask;
    const int               m_forbiddenSign; // -1,1 for unilaterals only

    // Time localization information.
    Real                    m_accuracyRelativeTimeWindow;
    double                  m_absoluteTimeWindow;

    // Set when the witnesses are being collected from within the triggers.
    EventWitnessIndex       m_witnessIndex;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_EVENT_WITNESS_H_
