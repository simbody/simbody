#ifndef SimTK_SIMBODY_ASSEMBLER_H_
#define SimTK_SIMBODY_ASSEMBLER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-15 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <set>
#include <map>
#include <cassert>
#include <cmath>

namespace SimTK {

SimTK_DEFINE_UNIQUE_INDEX_TYPE(AssemblyConditionIndex);

class AssemblyCondition;


//==============================================================================
//                               ASSEMBLER
//==============================================================================
/** This Study attempts to find a configuration (set of joint coordinates q) 
of a Simbody MultibodySystem that satisfies the System's position Constraints
plus optional additional assembly conditions. If successful, the final set
of q's will satisfy the constraints to within a specified tolerance. The
Assembler also supports high-performance repeated assembly (also known as
"inverse kinematics" or "tracking") where only small changes are expected 
between a series of observation frames.

The complete specification for an Assembly study consists of four elements:
  - The subset of q's which may be modified by the study.
  - Limits on the allowable range of values for each q.
  - A set of assembly error conditions that \e must be satisfied.
  - A set of weighted assembly goals that are to be achieved as best we can.

By default, all q's may be modified with no range restrictions. The q's whose
values are a prescribed function of time will be set to those values, while free
q's are available for satisfying the assembly error conditions and goals. The 
assembly error conditions are just the errors in the unconditional position 
(holonomic) constraints that are present in the MultibodySystem and currently 
enabled, as well as the conditional (unilateral) holonomic constraints that
are both enabled and active *on entry*. (Quaternion normalization constraints 
will also be satisfied, but do not generate assembly errors.) There are no 
default assembly goals. This is very similar in behavior to the System's 
project() method except that project() considers it an error if the constraints 
aren't already close to being satisfied initially, while Assembler will attempt 
to satisfy them regardless, and may take a series of increasingly desperate 
measures to do so. 

@bug Currently the %Assembler makes no attempt to activate or deactivate
conditional constraints; it will simply work with them as it finds them in
the initial state. If none are active on entry, they will be completely ignored.

<h2>Basic assembly:</h2>
This is the most common use of the Assembler: modify a System's given State
so that its configuration (set of generalized coordinates q) satisfies the 
position-affecting Motion and Constraint objects in the System that are
currently enabled in that State. This is done to a default tolerance if you 
don't provide one, and that tolerance is suitable for use with subsequent 
dynamic studies that are run at their default tolerances. Although the assembly 
begins from the configuration provided in the initial state, no assumption is 
made about how close this initial configuration is to one that satisfies the 
assembly conditions.
@code
  MultibodySystem system;
  // ... build system; get initial state

  Assembler assembler(system); // construct the Assembler study object
  try // modify state to satisfy Constraints
  {   assembler.assemble(state); }
  catch (const std::exception& exc)
  {   std::cout << "Assembly failed: " << exc.what() << std::endl; }
@endcode

<h2>Inverse kinematics (repeated assembly):</h2>
After the initial assembly done as above, an inverse kinematic study consists
of a series of assembly solutions for a sequence of small changes to the 
assembly conditions. A common example is the tracking of a time series of
marker observations. Thus each "tracking" assembly computation may assume:
  - There has been no change to the problem structure.
  - The internal State is already close to the desired solution and has not
    been changed since the last tracking frame solution.
  - The set of active conditional constraints has not changed.

Allowable (gradual) changes between tracking frames are:
  - Time, which may affect prescribed motion or other time-dependent 
    Constraints as well as time-dependent assembly conditions.
  - Weights on assembly goals. Continuous changes to goal weights are 
    permitted, but changes between 0 and nonzero or between finite and
    Infinity should not be made because these are structural changes
    to the form of the problem and require reinitialization.
  - Any changes allowed by the assembly conditions currently included in this
    Assembler, for example marker location observations for the Markers
    assembly condition.

The track() method performs assembly analysis under these assumptions and can
be much faster than assemble(). Here is an outline of code that performs
repeated tracking of data from a series of observation frames, each associated
with a frame time:
@code
  MultibodySystem system;
  // ... build system; get initial state
  Assembler assembler(system); // construct the Assembler Study object
  // ... set up assembly conditions; perform initial assemble() as above;
  //     assume assembled result is in State myState.

  try // track a series of small changes to the assembly conditions
  {   for (int i=0; i < numFrames; ++i) {
          // ... update assembly conditions for frame[i]
          assembler.track(frameTime[i]);
          assembler.updateFromInternalState(myState); // update time and qs
          // ... do something with the results in myState
      }
  }
  catch (const std::exception& exc)
  {   std::cout << "Tracking failed: " << exc.what() << std::endl; }
@endcode

<h2>Optional settings:</h2>
Optional settings include:
  - Locking particular mobilizers so that their q's can't be changed.
  - Setting bounds on the acceptable range of values for some of the q's. 
  - Defining additional assembly conditions.
  - Assigning assembly conditions to be errors ("needs") or weighted goals 
    ("wants").

Assembly errors are specified by giving an assembly condition a weight of 
Infinity. Anything with a lower weight is a goal and will be combined with all
the other goals into a single scalar objective. The built-in Constraints are
normally treated with infinite weight, but you can change them to goals instead
if you like; sometimes that can be useful as a step in getting a difficult-to-
assemble system assembled. 
**/
class SimTK_SIMBODY_EXPORT Assembler : public Study {
    typedef std::set<MobilizedBodyIndex>            LockedMobilizers;
    typedef std::set<MobilizerQIndex>               QSet;
    typedef std::map<MobilizedBodyIndex, QSet>      LockedQs;
    typedef std::map<MobilizerQIndex, Vec2>         QRanges;
    typedef std::map<MobilizedBodyIndex, QRanges>   RestrictedQs;
public:

/** Assembler::FreeQIndex is a unique integer type used for accessing
the subset of q's that the Assembler is permitted to change. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Assembler,FreeQIndex);

/** This class is the exception object thrown when a call to assemble()
is unable to achieve the required error tolerance. It is derived from 
std::exception so details can be obtained via the what() method. **/
class AssembleFailed;
/** This class is the exception object thrown when a call to track()
is unable to achieve the required error tolerance. It is derived from 
std::exception so details can be obtained via the what() method. **/
class TrackFailed;

/** @name             Construction and setup
By default, the only assembly condition is that any Simbody Constraints
in the System that are enabled, and conditional constraints that are both
enabled and active, must be satisifed to within the assembly
tolerance. (Prescribed motions specified with Motion objects are satisfied
exactly.) You can selectively enable and disable Constraints in the
state using the ordinary Constraint::disable() and enable() methods. You
can also apply an overall weighting to these Constraints here if you want; 
if the weight is zero they will be ignored; if Infinity they are treated
as must-satisfy assembly error conditions; any other number is used to
weight the RMS Constraint error into the scalar objective along
with the other goals. Additional assembly conditions may be specified with
these methods; some predefined conditions are available, most notably 
Markers for tracking observed marker locations. **/
/*@{*/
/** Create an Assembler study for the given MultibodySystem. The
Assembler's current state is set to the System's default state but with
Euler angles used instead of quaternions. **/
explicit Assembler(const MultibodySystem& system);

/** Set the assembly error tolerance. This value is tested against a norm
of all the assembly error conditions to determine whether an assemble() or
track() operation was successful. Note that assembly errors may have
arbitrary units (for example, built in Constraint errors may be distances
or angles); in general they must be scaled so that the same tolerance
value can be used for all of them. By default, tolerance is set to 
accuracy/10 if accuracy has been set, otherwise 1e-4; calling 
setErrorTolerance() with no argument or with zero restores it to its default 
behavior. **/
Assembler& setErrorTolerance(Real tolerance=0) {
    SimTK_ERRCHK1_ALWAYS(0 <= tolerance,
        "Assembler::setTolerance()", "The requested error tolerance %g"
        " is illegal; we require 0 <= tolerance, with 0 indicating that"
        " the default tolerance (accuracy/10) is to be used.", tolerance);
    this->tolerance = tolerance;
    return *this;
}
/** Obtain the tolerance setting that will be used during the next 
assemble() or track() call. Note that this may be an explicitly-set
tolerance or a default value calculated as accuracy/10 if accuracy has been
set, otherwise 1e-4. **/
Real getErrorToleranceInUse() const {   
    return tolerance > 0 ? tolerance 
           : (accuracy > 0 ? accuracy/10 : Real(0.1)/OODefaultAccuracy); 
}

/** Set the accuracy to which a solution should be pursued. This is a
unitless value that is roughly interpreted as a request for a certain
number of "correct" digits in the "answer", so that an accuracy of 0.001
means "1/10 of 1 percent" or roughly three digits, meaning that we would
like the assembly to stop when the goal is within 0.1% of its minimum.
However, if you don't say otherwise, this number is also used to set the 
absolute error tolerance used to determine whether the assembly succeeded 
or failed, by the following formula: error tolerance = accuracy/10. By
default, we set accuracy=1e-3 and tolerance=1e-4. **/
Assembler& setAccuracy(Real accuracy=0) {
    SimTK_ERRCHK2_ALWAYS(0 <= accuracy && accuracy < 1,
        "Assembler::setAccuracy()", "The requested accuracy %g is illegal;"
        " we require 0 <= accuracy < 1, with 0 indicating that the default"
        " accuracy (%g) is to be used.", Real(1)/OODefaultAccuracy, accuracy);
    this->accuracy = accuracy;
    return *this;
}
/** Obtain the accuracy setting that will be used during the next 
assemble() or track() call. The default is to use 1e-3, i.e., 1/10 of 1%. **/
Real getAccuracyInUse() const 
{   return accuracy > 0 ? accuracy : Real(1)/OODefaultAccuracy; }


/** Change how the System's enabled built-in Constraints are weighted as
compared to other assembly conditions. If this is Infinity (the default) then
the built-ins are treated as must-satisfy constraints; otherwise they are 
included in the assembly cost function with the given weight. If the weight is 
given as zero the built-in Constraints will be ignored altogether.
@see setAssemblyConditionWeight() **/
Assembler& setSystemConstraintsWeight(Real weight)
{   assert(systemConstraints.isValid());
    setAssemblyConditionWeight(systemConstraints,weight);
    return *this; }

/** Return the current weight being given to the System's built-in
Constraints; the default is Infinity. 
@see getAssemblyConditionWeight() **/
Real getSystemConstraintsWeight() const
{   assert(systemConstraints.isValid());
    return getAssemblyConditionWeight(systemConstraints); }

/** Set the weight to be used for this AssemblyCondition. If the weight is set
to 0, this condition will be disabled and will be ignored. If the weight is 
set to Infinity, the condition will be treated as an assembly error condition
that must be satisfied to tolerance. Otherwise (finite weight) the condition
will be treated as an assembly goal and the weight will be used to combine its 
cost function with that of the other assembly goals. **/
Assembler& setAssemblyConditionWeight(AssemblyConditionIndex condition, 
                                      Real                   weight) {
    SimTK_INDEXCHECK_ALWAYS(condition, conditions.size(),
        "Assembler::setAssemblyConditionWeight()");
    SimTK_ERRCHK1_ALWAYS(weight >= 0, "Assembler::setAssemblyConditionWeight()",
        "Illegal weight %g; weight must be nonnegative.", weight);
    uninitialize();
    weights[condition] = weight;
    return *this;
}

/** Return the weight currently in use for this AssemblyCondition. If the
returned value is 0, this condition is being ignored. If the weight is 
Infinity, then the condition is being treated as an assembly error condition
that must be satisfied to tolerance. Otherwise (finite weight) this is an
assembly goal and the weight is used to combine its cost function with that
of the other assembly goals. **/
Real getAssemblyConditionWeight(AssemblyConditionIndex condition) const {
    SimTK_INDEXCHECK_ALWAYS(condition, conditions.size(),
        "Assembler::getAssemblyConditionWeight()");
    return weights[condition];
}

/** Add an assembly error condition to this Assembler study, taking over 
ownership of the heap-allocated AssemblyCondition object. We will use the
calcErrors() method of this object to determine errors whose norm \e must be
driven below tolerance for an assembly to be considered successful. **/
AssemblyConditionIndex 
    adoptAssemblyError(AssemblyCondition* p);
/** Add an assembly goal to this Assembler study, taking over ownership
of the heap-allocated AssemblyCondition object. We will use normally use the 
calcGoal() method of this object to calculate its contribution to the 
assembly goal cost function. An optional weight can be provided that is used
when combining this cost with those of other goals to form the overall cost
function; the default weight is 1. If the weight is 0 the goal is ignored and 
not evaluated at all; if the weight is Infinity this is actually an assembly
constraint and we'll use its calcErrors() method instead. **/
AssemblyConditionIndex 
    adoptAssemblyGoal(AssemblyCondition* p, Real weight=1);


/** Set the Assembler's internal state from an existing state which must
be suitable for use with the Assembler's System as supplied at the time
the Assembler was constructed. All variables are copied, not just q's, so
the Assembler must be reinitialized after this call in case modeling 
options, instance variables, or time have changed. **/
Assembler& setInternalState(const State& state) {
    uninitialize();
    getMatterSubsystem().convertToEulerAngles(state, internalState);
    system.realizeModel(internalState);
    return *this;
}
/** Initialize the Assembler to prepare for performing assembly analysis.
This is normally called automatically when assemble() is called, but you 
can call it explicitly and then access methods that report on the 
properties of the system on which the analysis will be performed. The 
internal state should already have been set; if you want to provide the 
state now use initialize(State). **/
void initialize() const;
/** Set the internal State and initialize. See setInternalState() and 
initialize() methods for more information. **/
void initialize(const State& state)
{   setInternalState(state); initialize(); }
/*@}*/

/** @name                    Execution
These methods perform assembly or tracking analysis, determine how
successful they were, and obtain results. **/
/*@{*/

/** Starting with the current value of the internally-maintained State, 
modify the q's in it to satisfy all the assembly conditions to within a 
tolerance. The actual tolerance achieved is returned as the function value. 
@return The goal value actually achieved (not the error norm; that is 
guaranteed to be no greater than the error tolerance if this returns at
all.  **/
Real assemble();

/** Continue a series of assembly steps that is already in progress,
without restarting or reanalyzing the system, and optionally providing
a new frame time. This is designed for use with a series of assembly 
frames that are close together so that no heroic measures are needed to 
go from one to the next. For the first frame, and any time there might be 
a change to the problem structure or a major change to the state, use 
assemble() instead of track(). See the Assembler class documentation for 
more information and usage examples. **/
Real track(Real frameTime = -1);

/** Given an initial value for the State, modify the q's in it to satisfy
all the assembly conditions to within a tolerance. The actual tolerance 
achieved is returned as the function value. 
@param[in,out]      state    
    The initial and final State value. Only q's are modified.
@return The tolerance actually achieved.  **/
Real assemble(State& state) {
    setInternalState(state);
    Real achievedCost = assemble(); // throws if it fails
    updateFromInternalState(state);
    return achievedCost;
}


/** Return the goal value attained by the internal State's current settings
for the free q's; this is a weighted sum of the individual goal values for 
each assembly goal. Goal values are nonnegative scalars. **/
Real calcCurrentGoal() const;
/** This is the weighted norm of the assembly constraint errors directly
comparable with the assembly error tolerance setting. That is, if this 
number is less than or equal to tolerance (as returned by
getErrorToleranceInUse()), then the current state is a feasible
assembly solution (although it may not be optimal). Note that by default
we use the infinity norm (maximum absolute value of any error term) but
that you can specify use of an RMS norm instead via setUseRMSErrorNorm().
@see getErrorToleranceInUse(), setUseRMSErrorNorm() **/
Real calcCurrentErrorNorm() const;


/** Given an existing State that is suitable for the Assembler's System, 
update its q's from those found in the Assembler's internal State, leaving 
everything else unchanged. We will convert from Euler angles to quaternions
if the destination State is set to use quaternions. **/
void updateFromInternalState(State& state) const {
    system.realizeModel(state); // allocates q's if they haven't been yet
    if (!getMatterSubsystem().getUseEulerAngles(state)) {
        State tempState;
        getMatterSubsystem().convertToQuaternions(getInternalState(),
                                                  tempState);
        state.updQ() = tempState.getQ();
    } else 
        state.updQ() = getInternalState().getQ();
}
/*@}*/

/** @name                Parameter restrictions
These methods restrict which q's are allowed to be modified while trying
to assemble the system, or restrict the range within which the final q's
must lie. A prescribed mobilizer is always treated as locked; its q's
are set to their prescribed values and are not changed to satisfy assembly
conditions. **/
/*@{*/

/** Lock this mobilizer at its starting position. This overrides any 
individual q specifications, so even if a q was specifically unlocked it
will not move until the mobilizer as a whole is unlocked. **/
void lockMobilizer(MobilizedBodyIndex mbx)
{   uninitialize(); userLockedMobilizers.insert(mbx); }
/** Unlock this mobilizer as a whole; some of its q's may remain locked
if they were locked individually. It is OK if this mobilizer was already
unlocked; in that case this does nothing. Attempts to unlock a prescribed
mobilizer have no effect; you must disable the mobilizer's Motion object
instead. **/
void unlockMobilizer(MobilizedBodyIndex mbx) 
{   uninitialize(); userLockedMobilizers.erase(mbx); }

/** Lock one of this mobilizer's q's at its initial value. Be careful with 
this method because it requires that you understand the order of the 
generalized coordinates used by this particular mobilizer during assembly. 
In particular, the mobilizer will be modeled with Euler angles rather than 
quaternions and you must know the Euler sequence it uses in that case (that
is, body- or space-fixed, 2 or 3 axes, and the rotation order). It is 
preferable to use lockMobilizer() instead since that will lock all the q's 
however they are defined. Note that locking individual q's with this method
is independent of whole-mobilizer locking. If you unlock the mobilizer with
unlockMobilizer(), any q's which have been explicitly locked with lockQ() 
will remain locked. **/
void lockQ(MobilizedBodyIndex mbx, MobilizerQIndex qx)
{   uninitialize(); userLockedQs[mbx].insert(qx); }

/** Unlock one of this mobilizer's q's if it was locked. Note that this will
not take effect immediately if the mobilizer as a whole has been locked with 
lockMobilizer(); you have to unlockMobilizer() first. This has no effect on
a prescribed q. **/
void unlockQ(MobilizedBodyIndex mbx, MobilizerQIndex qx)
{   LockedQs::iterator p = userLockedQs.find(mbx);
    if (p == userLockedQs.end()) return;
    QSet& qs = p->second;
    if (qs.erase(qx)) { // returns 0 if nothing erased
        uninitialize();
        if (qs.empty())
            userLockedQs.erase(p); // remove the whole mobilized body
    }
}

/** Restrict a q to remain within a given range. Caution: this requires that
you understand the order of the generalized coordinates used by this 
particular mobilizer during assembly; see lockQ() for a discussion. You can
use -Infinity or Infinity to indicate that the q is not bounded in one 
direction. This has no effect on a prescribed q. **/
void restrictQ(MobilizedBodyIndex mbx, MobilizerQIndex qx,
               Real lowerBound, Real upperBound)
{   SimTK_ERRCHK2_ALWAYS(lowerBound <= upperBound, "Assembler::restrictQ()", 
        "The given range [%g,%g] is illegal because the lower bound is"
        " greater than the upper bound.", lowerBound, upperBound);
    if (lowerBound == -Infinity && upperBound == Infinity)
    {   unrestrictQ(mbx,qx); return; }
    uninitialize(); 
    userRestrictedQs[mbx][qx] = Vec2(lowerBound,upperBound); 
}


/** Unrestrict a particular generalized coordinate q if it was previously 
restricted. Note that this is independent of whether the q has been locked 
with lockMobilizer() or lockQ(); that is, the q may still be locked even 
though it is now unrestricted. This has no effect on a prescribed q. **/
void unrestrictQ(MobilizedBodyIndex mbx, MobilizerQIndex qx)
{   RestrictedQs::iterator p = userRestrictedQs.find(mbx);
    if (p == userRestrictedQs.end()) return;
    QRanges& qranges = p->second;
    if (qranges.erase(qx)) { // returns 0 if nothing erased
        uninitialize();
        if (qranges.empty())
            userRestrictedQs.erase(p); // remove the whole mobilized body
    }
}
/*@}*/



/** @name                    Statistics
The Assembler keeps counters of various internal operations it performs
during execution; these methods access those counters. These can be helpful
in evaluating the effects of various ways of structuring the assembly or 
tracking problem. Counters can also be reset to zero manually by calling 
resetStats(). **/
/*@{*/
/** Return the number of goal evaluations. **/
int getNumGoalEvals()  const;
/** Return the number of assembly error condition evaluations. **/
int getNumErrorEvals() const;
/** Return the number of goal gradient evaluations. **/
int getNumGoalGradientEvals()   const;
/** Return the number of assembly error condition Jacobian evaluations. **/
int getNumErrorJacobianEvals()   const;
/** Return the number of assembly steps; that is, the number of calls to
assemble() or track() since last initialization. **/
int getNumAssemblySteps() const;
/** Return the number of system initializations performed since this
Assembler was created or the most recent resetStats() call. **/
int getNumInitializations() const;
/** Reset all counters to zero; except for the number of initializations
counter this also happens whenever the assembler system is reinitialized 
either explicitly or due to system changes. **/
void resetStats() const;
/*@}*/


/** @name                Advanced options
These are primarily useful for debugging while developing new 
AssemblyCondition classes. **/
/*@{*/

/** This is useful for debugging but should not be used otherwise
since the analytic gradient is to be preferred. **/
void setForceNumericalGradient(bool yesno)
{   forceNumericalGradient = yesno; }
/** This is useful for debugging but should not be used otherwise
since the analytic Jacobian is to be preferred. **/
void setForceNumericalJacobian(bool yesno)
{   forceNumericalJacobian = yesno; }

/** Use an RMS norm for the assembly errors rather than the default
infinity norm (max absolute value). RMS is less stringent and defines
success based on on a good "average" case rather than a good worst case.
If there are n error terms ei, the default norm is e=max_i(abs(ei)) and
the RMS norm is e=sqrt(sum_i(ei^2)/n). **/
void setUseRMSErrorNorm(bool yesno)
{   useRMSErrorNorm = yesno; }
/** Determine whether we are currently using the RMS norm for constraint
errors; if not we're using the default infinity norm (max absolute value).
**/
bool isUsingRMSErrorNorm() const {return useRMSErrorNorm;}

/** Uninitialize the Assembler. After this call the Assembler must be
initialized again before an assembly study can be performed. Normally this
is called automatically when changes are made; you can call it explicitly
if you want. **/
void uninitialize() const;
/** Check whether the Assembler has been initialized since the last change
was made to its contents. **/
bool isInitialized() const {return alreadyInitialized;}

/** This provides read-only access to the Assembler's internal State; you
probably should use updateFromInternalState() to transfer just q's from
the internal state to your own State. Be aware that the internal state is
always maintained using Euler angles for rotations rather than quaternions,
while updateFromInternalState() will make sure you get the rotations in the
form you want. **/
const State& getInternalState() const {return internalState;}

/** Given a reference to an EventReporter, use this Reporter to provide 
progress reporting. The EventReporter object must be owned by someone
else and persist throughout the lifetime of this Assembler object. **/
void addReporter(const EventReporter& reporter) {
    reporters.push_back(&reporter);
}

/** Return the number of q's which are free to be changed by this 
already-initialized assembly analysis. The rest of the q's are locked
at their initial values or set to their prescribed values. **/
int getNumFreeQs() const 
{   return freeQ2Q.size(); }

/** Return the absolute q index associated with a free q. Every free q
is associated with a q so this will always return a valid index if the
free q index is in range. **/
QIndex getQIndexOfFreeQ(FreeQIndex freeQIndex) const
{   return freeQ2Q[freeQIndex]; }

/** A subset of the q's will be used as free q's for solving the assembly
problem. Given an absolute q index, this will return the corresponding
free q index if there is one; otherwise, the returned index will be 
invalid meaning that this q is currently locked. **/
FreeQIndex getFreeQIndexOfQ(QIndex qx) const 
{   return q2FreeQ[qx]; }

/** Return the allowable range for a particular free q. If this free q is
unrestricted the returned range will be [-Infinity,Infinity]. **/
Vec2 getFreeQBounds(FreeQIndex freeQIndex) const {
    if (!lower.size()) return Vec2(-Infinity, Infinity);
    else return Vec2(lower[freeQIndex], upper[freeQIndex]);
}

/** Return a reference to the MultibodySystem associated with this 
Assembler (that is, the System that was supplied in the Assembler's
constructor. **/
const MultibodySystem& getMultibodySystem() const 
{   return system; }
/** Return a reference to the SimbodyMatterSubsystem that is contained
in the MultibodySystem that is associated with this Assembler. **/
const SimbodyMatterSubsystem& getMatterSubsystem() const
{   return system.getMatterSubsystem(); }
/*@}*/

/** Destruct the Assembler objects and any Assembly Condition objects it
contains. **/
~Assembler();



//------------------------------------------------------------------------------
                           private: // methods
//------------------------------------------------------------------------------
// Note that the internalState is realized to Stage::Position on return.
void setInternalStateFromFreeQs(const Vector& freeQs) {
    assert(freeQs.size() == getNumFreeQs());
    Vector& q = internalState.updQ();
    for (FreeQIndex fx(0); fx < getNumFreeQs(); ++fx)
        q[getQIndexOfFreeQ(fx)] = freeQs[fx];
    system.realize(internalState, Stage::Position);
}

Vector getFreeQsFromInternalState() const {
    Vector freeQs(getNumFreeQs());
    const Vector& q = internalState.getQ();
    for (FreeQIndex fx(0); fx < getNumFreeQs(); ++fx)
        freeQs[fx] = q[getQIndexOfFreeQ(fx)];
    return freeQs;
}

void reinitializeWithExtraQsLocked
    (const Array_<QIndex>& toBeLocked) const;

// Implement the Study interface.
const System& getSystemVirtual() const override
{   return getSystem(); }
const State& getCurrentStateVirtual() const override
{   return getInternalState(); } 
const State& getInternalStateVirtual() const override
{   return getInternalState(); }  

State& updInternalStateVirtual() override {
    static State empty;
    SimTK_ERRCHK_ALWAYS(!"write to internal state", 
        "Study::updInternalState()", 
        "The Assembler Study does not provide write access to its "
        " internal state.");
    return empty;
}

Real getAccuracyInUseVirtual() const override
{   return getAccuracyInUse(); }
Real getConstraintToleranceInUseVirtual() const override
{   return getErrorToleranceInUse(); }

//------------------------------------------------------------------------------
                           private: // data members 
//------------------------------------------------------------------------------
const MultibodySystem&          system;
Array_<const EventReporter*>    reporters; // just references; don't delete

// These members affect the behavior of the assembly algorithm.
static const int OODefaultAccuracy = 1000; // 1/accuracy if acc==0
Real    accuracy;               // 0 means use 1/OODefaultAccuracy
Real    tolerance;              // 0 means use accuracy/10
bool    forceNumericalGradient; // ignore analytic gradient methods
bool    forceNumericalJacobian; // ignore analytic Jacobian methods
bool    useRMSErrorNorm;        // what norm defines success?

// Changes to any of these data members set isInitialized()=false.
State                           internalState;

// These are the mobilizers that were set in lockMobilizer(). They are
// separate from those involved in individually-locked q's.
LockedMobilizers                userLockedMobilizers;
// These are locks placed on individual q's; they are independent of the
// locked mobilizer settings.
LockedQs                        userLockedQs;
// These are range restrictions placed on individual q's.
RestrictedQs                    userRestrictedQs;

// These are (condition,weight) pairs with weight==Infinity meaning
// constraint; weight==0 meaning currently ignored; and any other
// positive weight meaning a goal.
Array_<AssemblyCondition*,AssemblyConditionIndex> 
                                        conditions;
Array_<Real,AssemblyConditionIndex>     weights;

// We always have an assembly condition for the Constraints which are
// enabled in the System; this is the index which can be used to 
// retrieve that condition. The default weight is Infinity.
AssemblyConditionIndex                  systemConstraints;


// These are filled in when the Assembler is initialized.
mutable bool                            alreadyInitialized;

// These are extra q's we removed for numerical reasons.
mutable Array_<QIndex>                  extraQsLocked;

// These represent restrictions on the independent variables (q's).
mutable std::set<QIndex>                lockedQs;
mutable Array_<FreeQIndex,QIndex>       q2FreeQ;    // nq of these
mutable Array_<QIndex,FreeQIndex>       freeQ2Q;    // nfreeQ of these
// 0 length if no bounds; otherwise, index by FreeQIndex.
mutable Vector                          lower, upper;

// These represent the active assembly conditions.
mutable Array_<AssemblyConditionIndex>  errors;
mutable Array_<int>                     nTermsPerError;
mutable Array_<AssemblyConditionIndex>  goals;

class AssemblerSystem; // local class
mutable AssemblerSystem* asmSys;
mutable Optimizer*       optimizer;

mutable int nAssemblySteps;   // count assemble() and track() calls
mutable int nInitializations; // # times we had to reinitialize

friend class AssemblerSystem;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLER_H_
