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
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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
value is a prescribed function of time will be set to that value, while free
q's are available for satisfying the assembly error conditions and goals. The 
assembly
error conditions are just the errors in the position (holonomic) constraints 
that are present in the MultibodySystem and currently enabled. (Quaternion 
normalization constraints will also be satisfied, but do not generate assembly
errors.) There are no default assembly goals. This is very similiar in behavior
to the System's project() method except that project() considers it an error if
the constraints aren't already close to being satisfied initially, while 
Assembler will attempt to satisfy them regardless, and may take a series of 
increasingly desperate measures to do so. 

<h2>Basic assembly:</h2>
This is the most common use of the Assembler: modify a System's given State
so that its configuration (set of generalized coordinates q) satisfies the 
position-affecting Motion and Constraint objects in the System that are 
currently enabled in that State. This 
is done to a default tolerance if you don't provide one, and that tolerance is
suitable for use with subsequent dynamic studies that are run at their default
tolerances. Although the assembly begins from the configuration provided in
the initial state, no assumption is made about how close this initial 
configuration is to one that satisfies the assembly conditions.
@code
  MultibodySystem system;
  // ... build system; get initial state

  Assembler assembler(system); // construct the Assembler study object
  try // modify state to satisfy Constraints
  {   assembler.assemble(state); }
  catch (std::exception exc)
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
  catch (std::exception exc)
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
in the System that are enabled must be satisifed to within the assembly
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



//------------------------------------------------------------------------------
//                            ASSEMBLY CONDITION
//------------------------------------------------------------------------------
/** Define an assembly condition consisting of a scalar goal and/or a 
related set of assembly error equations (that is, an objective and/or some 
constraints). Whether the goal or error is used depends on the weighting
assigned to this AssemblyCondition. A finite weight indicates that the
goal should be used (and combined with other goals); an infinite weighting 
means that each error must independently be satisfied to tolerance.  **/
class SimTK_SIMBODY_EXPORT AssemblyCondition {
public:

/** Base class constructor just takes the assembly condition name and 
saves it. **/
explicit AssemblyCondition(const String& name) 
:   name(name), assembler(0) {}

/** Destructor is virtual for use by derived classes. **/
virtual ~AssemblyCondition() {}

/** This is called whenever the Assembler is initialized in case this
assembly condition wants to do some internal work before getting started.
None of the other virtual methods will be called until this one has been,
except possibly the destructor. The set of free q's and the internal
State are valid at this point and can be retrieved from the Assembler
stored in the base class. **/
virtual int initializeCondition() const {return 0;}

/** This is called whenever the containing Assembler is uninitialized in
case this assembly condition has some cleanup to do. **/
virtual void uninitializeCondition() const {}

/** Calculate the amount by which this assembly condition is violated
by the q values in the given state, with one scalar error per assembly
equation returned in \a err. The functional return should be zero if
successful; negative values are reserved with -1 meaning "not implemented";
return a positive value if your implementation is unable to evaluate the 
error at the current state. If this method is not implemented then you must
implement calcGoal() and this assembly condition may only be used as a 
goal, not a requirement. **/
virtual int calcErrors(const State& state, Vector& err) const
{   return -1; }

/** Override to supply an analytic Jacobian for the assembly errors
returned by calcErrors(). The returned Jacobian must be nErr X nFreeQs; 
that is, if there is only one assembly error equation the returned matrix 
is a single row (that's the transpose of the gradient). The functional 
return should be zero if this succeeds; negative values are reserved with
the default implementation returning -1 which indicates that the
Jacobian must be calculated numerically using the calcErrors() method.
Return a positive value if your implementation is unable to evaluate the 
Jacobian at the current state. **/
virtual int calcErrorJacobian(const State& state, Matrix& jacobian) const
{   return -1; }

/** Override to supply an efficient method for determining how many errors
will be returned by calcErrors(). Otherwise the default implementation 
determines this by making a call to calcErrors() and returning the size
of the returned error vector. The functional return should be zero if this
succeeds; negative values are reserved; return a positive value if your
implementation of this method can't determine the number of errors with
the given state (unlikely!). **/
virtual int getNumErrors(const State& state) const 
{   Vector err;
    const int status = calcErrors(state, err);
    if (status == 0)
        return err.size();
    SimTK_ERRCHK1_ALWAYS(status != -1, "AssemblyCondition::getNumErrors()",
        "The default implementation of getNumErrors() depends on"
        " calcErrors() but that method was not implemented for assembly"
        " condition '%s'.", name.c_str());
    SimTK_ERRCHK2_ALWAYS(status == 0,  "AssemblyCondition::getNumErrors()",
        "The default implementation of getNumErrors() uses calcErrors()"
        " which returned status %d (assembly condition '%s').", 
        status, name.c_str());
    return -1; // NOTREACHED
}

/** Calculate the current contribution (>= 0) of this assembly condition to
the goal value that is being minimized. If this isn't overridden we'll 
generate it by combining the m errors returned by calcErrors() in a mean
sum of squares: goal = err^2/m. **/
virtual int calcGoal(const State& state, Real& goal) const
{   static Vector err;
    const int status = calcErrors(state, err);
    if (status == 0)
    {   goal = err.normSqr() / std::max(1,err.size());
        return 0; }
    SimTK_ERRCHK1_ALWAYS(status != -1, "AssemblyCondition::calcGoal()",
        "The default implementation of calcGoal() depends on calcErrors()"
        " but that method was not implemented for assembly condition '%s'.",
        name.c_str());
    SimTK_ERRCHK2_ALWAYS(status == 0,  "AssemblyCondition::calcGoal()",
        "The default implementation of calcGoal() uses calcErrors() which"
        " returned status %d (assembly condition '%s').", 
        status, name.c_str());
    return -1; // NOTREACHED
}

/** Override to supply an analytic gradient for this assembly condition's
goal. The returned gradient must be nFreeQ X 1; that is, it is a column
vector giving the partial derivative of the goal with respect to each of
the free q's in order. The functional return should be zero if this 
succeeds. The default implementation return -1 which indicates that the
gradient must be calculated numerically using the calcGoal() method. **/
virtual int calcGoalGradient(const State& state, Vector& gradient) const
{   return -1; }

/** Return the name assigned to this AssemblyCondition on construction. **/
const char* getName() const {return name.c_str();}

/** Test whether this AssemblyCondition has already been adopted by an 
Assembler. **/
bool isInAssembler() const {return assembler != 0;}
/** Return the Assembler that has adopted this AssemblyCondition. This will
throw an exception if there is no such Assembler; use isInAssembler() first
if you're not sure. **/
const Assembler& getAssembler() const 
{   assert(assembler); return *assembler;}
/** Return the AssemblyConditionIndex of this concrete AssemblyCondition
within the Assembler that has adopted it. This returned index will be
invalid if this AssemblyCondition has not yet been adopted. **/
AssemblyConditionIndex getAssemblyConditionIndex() const 
{   return myAssemblyConditionIndex; }

//------------------------------------------------------------------------------
                                 protected:
//------------------------------------------------------------------------------
// These are useful when writing concrete AssemblyConditions.

/** Ask the assembler how many free q's there are; only valid after
initialization but does not invoke initialization. **/
int getNumFreeQs() const {return getAssembler().getNumFreeQs();}
/** Ask the assembler where to find the actual q in the State that corresponds
to a given free q; only valid after initialization but does not invoke 
initialization. **/
QIndex getQIndexOfFreeQ(Assembler::FreeQIndex fx) const
{   return getAssembler().getQIndexOfFreeQ(fx); }
/** Ask the assembler where to find the free q (if any) that corresponds
to a given q in the State; only valid after initialization but does not invoke 
initialization. **/
Assembler::FreeQIndex getFreeQIndexOfQ(QIndex qx) const
{   return getAssembler().getFreeQIndexOfQ(qx); }
/** Ask the assembler for the MultibodySystem with which it is associated. **/ 
const MultibodySystem& getMultibodySystem() const
{   return getAssembler().getMultibodySystem(); }
/** Ask the assembler for the MultibodySystem with which it is associated
and extract the SimbodyMatterSubsystem contained therein. **/
const SimbodyMatterSubsystem& getMatterSubsystem() const
{   return getMultibodySystem().getMatterSubsystem(); }

/** Call this method before doing anything that logically requires the 
Assembler, or at least this AssemblyCondition, to have been initialized. **/
void initializeAssembler() const {
    // The Assembler will in turn invoke initializeCondition().
    if (isInAssembler()) getAssembler().initialize();
    else                 initializeCondition();
}

/** Call this when modifying any parameter of the concrete AssemblyCondition
that would require reinitialization of the Assembler or the 
AssemblyCondition. **/
void uninitializeAssembler() const {
    // The Assembler will in turn invoke uninitializeCondition().
    if (isInAssembler()) getAssembler().uninitialize();
    else                 uninitializeCondition();
}

//------------------------------------------------------------------------------
                                   private:
//------------------------------------------------------------------------------
// This method is used by the Assembler when the AssemblyCondition object 
// is adopted.
friend class Assembler;
void setAssembler(const Assembler& assembler, AssemblyConditionIndex acx) {
    assert(!this->assembler);
    this->assembler = &assembler;
    this->myAssemblyConditionIndex = acx;
}

String                  name; // assembly condition name
const Assembler*        assembler;
AssemblyConditionIndex  myAssemblyConditionIndex;
};



//------------------------------------------------------------------------------
//                                 Q VALUE
//------------------------------------------------------------------------------
/** This AssemblyCondition requests that a particular generalized coordinate
end up with a specified value. You can use this as a goal or a constraint,
depending on how serious you are about this requirement. **/
class QValue : public AssemblyCondition {
public:
    /** Construct an assembly condition that requests that the specified
    generalized coordinate be brought to the indicated value. The value 
    can be changed subsequently using setValue(). **/
    QValue(MobilizedBodyIndex mbx, MobilizerQIndex qx,
           Real value)
    :   AssemblyCondition("QValue"), 
        mobodIndex(mbx), qIndex(qx), value(value) {}

    /** Return the currently set value to be used for this generalized
    coordinate. **/
    Real getValue() const {return value;}
    /** Change the value to be used for this generalized coordinate; this
    can be done repeatedly during tracking to follow changing requirements. **/
    void setValue(Real newValue) {value=newValue;}

    // For constraint:
    int getNumEquations(const State&) const {return 1;}
    int calcErrors(const State& state, Vector& error) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        error.resize(1);
        error[0] = mobod.getOneQ(state, qIndex) - value;
        return 0;
    }
    // Error jacobian is a zero-row except for a 1 in this q's entry (if
    // this q is free).
    int calcErrorJacobian(const State& state, Matrix& J) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        J.resize(1, getNumFreeQs());
        J = 0; // will have at most one non-zero

        // Find the FreeQIndex corresponding to this q.
        const QIndex thisIx = QIndex(mobod.getFirstQIndex(state)+qIndex);
        const Assembler::FreeQIndex thisFreeIx = getFreeQIndexOfQ(thisIx);

        // If this q isn't free then there is no way to affect the error
        // so the Jacobian stays all-zero.
        if (thisFreeIx.isValid())
            J(0,thisFreeIx) = 1;

        return 0;
    }

    // For goal: goal = (q-value)^2 / 2 (the /2 is for gradient beauty)
    int calcGoal(const State& state, Real& goal) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        goal = square(mobod.getOneQ(state, qIndex) - value) / 2;
        return 0;
    }
    // Return a gradient with only this q's entry non-zero (if
    // this q is free).
    int calcGoalGradient(const State& state, Vector& grad) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        grad.resize(getNumFreeQs());
        grad = 0; // will have at most one non-zero

        // Find the FreeQIndex corresponding to this q.
        const QIndex thisIx = QIndex(mobod.getFirstQIndex(state)+qIndex);
        const Assembler::FreeQIndex thisFreeIx = getFreeQIndexOfQ(thisIx);

        // If this q isn't free then there is no way to affect the goal
        // so the gradient stays all-zero.
        if (thisFreeIx.isValid())
            grad[thisFreeIx] = mobod.getOneQ(state, qIndex) - value;

        return 0;
    }

private:
    MobilizedBodyIndex mobodIndex;
    MobilizerQIndex    qIndex;
    Real               value;
};



//------------------------------------------------------------------------------
//                                  MARKERS
//------------------------------------------------------------------------------
/** This AssemblyCondition specifies a correspondence between stations on
mobilized bodies ("markers") and fixed ground-frame locations ("observations").
The idea is to adjust the q's so that each marker is located close to its
corresponding observation. This is normally used as a goal since we don't
expect a perfect fit, but you can use these as a set of assembly error 
conditions if there are enough degrees of freedom to achieve a near-perfect 
solution. 

Markers are defined one at a time and assigned sequential marker index values
of type Markers::MarkerIx. They may optionally be given unique, case-sensitive
names, and we will keep a map from name to MarkerIx. A default name will be 
assigned if none is given. A weight is assigned to every marker, with default 
weight=1. We do not expect that all the markers will be used; markers with 
weights of zero will not be included in the study, nor will markers for which 
no observation is given.

Once specified, the marker definitions do not change during a series of inverse
kinematic (tracking) steps. The observations, on the other hand, are expected 
to come from a time series of experimental measurements of marker locations and
will be different at every step. They typically come from a file organized by 
"frame", meaning an observation time and a set of observed locations, one per 
marker, corresponding to that time. During initial setup, the number of 
observations per frame and their correspondence to the defined markers is 
specified. They can be in any order, may skip some markers, and may include 
data for markers that are not defined. However, once initialized each frame 
must supply the same information in the same order. Data for an unobserved 
marker can be provided as NaN in which case it will be ignored in that frame. 
The frame time is supplied to the track() method which initiates assembly for 
a frame.

Observation-marker correspondence maps a ObservationIx to a unique MarkerIx. 
By default, we'll expect to get an observation for each marker and that the 
observation order and the marker order are the same, i.e. 
ObservationIx==MarkerIx for every marker. However, you can instead define 
observation/marker correspondence yourself, (\e after all markers have been 
defined), via one of the defineObservationOrder() methods. This is done by 
supplying an array of MarkerIx values, or an array of Marker names, with the 
array elements ordered by ObservationIx. Any invalid marker index or 
unrecognized marker name means we will ignore values provide for that 
observation; similarly, any markers whose index or name is not specified at 
all will be ignored. 
**/
class SimTK_SIMBODY_EXPORT Markers : public AssemblyCondition {

// This is a private class used in the implementation below but not
// accessible through the API.
struct Marker {
    Marker(const String& name, MobilizedBodyIndex bodyB, 
           const Vec3& markerInB, Real weight = 1)
    :   name(name), bodyB(bodyB), markerInB(markerInB), weight(weight) 
    { assert(weight >= 0); }

    Marker(MobilizedBodyIndex bodyB, const Vec3& markerInB, Real weight=1)
    :   name(""), bodyB(bodyB), markerInB(markerInB), weight(weight) 
    { assert(weight >= 0); }

    String              name;
    MobilizedBodyIndex  bodyB;
    Vec3                markerInB;
    Real                weight; 
};

public:

/** Define the MarkerIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,MarkerIx);
/** Define the ObservationIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,ObservationIx);



//------------------------------------------------------------------------------
/** @name                Construction and setup
These methods are used as an extended construction phase for Markers
objects, defining the markers and observations that will be used in the
subsequent tracking steps. **/
/*@{*/

/** The default constructor creates an empty Markers AssemblyCondition
object that should be filled in with calls to addMarker() and optionally
defineObservationOrder(). **/
Markers() : AssemblyCondition("Markers") {}

/** Define a new marker attached to a particular MobilizedBody. Note that
a marker will be ignored unless an observation is provided for it.
@param[in]      name
    A unique name to be used to identify this marker. If the name is
    empty or blank, a default name will be supplied.
@param[in]      bodyB
    The MobilizedBody to which this marker is fixed. Markers on Ground
    are allowed but will be ignored.
@param[in]      markerInB
    This is the position vector of the marker in \a bodyB's local frame,
    also known as the marker's "station" on \a bodyB.
@param[in]      weight
    An optional weight for use in defining the objective function, which
    combines errors in this marker's position with errors in other markers'
    positions. If the weight is zero this marker is ignored.
@return The unique marker index number assigned to this marker. These are
assigned sequentially as the marker are added. 
@note Adding a marker invalidates any observation/marker correspondence; be
sure to call defineObservationOrder() \e after defining all your markers. **/
MarkerIx addMarker(const String& name, MobilizedBodyIndex bodyB, 
                   const Vec3& markerInB, Real weight=1)
{   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "Markers::addMarker()", "Illegal marker weight %g.", weight);
    uninitializeAssembler();
    // Forget any previously-established observation/marker correspondence.
    observation2marker.clear(); marker2observation.clear(); 
    observations.clear();
    const MarkerIx ix(markers.size());
    String nm = String::trimWhiteSpace(name);
    if (nm.empty())
        nm = String("_UNNAMED_") + String(ix);

    std::pair< std::map<String,MarkerIx>::iterator, bool >
        found = markersByName.insert(std::make_pair(nm,ix));
    SimTK_ERRCHK2_ALWAYS(found.second, // true if insertion was done
        "Markers::addMarker()",
        "Marker name '%s' was already use for Marker %d.",
        nm.c_str(), (int)found.first->second); 

    markers.push_back(Marker(nm,bodyB,markerInB,weight));
    return ix; 
}

/** Define an unnamed marker. A default name will be assigned; that name 
will be "_UNNAMED_XX" where XX is the MarkerIx assigned to that marker 
(don't use names of that form yourself).  
@see addMarker(name,...) for more information. **/
MarkerIx addMarker(MobilizedBodyIndex bodyB, const Vec3& markerInB,
                   Real weight=1)
{   return addMarker("", bodyB, markerInB, weight); }


/** Define the meaning of the observation data by giving the MarkerIx 
associated with each observation. The length of the array of marker indices 
defines the expected number of observations to be provided for each observation
frame. Any marker index that is supplied with an invalid value means that the
corresponding observation will be present in the supplied data but should be 
ignored.
@param[in]          observationOrder
    This is an array of marker index values, one per observation, that defines
    both the number of expected observations and the marker corresponding
    to each observation. Markers can be in any order; an invalid marker index
    means that observation will be provided but should be ignored; 
    markers whose indices are never listed are ignored. If \a observationOrder
    is supplied as a zero-length array, then we'll assume there are as
    many observations as markers and that their indices match.

@note If you don't call this method at all, a default correspondence will
be defined as described for a zero-length \a observationOrder array (that is,
same number of observations and markers with matching indices). Whenever you 
add a new marker, any previously defined observation order is forgotten so the 
default correspondence will be used unless you call this again. **/
void defineObservationOrder(const Array_<MarkerIx>& observationOrder) {
    uninitializeAssembler();
    if (observationOrder.empty()) {
        observation2marker.resize(markers.size());
        for (MarkerIx mx(0); mx < markers.size(); ++mx)
            observation2marker[ObservationIx(mx)] = mx;
    } else 
        observation2marker = observationOrder;
    marker2observation.clear(); 
    // We might need to grow this more, but this is an OK starting guess.
    marker2observation.resize(observation2marker.size()); // all invalid
    for (ObservationIx ox(0); ox < observation2marker.size(); ++ox) {
        const MarkerIx mx = observation2marker[ox];
        if (!mx.isValid()) continue;

        if (marker2observation.size() <= mx)
            marker2observation.resize(mx+1);
        SimTK_ERRCHK4_ALWAYS(!marker2observation[mx].isValid(),
            "Markers::defineObservationOrder()", 
            "An attempt was made to associate Marker %d (%s) with" 
            " Observations %d and %d; only one Observation per Marker"
            " is permitted.",
            (int)mx, getMarkerName(mx).c_str(), 
            (int)marker2observation[mx], (int)ox);

        marker2observation[mx] = ox;
    }
    // Make room for marker observations.
    observations.clear();
    observations.resize(observation2marker.size(),Vec3(NaN));
}

/** Define the meaning of the observations by giving the marker name 
corresponding to each observation, as a SimTK::Array_<String>. The length of 
the array of marker indices defines the expected number of observations. Any
marker name that is unrecognized or empty means that the corresponding 
observation will be present in the supplied data but should be ignored. **/
void defineObservationOrder(const Array_<String>& observationOrder) 
{   Array_<MarkerIx> markerIxs(observationOrder.size());
    for (ObservationIx ox(0); ox < observationOrder.size(); ++ox)
        markerIxs[ox] = getMarkerIx(observationOrder[ox]);
    defineObservationOrder(markerIxs); }

/** Define observation order using an std::vector of SimTK::String. */
// no copy required
void defineObservationOrder(const std::vector<String>& observationOrder)
{   defineObservationOrder(ArrayViewConst_<String>(observationOrder)); }


/** Define observation order using an Array_ of std::string. */
// must copy
void defineObservationOrder(const Array_<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using an std::vector of std::string. */
// must copy
void defineObservationOrder(const std::vector<std::string>& observationOrder) 
{   const Array_<String> observations(observationOrder); // copy
    defineObservationOrder(observations); }

/** Define observation order using a C array of const char* names. */
void defineObservationOrder(int n, const char* const observationOrder[]) 
{   Array_<MarkerIx> markerIxs(n);
    for (ObservationIx ox(0); ox < n; ++ox)
        markerIxs[ox] = getMarkerIx(String(observationOrder[ox]));
    defineObservationOrder(markerIxs); }
/*@}*/



//------------------------------------------------------------------------------
/** @name               Retrieve setup information
These methods are used to query information associated with the construction
and setup of this Markers object. This information does not normally change
during a marker-tracking study, although marker weights may be changed by 
some inverse kinematics methods. **/
/*@{*/

/** Return a count n of the number of currently-defined markers. Valid
marker index values (of type Markers::MarkerIx) are 0..n-1. **/
int getNumMarkers() const {return markers.size();}

/** Return the unique marker name assigned to the marker whose index
is provided. If the marker was defined without a name, this will return
the default name that was assigned to it. **/
const String& getMarkerName(MarkerIx ix) 
{   return markers[ix].name; }
/** Return the marker index associated with the given marker name. If the
name is not recognized the returned index will be invalid (test with
index.isValid()). **/
const MarkerIx getMarkerIx(const String& name) 
{   std::map<String,MarkerIx>::const_iterator p = markersByName.find(name);
    return p == markersByName.end() ? MarkerIx() : p->second; }

/** Get the weight currently in use for the specified marker; this can
be changed dynamically via changeMarkerWeight(). **/
Real getMarkerWeight(MarkerIx mx)
{   return markers[mx].weight; }

/** Get the MobilizedBodyIndex of the body associated with this marker. **/
MobilizedBodyIndex getMarkerBody(MarkerIx mx) const
{   return markers[mx].bodyB; }

/** Get the station (fixed location in its body frame) of the given marker. **/
const Vec3& getMarkerStation(MarkerIx mx) const
{   return markers[mx].markerInB; }

/** Return the number of observations that were defined via the last call to
defineObservationOrder(). These are not necessarily all being used. If 
defineObservationOrder() was never called, we'll expect the same number of
observations as markers although that won't be set up until the Assembler has
been initialized. **/
int getNumObservations() const {return observation2marker.size();}

/** Return the ObservationIx of the observation that is currently associated
with the given marker, or an invalid index if the marker doesn't have any
corresponding observation (in which case it is being ignored). An exception 
will be thrown if the given MarkerIx is not in the range 
0..getNumMarkers()-1. **/
ObservationIx getObservationIxForMarker(MarkerIx mx) const 
{ return marker2observation[mx]; }

/** Return true if the supplied marker is currently associated with an 
observation. @see getObservationIxForMarker() **/
bool hasObservation(MarkerIx mx) const 
{ return getObservationIxForMarker(mx).isValid(); }

/** Return the MarkerIx of the marker that is associated with the 
given observation, or an invalid index if the observation doesn't correspond
to any marker (in which case it is being ignored). An exception will be
thrown if the given ObservationIx is not in the range 
0..getNumObservations()-1. **/
MarkerIx getMarkerIxForObservation(ObservationIx ox) const 
{ return observation2marker[ox]; }

/** Return true if the supplied observation is currently associated with a 
marker. @see getMarkerIxForObservation() **/
bool hasMarker(ObservationIx ox) const 
{ return getMarkerIxForObservation(ox).isValid();}

/** The Markers assembly condition organizes the markers by body after
initialization; call this to get the list of markers on any particular body.
If necessary the Assembler will be initialized. It is an error if this 
assembly condition has not yet been adopted by an Assembler. **/
const Array_<MarkerIx>& getMarkersOnBody(MobilizedBodyIndex mbx) {
    static const Array_<MarkerIx> empty;
    SimTK_ERRCHK_ALWAYS(isInAssembler(), "Markers::getMarkersOnBody()",
        "This method can't be called until the Markers object has been"
        " adopted by an Assembler.");
    initializeAssembler();
    PerBodyMarkers::const_iterator bodyp = bodiesWithMarkers.find(mbx);
    return bodyp == bodiesWithMarkers.end() ? empty : bodyp->second;
}
/*@}*/



//------------------------------------------------------------------------------
/** @name                Execution methods
These methods can be called between tracking steps to make step-to-step
changes without reinitialization, and to access the current values of
step-to-step data including the resulting marker errors. **/
/*@{*/

/** Move a single marker's observed location without moving any of the others.
If the value contains a NaN, this marker/observation pair will be ignored the
next time the assembly goal cost function is calculated. **/
void moveOneObservation(ObservationIx ox, const Vec3& observation) 
{   SimTK_ERRCHK_ALWAYS(!observations.empty(), "Assembler::moveOneObservation()",
        "There are currently no observations defined. Either the Assembler"
        " needs to be initialized to get the default observation order, or you"
        " should call defineObservationOrder() explicitly.");
    SimTK_ERRCHK2_ALWAYS(ox.isValid() && ox < observations.size(),
        "Assembler::moveOneObservation()", "ObservationIx %d is invalid or"
        " out of range; there are %d observations currently defined. Use"
        " defineObservationOrder() to specify the set of observations and how"
        " they correspond to markers.", 
        (int)ox, (int)observations.size()); 
    observations[ox] = observation; 
}

/** Set the observed marker locations for a new observation frame. These are
the locations to which we will next attempt to move all the corresponding 
markers. Note that not all observations necessarily have corresponding markers
defined; locations of those markers must still be provided here but they will 
be ignored. The length of the \a allObservations array must be the same as the 
number of defined observations; you can obtain that using getNumObservations().
Any observations that contain a NaN will be ignored; that marker/observation 
pair will not be used in the next calculation of the assembly goal cost 
function. **/
void moveAllObservations(const Array_<Vec3>& observations) 
{   SimTK_ERRCHK2_ALWAYS((int)observations.size() == (int)observation2marker.size(),
        "Markers::moveAllObservations()",
        "Number of observations provided (%d) differs from the number of"
        " observations (%d) last defined with defineObservationOrder().",
        observations.size(), observation2marker.size());
    this->observations = observations; }

/** Change the weight associated with a particular marker. If this is just
a quantitative change (e.g., weight was 0.3 now it is 0.4) then this does
not require any reinitialization and will affect the goal calculation next
time it is done. If the weight changes to or from zero (a qualitative change)
then this will uninitialize the Assembler and all the internal data structures
will be changed to remove or add this marker from the list of active markers.
If you want to temporarily ignore a marker without reinitializing, you can
set its corresponding observation to NaN in which case it will simply be
skipped when the goal value is calculated. **/
void changeMarkerWeight(MarkerIx mx, Real weight) {
   SimTK_ERRCHK1_ALWAYS(isFinite(weight) && weight >= 0, 
        "Markers::changeMarkerWeight()", "Illegal marker weight %g.", weight);

    Marker& marker = markers[mx];
    if (marker.weight == weight)
        return;

    if (marker.weight == 0 || weight == 0)
        uninitializeAssembler(); // qualitative change

    marker.weight = weight;
}

/** Return the current value of the location for this observation. This
is where we will try to move the corresponding marker if there is one. 
The result might be NaN if there is no current value for this observation;
you can check using Vec3's isFinite() method. **/
const Vec3& getObservation(ObservationIx ox) const {return observations[ox];}
/** Return the current values of all the observed locations. This is where we 
will try to move the corresponding markers, for those observations that have 
corresponding markers defined. Some of the values may be NaN if there is
currently no corresponding observation. Note that these are indexed by
ObservationIx; use getObservationIxForMarker() to map a MarkerIx to its
corresponding ObservationIx. **/
const Array_<Vec3,ObservationIx>& getAllObservations() const
{   return observations; }

/** Using the current value of the internal state, calculate the ground
frame location of a particular marker. The difference between this location
and the corresponding observation is the current error for this marker. **/
Vec3 findCurrentMarkerLocation(MarkerIx mx) const;

/** Using the current value of the internal state, calculate the distance 
between the given marker's current location and its corresponding observed
location (unweighted). If the marker is not associated with an observation, 
or if the observed location is missing (indicated by a NaN value), then the 
error is reported as zero. 
@note If you actually want the square of the distance, you can save some
computation time by using findCurrentMarkerErrorSquared() which avoids the
square root needed to find the actual distance.
@see findCurrentMarkerErrorSquared() **/
Real findCurrentMarkerError(MarkerIx mx) const
{   return std::sqrt(findCurrentMarkerErrorSquared(mx)); }

/** Using the current value of the internal state, calculate the (unweighted)
square of the distance between the given marker's current location and its 
corresponding observed location (the squared distance is less expensive to 
compute than the distance). If the marker is not associated with an 
observation, or if the observed location is missing (indicated by a NaN 
value), then the error is reported as zero. 
@see findCurrentMarkerError() **/
Real findCurrentMarkerErrorSquared(MarkerIx mx) const {
    const ObservationIx ox = getObservationIxForMarker(mx);
    if (!ox.isValid()) return 0; // no observation for this marker
    const Vec3& loc = getObservation(ox);
    if (!loc.isFinite()) return 0; // NaN in observation; error is ignored
    return (findCurrentMarkerLocation(mx) - loc).normSqr();
}
/*@}*/



//------------------------------------------------------------------------------
/** @name              AssemblyCondition virtuals
These methods are the implementations of the AssemblyCondition virtuals. **/
/*@{*/
int calcErrors(const State& state, Vector& err) const;
int calcErrorJacobian(const State& state, Matrix& jacobian) const;
int getNumErrors(const State& state) const;
int calcGoal(const State& state, Real& goal) const;
int calcGoalGradient(const State& state, Vector& grad) const;
int initializeCondition() const;
void uninitializeCondition() const;
/*@}*/

//------------------------------------------------------------------------------
                                    private:
//------------------------------------------------------------------------------
const Marker& getMarker(MarkerIx i) const {return markers[i];}
Marker& updMarker(MarkerIx i) {uninitializeAssembler(); return markers[i];}

                                // data members                               
                               
// Marker definition. Any change here except a quantitative change to the
// marker's weight uninitializes the Assembler.
Array_<Marker,MarkerIx>         markers;
std::map<String,MarkerIx>       markersByName;

// Observation-marker corresondence specification. Any change here 
// uninitializes the Assembler.
Array_<MarkerIx,ObservationIx>  observation2marker;

// For convience in mapping from a marker to its corresponding observation.
// ObservationIx will be invalid if a particular marker has no associated
// observation.
Array_<ObservationIx,MarkerIx>  marker2observation;

// This is the current set of marker location observations, one per entry in 
// the observation2marker array. Changing the values here does not uninitialize
// the Assembler.            
Array_<Vec3,ObservationIx>      observations;

// After initialize, this groups the markers by body and weeds out
// any zero-weighted markers. TODO: skip low-weighted markers, at
// least at the start of the assembly.
typedef std::map<MobilizedBodyIndex,Array_<MarkerIx> > PerBodyMarkers;
mutable PerBodyMarkers          bodiesWithMarkers;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLER_H_
