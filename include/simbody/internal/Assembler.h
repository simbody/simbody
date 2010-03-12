#ifndef SimTK_SIMBODY_ASSEMBLER_H_
#define SimTK_SIMBODY_ASSEMBLER_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <set>
#include <map>

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

By default, all q's may be modified with no range restrictions. The assembly
error conditions are just the errors in the position (holonomic) constraints 
that are present in the MultibodySystem and currently enabled. (Quaternion 
normalization constraints will also be satisfied, but do not generate assembly
errors.) There are no default assembly goals. This is very similiar in behavior
to the System's project() method except that project() considers it an error if
the constraints aren't already close to being satisfied initially, while 
Assembler will attempt to satisfy them regardless, and make take a series of 
increasingly desperate measures to do so. 

<h2>Basic assembly:</h2>
This is the most common use of the Assembler: modify a System's given State
so that its configuration (set of generalized coordinates q) satisfies the 
System's built-in Constraints that are currently enabled in that State. This 
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
tolerance. You can selectively enable and disable Constraints in the
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
accuracy/10 if accuracy has been set, otherwise 1e-4; calling setTolerance() 
with no argument or with zero restores it to its default behavior. **/
void setErrorTolerance(Real tolerance=0) {
    SimTK_ERRCHK1_ALWAYS(0 <= tolerance,
        "Assembler::setTolerance()", "The requested error tolerance %g"
        " is illegal; we require 0 <= tolerance, with 0 indicating that"
        " the default tolerance (0.1*accuracy) is to be used.", tolerance);
    this->tolerance = tolerance;
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
void setAccuracy(Real accuracy=0) {
    SimTK_ERRCHK1_ALWAYS(0 <= accuracy && accuracy < 1,
        "Assembler::setAccuracy()", "The requested accuracy %g is illegal;"
        " we require 0 <= accuracy < 1, with 0 indicating that the default"
        " accuracy (10*tolerance) is to be used.", accuracy);
    this->accuracy = accuracy;
}
/** Obtain the accuracy setting that will be used during the next 
assemble() or track() call. The default is to use 10*tolerance if the 
error tolerance has been set, otherwise 1e-3. **/
Real getAccuracyInUse() const {
    return accuracy > 0 ? accuracy 
           : (tolerance > 0 ? 10*tolerance : Real(1)/OODefaultAccuracy); 
}


/** Change how the System's enabled built-in Constraints are weighted as
compared to other assembly conditions. If this is Infinity (the default) then
the built-ins are treated as must-satisfy constraints; otherwise they are 
included in the assembly cost function with the given weight. If the weight is 
given as zero the built-in Constraints will be ignored altogether.
@see setAssemblyConditionWeight() **/
void setSystemConstraintsWeight(Real weight)
{   assert(systemConstraints.isValid());
    setAssemblyConditionWeight(systemConstraints,weight); }

/** Return the current weight being given to the System's built-in
Constraints; the default is Infinity. 
@see getAssemblyConditionWeight() **/
Real getSystemConstraintsWeight()
{   assert(systemConstraints.isValid());
    return getAssemblyConditionWeight(systemConstraints); }

/** Set the weight to be used for this AssemblyCondition. If the weight is set
to 0, this condition will be disabled and will be ignored. If the weight is 
set to Infinity, the condition will be treated as an assembly error condition
that must be satisfied to tolerance. Otherwise (finite weight) the condition
will be treated as an assembly goal and the weight will be used to combine its 
cost function with that of the other assembly goals. **/
void setAssemblyConditionWeight(AssemblyConditionIndex condition, Real weight) {
    SimTK_INDEXCHECK_ALWAYS(condition, conditions.size(),
        "Assembler::setAssemblyConditionWeight()");
    SimTK_ERRCHK1_ALWAYS(weight >= 0, "Assembler::setAssemblyConditionWeight()",
        "Illegal weight %g; weight must be nonnegative.", weight);
    uninitialize();
    weights[condition] = weight;
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
of the heap-allocated AssemblyCondition object. We will use noramlly use the 
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
void setInternalState(const State& state) {
    uninitialize();
    getMatterSubsystem().convertToEulerAngles(state, internalState);
    system.realizeModel(internalState);
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
must lie. **/
/*@{*/

/** Lock this mobilizer at its starting position. This overrides any 
individual q specifications, so even if a q was specifically unlocked it
will not move until the mobilizer as a whole is unlocked. **/
void lockMobilizer(MobilizedBodyIndex mbx)
{   uninitialize(); userLockedMobilizers.insert(mbx); }
/** Unlock this mobilizer as a whole; some of its q's may remain locked
if they were locked individually. It is OK if this mobilizer was already
unlocked; in that case this does nothing. **/
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
lockMobilizer(); you have to unlockMobilizer() first. **/
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
direction. **/
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


/** Unrestrict a particular generalized coordinate q if it was
previously restricted. Note that this is independent of whether the q has
been locked with lockMobilizer() or lockQ(); that is, the q may still be
locked even though it is now unrestricted. **/
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
/** Reset all counters to zero; this also happens whenever the assembler
system is reinitialized either explicitly or due to system changes. **/
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
at their initial values. **/
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
static const int OODefaultAccuracy = 1000; // 1/accuracy if acc==tol==0
Real    accuracy;               // 0 means use 10*tolerance
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
/** Define an assembly condition consisting of a related set of assembly
error equations, with one scalar error returned per equation. When used as
an assembly constraint, each error will be satisfied to tolerance. When used
as an assembly goal, the norm of the errors will be weighted and combined with
other assembly goals. **/
class SimTK_SIMBODY_EXPORT AssemblyCondition {
public:

/** Base class constructor just takes the assembly condition name and 
saves it. **/
explicit AssemblyCondition(const String& name) 
:   name(name), assembler(0) {}

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
    QValue(MobilizedBodyIndex mbx, MobilizerQIndex qx,
           Real value)
    :   AssemblyCondition("QValue"), 
        mobodIndex(mbx), qIndex(qx), value(value) {}

    // For constraint:
    int getNumEquations(const State&) const {return 1;}
    int calcErrors(const State& state, Vector& error) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        error.resize(1);
        error[0] = mobod.getOneQ(state, qIndex) - value;
        return 0;
    }
    // int calcErrorJacobian(const State& state, Matrix& J) const;

    // For goal:
    int calcGoal(const State& state, Real& goal) const {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIndex);
        goal = square(mobod.getOneQ(state, qIndex) - value);
        return 0;
    }
    // int calcGoalGradient(const State& state, Vector& grad) const;

private:
    MobilizedBodyIndex mobodIndex;
    MobilizerQIndex    qIndex;
    Real               value;
};



//------------------------------------------------------------------------------
//                                  MARKERS
//------------------------------------------------------------------------------
/** This AssemblyCondition specifies a correspondence between stations on
mobilized bodies ("markers") and fixed ground-frame locations ("targets").
The idea is to adjust the q's so that each marker is located close to its
corresponding target. This is normally used as a goal since we don't
expect a perfect fit, but you can use these as a set of assembly error 
conditions if there are enough degrees of freedom to achieve a near-perfect 
solution. 

Markers are defined one at a time and assigned sequential marker index values
of type Markers::MarkerIx. They may optionally be given unique, case-sensitive
names, and we will keep a map from name to MarkerIx. A default name will be 
assigned if none is given. A weight is assigned to every marker, with default 
weight=1. We do not expect that all the markers will be used; markers with 
weights of zero will not be included in the study, nor will markers for which 
no target position is given.

Once specified, the marker definitions do not change during a series of inverse
kinematic (tracking) steps. The targets, on the other hand, are expected to 
come from a series of observations of the markers and will be different at 
every step. They typically come from a file organized by "frame", meaning an 
observation time and a set of locations, one per marker, corresponding to that
time. During initial setup, the number of targets per frame and their 
correspondence to the defined markers is specified. They can be in any order, 
may skip some markers, and may include data for markers that are not defined. 
However, once initialized each frame must supply the same information in the 
same order. Data for an unobserved target can be provided as NaN in which case
it will be ignored in that frame. The frame time is supplied to the track()
method which initiates assembly for a frame.

Target-marker correspondence maps a TargetIx to a unique MarkerIx. By default,
we'll expect to get a target for each marker and that the target order and
the marker order are the same, i.e. TargetIx==MarkerIx for every marker.
However, you can instead define target/marker correspondence yourself, 
(\e after all markers have been defined), via one of the defineTargetOrder() 
methods. This is done by supplying an array of MarkerIx values, or an array of
Marker names, with the array elements ordered by TargetIx. Any invalid marker 
index or unrecognized marker name means we will ignore values provide for that 
target; similarly, any markers whose index or name is not specified at all 
will be ignored. 
**/
class SimTK_SIMBODY_EXPORT Markers : public AssemblyCondition {
public:

Markers() : AssemblyCondition("Markers") {}

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

/** Define the MarkerIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,MarkerIx);
/** Define the TargetIx type which is just a uniquely-typed int. **/
SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,TargetIx);

/** Define a new marker attached to a particular MobilizedBody. Note that
a marker will be ignored unless a target position is provided for it.
@param[in]      name
    A unique name to be used to identify this marker. If the name is
    empty or blank, a default name will be supplied.
@param[in]      bodyB
    The MobilizedBody to which this marker is fixed. Markers on Ground
    are allowed but will be ignored.
@param[in]      markerInB
    This is the position vector of the marker in \a bodyB's local frame.
@param[in]      weight
    An optional weight for use in defining the objective function, which
    combines errors in this marker's position with errors in other markers'
    positions. If the weight is zero this marker is ignored.
@return The unique marker index number assigned to this marker. These are
assigned sequentially as the marker are added. 
@note Adding a marker invalidates any target/marker correspondence; be
sure to call defineTargetOrder() \e after defining all your markers. **/
MarkerIx addMarker(const String& name, MobilizedBodyIndex bodyB, 
                   const Vec3& markerInB, Real weight=1)
{   SimTK_ERRCHK1_ALWAYS(weight >= 0, "Markers::addMarker()",
        "Illegal marker weight %g.", weight);
    uninitializeAssembler();
    // Forget any previously-established target/marker correspondence.
    target2marker.clear(); marker2target.clear(); observations.clear();
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

/** Define the meaning of the target observation data by giving the 
MarkerIx associated with each target. The length of the array of marker 
indices defines the expected number of target observations to be provided
for each observation frame. Any marker index that is supplied with an
invalid value means that the corresponding target observation will be 
present in the supplied data but should be ignored.
@param[in]          targetOrder
    This is an array of marker index values, one per target, that defines
    both the number of expected target values and the marker corresponding
    to each target. Markers can be in any order; an invalid marker index
    means that target's value will be provided but should be ignored; 
    markers whose indices are never listed are ignored. If \a targetOrder
    is supplied as a zero-length array, then we'll assume there are as
    many targets as markers and that their indices match.
@note If you don't call this method at all, a default correspondence will
be defined as described for a zero-length \a targetOrder array (that is,
same number of targets and markers with matching indices). Whenever you add
a new marker, any previously defined target order is forgotten so the 
default correspondence will be used unless you call this again. **/
void defineTargetOrder(const Array_<MarkerIx>& targetOrder) {
    uninitializeAssembler();
    if (targetOrder.empty()) 
        for (MarkerIx mx(0); mx < markers.size(); ++mx)
            target2marker[TargetIx(mx)] = mx;
    else 
        target2marker = targetOrder;
    marker2target.clear(); 
    // We might need to grow this more, but this is a good starting guess.
    marker2target.resize(target2marker.size()); // all invalid
    for (TargetIx tx(0); tx < target2marker.size(); ++tx) {
        const MarkerIx mx = target2marker[tx];
        if (!mx.isValid()) continue;

        if (marker2target.size() <= mx)
            marker2target.resize(mx+1);
        SimTK_ERRCHK4_ALWAYS(!marker2target[mx].isValid(),
            "Markers::defineTargetOrder()", 
            "An attempt was made to associate Marker %d (%s) with" 
            " targets %d and %d; only one target is permitted.",
            (int)mx, getMarkerName(mx).c_str(), 
            (int)marker2target[mx], (int)tx);

        marker2target[mx] = tx;
    }
    // Make room for target observations.
    observations.clear();
    observations.resize(target2marker.size(),Vec3(NaN));
}

/** Define the target data by giving the marker name corresponding to
each target position, as a SimTK::Array_<String>. The length of the array 
of marker indices defines the expected number of target positions. Any 
marker name that is unrecognized or empty means that the corresponding 
target position will be present in the supplied data but should be ignored. **/
void defineTargetOrder(const Array_<String>& targetOrder) 
{   Array_<MarkerIx> markerIxs(targetOrder.size());
    for (TargetIx tx(0); tx < targetOrder.size(); ++tx)
        markerIxs[tx] = getMarkerIx(targetOrder[tx]);
    defineTargetOrder(markerIxs); }

/** Define target data using an std::vector of SimTK::String. */
void defineTargetOrder(const std::vector<String>& targetOrder)
{   defineTargetOrder(ArrayViewConst_<String>(targetOrder)); } // no copy

/** Define target data using an std::vector of std::string. */
void defineTargetOrder(const std::vector<std::string>& targetOrder) 
{   const Array_<String> targets(targetOrder); // copy
    defineTargetOrder(targets); }

/** Define target data using a C array of const char* names. */
void defineTargetOrder(int n, const char* const targetOrder[]) 
{   Array_<MarkerIx> markerIxs(n);
    for (TargetIx tx(0); tx < n; ++tx)
        markerIxs[tx] = getMarkerIx(String(targetOrder[tx]));
    defineTargetOrder(markerIxs); }

/** Move a single target's location without moving any of the others. If
the value contains a NaN, this marker/target pair will be ignored the
next time the assembly goal cost function is calculated. **/
void moveOneTarget(TargetIx tx, const Vec3& observation) 
{   observations[tx] = observation; }

/** Set the marker location for a new observation frame. These are the 
target locations to which we will next attempt to move all the 
corresponding markers. Note that not all targets necessarily have 
corresponding markers defined; locations of those targets must still be 
provided here but they will be ignored. The length of the \a allTargets 
array must be the same as the number of defined targets; you
can obtain that using getNumTargets(). Any observations that contain a
NaN will be ignored; that marker/target pair will not be used in the
next calculation of the assembly goal cost function. **/
void moveAllTargets(const Array_<Vec3>& observations) 
{   SimTK_ERRCHK2_ALWAYS(observations.size() == target2marker.size(),
        "Markers::setAllTargets()",
        "Number of observations provided (%d) differs from the number of"
        " targets (%d) last defined with defineTargetOrder().",
        observations.size(), target2marker.size());
    this->observations = observations; }

/** Return the number of targets that were defined via the last call to
defineTargetOrder(). These are not necessarily all being used. If 
defineTargetOrder() was never called, there will be the same number of
targets as markers although that won't be set up until the Assembler has
been initialized. **/
int getNumTargets() const {return target2marker.size();}
/** Return the current value of the target location for this target. This
is where we will try to move the corresponding marker if there is one. **/
const Vec3& getTargetObservation(TargetIx tx) {return observations[tx];}
/** Return the current values of all the target locations. This
is where we will try to move the corresponding markers, for those targets
that have corresponding markers defined. **/
const Array_<Vec3,TargetIx>& getAllTargetObservations() 
{   return observations; }

/** Return the TargetIx of the target that is associated with the 
given marker, or an invalid index if the marker doesn't have any
corresponding target (in which case it is being ignored). An exception will
be
thrown if the given MarkerIx is not in the range 0..getNumMarkers()-1. **/
TargetIx getTargetForMarker(MarkerIx mx) const {return marker2target[mx];}

/** Return true if the supplied marker is currently associated with a 
target. **/
bool hasTarget(MarkerIx mx) const {return getTargetForMarker(mx).isValid();}

/** Return the MarkerIx of the marker that is associated with the 
given target point, or an invalid index if the target doesn't correspond
to any marker (in which case it is being ignored). An exception will be
thrown if the given TargetIx is not in the range 0..getNumTargets()-1. **/
MarkerIx getMarkerForTarget(TargetIx tx) const {return target2marker[tx];}

/** Return true if the supplied target is currently associated with a 
marker. **/
bool hasMarker(TargetIx tx) const {return getMarkerForTarget(tx).isValid();}

const Marker& getMarker(MarkerIx i) const {return markers[i];}
Marker&       updMarker(MarkerIx i)       
{   uninitializeAssembler(); return markers[i]; }

const Array_<MarkerIx>& getMarkersOnBody(MobilizedBodyIndex mbx) {
    static const Array_<MarkerIx> empty;
    SimTK_ERRCHK_ALWAYS(isInAssembler(), "Markers::getMarkersOnBody()",
        "This method can't be called until the Markers have been"
        " adopted by an Assembler.");
    initializeAssembler();
    PerBodyMarkers::const_iterator bodyp = bodiesWithMarkers.find(mbx);
    return bodyp == bodiesWithMarkers.end() ? empty : bodyp->second;
}

/** TODO: should this be one scalar distance per marker/target pair,
or the three measure numbers of the vector separating them? **/
int calcErrors(const State& state, Vector& err) const;
int calcErrorJacobian(const State& state, Matrix& jacobian) const;
int getNumErrors(const State& state) const;
/** One half of the weighted sum of squared distances between 
corresponding points. **/
int calcGoal(const State& state, Real& goal) const;
int calcGoalGradient(const State& state, Vector& grad) const;
int initializeCondition() const;
void uninitializeCondition() const;

//------------------------------------------------------------------------------
                           private: // data members 
//------------------------------------------------------------------------------
// Marker definition. Any change here uninitializes the Assembler.
Array_<Marker,MarkerIx>         markers;
std::map<String,MarkerIx>       markersByName;

// Target-marker corresondence specification. Any change here 
// uninitializes the Assembler.
Array_<MarkerIx,TargetIx>       target2marker;

// For convience in mapping from a marker to its corresponding target.
// TargetIx will be invalid if a particular marker has no associated
// target.
Array_<TargetIx,MarkerIx>       marker2target;

// This is the current set of marker location observations (targets), one
// per entry in the target2marker array. Changing the values here does not
// uninitialize the Assembler.            
Array_<Vec3,TargetIx>           observations;

// After initialize, this groups the markers by body and weeds out
// any zero-weighted markers. TODO: skip low-weighted markers, at
// least at the start of the assembly.
typedef std::map<MobilizedBodyIndex,Array_<MarkerIx> > PerBodyMarkers;
mutable PerBodyMarkers          bodiesWithMarkers;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLER_H_
