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
of a Simbody MultibodySystem that satisfies the System's position constraints
plus optional additional assembly conditions, including limitations on which
joint coordinates may be changed and by how much. If successful, the final set
of q's will satisfy the constraints to within a specified tolerance. 

The complete specification for an Assembly study consists of four elements:
- The set of q's which may be modified by the study
- Limits on the allowable range of values for each q.
- A set of assembly constraints that \e must be met
- A set of weighted assembly goals that are to be achieved if possible

By default, all q's may be modified with no range restrictions. The assembly
constraints are just the position (holonomic) constraints that are present
in the MultibodySystem and currently enabled. Quaternion normalization 
constraints will also be satisfied if necessary, but they do not compete with
the other position constraints. (The Assembler's internal State is always
maintained with quaternions disabled.) There are no default assembly goals. 
This is very similiar in behavior to the System's project() method except that
project() considers it an error if the constraints aren't already close to 
being satisfied initially, while Assembler will attempt to satisfy them 
regardless, and make take a series of increasingly desperate measures to do so. 

Optional settings include:
- Locking particular mobilizers so that their q's can't be changed.
- Setting bounds on the acceptable range of values for some of the q's. 
- Defining additional assembly conditions.
- Assigning assembly conditions to be constraints or weighted goals.

Constraints are specified by giving an assembly condition a weight of 
Infinity. Anything with a lower weight is a goal and will be combined with all
the other goals into a single scalar objective. The built-in Constraints are
normally treated with infinite weight, but you can change them to goals instead
if you like; sometimes that can be useful as a step in getting a difficult-to-
assemble system assembled. 

<h2>Basic usage:</h2>
This is the most common use of the Assembler: modify a System's existing State
so that its configuration (set of generalized coordinates q) satisfies the 
System's built-in Constraints that are currently enabled in that State. This 
is done to a default tolerance if you don't provide one, and that tolerance is
suitable for use with subsequent dynamic studies that are run at their default
tolerances.
@code
    try {
    MultibodySystem system;
    // ... build system; get initial state

    Assembler assembler(system); // construct the Assembler study object
    try // modify state to satisfy Constraints
    {   assembler.assemble(state); }
    catch (std::exception exc)
    {   std::cout << "Assembly failed: " << exc.what() << std::endl; }
@endcode
**/
class SimTK_SIMBODY_EXPORT Assembler : public Study {
    typedef std::set<MobilizedBodyIndex>            LockedMobilizers;
    typedef std::set<MobilizerQIndex>               QSet;
    typedef std::map<MobilizedBodyIndex, QSet>      LockedQs;
    typedef std::map<MobilizerQIndex, Vec2>         QRanges;
    typedef std::map<MobilizedBodyIndex, QRanges>   RestrictedQs;
public:
    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Assembler,FreeQIndex);

    /** Create an Assembler study for the given MultibodySystem. The
    Assembler's current state is set to the System's default state but with
    Euler angles used instead of quaternions. **/
    explicit Assembler(const MultibodySystem& system);

    /** Destruct the Assembler objects and any Assembly Condition objects it
    contains. **/
    ~Assembler();

    /** Given a reference to an EventReporter, use this Reporter to provide 
    progress reporting. The EventReporter object must be owned by someone
    else and persist throughout the lifetime of this Assembler object. **/
    void addReporter(const EventReporter& reporter) {
        reporters.push_back(&reporter);
    }


    /** Return the goal value attained by the internal State's current settings
    for the free q's; this is a weighted sum of the individual goal values for 
    each assembly goal. Goal values are nonnegative scalars. **/
    Real calcCurrentGoal() const;
    /** This is the weighted norm of the assembly constraint errors directly
    comparable with the assembly tolerance setting. That is, if this number
    is less than or equal to tolerance, then the current state is a feasible
    assembly solution (although it may not be optimal). **/
    Real calcCurrentError() const;

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

    /** Initialize the Assembler to prepare for performing assembly analysis.
    This is normally called when needed but you can call it explicitly and then
    access methods that report on the properties of the system on which the
    analysis will be performed. The internal state should already have been
    set; if you want to provide the state use initialize(State). **/
    void initialize() const;
    /** Set the internal State and initialize. See setInternalState() and 
    initialize() methods for more information. **/
    void initialize(const State& state)
    {   setInternalState(state); initialize(); }
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
    while updateFromInternalState() will make sure you get the rotations in the form
    you want. **/
    const State& getInternalState() const {return internalState;}

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

    /** Starting with the current value of the internally-maintained State, 
    modify the q's in it to satisfy all the assembly conditions to within a 
    tolerance. The actual tolerance achieved is returned as the function value. 
    @param[in]          tol
        Optional tolerance; the default is 1e-4.
    @return The tolerance actually achieved.  **/
    Real assemble(Real tol = 1e-4);

    /** Continue a series of assembly steps that is already in progress,
    without restarting or reanalyzing the system. This is designed for
    use with a series of assembly frames that are close together so that
    no heroic measures are needed to go from one to the next. For the 
    first frame, and any time there might be a change to the problem 
    structure or a major change to the state, use assemble() instead of
    track(). **/
    Real track(Real tol = 1e-4);

    /** Given an initial value for the State, modify the q's in it to satisfy
    all the assembly conditions to within a tolerance. The actual tolerance 
    achieved is returned as the function value. 
    @param[in,out]      state    
        The initial and final State value. Only q's are modified.
    @param[in]          tol
        Optional tolerance; the default is 1e-4.
    @return The tolerance actually achieved.  **/
    Real assemble(State& state, Real tol = 1e-4) {
        setInternalState(state);
        Real achievedTol = assemble(tol);
        updateFromInternalState(state);
        return achievedTol;
    }

    /** @name              Parameter restrictions
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

    /** @name                Assembly conditions
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

    /** Change how the System's enabled built-in Constraints are weighted as
    compared to other assembly conditions. If this is Infinity the built-ins
    are treated as must-satisfy constraints; otherwise they are included in the
    objective function with the given weight. If the weight is given as zero
    the built-in Constraints are ignored altogether. **/
    void setSystemConstraintsWeight(Real weight)
    {   uninitialize(); assert(systemConstraints.isValid());
        weights[systemConstraints] = weight; }
    /** Return the current weight being given to the System's built-in
    Constraints; the default is Infinity. **/
    Real getSystemConstraintsWeight()
    {   assert(systemConstraints.isValid());
        return weights[systemConstraints]; }

    /** Add an assembly requirement to this Assembler study, taking over 
    ownership of the heap-allocated AssemblyCondition object. We will use the
    calcErrors() method of this object to determine errors which \e must be
    driven below tolerance for an assembly to be considered successful. **/
    AssemblyConditionIndex 
        adoptAssemblyConstraint(AssemblyCondition* p);
    /** Add an assembly goal to this Assembler study, taking over ownership
    of the heap-allocated AssemblyCondition object. An optional weight is used
    to combine this goal with others to form the assembly cost function. The
    default weight is 1; if the weight is 0 the goal is ignored and not
    evaluated at all; if the weight is Infinity this is actually an assembly
    constraint. **/
    AssemblyConditionIndex 
        adoptAssemblyGoal(AssemblyCondition* p, Real weight=1);

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
    /** Return the number of constraint error evaluations. **/
    int getNumConstraintEvals() const;
    /** Return the number of goal gradient evaluations. **/
    int getNumGoalGradientEvals()   const;
    /** Return the number of constraint error Jacobian evaluations. **/
    int getNumConstraintJacobianEvals()   const;
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

    /** This is useful for debugging but should not be used otherwise
    since the analytic gradient is to be preferred. **/
    void setForceNumericalGradient(bool yesno)
    {   forceNumericalGradient = yesno; }
    /** This is useful for debugging but should not be used otherwise
    since the analytic Jacobian is to be preferred. **/
    void setForceNumericalJacobian(bool yesno)
    {   forceNumericalJacobian = yesno; }
private:
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
private:
    const MultibodySystem&          system;
    Array_<const EventReporter*>    reporters; // just references; don't delete

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

    bool                                    forceNumericalGradient;
    bool                                    forceNumericalJacobian;

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
    mutable Array_<AssemblyConditionIndex>  constraints;
    mutable Array_<int>                     nErrorsPerConstraint;
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

    void setAssembler(const Assembler& assembler) {
        assert(!this->assembler);
        this->assembler = &assembler;
    }

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

    const char* getName() const {return name.c_str();}
    const Assembler& getAssembler() const 
    {   assert(assembler); return *assembler;}
    bool hasAssembler() const {return assembler != 0;}

    int getNumFreeQs() const {return getAssembler().getNumFreeQs();}
    QIndex getQIndexOfFreeQ(Assembler::FreeQIndex fx) const
    {   return getAssembler().getQIndexOfFreeQ(fx); }
    Assembler::FreeQIndex getFreeQIndexOfQ(QIndex qx) const
    {   return getAssembler().getFreeQIndexOfQ(qx); }
    const SimbodyMatterSubsystem& getMatterSubsystem() const
    {   return getAssembler().getMatterSubsystem(); }

    void initializeAssembler() const {
        // The Assembler will invoke initializeCondition().
        if (hasAssembler()) getAssembler().initialize();
        else                initializeCondition();
    }
    void uninitializeAssembler() const {
        // The Assembler will invoke uninitializeCondition().
        if (hasAssembler()) getAssembler().uninitialize();
        else                uninitializeCondition();
    }
private:
    String           name; // assembly condition name
    const Assembler* assembler;
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
The idea is to adjust the q's so that each markers is located close to its
corresponding target. This is normally used as a goal since we don't
expect a perfect fit, but you can use it as a constraint if there are enough
degrees of freedom to achieve a near-perfect solution. 

Markers are defined one at a time and assigned sequential marker index values
of type Markers::MarkerIx. They may optionally be given unique, case-sensitive
names, and we will keep a map from name to MarkerIx. A default name will be 
assigned if none is given; that name will be "_UNNAMED_XX" where XX is the 
MarkerIx assigned to that marker (don't use names of that form yourself). A 
weight is assigned to every marker, with default 1. We do not expect that all 
the markers will be used; markers with weights of zero will not be included in 
the study, nor will markers for which no target position is given.

Once specified, the marker definitions do not change during a series of inverse
kinematic steps. The targets, on the other hand, are expected to come from a
series of observations of the markers and will be different at every step.
They typically come from a file organized by "frame", meaning an observation
time and a set of locations, one per marker, corresponding to that time. During
initial setup, the number of targets per frame and their correspondence to the
defined markers is specified. They can be in any order, may skip some markers,
and may include data for markers that are not defined. However, once
initialized each frame must supply the same information in the same order.

Target-marker correspondence maps a TargetIx to a unique MarkerIx. This
can be specified directly by giving an array of MarkerIx values, or indirectly
by giving an array of unique Marker names, with the array elements ordered
by TargetIx. Any invalid marker index or unrecognized marker name means we
will ignore values provide for that target; similarly, any markers whose index
or name is not specfied at all will be ignored.
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

    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(Markers,MarkerIx);
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
    assigned sequentially as the marker are added. **/
    MarkerIx addMarker(const String& name, MobilizedBodyIndex bodyB, 
                       const Vec3& markerInB, Real weight=1)
    {   SimTK_ERRCHK1_ALWAYS(weight >= 0, "Markers::addMarker()",
            "Illegal marker weight %g.", weight);
        uninitializeAssembler();
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

        return ix; }

    /** Define an unnamed marker. 
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
    present in the supplied data but should be ignored. **/
    void defineTargetOrder(const Array_<MarkerIx>& targetOrder) {
        uninitializeAssembler();
        target2marker = targetOrder;
        marker2target.clear(); 
        // We might need to grow this more, but this is a good starting
        // guess.
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
        observations.resize(target2marker.size(),Vec3(NaN));
    }

    /** Define the target data by giving the marker name corresponding to
    each target position, as a SimTK::Array_<String>. The length of the array 
    of marker indices defines the expected number of target positions. Any 
    marker name that is unrecognized or empty means that the corresponding 
    target position will be present in the supplied data but should be ignored. **/
    void defineTargetOrder(const Array_<String>& targetOrder) {
        Array_<MarkerIx> markerIxs(targetOrder.size());
        for (TargetIx tx(0); tx < targetOrder.size(); ++tx)
            markerIxs[tx] = getMarkerIx(targetOrder[tx]);
        defineTargetOrder(markerIxs);
    }
    /** Define target data using an std::vector of SimTK::String. */
    void defineTargetOrder(const std::vector<String>& targetOrder)
    {   defineTargetOrder(ArrayViewConst_<String>(targetOrder)); } // no copy

    /** Define target data using an std::vector of std::string. */
    void defineTargetOrder(const std::vector<std::string>& targetOrder) {
        const Array_<String> targets(targetOrder); // copy
        defineTargetOrder(targets);
    }

    /** Define target data using a C array of const char* names. */
    void defineTargetOrder(int n, const char* const targetOrder[]) {
        Array_<MarkerIx> markerIxs(n);
        for (TargetIx tx(0); tx < n; ++tx)
            markerIxs[tx] = getMarkerIx(String(targetOrder[tx]));
        defineTargetOrder(markerIxs);
    }

    /** Set the time associated with the current set of marker location
    observations. We will perform the assembly at this time which may have an
    effect on the result if there are position constraints or prescribed
    motions that depend on time. **/
    void setFrameTime(Real time) {frameTime=time;}
    /** Return the currently-set frame time. This will be the value from 
    the initial State if no other frame time was supplied. **/
    Real getFrameTime() const {return frameTime;}

    /** Move a single target's location without moving any of the others and
    without changing the frame time. **/
    void moveOneTarget(TargetIx tx, const Vec3& observation) 
    {   observations[tx] = observation; }

    /** Set a new frame time and all the corresponding marker location 
    observations. These are the target locations to which we will next attempt
    to move all the corresponding markers. Note that not all targets 
    necessarily have corresponding markers defined; locations of those targets
    must still be provided here but they will be ignored. The length of the 
    \a allTargets array must be the same as the number of defined targets; you
    can obtain that using getNumTargets(). **/
    void moveAllTargets(Real time, const Array_<Vec3>& observations) {
        SimTK_ERRCHK2_ALWAYS(observations.size() == target2marker.size(),
            "Markers::setAllTargets()",
            "Number of observations provided (%d) differs from the number of"
            " targets (%d) last defined with defineTargetOrder().",
            observations.size(), target2marker.size());
        this->observations = observations;
    }

    /** Return the number of targets that were defined via the last call to
    defineTargetOrder(). These are not necessarily all being used. **/
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
    corresponding target (in which case it is being ignored). An exception will be
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
        SimTK_ERRCHK_ALWAYS(hasAssembler(), "Markers::getMarkersOnBody()",
            "This method can't be called until the Markers have been"
            " associated with an Assembler.");
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

private:
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
    // per entry in the target2marker array, and the time at which those
    // observations were made. Changing the values here does not
    // uninitialize the Assembler.
    Real                            frameTime;               
    Array_<Vec3,TargetIx>           observations;

    // After initialize, this groups the markers by body and weeds out
    // any zero-weighted markers. TODO: skip low-weighted markers, at
    // least at the start of the assembly.
    typedef std::map<MobilizedBodyIndex,Array_<MarkerIx> > PerBodyMarkers;
    mutable PerBodyMarkers          bodiesWithMarkers;

};

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLER_H_
