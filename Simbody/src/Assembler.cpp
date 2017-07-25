/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/AssemblyCondition.h"
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//------------------------------------------------------------------------------
//                               EXCEPTIONS
//------------------------------------------------------------------------------
class Assembler::AssembleFailed : public Exception::Base {
public:
    AssembleFailed
       (const char* fn, int ln, const char* why, 
        Real tolAchieved, Real tolRequired) 
       : Base(fn,ln)
    {
        setMessage(
            "Method Assembler::assemble() failed because:\n" + String(why)
            + "\nAssembly error tolerance achieved: "
            + String(tolAchieved) + " required: " + String(tolRequired)
            + ".");
    }
};
class Assembler::TrackFailed : public Exception::Base {
public:
    TrackFailed
       (const char* fn, int ln, const char* why, 
        Real tolAchieved, Real tolRequired) 
       : Base(fn,ln)
    {
        setMessage(
            "Method Assembler::track() failed because:\n" + String(why)
            + "\nAssembly error tolerance achieved: "
            + String(tolAchieved) + " required: " + String(tolRequired)
            + ".");
    }
};

//------------------------------------------------------------------------------
//                            BUILT IN CONSTRAINTS
//------------------------------------------------------------------------------
// This is the Assembly Condition representing the System's built-in
// Constraints. Only Constraints that are currently enabled are included.
// This class provides an efficient implementation for treating these
// Constraints either as an assembly requirement or an assembly goal.
namespace { // this class is local to this file
class BuiltInConstraints : public AssemblyCondition {
public:
    BuiltInConstraints() 
    :   AssemblyCondition("System Constraints") {}

    // Note that we have turned off quaternions so the number of q error
    // slots in the State includes only real holonomic constraint equations.
    int getNumErrors(const State& s) const override {return s.getNQErr();}

    // Return the system holonomic constraint errors as found in the State.
    int calcErrors(const State& state, Vector& err) const override {
        err = state.getQErr();
        return 0;
    }

    // The Jacobian of the holonomic constraint errors is Pq=PN^-1. We can get
    // that analytically from Simbody but we might have to strip out some
    // of the columns if we aren't using all the q's.
    int calcErrorJacobian(const State& state, Matrix& jacobian) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const int np = getNumFreeQs();
        const int nq = state.getNQ();

        jacobian.resize(state.getNQErr(), np); // ok cuz no quaternions

        if (np == nq) {
            // Jacobian is already the right shape.
            matter.calcPq(state, jacobian);
        } else {
            Matrix fullJac(state.getNQErr(), nq);
            matter.calcPq(state, fullJac);
            // Extract just the columns corresponding to free Qs
            for (Assembler::FreeQIndex fx(0); fx < np; ++fx)
                jacobian(fx) = fullJac(getQIndexOfFreeQ(fx));
        }
        return 0;
    }

    // Goal is qerr^2/2 (the /2 makes the gradient easier).
    int calcGoal(const State& state, Real& goal) const override {
        goal = (~state.getQErr() * state.getQErr()) / 2;
        return 0;
    }

    // Gradient is ~(d goal/dq) = ~(~qerr * dqerr/dq) = ~(~qerr*Pq)
    // = ~Pq qerr. This can be done in O(n+m) time since we can calculate
    // the matrix-vector product ~Pq*v in O(n+m) time, where
    // n=#q's and m=# constraint equations.
    int calcGoalGradient(const State& state, Vector& grad) const override {
        const SimbodyMatterSubsystem& matter = getMatterSubsystem();
        const int np = getNumFreeQs();
        const int nq = state.getNQ();

        grad.resize(np);

        if (np == nq) {
            // Nothing locked; analytic gradient is the right size
            matter.multiplyByPqTranspose(state, state.getQErr(), grad);
        } else {
            Vector fullGrad(nq);
            matter.multiplyByPqTranspose(state, state.getQErr(), fullGrad);
            // Extract just the entries corresponding to free Qs
            for (Assembler::FreeQIndex fx(0); fx < np; ++fx)
                grad[fx] = fullGrad[getQIndexOfFreeQ(fx)];
        }

        //cout << "built in grad=" << grad << endl;
        return 0;
    }

private:
};
} // end anonymous namespace


//------------------------------------------------------------------------------
//                            ASSEMBLER SYSTEM
//------------------------------------------------------------------------------

// This class defines the objective function which is passed to the Optimizer.
class Assembler::AssemblerSystem : public OptimizerSystem {
public:
    AssemblerSystem(Assembler& assembler)
    :   OptimizerSystem(assembler.getNumFreeQs()), assembler(assembler)
    {   getSystem().realize(getInternalState(), Stage::Time); 
        resetStats(); }

    // Convenient interface to objective function.
    Real calcCurrentGoal() const {
        Real val;
        const int status = objectiveFunc(getFreeQsFromInternalState(),true,val);
        SimTK_ERRCHK1_ALWAYS(status==0, 
            "AssemblerSystem::calcCurrentGoal()",
            "objectiveFunc() returned status %d.", status);
        return val;
    }

    // Convenient interface to gradient of objective function.
    Vector calcCurrentGradient() const {
        Vector grad(getNumFreeQs());
        const int status = 
            gradientFunc(getFreeQsFromInternalState(),true,grad);
        SimTK_ERRCHK1_ALWAYS(status==0, 
            "AssemblerSystem::calcCurrentGradient()",
            "gradientFunc() returned status %d.", status);
        return grad;
    }

    // Convenient interface to assembly constraint error function.
    Vector calcCurrentErrors() const {
        Vector errs(getNumEqualityConstraints());
        const int status = 
            constraintFunc(getFreeQsFromInternalState(), true, errs);
        SimTK_ERRCHK1_ALWAYS(status==0, 
            "AssemblerSystem::calcCurrentErrors()",
            "constraintFunc() returned status %d.", status);
        return errs;
    }

    // Convenient interface to Jacobian of assembly constraint error function.
    Matrix calcCurrentJacobian() const {
        Matrix jac(getNumEqualityConstraints(), getNumFreeQs());
        const int status = 
            constraintJacobian(getFreeQsFromInternalState(), true, jac);
        SimTK_ERRCHK1_ALWAYS(status==0, 
            "AssemblerSystem::calcCurrentJacobian()",
            "constraintJacobian() returned status %d.", status);
        return jac;
    }

    // Return the value of the objective to be minimized when the freeQs
    // have the values given by the parameters.
    int objectiveFunc(const Vector&     parameters, 
                      bool              new_parameters, 
                      Real&             objectiveValue) const override
    {   ++nEvalObjective;

        if (new_parameters)
            setInternalStateFromFreeQs(parameters);

        objectiveValue = 0;
        for (unsigned i=0; i < assembler.goals.size(); ++i) {
            AssemblyConditionIndex   goalIx = assembler.goals[i];
            const AssemblyCondition& cond   = *assembler.conditions[goalIx];
            Real goalValue;
            const int stat = cond.calcGoal(getInternalState(), goalValue);
            if (stat != 0)
                return stat;
            SimTK_ERRCHK2_ALWAYS(goalValue >= 0,
                "AssemblerSystem::objectiveFunc()",
                "calcGoal() method of assembly condition %s returned a"
                " negative or non-finite value %g.", cond.getName(), goalValue);
            objectiveValue += assembler.weights[goalIx] * goalValue;
        }

        //static int count = 0;
        //cout << "    " << count++ << " obj=" << objectiveValue << endl;
        return 0;
    }

    class NumGradientFunc : public Differentiator::GradientFunction {
    public:
        NumGradientFunc(Assembler& assembler,
                        const Array_<AssemblyConditionIndex>& numGoals) 
        :   Differentiator::GradientFunction(assembler.getNumFreeQs()),
            assembler(assembler), numGoals(numGoals) {}

        // This is the function that gets differentiated. We want it to
        // return fy = sum( w[i] * goal[i] ) for each of the goals that needs
        // a numerical gradient. Then we can calculate all of them at once.
        int f(const Vector& y, Real& fy) const override {
            assembler.setInternalStateFromFreeQs(y);
            fy = 0;
            for (unsigned i=0; i < numGoals.size(); ++i) {
                AssemblyConditionIndex goalIx = numGoals[i];
                const AssemblyCondition& cond = 
                    *assembler.conditions[goalIx];
                Real goalValue;
                const int stat = cond.calcGoal(assembler.getInternalState(), 
                                               goalValue);
                if (stat != 0)
                    return stat;
                fy += assembler.weights[goalIx] * goalValue;
            }
            return 0;
        }
    private:
        Assembler&                              assembler;
        const Array_<AssemblyConditionIndex>&   numGoals;
    };

    int gradientFunc(const Vector&     parameters, 
                     bool              new_parameters, 
                     Vector&           gradient) const override 
    {   SimTK_ASSERT2_ALWAYS(gradient.size() == getNumFreeQs(),
            "AssemblySystem::gradientFunc(): expected gradient vector of"
            " size %d but got %d; this is likely a problem with the optimizer"
            " which is required to allocate the right amount of space.",
            getNumFreeQs(), gradient.size());

        ++nEvalGradient;

        if (new_parameters)
            setInternalStateFromFreeQs(parameters);

        for (unsigned i=0; i < assembler.reporters.size(); ++i)
            assembler.reporters[i]->handleEvent(getInternalState());

        // This will record the indices of any goals we encounter that can't
        // provide their own gradients; we'll handle them all together at
        // the end.
        Array_<AssemblyConditionIndex> needNumericalGradient;

        gradient = 0;
        Vector tmpGradient(gradient.size());
        for (unsigned i=0; i < assembler.goals.size(); ++i) {
            AssemblyConditionIndex   goalIx = assembler.goals[i];
            const AssemblyCondition& cond   = *assembler.conditions[goalIx];
            const int stat = (assembler.forceNumericalGradient 
                              ? -1
                              : cond.calcGoalGradient(getInternalState(), 
                                                      tmpGradient));
            if (stat == -1) {
                needNumericalGradient.push_back(goalIx);
                continue;
            }
            if (stat != 0)
                return stat;

            gradient += assembler.weights[goalIx] * tmpGradient;
        }

        if (!needNumericalGradient.empty()) {
            //cout << "Need numerical gradient for " 
            //     << needNumericalGradient.size() << " goals." << endl;
            NumGradientFunc numGoals(assembler, needNumericalGradient);
            // Essential to use central difference here so that the
            // approximate gradient is actually zero at the optimum
            // solution, otherwise IpOpt won't converge.
            Differentiator gradNumGoals
               (numGoals,Differentiator::CentralDifference);
            // weights are already included here
            gradient += gradNumGoals.calcGradient(getFreeQsFromInternalState());

            nEvalObjective += gradNumGoals.getNumCallsToUserFunction();
        }

        //cout << "Grad=" << gradient << endl;

        return 0;
    }


    // Return the errors in the hard assembly error conditions.
    int constraintFunc(const Vector&    parameters, 
                       bool             new_parameters, 
                       Vector&          qerrs) const override
    {   ++nEvalConstraints;

        if (new_parameters)
            setInternalStateFromFreeQs(parameters);

        int nxtEqn = 0;
        for (unsigned i=0; i < assembler.errors.size(); ++i) {
            AssemblyConditionIndex   consIx = assembler.errors[i];
            const AssemblyCondition& cond   = *assembler.conditions[consIx];
            const int m = cond.getNumErrors(getInternalState());
            int stat = cond.calcErrors(getInternalState(), qerrs(nxtEqn,m));
            if (stat != 0)
                return stat;
            nxtEqn += m;
        }

        //cout << "    err=" << qerrs << endl;

        return 0;
    }

    class NumJacobianFunc : public Differentiator::JacobianFunction {
    public:
        NumJacobianFunc(Assembler& assembler,
                        const Array_<AssemblyConditionIndex>& numCons,
                        const Array_<int>& nErrorEqns,
                        int totalNEqns) 
        :   Differentiator::JacobianFunction
                (totalNEqns, assembler.getNumFreeQs()),
            assembler(assembler), numCons(numCons), nEqns(nErrorEqns), 
            totalNEqns(totalNEqns) 
        {   assert(numCons.size() == nEqns.size()); }

        // This is the function that gets differentiated. We want it to
        // return fy = [ err[i] ] for each of the assembly constraint 
        // conditions that needs a numerical gradient. Then we can calculate
        // all their Jacobians at once.
        int f(const Vector& y, Vector& fy) const override {
            assert(y.size() == assembler.getNumFreeQs());
            assert(fy.size() == totalNEqns);

            assembler.setInternalStateFromFreeQs(y);
            int nxtSlot = 0;
            for (unsigned i=0; i < numCons.size(); ++i) {
                AssemblyConditionIndex consIx = numCons[i];
                const AssemblyCondition& cond = 
                    *assembler.conditions[consIx];
                const int stat = cond.calcErrors
                   (assembler.getInternalState(), fy(nxtSlot, nEqns[i]));
                if (stat != 0)
                    return stat;
                nxtSlot += nEqns[i];
            }

            assert(nxtSlot == totalNEqns); // must use all slots
            return 0;
        }
    private:
        Assembler&                              assembler;
        const Array_<AssemblyConditionIndex>&   numCons;
        const Array_<int>&                      nEqns;
        const int                               totalNEqns;
    };

    int constraintJacobian(const Vector&    parameters, 
                           bool             new_parameters, 
                           Matrix&          J) const override 
    {   ++nEvalJacobian;

        if (new_parameters)
            setInternalStateFromFreeQs(parameters);
        for (unsigned i=0; i < assembler.reporters.size(); ++i)
            assembler.reporters[i]->handleEvent(getInternalState());

        assert(J.nrow() == getNumEqualityConstraints());
        assert(J.ncol() == getNumFreeQs());

        const int n = getNumFreeQs();

        // This will record the indices of any constraints we encounter that 
        // can't provide their own gradients; we'll handle them all together 
        // at the end.
        Array_<AssemblyConditionIndex> needNumericalJacobian;
        Array_<int>                    firstEqn;
        Array_<int>                    nEqns;
        int                            needy = 0;

        int nxtEqn = 0;
        for (unsigned i=0; i < assembler.errors.size(); ++i) {
            AssemblyConditionIndex   consIx = assembler.errors[i];
            const AssemblyCondition& cond   = *assembler.conditions[consIx];
            const int m = cond.getNumErrors(getInternalState());
            const int stat = (assembler.forceNumericalJacobian 
                              ? -1 
                              : cond.calcErrorJacobian(getInternalState(),
                                                       J(nxtEqn,0,m,n)));
            if (stat == -1) {
                needNumericalJacobian.push_back(consIx);
                firstEqn.push_back(nxtEqn);
                nEqns.push_back(m);
                needy += m;
            } else if (stat != 0)
                return stat;
            nxtEqn += m;
        }

        if (!needNumericalJacobian.empty()) {
            //cout << "Need numerical Jacobian for " 
            //     << needNumericalJacobian.size() << " constraints." << endl;
            NumJacobianFunc numCons(assembler, needNumericalJacobian, 
                                    nEqns, needy);
            // Forward difference should be fine here, unlike for the
            // gradient because we converge on the solution value 
            // rather than the derivative norm.
            Differentiator jacNumCons(numCons);
            Matrix numJ = jacNumCons.calcJacobian(getFreeQsFromInternalState());
            nEvalConstraints += jacNumCons.getNumCallsToUserFunction();

            // Fill in the missing rows.
            int nxtInNumJ = 0;
            for (unsigned i=0; i < needNumericalJacobian.size(); ++i) {
                J(firstEqn[i],0,nEqns[i],n) = numJ(nxtInNumJ,0,nEqns[i],n);
                nxtInNumJ += nEqns[i];
            }

        }

        //cout << "J=" << J;
        return 0;
    }

    int getNumObjectiveEvals()  const {return nEvalObjective;}
    int getNumConstraintEvals() const {return nEvalConstraints;}
    int getNumGradientEvals()   const {return nEvalGradient;}
    int getNumJacobianEvals()   const {return nEvalJacobian;}

    void resetStats() const { // stats are mutable
        nEvalObjective=nEvalConstraints=nEvalGradient=nEvalJacobian=0;
    }
private:
    const MultibodySystem& getSystem() const 
    {   return assembler.getMultibodySystem(); }
    const State& getInternalState() const 
    {   return assembler.getInternalState(); }
    int getNumFreeQs() const {return assembler.getNumFreeQs();}
    QIndex getQIndexOfFreeQ(FreeQIndex fx) const 
    {   return assembler.getQIndexOfFreeQ(fx);}
    FreeQIndex getFreeQIndexOfQ(QIndex qx) const
    {   return assembler.getFreeQIndexOfQ(qx); }
    Vector getFreeQsFromInternalState() const 
    {   return assembler.getFreeQsFromInternalState(); }
    void setInternalStateFromFreeQs(const Vector& freeQs) const 
    {   assembler.setInternalStateFromFreeQs(freeQs); }

    Assembler& assembler;

    mutable int nEvalObjective;
    mutable int nEvalConstraints;
    mutable int nEvalGradient;
    mutable int nEvalJacobian;
};



//------------------------------------------------------------------------------
//                                 ASSEMBLER
//------------------------------------------------------------------------------
Assembler::Assembler(const MultibodySystem& system)
:   system(system), accuracy(0), tolerance(0), // i.e., 1e-3, 1e-4
    forceNumericalGradient(false), forceNumericalJacobian(false), 
    useRMSErrorNorm(false), alreadyInitialized(false), 
    asmSys(0), optimizer(0), nAssemblySteps(0), nInitializations(0)
{
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    matter.convertToEulerAngles(system.getDefaultState(),
                                internalState);
    system.realizeModel(internalState);

    // Make sure the System's Constraints are always present; this sets the
    // weight to Infinity which makes us treat this as an assembly error
    // rather than merely a goal; that can be changed by the user.
    systemConstraints = adoptAssemblyError(new BuiltInConstraints());
}


Assembler::~Assembler() {
    uninitialize();
    // To be polite, and to show off, delete in reverse order of allocation 
    // (this is easier on the heap system).
    Array_<AssemblyCondition*,AssemblyConditionIndex>::reverse_iterator p;
    for (p = conditions.rbegin(); p != conditions.rend(); ++p)
        delete *p;
}


AssemblyConditionIndex Assembler::
adoptAssemblyError(AssemblyCondition* p) {
    return adoptAssemblyGoal(p, Infinity);
}

AssemblyConditionIndex Assembler::
adoptAssemblyGoal(AssemblyCondition* p, Real weight) {
    SimTK_ERRCHK_ALWAYS(p != 0, "Assembler::adoptAssemblyGoal()",
        "Null assembly condition pointer.");
    SimTK_ERRCHK1_ALWAYS(weight >= 0, "Assembler::adoptAssemblyGoal()",
        "Illegal assembly goal weight %g.", weight);

    uninitialize();

    const AssemblyConditionIndex acx(conditions.size());
    assert(conditions.size() == weights.size());
    p->setAssembler(*this, acx);
    conditions.push_back(p);
    weights.push_back(weight);
    return acx;
}


void Assembler::initialize() const {
    if (alreadyInitialized)
        return;

    ++nInitializations;

    Array_<QIndex> toBeLocked;
    reinitializeWithExtraQsLocked(toBeLocked);
    alreadyInitialized = true;
    return;

    /*NOTREACHED*/
    // TODO: This currently unused code would allow the Assembler to lock out 
    // variables that it thinks aren't worth bothering with. Needs real-world
    // testing and probably some override options. And should there be a
    // desperation mode where all variables are tried if we can't assemble
    // with some of them removed?
    Vector grad = abs(asmSys->calcCurrentGradient());
    Real maxGrad = 0;
    for (FreeQIndex fx(0); fx < grad.size(); ++fx)
        maxGrad = std::max(maxGrad, grad[fx]);
    if (maxGrad == 0) // no q does anything; probably no objective
        maxGrad = Infinity; // no q will be kept for gradient purposes

    Matrix jac = asmSys->calcCurrentJacobian();
    Vector colNorm(getNumFreeQs());
    Real maxJac = 0;
    for (FreeQIndex fx(0); fx < grad.size(); ++fx) {
        colNorm[fx] = jac(fx).norm();
        maxJac = std::max(maxJac, colNorm[fx]);
    }
    if (maxJac == 0) // no q does anything; probably no constraints
        maxJac = Infinity; // no q will be kept for Jacobian purposes

    const Real QTol = SqrtEps;
    const Real minGradAllowed = maxGrad*QTol;
    const Real minJacAllowed = maxJac*QTol;
    for (FreeQIndex fx(0); fx < grad.size(); ++fx)
        if (grad[fx] < minGradAllowed && colNorm[fx] < minJacAllowed)
            toBeLocked.push_back(getQIndexOfFreeQ(fx));

    if (toBeLocked.size()) {
        cout << "Reinitializing with these q's locked:\n";
        cout << toBeLocked << endl;
        reinitializeWithExtraQsLocked(toBeLocked);
        alreadyInitialized = true;
    }
}

void Assembler::reinitializeWithExtraQsLocked
   (const Array_<QIndex>& toBeLocked) const {
    uninitialize();

    const SimbodyMatterSubsystem& matter = getMatterSubsystem();

    system.realize(internalState, Stage::Instance);
    const int nq = internalState.getNQ();

    // Initialized locked q's to all those that the user locked, plus 
    // the extras.
    q2FreeQ.resize(nq); // no q has an associated freeQ at this point

    extraQsLocked = toBeLocked;
    lockedQs.insert(extraQsLocked.begin(), extraQsLocked.end());

    // Find all the mobilizers that have prescribed positions and lock
    // all their q's.
    for (MobodIndex mbx(0); mbx < matter.getNumBodies(); ++mbx) {
        const MobilizedBody& mobod  = matter.getMobilizedBody(mbx);
        if (mobod.getQMotionMethod(internalState) == Motion::Free)
            continue;
        const QIndex         q0     = mobod.getFirstQIndex(internalState);
        const int            nq     = mobod.getNumQ(internalState);
        for (int i=0; i<nq; ++i)
            lockedQs.insert(QIndex(q0+i));
    }

    // Lock all the q's for locked mobilizers.
    for (LockedMobilizers::const_iterator p = userLockedMobilizers.begin();
         p != userLockedMobilizers.end(); ++p)
    {
        const MobilizedBody& mobod  = matter.getMobilizedBody(*p);
        const QIndex         q0     = mobod.getFirstQIndex(internalState);
        const int            nq     = mobod.getNumQ(internalState);
        for (int i=0; i<nq; ++i)
            lockedQs.insert(QIndex(q0+i));
    }

    // Next add in all the q's that were individually locked.
    for (LockedQs::const_iterator p = userLockedQs.begin();
         p != userLockedQs.end(); ++p)
    {
        const MobilizedBodyIndex mbx = p->first;
        const MobilizedBody& mobod  = matter.getMobilizedBody(mbx);
        const QIndex         q0     = mobod.getFirstQIndex(internalState);
        const int            nq     = mobod.getNumQ(internalState);

        const QSet& qs = p->second;
        for (QSet::const_iterator qp = qs.begin(); qp != qs.end(); ++qp) {
            const MobilizerQIndex q = *qp;
            SimTK_ERRCHK3_ALWAYS(q < nq, "Assembler::initialize()", 
                "An attempt was made to lock q[%d] (numbering from 0) of"
                " mobilized body %d but that mobilizer has only %d q(s).",
                (int)q, (int)mbx, nq);
            lockedQs.insert(QIndex(q0 + q));
        }
    }

    const int nlockedq = (int)lockedQs.size();
    const int nfreeq   = std::max(nq - nlockedq, 0);

    // Find all the free q's and fill in the maps from FreeQIndex
    // to full QIndex, and from QIndex to FreeQIndex.
    freeQ2Q.resize(nfreeq);
    if (nlockedq) {
        FreeQIndex nxtFree(0);
        for (QIndex qx(0); qx < nq; ++qx) {
            if (lockedQs.find(qx) == lockedQs.end()) {
                q2FreeQ[qx] = nxtFree;
                freeQ2Q[nxtFree++] = qx;
            }
        }
    } else // all q's are free
        for (QIndex qx(0); qx < nq; ++qx) {
            q2FreeQ[qx] = FreeQIndex(qx);
            freeQ2Q[FreeQIndex(qx)] = qx;
        }

    // If *any* of the free q's have a restricted range, we'll provide a
    // (lower,upper) range value for *all* of them, using (-Inf,Inf) by
    // default. It might turn out that all the restricted q's are currently
    // locked so we'll hold off allocating the bounds arrays until we 
    // actually see a freeQ that's restricted.
    bool foundAnyRestricted = false;
    for (RestrictedQs::const_iterator p = userRestrictedQs.begin();
         p != userRestrictedQs.end(); ++p)
    {
        const MobilizedBodyIndex mbx    = p->first;
        const MobilizedBody&     mobod  = matter.getMobilizedBody(mbx);
        const QIndex             q0     = mobod.getFirstQIndex(internalState);
        const int                nq     = mobod.getNumQ(internalState);

        // Run through each of the q's that was restricted for this mobilizer.
        const QRanges& qranges = p->second;
        for (QRanges::const_iterator qr = qranges.begin(); 
             qr != qranges.end(); ++qr) 
        {
            const MobilizerQIndex q = qr->first;
            const Vec2&           r = qr->second;

            SimTK_ERRCHK3_ALWAYS(q < nq, "Assembler::initialize()", 
                "An attempt was made to restrict q[%d] (numbering from 0) of"
                " mobilized body %d but that mobilizer has only %d q(s).",
                (int)q, (int)mbx, nq);

            const QIndex     qx = QIndex(q0 + q);
            const FreeQIndex fx = q2FreeQ[qx];
            if (!fx.isValid())
                continue; // this q is locked; no need to restrict it

            if (!foundAnyRestricted) { // this is first one
                lower.resize(nfreeq); lower = -Infinity;    // allocate and initialize
                upper.resize(nfreeq); upper =  Infinity;
                foundAnyRestricted = true;
            }
            lower[fx] = r[0];
            upper[fx] = r[1];
        }
    }

    system.realize(internalState, Stage::Position);

    // Set up the lists of errors and goals based on the weights 
    // currently assigned to assembly conditions, and initialize the 
    // conditions as they are added.
    errors.clear(); nTermsPerError.clear(); goals.clear();
    assert(conditions.size() == weights.size());

    int nErrorTerms = 0;
    for (AssemblyConditionIndex acx(0); acx < conditions.size(); ++acx) {
        assert(conditions[acx] != 0 && weights[acx] >= 0);
        if (weights[acx] == 0) 
            continue;
        conditions[acx]->initializeCondition();
        if (weights[acx] == Infinity) {
            const int n = conditions[acx]->getNumErrors(internalState);
            if (n == 0)
                continue; // never mind; no constraint errors
            nErrorTerms += n;
            errors.push_back(acx);
            nTermsPerError.push_back(n);
        } else                        
            goals.push_back(acx);
    }

    // Allocate an AssemblerSystem which is in the form of an objective
    // function for the SimTK::Optimizer class.
    asmSys = new AssemblerSystem(*const_cast<Assembler*>(this));
    asmSys->setNumEqualityConstraints(nErrorTerms);
    if (lower.size())
        asmSys->setParameterLimits(lower, upper);


    // Optimizer will choose LBFGS for unconstrained (or just bounds-constrained)
    // problems, InteriorPoint for constrained problems.
    optimizer = new Optimizer(*asmSys
        //,InteriorPoint
        //,LBFGS
        //,LBFGSB
        );
    //optimizer->useNumericalGradient(true); // of goals
    //optimizer->useNumericalJacobian(true); // of errors

    // The size of the limited memory history affects the various optimizers
    // differently; I found this to be a good compromise. Smaller or larger
    // can both cause degraded performance.
    optimizer->setLimitedMemoryHistory(50);
    optimizer->setDiagnosticsLevel(0);
    optimizer->setMaxIterations(3000);
}

// Clean up all the mutable stuff; don't touch any user-set members.
void Assembler::uninitialize() const {
    if (!alreadyInitialized)
        return;

    alreadyInitialized = false;
    nAssemblySteps = 0;
    delete optimizer; optimizer = 0;
    delete asmSys; asmSys = 0;
    // Run through conditions in reverse order when uninitializing them; 
    // watch out: negative index not allowed so it is easier to use a reverse
    // iterators.
    Array_<AssemblyCondition*,AssemblyConditionIndex>::const_reverse_iterator p;
    for (p = conditions.crbegin(); p != conditions.crend(); ++p)
        (*p)->uninitializeCondition();
    goals.clear();
    nTermsPerError.clear();
    errors.clear();
    lower.clear(); upper.clear();
    freeQ2Q.clear();
    q2FreeQ.clear();
    lockedQs.clear();
    extraQsLocked.clear();
}

Real Assembler::calcCurrentGoal() const {
    initialize();
    return asmSys->calcCurrentGoal();
}

// Return norm of constraint errors, using the appropriate norm.
Real Assembler::calcCurrentErrorNorm() const {
    initialize();
    const int nc = asmSys->getNumEqualityConstraints();
    if (nc == 0) return 0;
    const Vector errs = asmSys->calcCurrentErrors();
    return useRMSErrorNorm
        ? std::sqrt(~errs*errs / errs.size())   // RMS
        : max(abs(errs));                       // infinity norm
}

Real Assembler::assemble() {
    initialize();
    ++nAssemblySteps;

    const Real initialErrorNorm = calcCurrentErrorNorm();
    const Real initialGoalValue = calcCurrentGoal(); // squared, >=0

    // Short circuit if this is already good enough.
    if (initialErrorNorm <= getErrorToleranceInUse()) {
        // Already feasible. Is the goal good enough too? Note that we don't
        // know much about the goal units, and "accuracy" is a unitless 
        // fraction. So we're going to use "tolerance" here. We do know that
        // the goal is squared but error tolerance is not. So this is likely
        // to be a very strict test if the error tolerance is tight!
        if (initialGoalValue <= square(getErrorToleranceInUse())) {
            for (unsigned i=0; i < reporters.size(); ++i)
                reporters[i]->handleEvent(internalState);
            return initialGoalValue;
        }
        // Not short circuiting.
    }


    const int nfreeq = getNumFreeQs();
    const int nqerr = internalState.getNQErr();


   // std::cout << "assemble(): initial tol/goal is " 
     //         << calcCurrentError() << "/" << calcCurrentGoal() << std::endl;

    for (unsigned i=0; i < reporters.size(); ++i)
        reporters[i]->handleEvent(internalState);

    // First step: satisfy prescribed motion (exactly).
    system.realize(internalState, Stage::Time);
    system.prescribeQ(internalState);

    // Optimize

    // Save the starting solution so we can restore it if the optimizer makes
    // it worse, which IpOpt has been observed to do.
    const Vector initialFreeQs = getFreeQsFromInternalState();
    Vector freeQs = initialFreeQs;
    // Use tolerance if there are any error conditions, else accuracy.
    optimizer->setConvergenceTolerance(getAccuracyInUse());
    optimizer->setConstraintTolerance(getErrorToleranceInUse());
    try
    {   optimizer->optimize(freeQs); }
    catch (const std::exception& e)
    {   setInternalStateFromFreeQs(freeQs); // realizes to Stage::Position

        // Sometimes the optimizer will throw an exception after it has
        // already achieved a winning solution. One message that comes up
        // is "Ipopt: Restoration failed (status -2)". We'll ignore that as
        // long as we have a good result. Otherwise we'll re-throw here.
        if (calcCurrentErrorNorm() > getErrorToleranceInUse()) {
            SimTK_THROW3(AssembleFailed, 
                (String("Optimizer failed with message: ") + e.what()).c_str(), 
                calcCurrentErrorNorm(), getErrorToleranceInUse());
        }
    }

    // This will ensure that the internalState has its q's set to match the
    // parameters.
    setInternalStateFromFreeQs(freeQs);

    for (unsigned i=0; i < reporters.size(); ++i)
        reporters[i]->handleEvent(internalState);

    Real tolAchieved = calcCurrentErrorNorm();
    Real goalAchieved = calcCurrentGoal();

    // See if we should just revert to the initial solution.
    if (   initialErrorNorm <= getErrorToleranceInUse() // started feasible
        && goalAchieved > initialGoalValue)             // objective got worse
    {
        setInternalStateFromFreeQs(initialFreeQs);
        tolAchieved = initialErrorNorm;
        goalAchieved = initialGoalValue;
    }

    if (tolAchieved > getErrorToleranceInUse())
        SimTK_THROW3(AssembleFailed, 
            "Unable to achieve required assembly error tolerance.",
            tolAchieved, getErrorToleranceInUse());

    //std::cout << "assemble(): final tol/goal is " 
    //          << calcCurrentError() << "/" << calcCurrentGoal() << std::endl;

    return goalAchieved;
}


Real Assembler::track(Real frameTime) {
    initialize();
    ++nAssemblySteps;

    if (frameTime >= 0 && internalState.getTime() != frameTime) {
        internalState.setTime(frameTime);
        system.realize(internalState, Stage::Time);
        // Satisfy prescribed motion (exactly).
        system.prescribeQ(internalState);
        system.realize(internalState, Stage::Position); 
    }

    const Real initialErrorNorm = calcCurrentErrorNorm();
    const Real initialGoalValue = calcCurrentGoal(); // squared!

    // Short circuit if this is already good enough.
    if (initialErrorNorm <= getErrorToleranceInUse()) {
        // See comments in assemble() for more info.
        if (initialGoalValue <= square(getErrorToleranceInUse())) {
            for (unsigned i=0; i < reporters.size(); ++i)
                reporters[i]->handleEvent(internalState);
            return initialGoalValue;
        }
        // Not short circuiting.
    }


    const int nfreeq = getNumFreeQs();
    const int nqerr = internalState.getNQErr();


    // std::cout << "track(): initial tol/goal is " 
    //         << calcCurrentError() << "/" << calcCurrentGoal() << std::endl;

    // Optimize
    Vector freeQs = getFreeQsFromInternalState();
    optimizer->setConvergenceTolerance(getAccuracyInUse());
    optimizer->setConstraintTolerance(getErrorToleranceInUse());
    try
    {   optimizer->optimize(freeQs); }
    catch (const std::exception& e)
    {   setInternalStateFromFreeQs(freeQs); // realizes to Stage::Position

        // Sometimes the optimizer will throw an exception after it has
        // already achieved a winning solution. One message that comes up
        // is "Ipopt: Restoration failed (status -2)". We'll ignore that as
        // long as we have a good result. Otherwise we'll re-throw here.
        if (calcCurrentErrorNorm() > getErrorToleranceInUse()) {
            SimTK_THROW3(TrackFailed, 
                (String("Optimizer failed with message: ") + e.what()).c_str(), 
                calcCurrentErrorNorm(), getErrorToleranceInUse());
        }
    }

    // This will ensure that the internalState has its q's set to match the
    // parameters.
    // This will ensure that the internalState has its q's set to match the
    // parameters.
    setInternalStateFromFreeQs(freeQs);

    for (unsigned i=0; i < reporters.size(); ++i)
        reporters[i]->handleEvent(internalState);

    const Real tolAchieved = calcCurrentErrorNorm();
    if (tolAchieved > getErrorToleranceInUse())
        SimTK_THROW3(TrackFailed, 
            "Unable to achieve required assembly error tolerance.",
            tolAchieved, getErrorToleranceInUse());

    //std::cout << "track(): final tol/goal is " 
    //          << calcCurrentError() << "/" << calcCurrentGoal() << std::endl;

    return calcCurrentGoal();
}

int Assembler::getNumGoalEvals()  const 
{   return asmSys ? asmSys->getNumObjectiveEvals() : 0;}
int Assembler::getNumErrorEvals() const
{   return asmSys ? asmSys->getNumConstraintEvals() : 0;}
int Assembler::getNumGoalGradientEvals()   const
{   return asmSys ? asmSys->getNumGradientEvals() : 0;}
int Assembler::getNumErrorJacobianEvals()   const
{   return asmSys ? asmSys->getNumJacobianEvals() : 0;}
int Assembler::getNumAssemblySteps() const
{   return nAssemblySteps; }
int Assembler::getNumInitializations() const
{   return nInitializations; }
void Assembler::resetStats() const
{   if (asmSys) asmSys->resetStats(); 
    nAssemblySteps = nInitializations = 0; }

