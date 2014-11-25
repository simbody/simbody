#ifndef SimTK_SIMBODY_ASSEMBLY_CONDITION_H_
#define SimTK_SIMBODY_ASSEMBLY_CONDITION_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Assembler.h"

#include <map>

namespace SimTK {

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

} // namespace SimTK

#endif // SimTK_SIMBODY_ASSEMBLY_CONDITION_H_
