#ifndef SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
#define SimTK_SimTKCOMMON_SYSTEM_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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

namespace SimTK {

class Subsystem;
class DecorativeGeometry;

//==============================================================================
//                             SYSTEM :: GUTS
//==============================================================================
/** This is the declaration for the System::Guts class, the abstract object to
which a System handle points. This is in a separate header file from System
because only the very few people who are extending the System class to make 
their own Systems need to be aware of the details. End users access only methods
from the System class and classes derived from System, never anything from
System::Guts or its derived classes. **/
class SimTK_SimTKCOMMON_EXPORT System::Guts {
public:

/** Destructor is virtual to permit cleanup of derived classes. **/
virtual ~Guts() {
    clearMyHandle();
    m_subsystems.clear();
    invalidateSystemTopologyCache();
}

// Copy constructor is protected; copy assignment is deleted.

const String& getName()    const {return m_systemName;}
const String& getVersion() const {return m_systemVersion;}

void setDefaultTimeScale(double tc) 
{   m_defaultTimeScale = tc; }
double getDefaultTimeScale() const {return m_defaultTimeScale;}

void setDefaultLengthScale(Real lc)
{   m_defaultLengthScale = lc; }
Real getDefaultLengthScale() const {return m_defaultLengthScale;}

void setUpDirection(const CoordinateDirection& up) 
{   m_defaultUpDirection = up; }
CoordinateDirection getUpDirection() const {return m_defaultUpDirection;}
void setUseUniformBackground(bool useUniform)
{   m_useUniformBackground = useUniform; }
bool getUseUniformBackground() const {return m_useUniformBackground;}


int getNumSubsystems() const {return (int)m_subsystems.size();}
const Subsystem& getSubsystem(SubsystemIndex i) const 
{   return m_subsystems[i]; }
Subsystem& updSubsystem(SubsystemIndex i)
{   return m_subsystems[i]; }

void setHasTimeAdvancedEvents(bool hasEm) 
{   m_hasTimeAdvancedEventsFlag = hasEm; }
bool hasTimeAdvancedEvents() const
{   return m_hasTimeAdvancedEventsFlag; }

    //////////////////////////////
    // EVALUATION (REALIZATION) //
    //////////////////////////////

// These are the routines to which the System class forwards requests.

const State& getDefaultState() const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS
        (systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::getDefaultState()");

    return m_defaultState;
}

State& updDefaultState() {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS
        (systemTopologyHasBeenRealized(),
        "System", getName(), "System::Guts::updDefaultState()");

    return m_defaultState;
}

void realize(const State& s, Stage g = Stage::HighestRuntime) const;

// Take over ownership from the Subsystem handle, allocate a new
// subsystem slot for it, and return the slot number. This is only 
// allowed if the supplied Subsystem already has Guts, but is
// NOT part of some other System.
SubsystemIndex adoptSubsystem(Subsystem& child);

// Obtain the owner handle for this System::Guts object.
const System& getSystem() const
{   assert(m_myHandle); return *m_myHandle; }
System& updSystem()
{   assert(m_myHandle); return *m_myHandle; }

void setOwnerHandle(System& handle) {m_myHandle = &handle;}
bool hasOwnerHandle() const {return m_myHandle != nullptr;}
void clearMyHandle() {m_myHandle=nullptr;}

bool systemTopologyHasBeenRealized() const 
{   return m_systemTopologyRealized; }

StageVersion getSystemTopologyCacheVersion() const
{   return m_topologyCacheVersion; }

// Use this cautiously if at all!
void setSystemTopologyCacheVersion(StageVersion version) const {
    assert(version>0); 
    auto mThis = const_cast<System::Guts*>(this);
    mThis->m_topologyCacheVersion = version; 
}

// Invalidating the System topology cache requires invalidating all
// Subsystem topology caches also so that the next realizeTopology()
// will cause them to request their State resources again so we can
// build up the defaultState.
void invalidateSystemTopologyCache() const {
    if (m_systemTopologyRealized) {
        auto mThis = const_cast<System::Guts*>(this);
        mThis->invalidateCache();
    }
}


// Wrap the cloneImpl virtual method.
System::Guts* clone() const {return cloneImpl();}

// These routines wrap the virtual realize...Impl() methods to ensure
// good behavior such as checking that stage requirements are met and
// updating the stage at the end. Note that these will do nothing if
// the System stage is already at or greater than the indicated stage.

const State& realizeTopology() const;
void realizeModel(State&) const;
void realizeInstance    (const State& s) const;
void realizeTime        (const State& s) const;
void realizePosition    (const State& s) const;
void realizeVelocity    (const State& s) const;
void realizeDynamics    (const State& s) const;
void realizeAcceleration(const State& s) const;
void realizeReport      (const State& s) const;

// These wrap the other virtual methods.
void multiplyByN(const State& state, const Vector& u, 
                    Vector& dq) const;
void multiplyByNTranspose(const State& state, const Vector& fq, 
                            Vector& fu) const;
void multiplyByNPInv(const State& state, const Vector& dq, 
                        Vector& u) const;
void multiplyByNPInvTranspose(const State& state, const Vector& fu, 
                                Vector& fq) const;

bool prescribeQ(State&) const;
bool prescribeU(State&) const;
void getFreeQIndex(const State&, Array_<SystemQIndex>& freeQs) const;
void getFreeUIndex(const State&, Array_<SystemUIndex>& freeUs) const;

void projectQ(State&, Vector& qErrEst, 
                const ProjectOptions& options, ProjectResults& results) const;
void projectU(State&, Vector& uErrEst, 
                const ProjectOptions& options, ProjectResults& results) const;

void calcDecorativeGeometryAndAppend(const State&, Stage, 
                                        Array_<DecorativeGeometry>&) const;


//------------------------------------------------------------------------------
                                  protected:
 
// The destructor is already virtual; see above.

virtual System::Guts* cloneImpl() const = 0;

// Override these to change the evaluation order of the Subsystems.
// Except for Subsystem(0), the SystemGlobalSubsystem, the default is to 
// evaluate subsystems in increasing order of SubsystemIndex.
// These methods should not be called directly; they are invoked by the
// above wrapper methods. Note: the wrappers *will not* call these
// routines if the system stage has already met the indicated stage level.
// If fact these routines will be called only when the system stage
// is at the level just prior to the one indicated here. For example,
// realizeVelocityImpl() will be called only if the passed-in State
// has been determined to have its system stage exactly Stage::Position.
// Your implementation does not need to repeat those checks.

virtual int realizeTopologyImpl(State& state)       const {return 0;}
virtual int realizeModelImpl   (State& state)       const {return 0;}
virtual int realizeInstanceImpl(const State& state) const {return 0;}
virtual int realizeTimeImpl    (const State& state) const {return 0;}
virtual int realizePositionImpl(const State& state) const {return 0;}
virtual int realizeVelocityImpl(const State& state) const {return 0;}
virtual int realizeDynamicsImpl(const State& state) const {return 0;}
virtual int realizeAccelerationImpl(const State& state) const {return 0;}
virtual int realizeReportImpl  (const State& state) const {return 0;}

virtual void multiplyByNImpl(const State& state, const Vector& u, 
                                Vector& dq) const;
virtual void multiplyByNTransposeImpl(const State& state, const Vector& fq, 
                                        Vector& fu) const;
virtual void multiplyByNPInvImpl(const State& state, const Vector& dq, 
                                    Vector& u) const;
virtual void multiplyByNPInvTransposeImpl(const State& state, 
                                            const Vector& fu, 
                                            Vector& fq) const;

// Defaults assume no prescribed motion; hence, no change made.
virtual bool prescribeQImpl(State&) const {return false;}
virtual bool prescribeUImpl(State&) const {return false;}


// Defaults assume no constraints and return success meaning "all 
// constraints satisfied".
virtual void projectQImpl(State& state, Vector& qErrEst, 
            const ProjectOptions& options, ProjectResults& results) const
{   results.clear(); results.setExitStatus(ProjectResults::Succeeded); }
virtual void projectUImpl(State& state, Vector& uErrEst, 
            const ProjectOptions& options, ProjectResults& results) const
{   results.clear(); results.setExitStatus(ProjectResults::Succeeded); }


// Default is that all the state variables are free.
virtual void getFreeQIndexImpl
    (const State& s, Array_<SystemQIndex>& freeQs) const {
    const unsigned nq = (unsigned)s.getNQ();
    freeQs.resize(nq);
    for (unsigned i=0; i<nq; ++i)
        freeQs[i] = SystemQIndex(i);
}
virtual void getFreeUIndexImpl
    (const State& s, Array_<SystemUIndex>& freeUs) const  {
    const unsigned nu = (unsigned)s.getNU();
    freeUs.resize(nu);
    for (unsigned i=0; i<nu; ++i)
        freeUs[i] = SystemUIndex(i);
}

// CAREFUL: can't use defaults here so make sure these are updated if
// you change or reorder any data members.

// Note that this serves as a default constructor since both arguments have 
// defaults.
explicit Guts(const String& name="<NONAME>", 
                const String& version="0.0.0") 
:   m_myHandle(nullptr), 
    m_systemName(name), m_systemVersion(version),
    m_defaultTimeScale(Real(0.1)), 
    m_defaultLengthScale(Real(1)),
    m_defaultUpDirection(YAxis), 
    m_useUniformBackground(false),
    m_hasTimeAdvancedEventsFlag(false),
    m_topologyCacheVersion(1) // zero isn't allowed
{
    clearCache();
    resetAllCounters();
}

// copies the base class; for use from derived class copy constructors
Guts(const Guts& src)
:   m_myHandle(nullptr),
    m_systemName(src.m_systemName), m_systemVersion(src.m_systemVersion),
    m_subsystems(src.m_subsystems),
    m_defaultTimeScale(src.m_defaultTimeScale),
    m_defaultLengthScale(src.m_defaultLengthScale),
    m_defaultUpDirection(src.m_defaultUpDirection), 
    m_useUniformBackground(src.m_useUniformBackground),
    m_hasTimeAdvancedEventsFlag(src.m_hasTimeAdvancedEventsFlag),
    m_topologyCacheVersion(src.m_topologyCacheVersion)
{
    clearCache();
    resetAllCounters();
}

Guts& operator=(const Guts&) = delete;


//------------------------------------------------------------------------------
                                    private:
friend class System;

void invalidateCache() {
    // Mark system topology invalid *first* so that the invalidate
    // subsystem calls below don't recurse back here!
    clearCache();
    for (SubsystemIndex i(0); i < (int)m_subsystems.size(); ++i)
        m_subsystems[i].invalidateSubsystemTopologyCache();
}

void clearCache() {
    m_systemTopologyRealized = false;
    ++m_topologyCacheVersion;
    m_defaultState.clear();
}
    
    // TOPOLOGY STAGE STATE //

ReferencePtr<System>    m_myHandle;     // the owner handle of these guts

String                  m_systemName, m_systemVersion;
StableArray<Subsystem>  m_subsystems;

// Scaling hints
double                  m_defaultTimeScale;       // units of time
Real                    m_defaultLengthScale;     // units of length

CoordinateDirection     m_defaultUpDirection;     // visualization hint
bool                    m_useUniformBackground;   // visualization hint

//TODO: should be in State as a Model variable
bool                    m_hasTimeAdvancedEventsFlag; 
   
    // TOPOLOGY STAGE CACHE //

// This should only be true when *all* subsystems have successfully
// completed realizeTopology(). Anything which invalidates topology for
// one of the contained subsystem must invalidate topology for the system
// as a whole also.
bool                    m_systemTopologyRealized;
StageVersion            m_topologyCacheVersion;

// This is only meaningful if systemTopologyRealized==true. In that case
// its Topology stage version will match the above. A State with a different
// Topology version cannot be used with this Subsystem.
State                   m_defaultState;

    // MUTABLE STATISTICS //

mutable int nRealizationsOfStage[Stage::NValid];
mutable int nRealizeCalls; // counts realizeTopology(), realizeModel(), 
                            // realize()

mutable int nPrescribeQCalls, nPrescribeUCalls;

mutable int nProjectQCalls, nProjectUCalls;
mutable int nFailedProjectQCalls, nFailedProjectUCalls;
mutable int nQProjections, nUProjections; // the ones that did something
mutable int nQErrEstProjections, nUErrEstProjections;

void resetAllCounters() {
    for (int i=0; i<Stage::NValid; ++i)
        nRealizationsOfStage[i] = 0;
    nRealizeCalls = nPrescribeQCalls = nPrescribeUCalls = 0;
    nProjectQCalls = nProjectUCalls = 0;
    nFailedProjectQCalls = nFailedProjectUCalls = 0;
    nQProjections = nUProjections = 0;
    nQErrEstProjections = nUErrEstProjections = 0;
}
};


//==============================================================================
//                             SYSTEM INLINES
//==============================================================================
// These had to wait for System::Guts to be defined.

// None here but see SystemGlobalSubsystem.h.

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_SYSTEM_GUTS_H_
