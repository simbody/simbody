#ifndef SimTK_SimTKCOMMON_STAGE_H_
#define SimTK_SimTKCOMMON_STAGE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Exception.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdarg>

namespace SimTK {
 
/** This is the type to use for Stage and state variable version numbers that
get incremented whenever a state value changes. Whenever time or any state 
variable is modified, we increment the stage version for any stage that gets
invalidated. We also increment a separate version number for the state variable
that changes, so that cache entries can have finer-grained dependencies than
just on whole stages. -1 means "unintialized". 0 is never used as 
a %StageVersion, but is allowed as a cache value which is guaranteed never to 
look valid. **/
typedef long long StageVersion;

/** This class is basically a glorified enumerated type, type-safe and range
checked but permitting convenient (if limited) arithmetic. Constants look like 
Stage::Position, and loops can be written like
@code
    for(Stage s = Stage::LowestValid; s <= Stage::HighestValid; ++s) {
        // ...
    }
@endcode
Stage constants (of type Stage::Level) are implicitly converted to type
Stage when necessary.

Default construction gives Stage::Empty which really means "invalid". **/    
class Stage  {
public:
    enum Level {
        Empty          =  0, ///< Lower than any legitimate Stage.
        Topology       =  1, ///< System topology realized.
        Model          =  2, ///< Modeling choices made.
        Instance       =  3, ///< Physical parameters set.
        Time           =  4, ///< A new time has been realized.
        Position       =  5, ///< Spatial configuration available.
        Velocity       =  6, ///< Spatial velocities available.
        Dynamics       =  7, ///< Forces calculated.
        Acceleration   =  8, ///< Accelerations and multipliers calculated.
        Report         =  9, ///< Report-only quantities evaluated.
        Infinity       = 10, ///< Higher than any legitimate Stage.

        LowestValid = Empty,    ///< For iterating over all stage values.
        HighestValid = Infinity,
        LowestRuntime = Model,  ///< For iterating over meaningful stage values.
        HighestRuntime = Report
    };

    enum {
        NValid = HighestValid-LowestValid+1,
        NRuntime = HighestRuntime-LowestRuntime+1
    };
    
    /** Default construction gives Stage::Empty. **/
    Stage() : level(Stage::Empty) {}
    /** This is an implicit conversion from Stage::Level to Stage. **/
    Stage(Level l) {
        assert(LowestValid <= l && l <= HighestValid);
        level = l;
    }
    /** You can explicitly create a Stage from an int if it is in range. **/
    explicit Stage(int l) {
        assert(LowestValid <= l && l <= HighestValid);
        level = Level(l);
    }
    /** Stage will implicitly convert to int so you can use it as an index. **/
    operator int() const {return level;}

    bool operator==(Level other) const {return level==other;}
    bool operator!=(Level other) const {return level!=other;}
    bool operator<(Level other) const {return level<other;}
    bool operator<=(Level other) const {return level<=other;}
    bool operator>(Level other) const {return level>other;}
    bool operator>=(Level other) const {return level>=other;}
    bool operator==(Stage other) const {return level==other.level;}
    bool operator!=(Stage other) const {return level!=other.level;}
    bool operator<(Stage other) const {return level<other.level;}
    bool operator<=(Stage other) const {return level<=other.level;}
    bool operator>(Stage other) const {return level>other.level;}
    bool operator>=(Stage other) const {return level>=other.level;}

    // Prefix operators
    const Stage& operator++()
    {   assert(level<HighestValid); level=Level(level+1); return *this; }
    const Stage& operator--() 
    {   assert(level>LowestValid);  level=Level(level-1); return *this;}
    // Postfix operators
    Stage operator++(int)     
    {   assert(level<HighestValid); level=Level(level+1); return prev(); }
    Stage operator--(int)     
    {   assert(level>LowestValid);  level=Level(level-1); return next(); }

    /** Return the Stage following this one, with Stage::Infinity returned
    if this Stage is already at its highest value, Stage::Report. An exception
    is thrown if this Stage is already Stage::Infinity. **/
    Stage next() const 
    {   assert(level<HighestValid); return Stage(Level(level+1)); }
    /** Return the Stage before this one, with Stage::Empty returned
    if this Stage is already at its lowest value, Stage::Topology. An exception
    is thrown if this Stage is already Stage::Empty. **/
    Stage prev() const 
    {   assert(level>LowestValid); return Stage(Level(level-1)); }

    /** Return a printable name corresponding to the stage level currently
    stored in this Stage. **/
    String getName() const {
        switch (level) {
        case Empty:         return "Empty";    break;
        case Topology:      return "Topology"; break;
        case Model:         return "Model";    break;
        case Instance:      return "Instance"; break;
        case Time:          return "Time";     break;
        case Position:      return "Position"; break;
        case Velocity:      return "Velocity"; break;
        case Dynamics:      return "Dynamics"; break;
        case Acceleration:  return "Acceleration"; break;
        case Report:        return "Report";   break;
        case Infinity:      return "Infinity"; break;
        default: assert(!"Stage::getName(): illegal level");
        }
        return String("INVALID STAGE LEVEL ") + String(level);
    }

    /** Set this Stage=min(stageNow, tooHigh-1). **/
    void invalidate(Stage tooHigh) {
        if (level >= tooHigh.level)
            *this = tooHigh.prev();
    }

    /** Return true if this Stage has one of the meaningful values between
    Stage::Topology and Stage::Report, rather than one of the end markers
    Stage::Empty or Stage::Infinity. **/
    bool isInRuntimeRange() const 
    {   return    Stage::LowestRuntime <= level 
               && level <= Stage::HighestRuntime; }

private:
    Level level;
};




namespace Exception {

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4996) // don't warn about sprintf, etc.
#endif

class RealizeTopologyMustBeCalledFirst : public Base {
public:
    RealizeTopologyMustBeCalledFirst(const char* fn, int ln,
       const char* objectType, // e.g., "System", "Subsystem"
       const char* objectName, const char* methodName) : Base(fn,ln)
    {
        setMessage(String(methodName) + ": " + String(objectType) + " " + String(objectName)
           + " topology has not been realized since the last topological change"
             " -- you must call realizeTopology() first.");
    }
    virtual ~RealizeTopologyMustBeCalledFirst() throw() { }
};

class StateAndSystemTopologyVersionsMustMatch : public Base {
public:
    StateAndSystemTopologyVersionsMustMatch(const char* fn, int ln,
       const char* objectType, // e.g., "System", "Subsystem"
       const char* objectName, const char* methodName,
       int sysTopoVersion,
       int stateTopoVersion) : Base(fn,ln)
    {
        setMessage(String(methodName) 
        + ": The given State's Topology stage version number (" 
        + String(stateTopoVersion)
        + ") doesn't match the current topology cache version number (" 
        + String(sysTopoVersion)
        + ") of " + String(objectType) + " " + String(objectName) + "."
        + " That means there has been a topology change to this System since this"
          " State was created so they are no longer compatible. You should create"
          " a new State from the System's default State."
          " (Loopholes exist for advanced users.)");
    }
    virtual ~StateAndSystemTopologyVersionsMustMatch() throw() { }
};



class StageTooLow : public Base {
public:
    StageTooLow(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be at least " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
    virtual ~StageTooLow() throw() { }
};

class StageIsWrong : public Base {
public:
    StageIsWrong(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
    virtual ~StageIsWrong() throw() { }
};

class StageTooHigh : public Base {
public:
    StageTooHigh(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be less than " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
    virtual ~StageTooHigh() throw() { }
};

class StageOutOfRange : public Base {
public:
    StageOutOfRange(const char* fn, int ln,
        Stage lower, Stage currentStage, Stage upper, const char* where) : Base(fn,ln)
    {
        setMessage("Expected (" + lower.getName() + " <= stage <= " + upper.getName() + ") in " + String(where)
           + " but stage was " + currentStage.getName());
    }
    virtual ~StageOutOfRange() throw() { }
};

class CacheEntryOutOfDate : public Base {
public:
    CacheEntryOutOfDate(const char* fn, int ln,
        Stage currentStage, Stage dependsOn, 
        StageVersion dependsOnVersion, StageVersion lastCalculatedVersion) 
    :   Base(fn,ln)
    {
        setMessage("State Cache entry was out of date at Stage " + currentStage.getName() 
           + ". This entry depends on version " + String(dependsOnVersion) 
           + " of Stage " + dependsOn.getName() 
           + " but was last updated at version " + String(lastCalculatedVersion) + ".");
    }
    virtual ~CacheEntryOutOfDate() throw() { }
};

// An attempt to realize a particular subsystem to a particular stage failed.
class RealizeCheckFailed : public Base {
public:
    RealizeCheckFailed(const char* fn, int ln, Stage g, 
                       int subsystemId, const char* subsystemName,
                       const char* fmt ...) : Base(fn,ln)
    {
        char buf[1024];
        va_list args;
        va_start(args, fmt);
        vsprintf(buf, fmt, args);
        setMessage("Couldn't realize subsystem " + String(subsystemId)
                   + "(" + String(subsystemName) + ") to Stage "
                   + g.getName() + ": " + String(buf) + ".");
        va_end(args);
    }
    virtual ~RealizeCheckFailed() throw() { }
};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

} // namespace Exception

inline std::ostream& operator<<(std::ostream& o, Stage g) 
{   o << g.getName(); return o; }    


} // namespace SimTK

    // REALIZECHECKs: these should be used to catch and report problems that
    // occur when realizing a subsystem.

#define SimTK_REALIZECHECK_ALWAYS(cond,stage,subsysIx,subsysName,msg)       \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysIx),(subsysName),(msg));                 \
    }while(false)
#define SimTK_REALIZECHECK1_ALWAYS(cond,stage,subsysIx,subsysName,msg,a1)   \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysIx),(subsysName),(msg),(a1));            \
    }while(false)
#define SimTK_REALIZECHECK2_ALWAYS(cond,stage,subsysIx,subsysName,msg,a1,a2)\
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysIx),(subsysName),(msg),(a1),(a2));       \
    }while(false)
#define SimTK_REALIZECHECK3_ALWAYS(cond,stage,subsysIx,subsysName,msg,a1,a2,a3)     \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::RealizeCheckFailed,                \
                    (stage),(subsysIx),(subsysName),(msg),(a1),(a2),(a3));          \
    }while(false)
#define SimTK_REALIZECHECK4_ALWAYS(cond,stage,subsysIx,subsysName,msg,a1,a2,a3,a4)  \
    do{if(!(cond))SimTK_THROW8(SimTK::Exception::RealizeCheckFailed,                \
                    (stage),(subsysIx),(subsysName),(msg),(a1),(a2),(a3),(a4));     \
    }while(false)
#define SimTK_REALIZECHECK5_ALWAYS(cond,stage,subsysIx,subsysName,msg,a1,a2,a3,a4,a5)   \
    do{if(!(cond))SimTK_THROW9(SimTK::Exception::RealizeCheckFailed,                    \
                    (stage),(subsysIx),(subsysName),(msg),(a1),(a2),(a3),(a4),(a5));    \
    }while(false)

    
#endif // SimTK_SimTKCOMMON_STAGE_H_
