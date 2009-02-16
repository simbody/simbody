#ifndef SimTK_SimTKCOMMON_STAGE_H_
#define SimTK_SimTKCOMMON_STAGE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Exception.h"
#include "SimTKcommon/internal/Enumeration.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdarg>

namespace SimTK {

/**
 * This class is basically a glorified enumerated type, type-safe and range
 * checked but permitting convenient (if limited) arithmetic.
 * Constants look like Stage::Position, and loops can be written like
 * 		for(Stage s=Stage::Lowest; s <= Stage::Highest; ++s) ...
 * Stage constants (of type Stage::Num) are implicitly converted to type
 * Stage when necessary.
 */	
class SimTK_SimTKCOMMON_EXPORT Stage : public Enumeration<Stage> {
public:
    // sherm 060720
    // Noun version          Verb version
    //  Invalid                InvalidType
    //  Empty                  Allocate    Allocated
    //  Topology               Build       Built
    //
    //  Model                  Model       Modeled
    //  Instance               Specify     Specified
    //  Time                   Time        Timed
    //  Position               Position    Positioned
    //  Velocity               Move        Moving
    //  Dynamics               ?
    //  Acceleration           Accelerate  Accelerated
    //
    //  Report                 Report      Reported
    
	enum Num {
        EmptyIndex          =  0, // TODO: Initialized, Unbuilt, Empty, Vacant?
        TopologyIndex       =  1, // TODO: Constructed, Finalized?
        ModelIndex          =  2, // TODO: Instantiated, Resourced, Provisioned, Specialized?
        InstanceIndex       =  3, // TODO: Instanced, Specified?
        TimeIndex           =  4,
        PositionIndex       =  5, // TODO: Positioned?
        VelocityIndex       =  6, // TODO: Velocity, Speed, Rate?
        DynamicsIndex       =  7, // forces, dynamic properties & operators available
        AccelerationIndex   =  8, // TODO: Accelerated?
        ReportIndex         =  9  // TODO: Output?
	};
    static const Stage Empty;
    static const Stage Topology;
    static const Stage Model;
    static const Stage Instance;
    static const Stage Time;
    static const Stage Position;
    static const Stage Velocity;
    static const Stage Dynamics;
    static const Stage Acceleration;
    static const Stage Report;

    static const Stage LowestValid;
    static const Stage HighestValid;
    static const int NValid = ReportIndex-EmptyIndex+1;

    // LowestRuntime->HighestRuntime cover the post-construction stages only.
	static const Stage LowestRuntime;
	static const Stage HighestRuntime;
	static const int NRuntime = ReportIndex-ModelIndex+1;
		
	Stage	next() const { return getValue(getIndex()+1); }
    Stage	prev() const { return getValue(getIndex()-1); }

    // Set stage=min(stage, tooHigh-1).
    void invalidate(Stage tooHigh) {
        if (getIndex() >= tooHigh.getIndex())
            *this = tooHigh.prev();
    }
	
    bool isInRuntimeRange() const {return Stage::LowestRuntime <= getIndex() && getIndex() <= Stage::HighestRuntime;}

private:
    Stage();
    Stage(const Stage& thisElement, int index, const char* name);
    static void initValues();
    friend class Enumeration<Stage>;

};

namespace Exception {

class RealizeTopologyMustBeCalledFirst : public Base {
public:
    RealizeTopologyMustBeCalledFirst(const char* fn, int ln,
       const char* objectType, // e.g., "System", "Subsystem"
       const char* objectName, const char* methodName) : Base(fn,ln)
    {
        setMessage(String(methodName) + ": " + String(objectType) + " " + String(objectName)
           + " topology has not yet been realized -- must call realizeTopology() first");
    }
};

class StageTooLow : public Base {
public:
    StageTooLow(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be at least " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
};

class StageIsWrong : public Base {
public:
    StageIsWrong(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
};

class StageTooHigh : public Base {
public:
    StageTooHigh(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be less than " + targetStage.getName() + " in " + String(where)
           + " but current stage was " + currentStage.getName());
    }
};

class StageOutOfRange : public Base {
public:
    StageOutOfRange(const char* fn, int ln,
        Stage lower, Stage currentStage, Stage upper, const char* where) : Base(fn,ln)
    {
        setMessage("Expected (" + lower.getName() + " <= stage <= " + upper.getName() + ") in " + String(where)
           + " but stage was " + currentStage.getName());
    }
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
private:
};
}

inline std::ostream& operator<<(std::ostream& o, Stage g) { o << g.getName(); return o; }	


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
