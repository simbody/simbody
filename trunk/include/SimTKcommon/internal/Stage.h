#ifndef SimTK_SimTKCOMMON_STAGE_H_
#define SimTK_SimTKCOMMON_STAGE_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Exception.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cstdarg>

namespace SimTK {

/**
 * This class is basically a glorified enumerated type, type-safe and range
 * checked but permitting convenient (if limited) arithmetic.
 * Constants look like Stage::Configure, and loops can be written like
 * 		for(Stage s=Stage::Lowest; s <= Stage::Highest; ++s) ...
 * Stage constants (of type Stage::Num) are implicitly converted to type
 * Stage when necessary.
 */	
class Stage {
public:
    // sherm 060720
    // Noun version          Verb version
    //  Invalid                Invalid
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
        Invalid        = -99,
        Empty          =  0, // TODO: Initialized, Unbuilt, Empty, Vacant?
        Topology       =  1, // TODO: Constructed, Finalized?
        Model          =  2, // TODO: Instantiated, Resourced, Provisioned, Specialized?
        Instance       =  3, // TODO: Instanced, Specified?
        Time           =  4,
        Position       =  5, // TODO: Positioned?
        Velocity       =  6, // TODO: Velocity, Speed, Rate?
        Dynamics       =  7, // forces, dynamic properties & operators available
        Acceleration   =  8, // TODO: Accelerated?
        Report         =  9  // TODO: Output?
	};
    static const Stage::Num Kinematics      = Velocity; // i.e., Position+Velocity

    static const Stage::Num LowestValid     = Empty;
    static const Stage::Num HighestValid    = Report;
    static const int        NValid          = HighestValid-LowestValid+1;

    // LowestRuntime->HighestRuntime cover the post-construction stages only.
	static const Stage::Num	LowestRuntime	= Model;
	static const Stage::Num	HighestRuntime	= Acceleration;
	static const int		NRuntime     	= HighestRuntime-LowestRuntime+1;
		
	Stage() : n(Stage::Invalid) { }
	Stage(Stage::Num nn) : n(nn) { assert(isAllowable(n)); }	// note implicit conversion	
	// default copy, destructor, assignment
	
	// implicit conversion from Stage to Stage::Num
	operator Stage::Num() const { return Stage::Num(n); }

	// Only prefix increment/decrement supported. Note that we have to
	// allow increment/decrement to go out of range LowestRuntime:HighestRuntime by 1.	
    Stage& operator++() { assert(n==Stage::Num(Stage::LowestRuntime-1) || isInRuntimeRange(n)); ++n; return *this; }
	Stage& operator--() { assert(isInRuntimeRange(n)); --n; return *this; } 

	SimTK_SimTKCOMMON_EXPORT String	name() const;
	Stage	next() const { return isValid(n+1) ? Stage::Num(n+1) : Stage::Invalid; }
    Stage	prev() const { return isValid(n-1) ? Stage::Num(n-1) : Stage::Invalid; }

    // Set stage=min(stage, tooHigh-1).
    void invalidate(Stage::Num tooHigh) {
        assert(tooHigh > Stage::Empty); // don't allow this to make stage Invalid
        if (n >= tooHigh)
            n=tooHigh-1;
    }
	
	// Use this for arrays of NRuntime elements indexed from 0.
	int index() const { assert(isInRuntimeRange(n)); return n-Stage::LowestRuntime; }
   
    static bool isInRuntimeRange(int n) {return Stage::LowestRuntime <= n && n <= Stage::HighestRuntime;}
    static bool isValid         (int n) {return Stage::LowestValid   <= n && n <= Stage::HighestValid;}
    static bool isAllowable     (int n) {return isValid(n) || n==Stage::Invalid;}

private:
	int n;

};
SimTK_LIST_SPECIALIZE(Stage);

namespace Exception {
class StageTooLow : public Base {
public:
    StageTooLow(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be at least " + targetStage.name() + " in " + String(where)
           + " but current stage was " + currentStage.name());
    }
};

class StageTooHigh : public Base {
public:
    StageTooHigh(const char* fn, int ln,
        Stage currentStage, Stage targetStage, const char* where) : Base(fn,ln)
    {
        setMessage("Expected stage to be less than " + targetStage.name() + " in " + String(where)
           + " but current stage was " + currentStage.name());
    }
};

class StageOutOfRange : public Base {
public:
    StageOutOfRange(const char* fn, int ln,
        Stage lower, Stage currentStage, Stage upper, const char* where) : Base(fn,ln)
    {
        setMessage("Expected (" + lower.name() + " <= stage <= " + upper.name() + ") in " + String(where)
           + " but stage was " + currentStage.name());
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
                   + g.name() + ": " + String(buf) + ".");
        va_end(args);
    }
private:
};
}

inline std::ostream& operator<<(std::ostream& o, Stage g) { o << g.name(); return o; }	


} // namespace SimTK

    // REALIZECHECKs: these should be used to catch and report problems that
    // occur when realizing a subsystem.

#define SimTK_REALIZECHECK_ALWAYS(cond,stage,subsysId,subsysName,msg)       \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysId),(subsysName),(msg));                 \
    }while(false)
#define SimTK_REALIZECHECK1_ALWAYS(cond,stage,subsysId,subsysName,msg,a1)   \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysId),(subsysName),(msg),(a1));            \
    }while(false)
#define SimTK_REALIZECHECK2_ALWAYS(cond,stage,subsysId,subsysName,msg,a1,a2)\
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::RealizeCheckFailed,        \
                    (stage),(subsysId),(subsysName),(msg),(a1),(a2));       \
    }while(false)
#define SimTK_REALIZECHECK3_ALWAYS(cond,stage,subsysId,subsysName,msg,a1,a2,a3)     \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::RealizeCheckFailed,                \
                    (stage),(subsysId),(subsysName),(msg),(a1),(a2),(a3));          \
    }while(false)
#define SimTK_REALIZECHECK4_ALWAYS(cond,stage,subsysId,subsysName,msg,a1,a2,a3,a4)  \
    do{if(!(cond))SimTK_THROW8(SimTK::Exception::RealizeCheckFailed,                \
                    (stage),(subsysId),(subsysName),(msg),(a1),(a2),(a3),(a4));     \
    }while(false)
#define SimTK_REALIZECHECK5_ALWAYS(cond,stage,subsysId,subsysName,msg,a1,a2,a3,a4,a5)   \
    do{if(!(cond))SimTK_THROW9(SimTK::Exception::RealizeCheckFailed,                    \
                    (stage),(subsysId),(subsysName),(msg),(a1),(a2),(a3),(a4),(a5));    \
    }while(false)

    
#endif // SimTK_SimTKCOMMON_STAGE_H_
