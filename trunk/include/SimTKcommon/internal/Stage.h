#ifndef _SimTK_STAGE_H_
#define _SimTK_STAGE_H_

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
    //  Initialized            Allocate    Allocated
    //  Topology               Build       Built
    //  Model                  Model       Modeled
    //  Specifics              Specify     Specified
    //  Time                   Time        Timed
    //  Position               Position    Positioned
    //  Velocity               Move        Moving
    //  Dynamics               ?
    //  Acceleration           Accelerate  Accelerated
    //  Report                 Report      Reported

	enum Num {
        Invalid        = -99,
        Allocated      =  0, // TODO: Initialized, Unbuilt, Empty?
        Built          =  1, // TODO: Constructed, Finalized?
        Modeled        =  2, // TODO: Instantiated, Resourced, Provisioned, Specialized?
        Parametrized   =  3,
        Timed          =  4,
        Configured     =  5, // TODO: Positioned?
        Moving         =  6, // TODO: Velocity, Speed, Rate?
        Dynamics       =  7, // dynamic properties & operators available
        Reacting       =  8  // TODO: Accelerated?
	};

    static const Stage::Num LowestValid     = Allocated;
    static const Stage::Num HighestValid    = Reacting;
    static const int        NValid          = HighestValid-LowestValid+1;

    // LowestRuntime->HighestRuntime cover the post-construction stages only.
	static const Stage::Num	LowestRuntime	= Modeled;
	static const Stage::Num	HighestRuntime	= Reacting;
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

	SimTK_SimTKCOMMON_API String	name() const;
	Stage	next() const { return isValid(n+1) ? Stage::Num(n+1) : Stage::Invalid; }
    Stage	prev() const { return isValid(n-1) ? Stage::Num(n-1) : Stage::Invalid; }
	
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
}

inline std::ostream& operator<<(std::ostream& o, Stage g) { o << g.name(); return o; }	


} // namespace SimTK

    
#endif //_SimTK_STAGE_H_
