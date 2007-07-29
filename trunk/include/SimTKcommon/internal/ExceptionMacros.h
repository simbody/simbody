#ifndef SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_
#define SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_

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

/** @file
 * This file contains macros which are convenient to use for 
 * sprinkling error checking around liberally in SimTK programs, a
 * practice which is highly encouraged. You can think of this as
 * a generalization of the standard assert() macro. By default, 
 * these macros evaporate completely in a release build, but are
 * present in any debug build. Macros are also provided which are
 * always present in cases where the error checking is not a 
 * performance problem, and those should be used in preference
 * to the disappearing ones when appropriate. Also, you can force
 * the disappearing macros to remain present on a file-by-file basis,
 * primarily for use in debugging those annoying problems which only
 * occur in release builds and won't reproduce in a debug build.
 *
 * Most macros have a similar structure, something like this:
 *    SimTK_MACRONAME[nargs][_ALWAYS](printfString [, args...])
 * for example
 *    SimTK_ASSERT3("expected %d < count < %d but count=%d",
 *                  lower, upper, count);
 * or
 *    SimTK_ASSERT3_ALWAYS("expected %d < count < %d but count=%d",
 *                  lower, upper, count);
 * In addition to disappearing without a trace, these macros will
 * also capture the current file name and line number for reporting
 * to developers when appropriate.
 *
 * Note that these are *global* symbols, so we use the reserved
 * SimTK_ name prefix (since we can't use the SimTK:: namespace 
 * for macros) to attempt to avoid pollution of user programs.
 *
 * We distinguish between macros which are used as internal 
 * "bugcatchers" and those which are used to report errors to
 * API users. The C++ exception mechanism is used in both circumstances
 * but the meaning and intended audience is quite different. Any
 * macro with 'ASSERT' in the name represents an internal error
 * which cannot be attributed to user misbehavior. Otherwise these
 * are for communicating with users. Those need to be carefully documented
 * so that users can selectively catch the exceptions when appropriate.
 */

// Here are the symbols you can define to force various checks to remain
// while in release mode.

// SimTK_KEEP_ASSERT
// SimTK_KEEP_RANGECHECK
// SimTK_KEEP_STAGECHECK
// SimTK_KEEP_APIARGCHECK
// ...

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Exception.h"

#include <string>
#include <iostream>
#include <exception>

// TODO: SHAPECHECK, DOMAINCHECK

    // RANGECHECKs: these exception are to be used for situations in which a
    // user of an API screws up by providing bad indices or dimensions. 
    //   INDEXCHECK: Note that we allow the index to be equal to the lower
    //     bound but it must be strictly less than the upper bound.
    //   SIZECHECK: A size or size expression must be >= 0 and less than OR EQUAL
    //     to the maximum size.
    //   SIZECHECK_NONNEG: An size argument must be >= 0.

// This is a rangecheck that is always present, even in Release mode.
#define SimTK_INDEXCHECK_ALWAYS(lb,ix,ub,where) \
    do{if(!((lb)<=(ix)&&(ix)<(ub)))SimTK_THROW5(SimTK::Exception::IndexOutOfRange,   \
                    #ix,(lb),(ix),(ub),(where));}while(false)

// This is a rangecheck that is always present, even in Release mode.
#define SimTK_SIZECHECK_ALWAYS(sz,maxsz,where) \
    do{if(!(0<=(sz)&&(sz)<=(maxsz)))SimTK_THROW4(SimTK::Exception::SizeOutOfRange,   \
                    #sz,(sz),(maxsz),(where));}while(false)

// This is a rangecheck that is always present, even in Release mode.
#define SimTK_SIZECHECK_NONNEG_ALWAYS(sz,where) \
    do{if((sz)<0)SimTK_THROW3(SimTK::Exception::SizeWasNegative,   \
                    #sz,(sz),(where));}while(false)

    // Similar checks for floating point values.

#define SimTK_VALUECHECK_ALWAYS(lb,val,ub,valName,where) \
    do{if(!(lb)<=(val)&&(val)<=(ub)))SimTK_THROW5(SimTK::Exception::ValueOutOfRange,   \
                    (valName),(lb),(val),(ub),(where));}while(false)


#define SimTK_VALUECHECK_NONNEG_ALWAYS(val,valName,where) \
    do{if((val)<0)SimTK_THROW3(SimTK::Exception::ValueWasNegative,   \
                    (valName),(val),(where));}while(false)


    // APIARGCHECKs: these should be used to catch all manner of problems with
    // a user's call to a method that is part of a SimTK API. Note that these
    // are intended for direct consumption by an application programmer using
    // a SimTK API, so should be wordy and helpful. These macros accept
    // printf-style format strings and arguments of whatever are the 
    // appropriate types for those formats.

#define SimTK_APIARGCHECK_ALWAYS(cond,className,methodName,msg)     \
    do{if(!(cond))SimTK_THROW3(SimTK::Exception::APIArgcheckFailed, \
                    (className),(methodName),(msg));                \
    }while(false)
#define SimTK_APIARGCHECK1_ALWAYS(cond,className,methodName,fmt,a1) \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::APIArgcheckFailed, \
                    (className),(methodName),(fmt),(a1));           \
    }while(false)
#define SimTK_APIARGCHECK2_ALWAYS(cond,className,methodName,fmt,a1,a2)      \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::APIArgcheckFailed,         \
                    (className),(methodName),(fmt),(a1),(a2));              \
    }while(false)
#define SimTK_APIARGCHECK3_ALWAYS(cond,className,methodName,fmt,a1,a2,a3)   \
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::APIArgcheckFailed,         \
                    (className),(methodName),(fmt),(a1),(a2),(a3));         \
    }while(false)
#define SimTK_APIARGCHECK4_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4)    \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::APIArgcheckFailed,             \
                    (className),(methodName),(fmt),(a1),(a2),(a3),(a4));        \
    }while(false)
#define SimTK_APIARGCHECK5_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4,a5) \
    do{if(!(cond))SimTK_THROW8(SimTK::Exception::APIArgcheckFailed,             \
                    (className),(methodName),(fmt),(a1),(a2),(a3),(a4),(a5));   \
    }while(false)


#if defined(NDEBUG) && !defined(SimTK_KEEP_RANGECHECK)
    #define SimTK_INDEXCHECK(lb,ix,ub,where)
    #define SimTK_SIZECHECK(sz,maxsz,where)
    #define SimTK_SIZECHECK_NONNEG(sz,where)
    #define SimTK_VALUECHECK(lb,val,ub,valName,where)
    #define SimTK_VALUECHECK_NONNEG(val,valName,where)
#else
    #define SimTK_INDEXCHECK(lb,ix,ub,where) SimTK_INDEXCHECK_ALWAYS(lb,ix,ub,where)
    #define SimTK_SIZECHECK(sz,maxsz,where)  SimTK_SIZECHECK_ALWAYS(sz,maxsz,where)
    #define SimTK_SIZECHECK_NONNEG(sz,where) SimTK_SIZECHECK_NONNEG_ALWAYS(sz,where)
    #define SimTK_VALUECHECK(lb,val,ub,valName,where)  SimTK_VALUECHECK_ALWAYS(lb,val,ub,valName,where)
    #define SimTK_VALUECHECK_NONNEG(val,valName,where) SimTK_VALUECHECK_NONNEG_ALWAYS(val,valName,where)
#endif

#if defined(NDEBUG) && !defined(SimTK_KEEP_APIARGCHECK)
    #define SimTK_APIARGCHECK(cond,className,methodName,msg)
    #define SimTK_APIARGCHECK1(cond,className,methodName,fmt,a1)
    #define SimTK_APIARGCHECK2(cond,className,methodName,fmt,a1,a2)
    #define SimTK_APIARGCHECK3(cond,className,methodName,fmt,a1,a2,a3)
    #define SimTK_APIARGCHECK4(cond,className,methodName,fmt,a1,a2,a3,a4)
    #define SimTK_APIARGCHECK5(cond,className,methodName,fmt,a1,a2,a3,a4,a5)
#else
    #define SimTK_APIARGCHECK(cond,className,methodName,msg)                       \
        SimTK_APIARGCHECK_ALWAYS(cond,className,methodName,msg)
    #define SimTK_APIARGCHECK1(cond,className,methodName,fmt,a1)                   \
        SimTK_APIARGCHECK1_ALWAYS(cond,className,methodName,fmt,a1)
    #define SimTK_APIARGCHECK2(cond,className,methodName,fmt,a1,a2)                \
        SimTK_APIARGCHECK2_ALWAYS(cond,className,methodName,fmt,a1,a2)
    #define SimTK_APIARGCHECK3(cond,className,methodName,fmt,a1,a2,a3)             \
        SimTK_APIARGCHECK3_ALWAYS(cond,className,methodName,fmt,a1,a2,a3)
    #define SimTK_APIARGCHECK4(cond,className,methodName,fmt,a1,a2,a3,a4)          \
        SimTK_APIARGCHECK4_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4)
    #define SimTK_APIARGCHECK5(cond,className,methodName,fmt,a1,a2,a3,a4,a5)       \
        SimTK_APIARGCHECK5_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4,a5)
#endif

    // STAGECHECKs: these exception are to be used for situations in which a
    // user of an API screws up by attempting to access something in the 
    // state before it has been realized to the appropriate stage.
    //
    //   STAGECHECK_EQ: Check that the current stage is == a particular stage.
    //   STAGECHECK_GE: Check that the current stage is >= a particular stage.
    //   STAGECHECK_LT: Check that the current stage is <  a particular stage.
    //   STAGECHECK_RANGE: Check that lower <= stage <= upper.

// These are stagechecks that is always present, even in Release mode.
#define SimTK_STAGECHECK_EQ_ALWAYS(currentStage,targetStage,where) \
    do{if((currentStage)!=(targetStage)) SimTK_THROW3(SimTK::Exception::StageIsWrong,   \
        (currentStage),(targetStage),(where));}while(false)
#define SimTK_STAGECHECK_GE_ALWAYS(currentStage,targetStage,where) \
    do{if(!((currentStage)>=(targetStage))) SimTK_THROW3(SimTK::Exception::StageTooLow,   \
        (currentStage),(targetStage),(where));}while(false)
#define SimTK_STAGECHECK_LT_ALWAYS(currentStage,targetStage,where) \
    do{if((currentStage)>=(targetStage)) SimTK_THROW3(SimTK::Exception::StageTooHigh,   \
        (currentStage),(targetStage),(where));}while(false)
#define SimTK_STAGECHECK_RANGE_ALWAYS(lower,current,upper,where) \
    do{if(!((lower)<=(current)&&(current)<=(upper))) SimTK_THROW4(SimTK::Exception::StageOutOfRange,   \
        (lower),(current),(upper),(where));}while(false)

// This one is present only in Debug mode or if SimTK_KEEP_STAGECHECK is explicitly defined.
#if defined(NDEBUG) && !defined(SimTK_KEEP_STAGECHECK)
    #define SimTK_STAGECHECK_EQ(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_GE(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_LE(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_RANGE(lower,current,upper,where)
#else
    #define SimTK_STAGECHECK_EQ(currentStage,targetStage,where) \
        SimTK_STAGECHECK_EQ_ALWAYS(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_GE(currentStage,targetStage,where) \
        SimTK_STAGECHECK_GE_ALWAYS(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_LE(currentStage,targetStage,where) \
        SimTK_STAGECHECK_LE_ALWAYS(currentStage,targetStage,where)
    #define SimTK_STAGECHECK_RANGE(lower,current,upper,where) \
        SimTK_STAGECHECK_RANGE_ALWAYS(lower,current,upper,where)
#endif



    // ASSERT: use this *only* for internal errors, that is, bugs. This must
    // not be used to catch usage errors by clients; if you want to catch
    // user errors use different exceptions.

// This is an assertion that is always active, even in Release mode.
#define SimTK_ASSERT_ALWAYS(cond,msg) \
    do{if(!(cond))SimTK_THROW2(SimTK::Exception::Assert,#cond,(msg));}while(false)
#define SimTK_ASSERT1_ALWAYS(cond,msg,a1) \
    do{if(!(cond))SimTK_THROW3(SimTK::Exception::Assert,#cond,(msg),(a1));}while(false)
#define SimTK_ASSERT2_ALWAYS(cond,msg,a1,a2) \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::Assert,#cond,(msg),(a1),(a2));}while(false)
#define SimTK_ASSERT3_ALWAYS(cond,msg,a1,a2,a3) \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::Assert,#cond,(msg),(a1),(a2),(a3));}while(false)
#define SimTK_ASSERT4_ALWAYS(cond,msg,a1,a2,a3,a4) \
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::Assert,#cond,(msg),(a1),(a2),(a3),(a4));}while(false)

// Note: unlike the system assert() we're putting ours within the header guards.
// So if you want to override NDEBUG do it at the *beginning* (that is, before
// the first #include or #ifdef) of whatever compilation unit you are fiddling with.
#if defined(NDEBUG) && !defined(SimTK_KEEP_ASSERT)
    #define SimTK_ASSERT(cond,msg)
    #define SimTK_ASSERT(cond,msg)
    #define SimTK_ASSERT1(cond,msg,a1)
    #define SimTK_ASSERT2(cond,msg,a1,a2)
    #define SimTK_ASSERT3(cond,msg,a1,a2,a3)
    #define SimTK_ASSERT4(cond,msg,a1,a2,a3,a4)
#else
    #define SimTK_ASSERT(cond,msg) SimTK_ASSERT_ALWAYS(cond,msg)
    #define SimTK_ASSERT1(cond,msg,a1) SimTK_ASSERT1_ALWAYS(cond,msg,a1)
    #define SimTK_ASSERT2(cond,msg,a1,a2) SimTK_ASSERT2_ALWAYS(cond,msg,a1,a2)
    #define SimTK_ASSERT3(cond,msg,a1,a2,a3) SimTK_ASSERT3_ALWAYS(cond,msg,a1,a2,a3)
    #define SimTK_ASSERT4(cond,msg,a1,a2,a3,a4) SimTK_ASSERT4_ALWAYS(cond,msg,a1,a2,a3,a4)
#endif


#endif // SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_



