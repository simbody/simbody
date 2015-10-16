#ifndef SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_
#define SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_

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
 * <pre>
 *    SimTK_MACRONAME[nargs][_ALWAYS](cond, printfString [, args...])
 * for example
 *    SimTK_ASSERT3(lower<count && count<upper,
 *                  "expected %d < count < %d but count=%d",
 *                  lower, upper, count);
 * or
 *    SimTK_ASSERT3_ALWAYS(lower<count && count<upper,
 *                  "expected %d < count < %d but count=%d",
 *                  lower, upper, count);
 * </pre>
 * To override the disappearance of the non-ALWAYS macros, your 
 * compilation should define a preprocessor symbols like
 * SimTK_KEEP_ASSERT, SimTK_KEEP_ERRCHK, etc.
 *
 * These macros will also capture the current file name and line 
 * number for reporting to developers when appropriate, and those
 * with a condition that failed will include the condition in the 
 * message.
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
 * which cannot be attributed to user misbehavior. Other macros
 * are for communicating with users. Those need to be carefully documented
 * so that users can selectively catch the exceptions when appropriate.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Exception.h"

#include <string>
#include <iostream>
#include <exception>

// --------------------------------- RANGECHECK --------------------------------
// These exceptions are to be used for situations in which a user of a SimTK
// API method screws up by providing bad indices or dimensions. These are special
// cases of the more general APIARGCHECK macros, providing "canned" error messages
// for several common situations. Although there are several different macro 
// names here, all are controlled by SimTK_KEEP_RANGECHECK to allow enabling of
// these index- and size-validating tests together in Release mode.
//
//   INDEXCHECK: Note that we allow the index to be equal to the lower
//     bound (zero) but it must be strictly less than the upper bound.
//   SIZECHECK: A size or size expression must be >= 0 and less than OR EQUAL
//     to the maximum size.
//   SIZECHECK_NONNEG: A size argument must be >= 0.
//   VALUECHECK: A floating point is required to be within a certain range.
//   VALUECHECK_NONNEG: A floating point argument must be non-negative.
//
// TODO: SHAPECHECK, DOMAINCHECK
// -----------------------------------------------------------------------------

// This is a rangecheck that is always present, even in Release mode. This may be
// applied both to signed and unsigned types (the latter are always nonnegative) so
// to avoid warnings we use the isIndexInRange() method which doesn't perform
// a nonnegativity check on unsigned quantities.
#define SimTK_INDEXCHECK_ALWAYS(ix,ub,where) \
    do{if(!isIndexInRange((ix),(ub)))SimTK_THROW5(SimTK::Exception::IndexOutOfRange,   \
                    #ix,0,(ix),(ub),(where));}while(false)

// This is a rangecheck that is always present, even in Release mode. This may be
// applied both to signed and unsigned types (the latter are always nonnegative) so
// to avoid warnings we use the isSizeInRange() method which doesn't perform
// a nonnegativity check on unsigned quantities.
#define SimTK_SIZECHECK_ALWAYS(sz,maxsz,where) \
    do{if(!isSizeInRange((sz),(maxsz)))SimTK_THROW4(SimTK::Exception::SizeOutOfRange,   \
                    #sz,(sz),(maxsz),(where));}while(false)

// This is a rangecheck that is always present, even in Release mode. Use
// isNonnegative() here in case sz is an unsigned type to avoid compiler
// warning.
#define SimTK_SIZECHECK_NONNEG_ALWAYS(sz,where) \
    do{if(!isNonnegative(sz))SimTK_THROW3(SimTK::Exception::SizeWasNegative,   \
                    #sz,(sz),(where));}while(false)

    // Similar checks for floating point values.

#define SimTK_VALUECHECK_ALWAYS(lb,val,ub,valName,where) \
    do{if(!((lb)<=(val)&&(val)<=(ub)))SimTK_THROW5(SimTK::Exception::ValueOutOfRange,   \
                    (valName),(lb),(val),(ub),(where));}while(false)


#define SimTK_VALUECHECK_NONNEG_ALWAYS(val,valName,where) \
    do{if((val)<0)SimTK_THROW3(SimTK::Exception::ValueWasNegative,   \
                    (valName),(val),(where));}while(false)



#if defined(NDEBUG) && !defined(SimTK_KEEP_RANGECHECK)
    #define SimTK_INDEXCHECK(ix,ub,where)
    #define SimTK_SIZECHECK(sz,maxsz,where)
    #define SimTK_SIZECHECK_NONNEG(sz,where)
    #define SimTK_VALUECHECK(lb,val,ub,valName,where)
    #define SimTK_VALUECHECK_NONNEG(val,valName,where)
#else
    #define SimTK_INDEXCHECK(ix,ub,where) SimTK_INDEXCHECK_ALWAYS(ix,ub,where)
    #define SimTK_SIZECHECK(sz,maxsz,where)  SimTK_SIZECHECK_ALWAYS(sz,maxsz,where)
    #define SimTK_SIZECHECK_NONNEG(sz,where) SimTK_SIZECHECK_NONNEG_ALWAYS(sz,where)
    #define SimTK_VALUECHECK(lb,val,ub,valName,where)  SimTK_VALUECHECK_ALWAYS(lb,val,ub,valName,where)
    #define SimTK_VALUECHECK_NONNEG(val,valName,where) SimTK_VALUECHECK_NONNEG_ALWAYS(val,valName,where)
#endif


// --------------------------------- STAGECHECK --------------------------------
// These exceptions are to be used for situations in which a
// user of an API screws up by attempting to access something in the 
// state before it has been realized to the appropriate stage.
//
//   STAGECHECK_TOPOLOGY_REALIZED: Check that realizeTopology() has been done
//     since the last topological change.
//   STAGECHECK_EQ: Check that the current stage is == a particular stage.
//   STAGECHECK_GE: Check that the current stage is >= a particular stage.
//   STAGECHECK_LT: Check that the current stage is <  a particular stage.
//   STAGECHECK_RANGE: Check that lower <= stage <= upper.
// -----------------------------------------------------------------------------

// These are stagechecks that is always present, even in Release mode.
#define SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(cond,objType,objName,methodNm)  \
    do{if(!(cond)) SimTK_THROW3(SimTK::Exception::RealizeTopologyMustBeCalledFirst, \
        (objType),(objName),(methodNm));}while(false)
#define SimTK_STAGECHECK_TOPOLOGY_VERSION_ALWAYS(sysTopoVersion,            \
                                stateTopoVersion,objType,objName,methodNm)  \
    do{if((stateTopoVersion)!=(sysTopoVersion))                             \
        SimTK_THROW5(SimTK::Exception::StateAndSystemTopologyVersionsMustMatch, \
        (objType),(objName),(methodNm),                     \
        (int)(sysTopoVersion),(int)(stateTopoVersion));}    \
    while(false)
#define SimTK_STAGECHECK_EQ_ALWAYS(currentStage,targetStage,methodNm) \
    do{if((currentStage)!=(targetStage)) SimTK_THROW3(SimTK::Exception::StageIsWrong,   \
        (currentStage),(targetStage),(methodNm));}while(false)
#define SimTK_STAGECHECK_GE_ALWAYS(currentStage,targetStage,methodNm) \
    do{if(!((currentStage)>=(targetStage))) SimTK_THROW3(SimTK::Exception::StageTooLow,   \
        (currentStage),(targetStage),(methodNm));}while(false)
#define SimTK_STAGECHECK_LT_ALWAYS(currentStage,targetStage,methodNm) \
    do{if((currentStage)>=(targetStage)) SimTK_THROW3(SimTK::Exception::StageTooHigh,   \
        (currentStage),(targetStage),(methodNm));}while(false)
#define SimTK_STAGECHECK_RANGE_ALWAYS(lower,current,upper,methodNm) \
    do{if(!((lower)<=(current)&&(current)<=(upper))) SimTK_THROW4(SimTK::Exception::StageOutOfRange,   \
        (lower),(current),(upper),(methodNm));}while(false)

// This one is present only in Debug mode or if SimTK_KEEP_STAGECHECK is explicitly defined.
#if defined(NDEBUG) && !defined(SimTK_KEEP_STAGECHECK)
    #define SimTK_STAGECHECK_TOPOLOGY_REALIZED(cond,objType,objName,methodName)
    #define SimTK_STAGECHECK_TOPOLOGY_VERSIONS(sysTopoVersion,stateTopoVersion,\
                                               objType,objName,methodNm)   
    #define SimTK_STAGECHECK_EQ(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_GE(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_LT(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_RANGE(lower,current,upper,methodNm)
#else
    #define SimTK_STAGECHECK_TOPOLOGY_REALIZED(cond,objType,objName,methodName) \
        SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(cond,objType,objName,methodName)
    #define SimTK_STAGECHECK_TOPOLOGY_VERSION(sysTopoVersion,stateTopoVersion, \
                                              objType,objName,methodNm)  \
        SimTK_STAGECHECK_TOPOLOGY_VERSION_ALWAYS(sysTopoVersion,stateTopoVersion,\
                                                 objType,objName,methodNm)
    #define SimTK_STAGECHECK_EQ(currentStage,targetStage,methodNm) \
        SimTK_STAGECHECK_EQ_ALWAYS(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_GE(currentStage,targetStage,methodNm) \
        SimTK_STAGECHECK_GE_ALWAYS(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_LT(currentStage,targetStage,methodNm) \
        SimTK_STAGECHECK_LE_ALWAYS(currentStage,targetStage,methodNm)
    #define SimTK_STAGECHECK_RANGE(lower,current,upper,methodNm) \
        SimTK_STAGECHECK_RANGE_ALWAYS(lower,current,upper,methodNm)
#endif

// -------------------------------- APIARGCHECK --------------------------------
// These should be used to catch all manner of problems with the arguments passed
// in an API user's call to a method that is part of a SimTK API. Note that these
// are intended for direct consumption by an application programmer using a SimTK 
// API, so should be wordy and helpful. These macros accept printf-style format 
// strings and arguments of whatever are the appropriate types for those formats.
// -----------------------------------------------------------------------------

#define SimTK_APIARGCHECK_ALWAYS(cond,className,methodName,msg)     \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::APIArgcheckFailed, \
              #cond,(className),(methodName),(msg));                \
    }while(false)
#define SimTK_APIARGCHECK1_ALWAYS(cond,className,methodName,fmt,a1) \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::APIArgcheckFailed, \
              #cond,(className),(methodName),(fmt),(a1));           \
    }while(false)
#define SimTK_APIARGCHECK2_ALWAYS(cond,className,methodName,fmt,a1,a2)      \
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::APIArgcheckFailed,         \
              #cond,(className),(methodName),(fmt),(a1),(a2));              \
    }while(false)
#define SimTK_APIARGCHECK3_ALWAYS(cond,className,methodName,fmt,a1,a2,a3)   \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::APIArgcheckFailed,         \
              #cond,(className),(methodName),(fmt),(a1),(a2),(a3));         \
    }while(false)
#define SimTK_APIARGCHECK4_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4)    \
    do{if(!(cond))SimTK_THROW8(SimTK::Exception::APIArgcheckFailed,             \
              #cond,(className),(methodName),(fmt),(a1),(a2),(a3),(a4));        \
    }while(false)
#define SimTK_APIARGCHECK5_ALWAYS(cond,className,methodName,fmt,a1,a2,a3,a4,a5) \
    do{if(!(cond))SimTK_THROW9(SimTK::Exception::APIArgcheckFailed,             \
              #cond,(className),(methodName),(fmt),(a1),(a2),(a3),(a4),(a5));   \
    }while(false)

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


// ----------------------------------- ERRCHK ----------------------------------
// ERRCHK: these should be used to catch all manner of problems that occur
// during execution of an API user's request by a method that is part of 
// a SimTK API. Note that these are intended for direct consumption by 
// an application programmer using a SimTK API, so should be wordy and 
// helpful. These macros accept printf-style format strings and arguments 
// of whatever are the appropriate types for those formats.
// -----------------------------------------------------------------------------

#define SimTK_ERRCHK_ALWAYS(cond,whereChecked,msg)              \
    do{if(!(cond))SimTK_THROW3(SimTK::Exception::ErrorCheck,    \
              #cond,(whereChecked),(msg));                      \
    }while(false)
#define SimTK_ERRCHK1_ALWAYS(cond,whereChecked,fmt,a1)          \
    do{if(!(cond))SimTK_THROW4(SimTK::Exception::ErrorCheck,    \
              #cond,(whereChecked),(fmt),(a1));                 \
    }while(false)
#define SimTK_ERRCHK2_ALWAYS(cond,whereChecked,fmt,a1,a2)       \
    do{if(!(cond))SimTK_THROW5(SimTK::Exception::ErrorCheck,    \
              #cond,(whereChecked),(fmt),(a1),(a2));            \
    }while(false)
#define SimTK_ERRCHK3_ALWAYS(cond,whereChecked,fmt,a1,a2,a3)    \
    do{if(!(cond))SimTK_THROW6(SimTK::Exception::ErrorCheck,    \
              #cond,(whereChecked),(fmt),(a1),(a2),(a3));       \
    }while(false)
#define SimTK_ERRCHK4_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4) \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::ErrorCheck,    \
              #cond,(whereChecked),(fmt),(a1),(a2),(a3),(a4));  \
    }while(false)
#define SimTK_ERRCHK5_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5)      \
    do{if(!(cond))SimTK_THROW8(SimTK::Exception::ErrorCheck,            \
              #cond,(whereChecked),(fmt),(a1),(a2),(a3),(a4),(a5));     \
    }while(false)
#define SimTK_ERRCHK6_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6)       \
    do{if(!(cond))SimTK_THROW9(SimTK::Exception::ErrorCheck,                \
              #cond,(whereChecked),(fmt),(a1),(a2),(a3),(a4),(a5),(a6));    \
    }while(false)
#define SimTK_ERRCHK7_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6,a7)        \
    do{if(!(cond))SimTK_THROW10(SimTK::Exception::ErrorCheck,                   \
              #cond,(whereChecked),(fmt),(a1),(a2),(a3),(a4),(a5),(a6),(a7));   \
    }while(false)

#if defined(NDEBUG) && !defined(SimTK_KEEP_ERRCHK)
    #define SimTK_ERRCHK(cond,whereChecked,msg)
    #define SimTK_ERRCHK1(cond,whereChecked,fmt,a1)
    #define SimTK_ERRCHK2(cond,whereChecked,fmt,a1,a2)
    #define SimTK_ERRCHK3(cond,whereChecked,fmt,a1,a2,a3)
    #define SimTK_ERRCHK4(cond,whereChecked,fmt,a1,a2,a3,a4)
    #define SimTK_ERRCHK5(cond,whereChecked,fmt,a1,a2,a3,a4,a5)
    #define SimTK_ERRCHK6(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6)
    #define SimTK_ERRCHK7(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6)
#else
    #define SimTK_ERRCHK(cond,whereChecked,msg)                       \
        SimTK_ERRCHK_ALWAYS(cond,whereChecked,msg)
    #define SimTK_ERRCHK1(cond,whereChecked,fmt,a1)                   \
        SimTK_ERRCHK1_ALWAYS(cond,whereChecked,fmt,a1)
    #define SimTK_ERRCHK2(cond,whereChecked,fmt,a1,a2)                \
        SimTK_ERRCHK2_ALWAYS(cond,whereChecked,fmt,a1,a2)
    #define SimTK_ERRCHK3(cond,whereChecked,fmt,a1,a2,a3)             \
        SimTK_ERRCHK3_ALWAYS(cond,whereChecked,fmt,a1,a2,a3)
    #define SimTK_ERRCHK4(cond,whereChecked,fmt,a1,a2,a3,a4)          \
        SimTK_ERRCHK4_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4)
    #define SimTK_ERRCHK5(cond,whereChecked,fmt,a1,a2,a3,a4,a5)       \
        SimTK_ERRCHK5_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5)
    #define SimTK_ERRCHK6(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6)    \
        SimTK_ERRCHK6_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6)
    #define SimTK_ERRCHK7(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6,a7) \
        SimTK_ERRCHK7_ALWAYS(cond,whereChecked,fmt,a1,a2,a3,a4,a5,a6,a7)
#endif

// ----------------------------------- ASSERT ----------------------------------
// ASSERT: use this *only* for internal errors, that is, bugs. This must
// not be used to catch usage errors by clients; if you want to catch
// user errors use different exceptions.
// -----------------------------------------------------------------------------

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
#define SimTK_ASSERT5_ALWAYS(cond,msg,a1,a2,a3,a4,a5) \
    do{if(!(cond))SimTK_THROW7(SimTK::Exception::Assert,#cond,(msg),(a1),(a2),(a3),(a4),(a5));}while(false)

// Note: unlike the system assert() we're putting ours within the header guards.
// So if you want to override NDEBUG do it at the *beginning* (that is, before
// the first #include or #ifdef) of whatever compilation unit you are fiddling with.
#if defined(NDEBUG) && !defined(SimTK_KEEP_ASSERT)
    #define SimTK_ASSERT(cond,msg)
    #define SimTK_ASSERT1(cond,msg,a1)
    #define SimTK_ASSERT2(cond,msg,a1,a2)
    #define SimTK_ASSERT3(cond,msg,a1,a2,a3)
    #define SimTK_ASSERT4(cond,msg,a1,a2,a3,a4)
    #define SimTK_ASSERT5(cond,msg,a1,a2,a3,a4,a5)
#else
    #define SimTK_ASSERT(cond,msg) SimTK_ASSERT_ALWAYS(cond,msg)
    #define SimTK_ASSERT1(cond,msg,a1) SimTK_ASSERT1_ALWAYS(cond,msg,a1)
    #define SimTK_ASSERT2(cond,msg,a1,a2) SimTK_ASSERT2_ALWAYS(cond,msg,a1,a2)
    #define SimTK_ASSERT3(cond,msg,a1,a2,a3) SimTK_ASSERT3_ALWAYS(cond,msg,a1,a2,a3)
    #define SimTK_ASSERT4(cond,msg,a1,a2,a3,a4) SimTK_ASSERT4_ALWAYS(cond,msg,a1,a2,a3,a4)
    #define SimTK_ASSERT5(cond,msg,a1,a2,a3,a4,a5) SimTK_ASSERT5_ALWAYS(cond,msg,a1,a2,a3,a4,a5)
#endif


#endif // SimTK_SimTKCOMMON_EXCEPTION_MACROS_H_



