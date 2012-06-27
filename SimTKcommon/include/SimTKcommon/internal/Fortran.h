#ifndef SimTK_SimTKCOMMON_FORTRAN_H_
#define SimTK_SimTKCOMMON_FORTRAN_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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

/**@file
 * This header defines a set of macros which are useful for 
 * making calls to Fortran routines from C++ programs. Please don't
 * take this as encouragement to use Fortran, but sometimes it
 * is necessary.
 */

// Although we are currently triggering this off the OS, these
// really are compiler dependencies. 
//      (1) The calling convention (__stdcall for Windows)
//      (2) Name capitalization (either all-lower or all-uppercase)
//      (3) Is a trailing underscore added to the name?
//      (4) And ugliest, Fortran passes string lengths as a hidden
//          value parameter. On some compilers, that length follows
//          the string immediately. On others, all the lengths
//          appear at the end of the argument list, in the same
//          order as the strings to which they correspond.
// Point (4) requires four ugly macros to be used, two in declarations
// and two in calls to Fortran routines. One macro appears immediately
// after each string, and the other appears at the end of the argument
// list, repeated as many times as necessary. One or the other will
// evaluate to nothing.

// These macros should be used for whatever the expected default
// Fortran behavior is for whatever Fortran is typically used in
// conjunction with the current C++ compiler.
#ifdef _WIN32
    #define SimTK_FORTRAN_STDCALL __stdcall
    #define SimTK_FORTRAN(x,X) X
    #define SimTK_FORTRAN_STRLEN_FOLLOWS_DECL       ,int
    #define SimTK_FORTRAN_STRLEN_FOLLOWS_CALL(n)    ,n
    #define SimTK_FORTRAN_STRLEN_ATEND_DECL         // nothing
    #define SimTK_FORTRAN_STRLEN_ATEND_CALL(n)
#else
    #define SimTK_FORTRAN_STDCALL
    #define SimTK_FORTRAN(x,X) x ## _
    #define SimTK_FORTRAN_STRLEN_FOLLOWS_DECL       // nothing
    #define SimTK_FORTRAN_STRLEN_FOLLOWS_CALL(n)
    #define SimTK_FORTRAN_STRLEN_ATEND_DECL         ,int
    #define SimTK_FORTRAN_STRLEN_ATEND_CALL(n)      ,n
#endif

// These macros should be used for whatever our chosen LAPACK and
// BLAS libraries will look like from here.
#ifdef SimTK_USE_ACML_LAPACK
  #ifdef _WIN32
    #define SimTK_LAPACK_STDCALL __stdcall
    #define SimTK_LAPACK(x,X) X
    #define SimTK_LAPACK_STRLEN_FOLLOWS_DECL       ,int
    #define SimTK_LAPACK_STRLEN_FOLLOWS_CALL(n)    ,n
    #define SimTK_LAPACK_STRLEN_ATEND_DECL         // nothing
    #define SimTK_LAPACK_STRLEN_ATEND_CALL(n)
  #else
    #define SimTK_LAPACK_STDCALL
    #define SimTK_LAPACK(x,X) x ## _
    #define SimTK_LAPACK_STRLEN_FOLLOWS_DECL       // nothing
    #define SimTK_LAPACK_STRLEN_FOLLOWS_CALL(n)
    #define SimTK_LAPACK_STRLEN_ATEND_DECL         ,int
    #define SimTK_LAPACK_STRLEN_ATEND_CALL(n)      ,n
  #endif
#else // default assumes we're using libSimTKlapack
    #define SimTK_LAPACK_STDCALL
    #define SimTK_LAPACK(x,X) x ## _
    #define SimTK_LAPACK_STRLEN_FOLLOWS_DECL       // nothing
    #define SimTK_LAPACK_STRLEN_FOLLOWS_CALL(n)
    #define SimTK_LAPACK_STRLEN_ATEND_DECL         ,int
    #define SimTK_LAPACK_STRLEN_ATEND_CALL(n)      ,n
#endif

// TODO: Currently this is unused and may not be needed anymore.
// Call these routines to intialize the GNU Fortran RTL.
// 

#ifdef USING_G77
    extern "C" {
        void f_setsig();
        void f_init();
    }
    #define SimTK_FORTRAN_INIT do {f_setsig(); f_init();} while(false)
#endif

#endif // SimTK_SimTKCOMMON_FORTRAN_H_
