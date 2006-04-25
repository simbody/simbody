#ifndef _SimTK_SIMMATRIX_COMMON_H_
#define _SimTK_SIMMATRIX_COMMON_H_

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

/** @file
 * Every Simmatrix header and source file should include this header before
 * any other Simmatrix header.
 */

#include "SimTKcommon.h"

// When building a shared library 'xyz', CMake defines a symbol 'xyz_EXPORTS'
// for use in distinguishing builds from client use of a header. The following
// is specific for the current 'simtk' library and doesn't affect other
// libraries even if they use this one.
#ifndef SimTK_SIMBODY_API
    #ifdef WIN32
        #ifdef simbody_EXPORTS
            #define SimTK_SIMBODY_API __declspec(dllexport)
        #elif defined(SimTK_OPTIMIZE_FOR_DYNAMIC_LIBRARY)
            #define SimTK_SIMBODY_API __declspec(dllimport)   // can't link with static lib now
        #else
            #define SimTK_SIMBODY_API // This works both for static & dynamic clients
        #endif
    #else
        #define SimTK_SIMBODY_API // Linux, Mac
    #endif
#endif

// Every SimTK Core library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_SIMBODY_API void SimTK_version_simmatrix(int* major, int* minor, int* build);
    SimTK_SIMBODY_API void SimTK_about_simmatrix(const char* key, int maxlen, char* value);
}

#endif //_SimTK_SIMMATRIX_COMMON_H_
