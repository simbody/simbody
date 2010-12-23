/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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
 * Defines the standard SimTK core "version" and "about" routines.
 */


#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <string>
#include <cstring>
#include <cctype>

#define STR(var) #var
#define MAKE_VERSION_STRING(maj,min,build)  STR(maj.min.build)
#define MAKE_COPYRIGHT_STRING(y,a) \
    "Copyright (c) " STR(y) " Stanford University, " STR(a)
#define MAKE_STRING(a) STR(a)

#define GET_VERSION_STRING  \
    MAKE_VERSION_STRING(SimTK_SIMBODY_MAJOR_VERSION,  \
                        SimTK_SIMBODY_MINOR_VERSION,  \
                        SimTK_SIMBODY_PATCH_VERSION)

#define GET_COPYRIGHT_STRING \
    MAKE_COPYRIGHT_STRING(SimTK_SIMBODY_COPYRIGHT_YEARS, \
                          SimTK_SIMBODY_AUTHORS)

#define GET_SVN_REVISION_STRING \
    MAKE_STRING(SimTK_SIMBODY_SVN_REVISION)

#define GET_AUTHORS_STRING \
    MAKE_STRING(SimTK_SIMBODY_AUTHORS)

#define GET_LIBRARY_STRING \
    MAKE_STRING(SimTK_SIMBODY_LIBRARY_NAME)

#if defined(SimTK_SIMBODY_BUILDING_SHARED_LIBRARY)
    #define GET_TYPE_STRING "shared"
#elif defined(SimTK_SIMBODY_BUILDING_STATIC_LIBRARY)
    #define GET_TYPE_STRING "static"
#else
    #define GET_TYPE_STRING "<unknown library type?!>"
#endif

#ifndef NDEBUG
    #define GET_DEBUG_STRING "debug"
#else
    #define GET_DEBUG_STRING "release"
#endif

extern "C" {

void SimTK_version_simbody(int* major, int* minor, int* patch) {
    static const char* l = "SimTK library="   GET_LIBRARY_STRING;
    static const char* t = "SimTK type="      GET_TYPE_STRING;
    static const char* d = "SimTK debug="     GET_DEBUG_STRING;
    static const char* v = "SimTK version="   GET_VERSION_STRING;
    static const char* r = "SimTK svn_revision=" GET_SVN_REVISION_STRING;
    static const char* c = "SimTK copyright=" GET_COPYRIGHT_STRING;

    if (major) *major = SimTK_SIMBODY_MAJOR_VERSION;
    if (minor) *minor = SimTK_SIMBODY_MINOR_VERSION;
    if (patch) *patch = SimTK_SIMBODY_PATCH_VERSION;

    // Force statics to be present in the binary (Release mode otherwise 
    // optimizes them away).
    volatile int i=0;
    if (i) { // never true, but compiler doesn't know ...
        *major = *l + *t + *d + *v + *r + *c;
    }
}

void SimTK_about_simbody(const char* key, int maxlen, char* value) {
    if (maxlen <= 0 || value==0) return;
    value[0] = '\0'; // in case we don't find a match
    if (key==0) return;

    // downshift the key
    std::string skey(key);
    for (size_t i=0; i<skey.size(); ++i)
        skey[i] = std::tolower(skey[i]);

    const char* v = 0;
    if      (skey == "version")   v = GET_VERSION_STRING;
    else if (skey == "library")   v = GET_LIBRARY_STRING;
    else if (skey == "type")      v = GET_TYPE_STRING;
    else if (skey == "copyright") v = GET_COPYRIGHT_STRING;
    else if (skey == "svn_revision") v = GET_SVN_REVISION_STRING;
    else if (skey == "authors")   v = GET_AUTHORS_STRING;
    else if (skey == "debug")     v = GET_DEBUG_STRING;

    if (v) {
        std::strncpy(value,v,maxlen-1);
        value[maxlen-1] = '\0'; // in case we ran out of room
    }
}

}
