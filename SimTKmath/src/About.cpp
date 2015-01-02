/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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
 * Defines the standard SimTK core "version" and "about" routines.
 */


#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <string>
#include <cstring>
#include <cctype>

#define STR(var) #var
#define MAKE_VERSION_STRING(maj,min,build)  STR(maj.min.build)
#define MAKE_COPYRIGHT_STRING(y,a) \
    "Copyright (c) " STR(y) " Stanford University, " STR(a)
#define MAKE_STRING(a) STR(a)

#define GET_VERSION_STRING  \
    MAKE_VERSION_STRING(SimTK_SIMMATH_MAJOR_VERSION,  \
                        SimTK_SIMMATH_MINOR_VERSION,  \
                        SimTK_SIMMATH_PATCH_VERSION)

#define GET_COPYRIGHT_STRING \
    MAKE_COPYRIGHT_STRING(SimTK_SIMMATH_COPYRIGHT_YEARS, \
                          SimTK_SIMMATH_AUTHORS)

#define GET_AUTHORS_STRING \
    MAKE_STRING(SimTK_SIMMATH_AUTHORS)

#define GET_LIBRARY_STRING \
    MAKE_STRING(SimTK_SIMMATH_LIBRARY_NAME)

#ifndef NDEBUG
    #define GET_DEBUG_STRING "debug"
#else
    #define GET_DEBUG_STRING "release"
#endif

#if defined(SimTK_SIMMATH_BUILDING_SHARED_LIBRARY)
    #define GET_TYPE_STRING "shared"
#elif defined(SimTK_SIMMATH_BUILDING_STATIC_LIBRARY)
    #define GET_TYPE_STRING "static"
#else
    #define GET_TYPE_STRING "<unknown library type?!>"
#endif

extern "C" {

void SimTK_version_simmath(int* major, int* minor, int* patch) {
    static const char* l = "SimTK library="   GET_LIBRARY_STRING;
    static const char* t = "SimTK type="      GET_TYPE_STRING;
    static const char* d = "SimTK debug="     GET_DEBUG_STRING;
    static const char* v = "SimTK version="   GET_VERSION_STRING;
    static const char* c = "SimTK copyright=" GET_COPYRIGHT_STRING;

    if (major) *major = SimTK_SIMMATH_MAJOR_VERSION;
    if (minor) *minor = SimTK_SIMMATH_MINOR_VERSION;
    if (patch) *patch = SimTK_SIMMATH_PATCH_VERSION;

    // Force statics to be present in the binary (Release mode otherwise 
    // optimizes them away).
    volatile int i=0;
    if (i) { // never true, but compiler doesn't know ...
        *major = *l + *t + *d + *v + *c;
    }
}

void SimTK_about_simmath(const char* key, int maxlen, char* value) {
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
    else if (skey == "authors")   v = GET_AUTHORS_STRING;
    else if (skey == "debug")     v = GET_DEBUG_STRING;

    if (v) {
        std::strncpy(value,v,maxlen-1);
        value[maxlen-1] = '\0'; // in case we ran out of room
    }
}

}
