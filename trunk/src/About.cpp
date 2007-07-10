/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
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

/** @file
 * Defines the standard SimTK core "version" and "about" routines.
 */


#include "SimTKcommon.h"

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
                        SimTK_SIMMATH_BUILD_VERSION)

#define GET_COPYRIGHT_STRING \
    MAKE_COPYRIGHT_STRING(SimTK_SIMMATH_COPYRIGHT_YEARS, \
                          SimTK_SIMMATH_AUTHORS)

#define GET_SVN_REVISION_STRING \
    MAKE_STRING(SimTK_SIMMATH_SVN_REVISION)

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

void SimTK_version_simmath(int* major, int* minor, int* build) {
    static const char* l = "SimTK library="   GET_LIBRARY_STRING;
    static const char* t = "SimTK type="      GET_TYPE_STRING;
    static const char* d = "SimTK debug="     GET_DEBUG_STRING;
    static const char* v = "SimTK version="   GET_VERSION_STRING;
    static const char* r = "SimTK svn_revision=" GET_SVN_REVISION_STRING;
    static const char* c = "SimTK copyright=" GET_COPYRIGHT_STRING;

    if (major) *major = SimTK_SIMMATH_MAJOR_VERSION;
    if (minor) *minor = SimTK_SIMMATH_MINOR_VERSION;
    if (build) *build = SimTK_SIMMATH_BUILD_VERSION;

    // Force statics to be present in the binary (Release mode otherwise 
    // optimizes them away).
    volatile int i=0;
    if (i) { // never true, but compiler doesn't know ...
        *major = *l + *t + *d + *v + *r + *c;
    }
}

void SimTK_about_simmath(const char* key, int maxlen, char* value) {
    const std::string version("version");
    const std::string library("library");
    const std::string type("type");
    const std::string copyright("copyright");
    const std::string svn_revision("svn_revision");
    const std::string authors("authors");
    const std::string debug("debug");

    if (maxlen <= 0 || value==0) return;
    value[0] = '\0'; // in case we don't find a match
    if (key==0) return;

    // downshift the key
    std::string skey(key);
    for (size_t i=0; i<skey.size(); ++i)
        skey[i] = std::tolower(skey[i]);

    char* v = 0;
    if      (skey == version)      std::string v(GET_VERSION_STRING);
    else if (skey == library)      std::string v(GET_LIBRARY_STRING);
    else if (skey == type)         std::string v(GET_TYPE_STRING);
    else if (skey == copyright)    std::string v(GET_COPYRIGHT_STRING);
    else if (skey == svn_revision) std::string v(GET_SVN_REVISION_STRING);
    else if (skey == authors)      std::string v(GET_AUTHORS_STRING);
    else if (skey == debug)        std::string v(GET_DEBUG_STRING);

    if (v) {
        std::strncpy(value,v,maxlen-1);
        value[maxlen-1] = '\0'; // in case we ran out of room
    }
}

}
