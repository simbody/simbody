#ifndef SimTK_SIMBODY_COMMON_H_
#define SimTK_SIMBODY_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
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
 * Every Simbody header and source file should include this header before
 * any other Simbody header.
 */

#include "SimTKcommon.h"

#include <cassert>
#include <vector>
#include <limits>


// Shared libraries are messy in Visual Studio. We have to distinguish three
// cases:
//   (1) this header is being used to build the simbody shared library (dllexport)
//   (2) this header is being used by a *client* of the simbody shared
//       library (dllimport)
//   (3) we are building the simbody static library, or the client is
//       being compiled with the expectation of linking with the
//       simbody static library (nothing special needed)
// In the CMake script for building this library, we define one of the symbols
//     SimTK_SIMBODY_BUILDING_{SHARED|STATIC}_LIBRARY
// Client code normally has no special symbol defined, in which case we'll
// assume it wants to use the shared library. However, if the client defines
// the symbol SimTK_USE_STATIC_LIBRARIES we'll suppress the dllimport so
// that the client code can be linked with static libraries. Note that
// the client symbol is not library dependent, while the library symbols
// affect only the simbody library, meaning that other libraries can
// be clients of this one. However, we are assuming all-static or all-shared.

#ifdef WIN32
    #if defined(SimTK_SIMBODY_BUILDING_SHARED_LIBRARY)
        #define SimTK_SIMBODY_EXPORT __declspec(dllexport)
    #elif defined(SimTK_SIMBODY_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_SIMBODY_EXPORT
    #else
        #define SimTK_SIMBODY_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define SimTK_SIMBODY_EXPORT // Linux, Mac
#endif

// Every SimTK Core library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_SIMBODY_EXPORT void SimTK_version_simbody(int* major, int* minor, int* build);
    SimTK_SIMBODY_EXPORT void SimTK_about_simbody(const char* key, int maxlen, char* value);
}

namespace SimTK {

// This is just a type-safe non-negative int, augmented with a "NaN" 
// value called InvalidMobilizedBodyId. For most uses it will behave like an int,
// and it has an implicit conversion *to* int. Importantly though,
// it has no implicit conversion *from* int so you can't pass some
// other kind of number as a MobilizedBodyId.
// 
// We also predefine GroundId which is always MobilizedBodyId(0).
//
SimTK_DEFINE_UNIQUE_ID_TYPE(MobilizedBodyId)
static const MobilizedBodyId GroundId(0);

// And similarly for other unique Id types.
SimTK_DEFINE_UNIQUE_ID_TYPE(ConstraintId)
SimTK_DEFINE_UNIQUE_ID_TYPE(ParticleId)


namespace Exception {


class APIMethodFailed : public Base {
public:
    APIMethodFailed(const char* fn, int ln, String method, String cause) : Base(fn,ln)
    {
        setMessage(method + " failed because:\n  " + cause);
    }
};


// This just reports rep-level bad things up to the API level with a helpful string.
class RepLevelException : public Base {
public:
    RepLevelException(const char* fn, int ln, String message) : Base(fn,ln)
    {
        setMessage(message);
    }
};

class MobilizerCantExactlyRepresentRequestedQuantity : public Base {
public:
    MobilizerCantExactlyRepresentRequestedQuantity(const char* fn, int ln, 
       String method, MobilizedBodyId body, String quantity) : Base(fn,ln)
    {
        setMessage(method + "(): the mobilizer for body " + String((int)body)
            + " can't represent the given " + quantity + " to machine precision");
    }
private:
};

class NewtonRaphsonFailure : public Base {
public:
    NewtonRaphsonFailure(const char* fn, int ln, String msg) : Base(fn,ln)
    {
        setMessage("NewtonRaphson failure: " + msg);
    }
private:
};


class LoopConstraintConstructionFailure : public Base {
public:
    LoopConstraintConstructionFailure(const char* fn, int ln, String msg) : Base(fn,ln)
    {
        setMessage("Loop constraint construction failure: " + msg);
    }
private:
};

} // namespace SimTK::Exception



} // namespace SimTK

#endif // SimTK_SIMBODY_COMMON_H_
