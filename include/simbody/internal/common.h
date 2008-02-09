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
    #ifdef _MSC_VER
    #pragma warning(disable:4231) // need to use 'extern' template explicit instantiation
    #endif
    #if defined(SimTK_SIMBODY_BUILDING_SHARED_LIBRARY)
        #define SimTK_SIMBODY_EXPORT __declspec(dllexport)
        // Keep MS VC++ quiet when it tries to instantiate incomplete template classes in a DLL.
        #ifdef _MSC_VER
        #pragma warning(disable:4661)
        #endif
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

    
// There are various arrays containing data about certain categories of objects. We want them
// indexable by simple ints for speed, but it would be a disaster to accidentally use the
// wrong int as an index! So we define unique index types here for accessing each category
// to help stay out of trouble. We use different index types for accessing the full array
// of all coordinates (for example) of a Matter Subsystem than we do for accessing subsets
// of those coordinates.
//
// A unique index type is just a type-safe non-negative int, augmented with a "NaN" 
// value called InvalidBLAH where BLAH is the type name. For most uses it will behave
// like an int, and it has an implicit conversion *to* int. Importantly though,
// it has no implicit conversion *from* int so you can't pass a plain int or any
// other Index type to an argument expecting a certain Index type.


    // MATTER SUBSYSTEM-GLOBAL INDEX TYPES

// This is for arrays indexed by MatterSubsystem-global MobilizedBodyIndex, assigned
// when a MobilizedBody is added to a MatterSubsystem. These will
// be used of course to hold MobilizedBodies but also to hold many different kinds of
// information that may be stored on a per-MobilizedBody basis.
// We also predefine GroundIndex here which is always MobilizedBodyIndex(0).
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizedBodyIndex)
static const MobilizedBodyIndex GroundIndex(0);

// This is for arrays indexed by MatterSubsystem-global ConstraintIndex, assigned when
// a Constraint is added to a Matter Subsystem.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstraintIndex)

// TODO: This is for arrays indexed by MatterSubsystem-global ParticleIndex, as yet to be defined.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ParticleIndex)

// This is for "q-like" arrays, that is, arrays which inherently have the same dimension as
// the totoal number of generalized coordinates for the whole Matter Subsystem.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QIndex)    // an index into generalized coordinates q

// This is for "u-like" arrays, that is, arrays which inherently have the same dimension as
// the total number of mobilities (generalized speeds) for the whole Matter Subsystem. This
// includes both u and udot.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UIndex)    // an index into generalized speeds u (and accelerations udot)

// This is for "u-squared" arrays, that is, arrays which allocate space for an nuXnu block
// for each MobilizedBody.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(USquaredIndex)

// This is for "quaternion information" arrays, which have total dimension equal to the
// number of quaternions currently in use as generalized coordinates for modeling the Matter
// Subsystem's MobilizedBodies. Primarily this is for storing the norm of quaternions so we need
// calculate them only once.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QuaternionPoolIndex)

// This is for "angle information" arrays, which have dimension equal to the total number
// of generalized coordinates currently in use which are just angles in radians. For those
// we want to precalculate some expensive things like sine and cosine so that we can just
// do that once.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(AnglePoolIndex)

    // PER-MOBILIZED BODY INDEX TYPES

// Each MobilizedBody, once modeled, has a specific number of generalized coordinates q
// (0-7) and generalized speeds (mobilities) u (0-6). These are the index types for the
// small arrays of MobilizedBody-local q's and u's.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizedBodyQIndex)
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizedBodyUIndex)

    // PER-CONSTRAINT INDEX TYPES
    
// This is the Constraint-specific index of the MobilizedBodies which are *directly* affected
// by a constraint. That is, the Constraint expects to apply constraint forces as body forces
// on these bodies or as mobility forces on these bodies' mobilizers.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedBodyIndex)

// This is the Constraint-specific index of a coordinate q which can be *directly* affected
// by this constraint through generation of a mobility force on a corresponding mobility. These
// are numbered in order of ConstrainedBodyIndex for the bodies for which these are the q's.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedQIndex)

// This is the Constraint-specific index of a mobility u which can be *directly* affected
// by this constraint through generation of a mobility force. These are numbered in order
// of ConstrainedBodyIndex for the bodies for which these are the u's.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedUIndex)

    // SUBTREE INDEX TYPES

// And similarly for other unique Index types.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeBodyIndex)
static const SubtreeBodyIndex SubtreeAncestorIndex(0);

SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeQIndex)
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeUIndex)




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
       String method, MobilizedBodyIndex body, String quantity) : Base(fn,ln)
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
