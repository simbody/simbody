#ifndef SimTK_SIMBODY_COMMON_H_
#define SimTK_SIMBODY_COMMON_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
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
Every Simbody header and source file should include this header before
any other Simbody header. **/

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

#if defined(_WIN32) && defined(_MSC_VER)
    #pragma warning(disable:4231) // need to use 'extern' template explicit instantiation
    #if defined(SimTK_SIMBODY_BUILDING_SHARED_LIBRARY)
        #define SimTK_SIMBODY_EXPORT __declspec(dllexport)
        // Keep MS VC++ quiet when it tries to instantiate incomplete template classes in a DLL.
        #pragma warning(disable:4661)
    #elif defined(SimTK_SIMBODY_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_SIMBODY_EXPORT
    #else
        #define SimTK_SIMBODY_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #define SimTK_SIMBODY_EXPORT // Linux, Mac, MinGW
#endif

// Every SimTK library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_SIMBODY_EXPORT void SimTK_version_simbody(int* major, int* minor, int* build);
    SimTK_SIMBODY_EXPORT void SimTK_about_simbody(const char* key, int maxlen, char* value);
}

namespace SimTK {
    
    // MATTER SUBSYSTEM-GLOBAL INDEX TYPES


/** @class SimTK::MobilizedBodyIndex    
This is for arrays indexed by mobilized body number within a subsystem
(typically the SimbodyMatterSubsystem).\ It is assigned when a MobilizedBody is
added to a subsystem.\ You can abbreviate this as MobodIndex if you prefer.
These will be used of course to index MobilizedBody objects but also to index 
many different kinds of information that may be stored on a per-MobilizedBody 
basis. @see SimTK::GroundIndex, SimTK::MobodIndex 
@ingroup UniqueIndexTypes **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizedBodyIndex);

/** This is the approved abbeviation for MobilizedBodyIndex.\ Feel free to
use it if you get tired of typing or seeing the full name. 
@ingroup UniqueIndexTypes **/
typedef MobilizedBodyIndex MobodIndex;

/** This is the MobilizedBodyIndex corresponding to the unique Ground body; 
its index is always zero. 
@ingroup UniqueIndexTypes **/
static const MobilizedBodyIndex GroundIndex(0);

/** @class SimTK::ConstraintIndex    
This is for arrays indexed by constraint number within a subsystem (typically
the SimbodyMatterSubsystem).\ It is assigned when a Constraint is added to the
subsystem.
@ingroup UniqueIndexTypes **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstraintIndex);

// TODO: experimental
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UnilateralContactIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(UnilateralSpeedConstraintIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(BoundedSpeedConstraintIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstraintLimitedFrictionIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(StateLimitedFrictionIndex);

// TODO: This is for arrays indexed by MatterSubsystem-global ParticleIndex, 
// as yet to be defined.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ParticleIndex);

// Constrained Bodies in constraints where the Ancestor body is not Ground (we 
// call these "Ancestor Constrained Bodies") require some additional cached 
// data, such as their orientations and velocities in the Ancestor frame, so 
// are each allocated a slot in pools of that data. Those pools are indexed by
// this type.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(AncestorConstrainedBodyPoolIndex);

// This is for "u-squared" arrays, that is, arrays which allocate space for an 
// nuXnu block for each MobilizedBody.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(USquaredIndex);

// This is for "quaternion information" arrays, which have total dimension 
// equal to the number of quaternions currently in use as generalized 
// coordinates for modeling the Matter Subsystem's MobilizedBodies. Primarily 
// this is for storing the norm of quaternions so we need calculate them only 
// once.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(QuaternionPoolIndex);

// This is for miscellaneous Real-valued position cache data that individual
// mobilizers ask us to hold for generalized coordinate q precalculations
// (e.g. sines and cosines).
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobodQPoolIndex);

// These are for indexing the pools of prescribed q's, u's, udots, and calculated forces
// needed to produce the udots. The arrays are allocated in order of MobilizedBodyIndex, and
// then in q and u order within the mobilizer. A mobilier with prescribed positions q gets
// slots in the u and udot pools also to hold derivatives, and similarly if it is the
// velocities u that are prescribed there will be slots in the udot pools. Note that
// the Q index can be used to index qdot and qdotdot arrays if needed. Note that a 
// prescribed force is produced whenever there is a udot that is not force driven; that
// includes prescribed udots but also zero and discrete ones.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(PresQPoolIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(PresUPoolIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(PresUDotPoolIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(PresForcePoolIndex);

    // PER-MOBILIZER INDEX TYPES

/** @class SimTK::MobilizerQIndex 
The Mobilizer associated with each MobilizedBody, once modeled, has a specific number
of generalized coordinates \e q (0-7) and generalized speeds (mobilities) 
\e u (0-6).\ This is the index type for the small array of Mobilizer-local \e q's. 
@see MobilizerUIndex, QIndex 
@ingroup UniqueIndexTypes **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizerQIndex);
/** @class SimTK::MobilizerUIndex 
The Mobilizer associated with each MobilizedBody, once modeled, has a specific number
of generalized coordinates \e q (0-7) and generalized speeds (mobilities) 
\e u (0-6).\ This is the index type for the small array of Mobilizer-local \e u's. 
@see MobilizerQIndex, UIndex 
@ingroup UniqueIndexTypes **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(MobilizerUIndex);

    // PER-CONSTRAINT INDEX TYPES
    
// This is the Constraint-specific index of the MobilizedBodies which are *directly* affected
// by a constraint, through body forces or body torques on these bodies.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedBodyIndex);

// This is the Constraint-specific index of the MobilizedBodies whose mobilizers' mobilities
// can appear explicitly in constraint equations, and which are acted upon by the Constraint
// through generation of generalized (mobility) forces. Note that for a multi-dof mobilizer
// we don't select individual mobilities; it is all or nothing so we can use the MobilizedBody
// to stand for its mobilizer.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedMobilizerIndex);

// This is the Constraint-specific index of a coordinate q which can be *directly* affected
// by this constraint through generation of a mobility force on a corresponding mobility. These
// are numbered in order of ConstrainedMobilizerIndex for the mobilizers for which these are the q's.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedQIndex);

// This is the Constraint-specific index of a mobility u which can be *directly* affected
// by this constraint through generation of a mobility force. These are numbered in order
// of ConstrainedMobilizerIndex for the bodies for which these are the u's.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ConstrainedUIndex);


// This is the Constraint-specific index of a coordinate q which can be involved in any
// constraint equation of this constraint, either directly through ConstrainedMobilizers
// or indirectly as a result of its effects on ConstrainedBodies (that is, this list
// includes all the ConstraintQIndex entries above, plus possibly many more). These are in sorted
// order by subsystem-wide QIndex, and each QIndex appears at most once.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ParticipatingQIndex);

// This is the Constraint-specific index of a coordinate u which can be involved in any
// constraint equation of this constraint, either directly through ConstrainedMobilizers
// or indirectly as a result of its effects on ConstrainedBodies (that is, this list
// includes all the ConstraintUIndex entries above, plus possibly many more). These are in sorted
// order by subsystem-wide UIndex, and each UIndex appears at most once.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ParticipatingUIndex);

    // SUBTREE INDEX TYPES

// And similarly for other unique Index types.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeBodyIndex);
static const SubtreeBodyIndex SubtreeAncestorIndex(0);

SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeQIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SubtreeUIndex);

    // INDEX TYPES FOR OTHER SUBSYSTEMS

/** @class SimTK::ForceIndex
This type represents the index of a Force element within its subsystem.
@relates SimTK::Force
@ingroup UniqueIndexTypes **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ForceIndex);

SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactSetIndex);


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
