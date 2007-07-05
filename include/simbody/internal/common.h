#ifndef SimTK_SIMBODY_COMMON_H_
#define SimTK_SIMBODY_COMMON_H_

/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * Every Simbody header and source file should include this header before
 * any other Simbody header.
 */

#include "SimTKcommon.h"

#include <cassert>
#include <vector>
#include <limits>

// TODO: move to SimTKcommon
#define SimTK_PIMPL_DOWNCAST(Derived, Parent)           \
    static bool           isInstanceOf(const Parent&);  \
    static const Derived& downcast(const Parent&);      \
    static Derived&       updDowncast(Parent&)

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

static const int invalidId = -1111111111;

/**
 * This is just a type-safe non-negative int, augmented with a "NaN" 
 * value called InvalidSubsystemId. For most uses it will behave like an int,
 * and it has an implicit conversion *to* int. Importantly though,
 * it has no implicit conversion *from* int so you can't pass some
 * other kind of number as a SubsystemId.
 */
class SubsystemId {
    int id;
public:
    inline SubsystemId();
    inline explicit SubsystemId(int i);
    operator int() const {return id;}
    bool isValid() const {return id>=0;}
    const SubsystemId& operator++() {assert(id>=0); ++id;return *this;}           // prefix
    SubsystemId operator++(int)     {assert(id>=0); ++id; return SubsystemId(id-1);} // postfix
    const SubsystemId& operator--() {assert(id>=1); --id;return *this;}           // prefix
    SubsystemId operator--(int)     {assert(id>=1); --id; return SubsystemId(id+1);} // postfix
};
static const SubsystemId InvalidSubsystemId(invalidId);
inline SubsystemId::SubsystemId() : id(InvalidSubsystemId) { }
inline SubsystemId::SubsystemId(int i) : id(i) {
    assert(i>=0 || i==invalidId);
}

/**
 * This is just a type-safe non-negative int, augmented with a "NaN" 
 * value called InvalidMobilizedBodyId. For most uses it will behave like an int,
 * and it has an implicit conversion *to* int. Importantly though,
 * it has no implicit conversion *from* int so you can't pass some
 * other kind of number as a MobilizedBodyId.
 * 
 * We also predefine GroundId which is always MobilizedBodyId(0).
 */
class MobilizedBodyId {
    int id;
public:
    inline MobilizedBodyId();
    inline explicit MobilizedBodyId(int i);
    operator int() const {return id;}
    bool isValid() const {return id>=0;}
    const MobilizedBodyId& operator++() {assert(id>=0); ++id;return *this;}           // prefix
    MobilizedBodyId operator++(int)     {assert(id>=0); ++id; return MobilizedBodyId(id-1);} // postfix
    const MobilizedBodyId& operator--() {assert(id>=1); --id;return *this;}           // prefix
    MobilizedBodyId operator--(int)     {assert(id>=1); --id; return MobilizedBodyId(id+1);} // postfix
};
static const MobilizedBodyId GroundId(0);
static const MobilizedBodyId InvalidMobilizedBodyId(invalidId);
inline MobilizedBodyId::MobilizedBodyId() : id(InvalidMobilizedBodyId) { }
inline MobilizedBodyId::MobilizedBodyId(int i) : id(i) {
    assert(i>=0 || i==invalidId);
}

/**
 * This is just a type-safe non-negative int, augmented with a "NaN" 
 * value called InvalidParticleId. For most uses it will behave like an int,
 * and it has an implicit conversion *to* int. Importantly though,
 * it has no implicit conversion *from* int so you can't pass some
 * other kind of number as  a ParticleId.
 */
class ParticleId {
    int id;
public:
    inline ParticleId();
    inline explicit ParticleId(int i);
    operator int() const {return id;}
    bool isValid() const {return id>=0;}
    const ParticleId& operator++() {assert(id>=0); ++id;return *this;}             // prefix
    ParticleId operator++(int)     {assert(id>=0); ++id; return ParticleId(id-1);} // postfix
    const ParticleId& operator--() {assert(id>=1); --id;return *this;}             // prefix
    ParticleId operator--(int)     {assert(id>=1); --id; return ParticleId(id+1);} // postfix
};
static const ParticleId InvalidParticleId(invalidId);
inline ParticleId::ParticleId() : id(InvalidParticleId) { }
inline ParticleId::ParticleId(int i) : id(i) {
    assert(i>=0 || i==invalidId);
}

/**
 * This is just a type-safe non-negative int, augmented with a "NaN" 
 * value called InvalidConstraintId. For most uses it will behave like an int,
 * and it has an implicit conversion *to* int. Importantly though,
 * it has no implicit conversion *from* int so you can't pass some
 * other kind of number as  a ConstraintId.
 */
class ConstraintId {
    int id;
public:
    inline ConstraintId();
    inline explicit ConstraintId(int i);
    operator int() const {return id;}
    bool isValid() const {return id>=0;}
    const ConstraintId& operator++() {assert(id>=0); ++id;return *this;}          // prefix
    ConstraintId operator++(int)     {assert(id>=0); ++id; return ConstraintId(id-1);}  // postfix
    const ConstraintId& operator--() {assert(id>=1); --id;return *this;}          // prefix
    ConstraintId operator--(int)     {assert(id>=1); --id; return ConstraintId(id+1);}  // postfix
};
static const ConstraintId InvalidConstraintId(invalidId);
inline ConstraintId::ConstraintId() : id(InvalidConstraintId) { }
inline ConstraintId::ConstraintId(int i) : id(i) {
    assert(i>=0 || i==invalidId);
}


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

/*
class FeatureAPIMethodFailed : public Base {
public:
    FeatureAPIMethodFailed(const char* fn, int ln, String fullFeatureName,
        String method, String argInfo, String cause) : Base(fn,ln)
    {
        setMessage("Feature('" + fullFeatureName + "')."
                   + method + "(" + argInfo + ") failed because:\n  " + cause);
    }
};

class PlacementAPIMethodFailed : public Base {
public:
    PlacementAPIMethodFailed(const char* fn, int ln,
        String method, String argInfo, String cause) : Base(fn,ln)
    {
        setMessage("Placement."
                   + method + "(" + argInfo + ") failed because:\n  " + cause);
    }
};

class AccessToInvalidPlacementValue : public Base {
public:
    AccessToInvalidPlacementValue(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("An attempt was made to access an invalid value for Feature '" + featureName 
            + "'s Placement. Features must be realized before their values can be obtained.");
    }
private:
};

class EmptySubsystemPathname : public Base {
public:
    EmptySubsystemPathname(const char* fn, int ln) : Base(fn,ln)
    {
        setMessage("Subsystem pathname was empty.");
    }
private:
};

class IllegalSubsystemPathname : public Base {
public:
    IllegalSubsystemPathname(const char* fn, int ln, String pathname, String badseg) : Base(fn,ln)
    {
        setMessage("Subsystem pathname '" + pathname + "' is illegal at segment '" + badseg + "'.");
    }
private:
};

class SubsystemNameNotFound : public Base {
public:
    SubsystemNameNotFound(const char* fn, int ln, String subname, String parentname) : Base(fn,ln)
    {
        setMessage("Can't find any subsystem named '" + subname + "' in '" + parentname + "'.");
    }
private:
};

class FeatureHasAlreadyBeenPlaced : public Base {
public:
    FeatureHasAlreadyBeenPlaced(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("Can't place Feature " + featureName + " because it already has a Placement.");
    }
private:
};

class FeatureCantBeUsedAsPlacement : public Base {
public:
    FeatureCantBeUsedAsPlacement(const char* fn, int ln, String featureName, String featureTypeName,
        String placementTypeNeeded) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
                   + "' as a " + placementTypeNeeded + " Placement.");
    }
private:
};

class PlacementCantBeConvertedToRightType : public Base {
public:
    PlacementCantBeConvertedToRightType(const char* fn, int ln, 
        String placementTypeNeeded, String srcType, String srcDescription) : Base(fn,ln)
    {
        setMessage("Can't use the supplied " + srcType + " Placement to create a "
            + placementTypeNeeded + " Placement. Supplied Placement was:\n"
            + srcDescription);
    }
private:
};

class OnlyFeaturesHavePlacements : public Base {
public:
    OnlyFeaturesHavePlacements(const char* fn, int ln, String subsysName) : Base(fn,ln)
    {
        setMessage("An attempt was made to access a Placement for Subsystem " + subsysName 
                   + " which is not a Feature. Only Features have Placements.");
    }
private:
};

class ExpectedFeatureIndexButGotSubsystem : public Base {
public:
    ExpectedFeatureIndexButGotSubsystem(const char* fn, int ln, String subsysName, int index) : Base(fn,ln)
    {
        setMessage("Child Subsystem " + String(index) + " of Subsystem " + subsysName 
                   + " is not a Feature, but this operation expects a Feature.");
    }
private:
};

class ExpectedFeatureButGotSubsystem : public Base {
public:
    ExpectedFeatureButGotSubsystem(const char* fn, int ln, String subsysName) : Base(fn,ln)
    {
        setMessage("Subsystem " + subsysName 
                   + " is not a Feature, but this operation expects a Feature.");
    }
private:
};

class FeatureIsNotPlacedOnABody : public Base {
public:
    FeatureIsNotPlacedOnABody(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("Feature " + featureName + "'s Placement is not owned by a Body.");
    }
private:
};

class ExpectedFeatureToHaveFeatureReferencePlacement : public Base {
public:
    ExpectedFeatureToHaveFeatureReferencePlacement(const char* fn, int ln, String featureName) : Base(fn,ln)
    {
        setMessage("Feature " + featureName + "'s Placement is not a Feature reference.");
    }
private:
};

class NotAFeatureReferencePlacement : public Base {
public:
    NotAFeatureReferencePlacement(const char* fn, int ln) : Base(fn,ln)
    {
        setMessage("Expected a Feature reference Placement.");
    }
private:
};

class ExpectedFeaturePrototypeButGotSubsystem : public Base {
public:
    ExpectedFeaturePrototypeButGotSubsystem(const char* fn, int ln, 
        String subsysName, String newName, String protoName) : Base(fn,ln)
    {
        setMessage("An attempt was made to add a Feature " + newName + " to Subsystem "
                   + subsysName + " using a prototype named " + protoName 
                   + " but the prototype was not a Feature.");
    }
private:
};

class NoFeatureLevelPlacementForThisKindOfFeature : public Base {
public:
    NoFeatureLevelPlacementForThisKindOfFeature(const char* fn, int ln, 
        String featureName, String featureTypeName) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
                   + "' as a Placement because this kind of Feature has no"
                     " Feature-level Placement. Did you mean to use one of its Subfeatures?");
    }
private:
};

class IndexOutOfRangeForFeaturePlacementReference : public Base {
public:
    IndexOutOfRangeForFeaturePlacementReference(const char* fn, int ln, 
        String featureName, String featureTypeName, int index) : Base(fn,ln)
    {
        setMessage("Index " + String(index)
                   + " is out of range for a reference to the Placement of "
                   + featureTypeName + " Feature '" + featureName + "'.");
    }
private:
};

class PlacementCantBeUsedForThisFeature : public Base {
public:
    PlacementCantBeUsedForThisFeature(const char* fn, int ln,
                                      String placementTypeName,
                                      String featureName, String featureTypeName) : Base(fn,ln)
    {
        setMessage("Can't use " + placementTypeName + "Placement for "
                   + featureTypeName + "Feature '" + featureName 
                   + "'.");
    }
private:
};

class AFeatureCantBePlacedOnItself : public Base {
public:
    AFeatureCantBePlacedOnItself(const char* fn, int ln,
                                 String featureName) : Base(fn,ln)
    {
        setMessage("Can't place Feature '" + featureName + "' on a Placement expression "
                   "which directly or indirectly depends on that Feature.");
    }
private:
};

class UnaryOperationNotAllowedForPlacementType : public Base {
public:
    UnaryOperationNotAllowedForPlacementType(const char* fn, int ln, 
                                             String opName, String placementTypeName) 
      : Base(fn,ln)
    {
        setMessage("Unary operator '" + opName + "' can't be used on a "
                   + placementTypeName + "Placement.");
    }
};

class InfixPlacementOperationNotAllowed : public Base {
public:
    InfixPlacementOperationNotAllowed(const char* fn, int ln, 
                                      String leftTypeName, String opName,
                                      String rightTypeName) 
      : Base(fn,ln)
    {
        setMessage("Operation (" + leftTypeName + "Placement " + opName + 
                   " " + rightTypeName + "Placement) is not supported.");
    }
};

class FeatureUsedAsFramePlacementMustBeOnFrame : public Base {
public:
    FeatureUsedAsFramePlacementMustBeOnFrame(const char* fn, int ln, 
        String featureName, String featureTypeName,
        String missingPlacementType) : Base(fn,ln)
    {
        setMessage("Can't use " + featureTypeName + " Feature '" + featureName 
         + "' as a Frame Placement because it doesn't have a parent Frame"
           " from which to inherit the " + missingPlacementType
           + " Placement.");
    }
private:
};

class PlacementMustBeLocal : public Base {
public:
    PlacementMustBeLocal(const char* fn, int ln, String method, 
                         String hostFeature, String offendingFeature) : Base(fn,ln)
    {
        setMessage(method + ": can't add placement expression to Feature '" 
            + hostFeature + "' because it references Feature '"
            + offendingFeature 
            + "' (and possibly others) which is not a descendent of '"
            + hostFeature + "'.");
    }
private:
};

class FeatureAndPlacementOnDifferentTrees : public Base {
public:
    FeatureAndPlacementOnDifferentTrees(const char* fn, int ln, 
                                        String hostFeature, String offendingFeature) : Base(fn,ln)
    {
        setMessage(
            "can't place Feature '" + hostFeature 
            + "' because the supplied placement references Feature '" + offendingFeature
            + "' (and possibly others) and there is no common ancestor.");
    }
private:
};
*/

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
