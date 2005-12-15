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

/**@file
 * Implementation of non-inline PlacementRep methods.
 */

#include "SimbodyCommon.h"
#include "Placement.h"
#include "Feature.h"
#include "PlacementRep.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::cout;
using std::endl;
using std::ostream;

namespace simtk {

    // PLACEMENT EXPR //

void PlacementExpr::exprRealize(/*State,*/ Stage g) const {
    for (size_t i=0; i < args.size(); ++i)
        args[i].realize(/*State,*/ g);
}

bool PlacementExpr::exprIsConstant() const {
    for (size_t i=0; i < args.size(); ++i)
        if (!args[i].isConstant()) return false;
    return true;
}

bool PlacementExpr::exprDependsOn(const Feature& f) const {
    for (size_t i=0; i < args.size(); ++i)
        if (args[i].dependsOn(f))
            return true;
    return false;
}

const Subsystem* 
PlacementExpr::exprFindAncestorSubsystem(const Subsystem& youngestAllowed) const {
    const Subsystem* ancestor = &youngestAllowed;
    for (size_t i=0; ancestor && i < args.size(); ++i)
        ancestor = args[i].getRep().findAncestorSubsystem(*ancestor);
    return ancestor;
}

const Subsystem* 
PlacementExpr::exprFindPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const {
    const Subsystem* valueOwner = &youngestAllowed;
    for (size_t i=0; valueOwner && i < args.size(); ++i)
        valueOwner = args[i].getRep().findPlacementValueOwnerSubsystem(*valueOwner);
    return valueOwner;
}

std::string
PlacementExpr::exprToString(const std::string& linePrefix) const {
    std::stringstream s;
    s << func->getOpName() << "(";
    for (size_t i=0; i<args.size(); ++i)
        s << (i>0?", ":"") 
            << args[i].getRep().toString(linePrefix);
    s << ")";
    return s.str();
}

bool PlacementExpr::exprIsLimitedToSubtree
    (const Subsystem& root, const Feature*& offender) const
{
    for (size_t i=0; i<args.size(); ++i)
        if (!args[i].getRep().isLimitedToSubtree(root,offender))
            return false;
    return true;
}

void PlacementExpr::exprRepairFeatureReferences
    (const Subsystem& oldRoot, const Subsystem& newRoot)
{
    for (size_t i=0; i<args.size(); ++i)
        args[i].updRep().repairFeatureReferences(oldRoot,newRoot);
}

    // FEATURE REFERENCE //

FeatureReference::FeatureReference(const Feature& f, int i) 
  : feature(&f), index(i)
{
    const PlacementType t = f.getRep().getRequiredPlacementType();
    const int nElements = PlacementRep::getNIndicesAllowed(t);

    // TODO: should this throw a nice message, or should we check higher up?
    assert(nElements > 0 // i.e., not void
        && (i == -1 || (0 <= i && i < nElements)));
}


void FeatureReference::refRealize(/*State,*/ Stage g) const {
    if (getReferencedFeature().hasPlacement())
        getReferencedFeature().getPlacement().realize(/*State,*/ g);
}

const Subsystem* 
FeatureReference::refFindAncestorSubsystem(const Subsystem& youngestAllowed) const {
    assert(feature);
    return SubsystemRep::findYoungestCommonAncestor(youngestAllowed, *feature);
}

// If the referenced Feature has no placement, we can't evaluate this placement so
// there is *no* acceptable owner for the value. If it does have a placement, then
// the youngest place for its value is the common ancestor of youngestAllowed and
// the placement's owner.
const Subsystem* 
FeatureReference::refFindPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const {
    assert(feature);
    if (!feature->hasPlacement())
        return 0;
    const PlacementRep& pr = feature->getPlacement().getRep();
    assert(pr.hasOwner());  // a feature can't be placed on an unowned placement
    const Subsystem* newYoungest = SubsystemRep::findYoungestCommonAncestor(youngestAllowed,
                                                                            pr.getOwner());
    assert(newYoungest);    // wrong tree???

    return feature->getPlacement().getRep().findPlacementValueOwnerSubsystem(*newYoungest);
}

// Check that this feature is on the feature subtree rooted by "ancestor". If
// not return a pointer to this feature for use in a friendly error message.
// If this is the right tree, we return true with offender==NULL.
bool FeatureReference::refIsLimitedToSubtree
    (const Subsystem& root, const Feature*& offender) const
{
    assert(feature);
    if (SubsystemRep::isSubsystemInSubsystemTree(root, *feature)) {
        offender = 0;
        return true;
    }
    offender = feature;
    return false;
}

void FeatureReference::refRepairFeatureReferences
    (const Subsystem& oldRoot, const Subsystem& newRoot)
{
    assert(feature);
    const Subsystem* corrFeature = SubsystemRep::findCorrespondingSubsystem(oldRoot,*feature,newRoot);
    assert(corrFeature && Feature::isInstanceOf(*corrFeature));
    feature = &Feature::downcast(*corrFeature);
}

std::string FeatureReference::refToString(const std::string&) const {
    std::stringstream s;
    const PlacementType t = feature 
        ? feature->getRep().getRequiredPlacementType()
        : InvalidPlacementType;

    s << "Ref<" << PlacementRep::getPlacementTypeName(t) << ">[";
    s << (feature ? feature->getFullName()
                    : std::string("NULL FEATURE"));
    s << "]"; 
    if (index != -1) s << "[" << index << "]";
    return s.str();
}

PlacementType FeatureReference::refGetPlacementType() const {
    assert(feature);
    const PlacementType whole = (*feature).getRep().getRequiredPlacementType();
    return index == -1 ? whole : PlacementRep::getIndexedPlacementType(whole, index);
}

    // PLACEMENT REP //


// We have just copied a Feature tree and this PlacementRep is the new copy. If
// it had a valueSlot, that valueSlot is still pointing into the old Feature tree
// and needs to be repaired to point to the corresponding valueSlot in the new tree.
void PlacementRep::repairValueReference(const Subsystem& oldRoot, const Subsystem& newRoot) {
    if (valueSlot) {
        valueSlot = const_cast<PlacementValue*>
                        (FeatureRep::findCorrespondingPlacementValue
                                    (oldRoot,*valueSlot,newRoot));
        if (valueSlot)
            valueSlot->updRep().setClientPlacement(getMyHandle());
    }
}

// These are the default implementations for the generic operators. Any 
// concrete class which thinks it knows how to perform one of these
// operations should override.

/*virtual*/ Placement PlacementRep::genericNegate() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "-", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericAbs() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "abs", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericSqrt() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "sqrt", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericExp() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "exp", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericLog() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "log", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericSin() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "sin", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericCos() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "cos", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericAsin() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "asin", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericAcos() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "acos", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericLength() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "length", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericNormalize() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "normalize", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericAdd(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "+", getPlacementTypeName(r.getRep().getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericSub(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "-", getPlacementTypeName(r.getRep().getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericMul(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "*", getPlacementTypeName(r.getRep().getPlacementType()));
}
/*virtual*/ Placement PlacementRep::genericDvd(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "/", getPlacementTypeName(r.getRep().getPlacementType()));
}

/*virtual*/ Placement PlacementRep::genericDistance(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "distance", getPlacementTypeName(r.getRep().getPlacementType()));
}

/*virtual*/ Placement PlacementRep::genericAngle(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "angle", getPlacementTypeName(r.getRep().getPlacementType()));
}

/*virtual*/ Placement PlacementRep::genericDotProduct(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "dot", getPlacementTypeName(r.getRep().getPlacementType()));
}

/*virtual*/ Placement PlacementRep::genericCrossProduct(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(getPlacementType()),
                "cross", getPlacementTypeName(r.getRep().getPlacementType()));
}

/*virtual*/ RealPlacement PlacementRep::castToRealPlacement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToRealPlacement()", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ Vec3Placement PlacementRep::castToVec3Placement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToVec3Placement()", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ StationPlacement PlacementRep::castToStationPlacement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToStationPlacement()", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ DirectionPlacement PlacementRep::castToDirectionPlacement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToDirectionPlacement()", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ OrientationPlacement PlacementRep::castToOrientationPlacement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToOrientationPlacement()", getPlacementTypeName(getPlacementType()));
}
/*virtual*/ FramePlacement PlacementRep::castToFramePlacement() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "castToFramePlacement()", getPlacementTypeName(getPlacementType()));
}

/*static*/ const char* 
PlacementRep::getPlacementTypeName(PlacementType t) {
    switch(t) {
    case InvalidPlacementType:      return "INVALID";
    case VoidPlacementType:         return "void";
    case BoolPlacementType:         return "bool";
    case IntPlacementType:          return "int";
    case RealPlacementType:         return "Real";
    case StationPlacementType:      return "Station";
    case DirectionPlacementType:    return "Direction";
    case OrientationPlacementType:  return "Orientation";
    case FramePlacementType:        return "Frame";
    case Vec2PlacementType:         return "Vec2";
    case Vec3PlacementType:         return "Vec3";
    case Mat33PlacementType:        return "Mat33";
    default: return "ILLEGAL PLACEMENT TYPE";
    };
}

/*static*/ int 
PlacementRep::getNIndicesAllowed(PlacementType t) {
    switch(t) {
    case VoidPlacementType:         return 0; // can't use at all

    case BoolPlacementType:
    case IntPlacementType:
    case RealPlacementType:         return 1; // no index or index==0 OK

    case Vec2PlacementType:         return 2; // 2 Reals

    case Vec3PlacementType:
    case StationPlacementType:
    case DirectionPlacementType:    return 3; // 3 Reals

    case Mat33PlacementType:        return 3; // 3 Vec3's (columns)
    case OrientationPlacementType:  return 3; // 3 Directions

    case FramePlacementType:        return 2; // Orientation, Station

    default: 
        assert(false);
    };
    //NOTREACHED
    return -1;
}

// If a PlacementType is indexed with this index, 
// what is the resulting PlacementType?
/*static*/ PlacementType 
PlacementRep::getIndexedPlacementType(PlacementType t, int i) {
    if (i == -1) 
        return t;   // -1 means not indexed, i.e., the whole thing

    assert(0 <= i && i <= getNIndicesAllowed(t));
    switch(t) {
    case StationPlacementType:
    case DirectionPlacementType: 
    case Vec2PlacementType:
    case Vec3PlacementType:         
        return RealPlacementType;  

    case Mat33PlacementType:
        return Vec3PlacementType;

    case OrientationPlacementType: 
        return DirectionPlacementType;

    case FramePlacementType:
        return i==0 ? OrientationPlacementType : StationPlacementType;

    default: assert(false);
        //NOTREACHED
    };
    //NOTREACHED
    return InvalidPlacementType;
}

void PlacementRep::checkPlacementConsistency(const Subsystem* expOwner, 
                                             int              expIndexInOwner,
                                             const Subsystem& expRoot) const
{
    cout << "CHECK PLACEMENT CONSISTENCY FOR PlacementRep@" << this << endl;
    if (!myHandle) 
        cout << "*** NO HANDLE ***" << endl;
    else if (myHandle->rep != this)
        cout << "*** Handle->rep=" << myHandle->rep << " which is *** WRONG ***" << endl;
    if (owner != expOwner)
        cout << "*** WRONG OWNER@" << owner << "; should have been " << expOwner << endl;
    if (indexInOwner != expIndexInOwner)
        cout << "*** WRONG INDEX " << indexInOwner << "; should have been " << expIndexInOwner << endl;

    if (expOwner == 0) {
        if (client)
          cout << "*** UNOWNED PLACEMENT HAD CLIENT " << client->getFullName() << " ***" << endl;
    } else {
        if (!client) 
            cout << "*** NO CLIENT ***" << endl;
        else if (!client->hasPlacement())
            cout << "*** CLIENT " << client->getFullName() << " HAS NO PLACEMENT??? ***" << endl;
        else if (&(client->getPlacement().getRep()) != this)
            cout << "*** CLIENT " << client->getFullName() << " HAS WRONG PLACEMENT@" 
                << &client->getPlacement().getRep() << endl;
    }

    const Feature* offender;
    if (!isLimitedToSubtree(expRoot, offender)) {
        cout << "*** Placement referenced Feature '" << offender->getFullName() << "' in wrong tree" << endl;
        cout << "*** Root should have been @" << &expRoot << ", not " << 
            &offender->getRep().findRootSubsystem() << endl;
    }

    if (hasValueSlot()) {
        std::string nm = client ? client->getFullName() : "(NO CLIENT)";
        if (!getValueSlot().hasOwner())
            cout << "*** Placement for " << nm << "'s value slot is unowned." << endl;
        else if (!getValueSlot().getOwner().getRep().findRootSubsystem().isSameSubsystem(expRoot))
            cout << "*** Placement for " << nm << "'s value slot is in wrong tree." << endl;
    }
}


    // REAL PLACEMENT REP //

// result = -real (result is always real)
Placement RealPlacementRep::genericNegate() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(-calcValue())
                     : (RealPlacementRep*)RealExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = abs(real) (result is always real)
Placement RealPlacementRep::genericAbs() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::abs(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::absOp(getMyHandle());
    return Placement(result);
}

// result = sqrt(real) (result is always real)
Placement RealPlacementRep::genericSqrt() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sqrt(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sqrtOp(getMyHandle());
    return Placement(result);
}

// result = exp(real) (result is always real)
Placement RealPlacementRep::genericExp() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::exp(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::expOp(getMyHandle());
    return Placement(result);
}

// result = log(real) (natural log, result is always real)
Placement RealPlacementRep::genericLog() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::log(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::logOp(getMyHandle());
    return Placement(result);
}

// result = sin(real) (result is always real)
Placement RealPlacementRep::genericSin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sin(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sinOp(getMyHandle());
    return Placement(result);
}

// result = cos(real) (result is always real)
Placement RealPlacementRep::genericCos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::cos(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::cosOp(getMyHandle());
    return Placement(result);
}

// result = asin(real) (result is always real)
Placement RealPlacementRep::genericAsin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::asin(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::asinOp(getMyHandle());
    return Placement(result);
}

// result = acos(real) (result is always real)
Placement RealPlacementRep::genericAcos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::acos(calcValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::acosOp(getMyHandle());
    return Placement(result);
}

// result = real + placement
// We support only 
//   real = real + real
Placement RealPlacementRep::genericAdd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        RealPlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcValue()+rp.getRep().calcValue())
                : (RealPlacementRep*)RealExprPlacementRep::addOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericAdd(r);    // die
}

// result = real - placement
// We support only 
//   real = real - real
Placement RealPlacementRep::genericSub(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        RealPlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcValue()-rp.getRep().calcValue())
                : (RealPlacementRep*)RealExprPlacementRep::subOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericSub(r);    // die
}

// result = real * placement
// We support
//    real = real * real
//    vec3 = real * vec3
//    vec3 = real * direction
//    vec3 = real * station
Placement RealPlacementRep::genericMul(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        RealPlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcValue()*rp.getRep().calcValue())
                : (RealPlacementRep*)RealExprPlacementRep::mulOp(getMyHandle(),rp);
        return Placement(result);
    }

    if (Vec3Placement::isInstanceOf(r))
        return Vec3Placement::downcast(r) * getMyHandle();   // punt to vec3*real

    if (DirectionPlacement::isInstanceOf(r) || StationPlacement::isInstanceOf(r))
        return r.getRep().castToVec3Placement() * getMyHandle(); // punt to vec3(r) * real

    return PlacementRep::genericMul(r);    // die
}

// result = real / placement
// We support only
//    real = real / real
Placement RealPlacementRep::genericDvd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        RealPlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcValue()/rp.getRep().calcValue())
                : (RealPlacementRep*)RealExprPlacementRep::dvdOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericDvd(r);    // die
}

    // REAL FEATURE PLACEMENT REP //
const Real& RealFeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();
    if (!isIndexed()) 
        return RealPlacementRep::downcast(p).getValue(/*State*/);

    // indexed
    if (Vec3PlacementRep::isA(p))
        return Vec3PlacementRep::downcast(p).getValue(/*State*/)
                [getPlacementIndex()];
    if (StationPlacementRep::isA(p))
        return StationPlacementRep::downcast(p).getValue(/*State*/)
                [getPlacementIndex()];
    if (DirectionPlacementRep::isA(p))
        return DirectionPlacementRep::downcast(p).getValue(/*State*/)
                [getPlacementIndex()];

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Real*>(0);
}

    // REAL EXPR PLACEMENT REP //

bool 
RealOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    // unary
    case Negate: case Abs: case Sqrt: case Exp:
    case Sin:    case Cos: case Asin: case Acos:
        return args.size()==1 && RealPlacement::isInstanceOf(args[0]); 
    case VectorLength:
        return args.size()==1 && Vec3Placement::isInstanceOf(args[0]);

    // binary
    case Add:      case Subtract:
    case Multiply: case Divide:
        return args.size()==2 && RealPlacement::isInstanceOf(args[0])
                              && RealPlacement::isInstanceOf(args[1]);
    //case DotProduct2:
    //   return args.size()==2 && Vec2Placement::isInstanceOf(args[0])
    //                         && Vec2Placement::isInstanceOf(args[1]);  
    case DotProduct3:
        return args.size()==2 && Vec3Placement::isInstanceOf(args[0])
                              && Vec3Placement::isInstanceOf(args[1]);  
    case PointDistance:
        return args.size()==2 && StationPlacement::isInstanceOf(args[0])
                              && StationPlacement::isInstanceOf(args[1]); 

    case AngleBetweenDirections:
        return args.size()==2 && DirectionPlacement::isInstanceOf(args[0])
                              && DirectionPlacement::isInstanceOf(args[1]); 

    default:
        assert(false);
    }
    //NOTREACHED
    return false;
}

Real RealOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    Real val = NTraits<Real>::getNaN();
    Real arg1, arg2;
    if (args.size() > 0 && RealPlacement::isInstanceOf(args[0]))
        arg1 = RealPlacement::downcast(args[0]).getRep().calcValue();
    if (args.size() > 1 && RealPlacement::isInstanceOf(args[1]))
        arg2 = RealPlacement::downcast(args[1]).getRep().calcValue();

    switch (op) {
    case Negate: val = -arg1; break;
    case Abs:    val = std::abs(arg1); break;
    case Sqrt:   val = std::sqrt(arg1); break;
    case Exp:    val = std::exp(arg1); break;
    case Log:    val = std::log(arg1); break;
    case Sin:    val = std::sin(arg1); break;
    case Cos:    val = std::cos(arg1); break;
    case Asin:   val = std::asin(arg1); break;
    case Acos:   val = std::acos(arg1); break;
    case VectorLength:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue().norm(); 
        break;

    case Add:      val = arg1+arg2; break;
    case Subtract: val = arg1-arg2; break;
    case Multiply: val = arg1*arg2; break;
    case Divide:   val = arg1/arg2; break; 

    case DotProduct3: {
        const Vec3 l = Vec3Placement::downcast(args[0]).getRep().calcValue();
        const Vec3 r = Vec3Placement::downcast(args[1]).getRep().calcValue();
        val = (~l)*r;
        break;
    }
    case DotProduct2: {
        assert(false);
        break;
    }
    case PointDistance: {
        const Vec3 head = StationPlacement::downcast(args[0]).getRep().calcValue();
        const Vec3 tail = StationPlacement::downcast(args[1]).getRep().calcValue();
        val = (head-tail).norm();
        break;
    }
    case AngleBetweenDirections: {
        const Vec3 v1 = DirectionPlacement::downcast(args[0]).getRep().calcValue();
        const Vec3 v2 = DirectionPlacement::downcast(args[1]).getRep().calcValue();
        Real dotprod = ~v1*v2;
        if (dotprod < -1.) dotprod = -1.;   // watch for roundoff
        else if (dotprod > 1.) dotprod = 1.;
        val = std::acos(dotprod);
        break;
    }
    default: assert(false);
    }
    return val;
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::unaryOp(RealOps::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new RealExprPlacementRep(RealOps(op), args);
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::binaryOp(RealOps::OpKind op, 
                               const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new RealExprPlacementRep(RealOps(op), args);
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::negateOp(const RealPlacement& rp)
  { return unaryOp(RealOps::Negate, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::absOp   (const RealPlacement& rp)
  { return unaryOp(RealOps::Abs, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::sqrtOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Sqrt, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::expOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Exp, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::logOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Log, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::sinOp   (const RealPlacement& rp)
  { return unaryOp(RealOps::Sin, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::cosOp   (const RealPlacement& rp)
  { return unaryOp(RealOps::Cos, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::asinOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Asin, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::acosOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Acos, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::lengthOp(const Vec3Placement& vp)
  { return unaryOp(RealOps::VectorLength, vp); }

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::addOp(const RealPlacement& l, const RealPlacement& r)
  { return binaryOp(RealOps::Add, l, r); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::subOp(const RealPlacement& l, const RealPlacement& r)
  { return binaryOp(RealOps::Subtract, l, r); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::mulOp(const RealPlacement& l, const RealPlacement& r)
  { return binaryOp(RealOps::Multiply, l, r); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::dvdOp(const RealPlacement& l, const RealPlacement& r)
  { return binaryOp(RealOps::Divide, l, r); }

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::distanceOp(const StationPlacement& l, const StationPlacement& r)
  { return binaryOp(RealOps::PointDistance, l, r); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::dot3Op    (const Vec3Placement& l,    const Vec3Placement& r)
  { return binaryOp(RealOps::DotProduct3, l, r); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::angleOp(const DirectionPlacement& l, const DirectionPlacement& r)
  { return binaryOp(RealOps::AngleBetweenDirections, l, r); }

    // VEC3 PLACEMENT REP //

DirectionPlacement Vec3PlacementRep::castToDirectionPlacement() const {
    DirectionPlacementRep* result = 
        isConstant() ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(calcValue())
                     : (DirectionPlacementRep*)DirectionExprPlacementRep::normalizeOp(getMyHandle()); 
    return DirectionPlacement(result);
}

StationPlacement Vec3PlacementRep::castToStationPlacement() const {
    StationPlacementRep* result = 
        isConstant() ? (StationPlacementRep*)new StationConstantPlacementRep(calcValue())
                     : (StationPlacementRep*)StationExprPlacementRep::recastVec3Op(getMyHandle()); 
    return StationPlacement(result);
}

// result = -vec3 (result is always vec3)
Placement Vec3PlacementRep::genericNegate() const {
    Vec3PlacementRep* result =
        isConstant() ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(-calcValue())
                     : (Vec3PlacementRep*)Vec3ExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = length(vec3) (result is always real)
Placement Vec3PlacementRep::genericLength() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(calcValue().norm())
                     : (RealPlacementRep*)RealExprPlacementRep::lengthOp(getMyHandle());
    return Placement(result);
}

// result = normalize(vec3) (result is always Direction)
Placement Vec3PlacementRep::genericNormalize() const {
    DirectionPlacementRep* result = 
        isConstant() ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(calcValue())
                     : (DirectionPlacementRep*)DirectionExprPlacementRep::normalizeOp(getMyHandle());
    return Placement(result);
}

// result = vec3 + placement
//   We support
//     vec3 = vec3 + vec3
//     station = vec3 + station
Placement Vec3PlacementRep::genericAdd(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcValue()+vp.getRep().calcValue())
                : (Vec3PlacementRep*)Vec3ExprPlacementRep::addOp(getMyHandle(),vp);
        return Placement(result);
    }

    // Let StationPlacementRep deal with v+s; reverse the args.
    if (StationPlacement::isInstanceOf(r))
        return StationPlacement::downcast(r) + getMyHandle();

    return PlacementRep::genericAdd(r);    // die
}

// result = vec3 - placement
//   We support only
//     vec3 = vec3 - vec3
Placement Vec3PlacementRep::genericSub(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcValue()-vp.getRep().calcValue())
                : (Vec3PlacementRep*)Vec3ExprPlacementRep::subOp(getMyHandle(),vp);
        return Placement(result);
    }

    return PlacementRep::genericSub(r);    // die
}

// result = vec3 * placement
//   We support only
//     vec3 = vec3 * real
// TODO: should we allow real = vec3*vec3 to be a dot product?
Placement Vec3PlacementRep::genericMul(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcValue()*rp.getRep().calcValue())
                : (Vec3PlacementRep*)Vec3ExprPlacementRep::smulOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericMul(r);    // die
}

// result = vec3 / placement
//   We support only
//     vec3 = vec3 / real
Placement Vec3PlacementRep::genericDvd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcValue()/rp.getRep().calcValue())
                : (Vec3PlacementRep*)Vec3ExprPlacementRep::sdvdOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericDvd(r);    // die
}

// result = dot(vec3, placement)
//   We support 
//     real = dot(vec3, vec3)
//     real = dot(vec3, direction)
Placement Vec3PlacementRep::genericDotProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        RealPlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (RealPlacementRep*)new RealConstantPlacementRep(~calcValue()*vp.getRep().calcValue())
                : (RealPlacementRep*)RealExprPlacementRep::dot3Op(getMyHandle(),vp);
        return Placement(result);
    }

    if (DirectionPlacement::isInstanceOf(r)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(r);
        return genericDotProduct(dp.getRep().castToVec3Placement()); // recursive call with improved argument
    }

    return PlacementRep::genericDotProduct(r);    // die
}

// result = cross(vec3, placement)
//   We support
//     vec3 = cross(vec3, vec3)
//     vec3 = cross(vec3, station)
//     vec3 = cross(vec3, direction)
Placement Vec3PlacementRep::genericCrossProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(cross(calcValue(),vp.getRep().calcValue()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),vp);
        return Placement(result);
    }

    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && sp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(cross(calcValue(),sp.getRep().calcValue()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),
                                                 sp.getRep().castToVec3Placement());
        return Placement(result);
    }

    if (DirectionPlacement::isInstanceOf(r)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && dp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(cross(calcValue(),dp.getRep().calcValue()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),
                                                 dp.getRep().castToVec3Placement());
        return Placement(result);
    }

    return PlacementRep::genericCrossProduct(r);    // die
}

// result = angle(vec3, placement)
//   We support 
//     real = angle(vec3, vec3)
//     real = angle(vec3, station)
//     real = angle(vec3, direction)
Placement Vec3PlacementRep::genericAngle(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        return angle(normalize(getMyHandle()), normalize(vp));
    }
    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        return angle(normalize(getMyHandle()), normalize(sp));
    }
    if (DirectionPlacement::isInstanceOf(r)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(r);
        return angle(normalize(getMyHandle()), dp);
    }

    return PlacementRep::genericAngle(r);    // die
}

    // VEC3 FEATURE PLACEMENT REP //
const Vec3& Vec3FeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();
    if (!isIndexed())
        return Vec3PlacementRep::downcast(p).getValue(/*State*/);

    // indexed
    if (OrientationPlacementRep::isA(p))
        return OrientationPlacementRep::downcast(p).getValue(/*State*/)
                (getPlacementIndex()); // round () to get column

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Vec3*>(0);
}

    // VEC3 EXPR PLACEMENT REP //
bool Vec3Ops::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    case RecastStation:
        return args.size()==1 && StationPlacement::isInstanceOf(args[0]);

    case RecastDirection:
        return args.size()==1 && DirectionPlacement::isInstanceOf(args[0]);

    case Negate:
        return args.size()==1 && Vec3Placement::isInstanceOf(args[0]);

    // v= v+v, but s+s not allowed
    case Add:
        return args.size()==2 && Vec3Placement::isInstanceOf(args[0])
                              && Vec3Placement::isInstanceOf(args[1]);
    
    // v=v-v
    case Subtract:
        return args.size()==2 && Vec3Placement::isInstanceOf(args[0])
                              && Vec3Placement::isInstanceOf(args[1]);

    // v=s-s
    case StationDifference:
        return args.size()==2 && StationPlacement::isInstanceOf(args[0])
                              && StationPlacement::isInstanceOf(args[1]);

    case ScalarMultiply:
    case ScalarDivide:
        return args.size()==2 && Vec3Placement::isInstanceOf(args[0])
                              && RealPlacement::isInstanceOf(args[1]);

    case CrossProduct:
        return args.size()==2 && Vec3Placement::isInstanceOf(args[0])
                              && Vec3Placement::isInstanceOf(args[1]);
    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

Vec3 Vec3Ops::apply(/*State,*/ const std::vector<Placement>& args) const {
    Vec3 val;
    switch(op) {
    case RecastStation:
        val = StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    case RecastDirection:
        val = DirectionPlacement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    case Negate:
        val = -Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    case Add:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/)
              + Vec3Placement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;
    
    case Subtract:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/)
              - Vec3Placement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    case StationDifference:
        val = StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/)
              - StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    // real is always on the right for scalar mul & dvd
    case ScalarMultiply:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/)
              * RealPlacement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    case ScalarDivide:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/)
              / RealPlacement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    case CrossProduct:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/)
              % Vec3Placement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    default: 
        assert(false);
    }
    return val;
}

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::unaryOp(Vec3Ops::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new Vec3ExprPlacementRep(Vec3Ops(op), args);
}

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::binaryOp(Vec3Ops::OpKind op, 
                               const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new Vec3ExprPlacementRep(Vec3Ops(op), args);
}

// Supported Vec3Expr-building operators

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::negateOp(const Vec3Placement& vp)
  { return unaryOp(Vec3Ops::Negate, vp); }
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::recastStationOp(const StationPlacement& sp)
  { return unaryOp(Vec3Ops::RecastStation, sp); }
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::recastDirectionOp(const DirectionPlacement& dp)
  { return unaryOp(Vec3Ops::RecastDirection, dp); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::addOp (const Vec3Placement& l, const Vec3Placement& r)
  { return binaryOp(Vec3Ops::Add, l, r); }
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::subOp (const Vec3Placement& l, const Vec3Placement& r)
  { return binaryOp(Vec3Ops::Subtract, l, r); }
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::stationSubOp (const StationPlacement& head, 
                                    const StationPlacement& tail)
  { return binaryOp(Vec3Ops::StationDifference, head, tail); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::smulOp(const Vec3Placement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarMultiply, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::smulOp(const StationPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarMultiply, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::smulOp(const DirectionPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarMultiply, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const Vec3Placement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const StationPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const DirectionPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::crossOp(const Vec3Placement& l, const Vec3Placement& r)
  { return binaryOp(Vec3Ops::CrossProduct, l, r); }

    // STATION PLACEMENT REP //

Vec3Placement StationPlacementRep::castToVec3Placement() const {
    Vec3PlacementRep* result = 
        isConstant()
          ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcValue())
          : (Vec3PlacementRep*)Vec3ExprPlacementRep::recastStationOp(getMyHandle()); 
    return Vec3Placement(result);
}

// result = -station (result is always vec3)
Placement StationPlacementRep::genericNegate() const {
    return -castToVec3Placement();
}

// result = length(station) (result is always real)
Placement StationPlacementRep::genericLength() const {
    return length(castToVec3Placement());
}

// result = normalize(station) (result is always Direction)
Placement StationPlacementRep::genericNormalize() const {
    return normalize(castToVec3Placement());
}

// result = station + placement
//   We support
//     station = station + vec3
Placement StationPlacementRep::genericAdd(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        StationPlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (StationPlacementRep*)new StationConstantPlacementRep
                                              (calcValue()+vp.getRep().calcValue())
                : (StationPlacementRep*)StationExprPlacementRep::addOp(getMyHandle(),vp);
        return Placement(result);
    }

    return PlacementRep::genericAdd(r);    // die
}

// result = station - placement
//   We support
//     station = station - vec3
//     vec3    = station - station
Placement StationPlacementRep::genericSub(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        StationPlacementRep* result = 
            isConstant() && vp.getRep().isConstant()
                ? (StationPlacementRep*)new StationConstantPlacementRep
                                              (calcValue()-vp.getRep().calcValue())
                : (StationPlacementRep*)StationExprPlacementRep::subOp(getMyHandle(),vp);
        return Placement(result);
    }

    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && sp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep
                                              (calcValue()-sp.getRep().calcValue())
                : (Vec3PlacementRep*)Vec3ExprPlacementRep::stationSubOp(getMyHandle(),sp);
        return Placement(result);
    }

    return PlacementRep::genericSub(r);    // die
}

// result = station * placement
//   We support only
//     vec3 = station * real
// TODO: should we allow real = station*vec3, station*direction to be a dot product?
Placement StationPlacementRep::genericMul(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        return castToVec3Placement() * r;
    }

    return PlacementRep::genericMul(r);    // die
}

// result = station / placement
//   We support only
//     vec3 = station / real
Placement StationPlacementRep::genericDvd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        return castToVec3Placement() / r;
    }

    return PlacementRep::genericDvd(r);    // die
}

// result = dot(station, placement)
//   We support 
//     real = dot(station, vec3)
//     real = dot(station, direction)
Placement StationPlacementRep::genericDotProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r) || DirectionPlacement::isInstanceOf(r)) {
        return dot(castToVec3Placement(), r);
    }

    return PlacementRep::genericDotProduct(r);    // die
}


// result = cross(station, placement)
//   We support 
//     vec3 = cross(station, vec3)
//     vec3 = cross(station, direction)
//     vec3 = cross(station, station)
Placement StationPlacementRep::genericCrossProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r) 
        || DirectionPlacement::isInstanceOf(r)
        || StationPlacement::isInstanceOf(r))
    {
        return cross(castToVec3Placement(), r);
    }

    return PlacementRep::genericCrossProduct(r);    // die
}

// result = distance(station, placement)
//   We support only
//     real = distance(station, station)
Placement StationPlacementRep::genericDistance(const Placement& r) const {
    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        return length(getMyHandle()-sp);
    }

    return PlacementRep::genericDistance(r);    // die
}

// result = angle(station, placement)
//   We support 
//     real = angle(station, vec3)
//     real = angle(station, station)
//     real = angle(station, direction)
Placement StationPlacementRep::genericAngle(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        return angle(normalize(getMyHandle()), normalize(vp));
    }
    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        return angle(normalize(getMyHandle()), normalize(sp));
    }
    if (DirectionPlacement::isInstanceOf(r)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(r);
        return angle(normalize(getMyHandle()), dp);
    }

    return PlacementRep::genericAngle(r);    // die
}

    // STATION FEATURE PLACEMENT REP //
const Vec3& StationFeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();
    if (!isIndexed())
        return StationPlacementRep::downcast(p).getValue(/*State*/);

    // indexed
    if (FramePlacementRep::isA(p) && getPlacementIndex()==1)
        return FramePlacementRep::downcast(p).getOriginValue(/*State*/);

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Vec3*>(0);
}

    // STATION EXPR PLACEMENT REP //

bool 
StationOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    case RecastVec3:
        return args.size()==1 && Vec3Placement::isInstanceOf(args[0]);

    // s=s+v, s=s-v, BUT v=s-s and s+s is not meaningful without cast
    case Add:
    case Subtract:
        return args.size()==2 && StationPlacement::isInstanceOf(args[0])
                              && Vec3Placement::isInstanceOf(args[1]);

    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

Vec3 StationOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    Vec3 val;
    switch (op) {
    case RecastVec3:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    case Add:
        val = StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/)
              + Vec3Placement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    case Subtract:
        val = StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/)
              - Vec3Placement::downcast(args[1]).getRep().calcValue(/*State*/);
        break;

    default: 
        assert(false);
    //NOTREACHED
    }
    return val;
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::unaryOp(StationOps::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new StationExprPlacementRep(StationOps(op), args);
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::binaryOp(StationOps::OpKind op, 
                               const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new StationExprPlacementRep(StationOps(op), args);
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::recastVec3Op(const Vec3Placement& p)
  { return unaryOp(StationOps::RecastVec3, p); }
/*static*/ StationExprPlacementRep*
StationExprPlacementRep::addOp(const StationPlacement& l, const Vec3Placement& r)
  { return binaryOp(StationOps::Add, l, r); }
/*static*/ StationExprPlacementRep*
StationExprPlacementRep::subOp(const StationPlacement& l, const Vec3Placement& r)
  { return binaryOp(StationOps::Subtract, l, r); }


    // DIRECTION PLACEMENT REP //


Vec3Placement DirectionPlacementRep::castToVec3Placement() const {
    Vec3PlacementRep* result = 
        isConstant()
          ? (Vec3PlacementRep*)
             new Vec3ConstantPlacementRep(calcValue())
          : (Vec3PlacementRep*)
             Vec3ExprPlacementRep::recastDirectionOp(getMyHandle()); 
    return Vec3Placement(result);
}

// result = -direction (result is always direction)
Placement DirectionPlacementRep::genericNegate() const {
    DirectionPlacementRep* result =
        isConstant() ? (DirectionPlacementRep*)
                        new DirectionConstantPlacementRep(-calcValue())
                     : (DirectionPlacementRep*)
                        DirectionExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = direction * placement
//   We support only
//     vec3 = direction * real
// TODO: should we allow real = direction*direction, direction*vec3 to be a dot product?
Placement DirectionPlacementRep::genericMul(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(calcValue()*rp.getRep().calcValue())
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::smulOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericMul(r);    // die
}

// result = direction / placement
//   We support only
//     vec3 = direction / real
Placement DirectionPlacementRep::genericDvd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        const RealPlacement& rp = RealPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && rp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(calcValue()/rp.getRep().calcValue())
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::sdvdOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericDvd(r);    // die
}

// result = dot(direction, placement)
//   We support 
//     real = dot(direction, vec3)
//     real = dot(direction, direction)
Placement DirectionPlacementRep::genericDotProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r) || DirectionPlacement::isInstanceOf(r)) {
        return dot(castToVec3Placement(), r);
    }

    return PlacementRep::genericDotProduct(r);    // die
}

// result = cross(direction, placement)
//   We support 
//     vec3 = cross(direction, vec3)
//     vec3 = cross(direction, direction)
//     vec3 = cross(direction, station)
Placement DirectionPlacementRep::genericCrossProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r) 
        || DirectionPlacement::isInstanceOf(r)
        || StationPlacement::isInstanceOf(r))
    {
        return cross(castToVec3Placement(), r);
    }

    return PlacementRep::genericCrossProduct(r);    // die
}

// result = angle(direction, placement)
//   We support 
//     real = angle(direction, vec3)
//     real = angle(direction, station)
//     real = angle(direction, direction)
// The last one invokes the actual real operator since it requires
// no normalization.
Placement DirectionPlacementRep::genericAngle(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r)) {
        const Vec3Placement& vp = Vec3Placement::downcast(r);
        return angle(getMyHandle(), normalize(vp));
    }
    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        return angle(getMyHandle(), normalize(sp));
    }

    if (DirectionPlacement::isInstanceOf(r)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(r);
        if (isConstant() && dp.getRep().isConstant()) {
            Real dotprod = ~calcValue() * dp.getRep().calcValue();
            if (dotprod < -1.) dotprod = -1.;
            else if (dotprod > 1.) dotprod = 1.;
            return RealPlacement(new RealConstantPlacementRep(std::acos(dotprod)));
        }
        return RealPlacement(RealExprPlacementRep::angleOp(getMyHandle(),dp));
    }

    return PlacementRep::genericAngle(r);    // die
}

    // DIRECTION FEATURE PLACEMENT REP //
const Vec3& DirectionFeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();

    if (!isIndexed())
        return DirectionPlacementRep::downcast(p).getValue(/*State*/);

    // indexed
    if (OrientationPlacementRep::isA(p))
        return OrientationPlacementRep::downcast(p).getValue(/*State*/)
            (getPlacementIndex()); // round brackets () to get column not row

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Vec3*>(0);
}

    // DIRECTION EXPR PLACEMENT REP //

bool DirectionOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {

    case Negate:
        return args.size()==1 && DirectionPlacement::isInstanceOf(args[0]);

    case NormalizeVec3:
        return args.size()==1 && Vec3Placement::isInstanceOf(args[0]);

    case NormalizeStation:
        return args.size()==1 && StationPlacement::isInstanceOf(args[0]);

    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

Vec3 DirectionOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    Vec3 val;
    switch (op) {

    case Negate:
        val = DirectionPlacement::downcast(args[0]).getRep().calcValue(/*State*/);
        break;

    case NormalizeVec3:
        val = Vec3Placement::downcast(args[0]).getRep().calcValue(/*State*/);
        val /= val.norm();
        break;

    case NormalizeStation:
        val = StationPlacement::downcast(args[0]).getRep().calcValue(/*State*/);
        val /= val.norm();
        break;

    default: 
        assert(false);
    //NOTREACHED
    }

    return val;
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::unaryOp(DirectionOps::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new DirectionExprPlacementRep(DirectionOps(op), args);
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::binaryOp(DirectionOps::OpKind op, 
                                    const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new DirectionExprPlacementRep(DirectionOps(op), args);
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::negateOp(const DirectionPlacement& p)
  { return unaryOp(DirectionOps::Negate, p); }
/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalizeOp(const StationPlacement& p)
  { return unaryOp(DirectionOps::NormalizeStation, p); }
/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalizeOp(const Vec3Placement& p)
  { return unaryOp(DirectionOps::NormalizeVec3, p); }

    // ORIENTATION PLACEMENT REP //

    // ORIENTATION FEATURE PLACEMENT REP //
const Mat33& OrientationFeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();

    if (!isIndexed())
        return OrientationPlacementRep::downcast(p).getValue(/*State*/);

    // indexed
    if (FramePlacementRep::isA(p) && getPlacementIndex()==0)
        return FramePlacementRep::downcast(p).getOrientationValue(/*State*/);

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Mat33*>(0);
}

    // ORIENTATION EXPR PLACEMENT REP //

bool OrientationOps::checkArgs(const std::vector<Placement>& args) const {
    assert(false); // none yet
    return false;
}

Mat33 OrientationOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    assert(false);
    return Mat33();
}

/*static*/ OrientationExprPlacementRep*
OrientationExprPlacementRep::unaryOp(OrientationOps::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new OrientationExprPlacementRep(OrientationOps(op), args);
}

/*static*/ OrientationExprPlacementRep*
OrientationExprPlacementRep::binaryOp(OrientationOps::OpKind op, 
                               const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new OrientationExprPlacementRep(OrientationOps(op), args);
}

    // FRAME PLACEMENT REP //


    // FRAME FEATURE PLACEMENT REP //
const Mat34& FrameFeaturePlacementRep::getReferencedValue(/*State*/) const {
    const PlacementRep& p = getReferencedFeature().getPlacement().getRep();

    if (!isIndexed())
        return FramePlacementRep::downcast(p).getValue(/*State*/);

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Mat34*>(0);
}

    // FRAME EXPR PLACEMENT REP //

bool FrameExprPlacementRep::isLimitedToSubtree
    (const Subsystem& root, const Feature*& offender) const
{
    if (!orientation.getRep().isLimitedToSubtree(root,offender))
        return false;
    return origin.getRep().isLimitedToSubtree(root,offender);
}

void FrameExprPlacementRep::repairFeatureReferences
    (const Subsystem& oldRoot, const Subsystem& newRoot)
{
    orientation.updRep().repairFeatureReferences(oldRoot,newRoot);
    origin.updRep().repairFeatureReferences(oldRoot,newRoot);
}


const Subsystem*
FrameExprPlacementRep::findAncestorSubsystem(const Subsystem& youngestAllowed) const {
    const Subsystem* ancestor = &youngestAllowed;
    ancestor = orientation.getRep().findAncestorSubsystem(*ancestor);

    if (ancestor)
        ancestor = origin.getRep().findAncestorSubsystem(*ancestor);

    return ancestor;
}


const Subsystem*
FrameExprPlacementRep::findPlacementValueOwnerSubsystem(const Subsystem& youngestAllowed) const {
    const Subsystem* valueOwner = &youngestAllowed;
    valueOwner = orientation.getRep().findPlacementValueOwnerSubsystem(*valueOwner);

    if (valueOwner)
        valueOwner = origin.getRep().findPlacementValueOwnerSubsystem(*valueOwner);

    return valueOwner;
}

    // FRAME EXPR PLACEMENT REP //

bool FrameOps::checkArgs(const std::vector<Placement>& args) const {
    assert(false); // none yet
    return false;
}

Mat34 FrameOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    assert(false);
    return Mat34();
}

} // namespace simtk
