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
using std::endl;
using std::ostream;

namespace simtk {

    // PLACEMENT EXPR //

const Feature* 
PlacementExpr::exprFindAncestorFeature(const Feature& root) const {
    const Feature* ancestor = 0;
    bool foundNonConst = false;
    for (size_t i=0; i < args.size(); ++i) {
        if (args[i].isConstant())
            continue;
        foundNonConst = true;
        const Feature* argAncestor = 
            args[i].getRep().findAncestorFeature(root);
        if (ancestor && argAncestor)
            ancestor = FeatureRep::findYoungestCommonAncestor(*ancestor,*argAncestor);
        else ancestor = argAncestor;
    }
    assert(foundNonConst);
    return ancestor;
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
    (const Feature& root, const Feature*& offender) const
{
    for (size_t i=0; i<args.size(); ++i)
        if (!args[i].getRep().isLimitedToSubtree(root,offender))
            return false;
    return true;
}

void PlacementExpr::exprRepairFeatureReferences
    (const Feature& oldRoot, const Feature& newRoot)
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

const Feature* 
FeatureReference::refFindAncestorFeature(const Feature& root) const {
    assert(feature);
    return FeatureRep::isFeatureInFeatureTree(root, *feature) 
            ? feature : 0;
}


// Check that this feature is on the feature subtree rooted by "ancestor". If
// not return a pointer to this feature for use in a friendly error message.
// If this is the right tree, we return true with offender==NULL.
bool FeatureReference::refIsLimitedToSubtree
    (const Feature& root, const Feature*& offender) const
{
    assert(feature);
    if (FeatureRep::isFeatureInFeatureTree(root, *feature)) {
        offender = 0;
        return true;
    }
    offender = feature;
    return false;
}

void FeatureReference::refRepairFeatureReferences
    (const Feature& oldRoot, const Feature& newRoot)
{
    assert(feature);
    const Feature* corrFeature = FeatureRep::findCorrespondingFeature(oldRoot,*feature,newRoot);
    assert(corrFeature);
    feature = corrFeature;
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
void PlacementRep::repairValueReference(const Feature& oldRoot, const Feature& newRoot) {
    if (valueSlot)
        valueSlot = FeatureRep::findCorrespondingPlacementValue
                        (oldRoot,*valueSlot,newRoot);
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

    // REAL PLACEMENT REP //

// result = -real (result is always real)
Placement RealPlacementRep::genericNegate() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(-getValue())
                     : (RealPlacementRep*)RealExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = abs(real) (result is always real)
Placement RealPlacementRep::genericAbs() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::abs(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::absOp(getMyHandle());
    return Placement(result);
}

// result = sqrt(real) (result is always real)
Placement RealPlacementRep::genericSqrt() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sqrt(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sqrtOp(getMyHandle());
    return Placement(result);
}

// result = exp(real) (result is always real)
Placement RealPlacementRep::genericExp() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::exp(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::expOp(getMyHandle());
    return Placement(result);
}

// result = log(real) (natural log, result is always real)
Placement RealPlacementRep::genericLog() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::log(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::logOp(getMyHandle());
    return Placement(result);
}

// result = sin(real) (result is always real)
Placement RealPlacementRep::genericSin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sin(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sinOp(getMyHandle());
    return Placement(result);
}

// result = cos(real) (result is always real)
Placement RealPlacementRep::genericCos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::cos(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::cosOp(getMyHandle());
    return Placement(result);
}

// result = asin(real) (result is always real)
Placement RealPlacementRep::genericAsin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::asin(getValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::asinOp(getMyHandle());
    return Placement(result);
}

// result = acos(real) (result is always real)
Placement RealPlacementRep::genericAcos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::acos(getValue()))
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(getValue()+rp.getRep().getValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(getValue()-rp.getRep().getValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(getValue()*rp.getRep().getValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(getValue()/rp.getRep().getValue())
                : (RealPlacementRep*)RealExprPlacementRep::dvdOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericDvd(r);    // die
}

    // REAL FEATURE PLACEMENT REP //
Real RealFeaturePlacementRep::getValue(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Real value = NTraits<Real>::getNaN();
    if (!isIndexed()) 
        value = RealPlacementRep::downcast(p).getValue(/*State*/);
    else if (Vec3PlacementRep::isA(p))
        value = Vec3PlacementRep::downcast(p).getValue(/*State*/)
                [getPlacementIndex()];
    else if (StationPlacementRep::isA(p))
        value = StationPlacementRep::downcast(p).getMeasureNumbers(/*State*/)
                [getPlacementIndex()];
    else if (DirectionPlacementRep::isA(p))
        value = DirectionPlacementRep::downcast(p).getMeasureNumbers(/*State*/)
                [getPlacementIndex()];
    else
        assert(false);

    return value;
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
    default:
        assert(false);
    }
    //NOTREACHED
    return false;
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

    // VEC3 PLACEMENT REP //

DirectionPlacement Vec3PlacementRep::castToDirectionPlacement() const {
    DirectionPlacementRep* result = 
        isConstant() ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(getValue())
                     : (DirectionPlacementRep*)DirectionExprPlacementRep::normalizeOp(getMyHandle()); 
    return DirectionPlacement(result);
}

StationPlacement Vec3PlacementRep::castToStationPlacement() const {
    StationPlacementRep* result = 
        isConstant() ? (StationPlacementRep*)new StationConstantPlacementRep(getValue())
                     : (StationPlacementRep*)StationExprPlacementRep::recastVec3Op(getMyHandle()); 
    return StationPlacement(result);
}

// result = -vec3 (result is always vec3)
Placement Vec3PlacementRep::genericNegate() const {
    Vec3PlacementRep* result =
        isConstant() ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(-getValue())
                     : (Vec3PlacementRep*)Vec3ExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = length(vec3) (result is always real)
Placement Vec3PlacementRep::genericLength() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(getValue().norm())
                     : (RealPlacementRep*)RealExprPlacementRep::lengthOp(getMyHandle());
    return Placement(result);
}

// result = normalize(vec3) (result is always Direction)
Placement Vec3PlacementRep::genericNormalize() const {
    DirectionPlacementRep* result = 
        isConstant() ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(getValue())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(getValue()+vp.getRep().getValue())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(getValue()-vp.getRep().getValue())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(getValue()*rp.getRep().getValue())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(getValue()/rp.getRep().getValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(~getValue()*vp.getRep().getValue())
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
                   new Vec3ConstantPlacementRep(cross(getValue(),vp.getRep().getValue()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),vp);
        return Placement(result);
    }

    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && sp.getRep().isConstant()
                ? (Vec3PlacementRep*)
                   new Vec3ConstantPlacementRep(cross(getValue(),sp.getRep().getMeasureNumbers()))
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
                   new Vec3ConstantPlacementRep(cross(getValue(),dp.getRep().getMeasureNumbers()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),
                                                 dp.getRep().castToVec3Placement());
        return Placement(result);
    }

    return PlacementRep::genericCrossProduct(r);    // die
}

    // VEC3 FEATURE PLACEMENT REP //
Vec3 Vec3FeaturePlacementRep::getValue(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Vec3 value;
    if (!isIndexed())
        value = Vec3PlacementRep::downcast(p).getValue(/*State*/);
    else if (OrientationPlacementRep::isA(p))
        value = OrientationPlacementRep::downcast(p).getMeasureNumbers(/*State*/)
                (getPlacementIndex()); // round () to get column
    else
        assert(false);

    return value;
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
          ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(getMeasureNumbers())
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
                                              (getMeasureNumbers()+vp.getRep().getValue())
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
                                              (getMeasureNumbers()-vp.getRep().getValue())
                : (StationPlacementRep*)StationExprPlacementRep::subOp(getMyHandle(),vp);
        return Placement(result);
    }

    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && sp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep
                                              (getMeasureNumbers()-sp.getRep().getMeasureNumbers())
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

    // STATION FEATURE PLACEMENT REP //
Vec3 StationFeaturePlacementRep::getMeasureNumbers(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Vec3 value = Vec3(NTraits<Real>::getNaN());
    if (!isIndexed())
        value = StationPlacementRep::downcast(p).getMeasureNumbers(/*State*/);
    else if (FramePlacementRep::isA(p) && getPlacementIndex()==1)
        value = FramePlacementRep::downcast(p).getOriginMeasureNumbers(/*State*/);
    else
        assert(false);

    return value;
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
             new Vec3ConstantPlacementRep(getMeasureNumbers())
          : (Vec3PlacementRep*)
             Vec3ExprPlacementRep::recastDirectionOp(getMyHandle()); 
    return Vec3Placement(result);
}

// result = -direction (result is always direction)
Placement DirectionPlacementRep::genericNegate() const {
    DirectionPlacementRep* result =
        isConstant() ? (DirectionPlacementRep*)
                        new DirectionConstantPlacementRep(-getMeasureNumbers())
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
                   new Vec3ConstantPlacementRep(getMeasureNumbers()*rp.getRep().getValue())
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
                   new Vec3ConstantPlacementRep(getMeasureNumbers()/rp.getRep().getValue())
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

    // DIRECTION FEATURE PLACEMENT REP //
Vec3 DirectionFeaturePlacementRep::getMeasureNumbers(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Vec3 value = Vec3(NTraits<Real>::getNaN());
    if (!isIndexed())
        value = DirectionPlacementRep::downcast(p).getMeasureNumbers(/*State*/);
    else if (OrientationPlacementRep::isA(p))
        value = OrientationPlacementRep::downcast(p).getMeasureNumbers(/*State*/)
            (getPlacementIndex()); // round brackets () to get column not row
    else
        assert(false);

    return value;
}

    // DIRECTION EXPR PLACEMENT REP //

bool DirectionOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {

    case Normalize:
        return args.size()==1 && (   Vec3Placement::isInstanceOf(args[0])
                                  || StationPlacement::isInstanceOf(args[0]));

    default: 
        assert(false);
        //NOTREACHED
    }
    //NOTREACHED
    return false;
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
  { return unaryOp(DirectionOps::Normalize, p); }
/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalizeOp(const Vec3Placement& p)
  { return unaryOp(DirectionOps::Normalize, p); }

    // ORIENTATION PLACEMENT REP //

    // ORIENTATION FEATURE PLACEMENT REP //
Mat33 OrientationFeaturePlacementRep::getMeasureNumbers(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Mat33 value = Mat33(NTraits<Real>::getNaN());
    if (!isIndexed())
        value = OrientationPlacementRep::downcast(p).getMeasureNumbers(/*State*/);
    else if (FramePlacementRep::isA(p) && getPlacementIndex()==0)
        value = FramePlacementRep::downcast(p).getOrientationMeasureNumbers(/*State*/);
    else
        assert(false);

    return value;
}

    // ORIENTATION EXPR PLACEMENT REP //

bool OrientationOps::checkArgs(const std::vector<Placement>& args) const {
    assert(false); // none yet
    return false;
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
Mat33 FrameFeaturePlacementRep::getOrientationMeasureNumbers(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Mat33 value = Mat33(NTraits<Real>::getNaN());
    if (!isIndexed())
        value = FramePlacementRep::downcast(p).getOrientationMeasureNumbers(/*State*/);
    else
        assert(false);

    return value;
}

Vec3 FrameFeaturePlacementRep::getOriginMeasureNumbers(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Vec3 value = Vec3(NTraits<Real>::getNaN());
    if (!isIndexed())
        value = FramePlacementRep::downcast(p).getOriginMeasureNumbers(/*State*/);
    else
        assert(false);

    return value;
}

    // FRAME EXPR PLACEMENT REP //

bool FrameExprPlacementRep::isLimitedToSubtree
    (const Feature& root, const Feature*& offender) const
{
    if (!orientation.getRep().isLimitedToSubtree(root,offender))
        return false;
    return origin.getRep().isLimitedToSubtree(root,offender);
}

void FrameExprPlacementRep::repairFeatureReferences
    (const Feature& oldRoot, const Feature& newRoot)
{
    orientation.updRep().repairFeatureReferences(oldRoot,newRoot);
    origin.updRep().repairFeatureReferences(oldRoot,newRoot);
}


const Feature*
FrameExprPlacementRep::findAncestorFeature(const Feature& root) const {
    assert(!isConstant()); // don't call!

    const Feature* ancestor = 0;
    if (!orientation.isConstant())
        ancestor =  orientation.getRep().findAncestorFeature(root);

    const Feature* tmpAncestor = 0;
    if (!origin.isConstant())
        tmpAncestor =  origin.getRep().findAncestorFeature(root);

    if (ancestor && tmpAncestor)
        ancestor = FeatureRep::findYoungestCommonAncestor(*ancestor,*tmpAncestor);
    else ancestor = tmpAncestor;

    return ancestor;
}

    // FRAME EXPR PLACEMENT REP //

bool FrameOps::checkArgs(const std::vector<Placement>& args) const {
    assert(false); // none yet
    return false;
}

/////////////////////////////////////////////////////////////
// PlacementValue definitions & instantiations.            //
// Note that every supported subtype must be instantiated. //
/////////////////////////////////////////////////////////////

    // PLACEMENT VALUE //

PlacementValue::PlacementValue(const PlacementValue& src) : rep(0) {
    if (src.rep) src.rep->cloneUnownedWithNewHandle(*this);
}

PlacementValue::PlacementValue(PlacementValueRep* r) : rep(r) {
    assert(r && !r->hasHandle());
    r->setMyHandle(*this);
}

PlacementValue&
PlacementValue::operator=(const PlacementValue& src) {
    if (this != &src) {
        if (rep && (&rep->getMyHandle() == this)) delete rep; 
        rep=0;
        if (src.rep) src.rep->cloneUnownedWithNewHandle(*this);
    }
    return *this;
}

PlacementValue::~PlacementValue() {
    // This will blow up if rep doesn't have a handle -- we shouldn't
    // be pointing to it in that case!
    if (rep && (&rep->getMyHandle() == this)) delete rep; 
    rep=0;
}

bool PlacementValue::isValid() const {
    return rep && rep->isValid();
}

bool PlacementValue::hasOwner() const {
    return rep && rep->hasOwner();
}
int PlacementValue::getIndexInOwner() const {
    assert(rep && rep->hasOwner());
    assert(&rep->getMyHandle() == this);
    return rep->getIndexInOwner();
}

const Feature& PlacementValue::getOwner() const {
    assert(rep && rep->hasOwner());
    assert(&rep->getMyHandle() == this);
    return rep->getOwner();
}

String PlacementValue::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Placement Value";
    if (!rep) {
        s << "at 0x" << this << " HAS NULL REP";
        return s.str();
    }
    if (&rep->getMyHandle() != this) {
        s << "at 0x" << this << " HAS MISMATCHED REP";
        return s.str();
    }
    if (hasOwner())
        s << getOwner().getFullName() << ":"
          << std::left << std::setw(2) << getIndexInOwner();
    else s << "NO OWNER";
    s << " " << rep->toString(linePrefix);
    return s.str();
}

    // PLACEMENT VALUE <T> //

template <class T>
PlacementValue_<T>::PlacementValue_<T>() : PlacementValue()
  { rep = new PlacementValueRep_<T>(); }

template <class T>
PlacementValue_<T>::PlacementValue_<T>(const T& v) : PlacementValue()
  { rep = new PlacementValueRep_<T>(v); }

template <class T>
PlacementValue_<T>& PlacementValue_<T>::operator=(const T& v) {
    if (rep) set(v);
    else rep = new PlacementValueRep_<T>(v);
    return *this;
}

template <class T>
const T& PlacementValue_<T>::get() const {
    return PlacementValueRep_<T>::downcast(getRep()).getValue();
}

template <class T>
PlacementValue_<T>::operator const T&() const { return get(); }

template <class T>
void PlacementValue_<T>::set(const T& v) { 
    PlacementValueRep_<T>::downcast(updRep()).setValue(v);
}

template <class T>
bool PlacementValue_<T>::isInstanceOf(const PlacementValue& pv) {
    return PlacementValueRep_<T>::isA(pv.getRep());
}

template <class T>
const PlacementValue_<T>& PlacementValue_<T>::downcast(const PlacementValue& pv) {
    assert(isInstanceOf(pv)); 
    return reinterpret_cast<const PlacementValue_<T>&>(pv);
}

template <class T>
PlacementValue_<T>& PlacementValue_<T>::downcast(PlacementValue& pv) { 
    assert(isInstanceOf(pv)); 
    return reinterpret_cast<PlacementValue_<T>&>(pv);
}

template class PlacementValue_<Real>;
template class PlacementValue_<Vec3>;
template class PlacementValue_<Mat33>;

} // namespace simtk
