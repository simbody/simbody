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

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Geometry.h"
#include "simbody/internal/Mechanics.h"
#include "simbody/internal/Placement.h"
#include "simbody/internal/Feature.h"

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

void PlacementExpr::exprRealize(Stage g) const {
    for (size_t i=0; i < args.size(); ++i)
        args[i].realize(g);
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
    const Placement& sample = f.getRep().getSamplePlacement();
    const int nElements = sample.getRep().getNIndicesAllowed();

    // TODO: should this throw a nice message, or should we check higher up?
    assert(nElements > 0 // i.e., not void
        && (i == -1 || (0 <= i && i < nElements)));
}


void FeatureReference::refRealize(/*State,*/ Stage g) const {
    if (refGetReferencedFeature().hasPlacement())
        refGetReferencedFeature().getRep().getPlacementSlot().realize(g);
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
    const PlacementSlot& ps = feature->getRep().getPlacementSlot();
    assert(ps.hasOwner());  // a feature can't be placed on an unowned placement
    const Subsystem* newYoungest = SubsystemRep::findYoungestCommonAncestor(youngestAllowed,
                                                                            ps.getOwner());
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

bool
FeatureReference::refDependsOn(const Feature& f) const {
    assert(feature); 
    return feature->isSameSubsystem(f)
        || feature->hasPlacement() && feature->getPlacement().dependsOn(f);
}

std::string FeatureReference::refToString(const std::string&) const {
    std::stringstream s;

    s << "Ref<";
    s << (feature ? feature->getRep().getSamplePlacement().getPlacementTypeName()
                  : "NONE");
    s << ">[";
    s << (feature ? feature->getFullName()
                    : std::string("NULL FEATURE"));
    s << "]"; 
    if (index != -1) s << "[" << index << "]";
    return s.str();
}

    // PLACEMENT REP //


// These are the default implementations for the generic operators. Any 
// concrete class which thinks it knows how to perform one of these
// operations should override.

/*virtual*/ Placement PlacementRep::genericNegate() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "-", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericAbs() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "abs", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericSqrt() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "sqrt", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericSquare() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "square", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericCube() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "cube", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericExp() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "exp", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericLog() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "log", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericSin() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "sin", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericCos() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "cos", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericAsin() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "asin", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericAcos() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "acos", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericLength() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "length", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericNormalize() const {
    SIMTK_THROW2(Exception::UnaryOperationNotAllowedForPlacementType,
                    "normalize", getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericAdd(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "+", r.getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericSub(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "-", r.getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericMul(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "*", r.getPlacementTypeName());
}
/*virtual*/ Placement PlacementRep::genericDvd(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "/", r.getPlacementTypeName());
}

/*virtual*/ Placement PlacementRep::genericDistance(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "distance", r.getPlacementTypeName());
}

/*virtual*/ Placement PlacementRep::genericAngle(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "angle", r.getPlacementTypeName());
}

/*virtual*/ Placement PlacementRep::genericDotProduct(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "dot", r.getPlacementTypeName());
}

/*virtual*/ Placement PlacementRep::genericCrossProduct(const Placement& r) const {
    SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                getPlacementTypeName(),
                "cross", r.getPlacementTypeName());
}

    // REAL PLACEMENT REP //

/*static*/ PlacementRep* 
RealPlacementRep::createRealPlacementFrom(const Placement& p, bool dontThrow) {
    if (RealPlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Real", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// result = -real (result is always real)
Placement RealPlacementRep::genericNegate() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(-calcRealValue())
                     : (RealPlacementRep*)RealExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = abs(real) (result is always real)
Placement RealPlacementRep::genericAbs() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::abs(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::absOp(getMyHandle());
    return Placement(result);
}

// result = sqrt(real) (result is always real)
Placement RealPlacementRep::genericSqrt() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sqrt(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sqrtOp(getMyHandle());
    return Placement(result);
}

// result = square(real) (result is always real)
Placement RealPlacementRep::genericSquare() const {
    if (isConstant()) {
        const Real r = calcRealValue();
        return Placement(new RealConstantPlacementRep(r*r));
    }
    return Placement(RealExprPlacementRep::squareOp(getMyHandle()));
}

// result = cube(real) (result is always real)
Placement RealPlacementRep::genericCube() const {
    if (isConstant()) {
        const Real r = calcRealValue();
        return Placement(new RealConstantPlacementRep(r*r*r));
    }
    return Placement(RealExprPlacementRep::cubeOp(getMyHandle()));
}

// result = exp(real) (result is always real)
Placement RealPlacementRep::genericExp() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::exp(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::expOp(getMyHandle());
    return Placement(result);
}

// result = log(real) (natural log, result is always real)
Placement RealPlacementRep::genericLog() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::log(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::logOp(getMyHandle());
    return Placement(result);
}

// result = sin(real) (result is always real)
Placement RealPlacementRep::genericSin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::sin(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::sinOp(getMyHandle());
    return Placement(result);
}

// result = cos(real) (result is always real)
Placement RealPlacementRep::genericCos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::cos(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::cosOp(getMyHandle());
    return Placement(result);
}

// result = asin(real) (result is always real)
Placement RealPlacementRep::genericAsin() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::asin(calcRealValue()))
                     : (RealPlacementRep*)RealExprPlacementRep::asinOp(getMyHandle());
    return Placement(result);
}

// result = acos(real) (result is always real)
Placement RealPlacementRep::genericAcos() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(std::acos(calcRealValue()))
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcRealValue()+rp.getRep().calcRealValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcRealValue()-rp.getRep().calcRealValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcRealValue()*rp.getRep().calcRealValue())
                : (RealPlacementRep*)RealExprPlacementRep::mulOp(getMyHandle(),rp);
        return Placement(result);
    }

    if (Vec3Placement::isInstanceOf(r))
        return Vec3Placement::downcast(r) * getMyHandle();   // punt to vec3*real

    if (DirectionPlacement::isInstanceOf(r) || StationPlacement::isInstanceOf(r))
        return Vec3Placement::convert(r) * getMyHandle(); // punt to vec3(r) * real

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
                ? (RealPlacementRep*)new RealConstantPlacementRep(calcRealValue()/rp.getRep().calcRealValue())
                : (RealPlacementRep*)RealExprPlacementRep::dvdOp(getMyHandle(),rp);
        return Placement(result);
    }

    return PlacementRep::genericDvd(r);    // die
}

    // REAL FEATURE PLACEMENT REP //
const Real& RealFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed()) 
        return PlacementValue_<Real>::downcast(ps.getValue()).get();

    // indexed
    const PlacementRep& p = ps.getPlacement().getRep();
    if (Vec3PlacementRep::isA(p) || StationPlacementRep::isA(p) || DirectionPlacementRep::isA(p))
        return PlacementValue_<Vec3>::downcast(ps.getValue()).get()[getPlacementIndex()];

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const Real*>(0);
}

    // REAL EXPR PLACEMENT REP //

bool 
RealOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    // unary
    case Negate: case Abs: case Sqrt: case Square: case Cube: case Exp:
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

Real RealOps::apply(const std::vector<Placement>& args) const {
    Real val = NTraits<Real>::getNaN();
    Real arg1, arg2;
    if (args.size() > 0 && RealPlacement::isInstanceOf(args[0]))
        arg1 = RealPlacement::downcast(args[0]).getRep().calcRealValue();
    if (args.size() > 1 && RealPlacement::isInstanceOf(args[1]))
        arg2 = RealPlacement::downcast(args[1]).getRep().calcRealValue();

    switch (op) {
    case Negate: val = -arg1;           break;
    case Abs:    val = std::abs(arg1);  break;
    case Sqrt:   val = std::sqrt(arg1); break;
    case Square: val = arg1*arg1;       break;
    case Cube:   val = arg1*arg1*arg1;  break;
    case Exp:    val = std::exp(arg1);  break;
    case Log:    val = std::log(arg1);  break;
    case Sin:    val = std::sin(arg1);  break;
    case Cos:    val = std::cos(arg1);  break;
    case Asin:   val = std::asin(arg1); break;
    case Acos:   val = std::acos(arg1); break;
    case VectorLength:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value().norm(); 
        break;

    case Add:      val = arg1+arg2; break;
    case Subtract: val = arg1-arg2; break;
    case Multiply: val = arg1*arg2; break;
    case Divide:   val = arg1/arg2; break; 

    case DotProduct3: {
        const Vec3 l = Vec3Placement::downcast(args[0]).getRep().calcVec3Value();
        const Vec3 r = Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
        val = (~l)*r;
        break;
    }
    case DotProduct2: {
        assert(false);
        break;
    }
    case PointDistance: {
        const Vec3 head = StationPlacement::downcast(args[0]).getRep().calcVec3Value();
        const Vec3 tail = StationPlacement::downcast(args[1]).getRep().calcVec3Value();
        val = (head-tail).norm();
        break;
    }
    case AngleBetweenDirections: {
        const UnitVec3 v1 = DirectionPlacement::downcast(args[0]).getRep().calcUnitVec3Value();
        const UnitVec3 v2 = DirectionPlacement::downcast(args[1]).getRep().calcUnitVec3Value();
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
RealExprPlacementRep::squareOp(const RealPlacement& rp)
  { return unaryOp(RealOps::Square, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::cubeOp  (const RealPlacement& rp)
  { return unaryOp(RealOps::Cube, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::expOp   (const RealPlacement& rp)
  { return unaryOp(RealOps::Exp, rp); }
/*static*/ RealExprPlacementRep*
RealExprPlacementRep::logOp   (const RealPlacement& rp)
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

/*static*/ PlacementRep* 
Vec3PlacementRep::createVec3PlacementFrom(const Placement& p, bool dontThrow) {
    if (Vec3Placement::isInstanceOf(p))
        return p.getRep().clone();

    if (StationPlacement::isInstanceOf(p)) {
        const StationPlacement& sp = StationPlacement::downcast(p);
        return sp.isConstant()
            ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(sp.getRep().calcVec3Value())
            : (Vec3PlacementRep*)Vec3ExprPlacementRep::recastStationOp(sp);
    } else if (DirectionPlacement::isInstanceOf(p)) {
        const DirectionPlacement& dp = DirectionPlacement::downcast(p);
        return dp.isConstant()
            ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(dp.getRep().calcUnitVec3Value().asVec3())
            : (Vec3PlacementRep*)Vec3ExprPlacementRep::recastDirectionOp(dp);
    }
    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Vec3", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// result = -vec3 (result is always vec3)
Placement Vec3PlacementRep::genericNegate() const {
    Vec3PlacementRep* result =
        isConstant() ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(-calcVec3Value())
                     : (Vec3PlacementRep*)Vec3ExprPlacementRep::negateOp(getMyHandle());
    return Placement(result);
}

// result = length(vec3) (result is always real)
Placement Vec3PlacementRep::genericLength() const {
    RealPlacementRep* result = 
        isConstant() ? (RealPlacementRep*)new RealConstantPlacementRep(calcVec3Value().norm())
                     : (RealPlacementRep*)RealExprPlacementRep::lengthOp(getMyHandle());
    return Placement(result);
}

// result = normalize(vec3) (result is always Direction)
Placement Vec3PlacementRep::genericNormalize() const {
    DirectionPlacementRep* result = 
        isConstant() ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(UnitVec3(calcVec3Value()))
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcVec3Value()+vp.getRep().calcVec3Value())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcVec3Value()-vp.getRep().calcVec3Value())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcVec3Value()*rp.getRep().calcRealValue())
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
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep(calcVec3Value()/rp.getRep().calcRealValue())
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
                ? (RealPlacementRep*)new RealConstantPlacementRep(~calcVec3Value()*vp.getRep().calcVec3Value())
                : (RealPlacementRep*)RealExprPlacementRep::dot3Op(getMyHandle(),vp);
        return Placement(result);
    }

    if (Vec3Placement::canConvert(r))
        return genericDotProduct(Vec3Placement::convert(r)); // recursive call with improved argument

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
                   new Vec3ConstantPlacementRep(cross(calcVec3Value(),vp.getRep().calcVec3Value()))
                : (Vec3PlacementRep*)
                   Vec3ExprPlacementRep::crossOp(getMyHandle(),vp);
        return Placement(result);
    }
    
    if (Vec3Placement::canConvert(r))
        return genericCrossProduct(Vec3Placement::convert(r));

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
    if (Vec3Placement::canConvert(r)) 
        return genericAngle(Vec3Placement::convert(r));

    return PlacementRep::genericAngle(r);    // die
}

    // VEC3 FEATURE PLACEMENT REP //
const Vec3& Vec3FeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed())
        return PlacementValue_<Vec3>::downcast(ps.getValue()).get();

    // indexed
    const PlacementRep& p = ps.getPlacement().getRep();
    if (OrientationPlacementRep::isA(p))
        return PlacementValue_<Mat33>::downcast(ps.getValue()).get()
            (getPlacementIndex()); // round brackets () to get column not row

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
        val = StationPlacement::downcast(args[0]).getRep().calcVec3Value();
        break;

    case RecastDirection:
        val = DirectionPlacement::downcast(args[0]).getRep().calcUnitVec3Value().asVec3();
        break;

    case Negate:
        val = -Vec3Placement::downcast(args[0]).getRep().calcVec3Value();
        break;

    case Add:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value()
              + Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
        break;
    
    case Subtract:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value()
              - Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
        break;

    case StationDifference:
        val = StationPlacement::downcast(args[0]).getRep().calcVec3Value()
              - StationPlacement::downcast(args[1]).getRep().calcVec3Value();
        break;

    // real is always on the right for scalar mul & dvd
    case ScalarMultiply:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value()
              * RealPlacement::downcast(args[1]).getRep().calcRealValue();
        break;

    case ScalarDivide:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value()
              / RealPlacement::downcast(args[1]).getRep().calcRealValue();
        break;

    case CrossProduct:
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value()
              % Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
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
  { return binaryOp(Vec3Ops::ScalarMultiply, Vec3Placement(l), r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::smulOp(const DirectionPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarMultiply, Vec3Placement(l), r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const Vec3Placement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, l, r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const StationPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, Vec3Placement(l), r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::sdvdOp(const DirectionPlacement& l, const RealPlacement& r)
  { return binaryOp(Vec3Ops::ScalarDivide, Vec3Placement(l), r); }

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::crossOp(const Vec3Placement& l, const Vec3Placement& r)
  { return binaryOp(Vec3Ops::CrossProduct, l, r); }

    // STATION PLACEMENT REP //

/*static*/ PlacementRep*
StationPlacementRep::createStationPlacementFrom(const Placement& p, bool dontThrow) {
    if (StationPlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (Vec3Placement::isInstanceOf(p)) {
        const Vec3Placement& vp = Vec3Placement::downcast(p);
        return vp.isConstant()
            ? (StationPlacementRep*)new StationConstantPlacementRep(vp.getRep().calcVec3Value())
            : (StationPlacementRep*)StationExprPlacementRep::recastVec3Op(vp);
    }
    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Station", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// result = -station (result is always vec3)
Placement StationPlacementRep::genericNegate() const {
    return -Vec3Placement::convert(getMyHandle());
}

// result = length(station) (result is always real)
Placement StationPlacementRep::genericLength() const {
    return length(Vec3Placement::convert(getMyHandle()));
}

// result = normalize(station) (result is always Direction)
Placement StationPlacementRep::genericNormalize() const {
    return normalize(Vec3Placement::convert(getMyHandle()));
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
                                              (calcVec3Value()+vp.getRep().calcVec3Value())
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
                                              (calcVec3Value()-vp.getRep().calcVec3Value())
                : (StationPlacementRep*)StationExprPlacementRep::subOp(getMyHandle(),vp);
        return Placement(result);
    }

    if (StationPlacement::isInstanceOf(r)) {
        const StationPlacement& sp = StationPlacement::downcast(r);
        Vec3PlacementRep* result = 
            isConstant() && sp.getRep().isConstant()
                ? (Vec3PlacementRep*)new Vec3ConstantPlacementRep
                                              (calcVec3Value()-sp.getRep().calcVec3Value())
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
        return Vec3Placement::convert(getMyHandle()) * r;
    }

    return PlacementRep::genericMul(r);    // die
}

// result = station / placement
//   We support only
//     vec3 = station / real
Placement StationPlacementRep::genericDvd(const Placement& r) const {
    if (RealPlacement::isInstanceOf(r)) {
        return Vec3Placement::convert(getMyHandle()) / r;
    }

    return PlacementRep::genericDvd(r);    // die
}

// result = dot(station, placement)
//   We support 
//     real = dot(station, vec3)
//     real = dot(station, direction)
Placement StationPlacementRep::genericDotProduct(const Placement& r) const {
    if (Vec3Placement::isInstanceOf(r) || DirectionPlacement::isInstanceOf(r)) {
        return dot(Vec3Placement::convert(getMyHandle()), r);
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
        return cross(Vec3Placement::convert(getMyHandle()), r);
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
const Vec3& StationFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();
    if (!isIndexed())
        return PlacementValue_<Vec3>::downcast(ps.getValue()).get();

    // indexed
    const PlacementRep& p = ps.getPlacement().getRep();
    if (FramePlacementRep::isA(p) && getPlacementIndex()==1)
        return PlacementValue_<TransformMat>::downcast(ps.getValue()).get().getTranslation();

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
        val = Vec3Placement::downcast(args[0]).getRep().calcVec3Value();
        break;

    case Add:
        val = StationPlacement::downcast(args[0]).getRep().calcVec3Value()
              + Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
        break;

    case Subtract:
        val = StationPlacement::downcast(args[0]).getRep().calcVec3Value()
              - Vec3Placement::downcast(args[1]).getRep().calcVec3Value();
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

/*static*/ PlacementRep*
DirectionPlacementRep::createDirectionPlacementFrom(const Placement& p, bool dontThrow) {
    if (DirectionPlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (Vec3Placement::isInstanceOf(p)) {
        const Vec3Placement& vp = Vec3Placement::downcast(p);
        return vp.isConstant()
            ? (DirectionPlacementRep*)new DirectionConstantPlacementRep(UnitVec3(vp.getRep().calcVec3Value()))
            : (DirectionPlacementRep*)DirectionExprPlacementRep::normalizeOp(vp);
    }
    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Direction", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// result = -direction (result is always direction)
Placement DirectionPlacementRep::genericNegate() const {
    DirectionPlacementRep* result =
        isConstant() ? (DirectionPlacementRep*)
                        new DirectionConstantPlacementRep(-calcUnitVec3Value())
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
                   new Vec3ConstantPlacementRep(calcUnitVec3Value()*rp.getRep().calcRealValue())
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
                   new Vec3ConstantPlacementRep(calcUnitVec3Value()/rp.getRep().calcRealValue())
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
        return dot(Vec3Placement::convert(getMyHandle()), r);
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
        return cross(Vec3Placement::convert(getMyHandle()), r);
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
            Real dotprod = ~calcUnitVec3Value().asVec3() * dp.getRep().calcUnitVec3Value().asVec3();
            if (dotprod < -1.) dotprod = -1.;
            else if (dotprod > 1.) dotprod = 1.;
            return RealPlacement(new RealConstantPlacementRep(std::acos(dotprod)));
        }
        return RealPlacement(RealExprPlacementRep::angleOp(getMyHandle(),dp));
    }

    return PlacementRep::genericAngle(r);    // die
}

    // DIRECTION FEATURE PLACEMENT REP //
const UnitVec3& DirectionFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed())
        return PlacementValue_<UnitVec3>::downcast(ps.getValue()).get();

    // indexed
    const PlacementRep& p = ps.getPlacement().getRep();
    if (OrientationPlacementRep::isA(p))
        return PlacementValue_<RotationMat>::downcast(ps.getValue()).get()
            (getPlacementIndex()); // round brackets () to get column not row

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const UnitVec3*>(0);
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

UnitVec3 DirectionOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    UnitVec3 val;
    switch (op) {

    case Negate:
        val = -DirectionPlacement::downcast(args[0]).getRep().calcUnitVec3Value();
        break;

    case NormalizeVec3:
        val = UnitVec3(Vec3Placement::downcast(args[0]).getRep().calcVec3Value());
        break;

    case NormalizeStation:
        val = UnitVec3(StationPlacement::downcast(args[0]).getRep().calcVec3Value());
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

/*static*/ PlacementRep* 
OrientationPlacementRep::createOrientationPlacementFrom(const Placement& p, bool dontThrow) {
    if (OrientationPlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Orientation", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// rotation = this->invert()
OrientationPlacement 
OrientationPlacementRep::invert() const
{
    OrientationPlacementRep* result = 
        isConstant() ? (OrientationPlacementRep*)new OrientationConstantPlacementRep(~calcMatRotationValue())
                     : (OrientationPlacementRep*)OrientationExprPlacementRep::invertOp(getMyHandle());
    return OrientationPlacement(result);
}

// rotation = createFromZAxis(z)
/*static*/ OrientationPlacement 
OrientationPlacementRep::createFromZAxis(const DirectionPlacement& z)
{
    const bool calcNow = z.isConstant();

    OrientationPlacementRep* result = 
        calcNow ? (OrientationPlacementRep*)new OrientationConstantPlacementRep(RotationMat(z.getRep().calcUnitVec3Value()))
                : (OrientationPlacementRep*)OrientationExprPlacementRep::createFromZOp(z);

    return OrientationPlacement(result);
}

    // ORIENTATION FEATURE PLACEMENT REP //
const RotationMat& OrientationFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed())
        return PlacementValue_<RotationMat>::downcast(ps.getValue()).get();

    // indexed
    const PlacementRep& p = ps.getPlacement().getRep();

    if (FramePlacementRep::isA(p) && getPlacementIndex()==0)
        return PlacementValue_<TransformMat>::downcast(ps.getValue()).get().getRotation();

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const RotationMat*>(0);
}

    // ORIENTATION EXPR PLACEMENT REP //

bool OrientationOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {

    // given the direction of the z axis, create a complete set of 3 axes
    case CreateFromZAxis:
        return args.size()==1 && DirectionPlacement::isInstanceOf(args[0]);

    // given a rotation matrix, return its inverse
    case Invert:
        return args.size()==1 && OrientationPlacement::isInstanceOf(args[0]);

    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

RotationMat OrientationOps::apply(const std::vector<Placement>& args) const {
    RotationMat val;
    switch(op) {

    case CreateFromZAxis:
        val = RotationMat(DirectionPlacement::downcast(args[0]).getRep().calcUnitVec3Value());
        break;

    case Invert:
        val = ~OrientationPlacement::downcast(args[0]).getRep().calcMatRotationValue();
        break;

    default:
        assert(false);
        //NOTREACHED
    }
    return val;
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

/*static*/ OrientationExprPlacementRep*
OrientationExprPlacementRep::createFromZOp(const DirectionPlacement& z)
  { return unaryOp(OrientationOps::CreateFromZAxis, z); }

/*static*/ OrientationExprPlacementRep*
OrientationExprPlacementRep::invertOp(const OrientationPlacement& r)
  { return unaryOp(OrientationOps::Invert, r); }

    // INERTIA PLACEMENT REP //

/*static*/ PlacementRep* 
InertiaPlacementRep::createInertiaPlacementFrom(const Placement& p, bool dontThrow) {
    if (InertiaPlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "InertiaMat", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

// result = inertia + placement
//   We support only
//     inertia = inertia + inertia
Placement InertiaPlacementRep::genericAdd(const Placement& r) const {
    if (InertiaPlacement::isInstanceOf(r)) {
        const InertiaPlacement& ip = InertiaPlacement::downcast(r);
        InertiaPlacementRep* result = 
            isConstant() && ip.getRep().isConstant()
                ? (InertiaPlacementRep*)new InertiaConstantPlacementRep(calcInertiaValue()+ip.getRep().calcInertiaValue())
                : (InertiaPlacementRep*)InertiaExprPlacementRep::addOp(getMyHandle(),ip);
        return Placement(result);
    }

    return PlacementRep::genericAdd(r);    // die
}

// result = inertia - placement
//   We support only
//     inertia = inertia - inertia
Placement InertiaPlacementRep::genericSub(const Placement& r) const {
    if (InertiaPlacement::isInstanceOf(r)) {
        const InertiaPlacement& ip = InertiaPlacement::downcast(r);
        InertiaPlacementRep* result = 
            isConstant() && ip.getRep().isConstant()
                ? (InertiaPlacementRep*)new InertiaConstantPlacementRep(calcInertiaValue()-ip.getRep().calcInertiaValue())
                : (InertiaPlacementRep*)InertiaExprPlacementRep::subOp(getMyHandle(),ip);
        return Placement(result);
    }

    return PlacementRep::genericSub(r);    // die
}

// inertia = changeAxes(inertia, axes)
InertiaPlacement 
InertiaPlacementRep::changeAxes(const OrientationPlacement& r) const
{
    const bool calcNow = isConstant() && r.isConstant();

    InertiaPlacementRep* result = 
        calcNow ? (InertiaPlacementRep*)new InertiaConstantPlacementRep(
                    calcInertiaValue().changeAxes(r.getRep().calcMatRotationValue()))
                : (InertiaPlacementRep*)InertiaExprPlacementRep::changeAxesOp
                                            (getMyHandle(),r);
    return InertiaPlacement(result);
}

// inertia = shiftFromCOM(inertia, to, totalMass)
InertiaPlacement 
InertiaPlacementRep::shiftFromCOM(const StationPlacement& to,
                                  const RealPlacement&    totalMass) const
{
    const bool calcNow = isConstant() && to.isConstant() && totalMass.isConstant();

    InertiaPlacementRep* result = 
        calcNow ? (InertiaPlacementRep*)new InertiaConstantPlacementRep(
                    calcInertiaValue().shiftFromCOM(to.getRep().calcVec3Value(),
                                                    totalMass.getRep().calcRealValue()))
                : (InertiaPlacementRep*)InertiaExprPlacementRep::shiftFromCOMOp
                                            (getMyHandle(),to,totalMass);
    return InertiaPlacement(result);
}

// inertia = shiftToCOM(inertia, com, totalMass)
InertiaPlacement 
InertiaPlacementRep::shiftToCOM(const StationPlacement& com,
                                const RealPlacement&    totalMass) const
{
    const bool calcNow = isConstant() && com.isConstant() && totalMass.isConstant();

    InertiaPlacementRep* result = 
        calcNow ? (InertiaPlacementRep*)new InertiaConstantPlacementRep(
                    calcInertiaValue().shiftToCOM(com.getRep().calcVec3Value(),
                                                  totalMass.getRep().calcRealValue()))
                : (InertiaPlacementRep*)InertiaExprPlacementRep::shiftToCOMOp
                                            (getMyHandle(),com,totalMass);
    return InertiaPlacement(result);
}

    // INERTIA FEATURE PLACEMENT REP //
const InertiaMat& InertiaFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed())
        return PlacementValue_<InertiaMat>::downcast(ps.getValue()).get();

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const InertiaMat*>(0);
}

    // INERTIA EXPR PLACEMENT REP //
bool InertiaOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {

    // i=i+i
    case Add:
        return args.size()==2 && InertiaPlacement::isInstanceOf(args[0])
                              && InertiaPlacement::isInstanceOf(args[1]);
    
    // i=i-i
    case Subtract:
        return args.size()==2 && InertiaPlacement::isInstanceOf(args[0])
                              && InertiaPlacement::isInstanceOf(args[1]);

    // i=transform(i, rotation)
    case ChangeAxes:
        return args.size()==2 && InertiaPlacement::isInstanceOf(args[0])
                              && OrientationPlacement::isInstanceOf(args[1]);

    // i=shift(i,pt,totalMass)
    case ShiftToCOM:
    case ShiftFromCOM:
        return args.size()==3 && InertiaPlacement::isInstanceOf(args[0])
                              && StationPlacement::isInstanceOf(args[1])
                              && RealPlacement::isInstanceOf(args[2]);

    // i=pointMass(loc, mass)
    case PointMass:
        return args.size()==2 && StationPlacement::isInstanceOf(args[0])
                              && RealPlacement::isInstanceOf(args[1]);

    case PrincipalMoments: {
        if (args.size() != 3) return false;
        for (size_t i=0; i<args.size(); ++i)
            if (!RealPlacement::isInstanceOf(args[i])) return false;
        return true;
    }

    case FullInertia: {
        if (args.size() != 6) return false;
        for (size_t i=0; i<args.size(); ++i)
            if (!RealPlacement::isInstanceOf(args[i])) return false;
        return true;
    }

    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

InertiaMat InertiaOps::apply(const std::vector<Placement>& args) const {
    InertiaMat val;
    switch(op) {
    case Add:
        val = InertiaPlacement::downcast(args[0]).getRep().calcInertiaValue()
              + InertiaPlacement::downcast(args[1]).getRep().calcInertiaValue();
        break;
    
    case Subtract:
        val = InertiaPlacement::downcast(args[0]).getRep().calcInertiaValue()
              - InertiaPlacement::downcast(args[1]).getRep().calcInertiaValue();
        break;

    case ChangeAxes:
        val = InertiaPlacement::downcast(args[0]).getRep().calcInertiaValue()
              .changeAxes(OrientationPlacement::downcast(args[1]).getRep().calcMatRotationValue());
        break;

    // i=shift(i,to,totalMass)
    case ShiftFromCOM:
        val = InertiaPlacement::downcast(args[0]).getRep().calcInertiaValue()
              .shiftFromCOM(StationPlacement::downcast(args[1]).getRep().calcVec3Value(),
                            RealPlacement::downcast   (args[2]).getRep().calcRealValue());
        break;

    // i=shift(i,to,totalMass)
    case ShiftToCOM:
        val = InertiaPlacement::downcast(args[0]).getRep().calcInertiaValue()
              .shiftToCOM(StationPlacement::downcast(args[1]).getRep().calcVec3Value(),
                          RealPlacement::downcast   (args[2]).getRep().calcRealValue());
        break;

    // i=pointMass(loc,mass)
    case PointMass:
        val = InertiaMat(StationPlacement::downcast(args[0]).getRep().calcVec3Value(),
                         RealPlacement::downcast(args[1]).getRep().calcRealValue());
        break;

    // i=principalMoments(Ixx,Iyy,Izz)
    case PrincipalMoments:
        val = InertiaMat(RealPlacement::downcast(args[0]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[1]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[2]).getRep().calcRealValue());
        break;

    // i=fullInertia(Ixx,Iyy,Izz,Ixy,Iyz,Ixz)
    case FullInertia:
        val = InertiaMat(RealPlacement::downcast(args[0]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[1]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[2]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[3]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[4]).getRep().calcRealValue(),
                         RealPlacement::downcast(args[5]).getRep().calcRealValue());
        break;

    default: 
        assert(false);
    }
    return val;
}

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::unaryOp(InertiaOps::OpKind op, const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new InertiaExprPlacementRep(InertiaOps(op), args);
}

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::binaryOp(InertiaOps::OpKind op, 
                               const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new InertiaExprPlacementRep(InertiaOps(op), args);
}

// Supported InertiaExpr-building operators

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::addOp (const InertiaPlacement& l, const InertiaPlacement& r)
  { return binaryOp(InertiaOps::Add, l, r); }
/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::subOp (const InertiaPlacement& l, const InertiaPlacement& r)
  { return binaryOp(InertiaOps::Subtract, l, r); }

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::changeAxesOp (const InertiaPlacement& i, const OrientationPlacement& r)
  { return binaryOp(InertiaOps::ChangeAxes, i, r); }

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::shiftFromCOMOp(const InertiaPlacement& i,
                                        const StationPlacement& to,
                                        const RealPlacement&    mtot)
{ 
    std::vector<const Placement*> args(3);
    args[0] = &i; args[1] = &to; args[2] = &mtot;
    return new InertiaExprPlacementRep(InertiaOps(InertiaOps::ShiftFromCOM), args);
}

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::shiftToCOMOp(const InertiaPlacement& i,
                                      const StationPlacement& com,
                                      const RealPlacement&    mtot)
{ 
    std::vector<const Placement*> args(3);
    args[0] = &i; args[1] = &com; args[2] = &mtot;
    return new InertiaExprPlacementRep(InertiaOps(InertiaOps::ShiftToCOM), args);
}

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::ptMassOp(const StationPlacement& loc,
                                  const RealPlacement&    mass)
  { return binaryOp(InertiaOps::PointMass, loc, mass); }

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::principalMomentsOp(const RealPlacement& Ixx,
                                            const RealPlacement& Iyy,
                                            const RealPlacement& Izz)
{
    std::vector<const Placement*> args(3);
    args[0] = &Ixx; args[1] = &Iyy; args[2] = &Izz;
    return new InertiaExprPlacementRep(InertiaOps(InertiaOps::PrincipalMoments), args);
}

/*static*/ InertiaExprPlacementRep*
InertiaExprPlacementRep::fullInertiaOp(const RealPlacement& Ixx, const RealPlacement& Iyy, const RealPlacement& Izz,
                                       const RealPlacement& Ixy, const RealPlacement& Iyz, const RealPlacement& Ixz)
{
    std::vector<const Placement*> args(6);
    args[0] = &Ixx; args[1] = &Iyy; args[2] = &Izz; args[3] = &Ixy; args[4] = &Iyz; args[5] = &Ixz;
    return new InertiaExprPlacementRep(InertiaOps(InertiaOps::FullInertia), args);
}

    // FRAME PLACEMENT REP //

// In addition to a straightforward Frame placement, a FrameFeature can be
// placed on a Station or Orientation Feature, if we can pick up the missing
// piece from that Feature's parent Frame. So we'll say yes if the proposed
// Placement is an unindexed Feature reference to either a Station or
// Orientation Feature, provided that Feature's parent is a FrameFeature.

/*static*/ PlacementRep* 
FramePlacementRep::createFramePlacementFrom(const Placement& p, bool dontThrow) {
    if (FramePlacement::isInstanceOf(p))
        return p.getRep().clone();

    if (p.isFeatureReference()) {
        const Feature& refFeature = p.getReferencedFeature();

        if (refFeature.hasParentSubsystem() 
            && FrameFeature::isInstanceOf(refFeature.getParentSubsystem()))
        {
            const FrameFeature& frame = FrameFeature::downcast(refFeature.getParentSubsystem());

            // Check that both the Placement and reference feature have the same type to
            // make sure there is no indexing.
            if (StationPlacement::isInstanceOf(p) && Station::isInstanceOf(refFeature))
                return new FrameExprPlacementRep(OrientationPlacement(frame.getOrientation()), 
                                                 StationPlacement::downcast(p));
            else if (OrientationPlacement::isInstanceOf(p) && Orientation::isInstanceOf(refFeature)) 
                return new FrameExprPlacementRep(OrientationPlacement::downcast(p), 
                                                 StationPlacement(frame.getOrigin()));
        }
    }
    if (!dontThrow) {
        SIMTK_THROW3(Exception::PlacementCantBeConvertedToRightType,
            "Frame", p.getPlacementTypeName(), p.toString());
        //NOTREACHED
    }
    return 0;
}

    // FRAME FEATURE PLACEMENT REP //
const TransformMat& FrameFeaturePlacementRep::getReferencedValue() const {
    const PlacementSlot& ps = getReferencedFeature().getRep().getPlacementSlot();

    if (!isIndexed())
        return PlacementValue_<TransformMat>::downcast(ps.getValue()).get();

    assert(false);
    //NOTREACHED

    return *reinterpret_cast<const TransformMat*>(0);
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

TransformMat FrameOps::apply(/*State,*/ const std::vector<Placement>& args) const {
    assert(false);
    return TransformMat();
}

} // namespace simtk
