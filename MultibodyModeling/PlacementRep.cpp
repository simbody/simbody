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

    // PLACEMENT REP //


    // FEATURE PLACEMENT REP //

const Feature* 
FeaturePlacementRep::findAncestorFeature(const Feature& root) const {
    assert(feature);
    return FeatureRep::isFeatureInFeatureTree(root, *feature) 
            ? feature : 0;
}

PlacementType FeaturePlacementRep::getPlacementType() const {
    assert(feature);
    const PlacementType whole = (*feature).getRep().getRequiredPlacementType();
    return index == -1 ? whole : getIndexedPlacementType(whole, index);
}

// Check that this feature is on the feature subtree rooted by "ancestor". If
// not return a pointer to this feature for use in a friendly error message.
// If this is the right tree, we return true with offender==NULL.
bool FeaturePlacementRep::isLimitedToSubtree
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

void FeaturePlacementRep::repairFeatureReferences
    (const Feature& oldRoot, const Feature& newRoot)
{
    assert(feature);
    const Feature* corrFeature = FeatureRep::findCorrespondingFeature(oldRoot,*feature,newRoot);
    assert(corrFeature);
    feature = corrFeature;
}

std::string FeaturePlacementRep::toString(const std::string&) const {
    std::stringstream s;
    s << "Feature[";
    s << (feature ? feature->getFullName()
                    : std::string("NULL FEATURE"));
    s << "]"; 
    if (index != -1) s << "[" << index << "]";
    return s.str();
}

    // FRAME PLACEMENT REP //
bool FramePlacementRep::isLimitedToSubtree
    (const Feature& root, const Feature*& offender) const
{
    assert(orientation && station);
    if (!FeatureRep::isFeatureInFeatureTree(root, *orientation)) {
        offender = orientation;
        return false;
    }
    if (!FeatureRep::isFeatureInFeatureTree(root, *station)) {
        offender = station;
        return false;
    }
    offender = 0;
    return true;
}

void FramePlacementRep::repairFeatureReferences
    (const Feature& oldRoot, const Feature& newRoot)
{
    assert(orientation && station);
    const Feature* corrOri = FeatureRep::findCorrespondingFeature(oldRoot,*orientation,newRoot);
    const Feature* corrSta = FeatureRep::findCorrespondingFeature(oldRoot,*station,newRoot);
    assert(corrOri && corrSta);
    orientation = &Orientation::downcast(*corrOri);
    station     = &Station::downcast(*corrSta);
}


const Feature*
FramePlacementRep::findAncestorFeature(const Feature& root) const {
    assert(orientation && station);
    const Feature* ancestor = 
        FeatureRep::findYoungestCommonAncestor(*orientation,*station);
    return (ancestor && FeatureRep::isFeatureInFeatureTree(root, *ancestor))
            ? ancestor : 0;
}

    // REAL PLACEMENT REP //

Placement
RealPlacementRep::mul(const Placement& r) const {
    switch(r.getRep().getPlacementType()) {
    case RealPlacementType: {
        RealPlacement x; // null rep
        RealExprPlacementRep::binop(x,RealOps::Times,getMyHandle(),r);
        return x;
    }
    case Vec3PlacementType: {
        Vec3Placement x;
        Vec3ExprPlacementRep::scale(x,getMyHandle(),r);
        return x;    
    }
    case StationPlacementType: 
    case DirectionPlacementType:
    {
        Vec3Placement x;
        Vec3ExprPlacementRep::scale(x,getMyHandle(),
                                      Vec3Placement::cast(r));
        return x; 
    }
    default:
        SIMTK_THROW3(Exception::InfixPlacementOperationNotAllowed,
                        getPlacementTypeName(getPlacementType()),
                        "*", getPlacementTypeName(r.getRep().getPlacementType()));
    };
    //NOTREACHED
    return Placement();
}

    // REAL EXPR PLACEMENT REP //

bool 
RealOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    case Negate: 
        return args.size()==1
            && args[0].getRep().getPlacementType()==RealPlacementType;
    case Plus:
    case Minus:
    case Times:
    case Divide:
        return args.size() == 2 
                && args[0].getRep().getPlacementType()==RealPlacementType
                && args[1].getRep().getPlacementType()==RealPlacementType;
    default:
        assert(false);
    }
    //NOTREACHED
    return false;
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::unop(RealPlacement& handle, RealOps::OpKind op,
                      const Placement& a) {
    std::vector<const Placement*> args(1);
    args[0] = &a;
    return new RealExprPlacementRep(handle, RealOps(op), args);
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::binop(RealPlacement& handle, RealOps::OpKind op,
                      const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new RealExprPlacementRep(handle, RealOps(op), args);
}

    // VEC3 PLACEMENT REP //

    // VEC3 EXPR PLACEMENT REP //
bool Vec3Ops::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    case Scale:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==RealPlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    case Cast:
        return args.size() == 1
            && (   args[0].getRep().getPlacementType()==DirectionPlacementType
                || args[0].getRep().getPlacementType()==StationPlacementType);
    case Plus:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==Vec3PlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    case Minus:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==Vec3PlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    default: 
        assert(false);
    }
    //NOTREACHED
    return false;
}

/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::scale(Vec3Placement& handle, 
                            const Placement& s, const Placement& v) {
    std::vector<const Placement*> args(2);
    args[0] = &s; args[1] = &v;
    return new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Scale), args);
}
/*static*/ Vec3ExprPlacementRep* 
Vec3ExprPlacementRep::plus(Vec3Placement& handle,
                           const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Plus), args);
}
/*static*/ Vec3ExprPlacementRep* 
Vec3ExprPlacementRep::minus(Vec3Placement& handle,
                            const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    return new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Minus), args);
}
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::cast(Vec3Placement& handle, const Placement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    return new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Cast), args);
}

    // STATION PLACEMENT REP //

    // STATION EXPR PLACEMENT REP //

bool 
StationOps::checkArgs(const std::vector<Placement>& args) const {

    switch (op) {
    case Plus:
    case Minus:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==StationPlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    case Cast:
        return args.size() == 1
            && args[0].getRep().getPlacementType()==Vec3PlacementType;
    default: 
        assert(false);
    }
    return false;
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::plus(StationPlacement& handle,
                              const StationPlacement& s, const Vec3Placement& d) {
    std::vector<const Placement*> args(2);
    args[0] = &s; args[1] = &d;
    return new StationExprPlacementRep(handle, StationOps(StationOps::Plus), args);
}
/*static*/ StationExprPlacementRep*
StationExprPlacementRep::minus(StationPlacement& handle,
                               const StationPlacement& s, const Vec3Placement& d) {
    std::vector<const Placement*> args(2);
    args[0] = &s; args[1] = &d;
    return new StationExprPlacementRep(handle, StationOps(StationOps::Minus), args);
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::cast(StationPlacement& handle, const Vec3Placement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    return new StationExprPlacementRep(handle, StationOps(StationOps::Cast), args);
}

    // DIRECTION PLACEMENT REP //

    // DIRECTION EXPR PLACEMENT REP //

bool DirectionOps::checkArgs(const std::vector<Placement>& args) const {
    switch (op) {
    case Normalize:
        return args.size() == 1
            && (   args[0].getRep().getPlacementType()==Vec3PlacementType
                || args[1].getRep().getPlacementType()==StationPlacementType);
    default: 
        assert(false);
    }
    return false;
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalize(DirectionPlacement& handle, const Vec3Placement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    return new DirectionExprPlacementRep(handle, DirectionOps(DirectionOps::Normalize), args);
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalize(DirectionPlacement& handle, const StationPlacement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    return new DirectionExprPlacementRep(handle, DirectionOps(DirectionOps::Normalize), args);
}

} // namespace simtk
