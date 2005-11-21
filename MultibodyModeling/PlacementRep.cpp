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
    case Negate: 
        return args.size()==1
            && args[0].getRep().getPlacementType()==RealPlacementType;
    case Length: 
        return args.size()==1
            && args[0].getRep().getPlacementType()==Vec3PlacementType;    
    case Plus:
    case Minus:
    case Times:
    case Divide:
        return args.size() == 2 
                && args[0].getRep().getPlacementType()==RealPlacementType
                && args[1].getRep().getPlacementType()==RealPlacementType;
    case Distance:
        return args.size() == 2 
                && args[0].getRep().getPlacementType()==StationPlacementType
                && args[1].getRep().getPlacementType()==StationPlacementType;    
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
    RealExprPlacementRep* rep = new RealExprPlacementRep(handle, RealOps(op), args);
    handle.setRep(rep);
    return rep;
}

/*static*/ RealExprPlacementRep*
RealExprPlacementRep::binop(RealPlacement& handle, RealOps::OpKind op,
                      const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    RealExprPlacementRep* rep = new RealExprPlacementRep(handle, RealOps(op), args);
    handle.setRep(rep);
    return rep;
}

    // VEC3 PLACEMENT REP //

    // VEC3 FEATURE PLACEMENT REP //
Vec3 Vec3FeaturePlacementRep::getValue(/*State*/) const {
    const PlacementRep& p = getReferencedPlacement().getRep();
    Vec3 value = Vec3(NTraits<Real>::getNaN());
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
    case Scale:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==RealPlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    case Cast:
        return args.size() == 1
            && (   args[0].getRep().getPlacementType()==DirectionPlacementType
                || args[0].getRep().getPlacementType()==StationPlacementType);
    // v= v+v, but s+s not allowed
    case Plus:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==Vec3PlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;
    
    // v= v-v, s-s OK
    case Minus:
        return args.size() == 2
            && (  (   args[0].getRep().getPlacementType()==Vec3PlacementType
                   && args[1].getRep().getPlacementType()==Vec3PlacementType)
               || (   args[0].getRep().getPlacementType()==StationPlacementType
                   && args[1].getRep().getPlacementType()==StationPlacementType)
               );
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
    Vec3ExprPlacementRep* rep = new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Scale), args);
    handle.setRep(rep); 
    return rep;
}
/*static*/ Vec3ExprPlacementRep* 
Vec3ExprPlacementRep::plus(Vec3Placement& handle,
                           const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    Vec3ExprPlacementRep* rep = new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Plus), args);
    handle.setRep(rep); 
    return rep;
}
/*static*/ Vec3ExprPlacementRep* 
Vec3ExprPlacementRep::minus(Vec3Placement& handle,
                            const Placement& l, const Placement& r) {
    std::vector<const Placement*> args(2);
    args[0] = &l; args[1] = &r;
    Vec3ExprPlacementRep* rep = new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Minus), args);
    handle.setRep(rep); 
    return rep;
}
/*static*/ Vec3ExprPlacementRep*
Vec3ExprPlacementRep::cast(Vec3Placement& handle, const Placement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    Vec3ExprPlacementRep* rep = new Vec3ExprPlacementRep(handle, Vec3Ops(Vec3Ops::Cast), args);
    handle.setRep(rep); 
    return rep;
}

    // STATION PLACEMENT REP //

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
    case Cast:
        return args.size() == 1
            && args[0].getRep().getPlacementType()==Vec3PlacementType;
    case Plus:
    case Minus:
        return args.size() == 2
            && args[0].getRep().getPlacementType()==StationPlacementType
            && args[1].getRep().getPlacementType()==Vec3PlacementType;

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
    StationExprPlacementRep* rep = new StationExprPlacementRep(handle, StationOps(StationOps::Plus), args);
    handle.setRep(rep); 
    return rep;
}
/*static*/ StationExprPlacementRep*
StationExprPlacementRep::minus(StationPlacement& handle,
                               const StationPlacement& s, const Vec3Placement& d) {
    std::vector<const Placement*> args(2);
    args[0] = &s; args[1] = &d;
    StationExprPlacementRep* rep = new StationExprPlacementRep(handle, StationOps(StationOps::Minus), args);
    handle.setRep(rep); 
    return rep;
}

/*static*/ StationExprPlacementRep*
StationExprPlacementRep::cast(StationPlacement& handle, const Vec3Placement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    StationExprPlacementRep* rep = new StationExprPlacementRep(handle, StationOps(StationOps::Cast), args);
    handle.setRep(rep); 
    return rep;
}

    // DIRECTION PLACEMENT REP //

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
    DirectionExprPlacementRep* rep = new DirectionExprPlacementRep(handle, DirectionOps(DirectionOps::Normalize), args);
    handle.setRep(rep); 
    return rep;
}

/*static*/ DirectionExprPlacementRep*
DirectionExprPlacementRep::normalize(DirectionPlacement& handle, const StationPlacement& v) {
    std::vector<const Placement*> args(1);
    args[0] = &v;
    DirectionExprPlacementRep* rep = new DirectionExprPlacementRep(handle, DirectionOps(DirectionOps::Normalize), args);
    handle.setRep(rep); 
    return rep;
}

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

} // namespace simtk
