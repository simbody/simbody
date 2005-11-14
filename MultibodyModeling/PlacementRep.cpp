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
    return (*feature).getRep().getRequiredPlacementType();
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

    // REAL EXPR PLACEMENT REP //

const Feature* 
RealExprPlacementRep::findAncestorFeature(const Feature& root) const {
    const Feature* ancestor = 0;
    bool foundNonConst = false;
    for (size_t i=0; i < args.size(); ++i) {
        assert(args[i]);
        if (args[i]->isConstant())
            continue;
        foundNonConst = true;
        const Feature* argAncestor = 
            args[i]->getRep().findAncestorFeature(root);
        if (ancestor && argAncestor)
            ancestor = FeatureRep::findYoungestCommonAncestor(*ancestor,*argAncestor);
        else ancestor = argAncestor;
    }
    assert(foundNonConst);
    return ancestor;
}

} // namespace simtk
