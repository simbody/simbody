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
 * Implementation of Placement handles.
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

    // PLACEMENT //
Placement::Placement(const Placement& src) : rep(0) { 
    if (src.rep) src.rep->clone(*this);
}
Placement& Placement::operator=(const Placement& src) {
    if (this != &src) {
        delete rep; rep=0;
        if (src.rep) src.rep->clone(*this);
    }
    return *this;
}
Placement::~Placement() {
    if (rep==0) return;
    assert(&rep->getMyHandle() == this);
    delete rep; rep=0;
}
bool Placement::hasOwner() const {
    return rep && rep->hasOwner();
}
int Placement::getIndexInOwner() const {
    assert(rep && rep->hasOwner());
    assert(&rep->getMyHandle() == this);
    return rep->getIndexInOwner();
}

const Feature& Placement::getOwner() const {
    assert(rep && rep->hasOwner());
    assert(&rep->getMyHandle() == this);
    return rep->getOwner();
}
String Placement::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Placement ";
    if (!rep) {
        s << "at 0x" << this << " HAS NULL REP";
        return s.str();
    }
    if (&rep->getMyHandle() != this) {
        s << "at 0x" << this << " HAS MISMATCHED REP";
        return s.str();
    }
    s << rep->toString(linePrefix);
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Placement& p) {
    return o << p.toString() << endl;
}

    // FEATURE PLACEMENT //
FeaturePlacement::FeaturePlacement(const Feature& f) {
    rep = new FeaturePlacementRep(*this,f);
}


    // REAL PLACEMENT //
RealPlacement::RealPlacement(const Real& r) {
    rep = new RealConstantPlacementRep(*this,r);
}

    // STATION PLACEMENT //
StationPlacement::StationPlacement(const Vec3& v) {
    rep = new StationConstantPlacementRep(*this,v);
}

    // DIRECTION PLACEMENT //
DirectionPlacement::DirectionPlacement(const Vec3& v) {
    rep = new DirectionConstantPlacementRep(*this,v);
}

    // ORIENTATION PLACEMENT //
OrientationPlacement::OrientationPlacement(const Mat33& m) {
    rep = new OrientationConstantPlacementRep(*this,m);
}

    // FRAME PLACEMENT //
FramePlacement::FramePlacement(const Orientation& o, const Station& s) {
    rep = new FramePlacementRep(*this,o,s);
}

FramePlacement::FramePlacement(const Station& s) {
    assert(s.hasParentFeature());
    assert(Frame::isInstanceOf(s.getParentFeature()));
    const Frame& parentFrame = Frame::downcast(s.getParentFeature());
    rep = new FramePlacementRep(*this,parentFrame.getOrientation(),s);
}

    // FEATURE PLACEMENT REP //
PlacementType FeaturePlacementRep::getPlacementType() const { 
    return FeatureRep::getRep(*feature)->getRequiredPlacementType();
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

} // namespace simtk
