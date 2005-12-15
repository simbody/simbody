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
 * Implementation of BasicPlacement handle classes.
 */

#include "SimbodyCommon.h"
#include "Placement.h"
#include "BasicPlacements.h"
#include "Feature.h"
#include "PlacementRep.h"
#include "FeatureRep.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace simtk {

    // REAL PLACEMENT //

RealPlacement::RealPlacement(const Real& r) {
    rep = new RealConstantPlacementRep(r);
    rep->setMyHandle(*this);
}

RealPlacement::RealPlacement(const RealParameter& rp) {
    rep = rp.getRep().useFeatureAsRealPlacement(*this);
}

RealPlacement::RealPlacement(const RealMeasure& rm) {
    rep = rm.getRep().useFeatureAsRealPlacement(*this);
}

RealPlacement::RealPlacement(const Feature& f) {
    rep = f.getRep().useFeatureAsRealPlacement(*this);
}


/*static*/ bool             
RealPlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return RealPlacementRep::isA(p.getRep());
}
/*static*/ const RealPlacement& 
RealPlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const RealPlacement&>(p);
}

/*static*/ RealPlacement&       
RealPlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<RealPlacement&>(p);
}

    // VEC3 PLACEMENT //
Vec3Placement::Vec3Placement(const Vec3& r) {
    rep = new Vec3ConstantPlacementRep(r);
    rep->setMyHandle(*this);
}

Vec3Placement::Vec3Placement(const Vec3Parameter& rp) {
    rep = rp.getRep().useFeatureAsVec3Placement(*this);
}

Vec3Placement::Vec3Placement(const Vec3Measure& rm) {
    rep = rm.getRep().useFeatureAsVec3Placement(*this);
}

/*static*/ bool             
Vec3Placement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return Vec3PlacementRep::isA(p.getRep());
}
/*static*/ const Vec3Placement& 
Vec3Placement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const Vec3Placement&>(p);
}

/*static*/ Vec3Placement&       
Vec3Placement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<Vec3Placement&>(p);
}

    // STATION PLACEMENT //

StationPlacement::StationPlacement(const Station& s) {
    rep = s.getRep().useFeatureAsStationPlacement(*this);
}
StationPlacement::StationPlacement(const StationMeasure& s) {
    rep = s.getRep().useFeatureAsStationPlacement(*this);
}
StationPlacement::StationPlacement(const StationParameter& s) {
    rep = s.getRep().useFeatureAsStationPlacement(*this);
}
StationPlacement::StationPlacement(const Vec3& v) {
    rep = new StationConstantPlacementRep(v);
    rep->setMyHandle(*this);
}
StationPlacement::StationPlacement(const Frame& f) {
    rep = f.getOrigin().getRep().useFeatureAsStationPlacement(*this);
}
StationPlacement::StationPlacement(const Feature& f) {
    rep = f.getRep().useFeatureAsStationPlacement(*this);
}


/*static*/ bool             
StationPlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return StationPlacementRep::isA(p.getRep());
}
/*static*/ const StationPlacement& 
StationPlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const StationPlacement&>(p);
}

/*static*/ StationPlacement&       
StationPlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<StationPlacement&>(p);
}

    // DIRECTION PLACEMENT //
DirectionPlacement::DirectionPlacement(const Direction& d) {
    rep = d.getRep().useFeatureAsDirectionPlacement(*this);
}
DirectionPlacement::DirectionPlacement(const DirectionMeasure& d) {
    rep = d.getRep().useFeatureAsDirectionPlacement(*this);
}
DirectionPlacement::DirectionPlacement(const Vec3& v) {
    rep = new DirectionConstantPlacementRep(v);
    rep->setMyHandle(*this);
}
DirectionPlacement::DirectionPlacement(const Feature& f) {
    rep = f.getRep().useFeatureAsDirectionPlacement(*this);
}
DirectionPlacement::DirectionPlacement(const Orientation& o, int i) {
    rep = o.getAxis(i).getRep().useFeatureAsDirectionPlacement(*this);
}
DirectionPlacement::DirectionPlacement(const Frame& f, int i) {
    rep = f.getAxis(i).getRep().useFeatureAsDirectionPlacement(*this);
}

/*static*/ bool             
DirectionPlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return DirectionPlacementRep::isA(p.getRep());
}
/*static*/ const DirectionPlacement& 
DirectionPlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const DirectionPlacement&>(p);
}

/*static*/ DirectionPlacement&       
DirectionPlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<DirectionPlacement&>(p);
}

    // ORIENTATION PLACEMENT //

OrientationPlacement::OrientationPlacement(const Orientation& o) {
    rep = o.getRep().useFeatureAsOrientationPlacement(*this);
}
OrientationPlacement::OrientationPlacement(const OrientationMeasure& om) {
    rep = om.getRep().useFeatureAsOrientationPlacement(*this);
}
OrientationPlacement::OrientationPlacement(const Mat33& m) {
    rep = new OrientationConstantPlacementRep(m);
    rep->setMyHandle(*this);
}
OrientationPlacement::OrientationPlacement(const Frame& f) {
    rep = f.getOrientation().getRep().useFeatureAsOrientationPlacement(*this);
}
OrientationPlacement::OrientationPlacement(const Feature& f) {
    rep = f.getRep().useFeatureAsOrientationPlacement(*this);
}
/*static*/ bool             
OrientationPlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return OrientationPlacementRep::isA(p.getRep());
}
/*static*/ const OrientationPlacement& 
OrientationPlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const OrientationPlacement&>(p);
}

/*static*/ OrientationPlacement&       
OrientationPlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<OrientationPlacement&>(p);
}

    // FRAME PLACEMENT //
FramePlacement::FramePlacement(const Orientation& o, const Station& s) {
    rep = new FrameExprPlacementRep(o,s);
    rep->setMyHandle(*this);
}

FramePlacement::FramePlacement(const Frame& f) {
    rep = f.getRep().useFeatureAsFramePlacement(*this);
}

FramePlacement::FramePlacement(const Station& s) {
    rep = s.getRep().useFeatureAsFramePlacement(*this);
}

FramePlacement::FramePlacement(const Orientation& o) {
    rep = o.getRep().useFeatureAsFramePlacement(*this);
}

/*static*/ bool             
FramePlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return FramePlacementRep::isA(p.getRep());
}
/*static*/ const FramePlacement& 
FramePlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const FramePlacement&>(p);
}

/*static*/ FramePlacement&       
FramePlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<FramePlacement&>(p);
}

} // namespace simtk
