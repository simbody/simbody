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

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Placement.h"
#include "simbody/internal/BasicPlacements.h"
#include "simbody/internal/Feature.h"

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
    rep = rp.getRep().createFeatureReference(*this);
}

RealPlacement::RealPlacement(const RealMeasure& rm) {
    rep = rm.getRep().createFeatureReference(*this);
}

RealPlacement::RealPlacement(const Feature& f) {
    try {
        rep = RealPlacementRep::createRealPlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "RealPlacement", f.getFullName(), exc.getMessageText());
    }
}

RealPlacement::RealPlacement(const Placement& p) {
    try {
        rep = RealPlacementRep::createRealPlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "RealPlacement", "Placement", exc.getMessageText());
    }
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
    rep = rp.getRep().createFeatureReference(*this);
}

Vec3Placement::Vec3Placement(const Vec3Measure& rm) {
    rep = rm.getRep().createFeatureReference(*this);
}

Vec3Placement::Vec3Placement(const Feature& f) {
    try {
        rep = Vec3PlacementRep::createVec3PlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "Vec3Placement", f.getFullName(), exc.getMessageText());
    }
}

Vec3Placement::Vec3Placement(const Placement& p) {
    try {
        rep = Vec3PlacementRep::createVec3PlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "Vec3Placement", "Placement", exc.getMessageText());
    }
}

/*static*/bool Vec3Placement::canConvert(const Placement& p) {
    if (!p.hasRep()) return false;
    PlacementRep* converted = Vec3PlacementRep::createVec3PlacementFrom(p, true);
    if (!converted) return false;
    delete converted;
    return true;
}

/*static*/Vec3Placement Vec3Placement::convert(const Placement& p) {
    PlacementRep* converted = Vec3PlacementRep::createVec3PlacementFrom(p);
    assert(converted);
    return Vec3Placement(reinterpret_cast<Vec3PlacementRep*>(converted));
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
    rep = s.getRep().createFeatureReference(*this);
}
StationPlacement::StationPlacement(const StationMeasure& s) {
    rep = s.getRep().createFeatureReference(*this);
}
StationPlacement::StationPlacement(const StationParameter& s) {
    rep = s.getRep().createFeatureReference(*this);
}
StationPlacement::StationPlacement(const Vec3& v) {
    rep = new StationConstantPlacementRep(v);
    rep->setMyHandle(*this);
}
StationPlacement::StationPlacement(const FrameFeature& f) {
    rep = new StationFeaturePlacementRep(f.getOrigin());
    rep->setMyHandle(*this);
}
StationPlacement::StationPlacement(const Feature& f) {
    try {
        rep = StationPlacementRep::createStationPlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "StationPlacement", f.getFullName(), exc.getMessageText());
    }
}

StationPlacement::StationPlacement(const Placement& p) {
    try {
        rep = StationPlacementRep::createStationPlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "StationPlacement", "Placement", exc.getMessageText());
    }
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
    rep = d.getRep().createFeatureReference(*this);
}
DirectionPlacement::DirectionPlacement(const DirectionMeasure& d) {
    rep = d.getRep().createFeatureReference(*this);
}
DirectionPlacement::DirectionPlacement(const UnitVec3& v) {
    rep = new DirectionConstantPlacementRep(v);
    rep->setMyHandle(*this);
}
DirectionPlacement::DirectionPlacement(const Vec3& v) {
    rep = new DirectionConstantPlacementRep(UnitVec3(v));
    rep->setMyHandle(*this);
}
DirectionPlacement::DirectionPlacement(const Feature& f) {
    try {
        rep = DirectionPlacementRep::createDirectionPlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "DirectionPlacement", f.getFullName(), exc.getMessageText());
    }
}
DirectionPlacement::DirectionPlacement(const Placement& p) {
    try {
        rep = DirectionPlacementRep::createDirectionPlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "DirectionPlacement", "Placement", exc.getMessageText());
    }
}
DirectionPlacement::DirectionPlacement(const Orientation& o, int i) {
    rep = new DirectionFeaturePlacementRep(o.getAxis(i));
    rep->setMyHandle(*this);
}
DirectionPlacement::DirectionPlacement(const FrameFeature& f, int i) {
    rep = new DirectionFeaturePlacementRep(f.getAxis(i));
    rep->setMyHandle(*this);
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
    rep = o.getRep().createFeatureReference(*this);
}
OrientationPlacement::OrientationPlacement(const OrientationMeasure& om) {
    rep = om.getRep().createFeatureReference(*this);
}
OrientationPlacement::OrientationPlacement(const RotationMat& m) {
    rep = new OrientationConstantPlacementRep(m);
    rep->setMyHandle(*this);
}
OrientationPlacement::OrientationPlacement(const FrameFeature& f) {
    rep = new OrientationFeaturePlacementRep(f.getOrientation());
    rep->setMyHandle(*this);
}
OrientationPlacement::OrientationPlacement(const Feature& f) {
    try {
        rep = OrientationPlacementRep::createOrientationPlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "OrientationPlacement", f.getFullName(), exc.getMessageText());
    }
}

OrientationPlacement::OrientationPlacement(const Placement& p) {
    try {
        rep = OrientationPlacementRep::createOrientationPlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "OrientationPlacement", "Placement", exc.getMessageText());
    }
}

OrientationPlacement
OrientationPlacement::invert() const
{
    try {
        return getRep().invert();
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "OrientationPlacement::invert", "", exc.getMessageText());
    }
    //NOTREACHED
    return OrientationPlacement();
}


/*static*/ OrientationPlacement
OrientationPlacement::createFromZAxis(const DirectionPlacement& z)
{
    try {
        return OrientationPlacementRep::createFromZAxis(z);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "OrientationPlacement::createFromZAxis", "z", exc.getMessageText());
    }
    //NOTREACHED
    return OrientationPlacement();
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

    // INERTIA PLACEMENT //

InertiaPlacement::InertiaPlacement(const InertiaMeasure& im) {
    rep = im.getRep().createFeatureReference(*this);
}

InertiaPlacement::InertiaPlacement(const StationPlacement& p, const RealPlacement& m) {
    rep = InertiaExprPlacementRep::ptMassOp(p,m);
    rep->setMyHandle(*this);
}

InertiaPlacement::InertiaPlacement(const RealPlacement& Ixx, const RealPlacement& Iyy, const RealPlacement& Izz)
{
    rep = InertiaExprPlacementRep::principalMomentsOp(Ixx,Iyy,Izz);
    rep->setMyHandle(*this);
}

InertiaPlacement::InertiaPlacement(const RealPlacement& Ixx, const RealPlacement& Iyy, const RealPlacement& Izz,
                                   const RealPlacement& Ixy, const RealPlacement& Ixz, const RealPlacement& Iyz)
{
    rep = InertiaExprPlacementRep::fullInertiaOp(Ixx,Iyy,Izz,Ixy,Ixz,Iyz);
    rep->setMyHandle(*this);
}

InertiaPlacement::InertiaPlacement(const InertiaMat& m) {
    rep = new InertiaConstantPlacementRep(m);
    rep->setMyHandle(*this);
}

InertiaPlacement::InertiaPlacement(const Feature& f) {
    try {
        rep = InertiaPlacementRep::createInertiaPlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "InertiaPlacement", f.getFullName(), exc.getMessageText());
    }
}

InertiaPlacement::InertiaPlacement(const Placement& p) {
    try {
        rep = InertiaPlacementRep::createInertiaPlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "InertiaPlacement", "Placement", exc.getMessageText());
    }
}

InertiaPlacement
InertiaPlacement::changeAxes(const OrientationPlacement& r) const
{
    try {
        return getRep().changeAxes(r);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "InertiaPlacement::changeAxes", "", exc.getMessageText());
    }
    //NOTREACHED
    return InertiaPlacement();
}

InertiaPlacement
InertiaPlacement::shiftFromCOM(const StationPlacement& to,
                               const RealPlacement&    totalMass) const
{
    try {
        return getRep().shiftFromCOM(to,totalMass);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "InertiaPlacement::shiftFromCOM", "", exc.getMessageText());
    }
    //NOTREACHED
    return InertiaPlacement();
}

InertiaPlacement
InertiaPlacement::shiftToCOM(const StationPlacement& com,
                             const RealPlacement&    totalMass) const
{
    try {
        return getRep().shiftToCOM(com,totalMass);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "InertiaPlacement::shiftToCOM", "", exc.getMessageText());
    }
    //NOTREACHED
    return InertiaPlacement();
}

/*static*/ bool             
InertiaPlacement::isInstanceOf(const Placement& p) {
    if (!p.hasRep()) return false;
    return InertiaPlacementRep::isA(p.getRep());
}
/*static*/ const InertiaPlacement& 
InertiaPlacement::downcast(const Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<const InertiaPlacement&>(p);
}

/*static*/ InertiaPlacement&       
InertiaPlacement::downcast(Placement& p) {
    assert(isInstanceOf(p));
    return reinterpret_cast<InertiaPlacement&>(p);
}

    // FRAME PLACEMENT //
FramePlacement::FramePlacement(const OrientationPlacement& o, const StationPlacement& s) {
    rep = new FrameExprPlacementRep(o,s);
    rep->setMyHandle(*this);
}

FramePlacement::FramePlacement(const FrameFeature& f) {
    rep = f.getRep().createFeatureReference(*this);
}

FramePlacement::FramePlacement(const Frame& f) {
    rep = new FrameConstantPlacementRep(f);
    rep->setMyHandle(*this);
}

FramePlacement::FramePlacement(const Station& s) {
    try {
        rep = FramePlacementRep::createFramePlacementFrom(Placement(s));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "FramePlacement", s.getFullName(), exc.getMessageText());
    }
}

FramePlacement::FramePlacement(const Orientation& o) {
    try {
        rep = FramePlacementRep::createFramePlacementFrom(Placement(o));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "FramePlacement", o.getFullName(), exc.getMessageText());
    }
}

FramePlacement::FramePlacement(const Feature& f) {
    try {
        rep = FramePlacementRep::createFramePlacementFrom(Placement(f));
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "FramePlacement", f.getFullName(), exc.getMessageText());
    }
}

FramePlacement::FramePlacement(const Placement& p) {
    try {
        rep = FramePlacementRep::createFramePlacementFrom(p);
        rep->setMyHandle(*this);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW3(Exception::PlacementAPIMethodFailed,
            "FramePlacement", "Placement", exc.getMessageText());
    }
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
