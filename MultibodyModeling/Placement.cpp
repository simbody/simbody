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
#include <iomanip>
#include <sstream>

namespace simtk {

    // PLACEMENT //

Placement::Placement(const Placement& src) : rep(0) { 
    if (src.rep) src.rep->cloneUnownedWithNewHandle(*this);
}
Placement& Placement::operator=(const Placement& src) {
    if (this != &src) {
        delete rep; rep=0;
        if (src.rep) src.rep->cloneUnownedWithNewHandle(*this);
    }
    return *this;
}
Placement::~Placement() {
    if (rep==0) return;
    assert(&rep->getMyHandle() == this);
    delete rep; rep=0;
}

Placement::Placement(const Feature& f) : rep(0) {
    rep = f.getRep().createFeatureReference(*this);
}
Placement::Placement(const Feature& f, int i) : rep(0) {
    rep = f.getRep().createFeatureReference(*this, i);
}
Placement::Placement(const Real& r) : rep(0) {
    rep = new RealConstantPlacementRep(reinterpret_cast<RealPlacement&>(*this), r);
}
Placement::Placement(const Vec3& v) : rep(0) {
    rep = new Vec3ConstantPlacementRep(reinterpret_cast<Vec3Placement&>(*this), v);
}
//Placement::Placement(const Mat33& m) : rep(0) {
//    rep = new Mat33ConstantPlacementRep(reinterpret_cast<RealPlacement&>(*this), r);
//}

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

bool Placement::isConstant() const {
    return rep && rep->isConstant();
}

bool Placement::dependsOn(const Feature& f) const {
    return rep && rep->dependsOn(f);
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
    if (hasOwner())
        s << getOwner().getFullName() << ":"
          << std::left << std::setw(2) << getIndexInOwner();
    else s << "NO OWNER";
    s << " " << rep->toString(linePrefix);
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Placement& p) {
    return o << p.toString() << std::endl;
}


// unary
Placement operator+(const Placement& f)          {return f;}
Placement operator-(const Placement& f)          {return f.getRep().negate();}
RealPlacement length(const Placement& f)         {return f.getRep().length();}
DirectionPlacement normalize(const Placement& f) {return f.getRep().normalize();}

// binary
Placement operator+(const Placement& l, const Placement& r) {return l.getRep().add(r);} 
Placement operator-(const Placement& l, const Placement& r) {return l.getRep().sub(r);}
Placement operator*(const Placement& l, const Placement& r) {return l.getRep().mul(r);}
Placement operator/(const Placement& l, const Placement& r) {return l.getRep().dvd(r);}

    // REAL PLACEMENT //

RealPlacement::RealPlacement(const Real& r) {
    rep = new RealConstantPlacementRep(*this,r);
}

RealPlacement::RealPlacement(const RealParameter& rp) {
    rep = rp.getRep().useFeatureAsRealPlacement(*this);
}

RealPlacement::RealPlacement(const RealMeasure& rm) {
    rep = rm.getRep().useFeatureAsRealPlacement(*this);
}

/*static*/ RealPlacement
RealPlacement::negate(const RealPlacement& a) {
    RealPlacement x; // null rep
    RealExprPlacementRep::unop(x,RealOps::Negate,a);
    return x;
}
/*static*/ RealPlacement
RealPlacement::plus(const RealPlacement& l, const RealPlacement& r) {
    RealPlacement x; // null rep
    RealExprPlacementRep::binop(x,RealOps::Plus,l,r);
    return x;
}
/*static*/ RealPlacement
RealPlacement::minus(const RealPlacement& l, const RealPlacement& r) {
    RealPlacement x; // null rep
    RealExprPlacementRep::binop(x,RealOps::Minus,l,r);
    return x;
}
/*static*/ RealPlacement
RealPlacement::times(const RealPlacement& l, const RealPlacement& r) {
    RealPlacement x; // null rep
    RealExprPlacementRep::binop(x,RealOps::Times,l,r);
    return x;
}
/*static*/ RealPlacement
RealPlacement::divide(const RealPlacement& l, const RealPlacement& r) {
    RealPlacement x; // null rep
    RealExprPlacementRep::binop(x,RealOps::Divide,l,r);
    return x;
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
    rep = new Vec3ConstantPlacementRep(*this,r);
}

Vec3Placement::Vec3Placement(const Vec3Parameter& rp) {
    rep = rp.getRep().useFeatureAsVec3Placement(*this);
}

Vec3Placement::Vec3Placement(const Vec3Measure& rm) {
    rep = rm.getRep().useFeatureAsVec3Placement(*this);
}

/*static*/ Vec3Placement
Vec3Placement::plus(const Placement& l, const Placement& r) {
    Vec3Placement x; // null rep
    Vec3ExprPlacementRep::plus(x,l,r);
    return x;
}
/*static*/ Vec3Placement
Vec3Placement::minus(const Placement& l, const Placement& r) {
    Vec3Placement x; // null rep
    Vec3ExprPlacementRep::minus(x,l,r);
    return x;
}
/*static*/ Vec3Placement
Vec3Placement::scale(const Placement& s, const Placement& v) {
    Vec3Placement x; // null rep
    Vec3ExprPlacementRep::scale(x,s,v);
    return x;
}

/*static*/ Vec3Placement
Vec3Placement::cast(const Placement& v) {
    Vec3Placement x; // null rep
    Vec3ExprPlacementRep::cast(x,v);
    return x;
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
    rep = new StationConstantPlacementRep(*this,v);
}
StationPlacement::StationPlacement(const Frame& f) {
    rep = f.getOrigin().getRep().useFeatureAsStationPlacement(*this);
}
StationPlacement::StationPlacement(const Feature& f) {
    rep = f.getRep().useFeatureAsStationPlacement(*this);
}
/*static*/ StationPlacement
StationPlacement::cast(const Vec3Placement& v) {
    StationPlacement x; // null rep
    StationExprPlacementRep::cast(x,v);
    return x;
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
    rep = new DirectionConstantPlacementRep(*this,v);
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

/*static*/ DirectionPlacement
DirectionPlacement::normalize(const Vec3Placement& v) {
    DirectionPlacement x; // null rep
    DirectionExprPlacementRep::normalize(x,v);
    return x;
}

/*static*/ DirectionPlacement
DirectionPlacement::normalize(const StationPlacement& v) {
    DirectionPlacement x; // null rep
    DirectionExprPlacementRep::normalize(x,v);
    return x;
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
    rep = new OrientationConstantPlacementRep(*this,m);
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
    rep = new FrameExprPlacementRep(*this,o,s);
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
