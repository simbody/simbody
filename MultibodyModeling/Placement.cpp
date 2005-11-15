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
    if (src.rep) src.rep->cloneWithNewHandle(*this);
}
Placement& Placement::operator=(const Placement& src) {
    if (this != &src) {
        delete rep; rep=0;
        if (src.rep) src.rep->cloneWithNewHandle(*this);
    }
    return *this;
}
Placement::~Placement() {
    if (rep==0) return;
    assert(&rep->getMyHandle() == this);
    delete rep; rep=0;
}

Placement::Placement(const Feature& f) : rep(0) {
    (void)new FeaturePlacementRep(
        reinterpret_cast<FeaturePlacement&>(*this),f);
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
        s << std::left << std::setw(2) << getIndexInOwner() 
          << " " << getOwner().getFullName();
    else s << "NO OWNER";
    s << ":" << rep->toString(linePrefix);
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Placement& p) {
    return o << p.toString() << std::endl;
}

    // FEATURE PLACEMENT //
FeaturePlacement::FeaturePlacement(const Feature& f) {
    (void)new FeaturePlacementRep(*this,f);
}

FeaturePlacement::FeaturePlacement(const Feature& f, int index) {
    (void)new FeaturePlacementRep(*this,f,index);
}


    // REAL PLACEMENT //
RealPlacement::RealPlacement(const Real& r) {
    (void)new RealConstantPlacementRep(*this,r);
}

RealPlacement::RealPlacement(const Feature& f) {
    f.getRep().useAsRealPlacement(*this);
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

    // STATION PLACEMENT //
StationPlacement::StationPlacement(const Vec3& v) {
    (void)new StationConstantPlacementRep(*this,v);
}

StationPlacement::StationPlacement(const Station& s) {
    s.getRep().useAsStationPlacement(*this);
}

StationPlacement::StationPlacement(const Feature& f) {
    f.getRep().useAsStationPlacement(*this);
}

    // DIRECTION PLACEMENT //
DirectionPlacement::DirectionPlacement(const Vec3& v) {
    (void)new DirectionConstantPlacementRep(*this,v);
}

    // ORIENTATION PLACEMENT //
OrientationPlacement::OrientationPlacement(const Mat33& m) {
    (void)new OrientationConstantPlacementRep(*this,m);
}

    // FRAME PLACEMENT //
FramePlacement::FramePlacement(const Orientation& o, const Station& s) {
    (void)new FramePlacementRep(*this,o,s);
}

FramePlacement::FramePlacement(const Station& s) {
    s.getRep().useAsFramePlacement(*this);
}

FramePlacement::FramePlacement(const Orientation& o) {
    o.getRep().useAsFramePlacement(*this);
}

FramePlacement::FramePlacement(const Feature& f) {
    f.getRep().useAsFramePlacement(*this);
}

} // namespace simtk
