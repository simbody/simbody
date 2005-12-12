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
 * Implementations of MassElement Features for Simbody.
 */

#include "SimbodyCommon.h"
#include "MassElement.h"
#include "MassElementRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::endl;
using std::ostream;

namespace simtk {

    // MASS ELEMENT //
const RealMeasure& MassElement::getMassMeasure() const {
    return MassElementRep::downcast(getRep()).getMassMeasure();
}
const StationMeasure& MassElement::getCentroidMeasure() const {
    return MassElementRep::downcast(getRep()).getCentroidMeasure();
}

/*static*/ bool             
MassElement::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return MassElementRep::isA(s.getRep());
}
/*static*/ const MassElement& 
MassElement::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MassElement&>(s);
}

/*static*/ MassElement&       
MassElement::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MassElement&>(s);
}

    // POINT MASS ELEMENT //

PointMassElement::PointMassElement(const String& nm) : MassElement() {
    rep = new PointMassElementRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
PointMassElement::PointMassElement(const PointMassElement& src)
  : MassElement(src) { }
PointMassElement& PointMassElement::operator=(const PointMassElement& src) {
    MassElement::operator=(src); return *this;
}
PointMassElement::~PointMassElement() { }

PointMassElement::PointMassElement(const String& nm, const Real& m) {
    rep = new PointMassElementRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
    setMass(m);
}
void PointMassElement::setMass(const Real& m) {
    PointMassElementRep::downcast(updRep()).setMass(m);
}

/*static*/ bool             
PointMassElement::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return PointMassElementRep::isA(s.getRep());
}
/*static*/ const PointMassElement& 
PointMassElement::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const PointMassElement&>(s);
}

/*static*/ PointMassElement&       
PointMassElement::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<PointMassElement&>(s);
}

    // SPHERE MASS ELEMENT //

    // CYLINDER MASS ELEMENT //

CylinderMassElement::CylinderMassElement(const String& nm) : MassElement() {
    rep = new CylinderMassElementRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
CylinderMassElement::CylinderMassElement(const CylinderMassElement& src)
  : MassElement(src) { }
CylinderMassElement& CylinderMassElement::operator=(const CylinderMassElement& src) {
    MassElement::operator=(src); return *this;
}
CylinderMassElement::~CylinderMassElement() { }

void CylinderMassElement::setMass(const Real& m) {
    CylinderMassElementRep::downcast(updRep()).setMass(m);
}
void CylinderMassElement::setRadius(const Real& r) {
    CylinderMassElementRep::downcast(updRep()).setRadius(r);
}
void CylinderMassElement::setHalfLength(const Real& h) {
    CylinderMassElementRep::downcast(updRep()).setHalfLength(h);
}
void CylinderMassElement::placeCenter(const Vec3& c) {
    CylinderMassElementRep::downcast(updRep()).placeCenter(c);
}
void CylinderMassElement::placeAxis(const Vec3& a) {
    CylinderMassElementRep::downcast(updRep()).placeAxis(a);
}

/*static*/ bool             
CylinderMassElement::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return CylinderMassElementRep::isA(s.getRep());
}
/*static*/ const CylinderMassElement& 
CylinderMassElement::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const CylinderMassElement&>(s);
}

/*static*/ CylinderMassElement&       
CylinderMassElement::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<CylinderMassElement&>(s);
}
} // namespace simtk
