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
    assert(rep);
    return MassElementRep::downcast(*rep).getMassMeasure();
}
const StationMeasure& MassElement::getCentroidMeasure() const {
    assert(rep);
    return MassElementRep::downcast(*rep).getCentroidMeasure();
}

/*static*/ bool             
MassElement::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return MassElementRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const MassElement& 
MassElement::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const MassElement&>(f);
}

/*static*/ MassElement&       
MassElement::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<MassElement&>(f);
}

    // POINT MASS ELEMENT //

PointMassElement::PointMassElement(const String& nm) {
    rep = new PointMassElementRep(*this, std::string(nm));
}
PointMassElement::PointMassElement(const PointMassElement& src)
  : MassElement(src) { }
PointMassElement& PointMassElement::operator=(const PointMassElement& src) {
    MassElement::operator=(src); return *this;
}
PointMassElement::~PointMassElement() { }

PointMassElement::PointMassElement(const String& nm, const Real& m) {
    rep = new PointMassElementRep(*this, std::string(nm));
    PointMassElementRep::downcast(*rep).setMass(m);
}
void PointMassElement::setMass(const Real& m) {
    assert(rep);
    PointMassElementRep::downcast(*rep).setMass(m);
}
void PointMassElement::placePoint(const Vec3& v) {
    assert(rep);
    PointMassElementRep::downcast(*rep).placePoint(v);
}

/*static*/ bool             
PointMassElement::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return PointMassElementRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const PointMassElement& 
PointMassElement::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const PointMassElement&>(f);
}

/*static*/ PointMassElement&       
PointMassElement::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<PointMassElement&>(f);
}

    // SPHERE MASS ELEMENT //

    // CYLINDER MASS ELEMENT //

CylinderMassElement::CylinderMassElement(const String& nm) {
    rep = new CylinderMassElementRep(*this, std::string(nm));
}
CylinderMassElement::CylinderMassElement(const CylinderMassElement& src)
  : MassElement(src) { }
CylinderMassElement& CylinderMassElement::operator=(const CylinderMassElement& src) {
    MassElement::operator=(src); return *this;
}
CylinderMassElement::~CylinderMassElement() { }

void CylinderMassElement::setMass(const Real& m) {
    assert(rep);
    CylinderMassElementRep::downcast(*rep).setMass(m);
}
void CylinderMassElement::setRadius(const Real& r) {
    assert(rep);
    CylinderMassElementRep::downcast(*rep).setRadius(r);
}
void CylinderMassElement::setHalfLength(const Real& h) {
    assert(rep);
    CylinderMassElementRep::downcast(*rep).setHalfLength(h);
}
void CylinderMassElement::placeCenter(const Vec3& c) {
    assert(rep);
    CylinderMassElementRep::downcast(*rep).placeCenter(c);
}
void CylinderMassElement::placeAxis(const Vec3& a) {
    assert(rep);
    CylinderMassElementRep::downcast(*rep).placeAxis(a);
}

/*static*/ bool             
CylinderMassElement::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return CylinderMassElementRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const CylinderMassElement& 
CylinderMassElement::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const CylinderMassElement&>(f);
}

/*static*/ CylinderMassElement&       
CylinderMassElement::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<CylinderMassElement&>(f);
}
} // namespace simtk
