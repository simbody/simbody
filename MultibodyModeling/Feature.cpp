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
 * Implementations of high level multibody modeling objects for Simbody.
 */

#include "SimbodyCommon.h"
#include "Feature.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::endl;
using std::ostream;

namespace simtk {

    // FEATURE //

Feature::Feature(const Feature& f) : rep(0) {
    if (f.rep) 
        rep = f.rep->cloneWithoutPlacement(*this);
}
Feature& Feature::operator=(const Feature& f) {
    if (this != &f) {
        delete rep; 
        rep = f.rep ? f.rep->cloneWithoutPlacement(*this) : 0;
    }
    return *this;
}
Feature::~Feature() { 
    delete rep;
}
static String featureHasNoRep(const Feature& f) {
    std::ostringstream s;
    s << "<FEATURE AT 0x" << &f << " WITH NULL REP>";
    return String(s.str());
}
String Feature::getName() const {
    return rep ? String(rep->getName()) : featureHasNoRep(*this);
}
String Feature::getFullName() const {
    return rep ? String(rep->getFullName()) : featureHasNoRep(*this);
}
String Feature::getFeatureTypeName() const {
    return rep ? String(rep->getFeatureTypeName()) 
               : featureHasNoRep(*this);
}


ostream& operator<<(ostream& o, const Feature& f) {
    const FeatureRep* r = FeatureRep::getRep(f);
    o << "Feature ";
    if (!r)
        return o << featureHasNoRep(f) << endl;

    o << r->getFeatureTypeName() << " " << r->getFullName() << ": " 
        << r->getPlacement() << endl;
    o << "  Child Features:" << endl;
    for (size_t i=0; i<r->getChildFeatures().size(); ++i)
        o << "    " << r->getChildFeatures()[i];
    o << endl;
    return o;
}

int Feature::getIndexInParent() const {
    assert(rep && rep->hasParentFeature());
    return rep->getIndexInParent();
}
const Feature& Feature::getParentFeature() const {
    assert(rep && rep->hasParentFeature());
    return rep->getParentFeature();
}
const Feature& Feature::getFeature(const String& n) const {
    assert(rep);
    const Feature* f = rep->getChildFeature(std::string(n));
    if (!f) SIMTK_THROW2(Exception::FeatureNameNotFound,"Feature::getFeature",n);
    return *f;
}
Feature& Feature::updFeature(const String& n) {
    assert(rep);
    Feature* p = rep->updChildFeature(std::string(n));
    if (!p) SIMTK_THROW2(Exception::FeatureNameNotFound,"Feature::updFeature",n);
    return *p;
}
const Parameter& Feature::getParameter(const String& n) const {
    return Parameter::downcast(getFeature(n));
}
Parameter& Feature::updParameter(const String& n) {
    return Parameter::downcast(updFeature(n));
}
Parameter& Feature::addParameter(const String& n) {
    assert(rep); return Parameter::downcast(rep->addFeatureLike(Parameter(n), n));
}
Measure& Feature::addMeasure(const String& n) {
    assert(rep); return Measure::downcast(rep->addFeatureLike(Measure(n), n));
}
Station& Feature::addStation(const String& n) {
    assert(rep); return Station::downcast(rep->addFeatureLike(Station(n), n));
}

    // PARAMETER //
Parameter::Parameter(const String& nm)
  { rep = new ParameterRep(*this, std::string(nm)); }
Parameter::Parameter(const Parameter& src) : Feature(src) { }
Parameter& Parameter::operator=(const Parameter& src)
  { Feature::operator=(src); return *this; }
Parameter::~Parameter() { }

    // MEASURE //
Measure::Measure(const String& nm)
  { rep = new MeasureRep(*this, std::string(nm)); }
Measure::Measure(const Measure& src) : Feature(src) { }
Measure& Measure::operator=(const Measure& src)
  { Feature::operator=(src); return *this; }
Measure::~Measure() { }

    // STATION //
Station::Station(const String& nm)
  { rep = new StationRep(*this, std::string(nm)); }
Station::Station(const Station& src) : Feature(src) { }
Station& Station::operator=(const Station& src)
  { Feature::operator=(src); return *this; }
Station::~Station() { }

    // DIRECTION //
Direction::Direction(const String& nm)
  { rep = new DirectionRep(*this, std::string(nm)); }
Direction::Direction(const Direction& src) : Feature(src) { }
Direction& Direction::operator=(const Direction& src)
  { Feature::operator=(src); return *this; }
Direction::~Direction() { }

    // ORIENTATION //
Orientation::Orientation(const String& nm)
  { rep = new OrientationRep(*this, std::string(nm)); }
Orientation::Orientation(const Orientation& src) : Feature(src) { }
Orientation& Orientation::operator=(const Orientation& src)
  { Feature::operator=(src); return *this; }
Orientation::~Orientation() { }

    // FRAME //
Frame::Frame(const String& nm)
  { rep = new FrameRep(*this, std::string(nm)); }
Frame::Frame(const Frame& src) : Feature(src) { }
Frame& Frame::operator=(const Frame& src)
  { Feature::operator=(src); return *this; }
Frame::~Frame() { }


} // namespace simtk
