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

void Feature::setPlacement(const Placement& p) {
    assert(rep);
    rep->setPlacement(p);
}


String Feature::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Feature ";
    if (!rep) {
        s << featureHasNoRep(*this);
        return String(s.str());
    }

    const FeatureRep& f = *rep;
    s << f.getFeatureTypeName() << " " << f.getFullName() << ": ";
    s << (f.hasPlacement() ? f.getPlacement().toString(linePrefix)
                           : String("NO PLACEMENT"));
    s << endl; 

    const size_t nChildren = f.getChildFeatures().size();
    const size_t nPlacement = f.getPlacementExpressions().size();

    const std::string nextIndent = linePrefix + "    ";

    s << linePrefix << "  Child Features (" << nChildren << ")";
    if (nChildren) s << ":";
    for (size_t i=0; i < nChildren; ++i)
        s  << endl << nextIndent << f.getChildFeature(i).toString(nextIndent);
    s << endl;
    s << linePrefix << "  Placement Expressions (" << nPlacement << ")";
    if (nPlacement) s << ":";
    for (size_t i=0; i < nPlacement; ++i)
        s  << endl << nextIndent << f.getPlacementExpression(i).toString(nextIndent);
    return s.str();
}

ostream& operator<<(ostream& o, const Feature& f) {
    return o << f.toString() << endl;
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

/*static*/ bool             
Parameter::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return ParameterRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Parameter& 
Parameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Parameter&>(f);
}

/*static*/ Parameter&       
Parameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Parameter&>(f);
}

    // MEASURE //
Measure::Measure(const String& nm)
  { rep = new MeasureRep(*this, std::string(nm)); }
Measure::Measure(const Measure& src) : Feature(src) { }
Measure& Measure::operator=(const Measure& src)
  { Feature::operator=(src); return *this; }
Measure::~Measure() { }

/*static*/ bool             
Measure::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return MeasureRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Measure& 
Measure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Measure&>(f);
}

/*static*/ Measure&       
Measure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Measure&>(f);
}

    // STATION //
Station::Station(const String& nm)
  { rep = new StationRep(*this, std::string(nm)); }
Station::Station(const Station& src) : Feature(src) { }
Station& Station::operator=(const Station& src)
  { Feature::operator=(src); return *this; }
Station::~Station() { }

/*static*/ bool             
Station::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return StationRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Station& 
Station::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Station&>(f);
}

/*static*/ Station&       
Station::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Station&>(f);
}

    // DIRECTION //
Direction::Direction(const String& nm)
  { rep = new DirectionRep(*this, std::string(nm)); }
Direction::Direction(const Direction& src) : Feature(src) { }
Direction& Direction::operator=(const Direction& src)
  { Feature::operator=(src); return *this; }
Direction::~Direction() { }

/*static*/ bool             
Direction::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return DirectionRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Direction& 
Direction::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Direction&>(f);
}

/*static*/ Direction&       
Direction::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Direction&>(f);
}

    // ORIENTATION //
Orientation::Orientation(const String& nm)
  { rep = new OrientationRep(*this, std::string(nm)); }
Orientation::Orientation(const Orientation& src) : Feature(src) { }
Orientation& Orientation::operator=(const Orientation& src)
  { Feature::operator=(src); return *this; }
Orientation::~Orientation() { }

/*static*/ bool             
Orientation::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return OrientationRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Orientation& 
Orientation::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Orientation&>(f);
}

/*static*/ Orientation&       
Orientation::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Orientation&>(f);
}

    // FRAME //
Frame::Frame(const String& nm)
  { rep = new FrameRep(*this, std::string(nm)); }
Frame::Frame(const Frame& src) : Feature(src) { }
Frame& Frame::operator=(const Frame& src)
  { Feature::operator=(src); return *this; }
Frame::~Frame() { }

/*static*/ bool             
Frame::isInstanceOf(const Feature& f) {
    if (!FeatureRep::getRep(f)) return false;
    return FrameRep::isA(*FeatureRep::getRep(f));
}
/*static*/ const Frame& 
Frame::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Frame&>(f);
}

/*static*/ Frame&       
Frame::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Frame&>(f);
}

} // namespace simtk
