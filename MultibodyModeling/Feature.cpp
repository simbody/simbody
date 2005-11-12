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
        f.rep->cloneWithoutExternalPlacements(*this);
}
Feature& Feature::operator=(const Feature& f) {
    if (this != &f) {
        delete rep; rep=0;
        if (f.rep) 
            f.rep->cloneWithoutExternalPlacements(*this);
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
    updRep().setPlacement(p);
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

    const size_t nChildren = f.getNChildFeatures();
    const size_t nPlacement = f.getNPlacementExpressions();
    const std::string nextIndent = linePrefix + "    ";

    if (nChildren) {
        s << endl << linePrefix << "  Child Features (" << nChildren << "):";
        for (size_t i=0; i < nChildren; ++i)
            s  << endl << nextIndent << f.getChildFeature(i).toString(nextIndent);
    }
    if (nPlacement) {
        s << endl << linePrefix << "  Placement Expressions (" << nPlacement << "):";
        for (size_t i=0; i < nPlacement; ++i)
            s  << endl << nextIndent << f.getPlacementExpression(i).toString(nextIndent);
    }
    return s.str();
}

ostream& operator<<(ostream& o, const Feature& f) {
    return o << f.toString() << endl;
}
bool Feature::hasParentFeature() const {
    return getRep().hasParentFeature();
}
int Feature::getIndexInParent() const {
    assert(getRep().hasParentFeature());
    return getRep().getIndexInParent();
}
const Feature& Feature::getParentFeature() const {
    assert(getRep().hasParentFeature());
    return getRep().getParentFeature();
}
bool Feature::dependsOn(const Feature& f) const {
    if (isSameFeature(f)) return true;
    return hasRep() && getRep().dependsOn(f);
}

const Feature& Feature::getFeature(const String& n) const {
    const Feature* f = getRep().getChildFeature(std::string(n));
    if (!f) SIMTK_THROW2(Exception::FeatureNameNotFound,"Feature::getFeature",n);
    return *f;
}
Feature& Feature::updFeature(const String& n) {
    Feature* p = updRep().updChildFeature(std::string(n));
    if (!p) SIMTK_THROW2(Exception::FeatureNameNotFound,"Feature::updFeature",n);
    return *p;
}

const RealParameter& Feature::getRealParameter(const String& n) const {
    return RealParameter::downcast(getFeature(n));
}
RealParameter& Feature::updRealParameter(const String& n) {
    return RealParameter::downcast(updFeature(n));
}
RealParameter& Feature::addRealParameter(const String& n, const RealPlacement& p) {
    RealParameter& r = RealParameter::downcast(
                            updRep().addFeatureLike(RealParameter(n), n));
    if (PlacementRep::getRep(p)) {
        const RealPlacement& rp = RealPlacement::downcast(updRep().addPlacementLike(p));
        r.setPlacement(rp);
    }
    return r;
}

const StationParameter& Feature::getStationParameter(const String& n) const {
    return StationParameter::downcast(getFeature(n));
}
StationParameter& Feature::updStationParameter(const String& n) {
    return StationParameter::downcast(updFeature(n));
}
StationParameter& Feature::addStationParameter(const String& n, const StationPlacement& p) {
    return StationParameter::downcast(updRep().addFeatureLike(StationParameter(n), n));
}

RealMeasure& Feature::addRealMeasure(const String& n, const RealPlacement& p) {
    return RealMeasure::downcast(updRep().addFeatureLike(RealMeasure(n), n));
}
StationMeasure& Feature::addStationMeasure(const String& n, const StationPlacement& p) {
    return StationMeasure::downcast(updRep().addFeatureLike(StationMeasure(n), n));
}

Station& Feature::addStation(const String& n, const StationPlacement& p) {
    return Station::downcast(updRep().addFeatureLike(Station(n), n));
}
const Station& Feature::getStation(const String& nm) const {
    return Station::downcast(getFeature(nm));
}
Station& Feature::updStation(const String& nm) {
    return Station::downcast(updFeature(nm));
}
Frame& Feature::addFrame(const String& n, const FramePlacement& p) {
    return Frame::downcast(updRep().addFeatureLike(Frame(n), n));
}
const Frame& Feature::getFrame(const String& nm) const {
    return Frame::downcast(getFeature(nm));
}
Frame& Feature::updFrame(const String& nm) {
    return Frame::downcast(updFeature(nm));
}

    // REAL PARAMETER //
RealParameter::RealParameter(const String& nm)
  { (void)new RealParameterRep(*this, std::string(nm)); }
RealParameter::RealParameter(const RealParameter& src) : RealMeasure(src) { }
RealParameter& RealParameter::operator=(const RealParameter& src)
  { RealMeasure::operator=(src); return *this; }
RealParameter::~RealParameter() { }

/*static*/ bool             
RealParameter::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return RealParameterRep::isA(f.getRep());
}
/*static*/ const RealParameter& 
RealParameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealParameter&>(f);
}

/*static*/ RealParameter&       
RealParameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealParameter&>(f);
}

    // STATION PARAMETER //
StationParameter::StationParameter(const String& nm)
  { (void)new StationParameterRep(*this, std::string(nm)); }
StationParameter::StationParameter(const StationParameter& src) : StationMeasure(src) { }
StationParameter& StationParameter::operator=(const StationParameter& src)
  { StationMeasure::operator=(src); return *this; }
StationParameter::~StationParameter() { }

/*static*/ bool             
StationParameter::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationParameterRep::isA(f.getRep());
}
/*static*/ const StationParameter& 
StationParameter::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationParameter&>(f);
}

/*static*/ StationParameter&       
StationParameter::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationParameter&>(f);
}

    // REAL MEASURE //
RealMeasure::RealMeasure(const String& nm)
  { (void)new RealMeasureRep(*this, std::string(nm)); }
RealMeasure::RealMeasure(const RealMeasure& src) : Feature(src) { }
RealMeasure& RealMeasure::operator=(const RealMeasure& src)
  { Feature::operator=(src); return *this; }
RealMeasure::~RealMeasure() { }

/*static*/ bool             
RealMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return RealMeasureRep::isA(f.getRep());
}
/*static*/ const RealMeasure& 
RealMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealMeasure&>(f);
}

/*static*/ RealMeasure&       
RealMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealMeasure&>(f);
}

    // STATION MEASURE //
StationMeasure::StationMeasure(const String& nm)
  { (void)new StationMeasureRep(*this, std::string(nm)); }
StationMeasure::StationMeasure(const StationMeasure& src) : Feature(src) { }
StationMeasure& StationMeasure::operator=(const StationMeasure& src)
  { Feature::operator=(src); return *this; }
StationMeasure::~StationMeasure() { }

/*static*/ bool             
StationMeasure::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationMeasureRep::isA(f.getRep());
}
/*static*/ const StationMeasure& 
StationMeasure::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationMeasure&>(f);
}

/*static*/ StationMeasure&       
StationMeasure::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationMeasure&>(f);
}

    // STATION //
Station::Station(const String& nm)
  { (void)new StationRep(*this, std::string(nm)); }
Station::Station(const Station& src) : Feature(src) { }
Station& Station::operator=(const Station& src)
  { Feature::operator=(src); return *this; }
Station::~Station() { }

/*static*/ bool             
Station::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return StationRep::isA(f.getRep());
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
  { (void)new DirectionRep(*this, std::string(nm)); }
Direction::Direction(const Direction& src) : Feature(src) { }
Direction& Direction::operator=(const Direction& src)
  { Feature::operator=(src); return *this; }
Direction::~Direction() { }

/*static*/ bool             
Direction::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return DirectionRep::isA(f.getRep());
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
  { (void)new OrientationRep(*this, std::string(nm)); }
Orientation::Orientation(const Orientation& src) : Feature(src) { }
Orientation& Orientation::operator=(const Orientation& src)
  { Feature::operator=(src); return *this; }
Orientation::~Orientation() { }

/*static*/ bool             
Orientation::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return OrientationRep::isA(f.getRep());
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
  { (void)new FrameRep(*this, std::string(nm)); }
Frame::Frame(const Frame& src) : Feature(src) { }
Frame& Frame::operator=(const Frame& src)
  { Feature::operator=(src); return *this; }
Frame::~Frame() { }

const Orientation& Frame::getOrientation() const {
    return FrameRep::downcast(getRep()).getOrientation();
}
const Station& Frame::getOrigin() const {
    return FrameRep::downcast(getRep()).getOrigin();
}

/*static*/ bool             
Frame::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return FrameRep::isA(f.getRep());
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
