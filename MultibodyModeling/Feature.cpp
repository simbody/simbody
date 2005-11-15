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
 * Implementations of API-level multibody modeling objects for Simbody.
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

bool Feature::hasPlacement() const {
    return getRep().hasPlacement();
}
void Feature::place(const Placement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Feature::place", exc.getMessage());
    }
}

const Placement& Feature::getPlacement() const {
    assert(hasPlacement());
    return getRep().getPlacement();
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

    const size_t nSubfeatures = f.getNSubfeatures();
    const size_t nPlacement = f.getNPlacementExpressions();
    const std::string nextIndent = linePrefix + "    ";

    if (nSubfeatures) {
        s << endl << linePrefix << "  Subfeatures (" << nSubfeatures << "):";
        for (size_t i=0; i < nSubfeatures; ++i)
            s  << endl << nextIndent << f.getSubfeature(i).toString(nextIndent);
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

int Feature::getNSubfeatures() const
  { return getRep().getNSubfeatures(); }
const Feature& Feature::getSubfeature(int i) const
  { return getRep().getSubfeature(i); }
Feature& Feature::updSubfeature(int i)
  { return updRep().updSubfeature(i); }

// getXXX() methods
const Feature& Feature::getSubfeature(const String& n) const 
  { return getRep().getSubfeature(n); }
const RealParameter& Feature::getRealParameter(const String& n) const
  { return RealParameter::downcast(getSubfeature(n)); }
const StationParameter& Feature::getStationParameter(const String& n) const
  { return StationParameter::downcast(getSubfeature(n)); }
const Station& Feature::getStation(const String& n) const
  { return Station::downcast(getSubfeature(n)); }

// updXXX() methods
Feature& Feature::updSubfeature(const String& n)
  { return updRep().updSubfeature(n); }
RealParameter& Feature::updRealParameter(const String& n)
  { return RealParameter::downcast(updSubfeature(n)); }
StationParameter& Feature::updStationParameter(const String& n)
  { return StationParameter::downcast(updSubfeature(n)); }
Station& Feature::updStation(const String& n)
  { return Station::downcast(updSubfeature(n)); }

// addXXX() methods
RealParameter& Feature::addRealParameter(const String& n, const RealPlacement& p) {
    RealParameter& rp = RealParameter::downcast(updRep().addSubfeatureLike(RealParameter(n), n));
    if (p.hasRep()) rp.place(p);
    return rp;
}

StationParameter& Feature::addStationParameter(const String& n, const StationPlacement& p) {
    StationParameter& sp = StationParameter::downcast(updRep().addSubfeatureLike(StationParameter(n), n));
    if (p.hasRep()) sp.place(p);
    return sp;
}

RealMeasure& Feature::addRealMeasure(const String& n, const RealPlacement& p) {
    RealMeasure& rm = RealMeasure::downcast(updRep().addSubfeatureLike(RealMeasure(n), n));
    if (p.hasRep()) rm.place(p);
    return rm;
}
StationMeasure& Feature::addStationMeasure(const String& n, const StationPlacement& p) {
    StationMeasure& sm = StationMeasure::downcast(updRep().addSubfeatureLike(StationMeasure(n), n));
    if (p.hasRep()) sm.place(p);
    return sm;
}

Station& Feature::addStation(const String& n, const StationPlacement& p) {
    Station& s = Station::downcast(updRep().addSubfeatureLike(Station(n), n));
    if (p.hasRep()) s.place(p);
    return s;
}

Frame& Feature::addFrame(const String& n, const FramePlacement& p) {
    Frame& f = Frame::downcast(updRep().addSubfeatureLike(Frame(n), n));
    if (p.hasRep()) f.place(p);
    return f;
}
const Frame& Feature::getFrame(const String& nm) const {
    return Frame::downcast(getSubfeature(nm));
}
Frame& Feature::updFrame(const String& nm) {
    return Frame::downcast(updSubfeature(nm));
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

void Station::place(const StationPlacement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Station::place", exc.getMessage());
    }
}

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

void Direction::place(const DirectionPlacement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Direction::place", exc.getMessage());
    }
}

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

void Orientation::place(const OrientationPlacement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Orientation::place", exc.getMessage());
    }
}

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

void Frame::place(const FramePlacement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Frame::place", exc.getMessage());
    }
}

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



static int caseInsensitiveCompare(const std::string& key, const std::string& test) {
    const size_t minlen = std::min(key.size(), test.size());
    for (size_t i=0; i < minlen; ++i) {
        const int k = tolower(key[i]), t = tolower(test[i]);
        if (k < t) return -1;
        else if (k > t) return 1;
    }
    // caution -- size() is unsigned, don't get clever here
    if (key.size() > minlen) return 1;
    else if (test.size() > minlen) return -1;
    return 0;
}
