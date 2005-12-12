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


namespace simtk {

    // SUBSYSTEM //


Subsystem::Subsystem(const Subsystem& f) : rep(0) {
    if (f.rep) 
        f.rep->cloneWithoutParentOrExternalPlacements(*this);
}
Subsystem& Subsystem::operator=(const Subsystem& f) {
    if (this != &f) {
        // This will blow up if rep doesn't have a handle -- we shouldn't
        // be pointing to it in that case!
        if (rep && (&rep->getMyHandle() == this)) delete rep; 
        rep = 0;
        if (f.rep) 
            f.rep->cloneWithoutParentOrExternalPlacements(*this);
    }
    return *this;
}
Subsystem::~Subsystem() {
    // This will blow up if rep doesn't have a handle -- we shouldn't
    // be pointing to it in that case!
    if (rep && (&rep->getMyHandle() == this)) delete rep; 
    rep = 0;
}

// It's the same Subsystem only if (1) they both have a rep, and
// (2) both reps point to the same address.
bool Subsystem::isSameSubsystem(const Subsystem& s) const {
    return rep && (rep == s.rep);
}

static String subsystemHasNoRep(const Subsystem& f) {
    std::ostringstream s;
    s << "<FEATURE AT 0x" << &f << " WITH NULL REP>";
    return String(s.str());
}
String Subsystem::getName() const {
    return rep ? String(rep->getName()) : subsystemHasNoRep(*this);
}
String Subsystem::getFullName() const {
    return rep ? String(rep->getFullName()) : subsystemHasNoRep(*this);
}

void Subsystem::realize(/*State,*/ Stage g) const {
    try {
        getRep().realize(g);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Subsystem::realize()", exc.getMessage());
    }
}

const Placement&
Subsystem::getPlacement() const {
    if (!Feature::isInstanceOf(*this)) {
        SIMTK_THROW1(Exception::OnlyFeaturesHavePlacements, getFullName());
        //NOTREACHED
    }
    return Feature::downcast(*this).getPlacement();
}

const PlacementValue&
Subsystem::getValue() const {
    return getPlacement().getValue();
}

void Subsystem::place(const Placement& p) {
    if (!Feature::isInstanceOf(*this)) {
        SIMTK_THROW1(Exception::OnlyFeaturesHavePlacements, getFullName());
        //NOTREACHED
    }
    Feature::downcast(*this).place(p);
}

String Subsystem::toString(const String& linePrefix) const {
    std::stringstream s;
    s << "Subsystem ";
    if (!rep) {
        s << subsystemHasNoRep(*this);
        return String(s.str());
    }

    const SubsystemRep& sr = *rep;
    s << sr.getFullName() << ": ";

    if (FeatureRep::isA(sr)) {
        const FeatureRep& fr = FeatureRep::downcast(sr);
        s << "Feature " << fr.getFeatureTypeName() << " ";
        s << (fr.hasPlacement() ? fr.getPlacement().toString(linePrefix)
                                : String("NO PLACEMENT"));
    }

    const size_t nSubsystems      = sr.getNSubsystems();
    const size_t nPlacement       = sr.getNPlacementExpressions();
    const size_t nPlacementValues = sr.getNPlacementValues();
    const std::string nextIndent  = linePrefix + "    ";

    if (nSubsystems) {
        s << std::endl << linePrefix << "  Subsystems (" << nSubsystems << "):";
        for (size_t i=0; i < nSubsystems; ++i)
            s  << std::endl << nextIndent << sr.getSubsystem(i).toString(nextIndent);
    }
    if (nPlacement) {
        s << std::endl << linePrefix << "  Placement Expressions (" << nPlacement << "):";
        for (size_t i=0; i < nPlacement; ++i)
            s  << std::endl << nextIndent << sr.getPlacementExpression(i).toString(nextIndent);
    }
    if (nPlacementValues) {
        s << std::endl << linePrefix << "  Placement Values (" << nPlacementValues << "):";
        for (size_t i=0; i < nPlacementValues; ++i)
            s  << std::endl << nextIndent << sr.getPlacementValue(i).toString(nextIndent);
    }
    return s.str();
}

std::ostream& operator<<(std::ostream& o, const Subsystem& s) {
    return o << s.toString() << std::endl;
}
bool Subsystem::hasParentSubsystem() const {
    return getRep().hasParentSubsystem();
}
int Subsystem::getIndexInParent() const {
    assert(getRep().hasParentSubsystem());
    return getRep().getIndexInParent();
}
const Subsystem& Subsystem::getParentSubsystem() const {
    assert(getRep().hasParentSubsystem());
    return getRep().getParentSubsystem();
}


int Subsystem::getNSubsystems() const
  { return getRep().getNSubsystems(); }
const Subsystem& Subsystem::getSubsystem(int i) const
  { return getRep().getSubsystem(i); }
Subsystem& Subsystem::updSubsystem(int i)
  { return updRep().updSubsystem(i); }

// getXXX() methods
const Subsystem& Subsystem::getSubsystem(const String& n) const 
  { return getRep().getSubsystem(n); }
const Feature& Subsystem::getFeature(const String& n) const 
  { return Feature::downcast(getSubsystem(n)); }
const RealParameter& Subsystem::getRealParameter(const String& n) const
  { return RealParameter::downcast(getSubsystem(n)); }
const Vec3Parameter& Subsystem::getVec3Parameter(const String& n) const
  { return Vec3Parameter::downcast(getSubsystem(n)); }
const StationParameter& Subsystem::getStationParameter(const String& n) const
  { return StationParameter::downcast(getSubsystem(n)); }
const RealMeasure& Subsystem::getRealMeasure(const String& n) const
  { return RealMeasure::downcast(getSubsystem(n)); }
const Vec3Measure& Subsystem::getVec3Measure(const String& n) const
  { return Vec3Measure::downcast(getSubsystem(n)); }
const StationMeasure& Subsystem::getStationMeasure(const String& n) const
  { return StationMeasure::downcast(getSubsystem(n)); }
const Station& Subsystem::getStation(const String& n) const
  { return Station::downcast(getSubsystem(n)); }
const Direction& Subsystem::getDirection(const String& n) const
  { return Direction::downcast(getSubsystem(n)); }
const Orientation& Subsystem::getOrientation(const String& n) const
  { return Orientation::downcast(getSubsystem(n)); }
const Frame& Subsystem::getFrame(const String& n) const
  { return Frame::downcast(getSubsystem(n)); }

// updXXX() methods
Subsystem& Subsystem::updSubsystem(const String& n)
  { return updRep().updSubsystem(n); }
Feature& Subsystem::updFeature(const String& n)
{ return Feature::downcast(updSubsystem(n)); }
RealParameter& Subsystem::updRealParameter(const String& n)
  { return RealParameter::downcast(updSubsystem(n)); }
Vec3Parameter& Subsystem::updVec3Parameter(const String& n)
  { return Vec3Parameter::downcast(updSubsystem(n)); }
StationParameter& Subsystem::updStationParameter(const String& n)
  { return StationParameter::downcast(updSubsystem(n)); }
RealMeasure& Subsystem::updRealMeasure(const String& n)
  { return RealMeasure::downcast(updSubsystem(n)); }
Vec3Measure& Subsystem::updVec3Measure(const String& n)
  { return Vec3Measure::downcast(updSubsystem(n)); }
StationMeasure& Subsystem::updStationMeasure(const String& n)
  { return StationMeasure::downcast(updSubsystem(n)); }
Station& Subsystem::updStation(const String& n)
  { return Station::downcast(updSubsystem(n)); }
Direction& Subsystem::updDirection(const String& n)
  { return Direction::downcast(updSubsystem(n)); }
Orientation& Subsystem::updOrientation(const String& n)
  { return Orientation::downcast(updSubsystem(n)); }
Frame& Subsystem::updFrame(const String& n)
  { return Frame::downcast(updSubsystem(n)); }

// addXXX() methods
RealParameter& Subsystem::addRealParameter(const String& n, const Placement& p) {
    RealParameter& rp = RealParameter::downcast(updRep().addSubsystemLike(RealParameter(n), n));
    if (p.hasRep()) rp.place(p);
    return rp;
}
RealMeasure& Subsystem::addRealMeasure(const String& n, const Placement& p) {
    RealMeasure& rm = RealMeasure::downcast(updRep().addSubsystemLike(RealMeasure(n), n));
    if (p.hasRep()) rm.place(p);
    return rm;
}
Vec3Parameter& Subsystem::addVec3Parameter(const String& n, const Placement& p) {
    Vec3Parameter& vp = Vec3Parameter::downcast(updRep().addSubsystemLike(Vec3Parameter(n), n));
    if (p.hasRep()) vp.place(p);
    return vp;
}
Vec3Measure& Subsystem::addVec3Measure(const String& n, const Placement& p) {
    Vec3Measure& vm = Vec3Measure::downcast(updRep().addSubsystemLike(Vec3Measure(n), n));
    if (p.hasRep()) vm.place(p);
    return vm;
}
StationParameter& Subsystem::addStationParameter(const String& n, const Placement& p) {
    StationParameter& sp = StationParameter::downcast(updRep().addSubsystemLike(StationParameter(n), n));
    if (p.hasRep()) sp.place(p);
    return sp;
}
StationMeasure& Subsystem::addStationMeasure(const String& n, const Placement& p) {
    StationMeasure& sm = StationMeasure::downcast(updRep().addSubsystemLike(StationMeasure(n), n));
    if (p.hasRep()) sm.place(p);
    return sm;
}
Station& Subsystem::addStation(const String& n, const Placement& p) {
    Station& s = Station::downcast(updRep().addSubsystemLike(Station(n), n));
    if (p.hasRep()) s.place(p);
    return s;
}
Direction& Subsystem::addDirection(const String& n, const Placement& p) {
    Direction& d = Direction::downcast(updRep().addSubsystemLike(Direction(n), n));
    if (p.hasRep()) d.place(p);
    return d;
}
Orientation& Subsystem::addOrientation(const String& n, const Placement& p) {
    Orientation& o = Orientation::downcast(updRep().addSubsystemLike(Orientation(n), n));
    if (p.hasRep()) o.place(p);
    return o;
}
Frame& Subsystem::addFrame(const String& n, const Placement& p) {
    Frame& f = Frame::downcast(updRep().addSubsystemLike(Frame(n), n));
    if (p.hasRep()) f.place(p);
    return f;
}

Subsystem& Subsystem::addSubsystemLike(const Subsystem& f, const String& n) {
    return updRep().addSubsystemLike(f,n);
}

Feature& Subsystem::addFeatureLike(const Subsystem& f, const String& n, const Placement& p) {
    Feature& fnew = updRep().addFeatureLike(f,n);
    if (p.hasRep()) fnew.place(p);
    return fnew;
}

void Subsystem::checkSubsystemConsistency(const Subsystem* expParent,
                                          int              expIndexInParent,
                                          const Subsystem& expRoot) const {
    if (!rep)
        std::cout << "checkSubsystemConsistency(): NO REP!!!" << std::endl;
    else
        getRep().checkSubsystemConsistency(expParent,expIndexInParent,expRoot);
}

    // FEATURE //
Feature::Feature(const Feature& src) : Subsystem(src) { }
Feature& Feature::operator=(const Feature& src)
  { Subsystem::operator=(src); return *this; }
Feature::~Feature() { }

/*static*/ bool             
Feature::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return FeatureRep::isA(f.getRep());
}
/*static*/ const Feature& 
Feature::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Feature&>(f);
}

/*static*/ Feature&       
Feature::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Feature&>(f);
}

bool Feature::dependsOn(const Feature& f) const {
    if (isSameSubsystem(f)) return true;
    return hasRep() && getRep().dependsOn(f);
}

String Feature::getFeatureTypeName() const {
    return hasRep() ? String(getRep().getFeatureTypeName()) 
               : subsystemHasNoRep(*this);
}

const PlacementValue& Feature::getValue(/*State*/) const {
    try {
        return getRep().getPlacement().getValue(/*State*/);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW4(Exception::FeatureAPIMethodFailed, getFullName(), 
            "getValue", "", exc.getMessageText());
    }
}

bool Feature::hasPlacement() const {
    return getRep().hasPlacement();
}
void Feature::place(const Placement& p) {
    try {
        updRep().place(p);
    }
    catch (const Exception::Base& exc) {
        SIMTK_THROW2(Exception::APIMethodFailed, "Feature::place()", exc.getMessage());
    }
}

const Placement& Feature::getPlacement() const {
    assert(hasPlacement());
    return getRep().getPlacement();
}
    // REAL PARAMETER //

RealParameter::RealParameter(const String& nm) : RealMeasure() {
    rep = new RealParameterRep(*this, std::string(nm)); 
    rep->initializeStandardSubfeatures();
}
RealParameter::RealParameter(const RealParameter& src) : RealMeasure(src) { }
RealParameter& RealParameter::operator=(const RealParameter& src)
  { RealMeasure::operator=(src); return *this; }
RealParameter::~RealParameter() { }

const RealPlacement& 
RealParameter::getPlacement() const {
    return RealPlacement::downcast(getRep().getPlacement());
}

const Real& 
RealParameter::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
RealParameter::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return RealParameterRep::isA(f.getRep());
}
/*static*/ const RealParameter& 
RealParameter::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealParameter&>(f);
}

/*static*/ RealParameter&       
RealParameter::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealParameter&>(f);
}

    // VEC3 PARAMETER //

Vec3Parameter::Vec3Parameter(const String& nm) : Vec3Measure() {
    rep = new Vec3ParameterRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Vec3Parameter::Vec3Parameter(const Vec3Parameter& src) : Vec3Measure(src) { }
Vec3Parameter& Vec3Parameter::operator=(const Vec3Parameter& src)
  { Vec3Measure::operator=(src); return *this; }
Vec3Parameter::~Vec3Parameter() { }

const Vec3Placement& 
Vec3Parameter::getPlacement() const {
    return Vec3Placement::downcast(getRep().getPlacement());
}

const Vec3& 
Vec3Parameter::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
Vec3Parameter::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return Vec3ParameterRep::isA(f.getRep());
}
/*static*/ const Vec3Parameter& 
Vec3Parameter::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Vec3Parameter&>(f);
}

/*static*/ Vec3Parameter&       
Vec3Parameter::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Vec3Parameter&>(f);
}

    // STATION PARAMETER //
StationParameter::StationParameter(const String& nm) : StationMeasure() {
    rep = new StationParameterRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
StationParameter::StationParameter(const StationParameter& src) : StationMeasure(src) { }
StationParameter& StationParameter::operator=(const StationParameter& src)
  { StationMeasure::operator=(src); return *this; }
StationParameter::~StationParameter() { }

const StationPlacement& 
StationParameter::getPlacement() const {
    return StationPlacement::downcast(getRep().getPlacement());
}

const Vec3& 
StationParameter::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
StationParameter::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return StationParameterRep::isA(f.getRep());
}
/*static*/ const StationParameter& 
StationParameter::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationParameter&>(f);
}

/*static*/ StationParameter&       
StationParameter::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationParameter&>(f);
}

    // REAL MEASURE //
RealMeasure::RealMeasure(const String& nm) : Feature() {
    rep = new RealMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
RealMeasure::RealMeasure(const RealMeasure& src) : Feature(src) { }
RealMeasure& RealMeasure::operator=(const RealMeasure& src)
  { Feature::operator=(src); return *this; }
RealMeasure::~RealMeasure() { }

const RealPlacement& 
RealMeasure::getPlacement() const {
    return RealPlacement::downcast(getRep().getPlacement());
}

const Real& 
RealMeasure::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
RealMeasure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return RealMeasureRep::isA(f.getRep());
}
/*static*/ const RealMeasure& 
RealMeasure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RealMeasure&>(f);
}

/*static*/ RealMeasure&       
RealMeasure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RealMeasure&>(f);
}

    // VEC3 MEASURE //
Vec3Measure::Vec3Measure(const String& nm) : Feature() {
    rep = new Vec3MeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Vec3Measure::Vec3Measure(const Vec3Measure& src) : Feature(src) { }
Vec3Measure& Vec3Measure::operator=(const Vec3Measure& src)
  { Feature::operator=(src); return *this; }
Vec3Measure::~Vec3Measure() { }

const Vec3Placement& 
Vec3Measure::getPlacement() const {
    return Vec3Placement::downcast(getRep().getPlacement());
}

const Vec3& 
Vec3Measure::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
Vec3Measure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return Vec3MeasureRep::isA(f.getRep());
}
/*static*/ const Vec3Measure& 
Vec3Measure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Vec3Measure&>(f);
}

/*static*/ Vec3Measure&       
Vec3Measure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Vec3Measure&>(f);
}

    // STATION MEASURE //
StationMeasure::StationMeasure(const String& nm) : Feature() {
    rep = new StationMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
StationMeasure::StationMeasure(const StationMeasure& src) : Feature(src) { }
StationMeasure& StationMeasure::operator=(const StationMeasure& src)
  { Feature::operator=(src); return *this; }
StationMeasure::~StationMeasure() { }

const StationPlacement& 
StationMeasure::getPlacement() const {
    return StationPlacement::downcast(getRep().getPlacement());
}

const Vec3& 
StationMeasure::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
StationMeasure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return StationMeasureRep::isA(f.getRep());
}
/*static*/ const StationMeasure& 
StationMeasure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const StationMeasure&>(f);
}

/*static*/ StationMeasure&       
StationMeasure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<StationMeasure&>(f);
}

    // STATION //
Station::Station(const String& nm) : Feature() {
    rep = new StationRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Station::Station(const Station& src) : Feature(src) { }
Station& Station::operator=(const Station& src)
  { Feature::operator=(src); return *this; }
Station::~Station() { }

const StationPlacement& 
Station::getPlacement() const {
    return StationPlacement::downcast(getRep().getPlacement());
}

const Vec3& 
Station::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
Station::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return StationRep::isA(f.getRep());
}
/*static*/ const Station& 
Station::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Station&>(f);
}

/*static*/ Station&       
Station::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Station&>(f);
}

    // DIRECTION MEASURE //
DirectionMeasure::DirectionMeasure(const String& nm) : Feature() {
    rep = new DirectionMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
DirectionMeasure::DirectionMeasure(const DirectionMeasure& src) : Feature(src) { }
DirectionMeasure& DirectionMeasure::operator=(const DirectionMeasure& src)
  { Feature::operator=(src); return *this; }
DirectionMeasure::~DirectionMeasure() { }

const DirectionPlacement& 
DirectionMeasure::getPlacement() const {
    return DirectionPlacement::downcast(getRep().getPlacement());
}

const Vec3& 
DirectionMeasure::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
DirectionMeasure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return DirectionMeasureRep::isA(f.getRep());
}
/*static*/ const DirectionMeasure& 
DirectionMeasure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const DirectionMeasure&>(f);
}

/*static*/ DirectionMeasure&       
DirectionMeasure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<DirectionMeasure&>(f);
}

    // DIRECTION //
Direction::Direction(const String& nm) : Feature() {
    rep = new DirectionRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Direction::Direction(const Direction& src) : Feature(src) { }
Direction& Direction::operator=(const Direction& src)
  { Feature::operator=(src); return *this; }
Direction::~Direction() { }

const DirectionPlacement& 
Direction::getPlacement() const {
    return DirectionPlacement::downcast(getRep().getPlacement());
}

const Vec3& 
Direction::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
Direction::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return DirectionRep::isA(f.getRep());
}
/*static*/ const Direction& 
Direction::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Direction&>(f);
}

/*static*/ Direction&       
Direction::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Direction&>(f);
}

    // ORIENTATION MEASURE //
OrientationMeasure::OrientationMeasure(const String& nm) : Feature() {
    rep = new OrientationMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
OrientationMeasure::OrientationMeasure(const OrientationMeasure& src) : Feature(src) { }
OrientationMeasure& OrientationMeasure::operator=(const OrientationMeasure& src)
  { Feature::operator=(src); return *this; }
OrientationMeasure::~OrientationMeasure() { }

const OrientationPlacement& 
OrientationMeasure::getPlacement() const {
    return OrientationPlacement::downcast(getRep().getPlacement());
}

const Mat33& 
OrientationMeasure::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

/*static*/ bool             
OrientationMeasure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return OrientationMeasureRep::isA(f.getRep());
}
/*static*/ const OrientationMeasure& 
OrientationMeasure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const OrientationMeasure&>(f);
}

/*static*/ OrientationMeasure&       
OrientationMeasure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<OrientationMeasure&>(f);
}

    // ORIENTATION //
Orientation::Orientation(const String& nm) : Feature() {
    rep = new OrientationRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Orientation::Orientation(const Orientation& src) : Feature(src) { }
Orientation& Orientation::operator=(const Orientation& src)
  { Feature::operator=(src); return *this; }
Orientation::~Orientation() { }

const OrientationPlacement& 
Orientation::getPlacement() const {
    return OrientationPlacement::downcast(getRep().getPlacement());
}

const Mat33& 
Orientation::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

const Direction& 
Orientation::getAxis(int i) const {
    return OrientationRep::downcast(getRep()).getAxis(i);
}

/*static*/ bool             
Orientation::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return OrientationRep::isA(f.getRep());
}
/*static*/ const Orientation& 
Orientation::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Orientation&>(f);
}

/*static*/ Orientation&       
Orientation::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Orientation&>(f);
}

    // FRAME //
Frame::Frame(const String& nm) : Feature() {
    rep = new FrameRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Frame::Frame(const Frame& src) : Feature(src) { }
Frame& Frame::operator=(const Frame& src)
  { Feature::operator=(src); return *this; }
Frame::~Frame() { }

const FramePlacement& 
Frame::getPlacement() const {
    return FramePlacement::downcast(getRep().getPlacement());
}

const Mat34& 
Frame::getValue(/*State*/) const {
    return getPlacement().getRep().getValue(/*State*/);
}

const Orientation& Frame::getOrientation() const {
    return FrameRep::downcast(getRep()).getOrientation();
}
const Station& Frame::getOrigin() const {
    return FrameRep::downcast(getRep()).getOrigin();
}

/*static*/ bool             
Frame::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return FrameRep::isA(f.getRep());
}
/*static*/ const Frame& 
Frame::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Frame&>(f);
}

/*static*/ Frame&       
Frame::downcast(Subsystem& f) {
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
