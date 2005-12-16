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
 * Implementations of API-level methods behind the client-side Subsystem.
 */

#include "SimbodyCommon.h"
#include "Feature.h"
#include "FeatureRep.h"

#include <string>
#include <iostream> 
#include <sstream>


namespace simtk {


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

static String subsystemHasNoRep(const Subsystem& s) {
    std::ostringstream o;
    o << "<SUBSYSTEM AT 0x" << &s << " WITH NULL REP>";
    return String(o.str());
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
    return Feature::downcast(*this).getRep().getPlacementSlot().getValueSlot().getValue();
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
        s << (fr.hasPlacement() ? fr.getPlacementSlot().toString(linePrefix)
                                : String("NO PLACEMENT"));
    }

    const size_t nSubsystems      = sr.getNSubsystems();
    const size_t nPlacement       = sr.getNPlacements();
    const size_t nPlacementValues = sr.getNPlacementValues();
    const std::string nextIndent  = linePrefix + "    ";

    if (nSubsystems) {
        s << std::endl << linePrefix << "  Subsystems (" << nSubsystems << "):";
        for (size_t i=0; i < nSubsystems; ++i)
            s  << std::endl << nextIndent << sr.getSubsystem(i).toString(nextIndent);
    }
    if (nPlacement) {
        s << std::endl << linePrefix << "  Placement Slots (" << nPlacement << "):";
        for (size_t i=0; i < nPlacement; ++i)
            s  << std::endl << nextIndent << sr.getPlacementSlot(i).toString(nextIndent);
    }
    if (nPlacementValues) {
        s << std::endl << linePrefix << "  PlacementValue Slots (" << nPlacementValues << "):";
        for (size_t i=0; i < nPlacementValues; ++i)
            s  << std::endl << nextIndent << sr.getPlacementValueSlot(i).toString(nextIndent);
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

} // namespace simtk
