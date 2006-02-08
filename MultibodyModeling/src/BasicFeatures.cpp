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
 * Implementations of API-level BasicFeature methods.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/internal/Subsystem.h"
#include "simbody/internal/Feature.h"
#include "simbody/internal/BasicFeatures.h"

#include "SubsystemRep.h"
#include "FeatureRep.h"
#include "BasicFeaturesRep.h"

#include <string>
#include <iostream> 
#include <sstream>


namespace simtk {

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
RealParameter::getValue() const {
    return PlacementValue_<Real>::downcast(getRep().getPlacementSlot().getValue());
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
Vec3Parameter::getValue() const {
    return PlacementValue_<Vec3>::downcast(getRep().getPlacementSlot().getValue());
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
StationParameter::getValue() const {
    return PlacementValue_<Vec3>::downcast(getRep().getPlacementSlot().getValue());
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
RealMeasure::getValue() const {
    return PlacementValue_<Real>::downcast(getRep().getPlacementSlot().getValue());
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
Vec3Measure::getValue() const {
    return PlacementValue_<Vec3>::downcast(getRep().getPlacementSlot().getValue());
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
StationMeasure::getValue() const {
    return PlacementValue_<Vec3>::downcast(getRep().getPlacementSlot().getValue());
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
Station::getValue() const {
    return PlacementValue_<Vec3>::downcast(getRep().getPlacementSlot().getValue());
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

const UnitVec3& 
DirectionMeasure::getValue() const {
    return PlacementValue_<UnitVec3>::downcast(getRep().getPlacementSlot().getValue());
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

const UnitVec3& 
Direction::getValue() const {
    return PlacementValue_<UnitVec3>::downcast(getRep().getPlacementSlot().getValue());
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

const RotationMat& 
OrientationMeasure::getValue() const {
    return PlacementValue_<RotationMat>::downcast(getRep().getPlacementSlot().getValue());
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

const RotationMat& 
Orientation::getValue() const {
    return PlacementValue_<RotationMat>::downcast(getRep().getPlacementSlot().getValue());
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


    // INERTIA MEASURE //
InertiaMeasure::InertiaMeasure(const String& nm) : Feature() {
    rep = new InertiaMeasureRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
InertiaMeasure::InertiaMeasure(const InertiaMeasure& src) : Feature(src) { }
InertiaMeasure& InertiaMeasure::operator=(const InertiaMeasure& src)
  { Feature::operator=(src); return *this; }
InertiaMeasure::~InertiaMeasure() { }

const InertiaPlacement& 
InertiaMeasure::getPlacement() const {
    return InertiaPlacement::downcast(getRep().getPlacement());
}

const InertiaMat& 
InertiaMeasure::getValue() const {
    return PlacementValue_<InertiaMat>::downcast(getRep().getPlacementSlot().getValue());
}

/*static*/ bool             
InertiaMeasure::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return InertiaMeasureRep::isA(f.getRep());
}
/*static*/ const InertiaMeasure& 
InertiaMeasure::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const InertiaMeasure&>(f);
}

/*static*/ InertiaMeasure&       
InertiaMeasure::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<InertiaMeasure&>(f);
}

    // FRAME //
FrameFeature::FrameFeature(const String& nm) : Feature() {
    rep = new FrameRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
FrameFeature::FrameFeature(const FrameFeature& src) : Feature(src) { }
FrameFeature& FrameFeature::operator=(const FrameFeature& src)
  { Feature::operator=(src); return *this; }
FrameFeature::~FrameFeature() { }

const FramePlacement& 
FrameFeature::getPlacement() const {
    return FramePlacement::downcast(getRep().getPlacement());
}

const Frame& 
FrameFeature::getValue() const {
    return PlacementValue_<Frame>::downcast(getRep().getPlacementSlot().getValue());
}

const Orientation& FrameFeature::getOrientation() const {
    return FrameRep::downcast(getRep()).getOrientation();
}
const Station& FrameFeature::getOrigin() const {
    return FrameRep::downcast(getRep()).getOrigin();
}

/*static*/ bool             
FrameFeature::isInstanceOf(const Subsystem& f) {
    if (!f.hasRep()) return false;
    return FrameRep::isA(f.getRep());
}
/*static*/ const FrameFeature& 
FrameFeature::downcast(const Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const FrameFeature&>(f);
}

/*static*/ FrameFeature&       
FrameFeature::downcast(Subsystem& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<FrameFeature&>(f);
}

} // namespace simtk
