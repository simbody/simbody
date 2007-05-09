/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 *
 * Private implementation of Mobilizer, and its included subclasses which
 * represent the built-in mobilizer types.
 */

#include "SimTKsimbody.h"
#include "simbody/internal/Mobilizer.h"

#include "MobilizerRep.h"

namespace SimTK {

    ///////////////
    // MOBILIZER //
    ///////////////

bool Mobilizer::isEmptyHandle() const {return rep==0;}
bool Mobilizer::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}

Mobilizer::~Mobilizer() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

Mobilizer::Mobilizer(const Mobilizer& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();	// create a new object
        rep->setMyHandle(*this);
    }
}

Mobilizer& Mobilizer::operator=(const Mobilizer& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
			rep = src.rep->clone();	// create a new object
			rep->setMyHandle(*this);
        }
    }
    return *this;
}

    ////////////////////
    // MOBILIZER::PIN //
    ////////////////////

Mobilizer::Pin::Pin() {
    rep = new PinRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Pin::isInstanceOf(const Mobilizer& s) {
    return PinRep::isA(s.getRep());
}
const Mobilizer::Pin& Mobilizer::Pin::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Pin&>(s);
}
Mobilizer::Pin& Mobilizer::Pin::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Pin&>(s);
}
const Mobilizer::Pin::PinRep& Mobilizer::Pin::getRep() const {
    return dynamic_cast<const PinRep&>(*rep);
}
Mobilizer::Pin::PinRep& Mobilizer::Pin::updRep() {
    return dynamic_cast<PinRep&>(*rep);
}

    ///////////////////////
    // MOBILIZER::SLIDER //
    ///////////////////////

Mobilizer::Slider::Slider() {
    rep = new SliderRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Slider::isInstanceOf(const Mobilizer& s) {
    return SliderRep::isA(s.getRep());
}
const Mobilizer::Slider& Mobilizer::Slider::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Slider&>(s);
}
Mobilizer::Slider& Mobilizer::Slider::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Slider&>(s);
}
const Mobilizer::Slider::SliderRep& Mobilizer::Slider::getRep() const {
    return dynamic_cast<const SliderRep&>(*rep);
}
Mobilizer::Slider::SliderRep& Mobilizer::Slider::updRep() {
    return dynamic_cast<SliderRep&>(*rep);
}

    //////////////////////////
    // MOBILIZER::UNIVERSAL //
    //////////////////////////

Mobilizer::Universal::Universal() {
    rep = new UniversalRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Universal::isInstanceOf(const Mobilizer& s) {
    return UniversalRep::isA(s.getRep());
}
const Mobilizer::Universal& Mobilizer::Universal::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Universal&>(s);
}
Mobilizer::Universal& Mobilizer::Universal::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Universal&>(s);
}
const Mobilizer::Universal::UniversalRep& Mobilizer::Universal::getRep() const {
    return dynamic_cast<const UniversalRep&>(*rep);
}
Mobilizer::Universal::UniversalRep& Mobilizer::Universal::updRep() {
    return dynamic_cast<UniversalRep&>(*rep);
}

    /////////////////////////
    // MOBILIZER::CYLINDER //
    /////////////////////////

Mobilizer::Cylinder::Cylinder() {
    rep = new CylinderRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Cylinder::isInstanceOf(const Mobilizer& s) {
    return CylinderRep::isA(s.getRep());
}
const Mobilizer::Cylinder& Mobilizer::Cylinder::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Cylinder&>(s);
}
Mobilizer::Cylinder& Mobilizer::Cylinder::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Cylinder&>(s);
}
const Mobilizer::Cylinder::CylinderRep& Mobilizer::Cylinder::getRep() const {
    return dynamic_cast<const CylinderRep&>(*rep);
}
Mobilizer::Cylinder::CylinderRep& Mobilizer::Cylinder::updRep() {
    return dynamic_cast<CylinderRep&>(*rep);
}

    /////////////////////////////
    // MOBILIZER::BEND STRETCH //
    /////////////////////////////

Mobilizer::BendStretch::BendStretch() {
    rep = new BendStretchRep(); rep->setMyHandle(*this);
}
bool Mobilizer::BendStretch::isInstanceOf(const Mobilizer& s) {
    return BendStretchRep::isA(s.getRep());
}
const Mobilizer::BendStretch& Mobilizer::BendStretch::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const BendStretch&>(s);
}
Mobilizer::BendStretch& Mobilizer::BendStretch::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<BendStretch&>(s);
}
const Mobilizer::BendStretch::BendStretchRep& Mobilizer::BendStretch::getRep() const {
    return dynamic_cast<const BendStretchRep&>(*rep);
}
Mobilizer::BendStretch::BendStretchRep& Mobilizer::BendStretch::updRep() {
    return dynamic_cast<BendStretchRep&>(*rep);
}

    ///////////////////////
    // MOBILIZER::PLANAR //
    ///////////////////////

Mobilizer::Planar::Planar() {
    rep = new PlanarRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Planar::isInstanceOf(const Mobilizer& s) {
    return PlanarRep::isA(s.getRep());
}
const Mobilizer::Planar& Mobilizer::Planar::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Planar&>(s);
}
Mobilizer::Planar& Mobilizer::Planar::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Planar&>(s);
}
const Mobilizer::Planar::PlanarRep& Mobilizer::Planar::getRep() const {
    return dynamic_cast<const PlanarRep&>(*rep);
}
Mobilizer::Planar::PlanarRep& Mobilizer::Planar::updRep() {
    return dynamic_cast<PlanarRep&>(*rep);
}

    ///////////////////////
    // MOBILIZER::GIMBAL //
    ///////////////////////

Mobilizer::Gimbal::Gimbal() {
    rep = new GimbalRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Gimbal::isInstanceOf(const Mobilizer& s) {
    return GimbalRep::isA(s.getRep());
}
const Mobilizer::Gimbal& Mobilizer::Gimbal::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Gimbal&>(s);
}
Mobilizer::Gimbal& Mobilizer::Gimbal::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Gimbal&>(s);
}
const Mobilizer::Gimbal::GimbalRep& Mobilizer::Gimbal::getRep() const {
    return dynamic_cast<const GimbalRep&>(*rep);
}
Mobilizer::Gimbal::GimbalRep& Mobilizer::Gimbal::updRep() {
    return dynamic_cast<GimbalRep&>(*rep);
}

    ///////////////////////////////////
    // MOBILIZER::ORIENTATION (BALL) //
    ///////////////////////////////////

Mobilizer::Orientation::Orientation() {
    rep = new OrientationRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Orientation::isInstanceOf(const Mobilizer& s) {
    return OrientationRep::isA(s.getRep());
}
const Mobilizer::Orientation& Mobilizer::Orientation::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Orientation&>(s);
}
Mobilizer::Orientation& Mobilizer::Orientation::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Orientation&>(s);
}
const Mobilizer::Orientation::OrientationRep& Mobilizer::Orientation::getRep() const {
    return dynamic_cast<const OrientationRep&>(*rep);
}
Mobilizer::Orientation::OrientationRep& Mobilizer::Orientation::updRep() {
    return dynamic_cast<OrientationRep&>(*rep);
}

    ////////////////////////////
    // MOBILIZER::TRANSLATION //
    ////////////////////////////

Mobilizer::Translation::Translation() {
    rep = new TranslationRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Translation::isInstanceOf(const Mobilizer& s) {
    return TranslationRep::isA(s.getRep());
}
const Mobilizer::Translation& Mobilizer::Translation::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Translation&>(s);
}
Mobilizer::Translation& Mobilizer::Translation::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Translation&>(s);
}
const Mobilizer::Translation::TranslationRep& Mobilizer::Translation::getRep() const {
    return dynamic_cast<const TranslationRep&>(*rep);
}
Mobilizer::Translation::TranslationRep& Mobilizer::Translation::updRep() {
    return dynamic_cast<TranslationRep&>(*rep);
}

    /////////////////////
    // MOBILIZER::FREE //
    /////////////////////

Mobilizer::Free::Free() {
    rep = new FreeRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Free::isInstanceOf(const Mobilizer& s) {
    return FreeRep::isA(s.getRep());
}
const Mobilizer::Free& Mobilizer::Free::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Free&>(s);
}
Mobilizer::Free& Mobilizer::Free::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Free&>(s);
}
const Mobilizer::Free::FreeRep& Mobilizer::Free::getRep() const {
    return dynamic_cast<const FreeRep&>(*rep);
}
Mobilizer::Free::FreeRep& Mobilizer::Free::updRep() {
    return dynamic_cast<FreeRep&>(*rep);
}

    /////////////////////////////////
    // MOBILIZER::LINE ORIENTATION //
    /////////////////////////////////

Mobilizer::LineOrientation::LineOrientation() {
    rep = new LineOrientationRep(); rep->setMyHandle(*this);
}
bool Mobilizer::LineOrientation::isInstanceOf(const Mobilizer& s) {
    return LineOrientationRep::isA(s.getRep());
}
const Mobilizer::LineOrientation& Mobilizer::LineOrientation::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const LineOrientation&>(s);
}
Mobilizer::LineOrientation& Mobilizer::LineOrientation::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<LineOrientation&>(s);
}
const Mobilizer::LineOrientation::LineOrientationRep& Mobilizer::LineOrientation::getRep() const {
    return dynamic_cast<const LineOrientationRep&>(*rep);
}
Mobilizer::LineOrientation::LineOrientationRep& Mobilizer::LineOrientation::updRep() {
    return dynamic_cast<LineOrientationRep&>(*rep);
}

    //////////////////////////
    // MOBILIZER::FREE LINE //
    //////////////////////////

Mobilizer::FreeLine::FreeLine() {
    rep = new FreeLineRep(); rep->setMyHandle(*this);
}
bool Mobilizer::FreeLine::isInstanceOf(const Mobilizer& s) {
    return FreeLineRep::isA(s.getRep());
}
const Mobilizer::FreeLine& Mobilizer::FreeLine::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const FreeLine&>(s);
}
Mobilizer::FreeLine& Mobilizer::FreeLine::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<FreeLine&>(s);
}
const Mobilizer::FreeLine::FreeLineRep& Mobilizer::FreeLine::getRep() const {
    return dynamic_cast<const FreeLineRep&>(*rep);
}
Mobilizer::FreeLine::FreeLineRep& Mobilizer::FreeLine::updRep() {
    return dynamic_cast<FreeLineRep&>(*rep);
}

    /////////////////////
    // MOBILIZER::WELD //
    /////////////////////

Mobilizer::Weld::Weld() {
    rep = new WeldRep(); rep->setMyHandle(*this);
}
bool Mobilizer::Weld::isInstanceOf(const Mobilizer& s) {
    return WeldRep::isA(s.getRep());
}
const Mobilizer::Weld& Mobilizer::Weld::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Weld&>(s);
}
Mobilizer::Weld& Mobilizer::Weld::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Weld&>(s);
}
const Mobilizer::Weld::WeldRep& Mobilizer::Weld::getRep() const {
    return dynamic_cast<const WeldRep&>(*rep);
}
Mobilizer::Weld::WeldRep& Mobilizer::Weld::updRep() {
    return dynamic_cast<WeldRep&>(*rep);
}

    //////////////////////
    // MOBILIZER::SCREW //
    //////////////////////

Mobilizer::Screw::Screw(Real pitch) {
    rep = new ScrewRep(pitch); rep->setMyHandle(*this);
}
bool Mobilizer::Screw::isInstanceOf(const Mobilizer& s) {
    return ScrewRep::isA(s.getRep());
}
const Mobilizer::Screw& Mobilizer::Screw::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Screw&>(s);
}
Mobilizer::Screw& Mobilizer::Screw::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Screw&>(s);
}
const Mobilizer::Screw::ScrewRep& Mobilizer::Screw::getRep() const {
    return dynamic_cast<const ScrewRep&>(*rep);
}
Mobilizer::Screw::ScrewRep& Mobilizer::Screw::updRep() {
    return dynamic_cast<ScrewRep&>(*rep);
}

    /////////////////////
    // MOBILIZER::USER //
    /////////////////////

Mobilizer::User::User(int nMobilities, int nCoordinates) {
    rep = new UserRep(nMobilities, nCoordinates); rep->setMyHandle(*this);
}
bool Mobilizer::User::isInstanceOf(const Mobilizer& s) {
    return UserRep::isA(s.getRep());
}
const Mobilizer::User& Mobilizer::User::downcast(const Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const User&>(s);
}
Mobilizer::User& Mobilizer::User::updDowncast(Mobilizer& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<User&>(s);
}
const Mobilizer::User::UserRep& Mobilizer::User::getRep() const {
    return dynamic_cast<const UserRep&>(*rep);
}
Mobilizer::User::UserRep& Mobilizer::User::updRep() {
    return dynamic_cast<UserRep&>(*rep);
}


} // namespace SimTK

