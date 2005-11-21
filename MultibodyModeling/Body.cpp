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
 * Implementations of Body Features for Simbody.
 */

#include "SimbodyCommon.h"
#include "Body.h"
#include "BodyRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::endl;
using std::ostream;

namespace simtk {

    // BODY //
const RealMeasure& Body::getMassMeasure() const {
    assert(rep);
    return BodyRep::downcast(getRep()).getMassMeasure();
}
const StationMeasure& Body::getCentroidMeasure() const {
    assert(rep);
    return BodyRep::downcast(getRep()).getCentroidMeasure();
}


/*static*/ bool             
Body::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return BodyRep::isA(f.getRep());
}
/*static*/ const Body& 
Body::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Body&>(f);
}

/*static*/ Body&       
Body::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Body&>(f);
}

    // RIGID BODY //

RigidBody::RigidBody(const String& nm) : Body() {
    rep = new RigidBodyRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
RigidBody::RigidBody(const RigidBody& src)
  : Body(src) { }
RigidBody& RigidBody::operator=(const RigidBody& src) {
    Body::operator=(src); return *this;
}
RigidBody::~RigidBody() { }

MassElement& RigidBody::addMassElementLike(const MassElement& me, const String& nm,
                                           const Placement& pl) {
    Placement& p = updRep().addPlacementLike(pl);
    MassElement& m = MassElement::downcast(updRep().addSubfeatureLike(me, nm));
    m.place(p);
    return m;
}
MassElement& RigidBody::addMassElementLike(const MassElement& me, const String& nm) {
    return MassElement::downcast(updRep().addSubfeatureLike(me, nm));
}

/*static*/ bool             
RigidBody::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return RigidBodyRep::isA(f.getRep());
}
/*static*/ const RigidBody& 
RigidBody::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const RigidBody&>(f);
}

/*static*/ RigidBody&       
RigidBody::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<RigidBody&>(f);
}

    // DEFORMABLE BODY //

    // MULTIBODY //
Multibody::Multibody(const String& nm) : Body() {
    rep = new MultibodyRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Multibody::Multibody(const Multibody& src)
  : Body(src) { }
Multibody& Multibody::operator=(const Multibody& src) {
    Body::operator=(src); return *this;
}
Multibody::~Multibody() { }

const Frame& Multibody::getGroundFrame() const {
    return getFrame("Ground");
}

RigidBody& Multibody::addGroundBody() {
    RigidBody& subBody = 
        RigidBody::downcast(updRep().addSubfeatureLike(RigidBody("Ground"), "Ground"));
    return subBody;
}

RigidBody& Multibody::addRigidBody(const String& nm) {
    RigidBody& subBody = 
        RigidBody::downcast(updRep().addSubfeatureLike(RigidBody(nm), nm));
    return subBody;
}

RigidBody& Multibody::addRigidBodyLike(const RigidBody& b, const String& nm) {
    RigidBody& subBody = RigidBody::downcast(updRep().addSubfeatureLike(b, nm));
    return subBody;
}

Body& Multibody::addBodyLike(const Body& b, const String& nm) {
    Body& subBody = Body::downcast(updRep().addSubfeatureLike(b, nm));
    return subBody;
}

Joint& Multibody::addJoint(JointType jt, const String& nm) {
    Joint& j = Joint::downcast(updRep().addSubfeatureLike(Joint(jt,nm), nm));
    return j;
}

Joint& Multibody::addJoint(JointType jt, const String& nm,
                           const Placement& reference,
                           const Placement& moving) {
    Joint& j = Joint::downcast(updRep().addSubfeatureLike(Joint(jt,nm), nm));
    j.updFrame("reference").place(reference);
    j.updFrame("moving").place(moving);
    return j;
}

/*static*/ bool             
Multibody::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return MultibodyRep::isA(f.getRep());
}
/*static*/ const Multibody& 
Multibody::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Multibody&>(f);
}

/*static*/ Multibody&       
Multibody::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Multibody&>(f);
}

    // JOINT //

Joint::Joint(JointType jt, const String& nm) : Feature() {
    rep = new JointRep(*this, jt, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Joint::Joint(const Joint& src)
  : Feature(src) { }
Joint& Joint::operator=(const Joint& src) {
    Feature::operator=(src); return *this;
}
Joint::~Joint() { }

/*static*/ bool             
Joint::isInstanceOf(const Feature& f) {
    if (!f.hasRep()) return false;
    return JointRep::isA(f.getRep());
}
/*static*/ const Joint& 
Joint::downcast(const Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<const Joint&>(f);
}

/*static*/ Joint&       
Joint::downcast(Feature& f) {
    assert(isInstanceOf(f));
    return reinterpret_cast<Joint&>(f);
}

} // namespace simtk
