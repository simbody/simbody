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
const RealMeasure& Body::getMass() const {
    assert(rep);
    return BodyRep::downcast(getRep()).getMass();
}
const StationMeasure& Body::getCentroid() const {
    assert(rep);
    return BodyRep::downcast(getRep()).getCentroid();
}
//const InertiaMeasure& Body::getCentroidalInertia() const {
//    assert(rep);
//    return BodyRep::downcast(getRep()).getCentroidalInertia();
//}

/*static*/ bool             
Body::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return BodyRep::isA(s.getRep());
}
/*static*/ const Body& 
Body::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Body&>(s);
}

/*static*/ Body&       
Body::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Body&>(s);
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
    MassElement& m = MassElement::downcast(updRep().addFeatureLike(me, nm));
    m.place(pl);
    return m;
}
MassElement& RigidBody::addMassElementLike(const MassElement& me, const String& nm) {
    return MassElement::downcast(updRep().addSubsystemLike(me, nm));
}

/*static*/ bool             
RigidBody::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return RigidBodyRep::isA(s.getRep());
}
/*static*/ const RigidBody& 
RigidBody::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const RigidBody&>(s);
}

/*static*/ RigidBody&       
RigidBody::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<RigidBody&>(s);
}

    // DEFORMABLE BODY //

    // MULTIBODY //
Multibody::Multibody(const String& nm) : Subsystem() {
    rep = new MultibodyRep(*this, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Multibody::Multibody(const Multibody& src)
  : Subsystem(src) { }
Multibody& Multibody::operator=(const Multibody& src) {
    Subsystem::operator=(src); return *this;
}
Multibody::~Multibody() { }

const FrameFeature& Multibody::getGroundFrame() const {
    return getFrame("Ground");
}

RigidBody& Multibody::addGroundBody() {
    RigidBody& subBody = 
        RigidBody::downcast(updRep().addSubsystemLike(RigidBody("Ground"), "Ground"));
    return subBody;
}

RigidBody& Multibody::addRigidBody(const String& nm) {
    RigidBody& subBody = 
        RigidBody::downcast(updRep().addSubsystemLike(RigidBody(nm), nm));
    return subBody;
}

RigidBody& Multibody::addRigidBodyLike(const RigidBody& b, const String& nm) {
    RigidBody& subBody = RigidBody::downcast(updRep().addSubsystemLike(b, nm));
    return subBody;
}

Body& Multibody::addBodyLike(const Body& b, const String& nm) {
    Body& subBody = Body::downcast(updRep().addSubsystemLike(b, nm));
    return subBody;
}

Joint& Multibody::addJoint(JointType jt, const String& nm) {
    Joint& j = Joint::downcast(updRep().addSubsystemLike(Joint(jt,nm), nm));
    return j;
}

Joint& Multibody::addJoint(JointType jt, const String& nm,
                           const Placement& reference,
                           const Placement& moving) {
    Joint& j = Joint::downcast(updRep().addSubsystemLike(Joint(jt,nm), nm));
    j.updFrame("reference").place(reference);
    j.updFrame("moving").place(moving);
    return j;
}

/*static*/ bool             
Multibody::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return MultibodyRep::isA(s.getRep());
}
/*static*/ const Multibody& 
Multibody::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Multibody&>(s);
}

/*static*/ Multibody&       
Multibody::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Multibody&>(s);
}

    // JOINT //

Joint::Joint(JointType jt, const String& nm) : Subsystem() {
    rep = new JointRep(*this, jt, std::string(nm));
    rep->initializeStandardSubfeatures();
}
Joint::Joint(const Joint& src)
  : Subsystem(src) { }
Joint& Joint::operator=(const Joint& src) {
    Subsystem::operator=(src); return *this;
}
Joint::~Joint() { }

/*static*/ bool             
Joint::isInstanceOf(const Subsystem& s) {
    if (!s.hasRep()) return false;
    return JointRep::isA(s.getRep());
}
/*static*/ const Joint& 
Joint::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Joint&>(s);
}

/*static*/ Joint&       
Joint::downcast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Joint&>(s);
}

} // namespace simtk
