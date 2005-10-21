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

#include "MultibodyModeling.h"
#include "MultibodyModelingRep.h"

#include <string>
#include <iostream> 
#include <sstream>
using std::endl;
using std::ostream;

namespace simtk {

    // PLACEMENT //

Placement::Placement(const Placement& p) 
  : rep(p.rep?p.rep->clone():0) { 
}
Placement& Placement::operator=(const Placement& p) {
    if (this!=&p) {
        delete rep;
        rep = p.rep ? p.rep->clone() : 0;
    }
    return *this;
}
Placement::~Placement() { 
    delete rep;
}
bool Placement::hasPlacement() const {
    return rep && rep->hasPlacement();
}
const Frame&
Placement::getParentFrame() const {
    assert(rep); return rep->getParentFrame();
}

ostream& operator<<(ostream& o, const Placement& p) {
    const PlacementRep* r = PlacementRep::getRep(p);
    if (!r || !r->hasPlacement()) {
        o << "Placement at 0x" << &p << " has NULL ";
        o << (r ? "placement." : "rep.");
    } else {
        o << r->toString();
    }
    return o;
}

    // FRAME //

Frame::Frame(const String& nm)
  { rep = new FrameRep(*this, nm.c_str()); }
Frame::~Frame() { delete rep; }
String Frame::getFullName() const {
    if (rep) return String(rep->getFullName());
    std::ostringstream s;
    s << "FRAME AT 0x" << this << " HAS NO REP";
    return String(s.str());
}


ostream& operator<<(ostream& o, const Frame& f) {
    const FrameRep* r = FrameRep::getRep(f);
    o << "Frame ";
    if (!r)
        return o << "at 0x" << &f << " has NULL rep." << endl;

    o << r->getFullName() << ": " << r->getPlacement() << endl;
    o << "  Features:" << endl;
    for (size_t i=0; i<r->getStations().size(); ++i)
        o << "    " << r->getStations()[i];
    for (size_t i=0; i<r->getDirections().size(); ++i)
        o << "    " << r->getDirections()[i];
    for (size_t i=0; i<r->getFrames().size(); ++i)
        o << "    " << r->getFrames()[i];
    o << endl;
    return o;
}

    // STATION //

Station::Station(const String& nm) {
    rep = new StationRep(*this, nm.c_str());
}
Station::Station(const Station& s) : rep(0) {
    if (s.rep) {
        rep = new StationRep(*s.rep);
        rep->removePlacement();
    }
}
Station& Station::operator=(const Station& s) {
    if (this == &s) return *this;
    delete rep; rep=0;
    if (s.rep) {
        rep = new StationRep(*s.rep);
        rep->removePlacement();
    }
    return *this;
}
Station::~Station() { delete rep; }
void Station::setPlacement(const Frame& f, const Vec3& v) {
    if (!rep) rep = new StationRep(*this);
    rep->setPlacement(f,v);
}
std::ostream& operator<<(std::ostream& o, const Station& s) {
    const StationRep* r = StationRep::getRep(s);
    o << "Station ";
    if (r) o << r->getName() << ": " << r->getPlacement();
    else o << "at 0x" << &s << " has NULL rep.";
    o << endl;
    return o;
}

    // DIRECTION //

Direction::Direction(const String& nm)
  { rep = new DirectionRep(*this, nm.c_str()); }
Direction::Direction(const Direction& d) : rep(0) {
    if (d.rep) {
        rep = new DirectionRep(*d.rep);
        rep->removePlacement();
    }
}
Direction& Direction::operator=(const Direction& d) {
    if (this == &d) return *this;
    delete rep; rep=0;
    if (d.rep) {
        rep = new DirectionRep(*d.rep);
        rep->removePlacement();
    }
    return *this;
}
Direction::~Direction() { delete rep; }
void Direction::setPlacement(const Frame& f, const Vec3& v) {
    if (!rep) rep = new DirectionRep(*this);
    rep->setPlacement(f,v);
}
std::ostream& operator<<(std::ostream& o, const Direction& d) {
    const DirectionRep* r = DirectionRep::getRep(d);
    o << "Direction ";
    if (r) o << r->getName() << ": " << r->getPlacement();
    else o << "at 0x" << &d << " has NULL rep.";
    o << endl;
    return o;
}

    // MASS ELEMENT //

MassElement::MassElement(const MassElement& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->removePlacement();
    }
}
MassElement& MassElement::operator=(const MassElement& src) {
    if (this == &src) return *this;
    delete rep; rep=0;
    if (src.rep) {
        rep = src.rep->clone();
        rep->removePlacement();
    }
    return *this;
}

MassElement::~MassElement() {
    delete rep; rep=0;
}

PointMassElement::PointMassElement() {
    MassElementRep::setRep(*this, new PointMassElementRep(*this));
}

PointMassElement::PointMassElement(const String& nm) {
    PointMassElementRep* p = new PointMassElementRep(*this);
    p->setName(nm.c_str());
    MassElementRep::setRep(*this, p);
}

PointMassElement::PointMassElement(const Real& m) {
    PointMassElementRep* p = new PointMassElementRep(*this);
    p->setMass(m);
    MassElementRep::setRep(*this, p);
}

PointMassElement::PointMassElement(const Real& m, const String& nm) {
    PointMassElementRep* p = new PointMassElementRep(*this);
    p->setName(nm.c_str());
    p->setMass(m);
    MassElementRep::setRep(*this, p);
}

    // JOINT //
Joint::Joint(const String& nm)
  { rep = new JointRep(*this, nm.c_str()); }
Joint::~Joint() { delete rep; }

    // MULTIBODY SYSTEM //
MultibodySystem::MultibodySystem(const String& nm)
  { rep = new MultibodySystemRep(*this, nm.c_str()); }
MultibodySystem::~MultibodySystem() { delete rep; }

    // BODY //
Body::Body(const String& nm) : Frame(nm) 
  { rep = new BodyRep(*this); }
Body::~Body() { delete rep; }

    // RIGIDBODY //
RigidBody::RigidBody(const String& nm) : Body(nm)
  { rep = new RigidBodyRep(*this); }
RigidBody::~RigidBody() { }

    // DEFORMABLE BODY //
DeformableBody::DeformableBody(const String& nm) : Body(nm)
  { rep = new DeformableBodyRep(*this); }
DeformableBody::~DeformableBody() { }

    // MULTIBODY //
Multibody::Multibody(const String& nm) : Body(nm)
  { rep = new MultibodyRep(*this); }
Multibody::~Multibody() { }


} // namespace simtk
