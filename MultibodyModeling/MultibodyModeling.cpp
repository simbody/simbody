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
using std::endl;
using std::ostream;

namespace simtk {

Placement::Placement() : rep(0) { }
Placement::~Placement() { delete rep; }
bool Placement::hasPlacement() const
  { return rep && rep->hasPlacement(); }
const Frame&
Placement::getParentFrame() const
  { assert(rep); return rep->getParentFrame(); }

Frame::Frame() : rep(0) { }
Frame::Frame(const String& nm)
  { rep = new FrameRep(*this, nm.c_str()); }
Frame::~Frame() { delete rep; }

ostream& operator<<(ostream& o, const Frame& f) {
    o << "Frame ";
    if (f.rep) o << f.rep->getFullName();
    else o << "at " << &f << "has NULL rep.";
    o << endl;
    return o;
}

Station::Station() : rep(0) { }
Station::Station(const String& nm)
  { rep = new StationRep(*this, nm.c_str()); }
Station::~Station() { delete rep; }

Direction::Direction() : rep(0) { }
Direction::Direction(const String& nm)
  { rep = new DirectionRep(*this, nm.c_str()); }
Direction::~Direction() { delete rep; }

MassElement::MassElement() : rep(0) { }
MassElement::MassElement(const String& nm)
  { rep = new MassElementRep(*this, nm.c_str()); }
MassElement::~MassElement() { delete rep; }

Joint::Joint() : rep(0) { }
Joint::Joint(const String& nm)
  { rep = new JointRep(*this, nm.c_str()); }
Joint::~Joint() { delete rep; }

MultibodySystem::MultibodySystem() : rep(0) { }
MultibodySystem::MultibodySystem(const String& nm)
  { rep = new MultibodySystemRep(*this, nm.c_str()); }
MultibodySystem::~MultibodySystem() { delete rep; }

Body::Body() : Frame(), rep(0) { }
Body::Body(const String& nm) : Frame(nm) 
  { rep = new BodyRep(*this); }
Body::~Body() { delete rep; }

RigidBody::RigidBody() : Body() { }
RigidBody::RigidBody(const String& nm) : Body(nm)
  { rep = new RigidBodyRep(*this); }
RigidBody::~RigidBody() { }

DeformableBody::DeformableBody() : Body() { }
DeformableBody::DeformableBody(const String& nm) : Body(nm)
  { rep = new DeformableBodyRep(*this); }
DeformableBody::~DeformableBody() { }

Multibody::Multibody() : Body() { }
Multibody::Multibody(const String& nm) : Body(nm)
  { rep = new MultibodyRep(*this); }
Multibody::~Multibody() { }


} // namespace simtk
