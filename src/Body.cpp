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
 * Private implementation of Body and its built-in subclasses.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"

#include "BodyRep.h"

namespace SimTK {

    //////////
    // BODY //
    //////////

bool Body::isEmptyHandle() const {return rep==0;}
bool Body::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}

Body::~Body() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

// Copy constructor creates a new deep copy of the source object.
Body::Body(const Body& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

// Assignment puts a deep copy of the src Body's rep into the current handle.
Body& Body::operator=(const Body& src) {
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

Body::Body(const MassProperties& m) : rep(0) {
    rep = new Body::Rigid::RigidRep(m);
    rep->setMyHandle(*this);
}


Body& Body::addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
    updRep().addDecoration(X_BD, g);
    return *this;
}

const MassProperties& Body::getDefaultRigidBodyMassProperties() const {
    return getRep().getDefaultRigidBodyMassProperties();
}

Body& Body::setDefaultRigidBodyMassProperties(const MassProperties& m) {
    updRep().setDefaultRigidBodyMassProperties(m);
    return *this;
}

    /////////////////
    // BODY::RIGID //
    /////////////////

Body::Rigid::Rigid() {
    rep = new RigidRep();
    rep->setMyHandle(*this);
}

Body::Rigid::Rigid(const MassProperties& m) {
    rep = new RigidRep(m);
    rep->setMyHandle(*this);
}

bool Body::Rigid::isInstanceOf(const Body& b) {
    return RigidRep::isA(b.getRep());
}
const Body::Rigid& Body::Rigid::downcast(const Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<const Rigid&>(b);
}
Body::Rigid& Body::Rigid::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<Rigid&>(b);
}
const Body::Rigid::RigidRep& Body::Rigid::getRep() const {
    return dynamic_cast<const RigidRep&>(*rep);
}
Body::Rigid::RigidRep& Body::Rigid::updRep() {
    return dynamic_cast<RigidRep&>(*rep);
}

    //////////////////
    // BODY::GROUND //
    //////////////////

Body::Ground::Ground() {
    rep = new GroundRep();
    rep->setMyHandle(*this);
}

bool Body::Ground::isInstanceOf(const Body& b) {
    return GroundRep::isA(b.getRep());
}
const Body::Ground& Body::Ground::downcast(const Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<const Ground&>(b);
}
Body::Ground& Body::Ground::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<Ground&>(b);
}
const Body::Ground::GroundRep& Body::Ground::getRep() const {
    return dynamic_cast<const GroundRep&>(*rep);
}
Body::Ground::GroundRep& Body::Ground::updRep() {
    return dynamic_cast<GroundRep&>(*rep);
}


    ////////////////////
    // BODY::MASSLESS //
    ////////////////////

Body::Massless::Massless() {
    rep = new MasslessRep();
    rep->setMyHandle(*this);
}

bool Body::Massless::isInstanceOf(const Body& b) {
    return MasslessRep::isA(b.getRep());
}
const Body::Massless& Body::Massless::downcast(const Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<const Massless&>(b);
}
Body::Massless& Body::Massless::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return reinterpret_cast<Massless&>(b);
}
const Body::Massless::MasslessRep& Body::Massless::getRep() const {
    return dynamic_cast<const MasslessRep&>(*rep);
}
Body::Massless::MasslessRep& Body::Massless::updRep() {
    return dynamic_cast<MasslessRep&>(*rep);
}

} // namespace SimTK

