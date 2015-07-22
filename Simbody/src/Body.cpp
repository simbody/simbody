/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**@file
Private implementation of Body and its built-in subclasses. **/

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"

#include "BodyRep.h"

namespace SimTK {

//==============================================================================
//                                    BODY
//==============================================================================
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
            rep = src.rep->clone(); // create a new object
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

Body::Body(const MassProperties& m) : rep(0) {
    rep = new Body::Rigid::RigidRep(m);
    rep->setMyHandle(*this);
}


int Body::addDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
    return updRep().addDecoration(X_BD, g);
}
int Body::getNumDecorations() const
{   return (int)getRep().decorations.size(); }
const DecorativeGeometry& Body::getDecoration(int i) const
{   return getRep().decorations[i]; }
// Allow writable access on const Body since just a decoration.
DecorativeGeometry& Body::updDecoration(int i) const
{   return const_cast<DecorativeGeometry&>(getDecoration(i)); }


int Body::addContactSurface(const Transform&      X_BS,
                            const ContactSurface& surface) {
    BodyRep& rep = updRep();
    const int nxt = (int)rep.surfaces.size();
    rep.surfaces.push_back(std::make_pair(X_BS, surface));
    rep.surfaces.back().second.setIndexOnBody(nxt);
    return nxt;
}
int Body::getNumContactSurfaces() const
{   return (int)getRep().surfaces.size(); }
const ContactSurface& Body::getContactSurface(int i) const
{   return getRep().surfaces[i].second; }
const Transform& Body::getContactSurfaceTransform(int i) const
{   return getRep().surfaces[i].first; }
ContactSurface& Body::updContactSurface(int i)
{   return updRep().surfaces[i].second; }
Transform& Body::updContactSurfaceTransform(int i)
{   return updRep().surfaces[i].first; }

const MassProperties& Body::getDefaultRigidBodyMassProperties() const {
    return getRep().getDefaultRigidBodyMassProperties();
}

Body& Body::setDefaultRigidBodyMassProperties(const MassProperties& m) {
    updRep().setDefaultRigidBodyMassProperties(m);
    return *this;
}



//==============================================================================
//                               BODY::RIGID
//==============================================================================
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
    return static_cast<const Rigid&>(b);
}
Body::Rigid& Body::Rigid::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return static_cast<Rigid&>(b);
}
const Body::Rigid::RigidRep& Body::Rigid::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const RigidRep&>(*rep);
}
Body::Rigid::RigidRep& Body::Rigid::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<RigidRep&>(*rep);
}



//==============================================================================
//                              BODY::GROUND
//==============================================================================
Body::Ground::Ground() {
    rep = new GroundRep();
    rep->setMyHandle(*this);
}

bool Body::Ground::isInstanceOf(const Body& b) {
    return GroundRep::isA(b.getRep());
}
const Body::Ground& Body::Ground::downcast(const Body& b) {
    assert(isInstanceOf(b));
    return static_cast<const Ground&>(b);
}
Body::Ground& Body::Ground::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return static_cast<Ground&>(b);
}
const Body::Ground::GroundRep& Body::Ground::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const GroundRep&>(*rep);
}
Body::Ground::GroundRep& Body::Ground::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<GroundRep&>(*rep);
}



//==============================================================================
//                              BODY::MASSLESS
//==============================================================================
Body::Massless::Massless() {
    rep = new MasslessRep();
    rep->setMyHandle(*this);
}

bool Body::Massless::isInstanceOf(const Body& b) {
    return MasslessRep::isA(b.getRep());
}
const Body::Massless& Body::Massless::downcast(const Body& b) {
    assert(isInstanceOf(b));
    return static_cast<const Massless&>(b);
}
Body::Massless& Body::Massless::updDowncast(Body& b) {
    assert(isInstanceOf(b));
    return static_cast<Massless&>(b);
}
const Body::Massless::MasslessRep& Body::Massless::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const MasslessRep&>(*rep);
}
Body::Massless::MasslessRep& Body::Massless::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<MasslessRep&>(*rep);
}

} // namespace SimTK

