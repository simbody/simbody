#ifndef SimTK_SIMBODY_BODY_H_
#define SimTK_SIMBODY_BODY_H_

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

/** @file
This defines the API for the Body base class and concrete Body types like
Body::Rigid that are derived from it. These are handle classes; the
implementation is hidden. **/

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/ContactSurface.h"

#include <cassert>

namespace SimTK {

class DecorativeGeometry;

//==============================================================================
//                                    BODY
//==============================================================================
/** The Body class represents a reference frame that can be used to describe 
mass properties and geometry. These are in turn used to build MobilizedBodies 
which combine a body and a specified mobilizer that defines how the reference 
frame can move with respect to other MobilizedBodies. Attached geometric 
objects can serve as the basis for a variety of force-generating elements or
other algorithms that act on bodies.

Body is an abstract base class handle, with concrete classes defined for each 
kind of body. There are a set of built-in body types, with Body::Rigid the
most common. **/
class SimTK_SIMBODY_EXPORT Body {
public:
/** Default constructor creates an empty Body handle. **/
Body() : rep(0) { }
/** Destroy the handle and the body if this is the owner. **/
~Body();
/** Copy constructor is a deep copy; the new Body is separate from the
source Body. **/
Body(const Body& source);
/** Copy assignment is a deep copy; the original object is deleted if this
is the owner, then replaced with a \e copy of the source. **/
Body& operator=(const Body& source);

/** This is a default conversion from MassProperties to Body. It will result 
in a rigid body (concrete class Body::Rigid) being created using these as its 
default MassProperties. This is what allows you to provide MassProperties 
instead of a Body in the MobilizedBody constructors. **/
Body(const MassProperties& massProps);

/** Every type of Body should provide an initial set of rigid body mass 
properties defined at Topology stage (i.e., in the System rather than
the State). This is thus a Topology-stage change which will require a new
realizeTopology() before use. **/ 
Body& setDefaultRigidBodyMassProperties(const MassProperties&);

/** Get the default (that is, Topology stage) mass properties for this Body.
This may be overridden in a State if this Body has variable mass 
properties. **/
const MassProperties& getDefaultRigidBodyMassProperties() const;

/** Add a piece of decorative geometry fixed at some pose on this Body. 
This can be used for visualization of the Body's motion. Returns a small 
integer that can be used to identify this decoration in any copy of this
%Body. The supplied DecorativeGeometry object is copied and transformed by the
given Transform so that the actual geometry is always stored relative to the
body's frame. **/
int addDecoration(const Transform& X_BD, const DecorativeGeometry& geometry);

/** Convenience method for when the decorative geometry is to be placed at the
body frame. This is the same as addDecoration(Transform(),geometry). **/
int addDecoration(const DecorativeGeometry& geometry)
{   return addDecoration(Transform(), geometry); }

/** Obtain a count nd of how many pieces of DecorativeGeometry have been
attached to this Body. The indices i for individual decorations will be numbered
0 <= i < nd. **/
int getNumDecorations() const;

/** Get a read-only reference to the i'th piece of DecorativeGeometry that 
was added to this Body, with 0 <= i < getNumDecorations(). The index i is
the small integer that was returned by addDecoration(). **/
const DecorativeGeometry& getDecoration(int i) const;

/** Get a writable reference to the i'th piece of DecorativeGeometry that 
was added to this Body, with 0 <= i < getNumDecorations().  The index i is
the small integer that was returned by addDecoration(). Note that we allow
writable access to decorations even on a const Body -- these are after all
just decorations. **/
DecorativeGeometry& updDecoration(int i) const;

/** Create a new contact surface on a body and place it using the indicated
Transform. Returns a small integer that can be used to identify this contact 
surface in any copy of this %Body. Unlike decorations, the transform is kept 
separately rather than used to transform the copied surface. You can obtain it 
later with getContactSurfaceTransform(). **/
int addContactSurface(const Transform&          X_BS,
                      const ContactSurface&     shape); 

/** Convenience method for when the contact surface is to be placed at the
body frame. This is the same as addContactSurface(Transform(),shape). **/
int addContactSurface(const ContactSurface& shape)
{   return addContactSurface(Transform(), shape); }

/** Obtain the number of contact surfaces ns attached to this Body. The valid
body-local index values i will be 0 <= i < ns. **/
int getNumContactSurfaces() const;
/** Get a reference to the i'th contact surface on this body; be sure to get
the Transform also. **/
const ContactSurface& getContactSurface(int i) const;
/** Get the transform specifying the placement of the i'th contact surface
on this Body. **/
const Transform& getContactSurfaceTransform(int i) const;
/** Get write access to the i'th unique contact surface owned by this Body. This
is a Topology-stage change that will require a new realizeTopology() call if 
this Body is part of a System. **/
ContactSurface& updContactSurface(int i);
/** Get a writable reference to the transform specifying the placement of the 
i'th contact surface on this Body.  This is a Topology-stage change that will 
require a new realizeTopology() call if this Body is part of a System. **/
Transform& updContactSurfaceTransform(int i);

// These are the built-in Body types.
class Ground;     // infinitely massive
class Massless;   // just a reference frame
class Particle;   // point mass only; com=0; no inertia
class Linear;     // point masses along a line; scalar inertia
class Rigid;      // general rigid body
class Deformable; // base class for bodies with internal deformation coords

bool isOwnerHandle() const;
bool isEmptyHandle() const;

// Internal use only
class BodyRep; // local subclass
explicit Body(class BodyRep* r) : rep(r) { }
bool           hasRep() const {return rep!=0;}
const BodyRep& getRep() const {assert(rep); return *rep;}
BodyRep&       updRep() const {assert(rep); return *rep;}
void           setRep(BodyRep& r) {assert(!rep); rep = &r;}

protected:
class BodyRep* rep;
};



//==============================================================================
//                               BODY::RIGID
//==============================================================================
/** A general rigid body. This can represent a body with mass properties that 
are full, linear, inertialess (e.g. a point), or massless. **/
class SimTK_SIMBODY_EXPORT Body::Rigid : public Body {
public:
    /** Construct a rigid body with default mass properties which are
    (1,Vec3(0),Inertia(1,1,1)) **/
    Rigid(); 
    /** Construct a rigid body with the given mass properties; any set of
    mass properties is allowed since this is a general rigid body. **/
    explicit Rigid(const MassProperties&);

    Rigid& setDefaultRigidBodyMassProperties(const MassProperties& m) {
        (void)Body::setDefaultRigidBodyMassProperties(m);
        return *this;
    }

    class RigidRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Rigid, Body);
private:
    RigidRep&       updRep();
    const RigidRep& getRep() const;
};




//==============================================================================
//                               BODY::LINEAR
//==============================================================================
/** This is a rigid body in the shape of a line, which is inherently 
inertialess about its axis. Its mass properties may be modified later, but
only in such a way that the Body remains inertialess about one axis. **/
class SimTK_SIMBODY_EXPORT Body::Linear : public Body {
public:
    Linear(); // default mass properties (1,Vec3(0),Inertia(1,1,0))
    explicit Linear(const MassProperties&);

    Linear& setDefaultRigidBodyMassProperties(const MassProperties& m) {
        (void)Body::setDefaultRigidBodyMassProperties(m);
        return *this;
    }

    class LinearRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Linear, Body);
private:
    LinearRep&       updRep();
    const LinearRep& getRep() const;
};


// Particles aren't yet supported, so hide this in Doxygen.
/** @cond **/
//==============================================================================
//                              BODY::PARTICLE
//==============================================================================
/** This kind of body can only represent an inertialess point mass, with
mass center at (0,0,0) in the local frame. You can change mass later (even
to make it massless) but you can't move the mass center or add any inertia. **/
class SimTK_SIMBODY_EXPORT Body::Particle : public Body {
public:
    Particle(); // default mass properties (1,Vec3(0),Inertia(0))
    explicit Particle(const Real& mass);

    Particle& setDefaultRigidBodyMassProperties(const MassProperties& m) {
        (void)Body::setDefaultRigidBodyMassProperties(m);
        return *this;
    }

    class ParticleRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Particle, Body);
private:
    ParticleRep&       updRep();
    const ParticleRep& getRep() const;
};
/** @endcond **/



//==============================================================================
//                              BODY::MASSLESS
//==============================================================================
/** This is a Body that is constitutively massless (and inertialess); meaning 
that no amount of fiddling with it will ever give it any mass or inertia. **/
class SimTK_SIMBODY_EXPORT Body::Massless : public Body {
public:
    Massless();

    class MasslessRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Massless, Body);
private:
    MasslessRep&       updRep();
    const MasslessRep& getRep() const;
};



//==============================================================================
//                              BODY::GROUND
//==============================================================================
/** This is a Body representing something immobile, of effectively infinite
mass and inertia, that cannot be modified to be anything else. **/
class SimTK_SIMBODY_EXPORT Body::Ground : public Body {
public:
    Ground();

    class GroundRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Ground, Body);
private:
    GroundRep&       updRep();
    const GroundRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_BODY_H_



