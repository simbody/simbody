#ifndef SimTK_SIMBODY_COMPLIANT_CONTACT_SUBSYSTEM_H_
#define SimTK_SIMBODY_COMPLIANT_CONTACT_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MultibodySystem;
class SimbodyMatterSubsystem;
class ContactTrackerSubsystem;
class ContactForceGenerator;
class Contact;
class ContactForce;
class ContactPatch;



//==============================================================================
//                        COMPLIANT CONTACT SUBSYSTEM
//==============================================================================
/** This is a force subsystem that implements a compliant contact model to
respond to Contact objects as detected by a ContactTrackerSubsystem. The
subsystem contains an extendable collection of ContactForceGenerator 
objects, one per type of Contact. For example, a point contact would be
handled by a different force generator than would a mesh contact. **/
class SimTK_SIMBODY_EXPORT CompliantContactSubsystem : public ForceSubsystem {
public:
/** Default constructor creates an empty handle. **/
CompliantContactSubsystem() {}

/** Add a new CompliantContactSubsystem to the indicated MultibodySystem,
specifying the ContactTrackerSubsystem that we will use to obtain the list
of ActiveContacts from which we will generate compliant forces. The
MultibodySystem takes over ownership of the new subsystem object, but the
handle we construct here retains a reference to it. **/
CompliantContactSubsystem(MultibodySystem&, 
                          const ContactTrackerSubsystem&);

/** Get the transition velocity (vt) of the friction model. **/
Real getTransitionVelocity() const;
/** Set the transition velocity (vt) of the friction model. **/
void setTransitionVelocity(Real vt);
/** Get a precalculated 1/vt to avoid expensive runtime divisions. **/
Real getOOTransitionVelocity() const;

/** Specify whether to track energy dissipated by contacts. This permits 
conservation of energy tests to be performed in the presence of dissipative
contact elements. If you enable this, a state variable (z) will be allocated
to integrate instantaneous dissipated power. By default energy dissipation is
not tracked and an exception is thrown if it is requested. This is a topological
change, meaning you'll have to call realizeTopology() and get a new State if
you change the setting. 
@see getDissipatedEnergy(),setDissipatedEnergy(),getTrackDissipatedEnergy() **/
void setTrackDissipatedEnergy(bool shouldTrack);
/** Obtain the current setting of the "track dissipated energy" flag. You
must not call getDissipatedEnergy() or setDissipatedEnergy() if this flag
has not been set.
@see getDissipatedEnergy(),setDissipatedEnergy(),setTrackDissipatedEnergy() **/
bool getTrackDissipatedEnergy() const;

/** Determine how many of the active Contacts are currently generating
contact forces. You can call this at Velocity stage or later; the contact
forces will be realized first if necessary before we report how many there 
are. **/
int getNumContactForces(const State& state) const;
/** For each active Contact, get a reference to the most recently calculated
force there; the ContactId that produced this force is available from the
referenced ContactForce object. The ContactForce object is measured and 
expressed in the Ground frame. You can request this response at Velocity stage 
or later; the contact forces will be realized first if necessary. 
@see getContactForceById() **/
const ContactForce& getContactForce(const State& state, int n) const;
/** Get a reference to the ContactForce currently being produced by a particular 
ContactId; if that Contact is not currently producing forces then a reference to 
an invalid ContactForce is returned (you can check with isValid()). The 
ContactForce object is measured and expressed in the Ground frame. You can call 
this at Velocity stage or later; the contact forces (not just this one) will be 
realized first if necessary.  
@see getContactForce(), calcContactPatchDetailsById() **/
const ContactForce& getContactForceById(const State& state, ContactId id) const;

/** Calculate detailed information about a particular active contact patch, 
including deformed geometric information, pressure and friction force 
distribution, and resultant forces and moments. This detailed information is 
only calculated when requested because it may be expensive for some 
ContactForceGenerators; for simulation purposes only the resultants are needed. 
Resultant forces may be obtained cheaply by calling getContactForceById(); don't 
call this method unless you really need these details. The returned information 
will overwrite anything that is already in the supplied ContactPatch object; the 
contact point and vectors are expressed in Ground. You can call this operator at 
Velocity stage or higher. The result is calculated here and not saved internally.

@pre \a state realized to Stage::Velocity
@param[in]      state   
    The state from whose active contacts we are selecting.
@param[in]      id      
    The ContactId of the Contact whose patch details are to be returned. This
    refers to a Contact within the ContactTrackerSubsystem associated with this
    ContactForceSubsystem.
@param[out]     patch   
    The object that will be overwritten with the result or set invalid
    if the given ContactId is not currently producing a patch. Results are
    returned in the Ground frame.
@return
    True if the indicated Contact is currently generating contact forces;
    false otherwise (in which case \a patch is also marked invalid).

@see getContactForceById() **/
bool calcContactPatchDetailsById(const State&   state,
                                 ContactId      id,
                                 ContactPatch&  patch) const;

/** Obtain the total amount of energy dissipated by all the contact responses
that were generated by this subsystem since some arbitrary starting point. This
information is available only if you have requested tracking by calling
setTrackDissipatedEnergy(). 

This is the time integral of all the power dissipated during any of the contacts
by material dissipative and friction forces. For a system whose only 
non-conservative forces are contacts, the sum of potential, kinetic, and 
dissipated energies should be conserved with an exception noted below. This is 
particularly useful for debugging new ContactForceGenerators. This is a 
State variable so you can obtain its value any time after it is allocated.
@pre \a state realized to Stage::Model
@param[in]          state    
    The State from which to obtain the current value of the dissipated energy.
@return
    The total dissipated energy (a nonnegative scalar). 

An error message will be thrown if you call this method without having 
turned on energy tracking via setTrackDissipatedEnergy().

@note The Hunt and Crossley dissipation model used by many contact force
generators can occasionally detect that a body is being "yanked" out of 
a contact (this is more likely with very large dissipation coefficients). 
To track all the energy in that case we would have to allow the
surfaces to stick together, which is unreasonable. Instead, some of the 
energy will be lost to unmodeled effects because the surface is unable to
transfer energy back to the bodies and will instead vibrate or ring until
the energy is dissipated. 

@see setTrackDissipatedEnergy(), setDissipatedEnergy() **/
Real getDissipatedEnergy(const State& state) const;

/** Set the accumulated dissipated energy to an arbitrary value. This is only
permitted if you have requested tracking of dissipated energy by calling
setTrackDissipatedEnergy(). 

Typically this method is used only to reset the dissipated energy to zero, but 
non-zero values can be useful if you are trying to match some existing data or
continuing a simulation. This is a State variable so you can set its 
value any time after it is allocated.
@pre \a state realized to Stage::Model
@param[in,out]      state    
    The State whose dissipated energy variable for this subsystem is to
    be modified.
@param[in]          energy   
    The new value for the accumulated dissipated energy (must be a 
    nonnegative scalar). 

An error message will be thrown if you call this method without having 
turned on energy tracking via setTrackDissipatedEnergy().    
    
@see getDissipatedEnergy(), setTrackDissipatedEnergy() **/
void setDissipatedEnergy(State& state, Real energy) const;


/** Attach a new generator to this subsystem as the responder to be used when
we see the kind of Contact type for which this generator is defined,
replacing the previous generator for this Contact type if there was one. The 
subsystem takes over ownership of the generator; don't delete it yourself. **/
void adoptForceGenerator(ContactForceGenerator* generator);

/** Attach a new generator to this subsystem as the responder to be used when
we see a Contact type for which no suitable generator has been defined,
replacing the previous default generator type if there was one. The 
subsystem takes over ownership of the generator; don't delete it yourself. **/
void adoptDefaultForceGenerator(ContactForceGenerator* generator);

/** Return true if this subsystem has a force generator registered that can
respond to this kind of Contact. **/
bool hasForceGenerator(ContactTypeId contact) const;

/** Return true if this subsystem has a force generator registered that can
be used for response to an unrecognized type of Contact (typically this will
be a "do nothing" or "throw an error" generator. **/
bool hasDefaultForceGenerator() const;

/** Return the force generator to be used for a Contact of the indicated
type. If no generator was registered for this type of contact, this will
be the default generator. **/
const ContactForceGenerator& 
    getContactForceGenerator(ContactTypeId contact) const; 

/** Return the force generator to be used for a Contact type for which no
suitable force generator has been registered. **/
const ContactForceGenerator& getDefaultForceGenerator() const; 

/** Get a read-only reference to the ContactTrackerSubsystem associated
with this CompliantContactSubsystem. This is the contact tracker that is
maintaining the list of contacts for which this subsystem will be providing
the response forces. **/
const ContactTrackerSubsystem& getContactTrackerSubsystem() const;

/** Every Subsystem is owned by a System; a CompliantContactSubsystem expects
to be owned by a MultibodySystem. This method returns a const reference
to the containing MultibodySystem and will throw an exception if there is
no containing System or it is not a MultibodySystem. **/
const MultibodySystem& getMultibodySystem() const;

/** @cond **/   // don't show in Doxygen docs
SimTK_PIMPL_DOWNCAST(CompliantContactSubsystem, ForceSubsystem);
/** @endcond **/

//--------------------------------------------------------------------------
                                 private:
class CompliantContactSubsystemImpl& updImpl();
const CompliantContactSubsystemImpl& getImpl() const;
};



//==============================================================================
//                               CONTACT FORCE
//==============================================================================
/** This is a simple class containing the basic force information for a 
single contact between deformable surfaces S1 and S2 mounted on rigid
bodies B1 and B2. Every contact interaction between two rigid bodies, however 
complex, can be expressed as a resultant that can be contained in this class 
and is sufficient for advancing a simulation. Optionally, you may be able to 
get more details about the deformed geometry and pressure distribution over the
patch but you have to ask for that separately because it can be expensive to
calculate or report.

The information stored here is:
  - A point in space at which equal and opposite forces will be applied to
    corresponding stations of the two interacting rigid bodies. This is called
    the "contact point".
  - The force vector to be applied there (to each body with opposite sign).
  - A moment vector to be applied to both bodies (with opposite signs).
  - The potential energy currently stored by the elasticity of the contacting
    materials.
  - The instantaneous power dissipation due to inelastic behavior such as 
    friction and internal material damping.

Points and vectors are measured and expressed in some assumed but unspecified
frame A, which might for example be the frame of one of the two surfaces, or
one of the two bodies, or Ground. Whenever a method returns a ContactForce
object it must document what frame is being used; typically that will be
Ground for end user use.

<h3>Definition of center of pressure</h3>

When the contact patch itself involves many distributed contact points, the
center of pressure might not be a particularly meaningful quantity but 
serves to provide enough information to proceed with a simulation, and to
display resultant contact forces, without requiring detailed patch geometry
information. We define the location r_c of the center of pressure like this:
@verbatim
           sum_i (r_i * |r_i X Fn_i|) 
     r_c = --------------------------
              sum_i |r_i X Fn_i|
@endverbatim
where r_i is the vector locating contact point i, F_i=Fn_i+Ft_i is the 
contact force at point i resolved into locally-normal and locally-tangential
components, and M_i is a pure moment if any generated by the contact force 
model as a result of the i'th contact point. Note that the locally-
tangent contact force Ft_i (presumably from friction) and pure moment M_i
do not contribute to the center of pressure calculation; we're choosing
the point that minimizes the net moment generated by normal ("pressure")
forces only. Note that the normal force Fn_i includes both stiffness and
dissipation contributions, so the center of pressure can be 
velocity-dependent. **/
class ContactForce {
public:
/** Default constructor has invalid contact id, other fields garbage. **/
ContactForce() {} // invalid

/** Construct with values for all fields. Contact point and force must
be in the same but unspecified frame (typically Ground), and force and 
moment must be as applied at the contact point; we can't check so don't
mess up. **/
ContactForce(ContactId id, const Vec3& contactPt,
             const SpatialVec& forceOnSurface2,
             Real potentialEnergy, Real powerLoss)
:   m_contactId(id), m_contactPt(contactPt), 
    m_forceOnSurface2(forceOnSurface2),
    m_potentialEnergy(potentialEnergy), m_powerLoss(powerLoss) {}

/** Return the ContactId of the Contact that generated this ContactForce. **/
ContactId getContactId() const {return m_contactId;}
/** This is the point at which this contact element applies equal and opposite
forces to the two bodies, chosen as the center of pressure for composite
contact patches. **/
const Vec3& getContactPoint() const {return m_contactPt;}
/** Get the total spatial force applied to body 2 at the contact point
(that is, a force and a moment); negate this to find the force applied to 
body 1 at the same point. **/
const SpatialVec& getForceOnSurface2() const {return m_forceOnSurface2;}
/** Get the amount of potential energy currently stored in the deformation of
this contact patch. **/
Real getPotentialEnergy() const {return m_potentialEnergy;}
/** Get the energy dissipation rate (power loss) due to the deformation rate 
and friction losses for this contact patch (this is a signed value with 
positive indicating dissipation). **/
Real getPowerDissipation() const {return m_powerLoss;}

/** Replace the current contents of this ContactForce object with the
given arguments. This provides values for all fields. Contact point and 
force must be in the same but unspecified frame (typically Ground), and 
force and moment must be as applied at the contact point; we can't check
so don't mess up. **/
void setTo(ContactId id, const Vec3& contactPt,
           const SpatialVec& forceOnSurface2,
           Real potentialEnergy, Real powerLoss)
{   m_contactId         = id; 
    m_contactPt         = contactPt;
    m_forceOnSurface2   = forceOnSurface2;
    m_potentialEnergy   = potentialEnergy;
    m_powerLoss         = powerLoss; }

/** Change the ContactId contained in this ContactForce object. **/
void setContactId(ContactId id) {m_contactId=id;}
/** Change the contact point contained in this ContactForce object. **/
void setContactPoint(const Vec3& contactPt) {m_contactPt=contactPt;}
/** Change the value stored in this ContactForce object for the spatial 
force being applied on surface 2. **/
void setForceOnSurface2(const SpatialVec& forceOnSurface2) 
{   m_forceOnSurface2=forceOnSurface2; }
/** Change the value stored for potential energy in this ContactForce object. **/
void setPotentialEnergy(Real potentialEnergy) 
{   m_potentialEnergy=potentialEnergy; }
/** Change the value stored for power loss in this ContactForce object. **/
void setPowerDissipation(Real powerLoss) {m_powerLoss=powerLoss;}

/** Restore the ContactForce object to its default-constructed state with
an invalid contact id and garbage for the other fields. **/
void clear() {m_contactId.invalidate();}
/** Return true if this contact force contains a valid ContactId. **/
bool isValid() const {return m_contactId.isValid();}

/** This object is currently in an assumed frame A; given a transform from
another frame B to A we'll re-measure and re-express this in B. Cost is
48 flops. Note that this doesn't change the moment because we're not moving
the force application point physically; just remeasuring it from B. **/
void changeFrameInPlace(const Transform& X_BA) {
    m_contactPt         = X_BA*m_contactPt;        // shift & reexpress in B
    m_forceOnSurface2   = X_BA.R()*m_forceOnSurface2;      // reexpress in B
}

private:
ContactId       m_contactId;            // Which Contact produced this force?
Vec3            m_contactPt;            // In some frame A
SpatialVec      m_forceOnSurface2;      // at contact pt, in A; neg. for Surf1
Real            m_potentialEnergy;      // > 0 when due to compression
Real            m_powerLoss;            // > 0 means dissipation
};

// For debugging.
inline std::ostream& operator<<(std::ostream& o, const ContactForce& f) {
    o << "ContactForce for ContactId " << f.getContactId() << " (ground frame):\n";
    o << "  contact point=" << f.getContactPoint() << "\n";
    o << "  force on surf2 =" << f.getForceOnSurface2() << "\n";
    o << "  pot. energy=" << f.getPotentialEnergy() 
      << "  powerLoss=" << f.getPowerDissipation();
    return o << "\n";
}

//==============================================================================
//                              CONTACT DETAIL
//==============================================================================
/** This provides deformed geometry and force details for one element of a
contact patch that may be composed of many elements. Every element generates 
equal and opposite forces on both surfaces so we can use a consistent 
patch-centered convention for reporting details regardless of the original 
source of the element.

<h3>Deformed patch geometry</h3>

Vector results are expressed in Ground. We return 
  - the contact point C where equal and opposite forces are applied to 
    <em>material points</em> of the two surfaces.
  - the contact normal direction n pointing \e away from surface1's exterior 
    and towards surface2's interior
  - the slip velocity v of surface 2 relative to surface 1, measured at C
    and in the plane perpendicular to n
  - the force applied to body 2 by body 1; that is, its normal component
    is in direction n and its tangential component opposes velocity v

It is important to note that we are interested in the motion (deformations) of 
<em>material points</em> here, \e not the motion of the <em>contact point</em>
which is a non-material concept and thus not directly involved in producing 
forces. The patch area is provided separately and its significance also depends
on the model.

<h3>Elasticity and dissipation</h3>

We report the combined deformation (>= 0) of the two surfaces at C, along the 
normal n. We do not attempt to report angular deformation of the element. Elastic 
deformation rate is also provided. This is the rate that the material points 
currently at C are compressing or relaxing along the normal. That is, we are not
including the fact that the contact point C may be moving along 
the surfaces, we just report what is happening to the material in its own 
frame. Note that we do not attempt to report angular deformation rate.

<h3>Slipping and friction</h3>

Slip velocity measures the relative velocity of the body stations that are
instantaneously coincident with C. This will be in the plane for which n is the 
normal, but it is expressed in Ground so is a 3d vector.

<h3>%Force and pressure</h3>

%Contact force is the force being 
applied to surface 2 by surface 1 at the contact point C. This includes
the normal force due to elasticity and dissipation, and
the tangential force due to friction and to tangential elasticity and
dissipation, if any.

Peak pressure is a scalar providing the worst-case pressure present somewhere
in the patch, if the model can provide that. Otherwise it will be the
normal force divided by the patch area to give the average pressure across
the patch. 

<h3>Energy and power</h3>

Potential energy is the amount of energy currently stored in the elastic
deformation of this element. Summing this over all the contact detail 
elements should yield the same value as is reported as the resultant potential
energy for this contact.

Power loss is the rate at which energy is being lost due to dissipation and
to friction, but not due to elastic deformation because that is contributing
to potential energy and we expect to get it back. Summing this over all the
contact detail elements should yield the same value as is reported as the
resultant power loss for this contact. **/
class ContactDetail {
public:
/** This is the point at which this contact element applies equal and opposite
forces to the two bodies. **/
const Vec3& getContactPoint() const {return m_contactPt;}
/** This is the normal direction for this contact element, pointing away
from body 1's exterior and towards body 2's interior, that is, in the
direction of the normal force applied to body 2 by body1. **/ 
const UnitVec3& getContactNormal() const {return m_patchNormal;}
/** Get the relative slip velocity between the bodies at the contact point,
as body 2's velocity in body 1. This is a vector in the contact plane, that is,
it is perpendicular to the contact normal. **/
const Vec3& getSlipVelocity() const {return m_slipVelocity;}
/** Get the total force applied to body 2 at the contact point by this 
contact element; negate this to find the force on body 1 at the same point. **/
const Vec3& getForceOnSurface2() const {return m_forceOnSurface2;}
/** Get the total normal material deformation at the contact point; this is the 
sum of the deformations of the two surfaces there. This is sometimes called the
"approach" of the two bodies and represents the amount of overlap of the
\e undeformed surfaces at this point; that is, they have to deform this much
in order not to overlap. **/
Real getDeformation() const {return m_deformation;}
/** Get the instantaneous rate at which the material at the contact point is
deforming; this is the material derivative of the deformation and does not
reflect that fact that the contact point itself is changing with time.
Energy dissipation in the material depends on the rate at which the \e material
is deforming, it does not depend on contact point changes. **/
Real getDeformationRate() const {return m_deformationRate;}
/** This is the surface area represented by this contact element. **/
Real getPatchArea() const {return m_patchArea;}
/** This is the peak pressure on this element; typically it is just the
normal force divided by the patch area since most contact models assume
constant force across a contact element's patch. **/
Real getPeakPressure() const {return m_peakPressure;}
/** Get the amount of potential energy currently stored in the deformation of
this contact element. **/
Real getPotentialEnergy() const {return m_potentialEnergy;}
/** Get the energy dissipation rate (power loss) due to the deformation rate and
friction losses for this element (this is a signed value with positive indicating
dissipation). **/
Real getPowerDissipation() const {return m_powerLoss;}

    
/** This object is currently in an assumed frame A; given a transform from
another frame B to A we'll re-measure and re-express this in B. Cost is
63 flops. **/
void changeFrameInPlace(const Transform& X_BA) {
    const Rotation& R_BA = X_BA.R();
    m_contactPt       = X_BA*m_contactPt;       // shift & reexpress in B (18 flops)
    m_patchNormal     = R_BA*m_patchNormal;     // reexpress only       (3*15 flops)
    m_slipVelocity    = R_BA*m_slipVelocity;    //      "
    m_forceOnSurface2 = R_BA*m_forceOnSurface2; //      "
}

/** Assuming that this object is currently reporting surface 2 information in 
frame A, here we want to both change the frame to B and swap which surface is to 
be considered as surface 2. Cost is 72 flops. **/
void changeFrameAndSwitchSurfacesInPlace(const Transform& X_BA) {
    const Rotation& R_BA = X_BA.R();
    m_contactPt       = X_BA*m_contactPt;       // shift & reexpress in B (18 flops)
    m_patchNormal     = R_BA*(-m_patchNormal);  // reverse & reexpress  (3*18 flops)
    m_slipVelocity    = R_BA*(-m_slipVelocity); //          "
    m_forceOnSurface2 = R_BA*(-m_forceOnSurface2); //       "
}

Vec3            m_contactPt;            // location of contact point C in A
UnitVec3        m_patchNormal;          // points outwards from body 1, exp. in A
Vec3            m_slipVelocity;         // material slip rate, perp. to normal, in A
Vec3            m_forceOnSurface2;      // applied at C, -force to surf1
Real            m_deformation;          // total normal compression (approach)
Real            m_deformationRate;      // d/dt deformation, w.r.t. A frame
Real            m_patchArea;            // >= 0
Real            m_peakPressure;         // > 0 in compression
Real            m_potentialEnergy;      // > 0 when due to compression
Real            m_powerLoss;            // > 0 means dissipation
};

// For debugging.
inline std::ostream& operator<<(std::ostream& o, const ContactDetail& d) {
    o << "ContactDetail (ground frame):\n";
    o << "  contact point=" << d.m_contactPt << "\n";
    o << "  contact normal=" << d.m_patchNormal << "\n";
    o << "  slip velocity=" << d.m_slipVelocity << "\n";
    o << "  force on surf2 =" << d.m_forceOnSurface2 << "\n";
    o << "  deformation=" << d.m_deformation 
      << "  deformation rate=" << d.m_deformationRate << "\n";
    o << "  patch area=" << d.m_patchArea 
      << "  peak pressure=" << d.m_peakPressure << "\n";
    o << "  pot. energy=" << d.m_potentialEnergy << "  powerLoss=" << d.m_powerLoss;
    return o << "\n";
}



//==============================================================================
//                               CONTACT PATCH
//==============================================================================
/** A ContactPatch is the description of the forces and the deformed shape of
the contact surfaces that result from compliant contact interactions. This
should not be confused with the Contact object that describes only the 
overlap in \e undeformed surface geometry and knows nothing of forces or
materials. Although there are several qualitatively different kinds of 
compliant contact models, we assume that each can be described by some number 
of contact "elements" and report the detailed results in terms 
of those elements. Depending on the model, the elements may be associated with
the surfaces or may be associated with the patch. Either way, every element 
applies equal and opposite forces to both surfaces. There is not necessarily 
any direct correspondence between elements on one surface with elements on the 
other.

A Hertz contact will have only a single element that belongs to the patch,
while an elastic foundation contact will have many, with each element
associated with one or the other surface. Only the elements that are currently 
participating in contact will have entries here; the surface and element id is 
stored with each piece of ContactDetail information. There is also some basic 
information needed to advance the simulation and that is common to all contact
patch types and stored as a ContactForce resultant. **/
class SimTK_SIMBODY_EXPORT ContactPatch {
public:
void clear() {m_resultant.clear(); m_elements.clear();}
bool isValid() const {return m_resultant.isValid();}
const ContactForce& getContactForce() const {return m_resultant;}
int getNumDetails() const {return (int)m_elements.size();}
const ContactDetail& getContactDetail(int n) const {return m_elements[n];}

/** This object is currently in an assumed frame A; given a transform from
another frame B to A we'll re-measure and re-express this in B. This is an
expensive operation. **/
void changeFrameInPlace(const Transform& X_BA) {
    m_resultant.changeFrameInPlace(X_BA);
    for (unsigned i=0; i<m_elements.size(); ++i)
        m_elements[i].changeFrameInPlace(X_BA);
}

ContactForce            m_resultant;
Array_<ContactDetail>   m_elements;
};



//==============================================================================
//                          CONTACT FORCE GENERATOR
//==============================================================================
/** A ContactForceGenerator implements an algorithm for responding to overlaps 
or potential overlaps between pairs of ContactSurface objects, as detected by
a ContactTrackerSubsystem. This class is used internally by 
CompliantContactSubsystem and there usually is no reason to access it directly.
The exception is if you are defining a new Contact subclass (very rare). In 
that case, you will also need to define one or more ContactForceGenerators to 
respond to Contacts with the new type, then register it with the 
CompliantContactSubsystem. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator {
public:
// Reasonably good physically-based compliant contact models.
class ElasticFoundation;        // for TriangleMeshContact
class HertzCircular;            // for CircularPointContact
class HertzElliptical;          // for EllipticalPointContact

// Penalty-based models enforcing non-penetration but without attempting
// to model the contacting materials physically.
class BrickHalfSpacePenalty;    // for BrickHalfSpaceContact

// These are for response to unknown ContactTypeIds.
class DoNothing;     // do nothing if called
class ThrowError;    // throw an error if called

/** Base class constructor for use by the concrete classes. **/
explicit ContactForceGenerator(ContactTypeId type): m_contactType(type) {}

/** Return the ContactTypeId handled by this force generator. ContactTypeId(0)
is reserved and is used here for fallback force generators that deal with
unrecognized ContactTypeIds. **/
ContactTypeId getContactTypeId() const {return m_contactType;}

const CompliantContactSubsystem& getCompliantContactSubsystem() const
{   assert(m_compliantContactSubsys); return *m_compliantContactSubsys; }
void setCompliantContactSubsystem(const CompliantContactSubsystem* sub)
{   m_compliantContactSubsys = sub; }

/** Base class destructor is virtual but does nothing. **/
virtual ~ContactForceGenerator() {}

/** The CompliantContactSubsystem will invoke this method on any 
active contact pair of the right Contact type for which there is overlapping 
undeformed geometry. The force generator is expected to calculate a point
in space where equal and opposite contact forces should be applied to the two
contacting rigid bodies, the potential energy currently stored in this 
contact, and the power (energy dissipation rate). State should be used
for instance info only; use position information from \a overlapping and
velocity information from the supplied arguments. That allows this method
to be used as an operator, for example to calculate potential energy when
velocities are not yet available. **/
virtual void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,  // relative surface velocity (S2 in S1)
    ContactForce&           contactForce) const = 0;

/** The CompliantContactSubsystem will invoke this method in response to a user
request for contact patch information; this returns force, potential energy, 
and power as above but may also require expensive computations that can be 
avoided in calcContactForce(). Don't use the state for position or
velocity information; the only allowed positions are in the Contact object
and the velocities are supplied explicitly. **/
virtual void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,  // relative surface velocity (S2 in S1)
    ContactPatch&           patch) const = 0;


//--------------------------------------------------------------------------
private:
    // This generator should be called only for Contact objects of the 
    // indicated type id.
    ContactTypeId                       m_contactType;
    // This is a reference to the owning CompliantContactSubsystem if any;
    // don't delete on destruction.
    const CompliantContactSubsystem*    m_compliantContactSubsys;
};




//==============================================================================
//                         HERTZ CIRCULAR GENERATOR
//==============================================================================

/** This ContactForceGenerator handles contact between non-conforming
objects that meet at a point and generate a circular contact patch; those
generate a CircularPointContact tracking object. Although this is just a special
case of elliptical contact we treat it separately so we can take advantage
of the significant simplifications afforded by circular contact. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::HertzCircular 
:   public ContactForceGenerator {
public:
HertzCircular() 
:   ContactForceGenerator(CircularPointContact::classTypeId()) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override;

void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override;
};



//==============================================================================
//                         HERTZ ELLIPTICAL GENERATOR
//==============================================================================

/** This ContactForceGenerator handles contact between non-conforming
objects that meet at a point and generate an elliptical contact patch; those
generate an EllipticalPointContact tracking object. For objects that are 
known to produce circular contact, use the specialized HertzCircular 
generator instead. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::HertzElliptical 
:   public ContactForceGenerator {
public:
HertzElliptical() 
:   ContactForceGenerator(EllipticalPointContact::classTypeId()) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override;

void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override;
};




//==============================================================================
//                         BRICK HALFSPACE GENERATOR
//==============================================================================

/** This ContactForceGenerator handles contact between a brick and a half-space.
TODO: generalize to convex mesh/half-space. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::BrickHalfSpacePenalty 
:   public ContactForceGenerator {
public:
BrickHalfSpacePenalty() 
:   ContactForceGenerator(BrickHalfSpaceContact::classTypeId()) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override;

void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override;
};



//==============================================================================
//                       ELASTIC FOUNDATION GENERATOR
//==============================================================================
/** This ContactForceGenerator handles contact between a TriangleMesh
and a variety of other geometric objects, all of which produce a
TriangleMeshContact tracking object. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::ElasticFoundation 
:   public ContactForceGenerator {
public:
ElasticFoundation() 
:   ContactForceGenerator(TriangleMeshContact::classTypeId()) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override;

void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override;

private:
void calcContactForceAndDetails
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce,
    Array_<ContactDetail>*  contactDetails) const;

void calcWeightedPatchCentroid
   (const ContactGeometry::TriangleMesh&    mesh,
    const std::set<int>&                    insideFaces,
    Vec3&                                   weightedPatchCentroid,
    Real&                                   patchArea) const;
                       
void processOneMesh
   (const State&                            state,
    const ContactGeometry::TriangleMesh&    mesh,
    const std::set<int>&                    insideFaces,
    const Transform&                        X_MO, 
    const SpatialVec&                       V_MO,
    const ContactGeometry&                  other,
    Real                                    meshDeformationFraction, // 0..1
    Real                                    areaScaleFactor, // >= 0
    Real k, Real c, Real us, Real ud, Real uv,
    const Vec3&                 resultantPt_M, // where to apply forces
    SpatialVec&                 resultantForceOnOther_M, // at resultant pt
    Real&                       potentialEnergy,
    Real&                       powerLoss,
    Vec3&                       weightedCenterOfPressure_M, // COP
    Real&                       sumOfAllPressureMoments,    // COP weight
    Array_<ContactDetail>*      contactDetails) const;
};




//==============================================================================
//                        DO NOTHING FORCE GENERATOR
//==============================================================================
/** This ContactForceGenerator silently does nothing. It can be used as a way
to explicitly ignore a certain ContactTypeId, or more commonly it is used as
the fallback generator for unrecognized ContactTypeIds. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::DoNothing 
:   public ContactForceGenerator {
public:
explicit DoNothing(ContactTypeId type = ContactTypeId(0)) 
:   ContactForceGenerator(type) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override
{   SimTK_ASSERT_ALWAYS(!"implemented",
        "ContactForceGenerator::DoNothing::calcContactForce() not implemented yet."); }
void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override
{   SimTK_ASSERT_ALWAYS(!"implemented",
        "ContactForceGenerator::DoNothing::calcContactPatch() not implemented yet."); }
};



//==============================================================================
//                       THROW ERROR FORCE GENERATOR
//==============================================================================
/** This ContactForceGenerator throws an error if it is every invoked. It can 
be used as a way to explicitly catch a certain ContactTypeId and complain, or 
more commonly it is used as the fallback generator for unrecognized 
ContactTypeIds. **/
class SimTK_SIMBODY_EXPORT ContactForceGenerator::ThrowError 
:   public ContactForceGenerator {
public:
explicit ThrowError(ContactTypeId type = ContactTypeId(0)) 
:   ContactForceGenerator(type) {}

void calcContactForce
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactForce&           contactForce) const override
{   SimTK_ASSERT_ALWAYS(!"implemented",
        "ContactForceGenerator::ThrowError::calcContactForce() not implemented yet."); }
void calcContactPatch
   (const State&            state,
    const Contact&          overlapping,
    const SpatialVec&       V_S1S2,
    ContactPatch&           patch) const override
{   SimTK_ASSERT_ALWAYS(!"implemented",
        "ContactForceGenerator::ThrowError::calcContactPatch() not implemented yet."); }
};

} // namespace SimTK

#endif // SimTK_SIMBODY_COMPLIANT_CONTACT_SUBSYSTEM_H_
