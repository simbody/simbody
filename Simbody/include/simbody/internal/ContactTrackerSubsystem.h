#ifndef SimTK_SIMBODY_CONTACT_TRACKER_SUBSYSTEM_H_
#define SimTK_SIMBODY_CONTACT_TRACKER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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
#include "simbody/internal/ContactSurface.h"

namespace SimTK {


class MultibodySystem;
class MobilizedBody;
class ContactTracker;
class Contact;
class ContactSnapshot;

/** This subsystem identifies and tracks potential contacts between the
mobilized bodies of a multibody system. It operates on the \e undeformed,
material-independent geometry of ContactSurface objects that have been
associated with those mobilized bodies, identifying and
characterizing pairwise geometric relationships that occur between surfaces,
and tracking the evolution of particular contacts through time. No physical
response to contact is generated here; this subsystem provides only a general
contact-tracking service that can be used by other parts of the MultibodySystem
to generate forces, impulses, visualizations, messages, noises, or whatever.
The goal here is to provide robust and extremely high performance
characterization of geometric interactions in a way that is useful for a
variety of contact response models.

The %ContactTrackerSubsystem maintains an evolving set of tracked Contact
objects throughout a simulation; for any given State of the system the subsystem
can calculate the current value of that set, which we call the "contact status"
of that state. The most recently known \e prior valid contact status is
maintained as part of the system state because the correct instantaneous
evaluation of contact status may require past information, and because we are
interested in discrete contact events such as initiation and breaking of contact
as well as ongoing interactions.

Each Contact being tracked represents the interaction between a unique pair of
ContactSurface objects; by definition there is at most one Contact between any
such pair (per %ContactTrackerSubsystem). The presence of a Contact in the set
does not necessarily mean its two surfaces are touching, just that their
proximity is interesting in some way. Each Contact is characterized as
impending, initiated, ongoing, broken, or separating. At each evaluation, the
last-known set of Contacts is used in conjuction with the current State to
determine a disposition for each of the tracked Contacts. Contacts that have
become boring are removed from the tracked set, and newly-interesting ones are
assigned a new ContactId and added to the tracked contacts set. The ContactId
persists as long as that Contact continues to remain in the set.

Ambiguous or impossible contact situations can occur during trial steps and
can be detected as long as we know the correct contact status in the immediate
past. Problems are more common when a simulation contains fast-moving objects,
small objects, or thin surfaces, and a relatively large time step size is
being attempted by the integrator. For example, we may find that we could have
missed a contact that may have occurred since the end of the last step
(pass-through), or can't determine whether a small object is deeply penetrated
through a thin surface or has simple "gone around the back". In such cases the
subsystem reports an error condition that will cause the integrator to reduce
the step size. At the start of a time stepping study, there is no past
information available to help disambiguate; in that case heuristics are used
to make a "best guess" at which surfaces are in contact; those can be
overridden by a knowledgable human.

As mentioned above, we track at most one Contact at a time between any pair
of ContactSurface objects. However, for some surface types a single Contact may
involve many geometric interactions; a mesh Contact, for example, may include
a list of all the mesh faces that overlap between the two surfaces, and these do
not have to be contiguous over the mesh surface. You may want to break up a
concave surface into several ContactSurface objects so that you can separately
track several Contacts between that surface and others.

A ContactSurface consists of a piece of surface geometry and a contact material
that describes the contact-related physical properties of the surface.
A body may have many ContactSurface objects attached to it, each with its own
associated contact material. Within an instance of the %ContactTrackerSubsystem,
each ContactSurface is assigned a unique ContactSurfaceIndex. If you have
more than one ContactTrackerSubsystem (not common), they will each have
independent sets of ContactSurface objects (and know nothing about one another),
so the (Subsystem,ContactSurfaceIndex) pair would be unique but not the index
alone. You can map between the assigned ContactSurfaceIndex and the actual
surface using methods of this class.

Within a %ContactTrackerSubsystem, all the contact surfaces are presumed to be
capable of interacting with one another unless they share membership in a
"contact clique". All surfaces attached to the same mobilized body are placed
together in a clique, so they will never interact. It is also common to create a
clique associated with each joint and place nearby contact surfaces on adjacent
mobilized bodies into that clique to avoid having to build excessively precise
geometry around joints.

Each concrete type of ContactGeometry object (whether provided as part of
Simbody or as an extension by the user) has a unique integer
ContactGeometryTypeId. Ordered pairs of these Ids are used as a key to select
a ContactTracker that is able to identify and manage contacts or potential
contacts between those two kinds of geometric objects. You can provide a default
tracking algorithm for unhandled pairs that will either ignore them or throw an
error. Note that a ContactTracker is invoked only for "narrow phase" contact;
the %ContactTrackerSubsytem handles the "broad phase" and weeds out all but a
few possible contacting surfaces that are then passed to ContactTracker for
final disposition. All the ContactTracker objects to be used with a given
%ContactTrackerSubsystem must be registered with that subsystem during
extended construction (Topology stage). Simbody provides default ContactTracker
implementations for interactions among most of its built-in ContactGeometry
types, such as Sphere-Sphere, Sphere-HalfSpace, Mesh-Sphere, Mesh-Mesh, etc.
These will be pre-registered in every %ContactTrackerSubsystem but you can
replace them with something else if you want.

@note ContactGeometryTypeIds are unique and persistent within any execution of
a program that uses them, but they are assigned at run time, possibly by
multiple asynchronous threads, and are likely to be different in different
programs and in separate runs of the same program.

The result of a ContactTracker when applied to a pair of contact surfaces, is
either a determination that the surfaces are not in contact, or a Contact object
describing their contact interaction. There are different types of these Contact
objects (for example, PointContact, LineContact, MeshContact) and the same
algorithm may result in different kinds of Contact under different
circumstances. At each evaluation, the subsystem passes in the previous Contact
object, if any, that was associated with two ContactSurface objects, then
receives an update from the algorithm. **/
//==============================================================================
//                          CONTACT TRACKER SUBSYSTEM
//==============================================================================
class SimTK_SIMBODY_EXPORT ContactTrackerSubsystem : public Subsystem {
public:
/** Create a new %ContactTrackerSubsystem and install it into the given
MultibodySystem. **/
explicit ContactTrackerSubsystem(MultibodySystem& system);

/**@name                     Find Active Contacts
Determine what contacts are occurring. **/
/**@{**/
/** Get the calculated value of the ContactSnapshot cache entry representing the
current set of Contact pairs for this system, as determined by the various
Tracker algorithms registered with this subsystem. You can call this at
Stage::Position or later; computation of contact status will be initiated if
needed. Only the past contact status and current positions are used.  This
cache entry value is precisely what will become the "previous active contacts"
state variable at the beginning of the next time step. An error will be thrown
if we have to calculate the contacts here but fail to do so; if you don't want
to deal with the possibility that an error might occur here, you can realize
contacts explicitly first (not common).
@see realizeActiveContacts() **/
const ContactSnapshot& getActiveContacts(const State& state) const;

/** (Advanced) Get an additional set of predicted Contacts that can be
anticipated from current velocity and acceleration information. You can call
this at Stage::Acceleration; computation will be initiated if needed. This
cache entry value is precisely what will become the "previous predicted
contacts" state variable at the beginning of the next time step. An error
will be thrown if we have to calculate the contacts here but fail to do so;
to avoid that you can realize them explicitly first (not common).
@see realizePredictedContacts()  **/
const ContactSnapshot& getPredictedContacts(const State& state) const;
/**@}**/


/**@name                    Contact Surface Identification
These methods map between the assigned ContactSurfaceIndex and the surfaces
that were attached to each MobilizedBody in the MultibodySystem. **/
/**@{**/

/** Get the number of surfaces being managed by this contact tracker subsystem.
These are identified by ContactSurfaceIndex values from 0 to
`getNumSurfaces()-1`. This is available after `realizeTopology()` and does not
change subsequently. You can find out which ContactSurfaceIndex got assigned
any particular ContactSurface using `getContactSurfaceIndex()`. **/
int getNumSurfaces() const;

/** Obtain the ContactSurfaceIndex that was assigned by this
%ContactTrackerSubsystem to a particular instance of a ContactSurface
on a MobilizedBody. A ContactSurface is attached to a body description using
`Body::addContactSurface()` and then that body description is used to create
an actual body instance in a MobilizedBody. The surface is identified here by
specifying the MobilizedBody of interest and the small integer ordinal that was
returned when `Body::addContactSurface()` was called. (You can provide either a
MobilizedBody or MobilizedBodyIndex here; there is an implicit conversion from
MobilizedBody to its MobilizedBodyIndex.) This method will throw an exception
if either argument is illegal or out of range.
@see Body::addContactSurface() **/
ContactSurfaceIndex getContactSurfaceIndex(MobilizedBodyIndex mobod,
                                           int contactSurfaceOrdinal) const;

/** Get the MobilizedBody associated with a particular contact surface. **/
const MobilizedBody& getMobilizedBody(ContactSurfaceIndex surfIx) const;

/** Get the ContactSurface object (detailed geometry and material properties)
that is associated with the given ContactSurfaceIndex. **/
const ContactSurface& getContactSurface(ContactSurfaceIndex surfIx) const;

/** Get the transform X_BS that gives the pose of the indicated contact
surface with respect to the body frame of the body to which it is attached. **/
const Transform& getContactSurfaceTransform(ContactSurfaceIndex surfIx) const;
/**@}**/


/**@name                     Contact Tracker management
Most users won't need to use these methods. **/
/**@{**/

/** Register the contact tracking algorithm to use for a particular pair of
ContactGeometry types, replacing the existing tracker if any. If the tracker
takes a pair (id1,id2), we will use it both for that pair and for (id2,id1) by
calling it with the arguments reversed. The subsystem takes over ownership of
the supplied heap-allocated object. **/
void adoptContactTracker(ContactTracker* tracker);

/** Return true if this subsystem has a contact tracker registered that can
deal with ineractions between surfaces using the indicated pair of geometry
types, in either order. **/
bool hasContactTracker(ContactGeometryTypeId surface1,
                       ContactGeometryTypeId surface2) const;

/** Return the contact tracker to be used for an interaction between the
indicated types of contact geometry. If the tracker requires the geometry
types to be in reverse order from the (surface1,surface2) order given here,
then the return argument \a reverseOrder will be set true, otherwise it will
be false. If no tracker was registered, this will be the default tracker. **/
const ContactTracker& getContactTracker(ContactGeometryTypeId surface1,
                                        ContactGeometryTypeId surface2,
                                        bool& reverseOrder) const;
/**@}**/

/**@name                     Advanced/Obscure
You probably don't want to call any of these methods. Some may be
unimplemented. **/
/**@{**/

/** (Advanced) Obtain the value of the ContactSnapshot state variable
representing the most recently known set of Contacts for this system. **/
const ContactSnapshot& getPreviousActiveContacts(const State& state) const;

/** (Advanced) Obtain the value of the ContactSnapshot state variable
representing the most recently predicted set of \e impending Contacts for this
system. **/
const ContactSnapshot& getPreviousPredictedContacts(const State& state) const;

/** (Advanced) Calculate the current ActiveContacts set at Position stage or
later if it hasn't already been done and return true if successful. If we can't
unambiguously determine the contact status, we'll return false and give the
caller a hint as to the latest time at which we think we could have succeeded.
If \a lastTry is true, then we throw an error on failure rather than
returning false. **/
bool realizeActiveContacts(const State& state,
                           bool         lastTry,
                           Real&        stepAdvice) const;

/** (Advanced) Calculate the set of anticipated Contacts set at Acceleration
stage if not already calculated and return true if successful. Otherwise,
problems are handled as for realizeActiveContacts(). **/
bool realizePredictedContacts(const State& state,
                              bool         lastTry,
                              Real&        stepAdvice) const;

/** (Advanced) Default constructor creates an empty Subsystem not associated
with any MultibodySystem; not very useful. **/
ContactTrackerSubsystem();
/**@}**/

SimTK_PIMPL_DOWNCAST(ContactTrackerSubsystem, Subsystem);

//--------------------------------------------------------------------------
private:
class ContactTrackerSubsystemImpl& updImpl();
const ContactTrackerSubsystemImpl& getImpl() const;
};



//==============================================================================
//                            CONTACT SNAPSHOT
//==============================================================================
/** Objects of this class represent collections of surface-pair interactions
that are being tracked at a particular instant during a simulation. These are
suitable for use as state variables for remembering past contact status and
as calculated cache entries containing the current contact status. Each
tracked surface pair has an integer ContactId that is persistent for as long
as a particular interaction is being tracked. We maintain a map providing
very fast access to individual Contact entries by ContactId. There is also
a map from ContactSurfaceIndex pairs to ContactId that can be used to see
whether we are already tracking a Contact between those surfaces; there can
be at most one Contact between a given surface pair at any given moment. **/
class SimTK_SIMBODY_EXPORT ContactSnapshot {
//TODO: replace with fast hash tables
typedef std::map<ContactId,int>     ContactMap;
// Note: we always order the key so that the first surface index is less than
// the second (they can't be equal!).
typedef std::map<std::pair<ContactSurfaceIndex,ContactSurfaceIndex>,
                 ContactId>         SurfaceMap;
public:
/** Default constructor sets timestamp to NaN. **/
ContactSnapshot() : m_time(NaN) {}

/** Restore to default-constructed condition. **/
void clear() {
    m_time = NaN;
    m_contacts.clear();
    m_id2contact.clear();
    m_surfPair2id.clear();
}

/** Set the time at which this snapshot was taken. **/
void setTimestamp(Real time) {m_time=time;}
/** At what simulation time was this contact snapshot taken? **/
Real getTimestamp() const {return m_time;}

/** Add this Contact object to this snapshot; this is a shallow,
reference-counted copy so the Contact object is not duplicated. The Contact
is assigned a slot in the array of Contact objects that can be used for
very fast access; however, that index may change if other Contact objects
are removed from the snapshot. **/
void adoptContact(Contact& contact) {
    const ContactId id = contact.getContactId();
    ContactSurfaceIndex surf1(contact.getSurface1());
    ContactSurfaceIndex surf2(contact.getSurface2());
    assert(id.isValid() && surf1.isValid() && surf2.isValid());

    // Surface pair must be ordered (low,high) for the map.
    if (surf1 > surf2) std::swap(surf1,surf2);

    assert(!hasContact(id));
    assert(!hasContact(surf1,surf2));

    const int indx = m_contacts.size();
    m_contacts.push_back(contact); // shallow copy
    m_id2contact[id] = indx;
    m_surfPair2id[std::make_pair(surf1,surf2)] = id;
}

/** Does this snapshot contain a Contact object with the given ContactId? **/
bool hasContact(ContactId id) const
{   return m_id2contact.find(id) != m_id2contact.end(); }
/** Does this snapshot contain a Contact object for the given surface pair
(in either order)? **/
bool hasContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2) const
{   if (surf1 > surf2) std::swap(surf1,surf2);
    return m_surfPair2id.find(std::make_pair(surf1,surf2))
        != m_surfPair2id.end(); }

/** Find out how many Contacts are in this snapshot. **/
int getNumContacts() const {return m_contacts.size();}
/** Get a reference to the n'th Contact in this snapshot; note that the
position of a given Contact in this array is not guaranteed to remain
unchanged when Contacts are removed from this snapshot. When in doubt,
look up the contact by ContactId instead. @see getContactById() **/
const Contact& getContact(int n) const {return m_contacts[n];}
/** If this snapshot contains a contact with the given id, return a reference
to it; otherwise, return a reference to an empty contact handle (you can check
with isEmpty()). **/
const Contact& getContactById(ContactId id) const
{   static Contact empty;
    ContactMap::const_iterator p = m_id2contact.find(id);
    return p==m_id2contact.end() ? empty : m_contacts[p->second]; }
/** If this snapshot contains a contact for the given pair of contact surfaces
(order doesn't matter), return its ContactId; otherwise, return an invalid
ContactId (you can check with isValid()). **/
ContactId getContactIdForSurfacePair(ContactSurfaceIndex surf1,
                                     ContactSurfaceIndex surf2) const
{   if (surf1 > surf2) std::swap(surf1,surf2);
    SurfaceMap::const_iterator p =
        m_surfPair2id.find(std::make_pair(surf1,surf2));
    return p==m_surfPair2id.end() ? ContactId() : p->second; }

//--------------------------------------------------------------------------
                                private:

// Remove a Contact occupying a particular slot in the Contact array. This
// will result in another Contact object being moved to occupy the now-empty
// slot to keep the array compact.
void removeContact(int n) {
    assert(0 <= n && n < m_contacts.size());
    if (n+1 == m_contacts.size()) {
        m_contacts.pop_back(); // this was the last one
        return;
    }
    // Move the last one to replace this one and update the map.
    m_contacts[n] = m_contacts.back();  // shallow copy
    m_contacts.pop_back();              // destruct
    m_id2contact[m_contacts[n].getContactId()] = n;
}


Real                m_time;         // when this snapshot was taken
Array_<Contact,int> m_contacts;     // all the contact pairs
ContactMap          m_id2contact;   // the contact pairs by contactId
SurfaceMap          m_surfPair2id;  // surfacepair -> contactId map
};

// for debugging
inline std::ostream& operator<<(std::ostream& o, const ContactSnapshot& cs) {
    o << "Contact snapshot: time=" << cs.getTimestamp()
        << " numContacts=" << cs.getNumContacts() << std::endl;
    return o;
}


} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_TRACKER_SUBSYSTEM_H_
