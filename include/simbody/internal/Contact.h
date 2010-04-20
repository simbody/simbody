#ifndef SimTK_SIMBODY_CONTACT_H_
#define SimTK_SIMBODY_CONTACT_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-10 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"

namespace SimTK {


/** This defines a unique index for all the contact surfaces being handled
either by a ContactTrackerSubsystem or within a single ContactSet of a
GeneralContactSubsystem. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactSurfaceIndex);

/** This is a unique integer Id assigned to each contact pair when we first
begin to track it. The Id persists for as long as we are tracking this pair
of surfaces; it is the basis for our notions of "the same contact" and
"continuing contact". After we stop tracking a contact pair, its Id will not 
refer to any contact pair and any given Id will not be reused for a very long 
time; hence, these will typically be large numbers even if there are only a 
small number of contacts at any given moment. The Id is unique only within an 
instance of ContactTrackerSubsystem. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactId);

/** This is a small integer that serves as the unique typeid for each type
of concrete Contact class. This is used to select an appropriate contact
response method for a given type of Contact. This Id is unique across all
threads of a given program execution but you can't expect to get the same Id
for different programs or different executions of the same program. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactTypeId);


class ContactImpl;
class UntrackedContactImpl;
class BrokenContactImpl;
class CircularPointContactImpl;
class TriangleMeshContactImpl;
class PointContactImpl;


//==============================================================================
//                                CONTACT
//==============================================================================
/** A Contact contains information about two surfaces that are in contact with
each other. It usually is created by a CollisionDetectionAlgorithm, and is 
retrieved by calling getContacts() on a GeneralContactSubsystem.

The base class records only the indices of the two surfaces that are in 
contact. CollisionDetectionAlgorithms which characterize contacts in more 
complex ways will typically define subclasses of Contact that provide
additional information. **/
class SimTK_SIMBODY_EXPORT Contact {
public:
    /** The Contact::Condition tracks the status of a Contact through its
    lifetime. **/
    enum Condition {
        Unknown,    ///< this is an illegal value
        Untracked,  ///< this pair not yet being tracked; might not contact
        Anticipated,///< we expect these to contact soon
        NewContact, ///< first time seen; needs a ContactId assigned
        Ongoing,    ///< was new or ongoing before; still in contact now
        Broken      ///< was new or ongoing before; no longer in contact
    };
    /** Returns a human-readable name corresponding to the given Condition; 
    useful for debugging. If the Condition is unrecognized the method will 
    return some text to that effect rather than crashing. **/
    static const char* nameOfCondition(Condition);

    /** The default constructor creates an empty handle. **/
    Contact() : impl(0) {}
    /** Copy constructor is shallow and reference-counted; this handle will 
    point to the same object as does the \a source. **/
    Contact(const Contact& source);
    /** Destructor clears the handle, deleting the referenced object if this
    was the last reference. **/
    ~Contact() {clear();}
    /** Copy assignment is shallow and reference-counted; this handle will 
    point to the same object as does the \a source. **/
    Contact& operator=(const Contact& source);
    /** Clear this handle, deleting the referenced object if this
    was the last reference.  **/
    void clear();
    /** See if this handle is empty. **/
    bool isEmpty() const {return impl==0;}

    /** Find out the current condition of this Contact object. **/
    Condition getCondition() const;
    /** Get the persistent ContactId that has been assigned to this Contact
    object if there is one (otherwise this will be invalid -- you can check
    with isValid(). **/
    ContactId getContactId() const;
    /** Get the first surface involved in the contact, specified by 
    its index within its contact set or ContactTrackerSubsystem. **/
    ContactSurfaceIndex getSurface1() const;
    /** Get the second surface involved in the contact, specified by 
    its index within its contact set or ContactTrackerSubsystem. **/
    ContactSurfaceIndex getSurface2() const;

    /** Set the current Condition. **/
    Contact& setCondition(Condition condition);
    /** Set the surfaces tracked by this Contact object. **/
    Contact& setSurfaces(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2);
    /** Set the ContactId for this Contact object. This must persist over the
    lifetime of a single contact event. **/
    Contact& setContactId(ContactId id);

    /** Return a unique small integer corresponding to the concrete type
    of Contact object being referenced by this handle. **/
    ContactTypeId getTypeId() const;

    /** This creates a new ContactId starting from 1 and increasing for a very
    long time (to a billion or so) before repeating. ContactId 0 is never
    returned and this call is thread-safe. **/
    static ContactId createNewContactId();

    const ContactImpl& getImpl() const {assert(impl); return *impl;}
    ContactImpl&       updImpl()       {assert(impl); return *impl;}
protected:
    explicit Contact(ContactImpl* impl);
private:
    ContactImpl* impl;
};

inline std::ostream& operator<<(std::ostream& o, const Contact& c) {
    o << "Contact id=" << c.getContactId() 
                       << " (typeId=" << c.getTypeId() << "):\n";
    o << "  surf1,surf2=" << c.getSurface1() << ","
                          << c.getSurface2() << "\n";
    o << "  condition=" << Contact::nameOfCondition(c.getCondition()) << "\n";
    return o;
}



//==============================================================================
//                            UNTRACKED CONTACT
//==============================================================================
/** This subclass of Contact represents a pair of contact surfaces that are
not yet being tracked; there is no ContactId for them. We don't yet know what 
kind of Contact object would be appropriate for them, so this is a placeholder
until a ContactTracker replaces it with something meaningful. The contact
condition for one of these is always "Untracked". **/
class SimTK_SIMBODY_EXPORT UntrackedContact : public Contact {
public:
    /** Default constructor creates an empty handle. **/
    UntrackedContact() {}
    /** Create an UntrackedContact object.
    @param surf1    the index of the first surface involved in the contact, 
                    specified by its index within its contact set
    @param surf2    the index of the second surface involved in the contact, 
                    specified by its index within its contact set **/
    UntrackedContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2); 

    /** Determine whether a Contact object is an UntrackedContact. **/
    static bool isInstance(const Contact& contact);
    /** Obtain the unique small-integer id for the UntrackedContact class. **/ 
    static ContactTypeId classTypeId();

private:
    const UntrackedContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const UntrackedContactImpl&>
                    (Contact::getImpl()); }
};



//==============================================================================
//                            BROKEN CONTACT
//==============================================================================
/** This subclass of Contact represents a pair of contact surfaces that were
in contact (meaning within cutoff range) but have now gone out of range. This
is the last time we will use this ContactId. The only parameters here are the
surfaces and the separation distance (> cutoff). If someone cares, the 
separation distance can be used to estimate the time at which contact was
broken. **/
class SimTK_SIMBODY_EXPORT BrokenContact : public Contact {
public:
    /** Create a BrokenContact object.
    @param surf1        The index of the first surface involved in the contact.
    @param surf2        The index of the second surface involved in the contact. 
    @param separation   The minimum distance between the surfaces, with 
                        separation > cutoff >= 0 always. **/
    BrokenContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
                  Real separation); 

    /** Get the separation (> cutoff >= 0) between the two surfaces at the time
    we decided the contact had been broken. Note that the sign convention is
    opposite from \e depth which is negative when separated. **/
    Real getSeparation() const;

    /** Determine whether a Contact object is a BrokenContact. **/
    static bool isInstance(const Contact& contact);
    /** Obtain the unique small-integer id for the BrokenContact class. **/ 
    static ContactTypeId classTypeId();

private:
    const BrokenContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const BrokenContactImpl&>(Contact::getImpl()); }
};



//==============================================================================
//                           CIRCULAR POINT CONTACT
//==============================================================================
/** This subclass of Contact represents a contact between two non-conforming
surfaces 1 and 2 that initially meet at a point where each surface has a 
uniform radius of curvature in all directions (R1 and R2), like a sphere (inside 
or outside) or a halfspace, resulting in a contact region with circular 
symmetry. One of the objects may be concave, with negative curvature, as long
as the two objects are not conforming.

The undeformed geometry is characterized here by the effective radius 
R=1/(1/R1+1/R2), a normal vector z defining the penetration direction, a scalar
penetration depth d (d>0 when surfaces overlap), and a patch origin point OP
located centered on the patch normal such that each undeformed surface is d/2
up or down the normal from OP. **/
class SimTK_SIMBODY_EXPORT CircularPointContact : public Contact {
public:
    /** Create a CircularPointContact object.
    @param surf1    the index of the first surface involved in the contact 
    @param radius1  surf1's uniform radius at the contact initiation point
    @param surf2    the index of the second surface involved in the contact 
    @param radius2  surf2's uniform radius at the contact initiation point
    @param radius   the effective combined radius to use
    @param depth    the penetration depth d (>0) or separation distance (<0); 
                    surfaces are at +/- d/2 from the origin, up and down the 
                    normal 
    @param origin   origin point for the contact patch frame, in G
    @param normal   the common normal at onset, pointing from surface1 to
                    surface2, expressed in G. This is the z axis of the patch 
                    frame. **/
    CircularPointContact
       (ContactSurfaceIndex surf1, Real radius1, 
        ContactSurfaceIndex surf2, Real radius2,
        Real radius, Real depth, const Vec3& origin, const UnitVec3& normal);

    /** Get the radius of surface1 at the contact point. **/
    Real getRadius1() const;
    /** Get the radius of surface2 at the contact point. **/
    Real getRadius2() const;
    /** Get the precalculated effective radius R at the contact point, where 
    R=1/(1/R1+1/R2). **/
    Real getEffectiveRadius() const;
    /** Get the penetration depth (>0) or separation distance (<0), also known 
    as the "approach". This is defined as the minimum distance you would need to
    translate surface2 along the normal vector to make the surfaces just touch
    at their contact points without overlap. **/
    Real getDepth() const;
    /** Get the origin OP of the contact patch frame P, in G. **/
    const Vec3& getOrigin() const;
    /** Get the z axis of the contact patch frame, which is the common surface 
    normal at the initial contact point, pointing outward from surface1 towards
    surface2 at initial contact. This is a unit vector expressed in G. **/
    const UnitVec3& getNormal() const;

    /** Determine whether a Contact object is a CircularPointContact. **/
    static bool isInstance(const Contact& contact);
    static const CircularPointContact& getAs(const Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<const CircularPointContact&>(contact); }
    static CircularPointContact& updAs(Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<CircularPointContact&>(contact); }

    /** Get the unique small-integer id for the CircularPointContact class. **/
    static ContactTypeId classTypeId();

private:
    const CircularPointContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const CircularPointContactImpl&>
                    (Contact::getImpl()); }
};



//==============================================================================
//                           TRIANGLE MESH CONTACT
//==============================================================================
/**
 * This subclass of Contact is used when one or both of the ContactGeometry 
 * objects is a TriangleMesh. It stores a list of every face on each object 
 * that is partly or completely inside the other one.
 */
class SimTK_SIMBODY_EXPORT TriangleMeshContact : public Contact {
public:
    /**
     * Create a TriangleMeshContact object.
     *
     * @param surf1    the index of the first surface involved in the contact, 
     *                 specified by its index within its contact set
     * @param surf2    the index of the second surface involved in the contact, 
     *                 specified by its index within its contact set
     * @param faces1   the indices of all faces in the first surface which are 
     *                 inside the second one
     * @param faces2   the indices of all faces in the second surface which are
     *                 inside the first one
     */
    TriangleMeshContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
                        const std::set<int>& faces1, 
                        const std::set<int>& faces2);
    /**
     * Get the indices of all faces of surface1 that are partly
     * or completely inside surface2. If surface1 
     * is not a TriangleMesh, this will return an empty set.
     */
    const std::set<int>& getSurface1Faces() const;
    /**
     * Get the indices of all faces of surface2 that are 
     * partly or completely inside surface1. If surface2
     * is not a TriangleMesh, this will return an empty set.
     */
    const std::set<int>& getSurface2Faces() const;
    /**
     * Determine whether a Contact object is a TriangleMeshContact.
     */
    static bool isInstance(const Contact& contact);
    /** 
     * Obtain the unique small-integer id for the TriangleMeshContact class. 
     */
    static ContactTypeId classTypeId();

private:
    const TriangleMeshContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const TriangleMeshContactImpl&>
                    (Contact::getImpl()); }
};




//==============================================================================
//                           POINT CONTACT (OBSOLETE)
//==============================================================================
// OBSOLETE: this incorrectly includes the patch radius as part of the 
// geometric information, but the returned value assumes a Hertz elastic
// contact model. Use CircularPointContact instead.
/**
 * This subclass of Contact represents a symmetric contact centered at a single
 * point, such as between two spheres or a sphere and a half space. It 
 * characterizes the contact by the center location and radius of the contact 
 * patch, the normal vector, and the penetration depth.
 */
class SimTK_SIMBODY_EXPORT PointContact : public Contact {
public:
    /**
     * Create a PointContact object.
     *
     * @param surf1    the index of the first surface involved in the contact, 
     *                 specified by its index within its contact set
     * @param surf2    the index of the second surface involved in the contact, 
     *                 specified by its index within its contact set
     * @param location the location where the two surfaces touch, specified in 
     *                 the ground frame
     * @param normal   the surface normal at the contact location. This is 
     *                 specified in the ground frame, and points outward 
     *                 from surface1 towards surface2
     * @param radius   the radius of the contact patch
     * @param depth    the penetration depth
     */
    PointContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
                 Vec3& location, Vec3& normal, Real radius, Real depth);
    /**
     * The location where the two surfaces touch, specified in the ground frame.
     * More precisely, the contact region is represented as a circular patch 
     * centered at this point and perpendicular to the normal vector.
     */
    Vec3 getLocation() const;
    /**
     * Get the surface normal at the contact location. This is specified in the 
     * ground frame, and points outward from surface1 towards surface2.
     */
    Vec3 getNormal() const;
    /**
     * Get the radius of the contact patch.
     */
    Real getRadius() const;
    /**
     * Get the penetration depth. This is defined as the minimum distance you 
     * would need to translate one surface along the normal vector to make the
     * surfaces no longer overlap.
     */
    Real getDepth() const;
    /**
     * Determine whether a Contact object is a PointContact.
     */
    static bool isInstance(const Contact& contact);
    /** 
     * Obtain the unique small-integer id for the PointContact class. 
     */
    static ContactTypeId classTypeId();

private:
    const PointContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const PointContactImpl&>(Contact::getImpl()); }
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_H_
