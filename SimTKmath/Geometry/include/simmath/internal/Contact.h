#ifndef SimTK_SIMMATH_CONTACT_H_
#define SimTK_SIMMATH_CONTACT_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <set>

namespace SimTK {


/** @class SimTK::ContactSurfaceIndex
This defines a unique index for all the contact surfaces being handled
either by a ContactTrackerSubsystem or within a single ContactSet of a
GeneralContactSubsystem. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactSurfaceIndex);

/** @class SimTK::ContactId
This is a unique integer Id assigned to each contact pair when we first
begin to track it. The Id persists for as long as we are tracking this pair
of surfaces; it is the basis for our notions of "the same contact" and
"continuing contact". After we stop tracking a contact pair, its Id will not 
refer to any contact pair and any given Id will not be reused for a very long 
time; hence, these will typically be large numbers even if there are only a 
small number of contacts at any given moment. The Id is unique only within an 
instance of ContactTrackerSubsystem. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactId);

/** @class SimTK::ContactTypeId
This is a small integer that serves as the unique typeid for each type
of concrete Contact class. This is used to select an appropriate contact
response method for a given type of Contact. This Id is shared by all
threads of a given program execution but you can't expect to get the same Id
for different programs or different executions of the same program. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ContactTypeId);


class ContactImpl;
class UntrackedContactImpl;
class BrokenContactImpl;
class CircularPointContactImpl;
class EllipticalPointContactImpl;
class BrickHalfSpaceContactImpl;
class TriangleMeshContactImpl;

class PointContactImpl; // deprecated


//==============================================================================
//                                CONTACT
//==============================================================================
/** A Contact contains information about the spatial relationship between two 
surfaces that are near, or in contact with, each other. It usually is created 
by a ContactTracker or CollisionDetectionAlgorithm, and is retrieved from a
ContactTrackerSubsystem or GeneralContactSubsystem.

The base class records only the indices of the two surfaces that are in 
contact, and the relative Transform between them at the time the Contact
was recorded. ContactTrackers or CollisionDetectionAlgorithms which 
characterize contacts in more complex ways will return objects that are 
subclasses of Contact that provide additional information. **/
class SimTK_SIMMATH_EXPORT Contact {
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

    /** Get the persistent ContactId that has been assigned to this Contact
    object if there is one (otherwise this will be invalid -- you can check
    with isValid(). **/
    ContactId getContactId() const;
    /** Find out the current condition of this Contact object. **/
    Condition getCondition() const;
    /** Get the first surface involved in the contact, specified by 
    its index within its contact set or ContactTrackerSubsystem. **/
    ContactSurfaceIndex getSurface1() const;
    /** Get the second surface involved in the contact, specified by 
    its index within its contact set or ContactTrackerSubsystem. **/
    ContactSurfaceIndex getSurface2() const;
    /** Return the transform X_S1S2 giving the pose of surface 2's frame 
    measured and expressed in surface 1's frame, recorded at the time this 
    Contact object was calculated. **/
    const Transform& getTransform() const;

    /** Set the ContactId for this Contact object. This must persist over the
    lifetime of a single contact event. **/
    Contact& setContactId(ContactId id);
    /** Set the current Condition. **/
    Contact& setCondition(Condition condition);
    /** Set the surfaces tracked by this Contact object. **/
    Contact& setSurfaces(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2);
    /** Set the surface-to-surface relative transform X_S1S2. **/
    Contact& setTransform(const Transform& X_S1S2);

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
class SimTK_SIMMATH_EXPORT UntrackedContact : public Contact {
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
class SimTK_SIMMATH_EXPORT BrokenContact : public Contact {
public:
    /** Create a BrokenContact object.
    @param surf1        The index of the first surface involved in the contact.
    @param surf2        The index of the second surface involved in the contact.
    @param X_S1S2       The surface-to-surface relative transform.
    @param separation   The minimum distance between the surfaces, with 
                        separation > cutoff >= 0 always. **/
    BrokenContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
                  const Transform& X_S1S2, Real separation); 

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
class SimTK_SIMMATH_EXPORT CircularPointContact : public Contact {
public:
    /** Create a CircularPointContact object.
    @param surf1        the index of the first surface involved in the contact 
    @param radius1      surf1's uniform radius at the contact initiation point
    @param surf2        the index of the second surface involved in the contact 
    @param radius2      surf2's uniform radius at the contact initiation point
    @param X_S1S2       the surface-to-surface relative transform
    @param radius       the effective combined radius to use
    @param depth        the penetration depth d (>0) or separation distance 
                        (<0); surfaces are at +/- d/2 from the origin, up and 
                        down the normal 
    @param origin_S1    origin point for the contact patch frame, in S1
    @param normal_S1    the common normal at onset, pointing from surface1 to
                        surface2, expressed in S1. This is the z axis of the 
                        patch frame. **/
    CircularPointContact
       (ContactSurfaceIndex surf1, Real radius1, 
        ContactSurfaceIndex surf2, Real radius2, 
        const Transform& X_S1S2, Real radius, Real depth, 
        const Vec3& origin_S1, const UnitVec3& normal_S1);

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
    /** Get the origin OP of the contact patch frame P, in S1. **/
    const Vec3& getOrigin() const;
    /** Get the z axis of the contact patch frame, which is the common surface 
    normal at the initial contact point, pointing outward from surface1 towards
    surface2 at initial contact. This is a unit vector expressed in S1. **/
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
//                           ELLIPTICAL POINT CONTACT
//==============================================================================
/** This subclass of Contact represents a contact between two non-conforming
surfaces 1 and 2 that initially meet at a point and where each surface has two
principal curvatures (maximum and minimum) in perpendicular directions. The 
prototypical example is ellipsoid-ellipsoid contact, but this includes a wide
range of contacts between smooth surfaces, where the surfaces have two 
continuous spatial derivatives at the contact point. The contact plane on
each surface is parameterized by its principal curvature directions x,y with the
surface's contact point at the origin, and the common normal as the z axis. The
surface is thus approximated by a paraboloid z=Ax^2+By^2 for which A=kx/2, 
B=ky/2 where kx,ky are the curvatures in the x,y directions. Here A>=0, A>=B, 
but B can be negative indicating a hyperbolic paraboloid (saddle). Each surface
is parameterized separately: the z axes are along the same line, but the x,y 
axes are relatively rotated about z by an angle theta, with 
cos(theta)=dot(x1,x2)=dot(y1,y2).

The surface of relative separation of the two surfaces will also be a 
paraboloid, sharing the z axis with the contact surfaces but having its own
relative principal curvatures and directions. The undeformed contact is 
ultimately characterized by this relative paraboloid and a penetration depth
d (d>0 when surfaces overlap, <0 when separated), in a contact frame where 
x,y are the relative principal curvature directions, z is the common normal 
oriented to point away from surface1, and the origin OP is centered such that 
the surface2 contact point is at O-(d/2)z and surface1 contact point is at
O+(d/2)z. 

<h3>References</h3>
    - Johnson, K.L. "Contact Mechanics", Cambridge University Press, 1985
      (corrected ed. 1987), sec. 4.1, pp. 84-88. 
    - Goldsmith, W. "Impact", Dover, 2001, sec. 4.2, pp. 83-85.
**/
class SimTK_SIMMATH_EXPORT EllipticalPointContact : public Contact {
public:
    /** Create a EllipticalPointContact object.
    @param surf1        the index of the first surface involved in the contact 
    @param surf2        the index of the second surface involved in the contact 
    @param X_S1S2       the surface-to-surface relative transform
    @param X_S1C        contact paraboloid coordinate frame C in S1 frame; x is
                        kmax direction, y is kmin direction, z points away from
                        surf1; origin OC is at midpoint between contact 
                        points on the two surfaces
    @param k            maximum and minimum curvatures kmax,kmin of the 
                        relative contact paraboloid
    @param depth        penetration depth d(>0) or separation (<0); surf1
                        contact pt at OC+(d/2)z, surf2 contact pt at OC-(d/2)z
    **/
    EllipticalPointContact
       (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
        const Transform& X_S1S2, 
        const Transform& X_S1C, const Vec2& k, Real depth); 

    /** Get the relative curvatures at the contact point, ordered kmax,kmin
    with kmax >= kmin. Note that it is possible that kmin < 0. **/
    const Vec2& getCurvatures() const;
    /** Get the frame C in which the contact paraboloid is expressed, as the
    transform X_S1C. The Cx axis is the direction of maximum relative
    curvature kmax, Cy is the direction of minimum curvature kmin, and Cz
    is the contact normal direct away from S1's surface. The origin OC is
    a point centered between the contact points on the two surfaces; those
    points are at +/- depth/2 along Cz away from OC. **/
    const Transform& getContactFrame() const;
    /** Get the penetration depth (>0) or separation distance (<0), also known 
    as the "approach". This is defined as the minimum distance you would need to
    translate surface2 along the normal vector to make the surfaces just touch
    at their contact points without overlap. **/
    Real getDepth() const;

    /** Determine whether a Contact object is an EllipticalPointContact. **/
    static bool isInstance(const Contact& contact);
    static const EllipticalPointContact& getAs(const Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<const EllipticalPointContact&>(contact); }
    static EllipticalPointContact& updAs(Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<EllipticalPointContact&>(contact); }

    /** Get the unique small-integer id for the CircularPointContact class. **/
    static ContactTypeId classTypeId();

private:
    const EllipticalPointContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const EllipticalPointContactImpl&>
                    (Contact::getImpl()); }
};



//==============================================================================
//                           BRICK HALFSPACE CONTACT
//==============================================================================
/** This subclass of Contact is used when one ContactGeometry object is a
half plane and the other is a Brick. This is a warmup for general convex
mesh contact. **/
class SimTK_SIMMATH_EXPORT BrickHalfSpaceContact : public Contact {
public:
    /** Create a BrickHalfSpaceContact object.
    @param halfSpace    the surface index of the halfspace
    @param brick        the surface index of the brick
    @param X_HB         the transform giving the brick's frame measured and 
                        expressed in the halfspace's frame
    @param lowestVertex which vertex of the brick is closest to (if separated)
                        or furthest in (if penetrating) the halfspace
    @param depth        the penetration depth (if depth>0) or separation
                        distance between the lowestVertex and halfspace surface
    **/
    BrickHalfSpaceContact(ContactSurfaceIndex     halfSpace, 
                          ContactSurfaceIndex     brick,
                          const Transform&        X_HB,
                          int                     lowestVertex,
                          Real                    depth);

    /** Get the vertex index (0-7) of the brick's vertex that is closest to or
    most penetrated into the halfspace. **/
    int getLowestVertex() const;

    /** Get the penetration depth (>0) or separation distance (<0) from the
    brick's lowest vertex to the halfspace surface. **/
    Real getDepth() const;

    /** Determine whether a Contact object is a BrickHalfSpaceContact. **/
    static bool isInstance(const Contact& contact);
    
    /** Recast a brick-halfspace contact given as a generic Contact object to a 
    const reference to a concrete BrickHalfSpaceContact object. **/
    static const BrickHalfSpaceContact& getAs(const Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<const BrickHalfSpaceContact&>(contact); }
        
    /** Recast a brick-halfspace contact given as a generic Contact object to a 
    writable reference to a concrete BrickHalfSpaceContact object. **/
    static BrickHalfSpaceContact& updAs(Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<BrickHalfSpaceContact&>(contact); }

    /** Obtain the unique small-integer id for the BrickHalfSpaceContact 
    class. **/
    static ContactTypeId classTypeId();

private:
    const BrickHalfSpaceContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const BrickHalfSpaceContactImpl&>
                    (Contact::getImpl()); }
};



//==============================================================================
//                           TRIANGLE MESH CONTACT
//==============================================================================
/** This subclass of Contact is used when one or both of the ContactGeometry 
objects is a TriangleMesh. It stores a list of every face on each object 
that is partly or completely inside the other one. **/
class SimTK_SIMMATH_EXPORT TriangleMeshContact : public Contact {
public:
    /** Create a TriangleMeshContact object.
    @param surf1    the index of the first surface involved in the contact, 
                    specified by its index within its contact set
    @param surf2    the index of the second surface involved in the contact, 
                    specified by its index within its contact set
    @param X_S1S2   the transform giving surf2's frame measured and expressed
                    in surf1's frame
    @param faces1   the indices of all faces in the first surface which are 
                    inside the second one
    @param faces2   the indices of all faces in the second surface which are
                    inside the first one **/
    TriangleMeshContact(ContactSurfaceIndex     surf1, 
                        ContactSurfaceIndex     surf2,
                        const Transform&        X_S1S2,
                        const std::set<int>&    faces1, 
                        const std::set<int>&    faces2);

    /** Get the indices of all faces of surface1 that are partly or completely 
    inside surface2. If surface1 is not a TriangleMesh, this will return an 
    empty set. **/
    const std::set<int>& getSurface1Faces() const;
    /** Get the indices of all faces of surface2 that are partly or completely
    inside surface1. If surface2 is not a TriangleMesh, this will return an 
    empty set. **/
    const std::set<int>& getSurface2Faces() const;

    /** Determine whether a Contact object is a TriangleMeshContact. **/
    static bool isInstance(const Contact& contact);
    /** Recast a triangle mesh given as a generic Contact object to a 
    const reference to a concrete TriangleMeshContact object. **/
    static const TriangleMeshContact& getAs(const Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<const TriangleMeshContact&>(contact); }
    /** Recast a triangle mesh given as a generic Contact object to a 
    writable reference to a concrete TriangleMeshContact object. **/
    static TriangleMeshContact& updAs(Contact& contact)
    {   assert(isInstance(contact)); 
        return static_cast<TriangleMeshContact&>(contact); }

    /** Obtain the unique small-integer id for the TriangleMeshContact 
    class. **/
    static ContactTypeId classTypeId();

private:
    const TriangleMeshContactImpl& getImpl() const 
    {   assert(isInstance(*this)); 
        return reinterpret_cast<const TriangleMeshContactImpl&>
                    (Contact::getImpl()); }
};




//==============================================================================
//                               POINT CONTACT
//==============================================================================
/** OBSOLETE -- use CircularPointContact or EllipticalPointContact.
 * This subclass of Contact represents a symmetric contact centered at a single
 * point, such as between two spheres or a sphere and a half space. It 
 * characterizes the contact by the center location and radius of the contact 
 * patch, the normal vector, and the penetration depth.
 */
class SimTK_SIMMATH_EXPORT PointContact : public Contact {
public:
    /**
     * Create a PointContact object representing a general (elliptical) contact.
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
     * @param radius1  the first principal relative radius of curvature of the contact surface
     * @param radius2  the second principal relative radius of curvature of the contact surface
     * @param depth    the penetration depth
     */
    PointContact(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
                 Vec3& location, Vec3& normal, Real radius1, Real radius2, Real depth);
    /**
     * Create a PointContact object representing a circularly symmetric contact.
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
     * @param radius   the relative radius of curvature of the contact surface
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
     * Get the first principal relative radius of curvature of the contact surface.
     */
    Real getRadiusOfCurvature1() const;
    /**
     * Get the second principal relative radius of curvature of the contact surface.
     */
    Real getRadiusOfCurvature2() const;
    /**
     * Get the effective relative radius of curvature of the contact surface.  This is equal to
     * sqrt(R1*R2), where R1 and R2 are the principal relative radii of curvature.
     */
    Real getEffectiveRadiusOfCurvature() const;
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

#endif // SimTK_SIMMATH_CONTACT_H_
