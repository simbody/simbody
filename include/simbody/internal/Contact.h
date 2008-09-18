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
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

class ContactImpl;
class PointContactImpl;
class TriangleMeshContactImpl;

/**
 * A Contact contains information about two bodies that are in contact with each other.
 * It usually is created by a CollisionDetectionAlgorithm, and is retrieved by calling getContacts()
 * on a GeneralContactSubsystem.
 *
 * The base class records only the indices of the two bodies that are in contact.  CollisionDetectionAlgorithms
 * which characterize contacts in more complex ways will typically define subclasses of Contact that provide
 * additional information.
 */

class SimTK_SIMBODY_EXPORT Contact {
public:
    /**
     * Create a Contact object.
     *
     * @param body1    the index of the first body involved in the contact, specified by its index within
     *                 its contact set
     * @param body2    the index of the second body involved in the contact, specified by its index within
     *                 its contact set
     */
    Contact(int body1, int body2);
    Contact(const Contact& copy);
    Contact(ContactImpl* impl);
    ~Contact();
    Contact& operator=(const Contact& copy);
    /**
     * Get the first body involved in the contact, specified by its index within its contact set.
     */
    int getFirstBody() const;
    /**
     * Get the second body involved in the contact, specified by its index within its contact set.
     */
    int getSecondBody() const;
    const ContactImpl& getImpl() const {
        return *impl;
    }
private:
    ContactImpl* impl;
};

/**
 * This subclass of Contact represents a symmetric contact centered at a single point, such as
 * between two spheres or a sphere and a half space.  It characterizes the contact by the center location
 * and radius of the contact patch, the normal vector, and the penetration depth.
 */

class SimTK_SIMBODY_EXPORT PointContact : public Contact {
public:
    /**
     * Create a PointContact object.
     *
     * @param body1    the index of the first body involved in the contact, specified by its index within
     *                 its contact set
     * @param body2    the index of the second body involved in the contact, specified by its index within
     *                 its contact set
     * @param location the location where the two bodies touch, specified in the ground frame
     * @param normal   the surface normal at the contact location.  This is specified in the ground frame,
     *                 and points outward from body1
     * @param radius   the radius of the contact patch
     * @param depth    the penetration depth
     */
    PointContact(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth);
    /**
     * The location where the two bodies touch, specified in the ground frame.  More precisely, the
     * contact region is represented as a circular patch centered at this point and perpendicular to
     * the normal vector.
     */
    Vec3 getLocation() const;
    /**
     * Get the surface normal at the contact location.  This is specified in the ground frame,
     * and points outward from the first body.
     */
    Vec3 getNormal() const;
    /**
     * Get the radius of the contact patch.
     */
    Real getRadius() const;
    /**
     * Get the penetration depth.  This is defined as the minimum distance you would need to translate
     * one body along the normal vector to make them no longer overlap.
     */
    Real getDepth() const;
    /**
     * Determine whether a Contact object is a PointContact.
     */
    static bool isInstance(const Contact& contact);
};

/**
 * This subclass of Contact is used when one or both of the ContactGeometry objects is a TriangleMesh.
 * It stores a list of every face on each object that is partly or completely inside the other one.
 */

class SimTK_SIMBODY_EXPORT TriangleMeshContact : public Contact {
public:
    /**
     * Create a TriangleMeshContact object.
     *
     * @param body1    the index of the first body involved in the contact, specified by its index within
     *                 its contact set
     * @param body2    the index of the second body involved in the contact, specified by its index within
     *                 its contact set
     * @param faces1   the indices of all faces in the first body which are inside the second one
     * @param faces2   the indices of all faces in the second body which are inside the first one
     */
    TriangleMeshContact(int body1, int body2, const std::set<int>& faces1, const std::set<int>& faces2);
    /**
     * Get the indices of all faces of the first body that are partly or completely inside the second one.  If the first body is
     * not a TriangleMesh, this will return an empty set.
     */
    const std::set<int>& getFirstBodyFaces() const;
    /**
     * Get the indices of all faces of the second body that are partly or completely inside the first one.  If the second body is
     * not a TriangleMesh, this will return an empty set.
     */
    const std::set<int>& getSecondBodyFaces() const;
    /**
     * Determine whether a Contact object is a TriangleMeshContact.
     */
    static bool isInstance(const Contact& contact);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_H_
