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

/**
 * A Contact contains information about two bodies that are in contact with each other.
 * It usually is created by a CollisionDetectionAlgorithm, and is retrieved by calling getContacts()
 * on a GeneralContactSubsystem.
 *
 * This class provides a set of numbers which are only fully meaningful for point contacts.  That is,
 * it assumes that a contact can be described with a single contact location, a normal vector, a
 * contact radius, and a penetration depth.  CollisionDetectionAlgorithms which characterize contacts
 * in more complex ways may define subclasses of Contact that provide additional information.
 */

class SimTK_SIMBODY_EXPORT Contact {
public:
    /**
     * Create a contact object.
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
    Contact(int body1, int body2, Vec3& location, Vec3& normal, Real radius, Real depth);
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
private:
    ContactImpl* impl;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONTACT_H_
