#ifndef SimTK_SIMBODY_GENERAL_CONTACT_SUBSYSTEM_H_
#define SimTK_SIMBODY_GENERAL_CONTACT_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

namespace SimTK {

class MultibodySystem;
class MobilizedBody;
class ContactGeometry;

/**
 * This class performs collision detection for use in contact modeling.  It manages sets
 * of bodies that can potentially interact with each other.  At each time step, it identifies
 * all the contacts between them.  It does not directly affect the behavior of the system
 * in any way.  Instead, it simply provides a service that can be used by other classes
 * to implement forces, events, contraints, etc. that are based on contacts between bodies.
 *
 * To use this class, first create a "contact set" by calling createContactSet().  A contact
 * set is a group of bodies which can interact with each other.  If you create multiple contact
 * sets, the bodies within each one will be checked for contacts with each other, but bodies in
 * different sets will not interact.
 *
 * Next, add bodies to the contact set.  Each one is represented by a ContactGeometry object
 * that describes its shape, the MoblizedBody it is attached to, and a Transform giving the location
 * of the geometry in the MobilizedBody's reference frame.
 *
 * Finally, call getContacts() to get a list of all contacts which exist between bodies in a
 * contact set.  Each Contact specifies two bodies that overlap, along with a description of the
 * contact point, such as its location and normal vector.
 */

class SimTK_SIMBODY_EXPORT GeneralContactSubsystem : public Subsystem {
public:
    GeneralContactSubsystem();
    explicit GeneralContactSubsystem(MultibodySystem&);
    /**
     * Create a new contact set.  The return value is an index which can be used to refer to the
     * contact set when calling other methods.
     */
    ContactSetIndex createContactSet();
    /**
     * Get the total number of contact sets that have been created.
     */
    int getNumContactSets() const;
    /**
     * Add a body to a contact set.
     *
     * @param index     the index of the contact set the body should be added to
     * @param body      the MobilizedBody whose reference frame the body is attached to
     * @param geom      a ContactGeometry describing the shape of the body
     * @param transform the location and orientation of the ContactGeometry in the MobilizedBody's reference frame
     */
    void addBody(ContactSetIndex index, const MobilizedBody& body, const ContactGeometry& geom, Transform transform);
    /**
     * Get the number of bodies in a contact set.
     */
    int getNumBodies(ContactSetIndex set) const;
    /**
     * Get the MobilizedBody in whose reference frame a body is defined.
     * 
     * @param set    the contact set the body belongs to
     * @param index  the index of the body within the contact set
     */
    const MobilizedBody& getBody(ContactSetIndex set, ContactSurfaceIndex index) const;
    /**
     * Get the ContactGeometry which defines the shape of a body.
     * 
     * @param set    the contact set the body belongs to
     * @param index  the index of the body within the contact set
     */
    const ContactGeometry& getBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index) const;
    /**
     * Get a mutable reference to the ContactGeometry which defines the shape of a body.
     * 
     * @param set    the contact set the body belongs to
     * @param index  the index of the body within the contact set
     */
    ContactGeometry& updBodyGeometry(ContactSetIndex set, ContactSurfaceIndex index);
    /**
     * Get the location and orientation of a body.
     * 
     * @param set    the contact set the body belongs to
     * @param index  the index of the body within the contact set
     */
    const Transform& getBodyTransform(ContactSetIndex set, ContactSurfaceIndex index) const;
    /**
     * Get a mutable reference to the location and orientation of a body.
     * 
     * @param set    the contact set the body belongs to
     * @param index  the index of the body within the contact set
     */
    Transform& updBodyTransform(ContactSetIndex set, ContactSurfaceIndex index);
    /**
     * Get a list of all contacts between bodies in a contact set.  Contacts are calculated at
     * Dynamics stage, so the state must have been realized to at least Dynamics stage.  This
     * subsystem is guaranteed to be realized before all ForceSubsystems, however, so a Force
     * may still invoke it to calculate forces based on contacts.
     */
    const Array_<Contact>& getContacts(const State& state, ContactSetIndex set) const;
    SimTK_PIMPL_DOWNCAST(GeneralContactSubsystem, Subsystem);
private:
    class GeneralContactSubsystemImpl& updImpl();
    const GeneralContactSubsystemImpl& getImpl() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_GENERAL_CONTACT_SUBSYSTEM_H_
