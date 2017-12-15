#ifndef SimTK_SIMMATH_CONTACT_IMPL_H_
#define SimTK_SIMMATH_CONTACT_IMPL_H_

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

#include "simmath/internal/common.h"
#include "simmath/internal/Contact.h"

#include <atomic>

namespace SimTK {


//==============================================================================
//                                CONTACT IMPL
//==============================================================================
/** This is the internal implementation base class for Contact. **/
class ContactImpl {
public:
    ContactImpl(ContactSurfaceIndex     surf1, 
                ContactSurfaceIndex     surf2,
                Contact::Condition      condition=Contact::Unknown) 
    :   m_referenceCount(0), m_condition(condition), 
        m_id(), m_surf1(surf1), m_surf2(surf2), m_X_S1S2() {}

    ContactImpl(ContactSurfaceIndex     surf1, 
                ContactSurfaceIndex     surf2,
                const Transform&        X_S1S2,
                Contact::Condition      condition=Contact::Unknown) 
    :   m_referenceCount(0), m_condition(condition), 
        m_id(), m_surf1(surf1), m_surf2(surf2), m_X_S1S2(X_S1S2) {}

    void setTransform(const Transform& X_S1S2) {m_X_S1S2 = X_S1S2;}
    const Transform& getTransform() const {return m_X_S1S2;}

    void setCondition(Contact::Condition cond) {m_condition=cond;}
    Contact::Condition getCondition() const {return m_condition;}

    void setContactId(ContactId id) {m_id=id;}
    ContactId getContactId() const {return m_id;}

    virtual ~ContactImpl() {
        assert(m_referenceCount == 0);
    }
    virtual ContactTypeId getTypeId() const = 0;

    /* Create a new ContactTypeId and return this unique small integer 
    (thread safe). Each distinct type of Contact should use this to
    initialize a static variable for that concrete class. */
    static ContactTypeId  createNewContactTypeId()
    {   static std::atomic<int> nextAvailableId(1);
        return ContactTypeId(nextAvailableId++); }


    /* Create a new ContactId and return this unique integer 
    (thread safe). Each distinct type of Contact should use this to
    initialize a static variable for that concrete class. This will
    roll over at approximately 1 billion. */
    static ContactId  createNewContactId()
    {   static std::atomic<int> nextAvailableId(1);
        const int MaxContactId = 999999999; // 1 billion-1
        const int id = nextAvailableId++;
        // Other threads might get a few more high-numbered ids here before
        // we reset the next available to 1, but since only one thread gets
        // exactly MaxContactId as its id, only one will execute the reset.
        if (id == MaxContactId)
            nextAvailableId = 1;
        return ContactId(id); }

protected:
friend class Contact;

    mutable int         m_referenceCount;
    Contact::Condition  m_condition;
    ContactId           m_id;
    ContactSurfaceIndex m_surf1,
                        m_surf2;
    Transform           m_X_S1S2;
};



//==============================================================================
//                          UNTRACKED CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for UntrackedContact. **/
class UntrackedContactImpl : public ContactImpl {
public:
    UntrackedContactImpl(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2) 
    :   ContactImpl(surf1, surf2, Contact::Untracked) {}

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }
};


//==============================================================================
//                           BROKEN CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for BrokenContact. **/
class BrokenContactImpl : public ContactImpl {
public:
    BrokenContactImpl
       (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
        const Transform& X_S1S2, Real separation) 
    :   ContactImpl(surf1, surf2, X_S1S2, Contact::Broken), 
        separation(separation) 
    {
    }

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class BrokenContact;
    Real        separation;
};



//==============================================================================
//                        CIRCULAR POINT CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for CircularPointContact. **/
class CircularPointContactImpl : public ContactImpl {
public:
    CircularPointContactImpl
       (ContactSurfaceIndex surf1, Real radius1, 
        ContactSurfaceIndex surf2 ,Real radius2, 
        const Transform& X_S1S2, Real radiusEff, Real depth, 
        const Vec3& origin_S1, const UnitVec3& normal_S1)
    :   ContactImpl(surf1, surf2, X_S1S2), 
        radius1(radius1), radius2(radius2), radiusEff(radiusEff), 
        depth(depth), origin_S1(origin_S1), normal_S1(normal_S1) {}

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class CircularPointContact;
    Real        radius1, radius2, radiusEff, depth;
    Vec3        origin_S1;
    UnitVec3    normal_S1;
};



//==============================================================================
//                        ELLIPTICAL POINT CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for EllipticalPointContact. **/
class EllipticalPointContactImpl : public ContactImpl {
public:
    EllipticalPointContactImpl
       (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
        const Transform& X_S1S2, const Transform& X_S1C, 
        const Vec2& k, Real depth)
    :   ContactImpl(surf1, surf2, X_S1S2), 
        X_S1C(X_S1C), k(k), depth(depth) {}

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class EllipticalPointContact;
    Transform   X_S1C;
    Vec2        k; // kmax, kmin
    Real        depth;
};



//==============================================================================
//                       BRICK HALFSPACE CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for BrickHalfSpaceContact. **/
class BrickHalfSpaceContactImpl : public ContactImpl {
public:
    BrickHalfSpaceContactImpl
       (ContactSurfaceIndex halfSpace, ContactSurfaceIndex brick, 
        const Transform& X_HB, int lowestVertex, Real depth)
    :   ContactImpl(halfSpace, brick, X_HB), 
        lowestVertex(lowestVertex), depth(depth) {}

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class BrickHalfSpaceContact;
    int     lowestVertex;
    Real    depth;
};



//==============================================================================
//                            TRIANGLE MESH IMPL
//==============================================================================
/** This is the internal implementation class for TriangleMeshContact. **/
class TriangleMeshContactImpl : public ContactImpl {
public:
    TriangleMeshContactImpl(ContactSurfaceIndex     surf1, 
                            ContactSurfaceIndex     surf2,
                            const Transform&        X_S1S2,
                            const std::set<int>&    faces1, 
                            const std::set<int>&    faces2);

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class TriangleMeshContact;

    const std::set<int> faces1;
    const std::set<int> faces2;
};



//==============================================================================
//                             POINT CONTACT IMPL
//==============================================================================
/** This is the internal implementation class for PointContact. **/
class PointContactImpl : public ContactImpl {
public:
    PointContactImpl(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
                     Vec3& location, Vec3& normal, Real radius1, Real radius2, Real depth);
    PointContactImpl(ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
                     Vec3& location, Vec3& normal, Real radius, Real depth);

    ContactTypeId getTypeId() const override {return classTypeId();}
    static ContactTypeId classTypeId() {
        static const ContactTypeId tid = createNewContactTypeId();
        return tid;
    }

private:
friend class PointContact;

    Vec3 location, normal;
    Real radius1, radius2, effectiveRadius, depth;
};



} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_IMPL_H_
