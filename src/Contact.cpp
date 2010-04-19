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


#include "simbody/internal/ContactImpl.h"
#include <set>

using namespace SimTK;
using std::set;



//==============================================================================
//                                 CONTACT
//==============================================================================

/*static*/ const char* Contact::
nameOfCondition(Condition cond) {
    switch(cond) {
    case Unknown:       return "Unknown";
    case Untracked:     return "Untracked";
    case Anticipated:   return "Untracked";
    case NewContact:    return "NewContact";
    case Ongoing:       return "Ongoing";
    case Broken:        return "Broken";
    default:            break;         
    }
    return "INVALID CONTACT CONDITION";
}

Contact::Contact(ContactImpl* impl) : impl(impl) {
    if (impl)
        impl->m_referenceCount++;
}

Contact::Contact(const Contact& src) : impl(src.impl) {
    if (impl)
        impl->m_referenceCount++;
}

void Contact::clear() {
    if (impl) {
        impl->m_referenceCount--;
        if (impl->m_referenceCount == 0)
            delete impl;
        impl = 0;
    }
}

Contact& Contact::operator=(const Contact& src) {
    clear();
    if (src.impl) {
        impl = src.impl;
        impl->m_referenceCount++;
    }
    return *this;
}

Contact::Condition  Contact::getCondition() const {return getImpl().m_condition;}
ContactId           Contact::getContactId() const {return getImpl().m_id;}
ContactSurfaceIndex Contact::getSurface1()  const {return getImpl().m_surf1;}
ContactSurfaceIndex Contact::getSurface2()  const {return getImpl().m_surf2;}
ContactTypeId       Contact::getTypeId()    const {return getImpl().getTypeId();}

Contact& Contact::setCondition(Condition condition)
{   updImpl().m_condition = condition; return *this; }
Contact& Contact::setSurfaces
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2)
{   updImpl().m_surf1 = surf1; updImpl().m_surf2 = surf2; return *this; }
Contact& Contact::setContactId(ContactId id)
{   updImpl().m_id = id; return *this; }

/*static*/ ContactId Contact::createNewContactId()
{   return ContactImpl::createNewContactId(); }



//==============================================================================
//                            UNTRACKED CONTACT
//==============================================================================
UntrackedContact::UntrackedContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2) 
:   Contact(new UntrackedContactImpl(surf1, surf2)) {}

bool UntrackedContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const UntrackedContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId UntrackedContact::classTypeId() 
{   return UntrackedContactImpl::classTypeId(); }



//==============================================================================
//                             BROKEN CONTACT
//==============================================================================
BrokenContact::BrokenContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, Real separation) 
:   Contact(new BrokenContactImpl(surf1, surf2, separation)) {}

bool BrokenContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const BrokenContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId BrokenContact::classTypeId() 
{   return BrokenContactImpl::classTypeId(); }

Real BrokenContact::getSeparation() const 
{   return getImpl().separation; }



//==============================================================================
//                          CIRCULAR POINT CONTACT
//==============================================================================
CircularPointContact::CircularPointContact
   (ContactSurfaceIndex surf1, Real radius1, 
    ContactSurfaceIndex surf2, Real radius2,
    Real radius, Real depth, const Vec3& origin, const UnitVec3& normal)
:   Contact(new CircularPointContactImpl(surf1,radius1,surf2,radius2,
                                         radius,depth,origin,normal)) {}

Real CircularPointContact::getRadius1() const
{   return getImpl().radius1; }
Real CircularPointContact::getRadius2() const
{   return getImpl().radius2; }
Real CircularPointContact::getEffectiveRadius() const
{   return getImpl().radiusEff; }
Real CircularPointContact::getDepth() const
{   return getImpl().depth; }
const Vec3& CircularPointContact::getOrigin() const
{   return getImpl().origin_G; }
const UnitVec3& CircularPointContact::getNormal() const
{   return getImpl().normal_G; }

bool CircularPointContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const CircularPointContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId CircularPointContact::classTypeId() 
{   return CircularPointContactImpl::classTypeId(); }




//==============================================================================
//                      TRIANGLE MESH CONTACT & IMPL
//==============================================================================
TriangleMeshContact::TriangleMeshContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    const std::set<int>& faces1, const std::set<int>& faces2) 
:   Contact(new TriangleMeshContactImpl(surf1, surf2, faces1, faces2)) {}

const set<int>& TriangleMeshContact::getSurface1Faces() const 
{   return getImpl().faces1; }
const set<int>& TriangleMeshContact::getSurface2Faces() const 
{   return getImpl().faces2; }

bool TriangleMeshContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const TriangleMeshContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId TriangleMeshContact::classTypeId() 
{   return TriangleMeshContactImpl::classTypeId(); }

TriangleMeshContactImpl::TriangleMeshContactImpl
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    const set<int>& faces1, const set<int>& faces2) 
:   ContactImpl(surf1, surf2), faces1(faces1), faces2(faces2) {}




//==============================================================================
//                      POINT CONTACT & IMPL (OBSOLETE)
//==============================================================================
PointContact::PointContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    Vec3& location, Vec3& normal, Real radius, Real depth) 
:   Contact(new PointContactImpl(surf1, surf2, location, normal, 
                                 radius, depth)) {}

Vec3 PointContact::getLocation() const
{   return getImpl().location; }
Vec3 PointContact::getNormal() const 
{   return getImpl().normal; }
Real PointContact::getRadius() const 
{   return getImpl().radius; }
Real PointContact::getDepth() const 
{   return getImpl().depth; }

bool PointContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const PointContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId PointContact::classTypeId() 
{   return PointContactImpl::classTypeId(); }

PointContactImpl::PointContactImpl
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    Vec3& location, Vec3& normal, Real radius, Real depth) 
:   ContactImpl(surf1, surf2), location(location), normal(normal), 
    radius(radius), depth(depth) {}


