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

ContactId           Contact::getContactId() const {return getImpl().m_id;}
Contact::Condition  Contact::getCondition() const {return getImpl().m_condition;}
ContactSurfaceIndex Contact::getSurface1()  const {return getImpl().m_surf1;}
ContactSurfaceIndex Contact::getSurface2()  const {return getImpl().m_surf2;}
const Transform&    Contact::getTransform() const {return getImpl().m_X_S1S2;}
ContactTypeId       Contact::getTypeId()    const {return getImpl().getTypeId();}

Contact& Contact::setContactId(ContactId id)
{   updImpl().m_id = id; return *this; }
Contact& Contact::setCondition(Condition condition)
{   updImpl().m_condition = condition; return *this; }
Contact& Contact::setSurfaces
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2)
{   updImpl().m_surf1 = surf1; updImpl().m_surf2 = surf2; return *this; }
Contact& Contact::setTransform(const Transform& X_S1S2)
{   updImpl().m_X_S1S2 = X_S1S2; return *this; }

/*static*/ ContactId Contact::createNewContactId()
{   return ContactImpl::createNewContactId(); }



//==============================================================================
//                            UNTRACKED CONTACT
//==============================================================================
UntrackedContact::UntrackedContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2) 
:   Contact(new UntrackedContactImpl(surf1, surf2)) {}

/*static*/ bool UntrackedContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const UntrackedContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId UntrackedContact::classTypeId() 
{   return UntrackedContactImpl::classTypeId(); }



//==============================================================================
//                             BROKEN CONTACT
//==============================================================================
BrokenContact::BrokenContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    const Transform& X_S1S2, Real separation) 
:   Contact(new BrokenContactImpl(surf1, surf2, X_S1S2, separation)) {}

/*static*/ bool BrokenContact::isInstance(const Contact& contact) {
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
    const Transform& X_S1S2, Real radius, Real depth, 
    const Vec3& origin_S1, const UnitVec3& normal_S1)
:   Contact(new CircularPointContactImpl(surf1,radius1,surf2,radius2,X_S1S2,
                                         radius,depth,origin_S1,normal_S1)) {}

Real CircularPointContact::getRadius1() const
{   return getImpl().radius1; }
Real CircularPointContact::getRadius2() const
{   return getImpl().radius2; }
Real CircularPointContact::getEffectiveRadius() const
{   return getImpl().radiusEff; }
Real CircularPointContact::getDepth() const
{   return getImpl().depth; }
const Vec3& CircularPointContact::getOrigin() const
{   return getImpl().origin_S1; }
const UnitVec3& CircularPointContact::getNormal() const
{   return getImpl().normal_S1; }

bool CircularPointContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const CircularPointContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId CircularPointContact::classTypeId() 
{   return CircularPointContactImpl::classTypeId(); }



//==============================================================================
//                          ELLIPTICAL POINT CONTACT
//==============================================================================
EllipticalPointContact::EllipticalPointContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    const Transform& X_S1S2, const Transform& X_S1C,
    const Vec2& k, Real depth)
:   Contact(new EllipticalPointContactImpl(surf1,surf2,X_S1S2,X_S1C,
                                           k,depth)) {}

const Vec2& EllipticalPointContact::getCurvatures() const
{   return getImpl().k; }
const Transform& EllipticalPointContact::getContactFrame() const
{   return getImpl().X_S1C; }
Real EllipticalPointContact::getDepth() const
{   return getImpl().depth; }

bool EllipticalPointContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const EllipticalPointContactImpl*>
        (&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId EllipticalPointContact::classTypeId() 
{   return EllipticalPointContactImpl::classTypeId(); }



//==============================================================================
//                      TRIANGLE MESH CONTACT & IMPL
//==============================================================================
TriangleMeshContact::TriangleMeshContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
    const Transform& X_S1S2,
    const std::set<int>& faces1, const std::set<int>& faces2) 
:   Contact(new TriangleMeshContactImpl(surf1, surf2, X_S1S2, 
                                        faces1, faces2)) {}

const set<int>& TriangleMeshContact::getSurface1Faces() const 
{   return getImpl().faces1; }
const set<int>& TriangleMeshContact::getSurface2Faces() const 
{   return getImpl().faces2; }

/*static*/ bool TriangleMeshContact::isInstance(const Contact& contact) 
{   return (dynamic_cast<const TriangleMeshContactImpl*>(&contact.getImpl())
            != 0); }

/*static*/ ContactTypeId TriangleMeshContact::classTypeId() 
{   return TriangleMeshContactImpl::classTypeId(); }

TriangleMeshContactImpl::TriangleMeshContactImpl
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
    const Transform& X_S1S2,
    const set<int>& faces1, const set<int>& faces2) 
:   ContactImpl(surf1, surf2, X_S1S2), faces1(faces1), faces2(faces2) {}




//==============================================================================
//                      POINT CONTACT & IMPL (OBSOLETE)
//==============================================================================
PointContact::PointContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    Vec3& location, Vec3& normal, Real radius1, Real radius2, Real depth)
:   Contact(new PointContactImpl(surf1, surf2, location, normal, radius1, radius2, depth)) {}
PointContact::PointContact
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
    Vec3& location, Vec3& normal, Real radius, Real depth)
:   Contact(new PointContactImpl(surf1, surf2, location, normal, radius, depth)) {}

Vec3 PointContact::getLocation() const
{   return getImpl().location; }
Vec3 PointContact::getNormal() const 
{   return getImpl().normal; }
Real PointContact::getRadiusOfCurvature1() const
{   return getImpl().radius1; }
Real PointContact::getRadiusOfCurvature2() const
{   return getImpl().radius2; }
Real PointContact::getEffectiveRadiusOfCurvature() const
{   return getImpl().effectiveRadius; }
Real PointContact::getDepth() const 
{   return getImpl().depth; }

/*static*/ bool PointContact::isInstance(const Contact& contact) {
    return (dynamic_cast<const PointContactImpl*>(&contact.getImpl()) != 0);
}

/*static*/ ContactTypeId PointContact::classTypeId() 
{   return PointContactImpl::classTypeId(); }

PointContactImpl::PointContactImpl
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2, 
    Vec3& location, Vec3& normal, Real radius1, Real radius2, Real depth)
:   ContactImpl(surf1, surf2), location(location), normal(normal), 
    radius1(radius1), radius2(radius2), effectiveRadius(std::sqrt(radius1*radius2)), depth(depth) {}

PointContactImpl::PointContactImpl
   (ContactSurfaceIndex surf1, ContactSurfaceIndex surf2,
    Vec3& location, Vec3& normal, Real radius, Real depth)
:   ContactImpl(surf1, surf2), location(location), normal(normal),
    radius1(radius), radius2(radius), effectiveRadius(radius), depth(depth) {}


