#ifndef SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_
#define SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
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
#include "simmath/internal/common.h"
#include "simmath/internal/Contact.h"

namespace SimTK {

//==============================================================================
//                              CONTACT TRACKER
//==============================================================================
/** A ContactTracker implements an algorithm for detecting overlaps or 
potential overlaps between pairs of ContactGeometry objects, and managing
Contact objects that track individual contacts as they evolve through
time. This class is used internally by ContractTrackerSubsystem and there 
usually is no reason to access it directly. The exception is if you are 
defining a new ContactGeometry subclass. In that case, you will also need to 
define one or more ContactTrackers to detect collisions with your new geometry
type, then register it with the ContactTrackerSubsystem. 

The result of a ContactTracker when applied to a pair of contact
surfaces, is either a determination that the surfaces are not in contact,
or a Contact object describing their contact interaction. There are different
types of these Contact objects (for example, PointContact, LineContact, 
MeshContact) and the same algorithm may result in different kinds of Contact 
under different circumstances. At each evaluation, the caller passes in the 
previous Contact object, if any, that was associated with two ContactSurfaces,
then receives an update from the algorithm.

Note that ContactTrackers that manage dissimilar geometry type pairs expect 
the two types in a particular order, e.g. (halfspace,sphere) rather than
(sphere,halfspace) but are used for all contacts involving that pair of 
types. It is up to the ContactTrackerSubsystem to ensure that the contact
surfaces are presented in the correct order regardless of how they are 
encountered. The Contact objects that are created and managed by trackers
always have their (surface1,surface2) pairs in the order required by the
tracker that handles those types. **/
class SimTK_SIMMATH_EXPORT ContactTracker {
public:
class HalfSpaceSphere;
class HalfSpaceEllipsoid;
class SphereSphere;
class HalfSpaceTriangleMesh;
class SphereTriangleMesh;
class TriangleMeshTriangleMesh;

/** Base class constructor for use by the concrete classes. **/
ContactTracker(ContactGeometryTypeId typeOfSurface1,
               ContactGeometryTypeId typeOfSurface2)
:   m_surfaceTypes(typeOfSurface1, typeOfSurface2)
{
}

/** Return the pair of contact geometry type ids handled by this tracker,
in the order that they must be presented to the tracker's methods. **/
const std::pair<ContactGeometryTypeId,ContactGeometryTypeId>&
getContactGeometryTypeIds() const {return m_surfaceTypes;}

virtual ~ContactTracker() {}

/** The ContactTrackerSubsystem will invoke this method for any pair of
contact surfaces that is already being tracked, or for which the static broad 
phase analysis indicated that they might be in contact now. Only position 
information is available. Note that the arguments and Contact object surfaces
must be ordered by geometry type id as required by this tracker. **/
virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const = 0;

/** The ContactTrackerSubsystem will invoke this method for any tracked pair of
contact surfaces that is still not in contact after trackContact() looked at
it, or any untracked pair for which the dynamic broad phase indicated that 
they might be in contact within the interval of interest. Position, velocity, 
and acceleration information may be used. Ordering must be correct as 
discussed for trackContact(). **/
virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const = 0;

/** At the beginning of a simulation we will have no past information to
help disambiguate tricky contact situations. This method may use current
position and velocity information in heuristics for guessing the contact
status between the indicated pair of surfaces. Ordering must be correct as 
discussed for trackContact(). **/
virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const = 0;

//--------------------------------------------------------------------------
                                private:
// This tracker should be called only for surfaces of these two types,
// in this order.
std::pair<ContactGeometryTypeId,ContactGeometryTypeId> m_surfaceTypes;
};



//==============================================================================
//                     HALFSPACE-SPHERE CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::Sphere, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceSphere 
:   public ContactTracker {
public:
HalfSpaceSphere() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::Sphere::classTypeId()) {}

virtual ~HalfSpaceSphere() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                     HALFSPACE-ELLIPSOID CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::Ellipsoid, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceEllipsoid 
:   public ContactTracker {
public:
HalfSpaceEllipsoid() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::Ellipsoid::classTypeId()) {}

virtual ~HalfSpaceEllipsoid() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                       SPHERE-SPHERE CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two ContactGeometry::Sphere
objects. **/
class SimTK_SIMMATH_EXPORT ContactTracker::SphereSphere 
:   public ContactTracker {
public:
SphereSphere() 
:   ContactTracker(ContactGeometry::Sphere::classTypeId(),
                   ContactGeometry::Sphere::classTypeId()) {}

virtual ~SphereSphere() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,
    const Transform& X_GS2, 
    const ContactGeometry& surface2,
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;
};



//==============================================================================
//                 HALFSPACE-TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::HalfSpace
and a ContactGeometry::TriangleMesh, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::HalfSpaceTriangleMesh
:   public ContactTracker {
public:
HalfSpaceTriangleMesh() 
:   ContactTracker(ContactGeometry::HalfSpace::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~HalfSpaceTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // the half space
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // the mesh
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void processBox(const ContactGeometry::TriangleMesh&              mesh, 
                const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
                const Transform& X_HM, const UnitVec3& hsNormal_M, 
                Real hsFaceHeight_M, std::set<int>& insideFaces) const;
void addAllTriangles(const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
                     std::set<int>& insideFaces) const; 
};



//==============================================================================
//                 SPHERE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between a ContactGeometry::Sphere
and a ContactGeometry::TriangleMesh, in that order. **/
class SimTK_SIMMATH_EXPORT ContactTracker::SphereTriangleMesh
:   public ContactTracker {
public:
SphereTriangleMesh() 
:   ContactTracker(ContactGeometry::Sphere::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~SphereTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // the sphere
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // the mesh
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Vec3& center_M, Real radius2,   
    std::set<int>& insideFaces) const ;
};



//==============================================================================
//             TRIANGLE MESH - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
/** This ContactTracker handles contacts between two 
ContactGeometry::TriangleMesh surfaces. **/
class SimTK_SIMMATH_EXPORT ContactTracker::TriangleMeshTriangleMesh
:   public ContactTracker {
public:
TriangleMeshTriangleMesh() 
:   ContactTracker(ContactGeometry::TriangleMesh::classTypeId(),
                   ContactGeometry::TriangleMesh::classTypeId()) {}

virtual ~TriangleMeshTriangleMesh() {}

virtual bool trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& surface1,    // mesh1
    const Transform& X_GS2, 
    const ContactGeometry& surface2,    // mesh2
    Real                   cutoff,
    Contact&               currentStatus) const;

virtual bool predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const;

virtual bool initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const;

private:
void findIntersectingFaces
   (const ContactGeometry::TriangleMesh&                mesh1, 
    const ContactGeometry::TriangleMesh&                mesh2,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node1, 
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node2, 
    const OrientedBoundingBox&                          node2Bounds_M1,
    const Transform&                                    X_M1M2, 
    std::set<int>&                                      insideFaces1, 
    std::set<int>&                                      insideFaces2) const; 
void findBuriedFaces
   (const ContactGeometry::TriangleMesh&    mesh,
    const ContactGeometry::TriangleMesh&    otherMesh,
    const Transform&                        X_OM, 
    std::set<int>&                          insideFaces) const;
void tagFaces(const ContactGeometry::TriangleMesh&   mesh, 
              Array_<int>&                           faceType,
              std::set<int>&                         triangles, 
              int                                    index,
              int                                    depth) const;
};

} // namespace SimTK

#endif // SimTK_SIMMATH_CONTACT_TRACKER_SUBSYSTEM_H_
