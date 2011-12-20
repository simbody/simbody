/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-11 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "SimTKmath.h"

#include <algorithm>
using std::pair; using std::make_pair;
#include <iostream>
using std::cout; using std::endl;
#include <set>


namespace SimTK {

//==============================================================================
//                     HALFSPACE-SPHERE CONTACT TRACKER
//==============================================================================
// Cost is 21 flops if no contact, 67 with contact.
bool ContactTracker::HalfSpaceSphere::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GS, 
    const ContactGeometry& geoSphere,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoSphere.getTypeId()==ContactGeometry::Sphere::classTypeId(),
       "ContactTracker::HalfSpaceSphere::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Sphere& sphere = 
        ContactGeometry::Sphere::getAs(geoSphere);

    const Rotation R_HG = ~X_GH.R(); // inverse rotation; no flops

    // p_HC is vector from H origin to S's center C
    const Vec3 p_HC = R_HG*(X_GS.p() - X_GH.p()); // 18 flops
    
    // Calculate depth of sphere center C given that the halfspace occupies 
    // all of x>0 space.
    const Real r = sphere.getRadius();
    const Real depth = p_HC[0] + r;   // 1 flop

    if (depth <= -cutoff) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // Calculate the rest of the X_HS transform as required by Contact.
    const Transform X_HS(R_HG*X_GS.R(), p_HC); // 45 flops
    const UnitVec3 normal_H(Vec3(-1,0,0), true); // 0 flops
    const Vec3     origin_H = Vec3(depth/2, p_HC[1], p_HC[2]); // 1 flop
    // The surfaces are contacting (or close enough to be interesting).
    // The sphere's radius is also the effective radius.
    currentStatus = CircularPointContact(priorStatus.getSurface1(), Infinity,
                                         priorStatus.getSurface2(), r,
                                         X_HS, r, depth, origin_H, normal_H);
    return true; // success
}

bool ContactTracker::HalfSpaceSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceSphere::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceSphere::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceSphere::initializeContact() not implemented yet."); 
    return false; }



//==============================================================================
//                     HALFSPACE-ELLIPSOID CONTACT TRACKER
//==============================================================================
// Cost is ~135 flops if no contact, ~425 with contact.
// The contact point on the ellipsoid must be the unique point that has its
// outward-facing normal in the opposite direction of the half space normal.
// We can find that point very fast and see how far it is from the half
// space surface. If it is close enough, we'll evaluate the curvatures at
// that point in preparation for generating forces with Hertz theory.
bool ContactTracker::HalfSpaceEllipsoid::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GE, 
    const ContactGeometry& geoEllipsoid,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoEllipsoid.getTypeId()==ContactGeometry::Ellipsoid::classTypeId(),
       "ContactTracker::HalfSpaceEllipsoid::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Ellipsoid& ellipsoid = 
        ContactGeometry::Ellipsoid::getAs(geoEllipsoid);

    // Our half space occupies the +x half so the normal is -x.
    const Transform X_HE = ~X_GH*X_GE; // 63 flops
    // Halfspace normal is -x, so the ellipsoid normal we're looking for is
    // in the half space's +x direction.
    const UnitVec3& n_E = (~X_HE.R()).x(); // halfspace normal in E
    const Vec3 Q_E = ellipsoid.findPointWithThisUnitNormal(n_E); // 40 flops
    const Vec3 Q_H = X_HE*Q_E; // Q measured from half space origin (18 flops)
    const Real depth = Q_H[0]; // x > 0 is penetrated

    if (depth <= -cutoff) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // The surfaces are contacting (or close enough to be interesting).
    // The ellipsoid's principal curvatures k at the contact point are also
    // the curvatures of the contact paraboloid since the half space doesn't
    // add anything interesting.
    Transform X_EQ; Vec2 k;
    ellipsoid.findParaboloidAtPointWithNormal(Q_E, n_E, X_EQ, k); // 220 flops

    // We have the contact paraboloid expressed in frame Q but Qz=n_E has the
    // wrong sign since we have to express it using the half space normal.
    // We have to end up with a right handed frame, so one of x or y has
    // to be negated too. (6 flops)
    Rotation& R_EQ = X_EQ.updR();
    R_EQ.setRotationColFromUnitVecTrustMe(ZAxis, -R_EQ.z()); // changing X_EQ
    R_EQ.setRotationColFromUnitVecTrustMe(XAxis, -R_EQ.x());

    // Now the frame is pointing in the right direction. Measure and express in 
    // half plane frame, then shift origin to half way between contact point 
    // Q on the undeformed ellipsoid and the corresponding contact point P 
    // on the undeformed half plane surface. It's easier to do this shift
    // in H since it is in the -Hx direction.
    Transform X_HC = X_HE*X_EQ; X_HC.updP()[0] -= depth/2; // 65 flops

    currentStatus = EllipticalPointContact(priorStatus.getSurface1(),
                                           priorStatus.getSurface2(),
                                           X_HE, X_HC, k, depth);
    return true; // success
}

bool ContactTracker::HalfSpaceEllipsoid::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceEllipsoid::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceEllipsoid::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceEllipsoid::initializeContact() not implemented yet."); 
    return false; }



//==============================================================================
//                       SPHERE-SPHERE CONTACT TRACKER
//==============================================================================
// Cost is 12 flops if no contact, 139 if contact
bool ContactTracker::SphereSphere::trackContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, 
    const ContactGeometry& geoSphere1,
    const Transform& X_GS2, 
    const ContactGeometry& geoSphere2,
    Real                   cutoff, // 0 for contact
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoSphere1.getTypeId()==ContactGeometry::Sphere::classTypeId()
        && geoSphere2.getTypeId()==ContactGeometry::Sphere::classTypeId(),
       "ContactTracker::SphereSphere::trackContact()");

    // No need for an expensive dynamic casts here; we know what we have.
    const ContactGeometry::Sphere& sphere1 = 
        ContactGeometry::Sphere::getAs(geoSphere1);
    const ContactGeometry::Sphere& sphere2 = 
        ContactGeometry::Sphere::getAs(geoSphere2);

    currentStatus.clear();

    // Find the vector from sphere center C1 to C2, expressed in G.
    const Vec3 p_12_G = X_GS2.p() - X_GS1.p(); // 3 flops
    const Real d2 = p_12_G.normSqr();          // 5 flops
    const Real r1 = sphere1.getRadius();
    const Real r2 = sphere2.getRadius();
    const Real rr = r1 + r2;                    // 1 flop

    // Quick check. If separated we can return nothing, unless we were
    // in contact last time in which case we have to return one last
    // Contact indicating that contact has been broken and by how much.
    if (d2 > square(rr+cutoff)) {       // 3 flops
        if (!priorStatus.getContactId().isValid())
            return true; // successful return: still separated

        const Real separation = std::sqrt(d2) - rr;   // > cutoff, ~25 flops
        const Transform X_S1S2(~X_GS1.R()*X_GS2.R(), 
                               ~X_GS1.R()*p_12_G);    // 60 flops
        currentStatus = BrokenContact(priorStatus.getSurface1(),
                                      priorStatus.getSurface2(),
                                      X_S1S2, separation);
        return true;
    }

    const Real d = std::sqrt(d2); // ~20 flops
    if (d < SignificantReal) {    // 1 flop
        // TODO: If the centers are coincident we should use past information
        // to determine the most likely normal. For now just fail.
        return false;
    }

    const Transform X_S1S2(~X_GS1.R()*X_GS2.R(), 
                           ~X_GS1.R()*p_12_G);// 60 flops
    const Vec3& p_12 = X_S1S2.p(); // center-to-center vector in S1

    const Real depth = rr - d; // >0 for penetration (1 flop)
    const Real r     = r1*r2/rr; // r=r1r2/(r1+r2)=1/(1/r1+1/r2) ~20 flops

    const UnitVec3 normal_S1(p_12/d, true); // 1/ + 3* = ~20 flops
    const Vec3     origin_S1 = (r1 - depth/2)*normal_S1; // 5 flops

    // The surfaces are contacting (or close enough to be interesting).
    currentStatus = CircularPointContact(priorStatus.getSurface1(), r1,
                                         priorStatus.getSurface2(), r2,
                                         X_S1S2, r, depth, 
                                         origin_S1, normal_S1);
    return true; // success
}

bool ContactTracker::SphereSphere::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereSphere::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::SphereSphere::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereSphere::initializeContact() not implemented yet."); 
    return false; }




//==============================================================================
//                  HALFSPACE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::HalfSpaceTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH, 
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GM, 
    const ContactGeometry& geoMesh,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::HalfSpace::isInstance(geoHalfSpace)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh),
       "ContactTracker::HalfSpaceTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::HalfSpaceTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::TriangleMesh& mesh = 
        ContactGeometry::TriangleMesh::getAs(geoMesh);

    // Transform giving mesh (S2) frame in the halfspace (S1) frame.
    const Transform X_HM = (~X_GH)*X_GM; 

    // Normal is halfspace -x direction; xdir is first column of R_MH.
    // That's a unit vector and -unitvec is also a unit vector so this
    // doesn't require normalization.
    const UnitVec3 hsNormal_M = -(~X_HM.R()).x();
    // Find the height of the halfspace face along the normal, measured
    // from the mesh origin.
    const Real hsFaceHeight_M = dot((~X_HM).p(), hsNormal_M);
    // Now collect all the faces that are all or partially below the 
    // halfspace surface.
    std::set<int> insideFaces;
    processBox(mesh, mesh.getOBBTreeNode(), X_HM, 
               hsNormal_M, hsFaceHeight_M, insideFaces);
    
    if (insideFaces.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_HM, 
                                        std::set<int>(), insideFaces);
    return true; // success
}

bool ContactTracker::HalfSpaceTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::HalfSpaceTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::HalfSpaceTriangleMesh::initializeContact() not implemented yet."); 
    return false; }



// Check a single OBB and its contents (recursively) against the halfspace,
// appending any penetrating faces to the insideFaces list.
void ContactTracker::HalfSpaceTriangleMesh::processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Transform& X_HM, const UnitVec3& hsNormal_M, Real hsFaceHeight_M, 
    std::set<int>& insideFaces) const 
{   // First check against the node's bounding box.
    
    const OrientedBoundingBox& bounds = node.getBounds();
    const Transform& X_MB = bounds.getTransform(); // box frame in mesh
    const Vec3 p_BC = bounds.getSize()/2; // from box origin corner to center
    // Express the half space normal in the box frame, then reflect it into
    // the first (+,+,+) quadrant where it is the normal of a different 
    // but symmetric and more convenient half space.
    const UnitVec3 octant1hsNormal_B = (~X_MB.R()*hsNormal_M).abs();
    // Dot our octant1 radius p_BC with our octant1 normal to get
    // the extent of the box from its center in the direction of the octant1
    // reflection of the halfspace.
    const Real extent = dot(p_BC, octant1hsNormal_B);
    // Compute the height of the box center over the mesh origin,
    // measured along the real halfspace normal.
    const Vec3 boxCenter_M       = X_MB*p_BC;
    const Real boxCenterHeight_M = dot(boxCenter_M, hsNormal_M);
    // Subtract the halfspace surface position to get the height of the 
    // box center over the halfspace.
    const Real boxCenterHeight = boxCenterHeight_M - hsFaceHeight_M;
    if (boxCenterHeight >= extent)
        return;                             // no penetration
    if (boxCenterHeight <= -extent) {
        addAllTriangles(node, insideFaces); // box is entirely in halfspace
        return;
    }
    
    // Box is partially penetrated into halfspace. If it is not a leaf node, 
    // check its children.
    if (!node.isLeafNode()) {
        processBox(mesh, node.getFirstChildNode(), X_HM, hsNormal_M, 
                   hsFaceHeight_M, insideFaces);
        processBox(mesh, node.getSecondChildNode(), X_HM, hsNormal_M, 
                   hsFaceHeight_M, insideFaces);
        return;
    }
    
    // This is a leaf OBB node that is penetrating, so some of its triangles
    // may be penetrating.
    const Array_<int>& triangles = node.getTriangles();
    for (int i = 0; i < (int) triangles.size(); i++) {
        for (int vx=0; vx < 3; ++vx) {
            const int   vertex         = mesh.getFaceVertex(triangles[i], vx);
            const Vec3& vertexPos      = mesh.getVertexPosition(vertex);
            const Real  vertexHeight_M = dot(vertexPos, hsNormal_M);
            if (vertexHeight_M < hsFaceHeight_M) {
                insideFaces.insert(triangles[i]);
                break; // done with this face
            }
        }
    }
}

void ContactTracker::HalfSpaceTriangleMesh::addAllTriangles
   (const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    std::set<int>& insideFaces) const 
{
    if (node.isLeafNode()) {
        const Array_<int>& triangles = node.getTriangles();
        for (int i = 0; i < (int) triangles.size(); i++)
            insideFaces.insert(triangles[i]);
    }
    else {
        addAllTriangles(node.getFirstChildNode(), insideFaces);
        addAllTriangles(node.getSecondChildNode(), insideFaces);
    }
}





//==============================================================================
//                  SPHERE - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::SphereTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GS, 
    const ContactGeometry& geoSphere,
    const Transform&       X_GM, 
    const ContactGeometry& geoMesh,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::Sphere::isInstance(geoSphere)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh),
       "ContactTracker::SphereTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::SphereTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::Sphere&          sphere = 
        ContactGeometry::Sphere::getAs(geoSphere);
    const ContactGeometry::TriangleMesh&    mesh = 
        ContactGeometry::TriangleMesh::getAs(geoMesh);

    // Transform giving mesh (M) frame in the sphere (S) frame.
    const Transform X_SM = ~X_GS*X_GM; 

    // Want the sphere center measured and expressed in the mesh frame.
    const Vec3 p_MC = (~X_SM).p();
    std::set<int> insideFaces;
    processBox(mesh, mesh.getOBBTreeNode(), p_MC, square(sphere.getRadius()), 
               insideFaces);
    
    if (insideFaces.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_SM, 
                                        std::set<int>(), insideFaces);
    return true; // success
}

// Check a single OBB and its contents (recursively) against the sphere
// whose center location in M and radius squared is given, appending any 
// penetrating faces to the insideFaces list.
void ContactTracker::SphereTriangleMesh::processBox
   (const ContactGeometry::TriangleMesh&              mesh, 
    const ContactGeometry::TriangleMesh::OBBTreeNode& node, 
    const Vec3& center_M, Real radius2, 
    std::set<int>& insideFaces) const 
{   // First check against the node's bounding box.

    const Vec3 nearest_M = node.getBounds().findNearestPoint(center_M);
    if ((nearest_M-center_M).normSqr() >= radius2)
        return; // no intersection possible
    
    // Bounding box is penetrating. If it's not a leaf node, check its children.
    if (!node.isLeafNode()) {
        processBox(mesh, node.getFirstChildNode(), center_M, radius2,
                   insideFaces);
        processBox(mesh, node.getSecondChildNode(), center_M, radius2,
                   insideFaces);
        return;
    }
    
    // This is a leaf node that may be penetrating; check the triangles.
    const Array_<int>& triangles = node.getTriangles();
    for (unsigned i = 0; i < triangles.size(); i++) {
        Vec2 uv;
        Vec3 nearest_M = mesh.findNearestPointToFace
                                    (center_M, triangles[i], uv);
        if ((nearest_M-center_M).normSqr() < radius2)
            insideFaces.insert(triangles[i]);
    }
}

bool ContactTracker::SphereTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::SphereTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::SphereTriangleMesh::initializeContact() not implemented yet."); 
    return false; }






//==============================================================================
//               TRIANGLE MESH - TRIANGLE MESH CONTACT TRACKER
//==============================================================================
// Cost is TODO
bool ContactTracker::TriangleMeshTriangleMesh::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GM1, 
    const ContactGeometry& geoMesh1,
    const Transform&       X_GM2, 
    const ContactGeometry& geoMesh2,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   ContactGeometry::TriangleMesh::isInstance(geoMesh1)
        && ContactGeometry::TriangleMesh::isInstance(geoMesh2),
       "ContactTracker::TriangleMeshTriangleMesh::trackContact()");

    // We can't handle a "proximity" test, only penetration. 
    SimTK_ASSERT_ALWAYS(cutoff==0,
       "ContactTracker::TriangleMeshTriangleMesh::trackContact()");

    // No need for an expensive dynamic cast here; we know what we have.
    const ContactGeometry::TriangleMesh& mesh1 = 
        ContactGeometry::TriangleMesh::getAs(geoMesh1);
    const ContactGeometry::TriangleMesh& mesh2 = 
        ContactGeometry::TriangleMesh::getAs(geoMesh2);

    // Transform giving mesh2 (M2) frame in the mesh1 (M1) frame.
    const Transform X_M1M2 = ~X_GM1*X_GM2; 
    std::set<int> insideFaces1, insideFaces2;

    // Get M2's bounding box in M1's frame.
    const OrientedBoundingBox 
        mesh2Bounds_M1 = X_M1M2*mesh2.getOBBTreeNode().getBounds();

    // Find the faces that are actually intersecting faces on the other
    // surface (this doesn't yet include faces that may be completely buried).
    findIntersectingFaces(mesh1, mesh2, 
                          mesh1.getOBBTreeNode(), mesh2.getOBBTreeNode(), 
                          mesh2Bounds_M1, X_M1M2, insideFaces1, insideFaces2);
    
    // It should never be the case that one set of faces is empty and the
    // other isn't, however it is conceivable that roundoff error could cause
    // this to happen so we'll check both lists.
    if (insideFaces1.empty() && insideFaces2.empty()) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }
    
    // There was an intersection. We now need to identify every triangle and 
    // vertex of each mesh that is inside the other mesh. We found the border
    // intersections above; now we have to fill in the buried faces.
    findBuriedFaces(mesh1, mesh2, ~X_M1M2, insideFaces1);
    findBuriedFaces(mesh2, mesh1,  X_M1M2, insideFaces2);

    currentStatus = TriangleMeshContact(priorStatus.getSurface1(), 
                                        priorStatus.getSurface2(), 
                                        X_M1M2, 
                                        insideFaces1, insideFaces2);
    return true; // success
}

void ContactTracker::TriangleMeshTriangleMesh::
findIntersectingFaces
   (const ContactGeometry::TriangleMesh&                mesh1, 
    const ContactGeometry::TriangleMesh&                mesh2,
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node1, 
    const ContactGeometry::TriangleMesh::OBBTreeNode&   node2, 
    const OrientedBoundingBox&                          node2Bounds_M1,
    const Transform&                                    X_M1M2, 
    std::set<int>&                                      triangles1, 
    std::set<int>&                                      triangles2) const 
{   // See if the bounding boxes intersect.
    
    if (!node1.getBounds().intersectsBox(node2Bounds_M1))
        return;
    
    // If either node is not a leaf node, process the children recursively.
    
    if (!node1.isLeafNode()) {
        if (!node2.isLeafNode()) {
            const OrientedBoundingBox firstChildBounds = 
                X_M1M2*node2.getFirstChildNode().getBounds();
            const OrientedBoundingBox secondChildBounds = 
                X_M1M2*node2.getSecondChildNode().getBounds();
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        }
        else {
            findIntersectingFaces(mesh1, mesh2, node1.getFirstChildNode(), node2, node2Bounds_M1, X_M1M2, triangles1, triangles2);
            findIntersectingFaces(mesh1, mesh2, node1.getSecondChildNode(), node2, node2Bounds_M1, X_M1M2, triangles1, triangles2);
        }
        return;
    }
    else if (!node2.isLeafNode()) {
        const OrientedBoundingBox firstChildBounds = 
            X_M1M2*node2.getFirstChildNode().getBounds();
        const OrientedBoundingBox secondChildBounds = 
            X_M1M2*node2.getSecondChildNode().getBounds();
        findIntersectingFaces(mesh1, mesh2, node1, node2.getFirstChildNode(), firstChildBounds, X_M1M2, triangles1, triangles2);
        findIntersectingFaces(mesh1, mesh2, node1, node2.getSecondChildNode(), secondChildBounds, X_M1M2, triangles1, triangles2);
        return;
    }
    
    // These are both leaf nodes, so check triangles for intersections.
    
    const Array_<int>& node1triangles = node1.getTriangles();
    const Array_<int>& node2triangles = node2.getTriangles();
    for (unsigned i = 0; i < node2triangles.size(); i++) {
        const int face2 = node2triangles[i];
        Vec3 a1 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 0));
        Vec3 a2 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 1));
        Vec3 a3 = X_M1M2*mesh2.getVertexPosition(mesh2.getFaceVertex(face2, 2));
        const Geo::Triangle A(a1,a2,a3);
        for (unsigned j = 0; j < node1triangles.size(); j++) {
            const int face1 = node1triangles[j];
            const Vec3& b1 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 0));
            const Vec3& b2 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 1));
            const Vec3& b3 = mesh1.getVertexPosition(mesh1.getFaceVertex(face1, 2));
            const Geo::Triangle B(b1,b2,b3);
            if (A.overlapsTriangle(B)) 
            {   // The triangles intersect.
                triangles1.insert(face1);
                triangles2.insert(face2);
            }
        }
    }
}

static const int Outside  = -1;
static const int Unknown  =  0;
static const int Boundary =  1;
static const int Inside   =  2;

void ContactTracker::TriangleMeshTriangleMesh::
findBuriedFaces(const ContactGeometry::TriangleMesh&    mesh,       // M 
                const ContactGeometry::TriangleMesh&    otherMesh,  // O
                const Transform&                        X_OM, 
                std::set<int>&                          insideFaces) const 
{  
    // Find which faces are inside.
    // We're passed in the list of Boundary faces, that is, those faces of
    // "mesh" that intersect faces of "otherMesh".
    Array_<int> faceType(mesh.getNumFaces(), Unknown);
    for (std::set<int>::iterator iter = insideFaces.begin(); 
                                 iter != insideFaces.end(); ++iter)
        faceType[*iter] = Boundary;

    for (int i = 0; i < (int) faceType.size(); i++) {
        if (faceType[i] == Unknown) {
            // Trace a ray from its center to determine whether it is inside.           
            const Vec3     origin_O    = X_OM    * mesh.findCentroid(i);
            const UnitVec3 direction_O = X_OM.R()* mesh.getFaceNormal(i);
            Real distance;
            int face;
            Vec2 uv;
            if (   otherMesh.intersectsRay(origin_O, direction_O, distance, 
                                           face, uv) 
                && ~direction_O*otherMesh.getFaceNormal(face) > 0) 
            {
                faceType[i] = Inside;
                insideFaces.insert(i);
            } else
                faceType[i] = Outside;
            
            // Recursively mark adjacent inside or outside Faces.           
            tagFaces(mesh, faceType, insideFaces, i, 0);
        }
    }
}

//TODO: the following method uses depth-first recursion to iterate through
//unmarked faces. For a large mesh this was observed to produce a stack
//overflow in OpenSim. Here we limit the recursion depth; after we get that
//deep we'll pop back out and do another expensive intersectsRay() test in
//the method above.
static const int MaxRecursionDepth = 500;

void ContactTracker::TriangleMeshTriangleMesh::
tagFaces(const ContactGeometry::TriangleMesh&   mesh, 
         Array_<int>&                           faceType,
         std::set<int>&                         triangles, 
         int                                    index,
         int                                    depth) const 
{
    for (int i = 0; i < 3; i++) {
        const int edge = mesh.getFaceEdge(index, i);
        const int face = (mesh.getEdgeFace(edge, 0) == index 
                            ? mesh.getEdgeFace(edge, 1) 
                            : mesh.getEdgeFace(edge, 0));
        if (faceType[face] == Unknown) {
            faceType[face] = faceType[index];
            if (faceType[index] > 0)
                triangles.insert(face);
            if (depth < MaxRecursionDepth)
                tagFaces(mesh, faceType, triangles, face, depth+1);
        }
    }
}


bool ContactTracker::TriangleMeshTriangleMesh::predictContact
   (const Contact&         priorStatus,
    const Transform& X_GS1, const SpatialVec& V_GS1, const SpatialVec& A_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2, const SpatialVec& A_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               predictedStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::TriangleMeshTriangleMesh::predictContact() not implemented yet."); 
    return false; }

bool ContactTracker::TriangleMeshTriangleMesh::initializeContact
   (const Transform& X_GS1, const SpatialVec& V_GS1,
    const ContactGeometry& surface1,
    const Transform& X_GS2, const SpatialVec& V_GS2,
    const ContactGeometry& surface2,
    Real                   cutoff,
    Real                   intervalOfInterest,
    Contact&               contactStatus) const
{   SimTK_ASSERT_ALWAYS(!"implemented",
    "ContactTracker::TriangleMeshTriangleMesh::initializeContact() not implemented yet."); 
    return false; }



} // namespace SimTK

