/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
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

#include "SimTKmath.h"

#include <algorithm>
using std::pair; using std::make_pair;
#include <iostream>
using std::cout; using std::endl;
#include <set>

// Define this if you want to see voluminous output from MPR (XenoCollide)
//#define MPR_DEBUG


namespace SimTK {

//==============================================================================
//                             CONTACT TRACKER
//==============================================================================


//------------------------------------------------------------------------------
//                 ESTIMATE IMPLICIT PAIR CONTACT USING MPR
//------------------------------------------------------------------------------
// Generate a rough guess at the contact points. Point P is returned in A's
// frame and point Q is in B's frame, but all the work here is done in frame A.
// See Snethen, G. "XenoCollide: Complex Collision Made Simple", Game
// Programming Gems 7, pp.165-178.
//
// Note that we intend to use this on smooth convex objects. That means there
// can be no guarantee about how many iterations this will take to converge; it
// may take a very long time for objects that are just barely touching. On the
// other hand we only need an approximate answer because we're going to polish
// the solution to machine precision using a Newton iteration afterwards.

namespace {

// Given a direction in shape A's frame, calculate the support point of the
// Minkowski difference shape A-B in that direction. To do that we find A's
// support point P in that direction, and B's support point Q in the opposite
// direction; i.e. what would be -B's support in the original direction. Then
// return the vector v=(P-Q) expressed in A's frame.
struct Support {
    Support(const ContactGeometry& shapeA, const ContactGeometry& shapeB,
            const Transform& X_AB, const UnitVec3& dirInA)
    :   shapeA(shapeA), shapeB(shapeB), X_AB(X_AB)
    {   computeSupport(dirInA); }

    // Calculate a new set of supports for the same shapes.
    void computeSupport(const UnitVec3& dirInA) {
        const UnitVec3 dirInB = ~X_AB.R() * dirInA; // 15 flops
        dir = dirInA;
        A = shapeA.calcSupportPoint( dirInA); // varies; 40 flops for ellipsoid
        B = shapeB.calcSupportPoint(-dirInB); //             "
        v = A - X_AB*B;                       // 21 flops
        depth = dot(v,dir);                   //  5 flops
    }

    Support& operator=(const Support& src) {
        dir = src.dir; A = src.A; B = src.B; v = src.v;
        return *this;
    }

    void swap(Support& other) {
        std::swap(dir, other.dir);
        std::swap(A,other.A); std::swap(B,other.B); std::swap(v,other.v);
    }

    void getResult(Vec3& pointP_A, Vec3& pointQ_B, UnitVec3& normalInA) const
    {   pointP_A = A; pointQ_B = B; normalInA = dir; }

    UnitVec3 dir; // A support direction, exp in A
    Vec3 A; // support point in direction dir on shapeA, expressed in A
    Vec3 B; // support point in direction -dir on shapeB, expressed in B
    Vec3 v; // support point A-B in Minkowski difference shape, exp. in A
    Real depth; // v . dir (positive when origin is below support plane)

private:
    const ContactGeometry& shapeA;
    const ContactGeometry& shapeB;
    const Transform&       X_AB;
};

const Real MPRAccuracy     = Real(.05);   // 5% of the surface dimensions
const Real SqrtMPRAccuracy = Real(.2236); // roughly sqrt(MPRAccuracy)
}


/*static*/ bool ContactTracker::
estimateConvexImplicitPairContactUsingMPR
   (const ContactGeometry& shapeA, const ContactGeometry& shapeB,
    const Transform& X_AB,
    Vec3& pointP, Vec3& pointQ, UnitVec3& dirInA, int& numIterations)
{
    const Rotation& R_AB = X_AB.R();
    numIterations = 0;

    // Get some cheap, rough scaling information.
    // TODO: ideally this would be done using local curvature information.
    // A reasonable alternative would be to use the smallest dimension of the
    // shape's oriented bounding box.
    // We're going to assume the smallest radius is 1/4 of the bounding
    // sphere radius.
    Vec3 cA, cB; Real rA, rB;
    shapeA.getBoundingSphere(cA,rA); shapeB.getBoundingSphere(cB,rB);
    const Real lengthScale = Real(0.25)*std::min(rA,rB);
    const Real areaScale   = Real(1.5)*square(lengthScale); // ~area of octant

    // If we determine the depth to within a small fraction of the scale,
    // or localize the contact area to a small fraction of the surface area,
    // that is good enough and we can stop.
    const Real depthGoal = MPRAccuracy*lengthScale;
    const Real areaGoal  = SqrtMPRAccuracy*areaScale;

    #ifdef MPR_DEBUG
    printf("\nMPR acc=%g: r=%g, depthGoal=%g areaGoal=%g\n",
        MPRAccuracy, std::min(rA,rB), depthGoal, areaGoal);
    #endif

    // Phase 1: Portal Discovery
    // -------------------------

    // Compute an interior point v0 that is known to be inside the Minkowski
    // difference, and a ray dir1 directed from that point to the origin.
    const Vec3 v0 = (  Support(shapeA, shapeB, X_AB, UnitVec3( XAxis)).v
                     + Support(shapeA, shapeB, X_AB, UnitVec3(-XAxis)).v)/2;

    if (v0 == 0) {
        // This is a pathological case: the two objects are directly on top of
        // each other with their centers at exactly the same place. Just
        // return *some* vaguely plausible contact.
        pointP = shapeA.calcSupportPoint(      UnitVec3( XAxis));
        pointQ = shapeB.calcSupportPoint(~R_AB*UnitVec3(-XAxis)); // in B
        dirInA = XAxis;
        return true;
    }


    // Support 1's direction is initially the "origin ray" that points from
    // interior point v0 to the origin.
    Support s1(shapeA, shapeB, X_AB, UnitVec3(-v0));

    // Test for NaN once and get out to avoid getting stuck in loops below.
    if (isNaN(s1.depth)) {
        pointP = pointQ = NaN;
        dirInA = UnitVec3();
        return false;
    }

    if (s1.depth <= 0) { // origin outside 1st support plane
        s1.getResult(pointP, pointQ, dirInA);
        return false;
    }

    if (s1.v % v0 == 0) {   // v0 perpendicular to support plane; origin inside
        s1.getResult(pointP, pointQ, dirInA);
        return true;
    }

    // Find support point perpendicular to plane containing origin, interior
    // point v0, and first support v1.
    Support s2(shapeA, shapeB, X_AB, UnitVec3(s1.v % v0));
    if (s2.depth <= 0) { // origin is outside 2nd support plane
        s2.getResult(pointP, pointQ, dirInA);
        return false;
    }

    // Find support perpendicular to plane containing interior point v0 and
    // first two support points v1 and v2. Make sure it is on the side that
    // is closer to the origin; fix point ordering if necessary.
    UnitVec3 d3 = UnitVec3((s1.v-v0)%(s2.v-v0));
    if (~d3*v0 > 0) { // oops -- picked the wrong side
        s1.swap(s2);
        d3 = -d3;
    }
    Support s3(shapeA, shapeB, X_AB, d3);
    if (s3.depth <= 0) {  // origin is outside 3rd support plane
        s3.getResult(pointP, pointQ, dirInA);
        return false;
    }

    // We now have a candidate portal (triangle v1,v2,v3). We have to refine it
    // until the origin ray -v0 intersects the candidate. Check against the
    // three planes of the tetrahedron that contain v0. By construction above
    // we know the origin is inside the v0,v1,v2 face already.

    // We should find a candidate portal very fast, probably in 1 or 2
    // iterations. We'll allow an absurdly large number of
    // iterations and then abort just to make sure we don't get stuck in an
    // infinite loop.
    const Real MaxCandidateIters = 100; // should never get anywhere near this
    int candidateIters = 0;
    while (true) {
        ++candidateIters;
        SimTK_ERRCHK_ALWAYS(candidateIters <= MaxCandidateIters,
            "ContactTracker::estimateConvexImplicitPairContactUsingMPR()",
            "Unable to find a candidate portal; should never happen.");

        if (~v0*(s1.v % s3.v) < -SignificantReal)
            s2 = s3; // origin outside v0,v1,v3 face; replace v2
        else if (~v0*(s3.v % s2.v) < -SignificantReal)
            s1 = s3; // origin outside v0,v2,v3 face; replace v1
        else
            break;

        // Choose new candidate. The keepers are in v1 and v2; get a new v3.
        s3.computeSupport(UnitVec3((s1.v-v0) % (s2.v-v0)));
    }

    // Phase 2: Portal Refinement
    // --------------------------

    // We have a portal (triangle v1,v2,v3) that the origin ray passes through.
    // Now we need to refine v1,v2,v3 until we have portal such that the origin
    // is inside the tetrahedron v0,v1,v2,v3.

    const int MinTriesToFindSeparatingPlane = 5;
    const int MaxTriesToImproveContact = 5;
    int triesToFindSeparatingPlane=0, triesToImproveContact=0;
    while (true) {
        ++numIterations;

        // Get the normal to the portal, oriented in the general direction of
        // the origin ray (i.e., the outward normal).
        const Vec3 portalVec = (s2.v-s1.v) % (s3.v-s1.v);
        const Real portalArea = portalVec.norm(), ooPortalArea = 1/portalArea;
        UnitVec3 portalDir(portalVec*ooPortalArea, true);
        if (~portalDir*v0 > 0)
            portalDir = -portalDir;

        // Any portal vertex is a vector from the origin to the portal plane.
        // Dot one of them with the portal outward normal to get the origin
        // depth (positive if inside).
        const Real depth = ~s1.v*portalDir;

        // Find new support in portal direction.
        Support s4(shapeA, shapeB, X_AB, portalDir);
        if (s4.depth <= 0) { // origin is outside new support plane
            s4.getResult(pointP, pointQ, dirInA);
            return false;
        }

        const Real depthChange = std::abs(s4.depth - depth);


        bool mustReturn=false, okToReturn=false;
        if (depth >= 0) {   // We found a contact.
            mustReturn = (++triesToImproveContact >= MaxTriesToImproveContact);
            okToReturn = true;
        } else {            // No contact yet.
            okToReturn =
               (++triesToFindSeparatingPlane >= MinTriesToFindSeparatingPlane);
            mustReturn = false;
        }
        bool accuracyAchieved =
            (depthChange <= depthGoal || portalArea <= areaGoal);

        #ifdef MPR_DEBUG
        printf("  depth=%g, change=%g area=%g changeFrac=%g areaFrac=%g\n",
            depth, depthChange, portalArea,
            depthChange/depthGoal, portalArea/areaGoal);
        printf("    accuracyAchieved=%d okToReturn=%d mustReturn=%d\n",
            accuracyAchieved, okToReturn, mustReturn);
        #endif

        if (mustReturn || (okToReturn && accuracyAchieved)) {
            // The origin is inside the portal, so we have an intersection.
            // Compute the barycentric coordinates of the origin ray's
            // intersection with the portal, and map back to the two surfaces.
            const Vec3 origin = v0+v0*(~portalDir*(s1.v-v0)/(~portalDir*v0));
            const Real area1 = ~portalDir*((s2.v-origin)%(s3.v-origin));
            const Real area2 = ~portalDir*((s3.v-origin)%(s1.v-origin));
            const Real u = area1*ooPortalArea;
            const Real v = area2*ooPortalArea;
            const Real w = 1-u-v;

            // Compute the contact points in their own shape's frame.
            pointP = u*s1.A + v*s2.A + w*s3.A;
            pointQ = u*s1.B + v*s2.B + w*s3.B;
            dirInA = portalDir;
            return true;
        }

        // We know the origin ray entered the (v1,v2,v3,v4) tetrahedron via the
        // (v1,v2,v3) face (the portal). New portal is the face that it exits.
        const Vec3 v4v0 = s4.v % v0;
        if (~s1.v * v4v0 > 0) {
            if (~s2.v * v4v0 > 0) s1 = s4; // v4,v2,v3 (new portal)
            else                  s3 = s4; // v1,v2,v4
        } else {
            if (~s3.v * v4v0 > 0) s2 = s4; // v1,v4,v3
            else                  s1 = s4; // v4,v2,v3 again
        }
    }
}


//------------------------------------------------------------------------------
//                            REFINE IMPLICIT PAIR
//------------------------------------------------------------------------------
// We have a rough estimate of the contact points. Use Newton iteration to
// refine them. If the surfaces are separated, this will find the points of
// closest approach. If contacting, this will find the points of maximium
// penetration.

// Returns true if the desired accuracy is achieved, regardless of whether the
// surfaces are separated or in contact.
/*static*/ bool ContactTracker::
refineImplicitPair
   (const ContactGeometry& shapeA, Vec3& pointP,    // in/out (in A)
    const ContactGeometry& shapeB, Vec3& pointQ,    // in/out (in B)
    const Transform& X_AB, Real accuracyRequested,
    Real& accuracyAchieved, int& numIterations)
{
    // If the initial guess is very bad, or the ellipsoids pathological we
    // may have to crawl along for a while at the beginning.
    const int MaxSlowIterations = 8;
    const int MaxIterations = MaxSlowIterations + 8;
    const Real MinStepFrac = Real(1e-6); // if we can't take at least this
                                         // fraction of Newton step, give it up

    Vec6 err = findImplicitPairError(shapeA, pointP, shapeB, pointQ, X_AB);
    accuracyAchieved = err.norm();
    Vector errVec(6, &err[0], true); // share space with err

    Mat66 J; // to hold the Jacobian
    Matrix JMat(6,6,6,&J(0,0)); // share space with J

    Vec6 delta;
    Vector deltaVec(6, &delta[0], true); // share space with delta

    numIterations = 0;
    while (   accuracyAchieved > accuracyRequested
           && numIterations < MaxIterations)
    {
        ++numIterations;
        J = calcImplicitPairJacobian(shapeA, pointP, shapeB, pointQ,
                                     X_AB, err);

        // Try to use LU factorization; fall back to QTZ if singular.
        FactorLU lu(JMat);
        if (!lu.isSingular())
            lu.solve(errVec, deltaVec);     // writes into delta also
        else {
            FactorQTZ qtz(JMat, SqrtEps);
            qtz.solve(errVec, deltaVec);    // writes into delta also
        }

        // Line search for safety in case starting guess bad. Don't accept
        // any move that makes things worse.
        Real f = 2; // scale back factor
        Vec3 oldP = pointP, oldQ = pointQ;
        Real oldAccuracyAchieved = accuracyAchieved;
        do {
            f /= 2;
            pointP = oldP - f*delta.getSubVec<3>(0);
            pointQ = oldQ - f*delta.getSubVec<3>(3);
            err = findImplicitPairError(shapeA, pointP, shapeB, pointQ, X_AB);
            accuracyAchieved = err.norm();
        } while (accuracyAchieved > oldAccuracyAchieved && f > MinStepFrac);

        const bool noProgressMade = (accuracyAchieved > oldAccuracyAchieved);

        if (noProgressMade) { // Restore best points and fall through.
            pointP = oldP;
            pointQ = oldQ;
        }

        if (    noProgressMade
            || (f < 1 && numIterations >= MaxSlowIterations)) // Too slow.
        {
            // We don't appear to be getting anywhere. Just project the
            // points onto the surfaces and then exit.
            bool inside;
            UnitVec3 normal;
            pointP = shapeA.findNearestPoint(pointP, inside, normal);
            pointQ = shapeB.findNearestPoint(pointQ, inside, normal);
            err = findImplicitPairError(shapeA, pointP, shapeB, pointQ, X_AB);
            accuracyAchieved = err.norm();
            break;
        }
    }

    return accuracyAchieved <= accuracyRequested;
}


//------------------------------------------------------------------------------
//                          FIND IMPLICIT PAIR ERROR
//------------------------------------------------------------------------------
// We're given two implicitly-defined shapes A and B and candidate surface
// contact points P (for A) and Q (for B). There are six equations that real
// contact points should satisfy. Here we return the errors in each of those
// equations; if all six terms are zero these are contact points.
//
// Per profiling, this gets called a lot at runtime, so keep it tight!
/*static*/ Vec6 ContactTracker::
findImplicitPairError
   (const ContactGeometry& shapeA, const Vec3& pointP,  // in A
    const ContactGeometry& shapeB, const Vec3& pointQ,  // in B
    const Transform& X_AB)
{
    // Compute the function value and normal vector for each object.

    const Function& fA = shapeA.getImplicitFunction();
    const Function& fB = shapeB.getImplicitFunction();

    // Avoid some heap allocations by using stack arrays.
    Vec3 xData;
    Vector x(3, &xData[0], true); // shares space with xdata

    int compData; // just one integer
    ArrayView_<int> components(&compData, &compData+1);

    Vec3 gradA, gradB;
    xData = pointP; // writes into Vector x also
    for (int i = 0; i < 3; i++) {
        components[0] = i;
        gradA[i] = fA.calcDerivative(components, x);
    }
    Real errorA = fA.calcValue(x);
    xData = pointQ; // writes into Vector x also
    for (int i = 0; i < 3; i++) {
        components[0] = i;
        gradB[i] = fB.calcDerivative(components, x);
    }
    Real errorB = fB.calcValue(x);

    // Construct a coordinate frame for each object.
    // TODO: this needs some work to make sure it is as stable as possible
    // so the perpendicularity errors are stable as the solution advances
    // and especially as the Jacobian is calculated by perturbation.

    UnitVec3 nA(-gradA);
    UnitVec3 nB(X_AB.R()*(-gradB));
    UnitVec3 uA(std::abs(nA[0])>Real(0.5)? nA%Vec3(0, 1, 0) : nA%Vec3(1, 0, 0));
    UnitVec3 uB(std::abs(nB[0])>Real(0.5)? nB%Vec3(0, 1, 0) : nB%Vec3(1, 0, 0));
    Vec3 vA = nA%uA; // Already a unit vector, so we don't need to normalize it.
    Vec3 vB = nB%uB;

    // Compute the error vector. The components indicate, in order, that nA
    // must be perpendicular to both tangents of object B, that the separation
    // vector should be zero or perpendicular to the tangents of object A, and
    // that both points should be on their respective surfaces.

    Vec3 delta = pointP-X_AB*pointQ;
    return Vec6(~nA*uB, ~nA*vB, ~delta*uA, ~delta*vA, errorA, errorB);
}


//------------------------------------------------------------------------------
//                        CALC IMPLICIT PAIR JACOBIAN
//------------------------------------------------------------------------------
// Differentiate the findImplicitPairError() with respect to changes in the
// locations of points P and Q in their own surface frames.

/*static*/ Mat66 ContactTracker::calcImplicitPairJacobian
   (const ContactGeometry& shapeA, const Vec3& pointP,
    const ContactGeometry& shapeB, const Vec3& pointQ,
    const Transform& X_AB, const Vec6& err0)
{
    Real dt = SqrtEps;

    Vec3 d1 = dt*Vec3(1, 0, 0);
    Vec3 d2 = dt*Vec3(0, 1, 0);
    Vec3 d3 = dt*Vec3(0, 0, 1);
    Vec6 err1 = findImplicitPairError(shapeA, pointP+d1, shapeB, pointQ, X_AB)
                - err0;
    Vec6 err2 = findImplicitPairError(shapeA, pointP+d2, shapeB, pointQ, X_AB)
                - err0;
    Vec6 err3 = findImplicitPairError(shapeA, pointP+d3, shapeB, pointQ, X_AB)
                - err0;
    Vec6 err4 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d1, X_AB)
                - err0;
    Vec6 err5 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d2, X_AB)
                - err0;
    Vec6 err6 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d3, X_AB)
                - err0;
    return Mat66(err1, err2, err3, err4, err5, err6) / dt;

    // This is a Central Difference derivative (you should use a somewhat
    // larger dt in this case). However, I haven't seen any evidence that this
    // helps, even for some very eccentric ellipsoids. (sherm 20130408)
    //Vec6 err1 = findImplicitPairError(shapeA, pointP+d1, shapeB, pointQ, X_AB)
    //            - findImplicitPairError(shapeA, pointP-d1, shapeB, pointQ, X_AB);
    //Vec6 err2 = findImplicitPairError(shapeA, pointP+d2, shapeB, pointQ, X_AB)
    //            - findImplicitPairError(shapeA, pointP-d2, shapeB, pointQ, X_AB);
    //Vec6 err3 = findImplicitPairError(shapeA, pointP+d3, shapeB, pointQ, X_AB)
    //            - findImplicitPairError(shapeA, pointP-d3, shapeB, pointQ, X_AB);
    //Vec6 err4 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d1, X_AB)
    //            - findImplicitPairError(shapeA, pointP, shapeB, pointQ-d1, X_AB);
    //Vec6 err5 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d2, X_AB)
    //            - findImplicitPairError(shapeA, pointP, shapeB, pointQ-d2, X_AB);
    //Vec6 err6 = findImplicitPairError(shapeA, pointP, shapeB, pointQ+d3, X_AB)
    //            - findImplicitPairError(shapeA, pointP, shapeB, pointQ-d3, X_AB);
    //return Mat66(err1, err2, err3, err4, err5, err6) / (2*dt);
}



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



//==============================================================================
//                     HALFSPACE-BRICK CONTACT TRACKER
//==============================================================================
// Cost is XX flops if no contact, YY with contact.
bool ContactTracker::HalfSpaceBrick::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH,
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GB,
    const ContactGeometry& geoBrick,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoBrick.getTypeId()==ContactGeometry::Brick::classTypeId(),
       "ContactTracker::HalfSpaceBrick::trackContact()");

    const ContactGeometry::HalfSpace& halfSpace =
        ContactGeometry::HalfSpace::getAs(geoHalfSpace);
    // This is the half-space outward normal in its own frame.
    const UnitVec3 n_H = halfSpace.getNormal();

    const ContactGeometry::Brick& brick =
        ContactGeometry::Brick::getAs(geoBrick);
    const Geo::Box box = brick.getGeoBox();

    const Transform X_HB = ~X_GH * X_GB; // 63 flops

    // Negative normal direction in box frame locates the bottommost vertex.
    const UnitVec3 nn_B = ~X_HB.R()*(-n_H);
    const int lowestVertex = box.findSupportVertex(nn_B);
    const Vec3 pt_B = box.getVertexPos(lowestVertex);

    const Vec3 pt_H = X_HB * pt_B; // find vertex location in H
    const Real height = dot(pt_H, n_H); // -ve means penetrated

    if (height >= cutoff) {  // 1 flop
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    currentStatus = BrickHalfSpaceContact(priorStatus.getSurface1(),
                                          priorStatus.getSurface2(),
                                          X_HB,
                                          lowestVertex, -height);
    return true; // success
}



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





//==============================================================================
//                   HALFSPACE-CONVEX IMPLICIT CONTACT TRACKER
//==============================================================================
// The contact point on the convex implicit surface must be the unique point on
// that surface that has its outward-facing normal in the opposite direction of
// the half space normal. We will use the calcSupportPoint() method using the
// negated half-space normal to find the desired point. You should only be
// using this tracker for convex surfaces that can provide a high-accuracy
// support point very fast. If the point is close enough, we'll evaluate the
// curvatures at that point using the calcSurfacePrincipalCurvatures() method,
// in preparation for generating forces with Hertz theory. This will return an
// elliptical point contact.
bool ContactTracker::HalfSpaceConvexImplicit::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GH,
    const ContactGeometry& geoHalfSpace,
    const Transform&       X_GS,
    const ContactGeometry& geoImplSurface,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT
       (   geoHalfSpace.getTypeId()==ContactGeometry::HalfSpace::classTypeId()
        && geoImplSurface.isConvex() && geoImplSurface.isSmooth(),
       "ContactTracker::HalfSpaceConvexImplicit::trackContact()");

    // Our half space occupies the +x half so the normal is -x.
    const Transform X_HS = ~X_GH*X_GS; // 63 flops
    // Halfspace normal is -x, so the surface normal we're looking for is
    // in the half space's +x direction.
    const UnitVec3& n_S = (~X_HS.R()).x(); // halfspace normal in S
    const Vec3 Q_S = geoImplSurface.calcSupportPoint(n_S); // cost depends on S
    const Vec3 Q_H = X_HS*Q_S; // Q measured from half space origin (18 flops)
    const Real depth = Q_H[0]; // x > 0 is penetrated

    if (depth <= -cutoff) {  // 2 flops
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // The surfaces are contacting (or close enough to be interesting).
    // The ellipsoid's principal curvatures k at the contact point are also
    // the curvatures of the contact paraboloid since the half space doesn't
    // add anything interesting.
    Transform X_SQ; Vec2 k;
    X_SQ.updP() = Q_S;
    Rotation& R_SQ = X_SQ.updR();
    geoImplSurface.calcSurfacePrincipalCurvatures(Q_S, k, R_SQ); // cost?

    // We have the contact paraboloid expressed in frame Q but Qz=n_E has the
    // wrong sign since we have to express it using the half space normal.
    // We have to end up with a right handed frame, so one of x or y has
    // to be negated too. (6 flops)
    R_SQ.setRotationColFromUnitVecTrustMe(ZAxis, -R_SQ.z()); // changing X_SQ
    R_SQ.setRotationColFromUnitVecTrustMe(XAxis, -R_SQ.x());

    // Now the frame is pointing in the right direction. Measure and express in
    // half plane frame, then shift origin to half way between contact point Q
    // on the undeformed implicit surface and the corresponding contact point P
    // on the undeformed half plane surface. It's easier to do this shift
    // in H since it is in the -Hx direction.
    Transform X_HC = X_HS*X_SQ; X_HC.updP()[0] -= depth/2; // 65 flops

    currentStatus = EllipticalPointContact(priorStatus.getSurface1(),
                                           priorStatus.getSurface2(),
                                           X_HS, X_HC, k, depth);
    return true; // success
}


//==============================================================================
//               CONVEX IMPLICIT SURFACE PAIR CONTACT TRACKER
//==============================================================================
// This will return an elliptical point contact.
bool ContactTracker::ConvexImplicitPair::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GA,
    const ContactGeometry& shapeA,
    const Transform&       X_GB,
    const ContactGeometry& shapeB,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   shapeA.isConvex() && shapeA.isSmooth()
        && shapeB.isConvex() && shapeB.isSmooth(),
       "ContactTracker::ConvexImplicitPair::trackContact()");

    // We'll work in the shape A frame.
    const Transform X_AB = ~X_GA*X_GB; // 63 flops
    const Rotation& R_AB = X_AB.R();

    // 1. Get a rough guess at the contact points P and Q and contact normal.
    Vec3 pointP_A, pointQ_B; // on A and B, resp.
    UnitVec3 norm_A;
    int numMPRIters;
    const bool mightBeContact = estimateConvexImplicitPairContactUsingMPR
                                   (shapeA, shapeB, X_AB,
                                    pointP_A, pointQ_B, norm_A, numMPRIters);

    #ifdef MPR_DEBUG
    std::cout << "MPR: " << (mightBeContact?"MAYBE":"NO") << std::endl;
    std::cout << "  P=" << X_GA*pointP_A << " Q=" << X_GB*pointQ_B << std::endl;
    std::cout << "  N=" << X_GA.R()*norm_A << std::endl;
    #endif

    if (!mightBeContact) {
        currentStatus.clear(); // definitely not touching
        return true; // successful return
    }

    // 2. Refine the contact points to near machine precision.
    const Real accuracyRequested = SignificantReal;
    Real accuracyAchieved; int numNewtonIters;
    bool converged = refineImplicitPair(shapeA, pointP_A, shapeB, pointQ_B,
        X_AB, accuracyRequested, accuracyAchieved, numNewtonIters);

    const Vec3 pointQ_A = X_AB*pointQ_B;  // Q on B, measured & expressed in A

    // 3. Compute the curvatures and surface normals of the two surfaces at
    //    P and Q. Once we have the first normal we can check whether there was
    //    actually any contact and duck out early if not.
    Rotation R_AP; Vec2 curvatureP;
    shapeA.calcCurvature(pointP_A, curvatureP, R_AP);

    // If the surfaces are in contact then the vector from Q on surface B
    // (supposedly inside A) to P on surface A (supposedly inside B) should be
    // aligned with the outward normal on A.
    const Real depth = dot(pointP_A-pointQ_A, R_AP.z());

    #ifdef MPR_DEBUG
    printf("MPR %2d iters, Newton %2d iters->accuracy=%g depth=%g\n",
        numMPRIters, numNewtonIters, accuracyAchieved, depth);
    #endif

    if (depth <= 0) {
        currentStatus.clear(); // not touching
        return true; // successful return
    }

    // The surfaces are in contact.

    Rotation R_BQ; Vec2 curvatureQ;
    shapeB.calcCurvature(pointQ_B, curvatureQ, R_BQ);
    const UnitVec3 maxDirB_A(R_AB*R_BQ.x()); // re-express in A

    // 4. Compute the effective contact frame C and corresponding relative
    //    curvatures.
    Transform X_AC; Vec2 curvatureC;

    // Define the contact frame origin to be at the midpoint of P and Q.
    X_AC.updP() = (pointP_A+pointQ_A)/2;

    // Determine the contact frame orientations and composite curvatures.
    ContactGeometry::combineParaboloids(R_AP, curvatureP,
                                        maxDirB_A, curvatureQ,
                                        X_AC.updR(), curvatureC);

    // 5. Return the elliptical point contact for force generation.
    currentStatus = EllipticalPointContact(priorStatus.getSurface1(),
                                           priorStatus.getSurface2(),
                                           X_AB, X_AC, curvatureC, depth);
    return true; // success
}



//==============================================================================
//                GENERAL IMPLICIT SURFACE PAIR CONTACT TRACKER
//==============================================================================
// This will return an elliptical point contact. TODO: should return a set
// of contacts.
//TODO: not implemented yet -- this is just the Convex-convex code here.
bool ContactTracker::GeneralImplicitPair::trackContact
   (const Contact&         priorStatus,
    const Transform&       X_GA,
    const ContactGeometry& shapeA,
    const Transform&       X_GB,
    const ContactGeometry& shapeB,
    Real                   cutoff,
    Contact&               currentStatus) const
{
    SimTK_ASSERT_ALWAYS
       (   shapeA.isSmooth() && shapeB.isSmooth(),
       "ContactTracker::GeneralImplicitPair::trackContact()");

    // TODO: this won't work unless the shapes are actually convex.
    return ConvexImplicitPair(shapeA.getTypeId(),shapeB.getTypeId())
            .trackContact(priorStatus, X_GA, shapeA, X_GB, shapeB,
                          cutoff, currentStatus);
}


} // namespace SimTK

