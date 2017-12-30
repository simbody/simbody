/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/** @file
Non-inline static methods from the Geo::Point_ class. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/LinearAlgebra.h"
#include "simmath/Optimizer.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Box.h"
#include "simmath/internal/Geo_Sphere.h"

#include <set>
#include <cstdio>
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

//==============================================================================
//                              GEO :: POINT
//==============================================================================

// These local helpers find the minimum or maximum of several values and 
// return which one was the extreme value.
template <class P>
inline static void minOf(P a, P b, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
}

template <class P>
inline static void minOf(P a, P b, P c, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
    if (c<minVal) minVal=c, which=2;
}

template <class P>
inline static void minOf(P a, P b, P c, P d, P& minVal, int& which) {
    minVal=a; which=0;
    if (b<minVal) minVal=b, which=1;
    if (c<minVal) minVal=c, which=2;
    if (d<minVal) minVal=d, which=3;
}
template <class P>
inline static void maxOf(P a, P b, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
}

template <class P>
inline static void maxOf(P a, P b, P c, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
    if (c>maxVal) maxVal=c, which=2;
}

template <class P>
inline static void maxOf(P a, P b, P c, P d, P& maxVal, int& which) {
    maxVal=a; which=0;
    if (b>maxVal) maxVal=b, which=1;
    if (c>maxVal) maxVal=c, which=2;
    if (d>maxVal) maxVal=d, which=3;
}


// This helper replaces the 0-based indices in "which" with corresponding
// indices in "map". We assume map is long enough.
inline static void fixWhich(const int* map, Array_<int>& which) {
    for (unsigned i=0; i<which.size(); ++i)
        which[i] = map[which[i]];
}

// Given an array of point locations, create an indirect array of pointers
// to those point locations.
template <class P> void
makeIndirect(const Array_< Vec<3,P> >&  points,
             Array_<const Vec<3,P>*>&   pointsIndirect) {
    pointsIndirect.resize(points.size());
    for (unsigned i=0; i < points.size(); ++i)
        pointsIndirect[i] = &points[i];
}


// Given a set of points, find the one that is the furthest in a given
// direction. There must be at least one point in the set.
template <class P> /*static*/ void Geo::Point_<P>::
findSupportPoint(const Array_<Vec3P>& points, const UnitVec3P& direction,
                 int& most, RealP& maxCoord) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", "findSupportPoint()",
        "There must be at least one point in the set.");
    most=0; maxCoord = dot(points[0],direction);
    for (int i=1; i < (int)points.size(); ++i) {
        const RealP coord = dot(points[i],direction);
        if (coord > maxCoord) most=i, maxCoord=coord;
    }
}

template <class P> /*static*/ void Geo::Point_<P>::
findSupportPointIndirect(const Array_<const Vec3P*>& points, 
                         const UnitVec3P& direction,
                         int& most, RealP& maxCoord) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", 
        "findSupportPointIndirect()",
        "There must be at least one point in the set.");
    most=0; maxCoord = dot(*points[0],direction);
    for (int i=1; i < (int)points.size(); ++i) {
        const RealP coord = dot(*points[i],direction);
        if (coord > maxCoord) most=i, maxCoord=coord;
    }
}

// Given a set of points, find the two points that are the most extreme along
// a given direction (not necessarily distinct). There must be at least one 
// point in the set.
template <class P> /*static*/ void Geo::Point_<P>::
findExtremePoints(const Array_<Vec3P>& points, const UnitVec3P& direction,
                  int& least, int& most, RealP& leastCoord, RealP& mostCoord) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", "findExtremePoints()",
        "There must be at least one point in the set.");
    least=most=0; 
    leastCoord = dot(points[0],direction); mostCoord=leastCoord;
    for (int i=1; i < (int)points.size(); ++i) {
        const Vec3P& p = points[i];
        const RealP coord = dot(p,direction);
        if (coord < leastCoord) least=i, leastCoord=coord;
        if (coord > mostCoord)  most=i,  mostCoord=coord;
    }
}

template <class P> /*static*/ void Geo::Point_<P>::
findExtremePointsIndirect(const Array_<const Vec3P*>& points, 
                          const UnitVec3P& direction,
                          int& least, int& most, 
                          RealP& leastCoord, RealP& mostCoord) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", 
        "findExtremePointsIndirect()",
        "There must be at least one point in the set.");
    least=most=0; 
    leastCoord = dot(*points[0],direction); mostCoord=leastCoord;
    for (int i=1; i < (int)points.size(); ++i) {
        const Vec3P& p = *points[i];
        const RealP coord = dot(p,direction);
        if (coord < leastCoord) least=i, leastCoord=coord;
        if (coord > mostCoord)  most=i,  mostCoord=coord;
    }
}

template <class P> /*static*/ Vec<3,P> Geo::Point_<P>::
calcCentroid(const Array_<Vec3P>& points) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", "calcCentroid()",
        "There must be at least one point in the set.");

    const int   n   = (int)points.size();
    const RealP oon = RealP(1)/n;   // ~10 flops

    Vec3P centroid(0);
    for (int i=0; i < n; ++i)
        centroid += points[i];      // 3 flops
    centroid *= oon;

    return centroid;
}

template <class P> /*static*/ Vec<3,P> Geo::Point_<P>::
calcCentroidIndirect(const Array_<const Vec3P*>& points) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", "calcCentroidIndirect()",
        "There must be at least one point in the set.");

    const int   n   = (int)points.size();
    const RealP oon = RealP(1)/n;   // ~10 flops

    Vec3P centroid(0);
    for (int i=0; i < n; ++i)
        centroid += *points[i];      // 3 flops
    centroid *= oon;

    return centroid;
}

template <class P> /*static*/ void Geo::Point_<P>::
calcCovariance(const Array_<Vec3P>& points_F,
               Vec3P& centroid, SymMat33P& covariance) {
    Array_<const Vec3P*> indirect_F;
    makeIndirect(points_F, indirect_F);
    calcCovarianceIndirect(indirect_F, centroid, covariance);
}

template <class P> /*static*/ void Geo::Point_<P>::
calcCovarianceIndirect(const Array_<const Vec3P*>& points_F,
                       Vec3P& centroid, SymMat33P& covariance) {
    SimTK_APIARGCHECK(!points_F.empty(),"Geo::Point_", 
        "calcCovarianceIndirect()",
        "There must be at least one point in the set.");

    const int   n   = (int)points_F.size();
    const RealP oon = RealP(1)/n;   // ~10 flops

    // Pass 1: find the centroid
    centroid = RealP(0);
    for (int i=0; i < n; ++i)
        centroid += *points_F[i];    // 3 flops
    centroid *= oon;

    // Pass 2: calculate the covariance matrix.
    covariance.setToZero();
    Vec3P& diag  = covariance.updDiag();
    Vec3P& lower = covariance.updLower(); // 1,0 2,0 2,1
    for (int i=0; i < n; ++i) {
        const Vec3P p = *points_F[i] - centroid;
        diag[0]  += p[0]*p[0]; diag[1]  += p[1]*p[1]; diag[2]  += p[2]*p[2];
        lower[0] += p[0]*p[1]; lower[1] += p[0]*p[2]; lower[2] += p[1]*p[2];
    }

    covariance *= oon;
}


template <class P> /*static*/ void Geo::Point_<P>::
calcPrincipalComponents(const Array_<Vec3P>& points_F,
                        TransformP&          X_FP)  {
    Array_<const Vec3P*> indirect_F;
    makeIndirect(points_F, indirect_F);
    calcPrincipalComponentsIndirect(indirect_F, X_FP);
}

template <class P> /*static*/ void Geo::Point_<P>::
calcPrincipalComponentsIndirect(const Array_<const Vec3P*>& points_F,
                                TransformP&                 X_FP) {
    SimTK_APIARGCHECK(!points_F.empty(),"Geo::Point_", 
        "calcPrincipalComponentsIndirect()",
        "There must be at least one point in the set.");

    Vec3P     centroid;
    SymMat33P covariance;
    calcCovarianceIndirect(points_F, X_FP.updP(), covariance);

    // Calculate eigenvalues and eigenvectors.
    const Mat33P cov33(covariance);
    Matrix_<P> cov(cov33);
    Vector_< std::complex<P> > evals;
    Matrix_< std::complex<P> > evecs;
    Eigen(cov).getAllEigenValuesAndVectors(evals, evecs);

    // Find the largest and smallest eigenvalues and corresponding vectors.
    const Vec3P vals(evals[0].real(), evals[1].real(), evals[2].real());
    int minIx, maxIx; RealP minVal, maxVal;
    minOf(vals[0], vals[1], vals[2], minVal, minIx);
    maxOf(vals[0], vals[1], vals[2], maxVal, maxIx);

    if (maxIx == minIx) // eigenvalues must all be the same
        maxIx = (minIx+1) % 3; // pick a different axis

    // Eigenvectors for a real symmetric matrix are perpendicular and
    // normalized already.
    UnitVec3P longAxis( Vec3P(evecs(0,maxIx).real(), evecs(1,maxIx).real(), 
                              evecs(2,maxIx).real()), 
                        true); // "trust me" to prevent normalizing
    UnitVec3P shortAxis(Vec3P(evecs(0,minIx).real(), evecs(1,minIx).real(), 
                              evecs(2,minIx).real()), 
                        true);
    X_FP.updR() = RotationP(longAxis, XAxis, shortAxis, YAxis);
}

//==============================================================================
//                       CALC AXIS ALIGNED BOUNDING BOX
//==============================================================================


template <class P> /*static*/ void Geo::Point_<P>::
findAxisAlignedExtremePoints(const Array_<Vec3P>& points,
                             int least[3], int most[3],
                             Vec3P& low, Vec3P& high) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", 
        "findAxisAlignedExtremePoints()",
        "There must be at least one point in the set.");

    low=points[0]; high=points[0];
    for (int i=0; i<3; i++) least[i]=most[i]=0; // point 0 most extreme so far
    
    for (int i=1; i < (int)points.size(); ++i) {
        const Vec3P& p = points[i];
        if (p[0] < low [0]) low [0]=p[0], least[0]=i;
        if (p[0] > high[0]) high[0]=p[0], most [0]=i;
        if (p[1] < low [1]) low [1]=p[1], least[1]=i;
        if (p[1] > high[1]) high[1]=p[1], most [1]=i;
        if (p[2] < low [2]) low [2]=p[2], least[2]=i;
        if (p[2] > high[2]) high[2]=p[2], most [2]=i;
    }
}

// This signature taking an array of pointers to points rather than the
// points themselves.

template <class P> /*static*/ void Geo::Point_<P>::
findAxisAlignedExtremePointsIndirect(const Array_<const Vec3P*>& points,
                                     int least[3], int most[3],
                                     Vec3P& low, Vec3P& high) {
    SimTK_APIARGCHECK(!points.empty(),"Geo::Point_", 
        "findAxisAlignedExtremePointsIndirect()",
        "There must be at least one point in the set.");

    low=(*points[0]); high=(*points[0]);
    for (int i=0; i<3; i++) least[i]=most[i]=0; // point 0 most extreme so far
    
    for (int i=1; i < (int)points.size(); ++i) {
        const Vec3P& p = (*points[i]);
        if (p[0] < low [0]) low [0]=p[0], least[0]=i;
        if (p[0] > high[0]) high[0]=p[0], most [0]=i;
        if (p[1] < low [1]) low [1]=p[1], least[1]=i;
        if (p[1] > high[1]) high[1]=p[1], most [1]=i;
        if (p[2] < low [2]) low [2]=p[2], least[2]=i;
        if (p[2] > high[2]) high[2]=p[2], most [2]=i;
    }
}


// This signature taking an array of pointers to points rather than the
// points themselves.
template <class P> /*static*/
Geo::AlignedBox_<P> Geo::Point_<P>::
calcAxisAlignedBoundingBoxIndirect(const Array_<const Vec<3,P>*>& points,
                                   Array_<int>&                   support) {
    const unsigned npoints = points.size();
    if (npoints == 0) {
        const P inf = NTraits<P>::getInfinity();
        support.clear();
        return Geo::AlignedBox_<P>(Vec3P(inf),Vec3P(0));
    }

    // Find the extreme points along each of the coordinate axes.
    int least[3], most[3];
    Vec3P low, high;
    findAxisAlignedExtremePointsIndirect(points, least, most, low, high);

    // Sort the support points and eliminate duplicates.
    std::set<int> supportSet;
    for (unsigned i=0; i<3; ++i) {
        supportSet.insert(least[i]);
        supportSet.insert(most[i]);
    }
    support.assign(supportSet.begin(), supportSet.end());

    const Vec3P ctr = (low+high)/2;
    // Make sure nothing gets left out due to roundoff.
    const Vec3P hdim(std::max(high[0]-ctr[0], ctr[0]-low[0]),
                     std::max(high[1]-ctr[1], ctr[1]-low[1]),
                     std::max(high[2]-ctr[2], ctr[2]-low[2]));

    return AlignedBox_<P>(ctr, hdim).stretchBoundary(); 
}


// This signature takes an array of points, creates an array of pointers to 
// those points and calls the other signature.
template <class P> /*static*/ 
Geo::AlignedBox_<P> Geo::Point_<P>::
calcAxisAlignedBoundingBox(const Array_< Vec<3,P> >& points,
                           Array_<int>&         support) {
    Array_<const Vec<3,P>*> indirect;
    makeIndirect(points, indirect);
    return calcAxisAlignedBoundingBoxIndirect(indirect, support);
}

//==============================================================================
//                       CALC ORIENTED BOUNDING BOX
//==============================================================================

template <class P> /*static*/ void Geo::Point_<P>::
findOrientedExtremePoints(const Array_<Vec3P>& points_F, 
                          const RotationP& R_FB,
                          int least[3], int most[3],
                          Vec3P& low_B, Vec3P& high_B) {
    SimTK_APIARGCHECK(!points_F.empty(),"Geo::Point_", 
        "findOrientedExtremePoints()",
        "There must be at least one point in the set.");
    const RotationP R_BF = ~R_FB;

    low_B=R_BF*points_F[0]; high_B=R_BF*points_F[0];
    for (int i=0; i<3; i++) least[i]=most[i]=0; // point 0 most extreme so far
    
    for (int i=1; i < (int)points_F.size(); ++i) {
        const Vec3P p = R_BF*points_F[i];
        if (p[0] < low_B [0]) low_B [0]=p[0], least[0]=i;
        if (p[0] > high_B[0]) high_B[0]=p[0], most [0]=i;
        if (p[1] < low_B [1]) low_B [1]=p[1], least[1]=i;
        if (p[1] > high_B[1]) high_B[1]=p[1], most [1]=i;
        if (p[2] < low_B [2]) low_B [2]=p[2], least[2]=i;
        if (p[2] > high_B[2]) high_B[2]=p[2], most [2]=i;
    }
}


template <class P> /*static*/ void Geo::Point_<P>::
findOrientedExtremePointsIndirect(const Array_<const Vec3P*>& points_F, 
                                  const RotationP& R_FB,
                                  int least[3], int most[3],
                                  Vec3P& low_B, Vec3P& high_B) {
    SimTK_APIARGCHECK(!points_F.empty(),"Geo::Point_", 
        "findOrientedExtremePointsIndirect()",
        "There must be at least one point in the set.");
    const RotationP R_BF = ~R_FB;

    low_B=R_BF*(*points_F[0]); high_B=R_BF*(*points_F[0]);
    for (int i=0; i<3; i++) least[i]=most[i]=0; // point 0 most extreme so far
    
    for (int i=1; i < (int)points_F.size(); ++i) {
        const Vec3P p = R_BF*(*points_F[i]);
        if (p[0] < low_B [0]) low_B [0]=p[0], least[0]=i;
        if (p[0] > high_B[0]) high_B[0]=p[0], most [0]=i;
        if (p[1] < low_B [1]) low_B [1]=p[1], least[1]=i;
        if (p[1] > high_B[1]) high_B[1]=p[1], most [1]=i;
        if (p[2] < low_B [2]) low_B [2]=p[2], least[2]=i;
        if (p[2] > high_B[2]) high_B[2]=p[2], most [2]=i;
    }
}

// Differentiate the above function being careful to note the low accuracy
// when running in single precision.
template <class P>
class VolumeGradient : public Differentiator::GradientFunction {
    typedef Vec<3,P> Vec3P;
public:
    VolumeGradient(const Array_<const Vec3P*>&  points, 
                   const Rotation_<P>&          R_FB0)
    :   Differentiator::GradientFunction(3,(Real)Geo::getDefaultTol<P>()), 
        p(&points), R_FB0(R_FB0) {}

    // This function calculates the volume as a function of rotation angles
    // relative to the starting frame B0.
    int f(const Vector& angles, Real& volume) const override {
        Vec3P a; a[0] = P(angles[0]); a[1] = P(angles[1]); a[2] = P(angles[2]);
        Rotation_<P> R_B0B(BodyRotationSequence,
                           a[0], XAxis, a[1], YAxis, a[2], ZAxis);
        int least[3], most[3];
        Vec<3,P> low_B, high_B;
        Geo::Point_<P>::findOrientedExtremePointsIndirect
                                (*p, R_FB0*R_B0B, least, most, low_B, high_B);
        const Vec<3,P> extent(high_B - low_B);
        volume = (Real)(extent[0]*extent[1]*extent[2]);
        return 0;
    }
private:
    const Array_<const Vec3P*>* p; // points to array of points
    const Rotation_<P>          R_FB0; // starting orientation
};


template <class P> /*static*/
Geo::OrientedBox_<P> Geo::Point_<P>::
calcOrientedBoundingBoxIndirect(const Array_<const Vec3P*>& points_F,
                                Array_<int>&                support,
                                bool                        optimize) {
    TransformP X_FB0;
    //TODO: this is not a good initial guess because it is sensitive to
    // point clustering and distribution of interior points. Should start
    // with a convex hull, get the mass properties of that and use principal
    // moments as directions.
    calcPrincipalComponentsIndirect(points_F, X_FB0);

    // We'll update these as we go.
    RotationP R_FB = X_FB0.R();
    Vec3P center_F = X_FB0.p();

    int least[3], most[3];
    Vec3P low_B, high_B;
    findOrientedExtremePointsIndirect(points_F, R_FB,
                                      least, most, low_B, high_B);

    // Initial guess at OBB.
    Vec3P extent_B = high_B - low_B;
    RealP volume = extent_B[0]*extent_B[1]*extent_B[2];
    center_F = R_FB*(high_B+low_B)/2;

    // This is a very abbreviated steepest-descent optimizer. Starting with
    // R_FB0 it uses three Euler angles a as parameters for an incremental
    // rotation R_B0B(a) and calculates volume of the box whose orientation
    // is R_FB(a)=R_FB0*R_B0B(a). It calculates the downhill gradient once, 
    // then does a line search along that gradient and quits as soon as it 
    // stops making significant progress.
    if (optimize) {
        VolumeGradient<P> grad(points_F, R_FB);
        Differentiator diff(grad);
        Vector g;
        diff.calcGradient(Vector(3, Real(0)), (Real)volume, g);
        Vec3P dir; dir[0]=P(g[0]); dir[1]=P(g[1]); dir[2]=P(g[2]);

        // Gradient has units of volume/radian.
        // Set initial step to attempt a 10% volume reduction.
        RealP dirNorm = dir.norm();
        RealP incr = volume/(10*dirNorm);
        const RealP MinImprovement = RealP(.001); // .1% or give up
        RealP minIncr = incr / 1000000;
        RealP step = 0;
        for (int i=0; i < 20; ++i) {
            step -= incr;
            Vec3P a(step*dir); // angles away from B0
            Rotation_<P> R_B0B(BodyRotationSequence,
                               a[0], XAxis, a[1], YAxis, a[2], ZAxis);
            Rotation_<P> tryR_FB = X_FB0.R()*R_B0B;
            int tryLeast[3], tryMost[3];
            Vec3P tryLow_B, tryHigh_B;
            findOrientedExtremePointsIndirect
                   (points_F, tryR_FB, tryLeast, tryMost, tryLow_B, tryHigh_B);
            Vec3P tryExtent_B = tryHigh_B - tryLow_B;
            RealP tryVol = tryExtent_B[0]*tryExtent_B[1]*tryExtent_B[2];

            if (tryVol < volume) {
                const RealP improvement = (volume-tryVol)/volume;
                for (int j=0; j<3; ++j) 
                    least[j]=tryLeast[j], most[j]=tryMost[j];
                R_FB = tryR_FB;
                extent_B = tryExtent_B;
                center_F = R_FB*(tryHigh_B+tryLow_B)/2;
                volume = tryVol;
                if (improvement < MinImprovement)
                    break;           
                incr *= RealP(1.5); // grow slowly
                continue;
            } 

            // Volume got worse.
            step += incr; // back to previous best 
            if (incr <= minIncr) 
                break;
            incr /= 10; // shrink fast
        }
    }

    // Sort the support points and eliminate duplicates.
    std::set<int> supportSet;
    for (unsigned i=0; i<3; ++i) {
        supportSet.insert(least[i]);
        supportSet.insert(most[i]);
    }
    support.assign(supportSet.begin(), supportSet.end());
    
    OrientedBox_<P> obb(TransformP(R_FB,center_F), extent_B/2);
    return obb.stretchBoundary();
}

// This signature takes an array of points rather than an array of pointers
// to points. Allocate a temporary array of pointers and then call the other
// signature.
template <class P> /*static*/
Geo::OrientedBox_<P> Geo::Point_<P>::
calcOrientedBoundingBox(const Array_<Vec3P>& points,
                        Array_<int>&         support,
                        bool                 optimize) {
    Array_<const Vec3P*> indirect;
    makeIndirect(points, indirect);
    return calcOrientedBoundingBoxIndirect(indirect,support,optimize);
}



//==============================================================================
//                       CALC BOUNDING SPHERE - 2 POINTS
//==============================================================================
/* I bet you think this is easy! Unfortunately if the points in question are
far from the origin (at 100, say) but close together, we can't compute the 
separation very accurately (around 100*eps for this example). So the resulting
sphere can be too big (annoying but OK) or too small (very bad). Also, if the
line length is <= tol, we want to consider these a single point rather than
two supporting points for the sphere; in that case we pick one of the points
and center the sphere around that with radius tol. About 45 flops. */
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcBoundingSphere(const Vec3P& p0, const Vec3P& p1, Array_<int>& which) {
    const RealP tol = getDefaultTol<P>();

    // Choose the tentative center, subject to scaled roundoff.
    const Vec3P ctr = (p0 + p1)/2; // 6 flops

    // Rather than using (p0-p1)/2 as the radius which can result in a sphere
    // that doesn't include the points by a large margin (I tried it), measure 
    // the actual distances from the center to each point and use the larger 
    // as the radius. If that radius is <= tol/2 (meaning line length <= tol) 
    // then we'll just move the center to the closer point, set the radius to 
    // tol or more (depending on expected roundoff) and return with 1 support 
    // point (i.e., we consider this a point not a line).
    const RealP p0d2 = (p0-ctr).normSqr();  // squared dist from center
    const RealP p1d2 = (p1-ctr).normSqr();  // (8 flops each)

    RealP rad;
    which.clear();
    if (p0d2 >= p1d2) {
        rad = std::sqrt(p0d2);          // 20 flops
        if (rad <= tol/2) {             //  2 flops
            which.push_back(1);
            return Sphere_<P>(p1, 0).stretchBoundary(); 
        }
    } else { // p1d2 > p0d2
        rad = std::sqrt(p1d2);          // same cost here
        if (rad <= tol/2) {
            which.push_back(0);
            return Sphere_<P>(p0, 0).stretchBoundary(); 
        }
    }
    // use both points
    which.push_back(0); which.push_back(1);
    return Sphere_<P>(ctr, rad).stretchBoundary(); 
}



//==============================================================================
//                       CALC BOUNDING SPHERE - 3 POINTS
//==============================================================================
/* (Note that the 3- and 4-point methods from Nicolas Capens compute only the
sphere that passes through all points (the circumsphere), which can be *way* 
too big.)

This method produces the minimum bounding sphere around three points (i.e.,
a triangle) for all cases except some very, very small triangles where it will
give a small, but perhaps not minimal, bounding sphere. It is 
modified from Christer Ericson's blog: "Minimum bounding circle (sphere) for 
a triangle (tetrahedron)" July 27, 2007. 
http://realtimecollisiondetection.net/blog/?p=20
Accessed 12/12/2011 and reimplemented by Sherm. The main change is handling of
singular triangles and roundoff problems.

Cost is about 110 flops for a typical case.

Implementation
--------------
We have a triangle with vertices a,b,c with no restrictions on their placement 
(i.e., they can be singular in several different ways). Algorithm:
    - Calculate the area of the triangle.
    - If the area is very small:
         - Find the longest edge from ab, ac, bc.
         - Drop the vertex that is not used in the longest edge.
         - Calculate the bounding sphere of the longest edge.
         - If the dropped vertex is outside, stretch the sphere to include it.
         - Return that sphere as the result (might not be perfect but will
           be very good).
    - Else (area is OK):
         - Calculate barycentric coordinates s,t, u=1-s-t of the 
           circumsphere center P (the point equidistant to all vertices).
         - If s,t, or u < 0, P lies outside the triangle and a 2-point sphere 
           will be smaller. Drop the corresponding vertex and return the 
           bounding sphere for the remaining edge (2 supporting points).
         - Otherwise return the circumsphere (3 supporting points).
         - There is an option to force use of the circumsphere regardless of
           where P ends up; sometimes that is needed by the 4-point method.

How to calculate the circumsphere:
    We express center P as P=a + s*(b-a) + t*(c-a). This point is to be 
    equidistant from all vertices, giving us two equations for s and t:
           (P-b)^2 = (P-a)^2, (P-c)^2 = (P-a)^2
    Substituting the definition of P gives us the following system of equations
    for s and t:
           [ab^2  abac] [s]   [ab^2/2]
           [abac  ac^2] [t] = [ac^2/2]   or   M x = b
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of three matrices: M, Ms, Mt where the latter are M with its 
    1st or 2nd column replaced by b. Let dm=det(M), ds=det(Ms), dt=det(Mt). 
    Then s=ds/dm, t=dt/dm. Determinant dm is related to the area of the 
    triangle by 2*a = sqrt(dm), so dm = 4 a^2. Note that dm can't be negative
    except possibly due to roundoff error. Below we'll actually calculate
    2*dm, 2*ds, etc. since that saves some multiplies.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcBoundingSphere(const Vec3P& a, const Vec3P& b, const Vec3P& c,
                   bool forceCircumsphere, Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();

    const Vec3P ab = b - a, ac = c - a;                     //  6 flops
    const RealP ab2 = ab.normSqr(), ac2 = ac.normSqr();     // 10 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta     5 flops
    const RealP ab2ac2 = ab2*ac2;                           //  1 flop

    // The expression below is equivalent to 
    //      dm2 = 2*dm = 2 * (|ab||ac|)^2 * (1-cos^2 theta)
    // where theta is the angle between side ab and ac. That quantity can't
    // be negative except for roundoff error. Note that the area of this 
    // triangle is area = 1/2 sqrt(ab2*ac2-square(abac)) so dm2 = 8*area^2
    const RealP dm2 = 2*(ab2ac2 - square(abac));             //  3 flops

    // We'll use one of three ways to find the center point and determine
    // the set of supporting vertices (in which):
    //  1) singular (2 or 1 support vertices)
    //  2) non-singular with 2 support vertices
    //  3) non-singular using circumsphere (3 support vertices)

    Vec3P ctr;

    // Consider the triangle singular if its area is less than tol.
    if (dm2 <= 8*square(tol)) {                               // 3 flops
        // Triangle is near singular or very small.
        const Vec3P bc = c - b;
        const RealP bc2 = bc.normSqr();
        // Find the longest edge.
        RealP maxLen2; int maxEdge;
        maxOf(ab2,ac2,bc2,maxLen2,maxEdge);
        // Now create a sphere around that edge using 2 or 1 support points,
        // but we're only going to keep the center.
        Sphere_<P> edgeSphere; int map[2];
        switch(maxEdge) {
        case 0: map[0]=0; map[1]=1; // (a,b) drop c
            edgeSphere=calcBoundingSphere(a,b,which); break;
        case 1: map[0]=0; map[1]=2; // (a,c) drop b
            edgeSphere=calcBoundingSphere(a,c,which); break;
        case 2: map[0]=1; map[1]=2; // (b,c) drop a
            edgeSphere=calcBoundingSphere(b,c,which); break;
        };
        fixWhich(map, which); // fix the indices
        ctr = edgeSphere.getCenter();
    } else {
        // Triangle is non-singular. It is still possible that we won't need
        // all 3 points to be on the sphere surface. If one or more of the
        // barycentric coordinates is negative, it means that the circumsphere
        // center is outside the triangle. Pick the most negative, then drop
        // the opposite vertex. The circumsphere around the remaining edge will 
        // include the dropped one for free. 

        // s controls P's height over ac, t over ab, u=1-s-t over bc.
        const RealP ds2 = ab2ac2 - ac2*abac;                 // 2 flops
        const RealP dt2 = ab2ac2 - ab2*abac;                 // 2 flops
        const RealP du2 = dm2-ds2-dt2;                       // 2 flops
        const RealP oodm2 = 1/dm2; // (> 0)                   ~10 flops
        const RealP s = ds2*oodm2, t = dt2*oodm2, u = du2*oodm2; //  3 flops

        RealP minBary; int minBaryIx;
        minOf(s,t,u,minBary,minBaryIx); // 2 flops

        if (!forceCircumsphere && minBary <= 0) { // 1 flop
            Sphere_<P> edgeSphere; int map[2];
            switch(minBaryIx) {
            case 0: // s is the most negative
                map[0]=0; map[1]=2; // sphere around ac includes b
                edgeSphere=calcBoundingSphere(a,c,which); 
                break;
            case 1: // t is the most negative
                map[0]=0; map[1]=1; // sphere around ab includes c
                edgeSphere=calcBoundingSphere(a,b,which); 
                break;
            case 2: // u=1-s-t is the most negative
                map[0]=1; map[1]=2; // sphere around bc includes a
                edgeSphere=calcBoundingSphere(b,c,which); 
                break;
            };
            fixWhich(map, which);
            ctr = edgeSphere.getCenter();
        } else { // s,t,u > 0 or we're forced to use circumsphere
            // All barycentric coordinates are positive. The circumsphere's 
            // center will be inside the triangle and thus of minimal size
            // (unless we were forced to use the circumsphere).
            which.clear(); 
            which.push_back(0); which.push_back(1); which.push_back(2);
            const Vec3P aToCenter = s*ab + t*ac;    //   9 flops
            ctr = a + aToCenter;                    //   3 flops
        }
    }

    // All methods lead here. At this point we have the support points
    // and center point chosen but still need to pick the radius.

    // This cleans up ugly roundoff errors that cause the center not to be
    // exactly equidistant from the points. This makes sure that the radius
    // touches the outermost point with the rest inside. This is 
    // expensive but worth it to ensure a trouble-free sphere.
    RealP rmax2; int rmaxIx;
    maxOf((a-ctr).normSqr(),(b-ctr).normSqr(),(c-ctr).normSqr(),
            rmax2, rmaxIx);                          // 27 flops
    const RealP rad = std::sqrt(rmax2);              // 20 flops

    // If the max radius point wasn't in the support set, we have to add it.
    if (   which.size() < 3 
        && std::find(which.begin(), which.end(), rmaxIx) == which.end())
        which.push_back(rmaxIx);

    return Sphere_<P>(ctr, rad).stretchBoundary();
}



//==============================================================================
//                       CALC BOUNDING SPHERE - 4 POINTS
//==============================================================================
/* This is a 4 point bounding sphere calculator based on Crister Ericson's
suggestion on his blog page -- see the 3 point routine above. The derivation
here was done by Sherm since Ericson didn't work out the 4-point case.

Cost is about 190 flops in a typical case where all 4 points are used.

Implementation
--------------
We have a tetrahedron with vertices a,b,c,d with no restrictions on their 
placement (i.e., they can be singular in a bunch of different ways). Algorithm:
    - Calculate the volume of the tetrahedron.
    - If the volume is very small:
         - Find the face of largest area from abc, abd, acd, bcd.
         - Drop the vertex that is not used in the largest face.
         - Calculate the bounding sphere of the largest face.
         - If the dropped vertex is outside, stretch the sphere to include it.
         - Return that sphere as the result (might not be perfect but will
           be very good).
    - Else (volume is OK):
         - Calculate barycentric coordinates s,t,u, v=1-s-t-u of the 
           circumsphere center P (the point equidistant to all vertices).
         - If s,t,u, or v < 0, P lies outside the tetrahedron and a 3-point 
           sphere will be smaller. Drop the corresponding vertex and return the
           bounding sphere for the remaining face. (3 support points).
         - Otherwise return the circumsphere (4 support points).
         - There is an option to force use of the circumsphere regardless of
           where P ends up.

How to calculate the circumsphere:
    We express center P as P=a + s*(b-a) + t*(c-a) + u*(d-a). This point is to 
    be equidistant from all vertices, giving us 3 equations for s, t, and u:
           (P-b)^2 = (P-a)^2, (P-c)^2 = (P-a)^2, (P-d)^2 = (P-a)^2
    Substituting the definition of P gives us the following system of equations
    for s, t, and u:
         [ab^2  abac  abad] [s]   [ab^2/2]
         [abac  ac^2  acad] [t] = [ac^2/2]   or   M x = b
         [abad  acad  ad^2] [u] = [ad^2/2]
    We'll use Cramer's rule to solve the equation, which requires computing 
    determinants of four matrices: M, Ms, Mt, Mu where the latter are M with 
    its 1st, 2nd, or 3rd column replaced by b. Let dm=det(M), ds=det(Ms), 
    dt=det(Mt), du=det(Mu). Then s=ds/dm, t=dt/dm, u=du/dm. There are going to 
    be lots of common subexpressions, and we can pick up face areas along the 
    way. Also, determinant dm is related to the volume of the tetrahedron by 
    6*v = sqrt(dm), so dm = 36 v^2. Note that dm can't be negative except
    possibly due to roundoff error. Below we'll work with 2*dm, 2*ds, etc.
    because it saves some multiplies.
*/
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcBoundingSphere(const Vec3P& a, const Vec3P& b, 
                   const Vec3P& c, const Vec3P& d,
                   bool forceCircumsphere, Array_<int>& which) {
    const RealP tol = Geo::getDefaultTol<P>();

    const Vec3P ab = b-a, ac = c-a, ad = d-a;               //  9 flops
    const RealP abac = dot(ab, ac); // == |ab||ac|cos theta   ( 5 flops)
    const RealP abad = dot(ab, ad); // == |ab||ad|cos theta   ( 5 flops)
    const RealP acad = dot(ac, ad); // == |ac||ad|cos theta   ( 5 flops)
    // 25 flops in the next two lines.
    const RealP ab2=ab.normSqr(), ac2=ac.normSqr(), ad2=ad.normSqr(); 
    const RealP abac2=square(abac), abad2=square(abad), acad2=square(acad);
    const RealP ab2ac2=ab2*ac2, ab2ad2=ab2*ad2, ac2ad2=ac2*ad2;
    const RealP ab2ac2ad2=ab2ac2*ad2;
   
    const RealP dm2 = 2*(  ab2ac2ad2          // 10 flops
                         - ab2*acad2
                         - ad2*abac2
                         - ac2*abad2
                         + 2*abac*abad*acad); // 72 v^2

    // We'll use one of three ways to find the center point and determine
    // the set of supporting vertices (in which):
    //  1) singular (3, 2, or 1 support vertices)
    //  2) non-singular with 3 support vertices
    //  3) non-singular using circumsphere (4 support vertices)

    Vec3P ctr;

    // Consider the tetrahedron singular if its volume is less than tol.
    if (dm2 <= 72*square(tol)) {                             // 3 flops
        // Tetrahedron is near singular or very small.
        const Vec3P bc = c-b, cd = d-c;
        const RealP bccd = dot(bc, cd), bccd2 = square(bccd);
        const RealP bc2=bc.normSqr(), cd2=cd.normSqr(), bc2cd2=bc2*cd2;
        // Find the largest face.The expressions below are equivalent to, for 
        // example: (|ab||ac|)^2 * (1-cos^2 theta)
        // where theta is the angle between side ab and ac. Note that the area
        // of this triangle is area = 1/2 sqrt(ab2*ac2-square(abac)).
        const RealP abc4Area2 = ab2ac2-abac2; // 4*area(abc)^2
        const RealP abd4Area2 = ab2ad2-abad2; // 4*area(abd)^2
        const RealP acd4Area2 = ac2ad2-acad2; // 4*area(acd)^2
        const RealP bcd4Area2 = bc2cd2-bccd2; // 4*area(bcd)^2
        RealP maxArea2; int maxFace=0;
        maxOf(abc4Area2,abd4Area2,acd4Area2,bcd4Area2,maxArea2,maxFace);
        // Now create a sphere around that edge using 3,2, or 1 support points,
        // but we're only going to keep the center.
        Sphere_<P> faceSphere; int map[3];
        switch(maxFace) {
        case 0: map[0]=0; map[1]=1; map[2]=2; //abc, drop d
            faceSphere=calcBoundingSphere(a,b,c,false,which); break;
        case 1: map[0]=0; map[1]=1; map[2]=3; //abd, drop c
            faceSphere=calcBoundingSphere(a,b,d,false,which); break;
        case 2: map[0]=0; map[1]=2; map[2]=3; //acd, drop b
            faceSphere=calcBoundingSphere(a,c,d,false,which); break;
        case 3: map[0]=1; map[1]=2; map[2]=3; //bcd, drop a
            faceSphere=calcBoundingSphere(b,c,d,false,which); break;
        };
        fixWhich(map, which); // fix the indices
        ctr = faceSphere.getCenter();
    } else {
        // Tetrahedron is non-singular. It is still possible that we won't need
        // all four points to be on the sphere surface. If one or more of the
        // barycentric coordinates is negative, it means that the circumsphere
        // center is outside the tetrahedron. Pick the most negative, then drop
        // the opposite vertex. The circumsphere around the remaining face will 
        // include the dropped one for free. 

        // s controls height over acd, t over abd, u over abc, 1-s-t-u over bcd.
        const RealP ds2 =  ab2ac2ad2            // 12 flops
                         + ac2*abad*acad
                         + ad2*abac*acad
                         - ab2*acad2
                         - ac2ad2*abac
                         - ac2ad2*abad;
        const RealP dt2 =  ab2ac2ad2            // 12 flops
                         + ab2*abad*acad
                         + ad2*abac*abad
                         - ac2*abad2
                         - ab2ad2*acad
                         - ab2ad2*abac;
        const RealP du2 =  ab2ac2ad2            // 12 flops
                         + ab2*abac*acad
                         + ac2*abac*abad
                         - ad2*abac2
                         - ab2ac2*acad
                         - ab2ac2*abad;    
        const RealP dv2 = dm2-ds2-dt2-du2;      //  3 flops  
        const RealP oodm2 = 1/dm2; // (> 0)       ~10 flops
        const RealP s=ds2*oodm2, t=dt2*oodm2, u=du2*oodm2, v=dv2*oodm2; 
                                                //  4 flops

        RealP minBary; int minBaryIx;
        minOf(s,t,u,v,minBary,minBaryIx);   // 3 flops

        if (!forceCircumsphere && minBary <= 0) { // 1 flop
            Sphere_<P> faceSphere; int map[3];
            switch(minBaryIx) {
            case 0: // s is the most negative
                map[0]=0; map[1]=2; map[2]=3; // sphere around acd includes b
                faceSphere=calcBoundingSphere(a,c,d,false,which);
                if (which.size()<3 && faceSphere.isPointOutside(b))
                    faceSphere=calcBoundingSphere(a,c,d,true,which);
                break;
            case 1: // t is the most negative
                map[0]=0; map[1]=1; map[2]=3; // sphere around abd includes c
                faceSphere=calcBoundingSphere(a,b,d,false,which); 
                if (which.size()<3 && faceSphere.isPointOutside(c))
                    faceSphere=calcBoundingSphere(a,b,d,true,which); 
                break;
            case 2: // u is the most negative
                map[0]=0; map[1]=1; map[2]=2; // sphere around abc includes d
                faceSphere=calcBoundingSphere(a,b,c,false,which); 
                if (which.size()<3 && faceSphere.isPointOutside(d))
                    faceSphere=calcBoundingSphere(a,b,c,true,which); 
                break;
            case 3: // v is the most negative
                map[0]=1; map[1]=2; map[2]=3; // sphere around bcd includes a
                faceSphere=calcBoundingSphere(b,c,d,false,which); 
                if (which.size()<3 && faceSphere.isPointOutside(d))
                    faceSphere=calcBoundingSphere(b,c,d,true,which); 
                break;
            };
            fixWhich(map, which);
            ctr = faceSphere.getCenter();
        } else { // s,t,u,v > 0 or we're forced to use circumsphere
            // All barycentric coordinates are positive. The circumsphere's 
            // center will be inside the tetrahedron and thus of minimal size
            // (unless we were forced to use the circumsphere).
            which.clear(); 
            which.push_back(0); which.push_back(1);  
            which.push_back(2); which.push_back(3);
            const Vec3P aToCenter = s*ab + t*ac + u*ad;     //  12 flops
            ctr = a + aToCenter;                            //   3 flops
        }
    }

    // All methods lead here. At this point we have the support points
    // and center point chosen but still need to pick the radius.

    // This cleans up ugly roundoff errors that cause the center not to be
    // exactly equidistant from the points. This makes sure that the radius
    // touches the outermost point with the rest inside. This is 
    // expensive but worth it to ensure a trouble-free sphere.
    RealP rmax2; int rmaxIx;
    maxOf((a-ctr).normSqr(),(b-ctr).normSqr(),
          (c-ctr).normSqr(),(d-ctr).normSqr(),
          rmax2, rmaxIx);                            // 35 flops
    const RealP rad = std::sqrt(rmax2);              // 20 flops

    // If the max radius point wasn't in the support set, we have to add it.
    if (   which.size() < 4 
        && std::find(which.begin(), which.end(), rmaxIx) == which.end())
        which.push_back(rmaxIx);

    return Sphere_<P>(ctr, rad).stretchBoundary();   //  6 flops
}



//==============================================================================
//                       CALC BOUNDING SPHERE - N POINTS
//==============================================================================
// This is called recursively to calculate the minimum bounding sphere for a
// set of points (like the vertices of a mesh). It uses an algorithm developed 
// by Emo Welzl which, despite appearances, has O(n) expected running time.
// We tried the implementation by Nicolas Capens at 
// http://www.flipcode.com/archives/Smallest_Enclosing_Spheres.shtml.
// As described there, the algorithm is highly susceptible to numerical 
// instabilities. Bernd Gartner describes an improved version in "Fast and 
// robust smallest enclosing balls", Proc. 7th Annual ACM European Symposium on 
// Algorithms, v. 1643 Lecture Notes in Computer Science, pp. 325-338, 1999.
// The implementation here follows Gartner's "pivoting" method which uses
// Welzl's recursive algorithm only for maximum sets of 5 points. But here the 
// primitives have been reworked so that they deal nicely with singularities 
// and roundoff. TODO: there are still rare problems that crop up in random-
// point tests. With about 5 million random sets of 1000 points (see TestGeo)
// I found a case that went into an infinite loop (now falls back to Ritter
// sphere after enough attempts). I have also seen the
// minimal sphere come out considerably larger than the crude Ritter sphere
// generated by calcApproxBoundingSphere(); this is also very rare but bears
// investigating. On average I see the Welzl sphere's volume about 20% smaller
// than Ritter's. (sherm 20111227)

// It is possible that we didn't need as support points all the bIn points
// we were given, but instead used a smaller number, bActual. In that case 
// we have to make sure that the *first* bActual are  the support points. 
// So we'll look at the entries [bActual..bIn-1]; if one is now part of the
// support set we'll swap it with a now-unused point in [0..bActual-1].
static void moveSupportToFront(int bIn, const Array_<int>& support,
                               Array_<int>& ix) {
    const int bActual = (int)support.size();
    for (int i=bActual; i < bIn; ++i) {
        const int* s = std::find(support.begin(), support.end(), ix[i]);
        if (s == support.end()) continue; // not being used
        // ix[i] is part of the support set
        bool swapped=false; // There has to be one to swap with!
        for (int j=0; j < bActual; ++j) {
            const int* ss = std::find(support.begin(), support.end(), ix[j]);
            if (ss == support.end()) {
                // ix[j] is not part of the support set; swap with i
                std::swap(ix[i], ix[j]); swapped=true;
                break;
            }
        }
        assert(swapped); // can't happen
    }
}

template <class P> static
Geo::Sphere_<P>
findWelzlSphere(const Array_<const Vec<3,P>*>& p, Array_<int>& ix,
                int bIn, Array_<int>& which, int recursionLevel) {

    // The first bIn points should be an independent support set, and we're 
    // hoping to calculate their circumsphere here. Although in theory the 
    // algorithm wouldn't have gotten here if all bIn points weren't 
    // independent, roundoff issues may cause the underlying primitive to 
    // disagree. So we might use only a subset of the available points, and
    // "which" will list the ones we used. If necessary we'll rearrange the
    // first few entries in "ix" to make sure that the indices of the new
    // support set come first.
    
    Geo::Sphere_<P> minSphere;

    switch(bIn) {
    // Create a bounding sphere for 0, 1, 2, 3, or 4 points. Note which
    // points were actually used, and how many. The original algorithm returned
    // when four points are used, but we still have to check that we have the
    // *right* four points because our primitives might not have used all
    // the points they were given earlier and we need to make sure the new
    // 4-point circumsphere includes the points that were dropped.
    case 0: 
        minSphere = Geo::Sphere_<P>(Vec<3,P>(NTraits<P>::getInfinity()),0);
        which.clear();
        break;
    case 1:
        minSphere = Geo::Point_<P>::calcBoundingSphere(*p[ix[0]],which);
        break;
    case 2:
        minSphere = Geo::Point_<P>::calcBoundingSphere
                            (*p[ix[0]],*p[ix[1]],which);
        break;
    case 3:
        minSphere = Geo::Point_<P>::calcBoundingSphere
                        (*p[ix[0]],*p[ix[1]],*p[ix[2]],false,which);
        break;
    case 4:
        minSphere = Geo::Point_<P>::calcBoundingSphere
                        (*p[ix[0]],*p[ix[1]],*p[ix[2]],*p[ix[3]],false,which);
        break;
    }

    // We have a sphere that surrounds all bIn points but might not have
    // used them all as support in "which".
    
    // Fix the indices in which (the primitives number from 0).
    fixWhich(ix.begin(), which);
    // Make sure the support set is at the front of ix.
    moveSupportToFront(bIn, which, ix);

    if (bIn==4) {
        // There can only be 4 or 5 points. If 4 we're done.
        if (ix.size()==4 || !minSphere.isPointOutside(*p[ix[4]])) 
            return minSphere; // all the points already enclosed
        // We failed with the current "left out" point in ix[4]. We'll try
        // the other four possibilities until one works.
        for (int i=0; i<4; ++i) {
            for (int j = 4; j > 0; --j) // rotate last point in
                std::swap(ix[j], ix[j-1]);
            minSphere = Geo::Point_<P>::calcBoundingSphere
                        (*p[ix[0]],*p[ix[1]],*p[ix[2]],*p[ix[3]],false,which);
            fixWhich(ix.begin(), which);
            if (!minSphere.isPointOutside(*p[ix[4]])) {
                moveSupportToFront(bIn, which, ix);
                return minSphere;
            }
        }
        // If we get here then no sphere around 4 points included the 5th.
        // That should be impossible and I've never seen it happen. But ...
        // I can't prove it would never happen given finite arithmetic so
        // just in case we'll generate a Ritter sphere, which will always work.
        Array_<const Vec<3,P>*> ritterPoints;
        for (int i=0; i < (int)ix.size(); ++i)
            ritterPoints.push_back(p[ix[i]]);
        minSphere = 
            Geo::Point_<P>::calcApproxBoundingSphereIndirect(ritterPoints);
        which.clear(); // we don't know
        return minSphere;
    }


    // The indices of the support points are within the first bIn entries in 
    // ix, although they may not be in the first bActual slots. The first bIn
    // may therefore contain some unused points we already processed
    // but didn't need. Now run through all subsequent points and update the
    // sphere to include them if necessary, each time ensuring that all 
    // previous points remain included.
  
    for (int i = bIn; i < (int)ix.size(); ++i) {
        // We expect the sphere already to have been enlarged to deal with
        // roundoff so we can do an exact test here.
        if (minSphere.isPointOutside(*p[ix[i]])) {
            // This point is outside the current bounding sphere.  
            // Move it to the start of the list *without* reordering.
            for (int j = i; j > 0; --j)
                std::swap(ix[j], ix[j-1]);
            
            // Update the bounding sphere, taking the new point into account
            // and ensuring that the resulting sphere also includes all the
            // previous points as well. No more than 4 support points are ever
            // needed.
            ArrayView_<int> toBoundIx(ix.begin(), &ix[i]+1);
            minSphere = findWelzlSphere<P>(p, toBoundIx, 
                                           std::min(bIn+1,4), 
                                           which, recursionLevel+1);

        }
    }

    return minSphere; // Already stretched for roundoff.
}

// This signature takes an array of points, creates an array of pointers to 
// those points and calls the other signature.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcBoundingSphere(const Array_<Vec3P>& points, Array_<int>& which) {
    Array_<const Vec3P*> indirect;
    makeIndirect(points, indirect);
    return calcBoundingSphereIndirect(indirect, which);
}

// This was the outer block for the basic Emo Welzl "move to front" algorithm
// (algorithm 1 in Gartner's paper). I found it unreliable and sometimes slow
// in practice and switched to Gartner's "pivoting" algorithm below.

//template <class P> /*static*/
//Geo::Sphere_<P> Geo::Point_<P>::
//calcBoundingSphere(const Array_<const Vec3P*>& points, Array_<int>& which) {
//    const unsigned npoints = points.size();
//
//    // Allocate and initialize an array of point indices. These will get 
//    // moved around during the computation.
//    Array_<int> ix(npoints); for (unsigned i=0; i<npoints; ++i) ix[i] = i; 
//
//    if (npoints < 10) {
//        // Not worth rearranging.
//        return findWelzlSphere<P>(points, ix, 1, which, 0);
//    }
//
//    // There are enough points that we'll try to improve the ordering so that
//    // the bounding sphere gets large quickly. This optimization helps *a lot* 
//    // for large numbers of points.
//
//    // Find the six points that have the most extreme
//    // x,y, and z coordinates (not necessarily six unique points) and move
//    // them to the front so they get processed first.
//    Vec3P lo=*points[0], hi=*points[0]; // initialize extremes
//    int   ilo[3], ihi[3]; for (int i=0; i<3; ++i) ilo[i]=ihi[i]=0;
//    for (unsigned i=0; i<points.size(); ++i) {
//        const Vec3P& p = *points[i];
//        if (p[0] > hi[0]) hi[0]=p[0], ihi[0]=i;
//        if (p[0] < lo[0]) lo[0]=p[0], ilo[0]=i;
//        if (p[1] > hi[1]) hi[1]=p[1], ihi[1]=i;
//        if (p[1] < lo[1]) lo[1]=p[1], ilo[1]=i;
//        if (p[2] > hi[2]) hi[2]=p[2], ihi[2]=i;
//        if (p[2] < lo[2]) lo[2]=p[2], ilo[2]=i;
//    }
//    // Find the nx <= 6 unique extreme points.
//    std::set<int> pending;
//    pending.insert(ilo, ilo+3); pending.insert(ihi, ihi+3);
//    const int nx = pending.size();
//    // Go through the first n points. If the point is already an extreme,
//    // remove it from the pending list. If not, swap it with one of the
//    // extreme points.
//    for (int i=0; i < nx; ++i) {
//        std::set<int>::iterator p = pending.find(ix[i]);
//        if (p != pending.end()) {
//            pending.erase(p);
//            continue;
//        }
//        p = pending.begin(); // first unmoved extreme
//        const int extremeIx = *p;
//        pending.erase(p);
//        std::swap(ix[i], ix[extremeIx]);
//    }
//
//    return findWelzlSphere<P>(points, ix, 1, which, 0);
//}

template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcBoundingSphereIndirect(const Array_<const Vec3P*>& points, 
                           Array_<int>& which) {
    const unsigned npoints = points.size();
    if (npoints == 0)
        return Geo::Sphere_<P>(Vec<3,P>(NTraits<P>::getInfinity()),0);

    // Allocate and initialize an array of point indices. These will get 
    // moved around during the computation.
    Array_<int> ix(npoints); for (unsigned i=0; i<npoints; ++i) ix[i] = i; 

    // Sherm 20111225: there are still very obscure cases where this 
    // code may get stuck in a cycle and iterate forever. It usually only
    // takes a few iterations to finish, so 100 means something is very wrong.
    // In that case we'll give up and return a Ritter sphere.
    const int MaxIters = 100;

    unsigned nxt=0; // "t" in Gartner's paper
    Geo::Sphere_<P> minSphere = calcBoundingSphere(*points[ix[nxt++]], which);

    // Loop until all the points are inside minSphere.
    int s=1; // number of support points
    for(int iters=0; ; ++iters) {
        // Find worst point.
        const Vec3P& center = minSphere.getCenter();
        RealP maxDist2=0; unsigned worst=0;
        for (unsigned i=nxt; i < ix.size(); ++i) {
            const RealP d2 = (*points[ix[i]] - center).normSqr();
            if (d2 > maxDist2) maxDist2=d2, worst=i;
        }
        if (maxDist2 <= square(minSphere.getRadius()))
            break;

        if (iters == MaxIters) {
            minSphere = calcApproxBoundingSphereIndirect(points);
            which.clear(); // we don't know
            break;
        }

        // Point worst (>=nxt) is outside the current bounding sphere.  
        // Move it to the start of the list *without* reordering.
        for (int j = worst; j > 0; --j)
            std::swap(ix[j], ix[j-1]);
        
        // Make a list of the initial s+1 points, s <= 4.
        ArrayView_<int> toBoundIx(ix.begin(), &ix[s]+1);

        nxt = s+1; // we're going to bound points 0..s
        s=std::min(s+1,4); // number of supports to try for next
        minSphere = findWelzlSphere<P>(points, toBoundIx, s, which, 0);

    }

    return minSphere;
}

// This signature takes an array of points, creates an array of pointers to 
// those points and calls the other signature.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcApproxBoundingSphere(const Array_<Vec3P>& points) {
    Array_<const Vec3P*> indirect;
    makeIndirect(points, indirect);
    return calcApproxBoundingSphereIndirect(indirect);
}

// Calculate a Ritter sphere using the method described by Christer Ericson
// in Real Time Collision Detection, Elsevier 2005, pp. 89-91. This method
// makes one pass to find some reasonably far apart points, then starts with
// a sphere around those two points. Then it makes another pass growing the
// sphere to just include the current sphere plus the new point.
template <class P> /*static*/
Geo::Sphere_<P> Geo::Point_<P>::
calcApproxBoundingSphereIndirect(const Array_<const Vec3P*>& points) {
    const unsigned npoints = points.size();
    if (npoints == 0)
        return Geo::Sphere_<P>(Vec<3,P>(0),0);

    // Find the most-separated pair of points along one of the coordinate axes.
    Vec3P lo=*points[0], hi=*points[0]; // initialize extremes
    int   ilo[3], ihi[3]; for (int i=0; i<3; ++i) ilo[i]=ihi[i]=0;
    for (unsigned i=0; i<points.size(); ++i) {
        const Vec3P& p = *points[i];
        if (p[0] > hi[0]) hi[0]=p[0], ihi[0]=i;
        if (p[0] < lo[0]) lo[0]=p[0], ilo[0]=i;
        if (p[1] > hi[1]) hi[1]=p[1], ihi[1]=i;
        if (p[1] < lo[1]) lo[1]=p[1], ilo[1]=i;
        if (p[2] > hi[2]) hi[2]=p[2], ihi[2]=i;
        if (p[2] < lo[2]) lo[2]=p[2], ilo[2]=i;
    }
    const RealP xdist2 = (*points[ihi[0]]-*points[ilo[0]]).normSqr();
    const RealP ydist2 = (*points[ihi[0]]-*points[ilo[0]]).normSqr();
    const RealP zdist2 = (*points[ihi[0]]-*points[ilo[0]]).normSqr();
    RealP maxVal; int which;
    maxOf(xdist2, ydist2, zdist2, maxVal, which);
    const Vec3P& pmin = *points[ilo[which]];
    const Vec3P& pmax = *points[ihi[which]];

    Geo::Sphere_<P> minSphere;
    Vec3P& ctr = minSphere.updCenter(); // aliases
    RealP& rad = minSphere.updRadius();
    ctr = (pmin+pmax)/2;
    // Calculating radius this way ensures that roundoff won't leave one of
    // the points outside. We'll do a final roundoff adjustment at the end.
    rad = std::sqrt(std::max((pmax-ctr).normSqr(),
                             (pmin-ctr).normSqr()));

    // Now run through all the points again and grow the sphere if necessary. 
    // This requires moving the center too so we don't grow more than needed.
    for (unsigned i=0; i<points.size(); ++i) {
        const Vec3P& p = *points[i];
        const Vec3P ctr2p = p-ctr; 
        const RealP dist2 = ctr2p.normSqr();
        if (dist2 > square(rad)) {
            const RealP dist   = std::sqrt(dist2);
            const RealP newrad = (dist + rad)/2;
            ctr += (newrad - rad)/dist * ctr2p; // has roundoff issues
            rad = newrad;
            // Make sure we didn't miss the new point due to roundoff.
            const RealP newdist2 = (p-ctr).normSqr();
            if (newdist2 > square(rad)) 
                rad = std::sqrt(newdist2);
        }
    }

    return minSphere.stretchBoundary();
}

// Explicit instantiations for float and double.
template class Geo::Point_<float>;
template class Geo::Point_<double>;


}  // End of namespace SimTK
