/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman, Ian Stavness                      *
 * Contributors: Andreas Scholz, Matthew Millard                              *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/Geo_BicubicBezierPatch.h"
#include "simmath/internal/BicubicSurface.h"
#include "simmath/internal/ParticleConSurfaceSystem.h"
#include "simmath/Differentiator.h"
#include "simmath/RungeKutta3Integrator.h"
#include "simmath/RungeKuttaMersonIntegrator.h"
#include "simmath/RungeKuttaFeldbergIntegrator.h"
#include "simmath/VerletIntegrator.h"
#include "simmath/CPodesIntegrator.h"
#include "simmath/TimeStepper.h"
#include "simmath/internal/ContactGeometry.h"
#include "simmath/internal/GeodesicIntegrator.h"

#include "ContactGeometryImpl.h"
#include "GeodesicEquations.h"

#include <iostream>
#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::cout; using std::endl;

#define USE_NEW_INTEGRATOR

//==============================================================================
//                            CONTACT GEOMETRY
//==============================================================================

ContactGeometry::ContactGeometry(ContactGeometryImpl* impl) : impl(impl) {
    assert(impl);
    impl->setMyHandle(*this);
}

ContactGeometry::~ContactGeometry() {
    if (isOwnerHandle())
        delete impl;
    impl = 0;
}

bool ContactGeometry::isOwnerHandle() const {
    return (impl == 0 || impl->getMyHandle() == this);
}

bool ContactGeometry::isEmptyHandle() const {
    return (impl == 0);
}

ContactGeometry::ContactGeometry(const ContactGeometry& src) : impl(0) {
    if (src.impl) {
        impl = src.impl->clone();
        impl->setMyHandle(*this);
    }
}

ContactGeometry& ContactGeometry::operator=(const ContactGeometry& src) {
    if (&src != this) {
        if (isOwnerHandle())
            delete impl;
        impl = 0;
        if (src.impl) {
            impl = src.impl->clone();
            impl->setMyHandle(*this);
        }
    }
    return *this;
}

ContactGeometryTypeId ContactGeometry::
getTypeId() const {return getImpl().getTypeId();}

DecorativeGeometry ContactGeometry::createDecorativeGeometry() const {
    return getImpl().createDecorativeGeometry();
}

bool ContactGeometry::intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

Vec3 ContactGeometry::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
}

Vec3 ContactGeometry::projectDownhillToNearestPoint(const Vec3& Q) const {
    return getImpl().projectDownhillToNearestPoint(Q);
}

bool ContactGeometry::
trackSeparationFromLine(const Vec3& pointOnLine,
                        const UnitVec3& directionOfLine,
                        const Vec3& startingGuessForClosestPoint,
                        Vec3& newClosestPointOnSurface,
                        Vec3& closestPointOnLine,
                        Real& height) const {
     return getImpl().trackSeparationFromLine(pointOnLine,directionOfLine,
                startingGuessForClosestPoint,newClosestPointOnSurface,
                closestPointOnLine,height); 
}


void ContactGeometry::getBoundingSphere(Vec3& center, Real& radius) const {
    getImpl().getBoundingSphere(center, radius);
}

bool ContactGeometry::isSmooth() const {return getImpl().isSmooth();}
bool ContactGeometry::isConvex() const {return getImpl().isConvex();}

void ContactGeometry::calcCurvature(const Vec3& point, Vec2& curvature, 
                                    Rotation& orientation) const 
{   getImpl().calcCurvature(point, curvature, orientation); }

const Function& ContactGeometry::getImplicitFunction() const 
{   return getImpl().getImplicitFunction(); }


Real ContactGeometry::calcSurfaceValue(const Vec3& point) const {
    return getImpl().calcSurfaceValue(point);
}

UnitVec3 ContactGeometry::calcSurfaceUnitNormal(const Vec3& point) const {
    return getImpl().calcSurfaceUnitNormal(point);
}

Vec3 ContactGeometry::calcSurfaceGradient(const Vec3& point) const {
    return getImpl().calcSurfaceGradient(point);
}

Mat33 ContactGeometry::calcSurfaceHessian(const Vec3& point) const {
    return getImpl().calcSurfaceHessian(point);
}

Real ContactGeometry::calcGaussianCurvature(const Vec3& gradient,
                                            const Mat33& Hessian) const {
    return getImpl().calcGaussianCurvature(gradient,Hessian);
}

Real ContactGeometry::calcSurfaceCurvatureInDirection(const Vec3& point, const UnitVec3& direction) const {
    return getImpl().calcSurfaceCurvatureInDirection(point, direction);
}

void ContactGeometry::calcSurfacePrincipalCurvatures(const Vec3& point, 
                                                     Vec2&       k,
                                                     Rotation&   R_SP) const
{
    getImpl().calcSurfacePrincipalCurvatures(point, k, R_SP);
}

Vec3 ContactGeometry::calcSupportPoint(UnitVec3 direction) const 
{   return getImpl().calcSupportPoint(direction); }


//------------------------------------------------------------------------------
//                        EVAL PARAMETRIC CURVATURE
//------------------------------------------------------------------------------
/*static*/Vec2 ContactGeometry::
evalParametricCurvature(const Vec3& P, const UnitVec3& nn,
                        const Vec3& dPdu, const Vec3& dPdv,
                        const Vec3& d2Pdu2, const Vec3& d2Pdv2, 
                        const Vec3& d2Pdudv,
                        Transform& X_EP)
{
    // All this is 42 flops
    Real E =  ~dPdu*dPdu,  F =  ~dPdu*dPdv,   G =  ~dPdv*dPdv;
    Real e =-(~d2Pdu2*nn), f =-(~d2Pdudv*nn), g =-(~d2Pdv2*nn);
    Real A = F*g-G*f, B = E*g-G*e, C = E*f-F*e;

    Real kmax, kmin;
    UnitVec3 dmax;
    if (std::abs(F) < SignificantReal) {
        Real ku = e/E, kv = g/G; // two divides ~20 flops
        if (ku < kv) {
            kmax=kv, kmin=ku;
            dmax=UnitVec3(dPdv); // normalizing, ~35 flops
        } else {
            kmax=ku, kmin=kv;
            dmax=UnitVec3(dPdu); // normalizing, ~35 flops
        }
    } else {
        // ~40 flops
        // t = (-b +/- sqrt(b^2-4ac)) / 2a
        // Discriminant must be nonnegative for real surfaces
        // but could be slightly negative due to numerical noise.
        Real sqrtd = std::sqrt(std::max(B*B - 4*A*C, Real(0)));
        Vec2 t = Vec2(sqrtd - B, -sqrtd - B) / (2*A);

        // Two divides + misc: ~30 flops
        Real kr = (e + f*t[0])/(E+F*t[0]); // Struik, eq. 6-4, pg 80
        Real ks = (e + f*t[1])/(E+F*t[1]); // (works only because these are extremes)
                                           // otherwise use eq. 6-3.

        if (kr < ks) {
            kmax=ks, kmin=kr;
            dmax = UnitVec3(t[1]*dPdv + dPdu); // Sdir, normalizing, ~40 flops
        } else {
            kmax=kr, kmin=ks;
            dmax = UnitVec3(t[0]*dPdv + dPdu); // Rdir, normalizing, ~40 flops
        }
    }

    // y=z%x ensures right handed; already unit vec (9 flops)
    UnitVec3 dmin = UnitVec3(nn % dmax, true);
    X_EP.updR().setRotationFromUnitVecsTrustMe(dmax, dmin, nn);
    X_EP.updP() = P; // the origin point

    return Vec2(kmax, kmin);
}


//------------------------------------------------------------------------------
//                          COMBINE PARABOLOIDS
//------------------------------------------------------------------------------
// See the documentation in the header file for a complete description of
// what's being calculated here. This comment adds implementation information
// that isn't relevant to the API user. 
// 
// Given two paraboloids P1 and P2 sharing a common normal z and origin,
// compute the paraboloid that represents their difference, and express that 
// paraboloid in a frame that has been rotated around z so that x and y 
// coincide with the principal curvature directions. P1 and P2 may be
// elliptic (kmax>=kmin>=0) or hyperbolic (kmax>=0>kmin). If the surfaces are 
// non-conforming, their difference will be elliptic with kmax>=kmin>0. 
//
// We assume the paraboloids represent surfaces and that each has its z axis 
// oriented away from the surface, pointing outside the "bowl" of the elliptic 
// paraboloid or away from the convex direction of a hyperbolic paraboloid. 
// That's the opposite sense from a standard paraboloid parameterization.
// The z axes are antiparallel. We will return the resulting difference 
// paraboloid in a frame whose z axis is coincident with P1's z axis, and thus 
// antiparallel to P2's z axis.
//
//     P1: z = -(kmax1/2 x1^2 + kmin1/2 y1^2)
//     P2: z =   kmax2/2 x2^2 + kmin2/2 y2^2
//      P: z = -( kmax/2  x^2 +  kmin/2  y^2)
// Thus the right-handed coordinate frames are:
//     P1: (x1,y1, z)
//     P2: (x2,y2,-z)
//      P: ( x, y, z)
// The above distinctions don't matter a whole lot for the implementation here,
// but still, I thought you might like to know anyway.
//
// Cost is about 70 flops to get the curvatures kmax,kmin. Then if you want
// the curvature directions too it costs another 150 flops. 

// This local static helper method calculates the curvatures and returns 
// intermediates necessary for calculating the directions, but doesn't actually
// calculate them. So we use only about 70 flops here.
static void combineParaboloidsHelper
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Real& cos2w, Real& sin2w, Real& kdiff1, Real& kdiff2, Vec2& k)
{
    const UnitVec3& x1 = R_SP1.x(); // P1 kmax direction (y is kmin direction)
    const UnitVec3& z  = R_SP1.z(); // P1, P, -P2 normal

    const Real ksum1  = k1[0]+k1[1], ksum2  = k2[0]+k2[1]; // 4 flops
    kdiff1 = k1[0]-k1[1], kdiff2 = k2[0]-k2[1];

    // w is angle between x1, x2 max curvature directions defined
    // using right hand rule rotation of x1 about z until it is
    // coincident with x2. But ... we want -90 <= w <= 90, meaning
    // cos(w) >= 0. If necessary we flip x2 180 degrees around z, 
    // since -x2 is an equally good max curvature direction.
    const Real dotx1x2 = dot(x1,x2);         // 5 flops
    const UnitVec3 x2p = dotx1x2 < 0 ? -x2 : x2;
    const Real cosw = std::abs(dotx1x2);
    const Real sinw = dot(cross(x1,x2p), z); // signed, 14 flops

    // We'll need cos(2w), sin(2w); luckily these are easy to get (5 flops).
    cos2w = 2*square(cosw) - 1; // double angle formulas
    sin2w = 2*sinw*cosw;

    // Compute min/max curvatures of the difference surface.
    // See KL Johnson 1987 Ch. 4 and Appendix 2, and J-F Antoine, et al. 
    // 2006 pg 661. ~35 flops
    const Real ksum = ksum1 + ksum2;
    const Real kdiff = std::sqrt(square(kdiff1) + square(kdiff2)
                                 + 2*kdiff1*kdiff2*cos2w);
    k = Vec2(ksum + kdiff, ksum - kdiff)/2; // kmax, kmin (4 flops)  
}

// This is the full version that calculates the curvatures for 70 flops
// then spends another 150 to get the curvature directions.
/*static*/ void ContactGeometry::combineParaboloids
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Rotation& R_SP, Vec2& k)
{
    Real cos2w, sin2w, kdiff1, kdiff2;
    combineParaboloidsHelper(R_SP1, k1, x2, k2,
                             cos2w, sin2w, kdiff1, kdiff2, k);

    // Now find the rotated coordinate system by solving for the
    // angle -90 <= alpha <= 90 by which we need to rotate the x1 axis
    // about z to align it with the x axis. See KL Johnson Appendix 2
    // again, noting that beta = theta-alpha, then solving for tan(2 alpha).
    // This is about 130 flops.
    const Real yy = kdiff2*sin2w, xx = kdiff2*cos2w + kdiff1;
    Real a = std::atan2(yy,xx) / 2; // yy==xx==0 -> a=0
    Real cosa = std::cos(a), sina = std::sin(a);

    // Perform the actual rotations of x1,y1 to get x,y (18 flops)
    const UnitVec3& x1 = R_SP1.x(); // P1 kmax direction
    const UnitVec3& y1 = R_SP1.y(); // P1 kmin direction
    const UnitVec3& z  = R_SP1.z(); // P1, P, -P2 normal
    R_SP.setRotationFromUnitVecsTrustMe(UnitVec3(cosa*x1 + sina*y1, true),
                                        UnitVec3(cosa*y1 - sina*x1, true),
                                        z);
}

// This is the abridged version that costs only 70 flops but gives you just
// the curvatures without the curvature directions.
/*static*/ void ContactGeometry::combineParaboloids
   (const Rotation& R_SP1, const Vec2& k1,
    const UnitVec3& x2, const Vec2& k2,
    Vec2& k)
{
    Real cos2w, sin2w, kdiff1, kdiff2; // unneeded
    combineParaboloidsHelper(R_SP1, k1, x2, k2,
                             cos2w, sin2w, kdiff1, kdiff2, k);     
}

//------------------------------------------------------------------------------
//                   PROJECT DOWNHILL TO NEAREST POINT
//------------------------------------------------------------------------------

// Find the "local" nearest point to p on this surface. The algorithm assumes that
// p is close to the surface and uses Newton's method to proceed downhill until a
// feasible surface point is found.
//
// The nearest point, x, is a point on the surface whose normal points toward the
// initial point, p. The nearest point, x, must satisfy the following conditions:
// (1) surf(x) = 0
// (2) ~r_px*t = 0
// (3) ~r_px*b = 0
// where surf(x) is the implicit surface function, r_px = (p-x), and
//  t and b are orthogonal vectors that span the tangent plane at x.
//
// Condition (1) ensures that x is on the surface, and conditions (2) and (3)
// ensure that the line between p and x is perpendicular to the surface
// tangent plane at x.
//
// At point p, we choose arbitrary basis vectors tP and bP that are perpendicular
// to the gradient at p. The tangent plane at x is spanned by t and b as:
// t := tP - n*(~tP*n)
// b := bP - n*(~bP*n)
// where n is the unit normal vector at x
//
// Newton's method solves err(x) = 0 iteratively as
//
// f(x + dx) ~= f(x) + J dx = 0
//
// The jacobian J has the following rows:
// J(1) = dsurf/dx = ~g
// J(2) = ~t*Rx + ~r*Tx
// J(3) = ~b*Rx + ~r*Bx
//
// n   := g/norm(g) (unit normal vector)
// H   := dg/dx (the Hessian of the surface at x)
// Rx  := dr/dx = -I
// Nx  := dn/dx = dn/dg*dg/dx = -(I - n*(~n))*H/norm(g) // negative due to sign convention of gradient
// Tx  := dt/dx = -Nx*(~tP*n) -n*(~tP*Nx)
// Bx  := db/dx = -Nx*(~bP*n) -n*(~bP*Nx)

Vec3 ContactGeometryImpl::
projectDownhillToNearestPoint(const Vec3& Q) const {

    // Newton solver settings
    const Real ftol = SignificantReal;

    // Check for immediate return.
    if (std::abs(calcSurfaceValue(Q)) <= ftol)
        return Q;

    // construct arbitrary frame at p
    UnitVec3 nP = calcSurfaceUnitNormal(Q);
    SimTK_ASSERT_ALWAYS(!nP.isNaN(), "gradient is undefined at the query point.");
    UnitVec3 tP = nP.perp();
    UnitVec3 bP(tP%nP);

    // Estimate a scale for the local neighborhood of this surface by 
    // using the larger curvature in the t or b direction. We want to take
    // conservative steps that never move by more than a fraction of the
    // scale to avoid jumping out of the local minimum.
    const Real kt = calcSurfaceCurvatureInDirection(Q, tP);
    const Real kb = calcSurfaceCurvatureInDirection(Q, bP);
    const Real maxK = std::max(std::abs(kt),std::abs(kb));
    const Real scale = std::min(1/maxK, Real(1000)); // keep scale reasonable
    const Real MaxMove = Real(.25); // Limit one move to 25% of smaller radius.

    const Real xtol = Real(1e-12);  // TODO: what should these be in single
    const Real minlam = Real(1e-9); //      precision?
    const int maxNewtonIterations = 30;

    Vec3 x(Q); // initialize to query point
    Vec3 f, dx, xold, g, r, t, b;
    Mat33 J, H, Nx, Tx, Bx;
    Real gNorm, rmsError, rmsErrorOld, xchg, lam;

    // initialized to query point, therefore first step is in gradient direction
    f[0] = calcSurfaceValue(x);
    f[1] = 0;
    f[2] = 0;

    rmsError = std::sqrt(f.normSqr()/3);
    //std::cout << "BEFORE Q=" << Q << ", f=" << f << ", frms=" << rmsError << std::endl;

    int cnt = 0;
    do {
        if (rmsError <= ftol) { // found solution
            //std::cout << "CONVERGED in " << cnt << " steps, frms=" << rmsError << std::endl;
            break;
        }

        g = calcSurfaceGradient(x); // non-normalized gradient vector
        H = calcSurfaceHessian(x);
        gNorm = g.norm();
        UnitVec3 n(-g/gNorm, true);

        // calculate frame at x
        r = Q-x;
        t = tP-n*(~tP*n); // project tP to tangent plane at x
        b = bP-n*(~bP*n); // project bP to tangent plane at x

        SimTK_ASSERT_ALWAYS(t.norm() > Real(1e-6), 
            "t is aligned with the normal vector at the current point.");
        SimTK_ASSERT_ALWAYS(b.norm() > Real(1e-6), 
            "b is aligned with the with normal vector at the current point.");

        // calculate error
        f[0] = calcSurfaceValue(x);
        f[1] = ~r*t;
        f[2] = ~r*b;

        // calculate derivatives
        Nx = -((Mat33(1) - n*(~n))/(gNorm))*H; // negative due to sign convention of gradient
        Tx = -Nx*(~tP*n) - n*(~tP*Nx);
        Bx = -Nx*(~bP*n) - n*(~bP*Nx);

        J[0] = ~g;
        J[1] = -(~t) + ~r*Tx;
        J[2] = -(~b) + ~r*Bx;

        dx = J.invert()*f;
        const Real dxrms = std::sqrt(dx.normSqr()/3);

        rmsErrorOld = rmsError;
        xold = x;

        // Backtracking. Limit the starting step size if dx is too big.
        lam = std::min(Real(1), MaxMove*(scale/dxrms));
        if (lam < 1) {
            //std::cout << "PROJECT: LIMITED STEP: iter=" << cnt 
            //          << " lam=" << lam << endl;
        }
        while (true) {
            x = xold - lam*dx;

            n = calcSurfaceUnitNormal(x);
            r = Q-x;
            t = tP-n*(~tP*n); // project tP to tangent plane at x
            b = bP-n*(~bP*n); // project bP to tangent plane at x

            f[0] = calcSurfaceValue(x);
            f[1] = ~r*t;
            f[2] = ~r*b;

            rmsError = std::sqrt(f.normSqr()/3);
            if (rmsError > rmsErrorOld && lam > minlam) {
                lam = lam / 2;
            } else {
                break;
            }
        }

        //std::cout << cnt << ": AFTER x-=" << lam << "*dx, x=" << x 
        //          << ", f=" << f  << ", frms=" << rmsError << std::endl;

        if (rmsError > ftol) {
            xchg = dxrms*lam; // roughly, how much we changed x
            if (xchg < xtol) { // check step size
                std::cout << "PROJECT: STALLED on step size, xchg=" << xchg 
                          << " frms=" << rmsError << std::endl;
                break;
            }
        }

        cnt++;
        if (cnt > maxNewtonIterations) {
//            SimTK_ASSERT_ALWAYS(false,"Newton solve did not converge, max iterations taken");
            std::cout << "PROJECT: MAX iterations taken" << std::endl;
            break; // Return whatever we got.
        }

        rmsErrorOld = rmsError;

    } while (true);

    return x;
}



//------------------------------------------------------------------------------
//                        TRACK SEPARATION FROM LINE
//------------------------------------------------------------------------------

class TrackLineSeparationJacobian : public Differentiator::JacobianFunction {
public:
    TrackLineSeparationJacobian(const ContactGeometryImpl& geom, 
                                const Vec3& p, const UnitVec3& e)
    :   Differentiator::JacobianFunction(3,3), geom(geom), p(p), e(e) { }

    int f(const Vector& xvec, Vector& f) const override {
        const Vec3& x = Vec3::getAs(&xvec[0]);
        UnitVec3 n; Vec3 closestPointOnLine; // not used
        const Vec3 eps = calcExtremePointError(x, n, closestPointOnLine);

        f[0] = eps[0]; f[1] = eps[1]; f[2] = eps[2];
        return 0;
    }

    // Given a point x, return the value of the "closest point" error of x.
    // We return the outward unit normal n at x, and the line's closest point R
    // as side effects in case you are interested.
    Vec3 calcExtremePointError(const Vec3& x, UnitVec3& n, Vec3& R) const {
        n = geom.calcSurfaceUnitNormal(x);

        Vec3 Q; // Q is the point of the normal line closest to L
        bool linesAreParallel;
        Geo::findClosestPointsOfTwoLines(x, n, p, e, 
            Q, R, linesAreParallel);

        if (linesAreParallel)           
            cout << "findClosest: PARALLEL!!!" << endl;

        Vec3 errf( ~n * e,   // normal and line should be perpendicular
                   ~(R-Q) * (n % e), // no separation along common perpendicular
                   geom.calcSurfaceValue(x) ); // x should be on the surface

        return errf;
    }

private:
    const ContactGeometryImpl& geom;
    const Vec3      p;  // a point on the line
    const UnitVec3  e;  // a unit vector along the line
};

// This method is required to search downhill only from the starting guess.
// The Newton iteration is throttled back accordingly if it gets too
// agressive.
bool ContactGeometryImpl::
trackSeparationFromLine(const Vec3& pointOnLine,
                        const UnitVec3& directionOfLine,
                        const Vec3& startingGuessForClosestPoint,
                        Vec3& x, // the new extreme point
                        Vec3& closestPointOnLine,
                        Real& height) const
{
    // Newton solver settings.

    // RMS error must reach sqrt of this value (squared for speed).
    const Real Ftol2 = square(SignificantReal);

    // Limit the number of Newton steps. We don't
    // count steps that we limited because we were nervous about the size of
    // the change, so this value is *very* generous. It should never take 
    // more than 7 full Newton iterations to solve to machine precision.
    const int MaxNewtonIterations = 20;
    // Make sure there is at least some limit on the total number of steps
    // including ones we throttled back on purpose, just as a sanity check.
    const int MaxTotalIterations = 50;
    // We won't take a step unless it reduces the error and we'll take a
    // fractional step if necessary. If we have to use a fraction smaller than
    // this we're probably at a local minimum and it's time to give up.
    const Real MinStepFrac = Real(1e-6);
    // If the norm of a change to X during a step isn't at least this fraction
    // of a full-scale change (see scaling below), we'll treat that as an
    // independent reason to give up (squared for speed).
    const Real Xtol2 = square(Real(1e-12)); // TODO: single precision?
    // We'll calculate a length scale for the local patch of this object and
    // then limit any moves we make to this fraction of that scale to avoid
    // jumping out of one local minimum into another. This can cause a series
    // of slowly-converging steps to be taken at the beginning.
    const Real MaxMove2 = square(Real(.25)); // Limit one move to 25% of smaller radius.

    // Create an object that can calculate the error function for this
    // surface against the given line.
    TrackLineSeparationJacobian extremePointJac
       (*this, pointOnLine, directionOfLine);

    // Initialize the extreme point to the given value.
    x = startingGuessForClosestPoint;
    UnitVec3 nX; // normal at x
    Vec3 f = extremePointJac.calcExtremePointError(x, 
                nX, closestPointOnLine);
    Real frms2 = f.normSqr(); // initial error

    if (frms2 <= Ftol2) {
        //cout << "TRACK: already at extreme point with frms=" 
        //     << std::sqrt(frms2) << endl;
        height = ~(closestPointOnLine - x) * nX;
        return true; // Success
    }

    // We are going to have to move the starting point to find the new
    // extreme point.

    // Estimate a scale for the local neighborhood of this surface by
    // sampling the curvature around x, taking the largest curvature we find.
    // We want to take conservative steps that never move by more than a 
    // fraction of the scale to avoid jumping out of the local minimum.
    // TODO: could calculate the actual max curvature here but it is 
    // more expensive.
    const UnitVec3 tX = nX.perp(); // any perpendicular to the normal at x
    const UnitVec3 bX(tX % nX, true);    // another tangent vector
    // using the larger curvature in the t or b direction. 
    const Real kt = calcSurfaceCurvatureInDirection(x, tX);
    const Real kb = calcSurfaceCurvatureInDirection(x, bX);
    const Real maxK = std::max(std::abs(kt),std::abs(kb));
    const Real scale2 = square(clamp(Real(0.1), 1/maxK, Real(1000))); // keep scale reasonable

    //cout << "TRACK START: line p0=" << pointOnLine << " d=" << directionOfLine << "\n";
    //cout << "  starting x=" << x << " nX=" << nX << " scale est=" << std::sqrt(scale2) << "\n";
    //cout << "  err=" << f << " rms=" << std::sqrt(frms2) 
    //     << " closest line pt=" << closestPointOnLine
    //     << " height=" << ~(closestPointOnLine - x) * nX << "\n";

    Differentiator diff(extremePointJac);

    int stepCount = 0, limitedStepCount = 0;
    bool succeeded = false;
    do {
        ++stepCount; // we're going to take a step now
        const Matrix Jnum = diff.calcJacobian((Vector)x,
                                        Differentiator::ForwardDifference);
        const Mat33& J = Mat33::getAs(&Jnum(0,0));

        // This is the full step that the Newton would like us to take. We
        // might not use all of it.
        const Vec3 dx = J.invert()*f;
        const Real dxrms2 = dx.normSqr()/3;

        //cout << "det(J)=" << det(J)
        //     << "full dxrms=" << std::sqrt(dxrms2) << " dx=" << dx << endl;

        const Vec3 xOld     = x;        // Save previous solution & its norm.
        const Real frms2Old = frms2;

        // Backtracking. Limit the starting step size if dx is too big.
        // Calculate the square of the step fraction.
        const Real stepFrac2 = std::min(Real(1), MaxMove2*(scale2/dxrms2));
        Real stepFrac = 1;
        if (stepFrac2 < 1) {
            stepFrac = std::sqrt(stepFrac2); // not done often
            //cout << "TRACK: LIMITED STEP: iter=" << stepCount 
            //          << " stepFrac=" << stepFrac << endl;
            ++limitedStepCount;
        }
        Real xchgrms2; // norm^2 of the actual change we make to X
        while (true) {
            x = xOld - stepFrac*dx;
            xchgrms2 = stepFrac2*dxrms2; // = |stepFrac*dx|^2 / 3

            f = extremePointJac.calcExtremePointError(x, nX, closestPointOnLine);
            frms2 = f.normSqr()/3;
            if (frms2 < frms2Old || stepFrac <= MinStepFrac)
                break;

            stepFrac /= 2;
        }

        //cout << stepCount << ": TRACK lam=" << stepFrac << " |lam*dx|=" << (stepFrac*dx).norm() 
        //            << " lam*dx=" << stepFrac*dx << "-> new x=" << x << "\n"; 
        //cout << "     |f|=" << std::sqrt(frms2) << " f=" << f  << "\n";

        if (frms2 <= Ftol2) { // found solution
            //cout << "TRACK CONVERGED in " << stepCount << " steps, frms=" 
            //     << std::sqrt(frms2) << endl;
            succeeded = true;
            break;
        }

        // This method is supposed to converge quickly -- if we're not making
        // substantial progress in each step, give it up.
        if (frms2 >= 0.999*frms2Old) {
            if (frms2 > frms2Old) { // Oops -- made it worse.
                x = xOld; // Repair the damage.
                f = extremePointJac.calcExtremePointError(x, nX, closestPointOnLine);
                frms2 = f.normSqr()/3;
            }
            //cout << "TRACK FAILED at " << stepCount << " steps, frms=" 
            //     << std::sqrt(frms2) << endl;
            break;
        }

        // We took a step and made an improvement but haven't converged yet.

        if (xchgrms2 < Xtol2) { // check step size
            //std::cout << "TRACK: STALLED on step size, xchg=" << std::sqrt(xchgrms2) 
            //            << " frms=" << std::sqrt(frms2) << std::endl;
            break;
        }

        if (stepCount-limitedStepCount > MaxNewtonIterations
            || stepCount > MaxTotalIterations) {
            //cout << "TRACK: MAX iterations taken" << endl;
            break; // Return whatever we got.
        }

    } while (true);

    if (!succeeded)
        x = projectDownhillToNearestPoint(x); // push to surface
    
    height = ~(closestPointOnLine - x) * nX;
    
    //cout << "TRACK END:  x=" << x << " nX=" << nX << "\n";
    //cout << "  err=" << f << " rms=" << std::sqrt(frms2) 
    //     << " closest line pt=" << closestPointOnLine 
    //     << " height=" << height << "\n";
    return succeeded;
}

//==============================================================================
//                  GEODESIC EVALUATORS in CONTACT GEOMETRY
//==============================================================================


void ContactGeometry::initGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& xSP, const GeodesicOptions& options, Geodesic& geod) const
{
    // TODO
}


// Given two points and previous geodesic curve close to the points, find
// a geodesic curve connecting the points that is close to the previous geodesic.
void ContactGeometry::
continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
                 const GeodesicOptions& options, Geodesic& geod) const {
    getImpl().continueGeodesic(xP, xQ, prevGeod, options, geod);
}

void ContactGeometry::
makeStraightLineGeodesic(const Vec3& xP, const Vec3& xQ,
        const UnitVec3& defaultDirectionIfNeeded,
        const GeodesicOptions& options, Geodesic& geod) const {
    getImpl().makeStraightLineGeodesic(xP, xQ, defaultDirectionIfNeeded, 
                                       options, geod);
}


// Utility method to find geodesic between P and Q
// with starting directions tPhint and tQhint
// XXX tangent basis should be formed from previous geodesic
void ContactGeometry::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {
    getImpl().calcGeodesic(xP, xQ, tPhint, tQhint, geod);
}

void ContactGeometry::calcGeodesicUsingOrthogonalMethod
   (const Vec3& xP, const Vec3& xQ,
    const Vec3& tPhint, Real lengthHint, Geodesic& geod) const {
    getImpl().calcGeodesicUsingOrthogonalMethod
       (xP, xQ, tPhint, lengthHint, geod);
}

// Compute a geodesic curve starting at the given point, starting in the given
// direction, and terminating at the given plane.
// XXX what to do if tP is not in the tangent plane at P -- project it?
// XXX what to do if we don't hit the plane
void ContactGeometry::
shootGeodesicInDirectionUntilPlaneHit(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {
    getImpl().shootGeodesicInDirectionUntilPlaneHit(xP, tP, terminatingPlane,
            options, geod);
}


// Compute a geodesic curve of the given length, starting at the given point and
// in the given direction.
// XXX what to do if tP is not in the tangent plane at P -- project it?
void ContactGeometry::
shootGeodesicInDirectionUntilLengthReached(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options,
        Geodesic& geod) const {
    getImpl().shootGeodesicInDirection(xP, tP, terminatingLength, options, geod);
}

void ContactGeometry::
calcGeodesicReverseSensitivity(Geodesic& geodesic, const Vec2& initSensitivity)
    const
{
    getImpl().calcGeodesicReverseSensitivity(geodesic, initSensitivity);
}


void ContactGeometry::
shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {
    getImpl().shootGeodesicInDirectionUntilLengthReachedAnalytical(xP, tP,
            terminatingLength, options, geod);
}

void ContactGeometry::
shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {
    getImpl().shootGeodesicInDirectionUntilPlaneHitAnalytical(xP, tP,
            terminatingPlane, options, geod);
}

void ContactGeometry::
calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
            const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {
    getImpl().calcGeodesicAnalytical(xP, xQ, tPhint, tQhint, geod);
}


Vec2 ContactGeometry::calcSplitGeodError(const Vec3& xP, const Vec3& xQ,
        const UnitVec3& tP, const UnitVec3& tQ,
        Geodesic* geod) const {
    return getImpl().calcSplitGeodError(xP, xQ, tP, tQ, geod);
}

Vec2 ContactGeometry::calcSplitGeodErrorAnalytical(const Vec3& xP, const Vec3& xQ,
        const UnitVec3& tP, const UnitVec3& tQ,
        Geodesic* geod) const {
    return getImpl().calcSplitGeodErrorAnalytical(xP, xQ, tP, tQ, geod);
}

const Plane& ContactGeometry::getPlane() const  { return getImpl().getPlane(); }
void ContactGeometry::setPlane(const Plane& plane) const { getImpl().setPlane(plane); }
const Geodesic& ContactGeometry::getGeodP() const { return getImpl().getGeodP(); }
const Geodesic& ContactGeometry::getGeodQ() const { return getImpl().getGeodQ(); }
const int ContactGeometry::getNumGeodesicsShot() const { return getImpl().getNumGeodesicsShot(); }
void ContactGeometry::addVizReporter(ScheduledEventReporter* reporter) const {
    getImpl().addVizReporter(reporter);
}



//==============================================================================
//                          CONTACT GEOMETRY IMPL
//==============================================================================



//------------------------------------------------------------------------------
//                           CALC SURFACE VALUE
//------------------------------------------------------------------------------
Real ContactGeometryImpl::
calcSurfaceValue(const Vec3& p) const {
    const Vector point(p); // required by Function
    return getImplicitFunction().calcValue(point);
}



//------------------------------------------------------------------------------
//                          CALC SURFACE GRADIENT
//------------------------------------------------------------------------------
Vec3 ContactGeometryImpl::
calcSurfaceGradient(const Vec3& p) const {
    const Vector point(p); // required by Function
    const Function& f = getImplicitFunction();

    // Arguments to get first derivative from the calcDerivative interface.
    // Avoid heap allocation by making ArrayViews of stack data.
    const int x = 0; const ArrayViewConst_<int> fx(&x, &x+1);
    const int y = 1; const ArrayViewConst_<int> fy(&y, &y+1);
    const int z = 2; const ArrayViewConst_<int> fz(&z, &z+1);

    Vec3 grad;
    // Note that the gradient may point inward or outward depending on the
    // sign convention used by the implicit surface function.
    grad[0] = f.calcDerivative(fx,point);
    grad[1] = f.calcDerivative(fy,point);
    grad[2] = f.calcDerivative(fz,point);
    return grad;
}



//------------------------------------------------------------------------------
//                         CALC SURFACE UNIT NORMAL
//------------------------------------------------------------------------------
UnitVec3 ContactGeometryImpl::
calcSurfaceUnitNormal(const Vec3& p) const {
    // Implicit surface functions may have singularities away from the surface,
    // such as a point along the central axis of a cylinder. This would 
    // produce a NaN unit normal; we'll instead move the point slightly to
    // return a valid nearby normal. This helps algorithms that are looking
    // for the surface to get there rather than blow up.
    Vec3 grad = calcSurfaceGradient(p);
    Real gradMag = grad.norm();

    if (gradMag < TinyReal) {
        // Try perturbing in x, y, or z and take the first one that has a 
        // non-zero gradient.
        for (int i=0; i < 3; ++i) {
            Vec3 phat(p); phat[i] += SqrtEps;
            grad=calcSurfaceGradient(phat); 
            gradMag = grad.norm();
            if (gradMag >= TinyReal)
                break;
        }
        if (gradMag < TinyReal) {
            // We're desperate now. Pull a normal out of our hat.
            grad=Vec3(Real(1),Real(1.1),Real(1.2)); gradMag=grad.norm();
        }
    }

    // Implicit surfaces are defined as positive inside and negative outside
    // therefore normal is the negative of the gradient.
    // TODO: this should be changed to the opposite convention.
    return UnitVec3(-grad/gradMag, true);
}



//------------------------------------------------------------------------------
//                          CALC SURFACE HESSIAN
//------------------------------------------------------------------------------
// Note: Hessian is symmetric, although we're filling in all 9 elements here.
// TODO: use a SymMat33.
Mat33 ContactGeometryImpl::
calcSurfaceHessian(const Vec3& p) const {
    const Vector point(p); // required by Function
    const Function& f = getImplicitFunction();

    // Arguments to get second derivatives from the calcDerivative interface.
    // No heap allocation here.
    const int xx[] = {0,0}; ArrayViewConst_<int> fxx(xx,xx+2);
    const int xy[] = {0,1}; ArrayViewConst_<int> fxy(xy,xy+2);
    const int xz[] = {0,2}; ArrayViewConst_<int> fxz(xz,xz+2);
    //const int yx[] = {1,0}; ArrayViewConst_<int> fyx(yx,yx+2);
    const int yy[] = {1,1}; ArrayViewConst_<int> fyy(yy,yy+2);
    const int yz[] = {1,2}; ArrayViewConst_<int> fyz(yz,yz+2);
    //const int zx[] = {2,0}; ArrayViewConst_<int> fzx(zx,zx+2);
    //const int zy[] = {2,1}; ArrayViewConst_<int> fzy(zy,zy+2);
    const int zz[] = {2,2}; ArrayViewConst_<int> fzz(zz,zz+2);

    Mat33 hess;

    hess(0,0) = f.calcDerivative(fxx,point);
    hess(0,1) = f.calcDerivative(fxy,point);
    hess(0,2) = f.calcDerivative(fxz,point);
    hess(1,0) = hess(0,1);
    hess(1,1) = f.calcDerivative(fyy,point);
    hess(1,2) = f.calcDerivative(fyz,point);
    hess(2,0) = hess(0,2);
    hess(2,1) = hess(1,2);
    hess(2,2) = f.calcDerivative(fzz,point);

    return hess;
}



//------------------------------------------------------------------------------
//                        CALC GAUSSIAN CURVATURE
//------------------------------------------------------------------------------
Real ContactGeometryImpl::
calcGaussianCurvature(const Vec3&  g, const Mat33& H) const {
    // Calculate the adjoint matrix. Watch the signs!
    Mat33 A;
    A(0,0) = H(1,1)*H(2,2) - square(H(1,2)); // fyy*fzz - fyz^2
    A(0,1) = H(0,2)*H(1,2) - H(0,1)*H(2,2);  // fxz*fyz - fxy*fzz
    A(0,2) = H(0,1)*H(1,2) - H(0,2)*H(1,1);  // fxy*fyz - fxz*fyy
    A(1,0) = A(0,1);
    A(1,1) = H(0,0)*H(2,2) - square(H(0,2)); // fxx*fzz - fxz^2
    A(1,2) = H(0,1)*H(0,2) - H(0,0)*H(1,2);  // fxy*fxz - fxx*fyz
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);
    A(2,2) = H(0,0)*H(1,1) - square(H(0,1)); // fxx*fyy - fxy^2

    Real Kg = ~g * (A*g) / square(g.normSqr()); // |g|^4
    return Kg;
}



//------------------------------------------------------------------------------
//                    CALC SURFACE CURVATURE IN A DIRECTION
//------------------------------------------------------------------------------
Real ContactGeometryImpl::
calcSurfaceCurvatureInDirection(const Vec3& point, 
                                const UnitVec3& direction) const 
{
    const UnitVec3 nn = calcSurfaceUnitNormal(point);
    const Vec3     g  = calcSurfaceGradient(point);
    const Mat33 H = calcSurfaceHessian(point);
    const Real  knum = ~direction*H*direction; // numerator
    if (std::abs(knum) < TinyReal)
        return 0; // don't want to return 0/0.
        
    const Real k = knum/(~g*nn);

    return k;
}

//------------------------------------------------------------------------------
//                    CALC SURFACE PRINCIPAL CURVATURES
//------------------------------------------------------------------------------
// The idea here is to calculate the principal curvatures and their directions
// at a given point on an implicit surface, using only the implicit function and
// its derivatives without knowning anything about the actual surface.
//
// The method is to form the Shape Operator (a.k.a. Weingarten Map), a 2x2
// matrix whose eigenvalues and eigenvectors are the quantities we want.
//
// References (thanks to Ian Stavness for sending these):
//  [1] Snyder, John. Deriving the Weingarten Map. Microsoft Research 2011
//      http://research.microsoft.com/en-us/um/people/johnsny/papers/weingarten.docx
//  [2] Zhihong, Mao. Curvature computing based on shape operator for 
//      implicit surfaces. (www.paper.edu.cn, 2013)
//      http://www.paper.edu.cn/index.php/default/en_releasepaper/content/4574762
// If the above links are dead, I have pdfs.
//
// I used Snyder's idea for forming the Weingarten Map (see eqn. 10) and 
// Zhihong's equations for getting its eigenvalues and eignenvectors. 
// (sherm 140826)
void ContactGeometryImpl::
calcSurfacePrincipalCurvatures(const Vec3&  point,
                               Vec2&        curvature,
                               Rotation&    R_SP) const
{
    const UnitVec3 nn = calcSurfaceUnitNormal(point);
    const Vec3     g  = calcSurfaceGradient(point);
    const Mat33    H  = calcSurfaceHessian(point);

    // 1/2 signed length of gradient. The 1/2 here simplifies the rest of
    // the code but makes it not quite match the papers.
    const Real     oow2 = Real(0.5)/(~g*nn); 

    // Create an arbitrary frame with z=nn.
    const Rotation R_SF(nn, ZAxis); //                 ~60 flops
    const UnitVec3& t1 = R_SF.x(); const UnitVec3& t2 = R_SF.y();

    const Vec3 Ht1(H*t1), Ht2(H*t2); //                 30 flops
    const Real a = (~t1 * Ht1)*oow2; // k11             16
    const Real b = (~t1 * Ht2)*oow2; // k12             16
    const Real c = (~t2 * Ht2)*oow2; // k22             16

    // Weingarten or Shape operator is a 2x2 sym matrix 2*[a b;b c].
    // Its eigenvalues are the principal curvatures and eigenvectors the
    // principal curvature directions in the t1,t2 basis.

    // If b is zero then the matrix is diagonal and its eigenvalues are
    // 2a and 2c in directions [1,0] (t1) and [0,1] (t2) resp. Must order 
    // correctly so x is kmax direction and y kmin direction.
    if (std::abs(b) < SignificantReal) {
        if (a >= c) {
            curvature = 2*Vec2(a,c);
            R_SP = R_SF; // t1 is vmax, t2 is vmin
        } else { // c > a so t2 is max direction; flip t1 so right handed
            curvature = 2*Vec2(c,a);
            R_SP.setRotationFromUnitVecsTrustMe(t2, -t1, nn);
        }
        return;
    }

    const Real d = std::sqrt(square(a-c)+4*b*b);    // discriminant    ~25 flops
    curvature = Vec2(a+c + d, a+c - d);             // kmax, kmin        4
    const Vec2 vmax2(a-c + d, 2*b);                 //                   3
    const UnitVec3 vmax(vmax2[0]*t1 + vmax2[1]*t2); //                 ~40
    const UnitVec3 vmin(nn % vmax, true); // already a unit vector       9
    R_SP.setRotationFromUnitVecsTrustMe(vmax, vmin, nn);
}

// used in numerical differentiation. TODO: what value for single precision?
static const Real estimatedGeodesicAccuracy = Real(1e-12); 
static const Real pauseBetweenGeodIterations = 0; // sec, used in newton solver


// This local class is used to Calcualte the split geodesic error
//  given scalar angles at P and Q
class ContactGeometryImpl::SplitGeodesicError: public Differentiator::JacobianFunction {

public:
    SplitGeodesicError(int nf, int ny, const ContactGeometryImpl& geom,
            const Vec3& xP, const Vec3& xQ, 
            const Vec3& tPhint, const Vec3& tQhint) 
    :   Differentiator::JacobianFunction(nf, ny),
                    geom(geom),
                    P(xP), Q(xQ),
                    R_SP(geom.calcTangentBasis(P, tPhint)),
                    R_SQ(geom.calcTangentBasis(Q, tQhint)) { }

    // x = ~[thetaP, thetaQ]
    int f(const Vector& x, Vector& fx) const override  {
        UnitVec3 tP = calcUnitTangentVec(x[0], R_SP);
        UnitVec3 tQ = calcUnitTangentVec(x[1], R_SQ);

        Vec2 geodErr = geom.calcSplitGeodError(P, Q, tP, tQ);

        // error between geodesic end points at plane
        fx[0] = geodErr[0];
        fx[1] = geodErr[1];
        return 0;
    }


private:
    const ContactGeometryImpl& geom;
    const Vec3& P;
    const Vec3& Q;
    const Rotation R_SP;
    const Rotation R_SQ;

}; // class SplitGeodesicError

// This local class is used to calculate the orthogonal geodesic error
// Jacobian given angle theta and length.
class ContactGeometryImpl::OrthoGeodesicError: public Differentiator::JacobianFunction {

public:
    OrthoGeodesicError(const ContactGeometryImpl& geom,
            const Vec3& xP, const Vec3& xQ) :
            Differentiator::JacobianFunction(2, 2),
                    geom(geom), P(xP), Q(xQ) { }

    // x = ~[thetaP, length]
    int f(const Vector& x, Vector& fx) const override  {
        Geodesic geod;
        Vec2 geodErr = geom.calcOrthogonalGeodError(P, Q, x[0], x[1], geod);

        // error between geodesic end points at plane
        fx[0] = geodErr[0];
        fx[1] = geodErr[1];
        return 0;
    }


private:
    const ContactGeometryImpl& geom;
    const Vec3 P;
    const Vec3 Q;
}; // class OrthoGeodesicError


//------------------------------------------------------------------------------
//                              INIT GEODESIC
//------------------------------------------------------------------------------

void ContactGeometryImpl::initGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& xSP, const GeodesicOptions& options, Geodesic& geod) const
{
    // TODO
}



//------------------------------------------------------------------------------
//                            CONTINUE GEODESIC
//------------------------------------------------------------------------------
// Given two points and a previous geodesic curve close to the points, find
// a geodesic curve connecting the points that is close to the previous 
// geodesic. See header or doxygen for algorithmic details.
// Note that the geodesic runs from P' to Q', which are closest-point
// projections of the initial points which might not be on the surface.
void ContactGeometryImpl::
continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
                 const GeodesicOptions& options, Geodesic& geod) const 
{
    const Real StraightLineGeoFrac = Real(1e-5); // a straight line

    const Vec3 P = projectDownhillToNearestPoint(xP);
    const Vec3 Q = projectDownhillToNearestPoint(xQ);

    //cout << "continueGeo(P=" << P << ",Q=" << Q << ")\n";

    // If P and Q are bit-identical to P and Q from the previous geodesic,
    // then the new one is just a copy of the previous one. This is especially
    // likely when doing numerical Jacobians since many value remain unchanged.
    if (prevGeod.getNumPoints() 
        && prevGeod.getPointP()==P && prevGeod.getPointQ()==Q) {
            geod = prevGeod;
            cout << "REUSING OLD GEODESIC of length=" 
                << prevGeod.getLength() << "\n";
            return;
    }
   
    const Vec3 PQ = Q-P;
    const Real PQlength = PQ.norm();
    const UnitVec3 PQdir =
        PQlength == 0 ? UnitVec3(XAxis) : UnitVec3(PQ/PQlength, true);
    cout << "  PQlen=" << PQlength << " PQdir=" << PQdir << "\n";

    // If the length is less than this fraction of the maximum radius of
    // curvature (1/kdP) then the geodesic is indistinguishable from a 
    // straight line.
    // TODO: if the previous geodesic was very long and came all the way
    // around this will incorrectly switch to a very short geodesic. Does
    // that matter?
    const Real kdP = std::abs(calcSurfaceCurvatureInDirection(P,PQdir));
    if (PQlength*kdP <= StraightLineGeoFrac) {
        cout << "STRAIGHT LINE GEO length=" << PQlength << "\n";
        makeStraightLineGeodesic(P, Q, PQdir, options, geod);
        return;
    }

    // If there is no previous geodesic the best we can do is to start in
    // direction PQdir, and guess length |PQ|.
    if (!prevGeod.getNumPoints()) {
        calcGeodesicUsingOrthogonalMethod(P, Q, PQdir, PQlength, geod);
        //calcGeodesicAnalytical(P, Q, PQ, PQ, geod);
        return;
    }

    // First classify the previous geodesic as direct or indirect. Direct is
    // a strict classification; only if the end tangents are aligned with the
    // PQ line to within an allowed cone angle is it direct.
    const Real CosMaxDirectAngle = Real(0.9); // about 25 degrees

    // Find maximum curvature of previous geodesic to use as a length scale.
    // TODO: should store this with the geodesic so we can use any intermediate
    // points also. TODO: use binormal curvature also?
    const Real maxK = std::max(std::abs(prevGeod.getCurvatureP()),
                               std::abs(prevGeod.getCurvatureQ()));

    // We consider a geodesic straight when the separation of its end points
    // is a small fraction of the radius of curvature because then the chord
    // and the arc are the same length to a MUCH smaller tolerance than that.
    const Vec3 prevPQ = prevGeod.getPointQ() - prevGeod.getPointP();
    const Real prevPQlen = prevPQ.norm();
    const bool isPrevStraightLine = prevPQlen*maxK <= StraightLineGeoFrac;

    const UnitVec3 eprevPQ = isPrevStraightLine
                                ? prevGeod.getTangentP()
                                : UnitVec3(prevPQ/prevPQlen, true);

    const Real cosConetP = dot(prevGeod.getTangentP(), eprevPQ);
    const Real cosConetQ = dot(prevGeod.getTangentQ(), eprevPQ);
    const bool isDirect = std::min(cosConetP,cosConetQ) >= CosMaxDirectAngle;

    UnitVec3 tPhint = prevGeod.getTangentP(); // might flip
    UnitVec3 tQhint = prevGeod.getTangentQ();
    Real     sHint  = std::max(prevGeod.getLength(), PQlength);

    if (isDirect) {
        if (~tPhint*PQdir < 0 && ~tQhint*PQdir < 0) {
            tPhint = PQdir;
            tQhint = PQdir;
            sHint = PQlength;
            cout << "GEODESIC FLIPPED. Prev len was " 
                 << prevGeod.getLength() << endl;
        }
    }

    calcGeodesicUsingOrthogonalMethod(P, Q, tPhint, sHint, geod);
    //calcGeodesicAnalytical(P, Q, tPhint, tQhint, geod);
}



//------------------------------------------------------------------------------
//                        MAKE STRAIGHT LINE GEODESIC
//------------------------------------------------------------------------------
void ContactGeometryImpl::
makeStraightLineGeodesic(const Vec3& xP, const Vec3& xQ,
        const UnitVec3& defaultDirectionIfNeeded,
        const GeodesicOptions& options, Geodesic& geod) const
{
    geod.clear();

    Vec3 Pprime = projectDownhillToNearestPoint(xP);
    Vec3 Qprime = projectDownhillToNearestPoint(xQ);

    const bool isZeroLength = 
        Geo::Point::pointsAreNumericallyCoincident(Pprime,Qprime);

    UnitVec3 d;
    Real     length;
    if (isZeroLength) {
        const Vec3 mid = Geo::Point::findMidpoint(Pprime, Qprime);
        Pprime = Qprime = projectDownhillToNearestPoint(mid);
        d = defaultDirectionIfNeeded;
        length = 0;
    } else {
        const Vec3 PQ = Qprime - Pprime;
        length = PQ.norm();
        d = UnitVec3(PQ/length, true);
    }

    const UnitVec3 nP = calcSurfaceUnitNormal(Pprime);
    const UnitVec3 nQ = calcSurfaceUnitNormal(Qprime);
    const Rotation RP(nP, ZAxis, d, YAxis);
    const Rotation RQ(nQ, ZAxis, d, YAxis);
    geod.addFrenetFrame(Transform(RP, Pprime));
    geod.addFrenetFrame(Transform(RQ, Qprime));
    geod.addArcLength(0);
    geod.addArcLength(length);
    geod.addCurvature(calcSurfaceCurvatureInDirection(Pprime, RP.y()));
    geod.addCurvature(calcSurfaceCurvatureInDirection(Qprime, RQ.y()));
    geod.addDirectionalSensitivityPtoQ(Vec2(0,1));
    geod.addDirectionalSensitivityPtoQ(Vec2(length,1));
    geod.addDirectionalSensitivityQtoP(Vec2(length,1));
    geod.addDirectionalSensitivityQtoP(Vec2(0,1));
    geod.addPositionalSensitivityPtoQ(Vec2(1,0));
    geod.addPositionalSensitivityPtoQ(Vec2(1,0));
    geod.addPositionalSensitivityQtoP(Vec2(1,0));
    geod.addPositionalSensitivityQtoP(Vec2(1,0));
    geod.setBinormalCurvatureAtP(calcSurfaceCurvatureInDirection(Pprime, RP.x()));
    geod.setBinormalCurvatureAtQ(calcSurfaceCurvatureInDirection(Qprime, RQ.x()));

    //TODO: We're estimating torsion here as the change in binormal per unit
    // of arc length, in the -n direction. If length is zero we'll just say
    // the torsion is zero although that is wrong.
    const Real ooLength = isZeroLength ? Real(0) : 1/length;
    const Vec3 bChg = geod.getBinormalQ()-geod.getBinormalP();
    geod.setTorsionAtP( (~nP * bChg) * ooLength );
    geod.setTorsionAtQ( (~nQ * bChg) * ooLength );

    geod.setIsConvex(geod.getCurvatureP() >= 0 && geod.getCurvatureQ() >= 0);
    geod.setIsShortest(true);
    geod.setInitialStepSizeHint(NaN); // no clue
    geod.setAchievedAccuracy(Geo::getDefaultTol<Real>());
}



//------------------------------------------------------------------------------
//                        SHOOT GEODESIC IN DIRECTION
//------------------------------------------------------------------------------
// Utility method used by shootGeodesicInDirectionUntilPlaneHit
// and shootGeodesicInDirectionUntilLengthReached.
// P is not necessarily on the surface and tP is not necessarily tangent
// to the surface, so the first thing we do is project P to the nearest
// surface point P', calculate the outward unit normal n' at P', then
// project tP to tP' by removing any component it has in the n' direction,
// then renormalizing.
static const Real IntegratorAccuracy = Real(1e-6); // TODO: how to choose?
static const Real IntegratorConstraintTol = Real(1e-10);
void ContactGeometryImpl::
shootGeodesicInDirection(const Vec3& P, const UnitVec3& tP,
        const Real& finalTime, const GeodesicOptions& options,
        Geodesic& geod) const {

    // integrator settings
    const Real startTime = 0;

    ++numGeodesicsShot;

    // Initialize state
    State sysState = ptOnSurfSys->getDefaultState();
    sysState.setTime(startTime);
    Vector& q = sysState.updQ();
    Vector& u = sysState.updU();
    q[0] = P[0]; q[1] = P[1]; q[2] = P[2];
    u[0] = tP[0]; u[1] = tP[1]; u[2] = tP[2];

    // Jacobi field states
    q[3] = 0; 
    u[3] = 1;

    // Setup integrator to integrate until terminatingLength
    //RungeKutta3Integrator integ(*ptOnSurfSys);
    RungeKuttaMersonIntegrator integ(*ptOnSurfSys);
    //RungeKuttaFeldbergIntegrator integ(*ptOnSurfSys);
    integ.setAccuracy(IntegratorAccuracy);
    integ.setConstraintTolerance(IntegratorConstraintTol);
    integ.setFinalTime(finalTime);
    integ.setReturnEveryInternalStep(true); // save geodesic knot points

    // Setup timestepper in order to handle event when geodesic hits the plane
    TimeStepper ts(*ptOnSurfSys, integ);
    ts.setReportAllSignificantStates(true);
    ts.initialize(sysState);

    // Simulate it, and record geodesic knot points after each step
    // Terminate when geodesic hits the plane
    Integrator::SuccessfulStepStatus status;
    int stepcnt = 0;
    geod.setIsConvex(true); // Set false if we see negative curvature anywhere.
    while (true) {
        // Final time is already reported by the time we see end of simulation;
        // don't duplicate the last step.
        if ((status=ts.stepTo(Infinity)) == Integrator::EndOfSimulation)
            break;

        // If we stopped just for a report, don't add that to the geodesic;
        // that may happen if a visualization reporter has been hung on this
        // system for watching the geodesic grow.
        if (status == Integrator::ReachedReportTime)
            continue;

        // Careful -- integrator has its own internal copies of the state and
        // the one returned won't always be the same one so you need to get a
        // fresh reference each time.
        const State& state = integ.getState();
        const Real s = state.getTime();
        geod.addArcLength(s);

        const Vec3& pt = Vec3::getAs(&state.getQ()[0]);
        const UnitVec3 n = calcSurfaceUnitNormal(pt);
        const Vec3& tangent = Vec3::getAs(&state.getU()[0]);
        // Rotation will orthogonalize so x direction we get may not be
        // exactly the same as what we supply here.
        const Transform frenetFrame(Rotation(n, ZAxis, tangent, YAxis), pt);
        geod.addFrenetFrame(frenetFrame);
        geod.addDirectionalSensitivityPtoQ(Vec2(state.getQ()[3],
                                                state.getU()[3]));
        geod.addPositionalSensitivityPtoQ(Vec2(NaN,NaN)); // XXX
        const Real kappa = calcSurfaceCurvatureInDirection(pt, frenetFrame.y());
        geod.addCurvature(kappa);
        if (kappa < 0) 
            geod.setIsConvex(false);

        ++stepcnt;
    }

    geod.setBinormalCurvatureAtP(
        calcSurfaceCurvatureInDirection(geod.getPointP(),geod.getBinormalP()));
    geod.setBinormalCurvatureAtQ(
        calcSurfaceCurvatureInDirection(geod.getPointQ(),geod.getBinormalQ()));

    //TODO: numerical torsion estimate for now
    const Array_<Transform>& frenet = geod.getFrenetFrames();
    const Array_<Real>& arcLen = geod.getArcLengths();
    const int last = geod.getNumPoints()-1;
    const Real lFirst = arcLen[1];
    const Real lLast  = arcLen[last] - arcLen[last-1];
    const Real tauP = lFirst==0 ? Real(0)
        : -dot(frenet[0].z(), 
              (frenet[1].x()-frenet[0].x())) / lFirst; // dbP/ds
    const Real tauQ = lLast==0 ? Real(0)
        : -dot(frenet[last].z(), 
              (frenet[last].x()-frenet[last-1].x())) / lLast;

    geod.setTorsionAtP(tauP); geod.setTorsionAtQ(tauQ);

    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(IntegratorAccuracy); // TODO: accuracy of length?
    // TODO: better to use something like the second-to-last step, or average
    // excluding initial and last steps, so that we don't have to start small.
    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken());
}

void ContactGeometryImpl::
shootGeodesicInDirection2(const Vec3& P, const UnitVec3& tP,
        const Real& finalArcLength, const GeodesicOptions& options,
        Geodesic& geod) const {

    // integrator settings
    const Real startArcLength = 0;

    GeodesicOnImplicitSurface eqns(*this);
    GeodesicIntegrator<GeodesicOnImplicitSurface> 
        integ(eqns,IntegratorAccuracy,IntegratorConstraintTol);
    static const int N = GeodesicOnImplicitSurface::N;

    ++numGeodesicsShot;

    integ.initialize(startArcLength, eqns.getInitialState(P,tP));
    // Aliases for the integrators internal time and state variables.
    const Vec<N>& y = integ.getY();
    const Real&   s = integ.getTime();  // arc length

    // Simulate it, and record geodesic knot points after each step
    int stepcnt = 0;
    geod.setIsConvex(true); // Set false if we see negative curvature anywhere.
    while (true) {
        // Record a knot point.
        geod.addArcLength(s);
        const Vec3& pt = eqns.getP(y);
        const UnitVec3 n = calcSurfaceUnitNormal(pt);
        const Vec3& tangent = eqns.getV(y);
        // Rotation will orthogonalize so x direction we get may not be
        // exactly the same as what we supply here.
        const Transform frenetFrame(Rotation(n, ZAxis, tangent, YAxis), pt);
        geod.addFrenetFrame(frenetFrame);
        geod.addDirectionalSensitivityPtoQ(Vec2(eqns.getJRot(y),
                                                eqns.getJRotDot(y)));
        geod.addPositionalSensitivityPtoQ(Vec2(eqns.getJTrans(y),
                                               eqns.getJTransDot(y)));
        const Real kappa = calcSurfaceCurvatureInDirection(pt, frenetFrame.y());
        geod.addCurvature(kappa);
        if (kappa < 0) geod.setIsConvex(false);

        if (s == finalArcLength)
            break;

        integ.takeOneStep(finalArcLength);
        ++stepcnt;
    }

    //printf("RKM acc=%g tol=%g: %d/%d steps, errtest=%d projfail=%d\n",
    //    integ.getRequiredAccuracy(), integ.getConstraintTolerance(),
    //    integ.getNumStepsTaken(), integ.getNumStepsAttempted(),
    //    integ.getNumErrorTestFailures(), integ.getNumProjectionFailures());

    geod.setBinormalCurvatureAtP(
        calcSurfaceCurvatureInDirection(geod.getPointP(),geod.getBinormalP()));
    geod.setBinormalCurvatureAtQ(
        calcSurfaceCurvatureInDirection(geod.getPointQ(),geod.getBinormalQ()));

    //TODO: numerical torsion estimate for now
    const Array_<Transform>& frenet = geod.getFrenetFrames();
    const Array_<Real>& arcLen = geod.getArcLengths();
    const int last = geod.getNumPoints()-1;
    const Real lFirst = arcLen[1];
    const Real lLast  = arcLen[last] - arcLen[last-1];
    const Real tauP = lFirst==0 ? Real(0)
        : -dot(frenet[0].z(), 
              (frenet[1].x()-frenet[0].x())) / lFirst; // dbP/ds
    const Real tauQ = lLast==0 ? Real(0)
        : -dot(frenet[last].z(), 
              (frenet[last].x()-frenet[last-1].x())) / lLast;

    geod.setTorsionAtP(tauP); geod.setTorsionAtQ(tauQ);

    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(IntegratorAccuracy); // TODO: accuracy of length?
    // TODO: better to use something like the second-to-last step, or average
    // excluding initial and last steps, so that we don't have to start small.
    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken());
}


//------------------------------------------------------------------------------
//                      CALC GEODESIC REVERSE SENSITIVITY
//------------------------------------------------------------------------------
// After a geodesic has been calculated, this method integrates backwards
// to fill in the missing reverse Jacobi term.
void ContactGeometryImpl::
calcGeodesicReverseSensitivity(Geodesic& geod, const Vec2& initJacobi) const {
    
    // Don't look for a plane.
    geodHitPlaneEvent->setEnabled(false);

    // integrator settings

    //RungeKutta3Integrator integ(ptOnSurfSys);
    RungeKuttaMersonIntegrator integ(*ptOnSurfSys);
    integ.setAccuracy(IntegratorAccuracy);
    integ.setConstraintTolerance(IntegratorConstraintTol);
    State sysState = ptOnSurfSys->getDefaultState();
    Vector& q = sysState.updQ();
    Vector& u = sysState.updU();

    // Initial Jacobi field states.
    q[3] = initJacobi[0]; 
    u[3] = initJacobi[1];

    Array_<Vec2>& jQ = geod.updDirectionalSensitivityQtoP();
    jQ.resize(geod.getNumPoints());
    jQ.back() = initJacobi;

    for (int step=geod.getNumPoints()-1; step >= 1; --step) {
        // Curve goes from P to Q. We have to integrate backwards from Q to P.
        const Transform& QFrenet = geod.getFrenetFrames()[step];
        const Transform& PFrenet = geod.getFrenetFrames()[step-1];
        const Real sQ = geod.getArcLengths()[step];
        const Real sP = geod.getArcLengths()[step-1];
        const Vec3&      Q = QFrenet.p();
        const UnitVec3&  tQ = QFrenet.x(); // we'll reverse this
        const Vec3&      P = PFrenet.p();
        const UnitVec3&  tP = PFrenet.x();

        // Initialize state
        sysState.setTime(0);
        q[0] = Q[0]; q[1] = Q[1]; q[2] = Q[2];
        u[0] = -tQ[0]; u[1] = -tQ[1]; u[2] = -tQ[2];
        q[3] = jQ[step][0]; 
        u[3] = jQ[step][1];

        const Real arcLength = sQ-sP; // how far to integrate

        integ.initialize(sysState);

        Integrator::SuccessfulStepStatus status;
        do {
            status = integ.stepTo(arcLength);
            if (status == Integrator::StartOfContinuousInterval)
                continue;
            if (integ.getTime() < arcLength) {
                printf("integ to %g returned early at %g with status=%s\n",
                    arcLength, integ.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            }
        } while (integ.getTime() < arcLength);

        // Save Jacobi field value.
        const State& state = integ.getState();
        jQ[step-1] = Vec2(state.getQ()[3],state.getU()[3]);
    }
}

void ContactGeometryImpl::
calcGeodesicReverseSensitivity2
   (Geodesic& geod, const Vec2& initJRot, const Vec2& initJTrans) const {

    GeodesicOnImplicitSurface eqns(*this);
    GeodesicIntegrator<GeodesicOnImplicitSurface> 
        integ(eqns,IntegratorAccuracy,IntegratorConstraintTol);
    static const int N = GeodesicOnImplicitSurface::N;

    integ.initialize(0, 
        eqns.getInitialState(geod.getPointQ(),-geod.getTangentQ()));
    // Aliases for the integrators internal time and state variables.
    const Vec<N>& y = integ.getY();
    const Real&   s = integ.getTime();  // arc length

    // These two arrays need to be filled in to complete the geodesic.
    Array_<Vec2>& jrP = geod.updDirectionalSensitivityQtoP();
    jrP.resize(geod.getNumPoints());
    jrP.back() = initJRot;

    Array_<Vec2>& jtP = geod.updPositionalSensitivityQtoP();
    jtP.resize(geod.getNumPoints());
    jtP.back() = initJTrans;

    for (int step=geod.getNumPoints()-1; step >= 1; --step) {
        // Curve goes from P to Q. We have to integrate backwards from Q to P.
        const Transform& QFrenet = geod.getFrenetFrames()[step];
        const Transform& PFrenet = geod.getFrenetFrames()[step-1];
        const Real sQ = geod.getArcLengths()[step];
        const Real sP = geod.getArcLengths()[step-1];
        const Vec3&      Q = QFrenet.p();
        const UnitVec3&  tQ = QFrenet.x(); // we'll reverse this
        const Vec3&      P = PFrenet.p();
        const UnitVec3&  tP = PFrenet.x();

        // Initialize state
        Vec<N> yInit;
        eqns.updP(yInit) = Q;
        eqns.updV(yInit) = -tQ.asVec3();
        eqns.updJRot(yInit) = jrP[step][0];
        eqns.updJRotDot(yInit) = jrP[step][1];
        eqns.updJTrans(yInit) = jtP[step][0];
        eqns.updJTransDot(yInit) = jtP[step][1];

        const Real arcLength = sQ-sP; // how far to integrate

        integ.setTimeAndState(0, yInit);
        integ.setNextStepSizeToTry(arcLength);

        do { // usually we'll be able to do this in one step
            integ.takeOneStep(arcLength);
            //if (s < arcLength) printf("integ to %g returned early at %g\n", arcLength, s);
        } while (s < arcLength);

        // Save Jacobi field values.
        jrP[step-1] = Vec2(eqns.getJRot(y), eqns.getJRotDot(y));
        jtP[step-1] = Vec2(eqns.getJTrans(y), eqns.getJTransDot(y));
    }

    //printf("REVERSE acc=%g tol=%g: %d/%d steps, errtest=%d projfail=%d\n",
    //    integ.getRequiredAccuracy(), integ.getConstraintTolerance(),
    //    integ.getNumStepsTaken(), integ.getNumStepsAttempted(),
    //    integ.getNumErrorTestFailures(), integ.getNumProjectionFailures());
}


// Compute a geodesic curve starting at the given point, starting in the given
// direction, and terminating at the given plane.
// XXX what to do if tP is not in the tangent plane at P -- project it?
// XXX what to do if we don't hit the plane
void ContactGeometryImpl::
shootGeodesicInDirectionUntilPlaneHit(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {
    geodHitPlaneEvent->setEnabled(true);
    geodHitPlaneEvent->setPlane(terminatingPlane);
    // TODO: need a reasonable max length
    const Real MaxLength = /*Infinity*/100;
    shootGeodesicInDirection(xP, tP, MaxLength, options, geod);
}


// Compute a geodesic curve of the given length, starting at the given point and
// in the given direction.
// XXX what to do if tP is not in the tangent plane at P -- project it?
void ContactGeometryImpl::
shootGeodesicInDirectionUntilLengthReached(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options,
        Geodesic& geod) const {
    geodHitPlaneEvent->setEnabled(false);
#ifdef USE_NEW_INTEGRATOR
    shootGeodesicInDirection2(xP, tP, terminatingLength, options, geod);
#else
    shootGeodesicInDirection(xP, tP, terminatingLength, options, geod);
#endif
}



static Real cleanUpH(Real hEst, Real y0) {
    volatile Real temp = y0+hEst;
    return temp-y0;
}

static Real maxabs(Vec2 x) {
    return std::max(std::abs(x[0]), std::abs(x[1]));
}

static Real maxabsdiff(Vec2 x, Vec2 xold) {
    return std::max(std::abs(x[0]-xold[0])/std::max(x[0],Real(1)),
                    std::abs(x[1]-xold[1])/std::max(x[1],Real(1)));
}



//------------------------------------------------------------------------------
//                              CALC GEODESIC
//------------------------------------------------------------------------------
// Utility method to find geodesic between P and Q
// with starting directions tPhint and tQhint
// XXX tangent basis should be formed from previous geodesic
void ContactGeometryImpl::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {

    // Newton solver settings. TODO: single precision?
    const Real ftol = Real(1e-9);
    const Real xtol = Real(1e-9);
    const Real minlam = Real(1e-9);
    const int maxNewtonIterations = 25;

    // reset counter
    numGeodesicsShot = 0;

    // define basis
    R_SP = calcTangentBasis(xP, tPhint);
    R_SQ = calcTangentBasis(xQ, tQhint);

    // calculate plane bisecting P and Q, and use as termination condition for integrator
    UnitVec3 normal(xQ - xP);
    Real offset = (~(xP+xQ)*normal)/2 ;
    geodHitPlaneEvent->setPlane(Plane(normal, offset));

    Mat22 J;
    Vec2 x, xold, dx, Fx;

    // initial conditions
    x[0] = Pi/2; // thetaP
    x[1] = -Pi/2; // thetaQ: shoot toward P, i.e. in the opposite direction of tQ
    Real f, fold, lam = 1;

//    splitGeodErr = new SplitGeodesicError(2, 2, *const_cast<ContactGeometry*>(this),
//            xP, xQ, tPhint, tQhint);
//    splitGeodErr->setEstimatedAccuracy(estimatedGeodesicAccuracy);
//    Differentiator diff( *const_cast<SplitGeodesicError*>(splitGeodErr));

//    splitGeodErr->f(x, Fx);
    Fx = calcSplitGeodError(xP, xQ, x[0], x[1]);
    if (vizReporter != NULL) {
        vizReporter->handleEvent(ptOnSurfSys->getDefaultState());
        sleepInSec(pauseBetweenGeodIterations);
    }

    f = std::sqrt(~Fx*Fx);

    for (int i = 0; i < maxNewtonIterations; ++i) {
        if (maxabs(Fx) < ftol) {
            //std::cout << "geodesic converged in " << i << " iterations" << std::endl;
//            std::cout << "err = " << Fx << std::endl;
            break;
        }
//        diff.calcJacobian(x,  Fx, J, Differentiator::ForwardDifference);
        J = calcSplitGeodErrorJacobian(xP, xQ, x[0], x[1],
                Differentiator::ForwardDifference);
        dx = J.invert()*Fx;

        fold = f;
        xold = x;

        // backtracking
        lam = 1;
        while (true) {
            x = xold - lam*dx;
//            splitGeodErr->f(x, Fx);
            Fx = calcSplitGeodError(xP, xQ, x[0], x[1]);
            f = std::sqrt(~Fx*Fx);
            if (f > fold && lam > minlam) {
                lam = lam / 2;
            } else {
                break;
            }
        }
        if (maxabsdiff(x,xold) < xtol) {
            std::cout << "converged on step size after " << i << " iterations" << std::endl;
            std::cout << "err = " << Fx << std::endl;
            break;
        }

        if (vizReporter != NULL) {
            vizReporter->handleEvent(ptOnSurfSys->getDefaultState());
            sleepInSec(pauseBetweenGeodIterations);
        }

    }

    // Finish each geodesic with reverse Jacobi field.
    calcGeodesicReverseSensitivity(geodP,
        geodQ.getDirectionalSensitivityPtoQ().back());
    calcGeodesicReverseSensitivity(geodQ,
        geodP.getDirectionalSensitivityPtoQ().back());

    mergeGeodesics(geodP, geodQ, geod);
}



//------------------------------------------------------------------------------
//                  CALC GEODESIC USING ORTHOGONAL METHOD
//------------------------------------------------------------------------------
void ContactGeometryImpl::calcGeodesicUsingOrthogonalMethod
   (const Vec3& xP, const Vec3& xQ,
    const Vec3& tPhint, Real lengthHint, Geodesic& geod) const 
{
    const Vec3 P = projectDownhillToNearestPoint(xP);
    const Vec3 Q = projectDownhillToNearestPoint(xQ);

    // Newton solver settings. TODO: single precision values?
    const Real ftol = Real(1e-9);
    const Real xtol = Real(1e-12);
    const Real minlam = Real(1e-3);
    const int MaxIterations = 25;

    bool useNewtonIteration = true;

    // reset counter
    numGeodesicsShot = 0;

    // Define basis. This will not change during the solution.
    R_SP = calcTangentBasis(P, tPhint);

    Mat22 J;
    Vec2 x, xold, dx, Fx;
    Matrix JMat;
    // initial conditions
    x[0] = Pi/2; // thetaP
    x[1] = lengthHint;
    Real f, fold, dist, lam = 1;

    Fx = calcOrthogonalGeodError(P, Q, x[0], x[1], geod);
    if (vizReporter != NULL) {
        vizReporter->handleEvent(ptOnSurfSys->getDefaultState());
        sleepInSec(pauseBetweenGeodIterations);
    }

    //OrthoGeodesicError orthoErr(*this, P, Q);
    //orthoErr.setEstimatedAccuracy(1e-12); // TODO
    //Differentiator diff(orthoErr);

    //cout << "Using " << (useNewtonIteration ? "NEWTON" : "FIXED POINT")
    //     << " iteration\n";

    for (int i = 0; i < MaxIterations; ++i) {
        f = std::sqrt(~Fx*Fx);
        const Vec3 r_QQhat = geod.getPointQ()-Q;
        dist = r_QQhat.norm();
        //std::cout << "ORTHO x= " << x << " err = " << Fx << " |err|=" << f 
        //          << " dist=" << dist << std::endl;
        if (f <= ftol) {
           // std::cout << "ORTHO geodesic converged in " 
            //          << i << " iterations with err=" << f << std::endl;
            break;
        }

        if (useNewtonIteration) {
            // This numerical Jacobian is very bad; CentralDifference is 
            // required in order to produce a reasonable one.
            //diff.calcJacobian(Vector(x),  Vector(Fx), JMat, 
            //                  Differentiator::CentralDifference);
            //J = Mat22::getAs(&JMat(0,0));

            //XXX: numerical calculation of j and jdot
            //Geodesic geod0, geod1;
            //calcOrthogonalGeodError(P, Q, x[0]-1e-5, x[1], geod0);
            //calcOrthogonalGeodError(P, Q, x[0]+1e-5, x[1], geod1);
            //Vec3 qdiff = geod1.getPointQ() - geod0.getPointQ();
            //Real num_j = dot(qdiff,geod0.getBinormalQ())/2e-5;
            //printf("Jacobi Q num=%g, analytic=%g\n", 
            //    num_j, geod0.getJacobiQ());

            //calcOrthogonalGeodError(P, Q, x[0]-1e-5, x[1]+1e-5, geod0);
            //calcOrthogonalGeodError(P, Q, x[0]+1e-5, x[1]+1e-5, geod1);
            //Vec3 qdiff2 = geod1.getPointQ() - geod0.getPointQ();
            //Real num_j2 = dot(qdiff2,geod0.getBinormalQ())/2e-5;

            //Real num_jd = (num_j2-num_j)/1e-5;
            //printf("Jacobi dot Q num=%g, analytic=%g\n",
            //    num_jd, geod0.getJacobiQDot());

            //XXX: numerical calculation of jt and jtdot
            //Geodesic geod0, geod1;
            //calcOrthogonalGeodError(P-1e-5*geod.getBinormalP(), Q, x[0], x[1], geod0);
            //calcOrthogonalGeodError(P+1e-5*geod.getBinormalP(), Q, x[0], x[1], geod1);
            //Vec3 qdiff = geod1.getPointQ() - geod0.getPointQ();
            //Real num_jt = dot(qdiff,geod0.getBinormalQ())/2e-5;
            //printf("Jacobi Trans Q num=%g, analytic=%g\n", 
            //    num_jt, geod0.getJacobiTransQ());

            //calcOrthogonalGeodError(P-1e-5*geod.getBinormalP(), Q, x[0], x[1]+1e-5, geod0);
            //calcOrthogonalGeodError(P+1e-5*geod.getBinormalP(), Q, x[0], x[1]+1e-5, geod1);
            //Vec3 qdiff2 = geod1.getPointQ() - geod0.getPointQ();
            //Real num_jt2 = dot(qdiff2,geod0.getBinormalQ())/2e-5;

            //Real num_jtd = (num_jt2-num_jt)/1e-5;
            //printf("Jacobi dot Q num=%g, analytic=%g\n",
            //    num_jtd, geod0.getJacobiTransQDot());

            const Real jt = geod.getJacobiTransQ();
            const Real jtd = geod.getJacobiTransQDot();

            const Real eh = ~r_QQhat*geod.getNormalQ();
            const Real es = ~r_QQhat*geod.getTangentQ();
            const Real eb = ~r_QQhat*geod.getBinormalQ();
            const Real j = geod.getJacobiQ();
            const Real jd = geod.getJacobiQDot();
            const Real tau = geod.getTorsionQ();
            const Real kappa = geod.getCurvatureQ();
            const Real mu = geod.getBinormalCurvatureQ();

           // printf("eh=%g es=%g, eb=%g, j=%g, tau=%g, kappa=%g, mu=%g\n",
            //  eh,es,eb,j,tau,kappa,mu);

            Mat22 newJ( -(jd*es+j*(mu*eh-1)),      -tau*eh,
                         -(j*tau*eh-jd*eb),      1-kappa*eh );

            //cout << "   J=" << J;
            //cout << "newJ=" << newJ;

            //dx = J.invert()*Fx; // Newton
            dx = newJ.invert()*Fx;
        } else {
            // fixed point -- feed error back to variables
            dx = Vec2(Fx[0]/geod.getJacobiQ(), Fx[1]); 
        }

        //cout << "f=" << f << "-> dx=" << dx << endl;

        fold = f;
        xold = x;

        // backtracking. Limit angle changes to around 22 degrees.
        const Real dtheta = std::abs(dx[0]);
        lam = 1;
        if (dtheta > Pi/8) {
            lam = (Pi/8)/dtheta; 
            //cout << "ORTHO: lam reduced to " << lam << "\n";
        }

        while (true) {
            x = xold - lam*dx;
            //cout << "at lam=" << lam << " x-lam*dx=" << x << endl;
            // Negative length means flip direction.
            Vec2 xadj = x;
            bool flipped = false;
            if (x[1] < 0) {
                cout << "NEGATIVE LENGTH; x=" << x << "; flipping\n";
                cout << "P=" << P << " Q=" << Q << " Q-P=" << Q-P
                     << "d=" << (Q-P).norm() << "\n";
                cout << "tP=" << R_SP.x() <<
                     "(Q-P).tP=" << ~(Q-P)*R_SP.x() << "\n";
                if (x[1] < -SqrtEps) {
                    xadj[0] -= Pi; xadj[1] = -xadj[1];  // flip
                    flipped = true;
                } else xadj[1]=0; // ignore
            }
            Fx = calcOrthogonalGeodError(P, Q, xadj[0], xadj[1],geod);
            if (flipped) Fx[1] = -Fx[1];
            f = std::sqrt(~Fx*Fx);
            dist = (geod.getPointQ()-Q).norm();
            //cout << "step size=" << lam << " errNorm=" << f << " dist=" << dist << endl;
            if (f > fold && lam > minlam) {
                lam = lam / 2;
            } else {
                break;
            }
        }
        if (f > ftol && maxabsdiff(x,xold) < xtol) {
            std::cout << "ORTHO terminated on too-small step size after " 
                      << i << " iterations" << std::endl;
            std::cout << "err = " << Fx << std::endl;
            break;
        }

        if (vizReporter != NULL) {
            vizReporter->handleEvent(ptOnSurfSys->getDefaultState());
            sleepInSec(pauseBetweenGeodIterations);
        }

    }
    if (f > ftol)
        std::cout << "### ORTHO geodesic DIVERGED in " 
                    << MaxIterations << " iterations with err=" << f << std::endl;

    // Finish each geodesic with reverse Jacobi field.
#ifdef USE_NEW_INTEGRATOR
    calcGeodesicReverseSensitivity2(geod, Vec2(0,1), Vec2(1,0));
#else
    calcGeodesicReverseSensitivity(geod, Vec2(0,1));
#endif

}


void ContactGeometryImpl::shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {
    std::cout << "warning: no analytical shootGeodesic for ContactGeometry base class, computing numerically." << std::endl;
    shootGeodesicInDirectionUntilLengthReached(xP, tP, terminatingLength, options, geod);
}

void ContactGeometryImpl::shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {
    std::cout << "warning: no analytical shootGeodesic for ContactGeometry base class, computing numerically." << std::endl;
    shootGeodesicInDirectionUntilPlaneHit(xP, tP, terminatingPlane, options, geod);
}

void ContactGeometryImpl::calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
            const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {
    std::cout << "warning: no analytical calcGeodesic for ContactGeometry base class, computing numerically." << std::endl;
    calcGeodesic(xP, xQ, tPhint, tQhint, geod);
}

// Calculate the "geodesic error" for thetaP and thetaQ, and return the
// resulting (kinked) geodesic if the supplied pointer is non-null.
Vec2 ContactGeometryImpl::
calcSplitGeodError(const Vec3& xP, const Vec3& xQ,
              Real thetaP, Real thetaQ,
              Geodesic* geodesic) const
{
    UnitVec3 tP = calcUnitTangentVec(thetaP, R_SP);
    UnitVec3 tQ = calcUnitTangentVec(thetaQ, R_SQ);
    return calcSplitGeodError(xP, xQ, tP, tQ, geodesic);
}

// Calculate the "geodesic error" for tP and tQ
Vec2 ContactGeometryImpl::
calcSplitGeodError(const Vec3& xP, const Vec3& xQ,
              const UnitVec3& tP, const UnitVec3& tQ,
              Geodesic* geodesic) const
{
    geodP.clear();
    geodQ.clear();

    GeodesicOptions opts;
    shootGeodesicInDirectionUntilPlaneHit(xP, tP,
            geodHitPlaneEvent->getPlane(), opts, geodP);
    shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
            geodHitPlaneEvent->getPlane(), opts, geodQ);

    // Finish each geodesic with reverse Jacobi field.
    calcGeodesicReverseSensitivity(geodP,
        geodQ.getDirectionalSensitivityPtoQ().back());
    calcGeodesicReverseSensitivity(geodQ,
        geodP.getDirectionalSensitivityPtoQ().back());

    if (geodesic)
        mergeGeodesics(geodP, geodQ, *geodesic);

    return calcError(geodP, geodQ);
}

//------------------------------------------------------------------------------
//                      CALC ORTHOGONAL GEOD ERROR
//------------------------------------------------------------------------------
// Calculate the "orthogonal error" for thetaP with given arc length.
// The error is
//           err(theta,s) = [ dot(Qhat-Q, bQhat) ]
//                          [ dot(Qhat-Q, tQhat) ]
Vec2 ContactGeometryImpl::
calcOrthogonalGeodError(const Vec3& xP, const Vec3& xQ,
                        Real thetaP, Real length,
                        Geodesic& geod) const
{
    const UnitVec3 tP = calcUnitTangentVec(thetaP, R_SP);

    geod.clear();

    GeodesicOptions opts;

    Vec3 Qhat;
    if (length < 1e-3) { // TODO
        Qhat = xP + length*tP;
        this->makeStraightLineGeodesic(xP, Qhat, tP, opts, geod);
    } else {
        shootGeodesicInDirectionUntilLengthReached(xP, tP, length, opts, geod);
        Qhat = geod.getPointQ();
    }
    const Vec3 r_QQhat = Qhat - xQ;

    const Real eb = ~r_QQhat * geod.getBinormalQ(); // length errors
    const Real es = ~r_QQhat * geod.getTangentQ();

    return Vec2(eb, es);
}

// Calculate the "geodesic error" for tP and tQ
Vec2 ContactGeometryImpl::
calcSplitGeodErrorAnalytical(const Vec3& xP, const Vec3& xQ,
              const UnitVec3& tP, const UnitVec3& tQ,
              Geodesic* geodesic) const
{
    geodP.clear();
    geodQ.clear();

    GeodesicOptions opts;
    shootGeodesicInDirectionUntilPlaneHitAnalytical(xP, tP,
            geodHitPlaneEvent->getPlane(), opts, geodP);
    shootGeodesicInDirectionUntilPlaneHitAnalytical(xQ, tQ,
            geodHitPlaneEvent->getPlane(), opts, geodQ);

    // Finish each geodesic with reverse Jacobi field.
    calcGeodesicReverseSensitivity(geodP,
        geodQ.getDirectionalSensitivityPtoQ().back());
    calcGeodesicReverseSensitivity(geodQ,
        geodP.getDirectionalSensitivityPtoQ().back());

    if (geodesic)
        mergeGeodesics(geodP, geodQ, *geodesic);

    return calcError(geodP, geodQ);
}



// Calculate the "geodesic jacobian" by numerical perturbation
Mat22 ContactGeometryImpl::
calcSplitGeodErrorJacobian(const Vec3& xP, const Vec3& xQ,
        const Real& thetaP, const Real& thetaQ, Differentiator::Method order) const {

//    UnitVec3 tP = calcUnitTangentVecGG(thetaP, R_SP);
//    UnitVec3 tQ = calcUnitTangentVecGG(thetaQ, R_SQ);
//
//    geodP.clear();
//    geodQ.clear();
//
//    GeodesicOptions opts;
//    shootGeodesicInDirectionUntilPlaneHit(xP, tP,
//            geodHitPlaneEvent->getPlane(), opts, geodP);
//    shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
//            geodHitPlaneEvent->getPlane(), opts, geodQ);

    GeodesicOptions opts;

    UnitVec3 tP, tQ;

    Vec2 fy0 = calcError(geodP, geodQ);

    Mat22 dfdy;
    Vec2 fyptmp, fymtmp;
    Real hEst, h, accFactor, YMin;

    Geodesic geodPtmp;
    Geodesic geodQtmp;

    if (order == 1) {
        accFactor = std::sqrt(estimatedGeodesicAccuracy);
    } else {
        accFactor = std::pow(estimatedGeodesicAccuracy, OneThird);
    }
    YMin = Real(0.1);

    Vec2 y0(thetaP, thetaQ);

    // perturb thetaP
    hEst = accFactor*std::max(std::abs(thetaP), YMin);
    h = cleanUpH(hEst, thetaP);

    // positive perturb
    tP = calcUnitTangentVec(thetaP+h, R_SP);
    shootGeodesicInDirectionUntilPlaneHit(xP, tP,
            geodHitPlaneEvent->getPlane(), opts, geodPtmp);
    fyptmp = calcError(geodPtmp, geodQ);

    if (order==1) {
        dfdy(0) = (fyptmp-fy0)/h;
    } else {

        geodPtmp.clear();
        tP = calcUnitTangentVec(thetaP-h, R_SP);
        shootGeodesicInDirectionUntilPlaneHit(xP, tP,
                geodHitPlaneEvent->getPlane(), opts, geodPtmp);
        fymtmp = calcError(geodPtmp, geodQ);

        dfdy(0) = (fyptmp-fymtmp)/(2*h);
    }

    // perturb thetaQ
    hEst = accFactor*std::max(std::abs(thetaQ), YMin);
    h = cleanUpH(hEst, thetaP);

    // positive perturb
    tQ = calcUnitTangentVec(thetaQ+h, R_SQ);
    shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
            geodHitPlaneEvent->getPlane(), opts, geodQtmp);
    fyptmp = calcError(geodP, geodQtmp);

    if (order==1) {
        dfdy(1) = (fyptmp-fy0)/h;
    } else {

        geodQtmp.clear();
        tQ = calcUnitTangentVec(thetaQ-h, R_SQ);
        shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
                geodHitPlaneEvent->getPlane(), opts, geodQtmp);
        fymtmp = calcError(geodP, geodQtmp);

        dfdy(1) = (fyptmp-fymtmp)/(2*h);
    }

    return dfdy;
}

Vec2  ContactGeometryImpl::
calcError(const Geodesic& geodP, const Geodesic& geodQ) const {
    const Transform& Fphat = geodP.getFrenetFrames().back();
    const Transform& Fqhat = geodQ.getFrenetFrames().back();
    const Vec3&      Phat = Fphat.p(); const Vec3&      Qhat = Fqhat.p();
    const UnitVec3& tPhat = Fphat.x(); const UnitVec3& tQhat = Fqhat.x();
    const UnitVec3& bPhat = Fphat.y(); const UnitVec3& bQhat = Fqhat.y();
    const UnitVec3& nPhat = Fphat.z(); const UnitVec3& nQhat = Fqhat.z();

    // Error is separation distance along mutual b direction, and angle by
    // which the curves fail to connect smoothly.
    Vec2 geodErr(~(bPhat - bQhat) * (Phat - Qhat), ~bPhat * tQhat);
    return geodErr;
}

