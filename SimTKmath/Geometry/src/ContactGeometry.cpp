/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman, Ian Stavness, Andreas Scholz      *
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

#include "ContactGeometryImpl.h"

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

Real ContactGeometry::calcGaussianCurvature(const Vec3& point) const {
    return getImpl().calcGaussianCurvature(point);
}

Real ContactGeometry::calcSurfaceCurvatureInDirection(const Vec3& point, const UnitVec3& direction) const {
	return getImpl().calcSurfaceCurvatureInDirection(point, direction);
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
    const UnitVec3& x1 = R_SP1.x(); // P1 kmax direction
    const UnitVec3& y1 = R_SP1.y(); // P1 kmin direction
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
// Nx  := dn/dx = dn/dg*dg/dx = H*(I- n*(~n))/norm(g)
// Tx  := dt/dx = -Nx*(~tP*n) -n*(~tP*Nx)
// Bx  := db/dx = -Nx*(~bP*n) -n*(~bP*Nx)

class ProjToNearestPointJacobianAnalytical {
public:
    ProjToNearestPointJacobianAnalytical(const ContactGeometryImpl& geom, const Vec3& p,
            const Vec3& tP, const Vec3& bP) :  geom(geom), p(p), tP(tP), bP(bP) { }
    void f (const Vec3& x, Vec3& f) {
        UnitVec3 n = geom.calcSurfaceUnitNormal(x);

        // calculate frame at x
        Vec3 r = p-x;
        Vec3 t = tP-n*(~tP*n); // project tP to tangent plane at x
        Vec3 b = bP-n*(~bP*n); // project bP to tangent plane at x

        SimTK_ASSERT_ALWAYS(t.norm() > 1e-6, "t is aligned with the normal vector at the current point.");
        SimTK_ASSERT_ALWAYS(b.norm() > 1e-6, "b is aligned with the with normal vector at the current point.");
        // Note that t and b are not perpendicular in this scheme, but do span the tangent space at x
//        SimTK_ASSERT_ALWAYS(~t*b < 1e-6, "t and b are not perpendicular.");

        // calculate error
        f[0] = geom.calcSurfaceValue(x);
        f[1] = ~r*t;
        f[2] = ~r*b;
    }

    void J (const Vec3& x, Mat33& J) {
        Vec3 g = geom.calcSurfaceGradient(x); // non-normalized outward facing "normal" vector
        Mat33 H = geom.calcSurfaceHessian(x);
        Real gNorm = g.norm();
        UnitVec3 n(-g/gNorm, true);

        // calculate frame at x
        Vec3 r = p-x;
        Vec3 t = tP-n*(~tP*n); // project tP to tangent plane at x
        Vec3 b = bP-n*(~bP*n); // project bP to tangent plane at x

        SimTK_ASSERT_ALWAYS(t.norm() > 1e-6, "t is aligned with the normal vector at the current point.");
        SimTK_ASSERT_ALWAYS(b.norm() > 1e-6, "b is aligned with the with normal vector at the current point.");
        // Note that t and b are not perpendicular in this scheme, but do span the tangent space at x
//        SimTK_ASSERT_ALWAYS(~t*b < 1e-6, "t and b are not perpendicular.");

        // calculate derivatives
        // TODO fix analytical jacobian
        Mat33 I(1);
        Mat33 Nx = (1/gNorm)*H*(I - n*(~n));
        Mat33 Tx = -Nx*(~tP*n) - n*(~tP*Nx);
        Mat33 Bx = -Nx*(~bP*n) - n*(~bP*Nx);

//        cout << "n=" << n  << endl;
//        cout << "r=" << r << endl;
//        cout << "t=" << t << endl;
//        cout << "b=" << b << endl;
//        cout << "H=" << H << endl;
//        cout << "Nx=" << Nx << endl;
//        cout << "Tx=" << Tx << endl;
//        cout << "Bx=" << Tx << endl;
//
        J[0] = ~g;
        J[1] = -(~t) + ~r*Tx;
        J[2] = -(~b) + ~r*Bx;

    }

    const ContactGeometryImpl& geom;
    const Vec3& p;
    const Vec3& tP;
    const Vec3& bP;
};

class ProjToNearestPointJacobian : public Differentiator::JacobianFunction {
public:
    ProjToNearestPointJacobian(const ContactGeometryImpl& geom, const Vec3& p, const Vec3& tP, const Vec3& bP)
        : Differentiator::JacobianFunction(3,3), geom(geom), p(p), tP(tP), bP(bP) { }

    int f(const Vector& xvec, Vector& f) const {
        const Vec3& x = Vec3::getAs(&xvec[0]);
        UnitVec3 n = geom.calcSurfaceUnitNormal(x);
        Vec3 r = p-x;
        Vec3 t = tP-n*(~tP*n); // project tP to tangent plane at x
        Vec3 b = bP-n*(~bP*n); // project bP to tangent plane at x

        SimTK_ASSERT_ALWAYS(t.norm() > 1e-6, "t is aligned with the normal vector at the current point.");
        SimTK_ASSERT_ALWAYS(b.norm() > 1e-6, "b is aligned with the with normal vector at the current point.");
        // Note that t and b are not perpendicular in this scheme, but do span the tangent space at x
//        SimTK_ASSERT_ALWAYS(~t*b < 1e-6, "t and b are not perpendicular.");


        // calculate error
        f[0] = geom.calcSurfaceValue(x);
        f[1] = ~r*t;
        f[2] = ~r*b;

        return 0;
    }

    const ContactGeometryImpl& geom;
    const Vec3& p;
    const Vec3& tP;
    const Vec3& bP;
};

Vec3 ContactGeometryImpl::
projectDownhillToNearestPoint(const Vec3& Q) const {

    // Newton solver settings
    const Real ftol = 1e-14;

    // Check for immediate return.
    if (std::abs(calcSurfaceValue(Q)) <= ftol)
        return Q;

    // construct arbitrary frame at p
    UnitVec3 nP = calcSurfaceUnitNormal(Q);
    UnitVec3 tP = nP.perp();
    UnitVec3 bP(tP%nP);

    // Estimate a scale for the local neighborhood of this surface by 
    // using the larger curvature in the t or b direction. We want to take
    // conservative steps that never move by more than a fraction of the
    // scale to avoid jumping out of the local minimum.
    const Real kt = calcSurfaceCurvatureInDirection(Q, tP);
    const Real kb = calcSurfaceCurvatureInDirection(Q, bP);
    const Real maxK = std::max(std::abs(kt),std::abs(kb));
    const Real scale = std::min(1/maxK, 1000.); // keep scale reasonable
    const Real MaxMove = .25; // Limit one move to 25% of smaller radius.

    const Real xtol = 1e-12;
    const Real minlam = 1e-9;
    const int maxNewtonIterations = 30;

    Vec3 x(Q); // initialize to query point
    Vec3 f, dx, xold;
    Mat33 J, Jp;
    Vector ftmp(3);

//    Vec3 r, t, b; // temporary
//    Mat33 Nx, Tx, Bx; // temporary
//    Mat33 I(1);

    Real rmsError, rmsErrorOld, xchg, lam;

//    ProjToNearestPointJacobianAnalytical nearestPointJacAn(*this), p, tP, bP);

    ProjToNearestPointJacobian nearestPointJac(*this, Q, tP, bP);
    Differentiator diff(nearestPointJac);

    nearestPointJac.f((Vector)x, ftmp);
    f = Vec3::getAs(&ftmp[0]);
    rmsError = std::sqrt(f.normSqr()/3);
    //std::cout << "BEFORE Q=" << Q << ", f=" << f << ", frms=" << rmsError << std::endl;

    int cnt = 0;
    do {
        if (rmsError <= ftol) { // found solution
            //std::cout << "CONVERGED in " << cnt << " steps, frms=" << rmsError << std::endl;
            break;
        }

        Matrix Jnum = diff.calcJacobian((Vector)x);
        J = Mat33::getAs(&Jnum(0,0));

//        nearestPointJacAn.J(x, J);

        dx = J.invert()*f;
        const Real dxrms = std::sqrt(dx.normSqr()/3);

        rmsErrorOld = rmsError;
        xold = x;

        // Backtracking. Limit the starting step size if dx is too big.
        lam = std::min(1., MaxMove*(scale/dxrms));
        if (lam < 1) {
            //std::cout << "PROJECT: LIMITED STEP: iter=" << cnt 
            //          << " lam=" << lam << endl;
        }
        while (true) {
            x = xold - lam*dx;
            nearestPointJac.f((Vector)x, ftmp);
            f = Vec3::getAs(&ftmp[0]);

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

    int f(const Vector& xvec, Vector& f) const OVERRIDE_11 {
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
    const Real Ftol2 = square(1e-14);

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
    const Real MinStepFrac = 1e-6;
    // If the norm of a change to X during a step isn't at least this fraction
    // of a full-scale change (see scaling below), we'll treat that as an
    // independent reason to give up (squared for speed).
    const Real Xtol2 = square(1e-12);
    // We'll calculate a length scale for the local patch of this object and
    // then limit any moves we make to this fraction of that scale to avoid
    // jumping out of one local minimum into another. This can cause a series
    // of slowly-converging steps to be taken at the beginning.
    const Real MaxMove2 = square(.25); // Limit one move to 25% of smaller radius.

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
    const Real scale2 = square(clamp(0.1, 1/maxK, 1000.)); // keep scale reasonable

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
        const Real stepFrac2 = std::min(1., MaxMove2*(scale2/dxrms2));
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
            grad=Vec3(1,1.1,1.2); gradMag=grad.norm();
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
calcGaussianCurvature(const Vec3& point) const {
    const Vec3  g = calcSurfaceGradient(point);
    const Mat33 H = calcSurfaceHessian(point);
    // Calculate the adjoint. TODO: use SymMat33
    Mat33 A;
    A(0,0) = det(H.dropRowCol(0,0));
    A(0,1) = det(H.dropRowCol(0,1));
    A(0,2) = det(H.dropRowCol(0,2));
    A(1,0) = A(0,1);
    A(1,1) = det(H.dropRowCol(1,1));
    A(1,2) = det(H.dropRowCol(1,2));
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);
    A(2,2) = det(H.dropRowCol(2,2));

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
    const UnitVec3 n = calcSurfaceUnitNormal(point);
    const Vec3     g = calcSurfaceGradient(point);
	const Mat33 H = calcSurfaceHessian(point);
	const Real  knum = ~direction*H*direction; // numerator
    if (std::abs(knum) < TinyReal)
        return 0; // don't want to return 0/0.
        
    const Real k = knum/(~g*n);

	return k;
}





static const Real estimatedGeodesicAccuracy = 1e-12; // used in numerical differentiation
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
    int f(const Vector& x, Vector& fx) const  {
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
    int f(const Vector& x, Vector& fx) const  {
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
    const Vec3 P = projectDownhillToNearestPoint(xP);
    const Vec3 Q = projectDownhillToNearestPoint(xQ);

    // TODO: If no previous geodesic just make a wild attempt.
    if (!prevGeod.getNumPoints()) {
        const Vec3 PQ = Q-P;
        const Real length = PQ.norm();
        if (length < /*SqrtEps*/1e-3) { // TODO: should depend on curvature
            const UnitVec3 d = length<TinyReal ? UnitVec3(XAxis)
                : UnitVec3(PQ/length, true);
            makeStraightLineGeodesic(P, Q, d, options, geod);
            return;
        }
        calcGeodesicUsingOrthogonalMethod(P, Q, PQ/length, length, geod);
        //calcGeodesicAnalytical(P, Q, PQ, PQ, geod);
        return;
    }


    // First classify the previous geodesic as direct or indirect. Direct is
    // a strict classification; only if the end tangents are aligned with the
    // PQ line to within an allowed cone angle is it direct.
    const Real CosMaxDirectAngle = 0.9; // about 25 degrees

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
    const bool isPrevStraightLine = prevPQlen*maxK <= SqrtEps;

    const UnitVec3 eprevPQ = isPrevStraightLine
                                ? prevGeod.getTangentP()
                                : UnitVec3(prevPQ/prevPQlen, true);

    const Real cosConetP = dot(prevGeod.getTangentP(), eprevPQ);
    const Real cosConetQ = dot(prevGeod.getTangentQ(), eprevPQ);
    const bool isDirect = std::min(cosConetP,cosConetQ) >= CosMaxDirectAngle;

    UnitVec3 tPhint = prevGeod.getTangentP(); // might flip
    UnitVec3 tQhint = prevGeod.getTangentQ();
    Real     sHint  = prevGeod.getLength();

    const Vec3 newPQ = Q - P;
    const Real newPQlen = newPQ.norm();
    if (isDirect) {
        const UnitVec3 enewPQ = newPQlen > 0 
            ? UnitVec3(newPQ/newPQlen, true)
            : UnitVec3(prevGeod.getTangentP()+prevGeod.getTangentQ());

        if (~tPhint*enewPQ < 0 && ~tQhint*enewPQ < 0) {
            tPhint = -tPhint;
            tQhint = -tQhint;
            cout << "GEODESIC FLIPPED. Prev len was " << sHint << endl;
        }
    }


    //calcGeodesicUsingOrthogonalMethod(P, Q, tPhint, sHint, geod);
    //calcGeodesic(P, Q, tPhint, tQhint, geod);
    if (newPQlen < /*SqrtEps*/1e-3) //TODO
        makeStraightLineGeodesic(P, Q, tPhint, options, geod);
    else 
        //calcGeodesicAnalytical(P, Q, tPhint, tQhint, geod);
        calcGeodesicUsingOrthogonalMethod(P, Q, tPhint, sHint, geod);
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
    //TODO:
    Vec3 Pprime = xP;
        // projectToNearestDownhillPointOnSurface(xP);
    Vec3 Qprime = xQ;
        // projectToNearestDownhillPointOnSurface(xQ);

    const bool isZeroLength = 
        Geo::Point::pointsAreNumericallyCoincident(Pprime,Qprime);

    UnitVec3 d;
    Real     length;
    if (isZeroLength) {
        Pprime = Qprime = Geo::Point::findMidpoint(Pprime, Qprime);
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
void ContactGeometryImpl::
shootGeodesicInDirection(const Vec3& P, const UnitVec3& tP,
        const Real& finalTime, const GeodesicOptions& options,
        Geodesic& geod) const {

    // integrator settings
    const Real startTime = 0;
    const Real integratorAccuracy = 1e-6; // << TODO
    const Real integratorConstraintTol = 1e-6;

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
    integ.setAccuracy(integratorAccuracy);
    integ.setConstraintTolerance(integratorConstraintTol);
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
    geod.setAchievedAccuracy(integratorAccuracy); // TODO: accuracy of length?
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
    const Real integratorAccuracy = 1e-6;
    const Real integratorConstraintTol = 1e-6;

    //RungeKutta3Integrator integ(ptOnSurfSys);
    RungeKuttaMersonIntegrator integ(*ptOnSurfSys);
    integ.setAccuracy(integratorAccuracy);
    integ.setConstraintTolerance(integratorConstraintTol);
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
    shootGeodesicInDirection(xP, tP, terminatingLength, options, geod);
}



static Real cleanUpH(Real hEst, Real y0) {
    volatile Real temp = y0+hEst;
    return temp-y0;
}

static Real maxabs(Vec2 x) {
    return std::max(std::abs(x[0]), std::abs(x[1]));
}

static Real maxabsdiff(Vec2 x, Vec2 xold) {
    return std::max(std::abs(x[0]-xold[0])/std::max(x[0],1.0),
                    std::abs(x[1]-xold[1])/std::max(x[1],1.0));
}



//------------------------------------------------------------------------------
//                              CALC GEODESIC
//------------------------------------------------------------------------------
// Utility method to find geodesic between P and Q
// with starting directions tPhint and tQhint
// XXX tangent basis should be formed from previous geodesic
void ContactGeometryImpl::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {

    // Newton solver settings
    const Real ftol = 1e-9;
    const Real xtol = 1e-9;
    const Real minlam = 1e-9;
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

    // Newton solver settings
    const Real ftol = 1e-9;
    const Real xtol = 1e-12;
    const Real minlam = 1e-3;
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
    //orthoErr.setEstimatedAccuracy(1e-16); // TODO
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
            //std::cout << "ORTHO geodesic converged in " 
             //         << i << " iterations with err=" << f << std::endl;
            break;
        }

        if (useNewtonIteration) {
            // This numerical Jacobian is very bad; CentralDifference is 
            // required in order to produce a reasonable one.
            //diff.calcJacobian(Vector(x),  Vector(Fx), JMat, 
            //                  Differentiator::CentralDifference);
            //J = Mat22::getAs(&JMat(0,0));

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
            cout << "ORTHO: lam reduced to " << lam << "\n";
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
    calcGeodesicReverseSensitivity(geod, Vec2(0,1));

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



//==============================================================================
//                             HALF SPACE & IMPL
//==============================================================================
ContactGeometry::HalfSpace::HalfSpace()
:   ContactGeometry(new HalfSpace::Impl()) {}

/*static*/ ContactGeometryTypeId ContactGeometry::HalfSpace::classTypeId() 
{   return ContactGeometry::HalfSpace::Impl::classTypeId(); }

DecorativeGeometry ContactGeometry::HalfSpace::Impl::createDecorativeGeometry() const {
    return DecorativeBrick(Vec3(0.01,1,1));
}

// Point position is given in the half space frame.
Vec3 ContactGeometry::HalfSpace::Impl::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    inside = (position[0] >= 0);
    normal = -UnitVec3(XAxis); // this does not require normalization
    return Vec3(0, position[1], position[2]);
}

bool ContactGeometry::HalfSpace::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    if (std::abs(direction[0]) < SignificantReal)
        return false; // ray is parallel to halfspace surface

    const Real t = origin[0]/direction[0];
    if (t > 0)
        return false; // ray points away from surface

    distance = -t;
    normal = -UnitVec3(XAxis); // cheap; no normalization required
    return true;
}

void ContactGeometry::HalfSpace::Impl::getBoundingSphere
   (Vec3& center, Real& radius) const 
{   center = Vec3(0);
    radius = Infinity; }

const ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
getImpl() const {
    assert(impl);
    return static_cast<const HalfSpace::Impl&>(*impl);
}

ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
updImpl() {
    assert(impl);
    return static_cast<HalfSpace::Impl&>(*impl);
}



//==============================================================================
//                            CYLINDER & IMPL
//==============================================================================

ContactGeometry::Cylinder::Cylinder(Real radius)
:   ContactGeometry(new Cylinder::Impl(radius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Cylinder::classTypeId()
{   return ContactGeometry::Cylinder::Impl::classTypeId(); }

Real ContactGeometry::Cylinder::getRadius() const {
    return getImpl().getRadius();
}

void ContactGeometry::Cylinder::setRadius(Real radius) {
    updImpl().setRadius(radius);
}

const ContactGeometry::Cylinder::Impl& ContactGeometry::Cylinder::getImpl() const {
    assert(impl);
    return static_cast<const Cylinder::Impl&>(*impl);
}

ContactGeometry::Cylinder::Impl& ContactGeometry::Cylinder::updImpl() {
    assert(impl);
    return static_cast<Cylinder::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Cylinder::Impl::createDecorativeGeometry() const {
    DecorativeCylinder cyl(radius, radius*2);
    // DecorativeCylinder's axis is defined as the y-axis,
    // whereas ContactGeometry::Cylinder axis is defined as the z-axis
    cyl.setTransform(Rotation(UnitVec3(0, 1, 0), ZAxis, Vec3(0, 0, 1), YAxis));
    return cyl;
}

Vec3 ContactGeometry::Cylinder::Impl::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {

    normal = calcSurfaceUnitNormal(position);

    // long axis is z-axis, project to x-y plane
    Vec2 xy_position(position(0), position(1));
    inside = (xy_position.normSqr() <= radius*radius);

    // nearestPoint = point_on_surface_in_xy_plane + height_in_z
    Vec3 nearestPoint = normal*radius + Vec3(0,0,position(2));

    return nearestPoint;
}

bool ContactGeometry::Cylinder::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    // cylinder axis is z-axis, project to x-y plane
    const Vec3 xy_vec(direction(0),direction(1),0);
    const Real xy_vec_norm = xy_vec.norm();
    const UnitVec3 xy_direction(xy_vec/xy_vec_norm, true); // don't renormalize
    const Vec3 xy_origin(origin(0),origin(1),0);
    Real xy_distance;
    Real b = -~xy_direction*xy_origin;
    Real c = xy_origin.normSqr() - radius*radius;
    if (c > 0) {
        // Ray origin is outside cylinder.

        if (b <= 0)
          return false;  // Ray points away from axis of cylinder.
        Real d = b*b - c;
        if (d < 0)
          return false;
        Real root = std::sqrt(d);
        xy_distance = b - root;
      }
    else {
        // Ray origin is inside cylinder.

        Real d = b*b - c;
        if (d < 0)
          return false;
        xy_distance = b + std::sqrt(d);
      }
    distance = xy_distance/xy_vec_norm;
    normal = UnitVec3(xy_origin+xy_distance*xy_direction);
    return true;
}

void ContactGeometry::Cylinder::Impl::getBoundingSphere
    (Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = Infinity;
}

void ContactGeometry::Cylinder::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {

    // long axis (direction of min curvature) points in <0,0,1>
    orientation = Rotation(calcSurfaceUnitNormal(point), ZAxis, Vec3(0, 0, 1), YAxis);
    curvature[0] = 1/radius;
    curvature[1] = 0;

}

// Sample geodesic between two points P and Q on a cylinder analytically.
static void setGeodesicToHelicalArc(Real R, Real phiP, Real angle, Real m, Real c, Geodesic& geod)
{
   	// Clear current geodesic.
	geod.clear();

    const Real sqrt1m2 = sqrt(1+m*m);   // Avoid repeated calculation.
    const Real kappa = 1 / (R*(1+m*m)); // Curvature in tangent direction.
    const Real kb = 1 / (R*(1+1/(m*m))); // Curvature in binormal direction
                                         //   (slope is 1/m).
	const Real tau   = m*kappa;         // Torsion (signed).

	// Arc length of the helix. Always
	const Real L = R * sqrt1m2 * std::abs(angle);

	// Orientation of helix. 
	const Real orientation = angle < 0 ? Real(-1) : Real(1);

	// TODO: Make this generic, so long geodesics are sampled more than short ones.
    const int numGeodesicSamples = 12;
	const Real deltaPhi = std::abs(angle / Real(numGeodesicSamples-1));

    for (int i = 0; i < numGeodesicSamples; ++i)
	{
		// Watch out: Angle phi has an offset phiP
        Real phi = Real(i)*angle/Real(numGeodesicSamples-1) + phiP;
        const Real sphi = sin(phi), cphi = cos(phi);

		// Evaluate helix.
        Vec3	 p( R*cphi, R*sphi, R*m*(phi - phiP) + c);

        // We'll normalize so UnitVec3 doesn't have to do it.
		UnitVec3 t((orientation/sqrt1m2)*Vec3(-sphi, cphi, m), true);
		UnitVec3 n(Vec3(cphi, sphi, 0), true);

        // Though not needed, we use an orthogonalizing constructor for the rotation.
        geod.addFrenetFrame(Transform(Rotation(n, ZAxis, t, YAxis), p));

		// Current arc length s.
		Real s = R * sqrt1m2 * (Real(i)*deltaPhi);
        geod.addArcLength(s);
		geod.addCurvature(kappa);		

		// Solve the scalar Jacobi equation
		//
		//        j''(s) + K(s)*j(s) = 0 ,                                     (1)
		//
		// where K is the Gaussian curvature and (.)' := d(.)/ds denotes differentiation
		// with respect to the arc length s. Then, j is the directional sensitivity and
		// we obtain the corresponding variational vector field by multiplying b*j. For
		// a cylinder, K = 0 and the solution of equation (1) becomes
		//
		//        j  = s				                                       (2)
		//		  j' = 1 ,							                           (3)
		//
		// so the Jacobi field increases linearly in s.

		// Forward directional sensitivity from P to Q
		Vec2 jPQ(s, 1);
		geod.addDirectionalSensitivityPtoQ(jPQ);

		// Backwards directional sensitivity from Q to P
		Vec2 jQP(L-s, 1);
		geod.addDirectionalSensitivityQtoP(jQP);
    }

	// Only compute torsion and binormal curvature at the end points.
	geod.setTorsionAtP(tau); geod.setTorsionAtQ(tau);
    geod.setBinormalCurvatureAtP(kb); geod.setBinormalCurvatureAtQ(kb);

    geod.setIsConvex(true); // Curve on cylinder is always convex.

    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(1e-15); // TODO: accuracy of length?
//    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken()); // TODO
}

// Compute geodesic between two points P and Q on a cylinder analytically. Since a geodesic on a
// cylinder is a helix it is parameterized by
//
//			       [ R * cos(phi)    ]
//        p(phi) = [ R * sin(phi)    ]
//				   [ R * m * phi + c ]
//
// where R is the radius of the cylinder, phi parameterizes the opening angle of the helix, m is  
// the slope and c is an offset. We define the geodesic from P to Q, hence c = Pz.
void ContactGeometry::Cylinder::Impl::
calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
                       const Vec3& tPhint, const Vec3& tQhint,
                       Geodesic& geod) const
{
	// Compute angle between P and Q. Save both the positive (right handed) and the negative
	// (left handed) angle.
	Real phiP = atan2(xP[1], xP[0]);
	Real phiQ = atan2(xQ[1], xQ[0]);

	Real temp = phiQ - phiP;
	Real angleRightHanded, angleLeftHanded;

	// Left-handed angle will always be negative, right-handed angle will be positive.
	if (temp >= 0) {
		angleRightHanded = temp;
		angleLeftHanded  = temp - 2*Pi;
	}

	else {
		angleLeftHanded  = temp;
		angleRightHanded = temp + 2*Pi;
	}

	// Compute "moment" of tPhint at P and tQhint at Q around z-Axis.
	// Make sure tPhint and tQhint are unit vectors, otherwise moments are scaled.
	Real MP = xP[0]*tPhint[1] - xP[1]*tPhint[0];
	Real MQ = xQ[0]*tQhint[1] - xQ[1]*tQhint[0];

	// Average moment.
	Real M = (MP + MQ) / 2;

	// Decide whether helix is right or left handed. The sign of angle stores the 
	// information about the orientation (right handed, if positive)
	Real angle;
	if (M >= 0)	{
		angle = angleRightHanded;
	}

	else {
		angle = angleLeftHanded;
	}

	// Offset and slope.
	Real c =  xP[2];
	Real m = (xQ[2] - xP[2]) / (angle * radius);

	setGeodesicToHelicalArc(radius, phiP, angle, m, c, geod);
}

void ContactGeometry::Cylinder::Impl::shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {

    //TODO for Andreas :)
}

void ContactGeometry::Cylinder::Impl::shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {

    //TODO for Andreas :)
}

Real CylinderImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(x[0]*x[0]+x[1]*x[1])/square(ownerp->getRadius());
}

Real CylinderImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    if (derivComponents.size() == 1 && derivComponents[0] < 2)
        return -2*x[derivComponents[0]]/square(ownerp->getRadius());
    if (derivComponents.size() == 2 &&
        derivComponents[0] == derivComponents[1] &&
        derivComponents[0] < 2 )
        return -2/square(ownerp->getRadius());
    return 0;
}


//==============================================================================
//                               SPHERE & IMPL
//==============================================================================

ContactGeometry::Sphere::Sphere(Real radius)
:   ContactGeometry(new Sphere::Impl(radius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Sphere::classTypeId()
{   return ContactGeometry::Sphere::Impl::classTypeId(); }

Real ContactGeometry::Sphere::getRadius() const {
    return getImpl().getRadius();
}

void ContactGeometry::Sphere::setRadius(Real radius) {
    updImpl().setRadius(radius);
}

const ContactGeometry::Sphere::Impl& ContactGeometry::Sphere::getImpl() const {
    assert(impl);
    return static_cast<const Sphere::Impl&>(*impl);
}

ContactGeometry::Sphere::Impl& ContactGeometry::Sphere::updImpl() {
    assert(impl);
    return static_cast<Sphere::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Sphere::Impl::createDecorativeGeometry() const {
    return DecorativeSphere(radius);
}

Vec3 ContactGeometry::Sphere::Impl::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    inside = (position.normSqr() <= radius*radius);
    normal = UnitVec3(position); // expensive -- normalizing
    return normal*radius;
}

bool ContactGeometry::Sphere::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    Real b = -~direction*origin;;
    Real c = origin.normSqr() - radius*radius;
    if (c > 0) {
        // Ray origin is outside sphere.

        if (b <= 0)
          return false;  // Ray points away from center of sphere.
        Real d = b*b - c;
        if (d < 0)
          return false;
        Real root = std::sqrt(d);
        distance = b - root;
      }
    else {
        // Ray origin is inside sphere.

        Real d = b*b - c;
        if (d < 0)
          return false;
        distance = b + std::sqrt(d);
      }
    normal = UnitVec3(origin+distance*direction);
    return true;
}

void ContactGeometry::Sphere::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = this->radius;
}

void ContactGeometry::Sphere::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    orientation = Rotation(UnitVec3(point), ZAxis, fabs(point[0]) > 0.5 ? Vec3(0, 1, 0) : Vec3(1, 0, 0), XAxis);
    curvature = 1/radius;
}

//TODO: just an axis-aligned leaf box for now
void ContactGeometry::Sphere::Impl::createOBBTree() {
    OBBNode& root = obbTree.updRoot();
    root.box.setHalfLengths(Vec3(radius));
    root.normal = UnitVec3(XAxis);  // doesn't matter
    root.coneHalfAngle = Pi;        // has all possible normals
    root.pointOnSurface = Vec3(radius,0,0); // doesn't matter
    root.children.clear(); // This is a leaf

    // Leaf contents.
    root.centerUW = Vec2(0,0);
    root.dims = Vec2(Pi, Pi/2); // u in [-Pi,Pi], v in [-Pi/2,Pi/2]
}


// Compute geodesic between two points P and Q on a sphere analytically. Since a geodesic on a
// sphere is a great circle it is parameterized by
//
//        p(phi) = R * (e1*cos(phi) + e2*sin(phi)) ,
//
// where R is the radius of the sphere and the angle phi parameterizes the great circle with
// respect to an orthonormal basis {e1, e2}. By definition P = p(0) and the geodesic goes from
// P to Q, where Q = p(angle). Make sure e1 . e2 = 0 and |e1| = |e2| = 1.
static void setGeodesicToArc(const UnitVec3& e1, const UnitVec3& e2,
                             double R, double angle, Geodesic& geod)
{
    // Check if e1 and e2 are orthogonal.
    assert(std::abs(~e1*e2) <= SignificantReal);

	// Clear current geodesic.
	geod.clear();

	// TODO: Make this generic, so long geodesics are sampled more than short ones.
    const int numGeodesicSamples = 12;

	// Total arc length and orientation.
	const Real orientation = sign(angle);
	const Real L = R*angle*orientation;

	// Increment of phi in loop.
	const Real deltaPhi = std::abs(angle / Real(numGeodesicSamples-1));

    const Real k = 1/R; // curvature
    for (int i = 0; i < numGeodesicSamples; ++i){
        Real phi = Real(i)*angle / Real(numGeodesicSamples-1);
        const Real sphi = sin(phi), cphi = cos(phi);

		// Trust me, this is already normalized by definition of the input.
		UnitVec3 n(e1*cphi + e2*sphi, true);

        Vec3 p = R*n;

		// t = dp/dphi, hence pointing into direction of increasing phi. 
        Vec3 t = (-e1*sphi + e2*cphi)*orientation;

        // Though not needed, we use an orthogonalizing constructor for the rotation.
        geod.addFrenetFrame(Transform(Rotation(n, ZAxis, t, YAxis), p));

		// Current arc length s.
		Real s = R*Real(i)*deltaPhi;
        geod.addArcLength(s);

		// Solve the scalar Jacobi equation
		//
		//        j''(s) + K(s)*j(s) = 0 ,                                     (1)
		//
		// where K is the Gaussian curvature and (.)' := d(.)/ds denotes differentiation
		// with respect to the arc length s. Then, j is the directional sensitivity and
		// we obtain the corresponding variational vector field by multiplying b*j. For
		// a sphere, K = R^(-2) and the solution of equation (1) becomes
		//
		//        j  = R * sin(1/R * s)                                        (2)
		//		  j' =     cos(1/R * s) ,                                      (3)
		//
		// where equation (2) is the standard solution of a non-damped oscillator. Its
		// period is 2*pi*R and its amplitude is R.

		// Forward directional sensitivity from P to Q
		Vec2 jPQ(R*sin(k * s), cos(k * s));
		geod.addDirectionalSensitivityPtoQ(jPQ);

		// Backwards directional sensitivity from Q to P
		Vec2 jQP(R*sin(k * (L-s)), cos(k * (L-s)));
		geod.addDirectionalSensitivityQtoP(jQP);

        geod.addCurvature(k);
    }
    geod.setTorsionAtP(0); geod.setTorsionAtQ(0);
    geod.setBinormalCurvatureAtP(k); geod.setBinormalCurvatureAtQ(k);

    geod.setIsConvex(true); // Curve on sphere is always convex.
    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(1e-15); // TODO: accuracy of length?
//    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken()); // TODO
}


void ContactGeometry::Sphere::Impl::
calcGeodesicAnalytical(const Vec3& xP, const Vec3& xQ,
                       const Vec3& tPhint, const Vec3& tQhint,
                       Geodesic& geod) const
{
	// Build an orthonormal basis {e1, e2, e3}.
    const UnitVec3 e1(xP), e_OQ(xQ);
    const Vec3 arcAxis = e1 % e_OQ;

    const Real sinAngle = arcAxis.norm();
    const Real cosAngle = ~e1*e_OQ;

    UnitVec3 e3(arcAxis/sinAngle, true);

	// Tangent vectors tP and tQ at P and Q corresponding to a positive rotation
	// of the arc around e3.
    UnitVec3 tP(e3 % e1,   true);
    UnitVec3 tQ(e3 % e_OQ, true);

	// Average moment of of hint vectors applied e3.
	Real MP = ~(e1   % tPhint)*e3;
	Real MQ = ~(e_OQ % tQhint)*e3;
	Real M  =  (MP + MQ) / 2;

	// Small angle between e_OP and e_OQ corresponding to a short geodesic.
    Real temp = atan2(sinAngle, cosAngle);
	Real angleRightHanded, angleLeftHanded;

	// Left-handed angle will always be negative, right-handed angle will be positive.
	if (temp >= 0) {
		angleRightHanded = temp;
		angleLeftHanded  = temp - 2*Pi;
	}

	else {
		angleLeftHanded  = temp;
		angleRightHanded = temp + 2*Pi;
	}

	// Orientation of arc. A negative angle means a left-handed rotation around e3. 
	Real angle;
	if (M >= 0)	{
		angle = angleRightHanded;
	}

	else {
		angle = angleLeftHanded;
	}

	// Create the last unit vector to form the orthonormal basis to describe the arc.
	UnitVec3 e2(e3 % e1, true);

    setGeodesicToArc(e1, e2, radius, angle, geod);
}

void ContactGeometry::Sphere::Impl::shootGeodesicInDirectionUntilLengthReachedAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options, Geodesic& geod) const {

    UnitVec3 e_OP(xP);
    Real angle = terminatingLength/radius;

    setGeodesicToArc(e_OP, tP, radius, angle, geod);
}

void ContactGeometry::Sphere::Impl::shootGeodesicInDirectionUntilPlaneHitAnalytical(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {

    UnitVec3 e_OP(xP);

    // solve ~( e_OP * cos(t) + tP * sin(t) - pt_on_plane )*plane_normal = 0
    // for sphere plane offset is zero, therefore pt_on_plane = 0
    Real a = ~e_OP*terminatingPlane.getNormal();
    Real b = ~tP*terminatingPlane.getNormal();
    Real alpha = std::atan2(a,b);
    Real angle = (alpha > 0 ? Pi-alpha : -alpha);
//    std::cout << "a=" << a << ", b=" << b << ", alpha = " << alpha << std::endl;

    setGeodesicToArc(e_OP, tP, radius, angle, geod);
}


Real SphereImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/square(ownerp->getRadius());
}

Real SphereImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    if (derivComponents.size() == 1)
        return -2*x[derivComponents[0]]/square(ownerp->getRadius());
    if (derivComponents[0] == derivComponents[1])
        return -2/square(ownerp->getRadius());
    return 0;
}



//==============================================================================
//                               ELLIPSOID & IMPL
//==============================================================================

ContactGeometry::Ellipsoid::Ellipsoid(const Vec3& radii)
:   ContactGeometry(new Ellipsoid::Impl(radii)) {}

void ContactGeometry::Ellipsoid::setRadii(const Vec3& radii) 
{   updImpl().setRadii(radii); }

/*static*/ ContactGeometryTypeId ContactGeometry::Ellipsoid::classTypeId()
{   return ContactGeometry::Ellipsoid::Impl::classTypeId(); }

const Vec3& ContactGeometry::Ellipsoid::getRadii() const 
{   return getImpl().getRadii(); }

const Vec3& ContactGeometry::Ellipsoid::getCurvatures() const 
{   return getImpl().getCurvatures(); }

UnitVec3 ContactGeometry::Ellipsoid::
findUnitNormalAtPoint(const Vec3& Q) const
{   return getImpl().findUnitNormalAtPoint(Q); }

Vec3 ContactGeometry::Ellipsoid::
findPointWithThisUnitNormal(const UnitVec3& nn) const
{   return getImpl().findPointWithThisUnitNormal(nn); }

Vec3 ContactGeometry::Ellipsoid::
findPointInSameDirection(const Vec3& Q) const
{   return getImpl().findPointInSameDirection(Q); }

void ContactGeometry::Ellipsoid::
findParaboloidAtPoint(const Vec3& Q, Transform& X_EP, Vec2& k) const
{   return getImpl().findParaboloidAtPoint(Q,X_EP,k); }

void ContactGeometry::Ellipsoid::
findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& nn,
                                Transform& X_EP, Vec2& k) const
{   return getImpl().findParaboloidAtPointWithNormal(Q,nn,X_EP,k); }


const ContactGeometry::Ellipsoid::Impl& ContactGeometry::Ellipsoid::
getImpl() const {
    assert(impl);
    return static_cast<const Ellipsoid::Impl&>(*impl);
}

ContactGeometry::Ellipsoid::Impl& ContactGeometry::Ellipsoid::
updImpl() {
    assert(impl);
    return static_cast<Ellipsoid::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Ellipsoid::Impl::createDecorativeGeometry() const {
    return DecorativeEllipsoid(radii);
}

// Given a point Q on an ellipsoid, with outward unit normal nn at Q: find the 
// principal curvatures at the point and their directions. The result is a 
// coordinate frame with origin Q, z axis the ellipsoid normal nn at Q, x axis 
// is the direction dmax of maximum curvature kmax, y axis the direction dmin 
// of minimum curvature kmin, such that [dmax dmin n] forms a right-handed set.
// This is equivalent to fitting an elliptic paraboloid 
// z = -kmax/2 x^2 -kmin/2 y^2 to the ellipsoid at point Q. Note that for
// an ellipsoid we have kmax>=kmin>0.
//
// We'll find the ellipse on the central plane perpendicular to the normal by 
// intersecting the plane equation with the ellipsoid equation but working in 
// the plane frame P=[u v n], where u and v are arbitrary axes in the plane.
// Our goal is to obtain an equation for the ellipse in P and then rotate the 
// P frame about its normal until we get the ellipse in standard form 
// Ru^2+Sv^2=1 in which case d/R and d/S are the ellipsoid curvatures (d is the
// distance from the point on the ellipsoid to the plane).
// ref: McArthur, Neil. "Principal radii of curvature at a point on an 
// ellipsoid", Mathematical Notes 24 pp. xvi-xvii, 1929.
//
// In its own frame E=[x y z] the ellipsoid surface is the set of points such 
// that
//    ~e * diag(A,B,C) * e = 1
// where e is a vector expressed in E. The plane is the set of points 
// satisfying ~e * n = 0. We can write rotation matrix R_EP=[u v n] where 
// u,v,n are expressed in E. Now we can put the ellipsoid in P:
//   ~(R_EP*p) * diag(A,B,C) * (R_EP*p) = 1
// We can intersect that with the plane just by dropping the n coordinate of 
// p so p=[u v 0] (u,v scalars here), and the intersection equation is
//    A(u*ux + v*vx)^2 + B(u*uy+v*vy)^2 + C(u*uz + v*vz)^2 = 1
// which is
//    R u^2 + S v^2 + T u*v = 1
// with
//    R =   A ux^2  + B uy^2  + C uz^2
//    S =   A vx^2  + B vy^2  + C vz^2
//    T = 2(A ux*vx + B uy*vy + C uz*vz)
//
// We want to find a rotation about n that eliminates the cross term Tuv, 
// leaving us with
//    R' u'^2 + S' v'^2 = 1
// for new constants R' and S' and new basis u' and v'.
//
// Method
// ------
// We'll calculate an angle theta where theta=0 would be along u and 
// theta=pi/2 would be along v. Then theta+pi/2 is a perpendicular direction 
// that has the other curvature extreme. Per "Dr Rob" at Mathforum.org 2000:
//   t2t = tan(2*theta) = T/(R-S)
//   theta = atan(t2t)/2, c = cos(theta), s = sin(theta)
//   R' = Rc^2 + Tsc + Ss^2   (theta direction)
//   S' = Rs^2 - Tsc + Sc^2   (theta+pi/2 direction)
// Directions are u' = c*u + s*v, v' = c*v - s*u; these are automatically unit
// vectors.
//
// Optimization
// ------------
// The above requires an atan() to get 2*theta then sin & cos(theta) at
// a cost of about 120 flops. We can use half angle formulas to work
// exclusively with 2*theta, but then we'll have to normalize u' and v' 
// at the end:
//   t2t = tan(2*theta) = T/(R-S)
//   c2t = cos(2*theta) = 1/sqrt(1 + t2t^2)
//   s2t = sin(2*theta) = t2t*cos2t;
//   2*R' = R+S + Rc2t - Sc2t + Ts2t
//   2*S' = R+S - Rc2t + Sc2t - Ts2t
// By multiplying the u',v' formulas above by 2*c we change the lengths
// but get expressions that are easily converted to double angles:
//   u' = normalize((1+c2t)*u + s2t*v)
//   v' = normalize((1+c2t)*v - s2t*u)
// (but actually v' is n X u' which is cheap). This saves about 30 
// flops over the straightforward method above.
//
// Cost: given a point and normalized normal
//    curvatures ~160 flops
//    directions ~ 60 flops more
//               ----
//               ~220 flops
//
// So: Given an ellipsoid in its own frame E, with equation Ax^2+By^2+Cz^2=1, a 
// point Q=(x,y,z) on its surface, and the unit outward normal vector nn at Q,
// return (kmax,kmin) the principal curvatures at Q, and a Transform with 
// x=dmax, y=dmin, z=nn, O=Q that gives the principal curvature directions. 
// (Note: A=1/a^2, B=1/b^2, C=1/c^2 where a,b,c are the ellipsoid radii.)
void ContactGeometry::Ellipsoid::Impl::
findParaboloidAtPointWithNormal(const Vec3& Q, const UnitVec3& nn,
                                Transform& X_EP, Vec2& k) const
{
    const Real A = square(curvatures[0]), B = square(curvatures[1]), 
               C = square(curvatures[2]);

    // Sanity checks in debug.
    SimTK_ERRCHK(std::abs(A*Q[0]*Q[0]+B*Q[1]*Q[1]+C*Q[2]*Q[2]-1) < SqrtEps,
        "ContactGeometry::Ellipsoid::findParaboloidAtPointWithNormal()",
        "The given point was not on the surface of the ellipsoid.");
    SimTK_ERRCHK((nn-findUnitNormalAtPoint(Q)).normSqr() < SqrtEps,
        "ContactGeometry::Ellipsoid::findParaboloidAtPointWithNormal()",
        "The given normal was not consistent with the given point.");

    UnitVec3 tu = nn.perp();    // ~40 flops
    UnitVec3 tv(nn % tu, true); // y = z X x for plane, already normalized (9 flops)
    
    // 27 flops to get R,S,T
    Real R=   A*square(tu[0]) + B*square(tu[1]) + C*square(tu[2]);
    Real S=   A*square(tv[0]) + B*square(tv[1]) + C*square(tv[2]);
    Real T=2*(A*tu[0]*tv[0]   + B*tu[1]*tv[1]   + C*tu[2]*tv[2]);

    // T will be zero for spheres (A=B=C) and for various "clean" points
    // on the ellipsoid where tu[i]*tv[i]==0, i=0,1,2. In that case we
    // already have the ellipse we're looking for with R,S.
    // R==S means curvature is the same in every direction (that's called
    // an "umbilic" point). In that case tu and tv are good directions.
    // I *believe* R==S -> T==0 but I don't have a proof.
    Real kmax2, kmin2; // squared curvatures of ellipse
    UnitVec3 dmax;
    if (std::abs(R-S) < SignificantReal*std::max(R,S)) {
        kmax2 = kmin2 = (R+S)/2;
        dmax = tu;
    } else if (std::abs(T) < SignificantReal) {
        if (R < S) kmax2=S, dmax=tv, kmin2=R;
        else       kmax2=R, dmax=tu, kmin2=S;
    } else { // T,R-S both nonzero
        Real tan2t = T/(R-S);       // ~20 flops
        Real cos2t = 1/std::sqrt(1 + square(tan2t)); // ~40 flops
        Real sin2t = tan2t*cos2t;   //   1 flop
        // 11 flops here
        Real term = R*cos2t-S*cos2t+T*sin2t;
        Real Rp = (R+S + term)/2;
        Real Sp = (R+S - term)/2;

        // Sort into kmax, kmin; at most one normalization done below
        if (Rp < Sp) {
            kmax2=Sp, kmin2=Rp;
            dmax = UnitVec3((1+cos2t)*tv - sin2t*tu); // Sdir, must normalize, ~50 flops
        } else {
            kmax2=Rp,kmin2=Sp;
            dmax = UnitVec3((1+cos2t)*tu + sin2t*tv); // Rdir, must normalize, ~50 flops
        }
    }

    Real d = ~Q * nn; // distance along normal from center to point on ellipsoid (5 flops)
    Real kmax = d * kmax2, kmin = d * kmin2; // surface curvatures (2 flops)

    X_EP.updP() = Q; // the origin point
    Rotation& R_EP = X_EP.updR();
    // 9 flops
    UnitVec3 dmin = UnitVec3(nn % dmax, true); // y=z%x ensures right handedness (already unit vector too)
    R_EP.setRotationFromUnitVecsTrustMe(dmax, dmin, nn);

    k = Vec2(kmax, kmin);
}


// Peter E. says he implemented this from David Eberly's web site
// http://www.geometrictools.com/Documentation/DistancePointToEllipsoid.pdf
// Eberly says he got it from John Hart's article in Graphics Gems 4, page
// 113 "Distance to an Ellipsoid". Both Eberly and Hart recommend using a
// Newton iteration to solve this problem because the largest root is directly
// downhill given appropriate starting points, which they provide. However,
// the implementation here uses a direct solution of the 6th-order polynomial
// then searches for the largest real root. That is likely to be *much* slower
// than the recommended approach, although that should be measured.
//
// I asked Peter and he said he did not try and reject the Newton approach;
// he just took the direct approach. I believe the Newton method would be
// *much* faster, but Eberly hints that there are special cases that can
// cause convergence troubles and must be dealt with carefully. If the
// existing routine turns out to be a bottleneck, it would be worth revisiting
// this implementation. -- Sherm 20110203.
//
// TODO: use faster method?
Vec3 ContactGeometry::Ellipsoid::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    Real a2 = radii[0]*radii[0];
    Real b2 = radii[1]*radii[1];
    Real c2 = radii[2]*radii[2];
    Real a4 = a2*a2;
    Real b4 = b2*b2;
    Real c4 = c2*c2;
    Real px2 = position[0]*position[0];
    Real py2 = position[1]*position[1];
    Real pz2 = position[2]*position[2];
    Real a2b2 = a2*b2;
    Real b2c2 = b2*c2;
    Real a2c2 = a2*c2;
    Real a2b2c2 = a2b2*c2;
    Vector coeff(7);
    coeff[0] = 1;
    coeff[1] = 2*(a2+b2+c2);
    coeff[2] = -(a2*px2+b2*py2+c2*pz2) + a4+b4+c4 + 4*(a2b2+b2c2+a2c2);
    coeff[3] = -2*((a2b2+a2c2)*px2+(a2b2+b2c2)*py2+(b2c2+a2c2)*pz2) + 2*(a4*(b2+c2)+b4*(a2+c2)+c4*(a2+b2)) + 8*a2b2c2;
    coeff[4] = -a2*(b4+4*b2c2+c4)*px2-b2*(a4+4*a2c2+c4)*py2-c2*(a4+4*a2b2+b4)*pz2 + 4*(a2+b2+c2)*a2b2c2 + a4*b4+a4*c4+b4*c4;
    coeff[5] = 2*a2b2c2*(-(b2+c2)*px2-(a2+c2)*py2-(a2+b2)*pz2 + a2b2+b2c2+a2c2);
    coeff[6] = a2b2c2*(-b2c2*px2-a2c2*py2-a2b2*pz2+a2b2c2);
    Vector_<complex<Real> > roots(6);
    PolynomialRootFinder::findRoots(coeff, roots);
    Real root = NTraits<Real>::getMostNegative();
    for (int i = 0; i < 6; i++)
        if (fabs(roots[i].imag()) < 1e-10 && (roots[i].real()) > (root))
            root = roots[i].real();
    Vec3 result(position[0]*a2/(root+a2), position[1]*b2/(root+b2), position[2]*c2/(root+c2));
    Vec3 ri2(1/a2, 1/b2, 1/c2);
    inside = (position[0]*position[0]*ri2[0] + position[1]*position[1]*ri2[1] + position[2]*position[2]*ri2[2] < 1.0);
    normal = UnitVec3(result[0]*ri2[0], result[1]*ri2[1], result[2]*ri2[2]);
    return result;
}

// Peter says he took this algorithm from Art of Illusion but can't remember
// where it came from. It is similar to an algorithm presented in this thread:
// http://www.ogre3d.org/forums/viewtopic.php?f=2&t=26442&start=0
// and is most likely a special case of the general ray-quadric intersection
// method presented by Cychosz and Waggenspack in Graphics Gems III, pg. 275,
// "Intersecting a ray with a quadric surface."
bool ContactGeometry::Ellipsoid::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    Real rx2 = radii[0]*radii[0];
    Real sy = rx2/(radii[1]*radii[1]);
    Real sz = rx2/(radii[2]*radii[2]);
    Vec3 scaledDir(direction[0], sy*direction[1], sz*direction[2]);
    Real b = -(~scaledDir*origin);
    Real c = origin[0]*origin[0] + sy*origin[1]*origin[1] + sz*origin[2]*origin[2] - rx2;
    if (c > 0) {
        // Ray origin is outside ellipsoid.

        if (b <= 0)
          return false;  // Ray points away from the ellipsoid.
        Real a = ~scaledDir*direction;;
        Real d = b*b - a*c;
        if (d < 0)
          return false;
        distance = (b - std::sqrt(d))/a;
    }
    else {
        // Ray origin is inside ellipsoid.

        Real a = ~scaledDir*direction;;
        Real d = b*b - a*c;
        if (d < 0)
          return false;
        distance = (b + std::sqrt(d))/a;
    }
    Vec3 pos = origin+distance*direction;
    normal = UnitVec3(pos[0], pos[1]*sy, pos[2]*sz);
    return true;
}

void ContactGeometry::Ellipsoid::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = max(radii);
}

void ContactGeometry::Ellipsoid::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    Transform transform;
    findParaboloidAtPoint(point, transform, curvature);
    orientation = transform.R();
}


//TODO: just an axis-aligned leaf box for now
void ContactGeometry::Ellipsoid::Impl::createOBBTree() {
    OBBNode& root = obbTree.updRoot();
    root.box.setHalfLengths(radii);
    root.normal = UnitVec3(XAxis);  // doesn't matter
    root.coneHalfAngle = Pi;        // has all possible normals
    root.pointOnSurface = Vec3(radii[0],0,0); // doesn't matter
    root.children.clear(); // This is a leaf

    // Leaf contents.
    root.centerUW = Vec2(0,0);
    root.dims = Vec2(Pi, Pi/2); // u in [-Pi,Pi], v in [-Pi/2,Pi/2]
}

Real EllipsoidImplicitFunction::
calcValue(const Vector& x) const {
    const Vec3& radii = ownerp->getRadii();
    return 1-x[0]*x[0]/(radii[0]*radii[0])-x[1]*x[1]/(radii[1]*radii[1])-x[2]*x[2]/(radii[2]*radii[2]);
}

Real EllipsoidImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    const Vec3& radii = ownerp->getRadii();
    if (derivComponents.size() == 1) {
        int c = derivComponents[0];
        return -2*x[c]/(radii[c]*radii[c]);
    }
    if (derivComponents[0] == derivComponents[1]) {
        int c = derivComponents[0];
        return -2/(radii[c]*radii[c]);
    }
    return 0;
}



//==============================================================================
//                          SMOOTH HEIGHT MAP & IMPL
//==============================================================================

ContactGeometry::SmoothHeightMap::
SmoothHeightMap(const BicubicSurface& surface) 
:   ContactGeometry(new SmoothHeightMap::Impl(surface)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::SmoothHeightMap::
classTypeId() 
{   return ContactGeometry::SmoothHeightMap::Impl::classTypeId(); }

const BicubicSurface& ContactGeometry::SmoothHeightMap::
getBicubicSurface() const {return getImpl().getBicubicSurface();}

const OBBTree& ContactGeometry::SmoothHeightMap::
getOBBTree() const {return getImpl().getOBBTree();}

const ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
getImpl() const {
    assert(impl);
    return static_cast<const SmoothHeightMap::Impl&>(*impl);
}

ContactGeometry::SmoothHeightMap::Impl& ContactGeometry::SmoothHeightMap::
updImpl() {
    assert(impl);
    return static_cast<SmoothHeightMap::Impl&>(*impl);
}

// This is the main constructor.
ContactGeometry::SmoothHeightMap::Impl::
Impl(const BicubicSurface& surface) 
:   surface(surface) { 
    implicitFunction.setOwner(*this); 

    createBoundingVolumes();

}

void ContactGeometry::SmoothHeightMap::Impl::
assignPatch(const Geo::BicubicBezierPatch& patch, 
            OBBNode& node, int depth,
            Array_<const Vec3*>* parentControlPoints) const 
{
    const Mat<4,4,Vec3>& nodeB = patch.getControlPoints();
    const Vec2& nodeB11 = nodeB(0,0).getSubVec<2>(0); // just x,y
    const Vec2& nodeB44 = nodeB(3,3).getSubVec<2>(0);
    node.centerUW = (nodeB11+nodeB44)/2;
    node.dims     = (nodeB11-nodeB44).abs()/2;

    // For now just split 4 ways; need to be done recursively based on
    // flatness of patch.
    node.children.resize(4);
    patch.split(0.5,0.5,node.children[0].patch, node.children[1].patch,
                        node.children[2].patch, node.children[3].patch);
    Array_<const Vec3*> myControlPoints;
    for (int c=0; c<4; ++c) {
        OBBNode& child = node.children[c];
        child.depth = depth+1;
        child.height = 0;
        child.box = child.patch.calcOrientedBoundingBox();
        const Mat<4,4,Vec3>& B = child.patch.getControlPoints();
        const Vec2& b11 = B(0,0).getSubVec<2>(0); // just x,y
        const Vec2& b44 = B(3,3).getSubVec<2>(0);
        child.centerUW = (b11+b44)/2;
        child.dims     = (b11-b44).abs()/2;
        for (int i=0; i<4; ++i) 
            for (int j=0; j<4; ++j) 
                myControlPoints.push_back(&B(i,j));
    }

    node.depth = depth;
    node.height = 1 + std::max(node.children[0].height,
                      std::max(node.children[1].height,
                      std::max(node.children[2].height,
                               node.children[3].height)));
    node.box = Geo::Point::calcOrientedBoundingBoxIndirect(myControlPoints);

    if (parentControlPoints)
        for (unsigned i=0; i<myControlPoints.size(); ++i)
            parentControlPoints->push_back(myControlPoints[i]);
}

void ContactGeometry::SmoothHeightMap::Impl::
splitPatches(int x0,int y0, int nx, int ny, 
             OBBNode& node, int depth,
             Array_<const Vec3*>* parentControlPoints) const {
    assert(nx>0 && ny>0 && depth>=0);


    node.x0=0; node.y0=0; node.nx=nx; node.ny=ny;
    if (nx==1 && ny==1) {
        assignPatch(surface.calcBezierPatch(x0,y0), node, depth,
                    parentControlPoints);
        return;
    } 

    // Add two children.
    node.children.resize(2);
    Array_<const Vec3*> myControlPoints;

    // Split on the long direction
    if (nx > ny) {
        splitPatches(x0,      y0, nx/2,    ny, node.children[0],
            depth+1, &myControlPoints);
        splitPatches(x0+nx/2, y0, nx-nx/2, ny, node.children[1],
            depth+1, &myControlPoints);
    } else {
        splitPatches(x0, y0,      nx, ny/2,    node.children[0],
            depth+1, &myControlPoints);
        splitPatches(x0, y0+ny/2, nx, ny-ny/2, node.children[1],
            depth+1, &myControlPoints);
    }
    node.depth = depth;
    node.height = 1 + std::max(node.children[0].height,
                               node.children[1].height);

    node.box = Geo::Point::calcOrientedBoundingBoxIndirect(myControlPoints);

    if (parentControlPoints) {
        for (unsigned i=0; i<myControlPoints.size(); ++i)
            parentControlPoints->push_back(myControlPoints[i]);
    }
}

void ContactGeometry::SmoothHeightMap::Impl::
createBoundingVolumes() {
    // Temporarily convert the surface into a set of Bezier patches (using
    // a lot more memory than the original).

    int nx,ny; surface.getNumPatches(nx,ny);
    OBBNode& root = obbTree.updRoot();
    splitPatches(0,0,nx,ny,root,0);


    // Create bounding sphere.
    // TODO: fake this using mesh; this needs to be done correctly instead
    // by the BicubicSurface itself. Using 5 subdivisions per patch.
    PolygonalMesh mesh = surface.createPolygonalMesh(5);

    // Collect all the vertices.
    const int n = mesh.getNumVertices();
    Array_<const Vec3*> points(n);
    for (int i=0; i<n; ++i)
        points[i] = &mesh.getVertexPosition(i);
    boundingSphere = Geo::Point::calcBoundingSphereIndirect(points);
    // Add 10% as a hack to make it less likely we'll miss part of the surface.
    boundingSphere.updRadius() *= 1.1;
}

DecorativeGeometry ContactGeometry::SmoothHeightMap::Impl::createDecorativeGeometry() const {
    return DecorativeMesh(surface.createPolygonalMesh());
}

Vec3 ContactGeometry::SmoothHeightMap::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    assert(false);
    return Vec3(NaN);
}

bool ContactGeometry::SmoothHeightMap::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, 
    Real& distance, UnitVec3& normal) const 
{
    assert(false);
    return true;
}

Real SmoothHeightMapImplicitFunction::
calcValue(const Vector& p) const {
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real z = surf.calcValue(Vec2(p[0],p[1]), hint);
    //TODO: this is negated from convention
    return z - p[2]; // negative outside, positive inside
}

// First deriv with respect to p[2] (z component) is -1 to match the above
// implicit function definition, all higher derivs are
// with respect to that component are 0.
Real SmoothHeightMapImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& p) const {
    if (derivComponents.empty()) return calcValue(p);
    if (derivComponents.size() == 1 && derivComponents[0]==2)
        return -1;
    for (unsigned i=0; i<derivComponents.size(); ++i)
        if (derivComponents[i]==2) return 0;

    // We're asking only for derivatives in x and y.
    const BicubicSurface&      surf = ownerp->getBicubicSurface();
    BicubicSurface::PatchHint& hint = ownerp->updHint();
    const Real d = surf.calcDerivative(derivComponents, Vec2(p[0],p[1]), hint);
    return d;
}




//==============================================================================
//                              TRIANGLE MESH
//==============================================================================

ContactGeometry::TriangleMesh::TriangleMesh
   (const ArrayViewConst_<Vec3>& vertices, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(vertices, faceIndices, smooth)) {}

ContactGeometry::TriangleMesh::TriangleMesh
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometry(new TriangleMesh::Impl(mesh, smooth)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::TriangleMesh::classTypeId() 
{   return ContactGeometry::TriangleMesh::Impl::classTypeId(); }


int ContactGeometry::TriangleMesh::getNumEdges() const {
    return getImpl().edges.size();
}

int ContactGeometry::TriangleMesh::getNumFaces() const {
    return getImpl().faces.size();
}

int ContactGeometry::TriangleMesh::getNumVertices() const {
    return getImpl().vertices.size();
}

const Vec3& ContactGeometry::TriangleMesh::getVertexPosition(int index) const {
    assert(index >= 0 && index < getNumVertices());
    return getImpl().vertices[index].pos;
}

int ContactGeometry::TriangleMesh::getFaceEdge(int face, int edge) const {
    assert(face >= 0 && face < getNumFaces());
    assert(edge >= 0 && edge < 3);
    return getImpl().faces[face].edges[edge];
}

int ContactGeometry::TriangleMesh::getFaceVertex(int face, int vertex) const {
    assert(face >= 0 && face < getNumFaces());
    assert(vertex >= 0 && vertex < 3);
    return getImpl().faces[face].vertices[vertex];
}

int ContactGeometry::TriangleMesh::getEdgeFace(int edge, int face) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(face >= 0 && face < 2);
    return getImpl().edges[edge].faces[face];
}

int ContactGeometry::TriangleMesh::getEdgeVertex(int edge, int vertex) const {
    assert(edge >= 0 && edge < getNumEdges());
    assert(vertex >= 0 && vertex < 2);
    return getImpl().edges[edge].vertices[vertex];
}

const UnitVec3& ContactGeometry::TriangleMesh::getFaceNormal(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].normal;
}

Real ContactGeometry::TriangleMesh::getFaceArea(int face) const {
    assert(face >= 0 && face < getNumFaces());
    return getImpl().faces[face].area;
}

void ContactGeometry::TriangleMesh::
findVertexEdges(int vertex, Array_<int>& edges) const {
    // Begin at an arbitrary edge which intersects the vertex.
    
    int firstEdge = getImpl().vertices[vertex].firstEdge;
    int previousEdge = firstEdge;
    int previousFace = getImpl().edges[firstEdge].faces[0];
    
    // Walk around the vertex, using each edge to find the next face and each 
    // face to find the next edge.
    
    do {
        edges.push_back(previousEdge);
        const ContactGeometry::TriangleMesh::Impl::Edge& 
            edge = getImpl().edges[previousEdge];
        int nextFace = (edge.faces[0] == previousFace ? edge.faces[1] 
                                                      : edge.faces[0]);
        const ContactGeometry::TriangleMesh::Impl::Face& 
            face = getImpl().faces[nextFace];
        int nextEdge;
        if (    face.edges[0] != previousEdge
            && (face.vertices[0] == vertex || face.vertices[1] == vertex))
            nextEdge = face.edges[0];
        else if (   face.edges[1] != previousEdge 
                 && (face.vertices[1] == vertex || face.vertices[2] == vertex))
            nextEdge = face.edges[1];
        else
            nextEdge = face.edges[2];
        previousEdge = nextEdge;
        previousFace = nextFace;
    } while (previousEdge != firstEdge);
}

Vec3 ContactGeometry::TriangleMesh::findPoint(int face, const Vec2& uv) const {
    return getImpl().findPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findCentroid(int face) const {
    return getImpl().findCentroid(face);
}

UnitVec3 ContactGeometry::TriangleMesh::
findNormalAtPoint(int face, const Vec2& uv) const {
    return getImpl().findNormalAtPoint(face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPoint
   (const Vec3& position, bool& inside, int& face, Vec2& uv) const {
    return getImpl().findNearestPoint(position, inside, face, uv);
}

Vec3 ContactGeometry::TriangleMesh::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    return getImpl().findNearestPointToFace(position, face, uv);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, 
    UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

bool ContactGeometry::TriangleMesh::intersectsRay
   (const Vec3& origin, const UnitVec3& direction, Real& distance, int& face, 
    Vec2& uv) const {
    return getImpl().intersectsRay(origin, direction, distance, face, uv);
}

ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::getOBBTreeNode() const {
    return OBBTreeNode(getImpl().obb);
}

PolygonalMesh ContactGeometry::TriangleMesh::createPolygonalMesh() const {
    PolygonalMesh mesh;
    getImpl().createPolygonalMesh(mesh);
    return mesh;
}

const ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::getImpl() const {
    assert(impl);
    return static_cast<const TriangleMesh::Impl&>(*impl);
}

ContactGeometry::TriangleMesh::Impl& 
ContactGeometry::TriangleMesh::updImpl() {
    assert(impl);
    return static_cast<TriangleMesh::Impl&>(*impl);
}



//==============================================================================
//                            TRIANGLE MESH IMPL
//==============================================================================

DecorativeGeometry ContactGeometry::TriangleMesh::Impl::createDecorativeGeometry() const {
    PolygonalMesh mesh;
    createPolygonalMesh(mesh);
    return DecorativeMesh(mesh);
}

Vec3 ContactGeometry::TriangleMesh::Impl::findPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    return             uv[0] * vertices[f.vertices[0]].pos
           +           uv[1] * vertices[f.vertices[1]].pos
           +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].pos;
}

// same as findPoint(face, (1/3,1/3)) but faster
Vec3 ContactGeometry::TriangleMesh::Impl::findCentroid(int face) const {
    const Face& f = faces[face];
    return (  vertices[f.vertices[0]].pos
            + vertices[f.vertices[1]].pos
            + vertices[f.vertices[2]].pos) / 3;
}

UnitVec3 ContactGeometry::TriangleMesh::Impl::findNormalAtPoint
   (int face, const Vec2& uv) const {
    const Face& f = faces[face];
    if (smooth)
        return UnitVec3(            uv[0] * vertices[f.vertices[0]].normal
                        +           uv[1] * vertices[f.vertices[1]].normal
                        +  (1-uv[0]-uv[1])* vertices[f.vertices[2]].normal);
    return f.normal;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    int face;
    Vec2 uv;
    Vec3 nearestPoint = findNearestPoint(position, inside, face, uv);
    normal = findNormalAtPoint(face, uv);
    return nearestPoint;
}

Vec3 ContactGeometry::TriangleMesh::Impl::
findNearestPoint(const Vec3& position, bool& inside, int& face, Vec2& uv) const 
{
    Real distance2;
    Vec3 nearestPoint = obb.findNearestPoint(*this, position, MostPositiveReal, distance2, face, uv);
    Vec3 delta = position-nearestPoint;
    inside = (~delta*faces[face].normal < 0);
    return nearestPoint;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              UnitVec3& normal) const {
    int face;
    Vec2 uv;
    if (!intersectsRay(origin, direction, distance, face, uv))
        return false;
    normal = findNormalAtPoint(face, uv);
    return true;
}

bool ContactGeometry::TriangleMesh::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    Real boundsDistance;
    if (!obb.bounds.intersectsRay(origin, direction, boundsDistance))
        return false;
    return obb.intersectsRay(*this, origin, direction, distance, face, uv);
}

void ContactGeometry::TriangleMesh::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = boundingSphereCenter;
    radius = boundingSphereRadius;
}

void ContactGeometry::TriangleMesh::Impl::
createPolygonalMesh(PolygonalMesh& mesh) const {
    for (unsigned vx=0; vx < vertices.size(); ++vx)
        mesh.addVertex(vertices[vx].pos);
    for (unsigned fx=0; fx < faces.size(); ++fx) {
        const Face& face = faces[fx];
        const ArrayViewConst_<int> verts(face.vertices, face.vertices+3);
        mesh.addFace(verts);
    }
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const ArrayViewConst_<Vec3>& vertexPositions, 
    const ArrayViewConst_<int>& faceIndices, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) {
    init(vertexPositions, faceIndices);
}

ContactGeometry::TriangleMesh::Impl::Impl
   (const PolygonalMesh& mesh, bool smooth) 
:   ContactGeometryImpl(), smooth(smooth) 
{   // Create the mesh, triangulating faces as necessary.
    Array_<Vec3>    vertexPositions;
    Array_<int>     faceIndices;
    for (int i = 0; i < mesh.getNumVertices(); i++)
        vertexPositions.push_back(mesh.getVertexPosition(i));
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        int numVert = mesh.getNumVerticesForFace(i);
        if (numVert < 3)
            continue; // Ignore it.
        if (numVert == 3) {
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
        }
        else if (numVert == 4) {
            // Split it into two triangles.
            
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
            faceIndices.push_back(mesh.getFaceVertex(i, 1));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 2));
            faceIndices.push_back(mesh.getFaceVertex(i, 3));
            faceIndices.push_back(mesh.getFaceVertex(i, 0));
        }
        else {
            // Add a vertex at the center, then split it into triangles.
            
            Vec3 center(0);
            for (int j = 0; j < numVert; j++)
                center += vertexPositions[mesh.getFaceVertex(i, j)];
            center /= numVert;
            vertexPositions.push_back(center);
            int newIndex = vertexPositions.size()-1;
            for (int j = 0; j < numVert-1; j++) {
                faceIndices.push_back(mesh.getFaceVertex(i, j));
                faceIndices.push_back(mesh.getFaceVertex(i, j+1));
                faceIndices.push_back(newIndex);
            }
        }
    }
    init(vertexPositions, faceIndices);
    
    // Make sure the mesh normals are oriented correctly.
    
    Vec3 origin(0);
    for (int i = 0; i < 3; i++)
        origin += vertices[faces[0].vertices[i]].pos;
    origin /= 3; // this is the face centroid

    const UnitVec3 direction = -faces[0].normal;
    // Calculate a ray origin that is guaranteed to be outside the
    // mesh. If the topology is right (face 0 normal points outward), we'll be
    // outside on the side containing face 0. If it is wrong, we'll be outside
    // on the opposite side of the mesh. Then we'll shoot a ray back along the
    // direction we came from (that is, towards the interior of the mesh from
    // outside). We'll hit *some* face. If the topology is right, the hit 
    // face's normal will be pointing back at us. If it is wrong, the face 
    // normal will also be pointing inwards, in roughly the same direction as 
    // the ray.
    origin -= max(obb.bounds.getSize())*direction;
    Real distance;
    int face;
    Vec2 uv;
    bool intersects = intersectsRay(origin, direction, distance, face, uv);
    assert(intersects);
    // Now dot the hit face normal with the ray direction; correct topology
    // will have them pointing in more-or-less opposite directions.
    if (dot(faces[face].normal, direction) > 0) {
        // We need to invert the mesh topology.
        
        for (int i = 0; i < (int) faces.size(); i++) {
            Face& f = faces[i];
            int temp = f.vertices[0];
            f.vertices[0] = f.vertices[1];
            f.vertices[1] = temp;
            temp = f.edges[1];
            f.edges[1] = f.edges[2];
            f.edges[2] = temp;
            f.normal *= -1;
        }
        for (int i = 0; i < (int) vertices.size(); i++)
            vertices[i].normal *= -1;
    }
}

void ContactGeometry::TriangleMesh::Impl::init
   (const Array_<Vec3>& vertexPositions, const Array_<int>& faceIndices) 
{   SimTK_APIARGCHECK_ALWAYS(faceIndices.size()%3 == 0, 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl", 
        "The number of indices must be a multiple of 3.");
    int numFaces = faceIndices.size()/3;
    
    // Create the vertices.
    
    for (int i = 0; i < (int) vertexPositions.size(); i++)
        vertices.push_back(Vertex(vertexPositions[i]));
    
    // Create the faces and build lists of all the edges.
    
    map<pair<int, int>, int> forwardEdges;
    map<pair<int, int>, int> backwardEdges;
    for (int i = 0; i < numFaces; i++) {
        int start = i*3;
        int v1 = faceIndices[start], v2 = faceIndices[start+1], 
            v3 = faceIndices[start+2];
        SimTK_APIARGCHECK1_ALWAYS
           (   v1 >= 0 && v1 < (int) vertices.size() 
            && v2 >= 0 && v2 < (int) vertices.size() 
            && v3 >= 0 && v3 < (int) vertices.size(),
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d contains a vertex with an illegal index.", i);
        Vec3 cross =   (vertexPositions[v2]-vertexPositions[v1])
                     % (vertexPositions[v3]-vertexPositions[v1]);
        Real norm = cross.norm();
        cross *= 1.0/norm;
        SimTK_APIARGCHECK1_ALWAYS(norm > 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Face %d is degenerate.", i);
        faces.push_back(Face(v1, v2, v3, cross, 0.5*norm));
        int edges[3][2] = {{v1, v2}, {v2, v3}, {v3, v1}};
        for (int j = 0; j < 3; j++) {
            SimTK_APIARGCHECK1_ALWAYS(edges[j][0] != edges[j][1], 
                "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                "Vertices %d appears twice in a single face.", edges[j][0]);
            if (edges[j][0] < edges[j][1]) {
                SimTK_APIARGCHECK2_ALWAYS
                   (forwardEdges.find(pair<int, int>(edges[j][0], edges[j][1])) 
                    == forwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][0], edges[j][1]);
                forwardEdges[pair<int, int>(edges[j][0], edges[j][1])] = i;
            }
            else {
                SimTK_APIARGCHECK2_ALWAYS
                   (backwardEdges.find(pair<int, int>(edges[j][1], edges[j][0]))
                    == backwardEdges.end(),
                    "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
                    "Multiple faces have an edge between vertices %d and %d"
                    " in the same order.", edges[j][1], edges[j][0]);
                backwardEdges[pair<int, int>(edges[j][1], edges[j][0])] = i;
            }
        }
    }
    
    // Create the edges.
    
    SimTK_APIARGCHECK_ALWAYS(forwardEdges.size() == backwardEdges.size(), 
        "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
        "Each edge must be shared by exactly two faces.");
    for (map<pair<int, int>, int>::iterator iter = forwardEdges.begin(); 
         iter != forwardEdges.end(); ++iter) {
        int vert1 = iter->first.first;
        int vert2 = iter->first.second;
        int face1 = iter->second;
        map<pair<int, int>, int>::iterator iter2 = 
            backwardEdges.find(pair<int, int>(vert1, vert2));
        SimTK_APIARGCHECK_ALWAYS(iter2 != backwardEdges.end(), 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Each edge must be shared by exactly two faces.");
        int face2 = iter2->second;
        edges.push_back(Edge(vert1, vert2, face1, face2));
    }
    
    // Record the edges for each face.
    
    for (int i = 0; i < (int) edges.size(); i++) {
        Edge& edge = edges[i];
        int f[2] = {edge.faces[0], edge.faces[1]};
        for (int j = 0; j < 2; j++) {
            Face& face = faces[f[j]];
            if ((edge.vertices[0] == face.vertices[0] || edge.vertices[0] == face.vertices[1]) &&
                    (edge.vertices[1] == face.vertices[0] || edge.vertices[1] == face.vertices[1]))
                face.edges[0] = i;
            else if ((edge.vertices[0] == face.vertices[1] || edge.vertices[0] == face.vertices[2]) &&
                    (edge.vertices[1] == face.vertices[1] || edge.vertices[1] == face.vertices[2]))
                face.edges[1] = i;
            else if ((edge.vertices[0] == face.vertices[2] || edge.vertices[0] == face.vertices[0]) &&
                    (edge.vertices[1] == face.vertices[2] || edge.vertices[1] == face.vertices[0]))
                face.edges[2] = i;
            else
                SimTK_ASSERT_ALWAYS(false, 
                    "Face and edge vertices are inconsistent.");
        }
    }
    
    // Record a single edge for each vertex.
    
    for (int i = 0; i < (int) edges.size(); i++) {
        vertices[edges[i].vertices[0]].firstEdge = i;
        vertices[edges[i].vertices[1]].firstEdge = i;
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        SimTK_APIARGCHECK1_ALWAYS(vertices[i].firstEdge >= 0, 
            "ContactGeometry::TriangleMesh::Impl", "TriangleMesh::Impl",
            "Vertex %d is not part of any face.", i);
    
    // Calculate a normal for each vertex.
    
    Vector_<Vec3> vertNorm(vertices.size(), Vec3(0));
    for (int i = 0; i < (int) faces.size(); i++) {
        const Face& f = faces[i];
        UnitVec3 edgeDir[3];
        for (int j = 0; j < 3; j++) {
            edgeDir[j] = UnitVec3(  vertices[f.vertices[(j+1)%3]].pos
                                  - vertices[f.vertices[j]].pos);
        }
        for (int j = 0; j < 3; j++) {
            Real angle = std::acos(~edgeDir[j]*edgeDir[(j+2)%3]);
            vertNorm[f.vertices[j]] += f.normal*angle;
        }
    }
    for (int i = 0; i < (int) vertices.size(); i++)
        vertices[i].normal = UnitVec3(vertNorm[i]);
    
    // Create the OBBTree.
    
    Array_<int> allFaces(faces.size());
    for (int i = 0; i < (int) allFaces.size(); i++)
        allFaces[i] = i;
    createObbTree(obb, allFaces);
    
    // Find the bounding sphere.
    Array_<const Vec3*> points(vertices.size());
    for (int i = 0; i < (int) vertices.size(); i++)
        points[i] = &vertices[i].pos;
    const Geo::Sphere bnd = Geo::Point::calcBoundingSphereIndirect(points);
    boundingSphereCenter = bnd.getCenter();
    boundingSphereRadius = bnd.getRadius();
}

void ContactGeometry::TriangleMesh::Impl::createObbTree
   (OBBTreeNodeImpl& node, const Array_<int>& faceIndices) 
{   // Find all vertices in the node and build the OrientedBoundingBox.
    node.numTriangles = faceIndices.size();
    set<int> vertexIndices;
    for (int i = 0; i < (int) faceIndices.size(); i++) 
        for (int j = 0; j < 3; j++)
            vertexIndices.insert(faces[faceIndices[i]].vertices[j]);
    Vector_<Vec3> points((int)vertexIndices.size());
    int index = 0;
    for (set<int>::iterator iter = vertexIndices.begin(); 
                            iter != vertexIndices.end(); ++iter)
        points[index++] = vertices[*iter].pos;
    node.bounds = OrientedBoundingBox(points);
    if (faceIndices.size() > 3) {

        // Order the axes by size.

        int axisOrder[3];
        const Vec3& size = node.bounds.getSize();
        if (size[0] > size[1]) {
            if (size[0] > size[2]) {
                axisOrder[0] = 0;
                if (size[1] > size[2]) {
                    axisOrder[1] = 1;
                    axisOrder[2] = 2;
                }
                else {
                    axisOrder[1] = 2;
                    axisOrder[2] = 1;
                }
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 0;
                axisOrder[2] = 1;
            }
        }
        else if (size[0] > size[2]) {
            axisOrder[0] = 1;
            axisOrder[1] = 0;
            axisOrder[2] = 2;
        }
        else {
            if (size[1] > size[2]) {
                axisOrder[0] = 1;
                axisOrder[1] = 2;
            }
            else {
                axisOrder[0] = 2;
                axisOrder[1] = 1;
            }
            axisOrder[2] = 0;
        }

        // Try splitting along each axis.

        for (int i = 0; i < 3; i++) {
            Array_<int> child1Indices, child2Indices;
            splitObbAxis(faceIndices, child1Indices, child2Indices, 
                         axisOrder[i]);
            if (child1Indices.size() > 0 && child2Indices.size() > 0) {
                // It was successfully split, so create the child nodes.

                node.child1 = new OBBTreeNodeImpl();
                node.child2 = new OBBTreeNodeImpl();
                createObbTree(*node.child1, child1Indices);
                createObbTree(*node.child2, child2Indices);
                return;
            }
        }
    }
    
    // This is a leaf node.
    
    node.triangles.insert(node.triangles.begin(), faceIndices.begin(), 
                          faceIndices.end());
}

void ContactGeometry::TriangleMesh::Impl::splitObbAxis
   (const Array_<int>& parentIndices, Array_<int>& child1Indices, 
    Array_<int>& child2Indices, int axis) 
{   // For each face, find its minimum and maximum extent along the axis.
    Vector minExtent(parentIndices.size());
    Vector maxExtent(parentIndices.size());
    for (int i = 0; i < (int) parentIndices.size(); i++) {
        int* vertexIndices = faces[parentIndices[i]].vertices;
        Real minVal = vertices[vertexIndices[0]].pos[axis];
        Real maxVal = vertices[vertexIndices[0]].pos[axis];
        minVal = std::min(minVal, vertices[vertexIndices[1]].pos[axis]);
        maxVal = std::max(maxVal, vertices[vertexIndices[1]].pos[axis]);
        minExtent[i] = std::min(minVal, vertices[vertexIndices[2]].pos[axis]);
        maxExtent[i] = std::max(maxVal, vertices[vertexIndices[2]].pos[axis]);
    }
    
    // Select a split point that tries to put as many faces as possible 
    // entirely on one side or the other.
    
    Real split = 0.5*(median(minExtent)+median(maxExtent));
    
    // Choose a side for each face.
    
    for (int i = 0; i < (int) parentIndices.size(); i++) {
        if (maxExtent[i] <= split)
            child1Indices.push_back(parentIndices[i]);
        else if (minExtent[i] >= split)
            child2Indices.push_back(parentIndices[i]);
        else if (0.5*(minExtent[i]+maxExtent[i]) <= split)
            child1Indices.push_back(parentIndices[i]);
        else
            child2Indices.push_back(parentIndices[i]);
    }
}

Vec3 ContactGeometry::TriangleMesh::Impl::findNearestPointToFace
   (const Vec3& position, int face, Vec2& uv) const {
    // Calculate the distance between a point in space and a face of the mesh.
    // This algorithm is based on a description by David Eberly found at 
    // http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf.
    
    const ContactGeometry::TriangleMesh::Impl::Face& fc = faces[face];
    const Vec3& vert1 = vertices[fc.vertices[0]].pos;
    const Vec3& vert2 = vertices[fc.vertices[1]].pos;
    const Vec3& vert3 = vertices[fc.vertices[2]].pos;
    const Vec3 e0 = vert2-vert1;
    const Vec3 e1 = vert3-vert1;
    const Vec3 delta = vert1-position;
    const Real a = e0.normSqr();
    const Real b = ~e0*e1;
    const Real c = e1.normSqr();
    const Real d = ~e0*delta;
    const Real e = ~e1*delta;
    const Real f = delta.normSqr();
    const Real det = a*c-b*b;
    Real s = b*e-c*d;
    Real t = b*d-a*e;
    if (s+t <= det) {
        if (s < 0) {
            if (t < 0) {
                // Region 4

                if (d < 0) {
                    s = (-d >= a ? 1 : -d/a);
                    t = 0;
                }
                else {
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            }
            else {
                // Region 3

                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 5

            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
            t = 0;
        }
        else {
            // Region 0

            const Real invDet = 1.0/det;
            s *= invDet;
            t *= invDet;
        }
    }
    else {
        if (s < 0) {
            // Region 2

            Real temp0 = b+d;
            Real temp1 = c+e;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1-s;
            }
            else {
                s = 0;
                t = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
            }
        }
        else if (t < 0) {
            // Region 6

            Real temp0 = b+e;
            Real temp1 = a+d;
            if (temp1 > temp0) {
                Real numer = temp1-temp0;
                Real denom = a-2*b+c;
                t = (numer >= denom ? 1 : numer/denom);
                s = 1-t;
            }
            else {
                s = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -d/a));
                t = 0;
            }
        }
        else {
            // Region 1

            const Real numer = c+e-b-d;
            if (numer <= 0)
                s = 0;
            else {
                const Real denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
            }
            t = 1-s;
        }
    }
    uv = Vec2(1-s-t, s);
    return vert1 + s*e0 + t*e1;
}


//==============================================================================
//                            OBB TREE NODE IMPL
//==============================================================================

OBBTreeNodeImpl::OBBTreeNodeImpl(const OBBTreeNodeImpl& copy) 
:   bounds(copy.bounds), triangles(copy.triangles), 
    numTriangles(copy.numTriangles) {
    if (copy.child1 == NULL) {
        child1 = NULL;
        child2 = NULL;
    }
    else {
        child1 = new OBBTreeNodeImpl(*copy.child1);
        child2 = new OBBTreeNodeImpl(*copy.child2);
    }
}

OBBTreeNodeImpl::~OBBTreeNodeImpl() {
    if (child1 != NULL)
        delete child1;
    if (child2 != NULL)
        delete child2;
}

Vec3 OBBTreeNodeImpl::findNearestPoint
   (const ContactGeometry::TriangleMesh::Impl& mesh, 
    const Vec3& position, Real cutoff2, 
    Real& distance2, int& face, Vec2& uv) const 
{
    Real tol = 100*Eps;
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance2 = MostPositiveReal, 
             child2distance2 = MostPositiveReal;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        Vec3 child1point, child2point;
        Real child1BoundsDist2 = 
            (child1->bounds.findNearestPoint(position)-position).normSqr();
        Real child2BoundsDist2 = 
            (child2->bounds.findNearestPoint(position)-position).normSqr();
        if (child1BoundsDist2 < child2BoundsDist2) {
            if (child1BoundsDist2 < cutoff2) {
                child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
                if (child2BoundsDist2 < child1distance2 && child2BoundsDist2 < cutoff2)
                    child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
            }
        }
        else {
            if (child2BoundsDist2 < cutoff2) {
                child2point = child2->findNearestPoint(mesh, position, cutoff2, child2distance2, child2face, child2uv);
                if (child1BoundsDist2 < child2distance2 && child1BoundsDist2 < cutoff2)
                    child1point = child1->findNearestPoint(mesh, position, cutoff2, child1distance2, child1face, child1uv);
            }
        }
        if (   child1distance2 <= child2distance2*(1+tol) 
            && child2distance2 <= child1distance2*(1+tol)) {
            // Decide based on angle which one to use.
            
            if (  std::abs(~(child1point-position)*mesh.faces[child1face].normal) 
                > std::abs(~(child2point-position)*mesh.faces[child2face].normal))
                child2distance2 = MostPositiveReal;
            else
                child1distance2 = MostPositiveReal;
        }
        if (child1distance2 < child2distance2) {
            distance2 = child1distance2;
            face = child1face;
            uv = child1uv;
            return child1point;
        }
        else {
            distance2 = child2distance2;
            face = child2face;
            uv = child2uv;
            return child2point;
        }
    }    
    // This is a leaf node, so check each triangle for its distance to the point.
    
    distance2 = MostPositiveReal;
    Vec3 nearestPoint;
    for (int i = 0; i < (int) triangles.size(); i++) {
        Vec2 triangleUV;
        Vec3 p = mesh.findNearestPointToFace(position, triangles[i], triangleUV);
        Vec3 offset = p-position;
        // TODO: volatile to work around compiler bug
        volatile Real d2 = offset.normSqr(); 
        if (d2 < distance2 || (d2 < distance2*(1+tol) && std::abs(~offset*mesh.faces[triangles[i]].normal) > std::abs(~offset*mesh.faces[face].normal))) {
            nearestPoint = p;
            distance2 = d2;
            face = triangles[i];
            uv = triangleUV;
        }
    }
    return nearestPoint;
}

bool OBBTreeNodeImpl::
intersectsRay(const ContactGeometry::TriangleMesh::Impl& mesh,
              const Vec3& origin, const UnitVec3& direction, Real& distance, 
              int& face, Vec2& uv) const {
    if (child1 != NULL) {
        // Recursively check the child nodes.
        
        Real child1distance, child2distance;
        int child1face, child2face;
        Vec2 child1uv, child2uv;
        bool child1intersects = child1->bounds.intersectsRay(origin, direction, child1distance);
        bool child2intersects = child2->bounds.intersectsRay(origin, direction, child2distance);
        if (child1intersects) {
            if (child2intersects) {
                // The ray intersects both child nodes.  First check the closer one.
                
                if (child1distance < child2distance) {
                    child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                    if (!child1intersects || child2distance < child1distance)
                        child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                }
                else {
                    child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
                    if (!child2intersects || child1distance < child2distance)
                        child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
                }
            }
            else
                child1intersects = child1->intersectsRay(mesh, origin,  direction, child1distance, child1face, child1uv);
        }
        else if (child2intersects)
            child2intersects = child2->intersectsRay(mesh, origin,  direction, child2distance, child2face, child2uv);
        
        // If either one had an intersection, return the closer one.
        
        if (child1intersects && (!child2intersects || child1distance < child2distance)) {
            distance = child1distance;
            face = child1face;
            uv = child1uv;
            return true;
        }
        if (child2intersects) {
            distance = child2distance;
            face = child2face;
            uv = child2uv;
            return true;
        }
        return false;
    }
    
    // This is a leaf node, so check each triangle for an intersection with the 
    // ray.
    
    bool foundIntersection = false;
    for (int i = 0; i < (int) triangles.size(); i++) {
        const UnitVec3& faceNormal = mesh.faces[triangles[i]].normal;
        double vd = ~faceNormal*direction;
        if (vd == 0.0)
            continue; // The ray is parallel to the plane.
        const Vec3& vert1 = mesh.vertices[mesh.faces[triangles[i]].vertices[0]].pos;
        double v0 = ~faceNormal*(vert1-origin);
        double t = v0/vd;
        if (t < 0.0)
            continue; // Ray points away from plane of triangle.
        if (foundIntersection && t >= distance)
            continue; // We already have a closer intersection.

        // Determine whether the intersection point is inside the triangle by projecting onto
        // a plane and computing the barycentric coordinates.

        Vec3 ri = origin+direction*t;
        const Vec3& vert2 = mesh.vertices[mesh.faces[triangles[i]].vertices[1]].pos;
        const Vec3& vert3 = mesh.vertices[mesh.faces[triangles[i]].vertices[2]].pos;
        int axis1, axis2;
        if (std::abs(faceNormal[1]) > std::abs(faceNormal[0])) {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[1])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 0;
                axis2 = 2;
            }
        }
        else {
            if (std::abs(faceNormal[2]) > std::abs(faceNormal[0])) {
                axis1 = 0;
                axis2 = 1;
            }
            else {
                axis1 = 1;
                axis2 = 2;
            }
        }
        Vec2 pos(ri[axis1]-vert1[axis1], ri[axis2]-vert1[axis2]);
        Vec2 edge1(vert1[axis1]-vert2[axis1], vert1[axis2]-vert2[axis2]);
        Vec2 edge2(vert1[axis1]-vert3[axis1], vert1[axis2]-vert3[axis2]);
        double denom = 1.0/(edge1%edge2);
        edge2 *= denom;
        double v = edge2%pos;
        if (v < 0.0 || v > 1.0)
            continue;
        edge1 *= denom;
        double w = pos%edge1;
        if (w < 0.0 || w > 1.0)
            continue;
        double u = 1.0-v-w;
        if (u < 0.0 || u > 1.0)
            continue;
        
        // It intersects.
        
        distance = t;
        face = triangles[i];
        uv = Vec2(u, v);
        foundIntersection = true;
    }
    return foundIntersection;
}




//==============================================================================
//                               OBB TREE NODE
//==============================================================================

ContactGeometry::TriangleMesh::OBBTreeNode::
OBBTreeNode(const OBBTreeNodeImpl& impl) : impl(&impl) {}

const OrientedBoundingBox& 
ContactGeometry::TriangleMesh::OBBTreeNode::getBounds() const {
    return impl->bounds;
}

bool ContactGeometry::TriangleMesh::OBBTreeNode::isLeafNode() const {
    return (impl->child1 == NULL);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getFirstChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child1, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child1);
}

const ContactGeometry::TriangleMesh::OBBTreeNode 
ContactGeometry::TriangleMesh::OBBTreeNode::getSecondChildNode() const {
    SimTK_ASSERT_ALWAYS(impl->child2, 
        "Called getFirstChildNode() on a leaf node");
    return OBBTreeNode(*impl->child2);
}

const Array_<int>& ContactGeometry::TriangleMesh::OBBTreeNode::
getTriangles() const {
    SimTK_ASSERT_ALWAYS(impl->child2 == NULL, 
        "Called getTriangles() on a non-leaf node");
    return impl->triangles;
}

int ContactGeometry::TriangleMesh::OBBTreeNode::getNumTriangles() const {
    return impl->numTriangles;
}



//==============================================================================
//                            TORUS & IMPL
//==============================================================================

ContactGeometry::Torus::Torus(Real torusRadius, Real tubeRadius)
:   ContactGeometry(new Torus::Impl(torusRadius, tubeRadius)) {}

/*static*/ ContactGeometryTypeId ContactGeometry::Torus::classTypeId()
{   return ContactGeometry::Torus::Impl::classTypeId(); }

Real ContactGeometry::Torus::getTorusRadius() const {
    return getImpl().getTorusRadius();
}

void ContactGeometry::Torus::setTorusRadius(Real radius) {
    updImpl().setTorusRadius(radius);
}

Real ContactGeometry::Torus::getTubeRadius() const {
    return getImpl().getTubeRadius();
}

void ContactGeometry::Torus::setTubeRadius(Real radius) {
    updImpl().setTubeRadius(radius);
}

const ContactGeometry::Torus::Impl& ContactGeometry::Torus::getImpl() const {
    assert(impl);
    return static_cast<const Torus::Impl&>(*impl);
}

ContactGeometry::Torus::Impl& ContactGeometry::Torus::updImpl() {
    assert(impl);
    return static_cast<Torus::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Torus::Impl::createDecorativeGeometry() const {
    PolygonalMesh mesh;
    createPolygonalMesh(mesh);
    return DecorativeMesh(mesh);
}

//TODO change to analytical solution
Vec3 ContactGeometry::Torus::Impl::
findNearestPoint(const Vec3& Q, bool& inside, UnitVec3& normal) const {
    // for now use local projection (TODO: not guaranteed to return the nearest
    // point)
    const Vec3 P = projectDownhillToNearestPoint(Q);
    inside = calcSurfaceValue(Q) > 0; // TODO: wrong sign convention
    normal = calcSurfaceUnitNormal(P);
    return P;
}

//TODO
bool ContactGeometry::Torus::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    std::cout << "ContactGeometry::Torus::Impl::intersectsRay unimplemented" << std::endl;
    return false;
}

void ContactGeometry::Torus::Impl::getBoundingSphere
    (Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = tubeRadius + torusRadius;
}

// Create a polygonal mesh for this torus using parameterization as follows:
// u = [0, 2*Pi] traces a circle in the x-y plane with radius torusRadius,
// which is the centroid of the torus. A point P on this circle is
// given by P = torusRadius*~[cos(u) sin(u) 0].
// v = [0, 2*Pi] traces a circle arond the cross-section (or tube) of the
// torus with radius tubeRadius, at a given u. A point Q on this circle
// is given by Q = (torusRadius + tubeRadius*cos(v))*e1 + tubeRadius*(~[0 0 1]*sin(v))
// where e1 = ~[sin(u) cos(u) 0]. The tube circle is in a plane spanned
// by e1 and the z-axis.
void ContactGeometry::Torus::Impl::createPolygonalMesh(PolygonalMesh& mesh) const {
    // TODO add resolution argument
    const int numSides = 12; //*resolution;
    const int numSlices = 36; //*resolution;   

    // add vertices 
    for (int i = 0; i < numSlices; ++i) {
      Real u = ((i*2*SimTK_PI)/numSlices);
      UnitVec3 e1(std::sin(u), std::cos(u), 0); // torus circle aligned with z-axis (z-axis through hole)
      for (int j = 0; j < numSides; ++j) {
        Real v = ((j*2*SimTK_PI)/numSides);
        Vec3 vtx = (torusRadius + tubeRadius*std::cos(v))*e1 + tubeRadius*std::sin(v)*Vec3(0,0,1); // use ZAXIS? 
        mesh.addVertex(vtx);  
      }
    }

    // add faces, be careful to wrap indices for the last slice
    int numVertices = mesh.getNumVertices();
    cout << "num verts = " << numVertices << endl;
    for (int i = 0; i < numVertices; ++i) {
//      cout << "v" << i << ": " << mesh.getVertexPosition(i) << endl;
      // define counter-clockwise quad faces
      Array_<int> faceIndices;
      faceIndices.push_back(i); // u_i,v_i
      faceIndices.push_back((i+1)%numVertices); // u_i, v_i+1
      faceIndices.push_back((i+1+numSides)%numVertices); // u_i+1, v_i+1
      faceIndices.push_back((i+numSides)%numVertices); // u_i+1, v_i
      mesh.addFace(faceIndices);	
    }

}

//TODO
void ContactGeometry::Torus::Impl::
calcCurvature(const Vec3& point, Vec2& curvature, Rotation& orientation) const {
    std::cout << "ContactGeometry::Torus::Impl::calcCurvature unimplemented" << std::endl;
}

//TODO
Vec3  ContactGeometry::Torus::Impl::
calcSupportPoint(UnitVec3 direction) const {
    std::cout << "ContactGeometry::Torus::Impl::calcSupportPoint unimplemented" << std::endl;
    return Vec3(0);
}

Real TorusImplicitFunction::
calcValue(const Vector& x) const {
    return 1-(square(ownerp->getTorusRadius()-std::sqrt(x[0]*x[0]+x[1]*x[1]))+x[2]*x[2])/
            square(ownerp->getTubeRadius());
}

Real TorusImplicitFunction::
calcDerivative(const Array_<int>& derivComponents, const Vector& x) const {
    // first derivatives
    if (derivComponents.size() == 1) {
        if (derivComponents[0]<2) {
            Real sqrt_xy = std::sqrt(x[0]*x[0] + x[1]*x[1]);
            return 2*x[derivComponents[0]]*(ownerp->getTorusRadius() - sqrt_xy)/
                    (square(ownerp->getTubeRadius())*sqrt_xy);
        }
        else
            return -2*x[2]/square(ownerp->getTubeRadius());
    }

    // second derivatives
    if (derivComponents.size() == 2) {
        if (derivComponents[0] < 2) { // fx_ fy_
            if (derivComponents[1] < 2) {
                Real tubeRadiusSq = square(ownerp->getTubeRadius());
                Real xy = x[0]*x[0] + x[1]*x[1];
                Real sqrt_xy = std::sqrt(xy);
                Real den = tubeRadiusSq*xy*sqrt_xy;
                if (derivComponents[0]==derivComponents[1]) { // fxx or fyy
                    int idx = derivComponents[1]==0; // if 0 then idx=1, if 1 then idx=0
                    Real num = 2*ownerp->getTorusRadius()*x[idx]*x[idx];
                    return num/den - 2/tubeRadiusSq;
                }
                else { // fxy or fyx
                    return - 2*ownerp->getTorusRadius()*x[0]*x[1]/den;
                }
            }
            else // fxz = fyz = 0
                return 0;
        }
        else { // fz_
            if (derivComponents[1] == 2) // fzz
                return -2/square(ownerp->getTubeRadius());
            else // fzx = fzy = 0
                return 0;
        }
    }

    //TODO higher order derivatives
    return 0;
}



//==============================================================================
//                                 GEODESIC
//==============================================================================

void Geodesic::dump(std::ostream& o) const {
    o << "Geodesic: " << getNumPoints() << " points, length=" 
                      << getLength() << "\n";
    bool hasQtoP = !directionalSensitivityQtoP.empty();
    if (!hasQtoP)
        o << "  QtoP Jacobian not available\n";
    for (int i=0; i < getNumPoints(); ++i) {
        o << "  Point at s=" << arcLengths[i] << ":\n";
        o << "    p=" << frenetFrames[i].p() 
          << " t=" << frenetFrames[i].x() << "\n";
        o << "    jP=" << directionalSensitivityPtoQ[i][0];
        if (hasQtoP)
          o << " jQ=" << directionalSensitivityQtoP[i][0];
        o << std::endl;
    }
}

//==============================================================================
//                      PARTICLE CON SURFACE SYSTEM GUTS
//==============================================================================


/*
 * This system is a 3d particle mass constrained to move along a surface
 * with no applied force (other than the constraint reaction force normal
 * to the surface). With no applied for the particle traces a geodesic along
 * the surface.
 *
 *
 * The DAE for a generic multibody system is:
 *       qdot = Nu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r be the 3d coordinates of the particle
 *       g(r) be the implicit surface function
 *       G be the gradient of the implicit surface function
 *       H be the hessian of the implicit surface function
 *
 * We will express implicit surface constraint as
 *                g(r) = 0    (perr)
 *                 Gr' = 0    (verr)
 *        Gr'' = -G'r' = - ~r'Hr' (aerr)
 *
 * So the matrix A = G and b = -r'Hr', and the
 * equations of motion are:
 *     [ M  ~G ] [ r'' ]   [     0    ]
 *     [ G   0 ] [  L  ] = [ - ~r'Hr' ]
 * where L (the Lagrange multiplier) is proportional to
 * the constraint force.
 *
 * solving for L,
 *        L = (GM\~G)\~r'Hr'
 *      r'' = - M\~G(GM\~G)\~r'Hr'
 *
 *      where ~A = transpose of A
 *        and A\ = inverse of A
 */
int ParticleConSurfaceSystemGuts::realizeTopologyImpl(State& s) const {
    // Generalized coordinates:
    //     q0, q1, q2, q3 = x,  y,  z,  j
    //     u0, u1, u2, u3 = x', y', z', j'
    const Vector init(4, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);
    return 0;
}
int ParticleConSurfaceSystemGuts::realizeModelImpl(State& s) const {
    return 0;
}
int ParticleConSurfaceSystemGuts::realizeInstanceImpl(const State& s) const {
    qerr0 = s.allocateQErr(subsysIndex, 1);
    uerr0 = s.allocateUErr(subsysIndex, 2);
    udoterr0 = s.allocateUDotErr(subsysIndex, 2); // and multiplier
//    event0 = s.allocateEvent(subsysIndex, Stage::Position, 3);
    return 0;
}
int ParticleConSurfaceSystemGuts::realizePositionImpl(const State& s) const {
    const Vector& q = s.getQ(subsysIndex);
    const Vec3&   point = Vec3::getAs(&q[0]);

    // This is the perr() equation that says the point must be on the surface.
    s.updQErr(subsysIndex)[0] = geom.calcSurfaceValue(point);
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    const Vec3&   p    = Vec3::getAs(&q[0]); // point
    const Vec3&   dpdt = Vec3::getAs(&u[0]); // point's velocity on surface

    Vector&       qdot = s.updQDot(subsysIndex);

    // Calculate qdots. They are just generalized speeds here.
    qdot = u;

    // These are the two verr() equations. The first is the derivative of
    // the point-on-surface holonomic constraint defined in realizePositionImpl
    // above. The second is a nonholonomic velocity constraint restricting
    // the velocity along the curve to be 1.
    s.updUErr(subsysIndex)[0]  = ~geom.calcSurfaceGradient(p)*dpdt;
    s.updUErr(subsysIndex)[1]  = dpdt.norm() - 1;
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeDynamicsImpl(const State& s) const {
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeAccelerationImpl(const State& s) const {

    // XXX assume unit mass

    const Vector& q    = s.getQ(subsysIndex);
    const Vec3&   p    = Vec3::getAs(&q[0]); // point
    const Real&   j    = q[3];                  // Jacobi field value

    const Vector& u    = s.getU(subsysIndex);
    const Vec3&   v    = Vec3::getAs(&u[0]); // point's velocity on surface

    Vector&       udot = s.updUDot(subsysIndex);
    Vec3&         a    = Vec3::updAs(&udot[0]);
    Real&         jdotdot = udot[3];            // Jacobi field 2nd derivative

    const Real  ep =  geom.calcSurfaceValue(p); // position constraint error
    const Vec3  GT =  geom.calcSurfaceGradient(p);
    const Real  ev =  ~GT*v;                    // d/dt ep
    const Mat33 H  =  geom.calcSurfaceHessian(p);
    const Real Gdotv = ~v*(H*v);
    const Real  L = Gdotv/(~GT*GT);
    a = GT*-L; // fills in udot

    // Now evaluate the Jacobi field.
    const Real Kg = geom.calcGaussianCurvature(p);
    jdotdot = -Kg*j;    // sets udot[3]

    // Qdotdots are just udots here.
    s.updQDotDot() = udot;

    s.updMultipliers(subsysIndex)[0] = L;

    // This is the aerr() equation.
    s.updUDotErr(subsysIndex)[0] = ~GT*a + Gdotv;
    s.updUDotErr(subsysIndex)[1] = (~v*a)/u.norm();
    return 0;
}

static Real wrms(const Vector& y, const Vector& w) {
    Real sumsq = 0;
    for (int i=0; i<y.size(); ++i)
        sumsq += square(y[i]*w[i]);
    return std::sqrt(sumsq/y.size());
}

/*
 * Here we want to remove any constraint errors from the current state,
 * and project out any component of the integrator's error estimate
 * perpendicular to the constraint manifold. We will do this sequentially
 * rather than handling position and velocity simultaneously.
 */
// qerrest is in/out
void ParticleConSurfaceSystemGuts::projectQImpl(State& s, Vector& qerrest,
                                const ProjectOptions& opts,
                                ProjectResults& results) const

{
    const Real consAccuracy = opts.getRequiredAccuracy();

    // These are convenient aliases for state entries. Note that they are
    // "live" since they refer directly into the state we are changing.

    const Vector& qerr = s.getQErr(subsysIndex);
    const Real&   ep   = qerr[0];

    Vector&       q    = s.updQ(subsysIndex);
    Vec3&         p    = Vec3::updAs(&q[0]);

    realize(s, Stage::Position); // recalc QErr (ep)
//    std::cout << "BEFORE wperr=" << ep << std::endl;

    Real qchg;
    Vec3 dp(0);
    int cnt = 0;
    do {
        Vec3 g = geom.calcSurfaceGradient(p);

        // dp = Pinv*ep, where Pinv = ~P/(P*~P) and P=~g
        Row3 P = ~g;
        dp = ~P/(P*~P)*ep;
//        dp = (g/(~g*g))*ep;

        qchg = std::sqrt(dp.normSqr()/3); // rms norm

        p -= dp; // updates the state

        s.invalidateAll(Stage::Position); // force realize position
        realize(s, Stage::Position); // recalc QErr (ep)

//        std::cout << cnt << ": AFTER q-=dq, q=" << q << ", wperr=" << ep << ", wqchg=" << qchg << std::endl;
        cnt++;

        //sleepInSec(0.5);
        if (cnt > 10) {
            results.setExitStatus(ProjectResults::FailedToConverge);
            return;
        }

    } while (std::abs(ep) > consAccuracy && qchg >= 0.01*consAccuracy);


    // Now do error estimates.
    if (qerrest.size()) {
        Vec3& eq = Vec3::updAs(&qerrest[0]);

        Vec3 g = geom.calcSurfaceGradient(p);

        // qperp = Pinv*P*eq, where Pinv = ~P/(P*~P) and P=~g
        Vec3 qperp = (g/(~g*g))*(~g*eq);

//        std::cout << "ERREST before=" << qerrest
//             << " wrms=" << std::sqrt(qerrest.normSqr()/q.size()) << std::endl;
//        std::cout << "P*eq=" << ~g*eq << std::endl;

        eq -= qperp;

//        std::cout << "ERREST after=" << qerrest
//                << " wrms=" << std::sqrt(qerrest.normSqr()/q.size()) << std::endl;
//        std::cout << "P*eq=" << ~g*eq << std::endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

// uerrest is in/out
void ParticleConSurfaceSystemGuts::projectUImpl(State& s, Vector& uerrest,
             const ProjectOptions& opts, ProjectResults& results) const
{
    const Real consAccuracy = opts.getRequiredAccuracy();

    const Vector& uerr = s.getUErr(subsysIndex); // set up aliases
    const Vec2&   ev   = Vec2::getAs(&uerr[0]);

    const Vector& q = s.getQ(subsysIndex);
    const Vec3&   p = Vec3::getAs(&q[0]);

    Vector&       u = s.updU(subsysIndex);
    Vec3&         v = Vec3::updAs(&u[0]);

    realize(s, Stage::Velocity); // calculate UErr (ev)
//    std::cout << "vBEFORE wperr=" << ep << std::endl;
//    std::cout << "vBEFORE wverr=" << ev << std::endl;

//    // Do velocity projection at current values of q, which should have
//    // been projected already.
    Vec3 g = geom.calcSurfaceGradient(p);

    // dv = Pinv*ev, where Pinv = ~P/(P*~P) and P=~g
//    Vec3 dv = (g/(~g*g))*ev[0];
    Row3 P = ~g;
    Vec3 dv = ~P/(P*~P)*ev[0];
//
//    // force unit speed
    const UnitVec3 newv(v - dv);
    v = Vec3(newv);
//    v -= dv;

    s.invalidateAll(Stage::Velocity); // force realize velocity
    realize(s, Stage::Velocity); // recalc UErr
//    std::cout << "vAFTER wverr=" << ev << std::endl;

// XXX combined u projection not working...

//    Real vchg;
//    Vec3 dv(0);
//    Mat23 V;
//    int cnt = 0;
//    V[0] = ~g; // Pnonholo = ~g;
//    do {
//        V[1] = ~v/v.norm(); // Vholo
//
//
//        // du = Vinv*ev, where Vinv = ~V/(V*~V) and V = ~[Pnonholo Vholo]
//        Mat22 VVT = V*~V;
//        VVT.invert();
//        dv = (~V)*VVT*ev;
//
////        std::cout << "V = " << V << std::endl;
////        std::cout << "v = " << v << ", dv = " << dv << std::endl;
//
////        dv = ~V/(V*~V)*ev;
//
//        vchg = std::sqrt(dv.normSqr()/3); // wrms norm
//
////        s.updU(subsysIndex)[0] -= du[0];
////        s.updU(subsysIndex)[1] -= du[1];
////        s.updU(subsysIndex)[2] -= du[2];
//
//        v -= dv;
//
//        s.invalidateAll(Stage::Velocity); // force realize velocity
//        realize(s, Stage::Velocity); // recalc UErr (ev)
//
//        std::cout << cnt << ": AFTER v-=dv verr=" << ev << " vchg=" << vchg << std::endl;
////        std::cout << "v = " << v << ", dv = " << dv << std::endl;
//
//        cnt++;
//
//        sleep(0.5);
//        if (cnt > 10) {
//            results.setExitStatus(ProjectResults::FailedToConverge);
//            return;
//        }
//
//
//    } while (std::max(std::abs(ev[0]), std::abs(ev[1])) > consAccuracy &&
//             vchg >= 0.01*consAccuracy);

    // Now do error estimates.
    if (uerrest.size()) {
        Vec3& eu = Vec3::updAs(&uerrest[0]);

        // uperp = Vinv*(V*eu), where V=P, Pinv = ~P/(P*~P) and P=~n
        Vec3 uperp = (g/(~g*g))*(~g*eu);

//        std::cout << "ERREST before=" << uerrest
//             << " wrms=" << std::sqrt(uerrest.normSqr()/u.size())  << std::endl;
//        std::cout << " VW*eu=" << ~n*eu << std::endl;

        eu -= uperp;

//        std::cout << "ERREST after=" << uerrest
//             << " wrms=" << std::sqrt(uerrest.normSqr()/u.size())  << std::endl;
//        std::cout << " VW*eu=" << ~n*eu << std::endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
//    std::cout << "norm(u) = " << s.getU(subsysIndex).norm() << std::endl;
}




//==============================================================================
//                      PARTICLE ON SURFACE SYSTEM GUTS
//==============================================================================
/*
 * This system is a 3d particle mass constrained to move along a surface
 * with no applied force (other than the constraint reaction force normal
 * to the surface). With no applied for the particle traces a geodesic along
 * the surface.
 *
 * We also integrate to find the directional sensitivity of the geodesic
 * endpoint with respect to an angular perturbation in the starting direction.
 * This is the amplitude of a particular Jacobi field along the curve defined
 * by this equation:
 *      j'' + Kg*j = 0,   with j(0)=0, j'(0)=1
 * where j=j(s), Kg=Kg(s) is the Gaussian curvature of the surface evaluated at
 * s along the curve.
 *
 * NOTE: these equations are only valid if the independent variable s is
 * arc length along the geodesic curve; don't mess with the parameterization!
 * Ref: Do Carmo, M.P. 1976 Differential Geometry of Curves and Surfaces,
 * Chapter 5-5 Jacobi Fields and Conjugate Points.
 * For real understanding, see Andreas Scholz' master's thesis (2012).
 *
 * The DAE for a generic multibody system is:
 *       qdot = Nu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r be the 3d coordinates of the particle
 *       g(r) be the implicit surface function
 *       G be the gradient of the implicit surface function
 *       H be the hessian of the implicit surface function
 *
 * We will express implicit surface constraint as
 *                g(r) = 0    (perr)
 *                 Gr' = 0    (verr)
 *        Gr'' = -G'r' = - ~r'Hr' (aerr)
 *
 * So the matrix A = G and b = -r'Hr', and the
 * equations of motion are:
 *     [ M  ~G ] [ r'' ]   [     0    ]
 *     [ G   0 ] [  L  ] = [ - ~r'Hr' ]
 * where L (the Lagrange multiplier) is proportional to
 * the constraint force.
 *
 * solving for L,
 *        L = (GM\~G)\~r'Hr'
 *      r'' = - M\~G(GM\~G)\~r'Hr'
 *
 *      where ~A = transpose of A
 *        and A\ = inverse of A
 */
int ParticleOnSurfaceSystemGuts::realizeTopologyImpl(State& s) const {
    // Generalized coordinates:
    //     q0, q1, q2, q3 = x,  y,  z,  j
    //     u0, u1, u2, u3 = x', y', z', j'        
    const Vector init(4, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeModelImpl(State& s) const {
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeInstanceImpl(const State& s) const {
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizePositionImpl(const State& s) const {
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot = u;
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeDynamicsImpl(const State& s) const {
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeAccelerationImpl(const State& s) const {

    // XXX assume unit mass

    const Vector& q    = s.getQ(subsysIndex);
    const Vec3&   p    = Vec3::getAs(&q[0]);    // point
    const Real&   j    = q[3];                  // Jacobi field value

    const Vector& u    = s.getU(subsysIndex);
    const Vec3&   v    = Vec3::getAs(&u[0]);    // point velocity in S

    Vector&       udot = s.updUDot(subsysIndex);
    Vec3&         a    = Vec3::updAs(&udot[0]); // point acceleration in S
    Real&         jdotdot = udot[3];            // Jacobi field 2nd derivative

    const Real  ep =  geom.calcSurfaceValue(p); // position constraint error
    const Vec3  GT =  geom.calcSurfaceGradient(p);
    const Real  ev =  ~GT*v;                    // d/dt ep
    const Mat33 H  =  geom.calcSurfaceHessian(p);
    const Real Gdotv = ~v*(H*v);
    const Real L = (Gdotv + beta*ev + alpha*ep)/(~GT*GT);
    a = GT*-L;          // sets udot[0..2]

    // Now evaluate the Jacobi field.
    const Real Kg = geom.calcGaussianCurvature(p);
    jdotdot = -Kg*j;    // sets udot[3]

    // qdotdot is just udot.
    s.updQDotDot() = udot;
    return 0;
}








