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
#include "simmath/TimeStepper.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;

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

bool ContactGeometry::intersectsRay(const Vec3& origin, const UnitVec3& direction, Real& distance, UnitVec3& normal) const {
    return getImpl().intersectsRay(origin, direction, distance, normal);
}

Vec3 ContactGeometry::findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    return getImpl().findNearestPoint(position, inside, normal);
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

Real ContactGeometry::
calcSurfaceValue(const Vector& point) const {
    return getImpl().getImplicitFunction().calcValue(point);
}

Vec3 ContactGeometry::
calcSurfaceNormal(const Vector& point) const {
    const Function& f = getImpl().getImplicitFunction();

    // arguments to get first derivative from the calcDerivative interface
    Array_<int> fx(1, 0);
    Array_<int> fy(1, 1);
    Array_<int> fz(1, 2);

    Vec3 normal(3);
    // implicit surfaces are defined as positive inside and negative outside
    // therefore normal is the negative of the gradient
    normal[0] = -f.calcDerivative(fx,point);
    normal[1] = -f.calcDerivative(fy,point);
    normal[2] = -f.calcDerivative(fz,point);
    return normal;
}

Mat33 ContactGeometry::
calcSurfaceHessian(const Vector& point) const {
    const Function& f = getImpl().getImplicitFunction();

    // arguments to get second derivatives from the calcDerivative interface
    const int xx[] = {0,0}; Array_<int> fxx(xx,xx+2);
    const int xy[] = {0,1}; Array_<int> fxy(xy,xy+2);
    const int xz[] = {0,2}; Array_<int> fxz(xz,xz+2);
    const int yx[] = {1,0}; Array_<int> fyx(yx,yx+2);
    const int yy[] = {1,1}; Array_<int> fyy(yy,yy+2);
    const int yz[] = {1,2}; Array_<int> fyz(yz,yz+2);
    const int zx[] = {2,0}; Array_<int> fzx(zx,zx+2);
    const int zy[] = {2,1}; Array_<int> fzy(zy,zy+2);
    const int zz[] = {2,2}; Array_<int> fzz(zz,zz+2);

    Mat33 hess;

    hess(0,0) = f.calcDerivative(fxx,point);
    hess(0,1) = f.calcDerivative(fxy,point);
    hess(0,2) = f.calcDerivative(fxz,point);
    hess(1,0) = f.calcDerivative(fyx,point);
    hess(1,1) = f.calcDerivative(fyy,point);
    hess(1,2) = f.calcDerivative(fyz,point);
    hess(2,0) = f.calcDerivative(fzx,point);
    hess(2,1) = f.calcDerivative(fzy,point);
    hess(2,2) = f.calcDerivative(fzz,point);

    return hess;
}

Vec3 ContactGeometry::calcSupportPoint(UnitVec3 direction) const 
{   return getImpl().calcSupportPoint(direction); }

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
void ContactGeometry::continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
        const GeodesicOptions& options, Geodesic& geod) {
    getImpl().continueGeodesic(xP, xQ, prevGeod, options, geod);
}


// Utility method to find geodesic between P and Q
// with starting directions tPhint and tQhint
// XXX tangent basis should be formed from previous geodesic
void ContactGeometry::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {
    getImpl().calcGeodesic(xP, xQ, tPhint, tQhint, geod);
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


Vec2 ContactGeometry::calcGeodError(const Vec3& xP, const Vec3& xQ,
        const UnitVec3& tP, const UnitVec3& tQ,
        Geodesic* geod) const {
    return getImpl().calcGeodError(xP, xQ, tP, tQ, geod);
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
//             GEODESIC EVALUATORS in CONTACT GEOMETRY IMPL
//==============================================================================

Real ContactGeometryImpl::calcSurfaceValue(const Vector& point) const {
    return myHandle->calcSurfaceValue(point);
}

Vec3 ContactGeometryImpl::calcSurfaceNormal(const Vector& point) const {
    return myHandle->calcSurfaceNormal(point);
}

Mat33 ContactGeometryImpl::calcSurfaceHessian(const Vector& point) const {
    return myHandle->calcSurfaceHessian(point);
}

const Real estimatedGeodesicAccuracy = 1e-12; // used in numerical differentiation
const Real pauseBetweenGeodIterations = 0; // sec, used in newton solver


// This local class is used to Calcualte the split geodesic error
//  given scalar angles at P and Q
class ContactGeometryImpl::SplitGeodesicError: public Differentiator::JacobianFunction {

public:
    SplitGeodesicError(int nf, int ny, ContactGeometry& geom,
            const Vec3& xP, const Vec3& xQ, const Vec3& tPhint, const Vec3& tQhint) :
            Differentiator::JacobianFunction(nf, ny),
                    geom(geom),
                    P(xP), Q(xQ),
                    R_SP(geom.calcTangentBasis(P, tPhint)),
                    R_SQ(geom.calcTangentBasis(Q, tQhint)) { }

    // x = ~[thetaP, thetaQ]
    int f(const Vector& x, Vector& fx) const  {
        UnitVec3 tP = calcUnitTangentVec(x[0], R_SP);
        UnitVec3 tQ = calcUnitTangentVec(x[1], R_SQ);

        Vec2 geodErr = geom.calcGeodError(P, Q, tP, tQ);

        // error between geodesic end points at plane
        fx[0] = geodErr[0];
        fx[1] = geodErr[1];
        return 0;
    }


private:
    ContactGeometry& geom;
    const Vec3& P;
    const Vec3& Q;
    const Rotation R_SP;
    const Rotation R_SQ;

}; // class SplitGeodesicError


void ContactGeometryImpl::initGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& xSP, const GeodesicOptions& options, Geodesic& geod) const
{
    // TODO
}


// Given two points and previous geodesic curve close to the points, find
// a geodesic curve connecting the points that is close to the previous geodesic.
void ContactGeometryImpl::continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
        const GeodesicOptions& options, Geodesic& geod) const {

    // XXX could also estimate P and Q based on prevGeod's contact point velocities

    // set tP and tQ hints based on previous geodesic's endpoint tangents
    Vec3 tPhint = prevGeod.getTangents()[0]; // prevGeod.getTangent(0);
    Vec3 tQhint = prevGeod.getTangents()[prevGeod.getTangents().size()-1]; // prevGeod.getTangent(1);

    calcGeodesic(xP, xQ, tPhint, tQhint, geod);
}



// Utility method to used by calcGeodesicInDirectionUntilPlaneHit
// and calcGeodesicInDirectionUntilLengthReached
void ContactGeometryImpl::
shootGeodesicInDirection(const Vec3& P, const UnitVec3& tP,
        const Real& finalTime, const GeodesicOptions& options,
        Geodesic& geod) const {

    // integrator settings
    const Real startTime = 0;
    const Real integratorAccuracy = 1e-6;
    const Real integratorConstraintTol = 1e-6;

    ++numGeodesicsShot;

    // Initialize state
    State sysState = ptOnSurfSys->getDefaultState();
    sysState.setTime(startTime);
    Vector& q = sysState.updQ();
    Vector& u = sysState.updU();
    q[0] = P[0]; q[1] = P[1]; q[2] = P[2];
    u[0] = tP[0]; u[1] = tP[1]; u[2] = tP[2];

    // Setup integrator to integrate until terminatingLength
    //RungeKutta3Integrator integ(ptOnSurfSys);
    RungeKuttaMersonIntegrator integ(*ptOnSurfSys);
    integ.setAccuracy(integratorAccuracy);
    integ.setConstraintTolerance(integratorConstraintTol);
    integ.setFinalTime(finalTime);
    integ.setReturnEveryInternalStep(true); // save geodesic knot points

    // Setup timestepper in order to handle event when geodesic hits the plane
    TimeStepper ts(*ptOnSurfSys, integ);
    ts.setReportAllSignificantStates(true);
    ts.initialize(sysState);
    const State& state = integ.getState(); // integrator state

    // Simulate it, and record geodesic knot points after each step
    // Terminate when geodesic hits the plane
    Integrator::SuccessfulStepStatus status;
    int stepcnt = 0;
    while (true) {
        // Final time is already reported by the time we see end of simulation;
        // don't duplicate the last step.
        if ((status=ts.stepTo(Infinity)) == Integrator::EndOfSimulation)
            break;

        const Real s = state.getTime();
        geod.addPoint(Vec3(&state.getQ()[0]));
        geod.addTangent(Vec3(&state.getU()[0]));
        geod.addArcLength(s);

        ++stepcnt;
    }
    geod.setIsConvex(false); // TODO
    geod.setIsShortest(false); // TODO
    geod.setAchievedAccuracy(integratorAccuracy); // TODO: accuracy of length?
    // TODO: better to use something like the second-to-last step, or average
    // excluding initial and last steps, so that we don't have to start small.
    geod.setInitialStepSizeHint(integ.getActualInitialStepSizeTaken());
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
    R_SP = myHandle->calcTangentBasis(xP, tPhint);
    R_SQ = myHandle->calcTangentBasis(xQ, tQhint);

    // calculate plane bisecting P and Q, and use as termination condition for integrator
    UnitVec3 normal(xQ - xP);
    Real offset = (~(xP+xQ)*normal)/2 ;
    geodHitPlaneEvent->setPlane(Plane(normal, offset));

    Mat22 J;
    Vec2 x, xold, dx, Fx;

//    Matrix J(2,2);
//    Vector x(2), xold(2), dx(2), Fx(2);


    // initial conditions
    x[0] = Pi/2; // thetaP
    x[1] = Pi/2; // thetaQ

    Real f, fold, lam = 1;

//    splitGeodErr = new SplitGeodesicError(2, 2, *const_cast<ContactGeometry*>(this),
//            xP, xQ, tPhint, tQhint);
//    splitGeodErr->setEstimatedAccuracy(estimatedGeodesicAccuracy);
//    Differentiator diff( *const_cast<SplitGeodesicError*>(splitGeodErr));

//    splitGeodErr->f(x, Fx);
    Fx = calcGeodError(xP, xQ, x[0], x[1]);
    if (vizReporter != NULL) {
        vizReporter->handleEvent(ptOnSurfSys->getDefaultState());
        sleepInSec(pauseBetweenGeodIterations);
    }

    f = std::sqrt(~Fx*Fx);

    for (int i = 0; i < maxNewtonIterations; ++i) {
        if (maxabs(Fx) < ftol) {
            std::cout << "geodesic converged in " << i << " iterations" << std::endl;
//            std::cout << "err = " << Fx << std::endl;
            break;
        }
//        diff.calcJacobian(x,  Fx, J, Differentiator::ForwardDifference);
        J = calcGeodErrorJacobian(xP, xQ, x[0], x[1],
                Differentiator::CentralDifference);
        dx = J.invert()*Fx;

        fold = f;
        xold = x;

        // backtracking
        lam = 1;
        while (true) {
            x = xold - lam*dx;
//            splitGeodErr->f(x, Fx);
            Fx = calcGeodError(xP, xQ, x[0], x[1]);
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
//    std::cout << "err = " << Fx << std::endl;

    mergeGeodesics(geodP, geodQ, geod);
}


// Calculate the "geodesic error" for thetaP and thetaQ, and return the
// resulting (kinked) geodesic if the supplied pointer is non-null.
Vec2 ContactGeometryImpl::
calcGeodError(const Vec3& xP, const Vec3& xQ,
              const Real thetaP, const Real thetaQ,
              Geodesic* geodesic) const
{
    UnitVec3 tP = calcUnitTangentVec(thetaP, R_SP);
    UnitVec3 tQ = calcUnitTangentVec(thetaQ, R_SQ);
    return calcGeodError(xP, xQ, tP, tQ, geodesic);
}

// Calculate the "geodesic error" for tP and tQ
Vec2 ContactGeometryImpl::
calcGeodError(const Vec3& xP, const Vec3& xQ,
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

    if (geodesic)
        mergeGeodesics(geodP, geodQ, *geodesic);

    return calcError(geodP, geodQ);
}


// Calculate the "geodesic jacobian" by numerical perturbation
Mat22 ContactGeometryImpl::
calcGeodErrorJacobian(const Vec3& xP, const Vec3& xQ,
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
    Vec3 Phat = geodP.getPoints().back();
    Vec3 Qhat = geodQ.getPoints().back();
    UnitVec3 tPhat(geodP.getTangents().back());
    UnitVec3 tQhat(geodQ.getTangents().back());

    UnitVec3 nPhat(calcSurfaceNormal((Vector) Phat));
    UnitVec3 nQhat(calcSurfaceNormal((Vector) Qhat));
    UnitVec3 bPhat(tPhat % nPhat);
    UnitVec3 bQhat(tQhat % nQhat);

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

// First deriv with respect to z (component 2) is -1, all higher derivs are 0.
// Higher partials involving z are zero.
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

    const Vector init(3, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int ParticleConSurfaceSystemGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int ParticleConSurfaceSystemGuts::realizeInstanceImpl(const State& s) const {
    qerr0 = s.allocateQErr(subsysIndex, 1);
    uerr0 = s.allocateUErr(subsysIndex, 2);
    udoterr0 = s.allocateUDotErr(subsysIndex, 2); // and multiplier
//    event0 = s.allocateEvent(subsysIndex, Stage::Position, 3);

    System::Guts::realizeInstanceImpl(s);
    return 0;
}
int ParticleConSurfaceSystemGuts::realizePositionImpl(const State& s) const {
    const Vector& q = s.getQ(subsysIndex);
    // This is the perr() equation.
    s.updQErr(subsysIndex)[0] = geom.calcSurfaceValue(q);

    System::Guts::realizePositionImpl(s);
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vec3&   u    = Vec3::getAs(&s.getU(subsysIndex)[0]);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];
    qdot[2] = u[2];

    // This is the verr() equation.
    s.updUErr(subsysIndex)[0]  = -(~geom.calcSurfaceNormal(q)*u);
    s.updUErr(subsysIndex)[1]  = u.norm() - 1;

    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeDynamicsImpl(const State& s) const {

    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int ParticleConSurfaceSystemGuts::realizeAccelerationImpl(const State& s) const {

    // XXX assume unit mass

    const Vector& q    = s.getQ(subsysIndex);
    const Vec3&   u    = Vec3::getAs(&s.getU(subsysIndex)[0]);
    Vector&       udot = s.updUDot(subsysIndex);
    Vec3 a(0);

    Real g = geom.calcSurfaceValue(q);
    Vec3 GT = geom.calcSurfaceNormal(q);
    Mat33 H = geom.calcSurfaceHessian(q);
    Real Gdotu = ~u*(H*u);
//    Real L = (Gdotu + beta*~GT*v + alpha*g)/(~GT*GT); // Baumgarte stabilization
    Real L = Gdotu;
    a = GT*-L;

    udot[0] = a[0];
    udot[1] = a[1];
    udot[2] = a[2];

    s.updQDotDot() = udot;

    s.updMultipliers(subsysIndex)[0] = L;

    // This is the aerr() equation.
    s.updUDotErr(subsysIndex)[0] = ~GT*a + Gdotu;
    s.updUDotErr(subsysIndex)[1] = (~u*a)/u.norm();

    System::Guts::realizeAccelerationImpl(s);
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
    const Vec3& q = Vec3::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    Real& ep = s.updQErr(subsysIndex)[0];

    realize(s, Stage::Position); // recalc QErr (ep)
//    std::cout << "BEFORE wperr=" << ep << std::endl;

    Real qchg;
    Vec3 dq(0);
    int cnt = 0;
    do {
        Vec3 n = geom.calcSurfaceNormal((Vector)q);

        // dq = Pinv*ep, where Pinv = ~P/(P*~P) and P=-(~n)
        dq = -(n/(~n*n))*ep;

        qchg = std::sqrt(dq.normSqr()/q.size()); // wrms norm

        s.updQ(subsysIndex)[0] -= dq[0];
        s.updQ(subsysIndex)[1] -= dq[1];
        s.updQ(subsysIndex)[2] -= dq[2];

        realize(s, Stage::Position); // recalc QErr (ep)

//        std::cout << cnt << ": AFTER q-=dq wperr=" << ep << " wqchg=" << qchg << std::endl;
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

        Vec3 n = geom.calcSurfaceNormal((Vector)q);

        // qperp = Pinv*P*eq, where Pinv = ~P/(P*~P) and P=-(~n)
        Vec3 qperp = -(n/(~n*n))*(~n*eq);

//        std::cout << "ERREST before=" << qerrest
//             << " wrms=" << std::sqrt(qerrest.normSqr()/q.size()) << std::endl;
//        std::cout << "P*eq=" << ~n*eq << std::endl;

        eq -= qperp;

//        std::cout << "ERREST after=" << qerrest
//                << " wrms=" << std::sqrt(qerrest.normSqr()/q.size()) << std::endl;
//        std::cout << "P*eq=" << ~n*eq << std::endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

// uerrest is in/out
void ParticleConSurfaceSystemGuts::projectUImpl(State& s, Vector& uerrest,
             const ProjectOptions& opts, ProjectResults& results) const
{
    const Real consAccuracy = opts.getRequiredAccuracy();

    const Vector& q = s.getQ(subsysIndex); // set up aliases
    const Vec3& u = Vec3::getAs(&s.getU(subsysIndex)[0]);
    Vec2& ev = Vec2::updAs(&s.updUErr(subsysIndex)[0]);

//    const Real& ep = s.getQErr(subsysIndex)[0];

    realize(s, Stage::Velocity); // calculate UErr (ev)
//    std::cout << "vBEFORE wperr=" << ep << std::endl;
//    std::cout << "vBEFORE wverr=" << ev << std::endl;

    // Do velocity projection at current values of q, which should have
    // been projected already.
    Vec3 n = geom.calcSurfaceNormal((Vector)q);

    // du = Pinv*ev, where Pinv = ~P/(P*~P) and P=-(~n)
//    Vec3 du = -(n/(~n*n))*ev[0];
    Row3 P = -(~n);
    Vec3 du = ~P/(P*~P)*ev[0];

//    s.updU(subsysIndex)[0] -= du[0];
//    s.updU(subsysIndex)[1] -= du[1];
//    s.updU(subsysIndex)[2] -= du[2];
//
//    // force unit speed
    UnitVec3 newu(u - du);
    s.updU(subsysIndex)[0] = newu[0];
    s.updU(subsysIndex)[1] = newu[1];
    s.updU(subsysIndex)[2] = newu[2];

    realize(s, Stage::Velocity); // recalc UErr
//    std::cout << "vAFTER wverr=" << ev << std::endl;

// XXX combined u projection not working...

//    Real uchg;
//    Vec3 du(0);
//    Mat23 V;
//    int cnt = 0;
//    V[0] = -(~n);
//    do {
//        V[1] = ~u/u.norm();
//
//        // du = Vinv*ev, where Vinv = ~V/(V*~V) and V = ~[Pnonholo Vholo]
//        Mat22 VVT = V*~V;
//        VVT.invert();
//        du = (~V)*VVT*ev;
//
////        Vec3 du = ~V/(V*~V)*ev;
//
//        uchg = std::sqrt(du.normSqr()/u.size()); // wrms norm
//
//        s.updU(subsysIndex)[0] -= du[0];
//        s.updU(subsysIndex)[1] -= du[1];
//        s.updU(subsysIndex)[2] -= du[2];
//
//        realize(s, Stage::Velocity); // recalc UErr (ev)
//
//        std::cout << cnt << ": AFTER u-=uq verr=" << ev << " uchg=" << uchg << std::endl;
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
//             uchg >= 0.01*consAccuracy);

    // Now do error estimates.
    if (uerrest.size()) {
        Vec3& eu = Vec3::updAs(&uerrest[0]);

        // uperp = Vinv*(V*eu), where V=P, Pinv = ~P/(P*~P) and P=-(~n)
        Vec3 uperp = -(n/(~n*n))*(~n*eu);

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

    const Vector init(3, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizeInstanceImpl(const State& s) const {

    System::Guts::realizeInstanceImpl(s);
    return 0;
}
int ParticleOnSurfaceSystemGuts::realizePositionImpl(const State& s) const {

    System::Guts::realizePositionImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];
    qdot[2] = u[2];

    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeDynamicsImpl(const State& s) const {

    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int ParticleOnSurfaceSystemGuts::realizeAccelerationImpl(const State& s) const {

    // XXX assume unit mass

    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       udot = s.updUDot(subsysIndex);

    Vec3 v(u[0], u[1], u[2]);
    Vec3 a(0);

    Real g = geom.calcSurfaceValue(q);
    Vec3 GT = -geom.calcSurfaceNormal(q);
    Mat33 H = geom.calcSurfaceHessian(q);
    Real Gdotu = ~v*(H*v);
    Real L = (Gdotu + beta*~GT*v + alpha*g)/(~GT*GT);
    a = GT*-L;

    udot[0] = a[0]; udot[1] = a[1]; udot[2] = a[2];
    s.updQDotDot() = udot;

    System::Guts::realizeAccelerationImpl(s);
    return 0;
}








