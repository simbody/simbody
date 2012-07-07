/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKmath                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness, Michael Sherman                                     *
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

/*
 * This is a temporary class for the geodesic specific
 * additions to ContactGeometry
 */

#include "simmath/internal/GeodesicGeometry.h"
#include "simmath/internal/ParticleConSurfaceSystem.h"

using namespace SimTK;
using std::cout;
using std::endl;


//==============================================================================
//                          GEODESIC GEOMETRY
//==============================================================================

// This local class is used to Calcualte the split geodesic error
//  given scalar angles at P and Q
class GeodesicGeometry::SplitGeodesicError: public Differentiator::JacobianFunction {

public:
    SplitGeodesicError(int nf, int ny, GeodesicGeometry& geodgeom,
            const Vec3& xP, const Vec3& xQ, const Vec3& tPhint, const Vec3& tQhint) :
            Differentiator::JacobianFunction(nf, ny),
                    gg(geodgeom),
                    P(xP), Q(xQ),
                    R_SP(calcTangentBasis(P, tPhint, gg.geom)),
                    R_SQ(calcTangentBasis(Q, tQhint, gg.geom)) { }

    // x = ~[thetaP, thetaQ]
    int f(const Vector& x, Vector& fx) const  {
        UnitVec3 tP = calcUnitTangentVec(x[0], R_SP);
        UnitVec3 tQ = calcUnitTangentVec(x[1], R_SQ);

        Vec2 geodErr = gg.calcGeodError(P, Q, tP, tQ);

        // error between geodesic end points at plane
        fx[0] = geodErr[0];
        fx[1] = geodErr[1];
        return 0;
    }


private:
    GeodesicGeometry& gg;
    const Vec3& P;
    const Vec3& Q;
    const Rotation R_SP;
    const Rotation R_SQ;

}; // class SplitGeodesicError




void GeodesicGeometry::initGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& xSP, const GeodesicOptions& options, Geodesic& geod) const
{
    // TODO
}


// Given two points and previous geodesic curve close to the points, find
// a geodesic curve connecting the points that is close to the previous geodesic.
void GeodesicGeometry::continueGeodesic(const Vec3& xP, const Vec3& xQ, const Geodesic& prevGeod,
        const GeodesicOptions& options, Geodesic& geod) {

    // XXX could also estimate P and Q based on prevGeod's contact point velocities

    // set tP and tQ hints based on previous geodesic's endpoint tangents
    Vec3 tPhint = prevGeod.getTangents()[0]; // prevGeod.getTangent(0);
    Vec3 tQhint = prevGeod.getTangents()[prevGeod.getTangents().size()-1]; // prevGeod.getTangent(1);

    calcGeodesic(xP, xQ, tPhint, tQhint, geod);
}



// Utility method to used by calcGeodesicInDirectionUntilPlaneHit
// and calcGeodesicInDirectionUntilLengthReached
void GeodesicGeometry::
shootGeodesicInDirection(const Vec3& P, const UnitVec3& tP,
        const Real& finalTime, const GeodesicOptions& options,
        Geodesic& geod) const {

    ++numGeodesicsShot;

    // Initialize state
    State sysState = ptOnSurfSys.getDefaultState();
    sysState.setTime(startTime);
    Vector& q = sysState.updQ();
    Vector& u = sysState.updU();
    q[0] = P[0]; q[1] = P[1]; q[2] = P[2];
    u[0] = tP[0]; u[1] = tP[1]; u[2] = tP[2];

    // Setup integrator to integrate until terminatingLength
    RungeKutta3Integrator integ(ptOnSurfSys);
    integ.setAccuracy(integratorAccuracy);
    integ.setConstraintTolerance(integratorConstraintTol);
    integ.setFinalTime(finalTime);
    integ.setReturnEveryInternalStep(true); // save geodesic knot points

    // Setup timestepper in order to handle event when geodesic hits the plane
    TimeStepper ts(ptOnSurfSys, integ);
    ts.setReportAllSignificantStates(true);
    ts.initialize(sysState);
    const State& state = integ.getState(); // integrator state

    // Simulate it, and record geodesic knot points after each step
    // Terminate when geodesic hits the plane
    Integrator::SuccessfulStepStatus status;
    int stepcnt = 0;
    while (true) {
        status = ts.stepTo(Infinity);
        const Real T = integ.getTime();
        geod.addPoint(Vec3(&state.getQ()[0]));
        geod.addTangent(Vec3(&state.getU()[0]));
        geod.addArcLength(T);

        if (status == Integrator::EndOfSimulation) {
            break;
        }
        ++stepcnt;
    }
}


// Compute a geodesic curve starting at the given point, starting in the given
// direction, and terminating at the given plane.
// XXX what to do if tP is not in the tangent plane at P -- project it?
// XXX what to do if we don't hit the plane
void GeodesicGeometry::
shootGeodesicInDirectionUntilPlaneHit(const Vec3& xP, const UnitVec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {
    geodHitPlaneEvent->setEnabled(true);
    geodHitPlaneEvent->setPlane(terminatingPlane);
    shootGeodesicInDirection(xP, tP, Infinity, options, geod);
}


// Compute a geodesic curve of the given length, starting at the given point and
// in the given direction.
// XXX what to do if tP is not in the tangent plane at P -- project it?
void GeodesicGeometry::
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
void GeodesicGeometry::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {

    // reset counter
    numGeodesicsShot = 0;

    // define basis
    R_SP = calcTangentBasis(xP, tPhint, geom);
    R_SQ = calcTangentBasis(xQ, tQhint, geom);

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

//    splitGeodErr = new SplitGeodesicError(2, 2, *const_cast<GeodesicGeometry*>(this),
//            xP, xQ, tPhint, tQhint);
//    splitGeodErr->setEstimatedAccuracy(estimatedGeodesicAccuracy);
//    Differentiator diff( *const_cast<SplitGeodesicError*>(splitGeodErr));

//    splitGeodErr->f(x, Fx);
    Fx = calcGeodError(xP, xQ, x[0], x[1]);
    if (vizReporter != NULL) {
        vizReporter->handleEvent(ptOnSurfSys.getDefaultState());
        usleep((useconds_t)(pauseBetweenGeodIterations*1000000));
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
            vizReporter->handleEvent(ptOnSurfSys.getDefaultState());
            usleep((useconds_t)(pauseBetweenGeodIterations*1000000));
        }

    }
//    std::cout << "err = " << Fx << std::endl;

    mergeGeodesics(geodP, geodQ, geod);
}


// Calculate the "geodesic error" for thetaP and thetaQ
Vec2 GeodesicGeometry::
calcGeodError(const Vec3& xP, const Vec3& xQ, const Real thetaP, const Real thetaQ) const {

    UnitVec3 tP = calcUnitTangentVec(thetaP, R_SP);
    UnitVec3 tQ = calcUnitTangentVec(thetaQ, R_SQ);
    return calcGeodError(xP, xQ, tP, tQ);
}

// Calculate the "geodesic error" for tP and tQ
Vec2 GeodesicGeometry::
calcGeodError(const Vec3& xP, const Vec3& xQ, const UnitVec3& tP, const UnitVec3& tQ) const {

    geodP.clear();
    geodQ.clear();

    GeodesicOptions opts;
    shootGeodesicInDirectionUntilPlaneHit(xP, tP,
            geodHitPlaneEvent->getPlane(), opts, geodP);
    shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
            geodHitPlaneEvent->getPlane(), opts, geodQ);

    return calcError(geom, geodP, geodQ);
}


// Calculate the "geodesic jacobian" by numerical perturbation
Mat22 GeodesicGeometry::
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

    Vec2 fy0 = calcError(geom, geodP, geodQ);

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
    fyptmp = calcError(geom, geodPtmp, geodQ);

    if (order==1) {
        dfdy(0) = (fyptmp-fy0)/h;
    } else {

        geodPtmp.clear();
        tP = calcUnitTangentVec(thetaP-h, R_SP);
        shootGeodesicInDirectionUntilPlaneHit(xP, tP,
                geodHitPlaneEvent->getPlane(), opts, geodPtmp);
        fymtmp = calcError(geom, geodPtmp, geodQ);

        dfdy(0) = (fyptmp-fymtmp)/(2*h);
    }

    // perturb thetaQ
    hEst = accFactor*std::max(std::abs(thetaQ), YMin);
    h = cleanUpH(hEst, thetaP);

    // positive perturb
    tQ = calcUnitTangentVec(thetaQ+h, R_SQ);
    shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
            geodHitPlaneEvent->getPlane(), opts, geodQtmp);
    fyptmp = calcError(geom, geodP, geodQtmp);

    if (order==1) {
        dfdy(1) = (fyptmp-fy0)/h;
    } else {

        geodQtmp.clear();
        tQ = calcUnitTangentVec(thetaQ-h, R_SQ);
        shootGeodesicInDirectionUntilPlaneHit(xQ, tQ,
                geodHitPlaneEvent->getPlane(), opts, geodQtmp);
        fymtmp = calcError(geom, geodP, geodQtmp);

        dfdy(1) = (fyptmp-fymtmp)/(2*h);
    }

    return dfdy;
}

/*static*/Vec2  GeodesicGeometry::
calcError(const ContactGeometry& geom, const Geodesic& geodP, const Geodesic& geodQ) {
    int numP = geodP.getPoints().size();
    int numQ = geodQ.getPoints().size();

    Vec3 Phat = geodP.getPoints()[numP - 1];
    Vec3 Qhat = geodQ.getPoints()[numQ - 1];
    UnitVec3 tPhat(geodP.getTangents()[numP - 1]);
    UnitVec3 tQhat(geodQ.getTangents()[numQ - 1]);

    UnitVec3 nPhat(geom.calcSurfaceNormal((Vector) Phat));
    UnitVec3 nQhat(geom.calcSurfaceNormal((Vector) Qhat));
    UnitVec3 bPhat(nPhat % tPhat);
    UnitVec3 bQhat(nQhat % tQhat);

    //    cout << "bP = " << ~bPhat << endl;
    //    cout << "tQ = " << ~tQhat << endl;

    Vec2 geodErr(~(bPhat - bQhat) * (Phat - Qhat), ~bPhat * tQhat);
    return geodErr;
}

// calculate rotation from tangent plane at R to the surface frame
// the tangent direction is given by the projection of dir onto the tangent plane
/*static*/ Rotation  GeodesicGeometry::
calcTangentBasis(const Vec3& R, const Vec3& dir,
        const ContactGeometry& geom) {
    UnitVec3 n(geom.calcSurfaceNormal((Vector)R));
    Rotation R_GS;
    R_GS.setRotationFromTwoAxes(n, ZAxis, dir, XAxis);
    return R_GS;
}

// calculate a tangent vector based on angle theta measured from the binormal axis
/*static*/ UnitVec3  GeodesicGeometry::
calcUnitTangentVec(const Real& theta, const Rotation& R_GS) {
    UnitVec3 tR_S(std::sin(theta), std::cos(theta), 0);
    return R_GS*tR_S;
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
 *       qdot = Qu
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

        sleep(0.5);
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
 *       qdot = Qu
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




