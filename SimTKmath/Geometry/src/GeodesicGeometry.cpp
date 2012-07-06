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
            const Vec3& xP, const Vec3& xQ, const Vec3& tPprev, const Vec3& tQprev) :
            Differentiator::JacobianFunction(nf, ny),
                    gg(geodgeom),
                    P(xP), Q(xQ),
                    R_SP(calcTangentBasis(P, tPprev, gg.geom)),
                    R_SQ(calcTangentBasis(Q, tQprev, gg.geom)) { }

    // x = ~[thetaP, thetaQ]
    int f(const Vector& x, Vector& fx) const  {
        Vec3 tP = calcTangentVec(x[0], R_SP);
        Vec3 tQ = calcTangentVec(x[1], R_SQ);

        Vec2 geodErr = GeodesicGeometry::calcGeodError(gg, P, Q, tP, tQ);

        // error between geodesic end points at plane
        fx[0] = geodErr[0];
        fx[1] = geodErr[1];
        return 0;
    }

    // calculate rotation from tangent plane at R to the surface frame
    // the tangent direction is given by the projection of dir onto the tangent plane
    static Rotation calcTangentBasis(const Vec3& R, const Vec3& dir,
            const ContactGeometry& geom) {
        UnitVec3 n(geom.calcSurfaceGradient((Vector)R));
        Rotation R_GS;
        R_GS.setRotationFromTwoAxes(n, ZAxis, dir, XAxis);
        return R_GS;
    }

    // calculate a tangent vector based on angle theta measured from the binormal axis
    static Vec3 calcTangentVec(const Real& theta, const Rotation& R_GS) {
        Vec3 tR_S(std::sin(theta), std::cos(theta), 0);
        return R_GS*tR_S;
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



// Compute a geodesic curve of the given length, starting at the given point and
// in the given direction.
// XXX what to do if tP is not in the tangent plane at P -- project it?
void GeodesicGeometry::
calcGeodesicInDirectionUntilLengthReached(const Vec3& P, const Vec3& tP,
        const Real& terminatingLength, const GeodesicOptions& options,
        Geodesic& geod) const {

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
    integ.setFinalTime(terminatingLength);
    integ.setReturnEveryInternalStep(true); // save geodesic knot points
    integ.initialize(sysState);
    const State& state = integ.getState(); // integrator state

    // Simulate it, and record geodesic knot points after each step
    Integrator::SuccessfulStepStatus status;
    int stepcnt = 0;
    while (true) {
        status = integ.stepTo(terminatingLength);
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
calcGeodesicInDirectionUntilPlaneHit(const Vec3& P, const Vec3& tP,
        const Plane& terminatingPlane, const GeodesicOptions& options,
        Geodesic& geod) const {

    // Initialize state
    geodHitPlaneEvent->setPlane(terminatingPlane);
    State sysState = ptOnSurfSys.getDefaultState();
    sysState.setTime(startTime);
    Vector& q = sysState.updQ();
    Vector& u = sysState.updU();
    q[0] = P[0]; q[1] = P[1]; q[2] = P[2];
    u[0] = tP[0]; u[1] = tP[1]; u[2] = tP[2];

    // Setup integrator to integrate until terminatingLength
    RungeKutta3Integrator integ(ptOnSurfSys);
    integ.setAccuracy(integratorAccuracy);
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


// Utility method to find geodesic between P and Q
// with starting directions tPhint and tQhint
// XXX tangent basis should be formed from previous geodesic
void GeodesicGeometry::calcGeodesic(const Vec3& xP, const Vec3& xQ,
        const Vec3& tPhint, const Vec3& tQhint, Geodesic& geod) const {

    // calculate plane bisecting P and Q, and use as termination condition for integrator
    UnitVec3 normal(xQ - xP);
    Real offset = (~(xP+xQ)*normal)/2 ;
    geodHitPlaneEvent->setPlane(Plane(normal, offset));

    Vector x(2), dx(2), Fx(2), xold(2);
    Matrix J(2,2);

//    Vec2 x, dx, fx;
//    Mat22 J;

    // initial conditions
    x[0] = Pi/2;
    x[1] = Pi/2;

    Real tol = 1e-12, fold, lam = 1;
    splitGeodErr = new SplitGeodesicError(2, 2, *const_cast<GeodesicGeometry*>(this),
            xP, xQ, tPhint, tQhint);
    splitGeodErr->setEstimatedAccuracy(tol);
    Differentiator diff( *const_cast<SplitGeodesicError*>(splitGeodErr));

    splitGeodErr->f(x, Fx);

    Real f = 0.5*~Fx*Fx;
    int maxIter = 40;
    for (int i = 0; i < maxIter; ++i) {
        if (std::sqrt(f) < tol) {
            std::cout << "geodesic converged in " << i << " iterations" << std::endl;
            break;
        }
//        std::cout << "err = " << Fx << ", x = " << x << std::endl;

        diff.calcJacobian(x, Fx, J, Differentiator::ForwardDifference);
        fold = f;
        xold = x;
//        cout << "J = " << J << endl;
        dx = J.invert()*Fx;

        // backtracking
        lam = 1;
        while (true) {
            x = xold - lam*dx;
            splitGeodErr->f(x, Fx);
            f = 0.5*~Fx*Fx;
            if (f > fold) {
                lam = lam / 2;
            } else {
                break;
            }
        }
    }
//    std::cout << "err = " << Fx << std::endl;

    mergeGeodesics(geodP, geodQ, geod);
}


// Calculate the "geodesic error"
/*static*/Vec2 GeodesicGeometry::calcGeodError(const GeodesicGeometry& gg,
        const Vec3& P, const Vec3& Q, const Vec3& tP, const Vec3& tQ) {

    gg.geodP.clear();
    gg.geodQ.clear();

    GeodesicOptions opts;
    gg.calcGeodesicInDirectionUntilPlaneHit(P, tP,
            gg.geodHitPlaneEvent->getPlane(), opts, gg.geodP);
    gg.calcGeodesicInDirectionUntilPlaneHit(Q, tQ,
            gg.geodHitPlaneEvent->getPlane(), opts, gg.geodQ);

    int numP = gg.geodP.getPoints().size();
    int numQ = gg.geodQ.getPoints().size();

    Vec3 Phat = gg.geodP.getPoints()[numP - 1];
    Vec3 Qhat = gg.geodQ.getPoints()[numQ - 1];
    Vec3 tPhat = gg.geodP.getTangents()[numP - 1];
    Vec3 tQhat = gg.geodQ.getTangents()[numQ - 1];

    Vec3 nPhat = gg.geom.calcSurfaceGradient((Vector) P);
    nPhat = nPhat.normalize();
    Vec3 bPhat;
    bPhat = cross(nPhat, tPhat);

    Vec2 geodErr(~bPhat * (Phat - Qhat), ~bPhat * tQhat);
    return geodErr;
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
    Vec3 GT = geom.calcSurfaceGradient(q);
    Mat33 H = geom.calcSurfaceHessian(q);
    Real Gdotu = ~v*(H*v);
    Real L = (Gdotu + beta*~GT*v + alpha*g)/(~GT*GT);
    a = GT*-L;

    udot[0] = a[0]; udot[1] = a[1]; udot[2] = a[2];
    s.updQDotDot() = udot;

    System::Guts::realizeAccelerationImpl(s);
    return 0;
}








