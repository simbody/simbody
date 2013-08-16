/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Ian Stavness, Michael Sherman                                     *
 * Contributors: Andreas Scholz                                               *
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
#include "simmath/internal/Geodesic.h"
#include "simmath/internal/GeodesicIntegrator.h"
#include "simmath/internal/ContactGeometry.h"

#include "GeodesicEquations.h"

#include <iostream>
#include <cmath>

using namespace SimTK;

// Make sure these compile correctly.
template class GeodesicIntegrator<GeodesicOnImplicitSurface>;
template class GeodesicIntegrator<GeodesicOnParametricSurface>;



//==============================================================================
//                                 GEODESIC
//==============================================================================

void Geodesic::dump(std::ostream& o) const {
    o << "Geodesic: " << getNumPoints() << " points, length=" 
                      << getLength() << "\n";
    bool hasQtoP = !directionalSensitivityQtoP.empty();
    if (!hasQtoP)
        o << "  QtoP Jacobi fields not available\n";
    for (int i=0; i < getNumPoints(); ++i) {
        o << "  Point at s=" << arcLengths[i] << ":\n";
        o << "    p=" << frenetFrames[i].p() << " t=" 
                      << frenetFrames[i].x() << "\n";
        o << "    jrP=" << directionalSensitivityPtoQ[i][0];
        o << "    jtP=" << positionalSensitivityPtoQ[i][0];
        if (hasQtoP) {
          o << " jrQ=" << directionalSensitivityQtoP[i][0];
          o << " jtQ=" << positionalSensitivityQtoP[i][0];
        }
        o << std::endl;
    }
}

//==============================================================================
//                      GEODESIC ON IMPLICIT SURFACE
//==============================================================================

// See class header for documentation.
void GeodesicOnImplicitSurface::
calcDerivs(Real t, const Vec<N>& y, Vec<N>& ydot) const {
    const Vec3& p = getP(y);        // rename state variables
    const Vec3& v = getV(y);
    const Real& jr = getJRot(y);
    const Real& jt = getJTrans(y);

    // Evaluate the surface at p.
    const Vec3  g = geom.calcSurfaceGradient(p);
    const Mat33 H = geom.calcSurfaceHessian(p);
    Real Kg = geom.calcGaussianCurvature(g,H);

    const Real Gdotv = ~v*(H*v);
    const Real L = Gdotv/(~g*g);    // Lagrange multiplier

    // We have qdot = u; that part is easy.
    updQ(ydot) = getU(y);

    // These together are the udots.
    Vec3& a     = updV(ydot);          // d/dt v
    Real& jrdd  = updJRotDot(ydot);    // d/dt jdr
    Real& jtdd  = updJTransDot(ydot);  // d/dt jdt

    a    = -L*g;
    jrdd = -Kg*jr;
    jtdd = -Kg*jt;
}


// See class header for documentation.
void GeodesicOnImplicitSurface::
calcConstraintErrors(Real t, const Vec<N>& y, Vec<NC>& cerr) const {
    const Vec3& p = getP(y);
    const Vec3& v = getV(y);
    // This is the perr() equation that says the point must be on the surface.
    cerr[0] = geom.calcSurfaceValue(p);
    // These are the two verr() equations. The first is the derivative of
    // the above point-on-surface holonomic constraint above. The second is 
    // a nonholonomic velocity constraint restricting the velocity along 
    // the curve to be 1.
    cerr[1] = ~geom.calcSurfaceGradient(p)*v;
    cerr[2] = v.norm() - 1;
}

// Given a state y drive the infinity norm of the position and velocity 
// constraint errors to consTol or below by adjusting y.
bool GeodesicOnImplicitSurface::
projectIfNeeded(Real consTol, Real t, Vec<N>& y) const {
    const int MaxIter = 10;         // should take *far* fewer
    const Real OvershootFac = Real(0.1);  // try to do better than consTol
        
    const Real tryTol = consTol * OvershootFac;
    Vec3& p = updP(y); // aliases for the state variables
    Vec3& v = updV(y);

    // Fix the position constraint first. This is a Newton interation
    // that modifies only the point location to make sure it remains on
    // the surface. No position projection is done if we're already at
    // tryTol, which is a little tighter than the requested consTol.

    // NOTE: (sherm) I don't think this is exactly the right projection.
    // Here we project down the gradient, but the final result won't 
    // be exactly the nearest point on the surface if the gradient changes
    // direction on the way down. For correcting small errors this is
    // probably completely irrelevant since the starting and final gradient
    // directions will be the same.

    Real perr, ptolAchieved;
    int piters=0; 
    while (true) {
        perr = geom.calcSurfaceValue(p);
        ptolAchieved = std::abs(perr); 
        if (ptolAchieved <= tryTol || piters==MaxIter)
            break;

        ++piters;
        // We want a least squares solution dp to ~g*dp=perr which we
        // get using the pseudoinverse: dp=pinv(~g)*perr, where
        // pinv(~g) = g*inv(~g*g).
        const Vec3 g = geom.calcSurfaceGradient(p);
        const Vec3 pinvgt = g/(~g*g);
        const Vec3 dp = pinvgt*perr;

        p -= dp; // updates the state
    }


    // Now the velocities. There are two velocity constraints that have
    // to be satisfied simultaneously. They are (1) the time derivative of 
    // the perr equation which we just solved, and (2) the requirement that 
    // the velocity magnitude be 1. So verr=~[ ~g*v, |v|-1 ]. You might 
    // think these need to be solved simultaneously to find the least 
    // squares dv, but dv can be determined by two orthogonal projections.
    // The allowable velocity vectors form a unit circle whose normal is
    // in the gradient direction. The least squares dv is the shortest
    // vector from the end point of v to that cicle. To find the closest
    // point on the unit circle, first project the vector v onto the 
    // circle's plane by the shortest path (remove the normal component). 
    // Then stretch the result to unit length.
    // First we solve the linear least squares problem ~g*(v+dv0)=0 for
    // dv0, and set v0=v+dv0. Then set vfinal = v0/|v0|, giving dv=vfinal-v.

    // We're going to project velocities unconditionally because we
    // would have to evaluate the constraint anyway to see if it is
    // violated and that is most of the computation we need to fix it.

    const Vec3 g = geom.calcSurfaceGradient(p);
    const Vec3 pinvgt = g/(~g*g);
    const Real perrdot = ~g*v;

    const Vec3 dv0 = pinvgt*perrdot;
    const Vec3 v0 = v - dv0;    // fix direction
    v = v0/v0.norm();           // fix length; updates state

    const bool success = (ptolAchieved <= consTol);
    return success;
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
    const Real Kg = geom.calcGaussianCurvature(GT,H);
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

/*
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
*/
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

/*
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
*/
    results.setExitStatus(ProjectResults::Succeeded);
//    std::cout << "norm(u) = " << s.getU(subsysIndex).norm() << std::endl;
}






