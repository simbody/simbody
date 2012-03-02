#ifndef SimTK_SIMMATH_PENDULUMSYSTEM_H_
#define SimTK_SIMMATH_PENDULUMSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

/**
 * This is a simple System that is used by various test cases.
 */

#include "SimTKcommon.h"
#include "SimTKcommon/internal/SystemGuts.h"

using namespace SimTK;

class PendulumSystem;
class PendulumSystemGuts: public System::Guts {
    friend class PendulumSystem;

    // TOPOLOGY STATE
    SubsystemIndex subsysIndex;

    // TOPOLOGY CACHE
    mutable DiscreteVariableIndex massIndex, lengthIndex, gravityIndex;
    mutable QIndex q0;
    mutable UIndex u0;
    mutable QErrIndex qerr0;
    mutable UErrIndex uerr0;
    mutable UDotErrIndex udoterr0;
    mutable EventTriggerByStageIndex event0;
    mutable CacheEntryIndex mgForceIndex; // a cache entry m*g calculated at Dynamics stage
public:
    PendulumSystemGuts() : Guts() {
        // Index types set themselves invalid on construction.
    }

    const PendulumSystem& getPendulumSystem() const {
        return reinterpret_cast<const PendulumSystem&>(getSystem());
    }
    
    SubsystemIndex getSubsysIndex() const {
        return subsysIndex;
    }

    /*virtual*/PendulumSystemGuts* cloneImpl() const {return new PendulumSystemGuts(*this);}

        /////////////////////////////////////////////////////////
        // Implementation of continuous DynamicSystem virtuals //
        /////////////////////////////////////////////////////////

    /*virtual*/int realizeTopologyImpl(State&) const;
    /*virtual*/int realizeModelImpl(State&) const;
    /*virtual*/int realizeInstanceImpl(const State&) const;
    /*virtual*/int realizePositionImpl(const State&) const;
    /*virtual*/int realizeVelocityImpl(const State&) const;
    /*virtual*/int realizeDynamicsImpl(const State&) const;
    /*virtual*/int realizeAccelerationImpl(const State&) const;

    // qdot==u here so these are just copies
    /*virtual*/void multiplyByNImpl(const State& state, const Vector& u, 
                                 Vector& dq) const {dq=u;}
    /*virtual*/void multiplyByNTransposeImpl(const State& state, const Vector& fq, 
                                          Vector& fu) const {fu=fq;}
    /*virtual*/void multiplyByNPInvImpl(const State& state, const Vector& dq, 
                                     Vector& u) const {u=dq;}
    /*virtual*/void multiplyByNPInvTransposeImpl(const State& state, const Vector& fu, 
                                              Vector& fq) const {fq=fu;}

    // No prescribed motion.
    /*virtual*/bool prescribeQImpl(State&) const {return false;}
    /*virtual*/bool prescribeUImpl(State&) const {return false;}

    /*virtual*/void projectQImpl(State&, Vector& qErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;
    /*virtual*/void projectUImpl(State&, Vector& uErrEst, 
             const ProjectOptions& options, ProjectResults& results) const;

};

class PendulumSystem: public System {
public:
    PendulumSystem() : System()
    { 
        adoptSystemGuts(new PendulumSystemGuts());
        DefaultSystemSubsystem defsub(*this);
        updGuts().subsysIndex = defsub.getMySubsystemIndex();

        setHasTimeAdvancedEvents(false);
    }

    const PendulumSystemGuts& getGuts() const {
        return dynamic_cast<const PendulumSystemGuts&>(getSystemGuts());
    }

    PendulumSystemGuts& updGuts() {
        return dynamic_cast<PendulumSystemGuts&>(updSystemGuts());
    }

    // Instance variables are written to our defaultState.
    void setDefaultMass(Real mass) {
        const PendulumSystemGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.massIndex) = Value<Real>(mass);
    }

    void setDefaultLength(Real length) {
        const PendulumSystemGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.lengthIndex) = Value<Real>(length);
    }

    void setDefaultGravity(Real gravity) {
        const PendulumSystemGuts& guts = getGuts();
        updDefaultState().updDiscreteVariable(guts.subsysIndex, guts.gravityIndex) = Value<Real>(gravity);
    }

    void setDefaultTimeAndState(Real t, const Vector& q, const Vector& u) {
        const PendulumSystemGuts& guts = getGuts();
        updDefaultState().updU(guts.subsysIndex) = u;
        updDefaultState().updQ(guts.subsysIndex) = q;
        updDefaultState().updTime() = t;
    }

    Real getMass(const State& s) const {
        const PendulumSystemGuts& guts = getGuts();
        const AbstractValue& m = s.getDiscreteVariable(guts.subsysIndex, guts.massIndex);
        return Value<Real>::downcast(m).get();
    }
    Real getDefaultMass() const {return getMass(getDefaultState());}

    Real getLength(const State& s) const {
        const PendulumSystemGuts& guts = getGuts();
        const AbstractValue& d = s.getDiscreteVariable(guts.subsysIndex, guts.lengthIndex);
        return Value<Real>::downcast(d).get();
    }
    Real getDefaultLength() const {return getLength(getDefaultState());}

    Real getGravity(const State& s) const {
        const PendulumSystemGuts& guts = getGuts();
        const AbstractValue& g = s.getDiscreteVariable(guts.subsysIndex, guts.gravityIndex);
        return Value<Real>::downcast(g).get();
    }
    Real getDefaultGravity() const {return getGravity(getDefaultState());}
};

/*
 * This system is a 2d pendulum swinging in gravity. It is modeled as
 * a point mass free in the plane, plus a distance constraint to model
 * the rod.
 *
 *    y       | g               O
 *    ^       v                  \  d
 *    |                           \
 *    |                            * m
 *     ------> x
 *
 * Gravity acts in the y direction, the rod is length d, mass m, pivot
 * location is the ground origin (0,0).
 *
 * The DAE for a generic multibody system is:
 *       qdot = Qu
 *       M udot = f - ~A lambda
 *       A udot = b
 *       perr(t,q) = 0
 *       verr(t,q,u) = 0
 *
 * Let   r^2 = x^2  + y^2
 *       v^2 = x'^2 + y'^2
 * We will express the "rod length=d" constraint as 
 *       (r^2 - d^2)/2 = 0    (perr)
 *           xx' + yy' = 0    (verr)
 *         xx'' + yy'' = -v^2 (aerr)
 *
 * So the matrix A = d perr/dq = [x y] and b = -v^2, and the
 * equations of motion are:
 *     [ m 0 x ] [ x'' ]   [  0  ]
 *     [ 0 m y ] [ y'' ] = [ -mg ]
 *     [ x y 0 ] [ L   ]   [-v^2 ]
 * where L (the Lagrange multiplier) is proportional to
 * the rod tension. You can solve this to get
 *     L   = (m*v^2 - mg*y)/(r^2)
 *     x'' = - x*L/m
 *     y'' = - y*L/m - g
 *               
 */ 
int PendulumSystemGuts::realizeTopologyImpl(State& s) const {
    // Instance variables mass, length, gravity
    massIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(1));
    lengthIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(1));
    gravityIndex = s.allocateDiscreteVariable(subsysIndex, Stage::Instance,
                                                      new Value<Real>(13.7503716373294544));
    const Vector init(2, Real(0));
    q0 = s.allocateQ(subsysIndex, init);
    u0 = s.allocateU(subsysIndex, init);

    mgForceIndex = s.allocateCacheEntry(subsysIndex, Stage::Dynamics,
                                             new Value<Real>());
    System::Guts::realizeTopologyImpl(s);
    return 0;
}
int PendulumSystemGuts::realizeModelImpl(State& s) const {
    System::Guts::realizeModelImpl(s);
    return 0;
}
int PendulumSystemGuts::realizeInstanceImpl(const State& s) const {
    qerr0 = s.allocateQErr(subsysIndex, 1);
    uerr0 = s.allocateUErr(subsysIndex, 1);
    udoterr0 = s.allocateUDotErr(subsysIndex, 1); // and multiplier
//    event0 = s.allocateEvent(subsysIndex, Stage::Position, 3);
    System::Guts::realizeInstanceImpl(s);
    return 0;
}
int PendulumSystemGuts::realizePositionImpl(const State& s) const {
    const Real    d = getPendulumSystem().getLength(s);
    const Vector& q = s.getQ(subsysIndex);
    // This is the perr() equation.
    s.updQErr(subsysIndex)[0] = (q[0]*q[0] + q[1]*q[1] - d*d)/2;
    
//    s.updEventsByStage(subsysIndex, Stage::Position)[0] = 100*q[0]-q[1];
//
//    s.updEventsByStage(subsysIndex, Stage::Position)[1] = 
//        s.getTime() > 1.49552 && s.getTime() < 12.28937;
//
//    s.updEventsByStage(subsysIndex, Stage::Position)[2] = s.getTime()-1.495508;
    System::Guts::realizePositionImpl(s);
    return 0;
}

int PendulumSystemGuts::realizeVelocityImpl(const State& s) const {
    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       qdot = s.updQDot(subsysIndex);

    qdot[0] = u[0]; // qdot=u
    qdot[1] = u[1];

    // This is the verr() equation.
    s.updUErr(subsysIndex)[0]  = q[0]*u[0] + q[1]*u[1];
    System::Guts::realizeVelocityImpl(s);
    return 0;
}

int PendulumSystemGuts::realizeDynamicsImpl(const State& s) const {
    const Real m  = getPendulumSystem().getMass(s);
    const Real g  = getPendulumSystem().getGravity(s);

    Real& mg = Value<Real>::downcast(s.updCacheEntry(subsysIndex, mgForceIndex)).upd();
    // Calculate the force due to gravity.
    mg = m*g;
    System::Guts::realizeDynamicsImpl(s);
    return 0;
}

int PendulumSystemGuts::realizeAccelerationImpl(const State& s) const {
    const Real m  = getPendulumSystem().getMass(s);
    const Real g  = getPendulumSystem().getGravity(s);
    // we're pretending we couldn't calculate this here!
    const Real mg = Value<Real>::downcast
                       (s.updCacheEntry(subsysIndex, mgForceIndex)).get();

    const Vector& q    = s.getQ(subsysIndex);
    const Vector& u    = s.getU(subsysIndex);
    Vector&       udot = s.updUDot(subsysIndex);

    const Real r2 = q[0]*q[0] + q[1]*q[1];
    const Real v2 = u[0]*u[0] + u[1]*u[1];
    const Real L  = (m*v2 - mg*q[1])/r2;
    udot[0] = - q[0]*L/m;
    udot[1] = - q[1]*L/m - g;
    s.updQDotDot() = udot;
    s.updMultipliers(subsysIndex)[0] = L;
    s.updUDotErr(subsysIndex)[0] = q[0]*udot[0] + q[1]*udot[1] + v2;
    System::Guts::realizeAccelerationImpl(s);
    return 0;
}

/*
 * Here we want to remove any constraint errors from the current state,
 * and project out any component of the integrator's error estimate
 * perpendicular to the constraint manifold. We will do this sequentially
 * rather than handling position and velocity simultaneously.
 *
 * For this system we have P = d perr/dq = V = d verr/du = [x y].
 * Weighted, we have PW=tp*[x/wx y/wy] VW=tv*[x/wxd y/wyd]. 
 * With pinv(A)=~A*(A*~A)^-1, we have:
 *
 *    pinv(P)  = ~[            x             y] /  (    x ^2+     y ^2)
 *    pinv(PW) = ~(1/tp)*[(wx *wy ^2)*x (wx ^2*wy) *y] / ((wy *x)^2+(wx *y)^2)
 *    pinv(VW) = ~(1/tv)*[(wxd*wyd^2)*x (wxd^2*wyd)*y] / ((wyd*x)^2+(wxd*y)^2)
 *      (the latter assuming x,y already projected on position manifold)
 *
 * We want to solve
 *    |perr(q0 - dq)|_TRMS <= accuracy, such that dq=min_WLS(dq)
 *    PW(q0) dq = Tp * perr(q0); q = q0-dq
 * Then
 *    |verr(q,u0 - du)|_TRMS <= accuracy, du=min_WLS(du)
 *    VW(q) du = Tv * verr(q,u0); u = u0-du
 *
 *
 * To remove the corresponding error estimates:
 *    PW(q) qperp = PW(q) qerrest; qerrest -= qperp
 *    VW(q) uperp = VW(q) uerrest; uerrest -= uperp
 * 
 *
 */
static Real wrms(const Vector& y, const Vector& w) {
    Real sumsq = 0;
    for (int i=0; i<y.size(); ++i)
        sumsq += square(y[i]*w[i]);
    return std::sqrt(sumsq/y.size());
}

//int PendulumSystemGuts::projectImpl(State& s, Real consAccuracy,
//                                const Vector& yweights, const Vector& ctols,
//                                Vector& yerrest, System::ProjectOptions opts) const // yerrest is in/out
//{
//    const Vec2& wq = Vec2::getAs(&yweights[0]);
//    const Vec2& wu = Vec2::getAs(&yweights[2]);
//    const Real& tp = ctols[0]; // inverse tolerances 1/ti
//    const Real& tv = ctols[1];
//
//    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
//    const Vec2& u = Vec2::getAs(&s.getU(subsysIndex)[0]);
//    Real& ep = s.updQErr(subsysIndex)[0];
//    Real& ev = s.updUErr(subsysIndex)[0];
//
//    //cout << "BEFORE wperr=" << tp*ep << endl;
//
//    Real wqchg;
//    if (opts.hasAnyPositionOptions()) {
//        do {
//            // Position projection
//            Real r2 = ~q*q; // x^2+y^2
//            Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
//            Row2 P(~q), PW(tp*q[0]/wq[0], tp*q[1]/wq[1]);
//            Vec2 Pinv(q/r2);
//            Vec2 PWinv = Vec2(square(wq[1])*wq[0]*q[0], 
//                              square(wq[0])*wq[1]*q[1]) / (tp*wqr2);
//            Vec2 dq  = Pinv*(ep);      //cout << "dq=" << dq << endl;
//            Vec2 wdq = PWinv*(tp*ep);  //cout << "wdq=" << wdq << endl;
//    
//            wqchg = std::sqrt(wdq.normSqr()/q.size()); // wrms norm
//    
//            s.updQ(subsysIndex)[0] -= wdq[0]/wq[0]; 
//            s.updQ(subsysIndex)[1] -= wdq[1]/wq[1]; 
//            realize(s, Stage::Position); // recalc QErr (ep)
//    
//            //cout << "AFTER q-=wdq/W wperr=" << tp*ep << " wqchg=" << wqchg << endl;
//        } while (std::abs(tp*ep) > consAccuracy && wqchg >= 0.01*consAccuracy);
//    }
//
//    // Do velocity projection at new values of q
//    Real r2 = ~q*q; // x^2+y^2
//    Real wur2 = square(wu[1]*q[0]) + square(wu[0]*q[1]);
//    Row2 V(~q), VW(tv*q[0]/wu[0], tv*q[1]/wu[1]);
//    Vec2 Vinv(q/r2);
//    Vec2 VWinv = Vec2(square(wu[1])*wu[0]*q[0], 
//                      square(wu[0])*wu[1]*q[1]) / (tv*wur2);
//    realize(s, Stage::Velocity); // calculate UErr (ev)
//
//    //cout << "BEFORE wverr=" << tv*ev << endl;
//    Vec2 du  = Vinv*(ev);      //cout << "du=" << du << endl;
//    Vec2 wdu = VWinv*(tv*ev);  //cout << "wdu=" << wdu << endl;
//
//    s.updU(subsysIndex)[0] -= wdu[0]/wu[0]; 
//    s.updU(subsysIndex)[1] -= wdu[1]/wu[1];
//
//    realize(s, Stage::Velocity); // recalc UErr
//    //cout << "AFTER u-=wdu wverr=" << tv*ev << endl;
//
//    // Now do error estimates.
//
//
//    if (yerrest.size()) {
//        Vec2& eq = Vec2::updAs(&yerrest[0]);
//        Vec2& eu = Vec2::updAs(&yerrest[2]);
//
//        // Recalc PW, PWInv:
//        const Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
//        const Row2 PW = Row2(tp*q[0]/wq[0], tp*q[1]/wq[1]);
//        const Vec2 PWinv = Vec2(wq[0]*square(wq[1])*q[0], 
//                                square(wq[0])*wq[1]*q[1]) / (tp*wqr2);
//
//        Vec2 qperp = PWinv*(PW*eq);
//        Vec2 uperp = VWinv*(VW*eu);
//
//        //cout << "ERREST before=" << yerrest 
//        //     << " wrms=" << wrms(yerrest,yweights) << endl;
//        //cout << "PW*eq=" << PW*eq << " VW*eu=" << VW*eu << endl;
//        eq -= qperp; eu -= uperp;
//
//        //cout << "ERREST after=" << yerrest 
//        //     << " wrms=" << wrms(yerrest,yweights) << endl;
//        //cout << "PW*eq=" << PW*eq << " VW*eu=" << VW*eu << endl;
//    }
//
//    return 0;
//}

// qerrest is in/out
void PendulumSystemGuts::projectQImpl(State& s, Vector& qerrest, 
                                const ProjectOptions& opts,
                                ProjectResults& results) const 
                                
{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getQErrWeights(subsysIndex);
    // Since qdot=u here we can use uweights directly as qweights.
    const Vec2& wq = Vec2::getAs(&uweights[0]);
    const Real& tp = ctols[0]; // inverse tolerances 1/ti

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    Real& ep = s.updQErr(subsysIndex)[0];

    //cout << "BEFORE wperr=" << tp*ep << endl;

    Real wqchg;
    do {
        // Position projection
        Real r2 = ~q*q; // x^2+y^2
        Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
        Row2 P(~q), PW(tp*q[0]/wq[0], tp*q[1]/wq[1]);
        Vec2 Pinv(q/r2);
        Vec2 PWinv = Vec2(square(wq[1])*wq[0]*q[0], 
                            square(wq[0])*wq[1]*q[1]) / (tp*wqr2);
        Vec2 dq  = Pinv*(ep);      //cout << "dq=" << dq << endl;
        Vec2 wdq = PWinv*(tp*ep);  //cout << "wdq=" << wdq << endl;
    
        wqchg = std::sqrt(wdq.normSqr()/q.size()); // wrms norm
    
        s.updQ(subsysIndex)[0] -= wdq[0]/wq[0]; 
        s.updQ(subsysIndex)[1] -= wdq[1]/wq[1]; 
        realize(s, Stage::Position); // recalc QErr (ep)
    
        //cout << "AFTER q-=wdq/W wperr=" << tp*ep << " wqchg=" << wqchg << endl;
    } while (std::abs(tp*ep) > consAccuracy && wqchg >= 0.01*consAccuracy);

    // Now do error estimates.

    if (qerrest.size()) {
        Vec2& eq = Vec2::updAs(&qerrest[0]);

        // Recalc PW, PWInv:
        const Real wqr2 = square(wq[1]*q[0]) + square(wq[0]*q[1]);
        const Row2 PW = Row2(tp*q[0]/wq[0], tp*q[1]/wq[1]);
        const Vec2 PWinv = Vec2(wq[0]*square(wq[1])*q[0], 
                                square(wq[0])*wq[1]*q[1]) / (tp*wqr2);

        Vec2 qperp = PWinv*(PW*eq);

        //cout << "ERREST before=" << yerrest 
        //     << " wrms=" << wrms(qerrest,qweights) << endl;
        //cout << "PW*eq=" << PW*eq << endl;
        eq -= qperp;

        //cout << "ERREST after=" << yerrest 
        //     << " wrms=" << wrms(qerrest,qweights) << endl;
        //cout << "PW*eq=" << PW*eq << endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

void PendulumSystemGuts::projectUImpl(State& s, Vector& uerrest, 
             const ProjectOptions& opts, ProjectResults& results) const
{
    const Real consAccuracy = opts.getRequiredAccuracy();
    const Vector& uweights = s.getUWeights(subsysIndex);
    const Vector& ctols = s.getUErrWeights(subsysIndex);

    const Vec2& wu = Vec2::getAs(&uweights[0]);
    const Real& tv = ctols[0];

    const Vec2& q = Vec2::getAs(&s.getQ(subsysIndex)[0]); // set up aliases
    const Vec2& u = Vec2::getAs(&s.getU(subsysIndex)[0]);
    Real& ev = s.updUErr(subsysIndex)[0];

    //cout << "BEFORE wperr=" << tp*ep << endl;

    // Do velocity projection at current values of q, which should have
    // been projected already.
    Real r2 = ~q*q; // x^2+y^2
    Real wur2 = square(wu[1]*q[0]) + square(wu[0]*q[1]);
    Row2 V(~q), VW(tv*q[0]/wu[0], tv*q[1]/wu[1]);
    Vec2 Vinv(q/r2);
    Vec2 VWinv = Vec2(square(wu[1])*wu[0]*q[0], 
                      square(wu[0])*wu[1]*q[1]) / (tv*wur2);
    realize(s, Stage::Velocity); // calculate UErr (ev)

    //cout << "BEFORE wverr=" << tv*ev << endl;
    Vec2 du  = Vinv*(ev);      //cout << "du=" << du << endl;
    Vec2 wdu = VWinv*(tv*ev);  //cout << "wdu=" << wdu << endl;

    s.updU(subsysIndex)[0] -= wdu[0]/wu[0]; 
    s.updU(subsysIndex)[1] -= wdu[1]/wu[1];

    realize(s, Stage::Velocity); // recalc UErr
    //cout << "AFTER u-=wdu wverr=" << tv*ev << endl;

    // Now do error estimates.


    if (uerrest.size()) {
        Vec2& eu = Vec2::updAs(&uerrest[0]);
        Vec2 uperp = VWinv*(VW*eu);

        //cout << "ERREST before=" << uerrest 
        //     << " wrms=" << wrms(uerrest,uweights) << endl;
        //cout << " VW*eu=" << VW*eu << endl;
        eu -= uperp;

        //cout << "ERREST after=" << yerrest 
        //     << " wrms=" << wrms(uerrest,uweights) << endl;
        //cout << " VW*eu=" << VW*eu << endl;
    }

    results.setExitStatus(ProjectResults::Succeeded);
}

#endif /*SimTK_SIMMATH_PENDULUMSYSTEM_H_*/
