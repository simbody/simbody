/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Tests of Simbody integration methods.
 */

#include "SimTKcommon.h"
#include "Simbody.h"

#include "simbody/internal/NumericalMethods.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

class SinCos : public MatterSubsystem {
public:
    SinCos(const Real& omega) : w(omega), y(2), yd(2), weights(2), units(2), t(CNT<Real>::getNaN()), 
        solnAcc(1e-3), consAcc(1e-3), timeScale(1e-3) 
    {
        units[0] = units[1] = 1.;
        relativeOK[0] = relativeOK[1] = true;
        ydValid = false;
    }

    int size() const {return 2;}
    void setAccuracy(const Real& s, const Real& c) {
        solnAcc=s; consAcc=c; calcWeights(solnAcc);
    }
    void setState(const Real& time, const Vector& state) {
        t = time;
        y = state;
        yd.setToNaN(); ydValid = false;
        calcWeights(solnAcc);
    }
    bool realize() const {yd[0]=w*y[1]; yd[1]=-w*y[0]; ydValid=true; return true;}
    bool realizeAndProject(bool& chg, bool) { realize(); chg=false; return true;}
    Real getTimescale() const {return timeScale;}
    const Real& getT() const {return t;}
    const Vector& getY() const {return y;}
    const Vector& getWeights() const {return weights;}
    const Vector& getYDot() const {assert(ydValid);return yd;}
    const Vector& getPositionError() const {assert(ydValid);return perr;}
    const Vector& getVelocityError() const {assert(ydValid);return verr;}
    const Vector& getAccelerationError() const {assert(ydValid);return aerr;}


private:
    void calcWeights(Real& acc) {
        for (int i=0; i<2; ++i) {
            Real okErr = acc*units[i];
            if (relativeOK[i]) okErr += acc*fabs(y[i]);
            weights[i] = 1./okErr;
        }
    }
    Real w;
    Vector units;
    bool relativeOK[2];
    mutable bool ydValid;
    mutable Vector yd, perr, verr, aerr;
    Real t; Vector y;
    Vector weights;
    Real solnAcc, consAcc, timeScale;
};

/**
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
 *     [ 0 m y ] [ y'' ] = [  mg ]
 *     [ x y 0 ] [ L   ]   [-v^2 ]
 * where L (the Lagrange multiplier) is proportional to
 * the rod tension. You can solve this to get
 *     L   = (m/r^2)*(y*g + v^2)
 *     x'' =   - x*L/m
 *     y'' = g - y*L/m
 *               
 */  
class PointMass2dPendulum : public MatterSubsystem {
    static const int NStates = 4;
    static const int NConstraints = 1;
public:
    PointMass2dPendulum(const Real& mass, const Real& length, const Real& gravity) 
        : m(mass), d(length), g(gravity)
    {
        y.resize(NStates); yd.resize(NStates);
        perr.resize(NConstraints); verr.resize(NConstraints); aerr.resize(NConstraints);
        weights.resize(NStates); yUnits.resize(NStates); 
        pUnits.resize(NConstraints); vUnits.resize(NConstraints); aUnits.resize(NConstraints);
        t = CNT<Real>::getNaN();
        solnAcc = consAcc =  1e-3;

        // Make the time scale 1/10 of a period.
        const Real pi = acos(-1.);
        const Real period = 2*pi*sqrt(fabs(length/gravity));
        timeScale = 0.1*period;
        pUnits = 1; vUnits = 1*timeScale; aUnits = 1*timeScale*timeScale/2;
        yUnits[0] = yUnits[1] = 1;
        yUnits[2] = yUnits[3] = 1*timeScale; 
        relativeOK[0] = relativeOK[1] = false;
        relativeOK[2] = relativeOK[3] = true;
        ydValid = false;
    }

    int size() const {return 4;}
    void setAccuracy(const Real& s, const Real& c) {
        solnAcc=s; consAcc=c; calcWeights(solnAcc);
    }
    void setState(const Real& time, const Vector& state) {
        t = time;
        y = state;
        yd.setToNaN(); ydValid = false;
        calcWeights(solnAcc);
    }
    bool realize() const {
        calcPerr(); 
        calcVerr();
        calcYDot();
        calcAerr();
        return aerrNorm <= consAcc;
    }
    bool realizeAndProject(bool& chg, bool force) { 
        chg=false;
        calcPerr();
        if (force || perrNorm > consAcc) {
            if (!projectPositions(consAcc)) return false;
            chg=true;
        } 
        // recalculate at new positions
        calcVerr();
        if (force || verrNorm > consAcc) {
            if (!projectVelocities(consAcc)) return false;
            chg=true;
        }
        calcYDot();
        return calcAerr() <= consAcc;
    }
    Real getTimescale() const {return timeScale;}
    const Real& getT() const {return t;}
    const Vector& getY() const {return y;}
    const Vector& getWeights() const {return weights;}
    const Vector& getYDot() const {assert(ydValid);return yd;}
    const Vector& getPositionError() const {assert(ydValid);return perr;}
    const Vector& getVelocityError() const {assert(ydValid);return verr;}
    const Vector& getAccelerationError() const {assert(ydValid);return aerr;}
    const Real& getPositionErrorNorm() const {assert(ydValid);return perrNorm;}
    const Real& getVelocityErrorNorm() const {assert(ydValid);return verrNorm;}
    const Real& getAccelerationErrorNorm() const {assert(ydValid);return aerrNorm;}

    const Real calcEnergy() {
        assert(ydValid);
        Real pe = (-y[1])*m*g;
        Real ke = m*(y[2]*y[2]+y[3]*y[3])/2.;
        return pe+ke;
    }

private:
    void calcYDot() const {
        yd[0]=y[2]; yd[1]=y[3]; // qdot=u
        const Real r2 = y[0]*y[0] + y[1]*y[1]; // x^2+y^2
        const Real v2 = y[2]*y[2] + y[3]*y[3]; // x'^2 + y'^2
        const Real lambda = (m/r2)*(y[1]*g + v2);
        yd[2] =   - y[0]*lambda/m;
        yd[3] = g - y[1]*lambda/m;
        ydValid=true; 
    }

    // (r^2 - d^2)/2 = 0    (perr)
    //     xx' + yy' = 0    (verr)
    //   xx'' + yy'' = -v^2 (aerr)
    Real calcPerr() const {
        const Real r2 = y[0]*y[0] + y[1]*y[1];
        perr[0] = (r2 - d*d) / 2;
        return perrNorm = fabs(perr[0]/pUnits[0]);
    }
    // Solve perr(q)=0
    bool projectPositions(const Real& tol) {
        // Want to solve [x y][dx] = -perr(x,y)
        //                    [dy]
        // for least squares change [dx dy].
        // Use pseudo inverse:
        //   ~[dx dy] = -pinv(A) * perr(x,y)
        // A=[x y], pseudoInv(A)= ~A * inv(A*~A)
        //
        // = [x] 1/(x^2+y^2)
        //   [y] 

        //printf("POSITION REPAIR ... ERR=%10.5g", perrNorm/tol);
        int iterationsLeft = 5;
        Real oldPerrNorm;
        do {
            oldPerrNorm = perrNorm;
            const Real r2 = y[0]*y[0] + y[1]*y[1]; // A*~A
            y[0] -= (y[0]/r2)*perr[0];
            y[1] -= (y[1]/r2)*perr[0];
            calcPerr();
            //printf(" --> %10.5g", perrNorm/tol);
            --iterationsLeft;
        } while (iterationsLeft && perrNorm > 0.1*tol && perrNorm < oldPerrNorm);

        //printf("\n");

        return perrNorm <= tol;
    }
    Real calcVerr() const {
        verr[0] = y[0]*y[2] + y[1]*y[3];
        return verrNorm = fabs(verr[0])/vUnits[0];
    }
    // Solve verr(q,u)=0
    bool projectVelocities(const Real& tol) {

        // Want to solve [x y][dxd] = -verr(x,y,xd,yd)
        //                    [dyd]
        // for least squares velocity change [dxd dyd].
        // Use pseudo inverse as above except no need to
        // iterate since verr is linear in xd, yd.

        //printf("velocity repair ... ERR=%10.5g", verrNorm/tol);

        const Real r2 = y[0]*y[0] + y[1]*y[1]; // A*~A
        y[2] -= (y[0]/r2)*verr[0];
        y[3] -= (y[1]/r2)*verr[0];
        calcVerr();
       // printf(" --> %g\n", verrNorm/tol);

        return verrNorm <= tol;
    }

    Real calcAerr() const {
        assert(ydValid);
        const Real v2 = y[2]*y[2] + y[3]*y[3];
        aerr[0] = y[0]*yd[2] + y[1]*yd[3] + v2;
        return aerrNorm = fabs(aerr[0])/aUnits[0];
    }

    void calcWeights(Real& acc) {
        for (int i=0; i<NStates; ++i) {
            Real okErr = acc*yUnits[i];
            if (relativeOK[i]) okErr += acc*fabs(y[i]);
            weights[i] = 1./okErr;
        }
    }
    Real m, d, g;
    Vector yUnits, pUnits, vUnits, aUnits;
    bool relativeOK[NStates];
    mutable bool ydValid;
    mutable Vector yd;
    mutable Vector perr, verr, aerr;
    mutable Real perrNorm, verrNorm, aerrNorm;
    Real t; Vector y;
    Vector weights;
    Real solnAcc, consAcc, timeScale;
};

int main() {


    try {
        State scState;
        MultibodySystem scmbs;
        SinCos sc(1.);
        scmbs.addMatterSubsystem(sc); 
        scmbs.addForceSubsystem(EmptyForcesSubsystem());
        scmbs.realize(scState, Stage::Modeled);
        scState.updTime() = 0;
        Vector y(2); y[0] = 0.; y[1] = 1.;
        scState.updY() = y;

        ExplicitEuler ee(scmbs, scState);
        ee.setInitialStepSize(0.00001);

        ee.initialize();
        //std::cout << "t=" << ee.getT() << " y=" << ee.getY() << std::endl;

        while (false && scState.getTime() < 10.) {
            ee.step(scState.getTime() + 0.1);
            std::cout << "t=" << scState.getTime() << " y=" << scState.getY() << std::endl;
        }

        const Real mass = 5.;
        const Real length = 10.;
        const Real gravity = -9.8;
        const Real pi = acos(-1.);
        const Real halfPeriod = pi*sqrt(fabs(length/gravity));
        printf("Period should be %gs\n", 2*halfPeriod);
        MultibodySystem pendmbs;
        PointMass2dPendulum p(mass,length,gravity);
        pendmbs.addMatterSubsystem(p); 
        pendmbs.addForceSubsystem(EmptyForcesSubsystem());
        State pendState;
        pendmbs.realize(pendState, Stage::Modeled);

        Vector yp(4); 
        //yp[0]=sqrt(50.); yp[1]=-sqrt(50.); // -45 degrees
        yp[0] = length*cos(-0.9*(pi/2)); yp[1]=length*sin(-0.9*(pi/2)); // 9 degrees from bottom
        //yp[0]=0; yp[1]=-10;                  // -90 degrees (straight down)
        yp[2]=yp[3]=0;

        pendState.updTime() = 0;
        pendState.updY() = yp;

        //ExplicitEuler eep(p);
        RungeKuttaMerson eep(pendmbs, pendState);
        //eep.setInitialStepSize(0.00001);
        eep.setStopTime(100.);
        const Real acc = 1e-8;
        eep.setAccuracy(acc);
        //eep.setConstraintTolerance(1.);
      //  eep.setProjectEveryStep(true);


        if (!eep.initialize()) {
            printf("**** CAN'T SET ICS\n");
            exit(1);
        }

        while (true) {

            std::cout << "t=" << pendState.getTime() 
                << " yp=" << pendState.getY();

            pendmbs.realize(pendState, Stage::Reacting);
            std::cout << "   perr=" << p.getPositionErrorNorm()/eep.getConstraintTolerance()
                << "  verr=" << p.getVelocityErrorNorm()/eep.getConstraintTolerance()
                << "  aerr=" << p.getAccelerationErrorNorm()/eep.getConstraintTolerance()
                << " E=" << p.calcEnergy();
            std::cout << " hnext=" << eep.getPredictedNextStep() << std::endl;

            if (pendState.getTime()  >= 100.)
                break;

            Real h = 10.;
            //if (fabs(eep.getT()-floor(eep.getT()/halfPeriod+0.5)*halfPeriod) < 0.2)
             //   h = 0.001;
            if (!eep.step(pendState.getTime()  + h)) {
                printf("**** STEP FAILED t=%g -> %g\n", pendState.getTime() , pendState.getTime() +h);
                exit(1);
            }
        }

        printf("steps attempted=%ld, taken=%ld, size changes=%ld\n",
            eep.getStepsAttempted(), eep.getStepsTaken(), eep.getStepSizeChanges());
        printf("  error test failures=%ld, realize failures=%ld (with projection %ld)\n",
            eep.getErrorTestFailures(), eep.getRealizeFailures(), eep.getProjectionFailures());


    }
    catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    return 0;
}
