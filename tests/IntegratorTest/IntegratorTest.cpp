/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
#include "Simmatrix.h"
//#include "Simbody.h"

#include "simbody/internal/NumericalMethods.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

class SinCos : public MechanicalDAESystem {
public:
    SinCos(const Real& omega) : w(omega), y(2), yd(2), weights(2), units(2), t(CNT<Real>::getNaN()), 
        solnAcc(1e-3), consAcc(1e-3), timeScale(1e-3) 
    {
        units[0] = units[1] = 1.;
        relativeOK[0] = relativeOK[1] = true;
        ydValid = false;
    }

    long size() const {return 2;}
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
    bool realizeAndProject(bool& chg) { realize(); chg=false; return true;}
    const Real& getTimescale() const {return timeScale;}
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
 *     x'' = -x*L/m
 *     y'' = g - y*L/m
 *               
 */  
class PointMass2dPendulum : public MechanicalDAESystem {
    static const int NStates = 4;
    static const int NConstraints = 1;
public:
    PointMass2dPendulum(const Real& mass, const Real& length, const Real& gravity) 
        : m(mass), d(length), g(gravity)
    {
        y.resize(NStates); yd.resize(NStates);
        perr.resize(NConstraints); verr.resize(NConstraints); aerr.resize(NConstraints);
        weights.resize(NStates); 
        yUnits.resize(NStates);
        cUnits.resize(NConstraints);
        t = CNT<Real>::getNaN();
        solnAcc = consAcc = timeScale = 1e-3;
        yUnits = 1; 
        cUnits = 1;
        yUnits[0] = yUnits[1] = 1;
        yUnits[2] = yUnits[3] = 1*timeScale; 
        cUnits[0] = 1.;
        relativeOK[0] = relativeOK[1] = false;
        relativeOK[2] = relativeOK[3] = true;
        ydValid = false;
    }

    long size() const {return 4;}
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
        calcPerr(); calcVerr();
        yd[0]=y[2]; yd[1]=y[3]; // qdot=u
        const Real r2 = y[0]*y[0] + y[1]*y[1]; // x^2+y^2
        const Real v2 = y[2]*y[2] + y[3]*y[3]; // x'^2 + y'^2
        const Real lambda = (m/r2)*(y[1]*g + v2);
        yd[2] =   - y[0]*lambda/m;
        yd[3] = g - y[1]*lambda/m;
        ydValid=true; 
        calcAerr();
        return true;
    }
    bool realizeAndProject(bool& chg) { 
        realize(); 
        chg=false;  // TODO: not yet
        return true;
    }
    const Real& getTimescale() const {return timeScale;}
    const Real& getT() const {return t;}
    const Vector& getY() const {return y;}
    const Vector& getWeights() const {return weights;}
    const Vector& getYDot() const {assert(ydValid);return yd;}
    const Vector& getPositionError() const {assert(ydValid);return perr;}
    const Vector& getVelocityError() const {assert(ydValid);return verr;}
    const Vector& getAccelerationError() const {assert(ydValid);return aerr;}

private:
    // (r^2 - d^2)/2 = 0    (perr)
    //     xx' + yy' = 0    (verr)
    //   xx'' + yy'' = -v^2 (aerr)
    void calcPerr() const {
        const Real r2 = y[0]*y[0] + y[1]*y[1];
        perr[0] = (r2 - d*d) / 2;
    }
    void calcVerr() const {
        verr[0] = y[0]*y[2] + y[1]*y[3];
    }
    void calcAerr() const {
        assert(ydValid);
        const Real v2 = y[2]*y[2] + y[3]*y[3];
        aerr[0] = y[0]*yd[2] + y[1]*yd[3] + v2;
    }

    void calcWeights(Real& acc) {
        for (int i=0; i<NStates; ++i) {
            Real okErr = acc*yUnits[i];
            if (relativeOK[i]) okErr += acc*fabs(y[i]);
            weights[i] = 1./okErr;
        }
    }
    Real m, d, g;
    Vector yUnits;
    Vector cUnits;
    bool relativeOK[NStates];
    mutable bool ydValid;
    mutable Vector yd;
    mutable Vector perr, verr, aerr;
    Real t; Vector y;
    Vector weights;
    Real solnAcc, consAcc, timeScale;
};

int main() {


    try {

        SinCos sc(1.);
        ExplicitEuler ee(sc);
        ee.setInitialStepSize(0.00001);
        Vector y(2); y[0] = 0.; y[1] = 1.;
        ee.setInitialConditions(0., y);
        //std::cout << "t=" << ee.getT() << " y=" << ee.getY() << std::endl;

        while (false && ee.getT() < 10.) {
            ee.step(ee.getT() + 0.1);
            std::cout << "t=" << ee.getT() << " y=" << ee.getY() << std::endl;
        }

        const Real mass = 5.;
        const Real length = 10.;
        const Real gravity = -9.8;
        const Real pi = acos(-1.);
        const Real halfPeriod = pi*sqrt(fabs(length/gravity));
        printf("Period should be %gs\n", 2*halfPeriod);
        PointMass2dPendulum p(mass,length,gravity);
        ExplicitEuler eep(p);
        eep.setInitialStepSize(0.0001);
        Vector yp(4); 
        yp[0]=sqrt(50.); yp[1]=-sqrt(50.); // -45 degrees
        //yp[0]=0; yp[1]=-10;                  // -90 degrees (straight down)
        yp[2]=yp[3]=0;
        eep.setInitialConditions(0., yp);

        while (true) {

            std::cout << "t=" << eep.getT() 
                << " yp=" << eep.getY();
            p.setState(eep.getT(), eep.getY());
            p.realize();
            std::cout << "   perr=" << p.getPositionError()
                << "  verr=" << p.getVelocityError()
                << "  aerr=" << p.getAccelerationError() 
                << std::endl;

            if (eep.getT() >= 10.)
                break;

            Real h = 0.1;
            if (eep.getT()-floor(eep.getT()/halfPeriod)*halfPeriod < 0.2)
                h = 0.001;
            eep.step(eep.getT() + h);
        }


    }
    catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
    }

    return 0;
}
