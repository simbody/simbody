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
 * The simple 2d pendulum example from the user's manual.
 */

#include "Simbody.h"
#include "simbody/internal/NumericalMethods.h"
#include "simbody/internal/VTKReporter.h"

#include "windows.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 1.5; // meters
static const Real initialTheta   = 20;             // degrees
static const Real expectedPeriod = 2*Pi*sqrt(d/g); // s

class MySimbodyPendulum : public MechanicalDAESystem {
public:
    MySimbodyPendulum() 
    {
        pendBodyNum =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                GroundBodyNum,           // parent body
                Transform(Vec3(1,1,1)),             // jt frame on parent             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes

        int pendBodyNum2 =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes


        int pendBodyNum3 =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes
       
    int theConstraint =
        pend.addConstantDistanceConstraint(0, Vec3(2,-3,0),
                                          pendBodyNum3, Vec3(0,-d/2,0),
                                           2);
        //pend.addCoincidentStationsConstraint(0, Vec3(2,0,0),
         //                                  pendBodyNum3, Vec3(0,-d/2,0));




        pend.realize(s, Stage::Built);
        pend.setUseEulerAngles(s, true);
        pend.realize(s, Stage::Modeled);
        nq = s.getQ().size();
        nu = s.getU().size();
        nz = 0; // these would be for controllers and other forces with states

        printf("System size nq=%d nu=%d nz=%d (ny=%d)\n", nq, nu, nz, size());
        y.resize(nq+nu+nz);       y.setToNaN();
        ydot.resize(nq+nu+nz);    ydot.setToNaN();
        weights.resize(nq+nu+nz); weights.setToNaN();
        solutionAccuracy = constraintAccuracy = 1e-3; // default is 0.1%

        // Set state defaults
        s.updTime() = 0.;
        setPendulumAngle(initialTheta);
        copyStateToY();

        calcWeights(solutionAccuracy);
        ydotReady = false;
    }

    const SimbodySubsystem& getSimbodySubsystem() const {return pend;}
    const State& getState() const {return s;}

    void setTime(Real t) {s.updTime() = t;}

    Real getPendulumAngle() const {
        return pend.getJointQ(s,pendBodyNum,0)/RadiansPerDegree;
    }

    void setPendulumAngle(Real angleInDegrees) {
        for (int i=0; i<3; ++i)
            pend.setJointQ(s,pendBodyNum,i,angleInDegrees*RadiansPerDegree);
        copyStateToY();
    }

    // Supply required virtual methods.

    int size() const {return nq+nu+nz;}

    void setAccuracy(const Real& solution, const Real& constraint) {
        solutionAccuracy = solution; constraintAccuracy = constraint;
        calcWeights(solutionAccuracy);
    }

    void setState(const Real& t, const Vector& newY) {
        s.updTime() = t;
        y = newY;   // TODO: shouldn't need this temporary copy
        copyYtoState();
        calcWeights(solutionAccuracy);
        ydotReady = false;
    }

    bool realize() const {
        pend.realize(s, Stage::Moving); // kinematics

        // calculate and apply forces
        pend.clearAppliedForces(s);
        pend.applyGravity(s, Vec3(0, -g, 0));

        // calculate Simbody derivatives
        pend.realize(s, Stage::Reacting);

        // calculate other derivatives if needed
        //    none here

        // group derivatives together (TODO: should handle this without copying)
        copyStateToYDot();

        ydotReady = true;
        return true;
    }

    bool realizeAndProject(bool& anyChange, bool force) {
        pend.realize(s, Stage::Configured);
        pend.enforceConfigurationConstraints(s);
        pend.realize(s, Stage::Moving);
        pend.enforceMotionConstraints(s);
        return realize();
    }

    Real          getTimescale() const {return 0.1;} // a tenth of a second
    const Real&   getT()         const {return s.getTime();}
    const Vector& getY()         const {return y;}
    const Vector& getWeights()   const {return weights;}
    const Vector& getYDot()      const {assert(ydotReady); return ydot;}

private:
    void calcWeights(Real acc) {
        // Assume units around 1 for everything, and relative tolerance OK
        for (int i=0; i<size(); ++i)
            weights[i] = 1./(acc*fabs(y[i]) + acc);
    }

    //TODO: these shouldn't be needed
    void copyStateToY() {
        y(0,nq)     = s.getQ();
        y(nq,nu)    = s.getU();
        //y(nq+nu,nz) = s.getZ();
    }

    void copyStateToYDot() const {
        ydot(0,nq)     = s.getQDot();
        ydot(nq,nu)    = s.getUDot();
        //ydot(nq+nu,nz) = s.getZDot();
    }

    void copyYtoState() {
        s.updQ() = y(0,nq);
        s.updU() = y(nq,nu);
        //s.updZ() = y(nq+nu,nz);
    }

private:
    SimbodySubsystem pend;

    Vector y, weights;
    int nq, nu, nz;
    int pendBodyNum;
    Real solutionAccuracy, constraintAccuracy;

    mutable State s;         // TODO: state shouldn't have to be mutable; needed
                             //   because we are (unnecessarily) using it to hold
                             //   forces.

    // These are "cache" items and must be mutable.
    mutable Vector ydot;
    mutable bool   ydotReady; // sanity check
};


int main(int argc, char** argv) {
    std::vector<State> saveEm;

    try { // If anything goes wrong, an exception will be thrown.
        Real start = initialTheta;
        if (argc > 1) sscanf(argv[1], "%lg", &start);
        printf("Pendulum starting at angle +%g degrees from vertical.\n", start);

        // Create a multibody system using Simbody.
        MySimbodyPendulum myPend;
        myPend.setPendulumAngle(start);

        // And a study using the Runge Kutta Merson integrator
        RungeKuttaMerson myStudy(myPend);
        myStudy.setAccuracy(1e-3);
        myStudy.setProjectEveryStep(true);

        MultibodySystem mbs(myPend.getSimbodySubsystem(), EmptyForcesSubsystem());

        VTKReporter display(mbs);
        //for (int i=1; i<myPend.getSimbodySubsystem().getNBodies(); ++i)
       //     display.addDecoration(i, Transform(Vec3(0,-d/2,0)), DecorativeSphere(0.2).setOpacity(0.3));

        //display.addDecoration(0,VTKReporter::Sphere(0.01),
        //    Transform(Vec3(0.5,-0.2,0.1)));

        // Run for 5 periods without output every dt seconds,
        // starting at theta=start degrees.

        const Real dt = 0.01; // output intervals

        printf("time  theta (deg)  (period should be %gs)\n", expectedPeriod);

        myStudy.setInitialConditions(myPend.getT(), myPend.getY());
        display.report(myPend.getState());
        for (;;) {
            printf("%5g %10.3g\n", myStudy.getT(), myPend.getPendulumAngle());
            display.report(myPend.getState());
            saveEm.push_back(myPend.getState());

           // if (myStudy.getT() >= 10*expectedPeriod)
             //   break;
    
            if (myStudy.getT() >= 10.)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(myStudy.getT() + dt);
        }

        for (int i=0; i < (int)saveEm.size(); ++i)
            display.report(saveEm[i]);
    } 
    catch (const exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);
    }
}

