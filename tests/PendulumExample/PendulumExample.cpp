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
#include "simbody/internal/DecorativeGeometry.h"
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
static const Real initialTheta   = 30;             // degrees
static const Real expectedPeriod = 2*Pi*sqrt(d/g); // s

class MySimbodyPendulum : public SimbodySubsystem {
public:
    MySimbodyPendulum() 
    {
        pendBodyNum =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                GroundBodyNum,           // parent body
                Transform(Vec3(1,1,1)),             // jt frame on parent             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes

        int pendBodyNum2 =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes


        int pendBodyNum3 =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                JointSpecification(JointSpecification::Ball, false)); // joint type; pin always aligns z axes
       
    int theConstraint =
       addConstantDistanceConstraint(0, Vec3(2,-3,0),
                                          pendBodyNum3, Vec3(0,-d/2,0),
                                           1.5);
       //pend.addCoincidentStationsConstraint(0, Vec3(2,0,0),
         //                                  pendBodyNum3, Vec3(0,-d/2,0));

    }

    Real getPendulumAngle(const State& s) const {
        const Vec4 aa = getMobilizerConfiguration(s,pendBodyNum).R().convertToAngleAxis();
        return aa[0]/RadiansPerDegree;
    }

    // Assume rotation around z
    void setPendulumAngle(State& s, Real angleInDegrees) {
        const Vec4 aa(angleInDegrees*RadiansPerDegree,0, 0, 1);
        Quaternion q; q.setToAngleAxis(aa);
        setMobilizerConfiguration(s,pendBodyNum,Transform(RotationMat(q)));
    }
private:
    int pendBodyNum;
};


int main(int argc, char** argv) {
    std::vector<State> saveEm;

    try { // If anything goes wrong, an exception will be thrown.
        Real start = initialTheta;
        if (argc > 1) sscanf(argv[1], "%lg", &start);
        printf("Pendulum starting at angle +%g degrees from vertical.\n", start);

        // Create a multibody system using Simbody.
        MySimbodyPendulum myPend;
        TwoPointSpringSubsystem forces(0,Vec3(0),1,Vec3(0),0.,1.);
        State s;
        MultibodySystem mbs(myPend,forces);
        mbs.realize(s, Stage::Built);
        //myPend.setUseEulerAngles(s,true);
        mbs.realize(s, Stage::Modeled);
        forces.updGravity(s) = Vec3(0, -g, 0);
        cout << "STATE AS MODELED: " << s;
       
        myPend.setPendulumAngle(s, start);

        // And a study using the Runge Kutta Merson integrator
        RungeKuttaMerson myStudy(mbs, s);
        myStudy.setAccuracy(1e-2);
        myStudy.setProjectEveryStep(true);

        VTKReporter display(mbs);
        DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
        display.addRubberBandLine(0, Vec3(2,-3,0), 3, Vec3(0,-d/2,0), rbProto);

        //display.addDecoration(2, Transform(), DecorativeCircle(1).setColor(Yellow).setLineThickness(5));
        //for (int i=1; i<myPend.getSimbodySubsystem().getNBodies(); ++i)
       //     display.addDecoration(i, Transform(Vec3(0,-d/2,0)), DecorativeSphere(0.2).setOpacity(0.3));

        //display.addDecoration(0,VTKReporter::Sphere(0.01),
        //    Transform(Vec3(0.5,-0.2,0.1)));

        // Run for 5 periods without output every dt seconds,
        // starting at theta=start degrees.

        const Real dt = 0.025; // output intervals

        printf("time  theta (deg)  (period should be %gs)\n", expectedPeriod);

        myStudy.initialize();
        s.updTime() = 0;
        display.report(s);
        for (;;) {
            printf("%5g %10.3g hNext=%g\n", s.getTime(), myPend.getPendulumAngle(s), myStudy.getPredictedNextStep());
            display.report(s);
            saveEm.push_back(s);

           // if (myStudy.getT() >= 10*expectedPeriod)
             //   break;
    
            if (s.getTime() >= 10)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(s.getTime() + dt);
        }

        for (int i=0; i < (int)saveEm.size(); ++i)
            display.report(saveEm[i]);
    } 
    catch (const exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);
    }
}

