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

static const Real Pi      = (Real)SimTK_PI, 
                  Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const int  GroundBodyNum = 0; // ground is always body 0
static const Transform GroundFrame;

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters
static const Real initialTheta   = 20;             // degrees


int main(int argc, char** argv) {
    std::vector<State> saveEm;

    try { // If anything goes wrong, an exception will be thrown.
        Real start = initialTheta;
        if (argc > 1) sscanf(argv[1], "%g", &start);
        printf("Pendulum starting at angle +%g degrees from vertical.\n", start);


        // Create a multibody system using Simbody.
        SimbodyMatterSubsystem pend;

        const Vec3 weightLocation(0, -d/2, 0); // in local frame of swinging body

        const int swinger = pend.addRigidBody(
            MassProperties(m, weightLocation, m*Inertia::pointMassAt(weightLocation)),
            Vec3(0, d/2, 0),    // inboard joint location
            GroundBodyNum, GroundFrame,
            Mobilizer::Pin); // rotates around common z axis

        GeneralForceElements forces;
        //forces.addGlobalEnergyDrain(1000);
       /* forces.addTwoPointLinearSpring(0, -attachPt,
                                       myRNA.getNBodies()-1, Vec3(0),
                                       1000.,  // stiffness
                                       1.);    // natural length
        */

        MultibodySystem mbs;
        mbs.setMatterSubsystem(pend);
        mbs.addForceSubsystem(forces);

        UniformGravitySubsystem ugs(Vec3(0, -g, 0));
        mbs.addForceSubsystem(ugs);

        State s;
        mbs.realize(s, Stage::Topology);
        //myRNA.setUseEulerAngles(s,true);
        mbs.realize(s, Stage::Model);

        //ugs.updGravity(s) *= 10;
        ugs.disableGravity(s);
        ugs.enableGravity(s);
       // ugs.updZeroHeight(s) = -0.8; // to change how PE is calculated
        //cout << "STATE AS MODELED: " << s;
       
        pend.setMobilizerQ(s, swinger, 0, start*Deg2Rad);

        // And a study using the Runge Kutta Merson integrator
        bool suppressProject = false;
        RungeKuttaMerson myStudy(mbs, s, suppressProject);
        myStudy.setAccuracy(1e-8);
        //myStudy.setConstraintTolerance(1e-3);
        //myStudy.setProjectEveryStep(false);

        VTKReporter display(mbs);
        display.setDefaultBodyColor(swinger, Cyan);
       // display.addDecoration(swinger, Transform(weightLocation), 
     //       DecorativeSphere(d/8).setColor(Blue).setOpacity(.2));

        //DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
        //display.addRubberBandLine(0, attachPt,myRNA.getNBodies()-1,Vec3(0), rbProto);

        const Real actualG = ugs.getGravity(s).norm() * ugs.isEnabled(s);
        printf("d=%g, g=%g -> Expected period: %g seconds\n", 
            d, actualG,
            2*Pi*std::sqrt(d/actualG));

        const Real dt = 0.01; // output intervals

        printf("time  nextStepSize\n");

        s.updTime() = 0;
        for (int i=0; i<100; ++i)
            saveEm.push_back(s);    // delay
        display.report(s);

        myStudy.initialize();
        saveEm.push_back(s);
        for (int i=0; i<100; ++i)
            saveEm.push_back(s);    // delay
        display.report(s);
        for (;;) {
            printf("%5g q=%10.4g u=%10.4g hNext=%g\n", s.getTime(), 
                pend.getMobilizerQ(s,swinger,0)*Rad2Deg, pend.getMobilizerU(s,swinger,0)*Rad2Deg,
                myStudy.getPredictedNextStep());
            printf("      E=%14.8g (pe=%10.4g ke=%10.4g)\n",
                mbs.getEnergy(s), mbs.getPotentialEnergy(s), mbs.getKineticEnergy(s));

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

