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

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi      = (Real)SimTK_PI, 
                  Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const Transform GroundFrame;

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters

int main(int argc, char** argv) {
try { // If anything goes wrong, an exception will be thrown.

    MultibodySystem         mbs;
    UniformGravitySubsystem gravity(Vec3(0, -g, 0));

    SimbodyMatterSubsystem  pend;
    const Vec3 weightLocation(0, -d/2, 0); // in local frame of swinging body

    const BodyId swinger = pend.addRigidBody(
        MassProperties(m, weightLocation, m*Inertia::pointMassAt(weightLocation)),
        Vec3(0, d/2, 0),    // inboard joint location
        GroundId, GroundFrame,
        Mobilizer::Pin); // rotates around common z axis

    // Put the subsystems into the system.
    mbs.setMatterSubsystem(pend);
    mbs.addForceSubsystem(gravity);

    State s;
    mbs.realize(s); // define appropriate states for this System


    // Create a study using the Runge Kutta Merson integrator
    RungeKuttaMerson myStudy(mbs, s);
    myStudy.setAccuracy(1e-2);
    //CPodesIntegrator myStudy(mbs, s);
    //myStudy.setAccuracy(1e-4);

    // Visualize with VTK.
    VTKReporter display(mbs);

    // Add a blue sphere around the weight.
    display.addDecoration(swinger, weightLocation, 
          DecorativeSphere(d/8).setColor(Blue).setOpacity(.2));

    const Real expectedPeriod = 2*Pi*std::sqrt(d/g);
    printf("Expected period: %g seconds\n", expectedPeriod);

    const Real dt = 0.01; // output intervals
    const Real finalTime = 1*expectedPeriod;

    for (Real startAngle = 10; startAngle <= 90; startAngle += 10) {
        printf("time  theta      energy           *************\n");
        s.updTime() = 0;
        pend.setMobilizerRotation(s, swinger, Rotation::aboutZ(startAngle*Deg2Rad));
        pend.setMobilizerVelocity(s, swinger, SpatialVec(Vec3(0),Vec3(0)));
        //pend.setMobilizerQ(s, swinger, 0, startAngle*Deg2Rad);
        //pend.setMobilizerU(s, swinger, 0, 0);
        myStudy.initialize();

        cout << "MassProperties in B=" << pend.calcBodyMassPropertiesInBody(s,swinger,swinger);
        cout << "MassProperties in G=" << pend.calcBodyMassPropertiesInBody(s,swinger,GroundId);
        cout << "Spatial Inertia    =" << pend.calcBodySpatialInertiaMatrixInGround(s,swinger);

        for (;;) {
            printf("%5g %10.4g %10.8g\n", s.getTime(), pend.getMobilizerQ(s,swinger,0)*Rad2Deg,
                mbs.getEnergy(s));

            display.report(s);
            if (s.getTime() >= finalTime)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(s.getTime() + dt);
        }
    }
} 
catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
}
}

