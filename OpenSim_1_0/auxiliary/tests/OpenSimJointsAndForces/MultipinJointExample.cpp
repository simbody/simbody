/* MultipinJointExample - Ajay Seth
 * Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
#include "UserLigamentForce.h"
#include "UserMuscleForce.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;
using namespace OpenSim;

static const Real Pi      = (Real)SimTK_PI, 
                  Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN,
                  Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;

static const Transform GroundFrame;
static const Transform SwingerFrame;

static const Real m = 1;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.4; // meters
const Vec3 com(0, 0, 0); // in local frame of swinging body

const MassProperties massless  = MassProperties(0, com, Inertia(0)); 

int main(int argc, char** argv) {
try { // If anything goes wrong, an exception will be thrown.

    MultibodySystem         mbs;
    UniformGravitySubsystem gravity(Vec3(0, -g, 0));

    SimbodyMatterSubsystem  pend;
	Transform aboutZ = Transform(Rotation(UnitVec3(0, 0, 1)), Vec3(0));  //rotation about Z
	Transform aboutX = Transform(Rotation(UnitVec3(1, 0, 0)), Vec3(0));  //rotation about X
	Transform aboutY = Transform(Rotation(UnitVec3(0, 1, 0)), Vec3(0));  //rotation about Y

	const BodyId base_1 = pend.addRigidBody(
		massless,
        aboutZ,    // inboard joint location
        GroundId, aboutZ,
        Mobilizer::Pin());

	const BodyId base_2 = pend.addRigidBody(
		massless,
        aboutX,    // inboard joint location
        base_1, aboutX,
        Mobilizer::Pin());

	const BodyId base = pend.addRigidBody(
        MassProperties(2*m, com, m*Inertia::cylinderAlongY(d/3, d)),
        aboutY,    // inboard joint location
        base_2, aboutY,
        Mobilizer::Pin());

    const BodyId swinger = pend.addRigidBody(
        MassProperties(m, com, m*Inertia::cylinderAlongY(d/5, d)),
        Vec3(0, d, 0),    // inboard joint location
        base, SwingerFrame,
        Mobilizer::Pin()); //

    // Put the subsystems into the system.
    mbs.setMatterSubsystem(pend);
    mbs.addForceSubsystem(gravity);

	 // Add ligamentous force sub-system to this multi-body system.
    GeneralForceElements userForceElements;
    mbs.addForceSubsystem( userForceElements );

    State s;

	mbs.realize(s);

    // Create a study using the Runge Kutta Merson integrator
    RungeKuttaMerson myStudy(mbs, s);
    myStudy.setAccuracy(1e-5);
    //CPodesIntegrator myStudy(mbs, s);
    //myStudy.setAccuracy(1e-4);

    // Visualize with VTK.
    VTKReporter display(mbs);

    // Add a blue sphere around the weight.
    display.addDecoration(swinger, com, 
          DecorativeSphere(d/20).setColor(Blue).setOpacity(.2));

    const Real dt = 0.001; // output intervals
    const Real finalTime = 100;

	Real startAngle = -80;
	Real startAnglularVel = 0; //degrees per second

    //for (Real startAngle = 30; startAngle <= 90; startAngle += 10) {
        printf("time  theta      energy           *************\n");
        s.updTime() = 0;
	    pend.setMobilizerQ(s, base, 0, 0);
        pend.setMobilizerU(s, base, 0, 0);
        pend.setMobilizerQ(s, swinger, 0, startAngle*Deg2Rad);
        pend.setMobilizerU(s, swinger, 0, startAnglularVel*Deg2Rad);
        myStudy.initialize();

        cout << "MassProperties in B=" << pend.calcBodyMassPropertiesInBody(s,swinger,swinger);
        cout << "MassProperties in G=" << pend.calcBodyMassPropertiesInBody(s,swinger,GroundId);
        cout << "Spatial Inertia    =" << pend.calcBodySpatialInertiaMatrixInGround(s,swinger);

		Vector acc;

        for (;;) {
			acc = pend.getUDot(s);
            printf("%8.4g %10.4g %10.8g %10.8g\n", s.getTime(), pend.getMobilizerQ(s,swinger,0)*Rad2Deg,
                acc[0], acc[1]); //mbs.getEnergy(s));

            display.report(s);
            if (s.getTime() >= finalTime)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(s.getTime() + dt);
        }
    //}
} 
catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
}
}

