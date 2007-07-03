/* KneeJointExample - Ajay Seth
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
#include "GCVSpline.h"
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
static const Real g = 0; //9.8; // meters/s^2; apply in –y direction
static const Real d = 1; // meters

int main(int argc, char** argv) {
try { // If anything goes wrong, an exception will be thrown.

    MultibodySystem         mbs;
    //UniformGravitySubsystem gravity(Vec3(0, -g, 0));
	int i = 0;

    SimbodyMatterSubsystem  pend;
    const Vec3 com(0, 0, 0); // in local frame of swinging body

	const BodyId base = pend.addRigidBody(
        MassProperties(m, com, m*Inertia::pointMassAt(com)),
        Vec3(0, d, 0),    // inboard joint location
        GroundId, GroundFrame,
        Mobilizer::Pin());

    // Put the subsystems into the system.
    mbs.setMatterSubsystem(pend);
    //mbs.addForceSubsystem(gravity);

	 // Add ligamentous force sub-system to this multi-body system.
    GeneralForceElements userForceElements;
    mbs.addForceSubsystem( userForceElements );

    State s;


    // Although "new" was used to allocate this UserForce, do not "delete" it.
    // This bug will be fixed in the next version of Simbody so it can take an object from the stack or heap.
    // For now, addUserForce takes ownership of the allocated item and takes care of deleting it at the end.

	// Constant point force applied at COM perpendicular to base pendulumn
	UserMuscleForce *muscle = new UserMuscleForce(10, Vec3(0, 0, 0), base);
	userForceElements.addUserForce( muscle );
	
	mbs.realize(s);

    // Create a study using the Runge Kutta Merson integrator
    RungeKuttaMerson myStudy(mbs, s);
    myStudy.setAccuracy(1e-5);
    //CPodesIntegrator myStudy(mbs, s);
    //myStudy.setAccuracy(1e-4);

    // Visualize with VTK.
    VTKReporter display(mbs);

    const Real dt = 0.001; // output intervals
    const Real finalTime = 100;

	Real startAngle = 0;
	Real startAngularVel = 0; //degrees per second

    //for (Real startAngle = 30; startAngle <= 90; startAngle += 10) {
        printf("time  theta    udot   energy           *************\n");
        s.updTime() = 0;
	    pend.setMobilizerQ(s, base, 0, startAngle);
        pend.setMobilizerU(s, base, 0, startAngularVel);
        myStudy.initialize();

        cout << "MassProperties in G=" << pend.calcBodyMassPropertiesInBody(s,base,GroundId);
        cout << "Spatial Inertia    =" << pend.calcBodySpatialInertiaMatrixInGround(s,base);

		Vector acc;

        for (;;) {
			acc = pend.getUDot(s);
            printf("%8.4g %10.4g %10.8g %10.8g\n", s.getTime(), pend.getMobilizerQ(s,base,0)*Rad2Deg,
                acc[0], mbs.getEnergy(s));

            display.report(s);
            if (s.getTime() >= finalTime)
                break;

            // TODO: should check for errors or have or teach RKM t5o throw. 
            myStudy.step(s.getTime() + dt);
        }
    //}
} 
catch (const exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
}
}

