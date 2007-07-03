/* EllipsoidTest - Ajay Seth
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
 * A test of the ellipsoid joint.
 */

#include "SimTKsimbody.h"
#include "GCVSpline.h"

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

static const Real m = 10;   // kg
static const Real g = 9.81; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters

bool useEulerAngs = true;

int main(int argc, char** argv) {
try { // If anything goes wrong, an exception will be thrown.

    MultibodySystem         mbs;
    UniformGravitySubsystem gravity(Vec3(0, 0, -g));
	int i = 0;

    SimbodyMatterSubsystem  pend;
    const Vec3 com(0, 0, 0); // in local frame of swinging body

	const BodyId base = pend.addRigidBody(
        MassProperties(m, com, m*Inertia::cylinderAlongZ(d/4, d)),
        Vec3(0, 0, d),    // inboard joint location
        GroundId, GroundFrame,
        Mobilizer::Ellipsoid(0.5, 0.25, 0.1)); //Ball()); // 

    // Put the subsystems into the system.
    mbs.setMatterSubsystem(pend);
    mbs.addForceSubsystem(gravity);

    State s;

	mbs.realize(s, Stage::Topology);
	pend.setUseEulerAngles(s,useEulerAngs);
    mbs.realize(s, Stage::Model);

    printf("# quaternions in use = %d\n", pend.getNQuaternionsInUse(s));
    for (BodyId i(0); i<pend.getNBodies(); ++i) {
        printf("body %2d: using quat? %s; quat index=%d\n",
            (int)i, pend.isUsingQuaternion(s,i) ? "true":"false", 
            pend.getQuaternionIndex(s,i));
    }

    // Create a study using the Runge Kutta Merson integrator
	    // And a study using the Runge Kutta Merson integrator
    bool suppressProject = false;
    RungeKuttaMerson myStudy(mbs, s, suppressProject);
    myStudy.setAccuracy(1e-8);
    myStudy.setConstraintTolerance(1e-5);
    myStudy.setProjectEveryStep(false);

    // Visualize with VTK.
    VTKReporter display(mbs);

    const Real dt = 0.001; // output intervals
    const Real finalTime = 100;

	Real xAngle = 90*Deg2Rad;
	Real yAngle = 45*Deg2Rad;
	Real startAngularVel = 0*Deg2Rad; //degrees per second

    //for (Real startAngle = 30; startAngle <= 90; startAngle += 10) {
        s.updTime() = 0;

		if(useEulerAngs) { //Using 1-2-3 Euler angles
			pend.setMobilizerQ(s, base, 0, xAngle);
			pend.setMobilizerQ(s, base, 1, yAngle);
			pend.setMobilizerQ(s, base, 2, 0);
		}
		else{ // Quaternions
			Rotation R_GB;
			Quaternion quat;

			R_GB = R_GB.aboutXThenNewY(xAngle, yAngle);
			quat = R_GB.convertToQuaternion();
			for(int i = 0; i < 4; i++)
				pend.setMobilizerQ(s, base, i, quat[i]);
		}

		pend.setMobilizerU(s, base, 0, 0);
		pend.setMobilizerU(s, base, 1, 0);
		pend.setMobilizerU(s, base, 2, startAngularVel);

		myStudy.initialize();

        cout << "MassProperties in G=" << pend.calcBodyMassPropertiesInBody(s,base,GroundId);
        cout << "Spatial Inertia    =" << pend.calcBodySpatialInertiaMatrixInGround(s,base);

		Vector acc;
		Vec3 angs;
		Vec3 pos;

		printf("time  local-Xrot  local-Yrot  local-Zrot(spin)   X		Y		Z		energy \n");
        for (;;) {
			mbs.realize(s);
			acc = pend.getUDot(s);
			angs = pend.getBodyRotation(s, base).convertToBodyFixed123()*Rad2Deg;
			pos = pend.getBodyTransform(s, base).T(); 

	        printf("%8.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.10f\n", s.getTime(), angs[0], angs[1], angs[2], 
				   pos[0], pos[1], pos[2], mbs.getEnergy(s));

	        //printf("%8.4g %10.4g %10.4g %10.4g %10.8g\n", s.getTime(), pend.getMobilizerQ(s,base,0)*Rad2Deg,
 			//	pend.getMobilizerQ(s,base,1)*Rad2Deg, pend.getMobilizerQ(s,base,2)*Rad2Deg, mbs.getEnergy(s));
	
			//printf("%8.4g %10.4g %10.4g %10.4g %10.8g\n", s.getTime(), acc[0],
			//	acc[1], acc[2], mbs.getEnergy(s));

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

