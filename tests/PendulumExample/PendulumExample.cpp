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
#include <cmath>
#include <cstdio>
#include <exception>

using namespace std;
using namespace SimTK;

static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real m = 5;  // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 0.5; // meters
static const Real initialTheta   = 10; // degrees
static const Real expectedPeriod = 2*Pi*sqrt(d/g); // s

int main() {
    try { // If anything goes wrong, an exception will be thrown.

        // First, create a MechanicalSubsystem using Simbody.
        SimbodySubsystem pend;
        const int pendBodyNum =
            pend.addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,0,0), 
                               Inertia(m, Vec3(0,-d/2,0))),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                GroundBodyNum,           // parent body
                Transform(),             // jt frame on parent (origin in this case)              
                Joint::Pin);             // joint type; pin always aligns z axes
        pend.finishConstruction();

        // Next, create a ForceSubsystem for use with the already-built MechanicalSubsystem.
        BasicForceElements forces(pend);
        forces.setGravity(Vec3(0,-g,0));
        forces.finishConstruction();

        // Combine the subsystems into a MultibodySystem and get the starting State.
        MultibodySystem mbs(pend,forces);
        State s = mbs.getDefaultState();

        // Create a Study; run for 2.5 periods without output every 0.01s,
        // starting at theta=initialTheta degrees.
        MultibodyDynamicsStudy study(mbs,s);
        printf("time  theta (deg) (period should be %gs)\n", expectedPeriod);

        pend.setJointQ(s,pendBodyNum,0,initialTheta*RadiansPerDegree);
        for (;;) {
            printf("%5g %10g\n", s.getT(), 
               pend.getJointQ(s,pendBodyNum,0)/RadiansPerDegree);

            if (s.getT() >= 2.5*expectedPeriod)
                break;
            study.advanceTimeBy(0.01);
        }
    } 
    catch (const exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);
    }
}
