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
 * A one-body pendulum, to test proper frame alignment and basic
 * functioning of Simbody.
 */

/* Sketch:
 *
 *     |           \           | g
 *     *--          *--        v
 *    / G          / Ji
 *
 *
 *   |           |
 *   *==---------*==---------W
 *  / J         / B         weight
 *   <--- L/2 ---|--- L/2 --->
 *
 *
 * The pendulum is a massless rod with origin frame
 * B, joint attachment frame J, and a point mass W.
 * The rod length is L, with the joint and mass
 * located in opposite directions along the B
 * frame X axis.
 *
 * There is a frame Ji on Ground which will connect
 * to J via a torsion joint around their mutual z axis.
 * Gravity is in the -y direction of the Ground frame.
 * Note that Ji may not be aligned with G, and J may
 * differ from B so the reference configuration may 
 * involve twisting the pendulum around somewhat.
 */

#include "simbody/Simbody.h"
#include "SimbodyTree.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;

using namespace simtk;

int main() {
try {
    SimbodyTree pend;

    Real L = 5.; 
    Real m = 3.;
    TransformMat groundFrame;
    TransformMat jointFrame(Vec3(-L/2,0,0));
    MassProperties mprops(m, Vec3(L/2,0,0), InertiaMat(Vec3(L/2,0,0), m)+InertiaMat(1.,0.,0.));
    cout << "mprops about body frame: " << mprops.getMass() << ", " 
        << mprops.getCOM() << ", " << mprops.getInertia() << endl;

    int theBody = 
      pend.addRigidBody(0, groundFrame, 
                        //JointSpecification(JointSpecification::Pin, false),
                        JointSpecification(JointSpecification::Ball, false),
                        jointFrame, mprops);
    pend.realizeConstruction();
    SBState s = pend.getInitialState();

    // set Modeling stuff (s)
    pend.setUseEulerAngles(s, false); // this is the default
    pend.realizeModeling(s);

    pend.realizeConfiguration(s);
    TransformMat bodyConfig = pend.getBodyConfiguration(s, theBody);
    cout << "body frame: " << bodyConfig;

    pend.clearAppliedForces(s);
    pend.applyGravity(s, Vec3(0.,-9.8,0.));
    pend.applyJointForce(s, 1, 0, 147);

    pend.setJointU(s, 1, 0, 10.);

    cout << "after applying gravity, body forces=" << pend.getAppliedBodyForces(s) << endl;
    cout << "   joint forces=" << pend.getAppliedJointForces(s) << endl;

    pend.realizeMotion(s);

    SpatialVec bodyVel = pend.getBodyVelocity(s, theBody);
    cout << "body vel: " << bodyVel << endl;

    cout << "wXwXr=" << bodyVel[0] % (bodyVel[0] % Vec3(2.5,0,0)) << endl;

    pend.realizeReaction(s);

    SpatialVec bodyAcc = pend.getBodyAcceleration(s, theBody);
    cout << "body acc: " << bodyAcc << endl;

    pend.updQ(s) = Vector(4, &Vec4(1.,0.,0.,0.)[0]);
    //pend.updQ(s)[0] = -1.5; // almost hanging straight down
    pend.updU(s) = 0;

    const Real h = 0.0001;
    const Real tstart = 0.;
    const Real tmax = 10.;
    for (int step=0; ; ++step) { 
        const Real t = tstart + step*h;
        if (t > tmax) break;

        pend.clearAppliedForces(s);
        pend.realizeConfiguration(s);
        pend.realizeMotion(s); 
        if (!(step % 100))
            cout << t << " " 
                 << pend.getQ(s) << " " << pend.getU(s) 
                 << endl;
        Vector qdot = pend.getQDot(s);
        pend.applyGravity(s,Vec3(0,-9.8,0));
        pend.realizeReaction(s);
        Vector udot = pend.getUDot(s);

        //cout << "qdot=" << qdot << "  udot=" << udot << endl;
        pend.updQ(s) += h*qdot;
        pend.updU(s) += h*udot;
    }

}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}

    return 0;
}
