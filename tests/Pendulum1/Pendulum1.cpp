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
 *   <--- h/2 ---|--- h/2 --->
 *
 *
 * The pendulum is a massless rod with origin frame
 * B, joint attachment frame J, and a point mass W.
 * The rod length is h, with the joint and mass
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
#include "simbody/IVMSimbodyInterface.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;

using namespace simtk;

static void tryPerp() {
    UnitVec3 u(1.,-2.,.3);
    UnitVec3 up(u.perp());
    cout << "u=" << u << " u.perp=" << up << " u*u.perp=" << ~u*up << endl;
}

int main() {
try {
    //////////////////////////////////
    // Create rigid body prototype  //
    //////////////////////////////////

    RigidBody rod("P");

    rod.addRealParameter("h");
    rod.addRealParameter("m"); // for the mass

    // Create stations at which we will locate the J frame and 
    rod.addStation("weightLocation", rod.x()*rod["h"]/2);
    rod.addStation("jointLocation", -rod.x()*rod["h"]/2);

    // Add the weight and place its mass on m.
    rod.addFeatureLike(PointMassElement("W"), "weight", rod["weightLocation"]);
    rod["weight/mass"].place(rod["m"]);

    // Add the joint frame. 
    
    // Here it is aligned with the rod frame (rod is horizontal as pictured).
    rod.addFrame("jointFrame",  
        MatRotation(),
        rod["jointLocation"]);

    // Here it is aligned so that the rod is hanging straight down at 0.
    const Mat33 jj(Vec3(0,1,0),Vec3(-1,0,0),Vec3(0,0,1));
    //rod.addFrame("jointFrame",  
    //    reinterpret_cast<const MatRotation&>(jj),
    //    rod["jointLocation"]);

    ////////////////////////////////////////////
    // Create an articulated multibody system //
    ////////////////////////////////////////////
	
    Multibody mbs("example1");

    mbs.addGroundBody();
    mbs.addBodyLike(rod, "rod");

    mbs["ground"].addStation("rodJointLocation", Vec3(10,0,0));
    mbs["ground"].addFrame("rodJointFrame",
        MatRotation(),
        mbs["ground/rodJointLocation"]);

    mbs.addJoint(Joint::Pin, "rodJoint", 
                 mbs["ground/rodJointFrame"],        //reference frame
                 mbs["rod/jointFrame"]);             //moving frame


    ///////////////////////////
    // Build a RigidBodyTree //
    ///////////////////////////

    mbs["rod/h"].place(5);
    mbs["rod/m"].place(3);

   // IVMRigidBodyModel pend(mbs);
   // DynamicInvestigation d(pend);


    mbs.realize(Stage::Startup);
    //cout << "MBS=" << mbs << endl;

    IVMSimbodyInterface instanceOld(mbs, true);  // old style
    IVMSimbodyInterface instanceNew(mbs, false); // new style
    State sOld = instanceOld.getDefaultState();
    State sNew = instanceNew.getDefaultState();
    sOld.updQ()[0] = -1.5; // almost hanging straight down
    sNew.updQ()[0] = -1.5; // almost hanging straight down

    const Real h = 0.0001;
    const Real tstart = 0.;
    const Real tmax = 10.;
    cout << "     OLD         NEW    " << endl;
    for (int step=0; ; ++step) { 
        const Real t = tstart + step*h;
        if (t > tmax) break;

        Vector_<SpatialVec>  bodyForces;
        Vector               hingeForces;
        instanceOld.clearForces(bodyForces,hingeForces);
        instanceOld.realizeParameters(sOld);    instanceNew.realizeParameters(sNew);
        instanceOld.realizeConfiguration(sOld); instanceNew.realizeConfiguration(sNew);
        if (!(step % 100))
            cout << t << " " 
                 << sOld.getQ()[0] << " " << sOld.getU()[0] 
                 << "        \t" << sNew.getQ()[0] << " " << sNew.getU()[0] 
                 << endl;
        instanceOld.applyGravity(sOld,Vec3(0,-9.8,0),bodyForces);
        instanceOld.realizeMotion(sOld); 
        instanceNew.realizeMotion(sNew);
        Vector udotOld = instanceOld.calcUDot(sOld,bodyForces,hingeForces);
        Vector udotNew = instanceNew.calcUDot(sNew,bodyForces,hingeForces);
        //cout << "udot=" << udot << endl;
        sOld.updQ() += h*sOld.getU();
        sOld.updU() += h*udotOld;

        sNew.updQ() += h*sNew.getU();
        sNew.updU() += h*udotNew;
    }

}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
