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
 *     |           \  
 *     *--          *--
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

    // Add the joint frame. For now it is aligned with the rod frame.
    rod.addFrame("jointFrame",  
        MatRotation(UnitVec3(0,1,0)),
        rod["jointLocation"]);


    ////////////////////////////////////////////
    // Create an articulated multibody system //
    ////////////////////////////////////////////
	
    Multibody mbs("example1");

    mbs.addGroundBody();
    mbs.addBodyLike(rod, "rod");

    mbs["ground"].addStation("rodJointLocation", Vec3(10,0,0));
    mbs["ground"].addFrame("rodJointFrame",
        MatRotation(UnitVec3(0,1,0)),
        mbs["ground/rodJointLocation"]);

    mbs.addJoint(Joint::Pin, "rodJoint", 
                 mbs["ground/rodJointFrame"],        //reference frame
                 mbs["rod/jointFrame"]);             //moving frame


    ///////////////////////////
    // Build a RigidBodyTree //
    ///////////////////////////

    mbs["rod/h"].place(5);
    mbs["rod/m"].place(3);

    IVMSimbodyInterface instance(mbs);

}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
