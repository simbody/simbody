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
 * Test of high level multibody modeling objects for Simbody.
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
    ////////////////////////////////////
    // Create mass element prototypes //
    ////////////////////////////////////

    // Create point masses with constant values for their mass parameter.
    // These will require a "station" placement.
    PointMassElement    blue  ("blue",  2.5), 
                        orange("orange",1.),
                        green ("green", 0.1);

    // Create a cylinder mass element prototype. We're going
    // to give it constant radius and density, but leave the half-height
    // unspecified. So this will require a station, direction, and 
    // scalar parameter placement.
    CylinderMassElement tubeProto("tube");
    tubeProto.setRadius(0.1);
    tubeProto.setMass(100.);
  
    //////////////////////////////////
    // Create rigid body prototypes //
    //////////////////////////////////

    RigidBody upper("U");

    // Create some named stations and place them on the upper body.
    upper.addStation("greenPt",         Vec3(1, 0,0));
    upper.addStation("orangePt",        Vec3(0, 0,1));
    upper.addStation("pinnedBluePt",    Vec3(0, 1,0));
    upper.addStation("otherBluePt",     Vec3(1, 0,1));
    upper.addStation("leftAttachPt",    Vec3(0,-1,1));
    upper.addStation("rightAttachPt",   Vec3(1,-1,0));

    // Instantiate 'atom' mass elements and place them on their stations.
    upper.addMassElementLike(green,  "green",  upper["greenPt"]);
    upper.addMassElementLike(orange, "orange", upper["orangePt"]);
    upper.addMassElementLike(blue,   "blue1",  upper["pinnedBluePt"]);
    upper.addMassElementLike(blue,   "blue2",  upper["otherBluePt"]);

    // Create some joint frames at the appropriate stations (by default
    // they are aligned with the body frame).
    upper.addFrame("pinFrame",          upper["pinnedBluePt"]);
    upper.addFrame("leftBallFrame",     upper["leftAttachPt"]);
    upper.addFrame("rightBallFrame",    upper["rightAttachPt"]);

    // Now build the prototype for the lower bodies.

    RigidBody lower("L");

    // Create a parameter for the cylinder half-height, and a station at
    // the top center of the cylinder where we can hook the ball joint.
    lower.addRealParameter("halfHeight");
    lower.addStation      ("ballAttachPt",         lower["halfHeight"] * lower.y());
    lower.addFrame        ("upperAttachmentFrame", lower["ballAttachPt"]);

    // Now instantiate a tube on the body prototype.
    MassElement& tube = lower.addMassElementLike(tubeProto, "tube");
  
    // Place the center and axis, but leave the halfHeight parameter unresolved
    // because we want to control both with a single parameter.
    tube["center"].place(lower.getOrigin());

    //TODO: this should require an explicit normalize()
    tube["axis"].place(lower["ballAttachPt"] - lower.getOrigin());

    tube["halfLength"].place(lower["halfHeight"]);

    ////////////////////////////////////////////
    // Create an articulated multibody system //
    ////////////////////////////////////////////
	
    Multibody mbs("example1");

    mbs.addGroundBody();
    mbs.addBodyLike(upper, "upper");
    mbs.addBodyLike(lower, "left");
    mbs.addBodyLike(lower, "right");

    // Create a single parameter of the multibody which can be used
    // to control the two halfHeights together.
    mbs.addRealParameter("halfHeight");
    mbs["left/halfHeight"].place(mbs["halfHeight"]);
    mbs["right/halfHeight"].place(mbs["halfHeight"]);

    mbs.addJoint(PinJoint, "base2ground", 
                 mbs.getGroundFrame(),                  //reference frame
                 mbs["upper"]);                         //moving frame
    mbs.addJoint(BallJoint, "leftHipJoint",
                 mbs["upper/leftBallFrame"],         //reference frame
                 mbs["left/upperAttachmentFrame"]);  //moving frame
    mbs.addJoint(BallJoint, "rightHipJoint",
                 mbs["upper/rightBallFrame"],        //reference frame
                 mbs["right/upperAttachmentFrame"]); //moving frame


    mbs["halfHeight"].place(13.111);
    mbs.realize(Stage::Startup);

    cout << "***BODY MASSES:" << endl;
    cout << "mbs/upper/mass=" << mbs["upper/mass"].getValue() << endl;
    cout << "mbs/left/mass =" << mbs["left/mass"].getValue() << endl;
    cout << "mbs/right/mass=" << mbs["right/mass"].getValue() << endl;
    cout << endl << "***BODY CENTROIDS:" << endl;
    cout << "mbs/upper/centroid=" << mbs["upper/centroid"].getValue() << endl;
    cout << "mbs/left/centroid =" << mbs["left/centroid"].getValue() << endl;
    cout << "mbs/right/centroid=" << mbs["right/centroid"].getValue() << endl;
    cout << endl;
    cout << endl << "***BODY INERTIAS:" << endl;
    cout << "mbs/upper/inertia =" << mbs["upper/inertia"].getValue() << endl;
    cout << "mbs/left/inertia  =" << mbs["left/inertia"].getValue() << endl;
    cout << "mbs/right/inertia =" << mbs["right/inertia"].getValue() << endl;
    cout << endl << "***BODY CENTRAL INERTIAS:" << endl;
    cout << "upper =" << mbs["upper/centralInertia"].getValue() << endl;
    cout << "left  =" << mbs["left/centralInertia"].getValue() << endl;
    cout << "right =" << mbs["right/centralInertia"].getValue() << endl;
    cout << endl;

    try {cout << "left/tube/axis=" << mbs["left/tube/axis"].getValue() << endl;}
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    //mbs.checkSubsystemConsistency(0,-1,mbs);    
    //std::cout << "***MULTIBODY SYSTEM***" << std::endl;
    //std::cout << mbs << std::endl; //let’s see what we’ve got
    //std::cout << "***END OF MULTIBODY SYSTEM***" << std::endl;

    cout << "*** JOINTS ***" << endl;
    for (int i=0; i<mbs.getNSubsystems(); ++i)
        if (Joint::isInstanceOf(mbs[i]))
            cout << mbs[i].getFullName() << ":" << endl 
                 << "  reference: " << mbs[i]["reference"].getPlacement() 
                 << "  moving:    " << mbs[i]["moving"].getPlacement() 
                 << endl;

    cout << "*** BODIES ***" << endl;
    for (int i=0; i<mbs.getNSubsystems(); ++i)
        if (Body::isInstanceOf(mbs[i])) {
            cout << mbs[i].getFullName() << ":" << endl;
            Real tmass = 0.;
            Vec3 tcom(0.);
            for (int j=0; j<mbs[i].getNSubsystems(); ++j) {
                if (MassElement::isInstanceOf(mbs[i][j])) {
                    const MassElement& me = MassElement::downcast(mbs[i][j]);
                    Real mass = me.getMassMeasure().getValue();
                    Vec3 com  = me.getCentroidMeasure().getValue();
                    cout << "  massElt " << mbs[i][j].getName() << "mass=" << mass << " com=" << com << endl;
                    tmass += mass;
                    tcom += mass*com;
                }
            }
            cout << "... total mass=" << tmass << " centroid=" << tcom/tmass << endl;
        }


    ///////////////////////////
    // Build a RigidBodyTree //
    ///////////////////////////


    IVMSimbodyInterface instance(mbs);

    cout << "MatRotation(z)=" << MatRotation(UnitVec3(0,0,1)) << endl;
    cout << "MatRotation(" << UnitVec3(.1,.2,.3) << ")="
        << MatRotation(UnitVec3(.1,.2,.3)) << endl; 

    UnitVec3 u(0.1,-3.2,7);
    cout << "u=" << u << " u.perp=" << u.perp() 
         << " ~u*u.perp=" << ~u*u.perp() << " u%u.perp(),1-norm=" 
         << u%u.perp() << "," << 1-(u%u.perp()).norm() << endl;
    cout << "u*~u=" << u*~u << endl;
    cout << "u*~u.perp=" << u*~u.perp() << endl;
}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
