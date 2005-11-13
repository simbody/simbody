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

#include "SimbodyCommon.h"
#include "Placement.h"
#include "Feature.h"
#include "MassElement.h"
#include "Body.h"

#include <string>
#include <iostream>
using std::cout;
using std::endl;

using namespace simtk;

int main() {
try {
    Frame f("frame1");
    cout << f;

    ////////////////////////////////////
    // Create mass element prototypes //
    ////////////////////////////////////

    // Create point masses with constant values for their mass parameter.
    // These will require a "station" placement.
    PointMassElement    blue  ("blue",  2.5), 
                        orange("orange",1.),
                        green ("green", 0.1);

    cout << "blue=" << blue;
    cout << "orange=" << orange;
    cout << "green=" << green;

    // Create a cylinder mass element prototype. We're going
    // to give it constant radius and density, but leave the half-height
    // unspecified. So this will require a station, direction, and 
    // scalar parameter placement.
    CylinderMassElement tubeProto("tube");
    tubeProto.setRadius(0.1);
    tubeProto.setMass(100.);

    cout << "tube=" << tubeProto;
  
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

    XXX no placements above


    // Instantiate 'atom' mass elements and place them on their stations.
    upper.addMassElementLike(green,  "green",  upper.getStation("greenPt"));
    upper.addMassElementLike(orange, "orange", upper.getStation("orangePt"));
    upper.addMassElementLike(blue,   "blue1",  upper.getStation("pinnedBluePt"));
    upper.addMassElementLike(blue,   "blue2",  upper.getStation("otherBluePt"));


    // Create some joint frames at the appropriate stations (by default
    // they are aligned with the body frame).
    upper.addFrame("pinFrame",          upper.getStation("pinnedBluePt"));
    upper.addFrame("leftBallFrame",     upper.getStation("leftAttachPt"));
    upper.addFrame("rightBallFrame",    upper.getStation("rightAttachPt"));
    cout << "U=" << upper;

    upper.updStation("pinnedBluePt").
        place(upper.getFrame("pinFrame").getOrigin());
    cout << "after moving pinnedBluePt to pinFrame.O, U=" << upper;

/*
    // Now build the prototype for the lower bodies.

    RigidBody lower("L");

    // Create a parameter for the cylinder half-height, and a station at
    // the top center of the cylinder where we can hook the ball joint.
    lower.addRealParameter("halfHeight");
    lower.addStation("ballAttachPt");

    lower.updStation("ballAttachPt").setPlacement(
                       StationPlacement(lower.y(), lower.getRealParameter("halfHeight")));

    // Now instantiate a tube on the body prototype.
    MassElement& tube = lower.addMassElementLike(tubeProto, "tube");
    cout << "L=" << lower;
  
    // Place the center and axis, but leave the halfHeight parameter unresolved
    // because we want to control both with a single parameter.
    tube.placeFeature("center", lower.getOrigin());
    tube.placeFeature("axis", DirectionPlacement::minus(lower.getStation("ballAttachPt"),
                                                          lower.getOrigin()));
  
    ////////////////////////////////////////////
    // Create an articulated multibody system //
    ////////////////////////////////////////////
	
    Multibody mbs("example1");
    mbs.addGroundBody();
    Body &upperBody = mbs.addBodyLike(upper, "upper"),
         &leftLeg   = mbs.addBodyLike(lower, "left"),
         &rightLeg  = mbs.addBodyLike(lower, "right");

    // Create a single parameter of the multibody which can be used
    // to control the two halfHeights together.
    Parameter& hh = mbs.addParameter("halfHeight");
    leftLeg .updFeature("tube").updParameter("halfHeight").set(hh);
    rightLeg.updFeature("tube").updParameter("halfHeight").set(hh);

    mbs.addJoint(PinJoint, "base2ground", 
                 mbs.getGroundFrame(),                      //reference frame
                 upperBody);                                //moving frame
    mbs.addJoint(BallJoint, "leftHipJoint",
                 upperBody,                                 //reference frame
                 leftLeg.getFrame("upperAttachmentFrame")); //moving frame
    mbs.addJoint(BallJoint, "rightHipJoint",
                 upperBody,                                 //reference frame
                 rightLeg.getFrame("upperAttachmentFrame"));//moving frame

    std::cout << mbs << std::endl; //let’s see what we’ve got
    */
}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
