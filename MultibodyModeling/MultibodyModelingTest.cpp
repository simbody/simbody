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
    f.realize(/*State,*/ Stage::Startup);
    //cout << f;

    f.addRealParameter("one").place(1.);
    f.addRealParameter("negone");
    f.addRealMeasure("zzz").place(f["one"]+f["negone"]);
    f.addRealMeasure("ttt").place(f["one"]+abs(f["negone"]));
    f.addRealMeasure("x2y").place(angle(f.x(),f.y()));
    f.addVec3Measure("YCrossX").place(cross(f.y(),f.x()));

    f.realize(/*State,*/ Stage::Startup);

    try {
        cout << "one=" << f["one"].getValue() << endl;
        cout << "zzz=" << f["zzz"].getValue() << endl;
        cout << "ttt=" << f["ttt"].getValue() << endl;
    }
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    f["negone"].place(-1);
    Placement p;
    p = 9.*(f["one"]+abs(f["negone"]));
    f.addRealMeasure("ttt9").place(p);
    f.realize(/*State,*/ Stage::Startup);

    try{cout << "p.getValue()=" << p.getValue() << endl;}
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    try {
        cout << "one=" << f["one"].getValue() << endl;
        cout << "zzz=" << f["zzz"].getValue() << endl;
        cout << "ttt=" << f["ttt"].getValue() << endl;
        cout << "ttt9=" << f["ttt9"].getValue() << endl;
        cout << "x2y=" << f["x2y"].getValue() << endl;
        cout << "y%x=" << f["YCrossX"].getValue() << endl;
    }
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    cout << f;



    ////////////////////////////////////
    // Create mass element prototypes //
    ////////////////////////////////////

    // Create point masses with constant values for their mass parameter.
    // These will require a "station" placement.
    PointMassElement    blue  ("blue",  2.5), 
                        orange("orange",1.),
                        green ("green", 0.1);
    //blue.place(Vec3(1,2,3)); // a pointless self-placement

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
    //cout << "U=" << upper;

    cout << "U's subfeatures:" << endl;
    for (int i=0; i<upper.getNSubsystems(); ++i)
        cout << upper[i].getFullName() << endl;

    // Now build the prototype for the lower bodies.


    RigidBody lower("L");

    // Create a parameter for the cylinder half-height, and a station at
    // the top center of the cylinder where we can hook the ball joint.
    lower.addRealParameter("halfHeight");
    lower.addStation      ("ballAttachPt",         lower["halfHeight"] * lower.y());
    lower.addFrame        ("upperAttachmentFrame", lower["ballAttachPt"]);

    // Now instantiate a tube on the body prototype.
    MassElement& tube = lower.addMassElementLike(tubeProto, "tube");

    lower.addRealMeasure  ("hh9", lower["halfHeight"] + 9. / lower["tube/mass"]);
  
    // Place the center and axis, but leave the halfHeight parameter unresolved
    // because we want to control both with a single parameter.
    tube["center"].place(lower.getOrigin());

    //TODO: this should require an explicit normalize()
    tube["axis"].place(lower["ballAttachPt"] - lower.getOrigin());

    tube["halfLength"].place(lower["halfHeight"]);

    //cout << "L=" << lower; 

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

    mbs.realize(/*State,*/ Stage::Startup);

    try {cout << "left/tube/axis=" << mbs["left/tube/axis"].getValue() << endl;}
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    mbs["halfHeight"].place(13.111);
    mbs.realize(/*State,*/ Stage::Startup);

    try {cout << "left/tube/axis=" << mbs["left/tube/axis"].getValue() << endl;}
    catch(const Exception::Base& e) {std::cout << e.getMessage() << std::endl;}

    mbs.checkSubsystemConsistency(0,-1,mbs);    
    //std::cout << "***MULTIBODY SYSTEM***" << std::endl;
   // std::cout << mbs << std::endl; //let’s see what we’ve got
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
            for (int j=0; j<mbs[i].getNSubsystems(); ++j) {
                if (MassElement::isInstanceOf(mbs[i][j])) {
                    Real mass = MassElement::downcast(mbs[i][j]).getMassMeasure().getValue();
                    cout << "  massElt " << mbs[i][j].getName() << " " << mass << endl;
                    tmass += mass;
                }
            }
            cout << "... total mass=" << tmass << endl;
        }

    //FeatureList hasMass = mbs["upper"].select(MassElement::Selector());
    //FeatureList masses = hasMass.getSubfeature("mass"); // elementwise indexing
    //FeatureList centroids = hasMass.getSubfeature("centroid");
    //Real upperMass = sum(masses);
    //Vec3 upperCOM  = sum(prod(masses,centroids))/hasMass.size();

    // Any leftover parameters need external placements. We'll make a RuntimeFeature
    // to hold them.
//    RuntimeFeature rt("rt");
 //   rt.addFeatureLike(mbs, "mbs");
//    rt["mbs/halfHeight"].place(5.);

//    rt.realize();
//    cout << "upper mass=" << rt["mbs/upper/massMeasure"].getValue() << endl;


    ///////////////////////////
    // Build a RigidBodyTree //
    ///////////////////////////

}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
