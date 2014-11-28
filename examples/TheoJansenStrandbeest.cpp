/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: TheoJansenStrandbeest                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */
#include "Simbody.h"
#include <iostream>
using namespace SimTK;
using std::cout; using std::endl; using std::cin;

/* This is a walking mechanism due to Theo Jansen. See strandbeest.com and
many YouTube videos. This is a full 3D simulation although each leg moves
only in 2D. The user controls the speed. I haven't yet figured out how to
steer one of these things though so it's not a very exciting drive!

Two different ways of building this model are demonstrated here based on the
compile-time flag below. One models all the links as they appear; the other
treats some links as massless resulting in halving the model size and better
than a 2X speedup.
*/

// Define this to use a simplified model that replaces some of the links
// with distance constraints. That reduces the model size by about half without
// changing the functionality at all.
//#define USE_MASSLESS_LINKS

// Undefine this to get more accurate CPU times.
#define ANIMATE

// Put local classes and definitions in the file-scope anonymous namespace.
namespace {

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input, mostly for speed control.
const int SpeedControlSlider = 1;
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, 
                     const Motion::Steady&  motor, 
                     Real                   interval); 
    void handleEvent(State& state, Real accuracy,
                     bool& shouldTerminate) const OVERRIDE_11;
private:
    Visualizer::InputSilo& m_silo;
    Motion::Steady         m_motor;
};

// Handy utility routine: given two vertices v1, v2 of a triangle in the x-y 
// plane and the lengths of the other two sides, find the location of the third 
// vertex, assuming v1-v2-v3 have counterclockwise ordering about the plane 
// normal. z coordinate is ignored on input and zero on output.
Vec3 findOtherVertex(const Vec3& v1, const Vec3& v2, Real s1, Real s2);

// Write interesting integrator info to stdout at end of simulation.
void dumpIntegratorStats(double startCPU, double startTime,
                         const Integrator& integ);

const Real LinkDepth = .01; // half depth of links
const Real LinkWidth = .02; // half width of links

const Real mu_s = 0.7;       // Friction coefficients.
const Real mu_d = 0.5;
const Real mu_v = 0;
const Real transitionVelocity = 1e-3; // slide->stick velocity

// Rubber for feet
const Real rubber_density = 1100.;  // kg/m^3
const Real rubber_young   = 0.01e9/10; // pascals (N/m)
const Real rubber_poisson = 0.5;    // ratio
const Real rubber_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
const Real rubber_dissipation = /*0.005*/10;

const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
                               mu_s,mu_d,mu_v);

// Concrete for ground
const Real concrete_density = 2300.;  // kg/m^3
const Real concrete_young   = 25e9;  // pascals (N/m)
const Real concrete_poisson = 0.15;    // ratio
const Real concrete_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
const Real concrete_dissipation = 0.005;

const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
                               mu_s,mu_d,mu_v);

// Add one leg to the given torso T with the leg frame's pose in T given.
void addOneLeg(Visualizer& viz, MobilizedBody& torso, const Transform& X_TL,
               MobilizedBody& crank, const Vec3& crankConnect);

const Rotation YtoX(-Pi/2,ZAxis); // some useful rotations
const Rotation YtoZ( Pi/2,XAxis);
}



//==============================================================================
//                                   MAIN
//==============================================================================
int main() {
    try { // catch errors if any
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    ContactTrackerSubsystem   tracker(system);
    CompliantContactSubsystem contact(system, tracker);
    contact.setTransitionVelocity(transitionVelocity);
    Force::Gravity(forces, matter, -YAxis, 9.81);

    // Set up visualization and ask for a frame every 1/30 second.
    Visualizer viz(system);
    viz.setShowSimTime(true); viz.setShowFrameRate(true);
    viz.addSlider("Speed", SpeedControlSlider, -10, 10, 0);
    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(silo);   
    #ifdef ANIMATE
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    #endif
    DecorativeText help("Any input to start; ESC to quit");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0),Vec3(0),help);
    matter.setShowDefaultGeometry(false);

    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.
    MobilizedBody& Ground = matter.updGround(); // Nicer name for Ground.
    Ground.updBody().addContactSurface(Transform(YtoX,Vec3(0)),
        ContactSurface(ContactGeometry::HalfSpace(),concrete));

    // Add some speed bumps.
    const int NBumps = 2; const Vec3 BumpShape(.8,0.2,2);
    for (int i=0; i < NBumps; ++i) {
        const Real x = -2*(i+1);
        Ground.updBody().addContactSurface(Vec3(x,0,0),
            ContactSurface(ContactGeometry::Ellipsoid(BumpShape), rubber));
        Ground.updBody().addDecoration(Vec3(x,0,0),
            DecorativeEllipsoid(BumpShape).setColor(Gray).setResolution(3));
    }

    // TORSO
    const Real TorsoHeight = 1.1;
    const Vec3 torsoHDims(1,.08,.8);
    const Real torsoVolume = 8*torsoHDims[0]*torsoHDims[1]*torsoHDims[2];
    const Real torsoMass = torsoVolume*rubber_density/10;
    const Vec3 torsoCOM(0,-.75,0); // put it low for stability
    Body::Rigid torsoInfo(MassProperties(torsoMass,torsoCOM,
        UnitInertia::brick(torsoHDims).shiftFromCentroidInPlace(-torsoCOM)));
    torsoInfo.addDecoration(Vec3(0),
        DecorativeBrick(torsoHDims).setColor(Cyan));

    // CRANK
    const Vec3 crankCenter(0,0,0); // in torso frame
    const Vec3 crankOffset(0,0,torsoHDims[2]+LinkDepth); // left/right offset
    const Real MLen=15/100.; // crank radius
    Body::Rigid crankInfo(MassProperties(.1,Vec3(0),
                            UnitInertia::cylinderAlongZ(MLen, LinkDepth)));
    crankInfo.addDecoration(Vec3(0),
        DecorativeBrick(Vec3(LinkWidth,LinkWidth,torsoHDims[2]))
        .setColor(Black));
    const Vec3 CrankConnect(MLen,0,0); // in crank frame

    // Add the torso and crank mobilized bodies.
    MobilizedBody::Free torso(Ground,Vec3(0,TorsoHeight,0), torsoInfo,Vec3(0));
    MobilizedBody::Pin crank(torso,crankCenter, crankInfo, Vec3(0));
    
    // Add the legs.
    for (int i=-1; i<=1; ++i) {
        const Vec3 offset = crankCenter + i*crankOffset;
        const Vec3 linkSpace(0,0,2*LinkDepth);
        const Rotation R_CP(i*2*Pi/3,ZAxis);
        // Add crank bars for looks.
        crank.addBodyDecoration(
            Transform(R_CP, offset+1.5*MLen/2*R_CP.x()+(i==0?linkSpace:Vec3(0))),
            DecorativeBrick(Vec3(1.5*MLen/2,LinkWidth,LinkDepth))
                        .setColor(Yellow));

        addOneLeg(viz, torso, offset + i*linkSpace, 
                  crank, R_CP*CrankConnect);
        addOneLeg(viz, torso, Transform(Rotation(Pi,YAxis), offset + i*linkSpace), 
                  crank, R_CP*CrankConnect);
    }

    // Add speed control.
    Motion::Steady motor(crank, 0); // User controls speed.
    system.addEventHandler
       (new UserInputHandler(*silo, motor, Real(0.1))); //check input every 100ms
  
    // Initialize the system and state.    
    State state = system.realizeTopology();
    system.realize(state);
    printf("Theo Jansen Strandbeest in 3D:\n");
    printf("%d bodies, %d mobilities, -%d constraint equations -%d motions\n",
        matter.getNumBodies(), state.getNU(), state.getNMultipliers(), 
        matter.getMotionMultipliers(state).size());

    viz.report(state);
    printf("Hit any key to assemble ...");
    silo->waitForAnyUserInput(); silo->clear();
    Assembler(system).assemble(state);
    printf("ASSEMBLED\n");

    // Simulate.
    SemiExplicitEuler2Integrator integ(system);
    integ.setAccuracy(0.1);
    integ.setConstraintTolerance(.001);
    integ.setMaximumStepSize(1./60);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    viz.report(ts.getState());
    printf("Hit ENTER to simulate ... (ESC to quit)\n");
    silo->waitForAnyUserInput(); silo->clear();

    const double startCPU  = cpuTime(), startTime = realTime();
    ts.stepTo(Infinity); // RUN
    dumpIntegratorStats(startCPU, startTime, integ);

    } catch (const std::exception& e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

namespace {

//==============================================================================
//                               ADD ONE LEG
//==============================================================================
// This function can build a leg two different ways, depending on whether
// USE_MASSLESS_LINKS is defined:
// (1) A straightforward model where all links have mass, producing 6 mobilities
//     and requiring 6 constraints per leg. 
// (2) Only the shoulder and foot bodies have mass, so the model requires only 3
//     mobilities and 3 distance constraints instead, for much improved 
//     performance.    
void addOneLeg(Visualizer& viz, MobilizedBody& torso, const Transform& X_TL,
               MobilizedBody& crank, const Vec3& crankConnect)
{
    // Segment lettering is from TJ's drawing here: 
    //        http://www.strandbeest.com/beests_leg.php
    // These dimensions are what TJ calls the "holy numbers". These are scaled
    // from cm to m so that the mechanism is about 1m tall.
    const Real ALen=38/100., LLen=7.8/100.; 
    const Real BLen=41.5/100., ELen=55.8/100., DLen=40.1/100.; // shoulder sides
    const Real ILen=49/100., HLen=65.7/100., GLen=36.7/100.;   // foot sides

    const Real CLen = 39.3/100.; // Link lengths, TJ's lettering
    const Real JLen = 50/100.;
    const Real FLen = 39.4/100.;
    const Real KLen = 61.9/100.;

    // The pivot point is where the shoulder is pinned to the torso.
    const Vec3 pivotPt = X_TL*Vec3(-ALen,-LLen,0);  // to torso frame

    // SHOULDER (a triangular body)
    // Put the shoulder origin at the pivot, one side in the +Y.
    Body::Rigid shoulderInfo(MassProperties(1,Vec3(0),UnitInertia(1)));
    const Vec3 shoulderPivot(0,0,0);
    const Vec3 shoulderUpper(0,BLen,0);
    const Vec3 shoulderSide = findOtherVertex(shoulderPivot,shoulderUpper,
                                              DLen,ELen);
    const UnitVec3 BDir(shoulderUpper-shoulderPivot);
    const UnitVec3 DDir(shoulderSide-shoulderPivot);
    const UnitVec3 EDir(shoulderUpper-shoulderSide);

    Rotation R_SB(BDir, XAxis, Vec3(0,0,1), ZAxis);
    Rotation R_SD(DDir, XAxis, Vec3(0,0,1), ZAxis);
    Rotation R_SE(EDir, XAxis, Vec3(0,0,1), ZAxis);
    shoulderInfo.addDecoration(Transform(R_SB, shoulderPivot+BLen/2*BDir),
        DecorativeBrick(Vec3(BLen/2,LinkWidth,LinkDepth)).setColor(Orange));
    shoulderInfo.addDecoration(Transform(R_SD,shoulderPivot+DLen/2*DDir),
        DecorativeBrick(Vec3(DLen/2,LinkWidth,LinkDepth)).setColor(Orange));
    shoulderInfo.addDecoration(Transform(R_SE, shoulderSide+ELen/2*EDir),
        DecorativeBrick(Vec3(ELen/2,LinkWidth,LinkDepth)).setColor(Orange));

    // FOOT (a triangular body)
    // Put the foot origin at the pivot, one side in the -Y.
    Body::Rigid footInfo(MassProperties(10,Vec3(0),UnitInertia(1)));
    const Vec3 footPivot(0,0,0);
    const Vec3 footLower(0,-ILen,0);
    const Vec3 footSide = findOtherVertex(footLower,footPivot,HLen,GLen);

    const UnitVec3 IDir(footLower-footPivot);
    const UnitVec3 GDir(footSide-footPivot);
    const UnitVec3 HDir(footLower-footSide);

    Rotation R_SI(IDir, XAxis, Vec3(0,0,1), ZAxis);
    Rotation R_SG(GDir, XAxis, Vec3(0,0,1), ZAxis);
    Rotation R_SH(HDir, XAxis, Vec3(0,0,1), ZAxis);
    footInfo.addDecoration(Transform(R_SI, footPivot+ILen/2*IDir),
        DecorativeBrick(Vec3(ILen/2,LinkWidth,LinkDepth)).setColor(Orange));
    footInfo.addDecoration(Transform(R_SG, footPivot+GLen/2*GDir),
        DecorativeBrick(Vec3(GLen/2,LinkWidth,LinkDepth)).setColor(Orange));
    footInfo.addDecoration(Transform(R_SH, footSide+HLen/2*HDir),
        DecorativeBrick(Vec3(HLen/2,LinkWidth,LinkDepth)).setColor(Orange));

    const Real FootRad = .05;
    footInfo.addContactSurface(footLower-Vec3(0,FootRad/2,0), 
        ContactSurface(ContactGeometry::Sphere(FootRad), rubber));
    footInfo.addDecoration(footLower-Vec3(0,FootRad/2,0), 
        DecorativeSphere(FootRad).setColor(Orange));

    // LINKS

    // Link C connects foot to pivot point. Start aligned with Y.
    const Vec3 CDims(LinkWidth,CLen/2,LinkDepth);
    Body::Rigid linkCInfo(MassProperties(.1,Vec3(0), UnitInertia::brick(CDims)));
    linkCInfo.addDecoration(Vec3(0), DecorativeBrick(CDims).setColor(Orange));

    // Link J connects upper shoulder to crank. Start aligned with X.
    const Vec3 JDims(JLen/2,LinkWidth,LinkDepth);
    Body::Rigid linkJInfo(MassProperties(.1,Vec3(0), UnitInertia::brick(JDims)));
    linkJInfo.addDecoration(Vec3(0), DecorativeBrick(JDims).setColor(Orange));

    // Link F connects shoulder to foot. Start aligned with Y.
    const Vec3 FDims(LinkWidth,FLen/2,LinkDepth);
    Body::Rigid linkFInfo(MassProperties(.1,Vec3(0), UnitInertia::brick(FDims)));
    linkFInfo.addDecoration(Vec3(0), DecorativeBrick(FDims).setColor(Orange));


    // Link K connects foot to crank. Start aligned with X.
    const Vec3 KDims(KLen/2,LinkWidth,LinkDepth);
    Body::Rigid linkKInfo(MassProperties(.1,Vec3(0), UnitInertia::brick(KDims)));
    linkKInfo.addDecoration(Vec3(0), DecorativeBrick(KDims).setColor(Orange));

    // Create the tree of mobilized bodies.
    MobilizedBody::Pin shoulder(torso, Transform(X_TL.R(), pivotPt), 
                                shoulderInfo, shoulderPivot);
    MobilizedBody::Pin linkC(torso,    Transform(X_TL.R(), pivotPt),  
                             linkCInfo, Vec3(0,CLen/2,0));
    MobilizedBody::Pin foot( linkC,     Vec3(0,-CLen/2,0), 
                             footInfo,  footPivot);

    #ifdef USE_MASSLESS_LINKS
    Vec3 crankAttach(crankConnect[0],crankConnect[1],X_TL.p()[2]);
    DecorativeLine line; line.setColor(Gray).setLineThickness(3);
    Constraint::Rod linkF(shoulder, shoulderSide, foot, footSide, FLen);
    viz.addRubberBandLine(shoulder, shoulderSide, foot, footSide, line);
    Constraint::Rod linkJ(shoulder, shoulderUpper, crank, crankAttach, JLen);
    viz.addRubberBandLine(shoulder, shoulderUpper, crank, crankAttach, line);
    Constraint::Rod linkK(foot, footPivot, crank, crankAttach, KLen);
    viz.addRubberBandLine(foot, footPivot, crank, crankAttach, line);
    #else
    MobilizedBody::Pin linkJ(shoulder,  shoulderUpper,
                             linkJInfo, Vec3(-JLen/2,0,0));
    MobilizedBody::Pin linkF(shoulder,  shoulderSide,
                             linkFInfo, Vec3(0,FLen/2,0));
    MobilizedBody::Pin linkK(foot,      footPivot, 
                             linkKInfo, Vec3(-KLen/2,0,0));
    linkJ.setDefaultAngle(-Pi/6); // set default angles to guide assembly
    linkK.setDefaultAngle( Pi/4);


    // Add 2d pin joint constraints (each pair of point-in-planes is a pin).
    Constraint::PointInPlane f2footx(foot, XAxis, footSide[0],
                                     linkF, Vec3(0,-FLen/2,0));
    Constraint::PointInPlane f2footy(foot, YAxis, footSide[1],
                                     linkF, Vec3(0,-FLen/2,0));

    Constraint::PointInPlane k2crankx(crank, XAxis, crankConnect[0],
                                      linkK, Vec3(KLen/2,0,0));
    Constraint::PointInPlane k2cranky(crank, YAxis, crankConnect[1],
                                      linkK, Vec3(KLen/2,0,0));

    Constraint::PointInPlane j2crankx(crank, XAxis, crankConnect[0],
                                      linkJ, Vec3(JLen/2,0,0));
    Constraint::PointInPlane j2cranky(crank, YAxis, crankConnect[1],
                                      linkJ, Vec3(JLen/2,0,0));
    #endif
}


//==============================================================================
//                           USER INPUT HANDLER
//==============================================================================
UserInputHandler::UserInputHandler(Visualizer::InputSilo& silo, 
                                   const Motion::Steady&  motor, 
                                   Real                   interval) 
:   PeriodicEventHandler(interval), m_silo(silo), m_motor(motor) {}

void UserInputHandler::handleEvent(State& state, Real accuracy,
                                   bool& shouldTerminate) const  {
    while (m_silo.isAnyUserInput()) {
        unsigned key, modifiers;
        while (m_silo.takeKeyHit(key,modifiers))
            if (key == Visualizer::InputListener::KeyEsc) {
                shouldTerminate = true;
                m_silo.clear();
                return;
            }

        int whichSlider; Real sliderValue;
        while (m_silo.takeSliderMove(whichSlider, sliderValue))
            if (whichSlider == SpeedControlSlider)
                m_motor.setRate(state, sliderValue);
    }  
}


//==============================================================================
//                           FIND OTHER VERTEX
//==============================================================================
// Given two vertices v1, v2 of a triangle in the x-y plane and the lengths of 
// the other two sides, find the location of the third vertex, assuming 
// v1-v2-v3 have counterclockwise ordering about the plane normal.
//
//                           v2
//                            * 
//                             \
//                           .  \ s0
//                               \
//                          s2    \
//                               a * v1
//                          .    .
//                             s1
//                           .
//                         *
//                         v3
//
// Our strategy will be to find the angle a and rotate the unit vector 
// v=(v2-v1)/|v2-v1| ccw by a, giving unit vector w along v1v3. Then v3=v1+s1*w.

// Ignore z component of vectors -- we're working in x-y plane.
Vec3 findOtherVertex(const Vec3& v1, const Vec3& v2, 
                     Real s1, Real s2)
{
    const Real s0 = (v2-v1).norm();
    const Real ca = (s0*s0 + s1*s1 - s2*s2) / (2*s0*s1);  // cos(a)
    const Real a = std::acos(ca);
    const Real sa = std::sin(a);
    const Mat22 R(ca, -sa,  // 2d rotation matrix
                  sa,  ca);
    const Vec2 v = (v2.drop1(2)-v1.drop1(2))/s0;
    const Vec2 w = R*v;
    const Vec2 v3 = v1.drop1(2) + s1*w;
    return v3.append1(0);
}

//==============================================================================
//                        DUMP INTEGRATOR STATS
//==============================================================================
void dumpIntegratorStats(double startCPU, double startTime, 
                         const Integrator& integ) {
    std::cout << "DONE: Simulated " << integ.getTime() << " seconds in " <<
        realTime()-startTime << " elapsed s, CPU="<< cpuTime()-startCPU << "s\n";
    #ifdef ANIMATE
    printf("***CAUTION: CPU time not accurate when animation is enabled.\n");
    #endif

    const int evals = integ.getNumRealizations();
    std::cout << "\nUsed "  << integ.getNumStepsTaken() << " steps, avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms " 
        << (1000*integ.getTime())/evals << "ms/eval\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n",  integ.getNumStepsTaken(), 
                                          integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n",     integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), 
                                          integ.getNumProjections());
}
}