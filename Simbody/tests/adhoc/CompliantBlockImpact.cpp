/* -------------------------------------------------------------------------- *
 *               Simbody(tm) Adhoc test: Compliant Block Impact               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Thomas Uchida                                    *
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

/* This adhoc test is for use in comparing compliant contact impact with
higher-performance approximations such as the PLUS method. The system is a
single block with contact spheres at each corner. The intent here is
to extract a lot of detailed information from the simulation during the
evolution of the impact so we'll run it slowly.
*/

#include "Simbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
using std::cout; using std::endl;

using namespace SimTK;

//#define USE_SHERM_PARAMETERS
#define USE_TOM_PARAMETERS
static const bool GenerateAnimation = true;


Array_<State> saveEm;

static const Real TimeScale = 1;
static const Real FrameRate = GenerateAnimation ? 6000 : 60;
static const Real ForceScale = .25;
static const Real MomentScale = .5;
const Vec3 BrickColor    = Blue;
const Vec3 SphereColor   = Red;

#ifdef USE_SHERM_PARAMETERS
    static const Real ReportInterval = TimeScale/FrameRate;
    static const Real VizReportInterval = .001;
#endif
#ifdef USE_TOM_PARAMETERS
    static const Real ReportInterval = GenerateAnimation
                                       ? TimeScale/FrameRate*10. : 1.e-4;
    static const Real VizReportInterval = GenerateAnimation ? 0.00005 : 1.e-4;
#endif


//==============================================================================
//                            FORCE ARROW GENERATOR
//==============================================================================
// Draws ephemeral geometry that is displayed for only one frame. Also prints
// the velocity of the first contact point to the console for analysis.
class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont,
                        const MobilizedBody& brick) 
    :   m_mbs(system), m_compliant(complCont), m_brick(brick),
        m_inContact(false), m_hasCompressionEnded(false) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) override {
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_mbs.realize(state, Stage::Velocity);
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const Real TextScale = 3./16.;
        m_mbs.realize(state, Stage::Dynamics);
        const Real KE=m_mbs.calcKineticEnergy(state), E=m_mbs.calcEnergy(state);
        const Real diss=m_compliant.getDissipatedEnergy(state);
        const Real PE=m_mbs.calcPotentialEnergy(state);
        DecorativeText txt; 
        //txt.setIsScreenText(true);
        //txt.setScale(TextScale);
        //txt.setColor(Orange);
        //txt.setText("KE/Diss/E: " + String(KE, "%.6f") + String(diss, "/%.6f")
        //            + String(E+diss, "/%.6f") );
        //const Vector qBrick = m_brick.getQAsVector(state);
        //geometry.push_back(DecorativeText(txt)
        //                   .setTransform(Vec3(qBrick[4],0,1)));
        //geometry.push_back(DecorativeText(txt).setTransform(Vec3(0,0,1)));

        txt.setColor(Black).setScale(0.8);
        txt.setText("Ed=" + String(diss, "%0.3f") );
        geometry.push_back(DecorativeText(txt).setTransform(Vec3(-6.75,-15,6.75)));
        txt.setText("Et=" + String(E+diss, "%0.3f") );
        geometry.push_back(DecorativeText(txt).setTransform(Vec3(-6.95,-15,5.75)));

        int nContPts = 0;
        const int ncont = m_compliant.getNumContactForces(state);

        if (ncont) {
            if (!m_inContact) {
                printf("\n\n--------- NEW CONTACT @%g ----------\n", 
                    state.getTime());
                m_inContact = true;
            }
        } else if (m_inContact) {
            printf("-------- END OF CONTACT --------\n\n");
            m_inContact = false;
            m_velLines.clear();
        }


        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            const Vec3&         pt    = force.getContactPoint();
            const Vec3& frc = force.getForceOnSurface2()[1];
            const Vec3& mom = force.getForceOnSurface2()[0];
            Real  frcMag = frc.norm(), momMag=mom.norm();
            const UnitVec3 frcDir(frc/frcMag, true);
            const UnitVec3 momDir(mom/momMag, true);
            const Vec3 offs = .1*frcDir; // shift up to clear ground
            //int frcThickness = 2, momThickness = 2;
            //Real frcScale = ForceScale, momScale = ForceScale;
            //while (frcMag > /*10*/1000000)
            //    frcThickness++, frcScale /= 10, frcMag /= 10;
            //while (momMag > /*10*/1000000)
            //    momThickness++, momScale /= 10, momMag /= 10;
            //geometry.push_back(DecorativePoint(pt)
            //                   .setScale(5).setColor(Yellow));
            //DecorativeLine frcLine(pt,      pt + std::log10(frcMag)*frcDir);
            //DecorativeLine momLine(pt+offs, pt+offs + std::log10(momMag)*momDir);
            //frcLine.setColor(Black);
            //momLine.setColor(Purple);
            //frcLine.setLineThickness(frcThickness);
            //momLine.setLineThickness(momThickness);
            //geometry.push_back(frcLine);
            //geometry.push_back(momLine);


            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                ++nContPts;
                const ContactDetail& detail = patch.getContactDetail(i);
                const Vec3& pt = detail.getContactPoint();
                //geometry.push_back(DecorativePoint(pt).setColor(Purple));

                const Vec3& force = detail.getForceOnSurface2();
                const Real forceMag = force.norm();
                //const UnitVec3 forceDir(force/forceMag, true);
                //DecorativeLine frcLine(pt, pt+std::log10(forceMag)*forceDir);
                //frcLine.setColor(Black);
                //geometry.push_back(frcLine);

                // Draw a line that extends from the contact
                // point in the direction of the slip velocity.
                const Vec3 v     = detail.getSlipVelocity();
                const Real vMag  = std::max(0., std::log10(v.norm()*1.e3));
                const Vec3 vDraw = v.normalize() * vMag;
                const Vec3 pt0   = Vec3(pt[0], pt[1], 5.e-3);
                const Vec3 pt1   = Vec3(pt[0]+2.*vDraw[0], pt[1]+2.*vDraw[1],
                                        5.e-3);
                Real colorFactor = clamp(0.0, m_velLines.size() / 32., 1.0);
                DecorativeLine slip(pt0, pt1);
                slip.setLineThickness(3)
                    .setColor(Vec3(1-colorFactor,0.,colorFactor));

                // Store line for displaying in subsequent frames, but only if
                // in compression.
                if (m_velLines.size() > 0 && v[ZAxis] > 0 && !m_hasCompressionEnded) {
                    m_hasCompressionEnded = true;
                    cout << m_velLines.size() << " lines drawn." << endl;
                }

                const bool inCompression = (v[ZAxis] <= 0)
                                           && !m_hasCompressionEnded;
                if (inCompression && vMag>0)
                    m_velLines.push_back(slip);

                for (int k=0; k<(int)m_velLines.size(); ++k)
                    geometry.push_back(m_velLines[k]);

                if (i==0 && inCompression && !GenerateAnimation) // REPORT ONLY FIRST CONTACT
                    printf("%8.4f %8.4f %8.4f\n", state.getTime(), v[0], v[1]);
            }
        }
        //txt.setText(String("Num contact points: ") + String(nContPts));
        //geometry.push_back(DecorativeText(txt)
        //                   .setTransform(Vec3(state.getQ()[4],0,.75)));
        txt.setText("t=" + String(state.getTime(), "%0.3f") + "s");
        txt.setColor(Black).setScale(0.8);
        geometry.push_back(DecorativeText(txt).setTransform(Vec3(-7.5,-15,7.75)));

        if (nContPts>0) {
            txt.setText(String("Contacting"));
            geometry.push_back(DecorativeText(txt).setTransform(Vec3(9,-15,7.75)));
        }

    }
private:
    const MultibodySystem&              m_mbs;
    const CompliantContactSubsystem&    m_compliant;
    const MobilizedBody&                m_brick;
    bool                                m_inContact;
    Array_<DecorativeLine>              m_velLines;
    bool                                m_hasCompressionEnded;
};


// These are the item numbers for the entries on the Run menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem=2, QuitItem=3;

//==============================================================================
//                           USER INPUT HANDLER
//==============================================================================
// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
// This one does nothing but look for the Run->Quit selection.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, Real interval) 
    :   PeriodicEventHandler(interval), m_silo(silo) {}

    virtual void handleEvent(State& state, Real accuracy, 
                             bool& shouldTerminate) const override 
    {
        int menuId, item;
        if (m_silo.takeMenuPick(menuId, item) && menuId==RunMenuId && item==QuitItem)
            shouldTerminate = true;
    }

private:
    Visualizer::InputSilo& m_silo;
};


//==============================================================================
//                                 BODY WATCHER
//==============================================================================
// Prior to rendering each frame, point the camera at the given body's origin.
// Adapted from TimsBox.cpp.
class BodyWatcher : public Visualizer::FrameController {
public:
    explicit BodyWatcher(const MobilizedBody& body) : m_body(body) {}

    void generateControls(const Visualizer&             viz,
                          const State&                  state,
                          Array_< DecorativeGeometry >& geometry) override
    {
        const Vec3 Bo = m_body.getBodyOriginLocation(state);
        //const Vec3 p_GC = Bo + Vec3(0,4,2);
        //const Rotation R1(-SimTK::Pi/3, XAxis);
        const Vec3 p_GC = Bo + Vec3(-1.,4.,0.5);
        const Rotation R1(-1.6, XAxis);
        const Rotation R2(SimTK::Pi, ZAxis);
        viz.setCameraTransform(Transform(R1*R2, p_GC));
    }
private:
    const MobilizedBody m_body;
};

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    system.setUpDirection(ZAxis);
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, -ZAxis, 9.81);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTrackDissipatedEnergy(true);
    contactForces.setTransitionVelocity(1e-3);

    const Vec3 hdim(.2,.3,.4); // Brick half dimensions
    const Real rad = .1;    // Contact sphere radius
    const Real brickMass = 2;

    #ifdef USE_SHERM_PARAMETERS
        const Real mu_d =.3; // dynamic friction
        const Real mu_s =.3; // static friction
        const Real mu_v = 0;  // viscous friction (1/v)
        const Real dissipation = .1;
        const Real fK = 1e6; // stiffness in pascals
        const Real simDuration = 5.;
    #endif
    #ifdef USE_TOM_PARAMETERS
        const Real mu_d =.3; // dynamic friction
        const Real mu_s =.3; // static friction
        const Real mu_v = 0;  // viscous friction (1/v)
        const Real dissipation = .1756;  //Second impact at 0.685 s.
        const Real fK = 1e6; // stiffness in pascals
        const Real simDuration = 0.5; //3.0; //0.8;
    #endif

    const ContactMaterial material(fK,dissipation,mu_s,mu_d,mu_v);

    // Halfspace floor
    const Rotation R_xdown(Pi/2,YAxis);
    matter.Ground().updBody().addContactSurface(
        Transform(R_xdown, Vec3(0,0,0)),
        ContactSurface(ContactGeometry::HalfSpace(), material));

    Body::Rigid brickBody(MassProperties(brickMass, Vec3(0), 
                            UnitInertia::brick(hdim)));
    brickBody.addDecoration(Transform(), 
        DecorativeBrick(hdim).setColor(BrickColor).setOpacity(.7));

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(hdim);
        brickBody.addContactSurface(pt,
            ContactSurface(ContactGeometry::Sphere(rad), material));
        brickBody.addDecoration(pt,
            DecorativeSphere(rad).setColor(SphereColor));
    }

    MobilizedBody::Free brick(matter.Ground(), Transform(),
                              brickBody, Transform());

    Visualizer viz(system);
    viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces,
                                                       brick));
    //viz.addFrameController(new BodyWatcher(brick));
    viz.addFrameController(new BodyWatcher(matter.Ground()));
    //viz.setShowSimTime(true);
    //viz.setShowFrameNumber(true);
    viz.setDesiredFrameRate(FrameRate);
    //viz.setShowFrameRate(true);

    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(silo);
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Go", GoItem));
    runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
    runMenuItems.push_back(std::make_pair("Quit", QuitItem));
    viz.addMenu("Run", RunMenuId, runMenuItems);

    Array_<std::pair<String,int> > helpMenuItems;
    helpMenuItems.push_back(std::make_pair("TBD - Sorry!", 1));
    viz.addMenu("Help", HelpMenuId, helpMenuItems);

    // Check for a Run->Quit menu pick every 1/4 second.
    //system.adoptEventHandler(new UserInputHandler(*silo, .25));

    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    // SET INITIAL CONDITIONS
    #ifdef USE_SHERM_PARAMETERS
        brick.setQToFitTranslation(state, Vec3(0,2,.8));
        brick.setQToFitRotation(state, Rotation(BodyRotationSequence, 
                                                Pi/4, XAxis, Pi/6, YAxis));
        brick.setUToFitLinearVelocity(state, Vec3(-5,0,0));
    #endif
    #ifdef USE_TOM_PARAMETERS
        Vector initQ = Vector(Vec7(1,0,0,0, 0,1,0.8));
        initQ(0,4) = Vector(Quaternion(Rotation(SimTK::Pi/4, XAxis) *
                                       Rotation(SimTK::Pi/6, YAxis))
                            .asVec4());
        Vector initU = Vector(Vec6(0,0,0, 0,0,6));
        initQ[6] = 1.5;
        initU[5] = -3.96;  //First impact at 0.181 s.
        initU[3] = -5.0;
        state.setQ(initQ);
        state.setU(initU);
    #endif

    saveEm.reserve(10000);

    viz.report(state);
    printf("Default state\n");
    cout << "t=" << state.getTime() 
         << " q=" << brick.getQAsVector(state)  
         << " u=" << brick.getUAsVector(state)  
         << endl;

    cout << "\nChoose 'Go' from Run menu to simulate:\n";
    int menuId, item;
    do { silo->waitForMenuPick(menuId, item);
         if (menuId != RunMenuId || item != GoItem) 
             cout << "\aDude ... follow instructions!\n";
    } while (menuId != RunMenuId || item != GoItem);

   
    // Simulate it.

    // The system as parameterized is very stiff (mostly due to friction)
    // and thus runs best with CPodes which is extremely stable for
    // stiff problems. To get reasonable performance out of the explicit
    // integrators (like the RKs) you'll have to run at a very loose
    // accuracy like 0.1, or reduce the friction coefficients and
    // maybe the stiffnesses.

    //SemiExplicitEuler2Integrator integ(system);
    //CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
    RungeKuttaMersonIntegrator integ(system);
    integ.setReturnEveryInternalStep(true);
    integ.setAllowInterpolation(false);
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    //integ.setAccuracy(1e-3); // minimum for CPodes
    integ.setAccuracy(1e-5);
    //integ.setAccuracy(.01);

    integ.initialize(state);
    double cpuStart = cpuTime();
    double realStart = realTime();
    Real lastReport = -Infinity;
    while (integ.getTime() < simDuration) {
        // Advance time by no more than ReportInterval. Might require multiple 
        // internal steps.
        integ.stepBy(ReportInterval);

        if (integ.getTime() >= lastReport + VizReportInterval) {
            // The state being used by the integrator.
            const State& s = integ.getState();
            viz.report(s);
            saveEm.push_back(s); // save state for playback
            lastReport = s.getTime();
        }
    }

    const double timeInSec = realTime() - realStart;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s elapsed for " << integ.getTime() << "s sim (avg step=" 
        << (1000*integ.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*integ.getTime())/evals << "ms/eval\n";
    cout << "  CPU time was " << cpuTime() - cpuStart << "s\n";

    printf("Using Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

    viz.dumpStats(std::cout);

    // Add as slider to control playback speed.
    viz.addSlider("Speed", 1, 0, 2, 1);
    viz.setMode(Visualizer::PassThrough);

    silo->clear(); // forget earlier input
    double speed = 1; // will change if slider moves
    while(true) {
        cout << "Choose Run/Replay to see that again ...\n";

        int menuId, item;
        silo->waitForMenuPick(menuId, item);


        if (menuId != RunMenuId) {
            cout << "\aUse the Run menu!\n";
            continue;
        }

        if (item == QuitItem)
            break;
        if (item != ReplayItem) {
            cout << "\aHuh? Try again.\n";
            continue;
        }

        for (double i=0; i < (int)saveEm.size(); i += speed ) {
            int slider; Real newValue;
            if (silo->takeSliderMove(slider,newValue)) {
                speed = newValue;
            }
            viz.report(saveEm[(int)i]);
        }
    }

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}

