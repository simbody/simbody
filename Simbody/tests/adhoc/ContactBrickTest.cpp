/* -------------------------------------------------------------------------- *
 *                Simbody(tm) Adhoc test: Contact Brick Test                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

/* This adhoc test is for playing with Brick contact elements.
*/

#include "Simbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real TimeScale = 1;
static const Real FrameRate = 60;
static const Real ReportInterval = TimeScale/FrameRate;
static const Real ForceScale = .25;
static const Real MomentScale = .5;


class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont)
    :   m_mbs(system), m_compliant(complCont) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) {
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_mbs.realize(state, Stage::Velocity);
        const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
        const Real TextScale = m_mbs.getDefaultLengthScale()/10; // was .1
        m_mbs.realize(state, Stage::Dynamics);
        const Real KE=m_mbs.calcKineticEnergy(state), E=m_mbs.calcEnergy(state);
        const Real diss=m_compliant.getDissipatedEnergy(state);
        DecorativeText txt; txt.setIsScreenText(true);
        txt.setText("KE/Diss/E: " + String(KE, "%.6f") + String(diss, "/%.6f")
                    + String(E+diss, "/%.6f") );
        geometry.push_back(txt);

        int nContPts = 0;
        const int ncont = m_compliant.getNumContactForces(state);
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
            int frcThickness = 2, momThickness = 2;
            Real frcScale = ForceScale, momScale = ForceScale;
            while (frcMag > /*10*/1000000)
                frcThickness++, frcScale /= 10, frcMag /= 10;
            while (momMag > /*10*/1000000)
                momThickness++, momScale /= 10, momMag /= 10;
            geometry.push_back(DecorativePoint(pt)
                               .setScale(5).setColor(Yellow));
            DecorativeLine frcLine(pt,      pt + std::log10(frcMag)*frcDir);
            DecorativeLine momLine(pt+offs, pt+offs + std::log10(momMag)*momDir);
            frcLine.setColor(Black);
            momLine.setColor(Purple);
            frcLine.setLineThickness(frcThickness);
            momLine.setLineThickness(momThickness);
            geometry.push_back(frcLine);
            geometry.push_back(momLine);

            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                ++nContPts;
                const ContactDetail& detail = patch.getContactDetail(i);
                const Vec3& pt = detail.getContactPoint();
                geometry.push_back(DecorativePoint(pt).setColor(Purple));

                const Vec3& force = detail.getForceOnSurface2();
                const Real forceMag = force.norm();
                const UnitVec3 forceDir(force/forceMag, true);
                DecorativeLine frcLine(pt, pt+std::log10(forceMag)*forceDir);
                frcLine.setColor(Black);
                geometry.push_back(frcLine);

                // Make a red line that extends from the contact
                // point in the direction of the slip velocity, of length 3*slipvel.
                DecorativeLine slip(pt,pt+3.*detail.getSlipVelocity());
                slip.setColor(Red);
                geometry.push_back(slip);
            }
        }
        txt.setText(String("Num contact points: ") + String(nContPts));
        geometry.push_back(txt);

    }
private:
    const MultibodySystem&              m_mbs;
    const CompliantContactSubsystem&    m_compliant;
};

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system,
               const CompliantContactSubsystem& complCont,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(complCont)
    {}

    ~MyReporter() {}

    void handleEvent(const State& state) const {
        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};

// These are the item numbers for the entries on the Run menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem=2, QuitItem=3;

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
// This one does nothing but look for the Run->Quit selection.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, Real interval)
    :   PeriodicEventHandler(interval), m_silo(silo) {}

    virtual void handleEvent(State& state, Real accuracy,
                             bool& shouldTerminate) const
    {
        int menuId, item;
        if (m_silo.takeMenuPick(menuId, item) && menuId==RunMenuId && item==QuitItem)
            shouldTerminate = true;
    }

private:
    Visualizer::InputSilo& m_silo;
};

int main() {
  try
  { // Create the system.

    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, UnitVec3(.1,-1,0), 9.81);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTrackDissipatedEnergy(true);
    contactForces.setTransitionVelocity(1e-3);

    const Vec3 hdim(1,2,3);

    const Real fFac =.15; // to turn off friction
    const Real fDis = .5; // to turn off dissipation
    const Real fVis =  .1; // to turn off viscous friction
    const Real fK = .1*1e6; // pascals

    // Halfspace floor
    const Rotation R_xdown(-Pi/2,ZAxis);
    matter.Ground().updBody().addDecoration(
        Transform(Vec3(0,-.5,0)),
        DecorativeBrick(Vec3(10,.5,20)).setColor(Green).setOpacity(.1));
    matter.Ground().updBody().addContactSurface(
        Transform(R_xdown, Vec3(0,0,0)),
        ContactSurface(ContactGeometry::HalfSpace(),
                       ContactMaterial(fK*.1,fDis*.9,
                           fFac*.8,fFac*.7,fVis*10)));

    const Real brickMass = 10;
    Body::Rigid brickBody(MassProperties(brickMass, Vec3(0),
                            UnitInertia::brick(hdim)));
    brickBody.addDecoration(Transform(),
                            DecorativeBrick(hdim).setColor(Cyan).setOpacity(.3));
    const int surfx = brickBody.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Brick(hdim),
                       ContactMaterial(fK,fDis,
                                       fFac*.8,fFac*.7,fVis))
                       );
    //brickBody.addContactSurface(Transform(),
    //    ContactSurface(ContactGeometry::Ellipsoid(hdim),
    //                   ContactMaterial(fK*.1,fDis*.9,
    //                                   .1*fFac*.8,.1*fFac*.7,fVis*1))
    //                   );
    const ContactSurface& surf = brickBody.getContactSurface(surfx);
    const ContactGeometry& cg = surf.getShape();
    const ContactGeometry::Brick& cgbrick = ContactGeometry::Brick::getAs(cg);
    cout << "cgbrick.hdim=" << cgbrick.getHalfLengths() << endl;
    const Geo::Box& box = cgbrick.getGeoBox();
    cout << "box.hdim=" << box.getHalfLengths() << endl;

    // Vertices
    for (int i=0; i<8; ++i) {
        const Vec3 vpos = box.getVertexPos(i);
        const UnitVec3 vn = box.getVertexNormal(i);
        brickBody.addDecoration
            (DecorativePoint(vpos).setColor(Orange));
        brickBody.addDecoration
            (DecorativeText(String(i)).setTransform(vpos).setColor(White)
             .setScale(.5));
        brickBody.addDecoration
           (DecorativeLine(vpos, vpos + 0.5*vn).setColor(Orange));

        printf("vertex %d:\n", i);
        int e[3],ew[3],f[3],fw[3];
        box.getVertexEdges(i,e,ew);
        box.getVertexFaces(i,f,fw);
        for (int ex=0; ex<3; ++ex) {
            int ev[2]; box.getEdgeVertices(e[ex], ev);
            printf("  e%2d(%d) ev=%d\n", e[ex], ew[ex], ev[ew[ex]]);
        }
        for (int fx=0; fx<3; ++fx) {
            int fv[4]; box.getFaceVertices(f[fx], fv);
            printf("  f%2d(%d) fv=%d\n", f[fx], fw[fx], fv[fw[fx]]);
        }
    }

    // Edges
    for (int i=0; i<12; ++i) {
        const UnitVec3 n = box.getEdgeNormal(i);
        const UnitVec3 d = box.getEdgeDirection(i);
        const Vec3 ctr = box.getEdgeCenter(i);
        const Real len = .75;
        brickBody.addDecoration
           (DecorativePoint(ctr).setColor(Green).setScale(2));
        brickBody.addDecoration
            (DecorativeText(String(i)).setTransform(ctr+len*n)
             .setColor(Green).setScale(.3));
        brickBody.addDecoration
           (DecorativeLine(ctr, ctr + len*n).setColor(Green));
        brickBody.addDecoration
           (DecorativeLine(ctr, ctr + len*d).setColor(Green));

        printf("edge %d:\n", i);
        int f[2],fw[2];
        box.getEdgeFaces(i,f,fw);
        for (int fx=0; fx<2; ++fx) {
            int fe[4]; box.getFaceEdges(f[fx], fe);
            printf("  f%2d(%d) fe=%d\n", f[fx], fw[fx], fe[fw[fx]]);
        }
    }

    // Faces
    for (int i=0; i<6; ++i) {
        int vertices[4]; box.getFaceVertices(i,vertices);
        const UnitVec3 n = box.getFaceNormal(i);
        const Vec3 ctr = box.getFaceCenter(i);
        brickBody.addDecoration
           (DecorativePoint(ctr).setColor(Magenta).setScale(3));
        brickBody.addDecoration
           (Transform(Rotation(n,ZAxis,Vec3(0,1,0),YAxis),ctr),
            DecorativeText(String(i)).setColor(Magenta)
             .setScale(.75).setFaceCamera(false));
        brickBody.addDecoration
           (DecorativeLine(ctr, ctr + 1.*n).setColor(Magenta));
    }

    MobilizedBody::Free brick(matter.Ground(), Transform(Vec3(0,3,0)),
        brickBody, Transform(Vec3(0)));

    Visualizer viz(system);
    viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
    viz.setShowShadows(true);
    viz.setShowSimTime(true);
    viz.setDesiredFrameRate(FrameRate);
    viz.setShowFrameRate(true);
    viz.setBackgroundType(Visualizer::SolidColor);
    viz.setBackgroundColor(White*.9);

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

    system.addEventReporter(new MyReporter(system,contactForces,ReportInterval));
    system.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

    // Check for a Run->Quit menu pick every 1/4 second.
    system.addEventHandler(new UserInputHandler(*silo, .25));

    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();

    brick.setQToFitRotation(state, Rotation(SpaceRotationSequence,
                                            .1, ZAxis, .05, XAxis));
    brick.setUToFitLinearVelocity(state, Vec3(2,0,0));

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
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    //integ.setAccuracy(1e-3); // minimum for CPodes
    integ.setAccuracy(1e-5);
    //integ.setAccuracy(.01);
    TimeStepper ts(system, integ);


    ts.initialize(state);
    double cpuStart = cpuTime();
    double realStart = realTime();

    ts.stepTo(20.0);

    const double timeInSec = realTime() - realStart;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s elapsed for " << ts.getTime() << "s sim (avg step="
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) "
        << (1000*ts.getTime())/evals << "ms/eval\n";
    cout << "  CPU time was " << cpuTime() - cpuStart << "s\n";

    printf("Using Integrator %s at accuracy %g:\n",
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

    viz.dumpStats(std::cout);

    // Add as slider to control playback speed.
    viz.addSlider("Speed", 1, 0, 4, 1);
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

