
/**
This is a stripped down version of the collision playground example to illustrate the
behavior of the collision detection system and a few geometry dependent odities.
I have removed the ground plane and collision cliques for the purposes of this
example.

try the following:
1)all spheres: collisions are detected forces are drawn but not applied. No collision interaction;
2)all cubes(size 0.4-0.5): Collisions are detected but notice that the moving cube "cuts" into the stationary one
    before collision is detected.  No forces are generated at this stage.  Collisions seem to occur
    on a spherical boundary within the cube. This is similar to the contact playground example.

*/
#include "SimTKsimbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <string.h>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real TimeScale = 0.21;
static const Real FrameRate = 30;
static const Real ReportInterval = TimeScale/FrameRate;
static const Real ForceScale = .25;
static const Real MomentScale = .5;


class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont)
        :   m_system(system), m_compliant(complCont) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) {
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_system.realize(state, Stage::Velocity);

        const int ncont = m_compliant.getNumContactForces(state);
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            const Vec3& frc = force.getForceOnSurface2()[1];

            const Vec3& mom = force.getForceOnSurface2()[0];
            Real  frcMag = frc.norm(), momMag=mom.norm();
            int frcThickness = 1, momThickness = 1;
            Real frcScale = ForceScale, momScale = ForceScale;
            while (frcMag > 10)
                frcThickness++, frcScale /= 10, frcMag /= 10;
            while (momMag > 10)
                momThickness++, momScale /= 10, momMag /= 10;
            DecorativeLine frcLine(force.getContactPoint(),
                                   force.getContactPoint() + frcScale*frc);
            DecorativeLine momLine(force.getContactPoint(),
                                   force.getContactPoint() + momScale*mom);
            frcLine.setColor(frcColors[id%3]);
            momLine.setColor(momColors[id%3]);
            frcLine.setLineThickness(2*frcThickness);
            momLine.setLineThickness(2*momThickness);
            geometry.push_back(frcLine);
            geometry.push_back(momLine);

            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                const ContactDetail& detail = patch.getContactDetail(i);
                const Real peakPressure = detail.getPeakPressure();
                // Make a black line from the element's contact point in the normal
                // direction, with length proportional to log(peak pressure)
                // on that element.

                DecorativeLine normal(detail.getContactPoint(),
                                      detail.getContactPoint()+ std::log10(peakPressure)
                                      * detail.getContactNormal());
                normal.setColor(Black);
                geometry.push_back(normal);

                // Make a red line that extends from the contact
                // point in the direction of the slip velocity, of length 3*slipvel.
                DecorativeLine slip(detail.getContactPoint(),
                                    detail.getContactPoint()+3*detail.getSlipVelocity());
                slip.setColor(Red);
                geometry.push_back(slip);
            }
        }
    }
private:
    const MultibodySystem&              m_system;
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
        m_system.realize(state, Stage::Dynamics);
        /*   cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
             << " Ediss=" << m_compliant.getDissipatedEnergy(state)
             << " E+Ediss=" << m_system.calcEnergy(state)
                +m_compliant.getDissipatedEnergy(state)
             << endl;*/
        const int ncont = m_compliant.getNumContactForces(state);
        if(ncont>0)
        {
            cout << "Num contacts: " << m_compliant.getNumContactForces(state) << endl;
            for (int i=0; i < ncont; ++i) {
                const ContactForce& force = m_compliant.getContactForce(state,i);
                cout << force;
            }
        }
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

        /// uncoment gravity to get some sort of collision interaction
        /// for cylinder mesh
        // Force::UniformGravity gravity(forces, matter,Vec3(0,0.001,0), 2);

        ContactTrackerSubsystem  tracker(system);
        //GeneralContactSubsystem contactsys(system);
        CompliantContactSubsystem contactForces(system, tracker);
        contactForces.setTrackDissipatedEnergy(true);

        for(SubsystemIndex i(0); i<system.getNumSubsystems(); ++i)
        {
            fprintf(stderr,"subsytem name %d %s\n", (int)i,
                system.getSubsystem((SubsystemIndex)i).getName().c_str());
        }

        const Real rad = .4;
        PolygonalMesh pyramidMesh1,pyramidMesh2;

        /// load cylinder forces drawn, but interaction depends on gravity???


        const Real fFac =1; // to turn off friction
        const Real fDis = .5*0.2; // to turn off dissipation
        const Real fVis =  .1*.1; // to turn off viscous friction
        const Real fK = 100*1e6; // pascals

        Body::Rigid pendulumBody3(MassProperties(100.0, Vec3(0), 100*Inertia(1)));
        PolygonalMesh body3contact = PolygonalMesh::createSphereMesh(rad, 2);
        ContactGeometry::TriangleMesh geo3(body3contact);

        const DecorativeMesh mesh3(geo3.createPolygonalMesh());
        pendulumBody3.addDecoration(Transform(),
                                    DecorativeMesh(mesh3).setOpacity(.2));
        pendulumBody3.addDecoration(Transform(),
                                    DecorativeMesh(mesh3).setColor(Gray)
                                    .setRepresentation(DecorativeGeometry::DrawWireframe)
                                    .setOpacity(.1));
        ContactSurface s1(geo3,
                          ContactMaterial(fK*.1,fDis*.9,fFac*.8,fFac*.7,fVis*10));
        s1.setThickness(1);
        s1.setShape(geo3);
        //ContactGeometry::Sphere geo3(rad);
        pendulumBody3.addContactSurface(Transform(),s1);
        /*
                std::ifstream meshFile1,meshFile2;
                meshFile1.open("cyl3.obj");
                pyramidMesh1.loadObjFile(meshFile1); meshFile1.close();
*/
        pyramidMesh1 = PolygonalMesh::createSphereMesh(rad, 2);
        ContactGeometry::TriangleMesh pyramid1(pyramidMesh1);

        DecorativeMesh showPyramid1(pyramid1.createPolygonalMesh());
        const Real ballMass = 200;
        Body::Rigid ballBody(MassProperties(ballMass, Vec3(0),
                                            ballMass*UnitInertia::sphere(1)));

        ballBody.addDecoration(Transform(),
                               showPyramid1.setColor(Cyan).setOpacity(.2));
        ballBody.addDecoration(Transform(),
                               showPyramid1.setColor(Gray)
                               .setRepresentation(DecorativeGeometry::DrawWireframe));

        ContactSurface s2(pyramid1,
                          ContactMaterial(fK*.1,fDis*.9,
                                          .1*fFac*.8,.1*fFac*.7,fVis*1));
        s2.setThickness(1);
        s2.setShape(pyramid1);
        ballBody.addContactSurface(Transform(),/*ContactSurface(ContactGeometry::Sphere(rad),ContactMaterial(fK*.1,fDis*.9,
                                                              .1*fFac*.8,.1*fFac*.7,fVis*1))*/  s2/*.joinClique(clique1)*/);

        /*   Body::Rigid d(MassProperties(1.0, Vec3(0),Inertia(1)));

        MobilizedBody::Pin dud(matter.Ground(),Transform(),d,Transform());
*/
        MobilizedBody::Free ball(matter.Ground(), Transform(Vec3(-2,-2,0)),
                                 ballBody, Transform(Vec3(0)));



        MobilizedBody::Free ball1(matter.Ground(), Transform(Vec3(0,0,0)),
                                  ballBody, Transform(Vec3(0)));
        /*
        MobilizedBody::Free ball2(matter.Ground(), Transform(Vec3(-4,0,0)),
                                 ballBody, Transform(Vec3(0)));
*/

        MobilizedBody::Free ball3(matter.Ground(), Transform(Vec3(-1,-2,0)),
                                  ballBody, Transform(Vec3(0)));


        MobilizedBody::Pin pendulum3(matter.Ground(), Transform(Vec3(-2,0,0)),
                                     pendulumBody3, Transform(Vec3(0, 2, 0)));

        ball.updBody();
        ball1.updBody();

        Visualizer viz(system);
        viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
        viz.setMode(Visualizer::RealTime);
        viz.setDesiredBufferLengthInSec(1);
        viz.setDesiredFrameRate(FrameRate);
        viz.setGroundHeight(-3);
        viz.setShowShadows(true);
        viz.setBackgroundType(Visualizer::SolidColor);
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
        //  system.addEventHandler(new TriggeredEventHandler(Stage::Model));
        // Initialize the system and state.

        system.realizeTopology();

        State state = system.getDefaultState();
        /*
        ball.setQToFitTransform(state, Transform(Rotation(Pi/2,XAxis),
                                                 Vec3(0,-1.8,0)));
*/
        //pendulum.setOneQ(state, 0, -Pi/12);
        pendulum3.setOneQ(state, 0, -Pi/2);
        pendulum3.setOneU(state, 0, Pi/4);
        // ball.setOneU(state, 1, 0.1);
        viz.report(state);
        matter.updAllParticleVelocities(state);
        printf("Default state\n");
        /* cout << "t=" << state.getTime()
         << " q=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state)
         << " u=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state)
         << endl;
*/
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

        //ExplicitEulerIntegrator integ(system);
        CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
        //RungeKuttaFeldbergIntegrator integ(system);
        //RungeKuttaMersonIntegrator integ(system);
        //RungeKutta3Integrator integ(system);
        //VerletIntegrator integ(system);
        //integ.setMaximumStepSize(1e-1);
        //integ.setAllowInterpolation(false);
        integ.setAccuracy(1e-3); // minimum for CPodes
        //integ.setAccuracy(.1);
        TimeStepper ts(system, integ);


        ts.initialize(state);
        double cpuStart = cpuTime();
        double realStart = realTime();

        ts.stepTo(2000.0);

        const double timeInSec = realTime() - realStart;
        const int evals = integ.getNumRealizations();
        /*  cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s elapsed for " << ts.getTime() << "s sim (avg step="
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) "
        << (1000*ts.getTime())/evals << "ms/eval\n";
    cout << "  CPU time was " << cpuTime() - cpuStart << "s\n";

    printf("Using Integrator %s at accuracy %g:\n",
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());
*/
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


