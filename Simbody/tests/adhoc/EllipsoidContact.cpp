/* -------------------------------------------------------------------------- *
 *                            Simbody(tm) Example                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/* This example is for experimenting with ellipsoid contact which was 
introduced in Simbody 2.2.
*/

#include "SimTKsimbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real TimeScale = 1;
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
        cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
             << " Ediss=" << m_compliant.getDissipatedEnergy(state)
             << " E+Ediss=" << m_system.calcEnergy(state)
                               +m_compliant.getDissipatedEnergy(state)
             << endl;
        const int ncont = m_compliant.getNumContactForces(state);
        cout << "Num contacts: " << m_compliant.getNumContactForces(state) << endl;
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            //cout << force;
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

    virtual void handleEvent(State& state, Real accuracy, const Vector& yWeights, 
                             const Vector& ooConstraintTols, Stage& lowestModified, 
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
    Force::Gravity   gravity(forces, matter, UnitVec3(0,-1,0), 9.8);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTransitionVelocity(1e-3);

    ContactCliqueId clique1 = ContactSurface::createNewContactClique();
    ContactCliqueId clique2 = ContactSurface::createNewContactClique();
    ContactCliqueId clique3 = ContactSurface::createNewContactClique();

    const Real fFac =1; // to turn off friction
    const Real fDis = .5*0.2; // to turn off dissipation
    const Real fVis =  .1*.1; // to turn off viscous friction
    const Real fK = 100*1e6; // pascals

    // Right hand wall
    //matter.Ground().updBody().addDecoration(Vec3(.25+.01,0,0),
    //    DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.2));
    //matter.Ground().updBody().addContactSurface(Vec3(.25,0,0),
    //    ContactSurface(ContactGeometry::HalfSpace(),
    //                   ContactMaterial(fK*.01,fDis*.9,fFac*.8,fFac*.7,fVis*10))
    //                   .joinClique(clique3));

    // Halfspace floor
    const Rotation R_xdown(-Pi/2,ZAxis);
    //matter.Ground().updBody().addDecoration(
    //    Transform(R_xdown, Vec3(0,-3-.01,0)),
    //    DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.1));
    matter.Ground().updBody().addContactSurface(
        Transform(R_xdown, Vec3(0,0,0)),
        ContactSurface(ContactGeometry::HalfSpace(),
                       ContactMaterial(fK*.1,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       //ContactMaterial(2e6,.01,.1,.05,.01))
                       .joinClique(clique1));



    const Real ellipsoidMass = 100;
    const Vec3 halfDims(10, 2, 1);
    Body::Rigid ellipsoidBody(MassProperties(ellipsoidMass, Vec3(0), 
                              ellipsoidMass*UnitInertia::ellipsoid(halfDims)));

    ellipsoidBody.addDecoration(Transform(), 
        DecorativeEllipsoid(halfDims).setColor(Cyan)
         .setOpacity(.5).setResolution(10));
    ellipsoidBody.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Ellipsoid(halfDims),
                       ContactMaterial(fK*.1,fDis*.9,
                                       .1*fFac*.8,.1*fFac*.7,fVis*1))
                       //ContactMaterial(2e6,.01,.1,.05,.01))
                       //.joinClique(clique2)
                       );
    MobilizedBody::Free ellipsoid(matter.Ground(), Transform(Vec3(0,0,0)),
        ellipsoidBody, Transform(Vec3(0)));


    Visualizer viz(system);
    viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
    viz.setMode(Visualizer::RealTime);
    viz.setDesiredBufferLengthInSec(1);
    viz.setDesiredFrameRate(FrameRate);
    viz.setShowShadows(true);

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
    ellipsoid.setQToFitTransform(state, Transform(Rotation(),
                                             Vec3(0,3,0)));

    viz.report(state);
    printf("Default state\n");

    cout << "\nChoose 'Go' from Run menu to simulate:\n";
    int menuId, item;
    do { silo->waitForMenuPick(menuId, item);
         if (menuId != RunMenuId || item != GoItem) 
             cout << "\aDude ... follow instructions!\n";
    } while (menuId != RunMenuId || item != GoItem);

    ellipsoid.setOneU(state, 1, .1);
    ellipsoid.setOneU(state, 0, .05); 
    
    // Simulate it.

    //ExplicitEulerIntegrator integ(system);
    //CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    integ.setAccuracy(1e-3); // minimum for CPodes
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
