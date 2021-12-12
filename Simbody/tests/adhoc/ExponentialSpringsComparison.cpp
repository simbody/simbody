/*-----------------------------------------------------------------------------
                               Simbody(tm)
-------------------------------------------------------------------------------
 Copyright (c) 2021 Authors.
 Authors: Frank C. Anderson
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ------------------------------------------------------------------------------*/

// ExponentialSpringsComparison.cpp
// This file contains the main() function. Program execution begins and ends there.
//

#include "Simbody.h"
#include <iostream>

using namespace SimTK;
using std::endl;
using std::cout;
using std::unique_ptr;
using std::make_unique;

Array_<State> saveEm;

// These are the item numbers for the entries on the Run menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem = 2, QuitItem = 3;


//=============================================================================
// Record the State at regular intervals.
class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, Real reportInterval)
        : PeriodicEventReporter(reportInterval) {
    }

    ~MyReporter() {}

    void handleEvent(const State& state) const override {
        saveEm.push_back(state);
    }
private:
};


//=============================================================================
/** This class impmlements an event reported that identifies the highest
position reached by a body.  If energy is conserved, which should
be the case if there are no damping components of the contact force,
then the height should be constant across multiple bounces. */
class MaxHeightReporter : public TriggeredEventReporter {
public:
    MaxHeightReporter(const MultibodySystem& system, const MobilizedBody& body)
        : TriggeredEventReporter(Stage::Velocity), system(system), body(body) {
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }
    Real getValue(const State& state) const {
        Vec3 vel = body.getBodyOriginVelocity(state);
        return vel[1];
    }
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Velocity);
        Vec3 pos = body.getBodyOriginLocation(state);
        Vec3 vel = body.getBodyOriginVelocity(state);
        //cout << endl;
        //cout << state.getTime() << "\tp = " << pos << "\tv = " << vel << endl;
        //cout << endl;
    }
private:
    const MultibodySystem& system;
    const MobilizedBody& body;
};

//=============================================================================
/** This class impmlements an event reported that access the force exerted
by an ExponentialSpringForce at regular intervals. */
class ExpSprForceReporter : public PeriodicEventReporter {
public:
    ExpSprForceReporter(const MultibodySystem& system,
        const ExponentialSpringForce& spr, Real reportInterval)
        : PeriodicEventReporter(reportInterval), system(system), spr(spr) {
    }
    void handleEvent(const State& state) const override {
        system.realize(state, Stage::Dynamics);
        Vec3 f_G = spr.getForce(state);
        Real sliding = spr.getSliding(state);
        Vec3 fric = f_G; fric[1] = 0.0;
        //cout << state.getTime() << "\tSliding= "<< sliding <<
        //    "  \tmu= "<< fric.norm() / f_G[1] <<
        //    "  \tf_G= " << f_G << endl;

    }
private:
    const MultibodySystem& system;
    const ExponentialSpringForce& spr;
};


//=============================================================================
int main() {
    std::cout << "Running a simulation in which contact with a floor is " <<
        "modeled with exponential springs." << endl;
    try {
        // I should probably make the following options command line
        // areguments, but since this code is part of the build environment,
        // it is easy to recomplile.

        // Use compliant Contact?
        bool CmpContactOn = true;

        // Use exponential Spring Contact?
        bool ExpContactOn = true;

        // Run with the vizualizer?
        bool VisOn = true;

        // Specify initial conditions.
        int condition = 4;
        // 0 = sitting still
        // 1 = dropped
        // 2 = sliding
        // 3 = spinning
        // 4 = spinning & sliding
        // 5 = spinning like a top
        // 6 = tumbling
        // 7 = non-vertical gravity (simple way to create a sideways force)
        
        // Create the system.
        Real hs = 0.1;
        MultibodySystem system; system.setUseUniformBackground(true);
        SimbodyMatterSubsystem matter(system);
        GeneralForceSubsystem forces(system);
        Vec3 g(0., -9.8, 0.);
        if(condition == 7) {
            // An angle of 90 deg is vertical. 80 is 10 deg short of vert.
            Real angle = convertDegreesToRadians(80.0);
            g = Vec3(9.8 * cos(angle), -9.8 * sin(angle), 0.0);
        }
        Force::UniformGravity gravity(forces, matter, g);
        Body::Rigid BodyPropsExp(MassProperties(10.0, Vec3(0), Inertia(1)));
        BodyPropsExp.addDecoration(Transform(),
            DecorativeBrick(Vec3(hs)).setColor(Blue));
        Body::Rigid BodyPropsCmp(MassProperties(10.0, Vec3(0), Inertia(1)));
        BodyPropsCmp.addDecoration(Transform(),
            DecorativeBrick(Vec3(hs)).setColor(Red));

     
        // Add a 6 dof mass that uses the Compliant Contact system.
        // Some constants
        const Vec3 hdim(hs, hs, hs);
        const Real mu_s = 0.7; // static coefficient of friction
        const Real mu_k = 0.5; // kinetic (or dynamic) coefficient of friction
        const Real fDis = .55; // to turn off dissipation
        const Real fVis = 0.0; // to turn off viscous friction
        const Real fK = .5 * 1e6; // pascals
        const Rotation R_xdown(-Pi / 2, ZAxis);
        //int surfx;
        matter.Ground().updBody().addDecoration(
            Transform(Vec3(0, -0.5, 0)),
            DecorativeBrick(Vec3(2, 0.5, 2)).setColor(Green).setOpacity(.1));
        // Compliant contact infrastructure -------
        unique_ptr<ContactTrackerSubsystem> tracker;
        unique_ptr<CompliantContactSubsystem> contactForces;
        unique_ptr<MobilizedBody::Free> blockCmp;
        if (CmpContactOn == true) { 
            tracker = make_unique<ContactTrackerSubsystem>(system);
            contactForces =
                make_unique<CompliantContactSubsystem>(system, *tracker);
            contactForces->setTrackDissipatedEnergy(true);
            contactForces->setTransitionVelocity(1.0e-3);
            // Ground
            matter.Ground().updBody().addContactSurface(
                Transform(R_xdown, Vec3(0, 0, 0)),
                ContactSurface(ContactGeometry::HalfSpace(),
                    ContactMaterial(fK*.1, fDis*.9, mu_s, mu_k, fVis*10)));
            // Free Body (cube)
            int surfx =
                BodyPropsCmp.addContactSurface(Transform(),
                    ContactSurface(ContactGeometry::Brick(hdim),
                        ContactMaterial(fK, fDis, mu_s, mu_k, fVis)));
            blockCmp = make_unique<MobilizedBody::Free>(
                matter.Ground(), BodyPropsCmp);
        }

        // Add a 6 dof mass that uses exponential springs for contact (new)
        unique_ptr<MobilizedBody::Free> blockExp;
        unique_ptr<ExponentialSpringForce> spr[8];
        int i;
        // Define the corners of the block
        Vec3 corner[8];
        corner[0] = Vec3(hs, -hs, hs);
        corner[1] = Vec3(hs, -hs, -hs);
        corner[2] = Vec3(-hs, -hs, -hs);
        corner[3] = Vec3(-hs, -hs, hs);
        corner[4] = Vec3(hs, hs, hs);
        corner[5] = Vec3(hs, hs, -hs);
        corner[6] = Vec3(-hs, hs, -hs);
        corner[7] = Vec3(-hs, hs, hs);
        if(ExpContactOn) {
            blockExp = make_unique<MobilizedBody::Free>(
                    matter.Ground(), BodyPropsExp);

            // Create a Transform representing the contact plane.
            // The ZAxis of contact plane should point up.
            Real angle = convertDegreesToRadians(90.0);
            Rotation floorRot(-angle, XAxis);
            Vec3 floorOrigin(0., -0.004, 0.);
            Transform floorXForm(floorRot, floorOrigin);
            // Modify the default parameters if desired
            ExponentialSpringParameters params;
            params.setElasticityAndViscosityForCriticalDamping(20000.0);
            // Add an exponential spring at each corner of the block.
            for(i = 0; i < 8; ++i) {
                spr[i] = make_unique<ExponentialSpringForce>(system,
                    floorXForm, *blockExp, corner[i], mu_s, mu_k, params);
            }
        }

        // Add reporting and visualization
        unique_ptr<Visualizer> viz;
        Visualizer::InputSilo* silo;
        if (VisOn) {
            viz = make_unique<Visualizer>(system);
            viz->setShowShadows(true);

            // ----- Run Menu ----
            silo = new Visualizer::InputSilo();
            viz->addInputListener(silo);
            Array_<std::pair<String, int> > runMenuItems;
            runMenuItems.push_back(std::make_pair("Go", GoItem));
            runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
            viz->addMenu("Run", RunMenuId, runMenuItems);

            system.addEventReporter(new MyReporter(system, 1.0/30.0));
            system.addEventReporter(new Visualizer::Reporter(*viz, 1./30.));
            if(ExpContactOn) {
                //system.addEventReporter(
                //new MaxHeightReporter(system, *blockExp));
                system.addEventReporter(
                    new ExpSprForceReporter(system, *spr[0], 0.01));
            }
        }

        // Initialize the system and state.
        system.realizeTopology();
        State state = system.getDefaultState();
        system.realizeModel(state);
        Rotation R(convertDegreesToRadians(0.0), Vec3(1, 0, 0));
        // Compliant Contact
        // I should probably haved used a switch block!
        if (CmpContactOn) {
            // Sitting Still
            if ((condition == 0)||(condition>6)) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 0.1, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Dropped
            if (condition == 1) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 4.0, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Sliding
            if (condition == 2) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 0.1, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 4.0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Spinning flat
            if (condition == 3) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 0.1, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 0.0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 8.0, 0.0));
            }
            // Spinning and Sliding
            if (condition == 4) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 0.1, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 2.0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0,12.0,0.0));
            }
            // Spinning top
            if (condition == 5) {
                R.setRotationFromAngleAboutNonUnitVector(
                    convertDegreesToRadians(54.74),
                    Vec3(1, 0, 1));
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-0.5, 0.2, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 6.0, 0.0));
            }
            // Tumbling
            if (condition == 6) {
                blockCmp->setQToFitRotation(state, R);
                blockCmp->setQToFitTranslation(state, Vec3(-1.5, 1.0, 0.0));
                blockCmp->setU(state, Vec6(0, 0, 0, 0.8, 0, 0));
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0,0.0,-10.0));
            }
        }
        // Exponential Contact
        // Same initial conditions as for Compliant Contact except shifted so
        // both blocks can be seen.
        if (ExpContactOn) {
            // Test modifying cooefficients of friction (states)
            // Static
            Real mus = mu_s;
            Real muk = mu_k;
            for(i = 0; i < 8; ++i) {
                spr[i]->setMuStatic(state, mus);
                spr[i]->setMuKinetic(state, muk);
                if((condition == 0) || (condition == 7))
                    spr[i]->setSliding(state, 0.0);
            }
            // Sitting Still or Non-Vertical Force
            if ((condition == 0) || (condition > 6)) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(-0.5, 0.0993, 0.5));
                blockExp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Dropped
            if (condition == 1) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(0.5, 4.0, 0.0));
                blockExp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Sliding
            if (condition == 2) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(-0.5, 0.1, 0.5));
                blockExp->setU(state, Vec6(0, 0, 0, 4.0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, 0.0));
            }
            // Spinning
            if (condition == 3) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(0.5, 0.1, 0.0));
                blockExp->setU(state, Vec6(0, 0, 0, 0.0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 8.0, 0.0));
            }
            // Spinning and Sliding
            if (condition == 4) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(0.5, 0.1, 0.0));
                blockExp->setU(state, Vec6(0, 0, 0, 2.0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0,12.0,0.0));
            }
            // Spinning top
            if (condition == 5) {
                R.setRotationFromAngleAboutNonUnitVector(
                    convertDegreesToRadians(54.74),Vec3(1,0,1));
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(0.5, 0.2, 0.0));
                blockExp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0,6.0,0.0));
            }
            // Tumbling
            if (condition == 6) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(-1.5, 1.0, 0.5));
                blockExp->setU(state, Vec6(0, 0, 0, 0.8, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0,0.0,-10.0));
            }
            // Reset the spring zeros
            // The spring zeros are set to the point on the Floor that
            // coincices with the Station that was specified on the Body.
            for(i = 0; i < 8; ++i) spr[i]->resetSpringZero(state);
        }

        // Run Menu and storage for Replay
        if (VisOn) {
            saveEm.reserve(50000);
            cout << "\nChoose 'Go' from Run menu to simulate:\n";
            int menuId, item;
            do {
                silo->waitForMenuPick(menuId, item);
                if (menuId != RunMenuId || item != GoItem)
                    cout << "\aDude ... follow instructions!\n";
            } while (menuId != RunMenuId || item != GoItem);
        }

        // Simulate it.
        RungeKuttaMersonIntegrator integ(system);
        integ.setAccuracy(1.0e-5);
        integ.setMaximumStepSize(0.05);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        double cpuStart = cpuTime();
        double realStart = realTime();
        ts.stepTo(10.0);
        double cpuDuration = cpuTime() - cpuStart;
        double realDuration = realTime() - realStart;
        // Report integrator performance
        cout << std::endl << "Integrator = " << integ.getMethodName() << endl;
        cout << "    Accuracy in use = " << integ.getAccuracyInUse() << endl;
        cout << "Num steps attempted = "<<integ.getNumStepsAttempted()<<endl;
        cout << "    Num steps taken = " << integ.getNumStepsTaken() << endl;
        cout << " cpu time taken = " << cpuDuration << " sec" << endl;
        cout << "real time taken = " << realDuration << " sec" << endl;

        // Replay Loop
        if (VisOn) {
            silo->clear(); // forget earlier input
            while (true) {
                cout << "Choose Replay to see that again ...\n";

                int menuId, item;
                silo->waitForMenuPick(menuId, item);

                for (double i = 0; i < (int)saveEm.size(); i++) {
                    viz->report(saveEm[(int)i]);
                }
            }
        }

    }
    catch (const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cout << "UNKNOWN EXCEPTION\n";
        return 1;
    }
    return 0;
}
