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

// TestExponentialSprings.cpp
// Authors: Frank C. Anderson
// This file contains the main() function. Program execution begins and ends there.
//

#include "Simbody.h"
#include <iostream>

using namespace SimTK;
using std::endl;
using std::cout;

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
        cout << endl;
        cout << state.getTime() << "\tp = " << pos << "\tv = " << vel << endl;
        cout << endl;
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
        const ExponentialSpringData& data = spr.getData(state);
        cout << state.getTime() << "\tf_G = " << data.f_G << endl;
    }
private:
    const MultibodySystem& system;
    const ExponentialSpringForce& spr;
};


//=============================================================================
int main() {
    std::cout << "Running a simulation in which contact with the ground is modeled with exponential springs.\n";
    try {
        // I should probably make the following options command line areguments,
        // but since this code is part of the build environment, it is easy to
        // recomplile.  I'll refine once the exponential spring classes pass
        // basic muster.

        // Use compliant Contact?
        bool CmpContactOn = true;

        // Use exponential Spring Contact?
        bool ExpContactOn = true;

        // Run with the vizualizer?
        bool VisOn = true;

        // Specify initial conditions.
        int condition = 2;
        // 0 = sitting still
        // 1 = dropped
        // 2 = sliding
        // 3 = spinning
        // 4 = spinning & sliding
        // 5 = spinning like a top
        // 6 = tumbling
        // 7 = titled plane 5.71 deg
        // 8 = titled plane 24.2 deg
        // 9 = titled plane 26.6 deg
        // 10 = titled plane 36.87 deg
        
        // Create the system.
        MultibodySystem system; system.setUseUniformBackground(true);
        SimbodyMatterSubsystem matter(system);
        GeneralForceSubsystem forces(system);
        // tilt = 0 deg
        if(condition < 7) Force::UniformGravity
            gravity(forces, matter, Vec3(0., -9.8, 0.));
        // tilt=5.71 deg, mu=0.10 slips
        if(condition == 7) Force::UniformGravity
            gravity(forces, matter, Vec3(0.9750, -9.7514, 0.));
        // tilt=24.2 deg, mu=0.45 slips
        if(condition == 8) Force::UniformGravity
            gravity(forces, matter, Vec3(4.0172, -8.9388, 0.));
        // tilt=26.6 deg, mu=0.5008 slips
        if(condition == 9) Force::UniformGravity
            gravity(forces, matter, Vec3(4.3880, -8.7627, 0.));
        // tilt=36.87 deg, mu=0.75 slips
        if(condition == 10) Force::UniformGravity
            gravity(forces, matter, Vec3(5.0385, -7.8400, 0.));
        // tilt = 31.0 deg, mu = 0.60
        //if (condition == 10) Force::UniformGravity
        //gravity(forces, matter, Vec3(5.0385, -8.4055, 0.));
        Body::Rigid BodyPropsExp(MassProperties(10.0, Vec3(0), Inertia(1)));
        BodyPropsExp.addDecoration(Transform(), DecorativeBrick(Vec3(0.1)).setColor(Blue));
        Body::Rigid BodyPropsCmp(MassProperties(10.0, Vec3(0), Inertia(1)));
        BodyPropsCmp.addDecoration(Transform(), DecorativeBrick(Vec3(0.1)).setColor(Red));

     
        // Add a 6 dof mass that uses the Compliant Contact system.
        // Some constants
        const Vec3 hdim(0.1, 0.1, 0.1);
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
        ContactTrackerSubsystem* tracker = NULL;
        CompliantContactSubsystem* contactForces = NULL;
        MobilizedBody::Free* blockCmp = NULL;
        if (CmpContactOn == true) { 
            tracker = new ContactTrackerSubsystem(system);
            contactForces = new CompliantContactSubsystem(system, *tracker);
            contactForces->setTrackDissipatedEnergy(true);
            contactForces->setTransitionVelocity(1.0e-3);
            // Ground
            matter.Ground().updBody().addContactSurface(Transform(R_xdown, Vec3(0, 0, 0)),
                ContactSurface(ContactGeometry::HalfSpace(),
                    ContactMaterial(fK * .1, fDis * .9, mu_s, mu_k, fVis * 10)));
            // Free Body (cube)
            int surfx =
                BodyPropsCmp.addContactSurface(Transform(),
                    ContactSurface(ContactGeometry::Brick(hdim),
                        ContactMaterial(fK, fDis, mu_s, mu_k, fVis)));
            blockCmp = new MobilizedBody::Free(matter.Ground(), BodyPropsCmp);
        }

        // Add a 6 dof mass that uses exponential springs for contact (new)
        MobilizedBody::Free* blockExp = NULL;
        ExponentialSpringForce* spr1 = NULL;
        ExponentialSpringForce* spr2 = NULL;
        ExponentialSpringForce* spr3 = NULL;
        ExponentialSpringForce* spr4 = NULL;
        ExponentialSpringForce* spr5 = NULL;
        ExponentialSpringForce* spr6 = NULL;
        ExponentialSpringForce* spr7 = NULL;
        ExponentialSpringForce* spr8 = NULL;
        if(ExpContactOn) {
            blockExp = new MobilizedBody::Free(matter.Ground(), BodyPropsExp);

            // Create a Transform representing the floor
            Real angle = convertDegreesToRadians(0.0);
            Rotation floorTilt(angle, ZAxis);
            Vec3 floorOrigin(0., -0.004, 0.);
            Transform floorXForm(floorTilt, floorOrigin);
            // Modify the default parameters if desired
            ExponentialSpringParameters params;
            //params.setElasticityAndComputeViscosity(100000.0);
            // Add an exponential spring at each corner of the block.
            spr1 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(0.1, -0.1, 0.1), mu_s, mu_k, params);
            spr2 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(0.1, -0.1, -0.1), mu_s, mu_k, params);
            spr3 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(-0.1, -0.1, 0.1), mu_s, mu_k, params);
            spr4 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(-0.1, -0.1, -0.1), mu_s, mu_k, params);
            spr5 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(0.1, 0.1, 0.1), mu_s, mu_k, params);
            spr6 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(0.1, 0.1, -0.1), mu_s, mu_k, params);
            spr7 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(-0.1, 0.1, 0.1), mu_s, mu_k, params);
            spr8 = new ExponentialSpringForce(system, floorXForm,
                *blockExp, Vec3(-0.1, 0.1, -0.1), mu_s, mu_k, params);
        }

        // Add reporting and visualization
        Visualizer* viz = NULL;
        Visualizer::InputSilo* silo = NULL;
        if (VisOn) {
            viz = new Visualizer(system);
            viz->setShowShadows(true);

            // ----- Run Menu ----
            silo = new Visualizer::InputSilo();
            viz->addInputListener(silo);
            Array_<std::pair<String, int> > runMenuItems;
            runMenuItems.push_back(std::make_pair("Go", GoItem));
            runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
            viz->addMenu("Run", RunMenuId, runMenuItems);

            system.addEventReporter(new MyReporter(system, 1./120.));
            system.addEventReporter(new Visualizer::Reporter(*viz, 1./120.));
            if(blockExp!=NULL)
                //system.addEventReporter(
                //new MaxHeightReporter(system, *blockExp));
            if (spr1 != NULL)
                system.addEventReporter(
                    new ExpSprForceReporter(system, *spr3, 0.1));
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
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 12.0, 0.0));
            }
            // Spinning top
            if (condition == 5) {
                R.setRotationFromAngleAboutNonUnitVector(convertDegreesToRadians(54.74),
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
                blockCmp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, -10.0));
            }
        }
        // Exponential Contact
        // Same initial conditions as for Compliant Contact except shifted so
        // both blocks can be seen.
        if (ExpContactOn) {
            // Test modifying cooefficients of friction (states)
            //spr4->setMuStatic(state, -3.5);
            // Sitting Still
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
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 12.0, 0.0));
            }
            // Spinning top
            if (condition == 5) {
                R.setRotationFromAngleAboutNonUnitVector(convertDegreesToRadians(54.74),Vec3(1,0,1));
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(0.5, 0.2, 0.0));
                blockExp->setU(state, Vec6(0, 0, 0, 0, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 6.0, 0.0));
            }
            // Tumbling
            if (condition == 6) {
                blockExp->setQToFitRotation(state, R);
                blockExp->setQToFitTranslation(state, Vec3(-1.5, 1.0, 0.5));
                blockExp->setU(state, Vec6(0, 0, 0, 0.8, 0, 0));
                blockExp->setUToFitAngularVelocity(state, Vec3(0.0, 0.0, -10.0));
            }
            // Reset the spring zeros
            // The spring zeros are set to the point on the Floor that coincices
            // with the Station that was specified on the Body.
            if (spr1 != NULL) spr1->resetSpringZero(state);
            if (spr2 != NULL) spr2->resetSpringZero(state);
            if (spr3 != NULL) spr3->resetSpringZero(state);
            if (spr4 != NULL) spr4->resetSpringZero(state);
            if (spr5 != NULL) spr5->resetSpringZero(state);
            if (spr6 != NULL) spr6->resetSpringZero(state);
            if (spr7 != NULL) spr7->resetSpringZero(state);
            if (spr8 != NULL) spr8->resetSpringZero(state);
        }

        // Run Menu and storage for Replay
        if (viz != NULL) {
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
        cout << "Num steps attempted = " << integ.getNumStepsAttempted() << endl;
        cout << "    Num steps taken = " << integ.getNumStepsTaken() << endl;
        cout << " cpu time taken = " << cpuDuration << " sec" << endl;
        cout << "real time taken = " << realDuration << " sec" << endl;

        // Replay Loop
        if (viz != NULL) {
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
