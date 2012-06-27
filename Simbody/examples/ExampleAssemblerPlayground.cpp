/* -------------------------------------------------------------------------- *
 *                 Simbody(tm) Example: Assembler Playground                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
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

/* This example is for experimenting with the Assembler study that is used
for initial assembly and for inverse kinematics (IK) such as repeated tracking
of marker positions. This toy example is easy to understand and captures all
the basics in a simplified model. For a realistic biomechanical model of about
the same size but considerably more complexity, see the AmysIKProblem example.
*/


#include "Simbody.h"

#include <cstdio>
#include <exception>

using std::cout; using std::endl;

using namespace SimTK;

static const int NBodies = 10;
static const int NBodies2 = 7;

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    const Real mass = 0.1;
    const Vec3 hdims(2,4,1);
    Body::Rigid pendulumBody = Body::Rigid
       (MassProperties(mass, Vec3(0), mass*UnitInertia::ellipsoid(hdims)))
        .addDecoration(Transform(),DecorativeEllipsoid(hdims)
                                    .setOpacity(.2));

    MobilizedBody midBody, finalBody;
    MobilizedBody::Free baseBody(matter.Ground(), Transform(Vec3(0,-hdims[1],0)), 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));

    Force::TwoPointLinearSpring(forces, matter.Ground(), Vec3(0),
        baseBody, Vec3(0), 2, 10);

    Force::GlobalDamper(forces, matter, 1);


    MobilizedBody parent = baseBody;
    for (int i=0; i < NBodies; ++i) {
        MobilizedBody::Ball child(parent, Transform(Vec3(0,-hdims[1],0)), 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));
        if (i==NBodies/2) midBody   = child;
        if (i==NBodies-1) finalBody = child;
        parent = child;
    }

    // Second chain
    MobilizedBody endOfSecondChain;
    parent = matter.Ground();
    Vec3 loc(-10,0,0);
    for (int i=0; i < NBodies2; ++i) {
        MobilizedBody::Ball child(parent, loc, 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));
        if (i==NBodies2-1) endOfSecondChain = child;
        parent = child;
        loc = Vec3(0,-hdims[1],0);
    }

    Force::TwoPointLinearSpring(forces, endOfSecondChain, Vec3(0),
                                midBody, Vec3(0), 2, 1);

    Constraint pres[5];
    for (int i=0; i<5; ++i) {
        pres[i] = Constraint::PrescribedMotion(matter, 
            new Function::Constant(.1),
            MobilizedBodyIndex(NBodies+2+i), MobilizerQIndex(2));
        //pres[i].setDisabledByDefault(true);
    }

    //Constraint::Ball chainConnection(endOfSecondChain, Vec3(0,-hdims[1],0),
    //    matter.updMobilizedBody(MobilizedBodyIndex(4)), Vec3(0));

    //Constraint::Ball ball(matter.Ground(), Vec3(3*NBodies,-NBodies,0), 
    //                      finalBody, Vec3(0,-hdims[1],0));

    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.1));


    // Initialize the system and state.
    
    system.realizeTopology();

    State state = system.getDefaultState();
    //for (int i=0; i<5; ++i) pres[i].disable(state);

    for (int i=0; i<5; ++i)
        cout << "state pres[" << i << "].isDisabled=" 
             << pres[i].isDisabled(state) << endl;

    const Vec3 midTarget(1.5*NBodies,-NBodies*0.5*hdims[1],3);
    const Vec3 finalTarget(2*NBodies, -NBodies*0.3*hdims[1],2);
    viz.addDecoration(midBody, Vec3(0), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(midBody, Vec3(1,2,3), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(GroundIndex, midTarget, 
        DecorativeSphere(1).setColor(Green));
    viz.addDecoration(finalBody, Vec3(0), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(finalBody, Vec3(1,2,3), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(GroundIndex, finalTarget, 
        DecorativeSphere(1).setColor(Green));

    // Show initial configuration
    viz.report(state);
    State tempState = state; 
    for (int i=0; i<5; ++i)
        cout << "tempState copy pres[" << i << "].isDisabled=" 
             << pres[i].isDisabled(tempState) << endl;

    system.realize(tempState, Stage::Position);
    const int nQuat = matter.getNumQuaternionsInUse(tempState);
    cout << "INITIAL CONFIGURATION\n"; 
    cout << tempState.getNU() << " dofs, " 
         << tempState.getNQErr()-nQuat << " constraints.\n";
    
    cout << "Type any character to continue:\n";
    getchar();


    Assembler ik(system);
    Markers& markers = *new Markers();
    ik.adoptAssemblyGoal(&markers);

    QValue& qvalue = *new QValue(endOfSecondChain, MobilizerQIndex(2),
        Pi/2);
    ik.adoptAssemblyGoal(&qvalue);

    // Set treatment of constraints.
    //ik.setSystemConstraintsWeight(1); // means penalty rather than constraint

    ik.lockMobilizer(MobilizedBodyIndex(2));
    ik.lockMobilizer(MobilizedBodyIndex(3));
    ik.lockMobilizer(MobilizedBodyIndex(4));

    ik.restrictQ(midBody, MobilizerQIndex(0), -Pi/10, Pi/10);

    markers.addMarker(midBody, Vec3(0));
    markers.addMarker(midBody, Vec3(1,2,3));
    markers.addMarker(finalBody, Vec3(0));
    markers.addMarker(finalBody, Vec3(1,2,3));

    // Manual observation/marker correspondence.
    //Array_<Markers::MarkerIx> observationOrder;
    //observationOrder.push_back(Markers::MarkerIx(2));
    //observationOrder.push_back(Markers::MarkerIx(3));
    //observationOrder.push_back(Markers::MarkerIx(0));
    //observationOrder.push_back(Markers::MarkerIx(1));
    //observationOrder.push_back(Markers::MarkerIx()); // unused
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //markers.defineTargetOrder(observationOrder);


    const Real Accuracy = 1e-3;
    ik.setAccuracy(Accuracy);

    ik.initialize(state);

    markers.moveOneObservation(Markers::ObservationIx(0), midTarget);
    markers.moveOneObservation(Markers::ObservationIx(1), midTarget);
    markers.moveOneObservation(Markers::ObservationIx(2), finalTarget);
    markers.moveOneObservation(Markers::ObservationIx(3), finalTarget);

    for (Markers::MarkerIx mx(0); mx < markers.getNumMarkers(); ++mx)
        printf("mx=%d ox=%d err=%g\n", 
            (int)mx, (int)markers.getObservationIxForMarker(mx),
            markers.findCurrentMarkerError(mx));


    ik.assemble(state);

    viz.report(state);
    cout << "ASSEMBLED CONFIGURATION\n"; 
    
    cout << "Type any character to continue:\n";
    getchar();

    for (Markers::MarkerIx mx(0); mx < markers.getNumMarkers(); ++mx)
        printf("mx=%d ox=%d err=%g\n", 
            (int)mx, (int)markers.getObservationIxForMarker(mx),
            markers.findCurrentMarkerError(mx));

    State startState = state;

    const double startCPU = cpuTime(), startReal = realTime();

    const int NSteps = 100;
    const int NToSkip = 4;
    for (int iters=0; iters <= NSteps; ++iters) {
        Vec3 newMidTarget = midTarget + NBodies/2*std::sin(2*Pi*iters/100);
        Vec3 newFinalTarget = finalTarget - NBodies/3*std::sin(2*Pi*iters/100);
        markers.moveOneObservation(Markers::ObservationIx(0), newMidTarget);
        markers.moveOneObservation(Markers::ObservationIx(1), newMidTarget);
        markers.moveOneObservation(Markers::ObservationIx(2), newFinalTarget);
        markers.moveOneObservation(Markers::ObservationIx(3), newFinalTarget);

        qvalue.setValue(qvalue.getValue() + Pi/10);
                                        
        ik.track();
        if (iters%NToSkip == 0) {
            ik.updateFromInternalState(state);
            viz.report(state);
        }
    }

    cout << "ASSEMBLED " << NSteps << " steps in " 
         << cpuTime()-startCPU   << " CPU s, " 
         << realTime()-startReal << " REAL s\n";

    cout << "DONE ASSEMBLING -- SIMULATE ...\n";
    viz.report(state);

    cout << "Type any character to continue:\n";
    getchar();
   
    // Simulate it.
    for (int i=0; i<5; ++i)
        cout << "state pres[" << i << "].isDisabled=" 
             << pres[i].isDisabled(state) << endl;

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);

    cout << "DONE SIMULATING.\n";
    cout << "Type any character to quit:\n";
    getchar();


  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}
