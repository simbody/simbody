/* -------------------------------------------------------------------------- *
 *                     Simbody(tm) Example: IMU Tracking                      *
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

/* This example shows a simple example using the Assembler for inverse
kinematics tracking where the input source is orientation observations from
an attached IMU (Inertial Measurement Unit).
*/


#include "Simbody.h"

#include <cstdio>
#include <exception>

using std::cout; using std::endl;

using namespace SimTK;

typedef SimTK::Markers::MarkerIx        MarkerIx;
typedef SimTK::Markers::ObservationIx   MarkerObsIx;

typedef SimTK::OrientationSensors::OSensorIx     IMUIx;
typedef SimTK::OrientationSensors::ObservationIx IMUObsIx;

int main() {
  try
  { // Create the system.

    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, -YAxis, 9.81);

    system.setUseUniformBackground(true);

    // An identity Rotation represents an IMU aligned with the body frame.
    const Rotation bodyAxesIMU;
    const Rotation tiltedIMU(Pi/4, ZAxis);

    const Rotation torsoIMUOri = bodyAxesIMU;
    const Rotation headIMUOri = tiltedIMU;


    const Real mass_torso = 20, mass_head = 2;
    const Vec3 hdims_torso(1,4,2), hdims_head(.3,1,.6);
    const Vec3 marker_head(hdims_head[0],hdims_head[1],0); // top front

    Body::Rigid torsoBody(MassProperties(mass_torso, Vec3(0),
                                            UnitInertia::brick(hdims_torso)));
    torsoBody.addDecoration(DecorativeBrick(hdims_torso).setOpacity(.4));
    torsoBody.addDecoration(torsoIMUOri,
        DecorativeFrame(1).setColor(Green).setLineThickness(5));

    Body::Rigid headBody(MassProperties(mass_head, Vec3(0),
                                        UnitInertia::brick(hdims_head)));
    headBody.addDecoration(DecorativeBrick(hdims_head)
                           .setOpacity(.4).setColor(Cyan));
    headBody.addDecoration(headIMUOri,
        DecorativeFrame(1).setColor(Green).setLineThickness(5));
    headBody.addDecoration(marker_head,
        DecorativePoint().setColor(Green));

    MobilizedBody::Free torso(matter.Ground(), Vec3(0,hdims_torso[1],0),
                              torsoBody, Vec3(0,-hdims_torso[1],0));

    MobilizedBody::Ball head(torso, Vec3(0,hdims_torso[1],0),
                             headBody, Vec3(0,-hdims_head[1],0));



    // Initialize the system and state.
    State state = system.realizeTopology();

    const Real Accuracy = 1e-5;
    Assembler ik(system);
    ik.setAccuracy(Accuracy);
    //ik.setForceNumericalGradient(true);


    Markers*            markers = new Markers();
    OrientationSensors* imus    = new OrientationSensors();
    ik.adoptAssemblyGoal(markers);
    ik.adoptAssemblyGoal(imus);

    const IMUIx    torsoIMU   = imus->addOSensor(torso, torsoIMUOri);
    const IMUIx    headIMU    = imus->addOSensor(head, headIMUOri);
    const MarkerIx headMarker = markers->addMarker(head, marker_head);

    ik.initialize(state);

    const IMUObsIx torsoObsIx = imus->getObservationIxForOSensor(torsoIMU);
    const IMUObsIx headObsIx  = imus->getObservationIxForOSensor(headIMU);
    const MarkerObsIx headMarkerObsIx =
        markers->getObservationIxForMarker(headMarker);

    // Try an initial assembly to an arbitrary pose.
    imus->moveOneObservation(torsoObsIx,
                             Rotation(SimTK::BodyRotationSequence,
                                      Pi/8, ZAxis, Pi/8, YAxis));
    imus->moveOneObservation(headObsIx,
                             Rotation()); // keep aligned with Ground
    markers->moveOneObservation(headMarkerObsIx,
                                Vec3(0, 12, 0));

    for (OrientationSensors::OSensorIx mx(0);
         mx < imus->getNumOSensors(); ++mx)
    {
        printf("mx=%d ox=%d err=%g\n",
            (int)mx, (int)imus->getObservationIxForOSensor(mx),
            imus->findCurrentOSensorError(mx));
    }

    Visualizer viz(system);
    // Show initial configuration
    viz.report(state);
    cout << "Initial state. Type any character to continue:\n";
    getchar();

    printf("Using accuracy=%g\n", ik.getAccuracyInUse());
    ik.assemble(state);

    for (OrientationSensors::OSensorIx mx(0);
         mx < imus->getNumOSensors(); ++mx)
    {
        printf("mx=%d ox=%d err=%g\n",
            (int)mx, (int)imus->getObservationIxForOSensor(mx),
            imus->findCurrentOSensorError(mx));
    }

    viz.report(state);
    cout << "ASSEMBLED CONFIGURATION\n";
    cout << "Type any character to start tracking:\n";
    getchar();

    const double startCPU = cpuTime(), startReal = realTime();

    const int NSteps = 200;
    for (int iters=0; iters <= NSteps; ++iters) {
        const Real slow = std::sin(2*Pi*iters/100);
        const Real fast = std::sin(10*Pi*iters/100);
        Rotation torsoObs(SpaceRotationSequence,
                          (Pi/4)*slow, ZAxis,
                          (Pi/2)*slow, YAxis);
        Rotation headObs((Pi/8)*fast, YAxis); // shake head
        Vec3 markerObs(Vec3(0,12,0) + slow*Vec3(5,0,0));

        imus->moveOneObservation(torsoObsIx, torsoObs);
        imus->moveOneObservation(headObsIx, headObs);
        markers->moveOneObservation(headMarkerObsIx, markerObs);

        ik.track();
        ik.updateFromInternalState(state);
        viz.report(state);
    }

    cout << "TRACKED " << NSteps << " steps in "
         << cpuTime()-startCPU   << " CPU s, "
         << realTime()-startReal << " REAL s\n";

    cout << "DONE TRACKING ...\n";
    viz.report(state);

    cout << "Type any character to continue:\n";
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
