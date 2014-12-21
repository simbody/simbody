/* -------------------------------------------------------------------------- *
 *                  Simbody(tm) Example: Bricard Mechanism                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Moonki Jung                                                  *
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
#include "shared/SimbodyExampleHelper.h"

#include <fstream>
#include <iostream>

using std::cout; using std::endl;


using namespace SimTK;

/* Note: this example is due to Moonki Jung from Seoul National University's
Human-Centered CAD Laboratory. See SimTK Core open-discussion forum startin
g 2/4/2010 with message 3292. Moonki kindly gave us permission to include this 
very difficult mechanism as a test case for Simbody and as a user example. The
Bricard mechanism is particularly difficult because it has six mobilities and 
six constraints but has one net degree of freedom. Thus one of the constraints
is redundant; but it is not any particular one. The mechanism must be perfectly
aligned in order to move, but in practice we can only satisfy the position
constraints approximately. This can hide the redundancy, causing the mechanism
to appear to have zero dofs and thus lock up. */


class EnergyReport : public PeriodicEventReporter {
public:
    EnergyReport(const MultibodySystem& system, Real interval) 
    :   PeriodicEventReporter(interval), system(system) {}

    void handleEvent(const State& s) const {
        cout << "\n*** t=" << s.getTime() 
             << " ke=" << system.calcKineticEnergy(s) 
             << " E=" << system.calcEnergy(s)
             << endl;
        cout << "    u=" << s.getU() << "\n";
        cout << " uerr=" << s.getUErr() << "\n\n";
    }
private:
    const MultibodySystem& system;
};


int main()
{
  try {
    cout << "This is Simbody example '" 
         << SimbodyExampleHelper::getExampleName() << "'\n";
    cout << "Working dir=" << Pathname::getCurrentWorkingDirectory() << endl;

    const std::string auxDir = 
        SimbodyExampleHelper::findAuxiliaryDirectoryContaining
        ("geometry/Bricard_EVEN_PART.obj");
    std::cout << "Getting geometry from '" << auxDir << "'\n";

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, UnitVec3(0, -1, 0), 9.8);

    const Real Mass = 20;
    const Vec3 EvenCOM(1.00000000, -0.16416667, -0.16416667);
    const Vec3 OddCOM(1.00000000, 0.16416667, -0.16416667);

    const Inertia EvenCentralInertia
       ( 3.33400000, 28.33366667, 28.33366667,  // xx, yy, zz
        -4.96666667, -1.60000000, -0.03333333); // xy, xz, yz

    const Inertia OddCentralInertia
       ( 3.33400000, 28.33366667, 28.33366667,  // xx, yy, zz 
         4.96666667, -1.60000000, 0.03333333);  // xy, xz, yz

    // Inertias must be given about the body origin.
    const Inertia EvenBodyInertia = 
        EvenCentralInertia.shiftFromMassCenter(-EvenCOM, Mass);
    const Inertia OddBodyInertia = 
        OddCentralInertia.shiftFromMassCenter(-OddCOM, Mass);

    Body::Rigid EVEN_PART_1(MassProperties(Mass, EvenCOM, EvenBodyInertia));
    Body::Rigid EVEN_PART_2(MassProperties(Mass, EvenCOM, EvenBodyInertia));
    Body::Rigid EVEN_PART_3(MassProperties(Mass, EvenCOM, EvenBodyInertia));


    Body::Rigid ODD_PART_1(MassProperties(Mass, OddCOM, OddBodyInertia));
    Body::Rigid ODD_PART_2(MassProperties(Mass, OddCOM, OddBodyInertia));

    // Split the last body and weld back together to close loop.
    Body::Rigid ODD_PART_3_HALF1(MassProperties(Mass/2, OddCOM, 
                                                OddBodyInertia/2));
    Body::Rigid ODD_PART_3_HALF2(MassProperties(Mass/2, OddCOM, 
                                                OddBodyInertia/2));

    PolygonalMesh Mesh1, Mesh2; 
    Mesh1.loadObjFile(auxDir + "geometry/Bricard_EVEN_PART.obj"); 
    Mesh2.loadObjFile(auxDir + "geometry/Bricard_ODD_PART.obj"); 

    EVEN_PART_1.addDecoration(Transform(), DecorativeMesh(Mesh1).setColor(Vec3(0.00000000, 1.00000000, 0.00000000)));
    EVEN_PART_2.addDecoration(Transform(), DecorativeMesh(Mesh1).setColor(Vec3(1.00000000, 0.00000000, 1.00000000)));
    EVEN_PART_3.addDecoration(Transform(), DecorativeMesh(Mesh1).setColor(Vec3(1.00000000, 1.00000000, 0.00000000)));
    ODD_PART_1.addDecoration(Transform(), DecorativeMesh(Mesh2).setColor(Vec3(1.00000000, 0.00000000, 0.00000000)));
    ODD_PART_2.addDecoration(Transform(), DecorativeMesh(Mesh2).setColor(Vec3(0.00000000, 0.00000000, 1.00000000)));
    ODD_PART_3_HALF1.addDecoration(Transform(), DecorativeMesh(Mesh2).setColor(Vec3(0.00000000, 1.00000000, 1.00000000)));
    
    MobilizedBody::Weld EVEN_PART_1_body(matter.updGround(), Transform(Rotation(Mat33(1,0,0,0,-1,0,0,0,-1)), Vec3(0, 0, 0))
        ,EVEN_PART_1, Transform());
    MobilizedBody::Pin ODD_PART_1_body(EVEN_PART_1_body, Transform(Rotation(Mat33(0,-1,0,1,0,0,0,0,1)), Vec3(0, 0, 0))
        ,ODD_PART_1, Transform(Rotation(Mat33(0,-1,0,1,0,0,0,0,1)), Vec3(0,0,0)));

    MobilizedBody::Pin EVEN_PART_2_body(ODD_PART_1_body, Transform(Rotation(Mat33(0,-1,0,0,0,1,-1,0,0)), Vec3(2, 0, 0))
        ,EVEN_PART_2, Transform(Rotation(Mat33(0,-1,0,-1,0,0,0,0,-1)), Vec3(0,0,0)));

    MobilizedBody::Pin ODD_PART_2_body(EVEN_PART_1_body, Transform(Rotation(Mat33(0,-1,0,0,0,1,-1,0,0)), Vec3(2, 0, 0))
        ,ODD_PART_2, Transform(Rotation(Mat33(0,-1,0,1,0,0,0,0,1)), Vec3(0,0,0)));

    MobilizedBody::Pin EVEN_PART_3_body(ODD_PART_2_body, Transform(Rotation(Mat33(0,-1,0,0,0,1,-1,0,0)), Vec3(2, 0, 0))
        ,EVEN_PART_3, Transform(Rotation(Mat33(0,-1,0,0,0,-1,1,0,0)), Vec3(2,0,0)));

    MobilizedBody::Pin ODD_PART_3_HALF1_body(EVEN_PART_3_body, Transform(Rotation(Mat33(0,-1,0,-1,0,0,0,0,-1)), Vec3(0, 0, 0))
        ,ODD_PART_3_HALF1, Transform(Rotation(Mat33(0,-1,0,0,0,1,-1,0,0)), Vec3(2,0,0)));

    MobilizedBody::Pin ODD_PART_3_HALF2_body(EVEN_PART_2_body, Transform(Rotation(Mat33(0,-1,0,0,0,1,-1,0,0)), Vec3(2, 0, 0))
        ,ODD_PART_3_HALF2, Transform(Rotation(Mat33(0,-1,0,1,0,0,0,0,1)), Vec3(0,0,0)));

    Constraint::Weld ODD_PART_3_UNION(ODD_PART_3_HALF1_body, Transform(), ODD_PART_3_HALF2_body, Transform());

    //Constraint::ConstantSpeed motion(EVEN_PART_3_body, -.1);
    //Force::MobilityLinearSpring frc(forces, EVEN_PART_3_body, 
       // MobilizerUIndex(0), 100, 0);

    Visualizer viz(system);
    viz.setCameraTransform(Vec3(0.5,0.5,0.5));
    viz.pointCameraAt(Vec3(0), Vec3(0,1,0));
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));

    system.addEventReporter(new EnergyReport(system, .01));
    system.realizeTopology();
    State state = system.getDefaultState();

    // Set initial states (Q's and U's)
    // Position
    ODD_PART_1_body.setOneQ(state, 0, 180.0*Pi/180.0);
    EVEN_PART_3_body.setOneQ(state, 0, 180.0*Pi/180.0);
    ODD_PART_3_HALF2_body.setOneQ(state, 0, 0.0*Pi/180.0);

    EVEN_PART_2_body.setOneQ(state, 0, -120.0*Pi/180.0);
    ODD_PART_2_body.setOneQ(state, 0, -120.0*Pi/180.0);
    ODD_PART_3_HALF1_body.setOneQ(state, 0, 120.0*Pi/180.0);

    // Velocity
    ODD_PART_1_body.setOneU(state,0, -11.2);

    //RungeKuttaMersonIntegrator integ(system);
    RungeKutta3Integrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //VerletIntegrator integ(system);
    //CPodesIntegrator integ(system);

    // Accuracy needs to be fairly tight to avoid lockup that would occur
    // if the constraints were allowed to drift.
    integ.setAccuracy(1e-5);
    //integ.setConstraintTolerance(1e-8);
    integ.initialize(state);


    const double startCPU = cpuTime(), startReal = realTime();

    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(20.0);    

    cout << "DONE. CPU=" << cpuTime()-startCPU 
         << "s, REAL=" << realTime()-startReal << "s\n";

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());


  } catch (const std::exception& e) {
      std::cout << "std::exception: " << e.what() << std::endl;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
  }

  return 0;
}

