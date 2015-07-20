/* -------------------------------------------------------------------------- *
 *            Simbody(tm) Adhoc Test: Cable Over Bicubic Surfaces             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Andreas Scholz                                   *
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

/*                     Simbody OpenSimPartyDemoCable
THIS DOESN'T WORK YET */

#include "Simbody.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// This gets called periodically to dump out interesting things about
// the cables and the system as a whole. It also saves states so that we
// can play back at the end.
static Array_<State> saveStates;
class ShowStuff : public PeriodicEventReporter {
public:
    ShowStuff(const MultibodySystem& mbs, 
              const CableSpring& cable1, Real interval) 
    :   PeriodicEventReporter(interval), 
        mbs(mbs), cable1(cable1) {}

    static void showHeading(std::ostream& o) {
        printf("%8s %10s %10s %10s %10s %10s %10s %10s %10s %12s\n",
            "time", "length", "rate", "integ-rate", "unitpow", "tension", "disswork",
            "KE", "PE", "KE+PE-W");
    }

    /** This is the implementation of the EventReporter virtual. **/ 
    void handleEvent(const State& state) const override {
        const CablePath& path1 = cable1.getCablePath();
        printf("%8g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %12.6g CPU=%g\n",
            state.getTime(),
            path1.getCableLength(state),
            path1.getCableLengthDot(state),
            path1.getIntegratedCableLengthDot(state),
            path1.calcCablePower(state, 1), // unit power
            cable1.getTension(state),
            cable1.getDissipatedEnergy(state),
            mbs.calcKineticEnergy(state),
            mbs.calcPotentialEnergy(state),
            mbs.calcEnergy(state)
                + cable1.getDissipatedEnergy(state),
            cpuTime());
        saveStates.push_back(state);
    }
private:
    const MultibodySystem&  mbs;
    CableSpring             cable1;
};

int main() {
  try {    
    // Create the system.   
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    matter.setShowDefaultGeometry(false);

    CableTrackerSubsystem cables(system);
    GeneralForceSubsystem forces(system);

    Force::Gravity gravity(forces, matter, -YAxis, 9.81);
    // Force::GlobalDamper(forces, matter, 5);

    system.setUseUniformBackground(true);    // no ground plane in display
    MobilizedBody Ground = matter.Ground(); // convenient abbreviation

    // Read in some bones.
    PolygonalMesh femur;
    PolygonalMesh tibia;

    femur.loadVtpFile("CableOverBicubicSurfaces-femur.vtp");
    tibia.loadVtpFile("CableOverBicubicSurfaces-tibia.vtp");
    femur.scaleMesh(30);
    tibia.scaleMesh(30);

    // Build a pendulum
    Body::Rigid pendulumBodyFemur(    MassProperties(1.0, Vec3(0, -5, 0), 
                                    UnitInertia(1).shiftFromCentroid(Vec3(0, 5, 0))));

    pendulumBodyFemur.addDecoration(Transform(), DecorativeMesh(femur).setColor(Vec3(0.8, 0.8, 0.8)));

    Body::Rigid pendulumBodyTibia(    MassProperties(1.0, Vec3(0, -5, 0), 
                                    UnitInertia(1).shiftFromCentroid(Vec3(0, 5, 0))));

    pendulumBodyTibia.addDecoration(Transform(), DecorativeMesh(tibia).setColor(Vec3(0.8, 0.8, 0.8)));

    Rotation z180(Pi, YAxis);

    MobilizedBody::Pin pendulumFemur(    matter.updGround(),
                                        Transform(Vec3(0, 0, 0)),
                                        pendulumBodyFemur,
                                        Transform(Vec3(0, 0, 0)) );

    Rotation rotZ45(-Pi/4, ZAxis);

    MobilizedBody::Pin pendulumTibia(   pendulumFemur,
                                        Transform(rotZ45, Vec3(0, -12, 0)),
                                        pendulumBodyTibia,
                                        Transform(Vec3(0, 0, 0)) );

    Real initialPendulumOffset = -0.25*Pi;

    Constraint::PrescribedMotion pres(matter, 
       new Function::Sinusoid(0.25*Pi, 0.2*Pi, 0*initialPendulumOffset), pendulumTibia, MobilizerQIndex(0));
               
    // Build a wrapping cable path
    CablePath path2(cables, Ground, Vec3(1, 3, 1),             // origin
                            pendulumTibia, Vec3(1, -4, 0));  // termination
    
    // Create a bicubic surface
    Vec3 patchOffset(0, -5, -1);
    Rotation rotZ90(0.5*Pi, ZAxis);
    Rotation rotX90(0.2*Pi, XAxis);

    Rotation patchRotation = rotZ90 * rotX90 * rotZ90;
    Transform patchTransform(patchRotation, patchOffset);

    Real patchScaleX = 2.0;
    Real patchScaleY = 2.0;
    Real patchScaleF = 0.75;

    const int Nx = 4, Ny = 4;
  
    const Real xData[Nx] = {  -2, -1, 1, 2 };
    const Real yData[Ny] = {  -2, -1, 1, 2 };

    const Real fData[Nx*Ny] = { 2,        3,        3,        1,
                                0,         1.5,  1.5,        0,
                                0,        1.5,  1.5,        0,
                                2,        3,        3,        1    };

    const Vector x_(Nx,        xData);
    const Vector y_(Ny,     yData);
    const Matrix f_(Nx, Ny, fData);

    Vector x = patchScaleX*x_;
    Vector y = patchScaleY*y_;
    Matrix f = patchScaleF*f_; 

    BicubicSurface patch(x, y, f, 0);

    Real highRes = 30;
    Real lowRes  = 1;

    PolygonalMesh highResPatchMesh = patch.createPolygonalMesh(highRes);
    PolygonalMesh lowResPatchMesh = patch.createPolygonalMesh(lowRes);

   
    pendulumFemur.addBodyDecoration(patchTransform,
        DecorativeMesh(highResPatchMesh).setColor(Cyan).setOpacity(.75));

    pendulumFemur.addBodyDecoration(patchTransform,
         DecorativeMesh(lowResPatchMesh).setRepresentation(DecorativeGeometry::DrawWireframe));

    Vec3 patchP(-0.5,-1,2), patchQ(-0.5,1,2);

    pendulumFemur.addBodyDecoration(patchTransform,
        DecorativePoint(patchP).setColor(Green).setScale(2));

    pendulumFemur.addBodyDecoration(patchTransform,
        DecorativePoint(patchQ).setColor(Red).setScale(2));

     CableObstacle::Surface patchObstacle(path2, pendulumFemur, patchTransform,
         ContactGeometry::SmoothHeightMap(patch));
        
      patchObstacle.setContactPointHints(patchP, patchQ);
    
      patchObstacle.setDisabledByDefault(true);

    // Sphere
    Real      sphRadius = 1.5;

    Vec3      sphOffset(0, -0.5, 0);
    Rotation  sphRotation(0*Pi, YAxis);
    Transform sphTransform(sphRotation, sphOffset);

    CableObstacle::Surface tibiaSphere(path2, pendulumTibia, sphTransform,
        ContactGeometry::Sphere(sphRadius));

    Vec3 sphP(1.5,-0.5,0), sphQ(1.5,0.5,0);
    tibiaSphere.setContactPointHints(sphP, sphQ);

    pendulumTibia.addBodyDecoration(sphTransform,
        DecorativeSphere(sphRadius).setColor(Red).setOpacity(0.5));

    // Make cable a spring
    CableSpring cable2(forces, path2, 50., 18., 0.1); 

    Visualizer viz(system);
    viz.setShowFrameNumber(true);
    system.adoptEventReporter(new Visualizer::Reporter(viz, 1./30));
    system.adoptEventReporter(new ShowStuff(system, cable2, 0.02));    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    system.realize(state, Stage::Position);
    viz.report(state);
    cout << "path2 init length=" << path2.getCableLength(state) << endl;
    cout << "Hit ENTER ...";
    getchar();

    // path1.setIntegratedCableLengthDot(state, path1.getCableLength(state));

    // Simulate it.
    saveStates.clear(); saveStates.reserve(2000);

    // RungeKutta3Integrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    // CPodesIntegrator integ(system);
    // integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-5);
    TimeStepper ts(integ);
    ts.initialize(state);
    ShowStuff::showHeading(cout);

    const Real finalTime = 10;
    const double startTime = realTime();
    ts.stepTo(finalTime);
    cout << "DONE with " << finalTime 
         << "s simulated in " << realTime()-startTime
         << "s elapsed.\n";

    while (true) {
        cout << "Hit ENTER FOR REPLAY, Q to quit ...";
        const char ch = getchar();
        if (ch=='q' || ch=='Q') break;
        for (unsigned i=0; i < saveStates.size(); ++i)
            viz.report(saveStates[i]);
    }

  } catch (const std::exception& e) {
    cout << "EXCEPTION: " << e.what() << "\n";
  }
}
