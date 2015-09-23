/* -------------------------------------------------------------------------- *
 *                       Simbody(tm) Example: Cable Path                      *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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

/*                      Simbody ExampleCablePath
This example shows how to use a CableTrackerSubsystem to follow the motion of
a cable that connects two bodies and passes around obstacles. We'll then
create a force element that generates spring forces that result from the
stretching and stretching rate of the cable. */

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
    CableTrackerSubsystem cables(system);
    GeneralForceSubsystem forces(system);

    system.setUseUniformBackground(true); // no ground plane in display

    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    //Force::GlobalDamper(forces, matter, 5);

    Body::Rigid someBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    const Real Rad = .25;
    someBody.addDecoration(Transform(),
        DecorativeSphere(Rad).setOpacity(.75).setResolution(4));

    Body::Rigid biggerBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    const Real BiggerRad = .5;
    biggerBody.addDecoration(Transform(),
        DecorativeSphere(BiggerRad).setOpacity(.75).setResolution(4));

    const Vec3 radii(.3, .25, .25);
    Body::Rigid ellipsoidBody(MassProperties(1.0, Vec3(0),
        UnitInertia::ellipsoid(radii)));
    ellipsoidBody.addDecoration(Transform(),
        DecorativeEllipsoid(radii).setOpacity(.75).setResolution(4)
                                  .setColor(Orange));

    const Real CylRad = .25, HalfLen = .5;
    Body::Rigid cylinderBody(MassProperties(1.0, Vec3(0),
        1.*UnitInertia::cylinderAlongX(Rad,HalfLen)));
    cylinderBody.addDecoration(Rotation(-Pi/2,ZAxis),
        DecorativeCylinder(CylRad,HalfLen).setOpacity(.75)
           .setResolution(4).setColor(Orange));

    Body::Rigid fancyBody = biggerBody; // NOT USING ELLIPSOID

    MobilizedBody Ground = matter.Ground();

    MobilizedBody::Free body1(Ground,           Transform(Vec3(0)),
                              someBody,         Transform(Vec3(0, 1, 0)));
    //MobilizedBody::Ball body2(body1,            Transform(Vec3(0)),
    //                          someBody,         Transform(Vec3(0, 1, 0)));
    //MobilizedBody::Ball body3(body2,            Transform(Vec3(0)),
    //                          someBody,         Transform(Vec3(0, 1, 0)));
    //MobilizedBody::Ball body4(body3,            Transform(Vec3(0)),
    //                          fancyBody,    Transform(Vec3(0, 1, 0)));
    //MobilizedBody::Ball body5(body4,            Transform(Vec3(0)),
    //                          someBody,         Transform(Vec3(0, 1, 0)));

    CablePath path1(cables, Ground, Vec3(0,0,0),   // origin
                            body1, Vec3(0,Rad,0));  // termination

    //CableObstacle::Surface obstacle(path1, Ground, Vec3(1,-1,0),
    //    ContactGeometry::Sphere(Rad));
    //obstacle.setContactPointHints(Rad*UnitVec3(1,1,0),Rad*UnitVec3(1,.5,0));
    //
    CableObstacle::Surface obstacle2(path1, Ground, Vec3(0,-1,0),
        //ContactGeometry::Sphere(CylRad));
        ContactGeometry::Ellipsoid(radii));
        //ContactGeometry::Torus(CylRad,.1));
        //ContactGeometry::Cylinder(CylRad));

    obstacle2.setContactPointHints(1.5*CylRad*UnitVec3(1,1,0),1.5*CylRad*UnitVec3(1,.5,0));
    obstacle2.setDisabledByDefault(true);

    ////CableObstacle::ViaPoint p4(path1, body4, Rad*UnitVec3(0,1,1));
    ////CableObstacle::ViaPoint p5(path1, body4, Rad*UnitVec3(1,0,1));
    //CableObstacle::Surface obs5(path1, body4,
    //    // Transform(), ContactGeometry::Ellipsoid(radii));
    //    //Rotation(Pi/2, YAxis), ContactGeometry::Cylinder(CylRad)); // along y
    //    //Transform(), ContactGeometry::Sphere(Rad));
    //    Transform(), ContactGeometry::Sphere(BiggerRad));
    ////obs5.setContactPointHints(Rad*UnitVec3(0,-1,-1),Rad*UnitVec3(0.1,-1,-1));
    //obs5.setContactPointHints(Rad*UnitVec3(.1,.125,-.2),
    //                          Rad*UnitVec3(0.1,-.1,-.2));

    CableSpring cable1(forces, path1, .5*100., 2, 0.1);

    //obs1.setPathPreferencePoint(Vec3(2,3,4));
    //obs1.setDecorativeGeometry(DecorativeSphere(0.25).setOpacity(.5));

    Visualizer viz(system);
    viz.setShowFrameNumber(true);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.1*1./30));
    system.addEventReporter(new ShowStuff(system, cable1, 0.1*0.1));
    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    //Random::Gaussian random;
    //for (int i = 0; i < state.getNQ(); ++i)
    //    state.updQ()[i] = random.getValue();
    //for (int i = 0; i < state.getNU(); ++i)
    //    state.updU()[i] = 0.1*random.getValue();
    body1.setQToFitTranslation(state, Vec3(2,-1,.5));

    system.realize(state, Stage::Position);
    viz.report(state);
    cout << "path1 init length=" << path1.getCableLength(state) << endl;
    cout << "Hit ENTER ...";
    getchar();

    path1.setIntegratedCableLengthDot(state, path1.getCableLength(state));


    // Simulate it.
    saveStates.clear(); saveStates.reserve(2000);

    //RungeKutta3Integrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //CPodesIntegrator integ(system);
    //integ.setAllowInterpolation(false);
    integ.setAccuracy(1e-3);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ShowStuff::showHeading(cout);

    const Real finalTime = 2;
    const double startTime = realTime(), startCPU = cpuTime();
    ts.stepTo(finalTime);
    cout << "DONE with " << finalTime
         << "s simulated in " << realTime()-startTime
         << "s elapsed, " << cpuTime()-startCPU << "s CPU.\n";


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
