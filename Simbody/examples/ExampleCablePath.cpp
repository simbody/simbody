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
#include "simbody/internal/CableTrackerSubsystem.h"
#include "simbody/internal/CablePath.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// This force element implements an elastic cable of a given nominal length,
// and a stiffness k that generates a k*x force opposing stretch beyond
// nominal. There is also a damping term c*xdot that applies only when the
// cable is stretched and is being extended (x>0 && xdot>0). We keep track
// of dissipated power here so we can use conservation of energy to check that
// the cable and force element aren't obviously broken.
class MyCableSpringImpl : public Force::Custom::Implementation {
public:
    MyCableSpringImpl(const GeneralForceSubsystem& forces, 
                      const CablePath& path, 
                      Real stiffness, Real nominal, Real damping) 
    :   forces(forces), path(path), k(stiffness), x0(nominal), c(damping)
    {   assert(stiffness >= 0 && nominal >= 0 && damping >= 0); }

    const CablePath& getCablePath() const {return path;}

    // Must be at stage Velocity. Evalutes tension if necessary.
    Real getTension(const State& state) const {
        ensureTensionCalculated(state);
        return Value<Real>::downcast(forces.getCacheEntry(state, tensionx));
    }

    // Must be at stage Velocity.
    Real getPowerDissipation(const State& state) const {
        const Real stretch = calcStretch(state);
        if (stretch == 0) return 0;
        const Real rate = path.getCableLengthDot(state);
        return k*stretch*std::max(c*rate, -1.)*rate;
    }

    // This integral is always available.
    Real getDissipatedEnergy(const State& state) const {
        return forces.getZ(state)[workx];
    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals

    // Ask the cable to apply body forces given the tension calculated here.
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11
    {   path.applyBodyForces(state, getTension(state), bodyForces); }

    // Return the potential energy currently stored by the stretch of the cable.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {
        const Real stretch = calcStretch(state);
        if (stretch == 0) return 0;
        return k*square(stretch)/2;
    }

    // Allocate the state variable for tracking dissipated energy, and a
    // cache entry to hold the calculated tension.
    void realizeTopology(State& state) const OVERRIDE_11 {
        Vector initWork(1, 0.);
        workx = forces.allocateZ(state, initWork);
        tensionx = forces.allocateLazyCacheEntry(state, Stage::Velocity,
                                             new Value<Real>(NaN));
    }

    // Report power dissipation as the derivative for the work variable.
    void realizeAcceleration(const State& state) const OVERRIDE_11 {
        Real& workDot = forces.updZDot(state)[workx];
        workDot = getPowerDissipation(state);
    }
    //--------------------------------------------------------------------------

private:
    // Return the amount by which the cable is stretched beyond its nominal
    // length or zero if the cable is slack. Must be at stage Position.
    Real calcStretch(const State& state) const {
        const Real stretch = path.getCableLength(state) - x0;
        return std::max(stretch, 0.);
    }

    // Must be at stage Velocity to calculate tension.
    Real calcTension(const State& state) const {
        const Real stretch = calcStretch(state);
        if (stretch == 0) return 0;
        const Real rate = path.getCableLengthDot(state);
        if (c*rate < -1)
            cout << "c*rate=" << c*rate << "; limited to -1\n";
        const Real tension = k*stretch*(1+std::max(c*rate,-1.));
        return tension;
    }

    // If state is at stage Velocity, we can calculate and store tension
    // in the cache if it hasn't already been calculated.
    void ensureTensionCalculated(const State& state) const {
        if (forces.isCacheValueRealized(state, tensionx))
            return;
        Value<Real>::updDowncast(forces.updCacheEntry(state, tensionx)) 
            = calcTension(state);
        forces.markCacheValueRealized(state, tensionx);
    }

    const GeneralForceSubsystem&    forces;
    CablePath                       path;
    Real                            k, x0, c;
    mutable ZIndex                  workx;
    mutable CacheEntryIndex         tensionx;
};

// A nice handle to hide most of the cable spring implementation. This defines
// a user's API.
class MyCableSpring : public Force::Custom {
public:
    MyCableSpring(GeneralForceSubsystem& forces, const CablePath& path, 
                  Real stiffness, Real nominal, Real damping) 
    :   Force::Custom(forces, new MyCableSpringImpl(forces,path,
                                                    stiffness,nominal,damping)) 
    {}
    
    // Expose some useful methods.
    const CablePath& getCablePath() const 
    {   return getImpl().getCablePath(); }
    Real getTension(const State& state) const
    {   return getImpl().getTension(state); }
    Real getPowerDissipation(const State& state) const
    {   return getImpl().getPowerDissipation(state); }
    Real getDissipatedEnergy(const State& state) const
    {   return getImpl().getDissipatedEnergy(state); }

private:
    const MyCableSpringImpl& getImpl() const
    {   return dynamic_cast<const MyCableSpringImpl&>(getImplementation()); }
};

static Array_<State> saveStates;
// This gets called periodically to dump out interesting things about
// the cables and the system as a whole.
class ShowStuff : public PeriodicEventReporter {
public:
    ShowStuff(const MultibodySystem& mbs, 
              const MyCableSpring& cable1, 
              const MyCableSpring& cable2, Real interval) 
    :   PeriodicEventReporter(interval), 
        mbs(mbs), cable1(cable1), cable2(cable2) {}

    static void showHeading(std::ostream& o) {
        printf("%8s %10s %10s %10s %10s %10s %10s %10s %10s %12s\n",
            "time", "length", "rate", "integ-rate", "unitpow", "tension", "disswork",
            "KE", "PE", "KE+PE-W");
    }

    /** This is the implementation of the EventReporter virtual. **/ 
    void handleEvent(const State& state) const OVERRIDE_11 {
        const CablePath& path1 = cable1.getCablePath();
        const CablePath& path2 = cable2.getCablePath();
        printf("%8g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g CPU=%g\n",
            state.getTime(),
            path1.getCableLength(state),
            path1.getCableLengthDot(state),
            path1.getIntegratedCableLengthDot(state),
            path1.calcCablePower(state, 1), // unit power
            cable1.getTension(state),
            cable1.getDissipatedEnergy(state), 
            cpuTime());
        printf("%8s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %12.6g\n",
            "",
            path2.getCableLength(state),
            path2.getCableLengthDot(state),
            path2.getIntegratedCableLengthDot(state),
            path2.calcCablePower(state, 1), // unit power
            cable2.getTension(state),
            cable2.getDissipatedEnergy(state),
            mbs.calcKineticEnergy(state),
            mbs.calcPotentialEnergy(state),
            mbs.calcEnergy(state)
                + cable1.getDissipatedEnergy(state)
                + cable2.getDissipatedEnergy(state));
        saveStates.push_back(state);
    }
private:
    const MultibodySystem&  mbs;
    MyCableSpring           cable1, cable2;
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

    const Vec3 radii(.4, .25, .15);
    Body::Rigid ellipsoidBody(MassProperties(1.0, Vec3(0), 
        1.*UnitInertia::ellipsoid(radii)));
    ellipsoidBody.addDecoration(Transform(), 
        DecorativeEllipsoid(radii).setOpacity(.75).setResolution(4)
                                  .setColor(Orange));

    const Real CylRad = .3, HalfLen = .5;
    Body::Rigid cylinderBody(MassProperties(1.0, Vec3(0), 
        1.*UnitInertia::cylinderAlongX(Rad,HalfLen)));
    cylinderBody.addDecoration(Rotation(-Pi/2,ZAxis), 
        DecorativeCylinder(CylRad,HalfLen).setOpacity(.75)
           .setResolution(4).setColor(Orange));

    Body::Rigid fancyBody = biggerBody; // NOT USING ELLIPSOID

    MobilizedBody Ground = matter.Ground();

    MobilizedBody::Ball body1(Ground,           Transform(Vec3(0)), 
                              someBody,         Transform(Vec3(0, 1, 0)));
    MobilizedBody::Ball body2(body1,            Transform(Vec3(0)), 
                              someBody,         Transform(Vec3(0, 1, 0)));
    MobilizedBody::Ball body3(body2,            Transform(Vec3(0)), 
                              someBody,         Transform(Vec3(0, 1, 0)));
    MobilizedBody::Ball body4(body3,            Transform(Vec3(0)), 
                              fancyBody,    Transform(Vec3(0, 1, 0)));
    MobilizedBody::Ball body5(body4,            Transform(Vec3(0)), 
                              someBody,         Transform(Vec3(0, 1, 0)));

    CablePath path1(cables, body1, Vec3(Rad,0,0),   // origin
                            body5, Vec3(0,0,Rad));  // termination

    CableObstacle::ViaPoint p1(path1, body2, Rad*UnitVec3(1,1,0));
    //CableObstacle::ViaPoint p2(path1, body3, Rad*UnitVec3(0,1,1));
    //CableObstacle::ViaPoint p3(path1, body3, Rad*UnitVec3(1,0,1));
    CableObstacle::Surface obs4(path1, body3, Transform(), 
        ContactGeometry::Sphere(Rad));
    //obs4.setContactPointHints(Rad*UnitVec3(-1,1,0),Rad*UnitVec3(-1,0,1));
    obs4.setContactPointHints(Rad*UnitVec3(-.25,.04,0.08),
                              Rad*UnitVec3(-.05,-.25,-.04));

    //CableObstacle::ViaPoint p4(path1, body4, Rad*UnitVec3(0,1,1));
    //CableObstacle::ViaPoint p5(path1, body4, Rad*UnitVec3(1,0,1));
    CableObstacle::Surface obs5(path1, body4, 
        // Transform(), ContactGeometry::Ellipsoid(radii));
        //Rotation(Pi/2, YAxis), ContactGeometry::Cylinder(CylRad)); // along y
        //Transform(), ContactGeometry::Sphere(Rad));
        Transform(), ContactGeometry::Sphere(BiggerRad));
    //obs5.setContactPointHints(Rad*UnitVec3(0,-1,-1),Rad*UnitVec3(0.1,-1,-1));
    obs5.setContactPointHints(Rad*UnitVec3(.1,.125,-.2),
                              Rad*UnitVec3(0.1,-.1,-.2));

    // NOTE: velocity-based force is disabled.
    MyCableSpring cable1(forces, path1, 100., 3.5, 0*0.1); 

    CablePath path2(cables, Ground, Vec3(-3,0,0),   // origin
                            Ground, Vec3(-2,1,0)); // termination
    CableObstacle::ViaPoint(path2, body3, 2*Rad*UnitVec3(1,1,1));
    CableObstacle::ViaPoint(path2, Ground, Vec3(-2.5,1,0));
    MyCableSpring cable2(forces, path2, 100., 8, 0.1); 


    //obs1.setPathPreferencePoint(Vec3(2,3,4));
    //obs1.setDecorativeGeometry(DecorativeSphere(0.25).setOpacity(.5));

    Visualizer viz(system);
    viz.setShowFrameNumber(true);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.1*1./30));
    system.addEventReporter(new ShowStuff(system, cable1, cable2, 0.1*0.1));    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue();
    for (int i = 0; i < state.getNU(); ++i)
        state.updU()[i] = 0.1*random.getValue(); 

    system.realize(state, Stage::Position);
    viz.report(state);
    cout << "path1 init length=" << path1.getCableLength(state) << endl;
    cout << "path2 init length=" << path2.getCableLength(state) << endl;
    cout << "Hit ENTER ...";
    getchar();

    path1.setIntegratedCableLengthDot(state, path1.getCableLength(state));
    path2.setIntegratedCableLengthDot(state, path2.getCableLength(state));


    // Simulate it.
    saveStates.clear(); saveStates.reserve(2000);

    RungeKuttaMersonIntegrator integ(system);
    //CPodesIntegrator integ(system);
    integ.setAccuracy(1e-3);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ShowStuff::showHeading(cout);

    const Real finalTime = 2;
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
