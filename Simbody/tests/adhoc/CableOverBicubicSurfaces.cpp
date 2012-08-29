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

/*                     Simbody CableOverBicubicSurfaces
This example shows how to use a CableTrackerSubsystem to follow the motion of
a cable that crosses bicubic surfaces. We'll then
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
              const MyCableSpring& cable1, Real interval) 
    :   PeriodicEventReporter(interval), 
        mbs(mbs), cable1(cable1) {}

    static void showHeading(std::ostream& o) {
        printf("%8s %10s %10s %10s %10s %10s %10s %10s %10s %12s\n",
            "time", "length", "rate", "integ-rate", "unitpow", "tension", "disswork",
            "KE", "PE", "KE+PE-W");
    }

    /** This is the implementation of the EventReporter virtual. **/ 
    void handleEvent(const State& state) const OVERRIDE_11 {
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
    MyCableSpring           cable1;
};

int main() {
  try {    
    // Create the system.   
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableTrackerSubsystem cables(system);
    GeneralForceSubsystem forces(system);

    Force::Gravity gravity(forces, matter, -YAxis, 9.81);
    //Force::GlobalDamper(forces, matter, 5);

    system.setUseUniformBackground(true); // no ground plane in display
    MobilizedBody Ground = matter.Ground(); // convenient abbreviation

    // Read in some bones.
    PolygonalMesh femur, tibia;
    femur.loadVtpFile("CableOverBicubicSurfaces-femur.vtp");
    tibia.loadVtpFile("CableOverBicubicSurfaces-tibia.vtp");
    femur.scaleMesh(20);
    tibia.scaleMesh(20);

    // Create some bicubic surfaces.
    const int Nx = 4, Ny = 5;
    const Real xData[Nx] = { .1, 1, 2, 4 };
    const Real yData[Ny] = { -3, -2, 0, 1, 3 };
    const Real fData[Nx*Ny] = { 1,   2,   3,   3,   2,
                                1.1, 2.1, 3.1, 3.1, 2.1,
                                1,   2,   7,   3,   2,
                                1.2, 2.2, 3.2, 3.2, 2.2 };
    const Vector x(Nx,   xData);
    const Vector y(Ny,   yData);
    const Matrix f(Nx,Ny, fData);
    BicubicSurface rough(x, y, f, 0);       // raw
    BicubicSurface smooth(x, y, f, 1);      // smoothed

    Vector xp(Vec2(.25,3.25)), yp(Vec2(.75,5.75));
    // Nice patch:
    //Matrix fp(Mat22(1, 1,
    //                1, 1));
    //Matrix fxp(.5*Mat22(-1,  0,
    //                    0, -1));
    //Matrix fyp(.5*Mat22(-1,  0,
    //                    0,  -1));
    //Matrix fxyp(0*Mat22(1, -1,
    //                    -1, 1));
    // One-hump patch:
    Matrix fp(Mat22(1, 1,
                    1, 1));
    Matrix fxp(Mat22(1,  1,
                     -1, -1));
    Matrix fyp(Mat22(1,  -1,
                     1,  -1));
    Matrix fxyp(0.5*Mat22(1, 3,
                        -3, 4));
    BicubicSurface patch(xp, yp, fp, fxp, fyp, fxyp);
    Rotation xm90(-Pi/2, XAxis);
    Transform patchPose(xm90, Vec3(4,2,0));
    
    // Ask the bicubic surfaces for some meshes we can use for display.
    Real resolution = 31;
    PolygonalMesh patchMesh = patch.createPolygonalMesh(resolution);
    PolygonalMesh roughMesh = rough.createPolygonalMesh(resolution);
    PolygonalMesh smoothMesh = smooth.createPolygonalMesh(resolution);

    const Vec3 SmoothOrigin(-3,-3,-3);
    // Shift the drawing slightly in the -z direction so that the path
    // shows better.
    Ground.addBodyDecoration(SmoothOrigin - Vec3(0,0,.01),
        DecorativeMesh(smoothMesh).setColor(Cyan).setOpacity(.75));

    // Not using these yet:
    Ground.addBodyDecoration(Vec3(5,-5,0),
        DecorativeMesh(patchMesh).setColor(Gray));
    Ground.addBodyDecoration(Vec3(5,0,0),
        DecorativeMesh(roughMesh).setColor(Red));
    Ground.addBodyDecoration(Vec3(5,0,0),
        DecorativeMesh(femur).setColor(Vec3(.8,.8,.8)));
    Ground.addBodyDecoration(Vec3(5,-4,0),
        DecorativeMesh(tibia).setColor(Vec3(.8,.8,.8)));


    Body::Rigid someBody(MassProperties(2.0, Vec3(0,-4,0), 
        UnitInertia::cylinderAlongY(1,4).shiftFromCentroid(Vec3(0,4,0))));

    someBody.addDecoration(Transform(Rotation(Pi,ZAxis),Vec3(0,-4,0)), 
        DecorativeCylinder(1,4).setColor(Yellow)
                            .setOpacity(.5).setResolution(4));
    someBody.addDecoration(Transform(), 
        DecorativeMesh(femur).setColor(Vec3(.8,.8,.8)));

    MobilizedBody::Free body1(Ground,    Transform(Vec3(0)), 
                              someBody,  Transform(Vec3(0,0,0)));

    CablePath path1(cables, Ground, Vec3(.5,-.5,0),   // origin
                            body1, Vec3(0,0,0));  // termination

    CableObstacle::Surface obstacle1(path1, Ground, SmoothOrigin, 
                                     ContactGeometry::SmoothHeightMap(smooth));
    // Provide an initial guess for P and Q (in frame of "smooth").
    Vec3 P1(1.5,1,3.75), Q1(1.5,-1,3.75);
    obstacle1.setContactPointHints(P1, Q1);
    Ground.addBodyDecoration(SmoothOrigin,
        DecorativePoint(P1).setColor(Green).setScale(2));
    Ground.addBodyDecoration(SmoothOrigin,
        DecorativePoint(Q1).setColor(Red).setScale(2));

    MyCableSpring cable1(forces, path1, 50., 8., 0.1); 

    //Force::TwoPointLinearSpring spring1(forces, body1, Vec3(0,-8,0),
    //    Ground, Vec3(0,0,0), 30., 12.);

    Visualizer viz(system);
    viz.setShowFrameNumber(true);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    system.addEventReporter(new ShowStuff(system, cable1, 0.02));    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    //Random::Gaussian random;
    //for (int i = 0; i < state.getNQ(); ++i)
    //    state.updQ()[i] = random.getValue();
    //for (int i = 0; i < state.getNU(); ++i)
    //    state.updU()[i] = 0.1*random.getValue(); 
    body1.setQToFitTranslation(state, Vec3(4,-10,-3));
    body1.setQToFitRotation(state, Rotation(-Pi, ZAxis));

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
