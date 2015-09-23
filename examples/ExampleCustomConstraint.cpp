/* -------------------------------------------------------------------------- *
 *                  Simbody(tm) Example: Custom Constraint                    *
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

// This is the Custom Constraint example from the Simbody Advanced Programming
// Guide.
//
// Here we'll make two pendulums hanging near one another, and connect their
// body origins together with the massless rod-like constraint given in
// the above example. This is a simplified version of Simbody's "Rod" (distance)
// constraint, which is not restricted to body origins.
//
// We're going to run and visualize a short simulation and show that energy
// is conserved by the constraint implementation.


#include "Simbody.h"

using namespace SimTK;
using std::cout; using std::endl;

// This constraint holds the body origins of two bodies apart at a given
// distance, without constraining any other motion. This is like connecting
// the points with a massless rod with ball joints at each end. We will
// use the following equation for this holonomic (position-level) constraint:
//       perr(q) = 0
// where perr(q) = (dot(r, r) - d^2)/2
// with r(q) the vector from body1's origin to body2's origin, and d a given
// constant.
// We must supply the first and second time derivatives of perr here as a
// velocity error function verr and acceleration error function aerr:
//       verr(q,u) = d/dt perr = dot(v, r)
//       aerr(q,u,udot) = d/dt verr = dot(a, r) + dot(v, v)
// where v = d/dt r is the relative velocity between the body origins and
// a = d/dt v is their relative acceleration.
//
class ExampleConstraint : public Constraint::Custom::Implementation {
public:
    // Constructor takes two bodies and the desired separation distance
    // between their body frame origins. Tell the base class that this
    // constraint generates 1 holonomic (position level), 0 nonholonomic
    // (velocity level), and 0 acceleration-only constraint equations.
    ExampleConstraint(MobilizedBody& b1, MobilizedBody& b2, Real distance)
    :   Implementation(b1.updMatterSubsystem(), 1, 0, 0), distance(distance) {
        body1 = addConstrainedBody(b1);
        body2 = addConstrainedBody(b2);
    }

    // Implement required pure virtual method.
    Implementation* clone () const override {return new ExampleConstraint(*this);}

    // Implement the Implementation virtuals required for a holonomic
    // (position level) constraint.

    // Simbody supplies position information in argument list; we calculate
    // the constraint error that represents here.
    void calcPositionErrors
       (const State&                                    state,
        const Array_<Transform,ConstrainedBodyIndex>&   X_AB,
        const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
        Array_<Real>&                                   perr) const override
    {
        Vec3 r1 = getBodyOriginLocation(X_AB, body1);
        Vec3 r2 = getBodyOriginLocation(X_AB, body2);
        Vec3 r = r2-r1;
        perr[0] = (dot(r,r)-distance*distance)/2;
    }

    // Simbody supplies velocity information in argument list; position info
    // is in the state. Return time derivative of position constraint error.
    void calcPositionDotErrors
       (const State&                                    state,
        const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB,
        const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
        Array_<Real>&                                   pverr) const override
    {
        Vec3 r1 = getBodyOriginLocationFromState(state, body1);
        Vec3 r2 = getBodyOriginLocationFromState(state, body2);
        Vec3 r = r2-r1;
        Vec3 v1 = getBodyOriginVelocity(V_AB, body1);
        Vec3 v2 = getBodyOriginVelocity(V_AB, body2);
        Vec3 v = v2-v1;
        pverr[0] = dot(v, r);
    }

    // Simbody supplies acceleration information in argument list; position and
    // velocity info is in the state. Return second time derivative of position
    // constraint error.
    void calcPositionDotDotErrors
       (const State&                                    state,
        const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB,
        const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
        Array_<Real>&                                   paerr) const override
    {
        Vec3 r1 = getBodyOriginLocationFromState(state, body1);
        Vec3 r2 = getBodyOriginLocationFromState(state, body2);
        Vec3 r = r2-r1;
        Vec3 v1 = getBodyOriginVelocityFromState(state, body1);
        Vec3 v2 = getBodyOriginVelocityFromState(state, body2);
        Vec3 v = v2-v1;
        Vec3 a1 = getBodyOriginAcceleration(A_AB, body1);
        Vec3 a2 = getBodyOriginAcceleration(A_AB, body2);
        Vec3 a = a2-a1;
        paerr[0] = dot(a, r) + dot(v, v);
    }

    // Simbody provides calculated constraint multiplier in argument list; we
    // turn that into forces here and apply them to the two bodies as point
    // forces at the origins.
    void addInPositionConstraintForces
       (const State&                                state,
        const Array_<Real>&                         multipliers,
        Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
        Array_<Real,      ConstrainedQIndex>&       qForces) const override
    {
        Vec3 r1 = getBodyOriginLocationFromState(state, body1);
        Vec3 r2 = getBodyOriginLocationFromState(state, body2);
        Vec3 r = r2-r1;
        Vec3 force = multipliers[0]*r;
        addInStationForce(state, body2, Vec3(0),  force, bodyForcesInA);
        addInStationForce(state, body1, Vec3(0), -force, bodyForcesInA);
    }

private:
    ConstrainedBodyIndex    body1, body2;
    Real                    distance;
};

// This will be used to report energy periodically. See TextDataEventReporter
// for more information.
class MyEvaluateEnergy : public TextDataEventReporter::UserFunction<Real> {
public:
    Real evaluate(const System& system, const State& state) override {
        const MultibodySystem& mbs = MultibodySystem::downcast(system);
        mbs.realize(state, Stage::Dynamics);
        return mbs.calcEnergy(state);
    }
};

int main() {
  try {
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Hint to Visualizer: don't show ground plane.
    system.setUseUniformBackground(true);

    // Add gravity as a force element.
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.81, 0));

    // Create the body and some artwork for it.
    const Vec3 halfLengths(.5, .1, .25); // half-size of brick (m)
    const Real mass = 2; // total mass of brick (kg)
    Body::Rigid pendulumBody(MassProperties(mass, Vec3(0),
                                mass*UnitInertia::brick(halfLengths)));
    pendulumBody.addDecoration(Transform(),
        DecorativeBrick(halfLengths).setColor(Red));

    // Add an instance of the body to the multibody system by connecting
    // it to Ground via a Ball mobilizer.
    MobilizedBody::Ball pendulum1(matter.updGround(), Transform(Vec3(-1,-1, 0)),
                                  pendulumBody,       Transform(Vec3( 0, 1, 0)));

    // Add a second instance of the pendulum nearby.
    MobilizedBody::Ball pendulum2(matter.updGround(), Transform(Vec3(1,-1, 0)),
                                  pendulumBody,       Transform(Vec3(0, 1, 0)));

    // Connect the origins of the two pendulum bodies together with our
    // rod-like custom constraint.
    const Real d = 1.5; // desired separation distance
    Constraint::Custom rod(new ExampleConstraint(pendulum1, pendulum2, d));

    // Visualize with default options.
    Visualizer viz(system);

    // Add a rubber band line connecting the origins of the two bodies to
    // represent the rod constraint.
    viz.addRubberBandLine(pendulum1, Vec3(0), pendulum2, Vec3(0),
        DecorativeLine().setColor(Blue).setLineThickness(3));

    // Ask for a report every 1/30 of a second to match the Visualizer's
    // default rate of 30 frames per second.
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    // Output total energy to the console once per second.
    system.addEventReporter(new TextDataEventReporter
                                   (system, new MyEvaluateEnergy(), 1.0));

    // Initialize the system and state.
    State state = system.realizeTopology();

    // Orient the two pendulums asymmetrically so they'll do something more
    // interesting than just hang there.
    pendulum1.setQToFitRotation(state, Rotation(Pi/4, ZAxis));
    pendulum2.setQToFitRotation(state, Rotation(BodyRotationSequence,
                                                Pi/4, ZAxis, Pi/4, YAxis));

    // Evaluate the system at the new state and draw one frame manually.
    system.realize(state);
    viz.report(state);

    // Simulate it.
    cout << "Hit ENTER to run a short simulation.\n";
    cout << "(Energy should be conserved to about four decimal places.)\n";
    getchar();

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-4); // ask for about 4 decimal places (default is 3)
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);

  } catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }

    return 0;
}
