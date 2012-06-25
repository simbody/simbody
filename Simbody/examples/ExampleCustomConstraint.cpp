/* -------------------------------------------------------------------------- *
 *                  Simbody(tm): Custom Constraint Example                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
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
    Implementation* clone () const {return new ExampleConstraint(*this);}

    // Implement the Implementation virtuals required for a holonomic
    // (position level) constraint.

    void calcPositionErrors     
       (const State&                                    state,
        const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
        const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
        Array_<Real>&                                   perr) const
    {
        Vec3 r1 = getBodyOriginLocation(X_AB, body1);
        Vec3 r2 = getBodyOriginLocation(X_AB, body2);
        Vec3 r = r2-r1;
        perr[0] = (dot(r,r)-distance*distance)/2;
    }

    void calcPositionDotErrors   
       (const State&                                    state,
        const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
        const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
        Array_<Real>&                                   pverr) const
    {
        Vec3 r1 = getBodyOriginLocationFromState(state, body1);
        Vec3 r2 = getBodyOriginLocationFromState(state, body2);
        Vec3 r = r2-r1;
        Vec3 v1 = getBodyOriginVelocity(V_AB, body1);
        Vec3 v2 = getBodyOriginVelocity(V_AB, body2);
        Vec3 v = v2-v1;
        pverr[0] = dot(v, r);
    }

    void calcPositionDotDotErrors
       (const State&                                    state,
        const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
        const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
        Array_<Real>&                                   paerr) const 
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

    void addInPositionConstraintForces
       (const State&                                state, 
        const Array_<Real>&                         multipliers,
        Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
        Array_<Real,      ConstrainedQIndex>&       qForces) const 
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

class MyTextDataUserFunction 
:   public TextDataEventReporter::UserFunction<Real> {
public:
    Real evaluate(const System& system, const State& state) {
        const MultibodySystem& mbs = MultibodySystem::downcast(system);
        return mbs.calcEnergy(state);
    }
};

int main() {
  try {   
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    system.setUseUniformBackground(true); // don't show ground in visualizer
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

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
    // it to Ground via a pin mobilizer.
    MobilizedBody::Ball pendulum1(matter.updGround(), 
                                Transform(Vec3(-1,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));

    // Add an instance of the body to the multibody system by connecting
    // it to Ground via a pin mobilizer.
    MobilizedBody::Ball pendulum2(matter.updGround(), 
                                Transform(Vec3(1,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));

    const Real d = 1.5; // desired separation distance
    Constraint::Custom rod(new ExampleConstraint(pendulum1, pendulum2, d));

    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(system);

    // Add a rubber band line connecting the origins of the two bodies to
    // represent the rod constraint.
    viz.addRubberBandLine(pendulum1, Vec3(0), pendulum2, Vec3(0),
        DecorativeLine().setColor(Blue).setLineThickness(3));

    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    system.addEventReporter(new TextDataEventReporter(system,
                                                new MyTextDataUserFunction(),
                                                .1));
    
    // Initialize the system and state.
    
    State state = system.realizeTopology();

    pendulum1.setQToFitRotation(state, Rotation(Pi/4, ZAxis));
    pendulum2.setQToFitRotation(state, Rotation(BodyRotationSequence,
                                                Pi/4, ZAxis, Pi/4, YAxis));

    system.realize(state);

    viz.report(state);
    // Simulate it.
    cout << "Hit ENTER to run a short simulation ...";
    getchar();

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-4);
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
