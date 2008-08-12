/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

#include "SimTKsimbody.h"
#include "../src/MobilizedBodyImpl.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2) {
    ASSERT(abs(val1-val2) < TOL);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < TOL);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2) {
    assertEqual(val1[0], val2[0]);
    assertEqual(val1[1], val2[1]);
}

template<>
void assertEqual(Transform val1, Transform val2) {
    assertEqual(val1.T(), val2.T());
    assertEqual(val1.R().convertRotationToBodyFixedXYZ(), val2.R().convertRotationToBodyFixedXYZ());
}

void testCalculationMethods() {
    
    // Create a system with two bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free b1(matter.Ground(), body);
    MobilizedBody::Free b2(matter.Ground(), body);
    
    // Set all the state variables to random values.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNY(); ++i)
        state.updY()[i] = random.getValue();
    system.realize(state, Stage::Acceleration);
    
    // Test the low level methods for transforming points and vectors.
    
    const Vec3 point(0.5, 1, -1.5);
    assertEqual(b1.findStationLocationInGround(state, Vec3(0)), b1.getBodyOriginLocation(state));
    assertEqual(b1.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)), point);
    assertEqual(b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, point)), b1.findStationLocationInAnotherBody(state, point, b2));
    assertEqual(b2.findStationAtGroundPoint(state, b1.findStationLocationInGround(state, Vec3(0))).norm(), (b1.getBodyOriginLocation(state)-b2.getBodyOriginLocation(state)).norm());
    assertEqual(b2.findMassCenterLocationInGround(state), b2.findStationLocationInGround(state, b2.getBodyMassCenterStation(state)));
    assertEqual(b1.expressVectorInGroundFrame(state, Vec3(0)), Vec3(0));
    assertEqual(b1.expressVectorInGroundFrame(state, point), b1.getBodyRotation(state)*point);
    assertEqual(b1.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)), point);
    assertEqual(b2.expressGroundVectorInBodyFrame(state, b1.expressVectorInGroundFrame(state, point)), b1.expressVectorInAnotherBodyFrame(state, point, b2));
    
    // Test the routines for mapping locations, velocities, and accelerations.
    
    Vec3 r, v, a;
    b1.findStationLocationVelocityAndAccelerationInGround(state, point, r, v, a);
    assertEqual(v, b1.findStationVelocityInGround(state, point));
    assertEqual(a, b1.findStationAccelerationInGround(state, point));
    {
        Vec3 r2, v2;
        b1.findStationLocationAndVelocityInGround(state, point, r2, v2);
        assertEqual(r, r2);
        assertEqual(v, v2);
    }
    assertEqual(b1.findStationVelocityInGround(state, Vec3(0)), b1.getBodyOriginVelocity(state));
    assertEqual(b1.findStationAccelerationInGround(state, Vec3(0)), b1.getBodyOriginAcceleration(state));
    assertEqual(b1.findStationVelocityInGround(state, point), b1.findStationVelocityInAnotherBody(state, point, matter.Ground()));
}

void testWeld() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -1, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    
    // Create two pendulums, each with two welded bodies.  One uses a Weld MobilizedBody,
    // and the other uses a Weld constraint.
    
    Transform inboard(Vec3(0.1, 0.5, -1));
    Transform outboard(Vec3(0.2, -0.2, 0));
    MobilizedBody::Ball p1(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Ball p2(matter.updGround(), Vec3(0), body, Vec3(0, 1, 0));
    MobilizedBody::Weld c1(p1, inboard, body, outboard);
    MobilizedBody::Free c2(p2, inboard, body, outboard);
    Constraint::Weld constraint(p2, inboard, c2, outboard);
    State state = system.realizeTopology();
    p1.setU(state, Vec3(1, 2, 3));
    p2.setU(state, Vec3(1, 2, 3));
    system.realize(state, Stage::Velocity);
    system.project(state, 1e-10, Vector(state.getNY()), Vector(state.getNYErr()), Vector());
    system.realize(state, Stage::Velocity);
    assertEqual(c1.getBodyTransform(state), c2.getBodyTransform(state));
    assertEqual(c1.getBodyVelocity(state), c2.getBodyVelocity(state));
    
    // Simulate it and see if both pendulums behave identically.
    
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5.0);
    system.realize(integ.getState(), Stage::Velocity);
    assertEqual(c1.getBodyTransform(integ.getState()), c2.getBodyTransform(integ.getState()));
    assertEqual(c1.getBodyVelocity(integ.getState()), c2.getBodyVelocity(integ.getState()));
}

int main() {
    try {
        testCalculationMethods();
        testWeld();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

