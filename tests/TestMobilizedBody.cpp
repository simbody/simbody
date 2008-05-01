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

int testCalcQDotDot() {
    
    // Build a system consisting of a chain of bodies with every possible mobilizer.
    
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)))
                                  .addDecoration(Transform(), DecorativeSphere(.1));
    Random::Uniform random(0.0, 2.0);
    MobilizedBody lastBody = MobilizedBody::Pin(matter.Ground(), Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Slider(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Universal(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Cylinder(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::BendStretch(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Planar(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Gimbal(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Ball(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Translation(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Free(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::LineOrientation(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::FreeLine(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
//    lastBody = MobilizedBody::Weld(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue()))); // TODO Not yet implemented
    lastBody = MobilizedBody::Screw(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())), 0.5);
    lastBody = MobilizedBody::Ellipsoid(lastBody, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    mbs.realizeTopology();
    State& s = mbs.updDefaultState();
    mbs.realizeModel(s);
    
    // Choose a random initial conformation.
    
    for (int i = 0; i < s.getNQ(); ++i)
        s.updQ()[i] = random.getValue();
    mbs.realize(s, Stage::Instance);
    mbs.project(s, 0.01, Vector(s.getNY()), Vector(s.getNYErr()), Vector(s.getNY())); // Normalize the quaternions
    mbs.realize(s, Stage::Position);
    
    // Choose a random udot vector.  Then calculate qdotdot twice: once by calling calcQDotDot() and once by
    // finite differences.
    
    Vector udot(s.getNU());
    for (int i = 0; i < udot.size(); ++i)
        udot[i] = random.getValue();
    Vector qdot1(s.getNQ()), qdot2(s.getNQ());
    Real delta = 1;
    matter.calcQDot(s, s.getU(), qdot1);
    matter.calcQDot(s, s.getU()+delta*udot, qdot2);
    Vector qdotdot1 = (qdot2-qdot1)/delta;
    Vector qdotdot2(s.getNQ());
    matter.calcQDotDot(s, udot, qdotdot2);
std::cout << qdotdot1 << std::endl;
std::cout << qdotdot2 << std::endl;
    std::cout << max(abs(qdotdot2-qdotdot1)) << std::endl;
    ASSERT(max(abs(qdotdot2-qdotdot1)) < sqrt(delta));
    
    // Convert to Euler angles and make sure the positions are all the same.
    
    State euler = s;
    matter.convertToEulerAngles(s, euler);
    mbs.realize(euler, Stage::Position);
    for (int i = 0; i < matter.getNBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (body.getBodyOriginLocation(euler)-body.getBodyOriginLocation(s)).norm();
        ASSERT(dist < 1e-5);
    }

    // Now convert back to quaternions and make sure the positions are still the same.
    
    State quaternions = s;
    matter.convertToQuaternions(euler, quaternions);
    mbs.realize(quaternions, Stage::Position);
    for (int i = 0; i < matter.getNBodies(); ++i) {
        const MobilizedBody& body = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (body.getBodyOriginLocation(quaternions)-body.getBodyOriginLocation(s)).norm();
        ASSERT(dist < 1e-5);
    }
    
    // Compare the state variables to see if they have been accurately reproduced.
    
    mbs.project(s, 0.01, Vector(s.getNY()), Vector(s.getNYErr()), Vector(s.getNY())); // Normalize the quaternions
    Real diff = std::sqrt((s.getQ()-quaternions.getQ()).normSqr()/s.getNQ());
    ASSERT(diff < 1e-5);
}

int main() {
    try {
        testCalculationMethods();
        testCalcQDotDot();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

