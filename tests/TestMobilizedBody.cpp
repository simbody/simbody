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
    assertEqual(b1.locateBodyPointOnGround(state, Vec3(0)), b1.getBodyOriginLocation(state));
    assertEqual(b1.locateGroundPointOnBody(state, b1.locateBodyPointOnGround(state, point)), point);
    assertEqual(b2.locateGroundPointOnBody(state, b1.locateBodyPointOnGround(state, point)), b1.locateBodyPointOnBody(state, point, b2));
    assertEqual(b2.locateGroundPointOnBody(state, b1.locateBodyPointOnGround(state, Vec3(0))).norm(), (b1.getBodyOriginLocation(state)-b2.getBodyOriginLocation(state)).norm());
    assertEqual(b2.locateBodyMassCenterOnGround(state), b2.locateBodyPointOnGround(state, b2.getBodyMassCenterStation(state)));
    assertEqual(b1.expressBodyVectorInGround(state, Vec3(0)), Vec3(0));
    assertEqual(b1.expressBodyVectorInGround(state, point), b1.getBodyRotation(state)*point);
    assertEqual(b1.expressGroundVectorInBody(state, b1.expressBodyVectorInGround(state, point)), point);
    assertEqual(b2.expressGroundVectorInBody(state, b1.expressBodyVectorInGround(state, point)), b1.expressBodyVectorInBody(state, point, b2));
    
    // Test the routines for mapping locations, velocities, and accelerations.
    
    Vec3 r, v, a;
    b1.calcBodyFixedPointLocationVelocityAndAccelerationInGround(state, point, r, v, a);
    assertEqual(v, b1.calcBodyFixedPointVelocityInGround(state, point));
    assertEqual(a, b1.calcBodyFixedPointAccelerationInGround(state, point));
    {
        Vec3 r2, v2;
        b1.calcBodyFixedPointLocationAndVelocityInGround(state, point, r2, v2);
        assertEqual(r, r2);
        assertEqual(v, v2);
    }
    assertEqual(b1.calcBodyFixedPointVelocityInGround(state, Vec3(0)), b1.getBodyOriginVelocity(state));
    assertEqual(b1.calcBodyFixedPointAccelerationInGround(state, Vec3(0)), b1.getBodyOriginAcceleration(state));
    assertEqual(b1.calcBodyFixedPointVelocityInGround(state, point), b1.calcStationVelocityInBody(state, point, matter.Ground()));
}

int main() {
    try {
        testCalculationMethods();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

