/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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
#include <string>

using std::cout;
using std::endl;
using std::string;

using namespace SimTK;

/**
 * This test measures the speed of various multibody calculations.  It executes a collection of operations
 * on various systems.  The CPU time required to perform each operation 1000 times is measured and
 * printed to the console.
 *
 * Each test system contains 256 identical bodies (plus ground), but they differ in the type of
 * bodies and their arrangement into a multibody tree.  The arrangements include 1) all bodies attached
 * directly to ground, 2) the bodies linked in a single chain, and 3) the bodies arranged to form
 * a binary tree.
 */

// The following routines define the operations to be profiled.

void doRealizePosition(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Time);
    system.realize(state, Stage::Position);
}

void doRealizeVelocity(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Time);
    system.realize(state, Stage::Velocity);
}

void doRealizeAcceleration(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Time);
    system.realize(state, Stage::Acceleration);
}

void doCalcMV(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector v(matter.getNumMobilities(), 1.0);
    Vector mv;
    matter.calcMV(state, v, mv);
}

void doCalcMInverseV(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector v(matter.getNumMobilities(), 1.0);
    Vector minvv;
    matter.calcMInverseV(state, v, minvv);
}

void doCalcResidualForceIgnoringConstraints(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector appliedMobilityForces(matter.getNumMobilities(), 1.0);
    Vector_<SpatialVec> appliedBodyForces(matter.getNumBodies(), SpatialVec(Vec3(1, 0, 0), Vec3(0, 1, 0)));
    Vector knownUdot, residualMobilityForces;
    matter.calcResidualForceIgnoringConstraints(state, appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);
}

void doCalcMobilizerReactionForces(MultibodySystem& system, State& state) {
    Vector_<SpatialVec> forces;
    system.getMatterSubsystem().calcMobilizerReactionForces(state, forces);
}

void doCalcInternalGradientFromSpatial(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec> dEdR(matter.getNumBodies(), SpatialVec(Vec3(1, 0, 0), Vec3(0, 1, 0)));
    Vector dEdQ;
    matter.calcInternalGradientFromSpatial(state, dEdR, dEdQ);
}

void doCalcCompositeBodyInertias(MultibodySystem& system, State& state) {
    Vector_<SpatialMat> r;
    system.getMatterSubsystem().calcCompositeBodyInertias(state, r);
}

/**
 * Time how long it takes to perform an operation 1000 times.  The test is repeated 5 times,
 * and the average is returned.  The return value represents CPU time, *not* clock time.
 */
void timeComputation(MultibodySystem& system, void function(MultibodySystem& system, State& state), const string& name, int iterations) {
    const int repeats = 5;
    Vector cpuTimes(repeats);
    Vector clockTimes(repeats);
    State state = system.getDefaultState();
    system.realize(state, Stage::Acceleration);

    // Repeatedly measure the CPU time for performing the operation 1000 times.

    for (int i = 0; i < repeats; i++) {
        Real startCpu = cpuTime();
        Real startClock = realTime();
        for (int j = 0; j < iterations; j++)
            function(system, state);
        Real endCpu = cpuTime();
        Real endClock = realTime();
        cpuTimes[i] = endCpu-startCpu;
        clockTimes[i] = endClock-startClock;
    }
    std::cout << name<<": "<<(mean(cpuTimes)*1000/iterations)<<" "<<(mean(clockTimes)*1000/iterations)<< std::endl;
}

/**
 * Time all the different calculations for one system.
 */
void runAllTests(MultibodySystem& system) {
    timeComputation(system, doRealizePosition, "realizePosition", 1000);
    timeComputation(system, doRealizeVelocity, "realizeVelocity", 1000);
    timeComputation(system, doRealizeAcceleration, "realizeAcceleration", 1000);
    timeComputation(system, doCalcMV, "calcMV", 5000);
    timeComputation(system, doCalcMInverseV, "calcMInverseV", 5000);
    timeComputation(system, doCalcResidualForceIgnoringConstraints, "calcResidualForceIgnoringConstraints", 5000);
    timeComputation(system, doCalcMobilizerReactionForces, "calcMobilizerReactionForces", 1000);
    timeComputation(system, doCalcInternalGradientFromSpatial, "calcInternalGradientFromSpatial", 5000);
    timeComputation(system, doCalcCompositeBodyInertias, "calcCompositeBodyInertias", 5000);
}

// The following routines create the systems to be profiled.

void createParticles(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    for (int i = 0; i < 256; i++)
        MobilizedBody::Translation next(matter.updGround(), Vec3(1, 0, 0), body, Vec3(0));
    system.realizeTopology();
}

void createFreeBodies(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    for (int i = 0; i < 256; i++)
        MobilizedBody::Free next(matter.updGround(), Vec3(1, 0, 0), body, Vec3(0));
    system.realizeTopology();
}

void createPinChain(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    MobilizedBody last = matter.updGround();
    for (int i = 0; i < 256; i++) {
        MobilizedBody::Pin next(last, Vec3(1, 0, 0), body, Vec3(0));
        last = next;
    }
    system.realizeTopology();
}

void createBallChain(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    MobilizedBody last = matter.updGround();
    for (int i = 0; i < 256; i++) {
        MobilizedBody::Ball next(last, Vec3(1, 0, 0), body, Vec3(0));
        last = next;
    }
    system.realizeTopology();
}

void createPinTree(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    for (int i = 0; i < 256; i++) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(i/2));
        MobilizedBody::Pin next(parent, Vec3(1, 0, 0), body, Vec3(0));
    }
    system.realizeTopology();
}

void createBallTree(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    for (int i = 0; i < 256; i++) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(i/2));
        MobilizedBody::Ball next(parent, Vec3(1, 0, 0), body, Vec3(0));
    }
    system.realizeTopology();
}

int main() {
    {
        std::cout << "\nParticles:\n" << std::endl;
        MultibodySystem system;
        createParticles(system);
        runAllTests(system);
    }
    {
        std::cout << "\nFree Bodies:\n" << std::endl;
        MultibodySystem system;
        createFreeBodies(system);
        runAllTests(system);
    }
    {
        std::cout << "\nPin Chain:\n" << std::endl;
        MultibodySystem system;
        createPinChain(system);
        runAllTests(system);
    }
    {
        std::cout << "\nBall Chain:\n" << std::endl;
        MultibodySystem system;
        createBallChain(system);
        runAllTests(system);
    }
    {
        std::cout << "\nPin Tree:\n" << std::endl;
        MultibodySystem system;
        createPinTree(system);
        runAllTests(system);
    }
    {
        std::cout << "\nBall Tree:\n" << std::endl;
        MultibodySystem system;
        createBallTree(system);
        runAllTests(system);
    }
}