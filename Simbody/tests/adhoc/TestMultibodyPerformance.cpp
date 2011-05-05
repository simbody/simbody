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
#include <ctime>

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
    Array_<SpatialInertia> r;
    system.getMatterSubsystem().calcCompositeBodyInertias(state, r);
}

static Real flopTimeInNs;

/**
 * Time how long it takes to perform an operation 1000 times.  The test is repeated 5 times,
 * and the average is returned.  The return value represents CPU time, *not* clock time.
 */
void timeComputation(MultibodySystem& system, void function(MultibodySystem& system, State& state), 
                     const string& name, int iterations, bool useEulerAngles) {
    const int repeats = 5;
    Vector cpuTimes(repeats);
    State state = system.getDefaultState();
    if (useEulerAngles) {
        system.getMatterSubsystem().setUseEulerAngles(state, true);
        system.realizeModel(state);
    }
    system.realize(state, Stage::Acceleration);

    const int ndof = system.getMatterSubsystem().getNumMobilities();
    const int nmovbod = system.getMatterSubsystem().getNumBodies()-1; // not Ground

    // Repeatedly measure the CPU time for performing the operation 1000 times.

    for (int i = 0; i < repeats; i++) {
        double startCpu = threadCpuTime();
        for (int j = 0; j < iterations; j++)
            function(system, state);
        double endCpu = threadCpuTime();
        cpuTimes[i] = endCpu-startCpu;
    }

    Real timePerIterUs = mean(cpuTimes)*1000000/iterations; // us
    Real flopTimeUs = flopTimeInNs / 1000;
    Real flopTimePerIter = timePerIterUs/flopTimeUs;
    std::printf("%40s:%6.4gus -> %4d flp/dof, %4d flp/bod\n",
        name.c_str(), timePerIterUs, (int)(flopTimePerIter/ndof),
        (int)(flopTimePerIter/nmovbod));
}

/**
 * Time all the different calculations for one system.
 */
void runAllTests(MultibodySystem& system, bool useEulerAngles=false) {
    std::cout << "# dofs=" << system.getMatterSubsystem().getNumMobilities() << "\n";
    timeComputation(system, doRealizePosition, "realizePosition", 5000, useEulerAngles);
    timeComputation(system, doRealizeVelocity, "realizeVelocity", 5000, useEulerAngles);
    timeComputation(system, doRealizeAcceleration, "realizeAcceleration", 2000, useEulerAngles);
    timeComputation(system, doCalcMV, "calcMV", 5000, useEulerAngles);
    timeComputation(system, doCalcMInverseV, "calcMInverseV", 5000, useEulerAngles);
    timeComputation(system, doCalcResidualForceIgnoringConstraints, "calcResidualForceIgnoringConstraints", 5000, useEulerAngles);
    timeComputation(system, doCalcMobilizerReactionForces, "calcMobilizerReactionForces", 1000, useEulerAngles);
    timeComputation(system, doCalcInternalGradientFromSpatial, "calcInternalGradientFromSpatial", 5000, useEulerAngles);
    timeComputation(system, doCalcCompositeBodyInertias, "calcCompositeBodyInertias", 5000, useEulerAngles);
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

void createSliderChain(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    MobilizedBody last = matter.updGround();
    for (int i = 0; i < 256; i++) {
        MobilizedBody::Slider next(last, Vec3(1, 0, 0), body, Vec3(0));
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


void createGimbalChain(MultibodySystem& system) {
    SimbodyMatterSubsystem matter(system);
    Body::Rigid body;
    MobilizedBody last = matter.updGround();
    for (int i = 0; i < 256; i++) {
        MobilizedBody::Gimbal next(last, Vec3(1, 0, 0), body, Vec3(0));
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

static int tenInts[10];
static Real tenReals[10];
// These should multiply out to about 1.
static Real tenMults[10] = 
    {0.501,0.2501,0.201,0.101,1.000000001,
    (1/1.000000002),(1/.101), (1/.201), (1/.2501), (1/.501)};
void testFunctions(Real& flopTime, bool flopTimeOnly=false) {
    Real addRes=1,subRes=1,mulRes=1,divRes=1,sqrtRes=1,oosqrtRes=1,
         sinRes=1,cosRes=1,atan2Res=1,logRes=1,expRes=1;
    int intAddRes=1;
    Random::Uniform rand; rand.setMin(-5); rand.setMax(5);
    for (int i=0; i<10; i++) tenInts[i] = rand.getIntValue();
    for (int i=0; i<10; i++) tenReals[i] = rand.getValue();

    Real tprev = threadCpuTime();
    for (int i = 0; i < 10*100000000; i++) {
        intAddRes += tenInts[0];
        intAddRes -= tenInts[1];
        intAddRes += tenInts[2];
        intAddRes -= tenInts[3];
        intAddRes += tenInts[4];
        intAddRes -= tenInts[5];
        intAddRes += tenInts[6];
        intAddRes -= tenInts[7];
        intAddRes += tenInts[8];
        intAddRes -= tenInts[9];
    }
    Real t = threadCpuTime(); Real intAddTime = (t-tprev)/10; // time for 1e9 ops
    printf("intAdd %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 5*100000000; i++) {
        addRes += 1.1;
        addRes += 1.2;
        addRes += 1.3;
        addRes += 1.4;
        addRes += 1.501;
        addRes += 1.6;
        addRes += 1.7;
        addRes += 1.8;
        addRes += 1.9;
        addRes += 2.007;
    }
    t = threadCpuTime(); Real addTime = (t-tprev)/5;
    printf("add %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 5*100000000; i++) {
        subRes -= 1.1;
        subRes -= 1.2;
        subRes -= 1.3;
        subRes -= 1.4;
        subRes -= 1.501;
        subRes -= 1.6;
        subRes -= 1.7;
        subRes -= 1.8;
        subRes -= 1.9;
        subRes -= 2.007;
    }
    t = threadCpuTime(); Real subTime = (t-tprev)/5;
    printf("sub %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 3*100000000; i++) {
        mulRes *= 0.501;
        mulRes *= 0.2501;
        mulRes *= 0.201;
        mulRes *= 0.101;
        mulRes *= 1.000000001;
        mulRes *= (1/1.000000002); // done at compile time
        mulRes *= (1/.101);
        mulRes *= (1/.201);
        mulRes *= (1/.2501);
        mulRes *= (1/.501);
    }
    t = threadCpuTime(); Real mulTime=(t-tprev)/3;
    printf("mul %gs\n", t-tprev);
    flopTime = (addTime+mulTime)/2;
    std::cout << "1 flop=avg(add,mul)=" << flopTime << "ns\n";
    if (flopTimeOnly)
        return;

    tprev = threadCpuTime();
    for (int i = 0; i < 100000000; i++) {
        divRes /= tenMults[7];
        divRes /= tenMults[9];
        divRes /= tenMults[8];
        divRes /= tenMults[6];
        divRes /= tenMults[2];
        divRes /= tenMults[3];
        divRes /= tenMults[1];
        divRes /= tenMults[0];
        divRes /= tenMults[4];
        divRes /= tenMults[5];
        // prevent clever optimization VC10 did to turn divides
        // into multiplies.
        tenMults[i%10]     *= 1.0000000000001; 
        tenMults[(i+5)%10] *= 0.9999999999999; 
    }
    t = threadCpuTime(); Real divTime=(t-tprev);
    printf("div %gs\n", t-tprev);
    tprev = threadCpuTime();

    for (int i = 0; i < 100000000/2; i++) {
        const Real ir = (Real)i;
        sqrtRes += std::sqrt(ir+0.001); // two adds
        sqrtRes += std::sqrt(ir+0.1);
        sqrtRes += std::sqrt(ir+0.2);
        sqrtRes += std::sqrt(ir+0.3);
        sqrtRes += std::sqrt(ir+0.4);
        sqrtRes += std::sqrt(ir+0.501);
        sqrtRes += std::sqrt(ir+0.6);
        sqrtRes += std::sqrt(ir+0.7);
        sqrtRes += std::sqrt(ir+0.8);
        sqrtRes += std::sqrt(ir+0.9);
    }
    t = threadCpuTime(); Real sqrtTime=2*(t-tprev)-2*addTime;
    printf("sqrt %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/4; i++) {
        const Real ir = (Real)i;
        oosqrtRes += 1/std::sqrt(ir+0.001); // two adds
        oosqrtRes += 1/std::sqrt(ir+0.1);
        oosqrtRes += 1/std::sqrt(ir+0.2);
        oosqrtRes += 1/std::sqrt(ir+0.3);
        oosqrtRes += 1/std::sqrt(ir+0.4);
        oosqrtRes += 1/std::sqrt(ir+0.501);
        oosqrtRes += 1/std::sqrt(ir+0.6);
        oosqrtRes += 1/std::sqrt(ir+0.7);
        oosqrtRes += 1/std::sqrt(ir+0.8);
        oosqrtRes += 1/std::sqrt(ir+0.9);
    }
    t = threadCpuTime(); Real oosqrtTime=4*(t-tprev)-2*addTime;
    printf("1/sqrt %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/5; i++) {
        const Real ir = (Real)i;
        logRes += std::log(ir+0.001); // two adds
        logRes += std::log(ir+0.1);
        logRes += std::log(ir+0.2);
        logRes += std::log(ir+0.3);
        logRes += std::log(ir+0.4);
        logRes += std::log(ir+0.501);
        logRes += std::log(ir+0.6);
        logRes += std::log(ir+0.7);
        logRes += std::log(ir+0.8);
        logRes += std::log(ir+0.9);
    }
    t = threadCpuTime(); Real logTime=(t-tprev)*5-2*addTime;
    printf("log %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/5; i++) {
        const Real ir = .000000001*(Real)i;
        expRes += std::exp(ir+0.001); // two adds
        expRes += std::exp(ir+0.1);
        expRes += std::exp(ir+0.2);
        expRes += std::exp(ir+0.3);
        expRes += std::exp(ir+0.4);
        expRes += std::exp(ir+0.501);
        expRes += std::exp(ir+0.6);
        expRes += std::exp(ir+0.7);
        expRes += std::exp(ir+0.8);
        expRes += std::exp(ir+0.9);
    }
    t = threadCpuTime(); Real expTime=(t-tprev)*5-2*addTime;
    printf("exp %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        sinRes += std::sin(ir+0.001); // two adds
        sinRes += std::sin(ir+0.1);
        sinRes += std::sin(ir+0.2);
        sinRes += std::sin(ir+0.3);
        sinRes += std::sin(ir+0.4);
        sinRes += std::sin(ir+0.501);
        sinRes += std::sin(ir+0.6);
        sinRes += std::sin(ir+0.7);
        sinRes += std::sin(ir+0.8);
        sinRes += std::sin(ir+0.9);
    }
    t = threadCpuTime(); Real sinTime=(t-tprev)*10-2*addTime;
    printf("sin %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        cosRes += std::cos(ir+0.001); // two adds
        cosRes += std::cos(ir+0.1);
        cosRes += std::cos(ir+0.2);
        cosRes += std::cos(ir+0.3);
        cosRes += std::cos(ir+0.4);
        cosRes += std::cos(ir+0.501);
        cosRes += std::cos(ir+0.6);
        cosRes += std::cos(ir+0.7);
        cosRes += std::cos(ir+0.8);
        cosRes += std::cos(ir+0.9);
    }
    t = threadCpuTime(); Real cosTime=(t-tprev)*10-2*addTime;
    printf("cos %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        atan2Res += std::atan2(ir+0.001,ir-0.001); // three adds
        atan2Res += std::atan2(ir+0.1,ir-0.1);
        atan2Res += std::atan2(ir+0.2,ir-0.2);
        atan2Res += std::atan2(ir+0.3,ir-0.3);
        atan2Res += std::atan2(ir+0.4,ir-0.4);
        atan2Res += std::atan2(ir+0.501,ir-0.501);
        atan2Res += std::atan2(ir+0.6,ir-0.6);
        atan2Res += std::atan2(ir+0.7,ir-0.7);
        atan2Res += std::atan2(ir+0.8,ir-0.8);
        atan2Res += std::atan2(ir+0.9,ir-0.9);
    }
    t = threadCpuTime(); Real atan2Time=(t-tprev)*10-3*addTime;
    printf("atan2 %gs\n", t-tprev);
    tprev = threadCpuTime();

    std::cout << std::setprecision(5);
    std::cout << "1 flop=avg(add,mul)=" << flopTime << "ns\n";
    printf("op\t t/10^9\t flops\t final result\n");
    std::cout << "int+:\t"  <<intAddTime<<"\t"<<intAddTime   /flopTime<<"\t"<<intAddRes<< "\n";
    std::cout << "+:\t"     <<addTime<<"\t"<<addTime   /flopTime<<"\t"<<addRes<< "\n";
    std::cout << "-:\t"     <<subTime<<"\t"<<subTime   /flopTime<<"\t"<<subRes<< "\n";
    std::cout << "*:\t"     <<mulTime<<"\t"<<mulTime   /flopTime<<"\t"<<mulRes<< "\n";
    std::cout << "/:\t"     <<divTime<<"\t"<<divTime   /flopTime<<"\t"<<divRes<< "\n";
    std::cout << "sqrt:\t"  <<sqrtTime<<"\t"<<sqrtTime  /flopTime<<"\t"<<sqrtRes<< "\n";
    std::cout << "1/sqrt:\t"<<oosqrtTime<<"\t"<<oosqrtTime/flopTime<<"\t"<<oosqrtRes<< "\n";
    std::cout << "log:\t"   <<logTime<<"\t"<<logTime   /flopTime<<"\t"<<logRes<< "\n";
    std::cout << "exp:\t"   <<expTime<<"\t"<<expTime   /flopTime<<"\t"<<expRes<< "\n";
    std::cout << "sin:\t"   <<sinTime<<"\t"<<sinTime   /flopTime<<"\t"<<sinRes<< "\n";
    std::cout << "cos:\t"   <<cosTime<<"\t"<<cosTime   /flopTime<<"\t"<<cosRes<< "\n";
    std::cout << "atan2:\t" <<atan2Time<<"\t"<<atan2Time /flopTime<<"\t"<<atan2Res<< "\n";
}


int main() {
    time_t now;
    time(&now);
    printf("Starting: %s\n", ctime(&now));
    {   std::cout << "\nCPU performance\n" << std::endl;
        testFunctions(flopTimeInNs, true /*flop time only*/);
    }
    Real startClock  = realTime();
    Real startCpu    = cpuTime();
    Real startThread = threadCpuTime();

    {
        std::cout << "\nParticles:\n" << std::endl;
        MultibodySystem system;
        createParticles(system);
        runAllTests(system);
    }
    {
        std::cout << "\nFree Bodies (Quaternions):\n" << std::endl;
        MultibodySystem system;
        createFreeBodies(system);
        runAllTests(system, false);
    }
    {
        std::cout << "\nFree Bodies (Euler angles):\n" << std::endl;
        MultibodySystem system;
        createFreeBodies(system);
        runAllTests(system, true);
    }
    {
        std::cout << "\nPin Chain:\n" << std::endl;
        MultibodySystem system;
        createPinChain(system);
        runAllTests(system);
    }
    {
        std::cout << "\nSlider Chain:\n" << std::endl;
        MultibodySystem system;
        createSliderChain(system);
        runAllTests(system);
    }
    {
        std::cout << "\nBall Chain (Quaternions):\n" << std::endl;
        MultibodySystem system;
        createBallChain(system);
        runAllTests(system, false);
    }
    {
        std::cout << "\nBall Chain (Euler angles):\n" << std::endl;
        MultibodySystem system;
        createBallChain(system);
        runAllTests(system, true);
    }
    {
        std::cout << "\nGimbal Chain:\n" << std::endl;
        MultibodySystem system;
        createGimbalChain(system);
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

    std::cout << "Total time:\n";
    std::cout << "  process CPU=" << cpuTime()-startCpu << "s\n";
    std::cout << "  thread CPU =" << threadCpuTime()-startThread << "s\n";
    std::cout << "  real time  =" << realTime()-startClock << "s\n";
}