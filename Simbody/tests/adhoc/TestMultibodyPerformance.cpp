/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKsimbody.h"
#include <string>
#include <ctime>

using std::cout;
using std::endl;
using std::string;

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strerror, sprintf, etc.
#endif

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

void doRealizeTime(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Time);
    system.realize(state, Stage::Time);
}

void doRealizePosition(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Position);
    system.realize(state, Stage::Position);
}

void doRealizeVelocity(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Velocity);
    system.realize(state, Stage::Velocity);
}

void doRealizeArticulatedBodyInertias(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    matter.invalidateArticulatedBodyInertias(state);
    matter.realizeArticulatedBodyInertias(state);
}

void doRealizeDynamics(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    system.realize(state, Stage::Dynamics);
}

void doRealizeAcceleration(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Acceleration);
    system.realize(state, Stage::Acceleration);
}

// Cost to re-evaluate accelerations after applying some new forces, but leaving
// the state variables unchanged.
void doRealizeDynamics2Acceleration(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    system.realize(state, Stage::Acceleration);
}

// Cost to re-evaluate accelerations after updating velocities, but leaving
// the positions unchanged (e.g. semi-implicit Euler iterating velocities).
void doRealizeVelocity2Acceleration(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Velocity);
    system.realize(state, Stage::Acceleration);
}

// Cost for a complete acceleration calculation at a new time and state.
// This includes the cost of articulated body inertias.
void doRealizeTime2Acceleration(MultibodySystem& system, State& state) {
    state.invalidateAllCacheAtOrAbove(Stage::Time);
    system.realize(state, Stage::Acceleration);
}

void doMultiplyByM(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector v(matter.getNumMobilities(), 1.0);
    Vector mv;
    matter.multiplyByM(state, v, mv);
}

void doMultiplyByMInv(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector v(matter.getNumMobilities(), 1.0);
    Vector minvv;
    matter.multiplyByMInv(state, v, minvv);
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

void doMultiplyBySystemJacobianTranspose(MultibodySystem& system, State& state) {
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
    Vector_<SpatialVec> dEdR(matter.getNumBodies(), SpatialVec(Vec3(1, 0, 0), Vec3(0, 1, 0)));
    Vector dEdQ;
    matter.multiplyBySystemJacobianTranspose(state, dEdR, dEdQ);
}

void doCalcCompositeBodyInertias(MultibodySystem& system, State& state) {
    Array_<SpatialInertia, MobilizedBodyIndex> r;
    system.getMatterSubsystem().calcCompositeBodyInertias(state, r);
}

static double flopTimeInNs;

/**
 * Time how long it takes to perform an operation 1000 times.  The test is repeated 5 times,
 * and the average is returned.  The return value represents CPU time, *not* clock time.
 */
void timeComputation(MultibodySystem& system, void function(MultibodySystem& system, State& state), 
                     const string& name, int iterations, bool useEulerAngles) {
    const int repeats = 3;
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
        cpuTimes[i] = Real(endCpu-startCpu);
    }

    double timePerIterUs = mean(cpuTimes)*1000000/iterations; // us
    double flopTimeUs = flopTimeInNs / 1000;
    double flopTimePerIter = timePerIterUs/flopTimeUs;
    std::printf("%40s:%6.4gus -> %4d flp/dof, %4d flp/bod\n",
        name.c_str(), timePerIterUs, (int)(flopTimePerIter/ndof),
        (int)(flopTimePerIter/nmovbod));
}

/**
 * Time all the different calculations for one system.
 */
void runAllTests(MultibodySystem& system, bool useEulerAngles=false) {
    std::cout << "# dofs=" << system.getMatterSubsystem().getNumMobilities() << "\n";
    timeComputation(system, doRealizeTime, "realizeTime", 5000, useEulerAngles);
    timeComputation(system, doRealizePosition, "realizePosition", 5000, useEulerAngles);
    timeComputation(system, doRealizeVelocity, "realizeVelocity", 5000, useEulerAngles);
    timeComputation(system, doRealizeArticulatedBodyInertias, "doRealizeArticulatedBodyInertias", 3000, useEulerAngles);
    timeComputation(system, doRealizeDynamics, "realizeDynamics", 5000, useEulerAngles);
    timeComputation(system, doRealizeAcceleration, "realizeAcceleration", 5000, useEulerAngles);
    timeComputation(system, doRealizeDynamics2Acceleration, "doRealizeDynamics2Acceleration", 5000, useEulerAngles);
    timeComputation(system, doRealizeVelocity2Acceleration, "realizeVelocity2Acceleration", 3000, useEulerAngles);
    timeComputation(system, doRealizeTime2Acceleration, "realizeTime2Acceleration", 2000, useEulerAngles);
    timeComputation(system, doMultiplyByM, "multiplyByM", 5000, useEulerAngles);
    timeComputation(system, doMultiplyByMInv, "multiplyByMInv", 5000, useEulerAngles);
    timeComputation(system, doCalcResidualForceIgnoringConstraints, "calcResidualForceIgnoringConstraints", 5000, useEulerAngles);
    timeComputation(system, doCalcMobilizerReactionForces, "calcMobilizerReactionForces", 1000, useEulerAngles);
    timeComputation(system, doMultiplyBySystemJacobianTranspose, "multiplyBySystemJacobianTranspose", 5000, useEulerAngles);
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
    {Real(0.501),Real(0.2501),Real(0.201),Real(0.101),Real(1.000000001),
    Real(1/1.000000002),Real(1/.101), Real(1/.201), Real(1/.2501), Real(1/.501)};
void testFunctions(double& flopTime, bool flopTimeOnly=false) {
    Real addRes=1,subRes=1,mulRes=1,divRes=1,sqrtRes=1,oosqrtRes=1,
         sinRes=1,cosRes=1,atan2Res=1,logRes=1,expRes=1;
    int intAddRes=1;
    Random::Uniform rand; rand.setMin(-5); rand.setMax(5);
    for (int i=0; i<10; i++) tenInts[i] = rand.getIntValue();
    for (int i=0; i<10; i++) tenReals[i] = rand.getValue();

    double tprev = threadCpuTime();
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
    double t = threadCpuTime(); double intAddTime = (t-tprev)/10; // time for 1e9 ops
    printf("intAdd %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 5*100000000; i++) {
        addRes += Real(1.1);
        addRes += Real(1.2);
        addRes += Real(1.3);
        addRes += Real(1.4);
        addRes += Real(1.501);
        addRes += Real(1.6);
        addRes += Real(1.7);
        addRes += Real(1.8);
        addRes += Real(1.9);
        addRes += Real(2.007);
    }
    t = threadCpuTime(); double addTime = (t-tprev)/5;
    printf("add %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 5*100000000; i++) {
        subRes -= Real(1.1);
        subRes -= Real(1.2);
        subRes -= Real(1.3);
        subRes -= Real(1.4);
        subRes -= Real(1.501);
        subRes -= Real(1.6);
        subRes -= Real(1.7);
        subRes -= Real(1.8);
        subRes -= Real(1.9);
        subRes -= Real(2.007);
    }
    t = threadCpuTime(); double subTime = (t-tprev)/5;
    printf("sub %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 3*100000000; i++) {
        mulRes *= Real(0.501);
        mulRes *= Real(0.2501);
        mulRes *= Real(0.201);
        mulRes *= Real(0.101);
        mulRes *= Real(1.000000001);
        mulRes *= Real(1/1.000000002); // done at compile time
        mulRes *= Real(1/.101);
        mulRes *= Real(1/.201);
        mulRes *= Real(1/.2501);
        mulRes *= Real(1/.501);
    }
    t = threadCpuTime(); double mulTime=(t-tprev)/3;
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
        tenMults[i%10]     = Real(tenMults[i%10]*1.0000000000001); 
        tenMults[(i+5)%10] = Real(tenMults[(i+5)%10]*0.9999999999999); 
    }
    t = threadCpuTime(); double divTime=(t-tprev);
    printf("div %gs\n", t-tprev);
    tprev = threadCpuTime();

    for (int i = 0; i < 100000000/2; i++) {
        const Real ir = (Real)i;
        sqrtRes += std::sqrt(ir+(Real)0.001); // two adds
        sqrtRes += std::sqrt(ir+(Real)0.1);
        sqrtRes += std::sqrt(ir+(Real)0.2);
        sqrtRes += std::sqrt(ir+(Real)0.3);
        sqrtRes += std::sqrt(ir+(Real)0.4);
        sqrtRes += std::sqrt(ir+(Real)0.501);
        sqrtRes += std::sqrt(ir+(Real)0.6);
        sqrtRes += std::sqrt(ir+(Real)0.7);
        sqrtRes += std::sqrt(ir+(Real)0.8);
        sqrtRes += std::sqrt(ir+(Real)0.9);
    }
    t = threadCpuTime(); double sqrtTime=2*(t-tprev)-2*addTime;
    printf("sqrt %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/4; i++) {
        const Real ir = (Real)i;
        oosqrtRes += 1/std::sqrt(ir+(Real)0.001); // two adds
        oosqrtRes += 1/std::sqrt(ir+(Real)0.1);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.2);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.3);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.4);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.501);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.6);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.7);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.8);
        oosqrtRes += 1/std::sqrt(ir+(Real)0.9);
    }
    t = threadCpuTime(); double oosqrtTime=4*(t-tprev)-2*addTime;
    printf("1/sqrt %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/5; i++) {
        const Real ir = (Real)i;
        logRes += std::log(ir+(Real)0.001); // two adds
        logRes += std::log(ir+(Real)0.1);
        logRes += std::log(ir+(Real)0.2);
        logRes += std::log(ir+(Real)0.3);
        logRes += std::log(ir+(Real)0.4);
        logRes += std::log(ir+(Real)0.501);
        logRes += std::log(ir+(Real)0.6);
        logRes += std::log(ir+(Real)0.7);
        logRes += std::log(ir+(Real)0.8);
        logRes += std::log(ir+(Real)0.9);
    }
    t = threadCpuTime(); double logTime=(t-tprev)*5-2*addTime;
    printf("log %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/5; i++) {
        const Real ir = (Real).000000001*(Real)i;
        expRes += std::exp(ir+(Real)0.001); // two adds
        expRes += std::exp(ir+(Real)0.1);
        expRes += std::exp(ir+(Real)0.2);
        expRes += std::exp(ir+(Real)0.3);
        expRes += std::exp(ir+(Real)0.4);
        expRes += std::exp(ir+(Real)0.501);
        expRes += std::exp(ir+(Real)0.6);
        expRes += std::exp(ir+(Real)0.7);
        expRes += std::exp(ir+(Real)0.8);
        expRes += std::exp(ir+(Real)0.9);
    }
    t = threadCpuTime(); double expTime=(t-tprev)*5-2*addTime;
    printf("exp %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        sinRes += std::sin(ir+(Real)0.001); // two adds
        sinRes += std::sin(ir+(Real)0.1);
        sinRes += std::sin(ir+(Real)0.2);
        sinRes += std::sin(ir+(Real)0.3);
        sinRes += std::sin(ir+(Real)0.4);
        sinRes += std::sin(ir+(Real)0.501);
        sinRes += std::sin(ir+(Real)0.6);
        sinRes += std::sin(ir+(Real)0.7);
        sinRes += std::sin(ir+(Real)0.8);
        sinRes += std::sin(ir+(Real)0.9);
    }
    t = threadCpuTime(); double sinTime=(t-tprev)*10-2*addTime;
    printf("sin %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        cosRes += std::cos(ir+(Real)0.001); // two adds
        cosRes += std::cos(ir+(Real)0.1);
        cosRes += std::cos(ir+(Real)0.2);
        cosRes += std::cos(ir+(Real)0.3);
        cosRes += std::cos(ir+(Real)0.4);
        cosRes += std::cos(ir+(Real)0.501);
        cosRes += std::cos(ir+(Real)0.6);
        cosRes += std::cos(ir+(Real)0.7);
        cosRes += std::cos(ir+(Real)0.8);
        cosRes += std::cos(ir+(Real)0.9);
    }
    t = threadCpuTime(); double cosTime=(t-tprev)*10-2*addTime;
    printf("cos %gs\n", t-tprev);
    tprev = threadCpuTime();
    for (int i = 0; i < 100000000/10; i++) {
        const Real ir = (Real)i;
        atan2Res += std::atan2(ir+(Real)0.001,ir-(Real)0.001); // three adds
        atan2Res += std::atan2(ir+(Real)0.1,ir-(Real)0.1);
        atan2Res += std::atan2(ir+(Real)0.2,ir-(Real)0.2);
        atan2Res += std::atan2(ir+(Real)0.3,ir-(Real)0.3);
        atan2Res += std::atan2(ir+(Real)0.4,ir-(Real)0.4);
        atan2Res += std::atan2(ir+(Real)0.501,ir-(Real)0.501);
        atan2Res += std::atan2(ir+(Real)0.6,ir-(Real)0.6);
        atan2Res += std::atan2(ir+(Real)0.7,ir-(Real)0.7);
        atan2Res += std::atan2(ir+(Real)0.8,ir-(Real)0.8);
        atan2Res += std::atan2(ir+(Real)0.9,ir-(Real)0.9);
    }
    t = threadCpuTime(); double atan2Time=(t-tprev)*10-3*addTime;
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
    double startClock  = realTime();
    double startCpu    = cpuTime();
    double startThread = threadCpuTime();

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