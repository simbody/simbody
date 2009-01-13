/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
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
#include "../src/ForceImpl.h"

using namespace SimTK;
using namespace std;

const int NUM_BODIES = 10;
const Real BOND_LENGTH = 0.5;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

void verifyForces(const Force& force, const State& state, Vector_<SpatialVec> bodyForces, Vector_<Vec3> particleForces, Vector mobilityForces) {
    Vector_<SpatialVec> actualBodyForces(bodyForces.size());
    Vector_<Vec3> actualParticleForces(particleForces.size());
    Vector actualMobilityForces(mobilityForces.size());
    actualBodyForces = SpatialVec(Vec3(0), Vec3(0));
    actualParticleForces = Vec3(0);
    actualMobilityForces = 0;
    force.getImpl().calcForce(state, actualBodyForces, actualParticleForces, actualMobilityForces);
    for (int i = 0; i < bodyForces.size(); ++i)
        ASSERT((bodyForces[i]-actualBodyForces[i]).norm() < 1e-10);
    for (int i = 0; i < particleForces.size(); ++i)
        ASSERT((particleForces[i]-actualParticleForces[i]).norm() < 1e-10);
    for (int i = 0; i < mobilityForces.size(); ++i)
        ASSERT(std::abs(mobilityForces[i]-actualMobilityForces[i]) < 1e-10);
}

class MyForceImpl : public Force::Custom::Implementation {
public:
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
        for (int i = 0; i < mobilityForces.size(); ++i)
            mobilityForces[i] += i;
    }
    Real calcPotentialEnergy(const State& state) const {
        return 0.0;
    }
};

/**
 * Test all of the standard Force subclasses, and make sure they generate correct forces.
 */

void testStandardForces() {
    
    // Create a system consisting of a chain of bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNBodies()-1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }
    
    // Add one of each type of force.
    
    MobilizedBody& body1 = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& body9 = matter.updMobilizedBody(MobilizedBodyIndex(9));
    Force::ConstantForce constantForce(forces, body1, Vec3(0), Vec3(1, 2, 3));
    Force::ConstantTorque constantTorque(forces, body1, Vec3(1, 2, 3));
    Force::GlobalDamper globalDamper(forces, matter, 2.0);
    Force::MobilityConstantForce mobilityConstantForce(forces, body1, 1, 2.0);
    Force::MobilityLinearDamper mobilityLinearDamper(forces, body1, 1, 2.0);
    Force::MobilityLinearSpring mobilityLinearSpring(forces, body1, 1, 2.0, 1.0);
    Force::TwoPointConstantForce twoPointConstantForce(forces, body1, Vec3(0), body9, Vec3(0), 2.0);
    Force::TwoPointLinearDamper twoPointLinearDamper(forces, body1, Vec3(0), body9, Vec3(0), 2.0);
    Force::TwoPointLinearSpring twoPointLinearSpring(forces, body1, Vec3(0), body9, Vec3(0), 2.0, 0.5);
    Force::UniformGravity uniformGravity(forces, matter, Vec3(0, -2.0, 0));
    Force::Custom custom(forces, new MyForceImpl());

    // Create a random state for it.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Uniform random;
    for (int i = 0; i < state.getNY(); ++i)
        state.updY()[i] = random.getValue();
    system.realize(state, Stage::Velocity);
    Vec3 pos1 = body1.getBodyOriginLocation(state);
    Vec3 pos9 = body9.getBodyOriginLocation(state);
    Vec3 delta19 = pos9-pos1;
    
    // Calculate each force component and see if it is correct.
    
    Vector_<SpatialVec> bodyForces(matter.getNBodies());
    Vector_<Vec3> particleForces(0);
    Vector mobilityForces(state.getNU());
    Real pe = 0;
    
    // Check ConstantForce
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    bodyForces[1][1] = Vec3(1, 2, 3);
    verifyForces(constantForce, state, bodyForces, particleForces, mobilityForces);
    
    // Check ConstantTorque
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    bodyForces[1][0] = Vec3(1, 2, 3);
    verifyForces(constantTorque, state, bodyForces, particleForces, mobilityForces);
    
    // Check GlobalDamper
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = -2.0*state.getU();
    verifyForces(globalDamper, state, bodyForces, particleForces, mobilityForces);
    
    // Check MobilityConstantForce
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    body1.updOneFromUPartition(state, 1, mobilityForces) = 2.0;
    verifyForces(mobilityConstantForce, state, bodyForces, particleForces, mobilityForces);
    
    // Check MobilityLinearDamper
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    body1.updOneFromUPartition(state, 1, mobilityForces) = -2.0*body1.getOneU(state, 1);
    verifyForces(mobilityLinearDamper, state, bodyForces, particleForces, mobilityForces);
    
    // Check MobilityLinearSpring
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    body1.updOneFromUPartition(state, 1, mobilityForces) = -2.0*(body1.getOneQ(state, 1)-1.0);
    verifyForces(mobilityLinearSpring, state, bodyForces, particleForces, mobilityForces);
    
    // Check TwoPointConstantForce
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    bodyForces[1][1] = -2.0*delta19.normalize();
    bodyForces[9][1] = 2.0*delta19.normalize();
    verifyForces(twoPointConstantForce, state, bodyForces, particleForces, mobilityForces);
    
    // Check TwoPointLinearDamper
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    Vec3 v19 = body9.getBodyOriginVelocity(state)-body1.getBodyOriginVelocity(state);
    Real twoPointLinearDamperForce = 2.0*dot(v19, delta19.normalize());
    bodyForces[1][1] = twoPointLinearDamperForce*delta19.normalize();
    bodyForces[9][1] = -twoPointLinearDamperForce*delta19.normalize();
    verifyForces(twoPointLinearDamper, state, bodyForces, particleForces, mobilityForces);
    
    // Check TwoPointLinearSpring
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    Real twoPointLinearSpringForce = 2.0*(delta19.norm()-0.5);
    bodyForces[1][1] = twoPointLinearSpringForce*delta19.normalize();
    bodyForces[9][1] = -twoPointLinearSpringForce*delta19.normalize();
    verifyForces(twoPointLinearSpring, state, bodyForces, particleForces, mobilityForces);
    
    // Check UniformGravity
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0, -2.0, 0));
    bodyForces[0] = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    mobilityForces = 0;
    verifyForces(uniformGravity, state, bodyForces, particleForces, mobilityForces);

    // Check Custom
    
    bodyForces = SpatialVec(Vec3(0), Vec3(0));
    particleForces = Vec3(0);
    for (int i = 0; i < mobilityForces.size(); ++i)
        mobilityForces[i] = i;
    verifyForces(custom, state, bodyForces, particleForces, mobilityForces);
}

/**
 * Test the standard conservative forces to make sure they really conserve energy.
 */

void testEnergyConservation() {
    
    // Create a system consisting of a chain of bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNBodies()-1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }
    
    // Add one of each type of conservative force.
    
    MobilizedBody& body1 = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& body9 = matter.updMobilizedBody(MobilizedBodyIndex(9));
    Force::MobilityLinearSpring mobilityLinearSpring(forces, body1, 1, 0.1, 1.0);
    Force::TwoPointLinearSpring twoPointLinearSpring(forces, body1, Vec3(0), body9, Vec3(0), 1.0, 4.0);
    Force::UniformGravity uniformGravity(forces, matter, Vec3(0, -1.0, 0));

    // Create a random initial state for it.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Uniform random;
    for (int i = 0; i < state.getNY(); ++i)
        state.updY()[i] = random.getValue();
    
    // Simulate it for a while and see if the energy changes.
    
    system.realize(state, Stage::Dynamics);
    Real initialEnergy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-4);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);
    system.realize(state, Stage::Dynamics);
    Real finalEnergy = system.calcEnergy(ts.getState());
    ASSERT(std::abs(initialEnergy/finalEnergy-1.0) < 0.005);
}

int main() {
    try {
        testStandardForces();
        testEnergyConservation();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
