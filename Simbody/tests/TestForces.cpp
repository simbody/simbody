/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const int NUM_BODIES = 10;
const Real BOND_LENGTH = 0.5;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#define ASSERT_EQUAL(val1, val2) {ASSERT(std::abs(val1-val2) < 1e-10);}

void verifyForces(const Force& force, const State& state, Vector_<SpatialVec> bodyForces, Vector_<Vec3> particleForces, Vector mobilityForces) {
    Vector_<SpatialVec> actualBodyForces(bodyForces.size());
    Vector_<Vec3> actualParticleForces(particleForces.size());
    Vector actualMobilityForces(mobilityForces.size());
    force.calcForceContribution(state, actualBodyForces, actualParticleForces, 
                                actualMobilityForces);

    for (int i = 0; i < bodyForces.size(); ++i)
        ASSERT((bodyForces[i]-actualBodyForces[i]).norm() < 1e-10);
    for (int i = 0; i < particleForces.size(); ++i)
        ASSERT((particleForces[i]-actualParticleForces[i]).norm() < 1e-10);
    for (int i = 0; i < mobilityForces.size(); ++i)
        ASSERT(std::abs(mobilityForces[i]-actualMobilityForces[i]) < 1e-10);
}

class MyForceImpl : public Force::Custom::Implementation {
public:
    mutable bool hasRealized[Stage::Report+1];
    MyForceImpl() {
        for (int i = 0; i < Stage::NValid; i++)
            hasRealized[i] = false;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const override {
        for (int i = 0; i < mobilityForces.size(); ++i)
            mobilityForces[i] += i;
    }
    Real calcPotentialEnergy(const State& state) const override {
        return 0.0;
    }
    void realizeTopology(State& state) const override {
        hasRealized[Stage::Topology] = true;
    }
    void realizeModel(State& state) const override {
        hasRealized[Stage::Model] = true;
    }
    void realizeInstance(const State& state) const override {
        hasRealized[Stage::Instance] = true;
    }
    void realizeTime(const State& state) const override {
        hasRealized[Stage::Time] = true;
    }
    void realizePosition(const State& state) const override {
        hasRealized[Stage::Position] = true;
    }
    void realizeVelocity(const State& state) const override {
        hasRealized[Stage::Velocity] = true;
    }
    void realizeDynamics(const State& state) const override {
        hasRealized[Stage::Dynamics] = true;
    }
    void realizeAcceleration(const State& state) const override {
        hasRealized[Stage::Acceleration] = true;
    }
    void realizeReport(const State& state) const override {
        hasRealized[Stage::Report] = true;
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
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNumBodies()-1));
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
    
    Vector_<SpatialVec> bodyForces(matter.getNumBodies());
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
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNumBodies()-1));
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

/**
 * Make sure that all the "realize" methods on a custom force actually get called.
 */

void testCustomRealization() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    MyForceImpl* impl = new MyForceImpl();
    Force::Custom custom(forces, impl);
    State state = system.realizeTopology();
    for (Stage j = Stage::Model; j <= Stage::Report; j++) {
        system.realize(state, j);
        for (Stage i = Stage::Topology; i <= Stage::Report; i++)
            ASSERT(impl->hasRealized[i] == (i <= j));
    }
}

/**
 * Test enabling and disabling forces.
 */

void testDisabling() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free body1(matter.updGround(), Vec3(0), body, Vec3(0));
    MobilizedBody::Free body2(matter.updGround(), Vec3(0), body, Vec3(0));
    Force::TwoPointLinearSpring spring(forces, body1, Vec3(0), body2, Vec3(0), 2.0, 0.5);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -2.0, 0));

    // Create an initial state.
    
    State state = system.realizeTopology();
    body1.setQToFitTranslation(state, Vec3(0, 1, 0));
    body2.setQToFitTranslation(state, Vec3(1, 1, 0));
    
    // These are the contribution of each force to the energy and to the force on body1.
    
    Real springEnergy = 0.5*2.0*0.5*0.5;
    SpatialVec springForce(Vec3(0), Vec3(2.0*0.5, 0, 0));
    Real gravityEnergy = 2*2.0;
    SpatialVec gravityForce(Vec3(0), Vec3(0, -2.0, 0));
    
    // Verify the force and energy for each combination of the forces being enabled or disabled.
    
    system.realize(state, Stage::Dynamics);
    ASSERT_EQUAL(springEnergy+gravityEnergy, system.calcEnergy(state));
    ASSERT((springForce+gravityForce-system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm() < 1e-10);
    ASSERT(!forces.isForceDisabled(state, gravity.getForceIndex()));
    ASSERT(!forces.isForceDisabled(state, spring.getForceIndex()));
    forces.setForceIsDisabled(state, spring.getForceIndex(), true);
    system.realize(state, Stage::Dynamics);
    ASSERT_EQUAL(gravityEnergy, system.calcEnergy(state));
    ASSERT((gravityForce-system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm() < 1e-10);
    ASSERT(!forces.isForceDisabled(state, gravity.getForceIndex()));
    ASSERT(forces.isForceDisabled(state, spring.getForceIndex()));
    forces.setForceIsDisabled(state, gravity.getForceIndex(), true);
    system.realize(state, Stage::Dynamics);
    ASSERT_EQUAL(0, system.calcEnergy(state));
    ASSERT((system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm() < 1e-10);
    ASSERT(forces.isForceDisabled(state, gravity.getForceIndex()));
    ASSERT(forces.isForceDisabled(state, spring.getForceIndex()));
    forces.setForceIsDisabled(state, spring.getForceIndex(), false);
    system.realize(state, Stage::Dynamics);
    ASSERT_EQUAL(springEnergy, system.calcEnergy(state));
    ASSERT((springForce-system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm() < 1e-10);
    ASSERT(forces.isForceDisabled(state, gravity.getForceIndex()));
    ASSERT(!forces.isForceDisabled(state, spring.getForceIndex()));
}

int main() {
    try {
        testStandardForces();
        testEnergyConservation();
        testCustomRealization();
        testDisabling();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
