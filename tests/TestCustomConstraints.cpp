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
#include "../src/ConstraintImpl.h"

using namespace SimTK;
using namespace std;

const int NUM_BODIES = 10;
const Real BOND_LENGTH = 0.5;
const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2, double tol = TOL) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, double tol) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2, double tol) {
    assertEqual(val1[0], val2[0], tol);
    assertEqual(val1[1], val2[1], tol);
}

template <>
void assertEqual(Vector val1, Vector val2, double tol) {
    ASSERT(val1.size() == val2.size());
    for (int i = 0; i < val1.size(); ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

/**
 * A Function that takes a single argument and returns it.
 */

class LinearFunction : public Function<1> {
public:
    Vec<1> calcValue(const Vector& x) const {
        return Vec1(x[0]);
    }
    Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        if (derivComponents.size() == 1)
            return Vec1(1);
        return Vec1(0);
    }
    int getArgumentSize() const {
        return 1;
    }
    int getMaxDerivativeOrder() const {
        return 100;
    }
};

/**
 * A Function that relates three different arguments.
 */

class CompoundFunction : public Function<1> {
public:
    Vec<1> calcValue(const Vector& x) const {
        return Vec1(x[0]+x[1]+x[2]);
//        return Vec1(x[0]*x[1]+x[2]);
    }
    Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        if (derivComponents.size() == 1) {
return Vec1(1);
            switch (derivComponents[0]) {
            case 0:
                return Vec1(x[1]);
            case 1:
                return Vec1(x[0]);
            default:
                return Vec1(1);
            }
        }
        if (derivComponents.size() == 2) {
return Vec1(0);
            int count[] = {0, 0, 0};
            count[derivComponents[0]]++;
            count[derivComponents[1]]++;
            if (count[0] == 1 && count[1] == 1)
                return Vec1(1);
            if (count[0] == 1 && count[2] == 1)
                return Vec1(x[1]+1);
            if (count[1] == 1 && count[2] == 1)
                return Vec1(x[0]+1);
        }
        return Vec1(0);
    }
    int getArgumentSize() const {
        return 3;
    }
    int getMaxDerivativeOrder() const {
        return 2;
    }
};


/**
 * Create a system consisting of a chain of bodies.
 */

void createSystem(MultibodySystem& system) {
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -1, 0), 0);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNBodies()-1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }
}

/**
 * Create a random state for the system.
 */

void createState(MultibodySystem& system, State& state, const Vector& y=Vector()) {
    system.realizeTopology();
    state = system.getDefaultState();
    if (y.size() > 0)
        state.updY() = y;
    else {
        Random::Uniform random;
        for (int i = 0; i < state.getNY(); ++i)
            state.updY()[i] = random.getValue();
    }
    system.realize(state, Stage::Velocity);
    system.project(state, TOL, Vector(state.getNY(), 1), Vector(state.getNYErr(), 1), Vector(state.getNY()));
    system.realize(state, Stage::Acceleration);
}

void testCoordinateCoupler1() {

    // Create a system using three CoordinateCouplers to fix the orientation of one body.
    
    MultibodySystem system1;
    SimbodyMatterSubsystem matter1(system1);
    createSystem(system1);
    MobilizedBody& first = matter1.updMobilizedBody(MobilizedBodyIndex(1));
    std::vector<MobilizedBodyIndex> bodies(1);
    std::vector<MobilizerQIndex> coordinates(1);
    bodies[0] = MobilizedBodyIndex(1);
    coordinates[0] = MobilizerQIndex(0);
    Constraint::CoordinateCoupler coupler1(matter1, new LinearFunction(), bodies, coordinates);
    coordinates[0] = MobilizerQIndex(1);
    Constraint::CoordinateCoupler coupler2(matter1, new LinearFunction(), bodies, coordinates);
    coordinates[0] = MobilizerQIndex(2);
    Constraint::CoordinateCoupler coupler3(matter1, new LinearFunction(), bodies, coordinates);
    State state1;
    createState(system1, state1);

    // Create a system using a ConstantOrientation constraint to do the same thing.
    
    MultibodySystem system2;
    SimbodyMatterSubsystem matter2(system2);
    createSystem(system2);
    Constraint::ConstantOrientation orient(matter2.updGround(), Rotation(), matter2.updMobilizedBody(MobilizedBodyIndex(1)), Rotation());
    State state2;
    createState(system2, state2, state1.getY());
    
    // Compare the results.
    
    assertEqual(state1.getQ(), state2.getQ());
    assertEqual(state1.getQDot(), state2.getQDot());
    assertEqual(state1.getQDotDot(), state2.getQDotDot());
    assertEqual(state1.getU(), state2.getU());
    assertEqual(state1.getUDot(), state2.getUDot());
}

void testCoordinateCoupler2() {
    
    // Create a system involving a constraint that affects three different bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    createSystem(system);
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    std::vector<MobilizedBodyIndex> bodies(1);
    std::vector<MobilizerQIndex> coordinates(1);
    bodies[0] = MobilizedBodyIndex(1);
    bodies[1] = MobilizedBodyIndex(3);
    bodies[2] = MobilizedBodyIndex(5);
    coordinates[0] = MobilizerQIndex(0);
    coordinates[1] = MobilizerQIndex(0);
    coordinates[2] = MobilizerQIndex(1);
    Function<1>* function = new CompoundFunction();
//    Constraint::CoordinateCoupler coupler(matter, function, bodies, coordinates);

    bodies[0] = MobilizedBodyIndex(1);
    coordinates[0] = MobilizerQIndex(0);
    Constraint::CoordinateCoupler coupler1(matter, new LinearFunction(), bodies, coordinates);
//    coordinates[0] = MobilizerQIndex(1);
//    Constraint::CoordinateCoupler coupler2(matter, new LinearFunction(), bodies, coordinates);
//    coordinates[0] = MobilizerQIndex(2);
//    Constraint::CoordinateCoupler coupler3(matter, new LinearFunction(), bodies, coordinates);

    State state;
    createState(system, state, /**/ Vector(60, 1.0));
    
    // Make sure the constraint is satisfied.
    
    Vector args(function->getArgumentSize());
    for (int i = 0; i < args.size(); ++i)
        args[i] = matter.getMobilizedBody(bodies[i]).getOneQ(state, coordinates[i]);
//    assertEqual(0.0, function->calcValue(args)[0]);
    
    // Simulate it and make sure the constraint is working correctly and energy is being conserved.
    
    Real energy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setReturnEveryInternalStep(true);
    integ.initialize(state);
std::cout << "energy: "<<energy << std::endl;
    while (integ.getTime() < 10.0) {
        integ.stepTo(10.0);
std::cout << integ.getTime() << std::endl;
//std::cout << state.getQErr() << std::endl;
//std::cout << state.getUErr() << std::endl;
//std::cout << state.getUDotErr() << std::endl;
        for (int i = 0; i < args.size(); ++i)
            args[i] = matter.getMobilizedBody(bodies[i]).getOneQ(integ.getState(), coordinates[i]);
std::cout << integ.getState().getQ() << std::endl;
std::cout << integ.getState().getU() << std::endl;
std::cout << integ.getState().getUDot() << std::endl;
//std::cout << function->calcValue(args)[0] << std::endl;
std::cout << system.calcEnergy(integ.getState()) << std::endl;
//        assertEqual(0.0, function->calcValue(args)[0], integ.getConstraintToleranceInUse());
//        assertEqual(energy, system.calcEnergy(integ.getState()));
    }
    delete function;
}

void testSpeedCoupler1() {

    // Create a system using a SpeedCoupler to fix one speed.
    
    MultibodySystem system1;
    SimbodyMatterSubsystem matter1(system1);
    createSystem(system1);
    MobilizedBody& first = matter1.updMobilizedBody(MobilizedBodyIndex(1));
    std::vector<MobilizedBodyIndex> bodies(1);
    std::vector<MobilizerUIndex> speeds(1);
    bodies[0] = MobilizedBodyIndex(1);
    speeds[0] = MobilizerUIndex(2);
    Constraint::SpeedCoupler coupler1(matter1, new LinearFunction(), bodies, speeds);
    State state1;
    createState(system1, state1);

    // Create a system using a ConstantSpeed constraint to do the same thing.
    
    MultibodySystem system2;
    SimbodyMatterSubsystem matter2(system2);
    createSystem(system2);
    Constraint::ConstantSpeed orient(matter2.updMobilizedBody(MobilizedBodyIndex(1)), MobilizerUIndex(2), 0);
    State state2;
    createState(system2, state2, state1.getY());
    
    // Compare the results.
    
    assertEqual(state1.getQ(), state2.getQ());
    assertEqual(state1.getQDot(), state2.getQDot());
    assertEqual(state1.getQDotDot(), state2.getQDotDot());
    assertEqual(state1.getU(), state2.getU());
    assertEqual(state1.getUDot(), state2.getUDot());
}

void testSpeedCoupler2() {
    
    // Create a system involving a constraint that affects three different bodies.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    createSystem(system);
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    std::vector<MobilizedBodyIndex> bodies(3);
    std::vector<MobilizerUIndex> speeds(3);
    bodies[0] = MobilizedBodyIndex(1);
    bodies[1] = MobilizedBodyIndex(3);
    bodies[2] = MobilizedBodyIndex(5);
    speeds[0] = MobilizerUIndex(0);
    speeds[1] = MobilizerUIndex(0);
    speeds[2] = MobilizerUIndex(1);
    Function<1>* function = new CompoundFunction();
    Constraint::SpeedCoupler coupler(matter, function, bodies, speeds);
    State state;
    createState(system, state);
    
    // Make sure the constraint is satisfied.
    
    Vector args(function->getArgumentSize());
    for (int i = 0; i < args.size(); ++i)
        args[i] = matter.getMobilizedBody(bodies[i]).getOneU(state, speeds[i]);
    assertEqual(0.0, function->calcValue(args)[0]);
    
    // Simulate it and make sure the constraint is working correctly and energy is being conserved.
    
    Real energy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setReturnEveryInternalStep(true);
    integ.initialize(state);
    while (integ.getTime() < 10.0) {
        integ.stepTo(10.0);
        for (int i = 0; i < args.size(); ++i)
            args[i] = matter.getMobilizedBody(bodies[i]).getOneU(integ.getState(), speeds[i]);
        assertEqual(0.0, function->calcValue(args)[0], integ.getConstraintToleranceInUse());
        assertEqual(energy, system.calcEnergy(integ.getState()), energy*0.05);
    }
    delete function;
}

int main() {
    try {
        testCoordinateCoupler1();
//        testCoordinateCoupler2();
        testSpeedCoupler1();
        testSpeedCoupler2();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
