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
#include "SimTKcommon/Testing.h"

#include "../src/ConstraintImpl.h"

using namespace SimTK;
using namespace std;

const int NUM_BODIES = 10;
const Real BOND_LENGTH = 0.5;

// Keep constraints satisfied to this tolerance during testing.
static const Real ConstraintTol = 1e-10;

// Compare two quantities that are expected to agree to constraint tolerance.
#define CONSTRAINT_TEST(a,b) SimTK_TEST_EQ_TOL(a,b, ConstraintTol)

// Compare two quantities that should have been calculated to machine tolerance
// given the problem size, which we'll characterize by the number of mobilities.
#define MACHINE_TEST(a,b) SimTK_TEST_EQ_SIZE(a,b, state.getNU())


/**
 * Create a system consisting of a chain of bodies.
 */

MultibodySystem& createSystem() {
    MultibodySystem* system = new MultibodySystem();
    SimbodyMatterSubsystem matter(*system);
    GeneralForceSubsystem forces(*system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNumBodies()-1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }
    return *system;
}

/**
 * Create a random state for the system.
 */

void createState(MultibodySystem& system, State& state, const Vector& qOverride=Vector()) {
    system.realizeTopology();
    state = system.getDefaultState();
    Random::Uniform random;
    for (int i = 0; i < state.getNY(); ++i)
        state.updY()[i] = random.getValue();
    if (qOverride.size())
        state.updQ() = qOverride;
    system.realize(state, Stage::Velocity);

    Vector temp(state.getNY());
    system.project(state, ConstraintTol, Vector(state.getNY(), 1), Vector(state.getNYErr(), 1), temp);
    system.realize(state, Stage::Acceleration);
}

void testBallConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::Ball constraint(first, last);
    createState(system, state);
    Vec3 r1 = first.getBodyOriginLocation(state);
    Vec3 r2 = last.getBodyOriginLocation(state);
    Vec3 dr = r1-r2;
    Vec3 v1 = first.getBodyOriginVelocity(state);
    Vec3 v2 = last.getBodyOriginVelocity(state);
    Vec3 dv = v1-v2;
    Vec3 a1 = first.getBodyOriginAcceleration(state);
    Vec3 a2 = last.getBodyOriginAcceleration(state);
    Vec3 da = a1-a2;
    CONSTRAINT_TEST(dr.norm(), 0.0);
    CONSTRAINT_TEST(dr, constraint.getPositionErrors(state));
    CONSTRAINT_TEST(dv.norm(), 0.0);
    CONSTRAINT_TEST(dv, constraint.getVelocityErrors(state));
    // Accelerations should be satisfied to machine tolerance times the
    // size of the problem.
    MACHINE_TEST(da.norm(), 0);
    MACHINE_TEST(da, constraint.getAccelerationErrors(state));
    delete &system;
}

void testConstantAngleConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::ConstantAngle constraint(first, UnitVec3(1, 0, 0), last, UnitVec3(0, 1, 0), 1.1);
    createState(system, state);
    Vec3 dir1 = first.getBodyRotation(state)*Vec3(1, 0, 0);
    Vec3 dir2 = last.getBodyRotation(state)*Vec3(0, 1, 0);
    Vec3 perpDir = dir1%dir2;
    Vec3 v1 = first.getBodyAngularVelocity(state);
    Vec3 v2 = last.getBodyAngularVelocity(state);
    CONSTRAINT_TEST(dot(dir1, dir2), cos(1.1));
    CONSTRAINT_TEST(dot(dir1, dir2)-cos(1.1), constraint.getPositionError(state));
    CONSTRAINT_TEST(dot(v1-v2, perpDir), 0.0);
    CONSTRAINT_TEST(dot(v1-v2, perpDir), constraint.getVelocityError(state));
    delete &system;
}

void testConstantOrientationConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Rotation r1(0.2, CoordinateAxis::XCoordinateAxis());
    Rotation r2(Pi/2, CoordinateAxis::YCoordinateAxis());
    Constraint::ConstantOrientation constraint(first, r1, last, r2);
    createState(system, state);
    Rotation R_G1 = first.getBodyRotation(state);
    Rotation R_G2 = last.getBodyRotation(state);
    Vec3 v1 = first.getBodyAngularVelocity(state);
    Vec3 v2 = last.getBodyAngularVelocity(state);
    Vec3 a1 = first.getBodyAngularAcceleration(state);
    Vec3 a2 = last.getBodyAngularAcceleration(state);
    
    // Extra constraints are required for assembly.  Without them, this constraint only guarantees
    // that the second body's X/Y/Z axis is perpendicular to the first body's Y/Z/X axis.
    
    // Careful: constraint is x2 perp y1, y2 perp z1, z2 perp x1; this isn't the
    // same if bodies are interchanged.
    CONSTRAINT_TEST(dot(R_G2*r2*Vec3(1, 0, 0), R_G1*r1*Vec3(0, 1, 0)), 0.0);
    CONSTRAINT_TEST(dot(R_G2*r2*Vec3(0, 1, 0), R_G1*r1*Vec3(0, 0, 1)), 0.0);
    CONSTRAINT_TEST(dot(R_G2*r2*Vec3(0, 0, 1), R_G1*r1*Vec3(1, 0, 0)), 0.0);
    CONSTRAINT_TEST(v1, v2);
    // Should match to machine precision for this size problem.
    MACHINE_TEST(a1, a2);
    delete &system;
}

void testConstantSpeedConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    Constraint::ConstantSpeed constraint(first, MobilizerUIndex(1), 0.8);
    createState(system, state);
    CONSTRAINT_TEST(first.getOneU(state, 1), 0.8);
    CONSTRAINT_TEST(first.getOneU(state, 1)-0.8, constraint.getVelocityError(state));
    MACHINE_TEST(first.getOneUDot(state, 1), 0.0);
    MACHINE_TEST(first.getOneUDot(state, 1), constraint.getAccelerationError(state));
    delete &system;
}

void testNoSlip1DConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Vec3 p(1, 0.5, -2.0);
    UnitVec3 n(0, 1, 0);
    Constraint::NoSlip1D constraint(matter.Ground(), p, n, first, last);
    createState(system, state);
    Vec3 p1 = first.findStationAtGroundPoint(state, p);
    Vec3 p2 = last.findStationAtGroundPoint(state, p);
    Vec3 v1 = first.findStationVelocityInGround(state, p1);
    Vec3 v2 = last.findStationVelocityInGround(state, p2);
    Vec3 a1 = first.findStationAccelerationInGround(state, p1);
    Vec3 a2 = last.findStationAccelerationInGround(state, p2);
    CONSTRAINT_TEST(dot(v1, n), dot(v2, n));
    CONSTRAINT_TEST(dot(v1-v2, n), constraint.getVelocityError(state));
    MACHINE_TEST(dot(a1-a2, n), 0.0);
    MACHINE_TEST(dot(a1-a2, n), constraint.getAccelerationError(state));
    delete &system;
}

void testPointInPlaneConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    UnitVec3 normal(1, 0.5, 0);
    Real height = 2.0;
    Vec3 p(1.0, 2.5, -3.0);
    Constraint::PointInPlane constraint(first, normal, height, last, p);
    createState(system, state);
    Vec3 p1 = last.findStationLocationInAnotherBody(state, p, first);
    Vec3 v1 = last.findStationVelocityInAnotherBody(state, p, first);
    CONSTRAINT_TEST(dot(p1, normal), height);
    CONSTRAINT_TEST(dot(p1, normal)-height, constraint.getPositionError(state));
    CONSTRAINT_TEST(dot(v1, normal), 0.0);
    CONSTRAINT_TEST(dot(v1, normal), constraint.getVelocityError(state));    
    delete &system;
}

void testPointOnLineConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    UnitVec3 dir(1, 0.5, 0);
    Vec3 base(0.5, -0.5, 2.0);
    Vec3 p(1.0, 2.5, -3.0);
    Constraint::PointOnLine constraint(first, dir, base, last, p);
    createState(system, state);
    Vec3 p1 = last.findStationLocationInAnotherBody(state, p, first);
    Vec3 v1 = last.findStationVelocityInAnotherBody(state, p, first);
    CONSTRAINT_TEST(cross(p1-base, dir).norm(), 0.0);
    CONSTRAINT_TEST(cross(v1, dir).norm(), 0.0);
    delete &system;
}

void testRodConstraint() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::Rod constraint(first, last, 3.0);
    createState(system, state);
    Vec3 r1 = first.getBodyOriginLocation(state);
    Vec3 r2 = last.getBodyOriginLocation(state);
    Vec3 dr = r1-r2;
    Vec3 v1 = first.getBodyOriginVelocity(state);
    Vec3 v2 = last.getBodyOriginVelocity(state);
    Vec3 dv = v1-v2;
    Vec3 a1 = first.getBodyOriginAcceleration(state);
    Vec3 a2 = last.getBodyOriginAcceleration(state);
    Vec3 da = a1-a2;
    CONSTRAINT_TEST(dr.norm(), 3.0);
    CONSTRAINT_TEST(dr.norm()-3.0, constraint.getPositionError(state));
    CONSTRAINT_TEST(dot(dv, dr), 0.0);
    CONSTRAINT_TEST(dot(dv, dr), constraint.getVelocityError(state));
    MACHINE_TEST(dot(da, dr), -dv.normSqr());
    MACHINE_TEST(dot(da, dr)+dv.normSqr(), constraint.getAccelerationError(state));
    delete &system;
}

void testWeldConstraint() {

    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::Weld constraint(first, last);
    createState(system, state);
    CONSTRAINT_TEST(first.getBodyOriginLocation(state), last.getBodyOriginLocation(state));
    CONSTRAINT_TEST(first.getBodyVelocity(state), last.getBodyVelocity(state));
    MACHINE_TEST(first.getBodyAcceleration(state), last.getBodyAcceleration(state));

    const Rotation& R_G1 = first.getBodyRotation(state);
    const Rotation& R_G2 = last.getBodyRotation(state);
    
    // Extra constraints are required for assembly.  Without them, this constraint only guarantees
    // that the second body's X/Y/Z axis is perpendicular to the first body's Y/Z/X axis.
    
    // Careful: constraint is x2 perp y1, y2 perp z1, z2 perp x1; this isn't the
    // same if bodies are interchanged.
    CONSTRAINT_TEST(dot(R_G2*Vec3(1, 0, 0), R_G1*Vec3(0, 1, 0)), 0.0);
    CONSTRAINT_TEST(dot(R_G2*Vec3(0, 1, 0), R_G1*Vec3(0, 0, 1)), 0.0);
    CONSTRAINT_TEST(dot(R_G2*Vec3(0, 0, 1), R_G1*Vec3(1, 0, 0)), 0.0);
    delete &system;
}

void testWeldConstraintWithPreAssembly() {
    
    // Different constraints are required for assembly.  Without them, this constraint only guarantees
    // that the second body's X/Y/Z axis is perpendicular to the first body's Y/Z/X axis. Here we'll
    // attempt to point all the axes in roughly the right direction prior to Weld-ing them.

    State assemblyState;
    MultibodySystem& assemblySystem = createSystem();
    SimbodyMatterSubsystem& assemblyMatter = assemblySystem.updMatterSubsystem();
    MobilizedBody& afirst = assemblyMatter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& alast = assemblyMatter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::ConstantAngle(afirst, UnitVec3(1,0,0), alast, UnitVec3(1,0,0), 0.5); // 30 degrees
    Constraint::ConstantAngle(afirst, UnitVec3(0,1,0), alast, UnitVec3(0,1,0), 0.5);
    Constraint::ConstantAngle(afirst, UnitVec3(0,0,1), alast, UnitVec3(0,0,1), 0.5);
    Constraint::Ball(afirst, alast); // take care of translation
    createState(assemblySystem, assemblyState);
    delete &assemblySystem;

    // Now rebuild the system using a Weld instead of the three angle constraints, but
    // transfer the q's from the State we calculated above to use as a starting guess.

    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::Weld constraint(first, last);
    createState(system, state, assemblyState.getQ());
    CONSTRAINT_TEST(first.getBodyOriginLocation(state), last.getBodyOriginLocation(state));
    CONSTRAINT_TEST(first.getBodyVelocity(state), last.getBodyVelocity(state));
    MACHINE_TEST(first.getBodyAcceleration(state), last.getBodyAcceleration(state));

    const Rotation& R_G1 = first.getBodyRotation(state);
    const Rotation& R_G2 = last.getBodyRotation(state);

    // This is a much more stringent requirement than the one we can ask for without
    // the preassembly step. Here we expect the frames to be perfectly aligned.
    CONSTRAINT_TEST(~R_G1.x()*R_G2.x(), 1.);
    CONSTRAINT_TEST(~R_G1.y()*R_G2.y(), 1.);
    CONSTRAINT_TEST(~R_G1.z()*R_G2.z(), 1.);

    // Just for fun -- compare Rotation matrices using pointing error.
   // ASSERT(R_G1.isSameRotationToWithinAngle(R_G2, TOL));
    CONSTRAINT_TEST(R_G1, R_G2);

    delete &system;
}

void testDisablingConstraints() {
    
    State state;
    MultibodySystem& system = createSystem();
    SimbodyMatterSubsystem& matter = system.updMatterSubsystem();
    MobilizedBody& first = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& last = matter.updMobilizedBody(MobilizedBodyIndex(NUM_BODIES));
    Constraint::Rod constraint(first, last, 3.0);
    createState(system, state);
    class DisableHandler : public ScheduledEventHandler {
    public:
        DisableHandler(Constraint& constraint) : constraint(constraint) {
        }
        void handleEvent(State& state, Real accuracy, const Vector& yWeights, const Vector& ooConstraintTols, Stage& lowestModified, bool& shouldTerminate) const {
            constraint.disable(state);
            lowestModified = Stage::Instance;
        }
        Real getNextEventTime(const State&, bool includeCurrentTime) const {
            return 4.9;
        }
        Constraint& constraint;
    };
    system.updDefaultSubsystem().addEventHandler(new DisableHandler(constraint));
    RungeKuttaMersonIntegrator integ(system);
    integ.setConstraintTolerance(ConstraintTol);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    for (int i = 0; i < 10; i++) {
        ts.stepTo(i+1);
        if (i < 4) {
            CONSTRAINT_TEST(ts.getState().getNQErr(), 1);
            CONSTRAINT_TEST(ts.getState().getNUErr(), 1);
            MACHINE_TEST(ts.getState().getNUDotErr(), 1);
            Vec3 r1 = first.getBodyOriginLocation(ts.getState());
            Vec3 r2 = last.getBodyOriginLocation(ts.getState());
            Vec3 dr = r1-r2;
            // Verify that the constraint is enforced while enabled.
            CONSTRAINT_TEST(dr.norm(), 3.0);
        }
        else {
            CONSTRAINT_TEST(ts.getState().getNQErr(), 0);
            CONSTRAINT_TEST(ts.getState().getNUErr(), 0);
            MACHINE_TEST(ts.getState().getNUDotErr(), 0);
            Vec3 r1 = first.getBodyOriginLocation(ts.getState());
            Vec3 r2 = last.getBodyOriginLocation(ts.getState());
            Vec3 dr = r1-r2;
            // Verify that the constraint is *not* being enforced any more.
            SimTK_TEST_NOTEQ_TOL(dr.norm(), Real(3), ConstraintTol);
        }
    }
    delete &system;
}

void testConstraintForces() {
    // Weld one body to ground, push on it, verify that it reacts to match the load.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    MobilizedBody::Weld welded(matter.Ground(), Transform(),
                               body, Transform());

    MobilizedBody::Free loose(matter.Ground(), Transform(),
                               body, Transform());
    Constraint::Weld glue(matter.Ground(), Transform(),
                          loose, Transform());

    // Apply forces to the body welded straight to ground.
    Force::ConstantForce(forces, welded, Vec3(0,0,0), Vec3(1,2,3));
    Force::ConstantTorque(forces, welded, Vec3(20,30,40));

    // Apply the same forces to the "free" body which is welded by constraint.
    Force::ConstantForce(forces, loose, Vec3(0,0,0), Vec3(1,2,3));
    Force::ConstantTorque(forces, loose, Vec3(20,30,40));

    State state = system.realizeTopology();
    system.realize(state, Stage::Acceleration);

    Vector_<SpatialVec> mobilizerReactions;
    matter.calcMobilizerReactionForces(state, mobilizerReactions);

    //NOT IMPLEMENTED YET:
    //cout << "Weld constraint reaction on Ground: " << glue.getWeldReactionOnBody1(state) << endl;
    //cout << "Weld constraint reaction on Body: " << glue.getWeldReactionOnBody2(state) << endl;


    // Note that constraint forces have opposite sign to applied forces, because
    // we calculate the multiplier lambda from M udot + ~G lambda = f_applied. We'll negate
    // the calculated multipliers to turn these into applied forces.
    const Vector mults = -state.getMultipliers();
    Vector_<SpatialVec> bodyForces;
    Vector mobilityForces;
    matter.calcConstraintForcesFromMultipliers(state, mults,
        bodyForces, mobilityForces);

    MACHINE_TEST(bodyForces[loose.getMobilizedBodyIndex()], mobilizerReactions[welded.getMobilizedBodyIndex()]);
}

int main() {
    SimTK_START_TEST("TestConstraints");
        SimTK_SUBTEST(testBallConstraint);
        SimTK_SUBTEST(testConstantAngleConstraint);
        SimTK_SUBTEST(testConstantOrientationConstraint);
        SimTK_SUBTEST(testConstantSpeedConstraint);
        SimTK_SUBTEST(testNoSlip1DConstraint);
        SimTK_SUBTEST(testPointInPlaneConstraint);
        SimTK_SUBTEST(testPointOnLineConstraint);
        SimTK_SUBTEST(testRodConstraint);
        SimTK_SUBTEST(testWeldConstraint);
        SimTK_SUBTEST(testWeldConstraintWithPreAssembly);
        SimTK_SUBTEST(testDisablingConstraints);
        SimTK_SUBTEST(testConstraintForces);
    SimTK_END_TEST();
}
