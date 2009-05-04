/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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


// Test the functioning of Simbody operators which involve the mass matrix.
// The O(N) operators like calcMV() and calcMInverseV() are supposed to behave
// *as though* they used the mass matrix, without actually forming it.

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

using namespace SimTK;
using namespace std;


class MyForceImpl : public Force::Custom::Implementation {
public:
	MyForceImpl() {}
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, 
				   Vector& mobilityForces) const 
	{
		SimTK_TEST( f.size() == 0 || f.size() == mobilityForces.size() );
		if (f.size() > 0)
			mobilityForces += f;
    }
    Real calcPotentialEnergy(const State&) const {return 0;}

	void setForce(const Vector& frc) {
		f = frc;
	}
private:
	Vector f;
};

void testSystem(const MultibodySystem& system, MyForceImpl* frcp) {
	const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

	State state = system.realizeTopology();
	const int nq = state.getNQ();
	const int nu = state.getNU();

	// Attainable accuracy drops with problem size.
	const Real Slop = nu*SignificantReal;

	system.realizeModel(state);
	// Randomize state.
    state.updQ() = Test::randVector(nq);
    state.updU() = Test::randVector(nu);

    Vector randVec = 100*Test::randVector(nu);
	Vector result1, result2;

	// result1 = M*v
	system.realize(state, Stage::Position);
	matter.calcMV(state, randVec, result1);
	SimTK_TEST_EQ(result1.size(), nu);

	// result2 = M^-1 * result1 == M^-1 * M * v == v
	system.realize(state, Stage::Dynamics);
	matter.calcMInverseV(state, result1, result2);
	SimTK_TEST_EQ(result2.size(), nu);

    SimTK_TEST_EQ_TOL(result2, randVec, Slop);

	Matrix M(nu,nu), MInv(nu,nu);

	Vector v(nu, Real(0));
	for (int j=0; j < nu; ++j) {
		v[j] = 1;
		matter.calcMV(state, v, M(j));
		matter.calcMInverseV(state, v, MInv(j));
		v[j] = 0;
	}

    Matrix identity(nu,nu); identity=1;
    SimTK_TEST_EQ_SIZE(M*MInv, identity, nu);
    SimTK_TEST_EQ_SIZE(MInv*M, identity, nu);

	//assertIsIdentity(eye);
	//assertIsIdentity(MInv*M);

	frcp->setForce(randVec);
	//cout << "f=" << randVec << endl;
	system.realize(state, Stage::Acceleration);
	Vector accel = state.getUDot();
	//cout << "v!=0, accel=" << accel << endl;

	matter.calcMInverseV(state, randVec, result1);
	//cout << "With velocities, |a - M^-1*f|=" << (accel-result1).norm() << endl;

    SimTK_TEST_NOTEQ(accel, result1); // because of the velocities
	//SimTK_TEST((accel-result1).norm() > SignificantReal); // because of velocities

	// With no velocities M^-1*f should match calculated acceleration.
	state.updU() = 0;
	system.realize(state, Stage::Acceleration);
	accel = state.getUDot();
	//cout << "v=0, accel=" << accel << endl;

	//cout << "With v=0, |a - M^-1*f|=" << (accel-result1).norm() << endl;

    SimTK_TEST_EQ(accel, result1); // because no velocities

	// And then M*a should = f.
	matter.calcMV(state, accel, result2);
	//cout << "v=0, M*accel=" << result2 << endl;
	//cout << "v=0, |M*accel-f|=" << (result2-randVec).norm() << endl;


    // Test forward and inverse dynamics operators.
    // Apply random forces and a random prescribed acceleration to
    // get back the residual generalized forces. Then applying those
    // should result in zero residual, and applying them. 

	// Randomize state.
    state.updQ() = Test::randVector(nq);
    state.updU() = Test::randVector(nu);


    // Inverse dynamics should require realization only to Velocity stage.
    system.realize(state, Stage::Velocity);

    // Randomize body forces.
    Vector_<SpatialVec> bodyForces(matter.getNumBodies());
    for (int i=0; i < matter.getNumBodies(); ++i)
        bodyForces[i] = Test::randSpatialVec();

    // Random mobility forces and known udots.
    Vector mobilityForces = Test::randVector(matter.getNumMobilities());
    Vector knownUdots = Test::randVector(matter.getNumMobilities());

    // Check self consistency: compute residual, apply it, should be no remaining residual.
    Vector residualForces, shouldBeZeroResidualForces;
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, bodyForces, knownUdots, residualForces);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces+residualForces, bodyForces, knownUdots, shouldBeZeroResidualForces);

    SimTK_TEST(shouldBeZeroResidualForces.norm() <= Slop);

    // Now apply these forces in forward dynamics and see if we get the desired
    // acceleration. State must be realized to Dynamics stage.
    system.realize(state, Stage::Dynamics);
    Vector udots;
    Vector_<SpatialVec> bodyAccels;
    matter.calcAccelerationIgnoringConstraints(state, 
        mobilityForces+residualForces, bodyForces, udots, bodyAccels);

    SimTK_TEST_EQ_TOL(udots, knownUdots, Slop);

    // Verify that leaving out arguments makes them act like zeroes.
    Vector residualForces1, residualForces2;
    matter.calcResidualForceIgnoringConstraints(state,
        0*mobilityForces, 0*bodyForces, 0*knownUdots, residualForces1);
    // no, the residual is not zero here because of the angular velocities
    matter.calcResidualForceIgnoringConstraints(state,
        Vector(), Vector_<SpatialVec>(), Vector(), residualForces2);

    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);

    // Same, but leave out combinations of arguments.
    matter.calcResidualForceIgnoringConstraints(state,
        0*mobilityForces, bodyForces, knownUdots, residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
        Vector(), bodyForces, knownUdots, residualForces2);
    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, 0*bodyForces, knownUdots, residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, Vector_<SpatialVec>(), knownUdots, residualForces2);
    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, bodyForces, 0*knownUdots, residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, bodyForces, Vector(), residualForces2);
    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);
    matter.calcResidualForceIgnoringConstraints(state,
        0*mobilityForces, bodyForces, 0*knownUdots, residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
        Vector(), bodyForces, Vector(), residualForces2);
    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, 0*bodyForces, 0*knownUdots, residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
        mobilityForces, Vector_<SpatialVec>(), Vector(), residualForces2);
    SimTK_TEST_EQ_TOL(residualForces2, residualForces1, Slop);

    // Check that we object to wrong-length arguments.
    SimTK_TEST_MUST_THROW(matter.calcResidualForceIgnoringConstraints(state,
        Vector(3,Zero), bodyForces, knownUdots, residualForces2));
    SimTK_TEST_MUST_THROW(matter.calcResidualForceIgnoringConstraints(state,
         mobilityForces, Vector_<SpatialVec>(5), knownUdots, residualForces2));
    SimTK_TEST_MUST_THROW(matter.calcResidualForceIgnoringConstraints(state,
         mobilityForces, bodyForces, Vector(2), residualForces2));

}

void testTreeSystem() {
	MultibodySystem			mbs;
    SimbodyMatterSubsystem  pend(mbs);
    GeneralForceSubsystem   forces(mbs);
	MyForceImpl* frcp = new MyForceImpl();
	Force::Custom(forces, frcp);

    const Real randomAngle1 = (Pi/2)*Test::randReal();
    const Real randomAngle2 = (Pi/2)*Test::randReal();
	Vector_<Vec3> randomVecs(10);
	for (int i=0; i<10; ++i) 
        randomVecs[i] = Test::randVec3();

	const Real mass = 2.3;
	const Vec3 com = randomVecs[5];
	const Inertia inertia = Inertia(3,4,5,.01,-.02,.04).shiftFromMassCenter(com, mass);
    Body::Rigid pendulumBody = Body::Rigid(
		MassProperties(mass, com, inertia));


    MobilizedBody::Ball
        pendBody1(	pend.Ground(),
						Transform(Rotation(randomAngle1, randomVecs[0]),
							      randomVecs[1]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[2]),
								  randomVecs[3]));
    MobilizedBody::Weld
        pendBody2(	pendBody1,
						Transform(Rotation(randomAngle1, randomVecs[4]),
							      randomVecs[5]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[6]),
								  randomVecs[7]));

    MobilizedBody::Pin
        pendBody3(	pendBody2,
						Transform(Rotation(randomAngle1, randomVecs[8]),
							      randomVecs[9]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[8]),
								  randomVecs[7]));
    MobilizedBody::Screw
        pendBody4(	pendBody3,
						Transform(Rotation(randomAngle2, randomVecs[6]),
							      randomVecs[5]),
					pendulumBody,
						Transform(Rotation(randomAngle1, randomVecs[4]),
								  randomVecs[3]),
					3); // pitch
    MobilizedBody::Translation
        pendBody5(	pendBody4,
						Transform(Rotation(randomAngle2, randomVecs[2]),
							      randomVecs[1]),
					pendulumBody,
						Transform(Rotation(randomAngle1, randomVecs[0]),
								  randomVecs[1]));

	// Now add some side branches.
    MobilizedBody::BendStretch
        pendBody1a(	pendBody1,
						Transform(Rotation(randomAngle2, randomVecs[2]),
							      randomVecs[3]),
					pendulumBody,
						Transform(Rotation(randomAngle1, randomVecs[4]),
								  randomVecs[5]));

    MobilizedBody::Slider
        pendBody2a(	pendBody2,
						Transform(Rotation(randomAngle1, randomVecs[6]),
							      randomVecs[7]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[8]),
								  randomVecs[9]));

    MobilizedBody::Universal
        pendBody2b(	pendBody2a,
						Transform(Rotation(randomAngle1, randomVecs[8]),
							      randomVecs[7]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[6]),
								  randomVecs[5]));

    MobilizedBody::Planar
        pendBody4a(	pendBody4,
						Transform(Rotation(randomAngle1, randomVecs[4]),
							      randomVecs[3]),
					pendulumBody,
						Transform(Rotation(randomAngle2, randomVecs[2]),
								  randomVecs[1]));

	testSystem(mbs, frcp);
}

void testCompositeInertia() {
	MultibodySystem			mbs;
    SimbodyMatterSubsystem  pend(mbs);

    Body::Rigid pointMass(MassProperties(3, Vec3(0), Inertia(0)));

    // Point mass at x=1.5 rotating about (0,0,0).
    MobilizedBody::Pin
        body1( pend.Ground(), Transform(), 
               pointMass, Vec3(1.5,0,0));

    // A second body 2 units further along x, rotating about the
    // first point mass origin.
    MobilizedBody::Pin
        body2( body1, Transform(), 
               pointMass, Vec3(2,0,0));

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Position);

    Vector_<SpatialMat> R(pend.getNumBodies());
    pend.calcCompositeBodyInertias(state, R);

    // Calculate expected inertias about the joint axes.
    Real expInertia2 = body2.getBodyMassProperties(state).getMass()*square(2);
    Real expInertia1 = body1.getBodyMassProperties(state).getMass()*square(1.5)
                           + body2.getBodyMassProperties(state).getMass()*square(3.5);

    // Should be able to recover these inertias by projecting the composite
    // body inertias onto the joint axes using H matrices.
    const SpatialVec H1 = body1.getHCol(state, UIndex(0));
    const SpatialVec H2 = body2.getHCol(state, UIndex(0));
    SimTK_TEST_EQ(~H2*R[2]*H2, expInertia2);
    SimTK_TEST_EQ(~H1*R[1]*H1, expInertia1);
}

int main() {
    SimTK_START_TEST("TestMassMatrix");
        SimTK_SUBTEST(testCompositeInertia);
        SimTK_SUBTEST(testTreeSystem);
    SimTK_END_TEST();
}

