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

// This is an imitation of SD/FAST's sdrel2cart() subroutine. We
// are given a station point S fixed to a body B. S is given by
// the constant vector p_BS from B's origin to point S, expressed 
// in B's frame. Denote the position of S in the ground frame G
// p_GS = p_GB + p_BS_G, where p_BS_G=R_GB*p_BS is the vector p_BS
// reexpressed in G. The velocity of S in G is v_GS = d/dt p_GS, taken
// in G. So v_GS = v_GB + w_GB X p_BS_G = v_GB - p_BS_G % w_GB.
//
// We would like to obtain the partial velocity of S
// with respect to each of the generalized speeds u, taken in the
// Ground frame, that is, d v_GS / du. We have a method that can
// calculate J=d V_GB / du where V_GB=[w_GB;v_GB] is the spatial 
// velocity of B. So we need to calculate
//    d v_GS   d v_GS   d V_GB
//    ------ = ------ * ------
//      du     d V_GB     du
// 
//           = [ -px ; eye(3) ] * J
//  where px is the cross product matrix of p_BS_G.
// 
void sbrel2cart(const State& state,
              const SimbodyMatterSubsystem& matter,
              MobilizedBodyIndex            bodyIx,
              const Vec3&                   p_BS, // point in body frame
              Vector_<Vec3>&                dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of body origin.
    Vector_<SpatialVec> J(nu);
    J = SpatialVec(Vec3(0), Vec3(0)); // or J.setToZero();

    Vector u(nu); u = 0;
    Vector_<SpatialVec> Ju(nu); // d allV / d ui 
    for (int i=0; i < nu; ++i) {
        u[i] = 1;
        matter.calcSpatialKinematicsFromInternal(state,u,Ju);
        u[i] = 0;
        J[i] = Ju[bodyIx]; // pick out the body of interest
    }

    Row<2,Mat33> dvdV( -crossMat(p_BS_G), Mat33(1) );
    dvdu.resize(nu);
    for (int i=0; i < nu; ++i)
        dvdu[i] = dvdV * J[i]; // or J[i][0] % p_BS_G + J[i][1]
}

void testRel2Cart() {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    // Pendulum of length 1, initially along -x like this:
    //     B ----------- * O
    //    -1,0,0          0,0,0
    // At q=0, partial(B)/partial(u) = 0,-1,0.
    // At q=pi/2, partial(B)/partial(u) = 1,0,0.
    // Then try this with station S=(0,1,0)_B:
    //      S
    //      |
    //      |
    //      B ------ * O
    // Now |OS| = sqrt(2). At q=Pi/4, S will be horizontal
    // so partial(S)/partial(u) = (0, -sqrt(2)/2, 0).
    // 
    MobilizedBody::Pin pinBody
       (matter.Ground(),                       Transform(),
        MassProperties(1,Vec3(0),Inertia(1)),  Vec3(1,0,0));
    State s = system.realizeTopology();

    Vector_<Vec3> dvdu;

    pinBody.setQ(s, 0);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(0,-1,0));

    pinBody.setQ(s, Pi/2);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(1,0,0));

    pinBody.setQ(s, Pi/4);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0,1,0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(0,-Sqrt2,0));
}
              

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

    Matrix MInvCalc(M);
    MInvCalc.invertInPlace();
    SimTK_TEST_EQ_SIZE(MInv, MInvCalc, nu);

    Matrix identity(nu,nu); identity=1;
    SimTK_TEST_EQ_SIZE(M*MInv, identity, nu);
    SimTK_TEST_EQ_SIZE(MInv*M, identity, nu);

    // Compare above-calculated values with values returned by the
    // calcM() and calcMInv() methods.
    Matrix MM, MMInv;
    matter.calcM(state,MM); matter.calcMInv(state,MMInv);
    SimTK_TEST_EQ_SIZE(MM, M, nu);
    SimTK_TEST_EQ_SIZE(MMInv, MInv, nu);

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

    // This should force realization of the composite body inertias.
    SpatialMat cbi = pend.getCompositeBodyInertia(state, body1);

    body2.setAngle(state, Pi/4);
    // This is not allowed until Position stage.
    SimTK_TEST_MUST_THROW(pend.getCompositeBodyInertia(state, body1));
    mbs.realize(state, Stage::Position);
    // Now it should be OK.
    cbi = pend.getCompositeBodyInertia(state, body1);

    mbs.realize(state, Stage::Acceleration);
    cout << "udots=" << state.getUDot() << endl;

    body1.setRate(state, 27);
    mbs.realize(state, Stage::Acceleration);
    cout << "udots=" << state.getUDot() << endl;
}

int main() {
    SimTK_START_TEST("TestMassMatrix");
        SimTK_SUBTEST(testRel2Cart);
        SimTK_SUBTEST(testCompositeInertia);
        SimTK_SUBTEST(testTreeSystem);
    SimTK_END_TEST();
}

