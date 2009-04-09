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

using namespace SimTK;
using namespace std;

const Real TOL = 1e-10;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2, Real tol) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, Real tol) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2, Real tol) {
    assertEqual(val1[0], val2[0], tol);
    assertEqual(val1[1], val2[1], tol);
}

template<>
void assertEqual(Transform val1, Transform val2, Real tol) {
    assertEqual(val1.p(), val2.p(), tol);
    assertEqual(val1.R().convertRotationToBodyFixedXYZ(), 
                val2.R().convertRotationToBodyFixedXYZ(), tol);
}


template<class T>
void assertEqual(Vector_<T> val1, Vector_<T> val2, Real tol) {
    for (int i=0; i < val1.size(); ++i)
        assertEqual(val1[i], val2[i], tol);
}

template <class T>
void assertEqual(T val1, T val2) {
    assertEqual(val1, val2, TOL);
}

// Require this to be an identity matrix to within n*eps where
// n is the size of the matrix.
void assertIsIdentity(const Matrix& m) {
	const int minD = std::min(m.nrow(), m.ncol());
	const Real Slop =  minD * SignificantReal;

	const RowVector colSums = m.sum();
	for (int j=0; j < minD; ++j) {
		ASSERT( std::fabs(m(j,j)-1) <= Slop );
		ASSERT( std::fabs(colSums[j]-1) <= Slop );
	}

}

class MyForceImpl : public Force::Custom::Implementation {
public:
	MyForceImpl() {}
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, 
				   Vector& mobilityForces) const 
	{
		ASSERT( f.size() == 0 || f.size() == mobilityForces.size() );
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
	Random::Uniform rand(-1,1);
	for (int i=0; i < nq; ++i)
		state.updQ()[i] = rand.getValue();
	for (int i=0; i < nu; ++i)
		state.updU()[i] = rand.getValue();

	Vector randVec(nu);
	for (int i=0; i<nu; ++i) randVec[i] = 100*rand.getValue();

	Vector result1, result2;

	// result1 = M*v
	system.realize(state, Stage::Position);
	matter.calcMV(state, randVec, result1);
	ASSERT(result1.size() == nu);

	// result2 = M^-1 * result1 == M^-1 * M * v == v
	system.realize(state, Stage::Dynamics);
	Vector_<SpatialVec> A_GB; // TODO: get rid of this
	matter.calcMInverseV(state, result1, result2, A_GB);
	ASSERT(result2.size() == nu);

	//cout << "|v - M^-1*M*v|=" << (randVec - result2).norm() << endl;

	assertEqual(result2, randVec, Slop);

	Matrix M(nu,nu), MInv(nu,nu);

	Vector v(nu, Real(0));
	for (int j=0; j < nu; ++j) {
		v[j] = 1;
		matter.calcMV(state, v, M(j));
		matter.calcMInverseV(state, v, MInv(j), A_GB);
		v[j] = 0;
	}

	Matrix eye = M*MInv;

	assertIsIdentity(eye);
	assertIsIdentity(MInv*M);

	frcp->setForce(randVec);
	//cout << "f=" << randVec << endl;
	system.realize(state, Stage::Acceleration);
	Vector accel = state.getUDot();
	//cout << "v!=0, accel=" << accel << endl;

	matter.calcMInverseV(state, randVec, result1, A_GB);
	//cout << "With velocities, |a - M^-1*f|=" << (accel-result1).norm() << endl;

	ASSERT((accel-result1).norm() > SignificantReal); // because of velocities

	// With no velocities M^-1*f should match calculated acceleration.
	state.updU() = 0;
	system.realize(state, Stage::Acceleration);
	accel = state.getUDot();
	//cout << "v=0, accel=" << accel << endl;

	//cout << "With v=0, |a - M^-1*f|=" << (accel-result1).norm() << endl;

	ASSERT((accel-result1).norm() <= Slop); // because no velocities

	// And then M*a should = f.
	matter.calcMV(state, accel, result2);
	//cout << "v=0, M*accel=" << result2 << endl;
	//cout << "v=0, |M*accel-f|=" << (result2-randVec).norm() << endl;

}

void testTreeSystem() {
	MultibodySystem			mbs;
    SimbodyMatterSubsystem  pend(mbs);
    GeneralForceSubsystem   forces(mbs);
	MyForceImpl* frcp = new MyForceImpl();
	Force::Custom(forces, frcp);

	Random::Uniform rand(-1,1);

	const Real randomAngle1 = (Pi/2)*(rand.getValue());
	const Real randomAngle2 = (Pi/2)*(rand.getValue());
	Vector_<Vec3> randomVecs(10);
	for (int i=0; i<10; ++i) 
		randomVecs[i] = 
			Vec3( rand.getValue(), rand.getValue(), rand.getValue() );

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



int main() {
    try {
        cout << "*** TEST MASS MATRIX CALCS  ***\n\n"; 
		testTreeSystem();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

