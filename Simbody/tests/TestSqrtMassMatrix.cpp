/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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


// Test the functioning of Simbody operators which involve the mass matrix,
// and other system matrices like the Jacobian (partial velocity matrix) that
// maps between generalized and spatial coordinates.
// The O(N) operators like multiplyByM() and multiplyByMInv() are supposed to 
// behave *as though* they used the mass matrix, without actually forming it.

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <iostream>

using namespace SimTK;
using std::cout; using std::endl;

// This will apply a constant set of mobility forces that can be set
// externally.
class MyForceImpl : public Force::Custom::Implementation {
public:
    MyForceImpl() {}
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, 
                   Vector& mobilityForces) const override 
    {
        SimTK_TEST( f.size() == 0 || f.size() == mobilityForces.size() );
        SimTK_TEST( F.size() == 0 || F.size() == bodyForces.size() );
        if (f.size())
            mobilityForces += f;
        if (F.size())
            bodyForces += F;
    }
    Real calcPotentialEnergy(const State&) const override {return 0;}

    void setMobilityForces(const Vector& mobFrc) {
        f = mobFrc;
    }
    void setBodyForces(const Vector_<SpatialVec>& bodFrc) {
        F = bodFrc;
    }
private:
    Vector              f;
    Vector_<SpatialVec> F;
};

// Compare two representations of the same matrix: one as an mXn matrix
// of SpatialVecs, the other as a 6mXn matrix of scalars. Note that this
// will also work if the first actual parameter is a Vector_<SpatialVec> 
// or RowVector_<SpatialVec> since those have implicit conversions to mX1 
// or 1Xn Matrix_<SpatialVec>, resp.
static void compareElementwise(const Matrix_<SpatialVec>& J,
                               const Matrix&              Jf) 
{
    const int m = J.nrow(), n = J.ncol();
    SimTK_TEST(Jf.nrow()==6*m && Jf.ncol()==n);

    for (int b=0; b<m; ++b) {
        const int r = 6*b; // row start for Jf
        for (int i=0; i<6; ++i)
            for (int j=0; j<n; ++j) {
                SimTK_TEST_EQ(J (b,  j)[i/3][i%3], 
                              Jf(r+i,j));
            }
    }
}

// Same thing but for comparing matrices where one has Vec3 elements.
static void compareElementwise(const Matrix_<Vec3>& JS,
                               const Matrix&        JSf) 
{
    const int m = JS.nrow(), n = JS.ncol();
    SimTK_TEST(JSf.nrow()==3*m && JSf.ncol()==n);

    for (int b=0; b<m; ++b) {
        const int r = 3*b; // row start for JSf
        for (int i=0; i<3; ++i)
            for (int j=0; j<n; ++j) {
                SimTK_TEST_EQ(JS (b,  j)[i], 
                              JSf(r+i,j));
            }
    }
}
              
void makeSystem(bool constrained, MultibodySystem& mbs, MyForceImpl*& frcp) {    
    SimbodyMatterSubsystem  pend(mbs);
    GeneralForceSubsystem   forces(mbs);
    frcp = new MyForceImpl();
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
        pendBody1(  pend.Ground(),
                        Transform(Rotation(randomAngle1, randomVecs[0]),
                                  randomVecs[1]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[2]),
                                  randomVecs[3]));
    MobilizedBody::Weld
        pendBody2(  pendBody1,
                        Transform(Rotation(randomAngle1, randomVecs[4]),
                                  randomVecs[5]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[6]),
                                  randomVecs[7]));

    MobilizedBody::Pin
        pendBody3(  pendBody2,
                        Transform(Rotation(randomAngle1, randomVecs[8]),
                                  randomVecs[9]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[8]),
                                  randomVecs[7]));
    MobilizedBody::Screw
        pendBody4(  pendBody3,
                        Transform(Rotation(randomAngle2, randomVecs[6]),
                                  randomVecs[5]),
                    pendulumBody,
                        Transform(Rotation(randomAngle1, randomVecs[4]),
                                  randomVecs[3]),
                    3); // pitch
    MobilizedBody::Translation
        pendBody5(  pendBody4,
                        Transform(Rotation(randomAngle2, randomVecs[2]),
                                  randomVecs[1]),
                    pendulumBody,
                        Transform(Rotation(randomAngle1, randomVecs[0]),
                                  randomVecs[1]));

    // Now add some side branches.
    MobilizedBody::BendStretch
        pendBody1a( pendBody1,
                        Transform(Rotation(randomAngle2, randomVecs[2]),
                                  randomVecs[3]),
                    pendulumBody,
                        Transform(Rotation(randomAngle1, randomVecs[4]),
                                  randomVecs[5]));

    MobilizedBody::Slider
        pendBody2a( pendBody2,
                        Transform(Rotation(randomAngle1, randomVecs[6]),
                                  randomVecs[7]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[8]),
                                  randomVecs[9]));

    MobilizedBody::Universal
        pendBody2b( pendBody2a,
                        Transform(Rotation(randomAngle1, randomVecs[8]),
                                  randomVecs[7]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[6]),
                                  randomVecs[5]));
    MobilizedBody::Slider
        pendBody2x( pendBody2b,
                        Transform(Rotation(randomAngle1, randomVecs[6]),
                                  randomVecs[7]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[8]),
                                  randomVecs[9]));

    MobilizedBody::Planar
        pendBody4a( pendBody4,
                        Transform(Rotation(randomAngle1, randomVecs[4]),
                                  randomVecs[3]),
                    pendulumBody,
                        Transform(Rotation(randomAngle2, randomVecs[2]),
                                  randomVecs[1]));


    if (constrained) { // probably can't be satisfied, but doesn't matter
        Constraint::Rod(pendBody4, pendBody2b, 1.); // holonomic
        Constraint::ConstantSpeed
            (pendBody2a, MobilizerUIndex(0), -3.); // nonholo
        Constraint::ConstantAcceleration
            (pendBody5, MobilizerUIndex(2), 0.01); // acc only
        Constraint::Weld(pendBody4a, Test::randTransform(),
                         pendBody4, Test::randTransform());
    }
}


void testUnconstrainedSystem() {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    // Attainable accuracy drops with problem size.
    const Real Slop = nu*SignificantReal;

    system.realizeModel(state);
    // Randomize state.
    state.updQ() = Test::randVector(nq);
    state.updU() = Test::randVector(nu);

    Vector randVec = 100*Test::randVector(nu);

    system.realize(state, Stage::Position);


    // Compute SqrtMInv - calcSqrtMInv not yet implemented
    Matrix MInv, SqrtMInv; 
    SqrtMInv.resize(nu,nu);
    Vector v(nu);
    v.setToZero();
    for (int i=0; i < nu; ++i) {
      state.invalidateAllCacheAtOrAbove(Stage::Position);
      system.realize(state, Stage::Position);
      v[i] = 1;
      matter.multiplyBySqrtMInv(state, v, SqrtMInv(i));
      v[i] = 0;
    }

    matter.calcMInv(state, MInv); 
    SimTK_TEST_EQ_SIZE((SqrtMInv * ~SqrtMInv), MInv, nu);
    SimTK_TEST_EQ_TOL((SqrtMInv * ~SqrtMInv), MInv, Slop);

    // multiplyBySqrtMInvTranspose
    // TODO
    // calcSqrtMInv
    // TODO
    // Check Maxwell-Boltzmann distribution
    // TODO
}



int main() {
    SimTK_START_TEST("TestSqrtMassMatrix");
        SimTK_SUBTEST(testUnconstrainedSystem);
    SimTK_END_TEST();
}

