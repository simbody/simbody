/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-11 Stanford University and the Authors.        *
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
                   Vector& mobilityForces) const 
    {
        SimTK_TEST( f.size() == 0 || f.size() == mobilityForces.size() );
        SimTK_TEST( F.size() == 0 || F.size() == bodyForces.size() );
        if (f.size())
            mobilityForces += f;
        if (F.size())
            bodyForces += F;
    }
    Real calcPotentialEnergy(const State&) const {return 0;}

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

// This is an imitation of SD/FAST's sdrel2cart() subroutine. We
// are given a station point S fixed to a body B. S is given by
// the constant vector p_BS from B's origin Bo to point S, expressed 
// in B's frame. Denote the position of S in the ground frame G
// p_GS = p_GB + p_BS_G, where p_BS_G=R_GB*p_BS is the vector p_BS
// reexpressed in G. The velocity of S in G is v_GS = d/dt p_GS, taken
// in G. So v_GS = v_GB + w_GB X p_BS_G = v_GB - p_BS_G % w_GB.
//
// We would like to obtain the partial velocity of S with respect to each of 
// the generalized speeds u, taken in the Ground frame, that is, 
// JS=d v_GS / du. (JS is a 3xnu matrix, or a single row of Vec3s.) We have 
// a method that can calculate J=d V_GB / du where V_GB=[w_GB;v_GB] is the 
// spatial velocity of B at its origin. So we need to calculate
//        d v_GS   d v_GS   d V_GB
//   JS = ------ = ------ * ------ = 
//          du     d V_GB     du
// 
//               = [ -px | eye(3) ] * J
// where px is the cross product matrix of p_BS_G and eye(3) is a 3x3
// identity matrix.
// 
// This function should produce the same result as the SimbodyMatterSubsystem
// method calcStationJacobian().
void sbrel2cart(const State& state,
                const SimbodyMatterSubsystem& matter,
                MobilizedBodyIndex            bodyIx,
                const Vec3&                   p_BS, // point in body frame
                RowVector_<Vec3>&             dvdu) // v is dS/dt in G
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
        matter.multiplyBySystemJacobian(state,u,Ju);
        u[i] = 0;
        J[i] = Ju[bodyIx]; // pick out the body of interest
    }

    Row<2,Mat33> dvdV( -crossMat(p_BS_G), Mat33(1) );
    dvdu.resize(nu);
    for (int i=0; i < nu; ++i)
        dvdu[i] = dvdV * J[i]; // or J[i][0] % p_BS_G + J[i][1]
}

// Another way to calculate exactly what sbrel2cart() does -- can we do it
// faster using f=J^T*F rather than V=J*u?
void sbrel2cart2(const State& state,
                 const SimbodyMatterSubsystem& matter,
                 MobilizedBodyIndex            bodyIx,
                 const Vec3&                   p_BS, // point in body frame
                 RowVector_<Vec3>&             dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of body origin.
    // (This is one row of J.)
    Matrix Jt(nu, 6); // a column of Jt but with scalar elements

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx]; // the only one we'll change
    for (int which=0; which < 2; ++which) { // moment, force   
        for (int i=0; i < 3; ++i) {
            Fb[which][i] = 1;
            VectorView col = Jt(3*which + i);
            matter.multiplyBySystemJacobianTranspose(state,F,col);
            Fb[which][i] = 0;
        }
    }

    Row<2,Mat33> dvdV( -crossMat(p_BS_G), Mat33(1) );
    dvdu.resize(nu);
    for (int i=0; i < nu; ++i) {
        const RowVectorView r = Jt[i]; 
        SpatialVec V(Vec3::getAs(&r[0]), Vec3::getAs(&r[3]));
        dvdu[i] = dvdV * V; // or J[i][0] % p_BS_G + J[i][1]
    }
}

// This is a further refinement that still calculates exactly what sbrel2cart()
// does. But now try it without the intermediate storage for a row of J^T and
// using only 3 J*v multiplies since only the translational (station) Jacobian
// is wanted.
void sbrel2cart3(const State& state,
                 const SimbodyMatterSubsystem& matter,
                 MobilizedBodyIndex            bodyIx,
                 const Vec3&                   p_BS, // point in body frame
                 RowVector_<Vec3>&             dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dvdu where v is linear velocity of p_BS.
    // (This is three rows of J.)
    dvdu.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F
        matter.multiplyBySystemJacobianTranspose(state,F,col);
        for (int r=0; r < nu; ++r) dvdu[r][i] = col[r]; 
        Fb[1][i] = 0;
    }
}

// Using the method of sbrel2cart3() but with 6 J*v multiplies, this gives
// the full 6xnu "Frame Jacobian" for one body for a specified frame on that
// body (only the origin, not the orientation, matters).
// This should produce the same result as the built-in calcFrameJacobian()
// method.
void sbrel2cart4(const State& state,
              const SimbodyMatterSubsystem& matter,
              MobilizedBodyIndex            bodyIx,
              const Vec3&                   p_BS, // point in body frame
              RowVector_<SpatialVec>&       dVdu) // V is [w,v] in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of p_BS.
    // (This is six rows of J.)
    dVdu.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    // Rotational part.
    for (int i=0; i < 3; ++i) {
        Fb[0][i] = 1;
        matter.multiplyBySystemJacobianTranspose(state,F,col);
        for (int r=0; r < nu; ++r) dVdu[r][0][i] = col[r]; 
        Fb[0][i] = 0;
    }
    // Translational part.
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F
        matter.multiplyBySystemJacobianTranspose(state,F,col);
        for (int r=0; r < nu; ++r) dVdu[r][1][i] = col[r]; 
        Fb[1][i] = 0;
    }
}

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
    // In all cases the partial angular velocity is (0,0,1).
    // 
    MobilizedBody::Pin pinBody
       (matter.Ground(),                       Transform(),
        MassProperties(1,Vec3(0),Inertia(1)),  Vec3(1,0,0));
    State s = system.realizeTopology();

    RowVector_<Vec3> dvdu, dvdu2, dvdu3;
    RowVector_<Vec3> JS;
    RowVector_<SpatialVec> dvdu4, JF;
    Matrix_<SpatialVec> J, Jn;
    Matrix JSf, JFf, Jf; // flat

    // We'll compute Jacobians for the body origin Bo in 2 configurations, 
    // q==0,pi/2 then for station S (0,1,0) at q==pi/4. Not
    // much of a test, I know, but at least we know the right answer.

    // q == 0; answer is JB == (0,-1,0)

    pinBody.setQ(s, 0);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(0,-1,0));

    sbrel2cart2(s, matter, pinBody, Vec3(0), dvdu2);
    SimTK_TEST_EQ(dvdu2[0], Vec3(0,-1,0));

    sbrel2cart3(s, matter, pinBody, Vec3(0), dvdu3);
    SimTK_TEST_EQ(dvdu3[0], Vec3(0,-1,0));

    matter.calcStationJacobian(s, pinBody, Vec3(0), JS);
    matter.calcStationJacobian(s, pinBody, Vec3(0), JSf);
    SimTK_TEST_EQ(JS, dvdu3); // == dvdu2 == dvdu
    compareElementwise(JS, JSf);

    sbrel2cart4(s, matter, pinBody, Vec3(0), dvdu4);
    SimTK_TEST_EQ(dvdu4[0][0], Vec3(0,0,1));
    SimTK_TEST_EQ(dvdu4[0][1], Vec3(0,-1,0));

    matter.calcFrameJacobian(s, pinBody, Vec3(0), JF);
    matter.calcFrameJacobian(s, pinBody, Vec3(0), JFf);
    SimTK_TEST_EQ(dvdu4, JF);
    compareElementwise(JF, JFf);

    // Calculate the whole system Jacobian at q==0 in two different 
    // representations and make sure they are the same.
    matter.calcSystemJacobian(s, J);
    matter.calcSystemJacobian(s, Jf);
    compareElementwise(J, Jf);

    // q == 90 degrees; answer is JB == (1,0,0)

    pinBody.setQ(s, Pi/2);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(1,0,0));

    sbrel2cart2(s, matter, pinBody, Vec3(0), dvdu2);
    SimTK_TEST_EQ(dvdu2[0], Vec3(1,0,0));

    sbrel2cart3(s, matter, pinBody, Vec3(0), dvdu3);
    SimTK_TEST_EQ(dvdu3[0], Vec3(1,0,0));

    matter.calcStationJacobian(s, pinBody, Vec3(0), JS);
    matter.calcStationJacobian(s, pinBody, Vec3(0), JSf);
    SimTK_TEST_EQ(JS, dvdu3); // == dvdu2 == dvdu
    compareElementwise(JS, JSf);

    sbrel2cart4(s, matter, pinBody, Vec3(0), dvdu4);
    SimTK_TEST_EQ(dvdu4[0][0], Vec3(0,0,1));
    SimTK_TEST_EQ(dvdu4[0][1], Vec3(1,0,0));

    matter.calcFrameJacobian(s, pinBody, Vec3(0), JF);
    matter.calcFrameJacobian(s, pinBody, Vec3(0), JFf);
    SimTK_TEST_EQ(dvdu4, JF);
    compareElementwise(JF, JFf);

    // Calculate the whole system Jacobian at q==pi/2 in two different 
    // representations and make sure they are the same.
    matter.calcSystemJacobian(s, J);
    matter.calcSystemJacobian(s, Jf);
    compareElementwise(J, Jf);

    // now station S, q == 45 degrees; answer is JS == (0,-sqrt(2),0)

    pinBody.setQ(s, Pi/4);
    system.realize(s, Stage::Position);
    sbrel2cart(s, matter, pinBody, Vec3(0,1,0), dvdu);
    SimTK_TEST_EQ(dvdu[0], Vec3(0,-Sqrt2,0));

    sbrel2cart2(s, matter, pinBody, Vec3(0,1,0), dvdu2);
    SimTK_TEST_EQ(dvdu2[0], Vec3(0,-Sqrt2,0));

    sbrel2cart3(s, matter, pinBody, Vec3(0,1,0), dvdu3);
    SimTK_TEST_EQ(dvdu3[0], Vec3(0,-Sqrt2,0));

    matter.calcStationJacobian(s, pinBody, Vec3(0,1,0), JS);
    matter.calcStationJacobian(s, pinBody, Vec3(0,1,0), JSf);
    SimTK_TEST_EQ(JS, dvdu3); // == dvdu2 == dvdu
    compareElementwise(JS, JSf);

    // Calculate station Jacobian JS by multiplication to test that the
    // multiplyByStationJacobian[Transpose] methods are working.
    Vec3 JSn; // 3xnu
    Vector u(1); u[0] = 1.;
    JSn = matter.multiplyByStationJacobian(s, pinBody, Vec3(0,1,0), u);
    SimTK_TEST_EQ(JSn, dvdu3[0]);

    Vec3 FS(0);
    Row3 JSnt;
    for (int i=0; i<3; ++i) {
        FS[i] = 1;
        matter.multiplyByStationJacobianTranspose(s,pinBody,Vec3(0,1,0),
            FS, u);
        FS[i] = 0;
        JSnt[i] = u[0];
    }
    SimTK_TEST_EQ(JSnt, ~dvdu3[0]);

    sbrel2cart4(s, matter, pinBody, Vec3(0,1,0), dvdu4);
    SimTK_TEST_EQ(dvdu4[0][0], Vec3(0,0,1));
    SimTK_TEST_EQ(dvdu4[0][1], Vec3(0,-Sqrt2,0));

    matter.calcFrameJacobian(s, pinBody, Vec3(0,1,0), JF);
    matter.calcFrameJacobian(s, pinBody, Vec3(0,1,0), JFf);
    SimTK_TEST_EQ(dvdu4, JF);
    compareElementwise(JF, JFf);

    // Calculate frame Jacobian JF by multiplication to test that the
    // multiplyByFrameJacobian[Transpose] methods are working.
    SpatialVec JFn; // 6xnu
    u[0] = 1.;
    JFn = matter.multiplyByFrameJacobian(s, pinBody, Vec3(0,1,0), u);
    SimTK_TEST_EQ(JFn, dvdu4[0]);

    SpatialVec FF(Vec3(0));
    SpatialRow JFnt;
    for (int i=0; i<6; ++i) {
        FF[i/3][i%3] = 1;
        matter.multiplyByFrameJacobianTranspose(s,pinBody,Vec3(0,1,0),
            FF, u);
        FF[i/3][i%3] = 0;
        JFnt[i/3][i%3] = u[0];
    }
    SimTK_TEST_EQ(JFnt, ~dvdu4[0]);

    // Calculate the whole system Jacobian at q==pi/4 in two different 
    // representations and make sure they are the same.
    matter.calcSystemJacobian(s, J);
    matter.calcSystemJacobian(s, Jf);
    compareElementwise(J, Jf);

    // Generate the whole system Jacobian one body at a time by multiplication.
    Jn.resize(matter.getNumBodies(), matter.getNumMobilities());
    u[0] = 1;
    for (MobodIndex i(0); i < matter.getNumBodies(); ++i)
        Jn[i] = matter.multiplyByFrameJacobian(s,i,Vec3(0),u);
    SimTK_TEST_EQ(Jn, J);
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
    Vector result1, result2;

    // result1 = M*v
    system.realize(state, Stage::Position);
    matter.multiplyByM(state, randVec, result1);
    SimTK_TEST_EQ(result1.size(), nu);

    // result2 = M^-1 * result1 == M^-1 * M * v == v
    system.realize(state, Stage::Dynamics);
    matter.multiplyByMInv(state, result1, result2);
    SimTK_TEST_EQ(result2.size(), nu);

    SimTK_TEST_EQ_TOL(result2, randVec, Slop);

    Matrix M(nu,nu), MInv(nu,nu);

    Vector v(nu, Real(0));
    for (int j=0; j < nu; ++j) {
        v[j] = 1;
        matter.multiplyByM(state, v, M(j));
        matter.multiplyByMInv(state, v, MInv(j));
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

    frcp->setMobilityForces(randVec);
    //cout << "f=" << randVec << endl;
    system.realize(state, Stage::Acceleration);
    Vector accel = state.getUDot();
    //cout << "v!=0, accel=" << accel << endl;

    matter.multiplyByMInv(state, randVec, result1);
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
    matter.multiplyByM(state, accel, result2);
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
    Vector_<SpatialVec> bodyForces(nb);
    for (int i=0; i < nb; ++i)
        bodyForces[i] = Test::randSpatialVec();

    // Random mobility forces and known udots.
    Vector mobilityForces = Test::randVector(nu);
    Vector knownUdots = Test::randVector(nu);

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

    // See if we get back the same body accelerations by feeding in 
    // these udots.
    Vector_<SpatialVec> A_GB, AC_GB;
    matter.calcBodyAccelerationFromUDot(state, udots, A_GB);
    SimTK_TEST_EQ_TOL(A_GB, bodyAccels, Slop);

    // Collect coriolis accelerations.
    AC_GB.resize(matter.getNumBodies());
    for (MobodIndex i(0); i<nb; ++i)
        AC_GB[i] = matter.getTotalCoriolisAcceleration(state, i);

    // Verify that either a zero-length or all-zero udot gives just
    // coriolis accelerations.
    matter.calcBodyAccelerationFromUDot(state, Vector(), A_GB);
    SimTK_TEST_EQ_TOL(A_GB, AC_GB, Slop);

    Vector allZeroUdot(matter.getNumMobilities(), Real(0));
    matter.calcBodyAccelerationFromUDot(state, allZeroUdot, A_GB);
    SimTK_TEST_EQ_TOL(A_GB, AC_GB, Slop);

    // Now let's test noncontiguous input and output vectors.
    Matrix MatUdot(3, nu); // use middle row
    MatUdot.setToNaN();
    MatUdot[1] = ~udots;
    Matrix_<SpatialRow> MatA_GB(3, nb); // use middle row
    MatA_GB.setToNaN();
    matter.calcBodyAccelerationFromUDot(state, ~MatUdot[1], ~MatA_GB[1]);
    SimTK_TEST_EQ_TOL(MatA_GB[1], ~bodyAccels, Slop);

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



void testConstrainedSystem() {
    MultibodySystem mbs;
    MyForceImpl* frcp;
    makeSystem(true, mbs, frcp);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Instance); // allocate multipliers, etc.

    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int m  = state.getNMultipliers();
    const int nb = matter.getNumBodies();

    // Attainable accuracy drops with problem size.
    const Real Slop = nu*SignificantReal;

    mbs.realizeModel(state);
    // Randomize state.
    state.updQ() = Test::randVector(nq);
    state.updU() = Test::randVector(nu);

    Vector randMobFrc = 100*Test::randVector(nu);
    Vector_<SpatialVec> randBodyFrc(nb);
    for (int i=0; i < nb; ++i)
        randBodyFrc[i] = Test::randSpatialVec();

    // Apply random mobility forces
    frcp->setMobilityForces(randMobFrc);

    mbs.realize(state); // calculate accelerations and multipliers
    Vector udot = state.getUDot();
    Vector lambda = state.getMultipliers();
    Vector residual;
    matter.calcResidualForce(state,randMobFrc,Vector_<SpatialVec>(),
                             udot, lambda, residual);

    // Residual should be zero since we accounted for everything.
    SimTK_TEST_EQ_TOL(residual, 0*randMobFrc, Slop);

    // Add in some body forces
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    frcp->setBodyForces(randBodyFrc);
    mbs.realize(state);
    udot = state.getUDot();
    lambda = state.getMultipliers();
    matter.calcResidualForce(state,randMobFrc,randBodyFrc,
                             udot, lambda, residual);
    SimTK_TEST_EQ_TOL(residual, 0*randMobFrc, Slop);

    // Try body forces only.
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    frcp->setMobilityForces(0*randMobFrc);
    mbs.realize(state);
    udot = state.getUDot();
    lambda = state.getMultipliers();
    matter.calcResidualForce(state,Vector(),randBodyFrc,
                             udot, lambda, residual);
    SimTK_TEST_EQ_TOL(residual, 0*randMobFrc, Slop);

    // Put vectors in noncontiguous storage.
    Matrix udotmat(3,nu); // rows are noncontig
    Matrix mobFrcMat(11,nu);
    Matrix lambdamat(5,m);
    Matrix_<SpatialRow> bodyFrcMat(3,nb);
    udotmat[2]    = ~udot;
    lambdamat[3]  = ~lambda;
    mobFrcMat[8] = ~randMobFrc;
    bodyFrcMat[2] = ~randBodyFrc;
    Matrix residmat(4,nu);

    // We last computed udot,lambda with no mobility forces. This time
    // will throw some in and then make sure the residual tries to cancel them.
    matter.calcResidualForce(state,~mobFrcMat[8],~bodyFrcMat[2],
        ~udotmat[2],~lambdamat[3],~residmat[2]);
    SimTK_TEST_EQ_TOL(residmat[2], -1*mobFrcMat[8], Slop);
}



void testCompositeInertia() {
    MultibodySystem         mbs;
    SimbodyMatterSubsystem  pend(mbs);

    Body::Rigid pointMass(MassProperties(3, Vec3(0), Inertia(0)));

    // Point mass at x=1.5 rotating about (0,0,0).
    MobilizedBody::Pin
        body1( pend.Ground(), Transform(), 
               pointMass, Vec3(1.5,0,0));
    const MobilizedBodyIndex body1x = body1.getMobilizedBodyIndex();

    // A second body 2 units further along x, rotating about the
    // first point mass origin.
    MobilizedBody::Pin
        body2( body1, Transform(), 
               pointMass, Vec3(2,0,0));
    const MobilizedBodyIndex body2x = body2.getMobilizedBodyIndex();

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Position);

    Array_<SpatialInertia, MobilizedBodyIndex> R(pend.getNumBodies());
    pend.calcCompositeBodyInertias(state, R);

    // Calculate expected inertias about the joint axes.
    Real expInertia2 = body2.getBodyMassProperties(state).getMass()*square(2);
    Real expInertia1 = body1.getBodyMassProperties(state).getMass()*square(1.5)
                           + body2.getBodyMassProperties(state).getMass()*square(3.5);

    // Should be able to recover these inertias by projecting the composite
    // body inertias onto the joint axes using H matrices.
    const SpatialVec H1 = body1.getHCol(state, MobilizerUIndex(0));
    const SpatialVec H2 = body2.getHCol(state, MobilizerUIndex(0));
    SimTK_TEST_EQ(~H2*(R[body2x]*H2), expInertia2);
    SimTK_TEST_EQ(~H1*(R[body1x]*H1), expInertia1);

    // This should force realization of the composite body inertias.
    SpatialInertia cbi = pend.getCompositeBodyInertia(state, body1);

    body2.setAngle(state, Pi/4);
    // This is not allowed until Position stage.
    SimTK_TEST_MUST_THROW(pend.getCompositeBodyInertia(state, body1));
    mbs.realize(state, Stage::Position);
    // Now it should be OK.
    cbi = pend.getCompositeBodyInertia(state, body1);

    mbs.realize(state, Stage::Acceleration);
    //cout << "udots=" << state.getUDot() << endl;

    body1.setRate(state, 27);
    mbs.realize(state, Stage::Acceleration);
    //cout << "udots=" << state.getUDot() << endl;
}

int main() {
    SimTK_START_TEST("TestMassMatrix");
        SimTK_SUBTEST(testRel2Cart);
        SimTK_SUBTEST(testCompositeInertia);
        SimTK_SUBTEST(testUnconstrainedSystem);
        SimTK_SUBTEST(testConstrainedSystem);
    SimTK_END_TEST();
}

