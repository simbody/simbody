/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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


// System is two free bodies A and B pinned together.
// The system is Ground (6dof) A (pin) B (6dof) Ground.
//   The forward system is Ground (6dof)-> A (pin)-> B.
//   The reverse system1 is A <-(pin) B <-(6dof) Ground.
//   The reverse system2 is A <-(rpin) B <-(6dof) Ground.
// Reverse system2's pin joint q and u should have the
// same meaning as the forward system's.
//
// Note: because HDot is zero for both forward & reverse
// pin joints, this test uses only H reversal and says
// nothing about HDot reversal!

const Real d = 1.5; // length (m)
const Real mass = 200000; // kg
// Mobilizer frame on A (used for both inboard and outboard)
const Transform X_AM(Rotation(Pi/3, Vec3(.1,-.3,.3)), Vec3(-4,-5,-1));
// Mobilizer frame on B (used for both inboard and outboard)
const Transform X_BM(Rotation(-Pi/10, Vec3(7,5,3)), Vec3(0,d,0));

void testPin() {
    MultibodySystem forward;
    SimbodyMatterSubsystem fwdMatter(forward);
    GeneralForceSubsystem fwdForces(forward);
    Force::UniformGravity(fwdForces, fwdMatter, Vec3(0, -1, 0));

    MultibodySystem reverse1;
    SimbodyMatterSubsystem rev1Matter(reverse1);
    GeneralForceSubsystem rev1Forces(reverse1);
    Force::UniformGravity(rev1Forces, rev1Matter, Vec3(0, -1, 0));

    MultibodySystem reverse2;
    SimbodyMatterSubsystem rev2Matter(reverse2);
    GeneralForceSubsystem rev2Forces(reverse2);
    Force::UniformGravity(rev2Forces, rev2Matter, Vec3(0, -1, 0));

    const Vec3 com(1,2,3);
    const UnitInertia centralGyration(1, 1.5, 2, .1, .2, .3);
    Body::Rigid body(MassProperties(mass, com, mass*centralGyration.shiftFromMassCenter(com, 1)));

    MobilizedBody::Free fwdA (fwdMatter.Ground(),  Vec3(0), body, X_AM);
    MobilizedBody::Free rev1B(rev1Matter.Ground(), Vec3(0), body, X_BM);
    MobilizedBody::Free rev2B(rev2Matter.Ground(), Vec3(0), body, X_BM);

    MobilizedBody::Pin fwdB (fwdA,  X_AM, body, X_BM);
    MobilizedBody::Pin rev1A(rev1B, X_BM, body, X_AM);
    MobilizedBody::Pin rev2A(rev2B, X_BM, body, X_AM, MobilizedBody::Reverse);

    Force::MobilityConstantForce(fwdForces,  fwdB,  0, -34.5);
    Force::MobilityConstantForce(rev1Forces, rev1A, 0,  34.5);
    Force::MobilityConstantForce(rev2Forces, rev2A, 0, -34.5); // reversed


    State fwdState  = forward.realizeTopology();
    State rev1State = reverse1.realizeTopology();
    State rev2State = reverse2.realizeTopology();

    // Put body A somewhere arbitrary.
    fwdA.setQToFitTransform(fwdState, Transform(Rotation(-1.9, Vec3(-3,2,4)),
                                                Vec3(-.33, .66, -.99)));

    fwdB.setAngle (fwdState,   Pi/4); // 45 degrees
    rev1A.setAngle(rev1State, -Pi/4);
    rev2A.setAngle(rev2State,  Pi/4);

    // Calculate where body B ended up in the forward system and then
    // set the reverse Free joint to match so that the whole system
    // is identically configured.

    forward.realize (fwdState,  Stage::Position);
    const Transform& X_GB = fwdB.getBodyTransform(fwdState);
    rev1B.setQToFitTransform(rev1State, X_GB*X_BM);
    rev2B.setQToFitTransform(rev2State, X_GB*X_BM);
    reverse1.realize(rev1State, Stage::Position);
    reverse2.realize(rev2State, Stage::Position);

    assertEqual(fwdB.getBodyTransform(fwdState), rev1B.getBodyTransform(rev1State));
    assertEqual(fwdB.getBodyTransform(fwdState), rev2B.getBodyTransform(rev2State));
    assertEqual(fwdA.getBodyTransform(fwdState), rev1A.getBodyTransform(rev1State));
    assertEqual(fwdA.getBodyTransform(fwdState), rev2A.getBodyTransform(rev2State));

    fwdB.setRate (fwdState,   3);
    forward.realize (fwdState,  Stage::Velocity);

    const SpatialVec V_G_BM( fwdB.getBodyAngularVelocity(fwdState),
                             fwdB.findStationVelocityInGround(fwdState, X_BM.p()));
    rev1B.setUToFitVelocity(rev1State, V_G_BM);
    rev2B.setUToFitVelocity(rev2State, V_G_BM);

    // Need to construct velocity of A in B from B in A.
    const Transform  X_AB = fwdB.getMobilizerTransform(fwdState);
    const Transform  X_BA = ~X_AB;
    const SpatialVec V_AB = fwdB.getMobilizerVelocity(fwdState);
    const SpatialVec V_BA( -X_BA.R()* V_AB[0],
                           -X_BA.R()*(V_AB[1] + X_AB.p() % V_AB[0]) );

    rev1A.setUToFitVelocity(rev1State, V_BA);
    rev2A.setUToFitVelocity(rev2State, V_BA); // KLUDGE; should be V_BA

    assertEqual(fwdB.getU(fwdState), rev2A.getU(rev2State));

    reverse1.realize(rev1State, Stage::Velocity);
    reverse2.realize(rev2State, Stage::Velocity);

    cout << "Vels B:\n";
    cout << "  Fwd:   " << fwdB.getBodyVelocity(fwdState) << endl;
    cout << "  Rev1: " << rev1B.getBodyVelocity(rev1State)<< endl;
    cout << "  Rev2: " << rev2B.getBodyVelocity(rev2State) << endl;

    assertEqual(fwdB.getBodyVelocity(fwdState), rev1B.getBodyVelocity(rev1State));
    assertEqual(fwdB.getBodyVelocity(fwdState), rev2B.getBodyVelocity(rev2State));

    cout << "Vels A:\n";
    cout << "  Fwd:   " << fwdA.getBodyVelocity(fwdState) << endl;
    cout << "  Rev1: " << rev1A.getBodyVelocity(rev1State)<< endl;
    cout << "  Rev2: " << rev2A.getBodyVelocity(rev2State) << endl;
    assertEqual(fwdA.getBodyVelocity(fwdState), rev1A.getBodyVelocity(rev1State));
    assertEqual(fwdA.getBodyVelocity(fwdState), rev2A.getBodyVelocity(rev2State));

    forward.realize (fwdState,  Stage::Acceleration);
    reverse1.realize(rev1State, Stage::Acceleration);
    reverse2.realize(rev2State, Stage::Acceleration);

    cout << "Accels B:\n";
    cout << "  Fwd:   " << fwdB.getBodyAcceleration(fwdState) << endl;
    cout << "  -Rev1: " << (fwdB.getBodyAcceleration(fwdState)-rev1B.getBodyAcceleration(rev1State)).norm() << endl;
    cout << "  -Rev2: " << (fwdB.getBodyAcceleration(fwdState)-rev2B.getBodyAcceleration(rev2State)).norm() << endl;

    assertEqual(fwdB.getBodyAcceleration(fwdState), rev1B.getBodyAcceleration(rev1State));
    assertEqual(fwdB.getBodyAcceleration(fwdState), rev2B.getBodyAcceleration(rev2State));
    assertEqual(fwdA.getBodyAcceleration(fwdState), rev1A.getBodyAcceleration(rev1State));
    assertEqual(fwdA.getBodyAcceleration(fwdState), rev2A.getBodyAcceleration(rev2State));

    Vector_<SpatialVec> fwdReac, rev1Reac, rev2Reac;
    fwdMatter.calcMobilizerReactionForces(fwdState, fwdReac);
    rev1Matter.calcMobilizerReactionForces(rev1State, rev1Reac);
    rev2Matter.calcMobilizerReactionForces(rev2State, rev2Reac);

    // We expect the AB reaction to be -BA reaction for a pin joint because
    // there is no translation (that is, the pin joint holds the origins
    // of the mobilizer frames coincident). This is not true for all joint
    // types!
    const SpatialVec reacBA = -fwdReac[fwdB.getMobilizedBodyIndex()];
    const SpatialVec reacBA1 = rev1Reac[rev1A.getMobilizedBodyIndex()];
    const SpatialVec reacBA2 = rev1Reac[rev1A.getMobilizedBodyIndex()];

    cout << "Reacs BA:\n";
    cout << "  Fwd:   " << reacBA << endl;
    cout << "  Rev1: "  << reacBA1 << endl;
    cout << "  Rev2: "  << reacBA2 << endl;


    assertEqual(reacBA, reacBA1, 1e-7);
    assertEqual(reacBA, reacBA2, 1e-7);
}


// System is two free bodies A and B connected by a planar joint.
// The system is Ground (6dof) A (planar) B (6dof) Ground.
//   The forward system is Ground (6dof)-> A (planar)-> B.
//   The reverse system2 is A <-(rplanar) B <-(6dof) Ground.
// Reverse system2's planar joint q's and u's should have the
// same meaning as the forward system's.
//
// Note: in this test, reversing the planar joint introduces
// sines and cosines into H where there were none before and
// hence produces a non-zero HDot where the forward direction
// has HDot==0. So this tests both H reversal and HDot reversal.
// HOWEVER: the rotational axis doesn't change so some of
// the HDot terms are not used.

void testPlanar() {
    MultibodySystem forward;
    SimbodyMatterSubsystem fwdMatter(forward);
    GeneralForceSubsystem fwdForces(forward);
    Force::UniformGravity(fwdForces, fwdMatter, Vec3(0, -1, 0));

    MultibodySystem reverse2;
    SimbodyMatterSubsystem rev2Matter(reverse2);
    GeneralForceSubsystem rev2Forces(reverse2);
    Force::UniformGravity(rev2Forces, rev2Matter, Vec3(0, -1, 0));

    const Vec3 com(1,2,3);
    const UnitInertia centralGyration(1, 1.5, 2, .1, .2, .3);
    Body::Rigid body(MassProperties(mass, com, mass*centralGyration.shiftFromMassCenter(com, 1)));

    MobilizedBody::Free fwdA (fwdMatter.Ground(),  Vec3(0), body, X_AM);
    MobilizedBody::Free rev2B(rev2Matter.Ground(), Vec3(0), body, X_BM);

    MobilizedBody::Planar fwdB (fwdA,  X_AM, body, X_BM);
    MobilizedBody::Planar rev2A(rev2B, X_BM, body, X_AM, MobilizedBody::Reverse);

    // The generalized speeds should mean the same things, so the generalized
    // forces should too.
    Force::MobilityConstantForce(fwdForces,  fwdB,  0, -3.45);
    Force::MobilityConstantForce(fwdForces,  fwdB,  1,  2);
    Force::MobilityConstantForce(fwdForces,  fwdB,  2, -3);
    Force::MobilityConstantForce(rev2Forces, rev2A, 0, -3.45); // reversed
    Force::MobilityConstantForce(rev2Forces, rev2A, 1,  2);    // reversed
    Force::MobilityConstantForce(rev2Forces, rev2A, 2, -3);    // reversed


    State fwdState  = forward.realizeTopology();
    State rev2State = reverse2.realizeTopology();

    // Put body A somewhere arbitrary.
    fwdA.setQToFitTransform(fwdState, Transform(Rotation(-1.9, Vec3(-3,2,4)),
                                                Vec3(-.33, .66, -.99)));

    // Set all the q's; meanings should be identical.
    fwdB.setQ(fwdState, Vec3(Pi/4, -7, 4));
    rev2A.setQ(rev2State, Vec3(Pi/4, -7, 4));

    // Calculate where body B ended up in the forward system and then
    // set the reverse Free joint to match so that the whole system
    // is identically configured.

    forward.realize (fwdState,  Stage::Position);
    const Transform& X_GB = fwdB.getBodyTransform(fwdState);
    const Transform X_GMb = X_GB*X_BM;
    rev2B.setQToFitTransform(rev2State, X_GMb);
    reverse2.realize(rev2State, Stage::Position);

    assertEqual(fwdB.getBodyTransform(fwdState), rev2B.getBodyTransform(rev2State));
    assertEqual(fwdA.getBodyTransform(fwdState), rev2A.getBodyTransform(rev2State));

    // Handle velocities similarly.

    fwdB.setU(fwdState, Vec3(3, .3, .4));

    forward.realize (fwdState,  Stage::Velocity);

    const SpatialVec V_G_BM( fwdB.getBodyAngularVelocity(fwdState),
                             fwdB.findStationVelocityInGround(fwdState, X_BM.p()));
    rev2B.setUToFitVelocity(rev2State, V_G_BM);

    // Need to construct velocity of A in B from B in A.
    const Transform  X_AB = fwdB.getMobilizerTransform(fwdState);
    const Transform  X_BA = ~X_AB;
    const SpatialVec V_AB = fwdB.getMobilizerVelocity(fwdState);
    const SpatialVec V_BA( -X_BA.R()* V_AB[0],
                           -X_BA.R()*(V_AB[1] + X_AB.p() % V_AB[0]) );

    rev2A.setUToFitVelocity(rev2State, V_BA);
    cout << "rev2A.getU()=" << rev2A.getU(rev2State) << endl;

    assertEqual(fwdB.getU(fwdState), rev2A.getU(rev2State));

    reverse2.realize(rev2State, Stage::Velocity);

    cout << "Vels B:\n";
    cout << "  Fwd:   " << fwdB.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyVelocity(rev2State) << endl;

    assertEqual(fwdB.getBodyVelocity(fwdState), rev2B.getBodyVelocity(rev2State));

    cout << "Vels A:\n";
    cout << "  Fwd:   " << fwdA.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyVelocity(rev2State) << endl;
    assertEqual(fwdA.getBodyVelocity(fwdState), rev2A.getBodyVelocity(rev2State));

    forward.realize (fwdState,  Stage::Acceleration);
    reverse2.realize(rev2State, Stage::Acceleration);

    cout << "Accels B:\n";
    cout << "  Fwd:  " << fwdB.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyAcceleration(rev2State) << endl;

    assertEqual(fwdB.getBodyAcceleration(fwdState), rev2B.getBodyAcceleration(rev2State));

    cout << "Accels A:\n";
    cout << "  Fwd:  " << fwdA.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyAcceleration(rev2State) << endl;
    assertEqual(fwdA.getBodyAcceleration(fwdState), rev2A.getBodyAcceleration(rev2State));

    Vector_<SpatialVec> fwdReac, rev2Reac;
    fwdMatter.calcMobilizerReactionForces(fwdState, fwdReac);
    rev2Matter.calcMobilizerReactionForces(rev2State, rev2Reac);

    // Reaction AB != -BA for a planar joint unless the translation is 0.
    // Instead, we have to shift the reaction F_AB at B's mobilizer frame
    // over to A's mobilizer frame, then negate. So
    //    F_BA = -(F_AB + [f_AB x p_BA_G; 0])
    const SpatialVec F_AB = fwdReac[fwdB.getMobilizedBodyIndex()];
    const SpatialVec reacBA = -(F_AB + SpatialVec(F_AB[1] % (X_GMb.R()*X_BA.p()), Vec3(0)));
    const SpatialVec reacBA2 = rev2Reac[rev2A.getMobilizedBodyIndex()];

    cout << "Reacs BA:\n";
    cout << "  Fwd AB:" << F_AB << " p_BA: " << X_BA.p() << endl;
    cout << "  Fwd:   " << reacBA << endl;
    cout << "  Rev2:  " << reacBA2 << endl;
    cout << "  Fwd-Rev2:" << reacBA-reacBA2 << endl;

    assertEqual(reacBA, reacBA2, 5e-6);
}


// System is two free bodies A and B connected by an ellipsoid joint.
// The system is Ground (6dof) A (ellipsoid) B (6dof) Ground.
//   The forward system is Ground (6dof)-> A (ellipsoid)-> B.
//   The reverse system2 is A <-(rellipsoid) B <-(6dof) Ground.
// Reverse system2's ellipsoid joint q's and u's should have the
// same meaning as the forward system's.
//
// H and HDot are both populated for the ellipsoid, so this
// tests the reverse HDot terms which involve the forward
// HDot terms (those don't matter if forward HDot==0).

void testEllipsoid() {
    MultibodySystem forward;
    SimbodyMatterSubsystem fwdMatter(forward);
    GeneralForceSubsystem fwdForces(forward);
    Force::UniformGravity(fwdForces, fwdMatter, Vec3(0, -1, 0));

    MultibodySystem reverse2;
    SimbodyMatterSubsystem rev2Matter(reverse2);
    GeneralForceSubsystem rev2Forces(reverse2);
    Force::UniformGravity(rev2Forces, rev2Matter, Vec3(0, -1, 0));

    const Vec3 com(1,2,3);
    const UnitInertia centralGyration(1, 1.5, 2, .1, .2, .3);
    Body::Rigid body(MassProperties(mass, com, mass*centralGyration.shiftFromMassCenter(com, 1)));

    //Transform X_AM, X_BM; // identity for now
    MobilizedBody::Free fwdA (fwdMatter.Ground(),  Vec3(0), body, X_AM);
    MobilizedBody::Free rev2B(rev2Matter.Ground(), Vec3(0), body, X_BM);

    MobilizedBody::Ellipsoid fwdB (fwdA,  X_AM, body, X_BM, Vec3(1,2,3));
    MobilizedBody::Ellipsoid rev2A(rev2B, X_BM, body, X_AM, Vec3(1,2,3), MobilizedBody::Reverse);

    // The generalized speeds should mean the same things, so the generalized
    // forces should too.
    Force::MobilityConstantForce(fwdForces,  fwdB,  0, -3.45);
    Force::MobilityConstantForce(fwdForces,  fwdB,  1,  2);
    Force::MobilityConstantForce(fwdForces,  fwdB,  2, -3);
    Force::MobilityConstantForce(rev2Forces, rev2A, 0, -3.45); // reversed
    Force::MobilityConstantForce(rev2Forces, rev2A, 1,  2);    // reversed
    Force::MobilityConstantForce(rev2Forces, rev2A, 2, -3);    // reversed


    State fwdState  = forward.realizeTopology();
    State rev2State = reverse2.realizeTopology();

    // Put body A somewhere arbitrary.
    fwdA.setQToFitTransform(fwdState, Transform(Rotation(-1.9, Vec3(-3,2,4)),
                                                Vec3(-.33, .66, -.99)));

    // Set all the q's; meanings should be identical.
    const Quaternion quat(Rotation(Pi/7, Vec3(1.1,-2.2,3.4)));
    fwdB.setQ(fwdState, quat);
    rev2A.setQ(rev2State, quat);

    // Calculate where body B ended up in the forward system and then
    // set the reverse Free joint to match so that the whole system
    // is identically configured.

    forward.realize (fwdState,  Stage::Position);
    const Transform& X_GB = fwdB.getBodyTransform(fwdState);
    const Transform X_GMb = X_GB*X_BM;
    rev2B.setQToFitTransform(rev2State, X_GMb);
    reverse2.realize(rev2State, Stage::Position);

    assertEqual(fwdB.getBodyTransform(fwdState), rev2B.getBodyTransform(rev2State));
    assertEqual(fwdA.getBodyTransform(fwdState), rev2A.getBodyTransform(rev2State));

    // Handle velocities similarly.

    fwdB.setU(fwdState, Vec3(3, .3, .4));

    forward.realize (fwdState,  Stage::Velocity);

    const SpatialVec V_G_BM( fwdB.getBodyAngularVelocity(fwdState),
                             fwdB.findStationVelocityInGround(fwdState, X_BM.p()));
    rev2B.setUToFitVelocity(rev2State, V_G_BM);

    // Need to construct velocity of A in B from B in A.
    const Transform  X_AB = fwdB.getMobilizerTransform(fwdState);
    const Transform  X_BA = ~X_AB;
    const SpatialVec V_AB = fwdB.getMobilizerVelocity(fwdState);
    const SpatialVec V_BA( -X_BA.R()* V_AB[0],
                           -X_BA.R()*(V_AB[1] + X_AB.p() % V_AB[0]) );

    rev2A.setUToFitAngularVelocity(rev2State, V_BA[0]);
    cout << "rev2A.getU()=" << rev2A.getU(rev2State) << endl;

    assertEqual(fwdB.getU(fwdState), rev2A.getU(rev2State));

    reverse2.realize(rev2State, Stage::Velocity);

    cout << "Vels B:\n";
    cout << "  Fwd:   " << fwdB.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyVelocity(rev2State) << endl;

    assertEqual(fwdB.getBodyVelocity(fwdState), rev2B.getBodyVelocity(rev2State));

    cout << "Vels A:\n";
    cout << "  Fwd:   " << fwdA.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyVelocity(rev2State) << endl;
    assertEqual(fwdA.getBodyVelocity(fwdState), rev2A.getBodyVelocity(rev2State));

    forward.realize (fwdState,  Stage::Acceleration);
    reverse2.realize(rev2State, Stage::Acceleration);

    cout << "Accels B:\n";
    cout << "  Fwd:  " << fwdB.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyAcceleration(rev2State) << endl;

    assertEqual(fwdB.getBodyAcceleration(fwdState), rev2B.getBodyAcceleration(rev2State));

    cout << "Accels A:\n";
    cout << "  Fwd:  " << fwdA.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyAcceleration(rev2State) << endl;
    assertEqual(fwdA.getBodyAcceleration(fwdState), rev2A.getBodyAcceleration(rev2State));

    Vector_<SpatialVec> fwdReac, rev2Reac;
    fwdMatter.calcMobilizerReactionForces(fwdState, fwdReac);
    rev2Matter.calcMobilizerReactionForces(rev2State, rev2Reac);

    // Reaction AB != -BA for a planar joint unless the translation is 0.
    // Instead, we have to shift the reaction F_AB at B's mobilizer frame
    // over to A's mobilizer frame, then negate. So
    //    F_BA = -(F_AB + [f_AB x p_BA_G; 0])
    const SpatialVec F_AB = fwdReac[fwdB.getMobilizedBodyIndex()];
    const SpatialVec reacBA = -(F_AB + SpatialVec(F_AB[1] % (X_GMb.R()*X_BA.p()), Vec3(0)));
    const SpatialVec reacBA2 = rev2Reac[rev2A.getMobilizedBodyIndex()];

    cout << "Reacs BA:\n";
    cout << "  Fwd AB:" << F_AB << " p_BA: " << X_BA.p() << endl;
    cout << "  Fwd:   " << reacBA << endl;
    cout << "  Rev2:  " << reacBA2 << endl;
    cout << "  Fwd-Rev2:" << reacBA-reacBA2 << endl;

    assertEqual(reacBA, reacBA2, reacBA.norm()*1e-8);
}

// This one should test all the mobilizer reversal equations.
// System is two free bodies A and B connected by a Free joint.
// The system is Ground (6dof) A (free) B (6dof) Ground.
//   The forward system is Ground (6dof)-> A (free)-> B.
//   The reverse system2 is A <-(rfree) B <-(6dof) Ground.
// Reverse system2's free joint q's and u's should have the
// same meaning as the forward system's.
// NOTE: 6dof and free are the same thing above, but the one named
// "free" and "rfree" is the mobilizer under test.
//
// Note: in this test, reversing the free joint introduces
// sines and cosines into H where there were none before and
// hence produces a non-zero HDot where the forward direction
// has HDot==0. So this tests both H reversal and HDot reversal.

void testFree() {
    MultibodySystem forward;
    SimbodyMatterSubsystem fwdMatter(forward);
    GeneralForceSubsystem fwdForces(forward);
    Force::UniformGravity(fwdForces, fwdMatter, Vec3(0, -1, 0));

    MultibodySystem reverse2;
    SimbodyMatterSubsystem rev2Matter(reverse2);
    GeneralForceSubsystem rev2Forces(reverse2);
    Force::UniformGravity(rev2Forces, rev2Matter, Vec3(0, -1, 0));

    const Vec3 com(1,2,3);
    const UnitInertia centralGyration(1, 1.5, 2, .1, .2, .3);
    Body::Rigid body(MassProperties(mass, com, mass*centralGyration.shiftFromMassCenter(com, 1)));

    MobilizedBody::Free fwdA (fwdMatter.Ground(),  Vec3(0), body, X_AM);
    MobilizedBody::Free rev2B(rev2Matter.Ground(), Vec3(0), body, X_BM);

    // These are the mobilizers under test.
    MobilizedBody::Free fwdB (fwdA,  X_AM, body, X_BM);
    MobilizedBody::Free rev2A(rev2B, X_BM, body, X_AM, MobilizedBody::Reverse);

    // The generalized speeds should mean the same things, so the generalized
    // forces should too.
    const Real genFrc[6] = {-3.45, 2, -3, .0232, 4.33, -2.4};
    for (int i=0; i<6; ++i) {
        // We're creating force elements here and adding them to the
        // force subsystems.
        Force::MobilityConstantForce(fwdForces,   fwdB,  i, 1000*genFrc[i]);
        Force::MobilityConstantForce(rev2Forces,  rev2A,  i, 1000*genFrc[i]);
    }

    State fwdState  = forward.realizeTopology();
    State rev2State = reverse2.realizeTopology();

    // Put body A somewhere arbitrary.
    fwdA.setQToFitTransform(fwdState, Transform(Rotation(-1.9, Vec3(-3,2,4)),
                                                Vec3(-.33, .66, -.99)));

    // Set all the q's; meanings should be identical.
    const Quaternion quat(Rotation(Pi/7, Vec3(1,2,-3)));
    const Vec3       trans(1.23, 2.34, -3.45);
    const Vec7       init(quat[0],quat[1],quat[2],quat[3],trans[0],trans[1],trans[2]);
    fwdB.setQ (fwdState,  init);
    rev2A.setQ(rev2State, init);

    // Calculate where body B ended up in the forward system and then
    // set the reverse Free joint to match so that the whole system
    // is identically configured.

    forward.realize (fwdState,  Stage::Position);
    const Transform& X_GB = fwdB.getBodyTransform(fwdState);
    const Transform X_GMb = X_GB*X_BM;
    rev2B.setQToFitTransform(rev2State, X_GMb);
    reverse2.realize(rev2State, Stage::Position);

    cout << "Xforms B\n";
    cout << "  Fwd:  " << fwdB.getBodyTransform(fwdState);
    cout << "  Rev2: " << rev2B.getBodyTransform(rev2State);
    assertEqual(fwdB.getBodyTransform(fwdState), rev2B.getBodyTransform(rev2State));

    cout << "Xforms A\n";
    cout << "  Fwd:  " << fwdA.getBodyTransform(fwdState);
    cout << "  Rev2: " << rev2A.getBodyTransform(rev2State);
    assertEqual(fwdA.getBodyTransform(fwdState), rev2A.getBodyTransform(rev2State));

    // Handle velocities similarly.

    // Want high velocity so cross terms involving HDot will be big enough
    // to notice if they are wrong.
    fwdB.setU(fwdState, 1000*Vec6(.1, .2, -.33, 3, .3, .4));

    forward.realize (fwdState,  Stage::Velocity);

    const SpatialVec V_G_BM( fwdB.getBodyAngularVelocity(fwdState),
                             fwdB.findStationVelocityInGround(fwdState, X_BM.p()));
    rev2B.setUToFitVelocity(rev2State, V_G_BM);

    // Need to construct velocity of A in B from B in A.
    const Transform  X_AB = fwdB.getMobilizerTransform(fwdState);
    const Transform  X_BA = ~X_AB;
    const SpatialVec V_AB = fwdB.getMobilizerVelocity(fwdState);
    const SpatialVec V_BA( -X_BA.R()* V_AB[0],
                           -X_BA.R()*(V_AB[1] + X_AB.p() % V_AB[0]) );

    rev2A.setUToFitVelocity(rev2State, V_BA);
    cout << "rev2A.getU()=" << rev2A.getU(rev2State) << endl;

    assertEqual(fwdB.getU(fwdState), rev2A.getU(rev2State));

    reverse2.realize(rev2State, Stage::Velocity);

    cout << "Vels B:\n";
    cout << "  Fwd:   " << fwdB.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyVelocity(rev2State) << endl;

    assertEqual(fwdB.getBodyVelocity(fwdState), rev2B.getBodyVelocity(rev2State));

    cout << "Vels A:\n";
    cout << "  Fwd:   " << fwdA.getBodyVelocity(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyVelocity(rev2State) << endl;
    assertEqual(fwdA.getBodyVelocity(fwdState), rev2A.getBodyVelocity(rev2State));

    forward.realize (fwdState,  Stage::Acceleration);
    reverse2.realize(rev2State, Stage::Acceleration);

    cout << "Accels B:\n";
    cout << "  Fwd:  " << fwdB.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2B.getBodyAcceleration(rev2State) << endl;
    cout << "  Fwd-Rev2: " << fwdB.getBodyAcceleration(fwdState)-rev2B.getBodyAcceleration(rev2State) << endl;

    assertEqual(fwdB.getBodyAcceleration(fwdState), rev2B.getBodyAcceleration(rev2State), 1e-7);

    cout << "Accels A:\n";
    cout << "  Fwd:  " << fwdA.getBodyAcceleration(fwdState) << endl;
    cout << "  Rev2: " << rev2A.getBodyAcceleration(rev2State) << endl;
    cout << "  Fwd-Rev2: " << fwdA.getBodyAcceleration(fwdState)-rev2A.getBodyAcceleration(rev2State) << endl;
    assertEqual(fwdA.getBodyAcceleration(fwdState), rev2A.getBodyAcceleration(rev2State), 1e-7);

    Vector_<SpatialVec> fwdReac, rev2Reac;
    fwdMatter.calcMobilizerReactionForces(fwdState, fwdReac);
    rev2Matter.calcMobilizerReactionForces(rev2State, rev2Reac);

    // Reaction AB != -BA for a free joint unless the translation is 0.
    // Instead, we have to shift the reaction F_AB at B's mobilizer frame
    // over to A's mobilizer frame, then negate. So
    //    F_BA = -(F_AB + [f_AB x p_BA_G; 0])
    const SpatialVec F_AB = fwdReac[fwdB.getMobilizedBodyIndex()];
    const SpatialVec reacBA = -(F_AB + SpatialVec(F_AB[1] % (X_GMb.R()*X_BA.p()), Vec3(0)));
    const SpatialVec reacBA2 = rev2Reac[rev2A.getMobilizedBodyIndex()];

    cout << "Reacs BA:\n";
    cout << "  Fwd AB:" << F_AB << " p_BA: " << X_BA.p() << endl;
    cout << "  Fwd:   " << reacBA << endl;
    cout << "  Rev2:  " << reacBA2 << endl;
    cout << "  Fwd-Rev2:" << reacBA-reacBA2 << endl;

    cout << "  Fwd: " << fwdB.findMobilizerReactionOnParentAtFInGround(fwdState) << endl;
    cout << "  Rev: " << rev2A.findMobilizerReactionOnBodyAtMInGround(rev2State) << endl;

    assertEqual(reacBA, reacBA2, reacBA.norm()*1e-5);
}

int main() {
    try {
        cout << "*** TEST PIN ***\n\n"; testPin();
        cout << "\n\n*** TEST PLANAR ***\n\n"; testPlanar();
        cout << "\n\n*** TEST ELLIPSOID ***\n\n"; testEllipsoid();
        cout << "\n\n*** TEST FREE ***\n\n"; testFree();
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

