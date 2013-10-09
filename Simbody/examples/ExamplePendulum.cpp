/* -------------------------------------------------------------------------- *
 *                       Simbody(tm) Example: Pendulum                        *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Ajay Seth                                                         *
 * Contributors: Michael Sherman                                              *
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


// Define two identical double pendulums, one modeled the easy way using
// two pin mobilizers, the other modeled with free mobilizers plus a ball
// constraint, plus two "constant angle" constraints to get rid of the extra 
// rotational degrees of freedom.
// We're going to show that the resulting reaction forces are identical.


#include "Simbody.h"

using namespace SimTK;
using std::cout; using std::endl;

int main() {
  try {   
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    // Add gravity as a force element.
    Rotation x45(Pi/4, XAxis);
    Rotation y45(Pi/4, YAxis);
    Rotation z45(Pi/4, ZAxis);
    Force::UniformGravity gravity(forces, matter, Vec3(10, Real(-9.8), 3));
    // Create the body and some artwork for it.
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), 
                               DecorativeSphere(Real(0.1)).setColor(Red));

    // Add an instance of the body to the multibody system by connecting
    // it to Ground via a pin mobilizer.
    MobilizedBody::Pin pendulum1(matter.updGround(), 
                                Transform(/*x45,*/Vec3(0,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum1b(pendulum1, 
                                Transform(/*x45,*/Vec3(0,-1,0)), 
                                pendulumBody, 
                                Transform(Vec3(0, 1, 0)));

    // Now make an identical system with pin joints faked up using a
    // free mobilizer + ball constraint + two angle constraints.
    MobilizedBody::Free pendulum2(matter.updGround(), 
                                  Transform(/*x45,*/Vec3(2,-1,0)),
                                  pendulumBody, 
                                  Transform(Vec3(0,1,0)));
    Constraint::Ball ballcons2(matter.updGround(), Vec3(2,-1,0),
                               pendulum2, Vec3(0,1,0));
    const Transform& X_GF2 = pendulum2.getDefaultInboardFrame();
    const Transform& X_P2M = pendulum2.getDefaultOutboardFrame();
    Constraint::ConstantAngle angx2(matter.Ground(), X_GF2.x(),
                              pendulum2, X_P2M.z());
    Constraint::ConstantAngle angy2(matter.Ground(), X_GF2.y(),
                              pendulum2, X_P2M.z());

    MobilizedBody::Free pendulum2b(pendulum2, 
                                   Transform(/*x45,*/Vec3(0,-1,0)),
                                   pendulumBody, 
                                   Transform(Vec3(0,1,0)));
    Constraint::Ball ballcons2b(pendulum2, Vec3(0,-1,0),
                                pendulum2b, Vec3(0,1,0));
    const Transform& X_GF2b = pendulum2b.getDefaultInboardFrame();
    const Transform& X_P2Mb = pendulum2b.getDefaultOutboardFrame();
    Constraint::ConstantAngle angx2b(pendulum2, X_GF2b.x(),
                              pendulum2b, X_P2Mb.z());
    Constraint::ConstantAngle angy2b(pendulum2, X_GF2b.y(),
                              pendulum2b, X_P2Mb.z());


    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, Real(1./30)));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum1.setOneQ(state, 0, Pi/4);
    //pendulum1.setOneU(state, 0, 1.0); // initial velocity 1 rad/sec

    //pendulum1b.setOneU(state, 0, 1.0); // initial velocity 1 rad/sec
    pendulum1b.setOneQ(state, 0, Pi/4);

    pendulum2.setQToFitRotation(state, Rotation(Pi/4, ZAxis));
    //pendulum2.setUToFitAngularVelocity(state, Vec3(0,0,1));
    pendulum2b.setQToFitRotation(state, Rotation(Pi/4, ZAxis));
    //pendulum2b.setUToFitAngularVelocity(state, Vec3(0,0,1));

    system.realize(state);
    const Vector lambda = state.getMultipliers();
    Vector_<SpatialVec> consBodyForcesInG;
    Vector              consMobForces;
    matter.calcConstraintForcesFromMultipliers(state, -lambda, consBodyForcesInG,
        consMobForces);
    const MobodIndex p2x = pendulum2.getMobilizedBodyIndex();
    const MobodIndex p2bx = pendulum2b.getMobilizedBodyIndex();
    const Rotation& R_G2 = pendulum2.getBodyTransform(state).R();
    //consBodyForcesInG[p2x] = shiftForceFromTo(consBodyForcesInG[p2x],
     //                                         Vec3(0), R_G2*Vec3(0,1,0));

    const int nb = matter.getNumBodies();
    Vector_<SpatialVec> forcesAtMInG;
    matter.calcMobilizerReactionForces(state, forcesAtMInG);

    // The above method returns reactions *on the child body* at the outboard
    // mobilizer frame M (that is, the frame fixed on the child body). 
    // Calculate the same reactions, but *on the parent body* and at the 
    // inboard mobilizer frame F (that is, the frame fixed on the parent body). 
    // This is done by shifting the reaction forces across the mobilizer from 
    // M to F, and negating.
    // Note that for the example here the mobilizers aren't translating so
    // the origins Mo and Fo are identical. Nevertheless we're using the
    // generic code here that should work for any system.

    Vector_<SpatialVec> forcesAtFInG(nb); // to hold the result
    forcesAtFInG[0] = -forcesAtMInG[0]; // Ground is "welded" at origin
    for (MobilizedBodyIndex i(1); i < matter.getNumBodies(); ++i) {
        const MobilizedBody& body   = matter.getMobilizedBody(i);
        const MobilizedBody& parent = body.getParentMobilizedBody();
        // Want to shift negated reaction by p_MF_G, the vector from M
        // to F across the mobilizer, expressed in Ground. We can get p_FM, 
        // then re-express in Ground for the shift and negate.
        const Vec3& p_FM = body.getMobilizerTransform(state).p();
        const Rotation& R_PF = body.getInboardFrame(state).R(); // In parent.
        const Rotation& R_GP = parent.getBodyTransform(state).R();
        Rotation R_GF   =   R_GP*R_PF;  // F frame orientation in Ground.
        Vec3     p_MF_G = -(R_GF*p_FM); // Re-express and negate shift vector.
        forcesAtFInG[i] = -shiftForceBy(forcesAtMInG[i], p_MF_G);
    }

    std::cout << "Reactions @M: " << forcesAtMInG << "\n";
    std::cout << "Reactions @F: " << forcesAtFInG << "\n";
    std::cout << "norm of difference: " << (forcesAtMInG+forcesAtFInG).norm() 
              << "\n";

    const MobodIndex p1x = pendulum1.getMobilizedBodyIndex();
    const MobodIndex p1bx = pendulum1b.getMobilizedBodyIndex();
    const Rotation& R_G1 = pendulum1.getBodyTransform(state).R();
    const Rotation& R_G1b = pendulum1b.getBodyTransform(state).R();

    // This time shift the child reactions from the mobilizer frames M to the
    // child body frames B so we can compare with the constraint forces that
    // are returned at the child body frame.
    Vector_<SpatialVec> forcesAtBInG(nb);
    for (MobodIndex i(0); i < nb; ++i) {
        const Mobod& body = matter.getMobilizedBody(i);
        const Vec3&  p_BM = body.getOutboardFrame(state).p();
        const Rotation& R_GB = body.getBodyTransform(state).R();
        forcesAtBInG[i] = shiftForceFromTo(forcesAtMInG[i],
                                           R_GB*p_BM, Vec3(0));
    }

    std::cout << "Pin mobilizer reaction forces:\n";
    std::cout << "FB_G=" << forcesAtBInG[p1x] 
              << " " << forcesAtBInG[p1bx] << "\n";

    std::cout << "Constraint reaction forces (should be the same):\n";
    cout << "FC_G=" << -(ballcons2.getConstrainedBodyForcesAsVector(state)
        + angx2.getConstrainedBodyForcesAsVector(state)
        + angy2.getConstrainedBodyForcesAsVector(state))[1] << " ";
    cout << -(ballcons2b.getConstrainedBodyForcesAsVector(state) 
        + angx2b.getConstrainedBodyForcesAsVector(state)
        + angy2b.getConstrainedBodyForcesAsVector(state))[1] << endl;

    viz.report(state);
    // Simulate it.
    cout << "Hit ENTER to run a short simulation ...";
    getchar();

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);
    state = integ.getState();
    system.realize(state);
    matter.calcMobilizerReactionForces(state, forcesAtMInG);
    const Transform& X_GP = pendulum1.getBodyTransform(state);
    //forcesAtMInG[1][1] = X_GP.R()*forcesAtMInG[1][1];
    std::cout << "FM_G=" << forcesAtMInG << "\n";
    ts.stepTo(Real(1.2));
    state = integ.getState();
    system.realize(state);
    matter.calcMobilizerReactionForces(state, forcesAtMInG);
    std::cout << "FM_G=" << forcesAtMInG << "\n";

  } catch (const std::exception& e) {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }

    return 0;
}
