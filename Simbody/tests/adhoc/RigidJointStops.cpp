/* -------------------------------------------------------------------------- *
 *             Simbody(tm) Example: SimplePlanarMechanism                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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
#include "Simbody.h"
#include <iostream>

using namespace SimTK;
using std::cout; using std::endl;

// This very simple example builds a 3-body planar mechanism that does 
// nothing but rock back and forth for 10 seconds. Note that Simbody always
// works in 3D; this mechanism is planar because of the alignment of its 
// joints not because it uses any special 2D features. The mechanism looks
// like this:
//                              @
//                     @--------+--------@
//     Y               |        |         \
//     |               |        |          \
//     |               |        |           \
//     /-----X         *        *            *
//    /
//    Z
//
// It consists of a central T-shaped body pinned to ground, and
// two pendulum bodies pinned to either side of the T. The @'s above represent
// pin joints rotating about the Z axes. Each body's mass is concentrated into
// point masses shown by *'s above. Gravity is in the -Y direction.
//
// We'll add rigid joint stops to the two arms as a very simple test of
// unilateral contact constraints.


//==============================================================================
//                              SHOW ENERGY
//==============================================================================
// Generate text in the scene that displays the total energy, which should be 
// conserved to roughly the number of decimal places corresponding to the 
// accuracy setting (i.e., acc=1e-5 -> 5 digits).
class ShowEnergy : public DecorationGenerator {
public:
    explicit ShowEnergy(const MultibodySystem& mbs) : m_mbs(mbs) {}
    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) override;
private:
    const MultibodySystem& m_mbs;
};

namespace {
const Real StopCoefRestLeft = 1;
const Real StopCoefRestRight = .5;
const Real RopeCoefRest = .1;
const Real BoxCoefRest = .925;
const Real BoxMuS = Infinity;
const Real BoxMuD = Infinity;
const Real BoxMuV = 0;
const Real Gravity = 9.81;
}

int main() {
  try { // catch errors if any

    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; 
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, Gravity);
    MobilizedBody& Ground = matter.updGround();
   // matter.setShowDefaultGeometry(false);

    // Describe a body with a point mass at (0, -3, 0) and draw a sphere there.
    Real mass = 3; Vec3 pos(0,-3,0);
    Body::Rigid bodyInfo(MassProperties(mass, pos, UnitInertia::pointMassAt(pos)));
    bodyInfo.addDecoration(pos, DecorativeSphere(.2).setOpacity(.5));

    Body::Rigid heavyInfo(MassProperties(10*mass, pos, UnitInertia::pointMassAt(pos)));
    heavyInfo.addDecoration(pos, DecorativeSphere(.3).setOpacity(.5));

    // Create the tree of mobilized bodies, reusing the above body description.
    MobilizedBody::Pin bodyT  (matter.Ground(), Vec3(0,4,0), heavyInfo, Vec3(0));
    MobilizedBody::Pin leftArm(bodyT, Vec3(-2, 0, 0),    bodyInfo, Vec3(0));
    MobilizedBody::Pin rtArm  (bodyT, Vec3(2, 0, 0),     bodyInfo, Vec3(0));

    DecorativeLine stop(Vec3(0), Vec3(0,-2,0));
    stop.setColor(Red).setLineThickness(1);
    const Real inner=.3, outer=1.;
    bodyT.addBodyDecoration(Transform(Rotation(-inner, ZAxis), Vec3(-2,0,0)),
                            stop);
    bodyT.addBodyDecoration(Transform(Rotation(-outer, ZAxis), Vec3(-2,0,0)),
                            stop);
    bodyT.addBodyDecoration(Transform(Rotation(inner, ZAxis), Vec3(2,0,0)),
                            stop);
    bodyT.addBodyDecoration(Transform(Rotation(outer, ZAxis), Vec3(2,0,0)),
                            stop);

    UnilateralContactIndex ll, ul, lr, ur, rope;
    ll=matter.adoptUnilateralContact(
        new HardStopLower(leftArm, MobilizerQIndex(0), -outer, StopCoefRestLeft));
    ul=matter.adoptUnilateralContact(
        new HardStopUpper(leftArm, MobilizerQIndex(0), -inner, StopCoefRestLeft));

    lr=matter.adoptUnilateralContact(
        new HardStopLower(rtArm, MobilizerQIndex(0), inner, StopCoefRestRight));
    ur=matter.adoptUnilateralContact(
        new HardStopUpper(rtArm, MobilizerQIndex(0), outer, StopCoefRestRight));

    // ROPE
    rope=matter.adoptUnilateralContact(
        new Rope(matter.Ground(), Vec3(-5.25,3,1),
                 leftArm, Vec3(0, -1.5, 0), 2.75, RopeCoefRest));

    // BOX
    const Vec3 hBox(1,.25,.5); // half dims of box
    //Ground.addBodyDecoration(Vec3(0,.05,0), DecorativeFrame(2).setColor(Green));
    DecorativeBrick drawBox(hBox); drawBox.setOpacity(0.5).setColor(Gray);

    Body::Rigid boxBody(MassProperties(1, Vec3(0), UnitInertia::brick(hBox))); 
    boxBody.addDecoration(Vec3(0),drawBox);
    boxBody.addDecoration(Vec3(0), DecorativePoint());

    MobilizedBody::Free box(Ground, Vec3(0,2,0), boxBody, Vec3(0));

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(0) + Vec3(i,j,k).elementwiseMultiply(hBox);
        PointPlaneContact* contact = new PointPlaneContact
           (Ground, YAxis, 0., box, pt, BoxCoefRest , BoxMuS, BoxMuD, BoxMuV);
        matter.adoptUnilateralContact(contact);
        //i=j=k=2; // break 3
    }

    // Ask for visualization every 1/30 second.
    //system.setUseUniformBackground(true); // turn off floor 
    Visualizer viz(system);
    viz.setShowSimTime(true);
    viz.setDesiredFrameRate(30.);
    viz.addDecorationGenerator(new ShowEnergy(system));
    system.adoptEventReporter(new Visualizer::Reporter(viz, 1./30));

    class ReportQ : public TextDataEventReporter::UserFunction<Vector> {
    public:
        Vector evaluate(const System& sys, const State& state) {
            return state.getQ();
        }
    };
    //system.adoptEventReporter(new TextDataEventReporter(system,
    //                                                    new ReportQ(),
    //                                                    1./30));
    
    // Initialize the system and state.    
    State state = system.realizeTopology();
    //bodyT.lock(state);
    //bodyT.lockAt(state, -.1, Motion::Velocity);
    //bodyT.setRate(state, -2);
    //leftArm.setAngle(state, -Pi/2);
    //rtArm.setAngle(state, Pi/2);
    leftArm.setAngle(state, -.9);
    rtArm.setAngle(state, .9);
    box.setQToFitRotation(state, Rotation(BodyRotationSequence,
                                          Pi/4, ZAxis, -Pi/10, XAxis));
    box.setUToFitLinearVelocity(state, Vec3(-2,0,0));
    box.setUToFitAngularVelocity(state, Vec3(0,5,2));

    const double SimTime = 20;

    //matter.getUnilateralContact(ul).setCondition(state, CondConstraint::Active);
    //matter.getUnilateralContact(ur).setCondition(state, CondConstraint::Active);
    matter.getUnilateralContact(rope).setCondition(state, CondConstraint::Active);

    printf("SHOWING UNASSEMBLED SYSTEM -- hit ENTER\n");
    viz.report(state);
    getchar();


    Assembler(system).setErrorTolerance(1e-6).assemble(state);
    printf("SHOWING ASSEMBLED SYSTEM -- hit ENTER\n");
    viz.report(state);
    getchar();

    // Simulate with acceleration-level time stepper.
    SemiExplicitEuler2Integrator integ(system);
    //RungeKutta2Integrator integ(system);
    //RungeKutta3Integrator integ(system);
    //RungeKuttaMersonIntegrator integ(system);
    //CPodesIntegrator integ(system);
    //integ.setAllowInterpolation(false);
    //integ.setAccuracy(0.01);
    integ.setAccuracy(1e-2);
    integ.setConstraintTolerance(1e-4);
    integ.setMaximumStepSize(.1);
    //integ.setPrecision(1e-6);
    TimeStepper ts(integ);

    const double startReal = realTime(), startCPU=cpuTime();
    ts.initialize(state);
    ts.stepTo(SimTime);
    const double timeInSec = realTime()-startReal, 
                 cpuInSec = cpuTime()-startCPU;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s for " << ts.getTime() << "s sim (avg step=" 
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*ts.getTime())/evals << "ms/eval\n";
    cout << "CPUtime " << cpuInSec << endl;

    printf("Used Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());
  
    //printf("# EVENT STEPS = %d\n", nStepsWithEvent);

    //// Simulate with velocity-level time stepper.
    //SemiExplicitEulerTimeStepper sxe(system);
    //sxe.initialize(state);
    //const Real h = .001;
    //const int DrawEvery = 33; // draw every nth step ~= 33ms
    //int nSteps = 0;
    //while (true) {
    //    const State& sxeState = sxe.getState();
    //    const double t = sxe.getTime();
    //    const bool isLastStep = t + h/2 > SimTime;

    //    if((nSteps%DrawEvery)==0 || isLastStep) {
    //        viz.report(sxeState);
    //        //printf("%7.4f %9.3g %9.3g\n", t,
    //        //       leftArm.getAngle(sxeState),
    //        //       rtArm.getAngle(sxeState));
    //    }

    //    if (isLastStep)
    //        break;

    //    sxe.stepTo(sxeState.getTime() + h);
    //    ++nSteps;
    //}

  } catch (const std::exception& e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return 1;
  }

    return 0;
}



void ShowEnergy::generateDecorations(const State&                state, 
                                     Array_<DecorativeGeometry>& geometry)
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(state, Stage::Dynamics);
    const Real E=m_mbs.calcEnergy(state);
    const SpatialVec mom = matter.calcSystemMomentumAboutGroundOrigin(state);
 
    DecorativeText energy;
    energy.setText("Energy: " + String(E, "%.6f"))
        .setIsScreenText(true);
    geometry.push_back(energy);

    for (UnilateralContactIndex ucx(0); 
         ucx < matter.getNumUnilateralContacts(); ++ucx)
    {
        const UnilateralContact& uni = matter.getUnilateralContact(ucx);
        CondConstraint::Condition cond = uni.getCondition(state);
        geometry.push_back(DecorativeText()
                           .setText(CondConstraint::toString(cond))
                           .setIsScreenText(true));
    }
}
