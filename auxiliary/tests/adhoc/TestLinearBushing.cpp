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

/**@file
 * Test the Force::LinearBushing force element.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"   // requires VTK

#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

class BushingReporter : public PeriodicEventReporter {
public:
    BushingReporter(const MultibodySystem& sys,
                    const Force::LinearBushing& frc,
                    const Force::LinearBushing& bush2,
                    Real dt) 
    :    PeriodicEventReporter(dt), system(sys), bushing(frc), bushing2(bush2) {}

    void handleEvent(const State& state) const {
        printf("BUSHING t=%g, stage %s, energy=%g CONSERVED=%g\n", state.getTime(),
                state.getSystemStage().getName().c_str(),
                system.calcEnergy(state),
                system.calcEnergy(state)+bushing.getDissipatedEnergy(state)
                        +bushing2.getDissipatedEnergy(state));
        cout << "q=" << bushing.getQ(state) << endl;
        cout << "qdot=" << bushing.getQDot(state) << endl;
        cout << "X_FM=" << bushing.getX_FM(state);
        cout << "V_FM=" << bushing.getV_FM(state) << endl;
        cout << "f=" << bushing.getF(state) << endl;
        cout << "F_GM=" << bushing.getF_GM(state) << endl;
        cout << "F_GF=" << bushing.getF_GF(state) << endl;
        cout << "pe=" << bushing.getPotentialEnergy(state) << endl;
        cout << "power=" << bushing.getPowerDissipation(state) << endl;
        cout << "e_dissipated=" << bushing.getDissipatedEnergy(state) << endl;
    }
private:
    const MultibodySystem&         system;
    const Force::LinearBushing&    bushing;
    const Force::LinearBushing&    bushing2;
};

class ReactionReporter : public PeriodicEventReporter {
public:
    ReactionReporter(std::string name,
                     const MultibodySystem& system,
                     const MobilizedBody& mobod,
                     Real  dt)
    :    PeriodicEventReporter(dt), name(name), system(system), mobod(mobod) {}

    void handleEvent(const State& state) const {
        const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();
        printf("MOBILIZER t=%g\n", state.getTime());
        system.realize(state, Stage::Acceleration);
        Vector_<SpatialVec> reactions;
        matter.calcMobilizerReactionForces(state, reactions);
        cout << "q=" << mobod.getQAsVector(state) << endl;
        cout << "qdot=" << mobod.getQDotAsVector(state) << endl;

        cout << "Reaction for " << name << ": " 
             << reactions[mobod.getMobilizedBodyIndex()] << endl;
        cout << "  X_FM=" << mobod.getMobilizerTransform(state);
        cout << "  X_GB=" << mobod.getBodyTransform(state);
    }
private:
    const std::string        name;
    const MultibodySystem&   system;
    const MobilizedBody&     mobod;
};


int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, 0*Vec3(0, -9.8, 0));

    const Real Mass = 1;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;
    const Transform BodyAttach(Rotation(), Vec3(HalfShape[0],0,0));
    Body::Rigid brickBody(MassProperties(Mass, Vec3(.1,.2,.3), 
                                Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    //Body::Rigid brickBody(MassProperties(Mass, Vec3(0), 
    //                        Mass*Gyration::ellipsoid(HalfShape)));
    brickBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    brickBody.addDecoration(BodyAttach,
                DecorativeFrame(0.5).setColor(Red));

    MobilizedBody::Free brick1(matter.Ground(), Transform(), 
                               brickBody,       BodyAttach);

    Body::Rigid massless(MassProperties(0, Vec3(0), Inertia(0)));
    const Rotation ZtoX(Pi/2, YAxis);
    const Rotation ZtoY(-Pi/2, XAxis);
    //MobilizedBody::Cartesian dummy1(matter.Ground(), Vec3(1,1,1)+Vec3(2,0,0),
    //                                massless, Transform());
    //MobilizedBody::Pin dummy2(dummy1,   ZtoX,
    //                          massless, ZtoX);
    //MobilizedBody::Pin dummy3(dummy2,   ZtoY,
    //                          massless, ZtoY);
    //MobilizedBody::Pin brick2(dummy3,   Transform(),    // about Z
    //                          brickBody, BodyAttach);

    MobilizedBody::Pin dummy1(matter.Ground(), Vec3(1,1,1)+Vec3(2,0,0),
                              massless, Transform(), MobilizedBody::Reverse);
    MobilizedBody::Pin dummy2(dummy1,   ZtoY,
                              massless, ZtoY, MobilizedBody::Reverse);
    MobilizedBody::Pin dummy3(dummy2,   ZtoX,
                              massless, ZtoX, MobilizedBody::Reverse);
    MobilizedBody::Cartesian brick2(dummy3,   Transform(),
                                    brickBody, BodyAttach, MobilizedBody::Reverse);

    const Vec6 k = 10*Vec6(1,1,1,1,1,1);
    const Vec6 c =  1*Vec6(1,1,1,1,1,1);
    Transform GroundAttach(Rotation(), Vec3(1,1,1));
    matter.Ground().updBody().addDecoration(GroundAttach,
                DecorativeFrame(0.5).setColor(Green));

    const Vec6 k1 = k;
    const Vec6 c1 = c;
    //Force::LinearBushing bushing
    //    (forces, matter.Ground(), GroundAttach,
    //     brick1, BodyAttach, k1, c1);
    Force::LinearBushing bushing
        (forces, brick1, BodyAttach, 
         matter.Ground(), GroundAttach,
         k1, c1);

    Force::LinearBushing bushing2
        (forces, brick1, brick2, k, c);

    Transform GroundAttach2(Rotation(), Vec3(1,1,1) + Vec3(2,0,0));
    matter.Ground().updBody().addDecoration(GroundAttach2,
                DecorativeFrame(0.5).setColor(Green));
    const Vec6 k2 = k;
    const Vec6 c2 = 0*c;
    //Force::MobilityLinearSpring kqx(forces, dummy2, 0, k2[0], 0);
    //Force::MobilityLinearDamper cqx(forces, dummy2, 0, c2[0]);
    //Force::MobilityLinearSpring kqy(forces, dummy3, 0, k2[1], 0);
    //Force::MobilityLinearDamper cqy(forces, dummy3, 0, c2[1]);
    //Force::MobilityLinearSpring kqz(forces, brick2, 0, k2[2], 0);
    //Force::MobilityLinearDamper cqz(forces, brick2, 0, c2[2]);
    //Force::MobilityLinearSpring kpx(forces, dummy1, 0, k2[3], 0);
    //Force::MobilityLinearDamper cpx(forces, dummy1, 0, c2[3]);
    //Force::MobilityLinearSpring kpy(forces, dummy1, 1, k2[4], 0);
    //Force::MobilityLinearDamper cpy(forces, dummy1, 1, c2[4]);
    //Force::MobilityLinearSpring kpz(forces, dummy1, 2, k2[5], 0);
    //Force::MobilityLinearDamper cpz(forces, dummy1, 2, c2[5]);

    Force::MobilityLinearSpring kqx(forces, dummy3, 0, k2[0], 0);
    Force::MobilityLinearDamper cqx(forces, dummy3, 0, c2[0]);
    Force::MobilityLinearSpring kqy(forces, dummy2, 0, k2[1], 0);
    Force::MobilityLinearDamper cqy(forces, dummy2, 0, c2[1]);
    Force::MobilityLinearSpring kqz(forces, dummy1, 0, k2[2], 0);
    Force::MobilityLinearDamper cqz(forces, dummy1, 0, c2[2]);
    Force::MobilityLinearSpring kpx(forces, brick2, 0, k2[3], 0);
    Force::MobilityLinearDamper cpx(forces, brick2, 0, c2[3]);
    Force::MobilityLinearSpring kpy(forces, brick2, 1, k2[4], 0);
    Force::MobilityLinearDamper cpy(forces, brick2, 1, c2[4]);
    Force::MobilityLinearSpring kpz(forces, brick2, 2, k2[5], 0);
    Force::MobilityLinearDamper cpz(forces, brick2, 2, c2[5]);

    VTKEventReporter* reporter = new VTKEventReporter(system, 0.01);
    system.updDefaultSubsystem().addEventReporter(reporter);

    BushingReporter* bushingReport = new BushingReporter(system, bushing, bushing2, .01);
    system.updDefaultSubsystem().addEventReporter(bushingReport);

    ReactionReporter* dummy1Reac = new ReactionReporter("dummy1", system, dummy1, .01);
    ReactionReporter* dummy2Reac = new ReactionReporter("dummy2", system, dummy2, .01);
    ReactionReporter* dummy3Reac = new ReactionReporter("dummy3", system, dummy3, .01);
    ReactionReporter* brick2Reac = new ReactionReporter("brick2", system, brick2, .01);

    system.updDefaultSubsystem().addEventReporter(dummy3Reac);
    system.updDefaultSubsystem().addEventReporter(dummy2Reac);
    system.updDefaultSubsystem().addEventReporter(dummy1Reac);
    system.updDefaultSubsystem().addEventReporter(brick2Reac);


    const VTKVisualizer& viz = reporter->getVisualizer();
   
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    bushing.setStiffness(state, 1*bushing.getDefaultStiffness());
    bushing.setDamping(state, 1*bushing.getDefaultDamping());

    viz.report(state);
    printf("Default state -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << brick1.getQAsVector(state) << brick2.getQAsVector(state) 
         << " u=" << brick1.getUAsVector(state) << brick2.getUAsVector(state) 
         << "\ndefK=" << bushing.getDefaultStiffness()
         << " k="     << bushing.getStiffness(state)
         << "\ndefC=" << bushing.getDefaultDamping()
         << " c="     << bushing.getDamping(state)
         << endl;
    char ch=getchar();

    state.setTime(0);
    Rotation R30(BodyRotationSequence, 30*(Pi/180), XAxis,
                                       30*(Pi/180), YAxis,
                                       30*(Pi/180), ZAxis);
    brick1.setQToFitTransform(state, R30); 

    Vec3 q30i;
    q30i = Rotation(~R30).convertRotationToBodyFixedXYZ();
    cout << "q30i=" << q30i << endl;

    //dummy1.setQToFitTransform(state, Vec3(-1,-1,-1));
    //dummy2.setOneQ(state, 0, 30*(Pi/180));
    //dummy3.setOneQ(state, 0, 30*(Pi/180));
    //brick2.setOneQ(state, 0, 30*(Pi/180));
    brick2.setQToFitTransform(state, ~R30*Vec3(-1,-1,-1));
    dummy3.setOneQ(state, 0, q30i[0]);
    dummy2.setOneQ(state, 0, q30i[1]);
    dummy1.setOneQ(state, 0, q30i[2]);

    

    RungeKuttaMersonIntegrator integ(system);
    //integ.setMinimumStepSize(1e-1);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    const State& istate = integ.getState();

    viz.report(integ.getState());
    printf("After initialize -- hit ENTER\n");
    cout << "t=" << integ.getTime() 
         << "\nE=" << system.calcEnergy(istate)
         << "\nq=" << brick1.getQAsVector(istate) << brick2.getQAsVector(istate) 
         << "\nu=" << brick1.getUAsVector(istate) << brick2.getUAsVector(istate) 
         << "\nudot=" << brick1.getUDotAsVector(istate) << brick2.getUDotAsVector(istate) 
         << "\ntau=" << brick1.getTauAsVector(istate) << brick2.getTauAsVector(istate) 
         << "\nbrick1=" << brick1.getBodyTransform(istate)
         << "brick2=" << brick2.getBodyTransform(istate)
         << endl;
    bushingReport->handleEvent(istate);
    dummy1Reac->handleEvent(istate);
    dummy2Reac->handleEvent(istate);
    dummy3Reac->handleEvent(istate);
    brick2Reac->handleEvent(istate);
    ch=getchar();
    // Simulate it.
    ts.stepTo(100.0);

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}
