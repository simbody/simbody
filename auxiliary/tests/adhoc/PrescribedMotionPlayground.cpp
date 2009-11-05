#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"

#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));

    MobilizedBody::Pin pendulum(matter.Ground(), Transform(Vec3(0)), 
                                pendulumBody,    Transform(Vec3(0, 1, 0)));
    Motion::Sinusoid(pendulum, Motion::Position, Pi/8, 2*Pi, Pi/4); // amp, rate, phase

    MobilizedBody::Pin pendulum2(pendulum, Transform(Vec3(0)), 
                                 pendulumBody,    Transform(Vec3(0, 1, 0)));
    //Motion::Sinusoid(pendulum2, Motion::Position, Pi/8, 2*Pi, Pi/4); // amp, rate, phase

    Force::MobilityLinearSpring(forces, pendulum2, MobilizerUIndex(0),
        200, 10*(Pi/180));

    VTKEventReporter* reporter = new VTKEventReporter(system, 0.01);
    system.updDefaultSubsystem().addEventReporter(reporter);

    const VTKVisualizer& viz = reporter->getVisualizer();
   
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    viz.report(state);
    printf("Default state -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state) 
         << " u=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state) 
         << endl;
    char c=getchar();

    state.setTime(0);
    system.realize(state, Stage::Time);
    if (matter.prescribe(state, Stage::Position)) cout << "Some PresQ\n";
    else cout << "NO PresQ\n";
    viz.report(state);
    printf("After prescribe(Position) -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state)
         << " u=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state) 
         << endl;
    c=getchar();

    if (matter.prescribe(state, Stage::Velocity)) cout << "Some PresU\n";
    else cout << "NO PresU\n";
    viz.report(state);
    printf("After prescribe(Velocity) -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state) 
         << " u=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state) 
         << endl;
    c=getchar();

    system.realize(state, Stage::Acceleration);
    printf("After realize(Acceleration) -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << "\nq=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state) 
         << "\nu=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state) 
         << "\nudot=" << pendulum.getUDotAsVector(state) << pendulum2.getUDotAsVector(state) 
         << "\ntau=" << pendulum.getTauAsVector(state) << pendulum2.getTauAsVector(state) 
         << endl;
    c=getchar();

    //pendulum.setOneU(state, 0, 1.0);

    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
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
