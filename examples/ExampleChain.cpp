#include "SimTKsimbody.h"

using namespace SimTK;

int main() {
    
    // Create the system.
    
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));
    MobilizedBodyIndex lastBody = matter.getGround().getMobilizedBodyIndex();
    for (int i = 0; i < 10; ++i) {
        MobilizedBody::Ball pendulum(matter.updMobilizedBody(lastBody), Transform(Vec3(0)), pendulumBody, Transform(Vec3(0, 1, 0)));
        lastBody = pendulum.getMobilizedBodyIndex();
    }

    Visualizer viz(system);
    system.updDefaultSubsystem().addEventReporter
        (new Visualizer::Reporter(viz, 0.02));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Gaussian random;
    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue();
    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(100.0);
}
