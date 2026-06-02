#include <Simbody.h>
using namespace SimTK;

int main() {
  try {
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));

    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));

    MobilizedBody::Ball spinner(matter.Ground(), Transform(Vec3(0,0,0)), pendulumBody, Transform(Vec3(-1, 0, 0)));

    // Set up visualization.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30.));
   
    // Initialize the system and state.
    system.realizeTopology();
    State state = system.getDefaultState();
    double a = 20.; // degrees from the vertical
    double v = 10.; // generalized axial velocity
    a = (90.-a)*3.14159265/180.;
    spinner.setQ(state, Vec4(cos(a/2), 0, 0, sin(a/2))); // rotate the ball mobilizer to the initial position
    spinner.setU(state, Vec3(v*cos(a),v*sin(a),0));      // apply initial spinning velocity
   
    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(50.0);
  } catch (const std::exception& e) {
    std::cout << "Error: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
