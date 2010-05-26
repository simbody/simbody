#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"

#include <cstdio>
#include <exception>
#include <ctime>

using std::cout; using std::cin; using std::endl;

using namespace SimTK;

static const int NBodies = 10;
static const int NBodies2 = 7;

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    const Real mass = 0.1;
    const Vec3 hdims(2,4,1);
    Body::Rigid pendulumBody = Body::Rigid
       (MassProperties(mass, Vec3(0), mass*Gyration::ellipsoid(hdims)))
        .addDecoration(Transform(),DecorativeEllipsoid(hdims)
                                    .setOpacity(.2));

    MobilizedBody midBody, finalBody;
    MobilizedBody::Free baseBody(matter.Ground(), Transform(Vec3(0,-hdims[1],0)), 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));

    Force::TwoPointLinearSpring(forces, matter.Ground(), Vec3(0),
        baseBody, Vec3(0), 2, 10);

    Force::GlobalDamper(forces, matter, 1);


    MobilizedBody parent = baseBody;
    for (int i=0; i < NBodies; ++i) {
        MobilizedBody::Ball child(parent, Transform(Vec3(0,-hdims[1],0)), 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));
        if (i==NBodies/2) midBody   = child;
        if (i==NBodies-1) finalBody = child;
        parent = child;
    }

    // Second chain
    MobilizedBody endOfSecondChain;
    parent = matter.Ground();
    Vec3 loc(-10,0,0);
    for (int i=0; i < NBodies2; ++i) {
        MobilizedBody::Ball child(parent, loc, 
                                  pendulumBody, Transform(Vec3(0,hdims[1],0)));
        if (i==NBodies2-1) endOfSecondChain = child;
        parent = child;
        loc = Vec3(0,-hdims[1],0);
    }

    Force::TwoPointLinearSpring(forces, endOfSecondChain, Vec3(0),
                                midBody, Vec3(0), 2, 1);

    //Constraint::Ball chainConnection(endOfSecondChain, Vec3(0,-hdims[1],0),
    //    matter.updMobilizedBody(MobilizedBodyIndex(4)), Vec3(0));

    //Constraint::Ball ball(matter.Ground(), Vec3(3*NBodies,-NBodies,0), 
    //                      finalBody, Vec3(0,-hdims[1],0));

	VTKEventReporter& vtkReporter = *new VTKEventReporter(system, 0.1);
    VTKVisualizer& viz = vtkReporter.updVisualizer();
    system.updDefaultSubsystem().addEventReporter(&vtkReporter);


    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    const Vec3 midTarget(1.5*NBodies,-NBodies*0.5*hdims[1],3);
    const Vec3 finalTarget(2*NBodies, -NBodies*0.3*hdims[1],2);
    viz.addDecoration(midBody, Vec3(0), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(midBody, Vec3(1,2,3), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(GroundIndex, midTarget, 
        DecorativeSphere(1).setColor(Green));
    viz.addDecoration(finalBody, Vec3(0), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(finalBody, Vec3(1,2,3), DecorativeSphere(1).setColor(Red));
    viz.addDecoration(GroundIndex, finalTarget, 
        DecorativeSphere(1).setColor(Green));

    char c;
    // Show initial configuration
    viz.report(state);
    State tempState = state; 
    matter.setUseEulerAngles(tempState, true);
    system.realizeModel(tempState);
    system.realize(tempState, Stage::Position);
    cout << "INITIAL CONFIGURATION\n"; 
    cout << tempState.getNU() << " dofs, " 
         << tempState.getNQErr() << " constraints.\n";
    
    cin >> c;

    Assembler ik(system);
    Markers& markers = *new Markers();
    ik.adoptAssemblyGoal(&markers);
    //ik.addReporter(vtkReporter);

    ik.lockMobilizer(MobilizedBodyIndex(2));
    ik.lockMobilizer(MobilizedBodyIndex(3));
    ik.lockMobilizer(MobilizedBodyIndex(4));

    ik.restrictQ(midBody, MobilizerQIndex(0), -Pi/10, Pi/10);

    markers.addMarker(midBody, Vec3(0));
    markers.addMarker(midBody, Vec3(1,2,3));
    markers.addMarker(finalBody, Vec3(0));
    markers.addMarker(finalBody, Vec3(1,2,3));

    // Manual observation/marker correspondence.
    //Array_<Markers::MarkerIx> observationOrder;
    //observationOrder.push_back(Markers::MarkerIx(2));
    //observationOrder.push_back(Markers::MarkerIx(3));
    //observationOrder.push_back(Markers::MarkerIx(0));
    //observationOrder.push_back(Markers::MarkerIx(1));
    //observationOrder.push_back(Markers::MarkerIx()); // unused
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //observationOrder.push_back(Markers::MarkerIx());
    //markers.defineTargetOrder(observationOrder);


    const Real Accuracy = 1e-3;
    ik.setAccuracy(Accuracy);
    ik.initialize();

    markers.moveOneObservation(Markers::ObservationIx(0), midTarget);
    markers.moveOneObservation(Markers::ObservationIx(1), midTarget);
    markers.moveOneObservation(Markers::ObservationIx(2), finalTarget);
    markers.moveOneObservation(Markers::ObservationIx(3), finalTarget);

    ik.assemble(state);
    viz.report(state);
    cout << "ASSEMBLED CONFIGURATION\n"; cin >> c;
    State startState = state;

    const clock_t start = clock();

    const int NSteps = 100;
    for (int iters=0; iters <= NSteps; ++iters) {
        Vec3 newMidTarget = midTarget + NBodies/2*std::sin(2*Pi*iters/100);
        Vec3 newFinalTarget = finalTarget - NBodies/3*std::sin(2*Pi*iters/100);
        markers.moveOneObservation(Markers::ObservationIx(0), newMidTarget);
        markers.moveOneObservation(Markers::ObservationIx(1), newMidTarget);
        markers.moveOneObservation(Markers::ObservationIx(2), newFinalTarget);
        markers.moveOneObservation(Markers::ObservationIx(3), newFinalTarget);
                                        
        ik.track();
        ik.updateFromInternalState(state);
        viz.report(state);
    }
    cout << "ASSEMBLED " << NSteps << " steps in " <<
        (double)(clock()-start)/CLOCKS_PER_SEC << "s\n";

    cout << "DONE ASSEMBLING -- SIMULATE ...\n";
    viz.report(state);

    cin >> c;
   
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
