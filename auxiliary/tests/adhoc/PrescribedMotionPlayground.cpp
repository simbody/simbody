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
 * Adhoc main program for playing with prescribed motion.
 */

#include "SimTKsimbody.h"

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
        100, 0*(Pi/180));

    Visualizer viz(system);

   
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

    system.updDefaultSubsystem().addEventReporter(new Visualizer::Reporter(viz, 0.01));

    RungeKuttaMersonIntegrator integ(system);
    //integ.setMinimumStepSize(1e-1);
    integ.setAccuracy(1e-3);
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
