/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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
 * Adhoc main program for playing with contact.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"   // requires VTK

#include "simbody/internal/ContactTrackerSubsystem.h"
#include "simbody/internal/CompliantContactSubsystem.h"

#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
               const CompliantContactSubsystem& compl,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(compl) {}

    ~MyReporter() {}
    void handleEvent(const State& state) const {
        m_system.realize(state, Stage::Velocity);
        cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
             << " Ediss=" << m_compliant.getDissipatedEnergy(state)
             << " E+Ediss=" << m_system.calcEnergy(state)
                               +m_compliant.getDissipatedEnergy(state)
             << endl;
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(2, -9.8, 0));

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);

    contactForces.setTransitionVelocity(1e-3);

    system.updDefaultSubsystem().addEventReporter
        (new MyReporter(system,contactForces,.01));

    // Create a triangle mesh in the shape of a pyramid, with the
    // square base having area 1 (split into two triangles).
    
    Array_<Vec3> vertices;
    vertices.push_back(Vec3(0, 0, 0));
    vertices.push_back(Vec3(1, 0, 0));
    vertices.push_back(Vec3(1, 0, 1));
    vertices.push_back(Vec3(0, 0, 1));
    vertices.push_back(Vec3(0.5, 1, 0.5));
    Array_<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4}, 
                       {2, 1, 4}, {3, 2, 4}, {0, 3, 4}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);

    ContactCliqueId clique1 = ContactSurface::createNewContactClique();
    ContactCliqueId clique2 = ContactSurface::createNewContactClique();

    const Real fFac = 1; // to turn off friction
    const Real fDis = 1; // to turn off dissipation
    matter.Ground().updBody().addDecoration(Vec3(.25+.01,0,0),
        DecorativeBrick(Vec3(.01,2,1)).setColor(Blue));
    matter.Ground().updBody().addContactSurface(Vec3(.25,0,0),
        ContactSurface(ContactGeometry::HalfSpace(),
                       ContactMaterial(100000,fDis*.9,fFac*.8,fFac*.7))
                       .joinClique(clique1));

    Body::Rigid pendulumBody1(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody1.addDecoration(Transform(), DecorativeSphere(0.1));
    pendulumBody1.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Sphere(.1),
                       ContactMaterial(10000,fDis*.9,fFac*.8,fFac*.7))
                       .joinClique(clique2));

    Body::Rigid pendulumBody2(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody2.addDecoration(Transform(), DecorativeSphere(0.1).setColor(Orange));
    pendulumBody2.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Sphere(.1),
                       ContactMaterial(100000,fDis*.9,fFac*.8,fFac*.7))
                       .joinClique(clique2));

    MobilizedBody::Pin pendulum(matter.Ground(), Transform(Vec3(0)), 
                                pendulumBody1,    Transform(Vec3(0, 1, 0)));

    MobilizedBody::Pin pendulum2(pendulum, Transform(Vec3(0)), 
                                 pendulumBody2,    Transform(Vec3(0, 1, 0)));

    Force::MobilityLinearSpring(forces, pendulum2, MobilizerUIndex(0),
        10, 0*(Pi/180));

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


    pendulum.setOneU(state, 0, 2.0);

    
    // Simulate it.

    RungeKuttaMersonIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    integ.setAccuracy(1e-4);
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
