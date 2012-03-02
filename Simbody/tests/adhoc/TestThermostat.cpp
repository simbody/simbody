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
 * Test the Force::Thermostat force element.
 */

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"


#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

class ThermoReporter : public PeriodicEventReporter {
public:
    ThermoReporter(const MultibodySystem& sys,
                    const Force::Thermostat& thermo,
                    const Force::LinearBushing& bushing1,
                    const Force::LinearBushing& bushing2,
                    Real dt) 
    :   PeriodicEventReporter(dt), system(sys), thermo(thermo),
        bushing1(bushing1), bushing2(bushing2) {}

    void handleEvent(const State& state) const {
        printf("THERMO t=%g, stage %s, KE+PE=%g Ebath=%g CONSERVED=%g\n", 
                state.getTime(),
                state.getSystemStage().getName().c_str(),
                system.calcEnergy(state),
                thermo.calcBathEnergy(state),
                system.calcEnergy(state)
                    + thermo.calcBathEnergy(state)
                    + bushing1.getDissipatedEnergy(state)
                    + bushing2.getDissipatedEnergy(state));
        printf("TEMP=%g N=%d Power=%g Work=%g\n", 
            thermo.getCurrentTemperature(state),
            thermo.getNumThermalDofs(state),
            thermo.getExternalPower(state),
            thermo.getExternalWork(state));
    }
private:
    const MultibodySystem&      system;
    const Force::Thermostat&    thermo;
    const Force::LinearBushing& bushing1;
    const Force::LinearBushing& bushing2;
};

void testConservationOfEnergy() {
    // Create the system.
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    const Real Mass = 1;
    const Vec3 HalfShape = Vec3(1,.5,.25)/2;
    const Transform BodyAttach(Rotation(), Vec3(HalfShape[0],0,0));
    Body::Rigid brickBody(MassProperties(Mass, Vec3(.1,.2,.3), 
                                Mass*Inertia(1,1.1,1.2,0.01,0.02,0.03)));
    //Body::Rigid brickBody(MassProperties(Mass, Vec3(0), 
    //                        Mass*UnitInertia::ellipsoid(HalfShape)));
    brickBody.addDecoration(Transform(), DecorativeEllipsoid(HalfShape)
                                            .setOpacity(0.25)
                                            .setColor(Blue));
    brickBody.addDecoration(BodyAttach,
                DecorativeFrame(0.5).setColor(Red));

    const int NBod=50;
    MobilizedBody::Free brick1(matter.Ground(), Transform(), 
                               brickBody,       BodyAttach);
    MobilizedBody::Free brick2(brick1, Transform(), 
                               brickBody,       BodyAttach);
    MobilizedBody prev=brick2;
    MobilizedBody body25;
    for (int i=0; i<NBod; ++i) {
        MobilizedBody::Ball next(prev, -1*BodyAttach.p(), 
                                 brickBody, BodyAttach);
        if (i==25) body25=next;
        //Force::TwoPointLinearSpring(forces,
        //    prev, Vec3(0), next, Vec3(0), 1000, 1);
        prev=next;
    }

    Constraint::Ball(matter.Ground(), Vec3(0,1,0)-2*NBod/3*BodyAttach.p(),
        prev, BodyAttach.p());
    Constraint::Ball(matter.Ground(), Vec3(0,1,0)-1*NBod/3*BodyAttach.p(),
        body25, BodyAttach.p());


    Vec6 k1(1,100,1,100,100,100), c1(0);
    Force::LinearBushing(forces, matter.Ground(), -2*NBod/3*BodyAttach.p(), 
                 prev, BodyAttach.p(), k1, c1);
    matter.Ground().addBodyDecoration(-2*NBod/3*BodyAttach.p(),
        DecorativeFrame().setColor(Green));

    Force::Thermostat thermo(forces, matter,
        SimTK_BOLTZMANN_CONSTANT_MD,
        5000,
        .1);

    Vec6 k(1,100,1,100,100,100), c(0);
    Force::LinearBushing bushing1(forces, matter.Ground(), -1*NBod/3*BodyAttach.p(),
        brick1, BodyAttach, k, c);
    Force::LinearBushing bushing2(forces, brick1, Transform(),
        brick2, BodyAttach, k, c);
    matter.Ground().addBodyDecoration(-1*NBod/3*BodyAttach.p(),
        DecorativeFrame().setColor(Green));
 

    Visualizer viz(system);
    Visualizer::Reporter* reporter = new Visualizer::Reporter(viz, 1./30);
    viz.setBackgroundType(Visualizer::SolidColor);
    system.addEventReporter(reporter);

    ThermoReporter* thermoReport = new ThermoReporter
        (system, thermo, bushing1, bushing2, 1./10);
    system.addEventReporter(thermoReport);
   
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    viz.report(state);
    printf("Default state -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << brick1.getQAsVector(state) << brick2.getQAsVector(state) 
         << " u=" << brick1.getUAsVector(state) << brick2.getUAsVector(state) 
         << "\nnChains=" << thermo.getNumChains(state)
         << " T="        << thermo.getBathTemperature(state)
         << "\nt_relax=" << thermo.getRelaxationTime(state)
         << " kB="       << thermo.getBoltzmannsConstant()
         << endl;
    getchar();

    state.setTime(0);
    system.realize(state, Stage::Acceleration);
    Vector initU(state.getNU());
    initU = Test::randVector(state.getNU());
    state.updU()=initU;
    

    RungeKuttaMersonIntegrator integ(system);
    //integ.setMinimumStepSize(1e-1);
    integ.setAccuracy(1e-2);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    const State& istate = integ.getState();

    viz.report(integ.getState());
    viz.zoomCameraToShowAllGeometry();
    printf("After initialize -- hit ENTER\n");
    cout << "t=" << integ.getTime() 
         << "\nE=" << system.calcEnergy(istate)
         << "\nEbath=" << thermo.calcBathEnergy(istate)
         << endl;
    thermoReport->handleEvent(istate);
    getchar();

    // Simulate it.
    ts.stepTo(20.0);

    viz.report(integ.getState());
    viz.zoomCameraToShowAllGeometry();
    printf("After simulation:\n");
    cout << "t=" << integ.getTime() 
         << "\nE=" << system.calcEnergy(istate)
         << "\nEbath=" << thermo.calcBathEnergy(istate)
         << "\nNsteps=" << integ.getNumStepsTaken()
         << endl;
    thermoReport->handleEvent(istate);
}

int main() {
    SimTK_START_TEST("TestThermostat");
        SimTK_SUBTEST(testConservationOfEnergy);
    SimTK_END_TEST();
}
