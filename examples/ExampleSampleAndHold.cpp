/* -------------------------------------------------------------------------- *
 *                    Simbody(tm) Example: Sample and Hold                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
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

/* This example demonstrates one way to implement sample and hold behavior
in Simbody. The outline is:
    1 define a discrete variable to hold the sampled signal
    2 define a periodic event handler that reads the signal and updates
      the variable
    3 define a system that has some dependence on the variable
Here the system is just a simple pendulum swinging in gravity. We'll 
sample the position of its end point periodically, and stretch a rubber
band line from the sampled location to the current location.
*/

#include "Simbody.h"

using namespace SimTK;

// This is a periodic event handler that interrupts the simulation on a regular
// basis to record the ground-frame location of a given station on a given
// mobilized body and record the result in a given variable.
class StationSampler : public PeriodicEventHandler {
public:
    StationSampler(const MultibodySystem&           system,
                   const MobilizedBody&             mobod, 
                   const Vec3&                      station,
                   const Measure_<Vec3>::Variable&  sample, 
                   Real                             interval) 
    :   PeriodicEventHandler(interval), m_system(system),
        m_mobod(mobod), m_station(station),
        m_sample(sample) {}

    // This method is called whenever this event occurs.
    virtual void handleEvent
       (State& state, Real, bool&) const override 
    {
        m_system.realize(state, Stage::Position);
        Vec3 location = m_mobod.findStationLocationInGround(state,m_station);
        m_sample.setValue(state, location);
    }

private:
    const MultibodySystem&      m_system;
    MobilizedBody               m_mobod;
    Vec3                        m_station;
    Measure_<Vec3>::Variable    m_sample;
};

// The only effect our example sample is going to have on this system is
// to change what's drawn. For that purpose we're defining a
// DecorationGenerator here that we can register with the Visualizer which then
// calls it each time it generates a frame to allow dynamically-generated 
// geometry. Here we'll add a mark at the last-sampled location, and draw a
// line from there to the current location.
class RubberBand : public DecorationGenerator {
public:
    RubberBand(const MultibodySystem&           system,
               const MobilizedBody&             mobod, 
               const Vec3&                      station,
               const Measure_<Vec3>::Variable&  sample) 
    :   m_system(system), m_mobod(mobod), m_station(station),
        m_sample(sample) {}

    virtual void generateDecorations(const State& state, 
                                     Array_<DecorativeGeometry>& geometry) override 
    {
        m_system.realize(state, Stage::Position);
        Vec3 current = m_mobod.findStationLocationInGround(state,m_station);
        Vec3 previous = m_sample.getValue(state);

        // Draw a small black circle at the sampled location.
        geometry.push_back(DecorativeCircle(0.05)
                           .setColor(Black)
                           .setTransform(previous));
        // Draw an orange line from the sampled to the current location.
        geometry.push_back(DecorativeLine(previous,current)
                            .setColor(Orange)
                            .setLineThickness(3));
    }

private:
    const MultibodySystem&      m_system;
    MobilizedBody               m_mobod;
    Vec3                        m_station;
    Measure_<Vec3>::Variable    m_sample;
};

int main() {
  try {
    // Create the system, with subsystems for the bodies and some forces.
    MultibodySystem system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);

    Measure_<Vec3>::Variable prevPos(matter, Stage::Position, Vec3(0));

    // Add gravity as a force element.
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));

    // Create the body and some artwork for it.
    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1).setColor(Red));

    // Add an instance of the body to the multibody system by connecting
    // it to Ground via a pin mobilizer.
    MobilizedBody::Pin pendulum(matter.updGround(), Transform(Vec3(0)), 
                                pendulumBody, Transform(Vec3(0, 1, 0)));

    // Once a second we'll record the pendulum's ground location.
    const Real Interval = 1;
    system.addEventHandler(new StationSampler(system,
                                              pendulum, Vec3(0),
                                              prevPos, Interval));

    // Visualize with default options; ask for a report every 1/30 of a second
    // to match the Visualizer's default 30 frames per second rate.
    Visualizer viz(system);
    viz.addDecorationGenerator(new RubberBand(system,
                                              pendulum, Vec3(0),
                                              prevPos));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    pendulum.setOneU(state, 0, 2.0); // initial velocity 2 rad/sec
    
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
