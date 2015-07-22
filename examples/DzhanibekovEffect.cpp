/* -------------------------------------------------------------------------- *
 *                      Simbody(tm) Example: Dzhanibekov Effect               *
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

/*                      Simbody Dzhanibekov Effect
This example demonstrates a non-intuitive behavior of a freely rotating
rigid body, known as the Dzhanibekov Effect. This is best seen in zero
gravity on the space station. Here are some cool real-world videos:
    https://www.youtube.com/watch?v=L2o9eBl_Gzw
    https://www.youtube.com/watch?v=JB0OAt4zQ1E

    And here is one generated from this example:
    https://youtu.be/L8B83DUKiiA
*/

#include "Simbody.h"
#include <iostream>

using namespace SimTK;

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


//==============================================================================
//                                  MAIN
//==============================================================================
int main() {
  try {
    // Create the system.
    MultibodySystem system; system.setUpDirection(ZAxis);
    SimbodyMatterSubsystem matter(system);
    // No gravity or other forces

    matter.setShowDefaultGeometry(false); // turn off frames and other junk

    // Construct a single rigid body by welding together a cylindrical shaft
    // and a rectangular bar.
    Rotation YtoX(-Pi/2, ZAxis);
    Body::Rigid shaftBody(MassProperties(1, Vec3(0),
                            UnitInertia::cylinderAlongX(.02, .05)));
    shaftBody.addDecoration(YtoX,
                            DecorativeCylinder(.02, .05).setColor(Red));

    const Vec3 halfLengths(.02,.04,.3);
    Body::Rigid barBody(MassProperties(2, Vec3(0),
                    UnitInertia::brick(halfLengths)));
    barBody.addDecoration(Transform(),
                          DecorativeBrick(halfLengths).setColor(Blue));

    MobilizedBody::Free shaft(matter.Ground(), Transform(),
                              shaftBody, Transform());
    MobilizedBody::Weld bar(shaft, Vec3(-.05,0,0),
                            barBody, Vec3(halfLengths[0],0,0));

    // Visualize a frame every 1/60 s, and include the energy.
    Visualizer viz(system); viz.setDesiredFrameRate(60);
    viz.addDecorationGenerator(new ShowEnergy(system));
    system.addEventReporter(new Visualizer::Reporter(viz, 1./60));

    // Initialize the system and state.
    State state = system.realizeTopology();

    // Set initial conditions. Need a slight perturbation of angular velocity
    // to trigger the instability.
    shaft.setQToFitTranslation(state, Vec3(0,0,.5));
    shaft.setUToFitAngularVelocity(state, Vec3(10,0,1e-10)); // 10 rad/s

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-5);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(100.0);

  } catch(const std::exception& e) {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    return 1;
  }
    return 0;
}


void ShowEnergy::generateDecorations(const State&                state,
                                     Array_<DecorativeGeometry>& geometry)
{
    m_mbs.realize(state, Stage::Dynamics);
    const Real E=m_mbs.calcEnergy(state);
    DecorativeText energy;
    energy.setTransform(Vec3(-.2,0,.5))
            .setText("Energy: " + String(E, "%.6f"))
            .setScale(.09)
            .setColor(Black);
    geometry.push_back(energy);
}
