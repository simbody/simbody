/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKsimbody.h"

// define VISUALIZE for visualisation (for debugging)
//#define VISUALIZE 1

#include <fstream>

using namespace SimTK;
using namespace std;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}


class HarmonicOscillator
{
public:

    class OscillatorReporter : public PeriodicEventReporter 
    {
    public:
        mutable int eventCount;
        mutable Real sumEnergy;
        mutable Real sumEnergySquared;
        mutable Real sumVelocity;
        mutable Real sumAbsVelocity;
        mutable Real sumVelocitySquared;
        mutable Real sumPosition;
        mutable Real sumRMSVelPos;

        OscillatorReporter(HarmonicOscillator& oscillator, Real reportInterval) 
            : PeriodicEventReporter(reportInterval), oscillator(oscillator),
              eventCount(0), sumEnergy(0.0), sumEnergySquared(0.0), 
              sumVelocity(0.0), sumAbsVelocity(0.0), sumVelocitySquared(0.0),
              sumPosition(0.0), sumRMSVelPos(0.0)
              // for Matlab
              //, phaseOut("phaseOut.txt")
        {}

        void handleEvent(const State& state) const 
        {
            // Equilibrate a bit before collecting data
            if (state.getTime() <= 1)
                return;

            ++eventCount;
            oscillator.updSystem().realize(state, Stage::Dynamics);

            //std::ostream& po = phaseOut;
            //po << state.getTime() << " " << state.getQ()[0] << " " << state.getU()[0] << std::endl;

            Real energy = oscillator.getSystem().calcKineticEnergy(state);
            sumEnergy += energy;
            sumEnergySquared += energy*energy;

            Real position = oscillator.getPosition(state);
            sumPosition += position;

            Real velocity = oscillator.getVelocity(state);
            sumVelocity += velocity;
            sumAbsVelocity += std::abs(velocity);
            sumVelocitySquared += velocity*velocity;

            sumRMSVelPos += std::sqrt(position*position * velocity*velocity);

            // oscillator.savedPositions.push_back(oscillator.getPosition(state));
            // oscillator.savedVelocities.push_back(oscillator.getVelocity(state));
        }
    private:
        HarmonicOscillator& oscillator;

        //mutable std::ofstream phaseOut;
    };

    HarmonicOscillator() :
        system(), matter(system), forces(system), mass(1.0)
    {
        Vec3 station(0.0);

        // Use a slider to constrain the oscillator to one dimension of motion
        MobilizedBody::Slider body ( 
            matter.updGround(),
            Body::Rigid(MassProperties(mass, station, Inertia(1))) );
        body.setDefaultLength(-2.0); // initial position at -2
        sliderIndex = body.getMobilizedBodyIndex();

        // Unit spring constant to match example in Frenkel and Smit
        Force::TwoPointLinearSpring(forces, 
            matter.getGround(), Vec3(0),
            body, Vec3(0),
            1.0, 0.0);
    }

    void simulate() {
        // View in Visualizer - for testing only
#ifdef VISUALIZE
        VisualizationReporter* vizrep = new VisualizationReporter(system, 0.2);
        vizrep->updVisualizer().setBackgroundType(Visualizer::SolidColor);
        system.updDefaultSubsystem().addEventReporter(vizrep);
#endif

        reporter = new OscillatorReporter(*this, 0.1);
        system.getDefaultSubsystem().addEventReporter(reporter);

        State state = system.realizeTopology();
        Random::Uniform rand(-1,1);
        //state.updQ()[0] = rand.getValue();
        //state.updU()[0] = rand.getValue();

        // Simulate it.
        
        VerletIntegrator integ(system);
        //RungeKuttaMersonIntegrator integ(system);
        //integ.setAccuracy(0.01);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        ts.stepTo(150.0);
    }

    void assertTemperature(Real temperature) const 
    {
        // ensure we collected some data
        ASSERT(reporter->eventCount > 100);

        int degreesOfFreedom = 1; 

        // Sanity checks

        // Mean position should be zero
        Real expectedMeanPosition = 0.0;
        Real measuredMeanPosition = reporter->sumPosition / reporter->eventCount;
        ASSERT(std::abs(expectedMeanPosition - measuredMeanPosition) < 0.2);

        // Mean velocity should be zero
        Real expectedMeanVelocity = 0.0;
        Real measuredMeanVelocity = reporter->sumVelocity / reporter->eventCount;
        ASSERT(std::abs(expectedMeanVelocity - measuredMeanVelocity) < 0.2);

        // Check temperature
        Real measuredMeanEnergy = reporter->sumEnergy/reporter->eventCount;
        Real expectedMeanEnergy = degreesOfFreedom * 0.5 * SimTK_BOLTZMANN_CONSTANT_MD * temperature; // kT/2 per degree of freedom
        ASSERT(std::abs(1.0 - measuredMeanEnergy/expectedMeanEnergy) < 0.2);

        // Boltzmann distribution stuff

        // Mean squared velocity should be dof*kT/mass
        // Boltzmann distribution
        Real expectedMeanVelocitySquared = 
            degreesOfFreedom * SimTK_BOLTZMANN_CONSTANT_MD * temperature / mass;
        Real measuredMeanVelocitySquared =
            reporter->sumVelocitySquared / reporter->eventCount;
        ASSERT(std::abs(1.0 - measuredMeanVelocitySquared/expectedMeanVelocitySquared) < 0.2);

        // TODO: check this formula
        // Mean absolute velocity should be (8*v2bar/3PI)^1/2
        Real expectedMeanAbsVelocity = std::sqrt(
            8.0 * expectedMeanVelocitySquared / (degreesOfFreedom * SimTK_PI) );
        Real measuredMeanAbsVelocity = 
            reporter->sumAbsVelocity / reporter->eventCount;
       // ASSERT(std::abs(1.0 - measuredMeanAbsVelocity/expectedMeanAbsVelocity) < 0.2);
    }

    MultibodySystem& updSystem() {return system;}
    const MultibodySystem& getSystem() const {return system;}

    SimbodyMatterSubsystem& updMatterSubsystem() {return matter;}
    const SimbodyMatterSubsystem& getMatterSubsystem() const {return matter;}

    GeneralForceSubsystem& updForceSubsystem() {return forces;}
    const GeneralForceSubsystem& getForceSubsystem() const {return forces;}

    // slider coordinate
    Real getPosition(const State& state) const {
        const MobilizedBody::Slider& slider = MobilizedBody::Slider::downcast( matter.getMobilizedBody(sliderIndex) );

        return slider.getLength(state);
    }

    Real getVelocity(const State& state) const {
        const MobilizedBody::Slider& slider = MobilizedBody::Slider::downcast( matter.getMobilizedBody(sliderIndex) );

        return slider.getRate(state);
    }

    Real getTime(const State& state) const {
        return state.getTime();
    }

    // std::vector<Real> savedVelocities;
    // std::vector<Real> savedPositions;

private:
    MultibodySystem system;
    SimbodyMatterSubsystem matter;
    GeneralForceSubsystem forces;

    Real mass;
    MobilizedBodyIndex sliderIndex;

    OscillatorReporter* reporter;
};


// Case study 12, page 155 in
// Understanding Molecular Simulation: From Algorithms to Applications
// Frenkel and Smit
void testHarmonicOscillatorNoThermostat() 
{
    HarmonicOscillator oscillator;
    oscillator.simulate();
}

void testNoseHooverConstructorSmoke()
{
    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat(forces, oscillator.getMatterSubsystem(), 
        SimTK_BOLTZMANN_CONSTANT_MD, 300, .1);
    oscillator.simulate();
}

void testOscillatorTemperature(Real temperature, int nChains=-1)
{
    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat nhc(forces, oscillator.getMatterSubsystem(), 
        SimTK_BOLTZMANN_CONSTANT_MD, temperature, 0.1);
    if (nChains > 0)
        nhc.setDefaultNumChains(nChains);
    oscillator.simulate();
    oscillator.assertTemperature(temperature);
}


int main() 
{

    // Several tests commented out to lessen burden on nightly build tests

    //cout << "testHarmonicOscillatorNoThermostat" << endl;
    //testHarmonicOscillatorNoThermostat();

    //cout << "testNoseHooverConstructorSmoke" << endl;
    //testNoseHooverConstructorSmoke();

    cout << "oscillator 100K" << endl;
    testOscillatorTemperature(100); // use default #chains

    //cout << "oscillator 300K" << endl;
    //testOscillatorTemperature(300.0);

    //cout << "oscillator 5000K" << endl;
    //testOscillatorTemperature(5000.0);

    //cout << "argon 10K" << endl;
    //testArgonTemperature(10.0);

    //cout << "argon 300K" << endl;
    //testArgonTemperature(300.0);

    //cout << "argon 5000K" << endl;
    //testArgonTemperature(5000.0);

    return 0;
}
