/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

const int NUM_MOLECULES = 10;
const Real BOND_LENGTH = 0.5;
const Real TEMPERATURE = 300.0;

/**
 * Create a force between every pair of bodies (including ground) with the potential k*(1/x^2 - 1/x).
 */

class BetweenBodyForce : public GeneralForceElements::CustomForce {
public:
    BetweenBodyForce() {
    }
    BetweenBodyForce* clone() const {
        return new BetweenBodyForce(*this);
    }
    void calc(const SimbodyMatterSubsystem& matter, const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const {
        for (int i = 0; i < matter.getNBodies(); ++i) {
            const MobilizedBody& body1 = matter.getMobilizedBody(MobilizedBodyId(i));
            const Vec3 pos1 = body1.getBodyOriginLocation(state);
            for (int j = i+1; j < matter.getNBodies(); ++j) {
                const MobilizedBody& body2 = matter.getMobilizedBody(MobilizedBodyId(j));
                const Vec3 pos2 = body2.getBodyOriginLocation(state);
                const Real dist = (pos2-pos1).norm();
                const Real invDist = 1.0/dist;
                const Real invDist2 = invDist*invDist;
                const Real k = 20.0;
                const Vec3 f = k*(dist-2.0)*(pos2-pos1)*invDist2*invDist2;
                pe -= k*(invDist-1)*invDist;
                body1.applyBodyForce(state, SpatialVec(Vec3(0),  f), bodyForces);
                body2.applyBodyForce(state, SpatialVec(Vec3(0), -f), bodyForces);
            }
        }
    }
};

class EnergyMonitor : public ScheduledEventReporter {
public:
    int eventCount;
    Real sumEnergy, sumEnergySquared;
    EnergyMonitor(MultibodySystem& system) : system(system), lastEventTime(1.0), sumEnergy(0.0), sumEnergySquared(0.0) {
    }
    Real getNextEventTime(const State&, bool includeCurrentTime) const {
        return lastEventTime+0.05;
    }
    void handleEvent(const State& state) {
        eventCount++;
        system.realize(state, Stage::Dynamics);
        Real energy = system.getKineticEnergy(state);
        sumEnergy += energy;
        sumEnergySquared += energy*energy;
        lastEventTime = state.getTime();
    }
private:
    MultibodySystem& system;
    Real lastEventTime;
};

int main() {
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceElements forces(mbs);
    
    // Create a gas of two atom molecules.
    
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)))
                                  .addDecoration(Transform(), DecorativeSphere(.1));
    Random::Uniform random(-10.0, 10.0);
    for (int i = 0; i < NUM_MOLECULES; ++i) {
        Vec3 pos(random.getValue(), random.getValue(), random.getValue());
        MobilizedBody* first = new MobilizedBody::Translation(matter.Ground(), Transform(Vec3(0, 0, 0)), body, Transform(pos));
        MobilizedBody* second = new MobilizedBody::LineOrientation(*first, Transform(Vec3(0, 0, 0)), body, Transform(Vec3(0, 0, BOND_LENGTH)));
    }
    forces.addCustomForce(BetweenBodyForce());
    mbs.updDefaultSubsystem().addEventHandler(new VelocityRescalingThermostat(mbs, SimTK_BOLTZMANN_CONSTANT_MD, TEMPERATURE, 1.0));
    EnergyMonitor* monitor = new EnergyMonitor(mbs);
    mbs.getDefaultSubsystem().addEventReporter(monitor);
    mbs.realizeTopology();
    State& s = mbs.updDefaultState();
    mbs.realizeModel(s);
    
    // Simulate it.
    
    RungeKuttaMersonIntegrator integ(mbs);
    TimeStepper ts(mbs, integ);
    ts.initialize(mbs.getDefaultState());
    ts.stepTo(30);
    
    // See if the temperature was being correctly maintained.
    
    ASSERT(monitor->eventCount > 100);
    Real meanEnergy = monitor->sumEnergy/monitor->eventCount;
    Real expectedEnergy = 5*NUM_MOLECULES*0.5*SimTK_BOLTZMANN_CONSTANT_MD*TEMPERATURE; // kT/2 per degree of freedom
    Real meanSqrEnergy = monitor->sumEnergySquared/monitor->eventCount;
    Real stdDevEnergy = std::sqrt(meanSqrEnergy-meanEnergy*meanEnergy);
    ASSERT(std::abs(meanEnergy/expectedEnergy-1.0) < 0.2);
    ASSERT(stdDevEnergy < 0.5*meanEnergy);
    std::cout << "Done" << std::endl;
}
