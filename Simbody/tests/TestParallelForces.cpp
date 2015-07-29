/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-15 Stanford University and the Authors.        *
 * Authors: Thomas Lau                                                    *
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

#include "SimTKsimbody.h"
#include <time.h>
#include <chrono>
#include <thread>
#include <iostream>

using namespace SimTK;
using namespace std;

class ParallelForceImpl : public Force::Custom::Implementation {
public:
    mutable bool hasRealized[Stage::Report+1];
    ParallelForceImpl() {
        for (int i = 0; i < Stage::NValid; i++)
            hasRealized[i] = false;
    }
    bool shouldBeParallelIfPossible() const override{
      return true;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
         std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    Real calcPotentialEnergy(const State& state) const {
        return 0.0;
    }
    void realizeTopology(State& state) const {
        hasRealized[Stage::Topology] = true;
    }
    void realizeModel(State& state) const {
        hasRealized[Stage::Model] = true;
    }
    void realizeInstance(const State& state) const {
        hasRealized[Stage::Instance] = true;
    }
    void realizeTime(const State& state) const {
        hasRealized[Stage::Time] = true;
    }
    void realizePosition(const State& state) const {
        hasRealized[Stage::Position] = true;
    }
    void realizeVelocity(const State& state) const {
        hasRealized[Stage::Velocity] = true;
    }
    void realizeDynamics(const State& state) const {
        hasRealized[Stage::Dynamics] = true;
    }
    void realizeAcceleration(const State& state) const {
        hasRealized[Stage::Acceleration] = true;
    }
    void realizeReport(const State& state) const {
        hasRealized[Stage::Report] = true;
    }
};

class NonParallelForceImpl : public Force::Custom::Implementation {
public:
    mutable bool hasRealized[Stage::Report+1];
    NonParallelForceImpl() {
        for (int i = 0; i < Stage::NValid; i++)
            hasRealized[i] = false;
    }
    bool shouldBeParallelIfPossible() const override{
      return false;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const {
         std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    Real calcPotentialEnergy(const State& state) const {
        return 0.0;
    }
    void realizeTopology(State& state) const {
        hasRealized[Stage::Topology] = true;
    }
    void realizeModel(State& state) const {
        hasRealized[Stage::Model] = true;
    }
    void realizeInstance(const State& state) const {
        hasRealized[Stage::Instance] = true;
    }
    void realizeTime(const State& state) const {
        hasRealized[Stage::Time] = true;
    }
    void realizePosition(const State& state) const {
        hasRealized[Stage::Position] = true;
    }
    void realizeVelocity(const State& state) const {
        hasRealized[Stage::Velocity] = true;
    }
    void realizeDynamics(const State& state) const {
        hasRealized[Stage::Dynamics] = true;
    }
    void realizeAcceleration(const State& state) const {
        hasRealized[Stage::Acceleration] = true;
    }
    void realizeReport(const State& state) const {
        hasRealized[Stage::Report] = true;
    }
};

void testParallelForce()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    
    for(int x = 0; x < 50; x++)
      Force::Custom custom(forces, new ParallelForceImpl());
    
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realize(state, Stage::Dynamics);
}

void testNonParallelForce()
{
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    
    for(int x = 0; x < 50; x++)
      Force::Custom custom(forces, new NonParallelForceImpl());
    
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realize(state, Stage::Dynamics);
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    testParallelForce();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedParallel = end-start;
    
    start = std::chrono::high_resolution_clock::now();
    testNonParallelForce();
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedNonParallel = end-start;
    
    cout << "Parallel Forces Time: " << elapsedParallel.count() << endl;
    cout << "NonParallel Forces Time: " << elapsedNonParallel.count() << endl;
    
    assert(elapsedNonParallel.count() < elapsedParallel.count());
    return 0;
}