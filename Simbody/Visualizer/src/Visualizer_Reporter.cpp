/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Visualizer_Reporter.h"

using namespace SimTK;

class Visualizer::Reporter::Impl {
public:
    explicit Impl(const Visualizer& viz) 
    :   handle(0), visualizer(viz) {}
    explicit Impl(const MultibodySystem& system) 
    :   handle(0), visualizer(system) {}

    const Visualizer& getVisualizer() const {
        return visualizer;
    }

    void handleEvent(const State& state) const {
        visualizer.getSystem().realize(state, Stage::Acceleration);
        visualizer.report(state);
    }

    Visualizer::Reporter*   handle;
    Visualizer              visualizer; // shallow copy
};

Visualizer::Reporter::Reporter(const Visualizer& viz, Real reportInterval) 
:   PeriodicEventReporter(reportInterval) {
    impl = new Impl(viz);
    updImpl().handle = this;
}

Visualizer::Reporter::Reporter(const MultibodySystem& sys, Real reportInterval) 
:   PeriodicEventReporter(reportInterval) {
    impl = new Impl(sys);
    updImpl().handle = this;
}

Visualizer::Reporter::~Reporter() {
    if (impl->handle == this)
        delete impl;
}

const Visualizer& Visualizer::Reporter::getVisualizer() const {
    return getImpl().getVisualizer();
}

void Visualizer::Reporter::handleEvent(const State& state) const {
    getImpl().handleEvent(state);
}
