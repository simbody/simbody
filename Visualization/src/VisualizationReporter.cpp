/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/VisualizationReporter.h"

using namespace SimTK;

class VisualizationReporter::VisualizationReporterRep {
public:
    VisualizationReporterRep(MultibodySystem& system) 
    :   system(system), visualizer(system) {}
    VisualizationReporterRep(MultibodySystem& system, const String& title) 
    :   system(system), visualizer(system, title) {}
    const Visualizer& getVisualizer() const {
        return visualizer;
    }
    Visualizer& updVisualizer() {
        return visualizer;
    }
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Acceleration);
        visualizer.report(state);
    }
    MultibodySystem& system;
    VisualizationReporter* handle;
    mutable Visualizer visualizer;
};

VisualizationReporter::VisualizationReporter(MultibodySystem& system, Real reportInterval) : PeriodicEventReporter(reportInterval) {
    rep = new VisualizationReporterRep(system);
    updRep().handle = this;
}

VisualizationReporter::VisualizationReporter(MultibodySystem& system, const String& title, Real reportInterval)
:   PeriodicEventReporter(reportInterval) {
    rep = new VisualizationReporterRep(system, title);
    updRep().handle = this;
}

VisualizationReporter::~VisualizationReporter() {
    if (rep->handle == this)
        delete rep;
}

Visualizer& VisualizationReporter::updVisualizer() {
    return updRep().updVisualizer();
}

const Visualizer& VisualizationReporter::getVisualizer() const {
    return getRep().getVisualizer();
}

void VisualizationReporter::handleEvent(const State& state) const {
    updRep().handleEvent(state);
}
