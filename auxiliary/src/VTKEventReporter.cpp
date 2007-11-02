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

#include "simbody/internal/VTKEventReporter.h"

using namespace SimTK;

class VTKEventReporter::VTKEventReporterRep {
public:
    VTKEventReporterRep(MultibodySystem& system, Real reportInterval, Real defaultScaleForAutoGeometry) : reporter(VTKReporter(system, defaultScaleForAutoGeometry)), reportInterval(reportInterval) {
        lastReportTime = -1;
    }
    VTKReporter& getReporter() {
        return reporter;
    }
    Real getReportInterval() const {
        return reportInterval;
    }
    void setReportInterval(Real interval) {
        reportInterval = interval;
    }
    Real getNextEventTime(const State&) const {
        if (lastReportTime == -1)
            return 0;
        return lastReportTime+reportInterval;
    }
    void handleEvent(const State& state) {
        reporter.report(state);
        lastReportTime = state.getTime();
    }
    VTKEventReporter* handle;
    VTKReporter reporter;
    Real reportInterval;
    Real lastReportTime;
};

VTKEventReporter::VTKEventReporter(MultibodySystem& system, Real reportInterval, Real defaultScaleForAutoGeometry) {
    rep = new VTKEventReporterRep(system, reportInterval, defaultScaleForAutoGeometry);
    updRep().handle = this;
}

VTKEventReporter::~VTKEventReporter() {
    if (rep->handle == this)
        delete rep;
}

/**
 * Get the VTKReporter which generates the images.  It may be used to configure the display.
 */

VTKReporter& VTKEventReporter::getReporter() {
    return updRep().getReporter();
}

/**
 * Get the time interval at which images should be displayed.
 */

Real VTKEventReporter::getReportInterval() const {
    return getRep().getReportInterval();
}

/**
 * Set the time interval at which images should be displayed.
 */

void VTKEventReporter::setReportInterval(Real interval) {
    updRep().setReportInterval(interval);
}

Real VTKEventReporter::getNextEventTime(const State& state) const {
    return getRep().getNextEventTime(state);
}
void VTKEventReporter::handleEvent(const State& state) {
    updRep().handleEvent(state);
}
