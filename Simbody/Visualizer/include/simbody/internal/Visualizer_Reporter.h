#ifndef SimTK_SIMBODY_VISUALIZER_REPORTER_H_
#define SimTK_SIMBODY_VISUALIZER_REPORTER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Visualizer.h"

namespace SimTK {

class MultibodySystem;

/** This is an EventReporter that makes it easy to generate on-screen movies of
any simulation. Use it like this:
@code
    MultibodySystem system;
    // ... build your system

    // Create a Visualizer object for communication with the visualizer.
    Visualizer viz(system);
    // ... set visualization options by calling methods on viz

    // Create a Reporter that will make periodic calls to the Visualizer's
    // report() method to render frames. Note that ownership of the Reporter
    // is taken by the System; don't delete it yourself.
    system.addEventReporter(new Visualizer::Reporter(viz, interval));
@endcode
**/
class SimTK_SIMBODY_EXPORT Visualizer::Reporter : public PeriodicEventReporter {
public:
    /** Create a Reporter for the given Visualizer \a viz, and call its
    report() method every \a reportInterval time units of \e simulation time
    (not necessarily measured in seconds). Note that if you want to run your
    simulation in real time and you aren't using seconds as time units, you
    should set the time scale via the Visualizer's setRealTimeScale() method
    and set the report interval here to TimeScale/FrameRate. **/
    explicit Reporter(const Visualizer& viz, Real reportInterval=Infinity);

    /** This constructor will create a Visualizer with all the default
    settings for the supplied system \a sys. This is an abbreviation for
    @code Reporter(Visualizer(system), reportInterval); @endcode. **/
    explicit Reporter(const MultibodySystem& sys, Real reportInterval=Infinity);

    /** Destructor will also destroy the contained Visualizer object if there
    are no other references to it. **/
    ~Reporter();

    /** Get the Visualizer which this Reporter is using to generate images. **/
    const Visualizer& getVisualizer() const;

    /** This satisfies the pure virtual method in EventReporter. **/
    virtual void handleEvent(const State& state) const;

protected:
    class Impl;
    Impl* impl;
    const Impl& getImpl() const {assert(impl); return *impl;}
    Impl&       updImpl()       {assert(impl); return *impl;}
};

/** OBSOLETE: This provides limited backwards compatibility with the old VTK
Visualizer that is no longer supported.\ Switch to Visualizer::Reporter instead. **/
typedef Visualizer::Reporter VTKEventReporter;

} // namespace SimTK

#endif // SimTK_SIMBODY_VISUALIZER_REPORTER_H_
