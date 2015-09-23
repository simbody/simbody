#ifndef SimTK_SimTKCOMMON_DECORATION_GENERATOR_H_
#define SimTK_SimTKCOMMON_DECORATION_GENERATOR_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon/basics.h"

namespace SimTK {

class State;

/**
 * A DecorationGenerator is used to define geometry that may change over the
 * course of a simulation.  Example include
 *  - Geometry whose position is not fixed relative to any single body.
 *  - Geometry which may appear or disappear during the simulation.
 *  - Geometry whose properties (color, size, etc.) may change during the
 *    simulation.
 *
 * To use it, define a concrete subclass that implements generateDecorations()
 * to generate whatever geometry is appropriate for a given State. It can then
 * be added to a DecorationSubsystem, or directly to a Visualizer.
 */
class DecorationGenerator {
public:
    /**
     * This will be called every time a new State is about to be visualized.
     * It should generate whatever decorations are appropriate for the State
     * and append them to the array.
     */
    virtual void generateDecorations(const State& state,
                                     Array_<DecorativeGeometry>& geometry) = 0;

    /** Destructor is virtual; be sure to override it if you have something
    to clean up at the end. **/
    virtual ~DecorationGenerator() {}
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_DECORATION_GENERATOR_H_
