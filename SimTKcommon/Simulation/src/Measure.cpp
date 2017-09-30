/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-13 Stanford University and the Authors.        *
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


/**@file
 *
 * Implementation of non-inline methods from the Measure family of classes.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/MeasureImplementation.h"

#include <cassert>
#include <algorithm>

#if defined(__clang__)
    #if __has_warning("-Winstantiation-after-specialization")
        // Avoid the harmless warning "explicit instantiation of 'Zero' that
        // occurs after an explicit specialization has no effect"
        #pragma clang diagnostic ignored "-Winstantiation-after-specialization"
    #endif
#endif

namespace SimTK {

// These are here just to make sure they compile.
template class Measure_<Real>::Constant;
template class Measure_<Real>::Zero;
template class Measure_<Real>::One;
template class Measure_<Real>::Differentiate;
template class Measure_<Real>::Integrate;
template class Measure_<Real>::Result;
template class Measure_<Real>::Variable;
template class Measure_<Real>::Extreme;
template class Measure_<Real>::Delay;

template class Measure_<Vec3>::Constant;
template class Measure_<Vec3>::Zero;
template class Measure_<Vec3>::One;
template class Measure_<Vec3>::Differentiate;
template class Measure_<Vec3>::Integrate;
template class Measure_<Vec3>::Result;
template class Measure_<Vec3>::Variable;
template class Measure_<Vec3>::Extreme;
template class Measure_<Vec3>::Delay;

template class Measure_<Vector>::Constant;
template class Measure_<Vector>::Zero;
template class Measure_<Vector>::One;
template class Measure_<Vector>::Differentiate;
template class Measure_<Vector>::Integrate;
template class Measure_<Vector>::Result;
template class Measure_<Vector>::Variable;
template class Measure_<Vector>::Extreme;
template class Measure_<Vector>::Delay;

} // namespace SimTK

