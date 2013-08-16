/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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
 * Definitions of external global constants of types CoordinateAxis and
 * CoordinateDirection.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/CoordinateAxis.h"

namespace SimTK {

// These constants are global external symbols exported by the library. See
// the CoordinateAxis.h header file for information.
SimTK_SimTKCOMMON_EXPORT const CoordinateAxis::XCoordinateAxis      XAxis;
SimTK_SimTKCOMMON_EXPORT const CoordinateAxis::YCoordinateAxis      YAxis;
SimTK_SimTKCOMMON_EXPORT const CoordinateAxis::ZCoordinateAxis      ZAxis;

SimTK_SimTKCOMMON_EXPORT const CoordinateDirection::NegXDirection   NegXAxis;
SimTK_SimTKCOMMON_EXPORT const CoordinateDirection::NegYDirection   NegYAxis;
SimTK_SimTKCOMMON_EXPORT const CoordinateDirection::NegZDirection   NegZAxis;


}  // End of namespace SimTK


