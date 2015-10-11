/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
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

/** @file
Non-inline static methods from the Geo class. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_LineSeg.h"
#include "simmath/internal/Geo_Box.h"
#include "simmath/internal/Geo_CubicHermiteCurve.h"
#include "simmath/internal/Geo_BicubicHermitePatch.h"
#include "simmath/internal/Geo_CubicBezierCurve.h"
#include "simmath/internal/Geo_BicubicBezierPatch.h"

namespace SimTK {

//==============================================================================
//                                   GEO
//==============================================================================

// Template instantiations for subclasses that don't have their own source
// files.

template class Geo::LineSeg_<float>;
template class Geo::LineSeg_<double>;

template class Geo::CubicHermiteCurve_<float>;
template class Geo::CubicHermiteCurve_<double>;

template class Geo::BicubicHermitePatch_<float>;
template class Geo::BicubicHermitePatch_<double>;

template class Geo::CubicBezierCurve_<float>;
template class Geo::CubicBezierCurve_<double>;

template class Geo::BicubicBezierPatch_<float>;
template class Geo::BicubicBezierPatch_<double>;


}  // End of namespace SimTK
