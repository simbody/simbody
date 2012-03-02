/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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