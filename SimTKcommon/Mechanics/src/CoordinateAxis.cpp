/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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


