//-----------------------------------------------------------------------------
// File:     Orientation.h
// Class:    None 
// Parent:   None
// Purpose:  Includes UnitVec3, Quaternion, Rotation, Transform, and related classes
//-----------------------------------------------------------------------------
#ifndef SimTK_SIMMATRIX_ORIENTATION_H_
#define SimTK_SIMMATRIX_ORIENTATION_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman and Paul Mitiguy                                  *
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
 *
 * These are numerical utility classes for dealing with the relative orientations
 * of geometric objects. These build on the basic arithmetic classes for small
 * vectors and matrices.
 */

//-----------------------------------------------------------------------------
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/UnitVec.h"
#include "SimTKcommon/internal/Quaternion.h"
#include "SimTKcommon/internal/Rotation.h"
#include "SimTKcommon/internal/Transform.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Some handy conversion constants. 
// Use templatized inline conversion routines instead whenever possible.
// These are defined so that you can multiply by them. For example, if you have
// an angle qRad in radians and want to convert to degrees, write qDeg = qRad*SimTK_RTD.
// Note that the expressions here will always be evaluated at compile time, yielding
// long double results which you can cast to smaller sizes if you want.
//-----------------------------------------------------------------------------
#define SimTK_RTD (180/SimTK_PI)
#define SimTK_DTR (SimTK_PI/180)

//-----------------------------------------------------------------------------
namespace SimTK {

inline static Real  convertRadiansToDegrees(const Real rad) { return rad*Real(SimTK_RTD); }
inline static Real  convertDegreesToRadians(const Real deg) { return deg*Real(SimTK_DTR); }


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_SIMMATRIX_ORIENTATION_H_
//--------------------------------------------------------------------------
