//-----------------------------------------------------------------------------
// File:     Orientation.h
// Class:    None 
// Parent:   None
// Purpose:  Includes UnitVec3, Quaternion, Rotation, Transform, and related classes
//-----------------------------------------------------------------------------
#ifndef SimTK_SIMMATRIX_ORIENTATION_H_
#define SimTK_SIMMATRIX_ORIENTATION_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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
 *
 * These are numerical utility classes for dealing with the relative orientations
 * of geometric objects. These build on the basic arithmetic classes for small
 * vectors and matrices.
 */

//-----------------------------------------------------------------------------
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/SmallMatrix.h"
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
