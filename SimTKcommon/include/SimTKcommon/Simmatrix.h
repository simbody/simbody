#ifndef SimTK_SIMMATRIX_H_
#define SimTK_SIMMATRIX_H_

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
 * This is the header which should be included in user programs that would
 * like to make use of all the Simmatrix facilities, but none of the other
 * parts of SimTKcommon.
 */

// Each of these is independently user-includable, with later ones including
// former ones.
#include "SimTKcommon/Scalar.h"         // self-contained
#include "SimTKcommon/SmallMatrix.h"    // includes Scalar.h
#include "SimTKcommon/Orientation.h"    // includes SmallMatrix.h
#include "SimTKcommon/Mechanics.h"      // includes Orientation.h

// Here we add the missing pieces that provide large matrix functionality,
// and some additional small matrix functionality that depends on having
// access to large matrix capabilities.
#include "SimTKcommon/internal/BigMatrix.h"
#include "SimTKcommon/internal/SmallDefsThatNeedBig.h"
#include "SimTKcommon/internal/VectorMath.h"

#endif // SimTK_SIMMATRIX_H_
