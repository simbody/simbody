#ifndef SimTK_SimTKCOMMON_BASICS_H_
#define SimTK_SimTKCOMMON_BASICS_H_

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

/**@file
 * Includes internal headers providing declarations for the basic SimTK
 * Core classes. This file can be included from ANSI C code as well as
 * C++, although only a small subset of the definitions will be done
 * in a C program. These include default precision and simple macros
 * such as those used for physical constants.
 */

/* NOTE: don't use "//" comments in this file! */

/* These two are safe for C programs. */
#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"

#if defined(__cplusplus)
#include "SimTKcommon/internal/Exception.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/ClonePtr.h"
#include "SimTKcommon/internal/ReferencePtr.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Serialize.h"
#include "SimTKcommon/internal/Fortran.h"
#include "SimTKcommon/internal/Array.h"
#include "SimTKcommon/internal/StableArray.h"
#include "SimTKcommon/internal/Value.h"
#include "SimTKcommon/internal/Stage.h"
#include "SimTKcommon/internal/CoordinateAxis.h"
#endif


#endif /* SimTK_SimTKCOMMON_BASICS_H_ */
