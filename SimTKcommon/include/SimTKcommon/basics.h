#ifndef SimTK_SimTKCOMMON_BASICS_H_
#define SimTK_SimTKCOMMON_BASICS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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
