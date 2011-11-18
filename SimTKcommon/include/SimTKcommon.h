#ifndef SimTK_SimTKCOMMON_H_
#define SimTK_SimTKCOMMON_H_

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
 * Core classes, including Simmatrix.
 */

/**@page SimTKcommon    Fundamental SimTK objects and utilities
SimTKcommon is the collection of low-level objects and utility functions that
provide the fundamental capabilities needed by Simbody and useful in
applications that use Simbody. 

@todo Add links to various categories of tools here.
**/

#include "SimTKcommon/basics.h"

#if defined(__cplusplus)
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/MeasureImplementation.h"
#include "SimTKcommon/internal/DecorativeGeometry.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/Study.h"
#include "SimTKcommon/internal/Function.h"
#include "SimTKcommon/internal/Random.h"
#include "SimTKcommon/internal/PolynomialRootFinder.h"
#include "SimTKcommon/internal/Enumeration.h"
#include "SimTKcommon/internal/PrivateImplementation.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include "SimTKcommon/internal/Parallel2DExecutor.h"
#include "SimTKcommon/internal/ParallelWorkQueue.h"
#include "SimTKcommon/internal/ThreadLocal.h"
#include "SimTKcommon/internal/AtomicInteger.h"
#include "SimTKcommon/internal/Plugin.h"
#include "SimTKcommon/internal/Timing.h"
#include "SimTKcommon/internal/Xml.h"
#include "SimTKcommon/Testing.h"
#endif


#endif /* SimTK_SimTKCOMMON_H_ */
