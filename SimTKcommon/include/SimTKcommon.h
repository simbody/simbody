#ifndef SimTK_SimTKCOMMON_H_
#define SimTK_SimTKCOMMON_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
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
 * Core classes, including Simmatrix.
 */

#include "SimTKcommon/basics.h"

#if defined(__cplusplus)
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/MeasureImplementation.h"
#include "SimTKcommon/internal/PolygonalMesh.h"
#include "SimTKcommon/internal/DecorativeGeometry.h"
#include "SimTKcommon/internal/DecorationGenerator.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/SystemGuts.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/Study.h"
#include "SimTKcommon/internal/Function.h"
#include "SimTKcommon/internal/Random.h"
#include "SimTKcommon/internal/PolynomialRootFinder.h"
#include "SimTKcommon/internal/PrivateImplementation.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"
#include "SimTKcommon/internal/ParallelExecutor.h"
#include "SimTKcommon/internal/Parallel2DExecutor.h"
#include "SimTKcommon/internal/ParallelWorkQueue.h"
#include "SimTKcommon/internal/Pathname.h"
#include "SimTKcommon/internal/Plugin.h"
#include "SimTKcommon/internal/Timing.h"
#include "SimTKcommon/internal/Xml.h"
#include "SimTKcommon/Testing.h"
#endif


#endif /* SimTK_SimTKCOMMON_H_ */
