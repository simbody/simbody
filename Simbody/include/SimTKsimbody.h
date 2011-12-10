#ifndef SimTK_SIMBODY_SimTKSIMBODY_H_
#define SimTK_SIMBODY_SimTKSIMBODY_H_

/* -------------------------------------------------------------------------- *
 *                             SimTK: Simbody(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-11 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
 * Contributors: Jack Middleton, Christopher Bruns, Paul Mitiguy,             *
 *               Charles Schwieters, Isaac Newton                             *
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
This header file includes all the Simbody header files that need to be 
visible to a compiler processing a Simbody-using compilation unit.\ However,
user programs should included only the top-level Simbody.h header (which 
will include this one). **/

// This should be kept self-contained for backwards compatibility since
// in releases prior to Simbody 2.2 users were told to include "SimTKsimbody.h"
// rather than the now-preferred "Simbody.h".

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/Motion.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/CollisionDetectionAlgorithm.h"
#include "simbody/internal/ElasticFoundationForce.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/Force_Gravity.h"
#include "simbody/internal/Force_LinearBushing.h"
#include "simbody/internal/Force_Thermostat.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"
#include "simbody/internal/GeneralContactSubsystem.h"
#include "simbody/internal/GeneralForceSubsystem.h"
#include "simbody/internal/HuntCrossleyContact.h"
#include "simbody/internal/HuntCrossleyForce.h"
#include "simbody/internal/DecorationSubsystem.h"
#include "simbody/internal/TextDataEventReporter.h"
#include "simbody/internal/ObservedPointFitter.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/LocalEnergyMinimizer.h"
#include "simbody/internal/ContactTrackerSubsystem.h"
#include "simbody/internal/CompliantContactSubsystem.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/Visualizer_InputListener.h"
#include "simbody/internal/Visualizer_Reporter.h"

#endif // SimTK_SIMBODY_SimTKSIMBODY_H_
