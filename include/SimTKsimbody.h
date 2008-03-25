#ifndef SimTK_SIMBODY_SimTKSIMBODY_H_
#define SimTK_SIMBODY_SimTKSIMBODY_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
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

/** @file
 * This is the header file that user code should include to pick up all
 * Simbody capabilities.
 */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
// #include "simbody/internal/MolecularMechanicsSystem.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"
#include "simbody/internal/GeneralForceSubsystem.h"
#include "simbody/internal/HuntCrossleyContact.h"
// #include "simbody/internal/DuMMForceFieldSubsystem.h"
#include "simbody/internal/NumericalMethods.h"
#include "simbody/internal/DecorationSubsystem.h"
#include "simbody/internal/TextDataEventReporter.h"
#include "simbody/internal/ObservedPointFitter.h"
#include "simbody/internal/VelocityRescalingThermostat.h"
#include "simbody/internal/LocalEnergyMinimizer.h"
// #include "simbody/internal/Element.h"
// #include "simbody/internal/Biotype.h"

#endif // SimTK_SIMBODY_SimTKSIMBODY_H_
