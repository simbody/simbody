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

/** @page Simbody          Simbody Reference Guide
 *
 * This is the Doxygen-generated reference material for Simbody objects. Until
 * we finish this main page you'll have to look up the objects by their class
 * names.
 *
 * @section main_intro Introduction
 * Simbody is an API for building multibody systems, whose main components
 * are:
 *  - @ref main_mobod
 *  - @ref main_forces
 *  - @ref main_constraints
 *  - @ref main_motions
 *
 * @section main_mobod Mobilized Bodies
 * A MobilizedBody is a rigid body plus its mobilizer, which is the joint that 
 * defines the body's mobility with respect to its parent body in the multibody
 * tree. Here are some of the built-in mobilizer types:
 * - @link SimTK::MobilizedBody::Pin    Pin (torsion, revolute) @endlink
 * - @link SimTK::MobilizedBody::Ball   Ball (spherical) @endlink
 * - @link SimTK::MobilizedBody::Slider Slider (prismatic) @endlink
 * - @link SimTK::MobilizedBody::Free   Free (6 dof) @endlink
 * @section main_forces Forces
 * tbd
 * @section main_constraints Constraints
 * tbd
 * @section main_motions Motions
 * tbd
 */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/Motion.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Contact.h"
#include "simbody/internal/ContactGeometry.h"
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
#include "simbody/internal/OrientedBoundingBox.h"
#include "simbody/internal/LocalEnergyMinimizer.h"
#include "simbody/internal/ContactTrackerSubsystem.h"
#include "simbody/internal/CompliantContactSubsystem.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/VisualizationEventListener.h"
#include "simbody/internal/VisualizationReporter.h"

#endif // SimTK_SIMBODY_SimTKSIMBODY_H_
