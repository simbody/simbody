#ifndef SimTK_SIMBODY_SimTKSIMBODY_H_
#define SimTK_SIMBODY_SimTKSIMBODY_H_
/* -------------------------------------------------------------------------- *
 *                                 Simbody(tm)                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Please cite:                                                               *
 *   Michael A. Sherman, Ajay Seth, Scott L. Delp, Simbody: multibody         *
 *   dynamics for biomedical research, Procedia IUTAM 2:241-261 (2011)        *
 *   http://dx.doi.org/10.1016/j.piutam.2011.04.023.                          *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman, Peter Eastman                                    *
 * Contributors: Jack Middleton, Christopher Bruns, Paul Mitiguy, Matthew     *
 *   Millard, Charles Schwieters, Abhinandan Jain, Isaac Newton               *
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
#include "simbody/internal/MobilizedBody_BuiltIns.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/ElasticFoundationForce.h"
#include "simbody/internal/Force.h"
#include "simbody/internal/Force_BuiltIns.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/ForceSubsystemGuts.h"
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
#include "simbody/internal/CableTrackerSubsystem.h"
#include "simbody/internal/CablePath.h"
#include "simbody/internal/CableSpring.h"
#include "simbody/internal/Visualizer.h"
#include "simbody/internal/Visualizer_InputListener.h"
#include "simbody/internal/Visualizer_Reporter.h"
#include "simbody/internal/ConditionalConstraint.h"
#include "simbody/internal/SemiExplicitEulerTimeStepper.h"

#endif // SimTK_SIMBODY_SimTKSIMBODY_H_
