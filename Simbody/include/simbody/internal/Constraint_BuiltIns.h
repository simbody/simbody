#ifndef SimTK_SIMBODY_CONSTRAINT_BUILTINS_H_
#define SimTK_SIMBODY_CONSTRAINT_BUILTINS_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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
Include the header files that define each of the built-in constraint
subclasses of class Constraint. **/

#include "simbody/internal/Constraint_Rod.h"
#include "simbody/internal/Constraint_Ball.h"
#include "simbody/internal/Constraint_Weld.h"
#include "simbody/internal/Constraint_PointInPlane.h"
#include "simbody/internal/Constraint_PointOnPlaneContact.h"
#include "simbody/internal/Constraint_SphereOnPlaneContact.h"
#include "simbody/internal/Constraint_SphereOnSphereContact.h"
#include "simbody/internal/Constraint_LineOnLineContact.h"

#endif // SimTK_SIMBODY_CONSTRAINT_BUILTINS_H_



