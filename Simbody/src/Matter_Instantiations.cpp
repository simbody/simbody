/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

// This suppresses the 'extern template' instantiations in Motion.h so that
// we can instantiate them for real here.
#define SimTK_SIMBODY_DEFINING_MOTION
#include "simbody/internal/Motion.h"
#include "MotionImpl.h"

// This suppresses the 'extern template' instantiations in MobilizedBody.h so that
// we can instantiate them for real here.
#define SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
#include "simbody/internal/MobilizedBody.h"
#include "MobilizedBodyImpl.h"

// This suppresses the 'extern template' instantiations in Constraint.h so that
// we can instantiate them for real here.
#define SimTK_SIMBODY_DEFINING_CONSTRAINT
#include "simbody/internal/Constraint.h"
#include "ConstraintImpl.h"


namespace SimTK {

    // Motion
template class PIMPLHandle<Motion, MotionImpl, true>;
template class PIMPLImplementation<Motion, MotionImpl>;

    // MobilizedBody
template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl, true>;
template class PIMPLImplementation<MobilizedBody, MobilizedBodyImpl>;
template class PIMPLHandle<MobilizedBody::Custom::Implementation, MobilizedBody::Custom::ImplementationImpl>;
template class PIMPLImplementation<MobilizedBody::Custom::Implementation, MobilizedBody::Custom::ImplementationImpl>;

    // Constraint
template class PIMPLHandle<Constraint, ConstraintImpl, true>;
template class PIMPLImplementation<Constraint, ConstraintImpl>;
template class PIMPLHandle<Constraint::Custom::Implementation, Constraint::Custom::ImplementationImpl>;
template class PIMPLImplementation<Constraint::Custom::Implementation, Constraint::Custom::ImplementationImpl>;

} // namespace SimTK

