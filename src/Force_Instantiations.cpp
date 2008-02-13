/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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


#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"


// This suppresses the 'extern template' instantiations in Force.h so that
// we can instantiate them for real here.
#define SimTK_SIMBODY_DEFINING_FORCE
#include "simbody/internal/Force.h"
#include "ForceImpl.h"

namespace SimTK {

template class PIMPLHandle<Force, ForceImpl>;
template class PIMPLImplementation<Force, ForceImpl>;
template class PIMPLDerivedHandle<Force::TwoPointLinearSpring, Force::TwoPointLinearSpringImpl, Force>;
template class PIMPLDerivedHandle<Force::TwoPointLinearDamper, Force::TwoPointLinearDamperImpl, Force>;
template class PIMPLDerivedHandle<Force::TwoPointConstantForce, Force::TwoPointConstantForceImpl, Force>;
template class PIMPLDerivedHandle<Force::MobilityLinearSpring, Force::MobilityLinearSpringImpl, Force>;
template class PIMPLDerivedHandle<Force::MobilityLinearDamper, Force::MobilityLinearDamperImpl, Force>;
template class PIMPLDerivedHandle<Force::MobilityConstantForce, Force::MobilityConstantForceImpl, Force>;
template class PIMPLDerivedHandle<Force::ConstantForce, Force::ConstantForceImpl, Force>;
template class PIMPLDerivedHandle<Force::ConstantTorque, Force::ConstantTorqueImpl, Force>;
template class PIMPLDerivedHandle<Force::GlobalDamper, Force::GlobalDamperImpl, Force>;
template class PIMPLDerivedHandle<Force::UniformGravity, Force::UniformGravityImpl, Force>;
template class PIMPLDerivedHandle<Force::Custom, Force::CustomImpl, Force>;
} // namespace SimTK

