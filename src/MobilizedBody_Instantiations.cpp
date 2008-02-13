/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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


#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"

// This suppresses the 'extern template' instantiations in MobilizedBody.h so that
// we can instantiate them for real here.
#define SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
#include "simbody/internal/MobilizedBody.h"
#include "MobilizedBodyImpl.h"

namespace SimTK {

template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl>;
template class PIMPLImplementation<MobilizedBody, MobilizedBodyImpl>;
template class PIMPLDerivedHandle<MobilizedBody::Ball, MobilizedBody::BallImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::BendStretch, MobilizedBody::BendStretchImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Custom, MobilizedBody::CustomImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Cylinder, MobilizedBody::CylinderImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Ellipsoid, MobilizedBody::EllipsoidImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Free, MobilizedBody::FreeImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::FreeLine, MobilizedBody::FreeLineImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Gimbal, MobilizedBody::GimbalImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Ground, MobilizedBody::GroundImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::LineOrientation, MobilizedBody::LineOrientationImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Pin, MobilizedBody::PinImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Planar, MobilizedBody::PlanarImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Screw, MobilizedBody::ScrewImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Slider, MobilizedBody::SliderImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Translation, MobilizedBody::TranslationImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Universal, MobilizedBody::UniversalImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Weld, MobilizedBody::WeldImpl, MobilizedBody>;

} // namespace SimTK

