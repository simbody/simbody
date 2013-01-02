/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman, Charles Schwieters                            *
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
 * This file contains implementations for any of the classes derived from
 * RigidBodyNodeSpec that aren't defined in the headers.
 */

#include "RigidBodyNodeSpec_Pin.h"
#include "RigidBodyNodeSpec_Slider.h"
#include "RigidBodyNodeSpec_Cylinder.h"
#include "RigidBodyNodeSpec_SphericalCoords.h"
#include "RigidBodyNodeSpec_Ball.h"
#include "RigidBodyNodeSpec_Ellipsoid.h"
#include "RigidBodyNodeSpec_Free.h"
#include "RigidBodyNodeSpec_Screw.h"
#include "RigidBodyNodeSpec_Universal.h"
#include "RigidBodyNodeSpec_PolarCoords.h"
#include "RigidBodyNodeSpec_Planar.h"
#include "RigidBodyNodeSpec_Gimbal.h"
#include "RigidBodyNodeSpec_Bushing.h"
#include "RigidBodyNodeSpec_FreeLine.h"
#include "RigidBodyNodeSpec_LineOrientation.h"
#include "RigidBodyNodeSpec_Custom.h"
// Note: _Translation is handled separately so we can special case
// a lone particle for speed if we find one.

// A macro for instantiating rigid body nodes.
#define INSTANTIATE(CLASS, ...) \
    bool noX_MB = (getDefaultOutboardFrame().p() == 0 && getDefaultOutboardFrame().R() == Mat33(1)); \
    bool noR_PF = (getDefaultInboardFrame().R() == Mat33(1)); \
    if (noX_MB) { \
        if (noR_PF) \
            return new CLASS<true, true> (__VA_ARGS__); \
        else \
            return new CLASS<true, false> (__VA_ARGS__); \
    } \
    else { \
        if (noR_PF) \
            return new CLASS<false, true> (__VA_ARGS__); \
        else \
            return new CLASS<false, false> (__VA_ARGS__); \
    }

    /////////////////////////////////////////////////////////
    // MoblizedBodyImpl::createRigidBodyNode() definitions //
    /////////////////////////////////////////////////////////


// These probably don't belong here.

#include "MobilizedBodyImpl.h"

RigidBodyNode* MobilizedBody::PinImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeTorsion,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}


RigidBodyNode* MobilizedBody::SliderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeSlider,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}


RigidBodyNode* MobilizedBody::BallImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeBall,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}


RigidBodyNode* MobilizedBody::FreeImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeFree,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}



RigidBodyNode* MobilizedBody::ScrewImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeScrew,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        getDefaultPitch(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}


RigidBodyNode* MobilizedBody::UniversalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeUJoint,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::CylinderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeCylinder, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::BendStretchImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeBendStretch, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::PlanarImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodePlanar, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::SphericalCoordsImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeSphericalCoords, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        az0, negAz, ze0, negZe, axisT, negT,
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::GimbalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeGimbal, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::BushingImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeBushing, getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::EllipsoidImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeEllipsoid,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        getDefaultRadii(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::LineOrientationImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeLineOrientation,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

RigidBodyNode* MobilizedBody::FreeLineImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    INSTANTIATE(RBNodeFreeLine,
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        isReversed(),
        nextUSlot,nextUSqSlot,nextQSlot)
}

#define INSTANTIATE_CUSTOM(DOF, ...) \
    if (noX_MB) { \
        if (noR_PF) \
            return new RBNodeCustom<DOF, true, true> (__VA_ARGS__); \
        else \
            return new RBNodeCustom<DOF, true, false> (__VA_ARGS__); \
    } \
    else { \
        if (noR_PF) \
            return new RBNodeCustom<DOF, false, true> (__VA_ARGS__); \
        else \
            return new RBNodeCustom<DOF, false, false> (__VA_ARGS__); \
    }

RigidBodyNode* MobilizedBody::CustomImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    bool noX_MB = (getDefaultOutboardFrame().p() == 0 && getDefaultOutboardFrame().R() == Mat33(1));
    bool noR_PF = (getDefaultInboardFrame().R() == Mat33(1));
    switch (getImplementation().getImpl().getNU()) {
    case 1:
        INSTANTIATE_CUSTOM(1, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    case 2:
        INSTANTIATE_CUSTOM(2, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    case 3:
        INSTANTIATE_CUSTOM(3, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    case 4:
        INSTANTIATE_CUSTOM(4, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    case 5:
        INSTANTIATE_CUSTOM(5, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    case 6:
        INSTANTIATE_CUSTOM(6, getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), isReversed(), nextUSlot, nextUSqSlot, nextQSlot)
    default:
        assert(!"Illegal number of degrees of freedom for custom MobilizedBody");
        return 0;
    }
}

