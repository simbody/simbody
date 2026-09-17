/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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
#include "RigidBodyNodeSpec_CantileverFreeBeam.h"
#include "RigidBodyNodeSpec_Custom.h"
#include "RigidBodyNodeSpec_Translation.h"

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
    return new RBNodeTorsion(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}


RigidBodyNode* MobilizedBody::SliderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeSlider(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}


RigidBodyNode* MobilizedBody::BallImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeBall(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}


RigidBodyNode* MobilizedBody::FreeImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeFree(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}


RigidBodyNode* MobilizedBody::ScrewImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeScrew(
        getDefaultRigidBodyMassProperties(),
        getDefaultPitch(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}


RigidBodyNode* MobilizedBody::UniversalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeUJoint(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::CylinderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeCylinder(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::BendStretchImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeBendStretch(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::PlanarImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodePlanar(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::SphericalCoordsImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeSphericalCoords(
        getDefaultRigidBodyMassProperties(),
        az0, negAz, ze0, negZe, axisT, negT,
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::GimbalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeGimbal(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::BushingImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeBushing(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::EllipsoidImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeEllipsoid(
        getDefaultRigidBodyMassProperties(),
        getDefaultRadii(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::LineOrientationImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeLineOrientation(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::FreeLineImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeFreeLine(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::CantileverFreeBeamImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeCantileverFreeBeam(
        getDefaultRigidBodyMassProperties(),
        getDefaultLength(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::TranslationImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeTranslate(
        getDefaultRigidBodyMassProperties(),
        isReversed(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::CustomImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    switch (getImplementation().getImpl().getNU()) {
    case 1:
        return new RBNodeCustom<1>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    case 2:
        return new RBNodeCustom<2>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    case 3:
        return new RBNodeCustom<3>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    case 4:
        return new RBNodeCustom<4>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    case 5:
        return new RBNodeCustom<5>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    case 6:
        return new RBNodeCustom<6>(getImplementation(),
            getDefaultRigidBodyMassProperties(),
            isReversed(), nextUSlot, nextUSqSlot, nextQSlot);
    default:
        assert(!"Illegal number of degrees of freedom for custom MobilizedBody");
        return 0;
    }
}

