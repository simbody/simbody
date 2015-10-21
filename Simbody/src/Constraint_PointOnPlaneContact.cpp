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

/* Implementation of non-inline methods of the handle class
Constraint::PointOnPlaneContact, and its implementation class
Constraint::PointOnPlaneContactImpl. */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Constraint_PointOnPlaneContact.h"

#include "Constraint_PointOnPlaneContactImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {


//==============================================================================
//                         POINT ON PLANE CONTACT
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::PointOnPlaneContact, 
                                        Constraint::PointOnPlaneContactImpl, 
                                        Constraint);

Constraint::PointOnPlaneContact::PointOnPlaneContact
   (MobilizedBody& planeBody,    const Transform& defPlaneFrame,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
  : Constraint(new PointOnPlaneContactImpl())
{
    SimTK_ASSERT_ALWAYS(planeBody.isInSubsystem()&&followerBody.isInSubsystem(),
        "Constraint::PointOnPlaneContact(): both bodies must already be "
        "in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(planeBody.isInSameSubsystem(followerBody),
        "Constraint::PointOnPlaneContact(): both bodies to be connected "
        "must be in the same SimbodyMatterSubsystem.");

    //rep = new PointInPlaneRep(); rep->setMyHandle(*this);
    planeBody.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_surfaceBody_S   = updImpl().addConstrainedBody(planeBody);
    updImpl().m_followerBody_B  = updImpl().addConstrainedBody(followerBody);
    updImpl().m_X_SP            = defPlaneFrame;
    updImpl().m_p_BF            = defFollowerPoint;
}

Constraint::PointOnPlaneContact& Constraint::PointOnPlaneContact::
setDefaultPlaneFrame(const Transform& X_SP) {
    getImpl().invalidateTopologyCache();
    updImpl().m_X_SP = X_SP;
    return *this;
}

Constraint::PointOnPlaneContact& Constraint::PointOnPlaneContact::
setDefaultFollowerPoint(const Vec3& p) {
    getImpl().invalidateTopologyCache();
    updImpl().m_p_BF = p;
    return *this;
}

MobilizedBodyIndex Constraint::PointOnPlaneContact::
getPlaneMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody
                                                (getImpl().m_surfaceBody_S);
}
MobilizedBodyIndex Constraint::PointOnPlaneContact::
getFollowerMobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody
                                                (getImpl().m_followerBody_B);
}
const Transform& Constraint::PointOnPlaneContact::
getDefaultPlaneFrame() const {
    return getImpl().m_X_SP;
}

const Vec3& Constraint::PointOnPlaneContact::
getDefaultFollowerPoint() const {
    return getImpl().m_p_BF;
}

Constraint::PointOnPlaneContact& Constraint::PointOnPlaneContact::
setPlaneDisplayHalfWidth(Real h) {
    updImpl().setPlaneDisplayHalfWidth(h);
    return *this;
}
Constraint::PointOnPlaneContact& Constraint::PointOnPlaneContact::
setPointDisplayRadius(Real r) {
    updImpl().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointOnPlaneContact::getPlaneDisplayHalfWidth() const {
    return getImpl().getPlaneDisplayHalfWidth();
}

Real Constraint::PointOnPlaneContact::getPointDisplayRadius() const {
    return getImpl().getPointDisplayRadius();
}

Real Constraint::PointOnPlaneContact::
getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Vec3 Constraint::PointOnPlaneContact::
getVelocityErrors(const State& s) const {
    Real verr[3];
    getImpl().getVelocityErrors(s, 3, verr);
    return Vec3(verr[1],verr[2],verr[0]); // switch to x,y,z order
}

Vec3 Constraint::PointOnPlaneContact::
getAccelerationErrors(const State& s) const {
    Real aerr[3];
    getImpl().getAccelerationErrors(s, 3, aerr);
    return Vec3(aerr[1],aerr[2],aerr[0]); // switch to x,y,z order
}

Vec3 Constraint::PointOnPlaneContact::
getMultipliers(const State& s) const {
    Real lambda[3];
    getImpl().getMultipliers(s, 3, lambda);
    return Vec3(lambda[1],lambda[2],lambda[0]); //switch to x,y,z order;
}

Vec3 Constraint::PointOnPlaneContact::
getForceOnFollowerPointInP(const State& s) const {
    return -getMultipliers(s);
}

//==============================================================================
//                      POINT ON PLANE CONTACT IMPL
//==============================================================================

void Constraint::PointOnPlaneContact::PointOnPlaneContactImpl::
calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the frame and follower
    // point location, which might not be until Instance stage.
    if (   stage == Stage::Instance 
        && getMyMatterSubsystemRep().getShowDefaultGeometry()) 
    {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();

        const MobilizedBodyIndex planeMBId = 
            getMobilizedBodyIndexOfConstrainedBody(m_surfaceBody_S);
        const MobilizedBodyIndex followerMBId = 
            getMobilizedBodyIndexOfConstrainedBody(m_followerBody_B);

        if (m_planeHalfWidth > 0 && m_pointRadius > 0) {
            // On the plane body, draw a gray transparent rectangle, outlined 
            // in black lines.
            geom.push_back(DecorativeBrick
               (Vec3(m_planeHalfWidth,m_planeHalfWidth,m_pointRadius/2))
                .setColor(Gray)
                .setRepresentation(DecorativeGeometry::DrawSurface)
                .setOpacity(Real(0.3))
                .setBodyId(planeMBId)
                .setTransform(m_X_SP));
            geom.push_back(DecorativeBrick
               (Vec3(m_planeHalfWidth,m_planeHalfWidth,m_pointRadius/2))
                .setColor(Black)
                .setRepresentation(DecorativeGeometry::DrawWireframe)
                .setBodyId(planeMBId)
                .setTransform(m_X_SP));

            // On follower body draw an orange mesh sphere using point radius.
            geom.push_back(DecorativeSphere(m_pointRadius)
                .setColor(Orange)
                .setRepresentation(DecorativeGeometry::DrawWireframe)
                .setResolution(Real(0.5))
                .setBodyId(followerMBId)
                .setTransform(m_p_BF));
        }
    }
}



} // namespace SimTK

