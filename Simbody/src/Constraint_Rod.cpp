/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-14 Stanford University and the Authors.        *
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

/* Implementation of non-inline methods of the handle class Constraint::Rod, and
its implementation class Constraint::RodImpl. */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Constraint_Rod.h"

#include "Constraint_RodImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {

//==============================================================================
//                                   ROD
//==============================================================================

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::Rod, Constraint::RodImpl, 
                                        Constraint);

Constraint::Rod::Rod(MobilizedBody& body1, MobilizedBody& body2, 
                     Real defaultRodLength)
:   Constraint(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultRodLength = defaultRodLength;
}

Constraint::Rod::Rod(MobilizedBody& body1, const Vec3& point1,
                     MobilizedBody& body2, const Vec3& point2, Real defaultRodLength)
  : Constraint(new RodImpl())
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    body1.updMatterSubsystem().adoptConstraint(*this);

    updImpl().B1 = updImpl().addConstrainedBody(body1);
    updImpl().B2 = updImpl().addConstrainedBody(body2);

    updImpl().defaultPoint1 = point1;
    updImpl().defaultPoint2 = point2;
    updImpl().defaultRodLength = defaultRodLength;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody1(const Vec3& p1) {
    updImpl().defaultPoint1 = p1;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody2(const Vec3& p2) {
    updImpl().defaultPoint2 = p2;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultRodLength(Real length) {
    updImpl().defaultRodLength = length;
    return *this;
}


MobilizedBodyIndex Constraint::Rod::getBody1MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B1);
}
MobilizedBodyIndex Constraint::Rod::getBody2MobilizedBodyIndex() const {
    return getImpl().getMobilizedBodyIndexOfConstrainedBody(getImpl().B1);
}
const Vec3& Constraint::Rod::getDefaultPointOnBody1() const {
    return getImpl().defaultPoint1;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody2() const {
    return getImpl().defaultPoint2;
}
Real Constraint::Rod::getDefaultRodLength() const {
    return getImpl().defaultRodLength;
}

Real Constraint::Rod::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Real Constraint::Rod::getVelocityError(const State& s) const {
    Real pverr;
    getImpl().getVelocityErrors(s, 1, &pverr);
    return pverr;
}

Real Constraint::Rod::getAccelerationError(const State& s) const {
    Real pvaerr;
    getImpl().getAccelerationErrors(s, 1, &pvaerr);
    return pvaerr;
}

Real Constraint::Rod::getMultiplier(const State& s) const {
    Real mult;
    getImpl().getMultipliers(s, 1, &mult);
    return mult;
}

//==============================================================================
//                                ROD IMPL
//==============================================================================


void Constraint::RodImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    if (!getMyMatterSubsystemRep().getShowDefaultGeometry())
        return;

    // We can't generate the endpoint artwork until we know the end point stations,
    // which could be as late as Stage::Instance.
    if (stage == Stage::Instance && pointRadius != 0) {
        // TODO: point stations and rod length should be instance-stage data 
        // from State rather than topological data
        const MobilizedBodyIndex body1 = getMobilizedBodyIndexOfConstrainedBody(B1);
        const MobilizedBodyIndex body2 = getMobilizedBodyIndexOfConstrainedBody(B2);

        const Real useRadius = pointRadius > 0 ? pointRadius 
            : std::max(defaultRodLength, Real(.1)) * Real(0.02); // 2% of the length by default

        // Draw a blue mesh sphere at the first point.
        geom.push_back(DecorativeSphere(useRadius)
                            .setColor(Blue)
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setResolution(Real(0.5))
                            .setBodyId(body1)
                            .setTransform(defaultPoint1));

        // On the follower body draw an purple mesh sphere at the point radius.
        geom.push_back(DecorativeSphere(useRadius)
                            .setColor(Purple)
                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                            .setResolution(Real(0.5))
                            .setBodyId(body2)
                            .setTransform(defaultPoint2));
    }

    // We can't generate the line artwork until we know the two end point locations,
    // which isn't until Position stage since the ends are on different bodies.
    if (stage == Stage::Position) {
        // TODO: point stations and rod length should be instance-stage data 
        // from State rather than topological data

        const Vec3 p_GP1 = getMobilizedBodyFromConstrainedBody(B1)
                              .findStationLocationInGround(s, defaultPoint1);
        const Vec3 p_GP2 = getMobilizedBodyFromConstrainedBody(B2)
                              .findStationLocationInGround(s, defaultPoint2);

        const Vec3 p_P1P2 = p_GP2 - p_GP1;
        const Real d = p_P1P2.norm();

        if (d >= SignificantReal) {
            const Vec3 endPoint = p_GP1 + defaultRodLength * p_P1P2/d;
            geom.push_back(DecorativeLine(p_GP1, endPoint)
                                            .setColor(Gray)
                                            .setLineThickness(3)
                                            .setBodyId(GroundIndex));
        }
    }

}


} // namespace SimTK

