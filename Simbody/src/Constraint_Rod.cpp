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

Constraint::Rod::Rod(MobilizedBody& mobod_F, MobilizedBody& mobod_B,
                     Real defaultRodLength)
:   Constraint(new RodImpl())
{
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSubsystem() && mobod_B.isInSubsystem(),
        "Constraint::Rod","Rod",
        "Both mobilized bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSameSubsystem(mobod_B),
        "Constraint::Rod","Rod",
        "The two mobilized bodies to be connected must be in the same "
        "SimbodyMatterSubsystem.");
    SimTK_APIARGCHECK1_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod","Rod",
        "The rod length (distance) was %g but must be greater than zero. "
        "Use Constraint::Ball to make two points coincident.",
        defaultRodLength);

    mobod_F.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_mobod_F         = updImpl().addConstrainedBody(mobod_F);
    updImpl().m_mobod_B         = updImpl().addConstrainedBody(mobod_B);
    updImpl().m_def_p_FSf       = Vec3(0); // Sf = Fo
    updImpl().m_def_p_BSb       = Vec3(0); // Sb = Bo
    updImpl().m_def_length      = defaultRodLength;
}

Constraint::Rod::Rod(MobilizedBody& mobod_F, const Vec3& station_F,
                     MobilizedBody& mobod_B, const Vec3& station_B,
                     Real defaultRodLength)
:   Constraint(new RodImpl())
{
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSubsystem() && mobod_B.isInSubsystem(),
        "Constraint::Rod","Rod",
        "Both mobilized bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSameSubsystem(mobod_B),
        "Constraint::Rod","Rod",
        "The two mobilized bodies to be connected must be in the same "
        "SimbodyMatterSubsystem.");
    SimTK_APIARGCHECK1_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod","Rod",
        "The rod length (distance) was %g but must be greater than zero. "
        "Use Constraint::Ball to make two points coincident.",
        defaultRodLength);

    mobod_F.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_mobod_F         = updImpl().addConstrainedBody(mobod_F);
    updImpl().m_mobod_B         = updImpl().addConstrainedBody(mobod_B);
    updImpl().m_def_p_FSf       = station_F;
    updImpl().m_def_p_BSb       = station_B;
    updImpl().m_def_length      = defaultRodLength;
}

Constraint::Rod& Constraint::Rod::
setDefaultPointOnBody1(const Vec3& p1) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_p_FSf = p1;
    return *this;
}

Constraint::Rod& Constraint::Rod::
setDefaultPointOnBody2(const Vec3& p2) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_p_BSb = p2;
    return *this;
}

Constraint::Rod& Constraint::Rod::
setDefaultRodLength(Real length) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_length = length;
    return *this;
}

const Constraint::Rod& Constraint::Rod::
setPointOnBody1(State& state, const Vec3& point) const {
    getImpl().updParameters(state).m_p_FSf = point;
    return *this;
}
const Constraint::Rod& Constraint::Rod::
setPointOnBody2(State& state, const Vec3& point) const {
    getImpl().updParameters(state).m_p_BSb = point;
    return *this;
}
const Constraint::Rod& Constraint::Rod::
setRodLength(State& state, Real rodLength) const {
    getImpl().updParameters(state).m_length = rodLength;
    return *this;
}

const Vec3& Constraint::Rod::getPointOnBody1(const State& state) const
{   return getImpl().getParameters(state).m_p_FSf; }
const Vec3& Constraint::Rod::getPointOnBody2(const State& state) const
{   return getImpl().getParameters(state).m_p_BSb; }
Real Constraint::Rod::getRodLength(const State& state) const
{   return getImpl().getParameters(state).m_length; }

const MobilizedBody& Constraint::Rod::getMobilizedBody1() const {
    const RodImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
}
const MobilizedBody& Constraint::Rod::getMobilizedBody2() const {
    const RodImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_B);
}
MobilizedBodyIndex Constraint::Rod::getBody1MobilizedBodyIndex() const {
    const RodImpl& impl = getImpl();
    return impl.getMobilizedBodyIndexOfConstrainedBody(impl.m_mobod_F);
}
MobilizedBodyIndex Constraint::Rod::getBody2MobilizedBodyIndex() const {
    const RodImpl& impl = getImpl();
    return impl.getMobilizedBodyIndexOfConstrainedBody(impl.m_mobod_B);
}

const Vec3& Constraint::Rod::getDefaultPointOnBody1() const {
    return getImpl().m_def_p_FSf;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody2() const {
    return getImpl().m_def_p_BSb;
}
Real Constraint::Rod::getDefaultRodLength() const {
    return getImpl().m_def_length;
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

// The multiplier is the force on body 2 along the line from point 1 to
// point 2, which would be positive in compression. But multipliers have the
// opposite sign from applied forces, so the multiplier is positive when the
// rod is in tension. Hence we can just return it as the tension.
Real Constraint::Rod::getRodTension(const State& s) const {
    return getMultiplier(s);
}

UnitVec3 Constraint::Rod::findRodOrientationInG(const State& s) const {
    const RodImpl& impl = getImpl();
    const RodImpl::PositionCache& pc =
        impl.ensurePositionCacheRealized(s);

    const MobilizedBody& mobod_A = impl.getAncestorMobilizedBody();
    if (mobod_A.isGround())
        return pc.Cz_A; // == Cz_G

    const Transform& X_GA = mobod_A.getBodyTransform(s);

    return X_GA.R() * pc.Cz_A;  // 15 flops
}

Real Constraint::Rod::findLengthViolation(const State& s) const {
    const RodImpl& impl = getImpl();

    const RodImpl::Parameters& params = impl.getParameters(s);
    const Vec3&       p_FSf  = params.m_p_FSf;
    const Vec3&       p_BSb  = params.m_p_BSb;
    const Real        length = params.m_length;

    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
    const MobilizedBody& bodyB =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_B);

    const Vec3 p_FSb = bodyB.findStationLocationInAnotherBody(s,p_BSb,bodyF);
    const Vec3 p_SfSb_F = p_FSb - p_FSf;
    return p_SfSb_F.norm() - length;
}

//==============================================================================
//                                ROD IMPL
//==============================================================================

// The default point location and rod length parameters may be overridden by
// setting a discrete variable in the state. Also, we want to cache expensive
// position- and velocity-level calculations since they will be reused
// repeatedly. We allocate the state resources here.
void Constraint::RodImpl::
realizeTopologyVirtual(State& state) const {
    m_parametersIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Position,
            new Value<Parameters>(Parameters(m_def_p_FSf,  m_def_p_BSb,
                                             m_def_length)));

    m_posCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Position,
            new Value<PositionCache>());

    m_velCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Velocity,
            new Value<VelocityCache>());
}

const Constraint::RodImpl::Parameters& Constraint::RodImpl::
getParameters(const State& state) const {
    return Value<Parameters>::downcast
       (getMyMatterSubsystemRep().getDiscreteVariable(state,m_parametersIx));
}

Constraint::RodImpl::Parameters& Constraint::RodImpl::
updParameters(State& state) const {
    return Value<Parameters>::updDowncast
       (getMyMatterSubsystemRep().updDiscreteVariable(state,m_parametersIx));
}

const Constraint::RodImpl::PositionCache& Constraint::RodImpl::
getPositionCache(const State& state) const {
    return Value<PositionCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_posCacheIx));
}

Constraint::RodImpl::PositionCache& Constraint::RodImpl::
updPositionCache(const State& state) const {
    return Value<PositionCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_posCacheIx));
}

const Constraint::RodImpl::VelocityCache& Constraint::RodImpl::
getVelocityCache(const State& state) const {
    return Value<VelocityCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_velCacheIx));
}

Constraint::RodImpl::VelocityCache& Constraint::RodImpl::
updVelocityCache(const State& state) const {
    return Value<VelocityCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_velCacheIx));
}

// This costs about 72 flops.
const Constraint::RodImpl::PositionCache& Constraint::RodImpl::
ensurePositionCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_posCacheIx))
        return getPositionCache(s);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf  = params.m_p_FSf;
    const Vec3&       p_BSb  = params.m_p_BSb;
    const Real        length = params.m_length;

    PositionCache& pc = updPositionCache(s);

    const Transform&  X_AF = getBodyTransformFromState(s, m_mobod_F);
    const Transform&  X_AB = getBodyTransformFromState(s, m_mobod_B);
    const Vec3& p_AF = X_AF.p();
    const Vec3& p_AB = X_AB.p();

    pc.p_FSf_A = X_AF.R() * p_FSf;            // exp. in A, 15 flops
    const Vec3 p_ASf = X_AF.p() + pc.p_FSf_A; // meas. from Ao, 3 flops

    pc.p_BSb_A = X_AB.R() * p_BSb;            // exp. in A, 15 flops
    const Vec3 p_ASb = X_AB.p() + pc.p_BSb_A; // meas. from Ao, 3 flops

    pc.p_SfSb_A = p_ASb - p_ASf;  // vec from Sf to Sb, exp. in A, 3 flops
    pc.r = pc.p_SfSb_A.norm();    // ~20 flops
    pc.oor = 1/pc.r;              // ~10 flops (might be Infinity)

    // Assume non-singular.
    pc.Cz_A = UnitVec3(pc.p_SfSb_A * pc.oor, true); // 3 flops
    pc.isSingular = false;
    if (pc.r < TinyReal) {
        pc.isSingular = true;
        pc.Cz_A = X_AF.z(); // arbitrary
    }

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_posCacheIx);
    return pc;
}

// This costs about 57 flops if position info has already been calculated,
// otherwise we also pay for ensurePositionCacheRealized().
const Constraint::RodImpl::VelocityCache&
Constraint::RodImpl::
ensureVelocityCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_velCacheIx))
        return getVelocityCache(s);

    const PositionCache& pc = ensurePositionCacheRealized(s);
    VelocityCache& vc = updVelocityCache(s);

    const UnitVec3& Cz_A = pc.Cz_A;

    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_mobod_F);
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_mobod_B);
    const Vec3&       w_AB = V_AB[0];
    const Vec3&       v_AB = V_AB[1];

    // These are d/dt_A p_FSf and d/dt_A p_BSb
    const Vec3 wX_p_FSf_A = w_AF % pc.p_FSf_A;      // 9 flops
    const Vec3 wX_p_BSb_A = w_AB % pc.p_BSb_A;      // 9 flops
    const Vec3 v_ASf = v_AF + wX_p_FSf_A;           // 3 flops
    const Vec3 v_ASb = v_AB + wX_p_BSb_A;           // 3 flops
    vc.pd_SfSb_A = v_ASb - v_ASf;                   // 3 flops

    // These are the Coriolis accelerations of Sf and Sb, needed later.
    vc.wXwX_p_FSf_A = w_AF % wX_p_FSf_A;            // 9 flops
    vc.wXwX_p_BSb_A = w_AB % wX_p_BSb_A;            // 9 flops

    // Calculate d/dt_A Cz.
    vc.Czd_A = pc.isSingular
        ? w_AF % Cz_A // rare
        : pc.oor*(vc.pd_SfSb_A - (~vc.pd_SfSb_A*Cz_A)*Cz_A); // 12 flops

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_velCacheIx);

    return vc;
}

void Constraint::RodImpl::calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    if (!getMyMatterSubsystemRep().getShowDefaultGeometry())
        return;

    // We can't generate the endpoint artwork until we know the end point
    // stations, which could be as late as Stage::Position.
    if (stage == Stage::Position) {
        const Parameters& params = getParameters(s);
        const Vec3&       p_FSf  = params.m_p_FSf;
        const Vec3&       p_BSb  = params.m_p_BSb;
        const Real        length = params.m_length;

        const MobilizedBodyIndex body1 =
            getMobilizedBodyIndexOfConstrainedBody(m_mobod_F);
        const MobilizedBodyIndex body2 =
            getMobilizedBodyIndexOfConstrainedBody(m_mobod_B);

        // Draw a blue point at the point on body F.
        geom.push_back(DecorativePoint(p_FSf)
                            .setColor(Green).setScaleFactors(Vec3(2))
                            .setLineThickness(3)
                            .setBodyId(body1));

        // Draw a purple point at the point on body B.
        geom.push_back(DecorativePoint(p_BSb)
                            .setColor(Orange).setScaleFactors(Vec3(2))
                            .setLineThickness(3)
                            .setBodyId(body2));

        const Vec3 p_GP1 = getMobilizedBodyFromConstrainedBody(m_mobod_F)
                              .findStationLocationInGround(s, p_FSf);
        const Vec3 p_GP2 = getMobilizedBodyFromConstrainedBody(m_mobod_B)
                              .findStationLocationInGround(s, p_BSb);

        const Vec3 p_P1P2 = p_GP2 - p_GP1;
        const Real d = p_P1P2.norm();

        if (d >= SignificantReal) {
            const Vec3 endPoint = p_GP1 + length * p_P1P2/d;
            geom.push_back(DecorativeLine(p_GP1, endPoint)
                                            .setColor(Gray)
                                            .setLineThickness(3)
                                            .setBodyId(GroundIndex));
        }
    }

}


} // namespace SimTK

