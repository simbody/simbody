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
Constraint::SphereOnSphereContact, and its implementation class
Constraint::SphereOnSphereContactImpl. */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Constraint_SphereOnSphereContact.h"

#include "Constraint_SphereOnSphereContactImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {


//==============================================================================
//                         SPHERE ON SPHERE CONTACT
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::SphereOnSphereContact, 
                                        Constraint::SphereOnSphereContactImpl, 
                                        Constraint);

Constraint::SphereOnSphereContact::SphereOnSphereContact
   (MobilizedBody&      mobod_F, 
    const Vec3&         defaultCenter_F, 
    Real                defaultRadius_F, 
    MobilizedBody&      mobod_B, 
    const Vec3&         defaultCenter_B,
    Real                defaultRadius_B,
    bool                enforceRolling)
:   Constraint(new SphereOnSphereContactImpl(enforceRolling))
{
    SimTK_ASSERT_ALWAYS(mobod_F.isInSubsystem() && mobod_B.isInSubsystem(),
        "Constraint::SphereOnSphereContact(): both bodies must already be in a "
        "SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(mobod_F.isInSameSubsystem(mobod_B),
        "Constraint::SphereOnSphereContact(): both bodies to be connected must be "
        "in the same SimbodyMatterSubsystem.");

    mobod_F.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_mobod_F         = updImpl().addConstrainedBody(mobod_F);
    updImpl().m_mobod_B         = updImpl().addConstrainedBody(mobod_B);
    updImpl().m_def_p_FSf       = defaultCenter_F;
    updImpl().m_def_radius_F    = defaultRadius_F;
    updImpl().m_def_p_BSb       = defaultCenter_B;
    updImpl().m_def_radius_B    = defaultRadius_B;
}

Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setDefaultCenterOnF(const Vec3& defaultCenter) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_p_FSf = defaultCenter;
    return *this;
}

Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setDefaultRadiusOnF(Real defaultRadius) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_radius_F = defaultRadius;
    return *this;
}

Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setDefaultCenterOnB(const Vec3& defaultCenter) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_p_BSb = defaultCenter;
    return *this;
}

Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setDefaultRadiusOnB(Real defaultRadius) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_radius_B = defaultRadius;
    return *this;
}

const MobilizedBody& Constraint::SphereOnSphereContact::
getMobilizedBodyF() const {
    const SphereOnSphereContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
}
const MobilizedBody& Constraint::SphereOnSphereContact::
getMobilizedBodyB() const {
    const SphereOnSphereContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_B);
}

bool Constraint::SphereOnSphereContact::isEnforcingRolling() const 
{   return getImpl().m_enforceRolling; }

const Vec3& Constraint::SphereOnSphereContact::
getDefaultCenterOnF() const {return getImpl().m_def_p_FSf;}

Real Constraint::SphereOnSphereContact::
getDefaultRadiusOnF() const {return getImpl().m_def_radius_F;}

const Vec3& Constraint::SphereOnSphereContact::
getDefaultCenterOnB() const {return getImpl().m_def_p_BSb;}

Real Constraint::SphereOnSphereContact::
getDefaultRadiusOnB() const {return getImpl().m_def_radius_B;}

const Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setCenterOnF(State& state, const Vec3& sphereCenter) const {
    getImpl().updParameters(state).m_p_FSf = sphereCenter;
    return *this;
}

const Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setRadiusOnF(State& state, Real sphereRadius) const {
    getImpl().updParameters(state).m_radius_F = sphereRadius;
    return *this;
}

const Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setCenterOnB(State& state, const Vec3& sphereCenter) const {
    getImpl().updParameters(state).m_p_BSb = sphereCenter;
    return *this;
}

const Constraint::SphereOnSphereContact& Constraint::SphereOnSphereContact::
setRadiusOnB(State& state, Real sphereRadius) const {
    getImpl().updParameters(state).m_radius_B = sphereRadius;
    return *this;
}

const Vec3& Constraint::SphereOnSphereContact::
getCenterOnF(const State& state) const
{   return getImpl().getParameters(state).m_p_FSf; }
Real Constraint::SphereOnSphereContact::
getRadiusOnF(const State& state) const
{   return getImpl().getParameters(state).m_radius_F; }
const Vec3& Constraint::SphereOnSphereContact::
getCenterOnB(const State& state) const
{   return getImpl().getParameters(state).m_p_BSb; }
Real Constraint::SphereOnSphereContact::
getRadiusOnB(const State& state) const
{   return getImpl().getParameters(state).m_radius_B; }

Real Constraint::SphereOnSphereContact::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Vec3 Constraint::SphereOnSphereContact::getVelocityErrors(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();
    Vec3 verr_PC; // result is velocity error in P frame 
    if (impl.m_enforceRolling) {
        Real verr[3];
        impl.getVelocityErrors(s, 3, verr);
        verr_PC = Vec3(verr[1],verr[2],verr[0]); // switch to x,y,z order
    } else {
        Real pverr;
        getImpl().getVelocityErrors(s, 1, &pverr);
        verr_PC = Vec3(0,0,pverr); // lone error is in z direction
    }
    return verr_PC;
}

Vec3 Constraint::SphereOnSphereContact::getAccelerationErrors(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();
    Vec3 aerr_PC; // result is acceleration error in P frame 
    if (impl.m_enforceRolling) {
        Real aerr[3];
        impl.getAccelerationErrors(s, 3, aerr);
        aerr_PC = Vec3(aerr[1],aerr[2],aerr[0]); // switch to x,y,z order
    } else {
        Real paerr;
        getImpl().getAccelerationErrors(s, 1, &paerr);
        aerr_PC = Vec3(0,0,paerr); // lone error is in z direction
    }
    return aerr_PC;
}

Vec3 Constraint::SphereOnSphereContact::getMultipliers(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();
    Vec3 lambda_PC; // result is -force on point F in P frame 
    if (impl.m_enforceRolling) {
        Real lambda[3];
        impl.getMultipliers(s, 3, lambda);
        lambda_PC = Vec3(lambda[1],lambda[2],lambda[0]); //switch to x,y,z order
    } else {
        Real lambda;
        getImpl().getMultipliers(s, 1, &lambda);
        lambda_PC = Vec3(0,0,lambda); // lone force is in z direction
    }
    return lambda_PC;
}

Vec3 Constraint::SphereOnSphereContact::
findForceOnSphereBInG(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();
    if (impl.isDisabled(s)) 
        return Vec3(0);

    const Transform X_GC = findContactFrameInG(s);

    const Vec3 f_C = -getMultipliers(s); // watch sign convention
    return X_GC.R()*f_C; // TODO: return f_GC
}

// The contact frame C is defined first in the F frame so that it is dependent
// only on the relative pose of F and B. The origin is given by
//      Co=Sf + rf/(rf+rb) * p_SbSf.
// The z direction Cz is the contact normal given by
//      Cz=p_SbSf/||p_SbSf||
// (arbitrary if the centers are coincident, which is pathological).
// The x-y directions are an arbitrary parameterization of the plane 
// perpendicular to Cz. They are calculated using the Rotation class algorithm
// for constructing a frame given only one axis (in the F frame). That will tend
// to align frame C (vaguely) with the F frame coordinate axes. 
Transform Constraint::SphereOnSphereContact::
findContactFrameInG(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();
    const SphereOnSphereContactImpl::PositionCache& pc = 
        impl.ensurePositionCacheRealized(s);

    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
    const Transform& X_GF = bodyF.getBodyTransform(s);

    return X_GF * pc.X_FC;
}

// The separation is the difference between the spheres's center-to-center
// distance and the sum of their radii. 
Real Constraint::SphereOnSphereContact::
findSeparation(const State& s) const {
    const SphereOnSphereContactImpl& impl = getImpl();

    const SphereOnSphereContactImpl::Parameters& params = impl.getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Real        rf    = params.m_radius_F;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rb    = params.m_radius_B;

    const MobilizedBody& bodyF =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
    const MobilizedBody& bodyB =
        impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_B);

    const Vec3 p_FSb = bodyB.findStationLocationInAnotherBody(s,p_BSb,bodyF);
    const Vec3 p_SfSb_F = p_FSb - p_FSf;
    return p_SfSb_F.norm() - (rf+rb);
}


//==============================================================================
//                     SPHERE ON SPHERE CONTACT IMPL
//==============================================================================

// The default plane and sphere parameters may be overridden by setting
// a discrete variable in the state. We allocate the state resources here.
void Constraint::SphereOnSphereContactImpl::
realizeTopologyVirtual(State& state) const {
    m_parametersIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Position, 
            new Value<Parameters>(Parameters(m_def_p_FSf, m_def_radius_F, 
                                             m_def_p_BSb, m_def_radius_B)));

    m_posCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Position, 
            new Value<PositionCache>());

    m_velCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Velocity, 
            new Value<VelocityCache>());
}

// Return the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Note that although
// these are used to define the position error, only the station on body 2
// is used to generate constraint forces; the point of body 1 that is 
// coincident with the body 2 point receives the equal and opposite force.
const Constraint::SphereOnSphereContactImpl::Parameters& 
Constraint::SphereOnSphereContactImpl::
getParameters(const State& state) const {
    return Value<Parameters>::downcast
       (getMyMatterSubsystemRep().getDiscreteVariable(state,m_parametersIx));
}

// Return a writable reference into the Instance-stage state variable 
// containing the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Calling this
// method invalidates the Instance stage and above in the given state.
Constraint::SphereOnSphereContactImpl::Parameters& 
Constraint::SphereOnSphereContactImpl::
updParameters(State& state) const {
    return Value<Parameters>::updDowncast
       (getMyMatterSubsystemRep().updDiscreteVariable(state,m_parametersIx));
}

const Constraint::SphereOnSphereContactImpl::PositionCache& 
Constraint::SphereOnSphereContactImpl::
getPositionCache(const State& state) const {
    return Value<PositionCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_posCacheIx));
}

Constraint::SphereOnSphereContactImpl::PositionCache& 
Constraint::SphereOnSphereContactImpl::
updPositionCache(const State& state) const {
    return Value<PositionCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_posCacheIx));
}

const Constraint::SphereOnSphereContactImpl::VelocityCache& 
Constraint::SphereOnSphereContactImpl::
getVelocityCache(const State& state) const {
    return Value<VelocityCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_velCacheIx));
}

Constraint::SphereOnSphereContactImpl::VelocityCache& 
Constraint::SphereOnSphereContactImpl::
updVelocityCache(const State& state) const {
    return Value<VelocityCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_velCacheIx));
}

// This costs about 164 flops.
const Constraint::SphereOnSphereContactImpl::PositionCache& 
Constraint::SphereOnSphereContactImpl::
ensurePositionCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_posCacheIx))
        return getPositionCache(s);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rf    = params.m_radius_F;
    const Real        rb    = params.m_radius_B;
    const Real        rtot  = rf+rb; // 1 flop

    PositionCache& pc = updPositionCache(s);

    const Transform&  X_AF = getBodyTransformFromState(s, m_mobod_F);
    const Transform&  X_AB = getBodyTransformFromState(s, m_mobod_B);

    pc.p_FSf_A = X_AF.R() * p_FSf;      // exp. in A, 15 flops
    const Vec3 p_ASf = X_AF.p() + pc.p_FSf_A; // meas. from Ao, 3 flops

    pc.p_BSb_A = X_AB.R() * p_BSb;      // exp. in A, 15 flops
    const Vec3 p_ASb = X_AB.p() + pc.p_BSb_A; // meas. from Ao, 3 flops

    pc.p_SfSb_A = p_ASb - p_ASf;  // vec from Sf to Sb, exp. in A, 3 flops
    pc.r = pc.p_SfSb_A.norm();    // ~20 flops
    pc.oor = 1/pc.r;              // ~10 flops (might be Infinity)
    if (pc.r < TinyReal) {
        pc.isSingular = true;
        pc.Cz_A = X_AF.z();
    } else {
        pc.isSingular = false;
        pc.Cz_A = UnitVec3(pc.p_SfSb_A * pc.oor, true); // 3 flops
    }
    pc.perr = pc.r - rtot;        // 1 flop

    // Now compute the contact frame C, in F.
    const UnitVec3 Cz_F = ~X_AF.R() * pc.Cz_A; // 15 flops

    // Place the contact point along the center-to-center line.
    pc.X_FC.updP() = Vec3(p_FSf + (rf/rtot)*pc.r * Cz_F); // ~16 flops
    pc.X_FC.updR().setRotationFromOneAxis(Cz_F, ZAxis);   // ~60 flops

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_posCacheIx);

    return pc;
}

// This costs about 57 flops if position info has already been calculated,
// otherwise we also pay for ensurePositionCacheRealized().
const Constraint::SphereOnSphereContactImpl::VelocityCache& 
Constraint::SphereOnSphereContactImpl::
ensureVelocityCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_velCacheIx))
        return getVelocityCache(s);

    const PositionCache& pc = ensurePositionCacheRealized(s);
    VelocityCache& vc = updVelocityCache(s);

    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_mobod_F);
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_mobod_B);
    const Vec3&       w_AB = V_AB[0];
    const Vec3&       v_AB = V_AB[1];

    // These are d/dt_A p_FSf and d/dt_A p_BSb
    const Vec3 wX_p_FSf_A = w_AF % pc.p_FSf_A;  // 9 flops
    const Vec3 wX_p_BSb_A = w_AB % pc.p_BSb_A;  // 9 flops
    const Vec3 v_ASf = v_AF + wX_p_FSf_A;       // 3 flops
    const Vec3 v_ASb = v_AB + wX_p_BSb_A;       // 3 flops

    // These are the Coriolis accelerations of Sf and Sb, needed later.
    vc.wXwX_p_FSf_A = w_AF % wX_p_FSf_A;        // 9 flops
    vc.wXwX_p_BSb_A = w_AB % wX_p_BSb_A;        // 9 flops

    vc.pd_SfSb_A = v_ASb - v_ASf;            // 3 flops
    vc.verr = ~vc.pd_SfSb_A * pc.Cz_A;       // 5 flops

    // Calculate d/dt_A Cz
    vc.Czd_A = pc.isSingular 
        ? w_AF % pc.Cz_A // rare
        : pc.oor*(vc.pd_SfSb_A - vc.verr*pc.Cz_A);   // 7 flops

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_velCacheIx);
    
    return vc;
}

void Constraint::SphereOnSphereContactImpl::
calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the spheres' centers and
    // radii, which might not be until Position stage.
    if (   stage == Stage::Position 
        && getMyMatterSubsystemRep().getShowDefaultGeometry()) 
    {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Parameters& params = getParameters(s);
        const Vec3&       p_FSf = params.m_p_FSf;
        const Real        rf    = params.m_radius_F;
        const Vec3&       p_BSb = params.m_p_BSb;
        const Real        rb    = params.m_radius_B;

        const MobilizedBodyIndex mobodFIx = 
            getMobilizedBodyIndexOfConstrainedBody(m_mobod_F);
        const MobilizedBodyIndex mobodBIx = 
            getMobilizedBodyIndexOfConstrainedBody(m_mobod_B);

        // On body F draw a green mesh sphere.
        geom.push_back(DecorativeSphere(rf)
            .setColor(Green)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setBodyId(mobodFIx)
            .setTransform(p_FSf));

        // On the ball body draw an orange mesh sphere.
        geom.push_back(DecorativeSphere(rb)
            .setColor(Orange)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setBodyId(mobodBIx)
            .setTransform(p_BSb));
    }
}


} // namespace SimTK

