/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Force.h"

#include "ForceImpl.h"

namespace SimTK {

class Force::LinearBushingImpl : public ForceImpl {
    struct PositionCache {
        Transform X_GF, X_GM, X_FM;
        Vec3      p_AF_G, p_BM_G, p_FM_G;
        Vec6      q;
    };
    struct VelocityCache {
        SpatialVec V_GF, V_GM, V_FM;
        Vec6       qdot;
    };
    struct ForceCache {
        SpatialVec F_GF, F_GM; // at Bushing frames
        SpatialVec F_GA, F_GB; // at Body frames
        Vec6       f;          // scalar generalized forces
    };
public:
    LinearBushingImpl(const MobilizedBody& bodyA, const Transform& frameOnA, 
                      const MobilizedBody& bodyB, const Transform& frameOnB, 
                      const Vec6& stiffness, const Vec6& damping);
    LinearBushingImpl* clone() const {
        return new LinearBushingImpl(*this);
    }
    bool dependsOnlyOnPositions() const {
        return false;
    }
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;

    // Allocate the position and velocity cache entries. These are all
    // lazy-evaluation entries - be sure to check whether they have already
    // been calculated; calculate them if not; and then mark them done.
    // They will be invalidated when the indicated stage has changed and
    // can be recalculated any time after that stage is realized.
    void realizeTopology(State& s) const {
        LinearBushingImpl* mThis = const_cast<LinearBushingImpl*>(this);
        mThis->positionCacheIx = getForceSubsystem().allocateCacheEntry(s,
            Stage::Position, Stage::Infinity, new Value<PositionCache>());
        mThis->potEnergyCacheIx = getForceSubsystem().allocateCacheEntry(s,
            Stage::Position, Stage::Infinity, new Value<Real>(NaN));
        mThis->velocityCacheIx = getForceSubsystem().allocateCacheEntry(s,
            Stage::Velocity, Stage::Infinity, new Value<VelocityCache>());
        mThis->forceCacheIx = getForceSubsystem().allocateCacheEntry(s,
            Stage::Velocity, Stage::Infinity, new Value<ForceCache>());
    }
private:
    const PositionCache& getPositionCache(const State& s) const
    {   return Value<PositionCache>::downcast
            (getForceSubsystem().getCacheEntry(s,positionCacheIx)); }
    const Real& getPotentialEnergyCache(const State& s) const
    {   return Value<Real>::downcast
            (getForceSubsystem().getCacheEntry(s,potEnergyCacheIx)); }
    const VelocityCache& getVelocityCache(const State& s) const
    {   return Value<VelocityCache>::downcast
            (getForceSubsystem().getCacheEntry(s,velocityCacheIx)); }
    const ForceCache& getForceCache(const State& s) const
    {   return Value<ForceCache>::downcast
            (getForceSubsystem().getCacheEntry(s,forceCacheIx)); }

    PositionCache& updPositionCache(const State& s) const
    {   return Value<PositionCache>::updDowncast
            (getForceSubsystem().updCacheEntry(s,positionCacheIx)); }
    Real& updPotentialEnergyCache(const State& s) const
    {   return Value<Real>::updDowncast
            (getForceSubsystem().updCacheEntry(s,potEnergyCacheIx)); }
    VelocityCache& updVelocityCache(const State& s) const
    {   return Value<VelocityCache>::updDowncast
            (getForceSubsystem().updCacheEntry(s,velocityCacheIx)); }
    ForceCache& updForceCache(const State& s) const
    {   return Value<ForceCache>::updDowncast
            (getForceSubsystem().updCacheEntry(s,forceCacheIx)); }

    bool isPositionCacheValid(const State& s) const
    {   return getForceSubsystem().isCacheValueCurrent(s,positionCacheIx); }
    bool isPotentialEnergyValid(const State& s) const
    {   return getForceSubsystem().isCacheValueCurrent(s,potEnergyCacheIx); }
    bool isVelocityCacheValid(const State& s) const
    {   return getForceSubsystem().isCacheValueCurrent(s,velocityCacheIx); }
    bool isForceCacheValid(const State& s) const
    {   return getForceSubsystem().isCacheValueCurrent(s,forceCacheIx); }

    void markPositionCacheValid(const State& s) const
    {   getForceSubsystem().markCacheValueRealized(s,positionCacheIx); }
    void markPotentialEnergyValid(const State& s) const
    {   getForceSubsystem().markCacheValueRealized(s,potEnergyCacheIx); }
    void markVelocityCacheValid(const State& s) const
    {   getForceSubsystem().markCacheValueRealized(s,velocityCacheIx); }
    void markForceCacheValid(const State& s) const
    {   getForceSubsystem().markCacheValueRealized(s,forceCacheIx); }

    void ensurePositionCacheValid(const State&) const;
    void ensurePotentialEnergyValid(const State&) const;
    void ensureVelocityCacheValid(const State&) const;
    void ensureForceCacheValid(const State&) const;

    // TOPOLOGY STATE
    const MobilizedBody& bodyA;
    const MobilizedBody& bodyB;
    Transform X_AF, X_BM;
    Vec6 k, c;

    // TOPOLOGY CACHE
    CacheEntryIndex positionCacheIx;
    CacheEntryIndex potEnergyCacheIx;
    CacheEntryIndex velocityCacheIx;
    CacheEntryIndex forceCacheIx;

friend class Force::LinearBushing;
friend std::ostream& operator<<(std::ostream&,const PositionCache&);
friend std::ostream& operator<<(std::ostream&,const VelocityCache&);
friend std::ostream& operator<<(std::ostream&,const ForceCache&);
};

// These are required by Value<T>.
inline std::ostream& operator<<
   (std::ostream& o, const Force::LinearBushingImpl::PositionCache& pc)
{   assert(!"implemented"); return o; }
inline std::ostream& operator<<
   (std::ostream& o, const Force::LinearBushingImpl::VelocityCache& vc)
{   assert(!"implemented"); return o; }
inline std::ostream& operator<<
   (std::ostream& o, const Force::LinearBushingImpl::ForceCache& fc)
{   assert(!"implemented"); return o; }


//------------------------------ LinearBushing ---------------------------------
//------------------------------------------------------------------------------

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Force::LinearBushing, 
                                        Force::LinearBushingImpl, Force);

Force::LinearBushing::LinearBushing
   (GeneralForceSubsystem& forces, 
    const MobilizedBody& bodyA, const Transform& frameOnA,
    const MobilizedBody& bodyB, const Transform& frameOnB,
    const Vec6& stiffness,  const Vec6&  damping) 
:   Force(new LinearBushingImpl(bodyA,frameOnA,bodyB,frameOnB,
                                stiffness,damping))
{
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}


Force::LinearBushing::LinearBushing
   (GeneralForceSubsystem& forces, 
    const MobilizedBody& bodyA, // assume body frames
    const MobilizedBody& bodyB,
    const Vec6& stiffness,  const Vec6&  damping) 
:   Force(new LinearBushingImpl(bodyA,Transform(),bodyB,Transform(),
                                stiffness,damping))
{
    updImpl().setForceSubsystem(forces, forces.adoptForce(*this));
}

const Vec6& Force::LinearBushing::
getQ(const State& s) const 
{   getImpl().ensurePositionCacheValid(s); 
    return getImpl().getPositionCache(s).q; }
const Vec6& Force::LinearBushing::
getQDot(const State& s) const 
{   getImpl().ensureVelocityCacheValid(s); 
    return getImpl().getVelocityCache(s).qdot; }

const Transform& Force::LinearBushing::
getX_FM(const State& s) const 
{   getImpl().ensurePositionCacheValid(s); 
    return getImpl().getPositionCache(s).X_FM; }
const SpatialVec& Force::LinearBushing::
getV_FM(const State& s) const 
{   getImpl().ensureVelocityCacheValid(s); 
    return getImpl().getVelocityCache(s).V_FM; }

const Vec6& Force::LinearBushing::
getF(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).f; }
const SpatialVec& Force::LinearBushing::
getF_GF(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).F_GF; }
const SpatialVec& Force::LinearBushing::
getF_GM(const State& s) const
{   getImpl().ensureForceCacheValid(s); 
    return getImpl().getForceCache(s).F_GM; }

const Real& Force::LinearBushing::
getPotentialEnergy(const State& s) const
{   getImpl().ensurePotentialEnergyValid(s); 
    return getImpl().getPotentialEnergyCache(s); }

Force::LinearBushingImpl::LinearBushingImpl
   (const MobilizedBody& bodyA, const Transform& frameOnA,
    const MobilizedBody& bodyB, const Transform& frameOnB,
    const Vec6& stiffness, const Vec6& damping)
:   bodyA(bodyA), X_AF(frameOnA),
    bodyB(bodyB), X_BM(frameOnB), 
    k(stiffness), c(damping) 
{
}

void Force::LinearBushingImpl::
ensurePositionCacheValid(const State& state) const {
    if (isPositionCacheValid(state)) return;

    PositionCache& pc = updPositionCache(state);

    const Transform& X_GA = bodyA.getBodyTransform(state);
    const Transform& X_GB = bodyB.getBodyTransform(state);
    pc.X_GF =     X_GA*   X_AF;   // 63 flops
    pc.X_GM =     X_GB*   X_BM;   // 63 flops
    pc.X_FM = ~pc.X_GF*pc.X_GM;   // 63 flops

    // Re-express local vectors in the Ground frame.
    pc.p_AF_G =    X_GA.R() *    X_AF.p();    // 15 flops
    pc.p_BM_G =    X_GB.R() *    X_BM.p();    // 15 flops
    pc.p_FM_G = pc.X_GF.R() * pc.X_FM.p();    // 15 flops

    // Calculate the 1-2-3 body B-fixed Euler angles; these are the
    // rotational coordinates.
    pc.q.updSubVec<3>(0) = pc.X_FM.R().convertRotationToBodyFixedXYZ();

    // The translation vector in X_FM contains our translational coordinates.
    pc.q.updSubVec<3>(3) = pc.X_FM.p();

    markPositionCacheValid(state);
}

void Force::LinearBushingImpl::
ensureVelocityCacheValid(const State& state) const {
    if (isVelocityCacheValid(state)) return;

    // We'll be needing this.
    ensurePositionCacheValid(state);
    const PositionCache& pc = getPositionCache(state);
    const Rotation& R_GF = pc.X_GF.R();
    const Rotation& R_FM = pc.X_FM.R();
    const Vec3&     q    = pc.q.getSubVec<3>(0);
    const Vec3&     p    = pc.q.getSubVec<3>(3);

    VelocityCache& vc = updVelocityCache(state);

    // Now do velocities.
    const SpatialVec& V_GA = bodyA.getBodyVelocity(state);
    const SpatialVec& V_GB = bodyB.getBodyVelocity(state);

    vc.V_GF = SpatialVec(V_GA[0], V_GA[1] + V_GA[0] % pc.p_AF_G);
    vc.V_GM = SpatialVec(V_GB[0], V_GB[1] + V_GB[0] % pc.p_BM_G);

    // This is the velocity of M in F, but with the time derivative
    // taken in the Ground frame.
    const SpatialVec V_FM_G = vc.V_GM - vc.V_GF;

    // To get derivative in F, we must remove the part due to the
    // angular velocity w_GF of F in G.
    vc.V_FM = ~R_GF * SpatialVec(V_FM_G[0], 
                                 V_FM_G[1] - vc.V_GF[0] % pc.p_FM_G);

    // Need angular velocity in M frame for conversion to qdot.
    const Vec3  w_FM_M = ~R_FM*vc.V_FM[0];
    const Mat33 N_FM   = Rotation::calcNForBodyXYZInBodyFrame(q);
    vc.qdot.updSubVec<3>(0) = N_FM * w_FM_M;
    vc.qdot.updSubVec<3>(3) = vc.V_FM[1];

    markVelocityCacheValid(state);
}

// This will also calculate potential energy since we can do it on the
// cheap simultaneously with the force.
void Force::LinearBushingImpl::
ensureForceCacheValid(const State& state) const {
    if (isForceCacheValid(state)) return;

    ForceCache& fc = updForceCache(state);

    ensurePositionCacheValid(state);
    const PositionCache& pc = getPositionCache(state);

    const Transform& X_GA = bodyA.getBodyTransform(state);
    const Rotation&  R_GA = X_GA.R();

    const Transform& X_GB = bodyB.getBodyTransform(state);
    const Rotation&  R_GB = X_GB.R();

    const Vec3&      p_AF = X_AF.p();
    const Vec3&      p_BM = X_BM.p();

    const Rotation&  R_GF = pc.X_GF.R();
    const Vec3&      p_GF = pc.X_GF.p();

    const Rotation&  R_GM = pc.X_GM.R();
    const Vec3&      p_GM = pc.X_GM.p();

    const Rotation&  R_FM = pc.X_FM.R();
    const Vec3&      p_FM = pc.X_FM.p();

    // Calculate stiffness generalized forces and potential
    // energy (cheap to do here).
    const Vec6& q = pc.q;
    Vec6 fk; Real pe2=0;
    for (int i=0; i<6; ++i) pe2 += (fk[i]=k[i]*q[i])*q[i]; 
    updPotentialEnergyCache(state) = pe2/2;
    markPotentialEnergyValid(state);

    ensureVelocityCacheValid(state);
    const VelocityCache& vc = getVelocityCache(state);

    const Vec6& qd = vc.qdot;
    Vec6 fv;
    for (int i=0; i<6; ++i) fv[i]=c[i]*qd[i];

    fc.f             = -(fk+fv); // generalized forces on body B
    const Vec3& fB_q = fc.f.getSubVec<3>(0); // in q basis
    const Vec3& fM_F = fc.f.getSubVec<3>(3); // acts at M, but exp. in F frame

    // Calculate the matrix relating q-space generalized forces to a real-space
    // moment vector. We know qforce = ~H * moment (where H is the
    // the hinge matrix for a mobilizer using qdots as generalized speeds).
    // In that case H would be N^-1, qforce = ~(N^-1)*moment so
    // moment = ~N*qforce. Caution: our N wants the moment in the outboard
    // body frame, in this case M.
    const Mat33 N_FM = Rotation::calcNForBodyXYZInBodyFrame(q.getSubVec<3>(0));
    const Vec3  mB_M = ~N_FM * fB_q; // moment acting on body B, exp. in M
    const Vec3  mB_G = R_GM*mB_M;    // moment on B, now exp. in G

    // Transform force from F frame to ground. This is the force to 
    // apply to body B at point OM; -f goes on body A at the same
    // spatial location. Here we actually apply it at OF, and since
    // we know the force acts along the line OF-OM the change in 
    // location does not generate a moment.
    const Vec3 fM_G = R_GF*fM_F;

    fc.F_GM = SpatialVec(mB_G, fM_G);
    fc.F_GF = -fc.F_GM; // see above for why force is OK here w/o shift

    // Shift forces to body origins.
    fc.F_GB = SpatialVec(fc.F_GM[0] + pc.p_BM_G % fc.F_GM[1], fc.F_GM[1]);
    fc.F_GA = SpatialVec(fc.F_GF[0] + pc.p_AF_G % fc.F_GF[1], fc.F_GF[1]);

    markForceCacheValid(state);
}

// This calculate is only performed if the PE is requested without
// already having calculated the force.
void Force::LinearBushingImpl::
ensurePotentialEnergyValid(const State& state) const {
    if (isPotentialEnergyValid(state)) return;

    ensurePositionCacheValid(state);
    const PositionCache& pc = getPositionCache(state);
    const Vec6& q = pc.q;

    Real pe2=0;
    for (int i=0; i<6; ++i) pe2 += k[i]*square(q[i]);

    updPotentialEnergyCache(state) = pe2 / 2;
    markPotentialEnergyValid(state);
}

void Force::LinearBushingImpl::
calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
          Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
{
    const MobilizedBodyIndex bodyAx = bodyA.getMobilizedBodyIndex();
    const MobilizedBodyIndex bodyBx = bodyB.getMobilizedBodyIndex();

    ensureForceCacheValid(state);
    const ForceCache& fc = getForceCache(state);
    bodyForces[bodyBx] +=  fc.F_GB;
    bodyForces[bodyAx] +=  fc.F_GA; // apply forces
}

// If the force was calculated, then the potential energy will already
// be valid. Otherwise we'll have to calculate it.
Real Force::LinearBushingImpl::
calcPotentialEnergy(const State& state) const {
    ensurePotentialEnergyValid(state);
    return getPotentialEnergyCache(state);
}


} // namespace SimTK

