#ifndef SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_

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

/**@file
Private implementation of Constraint::Rod. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintImpl.h"

class SimbodyMatterSubsystemRep;

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;

//==============================================================================
//                                  ROD IMPL
//==============================================================================
/* This constraint can be implemented using perr = (~p*p - d^2)/2 which results 
in very tidy and cheap equations, but suffers from poor scaling, especially when
combined with other constraints where perr is an actual distance. That's how 
Constraint::Rod was originally implemented but I reworked it to use 
perr=sqrt(~p*p)-d instead. This constraint is nearly identical to the normal 
part of the sphere-on-sphere contact constraint, which says that the distance 
between the two sphere centers must be d=r1+r2. The terminology (and code) used
here is borrowed from sphere-on-sphere; one body is called "F" and the other is
"B" and the connected stations are Sf and Sb. However, I left the user-side
terminology as it was for compatibility. (sherm 140505) 

For comparison, here are the lovely squared equations:
    perr = (p . p - d^2)/2          [NOT using these equatiosn]
    verr =  v . p
    aerr =  a . p + v . v
       f = lambda * p (on Sb), -lambda * p (on Sf)
where
    p = p_ASb-p_ASf, v = d/dt_A p = v_ASb-v_ASf, a = d/dt_A v = a_ASb-a_ASf

Instead, we're going to do it this way:

Let 
    P = Sf, Q = Sb
    p_PQ   = p_AQ - p_AP
    pd_PQ  = d/dt_A p_PQ  = v_AQ - v_AP
    pdd_PQ = d/dt_A pd_PQ = a_AQ - a_AP
    r(q)   = ||p_PQ||                       [center-center distance]

There is a pathological case in which the two connected points are coincident, 
so r(q)=0. In that case we can't use the point-to-point line as the
force direction Cz since it is zero length. This condition will be corrected 
eventually by assembly analysis; we just need some consistent return values 
until then. In this case we arbitrarily declare that Cz=Fz:

(1)   Cz = {p_PQ/r,    r != 0                 [force direction]
           {  Fz,      r == 0

Then our normal contact conditions can always be written like this:

(2)   perr = r(q)-d
(3)   verr = pd_PQ  . Cz
(4)   aerr = pdd_PQ . Cz + pd_PQ . d/dt_A Cz
         f = lambda * Cz (on Sb), -lambda * Cz (on Sf)

In case r==0, Cz is fixed in F (it is Fz). So Czd = d/dt_A Cz = w_AF x Cz and
      aerr = pdd_PQ . Cz +  pd_PQ . (w_AF x Cz)
           = pdd_PQ . Cz + (pd_PQ x w_AF) . Cz       [triple product identity]
(5)        = (pdd_PQ - w_AF x pd_PQ) . Cz            [if r==0]

Otherwise we need to calculate 
       Czd = d/dt_A Cz = d/dt_A p_PQ/r
           = pd_PQ * oor + p_PQ * oord
           = pd_PQ/r - p_PQ * (pd_PQ . p_PQ)/r^3
           = pd_PQ/r - p_PQ/r * (pd_PQ . p_PQ/r)/r
           = [pd_PQ - (pd_PQ . Cz)*Cz]/r
(6)        = (pd_PQ - verr * Cz) / r

where
    oor(q) = 1/r(q) = (p_PQ.p_PQ)^-1/2
    oord(q) = d/dt oor(q) = -pd_PQ . p_PQ/r^3

Substituting (6) into (4) gives
    aerr = pdd_PQ . Cz + pd_PQ . [pd_PQ - (pd_PQ . Cz)*Cz]/r 
         = pdd_PQ . Cz + [pd_PQ^2 - (pd_PQ . Cz)^2] / r
(7)      = pdd_PQ . Cz + [pd_PQ^2 - verr^2] / r

The scalars returned by perr,verr, and aerr are measure numbers along Cz (a 
unit vector normally aligned with p_SfSb). The multiplier will consequently act
along Cz.
*/
class Constraint::RodImpl : public ConstraintImpl {
public:

// An object of this local type is used as the value of a Position-stage
// discrete state variable that allows constraint parameters to be changed for
// a particular State.
struct Parameters {
    Parameters(const Vec3& p_FSf, const Vec3& p_BSb, Real d)
    :   m_p_FSf(p_FSf), m_p_BSb(p_BSb), m_length(d) {}
    Vec3 m_p_FSf;    // sphere center on F
    Vec3 m_p_BSb;    // sphere center on B
    Real m_length;   // the required distance between Sf and Sb
};

struct PositionCache {
    Vec3 p_FSf_A, p_BSb_A;  // stations expressed in A
    Vec3 p_SfSb_A;          // vector from Sf to Sb, exp. in A
    Real r, oor;            // r=||p_SfSb_A||; oor=1/r (might be Infinity)
    UnitVec3 Cz_A;          // force application direction, usually p_SfSb/r
    bool isSingular;        // if r is too small to be useful
};

struct VelocityCache {
    Vec3 wXwX_p_FSf_A;      // w_AF x (w_AF x p_FSf_A)  [coriolis accel.]
    Vec3 wXwX_p_BSb_A;      // w_AB x (w_AB x p_BSb_A)  [coriolis accel.]
    Vec3 pd_SfSb_A;         // v_ASb - vASf
    Vec3 Czd_A;             // d/dt_A Cz_A -- depends on isSingular
};

RodImpl() : ConstraintImpl(1, 0, 0), 
            m_def_p_FSf(0), m_def_p_BSb(0), m_def_length(NaN) {}

RodImpl* clone() const override
{   return new RodImpl(*this); }

// Draw some end points and a rubber band line.
void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    override;

// Allocates the discrete state variable for the parameters, and the cache
// entries.
void realizeTopologyVirtual(State& state) const override;

// Get the current value of the runtime-settable parameters from this state.
const Parameters& getParameters(const State& state) const;

// Get a writable reference to the value of the discrete variable; this 
// invalidated Position stage in the given state.
Parameters& updParameters(State& state) const;

const PositionCache& getPositionCache(const State& state) const;
PositionCache& updPositionCache(const State& state) const;
const PositionCache&  ensurePositionCacheRealized(const State& state) const;

const VelocityCache& getVelocityCache(const State& state) const;
VelocityCache& updVelocityCache(const State& state) const;
const VelocityCache& ensureVelocityCacheRealized(const State& state) const;


// Implementation of virtuals required for holonomic constraints.

/*  --------------------------------
    perr = ||p_SfSb|| - d
    --------------------------------
    ~60 flops

Note that the give State does not necessarily correspond to the given X_AB
and constrainedQ so we can't cache the results.

TODO: there ought to be a flag saying whether the X_AB and q are from the 
state (which is common) so that the results could be saved.
*/
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        d     = params.m_length;

    // 36 flops for these two together.
    const Vec3 p_ASf =  findStationLocation(allX_AB, m_mobod_F, p_FSf); 
    const Vec3 p_ASb =  findStationLocation(allX_AB, m_mobod_B, p_BSb); 

    const Vec3 p_SfSb_A = p_ASb - p_ASf;    //  3 flops
    const Real r = p_SfSb_A.norm();         // ~20 flops

    // Calculate this scalar using A-frame vectors; any frame would have done.
    perr[0] = r - d;  //  1 flop
}

/*  --------------------------------
    verr  =  pd_SfSb_A . Cz
    where pd_SfSb_A = v_ASb - v_ASf
    and   Cz = {p_SfSb_A / r,  r!=0
               {Fz_A,          r==0
    --------------------------------
    32 flops if position info is already available.

Note that we can't cache the velocity results here because we don't know that
V_AB and qdot were taken from the given State.
*/
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    const PositionCache& pc = ensurePositionCacheRealized(s);

    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3& w_AF = V_AF[0]; const Vec3& v_AF = V_AF[1];

    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];

    const Vec3 v_ASf = v_AF + w_AF % pc.p_FSf_A;        // 12 flops
    const Vec3 v_ASb = v_AB + w_AB % pc.p_BSb_A;        // 12 flops
    const Vec3 pd_SfSb_A = v_ASb - v_ASf;               //  3 flops
    pverr[0] = ~pd_SfSb_A * pc.Cz_A;                    //  5 flops
}

/*  ------------------------------------------------------------
    aerr = pdd_SfSb . Cz + pd_SfSb . Czd
    where pdd_SfSb_A = a_ASb - a_ASf
    and   Cz = {p_SfSb_A / r,  r!=0
               {Fz_A,          r==0
    and   Czd = d/dt_A Cz = {(pd_SfSb - Cz * (pd_SfSb . Cz))/r
                            {w_AF x Fz_A,     r==0
    ------------------------------------------------------------ 
    44 flops if position & velocity info already calculated
*/
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const
{
    assert(allA_AB.size()==2&&constrainedQDotDot.size()==0&&paerr.size()==1);

    // Ensuring the velocity cache is realized also ensures that the position
    // cache has been realized.
    const VelocityCache& vc = ensureVelocityCacheRealized(s);
    const PositionCache& pc = getPositionCache(s);

    // 30 flops for these two together.
    const Vec3 a_ASf = findStationInAAcceleration(s, allA_AB, m_mobod_F, 
                                                  pc.p_FSf_A, vc.wXwX_p_FSf_A);
    const Vec3 a_ASb = findStationInAAcceleration(s, allA_AB, m_mobod_B, 
                                                  pc.p_BSb_A, vc.wXwX_p_BSb_A);

    const Vec3 pdd_SfSb_A = a_ASb - a_ASf; // 2nd deriv in A,        3 flops
    paerr[0] = ~pdd_SfSb_A * pc.Cz_A + ~vc.pd_SfSb_A * vc.Czd_A; // 11 flops
}

// Because the normal force is applied along the line between the points, it
// does not matter where along that line we apply the force -- the force and
// resulting moment are the same anywhere along the line. Thus we are not 
// required to use the same point in space here -- for convenience we'll just
// apply forces to the defining end points, whose locations we have available.
//
// apply f=lambda*Cz to the end point on body B (B2),
//       -f          to the end point on body F (B1)
// 33 flops if position info already available.
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const
{
    assert(multipliers.size()==1&&bodyForcesInA.size()==2&&qForces.size()==0);
    const Real lambda = multipliers[0];

    const PositionCache& pc = ensurePositionCacheRealized(s);

    const Vec3 force_A = lambda * pc.Cz_A;     // 3 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BSb_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FSf_A, force_A, bodyForcesInA[m_mobod_F]);
}

SimTK_DOWNCAST(RodImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::Rod;

// These can't be changed after construction.
ConstrainedBodyIndex    m_mobod_F;          // F (B1)
ConstrainedBodyIndex    m_mobod_B;          // B (B2)

Vec3                    m_def_p_FSf;    // default for point on F
Vec3                    m_def_p_BSb;    // default for point on B
Real                    m_def_length;   // default rod length

// This Position-stage variable holds the constraint parameters to be used
// for a particular State.
mutable DiscreteVariableIndex   m_parametersIx;

// These cache entries hold precalculated position- and velocity-dependent
// computations. These are lazy evaluated by the first method to be called
// that requires the State to be at Position or Velocity stage.
mutable CacheEntryIndex         m_posCacheIx;
mutable CacheEntryIndex         m_velCacheIx;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_



