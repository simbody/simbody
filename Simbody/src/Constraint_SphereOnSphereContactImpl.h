#ifndef SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_IMPL_H_

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

/**@file
Private implementation of Constraint::SphereOnSphereContact. **/

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
//                      SPHERE ON SPHERE CONTACT IMPL
//==============================================================================
class Constraint::SphereOnSphereContactImpl : public ConstraintImpl {
public:

// An object of this local type is used as the value of a Position-stage
// discrete state variable that allows sphere parameters to be changed for a 
// particular State.
struct Parameters {
    Parameters(const Vec3& p_FSf, Real rad_F, const Vec3& p_BSb, Real rad_B)
    :   m_p_FSf(p_FSf), m_radius_F(rad_F), m_p_BSb(p_BSb), m_radius_B(rad_B) {}
    Vec3 m_p_FSf;    // sphere center on F
    Real m_radius_F; // radius for F's sphere
    Vec3 m_p_BSb;    // sphere center on B
    Real m_radius_B; // radius for B's sphere
};

struct PositionCache {
    Transform X_AC;         // contact frame in A
    Vec3 p_FSf_A, p_BSb_A;  // stations expressed in A
    Vec3 p_FCo_A, p_BCo_A;  // stations expressed in A
    Vec3 p_SfSb_A;          // vector from Sf to Sb, exp. in A
    Real r, oor;            // r=||p_SfSb_A||; oor=1/r (might be Infinity)
    Real kf;                // fraction along p_SfSb to find Co
    bool isSingular;        // if r is too small to be useful
};

struct VelocityCache {
    Vec3 wXwX_p_FSf_A;      // w_AF x (w_AF x p_FSf_A)  [coriolis accel.]
    Vec3 wXwX_p_BSb_A;      // w_AB x (w_AB x p_BSb_A)  [coriolis accel.]
    Vec3 pd_SfSb_A;         // v_ASb - vASf
    Vec3 vF_BCo_A;          // vA_BCo_A-vA_FCo_A, B/Co station vel in F
    Vec3 wXpd_FCo_A;        // w_AF % pd_FCo_A;
    Vec3 wXpd_BCo_A;        // w_AB % pd_BCo_A;
    Vec3 Czd_A;             // d/dt_A Cz_A -- depends on isSingular
    Vec3 Cxd_A, Cyd_A;      // d/dt_A Cx_A, Cy_A -- calculated from Czd_A
};

explicit SphereOnSphereContactImpl(bool enforceRolling)
:   ConstraintImpl(1, enforceRolling?2:0, 0), 
    m_enforceRolling(enforceRolling), 
    m_def_p_FSf(0), m_def_radius_F(NaN), m_def_p_BSb(0), m_def_radius_B(NaN)
{ }
SphereOnSphereContactImpl* clone() const OVERRIDE_11
{   return new SphereOnSphereContactImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    OVERRIDE_11;

// Allocates the discrete state variable for the parameters, and the cache
// entries.
void realizeTopologyVirtual(State& state) const OVERRIDE_11;

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

/* The normal contact condition is just a distance constraint enforcing that
the two center points Sf and Sb must exactly rf+rb apart. The scalar multiplier
acts along the center-center line, so the force application point doesn't
matter; it would be the same for any contact point along that line so we don't
need to calculate it explicitly for the normal contact constraint; it will
matter for tangential forces.
Let 
    P = Sf, Q = Sb
    p_PQ   = p_AQ - p_AP
    pd_PQ  = d/dt_A p_PQ  = v_AQ - v_AP
    pdd_PQ = d/dt_A pd_PQ = a_AQ - a_AP
    r(q)   = ||p_PQ||                       [center-center distance]

There is a pathological case in which the center points of the two spheres
are overlapping, r(q)=0. In that case we can't use the center-center line as the
contact normal Cz since it is zero length. This condition will be corrected 
eventually by assembly analysis; we just need some consistent return values 
until then. In this case we declare that Cz=Fz:

(1)   Cz = {p_PQ/r,    r != 0                 [contact normal]
           {  Fz,      r == 0

Then our normal contact conditions can always be written like this:

(2)   perr = r(q)-(rf+rb)
(3)   verr = pd_PQ  . Cz
(4)   aerr = pdd_PQ . Cz + pd_PQ . d/dt_A Cz

In case r==0, Cz is fixed in F (it is Fz). So Czd = d/dt_A Cz = w_AF x Cz and
      aerr = pdd_PQ . Cz +  pd_PQ . (w_AF x Cz)
           = pdd_PQ . Cz + (pd_PQ x w_AF) . Cz       [triple product identity]
(5)        = (pdd_PQ - w_AF x pd_PQ) . Cz            [r==0]

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
unit vector aligned with p_SfSb). The multiplier will consequently act along Cz.

(See comments below for friction.)
*/

/*  --------------------------------
    perr = ||p_SfSb|| - (rf+rb)
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
    const OVERRIDE_11
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rf    = params.m_radius_F;
    const Real        rb    = params.m_radius_B;

    // 36 flops for these two together.
    const Vec3 p_ASf =  findStationLocation(allX_AB, m_mobod_F, p_FSf); 
    const Vec3 p_ASb =  findStationLocation(allX_AB, m_mobod_B, p_BSb); 

    const Vec3 p_SfSb_A = p_ASb - p_ASf;    //  3 flops
    const Real r = p_SfSb_A.norm();         // ~20 flops

    // Calculate this scalar using A-frame vectors; any frame would have done.
    perr[0] = r - (rf+rb);  //  2 flops
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
    const OVERRIDE_11
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    const PositionCache& pc = ensurePositionCacheRealized(s);
    const UnitVec3& Cz_A = pc.X_AC.z();

    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3& w_AF = V_AF[0]; const Vec3& v_AF = V_AF[1];

    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];

    const Vec3 v_ASf = v_AF + w_AF % pc.p_FSf_A;        // 12 flops
    const Vec3 v_ASb = v_AB + w_AB % pc.p_BSb_A;        // 12 flops
    const Vec3 pd_SfSb_A = v_ASb - v_ASf;               //  3 flops
    pverr[0] = ~pd_SfSb_A * Cz_A;                       //  5 flops
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
    const OVERRIDE_11
{
    assert(allA_AB.size()==2&&constrainedQDotDot.size()==0&&paerr.size()==1);

    // Ensuring the velocity cache is realized also ensures that the position
    // cache has been realized.
    const VelocityCache& vc = ensureVelocityCacheRealized(s);
    const PositionCache& pc = getPositionCache(s);
    const UnitVec3& Cz_A = pc.X_AC.z();

    // 30 flops for these two together.
    const Vec3 a_ASf = findStationInAAcceleration(s, allA_AB, m_mobod_F, 
                                                  pc.p_FSf_A, vc.wXwX_p_FSf_A);
    const Vec3 a_ASb = findStationInAAcceleration(s, allA_AB, m_mobod_B, 
                                                  pc.p_BSb_A, vc.wXwX_p_BSb_A);

    const Vec3 pdd_SfSb_A = a_ASb - a_ASf; // 2nd deriv in A,     3 flops
    paerr[0] = ~pdd_SfSb_A * Cz_A + ~vc.pd_SfSb_A * vc.Czd_A; // 11 flops
}

// Because the normal force is applied along the line between the centers, it
// does not matter where along that line we apply the force -- the force and
// resulting moment are the same anywhere along the line. Thus we are not 
// required to use the same point in space here -- for convenience we'll just
// apply forces to the sphere centers whose locations we have available.
// We won't have this luxury for tangential forces -- for those the application
// point does matter.
//
// apply f=lambda*Cz to the sphere center on body B,
//       -f          to the sphere center on body F
// 33 flops if position info already available.
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const OVERRIDE_11
{
    assert(multipliers.size()==1&&bodyForcesInA.size()==2&&qForces.size()==0);
    const Real lambda = multipliers[0];

    const PositionCache& pc = ensurePositionCacheRealized(s);
    const UnitVec3& Cz_A = pc.X_AC.z();

    const Vec3 force_A = lambda * Cz_A;     // 3 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BSb_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FSf_A, force_A, bodyForcesInA[m_mobod_F]);
}


// Implementation of virtuals required for nonholonomic constraints.

/* The rolling friction constraint says that the velocity in F of the material 
point of B coincident with the contact point Co, must be zero in the contact 
frame's Cx and Cy directions. 

We have:
(1) p_FCo = p_FSf +   k   p_SfSb, k=rf/(rf+rb)
    p_BCo = p_BSb + (k-1) p_SfSb
(2)       = p_FCo - p_FB
Note that the contact point Co is not a station of either body but moves with
respect to the body frames F and B. However, for the velocity constraint
we want to enforce no relative tangential velocity of the two *material* points
at Co, which we'll call B/Co and F/Co and write
(3) vA_FCo = d/dt_A p_FCo 
           = v_AF + w_AF x p_FCo
(4) vA_BCo = d/dt A p_BCo
             v_AB + w_AB x p_BCo
Then the velocity of B/Co measured in the F frame is
(5) vF_BCo = vA_BCo - vA_FCo
and
(6) verr = [vF_BCo . Cx]
           [vF_BCo . Cy]
is our velocity constraint.

We need to calculate aerr = d/dt verr. The components of verr are two scalars
so we can pick any frame in which to calculate the derivatives; we'll use A. So
from the chain rule we have:
(7) aerr = [d/dt_A vF_BCo . Cx + vF_BCo . d/dt_A Cx]
           [d/dt_A vF_BCo . Cy + vF_BCo . d/dt_A Cy]
         = [vdF_BCo . Cx + vF_BCo . Cxd]
           [vdF_BCo . Cy + vF_BCo . Cyd]

We already have vF_BCo, so we need these three quantities:
    Cxd     = d/dt_A Cx = w_AC x Cx
    Cyd     = d/dt_A Cy = w_AC x Cy
    vFd_BCo = d/dt_A vF_BCo

We have to be careful to account for the fact that Co moves in B, that is
d/dt_B p_BCo != 0. We'll use "d" to mean time derivative in A:
    vdF_BCo = vdA_BCo - vdA_FCo                     from (5)
    vdA_FCo = a_AF + b_AF x p_FCo + w_AF x pd_FCo   from (3)
    vdA_BCo = a_AB + b_AB x p_BCo + w_AB x pd_BCo   from (4)
    pd_FCo  = w_AF x p_FSf + kf * pd_SfSb           from (1)
    pd_BCo  = pdFCo - pd_FB                         from (2)
    pd_FB   = v_AB - v_AF

You can find the velocity-dependent portions of the above calculations, 
including Cxd and Cyd, in ensureVelocityCacheRealized().
*/

/*
    -------------------------------------------
    verr = [vF_BCo_A . Cx_A]
           [vF_BCo_A . Cy_A]
    -------------------------------------------
    37 flops (if position info already calculated)
*/
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const OVERRIDE_11
{
    assert(allV_AB.size()==2 && constrainedU.size()==0 && verr.size()==2);

    const PositionCache& pc = ensurePositionCacheRealized(s);
    const UnitVec3& Cx_A  = pc.X_AC.x();
    const UnitVec3& Cy_A  = pc.X_AC.y();

    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3& w_AF = V_AF[0]; const Vec3& v_AF = V_AF[1];

    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];

    const Vec3 vA_BCo_A = v_AB + w_AB % pc.p_BCo_A; // 12 flops
    const Vec3 vA_FCo_A = v_AF + w_AF % pc.p_FCo_A; // 12 flops
    const Vec3 vF_BCo_A = vA_BCo_A - vA_FCo_A;      //  3 flops

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Cx Cy].
    verr[0] = ~vF_BCo_A * Cx_A;                     // 5 flops
    verr[1] = ~vF_BCo_A * Cy_A;                     // 5 flops
}

/*
    -------------------------------------------
     aerr = [vdF_BCo . Cx + vF_BCo . Cxd]
            [vdF_BCo . Cy + vF_BCo . Cyd]
    -------------------------------------------
    55 flops if position and velocity info already calculated.
*/

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const OVERRIDE_11
{
    assert(allA_AB.size()==2 && constrainedUDot.size()==0 && vaerr.size()==2);

    // Ensuring the velocity cache is realized also ensures that the position
    // cache has been realized.
    const VelocityCache& vc = ensureVelocityCacheRealized(s);
    const PositionCache& pc = getPositionCache(s);

    const UnitVec3& Cx_A  = pc.X_AC.x();
    const UnitVec3& Cy_A  = pc.X_AC.y();
    const Vec3&     Cxd_A = vc.Cxd_A;
    const Vec3&     Cyd_A = vc.Cyd_A;

    const SpatialVec& A_AF = allA_AB[m_mobod_F];
    const Vec3& b_AF = A_AF[0]; const Vec3& a_AF = A_AF[1];

    const SpatialVec& A_AB = allA_AB[m_mobod_B];
    const Vec3& b_AB = A_AB[0]; const Vec3& a_AB = A_AB[1];

    const Vec3 vd_FCo_A = a_AF + b_AF % pc.p_FCo_A + vc.wXpd_FCo_A; // 15 flops
    const Vec3 vd_BCo_A = a_AB + b_AB % pc.p_BCo_A + vc.wXpd_BCo_A; // 15 flops
    const Vec3 vdF_BCo_A = vd_BCo_A - vd_FCo_A;                     //  3 flops

    vaerr[0] = ~vdF_BCo_A*Cx_A + ~vc.vF_BCo_A*Cxd_A;                // 11 flops
    vaerr[1] = ~vdF_BCo_A*Cy_A + ~vc.vF_BCo_A*Cyd_A;                // 11 flops
}

// apply f=lambda0*Px + lambda1*Py to the bottom point C of ball on B
//      -f           to point C (coincident point) of body F
// 39 flops if position info already available
void addInVelocityConstraintForcesVirtual
   (const State&                                    s,      // Stage::Velocity
    const Array_<Real>&                             multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const OVERRIDE_11
{
    assert(multipliers.size()==2 && mobilityForces.size()==0 
           && bodyForcesInA.size()==2);
    const Real lambda0 = multipliers[0], lambda1 = multipliers[1];

    const PositionCache& pc = ensurePositionCacheRealized(s);
    const UnitVec3& Cx_A  = pc.X_AC.x();
    const UnitVec3& Cy_A  = pc.X_AC.y();

    const Vec3 force_A = lambda0*Cx_A + lambda1*Cy_A; // 9 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BCo_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FCo_A, force_A, bodyForcesInA[m_mobod_F]);
}


SimTK_DOWNCAST(SphereOnSphereContactImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::SphereOnSphereContact;

// These can't be changed after construction.
ConstrainedBodyIndex    m_mobod_F;          // F (B1)
ConstrainedBodyIndex    m_mobod_B;          // B (B2)
const bool              m_enforceRolling;   // if so, add 2 constraints

Vec3                    m_def_p_FSf;    // default sphere center on F
Real                    m_def_radius_F; // default radius for F's sphere
Vec3                    m_def_p_BSb;    // default sphere center on B
Real                    m_def_radius_B; // default radius for B's sphere

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

#endif // SimTK_SIMBODY_CONSTRAINT_SPHERE_ON_SPHERE_CONTACT_IMPL_H_



