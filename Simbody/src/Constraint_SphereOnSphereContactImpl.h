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
    Transform X_FC;         // contact frame in F
    UnitVec3  Cz_A;         // X_FC.z() re-expressed in A
    Vec3 p_FSf_A, p_BSb_A;  // stations expressed in A
    Vec3 p_SfSb_A;          // vector from Sf to Sb, exp. in A
    Real r;                 // ||p_SfSb_A||
    Real oor;               // 1/r (might be Infinity)
    Real perr;              // perr = r - (rf+rb)
    bool isSingular;        // if r is too small to be useful
};

struct VelocityCache {
    Vec3 wXwX_p_FSf_A;      // w_AF x (w_AF x p_FSf_A)  [coriolis accel.]
    Vec3 wXwX_p_BSb_A;      // w_AB x (w_AB x p_BSb_A)  [coriolis accel.]
    Vec3 pd_SfSb_A;         // v_ASb - vASf
    Vec3 Czd_A;             // d/dt_A Cz_A -- depends on isSingular
    Real verr;              // verr = p_SfSb_A_dot . Cz
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
    33 flops if position info is already available.

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

    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3& w_AF = V_AF[0]; const Vec3& v_AF = V_AF[1];

    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];

    const Vec3 v_ASf = v_AF + w_AF % pc.p_FSf_A;        // 12 flops
    const Vec3 v_ASb = v_AB + w_AB % pc.p_BSb_A;        // 12 flops
    const Vec3 pd_SfSb_A = v_ASb - v_ASf;               //  3 flops
    pverr[0] = pc.oor * (~pd_SfSb_A * pc.p_SfSb_A);     //  6 flops
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

    // 30 flops for these two together.
    const Vec3 a_ASf = findStationInAAcceleration(s, allA_AB, m_mobod_F, 
                                                  pc.p_FSf_A, vc.wXwX_p_FSf_A);
    const Vec3 a_ASb = findStationInAAcceleration(s, allA_AB, m_mobod_B, 
                                                  pc.p_BSb_A, vc.wXwX_p_BSb_A);

    const Vec3 pdd_SfSb_A = a_ASb - a_ASf; // 2nd deriv in A,        3 flops
    paerr[0] = ~pdd_SfSb_A * pc.Cz_A + ~vc.pd_SfSb_A * vc.Czd_A; // 11 flops
}

// apply f=lambda*Pz to the contact point C of ball on body B,
//       -f         to point C (coincident point) of body F
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

    const Vec3 force_A = lambda * pc.Cz_A;     // 3 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BSb_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FSf_A, force_A, bodyForcesInA[m_mobod_F]);
}


// Implementation of virtuals required for nonholonomic constraints.

/* The rolling friction constraint says that the velocity of the material 
point of B at the contact point C = O - r Pz 
(at the bottom of the ball with center O and radius r), measured 
with respect to the plane body F, must be zero in the plane's Px and Py 
directions. Note that the contact point is not a station of body B but moves 
in the B frame. We have
    p_BC = p_BO - r Pz
    p_FC = p_FO - r Pz
    verr = v_FC = d/dt_F p_FC 
         = v_FO - r w_FB x Pz
where   v_FO = (v_AO-v_AF) - w_AF x (p_AO-p_AF)
and     w_FB = (w_AB-w_AF)

You have to differentiate verr carefully. Because Pz is fixed in F, the 
result is not the acceleration of the material point at C, but rather of
the moving contact point C.

    aerr = d/dt_F verr = a_FO - r d/dt_F (w_FB x Pz)
    d/dt_F (w_FB x Pz) = (d/dt_F w_FB) x Pz + w_FB x (d/dt_F Pz)
                       = b_FB x Pz  [2nd term dropped out!]
so  aerr = a_FC = a_FO - r b_FB x Pz

where a_FO = d/dt_F v_FO
and   b_FB = d/dt_F w_FB = (b_AB-b_AF) - w_AF x (w_AB-w_AF)
           = (b_AB-b_AF) - w_AF x w_AB

d/dt_F v_FO = d/dt_A v_FO - w_AF x v_FO
d/dt_A v_FO = (a_AO-a_AF) - b_AF x (p_AO-p_AF) - w_AF x (v_AO-v_AF)
w_AF x v_FO = w_AF x (v_AO-v_AF) - w_AF x (w_AF x (p_AO-p_AF))

so a_FO = (a_AO-a_AF) - b_AF x (p_AO-p_AF) - 2 w_AF x (v_AO-v_AF)
          + w_AF x (w_AF x (p_AO-p_AF))
and
aerr = a_FC = a_FO - r [b_AB-b_AF - w_AF x w_AB] x Pz

    -------------------------------------------
    verr = [v_FC_A . Px_A]
           [v_FC_A . Py_A]
    -------------------------------------------
    ~ 110 flops TODO should precalculate where possible
*/
#ifdef NOTDEF
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const OVERRIDE_11
{
    assert(allV_AB.size()==2 && constrainedU.size()==0 && verr.size()==2);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Real        rf    = params.m_radius_F;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rb    = params.m_radius_B;

    const Transform&  X_AF = getBodyTransformFromState(s, m_mobod_F);
    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const Rotation    R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform&  X_AB = getBodyTransformFromState(s, m_mobod_B);
    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)
    const Vec3        p_BC_A = p_BO_A - r*R_AP.z(); // bottom of ball (6 flops)

    Vec3 p_AC, v_AC;
    findStationInALocationVelocity(X_AB, V_AB, p_BC_A, p_AC, v_AC); // 15 flops

    // Calculate relative velocity of C in F, expressed in A.
    const Vec3 p_FC_A     = p_AC - p_AF; // measure C from Fo, in A (3 flops)
    const Vec3 p_FC_A_dot = v_AC - v_AF; // derivative in A (3 flops)
    const Vec3 v_FC_A = p_FC_A_dot - w_AF % p_FC_A; // deriv in F (12 flops)

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    verr[0] = ~R_AP.x()*v_FC_A; // 5 flops
    verr[1] = ~R_AP.y()*v_FC_A; // 5 flops
}

/*
    -------------------------------------------
    a_FC_A = a_FO_A - r b_FB_A x Pz_A
    aerr = [a_FC_A . Px_A]
           [a_FC_A . Py_A]
    -------------------------------------------
    ~ 200 flops TODO should precalculate where possible
*/

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const OVERRIDE_11
{
    assert(allA_AB.size()==2 && constrainedUDot.size()==0 && vaerr.size()==2);

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Real        rf    = params.m_radius_F;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rb    = params.m_radius_B;

    const Transform&  X_AF = getBodyTransformFromState(s, m_mobod_F);
    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_mobod_F);
    const SpatialVec& A_AF = allA_AB[m_mobod_F];
    const Vec3&       p_AF = X_AF.p();
    const Vec3&       w_AF = V_AF[0];
    const Vec3&       v_AF = V_AF[1];
    const Vec3&       b_AF = A_AF[0];
    const Vec3&       a_AF = A_AF[1];
    const Rotation    R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform&  X_AB = getBodyTransformFromState(s, m_mobod_B);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_mobod_B);
    const SpatialVec& A_AB = allA_AB[m_mobod_B];
    const Vec3&       w_AB = V_AB[0];
    const Vec3&       b_AB = A_AB[0];
    const Vec3        p_BO_A = X_AB.R() * p_BO; // re-express in A (15 flops)

    Vec3 p_AO, v_AO, a_AO;
    findStationInALocationVelocityAcceleration(X_AB, V_AB, A_AB, p_BO_A,
                                               p_AO, v_AO, a_AO); // 39 flops

    const Vec3 p_FO_A        = p_AO - p_AF; // measure O from Fo, in A (3 flops)
    const Vec3 p_FO_A_dot    = v_AO - v_AF; // deriv in A (3 flops)
    const Vec3 p_FO_A_dotdot = a_AO - a_AF; // 2nd deriv in A (3 flops)

    const Vec3 v_FO_A     = p_FO_A_dot    - w_AF % p_FO_A; // in F (12 flops)
    const Vec3 v_FO_A_dot = p_FO_A_dotdot -(b_AF % p_FO_A + w_AF % p_FO_A_dot);
                                                      // deriv in A (24 flops)
    const Vec3 a_FO_A = v_FO_A_dot - w_AF % v_FO_A; // now in F (12 flops)

    const Vec3 w_FB_A_dot = b_AB - b_AF;                        // 3 flops
    const Vec3 b_FB_A = w_FB_A_dot - w_AF % w_AB;               // 12 flops
    const Vec3 a_FC_A = a_FO_A - b_FB_A % (r*R_AP.z());         // 15 flops

    // Calculate these scalars using A-frame vectors, but the results are
    // measure numbers in [Px Py].
    vaerr[0] = ~R_AP.x()*a_FC_A; // 5 flops
    vaerr[1] = ~R_AP.y()*a_FC_A; // 5 flops
}

// apply f=lambda0*Px + lambda1*Py to the bottom point C of ball on B
//      -f           to point C (coincident point) of body F
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

    const Parameters& params = getParameters(s);
    const Vec3&       p_FSf = params.m_p_FSf;
    const Real        rf    = params.m_radius_F;
    const Vec3&       p_BSb = params.m_p_BSb;
    const Real        rb    = params.m_radius_B;

    const Transform& X_AF = getBodyTransformFromState(s, m_mobod_F);
    const Vec3&      p_AF = X_AF.p(); // p_AoFo_A
    const Rotation   R_AP = X_AF.R() * X_FP.R(); //R_AP=[Px Py Pz]_A (45 flops)

    const Transform& X_AB = getBodyTransformFromState(s, m_mobod_B);
    const Vec3&      p_AB = X_AB.p(); // p_AoBo_A

    const Vec3 p_BO_A = X_AB.R() * p_BO;   // re-express in A        (15 flops)
    const Vec3 p_BC_A = p_BO_A - r*R_AP.z(); // bottom of ball       ( 6 flops)
    const Vec3 p_AC   = p_AB + p_BC_A;     // measure from Ao        ( 3 flops)
    const Vec3 p_FC_A = p_AC - p_AF;       // measure C from Fo, in A (3 flops)

    const Vec3 force_A = lambda0*R_AP.x() + lambda1*R_AP.y(); // 9 flops

    // 30 flops for the next two calls
    addInStationInAForce(p_BC_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(p_FC_A, force_A, bodyForcesInA[m_mobod_F]);
}
#endif

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



