#ifndef SimTK_SIMBODY_CONSTRAINT_LINE_ON_LINE_CONTACT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_LINE_ON_LINE_CONTACT_IMPL_H_

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
Private implementation of Constraint::LineOnLineContact. **/

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
//                        LINE ON LINE CONTACT IMPL
//==============================================================================
class Constraint::LineOnLineContactImpl : public ConstraintImpl {
public:

// An object of this local type is used as the value of a Position-stage
// discrete state variable that allows sphere parameters to be changed for a 
// particular State.
struct Parameters {
    Parameters(const Transform& X_FEf, Real hf, 
               const Transform& X_BEb, Real hb)
    :   X_FEf(X_FEf), hf(hf), X_BEb(X_BEb), hb(hb)  {}
    Transform   X_FEf;  // edge Ef's frame
    Real        hf;     // edge Ef's half-length
    Transform   X_BEb;  // edge Eb's frame
    Real        hb;     // edge Eb's half-length
};

struct PositionCache {
    Transform   X_AC;       // z=n_A, x=df_A, y=n_A x df_A, Co=contact pt
    Vec3        p_FPf_A;    // edge centers, exp in A
    Vec3        p_BPb_A;
    Vec3        p_PfPb_A;   // vector from Pf to Pb, exp. in A
    UnitVec3    df_A;       // Ef direction, exp. in A
    UnitVec3    db_A;       // Eb direction, exp. in A
    Vec3        w_A;        // w_A = sense*(df_A X db_A)
    Real        sense;      // needed for derivs (1 or -1)
    Real        oos;        // 1/||w||=1/sin(angle between edges)
    UnitVec3    n_A;        // unit normal w_A/s
    Vec3        pXn;        // p_PfPb_A x n_A
    Real        tf, tb;     // location of Qf and Qb as edge parameters
    Vec3        p_AQf;      // closest points meas. from and exp. in A
    Vec3        p_AQb;
    Vec3        p_FCo_A;    // contact point measured from Fo, Bo, but in A
    Vec3        p_BCo_A;    // (for applying forces)
    bool        edgesAreParallel; // invalid; using different equations
};

struct VelocityCache {
    Vec3 dp_PfPb_A;         // v_APb - vAPf
    Vec3 ddf_A, ddb_A;      // w_AF x df_A, w_AB x db_A
    Vec3 dw_A;              // sense*(ddf_A%db_A + df_A%ddb_A)
    Vec3 dn_A;              // d/dt_A n_A
    Vec3 vF_BCo_A;          // vel of B's station at Co, meas. in F, exp. in A

    Real doos;              // d/dt oos

    Vec3 wXdp_FCo_A;        // w_AF x d/dt p_FCo  (F's station at Co)
    Vec3 wXdp_BCo_A;        // w_AB x d/dt p_BCo  (B's station at Co)
    Vec3 ddfXddb2;          // 2*(ddf_A % ddb_A)
    Vec3 wXddf_A, wXddb_A;  // w_AF x (w_AF x df_A), and B.
    Vec3 wXwX_p_FPf_A;      // w_AF x (w_AF x p_FPf_A)  [coriolis accel.]
    Vec3 wXwX_p_BPb_A;      // w_AB x (w_AB x p_BPb_A)  [coriolis accel.]

    Vec3 dCx_A,dCy_A,dCz_A; // contact frame derivatives in A
    Vec3 dCo_A;             // contact point derivative in A
};

explicit LineOnLineContactImpl(bool enforceRolling)
:   ConstraintImpl(1, enforceRolling?2:0, 0), 
    m_enforceRolling(enforceRolling),
    m_def_X_FEf(), m_def_X_BEb(), m_def_hf(NaN), m_def_hb(NaN) {}

LineOnLineContactImpl* clone() const override
{   return new LineOnLineContactImpl(*this); }

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

// Fill a PositionCache struct using the supplied Transforms rather than
// ones from the State. The PositionCache struct can be a throwaway, however; 
// it doesn't have to belong to a State
void calcPositionInfo(const State& state,
                      const Transform& X_AF, const Transform& X_AB,
                      PositionCache& pc) const;

// Make sure that the position cache info in this State is up to date, and
// return a reference to it.
const PositionCache& ensurePositionCacheRealized(const State& state) const;

const VelocityCache& getVelocityCache(const State& state) const;
VelocityCache& updVelocityCache(const State& state) const;

// Fill a VelocityCache struct using the supplied spatial velocities rather than
// ones from the State. Position calculations will be performed first if
// necessary and the State cache updated. The VelocityCache struct can be a 
// throwaway, however; it doesn't have to belong to a State.
void calcVelocityInfo(const State& state,
                      const SpatialVec& V_AF, const SpatialVec& V_AB,
                      VelocityCache& vc) const;

// Make sure that the velocity cache info in this State is up to date, and
// return a reference to it. Updates position info first if necessary.
const VelocityCache& ensureVelocityCacheRealized(const State& state) const;


// Implementation of virtuals required for holonomic constraints.

/* A "crossed edges" constraint serves several functions:
    - should work as a bilateral constraint
    - should use a signed-distance convention useful for unilateral contact
    - should provide geometric information suitable for triggering an
      event when contact slides off one of the edges, or the edges become
      parallel

There are two bodies F and B, with edges Ef fixed to F and Eb fixed to B.
For the constraint to be meaningful in unilateral contact, at least one of
the edges must be formed by the convex intersection of two faces of the 
polyhedral surface of a solid object (which does not have to be convex overall).
We assume Ef is defined that way by two faces of solid Sf and require the 
outward normals of those faces. That allows us to determine whether the contact
point on Eb is inside or outside of Sf, so that we can follow a useful sign
convention in which a positive perr indicates separation (Eb outside Sf) and 
negative perr indicates penetration (Eb inside Sf). Eb may also be produced by 
a convex intersection of faces or it may just be the edge of a flat plate or
even a thin wire; we don't need to know. Note that a "crossed edges" constraint
like this one should never occur when either Ef or Eb is a concave edge because
in that case the end points would have contacted before the edges.

Each edge is defined by three parameters: a direction unit vector d, an edge
midpoint P, and a half length h. In addition, edge Ef is assumed to be the
convex intersection of two faces of a solid and we are given the outward normals
n1 and n2 of those two faces.
so that we can make the constraint errors have meaningful
sign conventions.


Then the lines Lf and Lb containing the
edges are given by
    Lf = Pf + u*df, u in R
    Lb = Pb + v*db, v in R
and the edges themselves are the line segments -hf <= u <= hf and
-hb <= v <= hb.

Denote the lines' (not necessarily edges') closest points as Qf and Qb:
    Qf = Pf + tf*df
    Qb = Pb + tb*db
with tf and tb unknown scalars. These points are unique if the lines are not
parallel.


The edges are crossed iff |tf|<=hf and |tb|<=hb. 

We know that the vector between the closest 
points will be perpendicular to both lines, that is, parallel to w = df X db. 
Thus
    (1) Qf + r*n = Qb
where r is the distance between the closest points and n=w/||w||. By dotting 
both sides of (1) with n, we eliminate the terms involving df and db so can 
solve for r:
    (2) r = (Pb - Pf) . n
To make r a useful signed-distance value for unilateral contact, the user
must provide an outward normal for the solid for which Lf is an edge; we use
that to orient w so that positive r means separation and negative r means
penetration.

The normal contact condition is that the distance r between the two lines'
closest points should be zero:
    perr = p_PfPb . n,       p_PfPb=Pb-Pf, n=w/s, w = sense*(df X db), s=||w||

We'll take derivatives in A for convenience. Could use any frame.

    verr = d/dt_A perr = d/dt_A p_PfPb . n + p_PfPb . d/dt_A n

    pd_PfPb = d/dt_A p_PfPb = v_APb - v_APf
    ddf = d/dt_A df = w_AF x df
    ddb = d/dt_A db = w_AB x db

    dw = d/dt_A w= sense*(d/dt_A df X db + df X d/dt_A db)
                 = sense*(ddf x db + df x ddb)

    dws = dw/s

    ds = d/dt_F s= d/dt_F (w . w)^1/2 
                 = (w . dw) / s
                 = n . dw
    d/dt_A s^-1  = -s^-2 ds
                 = -(n . dw) / s^2
                 = -(n . dws) / s
   
    dn = d/dt_A (w/s) = (d/dt_A w) / s + w * (d/dt_A s^-1)
                      = dw/s - w * (n.dws)/s
                      = dws - (n.dws)*n


So
    verr = pd_PfPb . n + p_PfPb . dn

    pdd_PfPb = d/dt_A pd_PfPb = a_APb - a_APf

    d2df = d/dt_A ddf = b_AF x df + w_AF x ddf
    d2db = d/dt_A ddb = b_AB x db + w_AB x ddb

    d2w = d/dt_A dw = sense*(d/dt_A ddf x db + 2(ddf x ddb) + df x d/dt_A ddb)
        = sense * [d2df x db + 2(ddf x ddb) + df x d2db]
    ddws = d2w/s - dw * (n . dws) / s
         = d2w/s - (n . dws) * dws
    d2n = d/dt_A dn = ddws - (n.dws)*dn - (n.ddws + dn.dws)*n

    aerr = pdd_PfPb . n + 2(pd_PfPb . dn) + p_PfPb . d2n

(See comments below for friction.)
*/

/*  --------------------------------
    perr = p_PfPb . n
    --------------------------------
    n = w/||w||
    w = sign*(df X db)

    ~XX flops

Note that the give State does not necessarily correspond to the given X_AB
and constrainedQ so we can't cache the results.

TODO: there ought to be a flag saying whether the X_AB and q are from the 
state (which is common) so that the results could be saved.
*/
void calcPositionErrorsVirtual      
   (const State&                                    state,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    PositionCache tpc; // unfortunately just a temp
    calcPositionInfo(state, allX_AB[m_mobod_F], allX_AB[m_mobod_B], tpc);

    if (tpc.edgesAreParallel) {
        perr[0] = 0;
        return;
    }

    const Real r = ~tpc.p_PfPb_A * tpc.n_A;       // 5 flops
    //printf("r=%g, sign=%g\n", r, sign);

    perr[0] = r;
}


/*  --------------------------------
    verr = v_PfPb . n + p_PfPb . dn
    --------------------------------
    94 flops if position info is already available.

Note that we can't cache the velocity results here because we don't know that
V_AB and qdot were taken from the given State.
*/
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    const PositionCache& pc = ensurePositionCacheRealized(s);

    if (pc.edgesAreParallel) {
        pverr[0] = 0;
        return;
    }

    const UnitVec3& df_A     = pc.df_A;
    const UnitVec3& db_A     = pc.db_A;
    const UnitVec3& n_A      = pc.n_A;
    const Vec3&     p_PfPb_A = pc.p_PfPb_A;

    const SpatialVec& V_AF = allV_AB[m_mobod_F];
    const Vec3& w_AF = V_AF[0]; const Vec3& v_AF = V_AF[1];

    const SpatialVec& V_AB = allV_AB[m_mobod_B];
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];

    const Vec3 v_APf = v_AF + w_AF % pc.p_FPf_A;                // 12 flops
    const Vec3 v_APb = v_AB + w_AB % pc.p_BPb_A;                // 12
    const Vec3 dp_PfPb_A = v_APb - v_APf;                       //  3

    const Vec3 ddf_A = w_AF % df_A;                             //  9 flops
    const Vec3 ddb_A = w_AB % db_A;                             //  9 
    const Vec3 dw_A  = pc.sense*(ddf_A % db_A + df_A % ddb_A);  // 24
    const Vec3 dn_A  = pc.oos * (dw_A - (~n_A*dw_A)*n_A);       // 11

    pverr[0] = ~dp_PfPb_A * n_A + ~p_PfPb_A * dn_A;             // 11 flops
}

/*  --------------------------------------------------
    verr = v_PfPb . n + p_PfPb . dn
    aerr = a_PfPb . n + 2 (v_PfPb . dn) + p_PfPb . d2n
    --------------------------------------------------
    139 flops if position & velocity info is already available.
*/
void calcPositionDotDotErrorsVirtual      
   (const State&                                    state,  // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2&&constrainedQDotDot.size()==0&&paerr.size()==1);

    // Ensuring the velocity cache is realized also ensures that the position
    // cache has been realized.
    const VelocityCache& vc = ensureVelocityCacheRealized(state);
    const PositionCache& pc = getPositionCache(state);

    if (pc.edgesAreParallel) {
        paerr[0] = 0;
        return;
    }

    const SpatialVec& A_AF = allA_AB[m_mobod_F];
    const Vec3& b_AF = A_AF[0]; const Vec3& a_AF = A_AF[1];

    const SpatialVec& A_AB = allA_AB[m_mobod_B];
    const Vec3& b_AB = A_AB[0]; const Vec3& a_AB = A_AB[1];

    const UnitVec3& df_A  = pc.df_A;
    const UnitVec3& db_A  = pc.db_A;
    const UnitVec3& n_A   = pc.n_A;
    const Vec3&     dw_A  = vc.dw_A;
    const Vec3&     dn_A  = vc.dn_A;

    // 30 flops for these two together.
    const Vec3 a_APf = findStationInAAcceleration(state, allA_AB, m_mobod_F, 
                                                  pc.p_FPf_A, vc.wXwX_p_FPf_A);
    const Vec3 a_APb = findStationInAAcceleration(state, allA_AB, m_mobod_B, 
                                                  pc.p_BPb_A, vc.wXwX_p_BPb_A);
    const Vec3 d2p_PfPb_A = a_APb - a_APf;          //  3 flops

    const Vec3 d2df_A = b_AF % df_A + vc.wXddf_A;   // 12
    const Vec3 d2db_A = b_AB % db_A + vc.wXddb_A;   // 12
    const Vec3 d2w_A  = pc.sense*(d2df_A % db_A + vc.ddfXddb2 + df_A % d2db_A);
                                                    // 27

    const Real ndw = ~n_A * dw_A;                   //  5 flops
    const Vec3 ddw = d2w_A - (pc.oos*ndw)*dw_A;     //  7 
    const Vec3 d2n_A = pc.oos*(ddw-ndw*dn_A-(~n_A*ddw + ~dn_A*dw_A)*n_A);
                                                    // 26

    paerr[0] =   ~d2p_PfPb_A * n_A + 2*(~vc.dp_PfPb_A * dn_A)
               + ~pc.p_PfPb_A * d2n_A;              // 18 flops
}

// apply  f=lambda*n to body B's station at the contact point
//       -f          to body F's station at the contact point
// 33 flops if position info already available.
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==1&&bodyForcesInA.size()==2&&qForces.size()==0);
    const Real lambda = multipliers[0];
    const PositionCache& pc = ensurePositionCacheRealized(s);

    if (pc.edgesAreParallel)
        return;

    const Vec3 force_A = lambda * pc.n_A;    // 3 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BCo_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FCo_A, force_A, bodyForcesInA[m_mobod_F]);
}

// Implementation of virtuals required for nonholonomic constraints.

/* The rolling friction constraint says that the velocity in F of the material 
point of B coincident with the contact point Co, must be zero in the contact 
frame's Cx and Cy directions. 

We have already calculated p_FCo and p_BCo, the station locations on bodies
F and B of stations that are coincident with the contact point Co.
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

We are going to need the time derivative of the contact frame, including the
contact point Co, which is defined like this:
    Co  = (Qf + Qb)/2
Differentiating in the A frame, we have
    dCo = d/dt_A Co = (dQf + dQb)/2
The closest points Qf and Qb are defined 
    Qf  = Pf + tf * df
    Qb  = Pb + tb * db
so  dQf = dPf + dtf * df + tf * ddf
    dQb = dPb + dtb * db + tb * ddb

The scalars tf and tb are defined
    tf = -sense * (p_PfPb . (n x db)) / s
    tb = -sense * (p_PfPb . (n x df)) / s

We calculated the derivative of oos=1/s earlier, call it doos.
    dtf = -sense * d/dt_A [p_PfPb.(nxdb)] * oos + [p_PfPb.(nxdb)] * doos
    d/dt_A p_PfPb . (n x db) = pd_PfPb . (n x db) + p_PfPb . d/dt_A (nxdb)
    d/dt_A n x db = dn x db + n x ddb
and similarly we can get dtb.

We have 
     p_FCo = Co - Fo,        p_BCo = Co - Bo
    dp_FCo = dCo - dFo,     dp_BCo = dCo - dBo
where dFo = v_AF, dBo = v_AB.

    vA_FCo = v_AF + w_AF x p_FCo
    vA_BCo = v_AB + w_AB x p_BCo
    dvA_FCo = a_AF + b_AF x p_FCo + w_AF x dp_FCo
    dvA_BCo = a_AB + b_AB x p_BCo + w_AB x dp_BCo

We need to calculate aerr = d/dt verr. The components of verr are two scalars
so we can pick any frame in which to calculate the derivatives; we'll use A. So
from the chain rule we have:
(7) aerr = [d/dt_A vF_BCo . Cx + vF_BCo . d/dt_A Cx]
           [d/dt_A vF_BCo . Cy + vF_BCo . d/dt_A Cy]
         = [dvF_BCo . Cx + vF_BCo . Cxd]
           [dvF_BCo . Cy + vF_BCo . Cyd]

We already have vF_BCo, so we need these three quantities:
    Cxd     = d/dt_A Cx = d/dt_A df       = ddf
    Cyd     = d/dt_A Cy = d/dt_A (n x df) = dn x df + n x ddf
    dvF_BCo = d/dt_A vF_BCo               = dvA_BCo - dvA_FCo

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
    const override
{
    assert(allV_AB.size()==2 && constrainedU.size()==0 && verr.size()==2);

    const PositionCache& pc = ensurePositionCacheRealized(s);
    if (pc.edgesAreParallel) {
        verr[0] = verr[1] = 0;
        return;
    }


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
    // measure numbers in [Cx Cy]. This is the slip velocity of B at Co 
    // measured in F and expressed in C.
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
    const override
{
    assert(allA_AB.size()==2 && constrainedUDot.size()==0 && vaerr.size()==2);

    // Ensuring the velocity cache is realized also ensures that the position
    // cache has been realized.
    const VelocityCache& vc = ensureVelocityCacheRealized(s);
    const PositionCache& pc = getPositionCache(s);

    if (pc.edgesAreParallel) {
        vaerr[0] = vaerr[1] = 0;
        return;
    }

    const UnitVec3& Cx_A  = pc.X_AC.x();
    const UnitVec3& Cy_A  = pc.X_AC.y();
    const Vec3&     dCx_A = vc.dCx_A;
    const Vec3&     dCy_A = vc.dCy_A;

    const SpatialVec& A_AF = allA_AB[m_mobod_F];
    const Vec3& b_AF = A_AF[0]; const Vec3& a_AF = A_AF[1];

    const SpatialVec& A_AB = allA_AB[m_mobod_B];
    const Vec3& b_AB = A_AB[0]; const Vec3& a_AB = A_AB[1];

    const Vec3 dv_FCo_A = a_AF + b_AF % pc.p_FCo_A + vc.wXdp_FCo_A; // 15 flops
    const Vec3 dv_BCo_A = a_AB + b_AB % pc.p_BCo_A + vc.wXdp_BCo_A; // 15
    const Vec3 dvF_BCo_A = dv_BCo_A - dv_FCo_A;                     //  3

    vaerr[0] = ~dvF_BCo_A*Cx_A + ~vc.vF_BCo_A*dCx_A;                // 11 flops
    vaerr[1] = ~dvF_BCo_A*Cy_A + ~vc.vF_BCo_A*dCy_A;                // 11
}

// apply f=lambda0*Cx + lambda1*Cy to material point B/Co on B
//      -f                         to material point F/Co on F
// 39 flops if position info already available
void addInVelocityConstraintForcesVirtual
   (const State&                                    s,      // Stage::Velocity
    const Array_<Real>&                             multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const override
{
    assert(multipliers.size()==2 && mobilityForces.size()==0 
           && bodyForcesInA.size()==2);
    const Real lambda0 = multipliers[0], lambda1 = multipliers[1];

    const PositionCache& pc = ensurePositionCacheRealized(s);

    if (pc.edgesAreParallel)
        return;

    const UnitVec3& Cx_A  = pc.X_AC.x();
    const UnitVec3& Cy_A  = pc.X_AC.y();

    const Vec3 force_A = lambda0*Cx_A + lambda1*Cy_A; // 9 flops

    // 30 flops for the next two calls
    addInStationInAForce(pc.p_BCo_A, force_A, bodyForcesInA[m_mobod_B]);
    subInStationInAForce(pc.p_FCo_A, force_A, bodyForcesInA[m_mobod_F]);
}


SimTK_DOWNCAST(LineOnLineContactImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::LineOnLineContact;

// These can't be changed after construction.
ConstrainedBodyIndex    m_mobod_F;          // F (B1)
ConstrainedBodyIndex    m_mobod_B;          // B (B2)
const bool              m_enforceRolling;   // if so, add 2 constraints

Transform               m_def_X_FEf;    // default edge frame for Ef
Transform               m_def_X_BEb;    // default edge frame for Eb
Real                    m_def_hf;       // default edge Ef half length
Real                    m_def_hb;       // default edge Eb half length

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



