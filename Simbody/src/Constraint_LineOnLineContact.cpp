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
Constraint::LineOnLineContact, and its implementation class
Constraint::LineOnLineContactImpl. */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/Constraint_LineOnLineContact.h"

#include "Constraint_LineOnLineContactImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {


//==============================================================================
//                         LINE ON LINE CONTACT
//==============================================================================
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(Constraint::LineOnLineContact,
                                        Constraint::LineOnLineContactImpl,
                                        Constraint);

Constraint::LineOnLineContact::LineOnLineContact
   (MobilizedBody&      mobod_F,
    const Transform&    defaultEdgeFrameF,
    Real                defaultHalfLengthF,
    MobilizedBody&      mobod_B,
    const Transform&    defaultEdgeFrameB,
    Real                defaultHalfLengthB,
    bool                enforceRolling)
:   Constraint(new LineOnLineContactImpl(enforceRolling))
{
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSubsystem() && mobod_B.isInSubsystem(),
        "Constraint::LineOnLineContact","LineOnLineContact",
        "Both mobilized bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_APIARGCHECK_ALWAYS(mobod_F.isInSameSubsystem(mobod_B),
        "Constraint::LineOnLineContact","LineOnLineContact",
        "The two mobilized bodies to be connected must be in the same "
        "SimbodyMatterSubsystem.");

    mobod_F.updMatterSubsystem().adoptConstraint(*this);

    updImpl().m_mobod_F     = updImpl().addConstrainedBody(mobod_F);
    updImpl().m_mobod_B     = updImpl().addConstrainedBody(mobod_B);
    updImpl().m_def_X_FEf   = defaultEdgeFrameF;
    updImpl().m_def_hf      = defaultHalfLengthF;
    updImpl().m_def_X_BEb   = defaultEdgeFrameB;
    updImpl().m_def_hb      = defaultHalfLengthB;
}

Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setDefaultEdgeFrameF(const Transform& defaultEdgeFrameF) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_X_FEf = defaultEdgeFrameF;
    return *this;
}

Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setDefaultHalfLengthF(Real defaultHalfLengthF) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_hf = defaultHalfLengthF;
    return *this;
}

Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setDefaultEdgeFrameB(const Transform& defaultEdgeFrameB) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_X_BEb = defaultEdgeFrameB;
    return *this;
}

Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setDefaultHalfLengthB(Real defaultHalfLengthB) {
    getImpl().invalidateTopologyCache();
    updImpl().m_def_hb = defaultHalfLengthB;
    return *this;
}

const MobilizedBody& Constraint::LineOnLineContact::
getMobilizedBodyF() const {
    const LineOnLineContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_F);
}
const MobilizedBody& Constraint::LineOnLineContact::
getMobilizedBodyB() const {
    const LineOnLineContactImpl& impl = getImpl();
    return impl.getMobilizedBodyFromConstrainedBody(impl.m_mobod_B);
}

bool Constraint::LineOnLineContact::isEnforcingRolling() const
{   return getImpl().m_enforceRolling; }

const Transform& Constraint::LineOnLineContact::
getDefaultEdgeFrameF() const {return getImpl().m_def_X_FEf;}

Real Constraint::LineOnLineContact::
getDefaultHalfLengthF() const {return getImpl().m_def_hf;}

const Transform& Constraint::LineOnLineContact::
getDefaultEdgeFrameB() const {return getImpl().m_def_X_BEb;}

Real Constraint::LineOnLineContact::
getDefaultHalfLengthB() const {return getImpl().m_def_hb;}

const Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setEdgeFrameF(State& state, const Transform& edgeFrameF) const {
    getImpl().updParameters(state).X_FEf = edgeFrameF;
    return *this;
}

const Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setHalfLengthF(State& state, Real halfLengthF) const {
    getImpl().updParameters(state).hf = halfLengthF;
    return *this;
}

const Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setEdgeFrameB(State& state, const Transform& edgeFrameB) const {
    getImpl().updParameters(state).X_BEb = edgeFrameB;
    return *this;
}

const Constraint::LineOnLineContact& Constraint::LineOnLineContact::
setHalfLengthB(State& state, Real halfLengthB) const {
    getImpl().updParameters(state).hb = halfLengthB;
    return *this;
}
const Transform& Constraint::LineOnLineContact::
getEdgeFrameF(const State& state) const
{   return getImpl().getParameters(state).X_FEf; }
Real Constraint::LineOnLineContact::
getHalfLengthF(const State& state) const
{   return getImpl().getParameters(state).hf; }
const Transform& Constraint::LineOnLineContact::
getEdgeFrameB(const State& state) const
{   return getImpl().getParameters(state).X_BEb; }
Real Constraint::LineOnLineContact::
getHalfLengthB(const State& state) const
{   return getImpl().getParameters(state).hb; }

Real Constraint::LineOnLineContact::getPositionError(const State& s) const {
    Real perr;
    getImpl().getPositionErrors(s, 1, &perr);
    return perr;
}

Vec3 Constraint::LineOnLineContact::getVelocityErrors(const State& s) const {
    const LineOnLineContactImpl& impl = getImpl();
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

Vec3 Constraint::LineOnLineContact::getAccelerationErrors(const State& s) const {
    const LineOnLineContactImpl& impl = getImpl();
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

Vec3 Constraint::LineOnLineContact::getMultipliers(const State& s) const {
    const LineOnLineContactImpl& impl = getImpl();
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

Vec3 Constraint::LineOnLineContact::
findForceOnBodyBInG(const State& s) const {
    const LineOnLineContactImpl& impl = getImpl();
    if (impl.isDisabled(s))
        return Vec3(0);

    const Transform X_GC = findContactFrameInG(s);

    const Vec3 f_C = -getMultipliers(s); // watch sign convention
    return X_GC.R()*f_C;
}

Transform Constraint::LineOnLineContact::
findContactFrameInG(const State& s) const {
    const LineOnLineContactImpl& impl = getImpl();
    const LineOnLineContactImpl::PositionCache& pc =
        impl.ensurePositionCacheRealized(s);

    const MobilizedBody& mobod_A = impl.getAncestorMobilizedBody();
    if (mobod_A.isGround())
        return pc.X_AC; // == X_GC

    const Transform& X_GA = mobod_A.getBodyTransform(s);
    return X_GA * pc.X_AC;  // 63 flops
}

void Constraint::LineOnLineContact::
findClosestPointsInG(const State& s, Vec3& Qf_G, Vec3& Qb_G,
                     bool& linesAreParallel) const {
    const LineOnLineContactImpl& impl = getImpl();
    const LineOnLineContactImpl::PositionCache& pc =
        impl.ensurePositionCacheRealized(s);

    const MobilizedBody& mobod_A = impl.getAncestorMobilizedBody();
    if (mobod_A.isGround()) {
        Qf_G = pc.p_AQf;
        Qb_G = pc.p_AQb;
    } else {
        const Transform& X_GA = mobod_A.getBodyTransform(s);
        Qf_G = X_GA * pc.p_AQf; // 18 flops
        Qb_G = X_GA * pc.p_AQb; // 18 flops
    }
}

// The separation is the signed distance between the lines' closest points.
// Same as perr routine, but works when constraint is disabled.
Real Constraint::LineOnLineContact::
findSeparation(const State& s) const {
    const LineOnLineContactImpl::PositionCache& pc =
        getImpl().ensurePositionCacheRealized(s);

    const Real r = ~pc.p_PfPb_A * pc.n_A; // 5 flops
    return r;
}


//==============================================================================
//                        LINE ON LINE CONTACT IMPL
//==============================================================================

// The default parameters may be overridden by setting a discrete variable in
// the state, and we need a couple of cache entries to hold expensive
// computations. We allocate the state resources here.
void Constraint::LineOnLineContactImpl::
realizeTopologyVirtual(State& state) const {
    m_parametersIx = getMyMatterSubsystemRep().
        allocateDiscreteVariable(state, Stage::Position,
            new Value<Parameters>(Parameters(m_def_X_FEf, m_def_hf,
                                             m_def_X_BEb, m_def_hb)));

    m_posCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Position,
            new Value<PositionCache>());

    m_velCacheIx = getMyMatterSubsystemRep().
        allocateLazyCacheEntry(state, Stage::Velocity,
            new Value<VelocityCache>());
}

const Constraint::LineOnLineContactImpl::Parameters&
Constraint::LineOnLineContactImpl::
getParameters(const State& state) const {
    return Value<Parameters>::downcast
       (getMyMatterSubsystemRep().getDiscreteVariable(state,m_parametersIx));
}

Constraint::LineOnLineContactImpl::Parameters&
Constraint::LineOnLineContactImpl::
updParameters(State& state) const {
    return Value<Parameters>::updDowncast
       (getMyMatterSubsystemRep().updDiscreteVariable(state,m_parametersIx));
}

const Constraint::LineOnLineContactImpl::PositionCache&
Constraint::LineOnLineContactImpl::
getPositionCache(const State& state) const {
    return Value<PositionCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_posCacheIx));
}

Constraint::LineOnLineContactImpl::PositionCache&
Constraint::LineOnLineContactImpl::
updPositionCache(const State& state) const {
    return Value<PositionCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_posCacheIx));
}

const Constraint::LineOnLineContactImpl::VelocityCache&
Constraint::LineOnLineContactImpl::
getVelocityCache(const State& state) const {
    return Value<VelocityCache>::downcast
       (getMyMatterSubsystemRep().getCacheEntry(state,m_velCacheIx));
}

Constraint::LineOnLineContactImpl::VelocityCache&
Constraint::LineOnLineContactImpl::
updVelocityCache(const State& state) const {
    return Value<VelocityCache>::updDowncast
       (getMyMatterSubsystemRep().updCacheEntry(state,m_velCacheIx));
}

// This costs about 213 flops.
void Constraint::LineOnLineContactImpl::
calcPositionInfo(const State& state,
                 const Transform& X_AF, const Transform& X_AB,
                 PositionCache& pc) const
{
    const Parameters& params = getParameters(state);
    const UnitVec3&     df_F   = params.X_FEf.x();
    const UnitVec3&     db_B   = params.X_BEb.x();
    const UnitVec3&     sf_F   = params.X_FEf.z(); // outward normal
    const UnitVec3&     sb_B   = params.X_BEb.z();
    const Vec3&         p_FPf  = params.X_FEf.p();
    const Vec3&         p_BPb  = params.X_BEb.p();

    pc.df_A     = X_AF.R() * df_F;                  // 15 flops
    pc.db_A     = X_AB.R() * db_B;                  // 15
    pc.p_FPf_A  = X_AF.R() * p_FPf;                 // 15
    pc.p_BPb_A  = X_AB.R() * p_BPb;                 // 15
    const Vec3 p_APf = X_AF.p() + pc.p_FPf_A;       //  3
    const Vec3 p_APb = X_AB.p() + pc.p_BPb_A;       //  3
    pc.p_PfPb_A = p_APb - p_APf;                    //  3

    const UnitVec3 sf_A = X_AF.R() * sf_F;          // 15 flops
    const UnitVec3 sb_A = X_AB.R() * sb_B;          // 15
    const Vec3 wraw_A = pc.df_A % pc.db_A;          //  9

    // Use whichever of the outward normal directions gives a clearer signal.
    // We want w point out from F or into B.          ~12 flops
    const Real wsf = ~wraw_A*sf_A, wsb = ~wraw_A*sb_A;
    if (std::abs(wsf) > std::abs(wsb)) pc.sense = (wsf >= 0 ?  1 : -1);
    else                               pc.sense = (wsb >= 0 ? -1 :  1);

    pc.w_A = pc.sense * wraw_A;                     //  3 flops
    const Real sinTheta = pc.w_A.norm();            //~20

    //TODO: should do something better than this for parallel edges. Not clear
    // that an exception is best though since this might occur prior to an
    // assembly analysis that would fix the problem. Also consider if we're
    // doing event-driven contact we may need to catch the event of edges
    // going from crossed to parallel to turn off the edge/edge contact.

    if (sinTheta < SignificantReal) {               // 1 flop
        pc.edgesAreParallel = true;
        pc.oos = 1/SignificantReal; // just to avoid NaN-ing if sinTheta==0
    } else {
        pc.edgesAreParallel = false;
        pc.oos = 1/sinTheta;                        //~10
    }

    pc.n_A = UnitVec3(pc.w_A * pc.oos, true);       //  3

    pc.pXn = pc.p_PfPb_A % pc.n_A;                  //  9 flops

    const Real f = -pc.sense*pc.oos;                //  2 flops
    pc.tf = f*dot(pc.db_A, pc.pXn);                 //  6
    pc.tb = f*dot(pc.df_A, pc.pXn);                 //  6

    // Create the contact frame C. The origin should be half way between
    // Qf and Qb. The z axis is the contact normal n. x is along edge Ef, in
    // direction df. y = z X x = n X df.
    pc.p_AQf = p_APf + pc.tf * pc.df_A;             //  6 flops
    pc.p_AQb = p_APb + pc.tb * pc.db_A;             //  6
    const Vec3 p_ACo = 0.5*(pc.p_AQf + pc.p_AQb);   //  6
    pc.X_AC.updP() = p_ACo;

    // Since unit vector n is perpendicular to df (and db), n X df
    // is already a unit vector without normalizing.
    const UnitVec3 nXdf(pc.n_A % pc.df_A, true);    //  9 flops
    pc.X_AC.updR().setRotationFromUnitVecsTrustMe(pc.df_A, nXdf, pc.n_A);

    // Vectors from F and B body origins to contact point, expressed in A.
    // These are the station locations at which forces are applied.
    pc.p_FCo_A = p_ACo - X_AF.p();                  //  3 flops
    pc.p_BCo_A = p_ACo - X_AB.p();                  //  3
}


const Constraint::LineOnLineContactImpl::PositionCache&
Constraint::LineOnLineContactImpl::
ensurePositionCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_posCacheIx))
        return getPositionCache(s);
    PositionCache& pc = updPositionCache(s);

    const Transform& X_AF = getBodyTransformFromState(s, m_mobod_F);
    const Transform& X_AB = getBodyTransformFromState(s, m_mobod_B);

    calcPositionInfo(s, X_AF, X_AB, pc);

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_posCacheIx);
    return pc;
}


// This costs about 340 flops if position info has already been calculated,
// otherwise we also pay for ensurePositionCacheRealized().
void Constraint::LineOnLineContactImpl::
calcVelocityInfo(const State& state,
                 const SpatialVec& V_AF, const SpatialVec& V_AB,
                 VelocityCache& vc) const
{
    const PositionCache& pc = ensurePositionCacheRealized(state);
    if (pc.edgesAreParallel)
        return;

    const Vec3& w_AF = V_AF[0]; // handy abbreviations
    const Vec3& v_AF = V_AF[1];
    const Vec3& w_AB = V_AB[0];
    const Vec3& v_AB = V_AB[1];

    // These are d/dt_A p_FPf and d/dt_A p_BPb
    const Vec3 wX_p_FPf_A = w_AF % pc.p_FPf_A;                  //  9 flops
    const Vec3 wX_p_BPb_A = w_AB % pc.p_BPb_A;                  //  9
    const Vec3 v_APf = v_AF + wX_p_FPf_A;                       //  3
    const Vec3 v_APb = v_AB + wX_p_BPb_A;                       //  3
    vc.dp_PfPb_A = v_APb - v_APf;                               //  3

    vc.ddf_A = w_AF % pc.df_A;                                  //  9 flops
    vc.ddb_A = w_AB % pc.db_A;                                  //  9
    vc.dw_A = pc.sense*(vc.ddf_A % pc.db_A + pc.df_A % vc.ddb_A);//24
    vc.dn_A = pc.oos * (vc.dw_A - (~pc.n_A*vc.dw_A)*pc.n_A);    // 14

    // Calculate the velocity of B's material point (station) at Co,
    // measured in the F frame and expressed in A.
    const Vec3 vA_BCo_A = v_AB + w_AB % pc.p_BCo_A;             // 12 flops
    const Vec3 vA_FCo_A = v_AF + w_AF % pc.p_FCo_A;             // 12
    vc.vF_BCo_A = vA_BCo_A - vA_FCo_A;                          //  3

    // We have s=||w||, oos=1/s. We want doos = d/dt oos.
    vc.doos = -square(pc.oos) * dot(pc.n_A, vc.dw_A);           //  8 flops

    const Vec3 nXdb = pc.n_A % pc.db_A;                         //  9 flops
    const Vec3 nXdf = pc.n_A % pc.df_A;                         //  9
    const Vec3 d_nXdb = vc.dn_A % pc.db_A + pc.n_A % vc.ddb_A;  // 21
    const Vec3 d_nXdf = vc.dn_A % pc.df_A + pc.n_A % vc.ddf_A;  // 21
    const Real dtf = -pc.sense * (                              // 20
          pc.oos  * (~vc.dp_PfPb_A*nXdb + ~pc.p_PfPb_A*d_nXdb)
        + vc.doos * (~pc.p_PfPb_A*nXdb) );
    const Real dtb = -pc.sense * (                              // 20
          pc.oos  * (~vc.dp_PfPb_A*nXdf + ~pc.p_PfPb_A*d_nXdf)
        + vc.doos * (~pc.p_PfPb_A*nXdf) );

    const Vec3 dQf = v_APf + dtf * pc.df_A + pc.tf * vc.ddf_A;  // 12 flops
    const Vec3 dQb = v_APb + dtb * pc.db_A + pc.tb * vc.ddb_A;  // 12
    const Vec3 dCo = 0.5*(dQf + dQb);                           //  6
    const Vec3 dp_FCo = dCo - v_AF;                             //  3
    const Vec3 dp_BCo = dCo - v_AB;                             //  3

    vc.wXdp_FCo_A = w_AF % dp_FCo;                              //  9 flops
    vc.wXdp_BCo_A = w_AB % dp_BCo;                              //  9
    vc.ddfXddb2   = 2.*(vc.ddf_A % vc.ddb_A);                   // 12
    vc.wXddf_A    = w_AF % vc.ddf_A;                            //  9
    vc.wXddb_A    = w_AB % vc.ddb_A;                            //  9

    // These are the Coriolis accelerations of Pf and Pb, needed later.
    vc.wXwX_p_FPf_A = w_AF % wX_p_FPf_A;                        //  9 flops
    vc.wXwX_p_BPb_A = w_AB % wX_p_BPb_A;                        //  9

    // Record derivative of the contact frame.
    // We have Cx=df, Cz=n, Cy=n x df. Want derivatives in A.
    vc.dCx_A = vc.ddf_A;
    vc.dCz_A = vc.dn_A;
    vc.dCy_A = d_nXdf;
    vc.dCo_A = dCo;
}

const Constraint::LineOnLineContactImpl::VelocityCache&
Constraint::LineOnLineContactImpl::
ensureVelocityCacheRealized(const State& s) const {
    if (getMyMatterSubsystemRep().isCacheValueRealized(s, m_velCacheIx))
        return getVelocityCache(s);
    VelocityCache& vc = updVelocityCache(s);
    const SpatialVec& V_AF = getBodyVelocityFromState(s, m_mobod_F);
    const SpatialVec& V_AB = getBodyVelocityFromState(s, m_mobod_B);

    calcVelocityInfo(s, V_AF, V_AB, vc);

    getMyMatterSubsystemRep().markCacheValueRealized(s, m_velCacheIx);
    return vc;
}

void Constraint::LineOnLineContactImpl::
calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the lines' placements,
    // which might not be until Position stage.
    if (   stage != Stage::Position
        || !getMyMatterSubsystemRep().getShowDefaultGeometry())
        return;

    const Parameters&       params = getParameters(s);
    const Real              hf     = params.hf;
    const Real              hb     = params.hb;

    const PositionCache&    pc = ensurePositionCacheRealized(s);
    const Transform&        X_AC = pc.X_AC;

    const MobilizedBody&    mobod_A = getAncestorMobilizedBody();
    const Transform&        X_GA    = mobod_A.getBodyTransform(s);
    const Rotation&         R_GA    = X_GA.R();

    const Transform  X_GC = X_GA * X_AC;

    // Convert interesting stuff from A to G.
    const UnitVec3   df_G  = R_GA * X_AC.x();
    const UnitVec3   db_G  = R_GA * pc.db_A;
    const Vec3       p_GQf = X_GA * pc.p_AQf;
    const Vec3       p_GQb = X_GA * pc.p_AQb;
    const Vec3       half_Lf = hf * df_G;
    const Vec3       half_Lb = hb * db_G;

    const MobilizedBody& bodyF = getMobilizedBodyFromConstrainedBody(m_mobod_F);
    const MobilizedBody& bodyB = getMobilizedBodyFromConstrainedBody(m_mobod_B);

    const Transform& X_GF  = bodyF.getBodyTransform(s);
    const Transform& X_GB  = bodyB.getBodyTransform(s);
    const Transform  X_GEf = X_GF * params.X_FEf;
    const Transform  X_GEb = X_GB * params.X_BEb;


    // On body F draw a green line segment around the orange closest point.
    geom.push_back(DecorativeLine(p_GQf-half_Lf, p_GQf+half_Lf)
        .setColor(Green));
    geom.push_back(DecorativeFrame().setTransform(X_GEf)
        .setColor(Green*.9).setLineThickness(1).setScale(0.5)); // F color
    geom.push_back(DecorativePoint(p_GQf)
        .setColor(Orange).setLineThickness(2)); // B color

    // On body B draw an orange line segment around the green closest point.
    geom.push_back(DecorativeLine(p_GQb-half_Lb, p_GQb+half_Lb)
        .setColor(Orange));
    geom.push_back(DecorativeFrame().setTransform(X_GEb)
        .setColor(Orange*.9).setLineThickness(1).setScale(0.5)); // B color
    geom.push_back(DecorativePoint(p_GQb)
        .setColor(Green).setLineThickness(2)); // F color

    // Show the contact frame in red.
    geom.push_back(DecorativeFrame().setTransform(X_GC)
                   .setColor(Red));
}


} // namespace SimTK

