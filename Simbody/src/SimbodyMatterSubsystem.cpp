/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-9 Stanford University and the Authors.         *
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

/**@file
 *
 * Implementation of SimbodyMatterSubsystem, a concrete Subsystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/MobilizedBody.h"

#include "MobilizedBodyImpl.h"
#include "SimbodyMatterSubsystemRep.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {


/*static*/ bool 
SimbodyMatterSubsystem::isInstanceOf(const Subsystem& s) {
    return SimbodyMatterSubsystemRep::isA(s.getSubsystemGuts());
}
/*static*/ const SimbodyMatterSubsystem&
SimbodyMatterSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const SimbodyMatterSubsystem&>(s);
}
/*static*/ SimbodyMatterSubsystem&
SimbodyMatterSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<SimbodyMatterSubsystem&>(s);
}

const SimbodyMatterSubsystemRep& 
SimbodyMatterSubsystem::getRep() const {
    return dynamic_cast<const SimbodyMatterSubsystemRep&>(getSubsystemGuts());
}
SimbodyMatterSubsystemRep&       
SimbodyMatterSubsystem::updRep() {
    return dynamic_cast<SimbodyMatterSubsystemRep&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much 
// use except for making std::vector's, which require a default constructor 
// to be available.
SimbodyMatterSubsystem::SimbodyMatterSubsystem() 
  : Subsystem()
{
    adoptSubsystemGuts(new SimbodyMatterSubsystemRep());
    updRep().createGroundBody(); //TODO: handle this differently
}

SimbodyMatterSubsystem::SimbodyMatterSubsystem(MultibodySystem& mbs) 
  : Subsystem()
{
    adoptSubsystemGuts(new SimbodyMatterSubsystemRep());
    updRep().createGroundBody(); //TODO: handle this differently
    mbs.setMatterSubsystem(*this);
}

MobilizedBodyIndex SimbodyMatterSubsystem::adoptMobilizedBody(MobilizedBodyIndex parent, MobilizedBody& child) {
    return updRep().adoptMobilizedBody(parent,child);
}
const MobilizedBody& SimbodyMatterSubsystem::getMobilizedBody(MobilizedBodyIndex id) const {
    return getRep().getMobilizedBody(id);
}
MobilizedBody& SimbodyMatterSubsystem::updMobilizedBody(MobilizedBodyIndex id) {
    return updRep().updMobilizedBody(id);
}
const MobilizedBody::Ground& SimbodyMatterSubsystem::getGround() const {
    return getRep().getGround();
}
MobilizedBody::Ground& SimbodyMatterSubsystem::updGround() {
    return updRep().updGround();
}

bool SimbodyMatterSubsystem::getShowDefaultGeometry() const {
    return getRep().getShowDefaultGeometry();
}

void SimbodyMatterSubsystem::setShowDefaultGeometry(bool show) {
    updRep().setShowDefaultGeometry(show);
}


ConstraintIndex SimbodyMatterSubsystem::adoptConstraint(Constraint& child) {
    return updRep().adoptConstraint(child);
}
const Constraint& SimbodyMatterSubsystem::getConstraint(ConstraintIndex id) const {
    return getRep().getConstraint(id);
}
Constraint& SimbodyMatterSubsystem::updConstraint(ConstraintIndex id) {
    return updRep().updConstraint(id);
}

//TODO: should allow zero-length force arrays to stand for zeroes.
void SimbodyMatterSubsystem::calcAcceleration
   (const State&                state,
    const Vector&               appliedMobilityForces,
    const Vector_<SpatialVec>&  appliedBodyForces,
    Vector&                     udot,
    Vector_<SpatialVec>&        A_GB) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcAcceleration",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcAcceleration",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());

    Vector_<Vec3> appliedParticleForces; // TODO

    // Create a dummy acceleration cache to hold the result.
    const SBModelCache&    mc = getRep().getModelCache(state);
    const SBInstanceCache& ic = getRep().getInstanceCache(state);
    SBTreeAccelerationCache        tac;
    SBConstrainedAccelerationCache cac;
    tac.allocate(getRep().topologyCache, mc, ic);
    cac.allocate(getRep().topologyCache, mc, ic);

    Vector qdotdot; // unwanted return value
    Vector multipliers(getNMultipliers(state)); // unwanted return value
    Vector udotErr(getNUDotErr(state)); // unwanted return value

    getRep().calcLoopForwardDynamicsOperator(state, 
        appliedMobilityForces, appliedParticleForces, appliedBodyForces,
        tac, cac, udot, qdotdot, multipliers, udotErr);

    A_GB = tac.bodyAccelerationInGround;
}

//TODO: should allow zero-length force arrays to stand for zeroes.
void SimbodyMatterSubsystem::calcAccelerationIgnoringConstraints
   (const State&                state,
    const Vector&               appliedMobilityForces,
    const Vector_<SpatialVec>&  appliedBodyForces,
    Vector&                     udot, // output only; returns pres. accels
    Vector_<SpatialVec>&        A_GB) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcAccelerationIgnoringConstraints",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcAccelerationIgnoringConstraints",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());

    Vector netHingeForces(getNumMobilities()); // unwanted side effects
    Vector tau;
    Vector qdotdot;

    getRep().calcTreeAccelerations(state,
        appliedMobilityForces, appliedBodyForces,
        netHingeForces, A_GB, udot, qdotdot, tau);
}


void SimbodyMatterSubsystem::calcMInverseV(const State& s,
    const Vector&        v,
    Vector&              MinvV) const
{
    Vector_<SpatialVec> A_GB;
    getRep().calcMInverseF(s,v, A_GB, MinvV);
}



//==============================================================================
//                  CALC RESIDUAL FORCE IGNORING CONSTRAINTS
//==============================================================================
// This is inverse dynamics.
void SimbodyMatterSubsystem::calcResidualForceIgnoringConstraints
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    Vector&                    residualMobilityForces) const
{
    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==0 || appliedMobilityForces.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), getNumMobilities());
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForces.size()==0 || appliedBodyForces.size()==getNumBodies(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForces.size(), getNumBodies());
    SimTK_APIARGCHECK2_ALWAYS(
        knownUdot.size()==0 || knownUdot.size()==getNumMobilities(),
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d knownUdots but there are %d mobilities.",
        knownUdot.size(), getNumMobilities());

    residualMobilityForces.resize(getNumMobilities());

    Vector_<SpatialVec> A_GB(getNumBodies());
    getRep().calcTreeResidualForces(state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        A_GB, residualMobilityForces);
}

void SimbodyMatterSubsystem::calcMV(const State& s, 
    const Vector& v, 
    Vector& MV) const
{
    Vector_<SpatialVec> A_GB;
    getRep().calcMV(s,v, A_GB, MV);
}

void SimbodyMatterSubsystem::calcM(const State& s, Matrix& M) const 
{   getRep().calcM(s, M); }

void SimbodyMatterSubsystem::calcMInv(const State& s, Matrix& MInv) const 
{   getRep().calcMInv(s, MInv); }


// Note: the implementation methods that generate matrices do *not* require 
// contiguous storage, so we can just forward to them with no preliminaries.
void SimbodyMatterSubsystem::calcProjectedMInv(const State&   s,
                                               Matrix&        GMInvGt) const
{   getRep().calcGMInvGt(s, GMInvGt); }
void SimbodyMatterSubsystem::calcG(const State& s, Matrix& G) const 
{   getRep().calcPVA(s, true, true, true, G); }
void SimbodyMatterSubsystem::calcGTranspose(const State& s, Matrix& Gt) const 
{   getRep().calcPVATranspose(s, true, true, true, Gt); }
void SimbodyMatterSubsystem::calcPq(const State& s, Matrix& Pq) const 
{   getRep().calcPq(s,Pq); }
void SimbodyMatterSubsystem::calcPqTranspose(const State& s, Matrix& Pqt) const 
{   getRep().calcPqTranspose(s,Pqt); }

// OBSOLETE
void SimbodyMatterSubsystem::
calcPNInv(const State& s, Matrix& PNInv) const {
    return getRep().calcHolonomicConstraintMatrixPNInv(s,PNInv);
}

void SimbodyMatterSubsystem::
calcP(const State& s, Matrix& P) const {
    return getRep().calcHolonomicVelocityConstraintMatrixP(s,P);
}

void SimbodyMatterSubsystem::
calcPt(const State& s, Matrix& Pt) const {
    return getRep().calcHolonomicVelocityConstraintMatrixPt(s,Pt);
}



//==============================================================================
//                          MULTIPLY BY G TRANSPOSE
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necessary, call the
// implementation method to calculate f = ~G*lambda.
void SimbodyMatterSubsystem::
multiplyByGTranspose(const State&  s,
                     const Vector& lambda,
                     Vector&       f) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const SBInstanceCache& ic = rep.getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int m  = mHolo+mNonholo+mAccOnly;
    const int nu = rep.getNU(s);

    SimTK_ERRCHK2_ALWAYS(lambda.size() == m,
        "SimbodyMatterSubsystem::multiplyByGTranspose()",
        "Argument 'lambda' had length %d but should have the same length"
        " as the total number of active constraint equations m=%d.", 
        lambda.size(), m);

    f.resize(nu);
    if (nu==0) return;
    if (m==0) {f.setToZero(); return;}

    // Assume at first that both Vectors are contiguous.
    const Vector* clambda = &lambda;
    Vector*       cf      = &f;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_lambda, contig_f;

    if (!lambda.hasContiguousData()) {
        contig_lambda.resize(m); // contiguous memory
        contig_lambda(0, m) = lambda; // copy, prevent reallocation
        clambda = (const Vector*)&contig_lambda;
    }

    if (!f.hasContiguousData()) {
        contig_f.resize(nu); // contiguous memory
        cf = (Vector*)&contig_f;
        needToCopyBack = true;
    }

    rep.multiplyByPVATranspose(s, true, true, true, *clambda, *cf);
    if (needToCopyBack)
        f = *cf;
}




//==============================================================================
//                          MULTIPLY BY Pq TRANSPOSE
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necessary, call the
// implementation method to calculate fq = ~Pq*lambdap.
void SimbodyMatterSubsystem::
multiplyByPqTranspose(const State&  s,
                      const Vector& lambdap,
                      Vector&       fq) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const SBInstanceCache& ic = rep.getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = rep.getNQ(s);

    SimTK_ERRCHK2_ALWAYS(lambdap.size() == mp,
        "SimbodyMatterSubsystem::multiplyByPqTranspose()",
        "Argument 'lambdap' had length %d but should have had the same length"
        " as the number of active position (holonomic) constraint equations"
        " mp=%d.", lambdap.size(), mp);

    fq.resize(nq);
    if (nq==0) return;
    if (mp==0) {fq.setToZero(); return;}

    // Assume at first that both Vectors are contiguous.
    const Vector* clambdap = &lambdap;
    Vector*       cfq      = &fq;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_lambdap, contig_fq;

    if (!lambdap.hasContiguousData()) {
        contig_lambdap.resize(mp); // contiguous memory
        contig_lambdap(0, mp) = lambdap; // copy, prevent reallocation
        clambdap = (const Vector*)&contig_lambdap;
    }

    if (!fq.hasContiguousData()) {
        contig_fq.resize(nq); // contiguous memory
        cfq = (Vector*)&contig_fq;
        needToCopyBack = true;
    }

    rep.multiplyByPqTranspose(s, *clambdap, *cfq);
    if (needToCopyBack)
        fq = *cfq;
}



//==============================================================================
//                             MULTIPLY BY G
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necessary, call the
// implementation method.
void SimbodyMatterSubsystem::
multiplyByG(const State&  s,
            const Vector& ulike,
            const Vector& bias,
            Vector&       Gulike) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const SBInstanceCache& ic = rep.getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int m  = mHolo+mNonholo+mAccOnly;
    const int nu = rep.getNU(s);

    SimTK_ERRCHK2_ALWAYS(ulike.size() == nu,
        "SimbodyMatterSubsystem::multiplyByG()",
        "Argument 'ulike' had length %d but should have the same length"
        " as the total number of mobilities nu=%d.", ulike.size(), nu);

    SimTK_ERRCHK2_ALWAYS(bias.size() == m,
        "SimbodyMatterSubsystem::multiplyByG()",
        "Argument 'bias' had length %d but should have the same length"
        " as the total number of constraint equations m=%d.", bias.size(), m);

    Gulike.resize(m);
    if (m==0) return;

    // Assume at first that all Vectors are contiguous.
    const Vector* culike  = &ulike;
    const Vector* cbias   = &bias;
    Vector*       cGulike = &Gulike;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_ulike, contig_bias, contig_Gulike;

    if (!ulike.hasContiguousData()) {
        contig_ulike.resize(nu); // contiguous memory
        contig_ulike(0, nu) = ulike; // copy, prevent reallocation
        culike = (const Vector*)&contig_ulike;
    }
    if (!bias.hasContiguousData()) {
        contig_bias.resize(m); // contiguous memory
        contig_bias(0, m) = bias; // copy, prevent reallocation
        cbias = (const Vector*)&contig_bias;
    }
    if (!Gulike.hasContiguousData()) {
        contig_Gulike.resize(m); // contiguous memory
        cGulike = (Vector*)&contig_Gulike;
        needToCopyBack = true;
    }

    rep.multiplyByPVA(s, true, true, true, *cbias, *culike, *cGulike);
    if (needToCopyBack)
        Gulike = *cGulike;
}



//==============================================================================
//                      CALC BIAS FOR MULTIPLY BY G
//==============================================================================
// Here we just make sure that we have a contiguous array for the result and
// then call the implementation method.
void SimbodyMatterSubsystem::
calcBiasForMultiplyByG(const State& state,
                       Vector&      bias) const
{
    const SBInstanceCache& ic = getRep().getInstanceCache(state);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int m        = mHolo+mNonholo+mAccOnly;

    bias.resize(m);
    if (m==0) return;

    if (bias.hasContiguousData()) {
        getRep().calcBiasForMultiplyByPVA(state, true, true, true, bias);
    } else {
        Vector tmpbias(m); // contiguous
        getRep().calcBiasForMultiplyByPVA(state, true, true, true, tmpbias);
        bias = tmpbias;
    }
}




//==============================================================================
//                            MULTIPLY BY Pq
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necssary, call the
// implementation method.
void SimbodyMatterSubsystem::
multiplyByPq(const State&   s,
             const Vector&  qlike,
             const Vector&  biasp,
             Vector&        PqXqlike) const
{   
    const SimbodyMatterSubsystemRep& rep = getRep();
    const SBInstanceCache& ic = rep.getInstanceCache(s);

    // Problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = rep.getNQ(s);

    SimTK_ERRCHK2_ALWAYS(qlike.size() == nq,
        "SimbodyMatterSubsystem::multiplyByPq()",
        "Argument 'qlike' had length %d but should have been the same length"
        " as the total number of generalized coordinates nq=%d.", 
        qlike.size(), nq);

    SimTK_ERRCHK2_ALWAYS(biasp.size() == mp,
        "SimbodyMatterSubsystem::multiplyByPq()",
        "Argument 'biasp' had length %d but should have been the same length"
        " as the total number of position (holonomic) constraint"
        " equations mp=%d.", biasp.size(), mp);

    PqXqlike.resize(mp);
    if (mp==0) return;

    // Assume at first that all Vectors are contiguous.
    const Vector* cqlike    = &qlike;
    const Vector* cbiasp    = &biasp;
    Vector*       cPqXqlike = &PqXqlike;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_qlike, contig_biasp, contig_PqXqlike;

    if (!qlike.hasContiguousData()) {
        contig_qlike.resize(nq); // contiguous memory
        contig_qlike(0, nq) = qlike; // copy, prevent reallocation
        cqlike = (const Vector*)&contig_qlike;
    }
    if (!biasp.hasContiguousData()) {
        contig_biasp.resize(mp); // contiguous memory
        contig_biasp(0, mp) = biasp; // copy, prevent reallocation
        cbiasp = (const Vector*)&contig_biasp;
    }
    if (!PqXqlike.hasContiguousData()) {
        contig_PqXqlike.resize(mp); // contiguous memory
        cPqXqlike = (Vector*)&contig_PqXqlike;
        needToCopyBack = true;
    }

    rep.multiplyByPq(s, *cbiasp, *cqlike, *cPqXqlike);
    if (needToCopyBack)
        PqXqlike = *cPqXqlike;
}



//==============================================================================
//                       CALC BIAS FOR MULTIPLY BY Pq
//==============================================================================
// The bias term is the same for P as for Pq because you can view the 
// position error first derivative as either 
//           pverr = Pq * qdot + Pt
//      or   pverr = P  * u    + Pt
// since Pq = P*N^-1. Either way the bias term is just Pt (or c(t,q)).
void SimbodyMatterSubsystem::
calcBiasForMultiplyByPq(const State& state,
                        Vector&      biasp) const
{      
    const SBInstanceCache& ic = getRep().getInstanceCache(state);

    // Problem dimension.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;

    biasp.resize(mp);
    if (mp==0) return;

    if (biasp.hasContiguousData()) {
        // Just ask for P's bias term.
        getRep().calcBiasForMultiplyByPVA(state, true, false, false, biasp);
    } else {
        Vector tmpbias(mp); // contiguous
        getRep().calcBiasForMultiplyByPVA(state, true, false, false, tmpbias);
        biasp = tmpbias;
    }
}




//==============================================================================
//             CALC Gt -- OBSOLETE, use calcGTranspose()
//==============================================================================
void SimbodyMatterSubsystem::
calcGt(const State& s, Matrix& Gt) const {
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int mHolo    = rep.getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = rep.getNumNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = rep.getNumAccelerationOnlyConstraintEquationsInUse(s);
    const int m  = mHolo+mNonholo+mAccOnly;
    const int nu = rep.getNU(s);

    Gt.resize(nu,m);

    if (m==0 || nu==0)
        return;

    // Fill in all the columns of Gt
    rep.calcHolonomicVelocityConstraintMatrixPt(s, Gt(0,     0,          nu, mHolo));
    rep.calcNonholonomicConstraintMatrixVt     (s, Gt(0,   mHolo,        nu, mNonholo));
    rep.calcAccelerationOnlyConstraintMatrixAt (s, Gt(0, mHolo+mNonholo, nu, mAccOnly));
}




//==============================================================================
//                      CALC BODY ACCELERATION FROM UDOT
//==============================================================================
// Here we implement the zero-length udot, which is interpreted as an all-zero
// udot meaning that only coriolis accelerations contribute. Otherwise, we
// arrange to have contiguous input and output vectors to work with if the
// supplied arguments won't do, then invoke the implementation method.
void SimbodyMatterSubsystem::
calcBodyAccelerationFromUDot(const State&         state,
                             const Vector&        knownUDot,
                             Vector_<SpatialVec>& A_GB) const
{  
    // Interpret 0-length knownUDot as nu all-zero udots.
    if (knownUDot.size() == 0) {
        // Acceleration is just the coriolis acceleration.
        const SBTreeVelocityCache& vc = getRep().getTreeVelocityCache(state);
        const Array_<SpatialVec>& tca = vc.totalCoriolisAcceleration;
        const Vector_<SpatialVec> 
            AC_GB(tca.size(), (const Real*)tca.begin(), true); // shallow ref
        A_GB = AC_GB;
        return;
    }

    const int nu = getNumMobilities();

    SimTK_ERRCHK2_ALWAYS(knownUDot.size() == nu,
        "SimbodyMatterSubsystem::calcBodyAccelerationFromUDot()",
        "Length of knownUDot argument was %d but should have been either"
        " zero or the same as the number of mobilities nu=%d.\n", 
        knownUDot.size(), nu);

    const int nb = getNumBodies();
    A_GB.resize(nb);

    // If the arguments use contiguous memory we'll work in place, otherwise
    // we'll work in contiguous temporaries and copy back.

    Vector              udotspace; // allocate only if we need to
    Vector_<SpatialVec> Aspace;

    const Vector*        udotp;
    Vector_<SpatialVec>* Ap;

    if (knownUDot.hasContiguousData()) {
        udotp = &knownUDot;
    } else {
        udotspace.resize(nu); // contiguous memory
        udotspace(0, nu) = knownUDot; // prevent reallocation
        udotp = (const Vector*)&udotspace;
    }

    bool needToCopyBack = false;
    if (A_GB.hasContiguousData()) {
        Ap = &A_GB;
    } else {
        Aspace.resize(nb); // contiguous memory
        Ap = &Aspace;
        needToCopyBack = true;
    }

    getRep().calcBodyAccelerationFromUDot(state, *udotp, *Ap);

    if (needToCopyBack)
        A_GB = *Ap;
}



//==============================================================================
//                        MULTIPLY BY N, NInv, NDot
//==============================================================================
// These methods arrange for contiguous Vectors if necessary, then call the
// implementation method.
void SimbodyMatterSubsystem::multiplyByN
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   
    const bool inIsContig=in.hasContiguousData();
    const bool outIsContig=out.hasContiguousData();

    if (inIsContig && outIsContig) {
        getRep().multiplyByN(s,matrixOnRight,in,out); 
        return;
    }

    Vector inSpace, outSpace; // allocate if needed
    const Vector* inp  = inIsContig  ? &in  : (const Vector*)&inSpace;
    Vector*       outp = outIsContig ? &out : &outSpace;
    if (!inIsContig) {
        inSpace.resize(in.size());
        inSpace(0, in.size()) = in; // prevent reallocation
    }

    getRep().multiplyByN(s,matrixOnRight,*inp,*outp);

    if (!outIsContig)
        out = *outp;
}

void SimbodyMatterSubsystem::multiplyByNInv
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   
    const bool inIsContig=in.hasContiguousData();
    const bool outIsContig=out.hasContiguousData();

    if (inIsContig && outIsContig) {
        getRep().multiplyByNInv(s,matrixOnRight,in,out); 
        return;
    }

    Vector inSpace, outSpace; // allocate if needed
    const Vector* inp  = inIsContig  ? &in  : (const Vector*)&inSpace;
    Vector*       outp = outIsContig ? &out : &outSpace;
    if (!inIsContig) {
        inSpace.resize(in.size());
        inSpace(0, in.size()) = in; // prevent reallocation
    }

    getRep().multiplyByNInv(s,matrixOnRight,*inp,*outp); 

    if (!outIsContig)
        out = *outp;
}

void SimbodyMatterSubsystem::multiplyByNDot
   (const State& s, bool matrixOnRight, const Vector& in, Vector& out) const
{   
    const bool inIsContig=in.hasContiguousData();
    const bool outIsContig=out.hasContiguousData();

    if (inIsContig && outIsContig) {
        getRep().multiplyByNDot(s,matrixOnRight,in,out); 
        return;
    }

    Vector inSpace, outSpace; // allocate if needed
    const Vector* inp  = inIsContig  ? &in  : (const Vector*)&inSpace;
    Vector*       outp = outIsContig ? &out : &outSpace;
    if (!inIsContig) {
        inSpace.resize(in.size());
        inSpace(0, in.size()) = in; // prevent reallocation
    }

    getRep().multiplyByNDot(s,matrixOnRight,*inp,*outp); 

    if (!outIsContig)
        out = *outp;
}




//==============================================================================
//                            JACOBIAN METHODS
//==============================================================================

void SimbodyMatterSubsystem::multiplyBySystemJacobian
   (const State& s, const Vector& u, Vector_<SpatialVec>& Ju) const
{   getRep().multiplyBySystemJacobian(s,u,Ju); }

void SimbodyMatterSubsystem::multiplyBySystemJacobianTranspose
   (const State& s, const Vector_<SpatialVec>& F_G, Vector& f) const
{   getRep().multiplyBySystemJacobianTranspose(s,F_G,f); }

// Calculate J as an nb X nu matrix of SpatialVecs, by repeated calls 
// to J*u with u=0 except one u[i]=1.
void SimbodyMatterSubsystem::calcSystemJacobian
   (const State&            state,
    Matrix_<SpatialVec>&    J_G) const 
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    J_G.resize(nb,nu);
    Vector u(nu, Real(0));
    for (int j=0; j<nu; ++j) { 
        VectorView_<SpatialVec> col = J_G(j);
        u[j] = 1; rep.multiplyBySystemJacobian(state,u,col); u[j] = 0;
    }
}

// Alternate signature that returns a system Jacobian as a 6*nbod X nu Matrix 
// rather than as an nbod X nu matrix of spatial vectors. Note that we
// don't know whether the output matrix has contiguous rows, columns or
// neither.
void SimbodyMatterSubsystem::calcSystemJacobian
   (const State&            state,
    Matrix&                 J_G) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    J_G.resize(6*nb,nu); // we don't know how this is stored
    Vector u(nu, Real(0));
    Vector_<SpatialVec> Ju(nb); // temp Ju=J_G*u
    for (int j=0; j<nu; ++j) {
        u[j] = 1; rep.multiplyBySystemJacobian(state,u,Ju); u[j] = 0;
        VectorView col = J_G(j); // 6*nb long; maybe not contiguous!
        int nxt = 0; // index into col
        for (int b=0; b<nb; ++b) {
            const SpatialVec& V = Ju[b];
            for (int k=0; k<3; ++k) col[nxt++] = V[0][k]; // w
            for (int k=0; k<3; ++k) col[nxt++] = V[1][k]; // v
        }
    }
}

// We want v_GS = J_GS*u, the linear velocity of station S in Ground induced
// by the given generalized speeds. Station S is on body B and is given by
// the vector B-frame vector p_BoS. We can easily calculate V_GBo = J_GB*u,
// the spatial velocity of body B's origin Bo. Then 
//      v_GS = v_GBo + w_GB X p_BS_G
// where p_BS_G is p_BS re-expressed in Ground.
//
// Cost is 27 + 12*(nb+nu) flops. If nb ~= nu this is 27+24*nu.
Vec3 SimbodyMatterSubsystem::multiplyByStationJacobian
   (const State&         state,
    MobilizedBodyIndex   onBodyB,
    const Vec3&          p_BS,
    const Vector&        u) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_INDEXCHECK_ALWAYS(onBodyB, nb,
        "SimbodyMatterSubsystem::multiplyByStationJacobian()");

    SimTK_ERRCHK2_ALWAYS(u.size() == nu,
        "SimbodyMatterSubsystem::multiplyByStationJacobian()",
        "The supplied u-space Vector had length %d; expected %d.", 
            u.size(), nu);

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS); //15flops

    Vector_<SpatialVec> Ju(nb); // temp Ju=J_G*u
    rep.multiplyBySystemJacobian(state,u,Ju); // 12*(nb+nu) flops

    const SpatialVec& V_GB = Ju[onBodyB];
    const SpatialVec  V_GS = shiftVelocityBy(V_GB, p_BS_G); // 12 flops

    return V_GS[1]; // return linear velocity only
}

// We want f = ~J_GS*f_GS, the generalized forces produced by applying a
// translational force vector f_GS to station S. Station S is on body B and is 
// given by the vector B-frame vector p_BS (==p_BoS). We can easily calculate
// fb = ~J_GB*F_GB, and F_GB = [p_BS_G X f_GS, f_GS]. 
// Cost is 24 + 18nb + 11nu. If nb ~= nu, that's 24 + 29nu.
void SimbodyMatterSubsystem::multiplyByStationJacobianTranspose
   (const State&         state,
    MobilizedBodyIndex   onBodyB,
    const Vec3&          p_BS,
    const Vec3&          f_GS,
    Vector&              f) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_INDEXCHECK_ALWAYS(onBodyB, nb,
        "SimbodyMatterSubsystem::multiplyByStationJacobianTranspose()");

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS); //15flops

    // Need an array putting a spatial force on each body.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    F_G[onBodyB] = SpatialVec(p_BS_G % f_GS, f_GS); // 9 flops
    rep.multiplyBySystemJacobianTranspose(state,F_G,f); // 18nb+11nu flops
}

// Cost is 15+ 3*(9 + 18nb + 11nu) flops. If nb ~= nu, this is 42+87nu flops.
// Each subsequent multiply by J_GS*u would be 6*nu-3 flops.
void SimbodyMatterSubsystem::calcStationJacobian
   (const State&       state,
    MobilizedBodyIndex onBodyB,
    const Vec3&        p_BS,    // location of station S on B
    RowVector_<Vec3>&  J_GS) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);//15flops

    // Calculate J=dvdu where v is linear velocity of p_BS.
    // (This is three rows of J.)
    J_GS.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[onBodyB]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F (9 flops)
        rep.multiplyBySystemJacobianTranspose(state,F,col);// 18nb+11nu flops
        for (int r=0; r < nu; ++r) J_GS[r][i] = col[r]; 
        Fb[1][i] = 0;
    }
}

// Alternate signature that returns a station Jacobian as a 3 x nu Matrix 
// rather than as a row vector of Vec3s.
void SimbodyMatterSubsystem::calcStationJacobian
   (const State&       state,
    MobilizedBodyIndex onBodyB,
    const Vec3&        p_BS,    // location of station S on B
    Matrix&            J_GS) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dvdu where v is linear velocity of p_BS.
    // (This is three rows of J.)
    J_GS.resize(3,nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[onBodyB]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F
        rep.multiplyBySystemJacobianTranspose(state,F,col);
        J_GS[i] = ~col; // copy (TODO: avoid if rows are contiguous)
        Fb[1][i] = 0;
    }
}

// We want V_GA = J_GA*u, the spatial velocity of frame A in 
// Ground induced by the given generalized speeds. Frame A is fixed on body B 
// and would be given by the transform X_BA, except the result depends only 
// on A's origin position p_BA (==p_BoAo) because angular velocity is the same 
// for all frames fixed to the same body. We can easily calculate V_GB = J_GB*u,
// the spatial velocity at body B's origin Bo. Then 
//      V_GA = [w_GA, v_GAo] = [w_GB, v_GBo + w_GB X p_BA_G]
// where p_BA_G is p_BA re-expressed in Ground.
//
// Cost is 27 + 12*(nb+nu) flops. If nb ~= nu this is 27+24*nu.
SpatialVec SimbodyMatterSubsystem::multiplyByFrameJacobian
   (const State&         state,
    MobilizedBodyIndex   onBodyB,
    const Vec3&          p_BA,
    const Vector&        u) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_INDEXCHECK_ALWAYS(onBodyB, nb,
        "SimbodyMatterSubsystem::multiplyByFrameJacobian()");

    SimTK_ERRCHK2_ALWAYS(u.size() == nu,
        "SimbodyMatterSubsystem::multiplyByFrameJacobian()",
        "The supplied u-space Vector had length %d; expected %d.", 
            u.size(), nu);

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BA_G = mobod.expressVectorInGroundFrame(state, p_BA); //15flops

    Vector_<SpatialVec> Ju(nb); // temp Ju=J_G*u
    rep.multiplyBySystemJacobian(state,u,Ju); // 12*(nb+nu) flops

    const SpatialVec& V_GB = Ju[onBodyB];
    const SpatialVec  V_GA = shiftVelocityBy(V_GB, p_BA_G); // 12 flops

    return V_GA;
}

// We want f = ~J_GA*F_GA, the generalized forces produced by applying a spatial
// force vector F_GA=[m_G,f_GAo] at Ao, the origin of frame A, which is fixed
// to some body B. Frame A would be given by transform X_BA, but the result
// depends only on A's origin location p_BA (==p_BoAo) since a torque is the 
// same wherever it is applied. We can easily calculate
// fb = ~J_GB*F_GB, and F_GB = [m_G + p_BS_G X f_GAo, f_GAo]. 
// Cost is 27 + 18nb + 11nu. If nb ~= nu, that's 27 + 29nu.
void SimbodyMatterSubsystem::multiplyByFrameJacobianTranspose
   (const State&        state,
    MobilizedBodyIndex  onBodyB,
    const Vec3&         p_BA,
    const SpatialVec&   F_GAo,
    Vector&             f) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_INDEXCHECK_ALWAYS(onBodyB, nb,
        "SimbodyMatterSubsystem::multiplyByFrameJacobianTranspose()");

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BA_G = mobod.expressVectorInGroundFrame(state, p_BA); //15flops

    // Need an array putting a spatial force on each body.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    F_G[onBodyB] = SpatialVec(F_GAo[0] + p_BA_G % F_GAo[1], // 12 flops 
                              F_GAo[1]);
    rep.multiplyBySystemJacobianTranspose(state,F_G,f); // 18nb+11nu flops
}

// Cost is 42+ 6*(18nb + 11nu) flops. If nb ~= nu, this is 42+174nu flops.
// Each subsequent multiply by J_GF*u would be 12*nu-6 flops.
void SimbodyMatterSubsystem::calcFrameJacobian
   (const State&             state,
    MobilizedBodyIndex       onBodyB,
    const Vec3&              p_BFo,
    RowVector_<SpatialVec>&  J_GF) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BFo_G = mobod.expressVectorInGroundFrame(state, p_BFo);//15flops

    // Calculate J=dVdu where V is spatial velocity of F.
    // (This is six rows of J.)
    J_GF.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[onBodyB]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    // Rotational part.
    for (int i=0; i < 3; ++i) {
        Fb[0][i] = 1;
        rep.multiplyBySystemJacobianTranspose(state,F,col); // 18nb+11nu flops
        for (int r=0; r < nu; ++r) J_GF[r][0][i] = col[r]; 
        Fb[0][i] = 0;
    }
    // Translational part.
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BFo_G % Fb[1]; // r X F (9 flops)
        rep.multiplyBySystemJacobianTranspose(state,F,col); // 18nb+11nu flops
        for (int r=0; r < nu; ++r) J_GF[r][1][i] = col[r]; 
        Fb[1][i] = 0;
    }
}

// Alternate signature that returns a frame Jacobian as a 6 x nu Matrix 
// rather than as a row vector of SpatialVecs.
void SimbodyMatterSubsystem::calcFrameJacobian
   (const State&             state,
    MobilizedBodyIndex       onBodyB,
    const Vec3&              p_BFo,
    Matrix&                  J_GF) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();

    const MobilizedBody& mobod = rep.getMobilizedBody(onBodyB);
    const Vec3 p_BFo_G = mobod.expressVectorInGroundFrame(state, p_BFo);

    // Calculate J=dVdu where V is spatial velocity of F.
    // (This is six rows of J.)
    J_GF.resize(6,nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[onBodyB]; // the only one we'll change
    Vector col(nu); // temporary to hold column of J^T
    // Rotational part.
    for (int i=0; i < 3; ++i) {
        Fb[0][i] = 1;
        rep.multiplyBySystemJacobianTranspose(state,F,col);
        J_GF[i] = ~col; // copy row (TODO: avoid if rows are contiguous)
        Fb[0][i] = 0;
    }
    // Translational part.
    for (int i=0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BFo_G % Fb[1]; // r X F
        rep.multiplyBySystemJacobianTranspose(state,F,col);
        J_GF[i+3] = ~col; // copy row (TODO: avoid if rows are contiguous)
        Fb[1][i] = 0;
    }
}


//==============================================================================
//                              MISC OPERATORS
//==============================================================================

void SimbodyMatterSubsystem::calcCompositeBodyInertias
   (const State& s, Array_<SpatialInertia,MobilizedBodyIndex>& R) const
{   getRep().calcCompositeBodyInertias(s,R); }

void SimbodyMatterSubsystem::calcTreeEquivalentMobilityForces
   (const State& s, const Vector_<SpatialVec>& bodyForces, 
    Vector& mobForces) const
{   getRep().calcTreeEquivalentMobilityForces(s,bodyForces,mobForces); }

Real SimbodyMatterSubsystem::calcKineticEnergy(const State& s) const 
{   return getRep().calcKineticEnergy(s); }

void SimbodyMatterSubsystem::calcMobilizerReactionForces
   (const State& s, Vector_<SpatialVec>& forces) const 
{   getRep().calcMobilizerReactionForces(s, forces); }

void SimbodyMatterSubsystem::calcConstraintForcesFromMultipliers
   (const State& s, const Vector& lambda,
    Vector_<SpatialVec>& bodyForcesInG, Vector& mobilityForces) const
{   getRep().calcConstraintForcesFromMultipliers
                (s,lambda,bodyForcesInG,mobilityForces); }

void SimbodyMatterSubsystem::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    getRep().calcQDot(s, u, qdot);
}

void SimbodyMatterSubsystem::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    getRep().calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodyMatterSubsystem::getNumBodies()        const {return getRep().getNumBodies();}
int SimbodyMatterSubsystem::getNumMobilities()    const {return getRep().getNumMobilities();}
int SimbodyMatterSubsystem::getNumConstraints()   const {return getRep().getNumConstraints();}
int SimbodyMatterSubsystem::getNumParticles()     const {return getRep().getNumParticles();}

int SimbodyMatterSubsystem::getTotalQAlloc()    const {return getRep().getTotalQAlloc();}

// Modeling info.
void SimbodyMatterSubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { getRep().setUseEulerAngles(s,useAngles); }
void SimbodyMatterSubsystem::setMobilizerIsPrescribed(State& s, MobilizedBodyIndex body, bool prescribed) const
  { getRep().setMobilizerIsPrescribed(s,body,prescribed); }
void SimbodyMatterSubsystem::setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disabled) const
  { getRep().setConstraintIsDisabled(s,constraint,disabled); }
bool SimbodyMatterSubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodyMatterSubsystem::isMobilizerPrescribed(const State& s, MobilizedBodyIndex body) const
  { return getRep().isMobilizerPrescribed(s,body); }
bool SimbodyMatterSubsystem::isConstraintDisabled(const State& s, ConstraintIndex constraint) const
  { return getRep().isConstraintDisabled(s,constraint); }
void SimbodyMatterSubsystem::convertToEulerAngles(const State& inputState, State& outputState) const
  { return getRep().convertToEulerAngles(inputState, outputState); }
void SimbodyMatterSubsystem::convertToQuaternions(const State& inputState, State& outputState) const
  { return getRep().convertToQuaternions(inputState, outputState); }

int SimbodyMatterSubsystem::getNumQuaternionsInUse(const State& s) const {
    return getRep().getNumQuaternionsInUse(s);
}
bool SimbodyMatterSubsystem::isUsingQuaternion(const State& s, MobilizedBodyIndex body) const {
    return getRep().isUsingQuaternion(s, body);
}
QuaternionPoolIndex SimbodyMatterSubsystem::getQuaternionPoolIndex(const State& s, MobilizedBodyIndex body) const {
    return getRep().getQuaternionPoolIndex(s, body);
}
const SpatialVec&
SimbodyMatterSubsystem::getMobilizerCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
    return getRep().getMobilizerCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
    return getRep().getTotalCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getGyroscopicForce(const State& s, MobilizedBodyIndex body) const {
    return getRep().getGyroscopicForce(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getMobilizerCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getMobilizerCentrifugalForces(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getTotalCentrifugalForces(s,body);
}
const Vector& 
SimbodyMatterSubsystem::getAllParticleMasses(const State& s) const { 
    return getRep().getAllParticleMasses(s); 
}
Vector& SimbodyMatterSubsystem::updAllParticleMasses(State& s) const {
    return getRep().updAllParticleMasses(s); 
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleLocations(const State& s) const { 
    return getRep().getAllParticleLocations(s); 
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleVelocities(const State& s) const {
    return getRep().getAllParticleVelocities(s);
}

const Vector_<Vec3>& 
SimbodyMatterSubsystem::getAllParticleAccelerations(const State& s) const {
    return getRep().getAllParticleAccelerations(s);
}

void SimbodyMatterSubsystem::addInStationForce(const State& s, MobilizedBodyIndex body, const Vec3& stationInB, 
                                               const Vec3& forceInG, Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNumBodies());
    const Rotation& R_GB = getRep().getBodyTransform(s,body).R();
    bodyForces[body] += SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void SimbodyMatterSubsystem::realizeCompositeBodyInertias(const State& s) const {
    getRep().realizeCompositeBodyInertias(s);
}

void SimbodyMatterSubsystem::realizeArticulatedBodyInertias(const State& s) const {
    getRep().realizeArticulatedBodyInertias(s);
}

const SpatialInertia&
SimbodyMatterSubsystem::getCompositeBodyInertia(const State& s, MobilizedBodyIndex mbx) const {
    return getRep().getCompositeBodyInertias(s)[mbx]; // will lazy-evaluate if necessary
}

const ArticulatedInertia&
SimbodyMatterSubsystem::getArticulatedBodyInertia(const State& s, MobilizedBodyIndex mbx) const {
    return getRep().getArticulatedBodyInertias(s)[mbx]; // will lazy-evaluate if necessary
}

void SimbodyMatterSubsystem::addInBodyTorque(const State& s, MobilizedBodyIndex body, const Vec3& torqueInG,
                                             Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getRep().getNumBodies());
    bodyForces[body][0] += torqueInG; // no force
}
void SimbodyMatterSubsystem::addInMobilityForce(const State& s, MobilizedBodyIndex body, MobilizerUIndex which, Real d,
                                                Vector& mobilityForces) const 
{ 
    assert(mobilityForces.size() == getRep().getNumMobilities());
    UIndex uStart; int nu; getRep().findMobilizerUs(s,body,uStart,nu);
    assert(0 <= which && which < nu);
    mobilityForces[uStart+which] += d;
}

Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleLocations(State& s) const {
    return getRep().updAllParticleLocations(s);
}
Vector_<Vec3>& SimbodyMatterSubsystem::updAllParticleVelocities(State& s) const {
    return getRep().updAllParticleVelocities(s);
}

bool SimbodyMatterSubsystem::prescribe(State& s, Stage g) const {
    return getRep().prescribe(s,g);
}

bool SimbodyMatterSubsystem::projectQConstraints(State& s, Real consAccuracy, const Vector& yWeights,
                                                 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{ 
    return getRep().projectQConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts); 
}
bool SimbodyMatterSubsystem::projectUConstraints(State& s, Real consAccuracy, const Vector& yWeights,
                                                 const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{ 
    return getRep().projectUConstraints(s, consAccuracy, yWeights, ooTols, yErrest, opts); 
}

/// Calculate the total system mass.
/// TODO: this should be precalculated.
Real SimbodyMatterSubsystem::calcSystemMass(const State& s) const {
    Real mass = 0;
    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b)
        mass += getMobilizedBody(b).getBodyMassProperties(s).getMass();
    return mass;
}


// Return the location r_OG_C of the system mass center C, measured from the ground
// origin OG, and expressed in Ground. 
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterLocationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Transform&      X_GB    = getMobilizedBody(b).getBodyTransform(s);
        const Real            mb      = MB_OB_B.getMass();
        const Vec3            r_OG_CB = X_GB * MB_OB_B.getMassCenter();
        mass += mb;
        com  += mb * r_OG_CB; // weighted by mass
    }

    if (mass != 0) 
        com /= mass;

    return com;
}


// Return total system mass, mass center location measured from the Ground origin,
// and system inertia taken about the Ground origin, expressed in Ground.
MassProperties SimbodyMatterSubsystem::calcSystemMassPropertiesInGround(const State& s) const {
    Real    mass = 0;
    Vec3    com  = Vec3(0);
    Inertia I    = Inertia(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Transform&      X_GB    = getMobilizedBody(b).getBodyTransform(s);
        const MassProperties  MB_OG_G = MB_OB_B.calcTransformedMassProps(~X_GB);
        const Real            mb      = MB_OG_G.getMass();
        mass += mb;
        com  += mb * MB_OG_G.getMassCenter();
        I    += mb * MB_OG_G.getUnitInertia();
    }

    if (mass != 0)
        com /= mass;

    return MassProperties(mass, com, I);
}

// Return the system inertia matrix taken about the system center of mass,
// expressed in Ground.
Inertia SimbodyMatterSubsystem::calcSystemCentralInertiaInGround(const State& s) const {
    const MassProperties M_OG_G = calcSystemMassPropertiesInGround(s);
    return M_OG_G.calcCentralInertia();
}


// Return the velocity V_G_C = d/dt r_OG_C of the system mass center C in the Ground frame G,
// expressed in G.
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterVelocityInGround(const State& s) const {
    Real    mass = 0;
    Vec3    comv = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 v_G_CB = getMobilizedBody(b).findStationVelocityInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        comv += mb * v_G_CB; // weighted by mass
    }

    if (mass != 0) 
        comv /= mass;

    return comv;
}

// Return the acceleration A_G_C = d^2/dt^2 r_OG_C of the system mass center C in
// the Ground frame G, expressed in G.
Vec3 SimbodyMatterSubsystem::calcSystemMassCenterAccelerationInGround(const State& s) const {
    Real    mass = 0;
    Vec3    coma = Vec3(0);

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MassProperties& MB_OB_B = getMobilizedBody(b).getBodyMassProperties(s);
        const Vec3 a_G_CB = getMobilizedBody(b).findStationAccelerationInGround(s, MB_OB_B.getMassCenter());
        const Real mb     = MB_OB_B.getMass();

        mass += mb;
        coma += mb * a_G_CB; // weighted by mass
    }

    if (mass != 0) 
        coma /= mass;

    return coma;
}

// Return the momentum of the system as a whole (angular, linear) measured
// in the ground frame, taken about the ground origin and expressed in ground.
// (The linear component is independent of the "about" point.)
SpatialVec SimbodyMatterSubsystem::calcSystemMomentumAboutGroundOrigin(const State& s) const {
    SpatialVec mom(Vec3(0), Vec3(0));
    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const SpatialVec mom_CB_G = getMobilizedBody(b).calcBodyMomentumAboutBodyMassCenterInGround(s);
        const Vec3&      Iw = mom_CB_G[0];
        const Vec3&      mv = mom_CB_G[1];
        const Vec3       r = getMobilizedBody(b).findMassCenterLocationInGround(s);
        mom[0] += (Iw + r % mv); // add central angular momentum plus contribution from mass center location
        mom[1] += mv;            // just add up central linear momenta
    }
    return mom;
}

// Return the momentum of the system as a whole (angular, linear) measured
// in the ground frame, taken about the current system center of mass
// location and expressed in ground.
// (The linear component is independent of the "about" point.)
SpatialVec SimbodyMatterSubsystem::calcSystemCentralMomentum(const State& s) const {
    SpatialVec mom(Vec3(0), Vec3(0));
    Real mtot = 0;  // system mass
    Vec3 com(0);    // system mass center

    for (MobilizedBodyIndex b(1); b < getNumBodies(); ++b) {
        const MobilizedBody& mobod = getMobilizedBody(b);
        const Real m    = mobod.getBodyMass(s);
        const Vec3 CB_G = mobod.findMassCenterLocationInGround(s);
        mtot += m;
        com  += m * CB_G;
        const SpatialVec mom_CB_G = mobod.calcBodyMomentumAboutBodyMassCenterInGround(s);
        const Vec3&      Iw = mom_CB_G[0];
        const Vec3&      mv = mom_CB_G[1];
        mom[0] += (Iw + CB_G % mv); // add central angular momentum plus contribution from mass center location
        mom[1] += mv;               // just add up central linear momenta
    }
    if (mtot != 0)
        com /= mtot;

    // Shift momentum from ground origin to system COM (only angular affected).
    mom[0] -= com % mom[1];
    return mom;
}

} // namespace SimTK

