/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
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
    return static_cast<const SimbodyMatterSubsystem&>(s);
}
/*static*/ SimbodyMatterSubsystem&
SimbodyMatterSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return static_cast<SimbodyMatterSubsystem&>(s);
}

const SimbodyMatterSubsystemRep& 
SimbodyMatterSubsystem::getRep() const {
    return SimTK_DYNAMIC_CAST_DEBUG<const SimbodyMatterSubsystemRep&>(getSubsystemGuts());
}
SimbodyMatterSubsystemRep&       
SimbodyMatterSubsystem::updRep() {
    return SimTK_DYNAMIC_CAST_DEBUG<SimbodyMatterSubsystemRep&>(updSubsystemGuts());
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


ConstraintIndex SimbodyMatterSubsystem::
adoptConstraint(Constraint& child) {return updRep().adoptConstraint(child);}
const Constraint& SimbodyMatterSubsystem::
getConstraint(ConstraintIndex id) const {return getRep().getConstraint(id);}
Constraint& SimbodyMatterSubsystem::
updConstraint(ConstraintIndex id) {return updRep().updConstraint(id);}


UnilateralContactIndex SimbodyMatterSubsystem::
adoptUnilateralContact(UnilateralContact* child)
{   return updRep().adoptUnilateralContact(child); }
int SimbodyMatterSubsystem::
getNumUnilateralContacts() const 
{   return getRep().getNumUnilateralContacts(); }
const UnilateralContact& SimbodyMatterSubsystem::
getUnilateralContact(UnilateralContactIndex ix) const
{   return getRep().getUnilateralContact(ix); }
UnilateralContact& SimbodyMatterSubsystem::
updUnilateralContact(UnilateralContactIndex ix)
{   return updRep().updUnilateralContact(ix); }


StateLimitedFrictionIndex SimbodyMatterSubsystem::
adoptStateLimitedFriction(StateLimitedFriction* child)
{   return updRep().adoptStateLimitedFriction(child); }
int SimbodyMatterSubsystem::
getNumStateLimitedFrictions() const
{   return getRep().getNumStateLimitedFrictions(); }
const StateLimitedFriction& SimbodyMatterSubsystem::
getStateLimitedFriction(StateLimitedFrictionIndex ix) const
{   return getRep().getStateLimitedFriction(ix); }
StateLimitedFriction& SimbodyMatterSubsystem::
updStateLimitedFriction(StateLimitedFrictionIndex ix)
{   return updRep().updStateLimitedFriction(ix); }


//==============================================================================
//                            CALC ACCELERATION
//==============================================================================
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



//==============================================================================
//                    CALC ACCELERATION IGNORING CONSTRAINTS
//==============================================================================
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
    Array_<SpatialVec,MobilizedBodyIndex> abForcesZ(getNumBodies());   
    Array_<SpatialVec,MobilizedBodyIndex> abForcesZPlus(getNumBodies());   
    Vector tau;
    Vector qdotdot;

    const SBDynamicsCache& dc = getRep().getDynamicsCache(state);

    getRep().calcTreeAccelerations(state,
        appliedMobilityForces, appliedBodyForces, dc.presUDotPool,
        netHingeForces, abForcesZ, abForcesZPlus, 
        A_GB, udot, qdotdot, tau);
}



//==============================================================================
//                  CALC RESIDUAL FORCE IGNORING CONSTRAINTS
//==============================================================================
// This is inverse dynamics.
// This just checks the arguments, arranges for contiguous vectors to work
// with if necessary, and then calls the implementation method.
void SimbodyMatterSubsystem::calcResidualForceIgnoringConstraints
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForcesInG,
    const Vector&              knownUdot,
    Vector&                    residualMobilityForces) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nu = rep.getNU(state);

    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==0 || appliedMobilityForces.size()==nu,
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), nu);
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForcesInG.size()==0 || appliedBodyForcesInG.size()==nb,
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForcesInG.size(), nb);
    SimTK_APIARGCHECK2_ALWAYS(
        knownUdot.size()==0 || knownUdot.size()==nu,
        "SimbodyMatterSubsystem", "calcResidualForceIgnoringConstraints",
        "Got %d knownUdots but there are %d mobilities.",
        knownUdot.size(), nu);

    residualMobilityForces.resize(nu);

    // Assume at first that all Vectors are contiguous.
    const Vector*               cmobForces  = &appliedMobilityForces;
    const Vector_<SpatialVec>*  cbodyForces = &appliedBodyForcesInG;
    const Vector*               cudot       = &knownUdot;
    Vector*                     cresid      = &residualMobilityForces;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_mobForces, contig_udot, contig_resid;
    Vector_<SpatialVec> contig_bodyForces;

    if (!appliedMobilityForces.hasContiguousData()) {
        contig_mobForces.resize(nu); // contiguous memory
        contig_mobForces(0, nu) = appliedMobilityForces; // copy, no reallocation
        cmobForces = (const Vector*)&contig_mobForces;
    }
    if (!appliedBodyForcesInG.hasContiguousData()) {
        contig_bodyForces.resize(nb); // contiguous memory
        contig_bodyForces(0, nb) = appliedBodyForcesInG; // copy, no reallocation
        cbodyForces = (const Vector_<SpatialVec>*)&contig_bodyForces;
    }
    if (!knownUdot.hasContiguousData()) {
        contig_udot.resize(nu); // contiguous memory
        contig_udot(0, nu) = knownUdot; // copy, no reallocation
        cudot = (const Vector*)&contig_udot;
    }
    if (!residualMobilityForces.hasContiguousData()) {
        contig_resid.resize(nu); // contiguous memory
        cresid = (Vector*)&contig_resid;
        needToCopyBack = true;
    }

    Vector_<SpatialVec> A_GB(nb); // temp for unwanted result
    rep.calcTreeResidualForces(state,
        *cmobForces, *cbodyForces, *cudot,
        A_GB, *cresid);

    if (needToCopyBack)
        residualMobilityForces = *cresid;
}



//==============================================================================
//                          CALC RESIDUAL FORCE 
//==============================================================================
// This is inverse dynamics with constraints.
void SimbodyMatterSubsystem::calcResidualForce
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForcesInG,
    const Vector&              knownUdot,
    const Vector&              knownLambda,
    Vector&                    residualMobilityForces) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nu = rep.getNU(state);
    const int m  = rep.getNMultipliers(state);

    SimTK_APIARGCHECK2_ALWAYS(
        appliedMobilityForces.size()==0 || appliedMobilityForces.size()==nu,
        "SimbodyMatterSubsystem", "calcResidualForce",
        "Got %d appliedMobilityForces but there are %d mobilities.",
        appliedMobilityForces.size(), nu);
    SimTK_APIARGCHECK2_ALWAYS(
        appliedBodyForcesInG.size()==0 || appliedBodyForcesInG.size()==nb,
        "SimbodyMatterSubsystem", "calcResidualForce",
        "Got %d appliedBodyForces but there are %d bodies (including Ground).",
        appliedBodyForcesInG.size(), nb);
    SimTK_APIARGCHECK2_ALWAYS(
        knownUdot.size()==0 || knownUdot.size()==nu,
        "SimbodyMatterSubsystem", "calcResidualForce",
        "Got %d knownUdots but there are %d mobilities.",
        knownUdot.size(), nu);
    SimTK_APIARGCHECK2_ALWAYS(
        knownLambda.size()==0 || knownLambda.size()==m,
        "SimbodyMatterSubsystem", "calcResidualForce",
        "Got %d knownLambdas but there are %d constraint equations.",
        knownUdot.size(), m);

    if (knownLambda.size() == 0) { // no constraint forces
        // Call above method instead.
        calcResidualForceIgnoringConstraints(state,
            appliedMobilityForces, appliedBodyForcesInG, knownUdot,
            residualMobilityForces);
        return; 
    }

    // There are some lambdas, so calculate the forces they produce, with
    // the result going in newly-allocated contiguous storage. We have to
    // negate lambda to make the constraint forces have the sign of applied
    // forces; we'll do that into a contiguous Vector also.
    Vector_<SpatialVec> bodyForcesInG(nb);
    Vector              mobilityForces(nu);
    Vector              negLambda = -knownLambda;
    rep.calcConstraintForcesFromMultipliers(state, negLambda,
        bodyForcesInG, mobilityForces);

    // Now add in the applied forces, and call the unconstrained routine.
    if (appliedBodyForcesInG.size())
        bodyForcesInG  += appliedBodyForcesInG;
    if (appliedMobilityForces.size())
        mobilityForces += appliedMobilityForces;

    calcResidualForceIgnoringConstraints(state,
        mobilityForces, bodyForcesInG, knownUdot,
        residualMobilityForces);
}



//==============================================================================
//                               MULTIPLY BY M
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necessary, call the
// implementation method to calculate f = M*a.
void SimbodyMatterSubsystem::multiplyByM(const State&  state, 
                                         const Vector& a, 
                                         Vector&       Ma) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nu = rep.getNU(state);

    SimTK_ERRCHK2_ALWAYS(a.size() == nu,
        "SimbodyMatterSubsystem::multiplyByM()",
        "Argument 'a' had length %d but should have the same length"
        " as the number of mobilities (generalized speeds u) %d.", 
        a.size(), nu);

    Ma.resize(nu);
    if (nu==0) return;

    // Assume at first that both Vectors are contiguous.
    const Vector* ca    = &a;
    Vector*       cMa   = &Ma;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_a, contig_Ma;

    if (!a.hasContiguousData()) {
        contig_a.resize(nu); // contiguous memory
        contig_a(0, nu) = a; // copy, prevent reallocation
        ca = (const Vector*)&contig_a;
    }

    if (!Ma.hasContiguousData()) {
        contig_Ma.resize(nu); // contiguous memory
        cMa = (Vector*)&contig_Ma;
        needToCopyBack = true;
    }

    rep.multiplyByM(state, *ca, *cMa);

    if (needToCopyBack)
        Ma = *cMa;
}



//==============================================================================
//                             MULTIPLY BY M INV
//==============================================================================
// Check arguments, copy in/out of contiguous Vectors if necessary, call the
// implementation method to calculate a = M^-1*v.
void SimbodyMatterSubsystem::multiplyByMInv(const State&    state,
                                            const Vector&   v,
                                            Vector&         MInvV) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nu = rep.getNU(state);

    SimTK_ERRCHK2_ALWAYS(v.size() == nu,
        "SimbodyMatterSubsystem::multiplyByMInv()",
        "Argument 'v' had length %d but should have the same length"
        " as the number of mobilities (generalized speeds u) %d.", 
        v.size(), nu);

    MInvV.resize(nu);
    if (nu==0) return;

    // Assume at first that both Vectors are contiguous.
    const Vector* cv    = &v;
    Vector*       cMInvV   = &MInvV;
    bool needToCopyBack = false;

    // We'll allocate these or not as needed.
    Vector contig_v, contig_MInvV;

    if (!v.hasContiguousData()) {
        contig_v.resize(nu); // contiguous memory
        contig_v(0, nu) = v; // copy, prevent reallocation
        cv = (const Vector*)&contig_v;
    }

    if (!MInvV.hasContiguousData()) {
        contig_MInvV.resize(nu); // contiguous memory
        cMInvV = (Vector*)&contig_MInvV;
        needToCopyBack = true;
    }

    rep.multiplyByMInv(state, *cv, *cMInvV);

    if (needToCopyBack)
        MInvV = *cMInvV;
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

void SimbodyMatterSubsystem::
solveForConstraintImpulses(const State&     state,
                           const Vector&    deltaV,
                           Vector&          impulse) const
{   getRep().solveForConstraintImpulses(state,deltaV,impulse); }


void SimbodyMatterSubsystem::calcG(const State& s, Matrix& G) const 
{   getRep().calcPVA(s, true, true, true, G); }
void SimbodyMatterSubsystem::calcGTranspose(const State& s, Matrix& Gt) const 
{   getRep().calcPVATranspose(s, true, true, true, Gt); }
void SimbodyMatterSubsystem::calcPq(const State& s, Matrix& Pq) const 
{   getRep().calcPq(s,Pq); }
void SimbodyMatterSubsystem::calcPqTranspose(const State& s, Matrix& Pqt) const 
{   getRep().calcPqTranspose(s,Pqt); }

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
//                    CALC BIAS FOR ACCELERATION CONSTRAINTS
//==============================================================================
// Here we just make sure that we have a contiguous array for the result and
// then call the implementation method.
void SimbodyMatterSubsystem::
calcBiasForAccelerationConstraints(const State& state,
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
        getRep().calcBiasForAccelerationConstraints(state,true,true,true,bias);
    } else {
        Vector tmpbias(m); // contiguous
        getRep().calcBiasForAccelerationConstraints(state,true,true,true,tmpbias);
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


//------------------------------------------------------------------------------
//                       MULTIPLY BY SYSTEM JACOBIAN
//------------------------------------------------------------------------------
// Ensure that we have contiguous storage and then call the underlying method.
void SimbodyMatterSubsystem::multiplyBySystemJacobian
   (const State& s, const Vector& u, Vector_<SpatialVec>& Ju) const
{   
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_ERRCHK2_ALWAYS(u.size() == nu,
        "SimbodyMatterSubsystem::multiplyBySystemJacobian()",
        "The supplied u-space Vector had length %d; expected %d.",u.size(),nu);

    const bool uIsContig  = u.hasContiguousData();
    const bool JuIsContig = Ju.hasContiguousData();

    Vector u_contig; Vector_<SpatialVec> Ju_contig; // allocate only if needed
    const Vector*        up  = uIsContig  ? &u  : (const Vector*)&u_contig;
    Vector_<SpatialVec>* Jup = JuIsContig ? &Ju : &Ju_contig;
    if (!uIsContig) {
        u_contig.resize(nu);
        u_contig(0, nu) = u; // prevent reallocation
    }

    rep.multiplyBySystemJacobian(s, *up, *Jup);

    if (!JuIsContig)
        Ju = Ju_contig;
}


//------------------------------------------------------------------------------
//                  MULTIPLY BY SYSTEM JACOBIAN TRANSPOSE
//------------------------------------------------------------------------------
void SimbodyMatterSubsystem::multiplyBySystemJacobianTranspose
   (const State& s, const Vector_<SpatialVec>& F_G, Vector& f) const
{   
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_ERRCHK2_ALWAYS(F_G.size() == nb,
        "SimbodyMatterSubsystem::multiplyBySystemJacobianTranspose()",
        "The supplied spatial forces vector had length %d; expected %d.",
        F_G.size(),nb);

    const bool F_GIsContig  = F_G.hasContiguousData();
    const bool fIsContig    = f.hasContiguousData();

    Vector_<SpatialVec> F_G_contig; Vector f_contig; // allocate only if needed
    const Vector_<SpatialVec>* F_Gp = F_GIsContig ? 
                                &F_G : (const Vector_<SpatialVec>*)&F_G_contig;
    Vector* fp   = fIsContig  ? &f  : &f_contig;
    if (!F_GIsContig) {
        F_G_contig.resize(nb);
        F_G_contig(0, nb) = F_G; // prevent reallocation
    }

    rep.multiplyBySystemJacobianTranspose(s, *F_Gp, *fp); 

    if (!fIsContig)
        f = f_contig;
}


//------------------------------------------------------------------------------
//                       CALC SYSTEM JACOBIAN (spatial)
//------------------------------------------------------------------------------
// Calculate J as an nb X n matrix of SpatialVecs, by repeated calls 
// to J*u with u=0 except one u[i]=1. Cost is 12*n*(nb+n).
// If the output matrix J_G isn't contiguous we have to allocate an nb-length 
// temporary and perform an extra copy from there into J_G.
void SimbodyMatterSubsystem::calcSystemJacobian
   (const State&            state,
    Matrix_<SpatialVec>&    J_G) const 
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    J_G.resize(nb,nu);
    Vector u(nu, Real(0));

    // If J_G is contiguous we can generate results directly into its columns.
    if (J_G.hasContiguousData()) {
        for (int j=0; j<nu; ++j) { 
            VectorView_<SpatialVec> col = J_G(j);
            u[j] = 1; rep.multiplyBySystemJacobian(state,u,col); u[j] = 0;
        }
        return;
    }

    // J_G is non-contiguous so generate results into a temporary column and
    // copy back.
    Vector_<SpatialVec> Ju(nb);
    for (int j=0; j<nu; ++j) { 
        u[j] = 1; rep.multiplyBySystemJacobian(state,u,Ju); u[j] = 0;
        VectorView_<SpatialVec> col = J_G(j);
        col = Ju;
    }
}


//------------------------------------------------------------------------------
//                       CALC SYSTEM JACOBIAN (scalars)
//------------------------------------------------------------------------------
// Alternate signature that returns a system Jacobian as a 6*nb X n Matrix 
// rather than as an nb X n matrix of spatial vectors. Note that we
// don't know whether the output matrix has contiguous rows, columns or
// neither; we'll always work with a contiguous temporary here and copy back.
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
        for (MobilizedBodyIndex mbx(0); mbx < nb; ++mbx) {
            const SpatialVec& V = Ju[mbx];
            for (int k=0; k<3; ++k) col[nxt++] = V[0][k]; // w
            for (int k=0; k<3; ++k) col[nxt++] = V[1][k]; // v
        }
    }
}


//------------------------------------------------------------------------------
//                 CALC BIAS FOR SYSTEM JACOBIAN (spatial)
//------------------------------------------------------------------------------
// This is just a synonym for getTotalCoriolisAcceleration(). No flops since
// we already computed this.
void SimbodyMatterSubsystem::calcBiasForSystemJacobian
   (const State&         state,
    Vector_<SpatialVec>& JDotu) const 
{
    // Just return the coriolis acceleration.
    const SBTreeVelocityCache& vc = getRep().getTreeVelocityCache(state);
    const Array_<SpatialVec,MobilizedBodyIndex>& tca = 
        vc.totalCoriolisAcceleration;
    const Vector_<SpatialVec> 
        AC_GB(tca.size(), (const Real*)tca.begin(), true); // shallow ref
    JDotu = AC_GB;
}


//------------------------------------------------------------------------------
//                  CALC BIAS FOR SYSTEM JACOBIAN (scalar)
//------------------------------------------------------------------------------
// Same as above but unpack into 6*nb vector rather nb spatial vecs.
void SimbodyMatterSubsystem::calcBiasForSystemJacobian
   (const State&    state,
    Vector&         JDotu) const 
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    const SBTreeVelocityCache& vc = rep.getTreeVelocityCache(state);
    const Array_<SpatialVec,MobilizedBodyIndex>& tca = 
        vc.totalCoriolisAcceleration;

    JDotu.resize(6*nb); // Might not be contiguous
    int nxt = 0; // index into JDotu
    for (MobilizedBodyIndex mbx(0); mbx < nb; ++mbx) {
        const SpatialVec& A = tca[mbx];
        for (int k=0; k<3; ++k) JDotu[nxt++] = A[0][k]; // b (angular accel)
        for (int k=0; k<3; ++k) JDotu[nxt++] = A[1][k]; // a
    }
}


//------------------------------------------------------------------------------
//                       MULTIPLY BY STATION JACOBIAN
//------------------------------------------------------------------------------
// We want v_GS = J_GS*u, the linear velocity of nt station tasks Si in Ground 
// induced by the given generalized speeds. Station Si is on body Bi and is 
// given by the Bi-frame vector p_BiS. We can easily calculate V_GB = J_GB*u,
// the spatial velocity of *all* the mobilized bodies. Then for each station
// task,
//      v_GSi = v_GBi + w_GBi X p_BiS_G
// where p_BiS_G is p_BiS re-expressed in Ground.
//
// Cost is 27*nt + 12*(nb+nu) flops.
//
// It is OK for the input u and output JSu vectors to be non-contiguous.
void SimbodyMatterSubsystem::multiplyByStationJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,    // task body
    const Array_<Vec3>&                 p_BS,       // task station
    const Vector&                       u,
    Vector_<Vec3>&                      JSu) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();

    SimTK_ERRCHK2_ALWAYS(u.size() == nu,
        "SimbodyMatterSubsystem::multiplyByStationJacobian()",
        "The supplied u-space Vector had length %d; expected %d.",u.size(),nu);

    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BS.size() == nt,
        "SimbodyMatterSubsystem::multiplyByStationJacobian()",
        "The given number of task bodies (%d) and station tasks (%d) must "
        "be the same.", nt, (int)p_BS.size());

    // First use the System Jacobian to obtain spatial velocities for *all*
    // mobilized body frames, at a cost of 12*(nb+nu) flops.
    Vector_<SpatialVec> Ju(nb); // temp Ju=J_G*u (contiguous)
    if (u.hasContiguousData())
        rep.multiplyBySystemJacobian(state,u,Ju); 
    else {
        Vector contig_u(nu); // contiguous data
        contig_u(0,nu) = u;  // no reallocation
        rep.multiplyBySystemJacobian(state,contig_u,Ju); 
    }

    // Then for each station task, determine its linear velocity at a cost of
    // 27 flops per task.
    JSu.resize(nt);
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::multiplyByStationJacobian()");

        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops
        const SpatialVec& V_GB = Ju[mobodx];
        const SpatialVec  V_GS = shiftVelocityBy(V_GB, p_BS_G);     // 12 flops
        JSu[task] = V_GS[1]; // return linear velocity only
    }
}


//------------------------------------------------------------------------------
//                   MULTIPLY BY STATION JACOBIAN TRANSPOSE
//------------------------------------------------------------------------------
// We want f = ~J_GS*f_GS, the generalized forces produced by applying
// translational task force vectors f_GSi to stations Si. Each station Si is on
// body Bi and is given by the Bi-frame vector p_BiSi. We can easily calculate
// f = ~J_GB*F_GB for task body origins, so we shift the task forces there
// with F_GBi=[p_BiSi_G X f_GSi, f_GSi]. 
// It is OK if f_GS and/or f are not contiguous.
// Cost is 30nt + 18nb + 11nu.
void SimbodyMatterSubsystem::multiplyByStationJacobianTranspose
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BS,
    const Vector_<Vec3>&                f_GS,
    Vector&                             f) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK3_ALWAYS(p_BS.size() == nt && f_GS.size() == nt,
        "SimbodyMatterSubsystem::multiplyByStationJacobianTranspose()",
        "The given number of task bodies (%d), task stations (%d), and "
        "applied task forces (%d) must all be the same.", 
        nt, (int)p_BS.size(), (int)f_GS.size());

    f.resize(nu); // might not be contiguous
    const bool fIsContig = f.hasContiguousData();
    Vector f_contig; // will get allocated only if used below
    Vector* fp = fIsContig ? &f  : &f_contig;

    // Need an array putting a spatial force on *every* body.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();

    // Collect the applied task forces into F_G.
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::multiplyByStationJacobianTranspose()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops
        F_G[mobodx] += SpatialVec(p_BS_G % f_GS[task], f_GS[task]); // 15 flops
    }

    rep.multiplyBySystemJacobianTranspose(state,F_G,*fp); // 18nb+11nu flops

    if (!fIsContig)
        f = f_contig; // copy result out
}


//------------------------------------------------------------------------------
//                       CALC STATION JACOBIAN (spatial)
//------------------------------------------------------------------------------
// Cost is 3*nt*(14 + 18nb + 11n) flops.
// Each subsequent multiply by JS_G*u would be 3*nt*(2n-1)~=6*nt*n flops.
void SimbodyMatterSubsystem::calcStationJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BS,
    Matrix_<Vec3>&                      JS_G) const // nt X nu Vec3s
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BS.size() == nt,
        "SimbodyMatterSubsystem::calcStationJacobian()",
        "The given number of task bodies (%d) and station tasks (%d) must "
        "be the same.", nt, (int)p_BS.size());

    // Calculate J=dvdu where v is linear velocity of task stations p_BS.
    // (This is nt half-rows of J.)
    JS_G.resize(nt,nu);

    // We're assuming that 3*nt << nu so that it is cheaper to calculate ~JS
    // than JS, using ~J*F rather than J*u.
    // TODO: check dimensions and use whichever method is cheaper.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    Vector col(nu); // temporary to hold column of ~J_G
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcStationJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops

        // Calculate the 3 rows of JS corresponding to this task.
        RowVectorView_<Vec3> row = JS_G[task];
        SpatialVec& Fb = F_G[mobodx]; // the only one we'll change
        for (int i=0; i < 3; ++i) {
            Fb[1][i] = 1;
            Fb[0] = p_BS_G % Fb[1]; // r X F (9 flops)
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            for (int r=0; r < nu; ++r) row[r][i] = col[r]; 
            Fb[1][i] = 0;
            Fb[0] = 0;
        }
    }
}


//------------------------------------------------------------------------------
//                       CALC STATION JACOBIAN (scalar)
//------------------------------------------------------------------------------
// Alternate signature that returns a station Jacobian as a 3*nt X nu Matrix 
// rather than as an nt X nu Matrix of Vec3s.
void SimbodyMatterSubsystem::calcStationJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BS,
    Matrix&                             JS_G) const // 3*nt X nu Vec3s
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BS.size() == nt,
        "SimbodyMatterSubsystem::calcStationJacobian()",
        "The specified number of task bodies (%d) and station tasks (%d) must "
        "be the same.", nt, (int)p_BS.size());

    // Calculate J=dvdu where v is linear velocity of p_BS.
    // (This is nt rows of J.)
    JS_G.resize(3*nt,nu);

    // We're assuming that 3*nt << nu so that it is cheaper to calculate ~JS
    // than JS, using ~J*F rather than J*u.
    // TODO: check dimensions and use whichever method is cheaper.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    Vector col(nu); // contiguous temporary to hold column of ~J_G
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcStationJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops

        // Calculate the 3 rows of JS corresponding to this task.
        SpatialVec& Fb = F_G[mobodx]; // the only one we'll change
        for (int i=0; i < 3; ++i) {
            Fb[1][i] = 1;
            Fb[0] = p_BS_G % Fb[1]; // r X F (9 flops)
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            JS_G[3*task + i] = ~col; 
            Fb[1][i] = 0;
            Fb[0] = 0;
        }
    }
}


//------------------------------------------------------------------------------
//                 CALC BIAS FOR STATION JACOBIAN (spatial)
//------------------------------------------------------------------------------
// Just get the total Coriolis acceleration for each task body and shift it to
// the task station S. Cost is 48*nt flops.
void SimbodyMatterSubsystem::calcBiasForStationJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,        // nt task bodies
    const Array_<Vec3>&                 p_BS,           // nt task stations
    Vector_<Vec3>&                      JSDotu) const   // nt of these
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BS.size() == nt,
        "SimbodyMatterSubsystem::calcBiasForStationJacobian()",
        "The specified number of task bodies (%d) and station tasks (%d) must "
        "be the same.", nt, (int)p_BS.size());

    const SBTreeVelocityCache& vc = rep.getTreeVelocityCache(state);

    JSDotu.resize(nt);
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcBiasForStationJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops
        const SpatialVec& A_GB = vc.totalCoriolisAcceleration[mobodx];
        const Vec3&       w_GB = mobod.getBodyAngularVelocity(state);
        const SpatialVec  A_GS = shiftAccelerationBy(A_GB, w_GB, p_BS_G); 
                                                                    // 33 flops 
        JSDotu[task] = A_GS[1]; // linear acceleration only
    }
}

//------------------------------------------------------------------------------
//                 CALC BIAS FOR STATION JACOBIAN (scalar)
//------------------------------------------------------------------------------
void SimbodyMatterSubsystem::calcBiasForStationJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,        // nt task bodies
    const Array_<Vec3>&                 p_BS,           // nt task stations
    Vector&                             JSDotu) const   // 3*nt  
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BS.size() == nt,
        "SimbodyMatterSubsystem::calcBiasForStationJacobian()",
        "The given number of task bodies (%d) and station tasks (%d) must "
        "be the same.", nt, (int)p_BS.size());

    const SBTreeVelocityCache& vc = rep.getTreeVelocityCache(state);

    JSDotu.resize(3*nt); // might not be contiguous
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcBiasForStationJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BS_G = 
            mobod.expressVectorInGroundFrame(state, p_BS[task]);    // 15 flops
        const SpatialVec& A_GB = vc.totalCoriolisAcceleration[mobodx];
        const Vec3&       w_GB = mobod.getBodyAngularVelocity(state);
        const SpatialVec  A_GS = shiftAccelerationBy(A_GB, w_GB, p_BS_G); 
                                                                    // 33 flops
        const Vec3&       a_GS = A_GS[1]; // linear acceleration only
        for (int k=0; k<3; ++k) JSDotu[3*task+k] = a_GS[k]; 
    }
}


//------------------------------------------------------------------------------
//                       MULTIPLY BY FRAME JACOBIAN
//------------------------------------------------------------------------------
// We want V_GA = J_GA*u, the spatial velocity of nt task frames Ai in 
// Ground induced by the given generalized speeds. Frames A are fixed on bodies
// B and would be given by the transform X_BA, except the result depends only 
// on A's origin position p_BA (==p_BoAo) because angular velocity is the same 
// for all frames fixed to the same body. We can easily calculate V_GB = J_GB*u,
// the spatial velocities at each body B's origin Bo. Then 
//      V_GAi = [w_GAi, v_GAi] = [w_GBi, v_GBi + w_GBi X p_BiAi_G]
// where p_BiAi_G is p_BiAi re-expressed in Ground.
//
// Cost is 27*nt + 12*(nb+nu) flops.
void SimbodyMatterSubsystem::multiplyByFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BA,
    const Vector&                       u,
    Vector_<SpatialVec>&                JFu) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(u.size() == nu,
        "SimbodyMatterSubsystem::multiplyByFrameJacobian()",
        "The supplied u-space Vector had length %d; expected %d.",u.size(),nu);

    SimTK_ERRCHK2_ALWAYS(p_BA.size() == nt,
        "SimbodyMatterSubsystem::multiplyByFrameJacobian()",
        "The given number of task bodies (%d) and frame tasks (%d) must "
        "be the same.", nt, (int)p_BA.size());

    // First use the System Jacobian to obtain spatial velocities for *all*
    // mobilized body frames, at a cost of 12*(nb+nu) flops.
    Vector_<SpatialVec> Ju(nb); // temp Ju=J_G*u (contiguous)
    if (u.hasContiguousData())
        rep.multiplyBySystemJacobian(state,u,Ju); 
    else {
        Vector contig_u(nu); // contiguous data
        contig_u(0,nu) = u;  // no reallocation
        rep.multiplyBySystemJacobian(state,contig_u,Ju); 
    }

    // Then for each frame task, determine its linear velocity at a cost of
    // 27 flops per task.
    JFu.resize(nt); // OK if not contiguous
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::multiplyByFrameJacobian()");

        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops
        const SpatialVec& V_GB = Ju[mobodx]; 
        JFu[task] = shiftVelocityBy(V_GB, p_BA_G);                  // 12 flops
    }
}

//------------------------------------------------------------------------------
//                    MULTIPLY BY FRAME JACOBIAN TRANSPOSE
//------------------------------------------------------------------------------
// We want f = ~J_GA*F_GA, the generalized forces produced by applying spatial
// force vectors F_GAi=[t_G,f_GAi] at Aio, the origin of task frame Ai, which 
// is fixed to some body Bi. Frame Ai would be given by transform X_BAi, but the
// result depends only on Ai's origin location p_BiAi (==p_BioAio) since a 
// torque is the same wherever it is applied. We can easily calculate
// fb = ~J_GB*F_GB so we shift each F_GAi to Bi via 
//      F_GBi = [t_G + p_BiAi_G X f_GAi, f_GAi]. 
// Cost is 33nt + 18nb + 11nu.
void SimbodyMatterSubsystem::multiplyByFrameJacobianTranspose
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BA,
    const Vector_<SpatialVec>&          F_GA,
    Vector&                             f) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(), nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK3_ALWAYS(p_BA.size() == nt && F_GA.size() == nt,
        "SimbodyMatterSubsystem::multiplyByFrameJacobianTranspose()",
        "The given number of task bodies (%d), task stations (%d), and "
        "applied task forces (%d) must all be the same.", 
        nt, (int)p_BA.size(), (int)F_GA.size());

    f.resize(nu); // might not be contiguous
    const bool fIsContig = f.hasContiguousData();
    Vector f_contig; // will get allocated only if used below
    Vector* fp = fIsContig ? &f  : &f_contig;

    // Need an array putting a spatial force on each body.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();

    // Collect the applied task forces into F_G.
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::multiplyByFrameJacobianTranspose()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops
        const SpatialVec& F_GAi = F_GA[task];
        F_G[mobodx] += SpatialVec(F_GAi[0] + p_BA_G % F_GAi[1],     // 18 flops 
                                  F_GAi[1]);
    }

    rep.multiplyBySystemJacobianTranspose(state,F_G,*fp); // 18nb+11nu flops

    if (!fIsContig)
        f = f_contig; // copy result out
}


//------------------------------------------------------------------------------
//                       CALC FRAME JACOBIAN (spatial)
//------------------------------------------------------------------------------
// Cost is 42*nt + 6*(18nb + 11nu) flops.
// Each subsequent multiply by JF_G*u would be 12*nu-6 flops.
void SimbodyMatterSubsystem::calcFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BA,
    Matrix_<SpatialVec>&                JF_G) const
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BA.size() == nt,
        "SimbodyMatterSubsystem::calcFrameJacobian()",
        "The given number of task bodies (%d) and frame tasks (%d) must "
        "be the same.", nt, (int)p_BA.size());

    // Calculate J=dVdu where V is spatial velocity of task frames A.
    // (This is nt rows of J.)
    JF_G.resize(nt,nu);

    // We're assuming that 6*nt << nu so that it is cheaper to calculate ~JF
    // than JF, using ~J*F rather than J*u.
    // TODO: check dimensions and use whichever method is cheaper.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    Vector col(nu); // temporary to hold column of ~J_G
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcFrameJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops

        // Calculate the 6 rows of JS corresponding to this task.
        RowVectorView_<SpatialVec> row = JF_G[task];
        SpatialVec& Fb = F_G[mobodx]; // the only one we'll change

        // Rotational part.
        for (int i=0; i < 3; ++i) {
            Fb[0][i] = 1;
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            for (int r=0; r < nu; ++r) row[r][0][i] = col[r]; 
            Fb[0][i] = 0;
        }

        // Translational part.
        for (int i=0; i < 3; ++i) {
            Fb[1][i] = 1;
            Fb[0] = p_BA_G % Fb[1]; // r X F (9 flops)
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            for (int r=0; r < nu; ++r) row[r][1][i] = col[r]; 
            Fb[1][i] = 0;
            Fb[0] = 0;
        }
    }
}


//------------------------------------------------------------------------------
//                        CALC FRAME JACOBIAN (scalar)
//------------------------------------------------------------------------------
// Alternate signature that returns a frame Jacobian as a 6*nt x n Matrix 
// rather than as a Matrix of SpatialVecs.
// Cost is 42*nt + 6*(18*nb + 11*n)
void SimbodyMatterSubsystem::calcFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 p_BA,
    Matrix&                             JF_G) const // 6*nt X n
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies(); // includes ground
    const int nu = rep.getNumMobilities();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BA.size() == nt,
        "SimbodyMatterSubsystem::calcFrameJacobian()",
        "The given number of task bodies (%d) and frame tasks (%d) must "
        "be the same.", nt, (int)p_BA.size());

    // Calculate J=dVdu where V is spatial velocity of task frames A.
    // (This is 6*nt rows of the scalar matrix form of J.)
    JF_G.resize(6*nt,nu);

    // We're assuming that 6*nt << nu so that it is cheaper to calculate ~JF
    // than JF, using ~J*F rather than J*u.
    // TODO: check dimensions and use whichever method is cheaper.
    Vector_<SpatialVec> F_G(nb); F_G.setToZero();
    Vector col(nu); // temporary to hold column of ~J_G
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcFrameJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops

        // Calculate the 6 rows of JS corresponding to this task.
        SpatialVec& Fb = F_G[mobodx]; // the only one we'll change

        // Rotational part.
        for (int i=0; i < 3; ++i) {
            Fb[0][i] = 1;
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            JF_G[6*task + i] = ~col;
            Fb[0][i] = 0;
        }

        // Translational part.
        for (int i=0; i < 3; ++i) {
            Fb[1][i] = 1;
            Fb[0] = p_BA_G % Fb[1]; // r X F (9 flops)
            rep.multiplyBySystemJacobianTranspose(state,F_G,col);// 18nb+11nu flops
            JF_G[6*task + 3 + i] = ~col;
            Fb[1][i] = 0;
            Fb[0] = 0;
        }
    }
}


//------------------------------------------------------------------------------
//                CALC BIAS FOR FRAME JACOBIAN (spatial)
//------------------------------------------------------------------------------
// Just get the total Coriolis acceleration for this body and shift it to
// the A frame origin Ao.
// Cost is 48*nt.
void SimbodyMatterSubsystem::calcBiasForFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,        // nt task bodies
    const Array_<Vec3>&                 p_BA,           // nt task frames
    Vector_<SpatialVec>&                JFDotu) const   // nt of these
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BA.size() == nt,
        "SimbodyMatterSubsystem::calcBiasForFrameJacobian()",
        "The given number of task bodies (%d) and frame tasks (%d) must "
        "be the same.", nt, (int)p_BA.size());

    const SBTreeVelocityCache& vc = rep.getTreeVelocityCache(state);

    JFDotu.resize(nt);
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcBiasForFrameJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops
        const SpatialVec& A_GB = vc.totalCoriolisAcceleration[mobodx];
        const Vec3&       w_GB = mobod.getBodyAngularVelocity(state);
        const SpatialVec  A_GA = shiftAccelerationBy(A_GB, w_GB, p_BA_G); 
                                                                    // 33 flops 
        JFDotu[task] = A_GA; // linear acceleration only
    }
}


//------------------------------------------------------------------------------
//                CALC BIAS FOR FRAME JACOBIAN (scalar)
//------------------------------------------------------------------------------
// Just get the total Coriolis acceleration for this body and shift it to
// the A frame origin Ao.
// Cost is 48*nt.
void SimbodyMatterSubsystem::calcBiasForFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,        // nt task bodies
    const Array_<Vec3>&                 p_BA,           // nt task frames
    Vector&                             JFDotu) const   // 6*nt
{
    const SimbodyMatterSubsystemRep& rep = getRep();
    const int nb = rep.getNumBodies();
    const int nt = (int)onBodyB.size(); // number of tasks

    SimTK_ERRCHK2_ALWAYS(p_BA.size() == nt,
        "SimbodyMatterSubsystem::calcBiasForFrameJacobian()",
        "The given number of task bodies (%d) and frame tasks (%d) must "
        "be the same.", nt, (int)p_BA.size());

    const SBTreeVelocityCache& vc = rep.getTreeVelocityCache(state);

    JFDotu.resize(6*nt); // might not be contiguous
    for (int task=0; task < nt; ++task) {
        const MobilizedBodyIndex mobodx = onBodyB[task];
        SimTK_INDEXCHECK(mobodx, nb,
            "SimbodyMatterSubsystem::calcBiasForFrameJacobian()");
        const MobilizedBody& mobod = rep.getMobilizedBody(mobodx);
        const Vec3 p_BA_G = 
            mobod.expressVectorInGroundFrame(state, p_BA[task]);    // 15 flops
        const SpatialVec& A_GB = vc.totalCoriolisAcceleration[mobodx];
        const Vec3&       w_GB = mobod.getBodyAngularVelocity(state);
        const SpatialVec  A_GA = shiftAccelerationBy(A_GB, w_GB, p_BA_G); 
                                                                    // 33 flops 
        for (int k=0; k<3; ++k) JFDotu[6*task+k]   = A_GA[0][k]; 
        for (int k=0; k<3; ++k) JFDotu[6*task+3+k] = A_GA[1][k]; 
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

const Vector& SimbodyMatterSubsystem::
getMotionMultipliers(const State& s) const 
{   return getRep().getMotionMultipliers(s); }

Vector SimbodyMatterSubsystem::
calcMotionErrors(const State& s, const Stage& stage) const
{   return getRep().calcMotionErrors(s,stage); }

void SimbodyMatterSubsystem::
findMotionForces(const State&         s,
                 Vector&              mobilityForces) const 
{   getRep().findMotionForces(s, mobilityForces); }

const Vector& SimbodyMatterSubsystem::
getConstraintMultipliers(const State& s) const
{   return getRep().getConstraintMultipliers(s); }

void SimbodyMatterSubsystem::
findConstraintForces(const State&         s, 
                     Vector_<SpatialVec>& bodyForcesInG,
                     Vector&              mobilityForces) const 
{   getRep().findConstraintForces(s, bodyForcesInG, mobilityForces); }

Real SimbodyMatterSubsystem::
calcMotionPower(const State& s) const
{   return getRep().calcMotionPower(s); }

Real SimbodyMatterSubsystem::
calcConstraintPower(const State& s) const
{   return getRep().calcConstraintPower(s); }

void SimbodyMatterSubsystem::
calcConstraintForcesFromMultipliers(const State&         s,
                                    const Vector&        lambda,
                                    Vector_<SpatialVec>& bodyForcesInG, 
                                    Vector&              mobilityForces) const
{   getRep().calcConstraintForcesFromMultipliers
                (s,lambda,bodyForcesInG,mobilityForces); }

void SimbodyMatterSubsystem::
calcMobilizerReactionForcesUsingFreebodyMethod
   (const State& s, Vector_<SpatialVec>& forces) const 
{   getRep().calcMobilizerReactionForcesUsingFreebodyMethod(s, forces); }


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
void SimbodyMatterSubsystem::setUseEulerAnglesByDefault(bool useAngles)
  { updRep().setUseEulerAnglesByDefault(useAngles); }
void SimbodyMatterSubsystem::setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disabled) const
  { getRep().setConstraintIsDisabled(s,constraint,disabled); }
bool SimbodyMatterSubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodyMatterSubsystem::getUseEulerAnglesByDefault() const
  { return getRep().getUseEulerAnglesByDefault(); }
bool SimbodyMatterSubsystem::isConstraintDisabled(const State& s, ConstraintIndex constraint) const
  { return getRep().isConstraintDisabled(s,constraint); }
void SimbodyMatterSubsystem::convertToEulerAngles(const State& inputState, State& outputState) const
  { return getRep().convertToEulerAngles(inputState, outputState); }
void SimbodyMatterSubsystem::convertToQuaternions(const State& inputState, State& outputState) const
  { return getRep().convertToQuaternions(inputState, outputState); }

void SimbodyMatterSubsystem::normalizeQuaternions(State& state) const {
    Vector dummy; // no error estimate to correct
    getRep().normalizeQuaternions(state, dummy);
}

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
SimbodyMatterSubsystem::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getTotalCentrifugalForces(s,body);
}

//TODO: user access to this quantity is deprecated in Simbody 3.6.
const SpatialVec&
SimbodyMatterSubsystem::getMobilizerCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
    return getRep().getArticulatedBodyCentrifugalForces(s,body);
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

void SimbodyMatterSubsystem::
realizePositionKinematics(const State& s) const {
    getRep().realizePositionKinematics(s);
}


void SimbodyMatterSubsystem::
realizeVelocityKinematics(const State& s) const {
    getRep().realizeVelocityKinematics(s);
}

void SimbodyMatterSubsystem::
realizeCompositeBodyInertias(const State& s) const {
    getRep().realizeCompositeBodyInertias(s);
}

void SimbodyMatterSubsystem::
realizeArticulatedBodyInertias(const State& s) const {
    getRep().realizeArticulatedBodyInertias(s);
}

void SimbodyMatterSubsystem::
realizeArticulatedBodyVelocity(const State& s) const {
    getRep().realizeArticulatedBodyVelocity(s);
}

void SimbodyMatterSubsystem::
invalidatePositionKinematics(const State& s) const {
    getRep().invalidatePositionKinematics(s);
}

void SimbodyMatterSubsystem::
invalidateVelocityKinematics(const State& s) const {
    getRep().invalidateVelocityKinematics(s);
}

void SimbodyMatterSubsystem::
invalidateCompositeBodyInertias(const State& s) const {
    getRep().invalidateCompositeBodyInertias(s);
}

void SimbodyMatterSubsystem::
invalidateArticulatedBodyInertias(const State& s) const {
    getRep().invalidateArticulatedBodyInertias(s);
}

void SimbodyMatterSubsystem::
invalidateArticulatedBodyVelocity(const State& s) const {
    getRep().invalidateArticulatedBodyVelocity(s);
}

bool SimbodyMatterSubsystem::
isPositionKinematicsRealized(const State& state) const
{   return getRep().isPositionKinematicsRealized(state); }
bool SimbodyMatterSubsystem::
isVelocityKinematicsRealized(const State& state) const
{   return getRep().isVelocityKinematicsRealized(state); }
bool SimbodyMatterSubsystem::
isCompositeBodyInertiasRealized(const State& state) const
{   return getRep().isCompositeBodyInertiasRealized(state); }
bool SimbodyMatterSubsystem::
isArticulatedBodyInertiasRealized(const State& state) const
{   return getRep().isArticulatedBodyInertiasRealized(state); }
bool SimbodyMatterSubsystem::
isArticulatedBodyVelocityRealized(const State& state) const
{   return getRep().isArticulatedBodyVelocityRealized(state); }

const Array_<QIndex>& SimbodyMatterSubsystem::
getFreeQIndex(const State& state) const
{   return getRep().getFreeQIndex(state); }
const Array_<UIndex>& SimbodyMatterSubsystem::
getFreeUIndex(const State& state) const
{   return getRep().getFreeUIndex(state); }
const Array_<UIndex>& SimbodyMatterSubsystem::
getFreeUDotIndex(const State& state) const
{   return getRep().getFreeUDotIndex(state); }
const Array_<UIndex>& SimbodyMatterSubsystem::
getKnownUDotIndex(const State& state) const
{   return getRep().getKnownUDotIndex(state); }
void SimbodyMatterSubsystem::
packFreeQ(const State& s, const Vector& allQ, Vector& packedFreeQ) const
{   getRep().packFreeQ(s,allQ,packedFreeQ); }
void SimbodyMatterSubsystem::
unpackFreeQ(const State& s, const Vector& packedFreeQ, Vector& unpackedFreeQ) const
{   getRep().unpackFreeQ(s,packedFreeQ,unpackedFreeQ); }
void SimbodyMatterSubsystem::
packFreeU(const State& s, const Vector& allU, Vector& packedFreeU) const
{   getRep().packFreeU(s,allU,packedFreeU); }
void SimbodyMatterSubsystem::
unpackFreeU(const State& s, const Vector& packedFreeU, Vector& unpackedFreeU) const
{   getRep().unpackFreeU(s,packedFreeU,unpackedFreeU); }

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

