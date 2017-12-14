#ifndef SimTK_SIMBODY_CONSTRAINT_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_IMPL_H_

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
Private implementation of Constraint and its included subclasses which
represent the built-in constraint types. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/SimbodyMatterSubtree.h"

#include "SimbodyTreeState.h"

#include <map>
#include <utility>  // std::pair
#include <iostream>
using std::cout; using std::endl;

class SimbodyMatterSubsystemRep;

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;

//==============================================================================
//                           CONSTRAINT IMPL
//==============================================================================
// This is what a Constraint handle points to.
class ConstraintImpl : public PIMPLImplementation<Constraint, ConstraintImpl> {
public:

ConstraintImpl()
    : myMatterSubsystemRep(0), 
    defaultMp(0), defaultMv(0), defaultMa(0), defaultDisabled(false),
    constraintIsConditional(false), myAncestorBodyIsNotGround(false)
{
}
virtual ~ConstraintImpl() { }
virtual ConstraintImpl* clone() const = 0;

ConstraintImpl(int mp, int mv, int ma)
    : myMatterSubsystemRep(0), 
    defaultMp(mp), defaultMv(mv), defaultMa(ma), defaultDisabled(false),
    constraintIsConditional(false), myAncestorBodyIsNotGround(false)
{
}

void setDefaultNumConstraintEquations(int mp, int mv, int ma) {
    assert(mp >= 0 && mv >= 0 && ma >= 0);
    invalidateTopologyCache();
    defaultMp = mp;
    defaultMv = mv;
    defaultMa = ma;
}

void getDefaultNumConstraintEquations(int& mp, int& mv, int& ma) const {
    mp = defaultMp;
    mv = defaultMv;
    ma = defaultMa;
}

void setDisabledByDefault(bool shouldBeDisabled) {
    invalidateTopologyCache();
    defaultDisabled = shouldBeDisabled;
}

bool isDisabledByDefault() const {
    return defaultDisabled;
}

void setDisabled(State& s, bool shouldBeDisabled) const ;
bool isDisabled(const State& s) const;

void setIsConditional(bool isConditional) {
    invalidateTopologyCache();
    constraintIsConditional = isConditional;
}

bool isConditional() const {return constraintIsConditional;}

typedef std::map<MobilizedBodyIndex,ConstrainedBodyIndex>       
    MobilizedBody2ConstrainedBodyMap;
typedef std::map<MobilizedBodyIndex,ConstrainedMobilizerIndex>  
    MobilizedBody2ConstrainedMobilizerMap;

// Call this during construction phase to add a body to the topological 
// structure of this Constraint. This body's mobilizer's mobilities are *not* 
// part of the constraint; mobilizers must be added separately. If this
// mobilized body has been seen as a constrained body before we'll return the
// same index as first assigned to it.
ConstrainedBodyIndex addConstrainedBody(const MobilizedBody&);

// Call this during construction phase to add a mobilizer to the topological 
// structure of this Constraint. All the coordinates q and mobilities u for 
// this mobilizer are added as "constrainable". We don't know how many
// coordinates and speeds there are until Stage::Model. If this
// mobilized body has been seen as a constrained mobilizer before we'll return
// the same index as first assigned to it.
ConstrainedMobilizerIndex addConstrainedMobilizer(const MobilizedBody&);

MobilizedBodyIndex getMobilizedBodyIndexOfConstrainedBody
   (ConstrainedBodyIndex c) const {
    assert(0 <= c && c < (int)myConstrainedBodies.size());
    return myConstrainedBodies[c];
}
MobilizedBodyIndex getMobilizedBodyIndexOfConstrainedMobilizer
   (ConstrainedMobilizerIndex c) const {
    assert(0 <= c && c < (int)myConstrainedMobilizers.size());
    return myConstrainedMobilizers[c];
}

//TODO: Constraint-local State allocation
//int allocateDiscreteVariable(State& s, Stage g, AbstractValue* v) const;
//int allocateCacheEntry(State& s, Stage g, AbstractValue* v) const;

QIndex getQIndexOfConstrainedQ(const State& s, ConstrainedQIndex cqx) const;
UIndex getUIndexOfConstrainedU(const State& s, ConstrainedUIndex cqx) const;

void convertQForcesToUForces(const State&                          s, 
                             const Array_<Real,ConstrainedQIndex>& qForces,
                             Array_<Real,ConstrainedUIndex>&       uForces) const;

// Given a position-stage state and an array of body velocities for all
// bodies (relative to Ground), select the short list of constrained bodies
// for this Constraint and ensure that their velocities are relative to
// the Ancestor body.
void convertBodyVelocityToConstrainedBodyVelocity
   (const State&                                    state,
    const Array_<SpatialVec, MobilizedBodyIndex>&   V_GB,
    Array_<SpatialVec, ConstrainedBodyIndex>&       V_AB) const;

// Given a velocity-stage state and an array of body accelerations for all
// bodies (relative to Ground), select the short list of constrained bodies
// for this Constraint and ensure that their accelerations are relative to
// the Ancestor body.
void convertBodyAccelToConstrainedBodyAccel
   (const State&                                    state,
    const Array_<SpatialVec, MobilizedBodyIndex>&   A_GB,
    Array_<SpatialVec, ConstrainedBodyIndex>&       A_AB) const;

void realizeTopology(State&)       const; // eventually calls realizeTopologyVirtual()
void realizeModel   (State&)       const; // eventually calls realizeModelVirtual() 
void realizeInstance(const State&) const; // eventually calls realizeInstanceVirtual() 

// These are called in loops over all the Constraints from the 
// SimbodyMatterSubsystem realize() methods, which will have already 
// deconstructed the State into an SBStateDigest object.
void realizeTime    (const SBStateDigest&) const; // eventually calls realizeTimeVirtual() 
void realizePosition(const SBStateDigest&) const; // eventually calls realizePositionVirtual() 
void realizeVelocity(const SBStateDigest&) const; // eventually calls realizeVelocityVirtual() 
void realizeDynamics(const SBStateDigest&) const; // eventually calls realizeDynamicsVirtual() 
void realizeAcceleration
                    (const SBStateDigest&) const; // eventually calls realizeAccelerationVirtual() 
void realizeReport  (const State&) const; // eventually calls realizeReportVirtual() 

// Given a state realized to Position stage, extract the position constraint 
// errors corresponding to this Constraint. The 'mp' argument is for sanity 
// checking -- it is an error if that isn't an exact match for the current 
// number of holonomic constraint equations generated by this Constraint. We 
// expect that perr points to an array of at least mp elements that we can 
// write on.
void getPositionErrors(const State& s, int mp, Real* perr) const;

// Given a State realized to Velocity stage, extract the velocity constraint 
// errors corresponding to this Constraint. This includes velocity constraints 
// which were produced by differentiation of holonomic (position) constraints, 
// and nonholonomic constraints which are introduced at the velocity level. The
// 'mpv' argument is for sanity checking -- it is an error if that isn't an 
// exact match for the current number of holonomic+nonholonomic (mp+mv) 
// constraint equations generated by this Constraint. We expect that pverr 
// points to an array of at least mp+mv elements that we can write on.
void getVelocityErrors(const State& s, int mpv, Real* pverr) const;

// Given a State realized to Acceleration stage, extract the accleration 
// constraint errors corresponding to this Constraint. This includes 
// acceleration constraints which were produced by twice differentiation of 
// holonomic (position) constraints, and differentiation of nonholonomic 
// (velocity) constraints, and acceleration-only constraints which are first 
// introduced at the acceleration level. The 'mpva' argument is for sanity 
// checking -- it is an error if that isn't an exact match for the current 
// number of holonomic+nonholonomic+accelerationOnly (mp+mv+ma) constraint
// equations generated by this Constraint. We expect that pvaerr points to an 
// array of at least mp+mv+ma elements that we can write on.
void getAccelerationErrors(const State& s, int mpva, Real* pvaerr) const;

// Given a State realized to Acceleration stage, extract the constraint 
// multipliers lambda corresponding to this constraint. This includes 
// multipliers for all the holonomic, nonholonomic, and acceleration-only 
// constraints (but not quaternion constraints which do not use multipliers). 
// The 'mpva' argument is for sanity checking -- it is an error if that isn't 
// an exact match for the current number (mp+mv+ma) of constraint equations 
// generated by this Constraint. We expect that lambda points to an array of at
// least mp+mv+ma elements that we can write on.
void getMultipliers(const State& s, int mpva, Real* lambda) const;

// Return a small, writable array directly referencing the segment of the longer 
// passed-in array that belongs to this constraint. State must be realized 
// through Instance stage.
ArrayView_<SpatialVec,ConstrainedBodyIndex>
updConstrainedBodyForces(const State&        state,
                         Array_<SpatialVec>& allConsBodyForces) const;
// Same, but for mobility forces.
ArrayView_<Real,ConstrainedUIndex>
updConstrainedMobilityForces(const State&  state,
                             Array_<Real>& allConsMobForces) const;

// Return a const reference to our segment. Can use above methods efficiently 
// since ArrayView is derived from ArrayViewConst; this just does a 
// shallow copy to fill in the ArrayViewConst handle; no heap is involved.
ArrayViewConst_<SpatialVec,ConstrainedBodyIndex>
getConstrainedBodyForces(const State&              state,
                         const Array_<SpatialVec>& allConsBodyForces) const 
{
    return updConstrainedBodyForces
       (state, const_cast<Array_<SpatialVec>&>(allConsBodyForces));
}
ArrayViewConst_<Real,ConstrainedUIndex>
getConstrainedMobilityForces(const State&        state,
                             const Array_<Real>& allConsMobForces) const 
{
    return updConstrainedMobilityForces
       (state, const_cast<Array_<Real>&>(allConsMobForces));
}

// Same as above but we use the matter subsystem's constrained acceleration
// cache as the source of the full-sized constrained body forces array,
// where the forces are expressed in Ground. That cache entry must have been 
// realized which occurs during Acceleration stage computations.
ArrayViewConst_<SpatialVec,ConstrainedBodyIndex>
getConstrainedBodyForcesInGFromState(const State& state) const
{
    const SBConstrainedAccelerationCache& cac = 
        getConstrainedAccelerationCache(state);
    return getConstrainedBodyForces(state, cac.constrainedBodyForcesInG);
}

ArrayView_<SpatialVec,ConstrainedBodyIndex>
updConstrainedBodyForcesInGFromState(const State& state) const
{
    SBConstrainedAccelerationCache& cac = 
        updConstrainedAccelerationCache(state);
    return updConstrainedBodyForces(state, cac.constrainedBodyForcesInG);
}

ArrayViewConst_<Real,ConstrainedUIndex>
getConstrainedMobilityForcesFromState(const State& state) const
{
    const SBConstrainedAccelerationCache& cac = 
        getConstrainedAccelerationCache(state);
    return getConstrainedMobilityForces(state, cac.constraintMobilityForces);
}

ArrayView_<Real,ConstrainedUIndex>
updConstrainedMobilityForcesFromState(const State& state) const
{
    SBConstrainedAccelerationCache& cac = 
        updConstrainedAccelerationCache(state);
    return getConstrainedMobilityForces(state, cac.constraintMobilityForces);
}

// Given a State and a set of m multipliers lambda, 
// calculate in O(m) time the constraint forces (body forces and torques and 
// u-space mobility forces) which would be generated by those multipliers. You can 
// restrict this to P,V,A subsets setting mp, mv, or ma to zero; in any case
// m = mp+mv+ma and any non-zero segments must match the corresponding number
// of constraint equations of that type.
//
// The State must be realized to at least Position stage, and if you ask for
// forces from holonomic (mv>0) or acceleration-only (ma>0) constraints then
// the State will have to be realized to Velocity stage if the corresponding
// constraint equations are velocity dependent.
void calcConstraintForcesFromMultipliers
   (const State& s,
    const Array_<Real>&                      lambdap, // 0 or mp of these
    const Array_<Real>&                      lambdav, // 0 or mv of these
    const Array_<Real>&                      lambdaa, // 0 or ma of these
    Array_<SpatialVec,ConstrainedBodyIndex>& bodyForcesInA, 
    Array_<Real,      ConstrainedUIndex>&    mobilityForces) const
{
    int actual_mp,actual_mv,actual_ma;
    getNumConstraintEquationsInUse(s, actual_mp, actual_mv, actual_ma);

    bodyForcesInA.resize(getNumConstrainedBodies()); 
    bodyForcesInA.fill(SpatialVec(Vec3(0), Vec3(0)));
    mobilityForces.resize(getNumConstrainedU(s));    
    mobilityForces.fill(Real(0));

    if (lambdap.size()) {
        assert(lambdap.size() == actual_mp);
        const int ncq = getNumConstrainedQ(s);
        const int ncu = getNumConstrainedU(s);
        // State need only be realized to Position stage
        Array_<Real, ConstrainedQIndex> qForces(ncq, Real(0));
        Array_<Real, ConstrainedUIndex> uForces(ncu);
        addInPositionConstraintForces(s, lambdap, 
                                      bodyForcesInA, qForces);
        convertQForcesToUForces(s, qForces, uForces);
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
            mobilityForces[cux] += uForces[cux];
    }
    if (lambdav.size()) {
        assert(lambdav.size() == actual_mv);
        // State may need to be realized to Velocity stage
        addInVelocityConstraintForces(s, lambdav, 
                                      bodyForcesInA, mobilityForces);
    }
    if (lambdaa.size()) {
        assert(lambdaa.size() == actual_ma);
        // State may need to be realized to Velocity stage
        addInAccelerationConstraintForces(s, lambdaa,
                                          bodyForcesInA, mobilityForces);
    }
}

// Given a State realized to Position stage, and a set of spatial forces applied
// to the constrained bodies and u-space generalized forces applied to the
// constrained mobilities, convert these to an equivalent set 
// of n generalized forces applied at each of the participating mobilities, in 
// O(n) time.
// TODO
void convertConstraintForcesToGeneralizedForces(const State& s,
    const Array_<SpatialVec,ConstrainedBodyIndex>& bodyForcesInA, 
    const Array_<Real,      ConstrainedUIndex>&    mobilityForces,
    Vector& generalizedForces) const
{
    // TODO
    assert(!"convertConstraintForcesToGeneralizedForces: not implemented yet");
}

// Calculate f = ~G*lambda in O(n+m) time. ~G=[~P ~V ~A] and you can work with 
// any subblock or combination by setting some of mp,mv,ma to zero. If nonzero
// they have to match the actual number of holonomic, nonholonomic, and 
// acceleration-only constraints. Vector lambda (typically Lagrange multipliers
// but not necessarily) is segmented lambda=[mp|mv|ma] where some of the 
// segments can be empty.
void calcGTransposeLambda
   (const State& s,
    const Array_<Real>&                      lambdap, // 0 or mp of these
    const Array_<Real>&                      lambdav, // 0 or mv of these
    const Array_<Real>&                      lambdaa, // 0 or ma of these
    Vector&                                  f) const
{
    Array_<SpatialVec,ConstrainedBodyIndex> bodyForcesInA;
    Array_<Real,      ConstrainedUIndex>    mobilityForces;

    calcConstraintForcesFromMultipliers
       (s, lambdap, lambdav, lambdaa, bodyForcesInA, mobilityForces);
    convertConstraintForcesToGeneralizedForces
       (s, bodyForcesInA, mobilityForces, f);
}

// Find the indicated cache in the passed-in State. The "get" methods require
// that the cache entry has already been marked valid.

const SBModelCache&             getModelCache(const State&) const;
const SBInstanceCache&          getInstanceCache(const State&) const;
const SBTreePositionCache&      getTreePositionCache(const State&) const;
const SBTreeVelocityCache&      getTreeVelocityCache(const State&) const;
const SBTreeAccelerationCache&  getTreeAccelerationCache(const State&) const;
const SBConstrainedAccelerationCache& 
    getConstrainedAccelerationCache(const State& s) const;
SBConstrainedAccelerationCache& 
    updConstrainedAccelerationCache(const State& s) const;

    // Methods for use with ConstrainedMobilizers.

Real getOneQFromState
   (const State&, ConstrainedMobilizerIndex, MobilizerQIndex) const;
Real getOneQDotFromState   
   (const State&, ConstrainedMobilizerIndex, MobilizerQIndex) const;
Real getOneQDotDotFromState
   (const State&, ConstrainedMobilizerIndex, MobilizerQIndex) const;

Real getOneUFromState
   (const State&, ConstrainedMobilizerIndex, MobilizerUIndex) const;
Real getOneUDotFromState
   (const State&, ConstrainedMobilizerIndex, MobilizerUIndex) const;

// Analogous methods for use when the generalized coordinate has been provided
// as an operand. The state is still necessary for modeling information.
Real getOneQ(const State& s,
             const Array_<Real,ConstrainedQIndex>&  cq,
             ConstrainedMobilizerIndex              M, 
             MobilizerQIndex                        whichQ) const
{
    assert(cq.size() == getNumConstrainedQ(s));
    assert(0 <= whichQ && whichQ < getNumConstrainedQ(s, M));
    return cq[getConstrainedQIndex(s,M,whichQ)];
}
Real getOneQDot(const State& s,
                const Array_<Real,ConstrainedQIndex>&  cqdot,
                ConstrainedMobilizerIndex              M, 
                MobilizerQIndex                        whichQ) const
{
    assert(cqdot.size() == getNumConstrainedQ(s));
    assert(0 <= whichQ && whichQ < getNumConstrainedQ(s, M));
    return cqdot[getConstrainedQIndex(s,M,whichQ)];
}
Real getOneQDotDot(const State& s,
                   const Array_<Real,ConstrainedQIndex>&  cqdotdot,
                   ConstrainedMobilizerIndex              M, 
                   MobilizerQIndex                        whichQ) const
{
    assert(cqdotdot.size() == getNumConstrainedQ(s));
    assert(0 <= whichQ && whichQ < getNumConstrainedQ(s, M));
    return cqdotdot[getConstrainedQIndex(s,M,whichQ)];
}

Real getOneU(const State& s,
             const Array_<Real,ConstrainedUIndex>&  cu,
             ConstrainedMobilizerIndex              M, 
             MobilizerUIndex                        whichU) const
{
    assert(cu.size() == getNumConstrainedU(s));
    assert(0 <= whichU && whichU < getNumConstrainedU(s, M));
    return cu[getConstrainedUIndex(s,M,whichU)];
}

Real getOneUDot(const State& s,
                const Array_<Real,ConstrainedUIndex>&  cudot,
                ConstrainedMobilizerIndex              M, 
                MobilizerUIndex                        whichU) const
{
    assert(cudot.size() == getNumConstrainedU(s));
    assert(0 <= whichU && whichU < getNumConstrainedU(s, M));
    return cudot[getConstrainedUIndex(s,M,whichU)];
}

// Apply a u-space (mobility) generalized force fu to a particular mobility of 
// the given constrained mobilizer, adding it in to the appropriate slot of the 
// mobilityForces vector which is of length numConstrainedU for this Constraint.
void addInOneMobilityForce
   (const State&                        s, 
    ConstrainedMobilizerIndex           cmx, 
    MobilizerUIndex                     whichU,
    Real                                fu, 
    Array_<Real,ConstrainedUIndex>&     mobilityForces) const 
{ 
    assert(mobilityForces.size() == getNumConstrainedU(s));
    assert(0 <= whichU && whichU < getNumConstrainedU(s, cmx));
    mobilityForces[getConstrainedUIndex(s,cmx,whichU)] += fu;
}

// Apply a q-space generalized force fq to a particular generalized coordinate q
// of the given constrained mobilizer, adding it in to the appropriate slot of 
// the qForces vector which is of length numConstrainedQ for this Constraint.
void addInOneQForce
   (const State&                        s, 
    ConstrainedMobilizerIndex           cmx, 
    MobilizerQIndex                     whichQ,
    Real                                fq, 
    Array_<Real,ConstrainedQIndex>&     qForces) const 
{ 
    assert(qForces.size() == getNumConstrainedQ(s));
    assert(0 <= whichQ && whichQ < getNumConstrainedQ(s, cmx));
    qForces[getConstrainedQIndex(s,cmx,whichQ)] += fq;
}

    // Methods for use with ConstrainedBodies.

// These are used to retrieve the indicated values from the State cache, with all values
// measured and expressed in the Ancestor (A) frame.
const Transform&  getBodyTransformFromState   
   (const State& s, ConstrainedBodyIndex B) const;      // X_AB
const SpatialVec& getBodyVelocityFromState   
   (const State& s, ConstrainedBodyIndex B) const;      // V_AB

// Extract just the rotational quantities from the spatial quantities above.
const Rotation& getBodyRotationFromState
   (const State& s, ConstrainedBodyIndex B)     const   // R_AB
    {return getBodyTransformFromState(s,B).R();}
const Vec3& getBodyAngularVelocityFromState
   (const State& s, ConstrainedBodyIndex B)     const   // w_AB
    {return getBodyVelocityFromState(s,B)[0];}

// Extract just the translational (linear) quantities from the spatial quantities above.
const Vec3& getBodyOriginLocationFromState
   (const State& s, ConstrainedBodyIndex B)     const   // p_AB
    {return getBodyTransformFromState(s,B).p();}  
const Vec3& getBodyOriginVelocityFromState
   (const State& s, ConstrainedBodyIndex B)     const   // v_AB
    {return getBodyVelocityFromState(s,B)[1];}     

// These are analogous methods for when you've been given X_AB, V_AB, or A_AB
// explicitly as an operand. These are all inline and trivial.
const Transform& getBodyTransform
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            B) const
{   return allX_AB[B]; }

const SpatialVec& getBodyVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            B) const
{   return allV_AB[B]; }

const SpatialVec& getBodyAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            B) const
{   return allA_AB[B]; }

const Rotation& getBodyRotation
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            B) const
{   return allX_AB[B].R(); }

const Vec3& getBodyAngularVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            B) const
{   return allV_AB[B][0]; }

const Vec3& getBodyAngularAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            B) const
{   return allA_AB[B][0]; }

const Vec3& getBodyOriginLocation
   (const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    ConstrainedBodyIndex                            B) const
{   return allX_AB[B].p(); }

const Vec3& getBodyOriginVelocity
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    ConstrainedBodyIndex                            B) const
{   return allV_AB[B][1]; }

const Vec3& getBodyOriginAcceleration
   (const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            B) const
{   return allA_AB[B][1]; }


// Given a station S on body B, specified by the vector p_BS from Bo to S
// expressed in B, find its location p_AS measured from Ao and expressed in the
// Ancestor's frame A. 18 flops.
Vec3 findStationLocationFromState
   (const State& s, ConstrainedBodyIndex B, const Vec3& p_BS) const {
    return getBodyTransformFromState(s,B) * p_BS; // re-measure and re-express
}

// Same, but we're given the constrained body poses as an operand (18 flops).
Vec3 findStationLocation
   (const Array_<Transform, ConstrainedBodyIndex>& allX_AB, 
    ConstrainedBodyIndex B, const Vec3& p_BS) const {
    const Transform& X_AB = allX_AB[B];
    return X_AB * p_BS; // re-measure and re-express (18 flops)
}

// Given a station S on body B, find its velocity measured and expressed in the
// Ancestor's frame A. 27 flops.
Vec3 findStationVelocityFromState
   (const State& s, ConstrainedBodyIndex B, const Vec3& p_BS) const {
    // p_BS rexpressed in A but not shifted to Ao
    const Vec3        p_BS_A = getBodyRotationFromState(s,B) * p_BS; // 15 flops
    const SpatialVec& V_AB   = getBodyVelocityFromState(s,B);
    return V_AB[1] + (V_AB[0] % p_BS_A);                             // 12 flops
}

// Same, but only configuration comes from state; velocities are an operand.
Vec3 findStationVelocity
   (const State& s,
    const Array_<SpatialVec, ConstrainedBodyIndex>& allV_AB, 
    ConstrainedBodyIndex B, const Vec3& p_BS) const 
{
    // p_BS rexpressed in A but not shifted to Ao
    const Vec3        p_BS_A  = getBodyRotationFromState(s,B) * p_BS; 
    const SpatialVec& V_AB = allV_AB[B];
    return V_AB[1] + (V_AB[0] % p_BS_A);
}

// Combo method is cheaper. Location comes from state, velocities from operand.
// NOTE: you must provide the p_BS vector expressed (but not measured) in A.
// 15 flops.
void findStationInALocationVelocity
   (const Transform&    X_AB,
    const SpatialVec&   V_AB,
    const Vec3&         p_BS_A,
    Vec3& p_AS, Vec3& v_AS) const 
{
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];
    const Vec3 pdot_BS_A = w_AB % p_BS_A;   //  9 flops

    p_AS = X_AB.p() + p_BS_A;               //  3 flops
    v_AS = v_AB + pdot_BS_A;                //  3 flops
}

// Given a station P on body F, and a station Q on body B, report the
// relative position p_PQ_A and relative velocity v_PQ_A=v_FQ_A. Note that
// the velocity is measured in the F frame, but expressed in the common A
// frame. 78 flops
void findRelativePositionVelocity
   (const Transform&    X_AF,
    const SpatialVec&   V_AF,
    const Vec3&         p_FP,
    const Transform&    X_AB,
    const SpatialVec&   V_AB,
    const Vec3&         p_BQ,
    Vec3& p_PQ_A, Vec3& v_FQ_A) const 
{
    const Vec3& w_AF = V_AF[0];

    // Express the point vectors in the A frame, but still measuring from the
    // body origins.
    const Vec3 p_FP_A = X_AF.R() * p_FP; // 15 flops
    const Vec3 p_BQ_A = X_AB.R() * p_BQ; // 15 flops
    
    Vec3 p_AP, v_AP;
    findStationInALocationVelocity(X_AF, V_AF, p_FP_A, p_AP, v_AP);//15 flops

    Vec3 p_AQ, v_AQ;
    findStationInALocationVelocity(X_AB, V_AB, p_BQ_A, p_AQ, v_AQ);//15 flops

    p_PQ_A = p_AQ - p_AP;                // 3 flops
    const Vec3 p_PQ_A_dot = v_AQ - v_AP; // derivative in A (3 flops)
    v_FQ_A = p_PQ_A_dot - w_AF % p_PQ_A; // derivative in F (12 flops)
}

// There is no findStationAccelerationFromState().

// Only configuration and velocity come from state; accelerations are an 
// operand (15 flops).
// p_BS_A      is p_BS rexpressed in A but not shifted to Ao.
// wXwX_p_BS_A is w_AB x (w_AB x p_BS_A)  (Coriolis acceleration)
Vec3 findStationInAAcceleration
   (const State&                                    s, 
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            B, 
    const Vec3&                                     p_BS_A,
    const Vec3&                                     wXwX_p_BS_A) const 
{   
    const SpatialVec& A_AB   = allA_AB[B];
    const Vec3& b_AB = A_AB[0]; const Vec3& a_AB = A_AB[1];

    // Result is a + b X r + w X (w X r).
    // ("b" is angular acceleration; w is angular velocity) 15 flops.
    const Vec3 a_AS = a_AB + (b_AB % p_BS_A) + wXwX_p_BS_A;
    return a_AS;
}

// Same as above but we only know the station in B. (48 flops)
Vec3 findStationAcceleration
   (const State&                                    s, 
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    ConstrainedBodyIndex                            B, 
    const Vec3&                                     p_BS) const 
{   // p_BS_A is p_BS rexpressed in A but not shifted to Ao
    const Rotation& R_AB   = getBodyRotationFromState(s,B);
    const Vec3&     w_AB   = getBodyAngularVelocityFromState(s,B);
    const Vec3 p_BS_A = R_AB * p_BS;                 // 15 flops
    const Vec3 wXwX_p_BS_A = w_AB % (w_AB % p_BS_A); // 18 flops
    return findStationInAAcceleration(s,allA_AB,B,p_BS_A,wXwX_p_BS_A);
}


// Combo method is cheaper. Location and velocity come from state, accelerations
// from operand. NOTE: you must provide the p_BS vector expressed (but not 
// measured) in A. 39 flops.
void findStationInALocationVelocityAcceleration
   (const Transform&                                X_AB,
    const SpatialVec&                               V_AB,
    const SpatialVec&                               A_AB,
    const Vec3&                                     p_BS_A,
    Vec3& p_AS, Vec3& v_AS, Vec3& a_AS) const 
{
    const Vec3& w_AB = V_AB[0]; const Vec3& v_AB = V_AB[1];
    const Vec3& b_AB = A_AB[0]; const Vec3& a_AB = A_AB[1];

    const Vec3 pdot_BS_A = w_AB % p_BS_A;   //  9 flops

    p_AS = X_AB.p() + p_BS_A;               //  3 flops
    v_AS = v_AB + pdot_BS_A;                //  3 flops

    // Result is a + b X r + w X (w X r).
    // ("b" is angular acceleration; w is angular velocity) 24 flops.
    a_AS = a_AB + (b_AB % p_BS_A) + (w_AB % pdot_BS_A);
}

// Apply an Ancestor-frame force to a B-frame station S, adding it to the 
// appropriate bodyForcesInA entry, where bodyForcesInA is *already* size 
// numConstrainedBodies for this Constraint. 30 flops.
void addInStationForce(const State& s, 
                       ConstrainedBodyIndex B, const Vec3& p_BS, 
                       const Vec3& forceInA, 
                       Array_<SpatialVec, ConstrainedBodyIndex>& bodyForcesInA) 
                       const 
{
    assert(bodyForcesInA.size() == getNumConstrainedBodies());
    const Rotation& R_AB = getBodyRotationFromState(s,B);
    const Vec3      p_BS_A = R_AB * p_BS;         // 15 flops
    bodyForcesInA[B] += SpatialVec(p_BS_A % forceInA, forceInA); // rXf, f
}

// If you already have the p_BS station vector re-expressed in A, use this
// faster method (15 flops).
void addInStationInAForce(const Vec3& p_BS_A, const Vec3& forceInA, 
                          SpatialVec& bodyForceOnBInA) 
                          const 
{
    bodyForceOnBInA += SpatialVec(p_BS_A % forceInA, forceInA); // rXf, f
}

// Same thing but subtract the force; this is just to save having to negate it.
void subInStationInAForce(const Vec3& p_BS_A, const Vec3& negForceInA, 
                          SpatialVec& bodyForceOnBInA) 
                          const 
{
    bodyForceOnBInA -= SpatialVec(p_BS_A % negForceInA, negForceInA); //-rXf,-f
}

// Apply an Ancestor-frame torque to body B, updating the appropriate 
// bodyForcesInA entry, where bodyForcesInA is *already* size 
// numConstrainedBodies for this Constraint. 3 flops.
void addInBodyTorque(const State& s, ConstrainedBodyIndex B,
                     const Vec3& torqueInA, 
                     Array_<SpatialVec, ConstrainedBodyIndex>& bodyForcesInA) 
                     const 
{
    assert(bodyForcesInA.size() == getNumConstrainedBodies());
    bodyForcesInA[B][0] += torqueInA; // no force
}

// After realizeTopology() we can look at the values of modeling variables in
// the State. A Constraint is free to use those in determining how many 
// constraint equations of each type to generate. The default implementation 
// doesn't look at the state but instead returns the default numbers of 
// equations supplied when the Constraint was constructed.
void calcNumConstraintEquationsInUse(const State& s, int& mp, int& mv, int& ma) const {
    calcNumConstraintEquationsInUseVirtual(s,mp,mv,ma);
}

// Abbreviated version of the above; returns number of holonomic constraint
// equations in use.
int calcNumPositionEquationsInUse(const State& s) const {
    int mp, mv, ma; calcNumConstraintEquationsInUse(s, mp, mv, ma);
    return mp;
}
// Abbreviated version of the above; returns number of nonholonomic constraint
// equations in use.
int calcNumVelocityEquationsInUse(const State& s) const {
    int mp, mv, ma; calcNumConstraintEquationsInUse(s, mp, mv, ma);
    return mv;
}
// Abbreviated version of the above; returns number of acceleration-only
// constraint equations in use.
int calcNumAccelerationEquationsInUse(const State& s) const {
    int mp, mv, ma; calcNumConstraintEquationsInUse(s, mp, mv, ma);
    return ma;
}

    // The next three methods are just interfaces to the constraint 
    // operators like calcPositionErrors(); they extract needed operands
    // from the supplied state and then call the operator.

// Calculate the mp position errors that would result from the configuration 
// present in the supplied state (that is, q's and body transforms). The state
// must be realized through Time stage and part way through realization of
// Position stage.
void calcPositionErrorsFromState(const State& s, Array_<Real>& perr) const;

// Calculate the mp velocity errors resulting from pdot equations, given a
// configuration and velocities in the supplied state which must be realized
// through Position stage and part way through realization of Velocity stage.
void calcPositionDotErrorsFromState(const State& s, Array_<Real>& pverr) const;

// Calculate the mv velocity errors resulting from the nonholonomic constraint
// equations here, taking the configuration and velocities (u, qdot, body
// spatial velocities) from the supplied state, which must be realized through
// Position stage and part way through realization of Velocity stage.
void calcVelocityErrorsFromState(const State& s, Array_<Real>& verr) const;

// Calculate position errors given pose of the constrained bodies and the
// values of the constrained q's. Pull t from state.
void calcPositionErrors     
   (const State&                                    s,     // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)  // mp of these
    const
{
    assert(X_AB.size()         == getNumConstrainedBodies());
    assert(constrainedQ.size() == getNumConstrainedQ(s));
    assert(perr.size()         == calcNumPositionEquationsInUse(s));

    calcPositionErrorsVirtual(s,X_AB,constrainedQ,perr);
}

// Calculate pdot errors given spatial velocity of the constrained bodies and 
// the values of the constrained qdot's. Pull t, X_AB and q from state.
void calcPositionDotErrors      
   (const State&                                    s, // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) // mp of these
    const
{
    assert(V_AB.size()            == getNumConstrainedBodies());
    assert(constrainedQDot.size() == getNumConstrainedQ(s));
    assert(pverr.size()           == calcNumPositionEquationsInUse(s));

    calcPositionDotErrorsVirtual(s,V_AB,constrainedQDot,pverr);
}

// Calculate pdotdot errors given spatial acceleration of the constrained 
// bodies and the values of the constrained qdotdot's. Pull t, X_AB, q, V_AB, 
// qdot from state.
void calcPositionDotDotErrors     
   (const State&                                    s, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) // mp of these
    const
{
    assert(A_AB.size()               == getNumConstrainedBodies());
    assert(constrainedQDotDot.size() == getNumConstrainedQ(s));
    assert(paerr.size()              == calcNumPositionEquationsInUse(s));

    calcPositionDotDotErrorsVirtual(s,A_AB,constrainedQDotDot,paerr);
}

// Given mp position constraint multipliers, generate the corresponding 
// constraint forces acting on the constrained bodies and the generalized
// coordinates q of the constrained mobilizers and add them in to the given 
// arrays. The state must be realized through Position stage.
void addInPositionConstraintForces
    (const State& s,
     const Array_<Real>&                       multipliers, // mp of these
     Array_<SpatialVec, ConstrainedBodyIndex>& bodyForcesInA,
     Array_<Real,       ConstrainedQIndex>&    qForces) const
{
    assert(multipliers.size()   == calcNumPositionEquationsInUse(s));
    assert(bodyForcesInA.size() == getNumConstrainedBodies());
    assert(qForces.size()       == getNumConstrainedQ(s));

    // Note that position constraints act on q (qdot,qdotdot) and for
    // convenience we allow them to generate generalized forces in qdot-space
    // rather than u-space. These will have to be mapped into u-space when
    // they are used, since Simbody only deals in u-space forces normally.
    // Since qdot=N*u, and we must have power ~f_qdot*qdot==~f_u*u, we have
    // f_u=~N * f_qdot.
    addInPositionConstraintForcesVirtual
       (s,multipliers,bodyForcesInA,qForces);
}

// Calculate velocity errors (errors produced by nonholonomic constraint
// equations) given spatial velocity of the constrained bodies and the
// values of the constrained u's. Pull t, X_AB, q from state.
void calcVelocityErrors     
   (const State&                                    s,    // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) // mv of these
    const
{
    assert(V_AB.size()         == getNumConstrainedBodies());
    assert(constrainedU.size() == getNumConstrainedU(s));
    assert(verr.size()         == calcNumVelocityEquationsInUse(s));

    calcVelocityErrorsVirtual(s,V_AB,constrainedU,verr);
}

// Calculate vdot errors (errors produced by nonholonomic constraint
// derivatives) given spatial acceleration of the constrained bodies and the
// values of the constrained udot's. Pull t, X_AB, q, V_AB, u from state.
void calcVelocityDotErrors     
   (const State&                                    s,     // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) // mv of these
    const
{
    assert(A_AB.size()            == getNumConstrainedBodies());
    assert(constrainedUDot.size() == getNumConstrainedU(s));
    assert(vaerr.size()           == calcNumVelocityEquationsInUse(s));

    calcVelocityDotErrorsVirtual(s,A_AB,constrainedUDot,vaerr);
}


// Given mv velocity constraint multipliers, generate the corresponding 
// constraint forces acting on the constrained bodies and the mobilities of the
// constrained mobilizers, and add them in to the given arrays. The state must 
// be realized through Velocity stage unless the V matrix is 
// velocity-independent in which case Position stage is enough.
void addInVelocityConstraintForces
    (const State& s,
     const Array_<Real>&                       multipliers, // mv of these
     Array_<SpatialVec, ConstrainedBodyIndex>& bodyForcesInA,
     Array_<Real,       ConstrainedUIndex>&    mobilityForces) const
{
    assert(multipliers.size()    == calcNumVelocityEquationsInUse(s));
    assert(bodyForcesInA.size()  == getNumConstrainedBodies());
    assert(mobilityForces.size() == getNumConstrainedU(s));

    addInVelocityConstraintForcesVirtual
       (s,multipliers,bodyForcesInA,mobilityForces);
}

// Calculate acceleration errors (errors produced by acceleration-only
// constraint equations) given spatial acceleration of the constrained bodies
// and the values of the constrained udot's. Pull t, X_AB, q, V_AB, u from state.
void calcAccelerationErrors      
   (const State&                                    s,    // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr) // ma of these
    const
{
    assert(A_AB.size()            == getNumConstrainedBodies());
    assert(constrainedUDot.size() == getNumConstrainedU(s));
    assert(aerr.size()            == calcNumAccelerationEquationsInUse(s));

    calcAccelerationErrorsVirtual(s,A_AB,constrainedUDot,aerr);
}


// Given ma acceleration constraint multipliers, generate the corresponding 
// constraint forces acting on the constrained bodies and the mobilities of the
// constrained mobilizers, and add them in to the given arrays. The state must 
// be realized through Velocity stage unless the A matrix is 
// velocity-independent in which case Position stage is enough.
void addInAccelerationConstraintForces
    (const State& s,
     const Array_<Real>&                       multipliers, // ma of these
     Array_<SpatialVec, ConstrainedBodyIndex>& bodyForcesInA,
     Array_<Real,       ConstrainedUIndex>&    mobilityForces) const
{
    assert(multipliers.size()    == calcNumAccelerationEquationsInUse(s));
    assert(bodyForcesInA.size()  == getNumConstrainedBodies());
    assert(mobilityForces.size() == getNumConstrainedU(s));

    addInAccelerationConstraintForcesVirtual
       (s,multipliers,bodyForcesInA,mobilityForces);
}

void calcDecorativeGeometryAndAppend
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // Let the individual constraint deal with any complicated stuff.
    calcDecorativeGeometryAndAppendVirtual(s,stage,geom);
}

    // Don't call these virtuals directly; use the provided interface
    // methods for safety (they evaporate in Release builds anyway).

virtual void calcNumConstraintEquationsInUseVirtual
   (const State&, int& mp, int& mv, int& ma) const 
{   mp = defaultMp; mv = defaultMv; ma = defaultMa; }

virtual void realizeTopologyVirtual     (State&)        const {}
virtual void realizeModelVirtual        (State&)        const {}
virtual void realizeInstanceVirtual     (const State&)  const {}
virtual void realizeTimeVirtual         (const State&)  const {}
virtual void realizePositionVirtual     (const State&)  const {}
virtual void realizeVelocityVirtual     (const State&)  const {}
virtual void realizeDynamicsVirtual     (const State&)  const {}
virtual void realizeAccelerationVirtual (const State&)  const {}
virtual void realizeReportVirtual       (const State&)  const {}

    // These must be defined if there are any position (holonomic) constraints.

// Pull t from state.
virtual void calcPositionErrorsVirtual      
   (const State&                                    state, // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)  // mp of these
    const;

// Pull t, X_AB and q from state.
virtual void calcPositionDotErrorsVirtual      
   (const State&                                    state, // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) // mp of these
    const;

// Pull t, X_AB, q, V_AB, qdot from state.
virtual void calcPositionDotDotErrorsVirtual      
   (const State&                                    state, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) // mp of these
    const;

// Pull t, X_AB and q from state.
virtual void addInPositionConstraintForcesVirtual
   (const State&                                    state, // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const;

    // These must be defined if there are velocity (nonholonomic) constraints.

// Pull t, X_AB, q from state.
virtual void calcVelocityErrorsVirtual      
   (const State&                                    state, // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) // mv of these
    const;

// Pull t, X_AB, q, V_AB, u from state.
virtual void calcVelocityDotErrorsVirtual      
   (const State&                                    state, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) // mv of these
    const;

// Pull t, X_AB, q, V_AB, u from state.
virtual void addInVelocityConstraintForcesVirtual
   (const State&                                    state, // Stage::Velocity
    const Array_<Real>&                             multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const;

    // These must be defined if there are any acceleration-only constraints.

// Pull t, X_AB, q, V_AB, u from state.
virtual void calcAccelerationErrorsVirtual      
   (const State&                                    state, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr) // ma of these
    const;

// Pull t, X_AB, q, V_AB, u from state.
virtual void addInAccelerationConstraintForcesVirtual
   (const State&                                    state, // Stage::Velocity
    const Array_<Real>&                             multipliers, // ma of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) 
    const;


virtual void calcDecorativeGeometryAndAppendVirtual
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const {}


void invalidateTopologyCache() const;
bool subsystemTopologyHasBeenRealized() const;

void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                          ConstraintIndex id);

const SimbodyMatterSubsystem& getMyMatterSubsystem() const;

bool isInSubsystem() const {
    return myMatterSubsystemRep != 0;
}

// Is the supplied body in the same subsystem as this Constraint? (Returns 
// false also if either the Constraint or the MobilizedBody is not in a 
// subsystem.)
bool isInSameSubsystem(const MobilizedBody& body) const;

int getNumConstrainedBodies() const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "Number of constrained bodies is not available until Topology stage has been realized.");
    return (int)myConstrainedBodies.size();
}
int getNumConstrainedMobilizers() const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "Number of constrained mobilizers is not available until Topology stage has been realized.");
    return (int)myConstrainedMobilizers.size();
}

const MobilizedBody& 
    getMobilizedBodyFromConstrainedMobilizer(ConstrainedMobilizerIndex) const;
const MobilizedBody& 
    getMobilizedBodyFromConstrainedBody(ConstrainedBodyIndex) const;

// Don't call this unless there is at least one Constrained Body.
const MobilizedBody& getAncestorMobilizedBody() const;

// After realizeTopology() we remember whether this constraint's Ancestor
// body is different from Ground.
bool isAncestorDifferentFromGround() const {
    SimTK_ASSERT(subsystemTopologyHasBeenRealized(),
        "isAncestorDifferentFromGround(): must call realizeTopology() first.");
    return myAncestorBodyIsNotGround;
}

// Find out how many holonomic (position), nonholonomic (velocity),
// and acceleration-only constraint equations are generated by this Constraint
// as currently modeled. State must be realized to Stage::Model.
void getNumConstraintEquationsInUse
   (const State&, int& mHolo, int& mNonholo, int& mAccOnly) const;

// Return the starting index within the multiplier or udot error arrays
// for each of the acceleration-level constraints produced by this Constraint.
// Let holo0, nonholo0, accOnly0 be the first index of the slots assigned to
// this Constraint's holonomic, nonholonomic, and acceleration-only constraints
// within each block for that category. Then the returns here are:
//     px = holo0
//     vx = totalNumHolo + nonholo0
//     ax = totalNumHolo+totalNumNonholo + accOnly0
// These are returned invalid if there are no constraint equations in that
// category.
void getIndexOfMultipliersInUse(const State& state,
                                MultiplierIndex& px0, 
                                MultiplierIndex& vx0, 
                                MultiplierIndex& ax0) const;

void setMyPartInConstraintSpaceVector(const State& state,
                                 const Vector& myPart,
                                 Vector& constraintSpace) const;

void getMyPartFromConstraintSpaceVector(const State& state,
                                   const Vector& constraintSpace,
                                   Vector& myPart) const;

int getNumConstrainedQ(const State&) const;
int getNumConstrainedU(const State&) const;
int getNumConstrainedQ(const State&, ConstrainedMobilizerIndex) const;
int getNumConstrainedU(const State&, ConstrainedMobilizerIndex) const;
ConstrainedQIndex getConstrainedQIndex
    (const State&, ConstrainedMobilizerIndex, MobilizerQIndex which) const;
ConstrainedUIndex getConstrainedUIndex
    (const State&, ConstrainedMobilizerIndex, MobilizerUIndex which) const;

const SimbodyMatterSubsystemRep& getMyMatterSubsystemRep() const {
    SimTK_ASSERT(myMatterSubsystemRep,
        "Operation illegal on a Constraint that is not in a Subsystem.");
    return *myMatterSubsystemRep;
}
SimbodyMatterSubsystemRep& updMyMatterSubsystemRep() {
    SimTK_ASSERT(myMatterSubsystemRep,
        "Operation illegal on a Constraint that is not in a Subsystem.");
    return *myMatterSubsystemRep;
}

ConstraintIndex getMyConstraintIndex() const {
    SimTK_ASSERT(myMatterSubsystemRep,
        "Operation illegal on a Constraint that is not in a Subsystem.");
    return myConstraintIndex;
}


// Calculate the transform X_AB of each ConstrainedBody in its Ancestor frame, 
// provided Ancestor!=Ground and A!=B. We expect a TreePositionCache in
// which the mobilizer- and ground-frame position kinematics results have
// already been calculated. We then fill in the missing ancestor-frame
// results back into that same TreePositionCache.
void calcConstrainedBodyTransformInAncestor(const SBInstanceVars&,  // in only
                                            SBTreePositionCache&    // in/out
                                            ) const;

// Similarly we calculate V_AB during realizeVelocity().
// Here we expect a StateDigest realized through Position stage, and a 
// partly-filled-in VelocityCache where we'll put V_AB for the 
// ConstrainedBodies.
void calcConstrainedBodyVelocityInAncestor(const SBInstanceVars&,   // in only
                                           const SBTreePositionCache&, // "
                                           SBTreeVelocityCache&     // in/out
                                           ) const;

// A_AB is not cached.

private:
friend class Constraint;

    // TOPOLOGY "STATE"

// These data members are filled in once the Constraint is added to
// a MatterSubsystem.
SimbodyMatterSubsystemRep* myMatterSubsystemRep;
ConstraintIndex            myConstraintIndex; // id within the matter subsystem

// We'll keep the constrained bodies and constrained mobilizers each in two 
// maps: one maps MobilizedBodyIndex->ConstrainedBody[Mobilizer]Index (O(log n)
// to look up), and the other maps ConstrainedBody[Mobilizer]Index->
// MobilizedBodyIndex (randomly addressable in constant time).
MobilizedBody2ConstrainedBodyMap        myMobilizedBody2ConstrainedBodyMap;
MobilizedBody2ConstrainedMobilizerMap   myMobilizedBody2ConstrainedMobilizerMap;

Array_<MobilizedBodyIndex> myConstrainedBodies;     // [ConstrainedBodyIndex]
Array_<MobilizedBodyIndex> myConstrainedMobilizers; // [ConstrainedMobilizerIndex]


// These are the defaults for the number of position (holonomic) constraint 
// equations, the number of velocity (nonholonomic) constraint equations, and 
// the number of acceleration-only constraint equations.
int defaultMp, defaultMv, defaultMa;

// This says whether the Model-stage "disabled" flag for this Constraint should
// be initially on or off. Most constraints are enabled by default.
bool defaultDisabled;

// ConditionalConstraints set this flag when they add a constraint to the
// system. It can be referenced by a time stepper to determine whether to
// treat a constraint as unconditional (in which case someone else has to 
// figure out whether it is active). This flag doesn't affect the operation of
// the Constraint object itself; it is just stored and reported.
bool constraintIsConditional;

    // TOPOLOGY "CACHE"

// When topology is realized we study the constrained bodies to identify the
// subtree of mobilized bodies which may be kinematically involved in 
// satisfaction of this Constraint. This requires finding the outmost common 
// ancestor of the constrained bodies. All mobilized bodies on the paths inward
// from the constrained bodies to the ancestor are included; nothing outboard 
// of the constrained bodies is included; and the ancestor is treated as ground
// so that its mobilities are *not* included. The Ancestor may be one of the
// Constrained Bodies or may be distinct.
mutable SimbodyMatterSubtree mySubtree;

// This is true only when (1) there are Constrained Bodies and (2) their 
// Ancestor is some MobilizedBody other than Ground.
mutable bool myAncestorBodyIsNotGround;

// When the Ancestor is not Ground, each of the Constrained Bodies (except
// the Ancestor) needs some additional precalculated data in the State cache
// so is assigned a MatterSubsystem-global AncestorConstrainedBodyPool slot.
// If Ancestor isn't Ground there is an entry here for each ConstrainedBody
// (including Ancestor if it is one), but the index is invalid for Ancestor's
// entry. When Ancestor is Ground we don't allocate this array.
mutable Array_<AncestorConstrainedBodyPoolIndex> myPoolIndex; 
                                            // index with ConstrainedBodyIndex
};



//==============================================================================
//                           POINT IN PLANE IMPL
//==============================================================================
class Constraint::PointInPlaneImpl : public ConstraintImpl {
public:
PointInPlaneImpl()
    : ConstraintImpl(1,0,0), defaultPlaneNormal(), defaultPlaneHeight(0), defaultFollowerPoint(0),
    planeHalfWidth(1), pointRadius(Real(0.05)) 
{ }
PointInPlaneImpl* clone() const override { return new PointInPlaneImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setPlaneDisplayHalfWidth(Real h) {
    // h <= 0 means don't display plane
    invalidateTopologyCache();
    planeHalfWidth = h > 0 ? h : 0;
}
Real getPlaneDisplayHalfWidth() const {return planeHalfWidth;}

void setPointDisplayRadius(Real r) {
    // r <= 0 means don't display point
    invalidateTopologyCache();
    pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return pointRadius;}

// Implementation of virtuals required for holonomic constraints.

// We have a point-in-plane connection between base body B, on which the plane 
// is fixed, and follower body F, on which the follower point S is fixed. All 
// forces will be applied at point S and the coincident material point C on B 
// which is instantaneously at the same spatial location as S. Then n is the 
// plane normal (a constant unit vector in B), h is the plane height measured 
// from the B origin along n (a scalar constant).Point C's location in B is 
// given by the vector p_BC from B's origin to the current location of S, and 
// expressed in B. That vector expressed in A is p_BC_A (= p_AS-p_AB). We will 
// express in the A frame but differentiate in the B frame.
//
// Derivation:
//   (1) Note that to take a derivative d/dt_B in a moving frame B, we can take
//       the derivative d/dt_A and then add in the contribution d_A/dt_B from 
//       A's rotation in B (which is the angular velocity of A in B, 
//       w_BA=-w_AB).
//   (2) p_CS = p_AS-p_AC = 0 by construction of C, but its derivative in A, 
//       v_CS_A = d/dt_A p_CS != 0.
//
//    perr = p_CS * n + constant 
//         = constant  (because p_CS==0 by construction)
//
//    verr = d/dt_B perr = d/dt_A perr + d_A/dt_B perr
//         = [v_CS_A*n + p_CS * (w_AB X n)] + [(w_BA X p_CS) * n + p_CS * (w_BA X n)]
//         = v_CS_A*n + p_CS * (w_AB X n) (because terms in 2nd [] cancel)
//         = v_CS_A * n  (because p_CS==0 by construction)
//
//    aerr = d/dt_B verr = d/dt_A verr + d_A/dt_B verr
//         = [a_CS_A*n + v_CS_A*(w_AB X n) + v_CS_A*(w_AB X n) + p_CS*(2 w_AB X(w_AB X n))]
//           + [w_BAXv_CS_A*n + v_CS_A*w_BAXn]
//         = (a_CS_A - 2 w_AB X v_CS_A) * n   (2nd bracket cancels, and p_CS==0)
//
// (The constant in perr is set so that S starts at the same height h as the 
// plane.)
//  
// Then, from examination of verr noting that v_CS_A=v_AS-v_AC:
//       ~v_AS*n                  (body F at point S) 
//     - ~v_AC*n                  (body B at point C)
// so we apply a forces lambda*n to F at S, -lambda*n to B at C.
//
//    --------------------------------
//    perr = ~p_BS*n - h
//    --------------------------------
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Vec3       p_AS = findStationLocation(allX_AB, followerBody, 
                                                defaultFollowerPoint);
    const Transform& X_AB = getBodyTransform(allX_AB, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                     // C is material pt of B coincident with S

    // We'll calculate this scalar using B-frame vectors, but any frame would 
    // have done.
    perr[0] = dot(p_BC, defaultPlaneNormal) - defaultPlaneHeight;
}

//    --------------------------------
//    verr = ~v_CS_A*n
//    --------------------------------
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(V_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    //TODO: should be able to get p info from State

    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3       v_AS = findStationVelocity(s, V_AB, followerBody, 
                                                defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocity(s, V_AB, planeBody, p_BC);

    // Calculate this scalar using A-frame vectors.
    pverr[0] = dot( v_AS-v_AC, n_A );
}

//    -------------------------------------
//    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
//    -------------------------------------
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(A_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size() == 1);
    //TODO: should be able to get p and v info from State
    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3&      w_AB = getBodyAngularVelocityFromState(s, planeBody);
    const Vec3       v_AS = findStationVelocityFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocityFromState(s, planeBody, p_BC);

    const Vec3       a_AS = findStationAcceleration(s, A_AB, followerBody, 
                                                    defaultFollowerPoint);;
    const Vec3       a_AC = findStationAcceleration(s, A_AB, planeBody, p_BC);

    paerr[0] = dot( (a_AS-a_AC) - 2.*w_AB % (v_AS-v_AC), n_A );
}

// apply f=lambda*n to the follower point S of body F,
//       -f         to point C (coincident point) of body B
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==1 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Real lambda = multipliers[0];

    //TODO: should be able to get p info from State
    const Vec3&      p_FS    = defaultFollowerPoint; // measured & expressed in F
    const Vec3       p_AS    = findStationLocationFromState(s, followerBody, 
                                                            defaultFollowerPoint);
    const Transform& X_AB    = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC    = ~X_AB * p_AS;         // measured & expressed in B
    const Vec3       force_A = X_AB.R()*(lambda*defaultPlaneNormal);

    addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
    addInStationForce(s, planeBody,    p_BC, -force_A, bodyForcesInA);
}

SimTK_DOWNCAST(PointInPlaneImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::PointInPlane;

ConstrainedBodyIndex    planeBody;    // B1
ConstrainedBodyIndex    followerBody; // B2

UnitVec3                defaultPlaneNormal;   // on body 1, exp. in B1 frame
Real                    defaultPlaneHeight;
Vec3                    defaultFollowerPoint; // on body 2, exp. in B2 frame

// These are just for visualization
Real                    planeHalfWidth;
Real                    pointRadius;
};



//==============================================================================
//                           POINT ON LINE IMPL
//==============================================================================
class Constraint::PointOnLineImpl : public ConstraintImpl {
public:
PointOnLineImpl()
:   ConstraintImpl(2,0,0), 
    defaultLineDirection(), defaultPointOnLine(), defaultFollowerPoint(0),
    lineHalfLength(1), pointRadius(Real(0.05)) 
{ }
PointOnLineImpl* clone() const override { return new PointOnLineImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setLineDisplayHalfLength(Real h) {
    // h <= 0 means don't display line
    invalidateTopologyCache();
    lineHalfLength = h > 0 ? h : 0;
}
Real getLineDisplayHalfLength() const {return lineHalfLength;}

void setPointDisplayRadius(Real r) {
    // r <= 0 means don't display point
    invalidateTopologyCache();
    pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return pointRadius;}

// Implementation of ContraintRep virtuals
void realizeTopologyVirtual(State& s) const override {
    x = defaultLineDirection.perp(); // x and y are mutable
    y = UnitVec3(defaultLineDirection % x);
}

// Implementation of virtuals required for holonomic constraints.

// We have a point-on-line connection between base body B, on which the line is
// fixed, and follower body F, on which the follower point S is fixed. All 
// forces will be applied at point S and the coincident material point C on B 
// which is instantaneously at the same spatial location as S. Then z is a unit
// vector in the direction of the line, and P is a point fixed to B that the 
// line passes through. We will enforce this using two point-on-plane 
// constraints, where the intersection of the two planes is the line. For that 
// we need two plane normals perpendicular to z. We'll use an arbitrary 
// perpendicular x, then use y=z X x as the other perpendicular. This 
// establishes a right handed coordinate system where the line is along the z 
// axis, and we'll apply constraint forces in the x-y plane.
//
// See the point-in-plane constraint for details; here we're just picking x and
// y as plane normals.

//    --------------------------------
//    perr = ~(p_BS-p_BP) * x
//           ~(p_BS-p_BP) * y
//    --------------------------------
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 2);

    const Transform& X_AB = getBodyTransform(allX_AB, lineBody);
    const Vec3       p_AS = findStationLocation(allX_AB, followerBody, 
                                                defaultFollowerPoint);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const Vec3       p_PC = p_BC - defaultPointOnLine;

    // We'll calculate these two scalars using B-frame vectors, but any frame 
    // would have done.
    perr[0] = ~p_PC * x;
    perr[1] = ~p_PC * y;
}

//    --------------------------------
//    verr = ~v_CS_A*n
//    --------------------------------
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 2);
    //TODO: should be able to get p info from State
    const Transform& X_AB = getBodyTransformFromState(s, lineBody);
    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Vec3       p_BC = ~X_AB * p_AS;
    const Vec3       p_PC = p_BC - defaultPointOnLine;

    const Vec3       v_AS = findStationVelocity(s, allV_AB, followerBody, 
                                                defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocity(s, allV_AB, lineBody, p_BC);

    const Vec3       v_CS_B = ~X_AB.R()*(v_AS-v_AC); // reexpress in B

    // Calculate these scalar using B-frame vectors, but any frame would 
    // have done.
    pverr[0] = ~v_CS_B * x;
    pverr[1] = ~v_CS_B * y;
}

//    -------------------------------------
//    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
//    -------------------------------------
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==2);
    //TODO: should be able to get p and v info from State
    const Transform& X_AB = getBodyTransformFromState(s, lineBody);
    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                 // C is material point of B coincident with S
    const Vec3       p_PC = p_BC - defaultPointOnLine;

    const Vec3&      w_AB = getBodyAngularVelocityFromState(s, lineBody);
    const Vec3       v_AS = findStationVelocityFromState(s, followerBody, defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocityFromState(s, lineBody, p_BC);

    const Vec3       a_AS = findStationAcceleration(s, allA_AB, followerBody, 
                                                    defaultFollowerPoint);
    const Vec3       a_AC = findStationAcceleration(s, allA_AB, lineBody, p_BC);
    const Vec3       a_CS_B = ~X_AB.R()*(a_AS-a_AC - 2 * w_AB % (v_AS-v_AC));

    // Calculate these scalar using B-frame vectors, but any frame would 
    // have done.
    paerr[0] = ~a_CS_B * x;
    paerr[1] = ~a_CS_B * y;
}

// apply f=lambda0*x + lambda1*y to the follower point S of body F,
//      -f                       to point C (coincident point) of body B
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==2 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Vec2 lambda = Vec2::getAs(&multipliers[0]);

    //TODO: should be able to get p info from State
    const Transform& X_AB    = getBodyTransformFromState(s, lineBody);
    const Vec3&      p_FS    = defaultFollowerPoint; // measured & expressed in F
    const Vec3       p_AS    = findStationLocationFromState(s, followerBody, 
                                                            defaultFollowerPoint);
    const Vec3       p_BC    = ~X_AB * p_AS;         // measured & expressed in B

    const Vec3       force_B = lambda[0] * x + lambda[1] * y;
    const Vec3       force_A = X_AB.R() * force_B;

    addInStationForce(s, followerBody, p_FS,  force_A, bodyForcesInA);
    addInStationForce(s, lineBody,     p_BC, -force_A, bodyForcesInA);
}

SimTK_DOWNCAST(PointOnLineImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::PointOnLine;

ConstrainedBodyIndex    lineBody;     // B
ConstrainedBodyIndex    followerBody; // F

UnitVec3                defaultLineDirection;   // z on B, exp. in B frame
Vec3                    defaultPointOnLine;     // P on B, meas&exp in B frame
Vec3                    defaultFollowerPoint;   // S on F, meas&exp in F frame

// These are just for visualization
Real                    lineHalfLength;
Real                    pointRadius;

// TOPOLOGY CACHE (that is, calculated from construction data)
mutable UnitVec3        x, y;
};



//==============================================================================
//                           CONSTANT ANGLE IMPL
//==============================================================================
class Constraint::ConstantAngleImpl : public ConstraintImpl {
public:
ConstantAngleImpl()
    : ConstraintImpl(1,0,0), defaultAxisB(), defaultAxisF(), defaultAngle(Pi/2),
    axisLength(1), axisThickness(1), cosineOfDefaultAngle(NaN)
{ }
ConstantAngleImpl* clone() const override { return new ConstantAngleImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setAxisLength(Real length) {
    // length <= 0 means don't display axis
    invalidateTopologyCache();
    axisLength = length > 0 ? length : 0;
}
Real getAxisLength() const {return axisLength;}

void setAxisThickness(Real t) {
    // t <= 0 means don't display axis
    invalidateTopologyCache();
    axisThickness = t > 0 ? t : 0;
}
Real getAxisThickness() const {return axisThickness;}

// Implementation of ContraintRep virtuals
void realizeTopologyVirtual(State& s) const override {
    cosineOfDefaultAngle = std::cos(defaultAngle);
}


// Implementation of virtuals required for holonomic constraints.

// Let B=B1 be the "base" body onto which unit vector b is fixed, and F=B2 the 
// "follower" body onto which unit vector f is fixed. The angle theta between 
// these vectors is given by cos(theta) = dot(b, f) with the axes expressed in 
// a common basis. This can range from 1 to -1, corresponding to angles 0 to 
// 180 respectively. We would like to enforce the constraint that cos(theta) is
// a constant. This can be done with a single constraint equation as long as 
// theta is sufficiently far away from 0 and 180, with the numerically best 
// performance at theta=90 degrees where cos(theta)==0.
//
// If you want to enforce that two axes are aligned with one another (that is, 
// the angle between them is 0 or 180), that takes *two* constraint equations 
// since the only remaining rotation is about the common axis.
//
// We will work in the A frame.
//
// ------------------------------
// perr = ~b_A * f_A - cos(theta)
// ------------------------------
//
// verr = d/dt perr (derivative taken in A)
//      = ~b_A * (w_AF % f_A) + ~f_A * (w_AB % b_A)
//      = ~w_AF * (f_A % b_A) - ~w_AB * (f_A % b_A)  (scalar triple product identity)
// => ------------------------------
// verr = ~(w_AF-w_AB) * (f_A % b_A)
// ---------------------------------
//
// aerr = d/dt verr (derivative taken in A)
//      = ~(b_AF-b_AB) * (f_A % b_A)
//        + (w_AF-w_AB) * ((w_AF%f_A) % b_A)
//        + (w_AF-w_AB) * (f_A % (w_AB%b_A))
//      =   ~(b_AF-b_AB) * (f_A % b_A)
//        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
// => -----------------------------------------------------------
// aerr =   ~(b_AF-b_AB) * (f_A % b_A)
//        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
// --------------------------------------------------------------
//
// Constraint torque can be determined by inspection of verr:
//    lambda * (f_A % b_A) applied to body F
//   -lambda * (f_A % b_A) applied to body B
//

// ------------------------------
// perr = ~b_A * f_A - cos(theta)
// ------------------------------
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Rotation& R_AB = getBodyRotation(allX_AB, B);
    const Rotation& R_AF = getBodyRotation(allX_AB, F);
    const UnitVec3  b_A  = R_AB * defaultAxisB;
    const UnitVec3  f_A  = R_AF * defaultAxisF;

    perr[0] = dot(b_A, f_A) - cosineOfDefaultAngle;
}

// ----------------------------------
// pverr = ~(w_AF-w_AB) * (f_A % b_A)
// ----------------------------------
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size()==1);
    //TODO: should be able to get p info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const UnitVec3  b_A  = R_AB * defaultAxisB;
    const UnitVec3  f_A  = R_AF * defaultAxisF;

    const Vec3&     w_AB = getBodyAngularVelocity(allV_AB, B);
    const Vec3&     w_AF = getBodyAngularVelocity(allV_AB, F);

    pverr[0] = dot( w_AF-w_AB,  f_A % b_A );
}

// --------------------------------------------------------------
// paerr =  ~(b_AF-b_AB) * (f_A % b_A)
//        + ~(w_AF-w_AB) * ((w_AF%f_A) % b_A) - (w_AB%b_A) % f_A)
// --------------------------------------------------------------
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==1);
    //TODO: should be able to get p and v info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const UnitVec3  b_A  = R_AB * defaultAxisB;
    const UnitVec3  f_A  = R_AF * defaultAxisF;
    const Vec3&     w_AB = getBodyAngularVelocityFromState(s, B);
    const Vec3&     w_AF = getBodyAngularVelocityFromState(s, F);

    const Vec3&     b_AB = getBodyAngularAcceleration(allA_AB, B);
    const Vec3&     b_AF = getBodyAngularAcceleration(allA_AB, F);

    paerr[0] =    dot( b_AF-b_AB, f_A % b_A )
                + dot( w_AF-w_AB, (w_AF%f_A) % b_A - (w_AB%b_A) % f_A);
}

//    lambda * (f_A % b_A) applied to body F
//   -lambda * (f_A % b_A) applied to body B
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==1 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Real lambda = multipliers[0];
    //TODO: should be able to get p info from State
    const Rotation&  R_AB = getBodyRotationFromState(s, B);
    const Rotation&  R_AF = getBodyRotationFromState(s, F);
    const UnitVec3   b_A = R_AB*defaultAxisB;
    const UnitVec3   f_A = R_AF*defaultAxisF;
    const Vec3       torque_F_A = lambda * (f_A % b_A); // on F, in A frame

    addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
    addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);
}

SimTK_DOWNCAST(ConstantAngleImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::ConstantAngle;

ConstrainedBodyIndex    B; // B1 is "base" body
ConstrainedBodyIndex    F; // B2 is "follower" body

UnitVec3                defaultAxisB; // fixed to B, expressed in B frame
UnitVec3                defaultAxisF; // fixed to F, expressed in F frame
Real                    defaultAngle; // required angle between axisB and axisF

// These are just for visualization
Real                    axisLength;
Real                    axisThickness;

// TOPOLOGY CACHE (that is, calculated from construction data)
mutable Real            cosineOfDefaultAngle;
};



//==============================================================================
//                                 BALL IMPL
//==============================================================================
class Constraint::BallImpl : public ConstraintImpl {
public:
BallImpl() : ConstraintImpl(3,0,0), defaultPoint1(0), defaultPoint2(0), 
             defaultRadius(Real(0.1)) { }
BallImpl* clone() const override { return new BallImpl(*this); }

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setDefaultRadius(Real r) {
    // r <= 0 means don't display
    invalidateTopologyCache();
    defaultRadius = r > 0 ? r : 0;
}
Real getDefaultRadius() const {return defaultRadius;}

// The default body stations may be overridden by setting instance variables
// in the state. We allocate the state resources here.
void realizeTopologyVirtual(State& state) const override;

// Return the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Note that although
// these are used to define the position error, only the station on body 2
// is used to generate constraint forces; the point of body 1 that is 
// coincident with the body 2 point receives the equal and opposite force.
const std::pair<Vec3,Vec3>& getBodyStations(const State& state) const;

// Return a writable reference into the Instance-stage state variable 
// containing the pair of constrained station points, with the first expressed 
// in the body 1 frame and the second in the body 2 frame. Calling this
// method invalidates the Instance stage and above in the given state.
std::pair<Vec3,Vec3>& updBodyStations(State& state) const;

// Implementation of virtuals required for holonomic constraints.

// We have a ball joint between base body B and follower body F, located at a 
// point P fixed to B and point S fixed on F. All forces will be applied at 
// point S and the coincident material point C on B which is instantaneously at
// the same spatial location as S. We will work in the A frame.
//
//  First, find the material point C of B that is coincident
//  in space with point S of F: p_BC = p_AS-p_AB. This vector
//  is *constant* in the B frame because it is a material point,
//  despite the fact that its definition involves a point which
//  moves with respect to B. The velocity constraint is then
//  very simple: the spatial velocity of point C of B should be
//  the same as the spatial velocity of point S of F:
//      verr = v_AS - v_AC = v_AS - (v_AB + w_AB X p_BC) = 0
//  Integrating to get perr, we get
//      perr = p_AS - p_AC + constant = 0
//  But p_AC=p_AS by construction, so perr=constant=0.
//  The constant is defined by the fact that we want material point
//  C of B to be in the same spatial location as point P of B,
//  so constant=p_BC-p_BP=0. Writing in the A frame we have:
//      perr = p_AS-(p_AB+R_AB*p_BP) = 0 (a constant)
//      verr = v_AS - (v_AB + w_AB X R_AB*p_BC)
//      aerr = a_AS - (a_AB + b_AB X R_AB*p_BC + w_AB X (w_AB X R_AB*p_BC))
//  apply +lambda to S of F, -lambda to C of B.
//      
//

void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 3);

    const std::pair<Vec3,Vec3>& pts = getBodyStations(s);

    const Vec3 p_AP = findStationLocation(allX_AB, B1, pts.first);
    const Vec3 p_AS = findStationLocation(allX_AB, B2, pts.second);

    // See above comments -- this is just the constant of integration; there is a missing (p_AS-p_AC)
    // term (always 0) here which is what we differentiate to get the verr equation.
    Vec3::updAs(&perr[0]) = p_AS - p_AP;
}
 
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size()==3);

    // Note that we're not using point1.
    const Vec3& point2 = getBodyStations(s).second;

    //TODO: should be able to get p info from State
    const Transform&  X_AB   = getBodyTransformFromState(s, B1);
    const Vec3        p_AS   = findStationLocationFromState(s, B2, point2);
    const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

    const Vec3        v_AS    = findStationVelocity(s, allV_AB, B2, 
                                                    point2);
    const Vec3        v_AC    = findStationVelocity(s, allV_AB, B1, p_BC);
    Vec3::updAs(&pverr[0]) = v_AS - v_AC;
}

void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==3);

    // Note that we're not using point1.
    const Vec3& point2 = getBodyStations(s).second;

    //TODO: should be able to get p and v info from State

    const Transform&  X_AB   = getBodyTransformFromState(s, B1);
    const Vec3        p_AS   = findStationLocationFromState(s, B2, point2);
    const Vec3        p_BC   = ~X_AB*p_AS; // C is a material point of body B

    const Vec3        a_AS    = findStationAcceleration(s, allA_AB, B2, 
                                                        point2);
    const Vec3        a_AC    = findStationAcceleration(s, allA_AB, B1, p_BC);
    Vec3::updAs(&paerr[0]) = a_AS - a_AC;
}

void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==3 && bodyForcesInA.size()==2 
           && qForces.size()==0);

    // Note that we're not using point1.
    const Vec3& point2 = getBodyStations(s).second;

    //TODO: should be able to get p info from State
    const Transform& X_AB  = getBodyTransformFromState(s,B1);
    const Vec3&      p_FS  = point2;
    const Vec3       p_AS  = findStationLocationFromState(s, B2, p_FS);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                  // C is material point of B coincident with S

    const Vec3 force_A = Vec3::getAs(&multipliers[0]);

    // Multipliers are force to be applied to S on F, but
    // apply the -force not to point P of B, but to the point "C" of B
    // coincident with S, which won't be exactly the same place
    // as P if the position-level constraint isn't met exactly.

    addInStationForce(s, B2, p_FS,  force_A, bodyForcesInA);
    addInStationForce(s, B1, p_BC, -force_A, bodyForcesInA);
}

SimTK_DOWNCAST(BallImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::Ball;

ConstrainedBodyIndex    B1;
ConstrainedBodyIndex    B2;

Vec3                    defaultPoint1; // on body 1, exp. in B1 frame
Vec3                    defaultPoint2; // on body 2, exp. in B2 frame
Real                    defaultRadius; // used for visualization only

// This Instance-stage variable holds the actual stations on B1 & B2.
mutable DiscreteVariableIndex   stationsIx;
};



//==============================================================================
//                          CONSTANT ORIENTATION IMPL
//==============================================================================
class Constraint::ConstantOrientationImpl : public ConstraintImpl {
public:
ConstantOrientationImpl()
    : ConstraintImpl(3,0,0), defaultRB(), defaultRF()
{ }
ConstantOrientationImpl* clone() const override { return new ConstantOrientationImpl(*this); }

//TODO: visualization?


// Implementation of virtuals required for holonomic constraints.

// Let B=B1 be the "base" body onto which rotation matrix RB is fixed, and F=B2
// the "follower" body onto which rotation matrix RF is fixed. We would like to
// enforce the constraint that RB==RF when both are expressed in a common basis.
// (Remember that a rotation matrix is just three axis vectors.)
// 
// Here the (redundant) assembly constraint is that all the axes are parallel, 
// that is RBx==RFx, RBy==RFy, and RBz==RFz. However, aligning two vectors 
// takes *two* constraints so that would be a total of 6 constraints, with only
// 3 independent. The independent runtime constraints just enforce 
// perpendicularity, but can be satisfied in cases where some of the axes are 
// antiparallel so are not suited for the initial assembly. The runtime 
// constraints are thus three "constant angle" constraints, where the angle
// is always 90 degrees:
//
//    ~RFx * RBy = 0
//    ~RFy * RBz = 0
//    ~RFz * RBx = 0
//
// We'll work in A. See the "constant angle" constraint for details.
//
// -----------------
// perr = ~RFx * RBy  (with all axes expressed in A)
//        ~RFy * RBz
//        ~RFz * RBx
// -----------------
//
// ---------------------------------
// verr = ~(w_AF-w_AB) * (RFx % RBy)
//      = ~(w_AF-w_AB) * (RFy % RBz)
//      = ~(w_AF-w_AB) * (RFz % RBx)
// ---------------------------------
//
// -----------------------------------------------------------------------
// aerr = ~(b_AF-b_AB) * (RFx % RBy)
//                 + ~(w_AF-w_AB) * ((w_AF%RFx) % RBy) - (w_AB%RBy) % RFx)
//        ~(b_AF-b_AB) * (RFy % RBz)
//                 + ~(w_AF-w_AB) * ((w_AF%RFy) % RBz) - (w_AB%RBz) % RFy)
//        ~(b_AF-b_AB) * (RFz % RBx)
//                 + ~(w_AF-w_AB) * ((w_AF%RFz) % RBx) - (w_AB%RBx) % RFz)
// -----------------------------------------------------------------------
//
// Constraint torque can be determined by inspection of verr:
//    t_F =   lambda_x * (RFx % RBy)   (applied to body F)
//          + lambda_y * (RFy % RBz)
//          + lambda_z * (RFz % RBx)
//    t_B = -t_F                       (applied to body B)
//

// -----------------
// perr = ~RFx * RBy  (with all axes expressed in A)
//        ~RFy * RBz
//        ~RFz * RBx
// -----------------
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size()==3);

    const Rotation& R_AB = getBodyRotation(allX_AB, B);
    const Rotation& R_AF = getBodyRotation(allX_AB, F);
    const Rotation  RB = R_AB * defaultRB; // now expressed in A
    const Rotation  RF = R_AF * defaultRF;

    Vec3::updAs(&perr[0]) = Vec3(~RF.x()*RB.y(),
                                 ~RF.y()*RB.z(),
                                 ~RF.z()*RB.x());
}

// ----------------------------------
// verr = ~(w_AF-w_AB) * (RFx % RBy)
//      = ~(w_AF-w_AB) * (RFy % RBz)
//      = ~(w_AF-w_AB) * (RFz % RBx)
// ----------------------------------
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size()==3);
    //TODO: should be able to get p info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultRB; // now expressed in A
    const Rotation  RF = R_AF * defaultRF;

    const Vec3&     w_AB = getBodyAngularVelocity(allV_AB, B);
    const Vec3&     w_AF = getBodyAngularVelocity(allV_AB, F);
    const Vec3      w_BF = w_AF-w_AB; // in A

    Vec3::updAs(&pverr[0]) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );
}

//------------------------------------------------------------------------
// aerr = ~(b_AF-b_AB) * (RFx % RBy)
//                 + ~(w_AF-w_AB) * ((w_AF%RFx) % RBy) - (w_AB%RBy) % RFx)
//        ~(b_AF-b_AB) * (RFy % RBz)
//                 + ~(w_AF-w_AB) * ((w_AF%RFy) % RBz) - (w_AB%RBz) % RFy)
//        ~(b_AF-b_AB) * (RFz % RBx)
//                 + ~(w_AF-w_AB) * ((w_AF%RFz) % RBx) - (w_AB%RBx) % RFz)
//------------------------------------------------------------------------
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==3);
    //TODO: should be able to get p and v info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultRB; // now expressed in A
    const Rotation  RF = R_AF * defaultRF;

    const Vec3&     w_AB = getBodyAngularVelocityFromState(s, B);
    const Vec3&     w_AF = getBodyAngularVelocityFromState(s, F);
    const Vec3      w_BF = w_AF-w_AB; // in A

    const Vec3&     b_AB = getBodyAngularAcceleration(allA_AB, B);
    const Vec3&     b_AF = getBodyAngularAcceleration(allA_AB, F);
    const Vec3      b_BF = b_AF-b_AB; // in A

    Vec3::updAs(&paerr[0]) = 
        Vec3( dot( b_BF, RF.x() % RB.y() )
                  + dot( w_BF, (w_AF%RF.x()) % RB.y() - (w_AB%RB.y()) % RF.x()),
              dot( b_BF, RF.y() % RB.z() )
                  + dot( w_BF, (w_AF%RF.y()) % RB.z() - (w_AB%RB.z()) % RF.y()),
              dot( b_BF, RF.z() % RB.x() )
                  + dot( w_BF, (w_AF%RF.z()) % RB.x() - (w_AB%RB.x()) % RF.z()));
}

//    t_F =   lambda_x * (RFx % RBy)   (applied to body F)
//          + lambda_y * (RFy % RBz)
//          + lambda_z * (RFz % RBx)
//    t_B = -t_F                       (applied to body B)
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==3 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Vec3& lambda = Vec3::getAs(&multipliers[0]);

    //TODO: should be able to get p info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultRB; // now expressed in A
    const Rotation  RF = R_AF * defaultRF;

    const Vec3 torque_F_A =   lambda[0] * (RF.x() % RB.y())
                            + lambda[1] * (RF.y() % RB.z())
                            + lambda[2] * (RF.z() % RB.x());

    addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
    addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);
}

SimTK_DOWNCAST(ConstantOrientationImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::ConstantOrientation;

ConstrainedBodyIndex    B; // B1 is "base" body
ConstrainedBodyIndex    F; // B2 is "follower" body

Rotation    defaultRB; // fixed to B, expressed in B frame; RB = R_B_RB
Rotation    defaultRF; // fixed to F, expressed in F frame; RF = R_F_RF
};



//==============================================================================
//                               WELD IMPL
//==============================================================================
class Constraint::WeldImpl : public ConstraintImpl {
static Real getDefaultAxisDisplayLength() {return 1;}
static Vec3 getDefaultFrameColor(int which) {
    return which==0 ? Blue : Purple;
}
public:
WeldImpl() 
    : ConstraintImpl(6,0,0), axisDisplayLength(-1), // means "use default axis length"
    frameBColor(-1), frameFColor(-1) // means "use default colors"
{   // default Transforms are identity, i.e. body frames
}
WeldImpl* clone() const override { return new WeldImpl(*this); }

// Draw the two frames.
void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setAxisDisplayLength(Real len) {
    // len == 0 means "don't display"
    // len < 0 means "use default"
    invalidateTopologyCache();
    axisDisplayLength = len >= 0 ? len : -1;
}
Real getAxisDisplayLength() const {
    return axisDisplayLength < 0 ? getDefaultAxisDisplayLength() : axisDisplayLength;
}

void setFrameColor(int which, const Vec3& color) {
    assert(which==0 || which==1);
    // color[0] < 0 means "use default color for this frame"
    invalidateTopologyCache();
    if (which==0) frameBColor = color[0] < 0 ? Vec3(-1) : color;
    else          frameFColor = color[0] < 0 ? Vec3(-1) : color;
}
Vec3 getFrameColor(int which) const {
    assert(which==0 || which==1);
    if (which==0) return frameBColor[0] < 0 ? getDefaultFrameColor(0) : frameBColor;
    else          return frameFColor[0] < 0 ? getDefaultFrameColor(1) : frameFColor;
}

// Implementation of virtuals required for holonomic constraints.

// For theory, look at the ConstantOrientation (1st 3 equations) and 
// Ball (last 3 equations) theory above. Otherwise just lay back and 
// enjoy the ride.

void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const override
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size()==6);

    const Rotation& R_AB = getBodyRotation(allX_AB, B);
    const Rotation& R_AF = getBodyRotation(allX_AB, F);
    const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
    const Rotation  RF = R_AF * defaultFrameF.R();

    // Orientation error
    Vec3::updAs(&perr[0]) = Vec3(~RF.x()*RB.y(),
                                 ~RF.y()*RB.z(),
                                 ~RF.z()*RB.x());

    const Vec3 p_AF1 = findStationLocation(allX_AB, B, defaultFrameB.p());
    const Vec3 p_AF2 = findStationLocation(allX_AB, F, defaultFrameF.p());

    // position error
    Vec3::updAs(&perr[3]) = p_AF2 - p_AF1;
}

void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const override 
{
    assert(allV_AB.size()==2 && constrainedQDot.size()==0 && pverr.size()==6);
    //TODO: should be able to get p info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
    const Rotation  RF = R_AF * defaultFrameF.R();

    const Vec3&     w_AB = getBodyAngularVelocity(allV_AB, B);
    const Vec3&     w_AF = getBodyAngularVelocity(allV_AB, F);
    const Vec3      w_BF = w_AF-w_AB; // in A

    // orientation error
    Vec3::updAs(&pverr[0]) = Vec3( ~w_BF * (RF.x() % RB.y()),
                                   ~w_BF * (RF.y() % RB.z()),
                                   ~w_BF * (RF.z() % RB.x()) );

    //TODO: should be able to get p info from State
    const Transform&  X_AB   = getBodyTransformFromState(s, B);
    const Vec3        p_AF2  = findStationLocationFromState(s, F, 
                                                            defaultFrameF.p());
    const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

    const Vec3        v_AF2  = findStationVelocity(s, allV_AB, F, 
                                                   defaultFrameF.p());
    const Vec3        v_AC   = findStationVelocity(s, allV_AB, B, p_BC);
 
    // position error
    Vec3::updAs(&pverr[3]) = v_AF2 - v_AC;
}

void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const override
{
    assert(allA_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==6);
    //TODO: should be able to get p and v info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
    const Rotation  RF = R_AF * defaultFrameF.R();

    const Vec3&     w_AB = getBodyAngularVelocityFromState(s, B);
    const Vec3&     w_AF = getBodyAngularVelocityFromState(s, F);
    const Vec3      w_BF = w_AF-w_AB; // in A

    const Vec3&     b_AB = getBodyAngularAcceleration(allA_AB, B);
    const Vec3&     b_AF = getBodyAngularAcceleration(allA_AB, F);
    const Vec3      b_BF = b_AF-b_AB; // in A

    // orientation error
    Vec3::updAs(&paerr[0]) = 
        Vec3( dot( b_BF, RF.x() % RB.y() )
                  + dot( w_BF, (w_AF%RF.x()) % RB.y() - (w_AB%RB.y()) % RF.x()),
              dot( b_BF, RF.y() % RB.z() )
                  + dot( w_BF, (w_AF%RF.y()) % RB.z() - (w_AB%RB.z()) % RF.y()),
              dot( b_BF, RF.z() % RB.x() )
                  + dot( w_BF, (w_AF%RF.z()) % RB.x() - (w_AB%RB.x()) % RF.z()));

    const Transform&  X_AB   = getBodyTransformFromState(s, B);
    const Vec3        p_AF2  = findStationLocationFromState(s, F, 
                                                            defaultFrameF.p());
    const Vec3        p_BC   = ~X_AB*p_AF2; // C is a material point of body B

    const Vec3        a_AF2  = findStationAcceleration(s, allA_AB, F, 
                                                       defaultFrameF.p());
    const Vec3        a_AC   = findStationAcceleration(s, allA_AB, B, p_BC);

    // position error
    Vec3::updAs(&paerr[3]) = a_AF2 - a_AC;
}

void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
    assert(multipliers.size()==6 && bodyForcesInA.size()==2 
           && qForces.size()==0);

    const Vec3& torques = Vec3::getAs(&multipliers[0]);
    const Vec3& force_A = Vec3::getAs(&multipliers[3]);

    //TODO: should be able to get p info from State
    const Rotation& R_AB = getBodyRotationFromState(s, B);
    const Rotation& R_AF = getBodyRotationFromState(s, F);
    const Rotation  RB = R_AB * defaultFrameB.R(); // now expressed in A
    const Rotation  RF = R_AF * defaultFrameF.R();

    const Vec3 torque_F_A =   torques[0] * (RF.x() % RB.y())
                            + torques[1] * (RF.y() % RB.z())
                            + torques[2] * (RF.z() % RB.x());

    addInBodyTorque(s, F,  torque_F_A, bodyForcesInA);
    addInBodyTorque(s, B, -torque_F_A, bodyForcesInA);

    const Transform& X_AB  = getBodyTransformFromState(s,B);
    const Vec3&      p_FF2 = defaultFrameF.p();
    const Vec3       p_AF2 = findStationLocationFromState(s, F, p_FF2);
    const Vec3       p_BC = ~X_AB * p_AF2;

    addInStationForce(s, F, p_FF2, force_A, bodyForcesInA);
    addInStationForce(s, B, p_BC, -force_A, bodyForcesInA);
}

SimTK_DOWNCAST(WeldImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::Weld;

ConstrainedBodyIndex    B; // aka "body 1"
ConstrainedBodyIndex    F; // aka "body 2"

Transform               defaultFrameB; // on body 1, relative to B frame
Transform               defaultFrameF; // on body 2, relative to F frame};

// These are for visualization control only.
Real                    axisDisplayLength; // for all 6 axes; <= 0 means "don't
Vec3                    frameBColor;       //       display"
Vec3                    frameFColor;
};



//==============================================================================
//                             NO SLIP 1D IMPL
//==============================================================================
class Constraint::NoSlip1DImpl : public ConstraintImpl {
public:
NoSlip1DImpl()
:   ConstraintImpl(0,1,0), defaultContactPoint(0), defaultNoSlipDirection(),
    directionLength(1), pointRadius(Real(0.05)) 
{ }
NoSlip1DImpl* clone() const override { return new NoSlip1DImpl(*this); }

// The default contact point and no-slip direction may be overridden by 
// setting an instance variable in the state. We allocate the state 
// resources here.
void realizeTopologyVirtual(State& state) const override;

// Return the contact point and no-slip direction, both expressed 
// in the Case body frame C.
const std::pair<Vec3,UnitVec3>& getContactInfo(const State& state) const;

// Return a writable reference into the Instance-stage state variable 
// containing the contact point and no-slip direction, both expressed 
// in the Case body frame C. Calling this
// method invalidates the Instance stage and above in the given state.
std::pair<Vec3,UnitVec3>& updContactInfo(State& state) const;

void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

void setDirectionDisplayLength(Real l) {
    // l <= 0 means don't display direction line
    invalidateTopologyCache();
    directionLength = l > 0 ? l : 0;
}
Real getDirectionDisplayLength() const {return directionLength;}

void setPointDisplayRadius(Real r) {
    // r <= 0 means don't display point
    invalidateTopologyCache();
    pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return pointRadius;}

// Implementation of virtuals required for nonholonomic constraints.

// One non-holonomic constraint equation. There is a contact point P and a 
// no-slip direction n fixed in a case body C. There are two moving bodies B0 
// and B1. The material point P0 of B0 and the material point P1 of B1 which 
// are each coincident with the contact point P must have identical velocities 
// in C, along the direction n. This can be used to implement simple rolling 
// contact between disks, such as occurs in gear trains.
//
// There is no perr equation here since this is a non-holonomic (velocity) 
// constraint. In the C frame, the constraint we want is
//    verr = ~(v_CP1 - v_CP0) * n_C
// that is, the two contact points have no relative velocity in C along the 
// normal. We can calculate this in A instead since the velocities in C of each
// point will differ from their velocities in A by a constant (because they are
// both in the same place in space). So:
//    verr = ~(v_AP1 - v_AP0) * n_A
// Differentiating material point velocities in A, we get the acceleration 
// error
//    aerr = ~(a_AP1 - a_AP0) * n_A + ~(v_AP1 - v_AP0) * (w_AC X n_A)
//         = ~(a_AP1 - a_AP0 - w_AC X (v_AP1 - v_AP0)) * n_A
// 
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const override
{
    // There ought to be at least 2 distinct bodies, up to 3.
    assert((allV_AB.size()==2 || allV_AB.size() == 3)
           && constrainedU.size()==0 && verr.size()==1);

    const std::pair<Vec3,UnitVec3>& info = getContactInfo(s);
    const Vec3&     P_C = info.first;
    const UnitVec3& n_C = info.second;

    //TODO: should be able to get p info from State
    const Transform& X_AC  = getBodyTransformFromState(s, caseBody);
    const Transform& X_AB0 = getBodyTransformFromState(s, movingBody0);
    const Transform& X_AB1 = getBodyTransformFromState(s, movingBody1);
    const Vec3       p_AP  =  X_AC * P_C;       // P's location in A
    const Vec3       p_P0  = ~X_AB0 * p_AP;     // P0's station in B0
    const Vec3       p_P1  = ~X_AB1 * p_AP;     // P1's station in B1
    const UnitVec3   n_A   = X_AC.R() * n_C;

    const Vec3       v_AP0 = findStationVelocity(s, allV_AB, movingBody0, p_P0);
    const Vec3       v_AP1 = findStationVelocity(s, allV_AB, movingBody1, p_P1);

    // Calculate this scalar using A-frame vectors.
    verr[0] = ~(v_AP1-v_AP0) * n_A;
}

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const override
{
    assert((allA_AB.size()==2 || allA_AB.size() == 3)
           && constrainedUDot.size()==0 && vaerr.size()==1);

    const std::pair<Vec3,UnitVec3>& info = getContactInfo(s);
    const Vec3&     P_C = info.first;
    const UnitVec3& n_C = info.second;

    //TODO: should be able to get p and v info from State
    const Transform& X_AC  = getBodyTransformFromState(s, caseBody);
    const Transform& X_AB0 = getBodyTransformFromState(s, movingBody0);
    const Transform& X_AB1 = getBodyTransformFromState(s, movingBody1);
    const Vec3       p_AP  =  X_AC * P_C;       // P's location in A
    const Vec3       p_P0  = ~X_AB0 * p_AP;     // P0's station in B0
    const Vec3       p_P1  = ~X_AB1 * p_AP;     // P1's station in B1
    const UnitVec3   n_A   = X_AC.R() * n_C;

    const Vec3  v_AP0 = findStationVelocityFromState(s, movingBody0, p_P0);
    const Vec3  v_AP1 = findStationVelocityFromState(s, movingBody1, p_P1);
    const Vec3& w_AC  = getBodyAngularVelocityFromState(s, caseBody);

    const Vec3  a_AP0 = findStationAcceleration(s, allA_AB, movingBody0, p_P0);
    const Vec3  a_AP1 = findStationAcceleration(s, allA_AB, movingBody1, p_P1);

    // Calculate this scalar using A-frame vectors.
    vaerr[0] = ~(a_AP1-a_AP0 - w_AC % (v_AP1-v_AP0)) * n_A;
}

// apply f=lambda*n to contact point P1 of B1,
//      -f          to contact point P2 of B2
void addInVelocityConstraintForcesVirtual
   (const State&                                    s,      // Stage::Velocity
    const Array_<Real>&                             multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) const override
{
    assert(multipliers.size()==1 && mobilityForces.size()==0 
           && (bodyForcesInA.size()==2 || bodyForcesInA.size()==3));

    const std::pair<Vec3,UnitVec3>& info = getContactInfo(s);
    const Vec3&     P_C = info.first;
    const UnitVec3& n_C = info.second;

    const Real lambda = multipliers[0];

    //TODO: should be able to get p info from State
    const Transform& X_AC  = getBodyTransformFromState(s, caseBody);
    const Transform& X_AB0 = getBodyTransformFromState(s, movingBody0);
    const Transform& X_AB1 = getBodyTransformFromState(s, movingBody1);
    const Vec3       p_AP  = X_AC * P_C; // P's location in A
    const Vec3       p_P0  = ~X_AB0 * p_AP;              // P0's station in B0
    const Vec3       p_P1  = ~X_AB1 * p_AP;              // P1's station in B1

    const Vec3       force_A = X_AC.R()*(lambda*n_C);

    addInStationForce(s, movingBody1, p_P1,  force_A, bodyForcesInA);
    addInStationForce(s, movingBody0, p_P0, -force_A, bodyForcesInA);
}

SimTK_DOWNCAST(NoSlip1DImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::NoSlip1D;

ConstrainedBodyIndex    caseBody;     // C
ConstrainedBodyIndex    movingBody0;  // B0
ConstrainedBodyIndex    movingBody1;  // B1

Vec3                    defaultContactPoint;      // on body C, exp. in C frame
UnitVec3                defaultNoSlipDirection;   // on body C, exp. in C frame

// These are just for visualization
Real                    directionLength;
Real                    pointRadius;

// This Instance-stage variable holds the actual contact point and no-slip
// direction on the Case body.
mutable DiscreteVariableIndex   contactInfoIx;
};


//==============================================================================
//                        CONSTANT COORDINATE IMPL
//==============================================================================
class Constraint::ConstantCoordinateImpl : public ConstraintImpl {
public:
ConstantCoordinateImpl()
:   ConstraintImpl(1,0,0), theMobilizer(), whichCoordinate(), 
    defaultPosition(NaN){}

ConstantCoordinateImpl* clone() const override 
{   return new ConstantCoordinateImpl(*this); }

// Allocate a state variable to hold the desired position.
void realizeTopologyVirtual(State& state) const override;
// Obtain the currently-set desired position from the state.
Real getPosition(const State& state) const;
// Get a reference to the desired position in the state; this 
// invalidates Position stage in the supplied state.
Real& updPosition(State& state) const;

// Implementation of virtuals required for holonomic constraints.

// One holonomic constraint equation.
//    perr = q - p
//    verr = qdot
//    aerr = qdotdot
// 
void calcPositionErrorsVirtual      
   (const State&                                    s, // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)  // mp of these
    const override
{
    // All the q's for a given constrained mobilizer are considered 
    // constrainedQ's, but we're just going to grab one of them.
    assert(X_AB.size()==0 && constrainedQ.size()>=1 && perr.size()==1);
    const Real q = getOneQ(s, constrainedQ, theMobilizer, whichCoordinate);
    perr[0] = q - getPosition(s);
}

void calcPositionDotErrorsVirtual      
   (const State&                                    s, // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) // mp of these
    const override
{
    // All the q's for a given constrained mobilizer are considered 
    // constrainedQDots's, but we're just going to grab one of them.
    assert(V_AB.size()==0 && constrainedQDot.size()>=1 && pverr.size()==1);
    const Real qdot = getOneQDot(s, constrainedQDot, 
                                 theMobilizer, whichCoordinate);
    pverr[0] = qdot;
}

// Pull t, X_AB, q, V_AB, qdot from state.
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s, // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) // mp of these
    const override
{
    // All the q's for a given constrained mobilizer are considered 
    // constrainedQDotDots's, but we're just going to grab one of them.
    assert(A_AB.size()==0 && constrainedQDotDot.size()>=1 && paerr.size()==1);
    const Real qdotdot = getOneQDotDot(s, constrainedQDotDot, 
                                       theMobilizer, whichCoordinate);
    paerr[0] = qdotdot;
}

// apply generalized force lambda to the mobility
void addInPositionConstraintForcesVirtual
   (const State&                                    s, // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const override
{
   // All the coordinates for a given constrained mobilizer have slots in
   // qForces, but we're just going to update one of them.
   assert(multipliers.size()==1 && bodyForcesInA.size()==0 
          && qForces.size()>=1);

    const Real lambda = multipliers[0];
    addInOneQForce(s, theMobilizer, whichCoordinate, lambda, qForces); 
}

SimTK_DOWNCAST(ConstantCoordinateImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::ConstantCoordinate;

// TOPOLOGY STATE
ConstrainedMobilizerIndex   theMobilizer;
MobilizerQIndex             whichCoordinate;
Real                        defaultPosition;

// TOPOLOGY CACHE
DiscreteVariableIndex       positionIx;
};




//==============================================================================
//                          CONSTANT SPEED IMPL
//==============================================================================
class Constraint::ConstantSpeedImpl : public ConstraintImpl {
public:
ConstantSpeedImpl()
:   ConstraintImpl(0,1,0), theMobilizer(), whichMobility(), 
    defaultSpeed(NaN){}

ConstantSpeedImpl* clone() const override 
{   return new ConstantSpeedImpl(*this); }

// Allocate a state variable to hold the desired speed.
void realizeTopologyVirtual(State& state) const override;
// Obtain the currently-set desired speed from the state.
Real getSpeed(const State& state) const;
// Get a reference to the desired speed in the state; this 
// invalidates Velocity stage in the supplied state.
Real& updSpeed(State& state) const;

// Implementation of virtuals required for nonholonomic constraints.

// One non-holonomic (well, velocity-level) constraint equation.
//    verr = u - s
//    aerr = udot
// 
void calcVelocityErrorsVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allV_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr)   // mv of these
    const override
{
    // All the u's for a given constrained mobilizer are considered 
    // constrainedU's, but we're just going to grab one of them.
    assert(allV_AB.size()==0 && constrainedU.size()>=1 && verr.size()==1);
    const Real u = getOneU(s, constrainedU, theMobilizer, whichMobility);
    verr[0] = u - getSpeed(s);
}

void calcVelocityDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr)  // mv of these
    const override
{
    // All the u's for a given constrained mobilizer are considered 
    // constrainedUDots's, but we're just going to grab one of them.
    assert(allA_AB.size()==0 && constrainedUDot.size()>=1 && vaerr.size()==1);
    const Real udot = getOneUDot(s, constrainedUDot, 
                                 theMobilizer, whichMobility);
    vaerr[0] = udot;
}

// apply generalized force lambda to the mobility
void addInVelocityConstraintForcesVirtual
   (const State&                              s,      // Stage::Velocity
    const Array_<Real>&                       multipliers, // mv of these
    Array_<SpatialVec,ConstrainedBodyIndex>&  bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&     mobilityForces) const override
{
   // All the mobilities for a given constrained mobilizer have slots in
   // mobilizedForces, but we're just going to update one of them.
   assert(multipliers.size()==1 && bodyForcesInA.size()==0
           && mobilityForces.size()>=1);

    const Real lambda = multipliers[0];
    addInOneMobilityForce(s, theMobilizer, whichMobility, lambda, 
                          mobilityForces); 
}

SimTK_DOWNCAST(ConstantSpeedImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::ConstantSpeed;

// TOPOLOGY STATE
ConstrainedMobilizerIndex   theMobilizer;
MobilizerUIndex             whichMobility;
Real                        defaultSpeed;

// TOPOLOGY CACHE
DiscreteVariableIndex       speedIx;
};



//==============================================================================
//                        CONSTANT ACCELERATION IMPL
//==============================================================================
class Constraint::ConstantAccelerationImpl : public ConstraintImpl {
public:
ConstantAccelerationImpl()
:   ConstraintImpl(0,0,1), theMobilizer(), whichMobility(), 
    defaultAcceleration(NaN) {}

ConstantAccelerationImpl* clone() const override 
{   return new ConstantAccelerationImpl(*this); }

// Allocate a state variable to hold the desired acceleration.
void realizeTopologyVirtual(State& state) const override;
// Obtain the currently-set desired acceleration from the state.
Real getAcceleration(const State& state) const;
// Get a reference to the desired acceleration in the state; this 
// invalidates Acceleration stage in the supplied state.
Real& updAcceleration(State& state) const;

// Implementation of virtuals required for acceleration-only constraints.

// One acceleration-only constraint equation.
//    aerr = udot - a
// 
void calcAccelerationErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  allA_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr) const override // ma of these
{
    // All the u's for a given constrained mobilizer are considered 
    // constrainedUDots's, but we're just going to grab one of them.
    assert(allA_AB.size()==0 && constrainedUDot.size()>=1 && aerr.size()==1);

    const Real udot = getOneUDot(s, constrainedUDot, 
                                 theMobilizer, whichMobility);
    aerr[0] = udot - getAcceleration(s);
}

// apply generalized force lambda to the mobility
void addInAccelerationConstraintForcesVirtual
   (const State&                                    state, // Stage::Velocity
    const Array_<Real>&                             multipliers, // ma of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedUIndex>&           mobilityForces) const override
{
   // All the mobilities for a given constrained mobilizer have slots in
   // mobilizedForces, but we're just going to update one of them.
    assert(multipliers.size()==1 && bodyForcesInA.size()==0
           && mobilityForces.size()>=1);

    const Real lambda = multipliers[0];
    addInOneMobilityForce(state, theMobilizer, whichMobility, lambda, 
                          mobilityForces); 
}

SimTK_DOWNCAST(ConstantAccelerationImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::ConstantAcceleration;

// TOPOLOGY STATE
ConstrainedMobilizerIndex   theMobilizer;
MobilizerUIndex             whichMobility;
Real                        defaultAcceleration;

// TOPOLOGY CACHE
DiscreteVariableIndex       accelIx;
};



//==============================================================================
//                         CUSTOM IMPLEMENTATION IMPL
//==============================================================================
// This class exists primarily to allow the Custom::Implementation class to keep
// a pointer to its handle class's CustomImpl class which is derived from 
// ConstraintImpl which has all the goodies that are needed for defining a 
// Constraint.
//
// At first this class is the owner of the CustomImpl. Then when this is put in
// a Custom handle, that handle takes over ownership of the CustomImpl and the 
// CustomImpl takes over ownership of this ImplementationImpl object.
class Constraint::Custom::ImplementationImpl 
:   public PIMPLImplementation<Implementation, ImplementationImpl>
{
public:
// no default constructor
explicit ImplementationImpl(CustomImpl* customImpl) : isOwner(true), builtInImpl(customImpl) { }
inline ~ImplementationImpl(); // see below -- have to wait for CustomImpl's definition

// Copying one of these just gives us a new one with a NULL CustomImpl pointer.
ImplementationImpl(const ImplementationImpl& src) : isOwner(false), builtInImpl(0) { }

ImplementationImpl* clone() const {return new ImplementationImpl(*this);}

bool isOwnerOfCustomImpl() const {return builtInImpl && isOwner;}
CustomImpl* removeOwnershipOfCustomImpl() {
    assert(isOwnerOfCustomImpl()); 
    isOwner=false; 
    return builtInImpl;
}

void setReferenceToCustomImpl(CustomImpl* cimpl) {
    assert(!builtInImpl); // you can only do this once
    isOwner=false;
    builtInImpl = cimpl;
}

bool hasCustomImpl() const {return builtInImpl != 0;}

const CustomImpl& getCustomImpl() const {
    assert(builtInImpl);
    return *builtInImpl;
}
CustomImpl& updCustomImpl() {
    assert(builtInImpl);
    return *builtInImpl;
}

//------------------------------------------------------------------------------
                                    private:
bool            isOwner;
CustomImpl*     builtInImpl; // just a reference; not owned

// suppress assignment
ImplementationImpl& operator=(const ImplementationImpl&);
};



//==============================================================================
//                                 CUSTOM IMPL
//==============================================================================
class Constraint::CustomImpl : public ConstraintImpl {
public:
CustomImpl() : implementation(0) { }
CustomImpl(int mp, int mv, int ma) : ConstraintImpl(mp,mv,ma), implementation(0) { }

void takeOwnershipOfImplementation(Custom::Implementation* userImpl);

explicit CustomImpl(Custom::Implementation* userImpl) : implementation(0) { 
    assert(userImpl);
    implementation = userImpl;
    implementation->updImpl().setReferenceToCustomImpl(this);
}    

// Copy constructor
CustomImpl(const CustomImpl& src) : implementation(0) {
    if (src.implementation) {
        implementation = src.implementation->clone();
        implementation->updImpl().setReferenceToCustomImpl(this);
    }
}
    
~CustomImpl() {
    delete implementation;
}
    
CustomImpl* clone() const override { return new CustomImpl(*this); }

const Custom::Implementation& getImplementation() const {
    assert(implementation);
    return *implementation;
}

Custom::Implementation& updImplementation() {
    assert(implementation);
    return *implementation;
}

// Forward all the virtuals to the Custom::Implementation virtuals.
void realizeTopologyVirtual(State& s) const override {getImplementation().realizeTopology(s);}
void realizeModelVirtual   (State& s) const override {getImplementation().realizeModel(s);}
void realizeInstanceVirtual(const State& s) const override {getImplementation().realizeInstance(s);}
void realizeTimeVirtual    (const State& s) const override {getImplementation().realizeTime(s);}
void realizePositionVirtual(const State& s) const override {getImplementation().realizePosition(s);}
void realizeVelocityVirtual(const State& s) const override {getImplementation().realizeVelocity(s);}
void realizeDynamicsVirtual(const State& s) const override {getImplementation().realizeDynamics(s);}
void realizeAccelerationVirtual(const State& s) const override {getImplementation().realizeAcceleration(s);}
void realizeReportVirtual  (const State& s) const override {getImplementation().realizeReport(s);}

void calcPositionErrorsVirtual     
   (const State&                                    state,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const override
{   getImplementation().calcPositionErrors(state,X_AB,constrainedQ,perr); }

void calcPositionDotErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const override
{   getImplementation().calcPositionDotErrors
                                        (state,V_AB,constrainedQDot,pverr); }

void calcPositionDotDotErrorsVirtual     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const override
{   getImplementation().calcPositionDotDotErrors
                                    (state,A_AB,constrainedQDotDot,paerr); }

void addInPositionConstraintForcesVirtual
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedQIndex>&             qForces) const override
{   getImplementation().addInPositionConstraintForces
        (s,multipliers,bodyForcesInA,qForces); }

void calcVelocityErrorsVirtual     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) const override
{   getImplementation().calcVelocityErrors(state,V_AB,constrainedU,verr); }

void calcVelocityDotErrorsVirtual     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) const override
{   getImplementation().calcVelocityDotErrors
                                        (state,A_AB,constrainedUDot,vaerr); }

void addInVelocityConstraintForcesVirtual
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const override
{   getImplementation().addInVelocityConstraintForces
        (s,multipliers,bodyForcesInA,mobilityForces); }

void calcAccelerationErrorsVirtual      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   aerr) const override
{   getImplementation().calcVelocityDotErrors
                                        (state,A_AB,constrainedUDot,aerr); }

void addInAccelerationConstraintForcesVirtual
   (const State&                                s, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForcesInA,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const override
{   getImplementation().addInAccelerationConstraintForces
        (s,multipliers,bodyForcesInA,mobilityForces); }

void calcDecorativeGeometryAndAppendVirtual
        (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override
    {getImplementation().calcDecorativeGeometryAndAppend(s,stage,geom);}

SimTK_DOWNCAST(CustomImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::Custom;

Custom::Implementation*     implementation;

CustomImpl& operator=(const CustomImpl&); // suppress assignment
};

// Need definition for CustomImpl here in case we have to delete it.
inline Constraint::Custom::ImplementationImpl::~ImplementationImpl() {
    if (isOwner) 
        delete builtInImpl; 
    builtInImpl=0;
}



//==============================================================================
//                          COORDINATE COUPLER IMPL
//==============================================================================
class Constraint::CoordinateCouplerImpl 
:   public Constraint::Custom::Implementation {
public:
CoordinateCouplerImpl(SimbodyMatterSubsystem&           matter, 
                      const Function*                   function, 
                      const Array_<MobilizedBodyIndex>& coordBody, 
                      const Array_<MobilizerQIndex>&    coordIndex);
    
~CoordinateCouplerImpl() {
    if (--referenceCount[0] == 0) {
        delete function;
        delete[] referenceCount;
    }
}
    
Implementation* clone() const override {
    referenceCount[0]++;
    CoordinateCouplerImpl* newCoupler = new CoordinateCouplerImpl(*this);
    return newCoupler;
}

void calcPositionErrors     
   (const State&                                    state,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const override;

void calcPositionDotErrors      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const override;

void calcPositionDotDotErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const override;

void addInPositionConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedQIndex>&             qForces) const override;

//------------------------------------------------------------------------------
                                    private:
friend class Constraint::CoordinateCoupler;

//  TOPOLOGY STATE
const Function*                     function;
Array_<ConstrainedMobilizerIndex>   coordBodies;
Array_<MobilizerQIndex>             coordIndices;

//  TOPOLOGY CACHE
//  None.

//  A reusable temporary variable allocated to the correct size
//  to hold all the Function arguments.
mutable Vector                      temp;

// This allows copies to be made of this constraint which share
// the function object.
int*                                referenceCount;

};



//==============================================================================
//                           SPEED COUPLER IMPL
//==============================================================================
class Constraint::SpeedCouplerImpl 
:   public Constraint::Custom::Implementation {
public:
SpeedCouplerImpl(SimbodyMatterSubsystem& matter, 
                 const Function*                    function, 
                 const Array_<MobilizedBodyIndex>&  speedBody, 
                 const Array_<MobilizerUIndex>&     speedIndex,
                 const Array_<MobilizedBodyIndex>&  coordBody, 
                 const Array_<MobilizerQIndex>&     coordIndex);
    
~SpeedCouplerImpl() {
    if (--referenceCount[0] == 0) {
        delete function;
        delete[] referenceCount;
    }
}
    
Implementation* clone() const override {
    referenceCount[0]++;
    return new SpeedCouplerImpl(*this);
}

void calcVelocityErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedU,
    Array_<Real>&                                   verr) const override;

void calcVelocityDotErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedUIndex>&     constrainedUDot,
    Array_<Real>&                                   vaerr) const override;

void addInVelocityConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedUIndex>&             mobilityForces) const override;

//------------------------------------------------------------------------------
                                    private:

const Function*                     function;
int*                                referenceCount;
Array_<ConstrainedMobilizerIndex>   speedBodies;
Array_<MobilizedBodyIndex>          coordBodies;
Array_<MobilizerUIndex>             speedIndices;
Array_<MobilizerQIndex>             coordIndices;
mutable Vector                      temp;
};



//==============================================================================
//                          PRESCRIBED MOTION IMPL
//==============================================================================
class Constraint::PrescribedMotionImpl 
:   public Constraint::Custom::Implementation {
public:
PrescribedMotionImpl(SimbodyMatterSubsystem&    matter, 
                     const Function*            function, 
                     MobilizedBodyIndex         coordBody, 
                     MobilizerQIndex            coordIndex);
    
~PrescribedMotionImpl() {
    if (--referenceCount[0] == 0) {
        delete function;
        delete[] referenceCount;
    }
}
    
Implementation* clone() const override {
    referenceCount[0]++;
    return new PrescribedMotionImpl(*this);
}

void calcPositionErrors     
   (const State&                                    state,
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr) const override;

void calcPositionDotErrors      
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr) const override;

void calcPositionDotDotErrors     
   (const State&                                    state,
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr) const override;

void addInPositionConstraintForces
   (const State&                                state, 
    const Array_<Real>&                         multipliers,
    Array_<SpatialVec,ConstrainedBodyIndex>&    bodyForces,
    Array_<Real,ConstrainedQIndex>&             qForces) const override;

//------------------------------------------------------------------------------
                                    private:
const Function*             function;
int*                        referenceCount;
ConstrainedMobilizerIndex   coordBody;
MobilizerQIndex             coordIndex;
mutable Vector              temp;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_IMPL_H_



