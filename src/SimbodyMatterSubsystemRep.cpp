/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Derived from NIH IVM code written by Charles Schwieters      *
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
 * Implementation of SimbodyMatterSubsystemRep.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"

#include "SimbodyMatterSubsystemRep.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "LengthConstraints.h"
#include "MultibodySystemRep.h"
#include "MobilizedBodyImpl.h"
#include "ConstraintRep.h"

#include <string>
#include <iostream>
using std::cout; using std::endl;

//#define USE_OLD_CONSTRAINTS

SimbodyMatterSubsystemRep::SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep& src)
  : SimTK::Subsystem::Guts("SimbodyMatterSubsystemRep", "X.X.X")
{
    assert(!"SimbodyMatterSubsystemRep copy constructor ... TODO!");
}


void SimbodyMatterSubsystemRep::clearTopologyState() {
    // Constraints are independent from one another, so any deletion order
    // is fine. However, they depend on bodies and not vice versa so we'll
    // delete them first just to be neat.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        delete constraints[cx];
    constraints.clear();

    // These are the owner handles, so this deletes the MobilizedBodyImpls also.
    // We'll delete from the terminal nodes inward just to be neat.
    for (MobilizedBodyIndex bx(mobilizedBodies.size()-1); bx >= 0; --bx)
        delete mobilizedBodies[bx];
    mobilizedBodies.clear();
}

void SimbodyMatterSubsystemRep::clearTopologyCache() {
    invalidateSubsystemTopologyCache();

    // TODO: state indices really shouldn't be dealt out until Stage::Model.
    // At the moment they are part of the topology.
    nextUSlot   = UIndex(0);
    nextUSqSlot = USquaredIndex(0);
    nextQSlot   = QIndex(0);

    DOFTotal=SqDOFTotal=maxNQTotal         = -1;
    topologyCache.clear();
    topologyCacheIndex = -1;

    delete lConstraints; lConstraints=0;

    for (int i=0; i<(int)distanceConstraints.size(); ++i)
        delete distanceConstraints[i];
    distanceConstraints.clear();

    // RigidBodyNodes themselves are owned by the MobilizedBodyImpls and will
    // be deleted when the MobilizedBodyImpl objects are.
    rbNodeLevels.clear();
    nodeNum2NodeMap.clear();
}

MobilizedBodyIndex SimbodyMatterSubsystemRep::adoptMobilizedBody
   (MobilizedBodyIndex parentIx, MobilizedBody& child) 
{
    invalidateSubsystemTopologyCache();

    const MobilizedBodyIndex ix(mobilizedBodies.size());
    assert(parentIx < ix);

    mobilizedBodies.push_back(new MobilizedBody()); // grow
    MobilizedBody& m = *mobilizedBodies.back(); // refer to the empty handle we just created

    child.disown(m); // transfer ownership to m

    // Now tell the MobilizedBody object its owning MatterSubsystem, id within
    // that Subsystem, and parent MobilizedBody object.
    m.updImpl().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), parentIx, ix);
    return ix;
}

ConstraintIndex SimbodyMatterSubsystemRep::adoptConstraint(Constraint& child) {
    invalidateSubsystemTopologyCache();

    const ConstraintIndex ix(constraints.size());

    constraints.push_back(new Constraint()); // grow
    Constraint& c = *constraints.back(); // refer to the empty handle we just created

    child.disown(c); // transfer ownership to c

    // Now tell the Constraint object its owning MatterSubsystem and id within
    // that Subsystem.
    c.updImpl().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), ix);
    return ix;
}


void SimbodyMatterSubsystemRep::createGroundBody() {
    assert(mobilizedBodies.empty());
    invalidateSubsystemTopologyCache();

    mobilizedBodies.push_back(new MobilizedBody::Ground());
    mobilizedBodies[0]->updImpl().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), 
                                                      MobilizedBodyIndex(), 
                                                      MobilizedBodyIndex(0));
}

// Add a distance constraint and assign it to use a particular slot in the
// qErr, uErr, and multiplier arrays.
// Return the assigned distance constraint index for caller's use.
// TODO: OBSOLETE
int SimbodyMatterSubsystemRep::addOneDistanceConstraintEquation(
    const RBStation& s1, const RBStation& s2, const Real& d)
{
    const int nxtIndex = (int)distanceConstraints.size();
    RBDistanceConstraint* dc = new RBDistanceConstraint(s1,s2,d);
    dc->setQErrIndex(nxtIndex);
    dc->setUErrIndex(nxtIndex);
    dc->setMultIndex(nxtIndex);
    dc->setDistanceConstraintNum(distanceConstraints.size());
    distanceConstraints.push_back(dc);
    return nxtIndex;
}

MobilizedBodyIndex SimbodyMatterSubsystemRep::getParent(MobilizedBodyIndex body) const { 
    return getRigidBodyNode(body).getParent()->getNodeNum();
}

std::vector<MobilizedBodyIndex> SimbodyMatterSubsystemRep::getChildren(MobilizedBodyIndex body) const {
    const RigidBodyNode& node = getRigidBodyNode(body);
    std::vector<MobilizedBodyIndex> children;
    for (MobilizedBodyIndex bx(0); bx < node.getNChildren(); ++bx)
        children.push_back(node.getChild(bx)->getNodeNum());
    return children;
}

const MassProperties&
SimbodyMatterSubsystemRep::getDefaultBodyMassProperties(MobilizedBodyIndex body) const
  { return getRigidBodyNode(body).getMassProperties_OB_B(); }

const Transform&
SimbodyMatterSubsystemRep::getDefaultMobilizerFrame(MobilizedBodyIndex body) const
  { return getRigidBodyNode(body).getX_BM(); }

const Transform&
SimbodyMatterSubsystemRep::getDefaultMobilizerFrameOnParent(MobilizedBodyIndex body) const
  { return getRigidBodyNode(body).getX_PF(); }

const Transform&
SimbodyMatterSubsystemRep::getBodyTransform(const State& s, const SBPositionCache& pc, MobilizedBodyIndex body) const
  { return getRigidBodyNode(body).getX_GB(pc); }

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyVelocity(const State& s, const SBVelocityCache& vc, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getV_GB(vc);
}

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyAcceleration(const State& s, const SBAccelerationCache& ac, MobilizedBodyIndex body) const {
    return getRigidBodyNode(body).getA_GB(ac);
}

const SpatialVec&
SimbodyMatterSubsystemRep::getCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getTotalCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getGyroscopicForce(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getGyroscopicForce(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getCentrifugalForces(getDynamicsCache(s));
}

const SpatialMat&
SimbodyMatterSubsystemRep::getArticulatedBodyInertia(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getArticulatedBodyInertia(getDynamicsCache(s));
}

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes we're going to need later for state
// variables and cache entries. We allocate and initialize all the
// Modeling variables here.
void SimbodyMatterSubsystemRep::endConstruction(State& s) {
    if (subsystemTopologyHasBeenRealized()) 
        return; // already done

    // This creates a RigidBodyNode owned by the the Topology cache of each MobilizedBody.
    // Each RigidBodyNode lists as its parent the RigidBodyNode contained in the 
    // MobilizedBody's parent.
    //
    // We simultaneously build up the computational version of the multibody tree,
    // based on RigidBodyNode objects rather than on MobilizedBody objects.
    nodeNum2NodeMap.clear();
    rbNodeLevels.clear();
    DOFTotal = SqDOFTotal = maxNQTotal = 0;

    // state allocation
    nextUSlot   = UIndex(0);
    nextUSqSlot = USquaredIndex(0);
    nextQSlot   = QIndex(0);

    //Must do these in order from lowest number (ground) to highest. 
    for (int i=0; i<getNumMobilizedBodies(); ++i) {
        // Create the RigidBodyNode properly linked to its parent.
        const MobilizedBodyImpl& mbr = getMobilizedBody(MobilizedBodyIndex(i)).getImpl();
        const RigidBodyNode& n = mbr.realizeTopology(nextUSlot,nextUSqSlot,nextQSlot);

        // Create the computational multibody tree data structures, organized by level.
        const int level = n.getLevel();
        if ((int)rbNodeLevels.size() <= level)
            rbNodeLevels.resize(level+1); // make room for the new level
        const int nodeIndexWithinLevel = rbNodeLevels[level].size();
        rbNodeLevels[level].push_back(&n);
        nodeNum2NodeMap.push_back(RigidBodyNodeIndex(level, nodeIndexWithinLevel));

        // Count up multibody tree totals.
        const int ndof = n.getDOF();
        DOFTotal += ndof; SqDOFTotal += ndof*ndof;
        maxNQTotal += n.getMaxNQ();
    }

    // Order doesn't matter for constraints as long as the bodies are already there.
    // Quaternion normalization constraints exist only at the 
    // position level, however they are not topological since modeling
    // choices affect whether we use them. See realizeModel() below.

//    groundAncestorConstraint.clear();
//    branches.clear(); 
//    if (rbNodeLevels.size() > 1)
//        branches.resize(rbNodeLevels[1].size()); // each level 1 body is a branch

    for (ConstraintIndex cx(0); cx<getNumConstraints(); ++cx) {
        const ConstraintRep& crep = getConstraint(cx).getImpl();
        crep.realizeTopology(s);

        // Create computational constraint data structure. This is organized by
        // the ancestor body's "branch" (ground or level 1 base body), then by 
        // ancestor level along that branch.
        /*
        const int level = crep.getAncestorLevel();
        if (crep.getAncestorLevel() == 0)
            groundAncestorConstraints.push_back(i);
        else {
            const int branch = crep.getAncestorBranch();
            if (branches[branch].size() <= level)
                branches[branch].resize(level+1);
            branches[branch][level].push_back(i);
        }
        */
    }
    
    // Now create the computational data structure for the length constraint
    // equations which currently are used to implement all the Constraints.
    // This is owned by the SimbodyMatterSubsystemRep.
    // TODO: OBSOLETE -- part of the old IVM constraint system
#ifdef USE_OLD_CONSTRAINTS
    lConstraints = new LengthConstraints(*this, 0);
    lConstraints->construct(distanceConstraints);
#endif
}

int SimbodyMatterSubsystemRep::realizeSubsystemTopologyImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "SimbodyMatterSubsystem::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this call*.
    // Afterwards they are truly const so we don't declare them mutable, but cheat
    // here instead.
    SimbodyMatterSubsystemRep* mutableThis = const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!subsystemTopologyHasBeenRealized()) mutableThis->endConstruction(s); // no more bodies after this!

    // Fill in the local copy of the topologyCache from the information
    // calculated in endConstruction(). Also ask the State for some room to
    // put Modeling variables & cache and remember the indices in our construction
    // cache.

    mutableThis->topologyCache.nBodies      = nodeNum2NodeMap.size();
    mutableThis->topologyCache.nParticles   = 0; // TODO
    mutableThis->topologyCache.nConstraints = constraints.size();
    mutableThis->topologyCache.nDOFs        = DOFTotal;
    mutableThis->topologyCache.maxNQs       = maxNQTotal;
    mutableThis->topologyCache.sumSqDOFs    = SqDOFTotal;
    mutableThis->topologyCache.nDistanceConstraints     = distanceConstraints.size();

    SBModelVars mvars;
    mvars.allocate(topologyCache);
    setDefaultModelValues(topologyCache, mvars);
    mutableThis->topologyCache.modelingVarsIndex  = 
        allocateDiscreteVariable(s,Stage::Model, new Value<SBModelVars>(mvars));

    mutableThis->topologyCache.modelingCacheIndex = 
        allocateCacheEntry(s,Stage::Model, new Value<SBModelCache>());

    mutableThis->topologyCache.valid = true;

    // Allocate a cache entry for the topologyCache, and save a copy there.
    mutableThis->topologyCacheIndex = 
        allocateCacheEntry(s,Stage::Topology, new Value<SBTopologyCache>(topologyCache));

    return 0;
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
int SimbodyMatterSubsystemRep::realizeSubsystemModelImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Topology, 
        "SimbodyMatterSubsystem::realizeModel()");

    const SBModelVars& mv = getModelVars(s);

    // Get the Model-stage cache and make sure it has been allocated and initialized if needed.
    // It is OK to hold a reference here because the discrete variables (and cache entries) in
    // the State are stable, that is, they don't change location even if more variables are added.
    SBModelCache& mc = updModelCache(s);
    mc.clear(); // forget any previous modeling information
    mc.allocate(topologyCache);


        // MOBILIZED BODY MODELING

    // Count quaternions, and assign a "quaternion pool" index to each MobilizedBody that
    // needs one, and do the same for the "angle pool". We can't do this until Model stage
    // because it is a Model stage variable which decides whether ball-like joints get
    // quaternions or Euler angles.
    mc.totalNQInUse = mc.totalNUInUse = mc.totalNQuaternionsInUse = mc.totalNAnglesInUse = 0;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            SBModelCache::PerMobilizedBodyModelInfo& 
                mbInfo = mc.updMobilizedBodyModelInfo(node.getNodeNum());

            // Assign q's.
            mbInfo.nQInUse = node.getNQInUse(mv);
            mbInfo.firstQIndex = QIndex(mc.totalNQInUse);
            mc.totalNQInUse += mbInfo.nQInUse;

            // Assign u's.
            mbInfo.nUInUse = node.getNUInUse(mv);
            mbInfo.firstUIndex = UIndex(mc.totalNUInUse);
            mc.totalNUInUse += mbInfo.nUInUse;

            // Assign quaternion pool slot.
            if (node.isUsingQuaternion(mv, mbInfo.startOfQuaternion)) {
                mbInfo.hasQuaternionInUse = true;
                mbInfo.quaternionPoolIndex = QuaternionPoolIndex(mc.totalNQuaternionsInUse);
                mc.totalNQuaternionsInUse++;
            }

            // Assign angle pool slots.
            if (node.isUsingAngles(mv, mbInfo.startOfAngles, mbInfo.nAnglesInUse)) {
                mbInfo.anglePoolIndex = AnglePoolIndex(mc.totalNAnglesInUse);
                mc.totalNAnglesInUse += mbInfo.nAnglesInUse;
            }
        }


    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModel(mv,mc); 


        // CONSTRAINT MODELING


    // Count position, velocity, and acceleration constraint equations generated by
    // each Constraint that has not been disabled. The State's QErr, UErr, UDotErr/Multiplier arrays
    // are laid out like this:
    //
    //            ---------------------- ---------------
    //    QErr   |      mHolonomic      |  mQuaternion  |
    //            ---------------------- ---------------
    //            ---------------------- -------------------
    //    UErr   |          "           |   mNonholonomic   |
    //            ---------------------- -------------------
    //            ---------------------- ------------------- ---------------------
    // UDotErr   |          "           |         "         |   mAccelerationOnly |
    //            ---------------------- ------------------- ---------------------
    //
    // Multipliers are allocated exactly as for UDotErr.
    //
    // Note that a Constraint with both holonomic and nonholonomic constraint equations
    // will get two disjoint segments in UErr (and UDotErr).
    //

    // Each Constraint's realizeModel() method will add its contribution to these.
    mc.totalNHolonomicConstraintEquationsInUse = mc.totalNNonholonomicConstraintEquationsInUse = 
        mc.totalNAccelerationOnlyConstraintEquationsInUse = 0;

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeModel(s);

    //cout << "-------------------------------\n";
    //cout << mc;
    //cout << "-------------------------------\n";

    // Build sets of kinematically coupled constraints. Kinematic coupling can be different at
    // position, velocity, and acceleration levels, with only holonomic constraints included
    // at the position level, holonomic+nonholonic at the velocity level, and
    // holonomic+nonholonomic+accelerationOnly coupled at the acceleration level.

    /*
    positionCoupledConstraints.clear();
    velocityCoupledConstraints.clear();
    accelerationCoupledConstraints.clear();

    for (int i=0; i < (int)constraints.size(); ++i) {
        if (!mc.mHolonomicEquationsInUse[i]) continue;

        positionCoupledConstraints.push_back(CoupledConstraintSet(*constraints[i]));
        CoupledConstraintSet& cset = positionCoupledConstraints.back();

    }
    */


    // quaternion errors are located after last holonomic constraint error; see diagram above
    mc.firstQuaternionQErrSlot = mc.totalNHolonomicConstraintEquationsInUse; 

    // Now allocate all remaining variables and cache entries. We can properly initialize only
    // the next stage, the parameters. Everything else gets some kind of meaningless initial
    // values but those could change as we realize the higher stages.
    SBInstanceVars iv;
    iv.allocate(topologyCache);
    setDefaultInstanceValues(mv, iv);

    mc.instanceVarsIndex = 
        allocateDiscreteVariable(s, Stage::Instance, new Value<SBInstanceVars>(iv));
    mc.instanceCacheIndex = 
        allocateCacheEntry(s, Stage::Instance, new Value<SBInstanceCache>());

    // No time vars or cache
    mc.timeVarsIndex = -1;
    mc.timeCacheIndex = -1;

    // Position variables are just q's, which the State knows how to deal with. 

    Vector qInit(maxNQTotal);
    setDefaultPositionValues(mv, qInit);

    // MobilizedBodies provide default values for their q's. The number and
    // values of these can depend on modeling variables, which are already
    // set in the State. Don't do this for Ground, which has no q's.
    for (MobilizedBodyIndex bx(1); bx < mobilizedBodies.size(); ++bx) {
        const MobilizedBodyImpl& mb = mobilizedBodies[bx]->getImpl();
        mb.copyOutDefaultQ(s, qInit);
    }

    mc.qIndex = s.allocateQ(getMySubsystemIndex(), qInit);
    mc.qVarsIndex = -1; // no config vars other than q
    mc.qCacheIndex = allocateCacheEntry(s, Stage::Position, new Value<SBPositionCache>());

    // We'll store the the physical constraint errors (which consist solely of distance
    // constraint equations at the moment), followed by the quaternion constraints.
    mc.qErrIndex = allocateQErr(s, mc.totalNHolonomicConstraintEquationsInUse + mc.totalNQuaternionsInUse);

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.

    Vector uInit(DOFTotal);
    setDefaultVelocityValues(mv, uInit);

    mc.uIndex = s.allocateU(getMySubsystemIndex(), uInit);
    mc.uVarsIndex = -1; // no velocity vars other than u
    mc.uCacheIndex = allocateCacheEntry(s, Stage::Velocity, new Value<SBVelocityCache>());
    // Note that qdots are automatically allocated in the Velocity stage cache.

    // Only physical constraints exist at the velocity and acceleration levels; 
    // the quaternion normalization constraints are gone.
    mc.uErrIndex    = allocateUErr   (s,   mc.totalNHolonomicConstraintEquationsInUse
                                         + mc.totalNNonholonomicConstraintEquationsInUse);
    mc.udotErrIndex = allocateUDotErr(s,   mc.totalNHolonomicConstraintEquationsInUse
                                         + mc.totalNNonholonomicConstraintEquationsInUse
                                         + mc.totalNAccelerationOnlyConstraintEquationsInUse);

    // no z's
    // We do have dynamic vars for now for forces & pres. accel. but those will
    // probably go away. TODO
    SBDynamicsVars dvars;
    dvars.allocate(topologyCache);
    setDefaultDynamicsValues(mv, dvars);
    mc.dynamicsVarsIndex = 
        allocateDiscreteVariable(s, Stage::Dynamics, new Value<SBDynamicsVars>(dvars));
    mc.dynamicsCacheIndex = 
        allocateCacheEntry(s, Stage::Dynamics, new Value<SBDynamicsCache>());

    // No reaction variables that I know of. But we can go through the
    // charade here anyway.
    SBAccelerationVars rvars;
    rvars.allocate(topologyCache);
    setDefaultAccelerationValues(mv, rvars);

    mc.accelerationVarsIndex = 
        allocateDiscreteVariable(s, Stage::Acceleration, new Value<SBAccelerationVars>(rvars));
    mc.accelerationCacheIndex = 
        allocateCacheEntry(s, Stage::Acceleration, new Value<SBAccelerationCache>());

    // Note that qdots, qdotdots, udots, zdots are automatically allocated by
    // the State when we advance the stage past modeling.

    return 0;
}

// Here we lock in parameterization of the model, such as body masses.
int SimbodyMatterSubsystemRep::realizeSubsystemInstanceImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "SimbodyMatterSubsystem::realizeInstance()");

    const SBModelVars&    mv = getModelVars(s);
    const SBInstanceVars& iv = getInstanceVars(s);

    // Get the Instance-stage cache and make sure it has been allocated and initialized if needed.
    SBInstanceCache& ic = updInstanceCache(s);
    ic.allocate(topologyCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeInstance(mv,iv,ic); 

    ic.totalMass = iv.particleMasses.sum();
    for (int i=0; i<getNBodies(); ++i)
        ic.totalMass += iv.bodyMassProperties[i].getMass();

    return 0;
}

int SimbodyMatterSubsystemRep::realizeSubsystemTimeImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "SimbodyMatterSubsystem::realizeTime()");

    // nothing yet

    return 0;
}

// Set generalized coordinates: sweep from base to tips.
int SimbodyMatterSubsystemRep::realizeSubsystemPositionImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "SimbodyMatterSubsystem::realizePosition()");

    // Set up StateDigest for calculating position information.
    SBStateDigest stateDigest(s, *this, Stage::Position);

    const SBModelVars&  mv   = getModelVars(s);
    const SBModelCache& mc   = getModelCache(s);
    const Vector&       q    = getQ(s);
    Vector&             qErr = updQErr(s);

    // Get the Position-stage cache and make sure it has been allocated and initialized if needed.
    SBPositionCache& pc = updPositionCache(s);
    pc.allocate(topologyCache);

    // Any body which is using quaternions should calculate the quaternion
    // constraint here and put it in the appropriate slot of qErr.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizePosition(mv,mc,q,qErr,pc); 


    //cout << "BEFORE qErr=" << qErr << endl;
#ifndef USE_OLD_CONSTRAINTS
    // Put position constraint equation errors in qErr
    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(c);
        const Segment& pseg = cInfo.holoErrSegment;
        if (pseg.length)
            constraints[c]->getImpl().realizePositionErrors(s, pc, pseg.length, &qErr[pseg.offset]);
    }
#else // USE_OLD_CONSTRAINTS
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(qErr,pc);
#endif

    //cout << "AFTER qErr=" << qErr << endl;
    return 0;
}

// Set generalized speeds: sweep from base to tip.
// realizePosition() must have been called already.
int SimbodyMatterSubsystemRep::realizeSubsystemVelocityImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "SimbodyMatterSubsystem::realizeVelocity()");

    const SBModelVars&     mv   = getModelVars(s);
    const SBModelCache&    mc   = getModelCache(s);
    const Vector&          q    = getQ(s);
    const SBPositionCache& pc   = getPositionCache(s);
    const Vector&          u    = getU(s);
    Vector&                uErr = updUErr(s);

    // Get the Motion-stage cache and make sure it has been allocated and initialized if needed.
    SBVelocityCache&       vc = updVelocityCache(s);
    vc.allocate(topologyCache);

    Vector& qdot = updQDot(s);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeVelocity(mv,q,pc,u,vc,qdot); 

#ifndef USE_OLD_CONSTRAINTS
    // Put velocity constraint equation errors in uErr
    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(c);

        const Segment& holoseg = cInfo.holoErrSegment; // for derivatives of holonomic constraints
        const Segment& nonholoseg = cInfo.nonholoErrSegment; // vseg includes holonomic+nonholonomic
        const int mHolo = holoseg.length, mNonholo = nonholoseg.length;
        if (mHolo)
            constraints[c]->getImpl().realizePositionDotErrors(s, vc, mHolo,    &uErr[holoseg.offset]);
        if (mNonholo)
            constraints[c]->getImpl().realizeVelocityErrors   (s, vc, mNonholo, &uErr[mc.totalNHolonomicConstraintEquationsInUse + nonholoseg.offset]);
    }
    //cout << "NEW UERR=" << uErr << endl;
#else // USE_OLD_CONSTRAINTS
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcVelInfo(pc,uErr,vc);
    //cout << "OLD UERR=" << uErr << endl;
#endif

    return 0;
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.
// Then go ask around to collect up all the applied forces from any
// force subsystems.

int SimbodyMatterSubsystemRep::realizeSubsystemDynamicsImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "SimbodyMatterSubsystem::realizeDynamics()");

    const SBModelVars&     mv = getModelVars(s);
    const SBPositionCache& pc = getPositionCache(s);
    const SBVelocityCache& vc = getVelocityCache(s);
    const Vector&          u  = getU(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    SBDynamicsCache& dc = updDynamicsCache(s);
    dc.allocate(topologyCache);

    // tip-to-base calculation
    calcArticulatedBodyInertias(s);

    // base-to-tip
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->realizeDynamics(mv,pc,u,vc,dc);

    // Update system kinetic energy
    const MultibodySystem& mbs = getMultibodySystem();  // owner of this subsystem
    mbs.getRep().updKineticEnergy(s, Stage::Dynamics) += calcKineticEnergy(s);

    return 0;
}

int SimbodyMatterSubsystemRep::realizeSubsystemAccelerationImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "SimbodyMatterSubsystem::realizeAcceleration()");

    // Get the Acceleration-stage cache and make sure it has been allocated
    // and initialized if needed.
    Vector&              udot    = updUDot(s);
    Vector&              qdotdot = updQDotDot(s);
    SBAccelerationCache& ac      = updAccelerationCache(s);
    ac.allocate(topologyCache);

    // We ask our containing MultibodySystem for a reference to the cached forces
    // accumulated from all the force subsystems. We use these to compute accelerations,
    // with all results going into the AccelerationCache.
    const MultibodySystem& mbs = getMultibodySystem();  // owner of this subsystem
    calcLoopForwardDynamics(s,
        mbs.getMobilityForces(s, Stage::Dynamics),
        mbs.getParticleForces(s, Stage::Dynamics),
        mbs.getRigidBodyForces(s, Stage::Dynamics));

    calcQDotDot(s, udot, qdotdot);

    return 0;
}


int SimbodyMatterSubsystemRep::realizeSubsystemReportImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "SimbodyMatterSubsystem::realizeReport()");

    // nothing yet

    return 0;
}

int SimbodyMatterSubsystemRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
{
    // Let the bodies and mobilizers have a chance to generate some geometry.
    for (MobilizedBodyIndex bx(0); bx<mobilizedBodies.size(); ++bx)
        mobilizedBodies[bx]->getImpl().calcDecorativeGeometryAndAppend(s,stage,geom);

    // Likewise for the constraints
    for (ConstraintIndex cx(0); cx<constraints.size(); ++cx)
        constraints[cx]->getImpl().calcDecorativeGeometryAndAppend(s,stage,geom);

    // Now add in any subsystem-level geometry.
    switch(stage) {
    case Stage::TopologyIndex: {
        assert(subsystemTopologyHasBeenRealized());
        // none yet
        break;
    }
    case Stage::PositionIndex: {
        assert(getStage(s) >= Stage::Position);
        //TODO: just to check control flow, put a ball at system COM
        //const Vec3 com = getMyMatterSubsystemHandle().calcSystemMassCenterLocationInGround(s);
        //geom.push_back(DecorativeSphere(0.02).setBodyId(GroundIndex).setTransform(com)
        //                .setColor(Green).setRepresentation(DecorativeGeometry::DrawPoints)
         //               .setResolution(1));
    }
    default: 
        assert(getStage(s) >= stage);
    }

    return 0;
}
int SimbodyMatterSubsystemRep::calcQUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNQ(s));
    weights = 1;
    /*
    Vec6 bounds;
    State tempState = s;
    getSystem().realize(tempState, Stage::Position);
    calcQUnitWeightsRecursively(s, tempState, weights, bounds, getRigidBodyNode(getGround().getMobilizedBodyIndex()));
    /**/
    return 0;
}
int SimbodyMatterSubsystemRep::calcUUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNU(s));
    weights = 1;
    /*
    Vec6 bounds;
    State tempState = s;
    tempState.updU() = 0;
    getSystem().realize(tempState, Stage::Position);
    calcUUnitWeightsRecursively(s, tempState, weights, bounds, getRigidBodyNode(getGround().getMobilizedBodyIndex()));
    /**/
    return 0;
}
void SimbodyMatterSubsystemRep::calcQUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const {
    bounds[0] = bounds[2] = bounds[4] = Infinity;
    bounds[1] = bounds[3] = bounds[5] = -Infinity;
    
    // Call this method recursively on the body's children, and build up a bounding box
    // for everything downstream of it.
    
    for (int i = 0; i < body.getNChildren(); ++i) {
        Vec6 childBounds;
        calcQUnitWeightsRecursively(s, tempState, weights, childBounds, *body.getChild(i));
        bounds[0] = std::min(bounds[0], childBounds[0]);
        bounds[2] = std::min(bounds[2], childBounds[2]);
        bounds[4] = std::min(bounds[4], childBounds[4]);
        bounds[1] = std::max(bounds[1], childBounds[1]);
        bounds[3] = std::max(bounds[3], childBounds[3]);
        bounds[5] = std::max(bounds[5], childBounds[5]);
    }
    const SBPositionCache& pc = getPositionCache(s);
    const SBModelVars& mv = getModelVars(s);
    Vec3 origin = body.getX_GB(pc).T();
    bounds[0] = std::min(bounds[0], origin[0]);
    bounds[2] = std::min(bounds[2], origin[1]);
    bounds[4] = std::min(bounds[4], origin[2]);
    bounds[1] = std::max(bounds[1], origin[0]);
    bounds[3] = std::max(bounds[3], origin[1]);
    bounds[5] = std::max(bounds[5], origin[2]);
    std::vector<Vec3> corners;
    corners.push_back(Vec3(bounds[0], bounds[2], bounds[4]));
    corners.push_back(Vec3(bounds[0], bounds[2], bounds[5]));
    corners.push_back(Vec3(bounds[0], bounds[3], bounds[4]));
    corners.push_back(Vec3(bounds[0], bounds[3], bounds[5]));
    corners.push_back(Vec3(bounds[1], bounds[2], bounds[4]));
    corners.push_back(Vec3(bounds[1], bounds[2], bounds[5]));
    corners.push_back(Vec3(bounds[1], bounds[3], bounds[4]));
    corners.push_back(Vec3(bounds[1], bounds[3], bounds[5]));
    
    // Try changing each Q, and see how much the position of each corner of the box is affected
    
    const Real delta = 1e-10;
    Transform t = ~body.getX_GB(pc);
    for (int i = 0; i < body.getNQInUse(mv); ++i) {
        int qindex = body.getQIndex()+i;
        tempState.updQ()[qindex] = s.getQ()[qindex]+delta;
        getSystem().realize(tempState, Stage::Position);
        const SBPositionCache& tempCache = getPositionCache(tempState);
        Transform deltaT = body.getX_GB(tempCache)*t;
        Real max = 0.0;
        for (int j = 0; j < (int) corners.size(); ++j) {
            Real dist = (corners[j]-deltaT*corners[j]).norm();
            max = std::max(max, dist);
        }
        weights[qindex] = std::max(max/delta, 0.1);
        tempState.updQ()[qindex] = s.getQ()[qindex];
    }
}
void SimbodyMatterSubsystemRep::calcUUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const {
    bounds[0] = bounds[2] = bounds[4] = Infinity;
    bounds[1] = bounds[3] = bounds[5] = -Infinity;
    
    // Call this method recursively on the body's children, and build up a bounding box
    // for everything downstream of it.
    
    for (int i = 0; i < body.getNChildren(); ++i) {
        Vec6 childBounds;
        calcUUnitWeightsRecursively(s, tempState, weights, childBounds, *body.getChild(i));
        bounds[0] = std::min(bounds[0], childBounds[0]);
        bounds[2] = std::min(bounds[2], childBounds[2]);
        bounds[4] = std::min(bounds[4], childBounds[4]);
        bounds[1] = std::max(bounds[1], childBounds[1]);
        bounds[3] = std::max(bounds[3], childBounds[3]);
        bounds[5] = std::max(bounds[5], childBounds[5]);
    }
    const SBPositionCache& pc = getPositionCache(s);
    const SBModelVars& mv = getModelVars(s);
    Vec3 origin = body.getX_GB(pc).T();
    bounds[0] = std::min(bounds[0], origin[0]);
    bounds[2] = std::min(bounds[2], origin[1]);
    bounds[4] = std::min(bounds[4], origin[2]);
    bounds[1] = std::max(bounds[1], origin[0]);
    bounds[3] = std::max(bounds[3], origin[1]);
    bounds[5] = std::max(bounds[5], origin[2]);
    std::vector<Vec3> corners;
    corners.push_back(Vec3(bounds[0], bounds[2], bounds[4]));
    corners.push_back(Vec3(bounds[0], bounds[2], bounds[5]));
    corners.push_back(Vec3(bounds[0], bounds[3], bounds[4]));
    corners.push_back(Vec3(bounds[0], bounds[3], bounds[5]));
    corners.push_back(Vec3(bounds[1], bounds[2], bounds[4]));
    corners.push_back(Vec3(bounds[1], bounds[2], bounds[5]));
    corners.push_back(Vec3(bounds[1], bounds[3], bounds[4]));
    corners.push_back(Vec3(bounds[1], bounds[3], bounds[5]));
    
    // Try changing each U, and see how much the position of each corner of the box is affected
    
    const Real delta = 1e-10;
    const Real timescale = getSystem().calcTimescale(s);
    Transform t = ~body.getX_GB(pc);
    for (int i = 0; i < body.getDOF(); ++i) {
        int uindex = body.getUIndex()+i;
        tempState.updU()[uindex] = 1.0;
        getSystem().realize(tempState, Stage::Velocity);
        tempState.updQ() = s.getQ()+tempState.getQDot()*(delta/timescale);
        getSystem().realize(tempState, Stage::Position);
        const SBPositionCache& tempCache = getPositionCache(tempState);
        Transform deltaT = body.getX_GB(tempCache)*t;
        Real max = 0.0;
        for (int j = 0; j < (int) corners.size(); ++j) {
            Real dist = (corners[j]-deltaT*corners[j]).norm();
            max = std::max(max, dist);
        }
        weights[uindex] = std::max(max/delta, 0.1);
        tempState.updU()[uindex] = 0.0;
        tempState.updQ() = s.getQ();
    }
}
int SimbodyMatterSubsystemRep::getQIndex(MobilizedBodyIndex body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getQIndex();
}
int SimbodyMatterSubsystemRep::getQAlloc(MobilizedBodyIndex body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getMaxNQ();
}
int SimbodyMatterSubsystemRep::getUIndex(MobilizedBodyIndex body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getUIndex();
}
int SimbodyMatterSubsystemRep::getDOF(MobilizedBodyIndex body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getDOF();
}

// We are in the process of realizingConstruction() when we need to make this call.
// We pass in the partially-completed Construction-stage cache, which must have all
// the dimensions properly filled in at this point.
void SimbodyMatterSubsystemRep::setDefaultModelValues(const SBTopologyCache& topologyCache, 
                                                      SBModelVars& modelVars) const 
{
    // Tree-level defaults
    modelVars.useEulerAngles = false;
    //modelVars.prescribed.assign(getNBodies(), false);
    for (int i = 0; i < (int)modelVars.prescribed.size(); ++i)
        modelVars.prescribed[i] = false;
    modelVars.prescribed[0] = true; // ground
    //modelVars.enabled.assign(getNConstraints(), false);
    for (int i = 0; i < (int)modelVars.disabled.size(); ++i)
        modelVars.disabled[i] = false;

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultModelValues(topologyCache, modelVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultInstanceValues(const SBModelVars& mv, 
                                                         SBInstanceVars& instanceVars) const 
{
    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultInstanceValues(mv, instanceVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultTimeValues(const SBModelVars& mv, 
                                         SBTimeVars& timeVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultTimeValues(mv, timeVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
{
    // Tree-level defaults (none)


    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultPositionValues(mv, q);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultVelocityValues(const SBModelVars& mv, Vector& u) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultVelocityValues(mv, u);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultDynamicsValues(const SBModelVars& mv, 
                                             SBDynamicsVars& dynamicsVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultDynamicsValues(mv, dynamicsVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultAccelerationValues(const SBModelVars& mv, 
                                             SBAccelerationVars& reactionVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultAccelerationValues(mv, reactionVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void SimbodyMatterSubsystemRep::setMobilizerIsPrescribed(State& s, MobilizedBodyIndex body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void SimbodyMatterSubsystemRep::setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disable) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.disabled[constraint] = disable;   
}

bool SimbodyMatterSubsystemRep::getUseEulerAngles(const State& s) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.useEulerAngles;
}
bool SimbodyMatterSubsystemRep::isMobilizerPrescribed(const State& s, MobilizedBodyIndex body) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.prescribed[body];
}
bool SimbodyMatterSubsystemRep::isConstraintDisabled(const State& s, ConstraintIndex constraint) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.disabled[constraint];
}

void SimbodyMatterSubsystemRep::convertToEulerAngles(const State& inputState, State& outputState) const {
    assert(!getUseEulerAngles(inputState));
    setUseEulerAngles(outputState, true);
    getMultibodySystem().realizeModel(outputState);
    outputState.updU() = inputState.getU();
    const Vector& inputQ = inputState.getQ();
    Vector& outputQ = outputState.updQ();
    for (int i = 0; i < getNumMobilizedBodies(); ++i) {
        const RigidBodyNode& node = getRigidBodyNode(i);
        node.convertToEulerAngles(inputQ, outputQ);
    }
}

void SimbodyMatterSubsystemRep::convertToQuaternions(const State& inputState, State& outputState) const {
    assert(getUseEulerAngles(inputState));
    setUseEulerAngles(outputState, false);
    getMultibodySystem().realizeModel(outputState);
    outputState.updU() = inputState.getU();
    const Vector& inputQ = inputState.getQ();
    Vector& outputQ = outputState.updQ();
    for (int i = 0; i < getNumMobilizedBodies(); ++i) {
        const RigidBodyNode& node = getRigidBodyNode(i);
        node.convertToQuaternions(inputQ, outputQ);
    }
}

bool SimbodyMatterSubsystemRep::isUsingQuaternion(const State& s, MobilizedBodyIndex body) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    MobilizerQIndex startOfQuaternion; // we don't need this information here
    return n.isUsingQuaternion(getModelVars(s), startOfQuaternion);
}

int SimbodyMatterSubsystemRep::getNQuaternionsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.totalNQuaternionsInUse;
}

QuaternionPoolIndex SimbodyMatterSubsystemRep::getQuaternionPoolIndex
   (const State& s, MobilizedBodyIndex body) const
{
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.getMobilizedBodyModelInfo(body).quaternionPoolIndex;
}

AnglePoolIndex SimbodyMatterSubsystemRep::getAnglePoolIndex
   (const State& s, MobilizedBodyIndex body) const
{
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.getMobilizedBodyModelInfo(body).anglePoolIndex;
}

int SimbodyMatterSubsystemRep::getNHolonomicConstraintEquationsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.totalNHolonomicConstraintEquationsInUse;
}
int SimbodyMatterSubsystemRep::getNNonholonomicConstraintEquationsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.totalNNonholonomicConstraintEquationsInUse;
}
int SimbodyMatterSubsystemRep::getNAccelerationOnlyConstraintEquationsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.totalNAccelerationOnlyConstraintEquationsInUse;
}

void SimbodyMatterSubsystemRep::calcHolonomicConstraintMatrixPQInverse(const State& s, Matrix& PQInv) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    const int mp = mc.totalNHolonomicConstraintEquationsInUse;
    const int nq = getNQ(s);

    PQInv.resize(mp,nq);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(c);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into qErr and mHolo (mp)

        PQInv(holoSeg.offset, 0, holoSeg.length, nq) = constraints[c]->calcPositionConstraintMatrixPQInverse(s);
    }
}
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixP(const State& s, Matrix& P) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mHolo = model.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    P.resize(mHolo,nu);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        P(holoSeg.offset, 0, holoSeg.length, nu) = constraints[c]->calcPositionConstraintMatrixP(s);
    }
}
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixPt(const State& s, Matrix& Pt) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mHolo = model.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Pt.resize(nu,mHolo);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        // Fill in columns of Pt
        Pt(0, holoSeg.offset, nu, holoSeg.length) = constraints[c]->calcPositionConstraintMatrixPt(s);
    }
}
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixV(const State& s, Matrix& V) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mNonholo = model.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    V.resize(mNonholo,nu);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        V(nonholoSeg.offset, 0, nonholoSeg.length, nu) = constraints[c]->calcVelocityConstraintMatrixV(s);
    }
}
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixVt(const State& s, Matrix& Vt) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mNonholo = model.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Vt.resize(nu,mNonholo);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        // Fill in columns of Vt
        Vt(0, nonholoSeg.offset, nu, nonholoSeg.length) = constraints[c]->calcVelocityConstraintMatrixVt(s);
    }
}
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixA (const State& s, Matrix& A) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mAccOnly = model.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    A.resize(mAccOnly,nu);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        A(accOnlySeg.offset, 0, accOnlySeg.length, nu) = constraints[c]->calcAccelerationConstraintMatrixA(s);
    }
}
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixAt(const State& s, Matrix& At) const {
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mAccOnly = model.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    At.resize(nu,mAccOnly);

    for (ConstraintIndex c(0); c < constraints.size(); ++c) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(c);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        At(0, accOnlySeg.offset, nu, accOnlySeg.length) = constraints[c]->calcAccelerationConstraintMatrixAt(s);
    }
}
// Must be realized to Stage::Position
void SimbodyMatterSubsystemRep::calcConstraintForcesFromMultipliers
  (const State& s, const Vector& lambda,
   Vector_<SpatialVec>& bodyForcesInG,
   Vector&              mobilityForces) const
{
    const SBModelCache& model = getModelCache(s); // must be >=Model stage
    const int mHolo    = model.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = model.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = model.totalNAccelerationOnlyConstraintEquationsInUse;
    const int ma = mHolo+mNonholo+mAccOnly;

    bodyForcesInG.resize(getNBodies()); bodyForcesInG.setToZero();
    mobilityForces.resize(getNU(s));    mobilityForces.setToZero();

    Vector_<SpatialVec> bodyF1;          // per constraint
    Vector              mobilityF1;
    Vector              lambda1; // multipliers for 1 constraint
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        const SBModelCache::PerConstraintModelInfo& cInfo = model.getConstraintModelInfo(cx);
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mh=holoSeg.length, mnh=nonholoSeg.length, mao=accOnlySeg.length;

        lambda1.resize(mh+mnh+mao);
        lambda1(0,mh)        = lambda(holoSeg.offset, mh);
        lambda1(mh, mnh)     = lambda(mHolo+nonholoSeg.offset, mnh);
        lambda1(mh+mnh, mao) = lambda(mHolo+mNonholo+accOnlySeg.offset, mao);

        const ConstraintRep& crep = constraints[cx]->getImpl();
 
        bodyF1.resize(crep.getNumConstrainedBodies());
        mobilityF1.resize(crep.getNumConstrainedU(s));
        crep.calcConstraintForcesFromMultipliers(s,mh,mnh,mao,&lambda1[0],bodyF1,mobilityF1);

        if (crep.getNumConstrainedBodies()) {
            const Rotation& R_GA = crep.getAncestorMobilizedBody().getBodyRotation(s);
            for (ConstrainedBodyIndex cb(0); cb < crep.getNumConstrainedBodies(); ++cb)
                bodyForcesInG[crep.getMobilizedBodyIndexOfConstrainedBody(cb)] += R_GA*bodyF1[cb];
        }

        for (ConstrainedUIndex cux(0); cux < crep.getNumConstrainedU(s); ++cux) 
            mobilityForces[crep.getUIndexOfConstrainedU(s,cux)] += mobilityF1[cux];
    }
}

static Real calcQErrestWeightedNorm(const SimbodyMatterSubsystemRep& matter, const State& s, const Vector& qErrest, const Vector& uWeights) {
    Vector qhatErrest(uWeights.size());
    matter.multiplyByQMatrixInverse(s, false, qErrest, qhatErrest); // qhatErrest = Q+ qErrest
    qhatErrest.rowScaleInPlace(uWeights);                    // qhatErrest = Wu Q+ qErrest
    return qhatErrest.normRMS();
}

void SimbodyMatterSubsystemRep::enforcePositionConstraints
   (State& s, Real consAccuracy, const Vector& yWeights,
    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{
    assert(getStage(s) >= Stage::Position-1);
    realizeSubsystemPosition(s);

    // First work only with the holonomic (position) constraints, which appear first in the QErr array.
    // Don't work on the quaternion constraints in this first section.
    const int mHolo  = getNHolonomicConstraintEquationsInUse(s);
    const int mQuats = getNQuaternionsInUse(s);
    const int nq     = getNQ(s);
    const int nu     = getNU(s);

    const VectorView uWeights   = yWeights(nq,nu);
    const Vector     ooUWeights = uWeights.elementwiseInvert();
    //Real oow[] = {1., 1., 1., .05, .05, .05};
   // Vector     ooUWeights(nu, oow);
    const VectorView ooPTols  = ooTols(0,mHolo);

    VectorView qErrest = yErrest.size() ? yErrest(0,nq) : yErrest(0,0);

    // This is a const view into the State; the contents it refers to will change though.
    const VectorView pErrs    = getQErr(s)(0,mHolo); // just leaving off quaternions

    bool anyChange = false;

    // Solve 
    //   (Tp P Wu^-1) dqhat_WLS = T perr, q -= Q*Wu^-1*dqhat_WLS
    // until perr(q)_TRMS <= 0.1*accuracy.
    //
#ifndef USE_OLD_CONSTRAINTS
    // This is a nonlinear least squares problem. This is a full Newton iteration since we
    // recalculate the iteration matrix each time around the loop. TODO: a solution could
    // be found using the same iteration matrix, since we are projecting from (presumably)
    // not too far away. Q1: will it be the same solution? Q2: if not, does it matter?
    Vector scaledPerrs = pErrs.rowScale(ooPTols);
    Real normAchievedTRMS = scaledPerrs.normRMS();

    /*
    cout << "!!!! initially @" << s.getTime() << ", perr TRMS=" << normAchievedTRMS << " consAcc=" << consAccuracy;
    if (qErrest.size())
        cout << " qErrest WRMS=" << calcQErrestWeightedNorm(*this,s,qErrest,uWeights);
    else cout << " NO Q ERROR ESTIMATE";
    cout << endl;
    cout << "!!!! PERR=" << pErrs << endl;
    */

    Real lastChangeMadeWRMS = 0; // size of last change in weighted dq
    int nItsUsed = 0;
    if (normAchievedTRMS > consAccuracy) {
        Vector saveQ = getQ(s);
        Matrix Pt(nu,mHolo);
        Vector dqhat_WLS(nu);
        Vector dq(nq); // = Q Wu^-1 dqhat_WLS
        FactorQTZ P_qtz;
        const int MaxIterations  = 20;
        do {
            //Matrix P(mHolo,nu);
            //calcHolonomicVelocityConstraintMatrixP(s, P); // mp X nu
            //cout << "  P(verr)=" << P;
            calcHolonomicVelocityConstraintMatrixPt(s, Pt); // nu X mp
            //cout << "  P(lambda)=" << ~Pt;
            //P.rowAndColScaleInPlace(ooPTols, ooUWeights); // P is now Tp P Wu^-1
            Pt.rowAndColScaleInPlace(ooUWeights, ooPTols); // Pt is now ~(TPW)=Wu^-1 Pt Tp (weights are symmetric)
            //cout << "TPW-1=" << P;
            //cout << "TPW-1=" << ~Pt;

            P_qtz.factor<Real>(~Pt, std::max(Pt.nrow(),Pt.ncol()) * SignificantReal); // this acts like a pseudoinverse
            P_qtz.solve(scaledPerrs, dqhat_WLS);
            lastChangeMadeWRMS = dqhat_WLS.normRMS(); // size of change in weighted norm
            //cout << "!!!! dqhat weighted=" << dqhat_WLS << endl;
            dqhat_WLS.rowScaleInPlace(ooUWeights);    // switch back to unweighted dqhat=Wu^-1*dqhat_WLS
            //cout << "!!!! dqhat unweighted=" << dqhat_WLS << endl;
            multiplyByQMatrix(s, false, dqhat_WLS, dq); // Q*(Wu^-1 dqhat_WLS)
            //cout << "!!!! dq unweighted=" << dq << endl;
            updQ(s) -= dq; // this is actually unweighted dq
            anyChange = true;

            // Now recalculate the position constraint errors at the new q.
            realizeSubsystemPosition(s);
            //cout << "!!!! PERR=" << pErrs << endl;
            scaledPerrs = pErrs.rowScale(ooPTols);
            normAchievedTRMS = scaledPerrs.normRMS();
            ++nItsUsed;
            //cout << "  !! iter " << nItsUsed << ": TRMS(perr)=" << normAchievedTRMS << ", WRMS(dq)=" << lastChangeMadeWRMS << endl;
        } while (normAchievedTRMS > 0.1*consAccuracy 
                 && nItsUsed < MaxIterations);

        if (nItsUsed >= MaxIterations && normAchievedTRMS > consAccuracy) {
            updQ(s) = saveQ;
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "maxIters exceeded in position projection");
        }


        // Next, if we projected out the position constraint errors, remove the corresponding
        // error from the integrator's error estimate.
        //   (Tp P Wu^-1) dqhat_WLS = (Tp P Wu^-1) Wu Q^+ qErrest, qErrest -= Q Wu^-1 dqhat_WLS
        // No iteration is required.
        if (qErrest.size()) {
            // Work in qhat W-norm
            Vector qhatErrest(nu), dqErrest(nq);
            multiplyByQMatrixInverse(s, false, qErrest, qhatErrest);
            qhatErrest.rowScaleInPlace(uWeights); // qbarErrest = Wu * Q^+ * qErrest

            P_qtz.solve(~Pt*qhatErrest, dqhat_WLS);
            const Real normOfAdjustment_WRMS = dqhat_WLS.normRMS();

            dqhat_WLS.rowScaleInPlace(ooUWeights); // unscale the result
            multiplyByQMatrix(s, false, dqhat_WLS, dqErrest);
            qErrest -= dqErrest; // unweighted
            //cout << "  !! FIXUP: now WRMS(qerrest)=" << calcQErrestWeightedNorm(*this,s,qErrest,uWeights) << " using WRMS(dq)=" << normOfAdjustment_WRMS << endl;
        }
    }

    //cout << "!!!! perr TRMS achieved " << normAchievedTRMS << " in " << nItsUsed << " iterations"  << endl;
    //cout << "!!!! ... PERR=" << pErrs << endl;
#else // USE_OLD_CONSTRAINTS
    // First, fix the position constraints produced by defined length constraints.
    if (lConstraints->enforcePositionConstraints(s, consAccuracy, 0.1*consAccuracy))
        anyChange = true;
#endif

    // By design, normalization of quaternions can't have any effect on the length
    // constraints we just fixed (because we normalize internally for calculations).
    // So now we can simply normalize the quaternions.
    const SBModelVars& mv = getModelVars(s);

    if (mQuats) {
        //cout << "!!!! QUAT START: errs=" << getQErr(s)(mHolo, mQuats) << " RMS(qErrest)=" << qErrest.normRMS() << endl;
        Vector& q  = updQ(s); //TODO: this invalidates q's already

        for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
            for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
                if (rbNodeLevels[i][j]->enforceQuaternionConstraints(mv,q,qErrest))
                    anyChange = true;

        //TODO: shouldn't need this
        realizeSubsystemPosition(s);

        //cout << "!!!! ... QUAT END: errs=" << getQErr(s)(mHolo, mQuats) << " RMS(qErrest)=" << qErrest.normRMS() << endl;
       //TODO: quaternion constraints shouldn't invalidate anything except
        // the qnorms, which will be all 1 now
    }

    if (anyChange)
        s.invalidateAll(Stage::Position);
}

void SimbodyMatterSubsystemRep::enforceVelocityConstraints
   (State& s, Real consAccuracy, const Vector& yWeights,
    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{
    assert(getStage(s) >= Stage::Velocity-1);
    realizeSubsystemVelocity(s);

    // Here we deal with the nonholonomic (velocity) constraints and the 
    // derivatives of the holonomic constraints.
    const int mHolo    = getNHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNNonholonomicConstraintEquationsInUse(s);
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    const VectorView uWeights   = yWeights(nq,nu); // skip the first nq weights
    const Vector     ooUWeights = uWeights.elementwiseInvert();
    const VectorView ooPVTols   = ooTols(0,mHolo+mNonholo);
    //TODO: scale holo part by time scale

    VectorView uErrest = yErrest.size() ? yErrest(nq,nu) : yErrest(0,0);

    // This is a const view into the State; the contents it refers to will change though.
    const Vector& vErrs = getUErr(s); // all the velocity constraint errors (mHolo+mNonholo)

    bool anyChange = false;

    // Solve 
    //   (Tpv [P;V] Wu^-1) du_WLS = Tpv uerr, u -= Wu^-1*du_WLS
    // Note that although this is a nonlinear least squares problem since uerr is a function
    // of u, we do not need to refactor the matrix since it does not depend on u.
    // TODO: Tp P Wu^-1 should already have been calculated for position projection (at least
    // if any position projection occurred)

#ifndef USE_OLD_CONSTRAINTS
    // This is a nonlinear least squares problem, but we only need to factor once since only
    // the RHS is dependent on u. 
    Vector scaledVerrs = vErrs.rowScale(ooPVTols);
    Real normAchievedTRMS = scaledVerrs.normRMS();
    
    /*
    cout << "!!!! initially @" << s.getTime() << ", verr TRMS=" << normAchievedTRMS << " consAcc=" << consAccuracy;
    if (uErrest.size())
        cout << " uErrest WRMS=" << uErrest.rowScale(uWeights).normRMS();
    else cout << " NO U ERROR ESTIMATE";
    cout << endl;
    cout << "!!!! VERR=" << vErrs << endl;
    */

    Real lastChangeMadeWRMS = 0;
    int nItsUsed = 0;
    if (normAchievedTRMS > consAccuracy) {
        const Vector saveU = getU(s);
        Matrix PVt(nu, mHolo+mNonholo);
        calcHolonomicVelocityConstraintMatrixPt(s, PVt(  0,  0,nu, mHolo));    // nu X mp
        calcNonholonomicConstraintMatrixVt     (s, PVt(0,mHolo,nu, mNonholo)); // nu X mv
        PVt.rowAndColScaleInPlace(ooUWeights, ooPVTols); // PVt is now Wu^-1 (Pt Vt) Tpv

        // Calculate pseudoinverse
        FactorQTZ PVqtz;
        PVqtz.factor<Real>(~PVt, std::max(PVt.nrow(), PVt.ncol()) * SignificantReal);

        Vector du_WLS(nu);
        const int MaxIterations  = 7;
        do {
            PVqtz.solve(scaledVerrs, du_WLS);
            lastChangeMadeWRMS = du_WLS.normRMS(); // size of change in weighted norm
            du_WLS.rowScaleInPlace(ooUWeights); // remove scaling: du=Wu^-1*du_WLS
            updU(s) -= du_WLS; // this is actually unweighted du
            anyChange = true;

            // Recalculate the constraint errors for the new u's.
            realizeSubsystemVelocity(s);
            scaledVerrs = vErrs.rowScale(ooPVTols);
            normAchievedTRMS=scaledVerrs.normRMS();
            ++nItsUsed;
            //cout << "  !! iter " << nItsUsed << ": TRMS(verr)=" << normAchievedTRMS << ", WRMS(du)=" << lastChangeMadeWRMS << endl;
        } while (normAchievedTRMS > 0.1*consAccuracy 
                 && nItsUsed < MaxIterations);

        if (nItsUsed >= MaxIterations && normAchievedTRMS > consAccuracy) {
            updU(s) = saveU;
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "maxIters exceeded in velocity projection");
        }

        // Next, if we projected out the velocity constraint errors, remove the corresponding
        // error from the integrator's error estimate.
        //   (Tpv [P;V] Wu^-1) du_WLS = (Tpv [P;V] Wu^-1) Wu*uErrEst, uErrEst = Wu^-1(Wu*uErrEst - du_WLS)
        // No iteration is required.
        if (uErrest.size()) {
            uErrest.rowScaleInPlace(uWeights); // uErrest = Wu*uErrEst
            PVqtz.solve(~PVt*uErrest, du_WLS);
            uErrest -= du_WLS; // still weighted
            //cout << "  !! U FIXUP: now WRMS(uerrest)=" << uErrest.normRMS() << " using WRMS(du)=" << du_WLS.normRMS() << endl;
            uErrest.rowScaleInPlace(ooUWeights); // back to unscaled error estimates
        }
    }

    
    //cout << "!!!! verr achieved " << normAchievedTRMS << " in " << nItsUsed << " iterations" << endl;
    //cout << "!!!! ... VERR=" << vErrs << endl;
#else // USE_OLD_CONSTRAINTS
    // Fix the velocity constraints produced by defined length constraints.
    if (lConstraints->enforceVelocityConstraints(s, consAccuracy, 0.1*consAccuracy))
        anyChange = true;
#endif

    if (anyChange)
        s.invalidateAll(Stage::Velocity);
}

// TREE FORWARD DYNAMICS OPERATOR
//
// Given a State realized through Stage::Dynamics, and a complete
// set of applied forces, calculate all acceleration results into
// the return arguments here. This routine *does not* affect the
// State cache -- it is an operator. In typical usage, the output
// arguments actually will be part of the state cache to effect
// a response, but this method can also be used to effect an
// operator.
//
// Note that although acceleration constraint errors will be
// calculated, the returned accelerations will not obey the
// constraints, unless the supplied forces already account
// for constraints. The argument list allows for some extra forces
// to be supplied, with the intent that these will be used to
// deal with internal forces generated by constraints. 
void SimbodyMatterSubsystemRep::calcTreeForwardDynamicsOperator(
    const State&               s,
    const Vector&              mobilityForces,
    const Vector_<Vec3>&       particleForces,
    const Vector_<SpatialVec>& bodyForces,
    const Vector*              extraMobilityForces,
    const Vector_<SpatialVec>* extraBodyForces,
    SBAccelerationCache&       ac,
    Vector&                    udot,
    Vector&                    udotErr) const
{
    const SBModelCache&    mc = getModelCache(s);
    const SBPositionCache& pc = getPositionCache(s);
    const SBVelocityCache& vc = getVelocityCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);

    // Ensure that output arguments have been allocated properly.
    ac.allocate(topologyCache);
    udot.resize(topologyCache.nDOFs);

    Vector              totalMobilityForces;
    Vector_<SpatialVec> totalBodyForces;

    // inputs

    // First assume we'll use the input forces as-is.
    const Vector*              mobilityForcesToUse  = &mobilityForces;
    const Vector_<SpatialVec>* bodyForcesToUse      = &bodyForces;

    if (extraMobilityForces) {
        totalMobilityForces = mobilityForces + *extraMobilityForces;
        mobilityForcesToUse = &totalMobilityForces;
    }

    if (extraBodyForces) {
        totalBodyForces = bodyForces + *extraBodyForces;
        bodyForcesToUse = &totalBodyForces;
    }

    // outputs
    Vector&              netHingeForces = ac.netHingeForces;
    Vector_<SpatialVec>& A_GB           = ac.bodyAccelerationInGround;

    calcTreeAccelerations(s, *mobilityForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot);


#ifndef USE_OLD_CONSTRAINTS
    // Put acceleration constraint equation errors in uErr
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        const SBModelCache::PerConstraintModelInfo& cInfo = mc.getConstraintModelInfo(cx);

        const Segment& holoseg    = cInfo.holoErrSegment;    // for 2nd derivatives of holonomic constraints
        const Segment& nonholoseg = cInfo.nonholoErrSegment; // for 1st derivatives of nonholonomic constraints
        const Segment& acconlyseg = cInfo.accOnlyErrSegment; // for acceleration-only constraints
        const int mHolo = holoseg.length, mNonholo = nonholoseg.length, mAccOnly = acconlyseg.length;
        if (mHolo)
            constraints[cx]->getImpl().realizePositionDotDotErrors(s, ac, mHolo,
                &udotErr[holoseg.offset]);
        if (mNonholo)
            constraints[cx]->getImpl().realizeVelocityDotErrors(s, ac, mNonholo, 
                &udotErr[  mc.totalNHolonomicConstraintEquationsInUse 
                         + nonholoseg.offset]);
        if (mAccOnly)
            constraints[cx]->getImpl().realizeAccelerationErrors(s, ac, mAccOnly, 
                &udotErr[  mc.totalNHolonomicConstraintEquationsInUse
                         + mc.totalNNonholonomicConstraintEquationsInUse 
                         + acconlyseg.offset]);
    }
    //cout << "Tree:NEW UDOT ERR=" << udotErr << endl;

#else // USE_OLD_CONSTRAINTS
    Vector oldUDotErr = udotErr; oldUDotErr = NaN;
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(pc,vc,oldUDotErr,ac);

    //cout << "Tree:OLD UDOT ERR=" << oldUDotErr << endl;
    udotErr=oldUDotErr;
#endif
}


// LOOP FORWARD DYNAMICS OPERATOR
// 
// Given a State realized through Stage::Dynamics, and a complete
// set of applied forces, calculate all acceleration results resulting
// from those forces AND enforcement of the acceleration constraints.
// The results go into the return arguments here. This routine *does not* affect the
// State cache -- it is an operator. In typical usage, the output
// arguments actually will be part of the state cache to effect
// a response, but this method can also be used to effect an
// operator.
void SimbodyMatterSubsystemRep::calcLoopForwardDynamicsOperator(const State& s, 
    const Vector&              mobilityForces,
    const Vector_<Vec3>&       particleForces,
    const Vector_<SpatialVec>& bodyForces,
    SBAccelerationCache&       ac,
    Vector&                    udot,
    Vector&                    multipliers,
    Vector&                    udotErr) const
{
    assert(getStage(s) >= Stage::Acceleration-1);


    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    0, 0, ac, udot, udotErr);

    const int mHolo    = getNHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = getNAccelerationOnlyConstraintEquationsInUse(s);
    const int ma = mHolo+mNonholo+mAccOnly;
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    if (ma==0 || nu==0) {
        multipliers.resize(0);
        return;
    }

    //cout << "---> BEFORE udotErr=" << udotErr << endl;


#ifndef USE_OLD_CONSTRAINTS
    Matrix Gt(nu,ma); // Gt==~P ~V ~A
    // Fill in all the columns of Gt
    calcHolonomicVelocityConstraintMatrixPt(s, Gt(0,     0,          nu, mHolo));
    calcNonholonomicConstraintMatrixVt     (s, Gt(0,   mHolo,        nu, mNonholo));
    calcAccelerationOnlyConstraintMatrixAt (s, Gt(0, mHolo+mNonholo, nu, mAccOnly));

    // Calculate multipliers lambda as
    //     (G M^-1 ~G) lambda = aerr
    // TODO: Optimally, we would calculate this mXm matrix in O(m^2) time. Then 
    // we'll factor it in O(m^3) time. I don't know how to calculate it that
    // fast, but using m calls to the M^-1*f and G*udot O(n) operators
    // we can calculate it in O(mn) time. As long as m << n, and
    // especially if m is a small constant independent of n, and even better
    // if we've partitioned it into little 
    // subblocks, this is all very reasonable. One slip up and you'll toss in a
    // factor of mn^2 or m^2n and screw this up -- be careful! (see below)
    Matrix MInvGt(nu, ma);
    Vector_<SpatialVec> A_GB(getNBodies()); // dummy
    for (int j=0; j<ma; ++j) // This is O(mn)
        calcMInverseF(s, Gt(j), A_GB, MInvGt(j));

    // TODO: Toldya! Check out this m^2n bit here ...
    Matrix GMInvGt = (~Gt)*MInvGt; // TODO: BAD!!! O(m^2n) -- Use m x G udot operators instead for O(mn)
    FactorQTZ qtz(GMInvGt, ma*SignificantReal); // specify 1/cond at which we declare rank deficiency
    Vector lambda(ma);
    qtz.solve(udotErr, lambda);
    /*cout << "qtz.getRank()=" << qtz.getRank() << endl;
    if (qtz.getRank() == 1) 
        cout << "OH MY GOD!! RANK 1" << endl;
    cout << "Solve GM^-1Gt lambda=rhs; GM^-1Gt=" << GMInvGt;
    cout << "  rhs=" << udotErr << endl;
    */
    //cout << "---> NEW lambda=" << lambda << endl;

    Vector_<SpatialVec> bodyF;
    Vector mobilityF;
    calcConstraintForcesFromMultipliers(s,lambda,bodyF,mobilityF);
    //TODO: this is wrong
        bodyF *= -1;

        //cout << "  NEW FORCES: " << bodyF << endl;

    multipliers = lambda;
    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    0, &bodyF, ac, udot, udotErr);
        //cout << "  NEW UDOTERR=" << udotErr << endl;
#else // USE_OLD_CONSTRAINTS
    Vector_<SpatialVec> cFrc(getNBodies()); 
    cFrc.setToZero();

    if (lConstraints->calcConstraintForces(s, udotErr, multipliers, ac)) {
        lConstraints->addInCorrectionForces(s, ac, cFrc);


        //cout << "  OLD FORCES: " << cFrc << endl;

        calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                        0, &cFrc, ac, udot, udotErr);
        //cout << "  OLD UDOTERR=" << udotErr << endl;
    }
#endif
}

// This is the response version of the above operator; that is, it uses
// the operator but puts the results in the State cache. Note that this
// only makes sense if the force arguments also come from the State
// somewhere else in the System that includes this Subsystem.
void SimbodyMatterSubsystemRep::calcTreeForwardDynamics(
    const State&               s,
    const Vector&              mobilityForces,
    const Vector_<Vec3>&       particleForces,
    const Vector_<SpatialVec>& bodyForces,
    const Vector*              extraMobilityForces,
    const Vector_<SpatialVec>* extraBodyForces) const
{
    // Output goes into State.
    SBAccelerationCache& ac      = updAccelerationCache(s);
    Vector&              udot    = updUDot(s);
    Vector&              udotErr = updUDotErr(s);

    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    extraMobilityForces, extraBodyForces,
                                    ac, udot, udotErr);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void SimbodyMatterSubsystemRep::calcLoopForwardDynamics(const State& s, 
    const Vector&               mobilityForces,
    const Vector_<Vec3>&        particleForces,
    const Vector_<SpatialVec>&  bodyForces) const 
{
    // Output goes into State.
    SBAccelerationCache& ac      = updAccelerationCache(s);
    Vector&              udot    = updUDot(s);
    Vector&              udotErr = updUDotErr(s);
    Vector&              multipliers = updMultipliers(s);

    calcLoopForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    ac, udot, multipliers, udotErr);
}

/* TODO:
// Calculate the position (holonomic) constraint matrix P for all mp position
// constraints in the system.
// The returned matrix is mp X n (where n is the number of mobilities and u's).
// The time complexity is also O(mp*n).
void SimbodyMatterSubsystemRep::calcPositionConstraintMatrix(const State& s,
    Matrix& P) const 
{
    const SBPositionCache& pc = getPositionCache(s);

    // The first time derivative of the position constraint error methods
    // contains the matrix we're interested in, like this (at the current configuration):
    //     positionDotError(u) = Pu - c(t) = 0
    // We can extract columns of P by setting a single u to 1 and the rest 0,
    // but first we need to evaluate the bias term -c(t), which we get when
    // u=0.

    // Evaluate pdot errors at u=0, save -c(t) 
    Vector pdotBias(getNumPositionConstraints());
    calcPositionDotBias(s, pdotBias);

    P.resize(getNumPositionConstraints(), getNumMobilities());
    for (int i=0; i < getNumMobilities(s); ++i) {
        calcPositionDotBiasedColumn(s, i, P(i)); // u[i]=1, all others 0
        P(i) -= pdotBias;
    }

    // ALTERNATIVE:
    // Use the applyPositionConstraintForces methods instead to calculate columns of ~P.

    P.resize(getNumPositionConstraints(), getNumMobilities());
    MatrixView Pt = ~P;

    Vector              multipliers(getNumPositionConstraints());
    Vector_<SpatialVec> bodyForces(getNBodies());
    Vector              directMobilityForces(getNumMobilities(s));

    multipliers.setToZero();
    for (int i=0; i < getNumPositionConstraints(); ++i) {
        multipliers[i] = 1;

        bodyForces.setToZero(); directMobilityForces.setToZero();
        applyPositionConstraintForces(s, multipliers.size(), &multipliers[0],
            bodyForces, directMobilityForces);

        calcInternalGradientFromSpatial(s, bodyForces, Pt(i));
        Pt(i) += directMobilityForces;

        multipliers[i] = 0;
    }

}
*/

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void SimbodyMatterSubsystemRep::calcArticulatedBodyInertias(const State& s) const {
    const SBPositionCache& pc = getPositionCache(s);
    SBDynamicsCache&       dc = updDynamicsCache(s);

    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcArticulatedBodyInertiasInward(pc,dc);
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
/*
void SimbodyMatterSubsystemRep::calcZ(const State& s, 
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces) const
{
    const SBPositionCache&    pc = getPositionCache(s);
    const SBDynamicsCache&    dc = getDynamicsCache(s);
    SBAccelerationCache&      ac = updAccelerationCache(s);

    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(pc,dc,mobilityForces,bodyForces,ac);
        }
}
*/
void SimbodyMatterSubsystemRep::calcZ(const State& s, 
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces) const
{
    const SBStateDigest sbs(s, *this, Stage::Acceleration);

    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(sbs,mobilityForces,bodyForces);
        }
}

// Y is used for length constraints: sweep from base to tip.
void SimbodyMatterSubsystemRep::calcY(const State& s) const {
    const SBPositionCache& pc = getPositionCache(s);
    SBDynamicsCache&       dc = updDynamicsCache(s);

    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcYOutward(pc,dc);
}

// Calc acceleration: sweep from base to tip.
void SimbodyMatterSubsystemRep::calcTreeAccel(const State& s) const {
    const SBModelVars&     mv      = getModelVars(s);
    const Vector&          q       = getQ(s);
    const SBPositionCache& pc      = getPositionCache(s);
    const Vector&          u       = getU(s);
    const SBDynamicsCache& dc      = getDynamicsCache(s);

    SBAccelerationCache&   ac      = updAccelerationCache(s);
    Vector&                udot    = updUDot(s);
    Vector&                qdotdot = updQDotDot(s);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel(mv,q,pc,u,dc,ac,udot,qdotdot);
}

void SimbodyMatterSubsystemRep::fixVel0(State& s, Vector& vel) const {
    lConstraints->fixVel0(s, vel);
}

Real SimbodyMatterSubsystemRep::calcKineticEnergy(const State& s) const {
    const SBPositionCache& pc = getPositionCache(s);
    const SBVelocityCache& vc = getVelocityCache(s);

    Real ke = 0.;

    // Skip ground level 0!
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            ke += rbNodeLevels[i][j]->calcKineticEnergy(pc,vc);

    return ke;
}

//
// Operator for open-loop dynamics.
//
void SimbodyMatterSubsystemRep::calcTreeAccelerations(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBStateDigest sbs(s, *this, Stage::Acceleration);
    const SBPositionCache& pc = sbs.getPositionCache();
    const SBDynamicsCache& dc = sbs.getDynamicsCache();

    assert(mobilityForces.size() == getTotalDOF());
    assert(bodyForces.size() == getNBodies());

    netHingeForces.resize(getTotalDOF());
    A_GB.resize(getNBodies());
    udot.resize(getTotalDOF());

    // Temporaries
    Vector_<SpatialVec> allZ(getNBodies());
    Vector_<SpatialVec> allGepsilon(getNBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass1Inward(pc,dc,
                mobilityForces, bodyForces, allZ, allGepsilon,
                netHingeForces);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass2Outward(pc,dc, netHingeForces, A_GB, udot);
        }
}

//
// Calculate udot = M^-1 f. We also get spatial accelerations A_GB for 
// each body as a side effect.
//
void SimbodyMatterSubsystemRep::calcMInverseF(const State& s,
    const Vector&              f,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBPositionCache& pc = getPositionCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);

    assert(f.size() == getTotalDOF());

    A_GB.resize(getNBodies());
    udot.resize(getTotalDOF());

    // Temporaries
    Vector              allEpsilon(getTotalDOF());
    Vector_<SpatialVec> allZ(getNBodies());
    Vector_<SpatialVec> allGepsilon(getNBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass1Inward(pc,dc,
                f, allZ, allGepsilon,
                allEpsilon);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass2Outward(pc,dc, allEpsilon, A_GB, udot);
        }
}

// q=Qu or u=~Qq
void SimbodyMatterSubsystemRep::multiplyByQMatrix(const State& s, bool transpose, const Vector& in, Vector& out) const
{
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next()); // i.e., we must be *done* with Stage::Position

    assert(in.size() == (transpose?getTotalQAlloc():getTotalDOF()));
    out.resize(transpose?getTotalDOF():getTotalQAlloc());

    //TODO: this shouldn't be necessary
    assert(in.size()  < 2 || &in[1]  == &in[0] +1); // for now must be contiguous in memory
    assert(out.size() < 2 || &out[1] == &out[0]+1); 

    const Real* qp   = sbState.getQ();
    const Real* inp  = in.size()  ? &in[0]  : 0;
    Real*       outp = out.size() ? &out[0] : 0;

    const bool useEulerAngles = sbState.getModelVars().useEulerAngles;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& rbn = *rbNodeLevels[i][j];

            // Find the right piece of the vectors to work with.
            const int qx = rbn.getQIndex();
            const int ux = rbn.getUIndex();
            const int inpx  = transpose ? qx : ux;
            const int outpx = transpose ? ux : qx;

            // TODO: kludge: for now q-like output may have an unused element because
            // we always allocate the max space. Set the last element to zero in case
            // it doesn't get written.
            if (!transpose) outp[outpx + rbn.getMaxNQ()-1] = 0;

            rbn.multiplyByQBlock(sbState, useEulerAngles, &qp[qx], 
                                 transpose, &inp[inpx], &outp[outpx]);
        }
}

// u= Qinv * q or q = ~QInv * u
void SimbodyMatterSubsystemRep::multiplyByQMatrixInverse(const State& s, bool transpose, const Vector& in, Vector& out) const
{
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next()); // i.e., we must be *done* with Stage::Position

    assert(in.size() == (transpose?getTotalDOF():getTotalQAlloc()));
    out.resize(transpose?getTotalQAlloc():getTotalDOF());

    //TODO: this shouldn't be necessary
    assert(in.size()  < 2 || &in[1]  == &in[0] +1); // for now must be contiguous in memory
    assert(out.size() < 2 || &out[1] == &out[0]+1); 

    const Real* qp   = sbState.getQ();
    const Real* inp  = in.size()  ? &in[0]  : 0;
    Real*       outp = out.size() ? &out[0] : 0;

    const bool useEulerAngles = sbState.getModelVars().useEulerAngles;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& rbn = *rbNodeLevels[i][j];

            // Find the right piece of the vectors to work with.
            const int qx = rbn.getQIndex();
            const int ux = rbn.getUIndex();
            const int inpx  = transpose ? ux : qx;
            const int outpx = transpose ? qx : ux;

            // TODO: kludge: for now q-like output may have an unused element because
            // we always allocate the max space. Set the last element to zero in case
            // it doesn't get written.
            if (transpose) outp[outpx + rbn.getMaxNQ()-1] = 0;

            rbn.multiplyByQInvBlock(sbState, useEulerAngles, &qp[qx],
                                    transpose, &inp[inpx], &outp[outpx]);
        }
}


// Must be in ConfigurationStage to calculate qdot = Q*u.
void SimbodyMatterSubsystemRep::calcQDot(const State& s, const Vector& u, Vector& qdot) const {
    const SBModelVars&     mv = getModelVars(s);
    const Vector&          q  = getQ(s);
    const SBPositionCache& pc = getPositionCache(s);

    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(mv,q,pc, u, qdot);
}

// Must be in Stage::Velocity to calculate qdotdot = Qdot*u + Q*udot.
void SimbodyMatterSubsystemRep::calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
    const SBModelVars&     mv = getModelVars(s);
    const Vector&          q  = getQ(s);
    const SBPositionCache& pc = getPositionCache(s);
    const Vector&          u  = getU(s);

    assert(udot.size() == getTotalDOF());
    qdotdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDotDot(mv,q,pc,u, udot, qdotdot);
}

// State must be in Stage::Position.
void SimbodyMatterSubsystemRep::
calcMobilizerQDotFromU(const State& s, MobilizedBodyIndex mb, int nu, const Real* u, 
                       int nq, Real* qdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Position);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQInUse(sbState.getModelVars()));

    n.calcLocalQDotFromLocalU(sbState, u, qdot);
}

// State must be realized to Stage::Velocity, so that we can extract Q(q), QDot(q,u), and u from it to calculate
// qdotdot=Q(q)*udot + QDot(q,u)*u for this mobilizer.
void SimbodyMatterSubsystemRep::
calcMobilizerQDotDotFromUDot(const State& s, MobilizedBodyIndex mb, int nu, const Real* udot, 
                                  int nq, Real* qdotdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Velocity);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQInUse(sbState.getModelVars()));

    n.calcLocalQDotDotFromLocalUDot(sbState, udot, qdotdot);
}

// State must be realized through Stage::Instance. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of q's is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized coordinates.
// Returns X_FM(q).
Transform SimbodyMatterSubsystemRep::
calcMobilizerTransformFromQ(const State& s, MobilizedBodyIndex mb, int nq, const Real* q) const {
    const SBStateDigest sbState(s, *this, Stage::Instance);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nq == n.getNQInUse(sbState.getModelVars()));

    return n.calcMobilizerTransformFromQ(sbState, q);
}

// State must be realized through Stage::Position. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of u's is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized speeds.
// Returns V_FM(q,u)=H_FM(q)*u, where the q dependency is extracted from the State via
// the hinge transition matrix H_FM(q).
SpatialVec SimbodyMatterSubsystemRep::
calcMobilizerVelocityFromU(const State& s, MobilizedBodyIndex mb, int nu, const Real* u) const {
    const SBStateDigest sbState(s, *this, Stage::Position);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());

    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}

// State must be realized through Stage::Velocity. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of u's (and udot's) is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized accelerations.
// Returns A_FM(q,u,udot)=H_FM(q)*udot + HDot_FM(q,u)*u where the q and u dependencies
// are extracted from the State via H_FM(q), and HDot_FM(q,u).
SpatialVec SimbodyMatterSubsystemRep::
calcMobilizerAccelerationFromUDot(const State& s, MobilizedBodyIndex mb, int nu, const Real* udot) const{
    const SBStateDigest sbState(s, *this, Stage::Velocity);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());

    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}

// These perform the same computations as above but then transform the results so that they
// relate the child body's frame B to its parent body's frame P, rather than the M and F frames
// which are attached to B and P respectively but differ by a constant transform.
Transform SimbodyMatterSubsystemRep::
calcParentToChildTransformFromQ(const State& s, MobilizedBodyIndex mb, int nq, const Real* q) const {
    const Transform& X_PF = getMobilizerFrameOnParent(s,mb);
    const Transform& X_BM = getMobilizerFrame(s,mb);

    const Transform X_FM = calcMobilizerTransformFromQ(s,mb,nq,q);
    return X_PF * X_FM * ~X_BM; // X_PB
}

SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildVelocityFromU(const State& s, MobilizedBodyIndex mb, int nu, const Real* u) const {
    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}
SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildAccelerationFromUDot(const State& s, MobilizedBodyIndex mb, int nu, const Real* udot) const {
    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}




//
// We have V_GB = J u where J=~Phi*~H is the kinematic Jacobian (partial velocity matrix)
// that maps generalized speeds to spatial velocities. 
//
void SimbodyMatterSubsystemRep::calcSpatialKinematicsFromInternal(const State& s,
    const Vector&              v,
    Vector_<SpatialVec>&       Jv) const 
{
    const SBPositionCache& pc = getPositionCache(s);

    assert(v.size() == getTotalDOF());

    Jv.resize(getNBodies());

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcSpatialKinematicsFromInternal(pc,v, Jv);
        }
}

// If V is a spatial velocity, and you have an X=d(something)/dV (one per body)
// this routine will return d(something)/du for internal generalized speeds u. If
// instead you have d(something)/dR where R is a spatial configuration, this routine
// returns d(something)/dq PROVIDED that dq/dt = u for all q's. That's not true for
// quaternions, so be careful how you use this routine.
// In Kane's terminology, we are calculating the product of a (generalized)
// partial velocity with some vector.
void SimbodyMatterSubsystemRep::calcInternalGradientFromSpatial(const State& s, 
                                                    const Vector_<SpatialVec>& X,
                                                    Vector& JX) const
{
    assert(X.size() == getNBodies());

    const SBPositionCache& pc = getPositionCache(s);

    Vector_<SpatialVec> zTemp(getNBodies()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(pc, zTemp, X, JX);
        }
}

// This routine does the same thing as the above but accounts for centrifugal
// forces induced by velocities. The equivalent joint forces returned include
// both the applied forces and the centrifugal ones. Constraints are ignored.
void SimbodyMatterSubsystemRep::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const
{
    const SBPositionCache& pc = getPositionCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);

    assert(bodyForces.size() == getNBodies());
    mobilityForces.resize(getTotalDOF());

    Vector_<SpatialVec> allZ(getNBodies());

    // Don't do ground's level since ground has no inboard joint.
    for (int i=rbNodeLevels.size()-1 ; i>0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcEquivalentJointForces(pc,dc,
                bodyForces, allZ,
                mobilityForces);
        }
}

// Pass in a set of internal forces in T; we'll modify them here.
void SimbodyMatterSubsystemRep::calcConstraintCorrectedInternalForces(const State& s, Vector& T) {
    lConstraints->projectUVecOntoMotionConstraints(s, T);
}

bool SimbodyMatterSubsystemRep::getShowDefaultGeometry() const {
    return showDefaultGeometry;
}

void SimbodyMatterSubsystemRep::setShowDefaultGeometry(bool show) {
    showDefaultGeometry = show;
}

std::ostream& operator<<(std::ostream& o, const SimbodyMatterSubsystemRep& tree) {
    o << "SimbodyMatterSubsystemRep has " << tree.getNBodies() << " bodies (incl. G) in "
      << tree.rbNodeLevels.size() << " levels." << std::endl;
    o << "NodeNum->level,offset;stored nodeNum,level (stateOffset:dim)" << std::endl;
    for (int i=0; i < tree.getNBodies(); ++i) {
        o << i << "->" << tree.nodeNum2NodeMap[i].level << "," 
                       << tree.nodeNum2NodeMap[i].offset << ";";
        const RigidBodyNode& n = tree.getRigidBodyNode(i);
        o << n.getNodeNum() << "," << n.getLevel() 
          <<"(u"<< n.getUIndex()<<":"<<n.getDOF() 
          <<",q"<< n.getQIndex()<<":"<<n.getMaxNQ()<<")"<< std::endl;
    }

    return o;
}

    /////////////////////
    // SB STATE DIGEST //
    /////////////////////

void SBStateDigest::fillThroughStage(const SimbodyMatterSubsystemRep& matter, Stage g) {
    assert((int)g <= (int)matter.getStage(state) + 1);
    clear();

    if (g >= Stage::Model) {
        mv = &matter.getModelVars(state);
        mc = &matter.updModelCache(state);
    }
    if (g >= Stage::Instance) {
        if (mc->instanceVarsIndex >= 0)
            iv = &matter.getInstanceVars(state);
        if (mc->instanceCacheIndex >= 0)
            ic = &matter.updInstanceCache(state);
    }
    if (g >= Stage::Time) {
        if (mc->timeVarsIndex >= 0)
            tv = &matter.getTimeVars(state);
        if (mc->timeCacheIndex >= 0)
            tc = &matter.updTimeCache(state);
    }
    if (g >= Stage::Position) {
        if (matter.getNQ(state))
            q = &matter.getQ(state)[0];
        if (mc->qVarsIndex >= 0)
            pv = &matter.getPositionVars(state);

        if (matter.getNQErr(state))
            qErr = &matter.updQErr(state)[0];
        if (mc->qCacheIndex >= 0)
            pc = &matter.updPositionCache(state);
    }
    if (g >= Stage::Velocity) {
        if (matter.getNU(state))
            u = &matter.getU(state)[0];
        if (mc->uVarsIndex >= 0)
            vv = &matter.getVelocityVars(state);

        if (matter.getNQ(state))
            qdot = &matter.updQDot(state)[0];
        if (matter.getNUErr(state))
            uErr = &matter.updUErr(state)[0];
        if (mc->uCacheIndex >= 0)
            vc = &matter.updVelocityCache(state);
    }
    if (g >= Stage::Dynamics) {
        if (mc->dynamicsVarsIndex >= 0)
            dv = &matter.getDynamicsVars(state);
        if (mc->dynamicsCacheIndex >= 0)
            dc = &matter.updDynamicsCache(state);
    }
    if (g >= Stage::Acceleration) {
        if (mc->accelerationVarsIndex >= 0)
            av = &matter.getAccelerationVars(state);

        if (matter.getNU(state))
            udot = &matter.updUDot(state)[0];
        if (matter.getNQ(state))
            qdotdot = &matter.updQDotDot(state)[0];
        if (matter.getNUDotErr(state))
            udotErr = &matter.updUDotErr(state)[0];
        if (mc->accelerationCacheIndex >= 0)
            ac = &matter.updAccelerationCache(state);
    }

    stage = g;
}

