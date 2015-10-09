/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Derived from IVM code written by Charles Schwieters          *
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

/* Implementation of SimbodyMatterSubsystemRep. */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/ConditionalConstraint.h"

#include "SimbodyMatterSubsystemRep.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "MultibodySystemRep.h"
#include "MobilizedBodyImpl.h"
#include "ConstraintImpl.h"
#include "ContactEventActions.h"

using namespace SimTK;

#include <string>
#include <iostream>
using std::cout; using std::endl;

SimbodyMatterSubsystemRep::SimbodyMatterSubsystemRep
   (const SimbodyMatterSubsystemRep& src)
:   SimTK::Subsystem::Guts("SimbodyMatterSubsystemRep", "X.X.X")
{
    assert(!"SimbodyMatterSubsystemRep copy constructor ... TODO!");
}


void SimbodyMatterSubsystemRep::clearTopologyState() {
    // Unilateral constraints reference Constraints but not vice versa,
    // so delete the conditional constraints first.

    for (UnilateralContactIndex ucx(0); ucx < uniContacts.size(); ++ucx)
        delete uniContacts[ucx];
    uniContacts.clear();

    for (StateLimitedFrictionIndex fx(0); fx < stateLtdFriction.size(); ++fx)
        delete stateLtdFriction[fx];
    stateLtdFriction.clear();
    //TODO: more conditional constraints to come.

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

    nextAncestorConstrainedBodyPoolSlot = AncestorConstrainedBodyPoolIndex(0);

    DOFTotal=SqDOFTotal=maxNQTotal             = -1;
    nextQErrSlot = nextUErrSlot = nextMultSlot = 0; // TODO: are these being used?

    topologyCache.clear();
    topologyCacheIndex.invalidate();

    // New constraint fields (TODO not used yet)
    branches.clear();
    positionCoupledConstraints.clear();
    velocityCoupledConstraints.clear();
    accelerationCoupledConstraints.clear();
    dynamicallyCoupledConstraints.clear();

    // RigidBodyNodes themselves are owned by the MobilizedBodyImpls and will
    // be deleted when the MobilizedBodyImpl objects are.
    rbNodeLevels.clear();
    nodeNum2NodeMap.clear();

    m_showDefaultGeometry = true;
    m_useEulerAnglesByDefault = false;
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


UnilateralContactIndex SimbodyMatterSubsystemRep::
adoptUnilateralContact(UnilateralContact* child) {
    assert(child);
    invalidateSubsystemTopologyCache();

    const UnilateralContactIndex ucx(uniContacts.size());
    uniContacts.push_back(child); // grow

    // Tell the contact object its index within the matter subsystem.
    child->setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), ucx);
    return ucx;
}

StateLimitedFrictionIndex SimbodyMatterSubsystemRep::
adoptStateLimitedFriction(StateLimitedFriction* child) {
    assert(child);
    invalidateSubsystemTopologyCache();

    const StateLimitedFrictionIndex fx(stateLtdFriction.size());
    stateLtdFriction.push_back(child); // grow

    // Tell the friction object its index within the matter subsystem.
    child->setMyIndex(fx);
    return fx;
}


void SimbodyMatterSubsystemRep::createGroundBody() {
    assert(mobilizedBodies.empty());
    invalidateSubsystemTopologyCache();

    mobilizedBodies.push_back(new MobilizedBody::Ground());
    mobilizedBodies[GroundIndex]->updImpl()
        .setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), 
                              MobilizedBodyIndex(),    // no parent 
                              GroundIndex); //== MobilizedBodyIndex(0)
}

MobilizedBodyIndex SimbodyMatterSubsystemRep::getParent(MobilizedBodyIndex body) const { 
    return getRigidBodyNode(body).getParent()->getNodeNum();
}

Array_<MobilizedBodyIndex>
SimbodyMatterSubsystemRep::getChildren(MobilizedBodyIndex body) const {
    const RigidBodyNode& node = getRigidBodyNode(body);
    Array_<MobilizedBodyIndex> children;
    for (MobilizedBodyIndex bx(0); bx < node.getNumChildren(); ++bx)
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
SimbodyMatterSubsystemRep::getBodyTransform(const State& s, MobilizedBodyIndex body) const { 
    return getRigidBodyNode(body).getX_GB(getTreePositionCache(s));
}

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyVelocity(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getV_GB(getTreeVelocityCache(s));
}

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyAcceleration(const State& s, MobilizedBodyIndex body) const {
    return getRigidBodyNode(body).getA_GB(getTreeAccelerationCache(s));
}

const SpatialVec&
SimbodyMatterSubsystemRep::getMobilizerCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getMobilizerCoriolisAcceleration(getTreeVelocityCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getTotalCoriolisAcceleration(getTreeVelocityCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getGyroscopicForce(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getGyroscopicForce(getTreeVelocityCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getTotalCentrifugalForces(getTreeVelocityCache(s));
}
const SpatialVec& SimbodyMatterSubsystemRep::
getArticulatedBodyCentrifugalForces(const State& s, MobodIndex body) const {
  return getRigidBodyNode(body).getArticulatedBodyCentrifugalForces
                                        (getArticulatedBodyVelocityCache(s));
}



//==============================================================================
//                     ACQUIRE SYSTEM RESOURCES IMPL
//==============================================================================
// This method is invoked after the matter subsystem has been added to a
// MultibodySystem. We can now create system-level objects.
void SimbodyMatterSubsystemRep::acquireSystemResourcesImpl() {
    auto& mbs = updMultibodySystem();

    // Add impact and contact Actions.
    ImpactEvent& impact = mbs.updImpactEvent();
    ContactEvent& contact = mbs.updContactEvent();

    impact.adoptEventAction(new ImpactEventAction(*this));

    contact.adoptEventAction(new ContactEventAction(*this));
}


//==============================================================================
//                               REALIZE TOPOLOGY
//==============================================================================

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes we're going to need later for state
// variables and cache entries. We allocate and initialize all the
// Modeling variables here.
void SimbodyMatterSubsystemRep::endConstruction(State& s) {
    if (subsystemTopologyHasBeenRealized()) 
        return; // already done

    // This creates a RigidBodyNode owned by the the Topology cache of each 
    // MobilizedBody. Each RigidBodyNode lists as its parent the RigidBodyNode
    // contained in the MobilizedBody's parent. We simultaneously build up the 
    // computational version of the multibody tree, based on RigidBodyNode 
    // objects rather than on MobilizedBody objects.
    nodeNum2NodeMap.clear();
    rbNodeLevels.clear();
    DOFTotal = SqDOFTotal = maxNQTotal = 0;

    // state allocation
    nextUSlot   = UIndex(0);
    nextUSqSlot = USquaredIndex(0);
    nextQSlot   = QIndex(0);

    //Must do these in order from lowest number (ground) to highest. 
    for (MobilizedBodyIndex mbx(0); mbx<getNumMobilizedBodies(); ++mbx) {
        // Create the RigidBodyNode properly linked to its parent.
        const MobilizedBodyImpl& mbr = getMobilizedBody(mbx).getImpl();
        const RigidBodyNode& n = mbr.realizeTopology(s,nextUSlot,nextUSqSlot,nextQSlot);

        // Create the computational multibody tree data structures, organized 
        // by level.
        const int level = n.getLevel();
        if ((int)rbNodeLevels.size() <= level)
            rbNodeLevels.resize(level+1); // make room for the new level
        const int nodeIndexWithinLevel = rbNodeLevels[level].size();
        rbNodeLevels[level].push_back(&n);
        nodeNum2NodeMap.push_back(RigidBodyNodeId(level, nodeIndexWithinLevel));

        // Count up multibody tree totals.
        const int ndof = n.getDOF();
        DOFTotal += ndof; SqDOFTotal += ndof*ndof;
        maxNQTotal += n.getMaxNQ();
    }
    
    // Order doesn't matter for constraints as long as the bodies are already 
    // there. Quaternion normalization constraints exist only at the 
    // position level, however they are not topological since modeling
    // choices affect whether we use them. See realizeModel() below.

//    groundAncestorConstraint.clear();
//    branches.clear(); 
//    if (rbNodeLevels.size() > 1)
//        branches.resize(rbNodeLevels[1].size()); // each level 1 body is a branch

    nextAncestorConstrainedBodyPoolSlot = AncestorConstrainedBodyPoolIndex(0);

    for (ConstraintIndex cx(0); cx<getNumConstraints(); ++cx) {
        // Note: currently there is no such thing as a disabled constraint at
        // the Topology stage. Even a constraint which is disabled by default 
        // is not actually disabled until Model stage, and can be re-enabled 
        // at Model stage.
        const ConstraintImpl& crep = getConstraint(cx).getImpl();
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

    for (UnilateralContactIndex ucx(0); ucx < getNumUnilateralContacts(); ++ucx) {
        const UnilateralContact& uc = getUnilateralContact(ucx);
        uc.realizeTopology(s);
    }

}

int SimbodyMatterSubsystemRep::realizeSubsystemTopologyImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "SimbodyMatterSubsystem::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this 
    // call*. Afterwards they are truly const so we don't declare them mutable,
    // but cheat here instead.
    SimbodyMatterSubsystemRep* mThis = 
        const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!subsystemTopologyHasBeenRealized()) 
        mThis->endConstruction(s); // no more bodies after this!

    // Fill in the local copy of the topologyCache from the information
    // calculated in endConstruction(). Also ask the State for some room to
    // put Modeling variables & cache and remember the indices in our 
    // construction cache.
    SBTopologyCache& tc = mThis->topologyCache;

    tc.nBodies      = nodeNum2NodeMap.size();
    tc.nConstraints = constraints.size();
    tc.nParticles   = 0; // TODO
    tc.nAncestorConstrainedBodies = (int)nextAncestorConstrainedBodyPoolSlot;

    tc.nDOFs        = DOFTotal;
    tc.maxNQs       = maxNQTotal;
    tc.sumSqDOFs    = SqDOFTotal;

    SBModelVars mvars;
    mvars.allocate(topologyCache);
    setDefaultModelValues(topologyCache, m_useEulerAnglesByDefault, mvars);
    tc.modelingVarsIndex  = 
        allocateDiscreteVariable(s,Stage::Model, new Value<SBModelVars>(mvars));

    tc.modelingCacheIndex = 
        allocateCacheEntry(s,Stage::Model, new Value<SBModelCache>());

    SBInstanceVars iv;
    iv.allocate(topologyCache);
    setDefaultInstanceValues(mvars, iv); // sets lock-by-default, but not q or u
    tc.topoInstanceVarsIndex = 
        allocateDiscreteVariable(s, Stage::Instance, 
                                 new Value<SBInstanceVars>(iv));
    tc.instanceCacheIndex = 
        allocateCacheEntry(s, Stage::Instance, new Value<SBInstanceCache>());

    // Allocate the rest of the cache entries now although they won't get
    // any interesting content until later.

    tc.timeCacheIndex = 
        allocateCacheEntry(s, Stage::Time, new Value<SBTimeCache>());

    // Basic tree position kinematics can be calculated any time after Instance
    // stage and should be filled in first during realizePosition() and then
    // marked valid so later computations during the same realization can
    // access these quantities. This has an extra dependency on the q-version
    // so that it is invalidated whenever a q changes. This can be evaluated
    // earlier using realizePositionKinematics(). We guarantee this has been
    // realized by Stage::Position.
    tc.treePositionCacheIndex = s.allocateCacheEntryWithPrerequisites
       (getMySubsystemIndex(), Stage::Instance, Stage::Position,
        true /*q*/, false /*u*/, false /*z*/, {} /*dv*/, {} /*ce*/,
        new Value<SBTreePositionCache>());

    // Here is where later computations during realizePosition() go; these
    // will assume that the TreePositionCache is available. So you can 
    // calculate these prior to Position stage's completion but not until
    // the TreePositionCache has been marked valid.
    tc.constrainedPositionCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Time,
                               new Value<SBConstrainedPositionCache>());

    // Composite body inertias *can* be calculated any time after 
    // PositionKinematics are available, but they aren't ever needed internally
    // so we won't compute them unless explicitly requested.
    tc.compositeBodyInertiaCacheIndex = s.allocateCacheEntryWithPrerequisites
       (getMySubsystemIndex(), Stage::Instance, Stage::Infinity,
        false /*q*/, false /*u*/, false /*z*/, {} /*dv*/, 
        {CacheEntryKey(getMySubsystemIndex(), tc.treePositionCacheIndex)},
        new Value<SBCompositeBodyInertiaCache>());

    // Articulated body inertias *can* be calculated any time after 
    // PositionKinematics are available but we want to put them off until 
    // Acceleration stage if possible.
    tc.articulatedBodyInertiaCacheIndex = s.allocateCacheEntryWithPrerequisites
       (getMySubsystemIndex(), Stage::Instance, Stage::Acceleration,
        false /*q*/, false /*u*/, false /*z*/, {} /*dv*/, 
        {CacheEntryKey(getMySubsystemIndex(), tc.treePositionCacheIndex)},
        new Value<SBArticulatedBodyInertiaCache>());

    // Basic tree velocity kinematics can be calculated any time after Instance
    // stage, provided PositionKinematics have been realized, or unconditionally
    // after stage Position. These should be filled in first during 
    // realizeVelocity() and then marked valid so later computations during the 
    // same realization can access these quantities. Note that qdots are 
    // automatically allocated in the State's Velocity-stage cache. This has 
    // extra dependencies on the position kinematics cache entry and u state
    // variables, so it is invalidated whenever a q or u changes. This can be 
    // evaluated earlier using realizeVelocityKinematics(). We guarantee this 
    // has been realized by Stage::Velociy.
    tc.treeVelocityCacheIndex = s.allocateCacheEntryWithPrerequisites
       (getMySubsystemIndex(), Stage::Instance, Stage::Velocity,
        false /*q*/, true /*u*/, false /*z*/, {} /*dv*/,
        {CacheEntryKey(getMySubsystemIndex(), tc.treePositionCacheIndex)},
        new Value<SBTreeVelocityCache>());

    // Here is where later computations during realizeVelocity() go; these
    // will assume that the TreeVelocityCache is available. So you can 
    // calculate these prior to Velocity stage's completion but not until
    // the TreeVelocityCache has been marked valid.
    tc.constrainedVelocityCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Position,
                               new Value<SBConstrainedVelocityCache>());

    // Articulated body velocity calculations *can* be calculated any time after 
    // VelocityKinematics and articulated body inertias are available but we 
    // want to put them off until Acceleration stage if possible.
    tc.articulatedBodyVelocityCacheIndex = s.allocateCacheEntryWithPrerequisites
       (getMySubsystemIndex(), Stage::Instance, Stage::Acceleration,
        false /*q*/, false /*u*/, false /*z*/, {} /*dv*/, 
        {CacheEntryKey(getMySubsystemIndex(), 
                       tc.treeVelocityCacheIndex),
         CacheEntryKey(getMySubsystemIndex(), 
                       tc.articulatedBodyInertiaCacheIndex)},
        new Value<SBArticulatedBodyVelocityCache>());

    tc.dynamicsCacheIndex = 
        allocateCacheEntry(s, Stage::Dynamics, 
                           new Value<SBDynamicsCache>());

    // Tree acceleration kinematics can be calculated any time after Dynamics
    // stage and should be filled in first during realizeAcceleration() and then
    // marked valid so later computations during the same realization can
    // access these quantities.
    // Note that qdotdots, udots, zdots are automatically allocated by
    // the State when we advance the stage past Model.
    tc.treeAccelerationCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Dynamics,
                               new Value<SBTreeAccelerationCache>());

    // Here is where later computations during realizeAcceleration() go; these
    // will assume that the TreeAccelerationCache is available. So you can 
    // calculate these prior to Acceleration stage's completion but not until
    // the TreeAccelerationCache has been marked valid.
    tc.constrainedAccelerationCacheIndex =
        allocateLazyCacheEntry(s, Stage::Dynamics,
                               new Value<SBConstrainedAccelerationCache>());

    tc.valid = true;

    // Allocate a cache entry for the topologyCache, and save a copy there.
    mThis->topologyCacheIndex = 
        allocateCacheEntry(s,Stage::Topology, new Value<SBTopologyCache>(tc));
    return 0;
}



//==============================================================================
//                                 REALIZE MODEL
//==============================================================================
// Here we lock in modeling choices as conveyed by the values of Model-stage
// state variables which now all have values. These choices determine the number 
// and types of state variables we're going to use to represent the changeable 
// properties of this matter subsystem. This is the last realization stage at
// which we are given a writable State. That means all the state variables 
// we'll ever need must be allocated here (although cache entries can be added
// later since they are mutable.)
//
// Issues we'll settle here include:
//  - whether to use 4 quaternions or 3 Euler angles to represent orientation
//  - the default values of all the state variables
//  - the locked q or u values for lock-by-default mobilizers
//
// We allocate and fill in the Model-stage cache with information that can be
// calculated now that we have values for the Model-stage state variables.

int SimbodyMatterSubsystemRep::realizeSubsystemModelImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Topology, 
        "SimbodyMatterSubsystem::realizeModel()");

    SBStateDigest sbs(s, *this, Stage::Model);
    const SBModelVars& mv = sbs.getModelVars();

    // We're going to finish initializing InstanceVars below.
    SBInstanceVars& iv = updInstanceVars(s);

    // Get the Model-stage cache and make sure it has been allocated and 
    // initialized if needed. It is OK to hold a reference here because the 
    // discrete variables (and cache entries) in the State are stable, that is, 
    // they don't change location even if more variables are added.
    SBModelCache& mc = updModelCache(s);
    mc.clear(); // forget any previous modeling information
    mc.allocate(topologyCache);


        // MOBILIZED BODY MODELING

    // Count quaternions, and assign a "quaternion pool" index to each 
    // MobilizedBody that needs one, and allow mobilizers to reserve some
    // position- and velocity-cache space for their own purposes. We 
    // can't do this until Model stage because it is a Model stage variable 
    // which decides whether ball-like joints get quaternions or Euler angles.
    mc.totalNQInUse = mc.totalNUInUse = mc.totalNQuaternionsInUse = 0;
    mc.totalNQPoolInUse = 0;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j) {
            const RigidBodyNode& node  = *rbNodeLevels[i][j];
            SBModelPerMobodInfo& mbInfo = 
                mc.updMobodModelInfo(node.getNodeNum());

            // Assign q's.
            mbInfo.nQInUse     = node.getNQInUse(mv);
            //KLUDGE: currently the Q slots are assigned at topology stage.
            //mbInfo.firstQIndex = QIndex(mc.totalNQInUse);
            mbInfo.firstQIndex = node.getQIndex(); // TODO
            mc.totalNQInUse   += mbInfo.nQInUse;

            // Assign u's.
            mbInfo.nUInUse     = node.getNUInUse(mv);
            //KLUDGE: currently the U slots are assigned at topology stage.
            //mbInfo.firstUIndex = UIndex(mc.totalNUInUse);
            mbInfo.firstUIndex = node.getUIndex(); // TODO
            mc.totalNUInUse   += mbInfo.nUInUse;

            // Assign quaternion pool slot.
            if (node.isUsingQuaternion(sbs, mbInfo.startOfQuaternion)) {
                mbInfo.hasQuaternionInUse  = true;
                mbInfo.quaternionPoolIndex = 
                    QuaternionPoolIndex(mc.totalNQuaternionsInUse);
                mc.totalNQuaternionsInUse++;
            }

            // Assign misc. cache data slots for q's of this mobilizer.
            if ((mbInfo.nQPoolInUse=node.calcQPoolSize(mv)) != 0) {
                mbInfo.startInQPool = 
                    MobodQPoolIndex(mc.totalNQPoolInUse);
                mc.totalNQPoolInUse += mbInfo.nQPoolInUse;
            } else mbInfo.startInQPool.invalidate();
        }


    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeModel(sbs); 

        // CONSTRAINT MODELING

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeModel(s);

        // STATE RESOURCE ALLOCATION

    // Now allocate all remaining variables and cache entries and record the 
    // state resource index numbers in the ModelCache. Although we allocate 
    // all the resources now, we can only initialize those that depend only on
    // Model-stage variables; initialization of the rest will be performed at 
    // Instance stage.

    // SBInstanceVars were mostly allocated at topology stage but the lockedQs
    // and lockedUs values can't be initialized until now. We're setting all
    // of them although only the ones that are lock-by-default will actually
    // get used. Note that lockedUs does double duty since it holds both u
    // and udot values; the udot ones must be initialized to zero.

    iv.lockedQs.resize(maxNQTotal); 
    setDefaultPositionValues(mv, iv.lockedQs); // set locked q's to init values

    // MobilizedBodies provide default values for their q's. The number and
    // values of these can depend on modeling variables, which are already
    // set in the State. Don't do this for Ground, which has no q's.
    for (MobilizedBodyIndex mbx(1); mbx < mobilizedBodies.size(); ++mbx) {
        const MobilizedBodyImpl& mb = mobilizedBodies[mbx]->getImpl();
        mb.copyOutDefaultQ(s, iv.lockedQs);
    }

    // Velocity variables are just the generalized speeds u, which the State 
    // knows how to deal with. Zero is always a reasonable value for velocity,
    // so we'll initialize it here.
    iv.lockedUs.resize(DOFTotal); // also holds locked udots
    setDefaultVelocityValues(mv, iv.lockedUs); // set locked Us to initial value
    // We'll overwrite the locked udots to zero below after we've used the
    // current setting to initialize the u's.


    mc.timeVarsIndex = 
        allocateDiscreteVariable(s, Stage::Time, new Value<SBTimeVars>());


    // Position variables are just q's, which the State knows how to deal with. 

    // Initialize state's q values to the same values we put into lockedQs.
    mc.qIndex = allocateQ(s, iv.lockedQs);
    mc.qVarsIndex.invalidate(); // no position-stage vars other than q


    // Initialize state's u values to the same values we put into lockedUs.
    mc.uIndex = allocateU(s, iv.lockedUs); // set state to same initial value

    // Done with the default u's stored in iv.lockedUs; now overwrite any
    // that were locked by default at Motion::Acceleration; those must be zero.
    for (MobilizedBodyIndex mbx(1); mbx < mobilizedBodies.size(); ++mbx) {
        if (iv.mobilizerLockLevel[mbx] != Motion::Acceleration)
            continue;

        const SBModelPerMobodInfo& mbInfo = mc.updMobodModelInfo(mbx);
        const int    nu     = mbInfo.nUInUse;
        const UIndex uStart = mbInfo.firstUIndex;
        for (int i=0; i < nu; ++i)
            iv.lockedUs[UIndex(uStart+i)] = 0;
    }

    mc.uVarsIndex.invalidate(); // no velocity-stage vars other than u


    // no z's
    // Probably no dynamics-stage variables but we'll allocate anyway.
    SBDynamicsVars dvars;
    dvars.allocate(topologyCache);
    setDefaultDynamicsValues(mv, dvars);

    mc.dynamicsVarsIndex = 
        allocateDiscreteVariable(s, Stage::Dynamics, 
                                 new Value<SBDynamicsVars>(dvars));


    // No Acceleration variables that I know of. But we can go through the
    // charade here anyway.
    SBAccelerationVars rvars;
    rvars.allocate(topologyCache);
    setDefaultAccelerationValues(mv, rvars);
    mc.accelerationVarsIndex = 
        allocateDiscreteVariable(s, Stage::Acceleration, 
                                 new Value<SBAccelerationVars>(rvars));


    return 0;
}



//==============================================================================
//                               REALIZE INSTANCE
//==============================================================================
// Here we lock in parameterization of ("instantiate") the model, including
//  - how the motion of each mobilizer is to be treated
//  - the total number of constraint equations
//  - physical parameters like mass and geometry. 
// All cache entries should be allocated at this stage although they can be 
// written into at any stage.
//
// This is the last stage that doesn't change during time stepping, so it is 
// important to calculate as much as possible now to avoid unnecessary work 
// later. The Instance-stage cache is fully calculated and filled in here. Any 
// values in higher-level caches that can be calculated now should be.

int SimbodyMatterSubsystemRep::
realizeSubsystemInstanceImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "SimbodyMatterSubsystem::realizeInstance()");

    const SBModelCache&   mc = getModelCache(s);
    const SBInstanceVars& iv = getInstanceVars(s);

    // Get the Instance-stage cache and make sure it has been allocated and 
    // initialized if needed. We'll fill it in here and then allocate all
    // the rest of the cache entries after that.
    SBInstanceCache& ic = updInstanceCache(s);
    ic.allocate(topologyCache, mc);


    // MOBILIZED BODY INSTANCE
    // Here we need to instantiate the Body, Mobilizer, Motion, and the 
    // implementing RigidBodyNode.

    // Body mass properties are now available in the InstanceVars.
    ic.totalMass = iv.particleMasses.sum();
    for (MobilizedBodyIndex i(1); i<getNumBodies(); ++i) // not Ground!
        ic.totalMass += iv.bodyMassProperties[i].getMass();
    // TODO: central and principal inertias

    // Mobilizer geometry is now available in the InstanceVars.
    // TODO: reference configuration

    // Count position-, velocity-, and acceleration- prescribed motions 
    // generated by Motion objects associated with MobilizedBodies and allocate
    // pools to hold the associated values in the state cache. When position is
    // prescribed (by specifying q(t)), the corresponding qdot and qdotdot are 
    // also prescribed and we use them to set u and udot (via u=N^-1 qdot and 
    // udot = N^-1(qdotdot-NDot*u)). Each prescribed udot will have a 
    // corresponding force calculated, and other known udots (zero or discrete;
    // anything but free) will also need force slots although they don't get 
    // UDotPool slots.
    //
    // There is no built-in support in the State for these pools, so we 
    // allocate them in Simbody's cache entries at the appropriate stages. The 
    // prescribed q pool is in the TimeCache, prescribed u (dependent on the 
    // TreePositionCache) is written into the ConstrainedPositionCache, 
    // prescribed udots are written to the DynamicsCache,
    // and the prescribed forces tau are calculated at the same time as the 
    // free (non-prescribed) udots and are thus in the TreeAccelerationCache.
    //
    // NOTE: despite appearances here, each pool is in MobilizedBodyIndex order, 
    // meaning that the prescribed position, velocity, and acceleration, and the
    // other known udot entries, will be intermingled in the ForcePool rather than 
    // neatly lined up as I've drawn them.
    //
    //            --------------------
    //     QPool |       nPresQ       |      NOTE: not really ordered like this
    //            \------------------/ 
    //             \----------------/-------------
    //     UPool   |              nPresU          |
    //             |----------------|-------------|
    //             |------------------------------|--------------
    //  UDotPool   |                    nPresUDot                |
    //             |------------------------------|--------------|
    //             |---------------------------------------------|----------
    // ForcePool   |                       nPresForces                      |
    //              ---------------------------------------------|----------
    //
    // Note that there are no slots allocated for q's, u's, or udots that are 
    // known to be zero; only the explicitly prescribed ones get a slot. And 
    // udots known for any reason get a slot in the ForcePool.

    // Motion options for all mobilizers are now available in the InstanceVars.
    // We need to figure out what cache resources are required, including
    // slots for the generalized forces that implement prescribed motion.
    ic.presQ.clear();       ic.zeroQ.clear();       ic.freeQ.clear();
    ic.presU.clear();       ic.zeroU.clear();       ic.freeU.clear();
    ic.presUDot.clear();    ic.zeroUDot.clear();    ic.freeUDot.clear();
    ic.presForce.clear();
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx) {
        const MobilizedBody&       mobod        = getMobilizedBody(mbx);
        const SBModelPerMobodInfo& modelInfo    = mc.getMobodModelInfo(mbx);
        SBInstancePerMobodInfo&    instanceInfo = ic.updMobodInstanceInfo(mbx);

        const int nq = modelInfo.nQInUse;
        const int nu = modelInfo.nUInUse;

        const QIndex qx = modelInfo.firstQIndex;
        const UIndex ux = modelInfo.firstUIndex;

        instanceInfo.clear(); // all motion is tentatively free

        // Treat Ground or Weld as prescribed to zero.
        if (mbx == GroundIndex || nq==0) {
            instanceInfo.qMethod = instanceInfo.uMethod = 
                instanceInfo.udotMethod = Motion::Zero;
            continue;
        }

        // Not Ground or a Weld. If locked, use the lock level to determine
        // how motion is calculated. Otherwise, if this mobilizer has a
        // non-disabled Motion, let the Motion decide.

        const Motion::Level lockLevel = iv.mobilizerLockLevel[mbx];
        if (lockLevel != Motion::NoLevel) { // locked
            if (lockLevel == Motion::Acceleration) {
                // If all the udots are zero, the method is Motion::Zero;
                // otherwise, it is Motion::Prescribed.
                instanceInfo.udotMethod = Motion::Zero;
                for (int i=0; i<nu; ++i)
                    if (iv.lockedUs[UIndex(ux+i)] != 0) {
                        instanceInfo.udotMethod = Motion::Prescribed;
                        break;
                    }
            } else if (lockLevel == Motion::Velocity) {
                // If all the u's are zero, the method is Motion::Zero;
                // otherwise, it is Motion::Prescribed.
                instanceInfo.uMethod = Motion::Zero;
                for (int i=0; i<nu; ++i)
                    if (iv.lockedUs[UIndex(ux+i)] != 0) {
                        instanceInfo.uMethod = Motion::Prescribed;
                        break;
                    }
                instanceInfo.udotMethod = Motion::Zero;
            } else if (lockLevel == Motion::Position) {
                //TODO: Motion::Zero should mean identity transform; not the
                //same as q==0. Always use Prescribed for positions for now.
                instanceInfo.qMethod    = Motion::Prescribed;
                instanceInfo.uMethod    = Motion::Zero;
                instanceInfo.udotMethod = Motion::Zero;
            }
        } else if (mobod.hasMotion() && !iv.prescribedMotionIsDisabled[mbx]) {
            // Not locked, but has an active Motion.
            const Motion& motion = mobod.getMotion();
            motion.calcAllMethods(s, instanceInfo.qMethod, 
                                     instanceInfo.uMethod,
                                     instanceInfo.udotMethod);
        }

        // Assign q's to appropriate index vectors for convenient
        // manipulation. Prescribed q's also need pool allocation.
        switch(instanceInfo.qMethod) {
        case Motion::Prescribed:
            instanceInfo.firstPresQ = PresQPoolIndex(ic.presQ.size());
            for (int i=0; i < nq; ++i)
                ic.presQ.push_back(QIndex(qx+i));
            break;
        case Motion::Zero:
            for (int i=0; i < nq; ++i)
                ic.zeroQ.push_back(QIndex(qx+i));
            break;
        case Motion::Free:
            for (int i=0; i < nq; ++i)
                ic.freeQ.push_back(QIndex(qx+i));
            break;
        case Motion::Discrete:
        case Motion::Fast:
            break; // nothing to do for these here
        default:
            SimTK_ASSERT1_ALWAYS(!"qMethod",
                "SimbodyMatterSubsystemRep::realizeSubsystemInstanceImpl(): "
                "Unrecognized qMethod %d.", (int)instanceInfo.qMethod);
        }

        // Assign u's to appropriate index vectors for convenient
        // manipulation. Prescribed u's also need pool allocation.
        switch(instanceInfo.uMethod) {
        case Motion::Prescribed:
            instanceInfo.firstPresU = PresUPoolIndex(ic.presU.size());
            for (int i=0; i < nu; ++i)
                ic.presU.push_back(UIndex(ux+i));
            break;
        case Motion::Zero:
            for (int i=0; i < nu; ++i)
                ic.zeroU.push_back(UIndex(ux+i));
            break;
        case Motion::Free:
            for (int i=0; i < nu; ++i)
                ic.freeU.push_back(UIndex(ux+i));
            break;
        case Motion::Discrete:
        case Motion::Fast:
            break; // nothing to do for these here
        default:
            SimTK_ASSERT1_ALWAYS(!"uMethod",
                "SimbodyMatterSubsystemRep::realizeSubsystemInstanceImpl(): "
                "Unrecognized uMethod %d.", (int)instanceInfo.uMethod);
        }

        // Assign udots to appropriate index vectors for convenient
        // manipulation. Prescribed udots also need pool allocation.

        switch(instanceInfo.udotMethod) {
        case Motion::Prescribed:
            instanceInfo.firstPresUDot = PresUDotPoolIndex(ic.presUDot.size());
            for (int i=0; i < nu; ++i)
                ic.presUDot.push_back(UIndex(ux+i));
            break;
        case Motion::Zero:
            for (int i=0; i < nu; ++i)
                ic.zeroUDot.push_back(UIndex(ux+i));
            break;
        case Motion::Free:
            for (int i=0; i < nu; ++i)
                ic.freeUDot.push_back(UIndex(ux+i));
            break;
        case Motion::Discrete:
        case Motion::Fast:
            break; // nothing to do for these here
        default:
            SimTK_ASSERT1_ALWAYS(!"udotMethod",
                "SimbodyMatterSubsystemRep::realizeSubsystemInstanceImpl(): "
                "Unrecognized udotMethod %d.", (int)instanceInfo.udotMethod);
        }

        // Any non-Free udot needs a prescribed force (tau) slot.
        // Count mobilities that need a slot to hold the calculated force due
        // to a known udot, whether prescribed or known for some other reason.
        if (instanceInfo.udotMethod != Motion::Free) {
            instanceInfo.firstPresForce = 
                PresForcePoolIndex(ic.presForce.size());
            for (int i=0; i < nu; ++i)
                ic.presForce.push_back(UIndex(ux+i));
        }
    }

    
    // CONSTRAINT INSTANCE


    // Count position, velocity, and acceleration constraint equations
    // generated by each Constraint that has not been disabled. The State's 
    // QErr, UErr, UDotErr/Multiplier arrays are laid out like this:
    //
    //           ------------------- -------------
    //    QErr  |    mHolonomic     | mQuaternion |
    //           ------------------- -------------
    //           ------------------- -------------------
    //    UErr  |         "         |   mNonholonomic   |
    //           ------------------- -------------------
    //           ------------------- ------------------- ---------------------
    // UDotErr  |         "         |         "         |  mAccelerationOnly  |
    //           ------------------- ------------------- ---------------------
    //
    // Multipliers are allocated exactly as for UDotErr.
    //
    // Note that a Constraint with both holonomic and nonholonomic constraint 
    // equations will get two disjoint segments in UErr (and UDotErr).

    // Each Constraint's realizeInstance() method will add its contribution to
    // these. "InUse" here means we only add up contributions from enabled
    // constraints and ignore disabled ones.
    ic.totalNHolonomicConstraintEquationsInUse        = 0;
    ic.totalNNonholonomicConstraintEquationsInUse     = 0;
    ic.totalNAccelerationOnlyConstraintEquationsInUse = 0;

    ic.totalNConstrainedBodiesInUse = 0;
    ic.totalNConstrainedMobilizersInUse = 0;
    ic.totalNConstrainedQInUse = 0; // q,u from the constrained mobilizers
    ic.totalNConstrainedUInUse = 0; 


    // Build sets of kinematically coupled constraints. Kinematic coupling can 
    // be different at position, velocity, and acceleration levels, with only 
    // holonomic constraints included at the position level, 
    // holonomic+nonholonic at the velocity level, and
    // holonomic+nonholonomic+accelerationOnly coupled at the acceleration 
    // level.

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

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeInstance(s);


    // Quaternion errors are located after last holonomic constraint error; 
    // see diagram above.
    ic.firstQuaternionQErrSlot = ic.totalNHolonomicConstraintEquationsInUse; 

    // We'll store the the physical constraint errors, followed by the 
    // quaternion constraints.
    ic.qErrIndex = allocateQErr(s, ic.totalNHolonomicConstraintEquationsInUse 
                                 + mc.totalNQuaternionsInUse);

    // Only physical constraints exist at the velocity and acceleration levels; 
    // the quaternion normalization constraints are gone.
    ic.uErrIndex    = allocateUErr   (s,   ic.totalNHolonomicConstraintEquationsInUse
                                         + ic.totalNNonholonomicConstraintEquationsInUse);
    ic.udotErrIndex = allocateUDotErr(s,   ic.totalNHolonomicConstraintEquationsInUse
                                         + ic.totalNNonholonomicConstraintEquationsInUse
                                         + ic.totalNAccelerationOnlyConstraintEquationsInUse);

    // ALLOCATE REMAINING CACHE ENTRIES
    // Now that we know all the Instance-stage info, we can allocate (or 
    // reallocate) the rest of the cache entries.
    updTimeCache(s).allocate(topologyCache, mc, ic);
    updTreePositionCache(s).allocate(topologyCache, mc, ic);
    updConstrainedPositionCache(s).allocate(topologyCache, mc, ic);
    updCompositeBodyInertiaCache(s).allocate(topologyCache, mc, ic);
    updArticulatedBodyInertiaCache(s).allocate(topologyCache, mc, ic);
    updTreeVelocityCache(s).allocate(topologyCache, mc, ic);
    updConstrainedVelocityCache(s).allocate(topologyCache, mc, ic);
    updArticulatedBodyVelocityCache(s).allocate(topologyCache, mc, ic);
    updDynamicsCache(s).allocate(topologyCache, mc, ic);
    updTreeAccelerationCache(s).allocate(topologyCache, mc, ic);
    updConstrainedAccelerationCache(s).allocate(topologyCache, mc, ic);

    // Now let the implementing RigidBodyNodes do their realization.
    SBStateDigest stateDigest(s, *this, Stage::Instance);
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeInstance(stateDigest); 
    
    return 0;
}



//==============================================================================
//                                REALIZE TIME
//==============================================================================
int SimbodyMatterSubsystemRep::realizeSubsystemTimeImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "SimbodyMatterSubsystem::realizeTime()");

    const SBStateDigest stateDigest(s, *this, Stage::Time);

    // the multibody tree cannot have time dependence

    // Now realize the MobilizedBodies, which will realize prescribed positions
    // if there are any; those will go into the TimeCache.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeTime(stateDigest);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeTime(stateDigest);

    // We're done with the TimeCache now.
    markCacheValueRealized(s, topologyCache.timeCacheIndex);
    return 0;
}



//==============================================================================
//                             REALIZE POSITION
//==============================================================================
// The goals here are:
// (1) fill in the TreePositionCache and mark it valid
// (2) realize the remaining position dependencies (which may depend on
//     step (1)) with the result going into ConstrainedPositionCache & QErr
//
// In step (1) we take the q's from the State and sweep outward from Ground
// through the multibody tree, calculating all position kinematics. Note that
// we *do not* look at the prescribed q's in the TimeCache except to calculate
// errors; unless someone has invoked a solver to update the State q's from
// the prescribed values they will not be the same.
//
// In step (2) we calculate the matter subsystem's other position dependencies
// which are:
//      - prescribed velocities (u's) --> PresUPool
//      - position constraint errors  --> State's QErr array
// These calculations, which may be user-written, depend on values from step (1)
// being valid already, even though the stage won't yet have been advanced to
// stage Position. So we explicitly mark the TreePositionCache valid as soon as
// it is known.

int SimbodyMatterSubsystemRep::realizeSubsystemPositionImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "SimbodyMatterSubsystem::realizePosition()");

    // We promise that position kinematics is available after Stage::Position,
    // so force that calculation now if it hasn't already been done.
    realizePositionKinematics(s);

    // Set up StateDigest for calculating remaining position information.
    const SBStateDigest     stateDigest(s, *this, Stage::Position);
    const SBInstanceCache&  ic = stateDigest.getInstanceCache();

    // MobilizedBodies
    // This will include writing the prescribed velocities (u's) into
    // the PresUPool in the ConstrainedPositionCache.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizePosition(stateDigest);


    // Put position constraint equation errors in qErr
    Vector& qErr = stateDigest.updQErr();
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& pseg = cInfo.holoErrSegment;
        if (pseg.length) {
            Real* perrp = &qErr[pseg.offset];
            ArrayView_<Real> perr(perrp, perrp+pseg.length);
            constraints[cx]->getImpl().calcPositionErrorsFromState(s, perr);
        }
    }

    // Now we're done with the ConstrainedPositionCache.
    markCacheValueRealized(s, topologyCache.constrainedPositionCacheIndex);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizePosition(stateDigest);
    return 0;
}


//==============================================================================
//                       REALIZE POSITION KINEMATICS
//==============================================================================
void SimbodyMatterSubsystemRep::
realizePositionKinematics(const State& state) const {
    const CacheEntryIndex tpcx = topologyCache.treePositionCacheIndex;

    if (isCacheValueRealized(state, tpcx))
        return; // already realized

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Instance, 
        "SimbodyMatterSubsystem::realizePositionKinematics()");

    const SBStateDigest     stateDigest(state, *this, Stage::Time);
    const SBInstanceVars&   iv = stateDigest.getInstanceVars();
    SBTreePositionCache&    tpc = stateDigest.updTreePositionCache();

    // realize tree positions (kinematics)
    // This includes all local cross-mobilizer kinematics (M in F, B in P)
    // and all global kinematics relative to Ground (G).

    // Any body which is using quaternions should calculate the quaternion
    // constraint here and put it in the appropriate slot of qErr.
    // Set generalized coordinates: sweep from base to tips.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizePosition(stateDigest); 

    // Ask the constraints to calculate ancestor-relative kinematics (still 
    // goes in TreePositionCache).
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl()
                            .calcConstrainedBodyTransformInAncestor(iv, tpc);

    markCacheValueRealized(state, tpcx);
}

// Position kinematics is realized only if 
//  - we are currently at Stage::Position or later
//      OR
//  - position kinematics was previously realized since the most recent change 
//    to Stage::Instance and the generalized coordinates q.
bool SimbodyMatterSubsystemRep::
isPositionKinematicsRealized(const State&  state) const {
    const CacheEntryIndex tpcx = topologyCache.treePositionCacheIndex;
    return isCacheValueRealized(state, tpcx);
}

void SimbodyMatterSubsystemRep::
invalidatePositionKinematics(const State& state) const {
    // Position kinematics is assumed calculated at Position stage, regardless 
    // of stage versions in the cache entry.
    state.invalidateAllCacheAtOrAbove(Stage::Position);
    const CacheEntryIndex tpcx = topologyCache.treePositionCacheIndex;
    markCacheValueNotRealized(state, tpcx);
}

//==============================================================================
//                      REALIZE COMPOSITE BODY INERTIAS
//==============================================================================
void SimbodyMatterSubsystemRep::
realizeCompositeBodyInertias(const State& state) const {
    const CacheEntryIndex cbx = topologyCache.compositeBodyInertiaCacheIndex;

    if (isCacheValueRealized(state, cbx))
        return; // already realized

    // Composite body inertias have not been realized. We don't want to realize 
    // position kinematics implicitly here so we'll fail if that hasn't been 
    // done already (typically by realize(Position)).

    SimTK_ERRCHK_ALWAYS(isPositionKinematicsRealized(state), 
        "SimbodyMatterSubsystem::realizeCompositeBodyInertias()",
        "Composite body inertias cannot be realized unless the state has been "
        "realized to Stage::Position or realizePositionKinematics() has been "
        "called explicitly.");

    // OK, we have everything we'll need to compute CBIs.

    SBCompositeBodyInertiaCache& cbc = Value<SBCompositeBodyInertiaCache>::
                                        updDowncast(updCacheEntry(state, cbx));

    calcCompositeBodyInertias(state, cbc.compositeBodyInertia);
    markCacheValueRealized(state, cbx);
}

// Composite body inertias are realized only if
//  - they were last realized since any change to Stage::Instance,
//    PositionKinematics, and generalized coordinates q.
bool SimbodyMatterSubsystemRep::
isCompositeBodyInertiasRealized(const State& state) const {
    const CacheEntryIndex cbx = topologyCache.compositeBodyInertiaCacheIndex;
    return isCacheValueRealized(state, cbx);
}

void SimbodyMatterSubsystemRep::
invalidateCompositeBodyInertias(const State& state) const {
    // There is no stage at which CBIs are guaranteed to have been computed
    // so we don't need to invalidate a stage here.
    const CacheEntryIndex cbx = topologyCache.compositeBodyInertiaCacheIndex;
    markCacheValueNotRealized(state, cbx);
}



//==============================================================================
//                     REALIZE ARTICULATED BODY INERTIAS
//==============================================================================
void SimbodyMatterSubsystemRep::
realizeArticulatedBodyInertias(const State& state) const {
    const CacheEntryIndex abx = topologyCache.articulatedBodyInertiaCacheIndex;

    if (isCacheValueRealized(state, abx))
        return; // already realized

    // Articulated body inertias have not been realized. We don't want to 
    // realize position kinematics implicitly here so we'll fail if that hasn't 
    // been done already (typically by realize(Position)).

    SimTK_ERRCHK_ALWAYS(isPositionKinematicsRealized(state), 
    "SimbodyMatterSubsystem::realizeArticulatedBodyInertias()",
    "Articulated body inertias cannot be realized unless the state has been "
    "realized to Stage::Position or realizePositionKinematics() has been "
    "called explicitly.");

    // OK, we have everything we'll need to compute ABIs.

    const SBInstanceCache&          ic  = getInstanceCache(state);
    const SBTreePositionCache&      tpc = getTreePositionCache(state);
    SBArticulatedBodyInertiaCache&  abc = updArticulatedBodyInertiaCache(state);

    // tip-to-base sweep
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; --i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeArticulatedBodyInertiasInward(ic,tpc,abc);

    markCacheValueRealized(state, abx);
}

// Articulated body inertias are realized only if
//  - we're already at Stage::Acceleration
//          OR
//  - they were last realized since any change to Stage::Instance,
//    PositionKinematics, and generalized coordinates q.
bool SimbodyMatterSubsystemRep::
isArticulatedBodyInertiasRealized(const State& state) const {
    const CacheEntryIndex abx = 
        topologyCache.articulatedBodyInertiaCacheIndex;
    return isCacheValueRealized(state, abx);
}

void SimbodyMatterSubsystemRep::
invalidateArticulatedBodyInertias(const State& state) const {
    // ABIs are assumed calculated at Acceleration stage, regardless of the 
    // flag in the cache entry.
    state.invalidateAllCacheAtOrAbove(Stage::Acceleration);
    const CacheEntryIndex abx = topologyCache.articulatedBodyInertiaCacheIndex;
    markCacheValueNotRealized(state, abx);
}



//==============================================================================
//                               REALIZE VELOCITY
//==============================================================================
// The goals here are:
// (1) fill in the TreeVelocityCache and mark it valid
// (2) realize the remaining velocity dependencies (which may depend on
//     step (1)) with the result going into ConstrainedVelocityCache & UErr
//
// In step (1) we take the u's from the State and sweep outward from Ground
// through the multibody tree, calculating all velocity kinematics. Note that
// we *do not* look at the prescribed u's in the ConstrainedPositionCache except 
// to calculate errors; unless someone has invoked a solver to update the State 
// u's from the prescribed values they will not be the same.
//
// In step (2) we calculate the matter subsystem's other velocity dependencies
// which are:
//      - (no need to do prescribed accelerations yet)
//      - velocity constraint errors  --> State's UErr array
// These calculations, which may be user-written, depend on values from step (1)
// being valid already, even though the stage won't yet have been advanced to
// stage Velocity. So we explicitly mark the TreeVelocityCache valid as soon as
// it is known.
int SimbodyMatterSubsystemRep::realizeSubsystemVelocityImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "SimbodyMatterSubsystem::realizeVelocity()");

    // We promise that velocity kinematics is available after Stage::Velocity,
    // so force that calculation now if it hasn't already been done.
    realizeVelocityKinematics(s);

    const SBStateDigest stateDigest(s, *this, Stage::Velocity);
    const SBInstanceCache& ic = stateDigest.getInstanceCache();

    // MobilizedBodies
    // probably not much to do here
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeVelocity(stateDigest);

    // Put velocity constraint equation errors in uErr
    Vector& uErr = stateDigest.updUErr();
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);

        const Segment& holoseg    = cInfo.holoErrSegment; // for derivs of holo constraints
        const Segment& nonholoseg = cInfo.nonholoErrSegment; // includes holo+nonholo
        const int mHolo = holoseg.length, mNonholo = nonholoseg.length;
        if (mHolo) {
            Real* pverrp = &uErr[holoseg.offset];
            ArrayView_<Real> pverr(pverrp, pverrp+mHolo);
            constraints[cx]->getImpl().calcPositionDotErrorsFromState(s, pverr);
        }
        if (mNonholo) {
            Real* verrp = &uErr[ic.totalNHolonomicConstraintEquationsInUse 
                                + nonholoseg.offset];
            ArrayView_<Real> verr(verrp, verrp+mNonholo);
            constraints[cx]->getImpl().calcVelocityErrorsFromState(s, verr);
        }
    }

    // Now we're done with the ConstrainedVelocityCache.
    markCacheValueRealized(s, topologyCache.constrainedVelocityCacheIndex);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeVelocity(stateDigest);
    return 0;
}


//==============================================================================
//                       REALIZE VELOCITY KINEMATICS
//==============================================================================
void SimbodyMatterSubsystemRep::
realizeVelocityKinematics(const State& state) const {
    const CacheEntryIndex velx = topologyCache.treeVelocityCacheIndex;
    if (isCacheValueRealized(state, velx))
        return; // already realized

    // Velocity kinematics is not realized. We don't want to realize position
    // kinematics implicitly so we'll fail if that hasn't been done.

    SimTK_ERRCHK_ALWAYS(isPositionKinematicsRealized(state), 
        "SimbodyMatterSubsystem::realizeVelocityKinematics()",
        "Velocity kinematics cannot be realized unless the state has been "
        "realized to Stage::Position or realizePositionKinematics() has been "
        "called explicitly.");

    // Here we know that position kinematics has already been realized but
    // velocity kinematics is out of date with respect to the position
    // kinematics version number, or because a u changed since last computed.
    // We know this subsystem's stage is at least Instance since otherwise
    // position kinematics couldn't have been valid.

    const SBStateDigest stateDigest(state, *this, Stage::Time);
    const SBInstanceVars&      iv  = stateDigest.getInstanceVars();
    const SBTreePositionCache& tpc = stateDigest.getTreePositionCache();
    SBTreeVelocityCache&       tvc = stateDigest.updTreeVelocityCache();

    // realize tree velocity kinematics
    // This includes all local cross-mobilizer velocities (M in F, B in P)
    // and all global velocities relative to Ground (G). Also computes qdots.

    // Set generalized speeds: sweep from base to tips.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeVelocity(stateDigest); 

    // Ask the constraints to calculate ancestor-relative velocity kinematics 
    // (still goes in TreeVelocityCache).
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl()
            .calcConstrainedBodyVelocityInAncestor(iv, tpc, tvc);

    // Velocity cache is now up to date.
    markCacheValueRealized(state, velx);
}


// Velocity kinematics is realized only if
//  - we're at Stage::Velocity or later
//          OR
//  - velocity kinematics was last realized since any change to Stage::Instance,
//    PositionKinematics, and generalized speeds u.
bool SimbodyMatterSubsystemRep::
isVelocityKinematicsRealized(const State& state) const {
    const CacheEntryIndex velx = topologyCache.treeVelocityCacheIndex;
    return isCacheValueRealized(state, velx);
}

void SimbodyMatterSubsystemRep::
invalidateVelocityKinematics(const State& state) const {
    // Velocity kinematics is assumed calculated at Velocity stage, regardless 
    // of the flag in the cache entry.
    state.invalidateAllCacheAtOrAbove(Stage::Velocity);
    const CacheEntryIndex tvcx = topologyCache.treeVelocityCacheIndex;
    markCacheValueNotRealized(state, tvcx);
}




//==============================================================================
//                 REALIZE ARTICULATED BODY VELOCITY CACHE
//==============================================================================
void SimbodyMatterSubsystemRep::
realizeArticulatedBodyVelocity(const State& state) const {
    const CacheEntryIndex abvx = 
        topologyCache.articulatedBodyVelocityCacheIndex;

    if (isCacheValueRealized(state, abvx))
        return; // already realized

    // Articulated body velocity computations have not been realized. We don't 
    // want to realize velocity kinematics or articulated body inertias 
    // implicitly here so we'll fail if they haven't been realized already.
    // Although velocity kinematics is available after Stage::Velocity, usually
    // ABIs are not calculated until Stage::Acceleration.

    SimTK_ERRCHK_ALWAYS(isVelocityKinematicsRealized(state), 
    "SimbodyMatterSubsystem::realizeArticulatedBodyVelocity()",
    "Articulated body velocity calculations cannot be realized unless "
    "velocity kinematics have already been realized, either by realizing the "
    "state to Stage::Velocity or by calling realizeVelocityKinematics() "
    "explicitly.");

    SimTK_ERRCHK_ALWAYS(isArticulatedBodyInertiasRealized(state), 
    "SimbodyMatterSubsystem::realizeArticulatedBodyVelocity()",
    "Articulated body velocity calculations cannot be realized unless "
    "articulated body inertias have already been realized, implicitly "
    "or by calling realizeArticulatedBodyInertias() explicitly.");

    // OK, we have everything we'll need to compute ABIs.

    const SBTreePositionCache&  tpc = getTreePositionCache(state);
    const SBTreeVelocityCache&  tvc = getTreeVelocityCache(state);
    const SBArticulatedBodyInertiaCache& 
                                abc = getArticulatedBodyInertiaCache(state);
    SBArticulatedBodyVelocityCache& 
                                abvc = updArticulatedBodyVelocityCache(state);

    // Order doesn't matter for this calculation. Ground's entries are
    // precalculated so start at level 1.
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeArticulatedBodyVelocityCache
                                                            (tpc,tvc,abc,abvc);

    markCacheValueRealized(state, abvx);
}

// Articulated body velocity computations are realized only if
//  - we're already at Stage::Acceleration
//          OR
//  - they were last realized since any change to Stage::Instance,
//    PositionKinematics, VelocityKinematics, ArticulatedBodyInertias,
//    and generalized coordinates and speeds q and u.
bool SimbodyMatterSubsystemRep::
isArticulatedBodyVelocityRealized(const State& state) const {
    const CacheEntryIndex abvx = 
        topologyCache.articulatedBodyVelocityCacheIndex;
    return isCacheValueRealized(state, abvx);
}

void SimbodyMatterSubsystemRep::
invalidateArticulatedBodyVelocity(const State& state) const {
    // AB velocity computations are assumed done at Acceleration stage, 
    // regardless of the flag in the cache entry.
    state.invalidateAllCacheAtOrAbove(Stage::Acceleration);
    const CacheEntryIndex abvx = 
        topologyCache.articulatedBodyVelocityCacheIndex;
    markCacheValueNotRealized(state, abvx);
}


//==============================================================================
//                               REALIZE DYNAMICS
//==============================================================================
// Prepare for dynamics by calculating position-dependent quantities like the 
// articulated body inertias P, and velocity-dependent quantities like the 
// Coriolis acceleration.

int SimbodyMatterSubsystemRep::realizeSubsystemDynamicsImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "SimbodyMatterSubsystem::realizeDynamics()");

    SBStateDigest stateDigest(s, *this, Stage::Dynamics);

    // Get the Dynamics-stage cache; it was already allocated at Instance stage.
    SBDynamicsCache& dc = stateDigest.updDynamicsCache();

    // Probably nothing to do here.
    for (int i=0; i < (int)rbNodeLevels.size(); ++i)
        for (int j=0; j < (int)rbNodeLevels[i].size(); ++j)
            rbNodeLevels[i][j]->realizeDynamics(stateDigest);

    // MobilizedBodies
    // This will include writing the prescribed accelerations into
    // the udot array in the State.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeDynamics(stateDigest);

    // Realize Constraint dynamics.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeDynamics(stateDigest);

    return 0;
}



//==============================================================================
//                             REALIZE ACCELERATION
//==============================================================================
int SimbodyMatterSubsystemRep::realizeSubsystemAccelerationImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "SimbodyMatterSubsystem::realizeAcceleration()");

    // Articulated body inertias and velocity will be realized if necessary
    // when realizeLoopForwardDynamics() is called, fulfilling our promise that
    // those will be calculated by Stage::Accleration.

    // We ask our containing MultibodySystem for a reference to the cached 
    // forces accumulated from all the force subsystems. We use these to 
    // compute accelerations, with all results going into the AccelerationCache.
    const MultibodySystem& mbs = getMultibodySystem(); // owner of this subsystem
    realizeLoopForwardDynamics(s,
        mbs.getMobilityForces(s, Stage::Dynamics),
        mbs.getParticleForces(s, Stage::Dynamics),
        mbs.getRigidBodyForces(s, Stage::Dynamics));

    SBStateDigest stateDigest(s, *this, Stage::Acceleration);

    // MobilizedBodies
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeAcceleration(stateDigest);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeAcceleration(stateDigest);

    return 0;
}



//==============================================================================
//                               REALIZE REPORT
//==============================================================================
int SimbodyMatterSubsystemRep::realizeSubsystemReportImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "SimbodyMatterSubsystem::realizeReport()");

    // realize RB nodes report
    SBStateDigest stateDigest(s, *this, Stage::Report);
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeReport(stateDigest);

    // MobilizedBodies
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeReport(stateDigest);

    // realize Constraint report
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeReport(s);

    return 0;
}



//==============================================================================
//                    CALC DECORATIVE GEOMETRY AND APPEND
//==============================================================================
int SimbodyMatterSubsystemRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // Let the bodies and mobilizers have a chance to generate some geometry.
    for (MobilizedBodyIndex bx(0); bx<mobilizedBodies.size(); ++bx)
        mobilizedBodies[bx]->getImpl().calcDecorativeGeometryAndAppend(s,stage,geom);

    // Likewise for the constraints
    for (ConstraintIndex cx(0); cx<constraints.size(); ++cx)
        constraints[cx]->getImpl().calcDecorativeGeometryAndAppend(s,stage,geom);

    // Now add in any subsystem-level geometry.
    switch(stage) {
    case Stage::Topology: {
        assert(subsystemTopologyHasBeenRealized());
        // none yet
        break;
    }
    case Stage::Position: {
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

// We are in the process of realizeTopology() when we need to make this call.
// We pass in the partially-completed Topology-stage cache, which must have all
// the dimensions properly filled in at this point.
void SimbodyMatterSubsystemRep::
setDefaultModelValues(const SBTopologyCache& topologyCache,
                      bool useEulerAnglesByDefault,
                      SBModelVars& modelVars) const 
{
    // Tree-level defaults
    modelVars.useEulerAngles = useEulerAnglesByDefault;

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultModelValues(topologyCache, modelVars);

}

void SimbodyMatterSubsystemRep::setDefaultInstanceValues(const SBModelVars& mv, 
                                                         SBInstanceVars& iv) const 
{
    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultInstanceValues(mv, iv);

    assert((int)iv.mobilizerLockLevel.size() == getNumBodies());
    assert((int)iv.prescribedMotionIsDisabled.size() == getNumBodies());
    for (MobilizedBodyIndex mbx(0); mbx < getNumBodies(); ++mbx) {
        const MobilizedBody& mobod = getMobilizedBody(mbx);
        iv.mobilizerLockLevel[mbx] = mobod.getLockByDefaultLevel();
        iv.prescribedMotionIsDisabled[mbx] =
            mobod.hasMotion() ? mobod.getMotion().isDisabledByDefault()
                              : false;
    }

    assert((int)iv.constraintIsDisabled.size() == getNumConstraints());
    for (ConstraintIndex cx(0); cx < getNumConstraints(); ++cx)
        iv.constraintIsDisabled[cx] =  getConstraint(cx).isDisabledByDefault();

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
                                             SBAccelerationVars& accVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultAccelerationValues(mv, accVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::
setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}

void SimbodyMatterSubsystemRep::
setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disable) const {
    SBInstanceVars& instanceVars = updInstanceVars(s); // check/adjust stage
    instanceVars.constraintIsDisabled[constraint] = disable;   
}

bool SimbodyMatterSubsystemRep::getUseEulerAngles(const State& s) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.useEulerAngles;
}

bool SimbodyMatterSubsystemRep::isConstraintDisabled(const State& s, ConstraintIndex constraint) const {
    const SBInstanceVars& instanceVars = getInstanceVars(s); // check stage
    return instanceVars.constraintIsDisabled[constraint];
}

void SimbodyMatterSubsystemRep::
convertToEulerAngles(const State& inputState, State& outputState) const {
    outputState = inputState;
    if (!getUseEulerAngles(inputState)) {
        setUseEulerAngles(outputState, true);
        getMultibodySystem().realizeModel(outputState);
        const Vector& inputQ  = inputState.getQ();
        Vector&       outputQ = outputState.updQ();
        for (MobilizedBodyIndex i(0); i < getNumMobilizedBodies(); ++i) {
            const RigidBodyNode& node = getRigidBodyNode(i);
            node.convertToEulerAngles(inputQ, outputQ);
        }
    }
}

void SimbodyMatterSubsystemRep::
convertToQuaternions(const State& inputState, State& outputState) const {
    outputState = inputState;
    if (getUseEulerAngles(inputState)) {
        setUseEulerAngles(outputState, false);
        getMultibodySystem().realizeModel(outputState);
        const Vector& inputQ  = inputState.getQ();
        Vector&       outputQ = outputState.updQ();
        for (MobilizedBodyIndex i(0); i < getNumMobilizedBodies(); ++i) {
            const RigidBodyNode& node = getRigidBodyNode(i);
            node.convertToQuaternions(inputQ, outputQ);
        }
    }
}


bool SimbodyMatterSubsystemRep::
isUsingQuaternion(const State& s, MobilizedBodyIndex body) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    MobilizerQIndex startOfQuaternion; // we don't need this information here
    SBStateDigest sbs(s, *this, Stage::Model);
    return n.isUsingQuaternion(sbs, startOfQuaternion);
}

int SimbodyMatterSubsystemRep::getNumQuaternionsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.totalNQuaternionsInUse;
}

QuaternionPoolIndex SimbodyMatterSubsystemRep::getQuaternionPoolIndex
   (const State& s, MobilizedBodyIndex body) const
{
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.getMobodModelInfo(body).quaternionPoolIndex;
}

int SimbodyMatterSubsystemRep::getNumHolonomicConstraintEquationsInUse(const State& s) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    return ic.totalNHolonomicConstraintEquationsInUse;
}
int SimbodyMatterSubsystemRep::getNumNonholonomicConstraintEquationsInUse(const State& s) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    return ic.totalNNonholonomicConstraintEquationsInUse;
}
int SimbodyMatterSubsystemRep::getNumAccelerationOnlyConstraintEquationsInUse(const State& s) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    return ic.totalNAccelerationOnlyConstraintEquationsInUse;
}


const Array_<QIndex>& SimbodyMatterSubsystemRep::
getFreeQIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.freeQ;
}

const Array_<QIndex>& SimbodyMatterSubsystemRep::
getPresQIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.presQ;
}
const Array_<QIndex>& SimbodyMatterSubsystemRep::
getZeroQIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.zeroQ;
}

const Array_<UIndex>& SimbodyMatterSubsystemRep::
getFreeUIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.freeU;
}

const Array_<UIndex>& SimbodyMatterSubsystemRep::
getPresUIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.presU;
}
const Array_<UIndex>& SimbodyMatterSubsystemRep::
getZeroUIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.zeroU;
}

const Array_<UIndex>& SimbodyMatterSubsystemRep::
getFreeUDotIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.freeUDot;
}

const Array_<UIndex>& SimbodyMatterSubsystemRep::
getKnownUDotIndex(const State& state) const {
    const SBInstanceCache& ic = getInstanceCache(state);
    return ic.presForce;
}



//==============================================================================
//                      PACK/UNPACK FREE Q/U, ZERO KNOWN Q/U
//==============================================================================
void SimbodyMatterSubsystemRep::
packFreeQ(const State& s, const Vector& allQ, 
          Vector& packedFreeQ) const 
{
    const Array_<QIndex>& freeQX = getFreeQIndex(s);
    const int nq = getNQ(s);
    const int nfq = freeQX.size();
    assert(allQ.size() == nq);
    assert(packedFreeQ.size() == nfq);

    if (nfq==0) return;                         // all prescribed
    if (nfq==nq) {packedFreeQ=allQ; return;}    // none prescribed

    if (allQ.hasContiguousData() && packedFreeQ.hasContiguousData()) {
        const Real* allQp        = &allQ[0];
        Real*       packedFreeQp = &packedFreeQ[0];
        for (int i=0; i < nfq; ++i)
            packedFreeQp[i] = allQp[freeQX[i]];
        return;
    } 

    // Slower copy for noncontiguous data.
    for (int i=0; i < nfq; ++i)
        packedFreeQ[i] = allQ[freeQX[i]];
}

void SimbodyMatterSubsystemRep::
unpackFreeQ(const State& s, const Vector& packedFreeQ, 
            Vector& unpackedFreeQ) const 
{
    const Array_<QIndex>& freeQX = getFreeQIndex(s);
    const int nq = getNQ(s);
    const int nfq = freeQX.size();
    assert(packedFreeQ.size()   == nfq);
    assert(unpackedFreeQ.size() == nq);

    if (nfq==0) return;                                 // all prescribed
    if (nfq==nq) {unpackedFreeQ=packedFreeQ; return;}   // none prescribed

    if (packedFreeQ.hasContiguousData() && unpackedFreeQ.hasContiguousData()) {
        const Real* packedFreeQp   = &packedFreeQ[0];
        Real*       unpackedFreeQp = &unpackedFreeQ[0];
        for (int i=0; i < nfq; ++i)
            unpackedFreeQp[freeQX[i]] = packedFreeQp[i];
        return;
    } 

    // Slower copy for noncontiguous data.
    for (int i=0; i < nfq; ++i)
        unpackedFreeQ[freeQX[i]] = packedFreeQ[i];
}


void SimbodyMatterSubsystemRep::
zeroKnownQ(const State& s, Vector& qlike) const
{   // known = prescribed + zero
    const int nq  = getNQ(s);
    assert(qlike.size() == nq);

    const Array_<QIndex>& presQX = getPresQIndex(s);
    const Array_<QIndex>& zeroQX = getZeroQIndex(s);
    const int npq = presQX.size();
    const int nzq = zeroQX.size();

    if (npq==0 && nzq==0) return; // all q's are free

    if (qlike.hasContiguousData()) {
        Real* qp = &qlike[0];
        for (int i=0; i < npq; ++i)
            qp[presQX[i]] = 0;
        for (int i=0; i < nzq; ++i)
            qp[zeroQX[i]] = 0;
        return;
    } 

    // Slower zeroing for noncontiguous data.
    for (int i=0; i < npq; ++i)
        qlike[presQX[i]] = 0;
    for (int i=0; i < nzq; ++i)
        qlike[zeroQX[i]] = 0;
}

void SimbodyMatterSubsystemRep::
packFreeU(const State& s, const Vector& allU, 
          Vector& packedFreeU) const 
{
    const Array_<UIndex>& freeUX = getFreeUIndex(s);
    const int nu = getNU(s);
    const int nfu = freeUX.size();
    assert(allU.size() == nu);
    assert(packedFreeU.size() == nfu);

    if (nfu==0) return;                         // all prescribed
    if (nfu==nu) {packedFreeU=allU; return;}    // none prescribed

    if (allU.hasContiguousData() && packedFreeU.hasContiguousData()) {
        const Real* allUp        = &allU[0];
        Real*       packedFreeUp = &packedFreeU[0];
        for (int i=0; i < nfu; ++i)
            packedFreeUp[i] = allUp[freeUX[i]];
        return;
    } 

    // Slower copy for noncontiguous data.
    for (int i=0; i < nfu; ++i)
        packedFreeU[i] = allU[freeUX[i]];
}

void SimbodyMatterSubsystemRep::
unpackFreeU(const State& s, const Vector& packedFreeU, 
            Vector& unpackedFreeU) const 
{
    const Array_<UIndex>& freeUX = getFreeUIndex(s);
    const int nu = getNU(s);
    const int nfu = freeUX.size();
    assert(packedFreeU.size()   == nfu);
    assert(unpackedFreeU.size() == nu);

    if (nfu==0) return;                                 // all prescribed
    if (nfu==nu) {unpackedFreeU=packedFreeU; return;}   // none prescribed

    if (packedFreeU.hasContiguousData() && unpackedFreeU.hasContiguousData()) {
        const Real* packedFreeUp   = &packedFreeU[0];
        Real*       unpackedFreeUp = &unpackedFreeU[0];
        for (int i=0; i < nfu; ++i)
            unpackedFreeUp[freeUX[i]] = packedFreeUp[i];
        return;
    } 

    // Slower copy for noncontiguous data.
    for (int i=0; i < nfu; ++i)
        unpackedFreeU[freeUX[i]] = packedFreeU[i];
}


void SimbodyMatterSubsystemRep::
zeroKnownU(const State& s, Vector& ulike) const
{   // known = prescribed + zero
    const int nu  = getNU(s);
    assert(ulike.size() == nu);

    const Array_<UIndex>& presUX = getPresUIndex(s);
    const Array_<UIndex>& zeroUX = getZeroUIndex(s);
    const int npu = presUX.size();
    const int nzu = zeroUX.size();

    if (npu==0 && nzu==0) return; // all u's are free

    if (ulike.hasContiguousData()) {
        Real* up = &ulike[0];
        for (int i=0; i < npu; ++i)
            up[presUX[i]] = 0;
        for (int i=0; i < nzu; ++i)
            up[zeroUX[i]] = 0;
        return;
    } 

    // Slower zeroing for noncontiguous data.
    for (int i=0; i < npu; ++i)
        ulike[presUX[i]] = 0;
    for (int i=0; i < nzu; ++i)
        ulike[zeroUX[i]] = 0;
}

//==============================================================================
//                     FIND ENFORCED POS CONSTRAINT EQNS
//==============================================================================
/* Currently any enabled, unconditional constraint's holonomic constraint
equations are enforced (meaning they should be projected) as should any
active conditional constraint's holonomic equations.
*/
Array_<MultiplierIndex> SimbodyMatterSubsystemRep::
findEnforcedPosConstraintEqns(const State& state) const {
    Array_<MultiplierIndex> eqns;

    // Unconditional constraints.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        const ConstraintImpl& cons = constraints[cx]->getImpl();
        if (cons.isDisabled(state) || cons.isConditional())
            continue;
        int mp, mv, ma;
        cons.getNumConstraintEquationsInUse(state, mp, mv, ma);
        if (mp==0) continue;
        MultiplierIndex px0, vx0, ax0;
        cons.getIndexOfMultipliersInUse(state, px0, vx0, ax0);
        for (int i=0; i<mp; ++i)
            eqns.push_back(MultiplierIndex(px0+i));
    }

    // TODO: all conditional constraints (now just uni contact).
    for (UnilateralContactIndex ucx(0); ucx < uniContacts.size(); ++ucx) {
        const UnilateralContact& ucont = *uniContacts[ucx];
        if (ucont.getCondition(state) <= CondConstraint::Off)
            continue;
        eqns.push_back(ucont.getContactMultiplierIndex(state));
    }

    return eqns;
}


//==============================================================================
//                     FIND ENFORCED VEL CONSTRAINT EQNS
//==============================================================================
Array_<MultiplierIndex> SimbodyMatterSubsystemRep::
findEnforcedVelConstraintEqns(const State& state) const {
    Array_<MultiplierIndex> eqns;

    // Unconditional constraints.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        const ConstraintImpl& cons = constraints[cx]->getImpl();
        if (cons.isDisabled(state) || cons.isConditional())
            continue;
        int mp, mv, ma;
        cons.getNumConstraintEquationsInUse(state, mp, mv, ma);
        if (mv==0) continue;
        MultiplierIndex px0, vx0, ax0;
        cons.getIndexOfMultipliersInUse(state, px0, vx0, ax0);
        for (int i=0; i<mv; ++i)
            eqns.push_back(MultiplierIndex(vx0+i));
    }

    // TODO: all conditional constraints (now just uni contact).
    for (UnilateralContactIndex ucx(0); ucx < uniContacts.size(); ++ucx) {
        const UnilateralContact& ucont = *uniContacts[ucx];
        if (!ucont.hasFriction(state))
            continue;
        const CondConstraint::Condition cond = ucont.getCondition(state);
        if (cond != CondConstraint::Active) // i.e., rolling
            continue;
        MultiplierIndex ix_x, ix_y;
        ucont.getFrictionMultiplierIndices(state, ix_x, ix_y);
        eqns.push_back(ix_x); 
        eqns.push_back(ix_y);
    }

    return eqns;
}

//==============================================================================
//                         FIND ACTIVE MULTIPLIERS
//==============================================================================

// Find the multiplier index of each active constraint equation from among
// the enabled constraints, in ascending order.
// TODO: precalculate this when active set changes
Array_<MultiplierIndex> SimbodyMatterSubsystemRep::
findActiveMultipliers(const State& s) const {
    // These are the enabled constraint equations; they may not all be active.
    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = getNumAccelerationOnlyConstraintEquationsInUse(s);
    const int m        = mHolo+mNonholo+mAccOnly;

    // Assume all constraints are active.
    Array_<bool> isActive(m, true);
    int mActive = m; // may be reduced below

    //TODO: generalize this
    const int nUniContacts  = getNumUnilateralContacts();
    for (UnilateralContactIndex ux(0); ux < nUniContacts; ++ux) {
        const UnilateralContact& contact = getUnilateralContact(ux);
        if (contact.isActive(s))
            continue;
        isActive[contact.getContactMultiplierIndex(s)] = false;
        --mActive;
        if (contact.hasFriction(s)) {
            MultiplierIndex xIx, yIx;
            contact.getFrictionMultiplierIndices(s,xIx,yIx);
            isActive[xIx] = isActive[yIx] = false;
            mActive -= 2;
        }
    }

    assert(mActive >= 0);

    // Use the isActive array to enumerate the active constraint equations.
    Array_<MultiplierIndex> active;
    active.reserve(mActive);
    for (MultiplierIndex mx(0); mx < m; ++mx)
        if (isActive[mx]) active.push_back(mx);

    return active; // counting on return value optimization
}

// Output array will be sized to active.size() since it does not need to be
// pre-initialized.
void SimbodyMatterSubsystemRep::
packActiveMultipliers(const Array_<MultiplierIndex>& active,
                      const Vector& multipliers,
                      Vector& activeMultipliers) const
{
    const int mActive = (int)active.size();
    assert(multipliers.size() >= mActive);
    activeMultipliers.resize(mActive);

    for (int a=0; a < mActive; ++a)
        activeMultipliers[a] = multipliers[active[a]];
}

// Output array will be sized to active1.size()+active2.size() since it does 
// not need to be pre-initialized.
void SimbodyMatterSubsystemRep::
packActiveMultipliers2(const Array_<MultiplierIndex>& active1,
                       const Array_<MultiplierIndex>& active2,
                       const Vector& multipliers,
                       Vector& activeMultipliers) const
{
    const int mActive1 = active1.size(), mActive2 = active2.size();
    const int mActive = mActive1 + mActive2;
    assert(multipliers.size() >= mActive);
    activeMultipliers.resize(mActive);

    for (int a=0; a < mActive1; ++a)
        activeMultipliers[a] = multipliers[active1[a]];
    for (int a=0; a < mActive2; ++a)
        activeMultipliers[mActive1+a] = multipliers[active2[a]];
}

// Note that the output array must already be sized and initialized (usually
// to zero).
void SimbodyMatterSubsystemRep::
unpackActiveMultipliers(const Array_<MultiplierIndex>&  active,
                        const Vector&                   activeMultipliers,
                        Vector&                         multipliers) const 
{
    const int mActive = (int)active.size();
    assert(activeMultipliers.size() == mActive);
    assert(multipliers.size() >= mActive);

    for (int a=0; a < mActive; ++a)
        multipliers[active[a]] = activeMultipliers[a];
}

void SimbodyMatterSubsystemRep::
formActiveSubmatrix(const Array_<MultiplierIndex>&  active,
                    const Matrix&                   fullMat,
                    Matrix&                         subMat) const
{
    const int mActive = (int)active.size();
    assert(fullMat.nrow() == fullMat.ncol() && fullMat.nrow() >= mActive);

    subMat.resize(mActive, mActive);
    for (int col=0; col < mActive; ++col)
        for (int row=0; row < mActive; ++row)
            subMat(row,col) = fullMat(active[row],active[col]);
}


//==============================================================================
//                   OBSOLETE -- see calcPq()
//==============================================================================
void SimbodyMatterSubsystemRep::calcHolonomicConstraintMatrixPNInv(const State& s, Matrix& PNInv) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = getNQ(s);

    PNInv.resize(mp,nq);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into qErr and mHolo (mp)

        PNInv(holoSeg.offset, 0, holoSeg.length, nq) = 
            constraints[cx]->calcPositionConstraintMatrixPNInv(s);
    }
}

//==============================================================================
//                  CALC HOLONOMIC VELOCITY CONSTRAINT MATRIX P
//==============================================================================
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixP(const State& s, Matrix& P) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mHolo = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    P.resize(mHolo,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        P(holoSeg.offset, 0, holoSeg.length, nu) = 
            constraints[cx]->calcPositionConstraintMatrixP(s);
    }
}
//==============================================================================
//                 CALC HOLONOMIC VELOCITY CONSTRAINT MATRIX P^T
//==============================================================================
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixPt(const State& s, Matrix& Pt) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mHolo = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Pt.resize(nu,mHolo);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        // Fill in columns of Pt
        Pt(0, holoSeg.offset, nu, holoSeg.length) = 
            constraints[cx]->calcPositionConstraintMatrixPt(s);
    }
}

//==============================================================================
//                  CALC NONHOLONOMIC CONSTRAINT MATRIX V
//==============================================================================
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixV(const State& s, Matrix& V) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    V.resize(mNonholo,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        V(nonholoSeg.offset, 0, nonholoSeg.length, nu) = 
            constraints[cx]->calcVelocityConstraintMatrixV(s);
    }
}

//==============================================================================
//                  CALC NONHOLONOMIC CONSTRAINT MATRIX V^T
//==============================================================================
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixVt(const State& s, Matrix& Vt) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Vt.resize(nu,mNonholo);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        // Fill in columns of Vt
        Vt(0, nonholoSeg.offset, nu, nonholoSeg.length) = 
            constraints[cx]->calcVelocityConstraintMatrixVt(s);
    }
}
//==============================================================================
//                  CALC ACCELERATION ONLY CONSTRAINT MATRIX A
//==============================================================================
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixA (const State& s, Matrix& A) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    A.resize(mAccOnly,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        A(accOnlySeg.offset, 0, accOnlySeg.length, nu) = 
            constraints[cx]->calcAccelerationConstraintMatrixA(s);
    }
}
//==============================================================================
//                  CALC ACCELERATION ONLY CONSTRAINT MATRIX A^T
//==============================================================================
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixAt(const State& s, Matrix& At) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    At.resize(nu,mAccOnly);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstancePerConstraintInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        At(0, accOnlySeg.offset, nu, accOnlySeg.length) = 
            constraints[cx]->calcAccelerationConstraintMatrixAt(s);
    }
}



//==============================================================================
//                  CALC CONSTRAINT FORCES FROM MULTIPLIERS
//==============================================================================
// Must be realized to Stage::Velocity.
// Note that constraint forces have the opposite sign from applied forces,
// because we calculate multipliers from
//    M udot + ~G lambda = f_applied
// If you want to view the constraint generated forces as though they
// were applied forces, negate lambda before the call here. This method
// returns the individual constraint force contributions as well as returning
// the combined forces.
void SimbodyMatterSubsystemRep::
calcConstraintForcesFromMultipliers
   (const State&         s, 
    const Vector&        lambda,
    Vector_<SpatialVec>& bodyForcesInG,
    Vector&              mobilityForces,
    Array_<SpatialVec>&  consBodyForcesInG,
    Array_<Real>&        consMobilityForces) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int m        = mHolo+mNonholo+mAccOnly;

    assert(lambda.size() == m);

    // These must be zeroed because we'll be adding in constraint force
    // contributions and multiple constraints may generate forces on the
    // same body or mobility.
    bodyForcesInG.resize(getNumBodies()); bodyForcesInG.setToZero();
    mobilityForces.resize(getNU(s));      mobilityForces.setToZero();

    // These Arrays are for one constraint at a time.
    Array_<Real> lambdap, lambdav, lambdaa; // multipliers

    // Loop over all enabled constraints, ask them to generate forces, and
    // accumulate the results in the global problem return vectors.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const ConstraintImpl& crep = constraints[cx]->getImpl();

        // No heap allocation is being done here. These are views directly
        // into the proper segment of the longer array.
        ArrayView_<SpatialVec,ConstrainedBodyIndex> bodyF1_G = 
            crep.updConstrainedBodyForces(s, consBodyForcesInG);
        ArrayView_<Real,ConstrainedUIndex>          mobilityF1 = 
            crep.updConstrainedMobilityForces(s, consMobilityForces);

        const int ncb = bodyF1_G.size();
        const int ncu = mobilityF1.size();

        // These have to be zeroed because a Constraint force method is not
        // *required* to apply forces to all its bodies and mobilities.
        bodyF1_G.fill(SpatialVec(Vec3(0), Vec3(0)));
        mobilityF1.fill(Real(0));

        const SBInstancePerConstraintInfo& 
                              cInfo = ic.getConstraintInstanceInfo(cx);

        // Find this Constraint's multipliers within the global array.
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp=holoSeg.length, mv=nonholoSeg.length, 
                  ma=accOnlySeg.length;

        // Pack the multipliers into small arrays lambdap for holonomic 2nd 
        // derivs, labmdav for nonholonomic 1st derivs, and lambda for
        // acceleration-only.
        // Note: these lengths are *very* small integers!
        lambdap.resize(mp); lambdav.resize(mv); lambdaa.resize(ma);
        for (int i=0; i<mp; ++i) 
            lambdap[i] = lambda[                 holoSeg.offset    + i];
        for (int i=0; i<mv; ++i) 
            lambdav[i] = lambda[mHolo          + nonholoSeg.offset + i];
        for (int i=0; i<ma; ++i) 
            lambdaa[i] = lambda[mHolo+mNonholo + accOnlySeg.offset + i];

        // Generate forces for this Constraint. Body forces will come back
        // in the A frame; if that's not Ground then we have to re-express
        // them in Ground before moving on.
        crep.calcConstraintForcesFromMultipliers
                        (s, lambdap, lambdav, lambdaa, bodyF1_G, mobilityF1);
        if (crep.isAncestorDifferentFromGround()) {
            const Rotation& R_GA = 
                crep.getAncestorMobilizedBody().getBodyRotation(s);
            for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
                bodyF1_G[cbx] = R_GA*bodyF1_G[cbx];  // 30 flops
        }

        // Unpack constrained body forces and add them to the proper slots 
        // in the global body forces array. They are already expressed in
        // the Ground frame.
        for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
            bodyForcesInG[crep.getMobilizedBodyIndexOfConstrainedBody(cbx)] 
                += bodyF1_G[cbx];       // 6 flops

        // Unpack constrained mobility forces and add them into global array.
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux) 
            mobilityForces[cInfo.getUIndexFromConstrainedU(cux)] 
                += mobilityF1[cux];     // 1 flop
    }
}



//==============================================================================
//                         MULTIPLY BY PVA TRANSPOSE
//==============================================================================
// We have these constraint equations available:
// (1)  [Fp,fp] = pforces(t,q;   lambdap)       holonomic (position)
// (2)  [Fv,fv] = vforces(t,q,u; lambdav)       non-holonomic (velocity)
// (3)  [Fa,fa] = aforces(t,q,u; lambdaa)       acceleration-only
// where F denotes body spatial forces and f denotes u-space generalized forces.
//
// such that ~J*Fp+fp = ~P*lambdap
//           ~J*Fv+fv = ~V*lambdav
//       and ~J*Fa+fa = ~A*lambdaa
//
// with P=P(t,q), V=V(t,q,u), A=A(t,q,u). Note that P is the u-space matrix
// P=Dperrdot/Du, *not* the q-space matrix Pq=Dperr/Dq=P*N^-1. 
// (See multiplyByPqTranspose() to work conveniently with Pq.)
//
// Here we will use those equations to perform the multiplications by the
// matrices ~P,~V, and/or ~A times a multiplier-like vector: 
//                                    [lambdap]                       
// (4)  fu = ~G*lambda = [~P ~V ~A] * [lambdav] = ~J*(Fp+Fv+Fa) + (fp+fv+fa).
//                                    [lambdaa]
//
// Individual constraint force equations are calculated in constant time, so 
// the whole multiplication can be done in O(m+n) time where m=mp+mv+ma is the 
// total number of active constraint equations, and n the number of mobilities.
//
// In general the state must be realized through Velocity stage, but if the 
// system contains only holonomic constraints, or if only ~P is included, then 
// the result is only time- and position-dependent since we just need to use 
// equation (1). In that case we require only that the state be realized to 
// stage Position.
//
// All of the Vector arguments must use contiguous storage.
// Complexity is O(m+n).
void SimbodyMatterSubsystemRep::
multiplyByPVATranspose( const State&     s,
                        bool             includeP,
                        bool             includeV,
                        bool             includeA,
                        const Vector&    lambda,
                        Vector&          allfuVector) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;
    const int m = mHolo+mNonholo+mAccOnly;

    const int nu       = getNU(s);
    const int nb       = getNumBodies();

    assert(lambda.size() == m);
    assert(lambda.hasContiguousData());

    allfuVector.resize(nu);
    assert(allfuVector.hasContiguousData());
    if (nu==0) return;
    if (m==0) {allfuVector.setToZero(); return;}

    // Allocate a temporary body forces vector here. We'll map these to 
    // generalized forces as the penultimate step, then add those into 
    // the output argument allfuVector which will have already accumulated 
    // all directly-generated mobility forces.
    Vector_<SpatialVec> allF_GVector(nb);

    // We'll be accumulating constraint forces into these Vectors so zero 
    // them now. Multiple constraints may contribute to forces on the same 
    // body or mobility. 
    allF_GVector.setToZero();
    allfuVector.setToZero();

    // Overlay arrays on the contiguous Vector memory for faster elementwise
    // access. Any of the lambda segments may be empty.
    const Real* first = &lambda[0];
    ArrayViewConst_<Real>  allLambdap(first, first+mHolo);
    ArrayViewConst_<Real>  allLambdav(first+mHolo,
                                      first+mHolo+mNonholo);
    ArrayViewConst_<Real>  allLambdaa(first+mHolo+mNonholo, 
                                      first+mHolo+mNonholo+mAccOnly);

    ArrayView_<SpatialVec> allF_G(&allF_GVector[0], &allF_GVector[0] + nb);
    ArrayView_<Real>       allfu (&allfuVector[0],  &allfuVector[0]  + nu);

    // These Arrays are for one constraint at a time. We need separate 
    // memory for these because constrained bodies and constrained u's are
    // not ordered the same as the global ones, nor are they necessarily
    // contiguous in the global arrays. We're declaring these arrays
    // outside the loop to avoid heap allocation -- they will grow to the
    // max size needed by any constraint, then get resized as needed without
    // further heap allocation.
    Array_<SpatialVec,ConstrainedBodyIndex> oneF_G; // body spatial forces
    Array_<Real,      ConstrainedUIndex>    onefu;  // u-space generalized forces     
    Array_<Real,      ConstrainedQIndex>    onefq;  // q-space generalized forces     

    // Loop over all enabled constraints, ask them to generate forces, and
    // accumulate the results in the global problem arrays (allF_G,allfu).
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const ConstraintImpl& crep = constraints[cx]->getImpl();
        const SBInstancePerConstraintInfo& 
                              cInfo = ic.getConstraintInstanceInfo(cx);
        const int ncb = crep.getNumConstrainedBodies();
        const int ncu = cInfo.getNumConstrainedU();

        // These have to be zeroed because we're accumulating forces from
        // each of the three possible constraint levels below.
        oneF_G.resize(ncb);                onefu.resize(ncu);
        oneF_G.fill(SpatialVec(Vec3(0)));  onefu.fill(Real(0));

        // Find this Constraint's multiplier segments within the global array.
        // (These match the acceleration error segments.)
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp = includeP ? holoSeg.length    : 0;
        const int mv = includeV ? nonholoSeg.length : 0;
        const int ma = includeA ? accOnlySeg.length : 0;

        // Now generate forces. Body forces will come back in the A frame; 
        // if that's not Ground then we have to re-express them in Ground 
        // before moving on.
        if (mp) {
            const int ncq = cInfo.getNumConstrainedQ();
            onefq.resize(ncq); onefq.fill(Real(0));
            ArrayViewConst_<Real> lambdap(&allLambdap[holoSeg.offset],
                                          &allLambdap[holoSeg.offset]+mp);
            crep.addInPositionConstraintForces(s, lambdap, oneF_G, onefq);
            // OK just to write on onefu here because it is zero.
            crep.convertQForcesToUForces(s, onefq, onefu);
        }
        if (mv) {
            ArrayViewConst_<Real> lambdav(&allLambdav[nonholoSeg.offset],
                                          &allLambdav[nonholoSeg.offset]+mv);
            crep.addInVelocityConstraintForces(s, lambdav, oneF_G, onefu);                                       
        }
        if (ma) {
            ArrayViewConst_<Real> lambdaa(&allLambdaa[accOnlySeg.offset],
                                          &allLambdaa[accOnlySeg.offset]+ma);
            crep.addInAccelerationConstraintForces(s, lambdaa, oneF_G, onefu);                                       
        }

        // Fix expressed-in frame for body forces if necessary.
        if (crep.isAncestorDifferentFromGround()) {
            const Rotation& R_GA = 
                crep.getAncestorMobilizedBody().getBodyRotation(s);
            for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
                oneF_G[cbx] = R_GA*oneF_G[cbx];  // 30 flops
        }

        // Unpack constrained body forces and add them to the proper slots 
        // in the global body forces array.
        for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
            allF_G[crep.getMobilizedBodyIndexOfConstrainedBody(cbx)] 
                += oneF_G[cbx];       // 6 flops per constrained body

        // Unpack constrained mobility forces and add them into global array.
        // (1 flop per constrained mobility).
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux) 
            allfu[cInfo.getUIndexFromConstrainedU(cux)] += onefu[cux]; 
    }


    // Map the body forces into u-space generalized forces.
    // 12*nu + 18*nb flops.
    Vector ftmp(nu);
    multiplyBySystemJacobianTranspose(s, allF_GVector, ftmp);
    allfuVector += ftmp;
}



//==============================================================================
//                             CALC PVA TRANSPOSE
//==============================================================================
// Arranges for contiguous workspace if necessary, then makes repeated calls
// to multiplyByPVATranspose() to compute one column at a time of ~G.
// Complexity is O(m^2 + m*n) = O(m*n).
void SimbodyMatterSubsystemRep::
calcPVATranspose(   const State&     s,
                    bool             includeP,
                    bool             includeV,
                    bool             includeA,
                    Matrix&          PVAt) const
    {
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;
    const int m = mHolo+mNonholo+mAccOnly;

    const int nu = getNU(s);

    PVAt.resize(nu,m);
    if (m==0 || nu==0)
        return;

    Vector lambda(m, Real(0));

    // If PVAt's columns are contiguous we can avoid copying.
    const bool isContiguous = PVAt(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : m);
    for (int j=0; j < m; ++j) {
        lambda[j] = 1; // column we're working on
        if (isContiguous) {
            multiplyByPVATranspose(s, includeP, includeV, includeA, lambda, 
                                   PVAt(j));
        } else {
            multiplyByPVATranspose(s, includeP, includeV, includeA, lambda, 
                                   contig_col);
            PVAt(j) = contig_col;
        }
        lambda[j] = 0;
    }
}



//==============================================================================
//                          MULTIPLY BY Pq TRANSPOSE
//==============================================================================
// First we calculate the u-space result fu=~P*lambdap, then we map that
// result into q-space via fq=N^-T*fu.
// See multiplyByPVATranspose() above for an explanation.
// Vectors lambdap and fq must be using contiguous storage.
// Complexity is O(mp + n).
void SimbodyMatterSubsystemRep::
multiplyByPqTranspose(  const State&     state,
                        const Vector&    lambdap,
                        Vector&          fq) const
{
    Vector fu;
    // Calculate fu = ~P*lambdap.
    multiplyByPVATranspose(state, true, false, false, lambdap, fu);
    // Calculate fq = ~(N^-1) * fu = ~(~fu * N^-1)
    multiplyByNInv(state, true/*transpose*/,fu,fq);
}



//==============================================================================
//                             CALC Pq TRANSPOSE
//==============================================================================
// Arranges for contiguous workspace if necessary, then makes repeated calls
// to multiplyByPqTranspose() to compute one column at a time of ~Pq.
// Complexity is O(mp*mp + mp*n) = O(mp*n).
void SimbodyMatterSubsystemRep::
calcPqTranspose(const State& s, Matrix& Pqt) const {
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = getNQ(s);

    Pqt.resize(nq, mp);
    if (mp==0 || nq==0)
        return;

    Vector lambdap(mp, Real(0));

    // If Pqt's columns are contiguous we can avoid copying.
    const bool isContiguous = Pqt(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : nq);

    for (int i=0; i < mp; ++i) {
        lambdap[i] = 1; // column we're working on
        if (isContiguous) {
            multiplyByPqTranspose(s, lambdap, Pqt(i));
        } else {
            multiplyByPqTranspose(s, lambdap, contig_col);
            Pqt(i) = contig_col;
        }
        lambdap[i] = 0;
    }
}



//==============================================================================
//                 CALC WEIGHTED Pq_r TRANSPOSE (private)
//==============================================================================
/* The full Pq matrix for all mp enabled position constraint equatiosn and all 
nq generalized coordinates is mp X nq. We want the mfp X nfq submatrix of Pq 
that retains only 
  - the rows corresponding to the mfp *enForced* position constraint equatiosn, 
    and
  - the columns that correspond to the nfq *Free* (not prescribed) q's.
Call that submatrix Pq_r for "retained". Also, we want the retained rows scaled 
by 1/constraint tolerances and retained columns scaled by 1/q weights; call that
Pqw_r. And we're actually going to compute the nfq X mfp transpose ~Pqw_r, which
we'll call Pqw_rt.

  Pq   = P N\                   N\, Wq\ are pseudoinverses, Wu\ is inverse
  Pqw  = Tp Pq     Wq\ 
       = Tp Pq  (N Wu\ N\)
       = (Tp P Wu\) N\          because N\ N = I
       =     Pw     N\
  Pqwt = ~Pqw = ~N\ ~Pw         ~Pw = Wu\ ~P Tp because weights are symmetric

  Pqw_rt = Pqwt submatrix retaining only rows corresponding to free generalized
           coordinates and columns corresponding to enforced position 
           constraint equations.

We calculate one column at a time to avoid any matrix ops. We can do the column 
scaling of ~P by Tp for free, but the row scaling requires nq*mfp flops. Return 
matrix Pqw_rt must have contiguous-data columns. */
void SimbodyMatterSubsystemRep::
calcWeightedPqrTranspose( 
    const State&                    s,
    const Array_<MultiplierIndex>&  enforcedPosConsEqns,
    const Vector&                   Tp,   // 1/perr tols (mp)
    const Vector&                   ooWu, // 1/u weights (Wu\, nu)
    Matrix&                         Pqw_rt) const // nfq X mfp
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);
    assert(Tp.size() == mp && ooWu.size() == nu);

    const int nq = getNQ(s);
    const int nfq = ic.getTotalNumFreeQ();
    assert(nfq <= nq);
    const bool mustPack = (nfq != nq);

    const int mfp = (int)enforcedPosConsEqns.size();

    Pqw_rt.resize(nfq, mfp);
    if (mfp==0 || nfq==0)
        return;

    assert(Pqw_rt(0).hasContiguousData());

    Vector Ptcol(nu), PNInvtcol;
    Vector lambdap(mp, Real(0)); // full size

    int nxtj = 0; // next available enforced constraint column
    for (auto j : enforcedPosConsEqns) {
        lambdap[j] = Tp[j]; // this gives column j of ~P scaled by Tp[j]
        multiplyByPVATranspose(s, true,false,false, lambdap, Ptcol);
        lambdap[j] = 0;
        Ptcol.rowScaleInPlace(ooWu); // now (Wu\ ~P Tp)_i
        // Calculate (~Pqw)(j) = ~N\ (~Pw)(j)
        if (mustPack) {
            multiplyByNInv(s, true/*transpose*/,Ptcol,PNInvtcol);
            packFreeQ(s, PNInvtcol, Pqw_rt(nxtj));
        } else
            multiplyByNInv(s, true/*transpose*/,Ptcol,Pqw_rt(nxtj));
        ++nxtj;
    }
}



//==============================================================================
//                  CALC WEIGHTED PV_r TRANSPOSE (private)
//==============================================================================
/* The full P;V matrix PV is mpv X nu, where mpv=(mp+mv), the number of enabled
position and velocity constraint equations. Here we want the mfpv X nfu 
submatrix that retains only the rows corresponding to the *enForced* position
and velocity constraint equations, and retains only the columns that correspond 
to *Free* (not prescribed) generalized speeds u; call that PV_r (for 
"retained"). Also, we want the retained rows scaled by 1/constraint tolerances 
and retained columns scaled by 1/u weights; call that PVw_r. And we're actually 
going to compute the nfu X mfpv transpose ~PVw_r, which we'll call PVw_rt.

  PV   = [P]
         [V]
  PVw  = Tpv PV Wu\             Wu\ is Wu^-1
  PVwt = ~PVw = Wu\ ~PV Tpv     because weights are symmetric

  PVw_rt = PVwt submatrix retaining only rows corresponding to free generalized
           speeds and columns corresponding to enforced position and velocity 
           constraint equations.

We calculate one column at a time to avoid any matrix ops. We can do the
column scaling of ~PV by Tpv for free, but the row scaling requires nu*mfpv 
flops. Return matrix PVw_rt must have contiguous-data columns. */
void SimbodyMatterSubsystemRep::
calcWeightedPVrTranspose(
    const State&                    s,
    const Array_<MultiplierIndex>&  enforcedPosConsEqns,
    const Array_<MultiplierIndex>&  enforcedVelConsEqns,
    const Vector&                   Tpv,  // 1/verr tols (mp+mv)
    const Vector&                   ooWu, // 1/u weights (Wu\, nu)
    Matrix&                         PVw_rt) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int mv = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mpv = mp+mv;
    const int nu = getNU(s);
    assert(Tpv.size() == mpv && ooWu.size() == nu);

    const int nfu = ic.getTotalNumFreeU();
    assert(nfu <= nu);
    const bool mustPack = (nfu != nu);


    const int mfp = (int)enforcedPosConsEqns.size();
    const int mfv = (int)enforcedVelConsEqns.size();
    const int mfpv = mfp + mfv;

    PVw_rt.resize(nfu,mfpv);
    if (mfpv==0 || nfu==0)
        return;
    assert(PVw_rt.hasContiguousData());

    Vector lambdapv(mpv, Real(0)); // full size

    int nxtj = 0; // next available enforced constraint column
    if (mustPack) {
        Vector PVtcol(nu);
        for (auto list : {&enforcedPosConsEqns, &enforcedVelConsEqns})
            for (auto j : *list) {
                lambdapv[j] = Tpv[j]; // gives column j of ~PV scaled by Tpv[j]
                multiplyByPVATranspose(s, true,true,false, lambdapv, PVtcol);
                lambdapv[j] = 0;
                PVtcol.rowScaleInPlace(ooWu); // now (Wu\ ~PV Tpv)_i
                packFreeU(s, PVtcol, PVw_rt(nxtj));
                ++nxtj;
            }
    } else { // No prescribed motion so we can work in place.
        for (auto list : {&enforcedPosConsEqns, &enforcedVelConsEqns})
            for (auto j : *list) {
                VectorView PVw_rt_j = PVw_rt(nxtj);
                lambdapv[j] = Tpv[j]; // gives column j of ~PV scaled by Tpv[j]
                multiplyByPVATranspose(s, true,true,false, lambdapv, PVw_rt_j);
                lambdapv[j] = 0;
                PVw_rt_j.rowScaleInPlace(ooWu); // now (Wu\ ~PV Tpv)_i
                ++nxtj;
            }
    }
}



//==============================================================================
//                      CALC BIAS FOR MULTIPLY BY PVA
//==============================================================================
// We have these constraint equations available:
// (1)  pverr = Pq*qdot - Pt            (Pq==P*N^+, Pt==c(t,q))
// (2)  vaerr = V*udot  - b_v(t,q,u)
// (3)  aerr  = A*udot  - b_a(t,q,u)
// with P=P(t,q), N=N(q), V=V(t,q,u), A=A(t,q,u). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(m) time where m is the total number of constraint equations.
//
// Our plan is to use those equations to perform the multiplications by the
// matrices P,V, and/or A times a u-like vector. But first we need to calculate
// those extra terms that don't involve the matrices so we can subtract them
// off later. So we're going to compute:
// (4)  bias=[  bias_p,   bias_v,      bias_a    ]
//          =[   -Pt,   -b_v(t,q,u), -b_a(t,q,u) ].
// which we can get by using equations (1)-(3) with zero qdot or udot.
//
// In general the state must be realized through Velocity stage, but if the 
// system contains only holonomic constraints, or if only P is requested, then 
// bias is just bias_p and only time- and position-dependent since we just need
// to use eq. (1). In that case we require only stage Position.
//
// Note that bias_p is correct for use when multiplying a u-like vector by P
// or a q-like vector by Pq (==P*N^-1).
//
// The output vector must use contiguous storage.
// Complexity is O(m).
void SimbodyMatterSubsystemRep::
calcBiasForMultiplyByPVA(const State& s,
                         bool         includeP,
                         bool         includeV,
                         bool         includeA,
                         Vector&      bias) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;
    const int m = mHolo+mNonholo+mAccOnly;

    bias.resize(m);
    if (m == 0) return;

    assert(bias.hasContiguousData());

    // Overlay this Array on the bias Vector's data so that we can manipulate
    // small chunks of it repeatedly with no heap activity or virtual method
    // calls.
    ArrayView_<Real> biasArray(&bias[0], &bias[0] + m);

    // Except for holonomic constraint equations where we can work at the
    // velocity level, we'll need to supply body accelerations to the constraint
    // acceleration error routines. Because the udots are zero, the body 
    // accelerations include velocity-dependent terms only, i.e. the coriolis 
    // accelerations. Those have already been calculated in the state, but they
    // are AC_GB, the coriolis accelerations in Ground. The constraint methods
    // want those relative to their Ancestor frames, which might not be Ground.
    const Array_<SpatialVec,MobilizedBodyIndex>* allAC_GB = 0;
    if (mNonholo || mAccOnly)
        allAC_GB = &getTreeVelocityCache(s).totalCoriolisAcceleration;

    // This array will be resized and filled with the Ancestor-relative
    // coriolis accelerations for the constrained bodies of each velocity
    // or acceleration-only Constraint in turn; we're declaring it outside the 
    // loop to minimize heap allocation (resizing down doesn't normally free 
    // heap space). This won't be used if we have only holonomic constraints.
    Array_<SpatialVec,ConstrainedBodyIndex> AC_AB;

    // Subarrays of these all-zero arrays will be used to supply zero body
    // velocities and qdots (holonomic) or zero udots (nonholonomic and
    // acceleration-only) for each Constraint in turn; we're declaring 
    // them outside the loop to minimize heap allocation. They'll grow until 
    // they hit the maximum size needed by any Constraint.
    Array_<SpatialVec,ConstrainedBodyIndex> zeroV_AB;
    Array_<Real,      ConstrainedQIndex>    zeroQDot;
    Array_<Real,      ConstrainedUIndex>    zeroUDot;

    // Loop over all enabled constraints, ask them to generate constraint
    // errors, and collect those in the output bias vector.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        // Find this Constraint's err segments within the global array.
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp = includeP ? holoSeg.length    : 0;
        const int mv = includeV ? nonholoSeg.length : 0;
        const int ma = includeA ? accOnlySeg.length : 0;

        const ConstraintImpl& crep = constraints[cx]->getImpl();
        const int ncb = crep.getNumConstrainedBodies();

        if (mp) { // holonomic -- use velocity equations
            const int ncq = cInfo.getNumConstrainedQ();
            // Make sure we have enough zeroes.
            if (zeroV_AB.size() < ncb) zeroV_AB.resize(ncb, SpatialVec(Vec3(0)));
            if (zeroQDot.size() < ncq) zeroQDot.resize(ncq, Real(0));

            // Make subarrays; this does not require heap allocation.
            const ArrayViewConst_<SpatialVec,ConstrainedBodyIndex> 
                V0_AB = zeroV_AB(ConstrainedBodyIndex(0), ncb);
            const ArrayViewConst_<Real,ConstrainedQIndex>    
                qdot0 = zeroQDot(ConstrainedQIndex(0), ncq);

            // The holonomic error slots start at beginning of bias array.
            ArrayView_<Real> pverr = biasArray(holoSeg.offset, mp);

            // Write errors into pverr, which is a segment of the bias argument.
            crep.calcPositionDotErrors(s, V0_AB, qdot0, pverr);
        }

        if (!(mv || ma))
            continue; // nothing else to do here

        const int ncu = cInfo.getNumConstrainedU();
        // Make sure we have enough zeroes for udots.
        if (zeroUDot.size() < ncu) zeroUDot.resize(ncu, Real(0));
        // Make a subarray of the right size.
        const ArrayViewConst_<Real,ConstrainedUIndex>    
            udot0 = zeroUDot(ConstrainedUIndex(0), ncu);

        // Now fill in coriolis accelerations. If the Ancestor is Ground
        // we're just reordering. If it isn't Ground we have to transform
        // the coriolis accelerations from Ground to Ancestor, at a cost
        // of 105 flops/constrained body (not just re-expressing).
        crep.convertBodyAccelToConstrainedBodyAccel(s, *allAC_GB, AC_AB);

        // At this point AC_AB holds the coriolis accelerations of each
        // constrained body in A.

        if (mv) {   // non-holonomic constraints
            // The error slots begin after skipping the holonomic part of
            // the bias array. (That could be empty if P wasn't included.)
            const int start = mHolo+nonholoSeg.offset;
            ArrayView_<Real> vaerr = biasArray(start, mv);
            crep.calcVelocityDotErrors(s, AC_AB, udot0, vaerr);
        }
        if (ma) {   // acceleration-only constraints
            // The error slots begin after skipping the holonomic and 
            // non-holonomic parts of the bias array (those could be empty
            // if P or V weren't included).
            const int start = mHolo+mNonholo+accOnlySeg.offset;
            ArrayView_<Real> aerr = biasArray(start, ma);
            crep.calcAccelerationErrors(s, AC_AB, udot0, aerr);
        }
    }
}





//==============================================================================
//                    CALC BIAS FOR ACCELERATION CONSTRAINTS
//==============================================================================
// We have these constraint equations available:
// (1)  paerr = Pq*qdotdot - b_p(t,q,u)            (Pq==P*N^+, Pt==c(t,q))
// (2)  vaerr =  V*udot    - b_v(t,q,u)
// (3)  aerr  =  A*udot    - b_a(t,q,u)
// with P=P(t,q), N=N(q), V=V(t,q,u), A=A(t,q,u). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(m) time where m is the total number of constraint equations.
//
// We want to calculate those extra terms that don't involve the matrices so 
// we can subtract them off later. So we're going to compute:
// (4)  bias=[  bias_p,      bias_v,      bias_a    ]
//          =[-b_p(t,q,u), -b_v(t,q,u), -b_a(t,q,u) ].
// which we can get by using equations (1)-(3) with zero qdotdot or udot.
//
// The state must be realized through Velocity stage.
//
// Note that bias_p is correct for use when multiplying a u-like vector by P
// or a q-like vector by Pq (==P*N^-1).
//
// The output vector must use contiguous storage.
// Complexity is O(m).
void SimbodyMatterSubsystemRep::
calcBiasForAccelerationConstraints(const State& s,
                                   bool         includeP,
                                   bool         includeV,
                                   bool         includeA,
                                   Vector&      bias) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;
    const int m = mHolo+mNonholo+mAccOnly;

    bias.resize(m);
    if (m == 0) return;

    assert(bias.hasContiguousData());

    // Overlay this Array on the bias Vector's data so that we can manipulate
    // small chunks of it repeatedly with no heap activity or virtual method
    // calls.
    ArrayView_<Real> biasArray(&bias[0], &bias[0] + m);

    // We'll need to supply body accelerations to the constraint
    // acceleration error routines. Because the udots are zero, the body 
    // accelerations include velocity-dependent terms only, i.e. the coriolis 
    // accelerations. Those have already been calculated in the state, but they
    // are AC_GB, the coriolis accelerations in Ground. The constraint methods
    // want those relative to their Ancestor frames, which might not be Ground.
    const Array_<SpatialVec,MobilizedBodyIndex>&
        allAC_GB = getTreeVelocityCache(s).totalCoriolisAcceleration;

    // This array will be resized and filled with the Ancestor-relative
    // coriolis accelerations for the constrained bodies of each Constraint in 
    // turn; we're declaring it outside the 
    // loop to minimize heap allocation (resizing down doesn't normally free 
    // heap space).
    Array_<SpatialVec,ConstrainedBodyIndex> AC_AB;

    // Subarrays of these all-zero arrays will be used to supply zero qdotdots
    // (holonomic) or zero udots (nonholonomic and
    // acceleration-only) for each Constraint in turn; we're declaring 
    // them outside the loop to minimize heap allocation. They'll grow until 
    // they hit the maximum size needed by any Constraint.
    Array_<Real,      ConstrainedQIndex>    zeroQDotDot;
    Array_<Real,      ConstrainedUIndex>    zeroUDot;

    // Loop over all enabled constraints, ask them to generate constraint
    // errors, and collect those in the output bias vector.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        // Find this Constraint's err segments within the global array.
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp = includeP ? holoSeg.length    : 0;
        const int mv = includeV ? nonholoSeg.length : 0;
        const int ma = includeA ? accOnlySeg.length : 0;

        const ConstraintImpl& crep = constraints[cx]->getImpl();
        const int ncb = crep.getNumConstrainedBodies();

        // Now fill in coriolis accelerations. If the Ancestor is Ground
        // we're just reordering. If it isn't Ground we have to transform
        // the coriolis accelerations from Ground to Ancestor, at a cost
        // of 105 flops/constrained body (not just re-expressing).
        crep.convertBodyAccelToConstrainedBodyAccel(s, allAC_GB, AC_AB);

        // At this point AC_AB holds the coriolis accelerations of each
        // constrained body in A.

        if (mp) { // holonomic -- use velocity equations
            const int ncq = cInfo.getNumConstrainedQ();
            // Make sure we have enough zeroes.
            if (zeroQDotDot.size() < ncq) zeroQDotDot.resize(ncq, Real(0));

            // Make subarray; this does not require heap allocation.
            const ArrayViewConst_<Real,ConstrainedQIndex>    
                qdotdot0 = zeroQDotDot(ConstrainedQIndex(0), ncq);

            // The holonomic error slots start at beginning of bias array.
            ArrayView_<Real> paerr = biasArray(holoSeg.offset, mp);

            // Write errors into paerr, which is a segment of the bias argument.
            crep.calcPositionDotDotErrors(s, AC_AB, qdotdot0, paerr);
        }

        if (!(mv || ma))
            continue; // nothing else to do here

        const int ncu = cInfo.getNumConstrainedU();
        // Make sure we have enough zeroes for udots.
        if (zeroUDot.size() < ncu) zeroUDot.resize(ncu, Real(0));
        // Make a subarray of the right size.
        const ArrayViewConst_<Real,ConstrainedUIndex>    
            udot0 = zeroUDot(ConstrainedUIndex(0), ncu);

        if (mv) {   // non-holonomic constraints
            // The error slots begin after skipping the holonomic part of
            // the bias array. (That could be empty if P wasn't included.)
            const int start = mHolo+nonholoSeg.offset;
            ArrayView_<Real> vaerr = biasArray(start, mv);
            crep.calcVelocityDotErrors(s, AC_AB, udot0, vaerr);
        }
        if (ma) {   // acceleration-only constraints
            // The error slots begin after skipping the holonomic and 
            // non-holonomic parts of the bias array (those could be empty
            // if P or V weren't included).
            const int start = mHolo+mNonholo+accOnlySeg.offset;
            ArrayView_<Real> aerr = biasArray(start, ma);
            crep.calcAccelerationErrors(s, AC_AB, udot0, aerr);
        }
    }
}



//==============================================================================
//                              MULTIPLY BY Pq
//==============================================================================
// We have these mp constraint equations available:
// (1)  pverr(t,q;qdot) = Pq*qdot - Pt     (where Pq=P*N^+, Pt=c(t,q))      
//                      = P *u    - Pt
// with P=P(t,q), N=N(q), and Pt=Pt(t,q). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(nq+mp) time where mp is the total number of holonomic
// constraint equations. We expect to be given bias_p=-Pt as a 
// precalculated argument. (See calcBiasForMultiplyByPVA().)
//
// Given a q-like vector we can calculate
//      PqXqlike = Pq*qlike = P*N^+*qlike = pverr(t,q;qlike) - bias_p.
//
// The state must be realized to stage Position.
// All vectors must be using contiguous storage.
// Complexity is O(mp + n).
void SimbodyMatterSubsystemRep::
multiplyByPq(const State&   s,
             const Vector&  bias_p,
             const Vector&  qlike,
             Vector&        PqXqlike) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo = ic.totalNHolonomicConstraintEquationsInUse;
    const int m  = mHolo;
    const int nq = getNQ(s);
    const int nu = getNU(s);
    const int nb = getNumBodies();

    assert(bias_p.size() == m);
    assert(qlike.size() == nq);

    PqXqlike.resize(m);
    if (m == 0) return;

    if (nq == 0) { // not likely!
        PqXqlike.setToZero();
        return;
    }

    assert(bias_p.hasContiguousData() && qlike.hasContiguousData());
    assert(PqXqlike.hasContiguousData());

    // Generate a u-like Vector via ulike = N^-1 * qlike. Then use that to
    // calculate spatial velocities V_GB = J * ulike.
    Vector ulike(nu);
    Vector_<SpatialVec> V_GB(nb);
    multiplyByNInv(s, false, qlike, ulike);   // cheap
    multiplyBySystemJacobian(s, ulike, V_GB); // 12*(nu+nb) flops

    // Overlay Arrays on the Vectors' data so that we can manipulate small 
    // chunks of them repeatedly with no heap activity or virtual method calls.
    const ArrayViewConst_<SpatialVec,MobilizedBodyIndex> 
                                       V_GBArray  (&V_GB[0],    &V_GB[0]    +nb);
    const ArrayViewConst_<Real,QIndex> qArray     (&qlike[0],   &qlike[0]   +nq);
    const ArrayViewConst_<Real>        biasArray  (&bias_p[0],  &bias_p[0]  +m);
    ArrayView_<Real>                   PNInvqArray(&PqXqlike[0],&PqXqlike[0]+m);

    // This array will be resized and filled with the Ancestor-relative
    // velocities for the constrained bodies of each Constraint in turn; 
    // we're declaring it outside the loop to minimize heap allocation 
    // (resizing down doesn't normally free heap space).
    Array_<SpatialVec,ConstrainedBodyIndex> V_AB;
    // Same, but for each constraint's qdot subset.
    Array_<Real,ConstrainedQIndex> qdot;

    // Loop over all enabled constraints, ask them to generate constraint
    // errors, and collect those in the output vector, subtracting off the bias
    // as we go so we can go through the memory just once.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        // Find this Constraint's perr segment within the global array.
        const Segment& holoSeg = cInfo.holoErrSegment;
        const int mp = holoSeg.length;

        if (!mp)
            continue;

        const ConstraintImpl& crep  = constraints[cx]->getImpl();
        const int ncq = cInfo.getNumConstrainedQ();

        // Need body velocities relative to Ancestor body.
        crep.convertBodyVelocityToConstrainedBodyVelocity(s, V_GBArray, V_AB);

        qdot.resize(ncq);
        for (ConstrainedQIndex cqx(0); cqx < ncq; ++cqx)
            qdot[cqx] = qArray[cInfo.getQIndexFromConstrainedQ(cqx)];

        const int start = holoSeg.offset;
        const ArrayViewConst_<Real> bias  = biasArray(start, mp);
        ArrayView_<Real>            pverr = PNInvqArray(start, mp);
        crep.calcPositionDotErrors(s, V_AB, qdot, pverr);
        for (int i=0; i < mp; ++i)
            pverr[i] -= bias[i];
    }
}



//==============================================================================
//                                 CALC Pq
//==============================================================================
// Arranges for contiguous workspace if necessary, then makes repeated calls
// to multiplyByPq() to compute one column at a time of Pq (=P*N^-1).
// Complexity is O(n*mp + n*n) = O(n^2) <-- EXPENSIVE! Use transpose instead.
void SimbodyMatterSubsystemRep::
calcPq(const State& s, Matrix& Pq) const 
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = getNQ(s);

    Pq.resize(mp,nq);
    if (mp==0 || nq==0)
        return;

    Vector biasp(mp);
    calcBiasForMultiplyByPVA(s, true, false, false, biasp);
    Vector qlike(nq, Real(0));

    // If Pq's columns are contiguous we can avoid copying.
    const bool isContiguous = Pq(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : mp);
    for (int j=0; j < nq; ++j) {
        qlike[j] = 1; // column we're working on
        if (isContiguous)
            multiplyByPq(s, biasp, qlike, Pq(j));
        else {
            multiplyByPq(s, biasp, qlike, contig_col);
            Pq(j) = contig_col;
        }
        qlike[j] = 0;
    }
}



//==============================================================================
//                              MULTIPLY BY PVA
//==============================================================================
// We have these constraint equations available:
// (1)  pverr = Pq*qdot - Pt       (Pq==P*N^+, Pt==c(t,q))
// (2)  vaerr = V*udot  - b_v(t,q,u)
// (3)  aerr  = A*udot  - b_a(t,q,u)
// with P=P(t,q), N=N(q), V=V(t,q,u), A=A(t,q,u). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(m) time where m is the total number of constraint equations.
//
// Here we will use those equations to perform the multiplications by the
// matrices P,V, and/or A times a u-like vector: 
//             [P]           [ pverr(N*ulike) ]
// (4)  PVAu = [V] * ulike = [  vaerr(ulike)  ] - bias.
//             [A]           [   aerr(ulike)  ]
//
// We expect to be supplied as a precalculated argument "bias" the terms in 
// equations (1)-(3) that don't involve P,V, or A:
// (5)  bias=[  bias_p,   bias_v,      bias_a    ]
//          =[   -Pt,   -b_v(t,q,u), -b_a(t,q,u) ].
// See calcBiasForMultiplyByPVA() for how to get the bias terms.
//
// In general the state must be realized through Velocity stage, but if the 
// system contains only holonomic constraints, or if only P is included, then 
// bias is just bias_p and the result is only time- and position-dependent 
// since we just need to use eq. (1). In that case we require only that the
// state be realized to stage Position.
//
// All of the Vector arguments must use contiguous storage.
// Complexity is O(m+n).
void SimbodyMatterSubsystemRep::
multiplyByPVA(  const State&     s,
                bool             includeP,
                bool             includeV,
                bool             includeA,
                const Vector&    bias,
                const Vector&    ulike,
                Vector&          PVAu) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;

    const int m  = mHolo+mNonholo+mAccOnly;
    const int nq = getNQ(s);
    const int nu = getNU(s);
    const int nb = getNumBodies();

    assert(bias.size() == m);
    assert(ulike.size() == nu);

    PVAu.resize(m);
    if (m == 0) return;

    if (nu == 0) { // not likely!
        PVAu.setToZero();
        return;
    }

    assert(bias.hasContiguousData() && ulike.hasContiguousData());
    assert(PVAu.hasContiguousData());

    // Generate body spatial velocities V=J*u or first term of the
    // body spatial accelerations A=J*udot + Jdot*u, depending on how we're
    // interpreting the ulike argument (as a u for holonomic constraints,
    // and as udot for everything else).
    Vector_<SpatialVec> Julike(nb);
    multiplyBySystemJacobian(s, ulike, Julike); // 12*(nu+nb) flops

    // Julike serves as V_GB when we're interpreting ulike as u.
    const ArrayViewConst_<SpatialVec,MobilizedBodyIndex> 
        allV_GB(&Julike[0], &Julike[0] + nb);

    // If we're doing any nonholonomic or acceleration-only constraints, we'll 
    // finish calculating body spatial accelerations and put them here.
    Array_<SpatialVec,MobilizedBodyIndex> allA_GB;
    if (mNonholo || mAccOnly) {
        allA_GB.resize(nb);
        const Array_<SpatialVec>& 
            allAC_GB = getTreeVelocityCache(s).totalCoriolisAcceleration;
        for (MobilizedBodyIndex b(0); b < nb; ++b)
            allA_GB[b] = allV_GB[b] + allAC_GB[b]; // i.e., J*udot + Jdot*u
    }

    // If we're going to be dealing with holonomic (position) constraints,
    // generate a q-like Vector via qlike = N * ulike since the position
    // error derivative routine wants qdots.
    Vector qlike(nq);
    if (mHolo)
        multiplyByN(s, false, ulike, qlike);   // cheap

    // Overlay Arrays on the Vectors' data so that we can manipulate small 
    // chunks of them repeatedly with no heap activity or virtual method calls.
    const ArrayViewConst_<Real,UIndex>  uArray   (&ulike[0], &ulike[0] + nu);
    const ArrayViewConst_<Real,QIndex>  qArray   (&qlike[0], &qlike[0] + nq);
    const ArrayViewConst_<Real>         biasArray(&bias[0],  &bias[0]  + m );
    ArrayView_<Real>                    PVAuArray(&PVAu[0],  &PVAu[0]  + m );

    // This array will be resized and filled with the Ancestor-relative
    // velocities for the constrained bodies of each holonomic Constraint in 
    // turn; we're declaring it outside the loop to minimize heap allocation 
    // (resizing down doesn't normally free heap space). This won't be used 
    // if we aren't processing holonomic constraints.
    Array_<SpatialVec,ConstrainedBodyIndex> V_AB;
    // Same, but for each holonomic constraint's qdot subset.
    Array_<Real,ConstrainedQIndex> qdot;

    // This array will be resized and filled with the Ancestor-relative
    // accelerations for the constrained bodies of each velocity
    // or acceleration-only Constraint in turn. This won't be used if we have 
    // only holonomic constraints.
    Array_<SpatialVec,ConstrainedBodyIndex> A_AB;
    // Same, but for each nonholonomic/acconly constraint's udot subset.
    Array_<Real,ConstrainedUIndex> udot;

    // Loop over all enabled constraints, ask them to generate constraint
    // errors, and collect those in the output argument PVAu. Remove bias
    // as we go so we only have to touch the memory once.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        // Find this Constraint's err segments within the global array.
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp = includeP ? holoSeg.length    : 0;
        const int mv = includeV ? nonholoSeg.length : 0;
        const int ma = includeA ? accOnlySeg.length : 0;

        const ConstraintImpl& crep = constraints[cx]->getImpl();

        if (mp) { // holonomic -- use velocity equations
            const int ncq = cInfo.getNumConstrainedQ();
            // Need body velocities in ancestor frame.
            crep.convertBodyVelocityToConstrainedBodyVelocity(s, allV_GB, V_AB);
            qdot.resize(ncq);
            for (ConstrainedQIndex cqx(0); cqx < ncq; ++cqx)
                qdot[cqx] = qArray[cInfo.getQIndexFromConstrainedQ(cqx)];

            // The error slots start at the beginning of the bias array.
            const int start = holoSeg.offset;
            const ArrayViewConst_<Real> bias  = biasArray(start, mp);
            ArrayView_<Real>            pverr = PVAuArray(start, mp);
            crep.calcPositionDotErrors(s, V_AB, qdot, pverr);
            for (int i=0; i < mp; ++i)
                pverr[i] -= bias[i];
        }

        if (!(mv || ma))
            continue; // nothing else to do here

        const int ncu = cInfo.getNumConstrainedU();

        // Now fill in accelerations. If the Ancestor is Ground
        // we're just reordering. If it isn't Ground we have to transform
        // the accelerations from Ground to Ancestor, at a cost
        // of 105 flops/constrained body (not just re-expressing).
        crep.convertBodyAccelToConstrainedBodyAccel(s, allA_GB, A_AB);
        // At this point A_AB holds the accelerations of each
        // constrained body in A.

        // Now pack together the appropriate udots.
        udot.resize(ncu);
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
            udot[cux] = uArray[cInfo.getUIndexFromConstrainedU(cux)];

        if (mv) {   // non-holonomic constraints
            // The error slots begin after skipping the holonomic part of
            // the arrays. (Those could be empty if P wasn't included.)
            const int start = mHolo + nonholoSeg.offset;
            const ArrayViewConst_<Real> bias  = biasArray(start, mv);
            ArrayView_<Real>            vaerr = PVAuArray(start, mv);
            crep.calcVelocityDotErrors(s, A_AB, udot, vaerr);
            for (int i=0; i < mv; ++i)
                vaerr[i] -= bias[i];
        }

        if (ma) {   // acceleration-only constraints
            // The error slots begin after skipping the holonomic and 
            // non-holonomic parts of the arrays (those could be empty
            // if P or V weren't included).
            const int start = mHolo+mNonholo+accOnlySeg.offset;
            const ArrayViewConst_<Real> bias = biasArray(start, ma);
            ArrayView_<Real>            aerr = PVAuArray(start, ma);
            crep.calcAccelerationErrors(s, A_AB, udot, aerr);
            for (int i=0; i < ma; ++i)
                aerr[i] -= bias[i];
        }
    }
}



//==============================================================================
//                                CALC PVA
//==============================================================================
// Arranges for contiguous workspace if necessary, then makes repeated calls
// to multiplyByPVA() to compute one column at a time of G=PVA or a submatrix.
// This is particularly useful for computing [P;V] which is the velocity-level
// constraint projection matrix.
// Complexity is O(n*m + n*n) = O(n^2) <-- EXPENSIVE! Use transpose instead.
void SimbodyMatterSubsystemRep::
calcPVA(const State&     s,
        bool             includeP,
        bool             includeV,
        bool             includeA,
        Matrix&          PVA) const
{
    const SBInstanceCache& ic  = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = includeP ? 
        ic.totalNHolonomicConstraintEquationsInUse : 0;
    const int mNonholo = includeV ? 
        ic.totalNNonholonomicConstraintEquationsInUse : 0;
    const int mAccOnly = includeA ? 
        ic.totalNAccelerationOnlyConstraintEquationsInUse : 0;

    const int m  = mHolo+mNonholo+mAccOnly;
    const int nu = getNU(s);

    PVA.resize(m,nu);
    if (m==0 || nu==0)
        return;

    Vector bias(m);
    calcBiasForMultiplyByPVA(s, includeP, includeV, includeA, bias);
    Vector ulike(nu, Real(0));

    // If PVA's columns are contiguous we can avoid copying.
    const bool isContiguous = PVA(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : m); // temp space if needed
    for (int j=0; j < nu; ++j) {
        ulike[j] = 1; // column we're working on
        if (isContiguous) {
            multiplyByPVA(s, includeP, includeV, includeA, bias, ulike, 
                          PVA(j));
        } else {
            multiplyByPVA(s, includeP, includeV, includeA, bias, ulike, 
                          contig_col);
            PVA(j) = contig_col;
        }
        ulike[j] = 0;
    }
}



// =============================================================================
//                            CALC G MInv G^T
// =============================================================================
/* We will later calculate multipliers lambda as
    C lambda = aerr
where C(t,q,u) = G M^-1 ~G is the constraint compliance matrix. We need a fast 
way to explicitly calculate the symmetric mXm matrix C, provided here.

Optimally, we would calculate it in O(m^2) time. I don't know how to calculate 
it that fast, but using m calls to operator sequence:
           Gt_j      = Gt lambda_j        O(m+n)
           MInvGt_j  = M^-1 Gt_j          O(n)
    C(j) = GMInvGt_j = G MInvGt_j         O(m+n)
we can calculate it in O(m^2 + mn) time (lambda_j is a unit vector with element
j=1). As long as m << n, and especially if m is a small constant independent 
of n, and even better if we've partitioned it into little subblocks, this is 
all very reasonable. One slip up and you'll toss in a factor of m*n^2 or m^2*n 
and screw this up -- be careful!

When there is prescribed motion in the system the matrix we want is
Gr Mrr^-1 ~Gr. That is still an mXm matrix and we are able to produce it with no
visible effort due to the definition of our a=M^-1*f operator. It ignores the 
prescribed part fp of f (here a column of ~Gp), and returns zeroes in the 
prescribed part ap of a. Those zeroes have the effect of removing the Gp columns
of G in the final operation. The resulting matrix is *not* a submatrix of C!

Note that we do not require contiguous storage for GMInvGt's columns, although 
we'll take advantage of it if they are. If not, we'll work in a contiguous temp 
and then copy back. This is because we want to allow any matrix at the Simbody 
API level and we don't want to force the API method to have to allocate a whole 
new mXm matrix when all we need for a temporary here is an m-length temporary.

Complexity is O(m^2 + m*n) = O(m*n) when m < n which is almost always.

TODO: the resulting matrix C is symmetric but we are not taking advantage of 
that here. */
void SimbodyMatterSubsystemRep::
calcGMInvGt(const State&   s,
            Matrix&        GMInvGt) const
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int m        = mHolo+mNonholo+mAccOnly;  
    const int nu       = getNU(s);

    GMInvGt.resize(m,m);
    if (m==0) return;

    // If the output matrix doesn't have columns in contiguous memory, we'll
    // allocate a contiguous-memory temp that we can use to hold one column
    // at a time as we compute them.
    const bool columnsAreContiguous = GMInvGt(0).hasContiguousData();
    Vector GMInvGt_j(columnsAreContiguous ? 0 : m);

    // These two temporaries are always needed to hold one column of Gt,
    // then one column of M^-1 * Gt.
    Vector Gtcol(nu), MInvGtcol(nu);

    // Precalculate bias so we can perform multiplication by G efficiently.
    Vector bias(m);
    calcBiasForMultiplyByPVA(s,true,true,true,bias);
   
    // Lambda is used to pluck out one column at a time of Gt. Exactly one
    // element at a time of lambda will be 1, the rest are 0.
    Vector lambda(m, Real(0));

    for (int j=0; j < m; ++j) {
        lambda[j] = 1;
        multiplyByPVATranspose(s, true, true, true, lambda, Gtcol);
        lambda[j] = 0;
        multiplyByMInv(s, Gtcol, MInvGtcol);
        if (columnsAreContiguous)
            multiplyByPVA(s, true, true, true, bias, MInvGtcol, GMInvGt(j));
        else {
            multiplyByPVA(s, true, true, true, bias, MInvGtcol, GMInvGt_j);
            GMInvGt(j) = GMInvGt_j;
        }
    }
} 



// =============================================================================
//                     SOLVE FOR CONSTRAINT IMPULSES
// =============================================================================
// Current implementation computes G*M^-1*~G, factors it, and does a single
// solve all at great expense.
// TODO: should realize factored matrix if needed and reuse if possible.
void SimbodyMatterSubsystemRep::
solveForConstraintImpulses(const State&     state,
                           const Vector&    deltaV,
                           Vector&          impulse) const
{
    Matrix GMInvGt;
    calcGMInvGt(state, GMInvGt);
    // MUST DUPLICATE SIMBODY'S METHOD HERE:
    const Real conditioningTol = GMInvGt.nrow() * Eps34;
    FactorQTZ qtz(GMInvGt, conditioningTol); 
    qtz.solve(deltaV, impulse);
}



// =============================================================================
//                     CALC BODY ACCELERATION FROM UDOT
// =============================================================================
// Input and output vectors must use contiguous storage.
// The knownUDot argument must be of length nu; the output argument will be
// resized if necessary to length nb.
// Cost is 12*nu + 18*nb flops.
void SimbodyMatterSubsystemRep::
calcBodyAccelerationFromUDot(const State&           s,
                             const Vector&          knownUDot,
                             Vector_<SpatialVec>&   A_GB) const {
    const int nu = getNU(s);
    const int nb = getNumBodies();
    A_GB.resize(nb); // always at least Ground

    // Should have been checked before we got here.
    assert(knownUDot.size() == nu);
    if (nu == 0) {
        A_GB.setToZero();
        return;
    }

    assert(knownUDot.hasContiguousData() && A_GB.hasContiguousData());
    const Real* knownUdotPtr = &knownUDot[0];
    SpatialVec* aPtr         = &A_GB[0];

    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const SBTreeVelocityCache& tvc = getTreeVelocityCache(s);

    // Note: this is equivalent to calculating J*u with 
    // multiplyBySystemJacobian() and adding in the totalCoriolisAcceleration
    // for each body (which is Jdot*u). (sherm 110829: I tried it both ways)

    // Sweep outward and delegate to RB nodes.
    for (int i=0; i<(int)rbNodeLevels.size(); i++)
        for (int j=0; j<(int)rbNodeLevels[i].size(); j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcBodyAccelerationsFromUdotOutward
               (tpc,tvc,knownUdotPtr,aPtr);
        }
}



//==============================================================================
//                   CALC CONSTRAINT ACCELERATION ERRORS
//==============================================================================
// Note: similar to multiplyByPVA() except:
//  - we're using the holonomic *second* derivative here
//  - we don't want to get rid of the bias terms
//  - state must be realized through Velocity stage even if there are only
//    holonomic constraints
//  - no option to work with only a subset of the constraints
//
// We have these constraint equations available:
// (1)  paerr = Pq*qdotdot - b_p(t,q,u)     (Pq = P*N^-1)
// (2)  vaerr = V*udot     - b_v(t,q,u)
// (3)  aerr  = A*udot     - b_a(t,q,u)
// with P=P(t,q), N=N(q), V=V(t,q,u), A=A(t,q,u). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(m) time where m is the total number of constraint equations.
//
// All of the Vector arguments must use contiguous storage. Output is resized
// to m=mp+mv+ma.  Complexity is O(m+n).
void SimbodyMatterSubsystemRep::
calcConstraintAccelerationErrors
       (const State&                s,
        const Vector_<SpatialVec>&  A_GB,
        const Vector_<Real>&        udot,
        const Vector_<Real>&        qdotdot,
        Vector&                     pvaerr) const
{
    const SBInstanceCache&     ic  = getInstanceCache(s);

    // Global problem dimensions.
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;

    const int m  = mHolo+mNonholo+mAccOnly;
    const int nq = getNQ(s);
    const int nu = getNU(s);
    const int nb = getNumBodies();

    assert(A_GB.size()    == nb);
    assert(udot.size()    == nu);
    assert(qdotdot.size() == nq);

    pvaerr.resize(m);
    if (m == 0) return;

    if (nu == 0) { // not likely!
        pvaerr.setToZero();
        return;
    }

    assert(udot.hasContiguousData() && qdotdot.hasContiguousData());
    assert(A_GB.hasContiguousData() && pvaerr.hasContiguousData());

    // Overlay Arrays on the Vectors' data so that we can manipulate small 
    // chunks of them repeatedly with no heap activity or virtual method calls.
    const ArrayViewConst_<SpatialVec, MobilizedBodyIndex> 
                                        allA_GB (&A_GB[0],    &A_GB[0]    + nb);
    const ArrayViewConst_<Real,UIndex>  udArray (&udot[0],    &udot[0]    + nu);
    const ArrayViewConst_<Real,QIndex>  qddArray(&qdotdot[0], &qdotdot[0] + nq);
    ArrayView_<Real>                    allAerr (&pvaerr[0],  &pvaerr[0]  + m );

    // These arrays will be resized and filled with the input needs of each 
    // Constraint in turn. We're declaring them outside the loop to minimize 
    // heap allocation (resizing down doesn't normally free heap space). 
    Array_<SpatialVec,ConstrainedBodyIndex> A_AB;
    Array_<Real,ConstrainedQIndex> qdd; // holonomic only
    Array_<Real,ConstrainedUIndex> ud;  // nonholonomic or acc-only

    // Loop over all enabled constraints, ask them to generate constraint
    // errors, and collect those in the output argument pvaerr.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstancePerConstraintInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        // Find this Constraint's err segments within the global array.
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mp = holoSeg.length;
        const int mv = nonholoSeg.length;
        const int ma = accOnlySeg.length;

        const ConstraintImpl& crep = constraints[cx]->getImpl();

        // Now fill in accelerations. If the Ancestor is Ground
        // we're just reordering. If it isn't Ground we have to transform
        // the accelerations from Ground to Ancestor, at a cost
        // of 105 flops/constrained body (not just re-expressing).
        crep.convertBodyAccelToConstrainedBodyAccel(s, allA_GB, A_AB);

        // At this point A_AB holds the accelerations of each
        // constrained body in A.


        if (mp) { // holonomic
            // Now pack together the appropriate qdotdots.
            const int ncq = cInfo.getNumConstrainedQ();
            qdd.resize(ncq);
            for (ConstrainedQIndex cqx(0); cqx < ncq; ++cqx)
                qdd[cqx] = qddArray[cInfo.getQIndexFromConstrainedQ(cqx)];

            // The error slots start at the beginning of the pvaerr array.
            const int start = holoSeg.offset;
            ArrayView_<Real>  paerr = allAerr(start, mp);
            crep.calcPositionDotDotErrors(s, A_AB, qdd, paerr);
        }

        if (!(mv || ma))
            continue; // nothing else to do here

        // Now pack together the appropriate udots.
        const int ncu = cInfo.getNumConstrainedU();
        ud.resize(ncu);
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux)
            ud[cux] = udArray[cInfo.getUIndexFromConstrainedU(cux)];

        if (mv) {   // non-holonomic constraints
            // The error slots begin after skipping the holonomic part of
            // the arrays.
            const int start = mHolo + nonholoSeg.offset;
            ArrayView_<Real> vaerr = allAerr(start, mv);
            crep.calcVelocityDotErrors(s, A_AB, ud, vaerr);
        }

        if (ma) {   // acceleration-only constraints
            // The error slots begin after skipping the holonomic and 
            // non-holonomic parts of the arrays.
            const int start = mHolo+mNonholo+accOnlySeg.offset;
            ArrayView_<Real> aerr = allAerr(start, ma);
            crep.calcAccelerationErrors(s, A_AB, ud, aerr);
        }
    }
}




// =============================================================================
//                          PRESCRIBE Q, PRESCRIBE U
// =============================================================================
// These are solvers that set continuous state variables q or u to their 
// prescribed values q(t) or u(t,q) that will already have been computed. 
// Note that prescribed udot=udot(t,q,u) is not dealt with here because it does 
// not involve a state change.
bool SimbodyMatterSubsystemRep::prescribeQ(State& s) const {
    const SBModelCache&    mc = getModelCache(s);
    const SBInstanceCache& ic = getInstanceCache(s);

    const int npq = ic.getTotalNumPresQ();
    const int nzq = ic.getTotalNumZeroQ();
    if (npq==0 && nzq==0) return false; // don't invalidate positions

    // copy prescribed q's from cache to state
    // set known-zero q's to zero (or reference configuration)
    const SBTimeCache& tc = getTimeCache(s);
    Vector& q = updQ(s); // this Subsystem's q's, now invalidated

    for (int i=0; i < npq; ++i)
        q[ic.presQ[i]] = tc.presQPool[i];

    //TODO: this isn't right -- need to use reference config
    //q's which will be 1000 for quaternion.
    for (int i=0; i < nzq; ++i)
        q[ic.zeroQ[i]] = 0;

    return true;
}

bool SimbodyMatterSubsystemRep::prescribeU(State& s) const {
    const SBModelCache&    mc = getModelCache(s);
    const SBInstanceCache& ic = getInstanceCache(s);

    const int npu = ic.getTotalNumPresU();
    const int nzu = ic.getTotalNumZeroU();
    if (npu==0 && nzu==0) return false; // don't invalidate velocities

    // copy prescribed u's from cache to state
    // set known-zero u's to zero
    const SBConstrainedPositionCache& cpc = getConstrainedPositionCache(s);
    Vector& u = updU(s); // this Subsystem's u's, now invalidated

    for (int i=0; i < npu; ++i)
        u[ic.presU[i]] = cpc.presUPool[i];

    for (int i=0; i < nzu; ++i)
        u[ic.zeroU[i]] = 0;

    return true;
}
//............................. PRESCRIBE Q, U .................................



//==============================================================================
//                                  PROJECT Q
//==============================================================================
/* A note on state variable weights Wq and Wu:
  - q and u weights are not independent.
  - We have    qdot=N  u 
                  u=N\ qdot
    where N\ is the left pseudoinverse of N, such that N\ N = I (N N\ != I).
  - We consider u weights Wu primary and want the weighted variables to be 
    related like the unweighted ones: qdot=N*u so qdotw=N*uw, where
    uw = Wu*u and we should have qdotw=Wq*qdot for some matrix Wq.
        qdotw = Wq qdot
              = N Wu u
              = N Wu N\ qdot
        =>     Wq  = N Wu N\
  - Wu is diagonal with full rank nu, but Wq is block diagonal and rank
    deficient since it is nqXnq but has rank nu.

The pseudoinverse of Wq, Wq\ = N Wu\ N\ (verified in Matlab). But since Wq has
neither full row nor full column rank we have Wq Wq\ = Wq\ Wq = N N\ != I.

We have fast operators for multiplying Wu, N, and N\ by columns, but not for 
producing Wq so we just create it operationally as we go. */
int SimbodyMatterSubsystemRep::projectQ
   (State&                  s, 
    Vector&                 qErrest, // q error estimate or empty 
    const ProjectOptions&   opts, 
    ProjectResults&         results) const
{
    SimTK_STAGECHECK_GE(getStage(s), Stage::Position,
        "SimbodyMatterSubsystemRep::projectQ()");

    results.clear();
    const Real consAccuracy = opts.getRequiredAccuracy();
    // Normally we'll use an RMS norm for the perrs.
    const bool useNormInf = opts.isOptionSet(ProjectOptions::UseInfinityNorm);
    // Force projection even if accuracy ok on entry.
    const bool forceOneIter = opts.isOptionSet(ProjectOptions::ForceProjection);
    // Normally we'll throw an exception with a helpful message. If this is
    // set we'll quietly return status instead.
    const bool dontThrow = opts.isOptionSet(ProjectOptions::DontThrow);

    const int mHolo  = getNumHolonomicConstraintEquationsInUse(s);
    const int mQuats = getNumQuaternionsInUse(s);

    // Thise array contains a unique MultiplierIndex for each perr 
    // that is to be enforced by projection. Baumgarte-stabilized constraint
    // equations don't count since then position constraints are 
    // softly enforced through the acceleration error equation instead.
    const Array_<MultiplierIndex> enforcedPosConsEqns =
                                            findEnforcedPosConstraintEqns(s);
    const int mfp = (int)enforcedPosConsEqns.size();

    // This is a const view into the State; the contents it refers to will 
    // change though.
    const VectorView pErrs = getQErr(s)(0,mHolo); // just leave off quaternions
    const VectorView quatErrs = getQErr(s)(mHolo,mQuats); // quaternions

    // We don't weight the quaternion errors.
    const VectorView perrWeights = getQErrWeights(s)(0,mHolo); // 1/unit error (Tp)

    // Determine norms on entry.
    int worstPerr, worstQuatErr;
    Vector scaledPerrs(mfp); // scaled and packed
    packActiveMultipliers(enforcedPosConsEqns,
                          pErrs.rowScale(perrWeights),
                          scaledPerrs);
    const Real perrNormOnEntry = useNormInf ? scaledPerrs.normInf(&worstPerr)
                                            : scaledPerrs.normRMS(&worstPerr);
    // Map worstPerr back to real index.
    if (worstPerr >= 0) // will be -1 if mfp==0
        worstPerr = enforcedPosConsEqns[worstPerr];
    const Real quatNormOnEntry = useNormInf ? quatErrs.normInf(&worstQuatErr)
                                            : quatErrs.normRMS(&worstQuatErr);
    
    Real normOnEntry;
    if (perrNormOnEntry >= quatNormOnEntry) {
        results.setNormOnEntrance(perrNormOnEntry, worstPerr);
        normOnEntry = perrNormOnEntry;
    }
    else {
        results.setNormOnEntrance(quatNormOnEntry, mHolo + worstQuatErr);
        normOnEntry = quatNormOnEntry;
    }

    if (normOnEntry > opts.getProjectionLimit()) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        return 1;
    }


    // Return quickly if (a) constraint norm is zero (probably because there
    // aren't any), or (b) constraints are already satisfied and we're not
    // being forced to go ahead anyway.
    if (    perrNormOnEntry == 0 
        || (perrNormOnEntry <= consAccuracy && !forceOneIter)) {
        // Perrs are good enough already. Might still need to project
        // quaternions, but that doesn't take long. The only way this can
        // fail is if some of the quaternions are prescribed but prescribedQ()
        // wasn't called earlier as it should have been..
        if (quatNormOnEntry > consAccuracy || forceOneIter) {
            const bool anyQuatChange = normalizeQuaternions(s,qErrest);
            results.setAnyChangeMade(anyQuatChange);
            const Real quatNorm = useNormInf ? quatErrs.normInf(&worstQuatErr)
                                             : quatErrs.normRMS(&worstQuatErr);
            results.setNormOnExit(quatNorm);
            if (quatNorm > consAccuracy) {
                results.setExitStatus(ProjectResults::FailedToAchieveAccuracy);
                if (!dontThrow) {
                    SimTK_ERRCHK2_ALWAYS(quatNorm <= consAccuracy, 
                         "SimbodyMatterSubsystem::projectQ()",
                         "Failed to normalize quaternions. Norm achieved=%g"
                         " but required norm=%g. Did you forget to call"
                         " prescribeQ()?", quatNorm, consAccuracy);
                }
                return 1;
            }
        } else {    // both perrs and quatErrs were good on entry
            // numIterations==0.
            results.setAnyChangeMade(false);
            results.setNormOnExit(normOnEntry);
        }
        results.setExitStatus(ProjectResults::Succeeded);
        return 0;
    }


    // We're going to have to project constraints. Get the remaining options.


    // This is the factor by which we try to achieve a tighter accuracy
    // than requested. E.g. if overshootFactor=0.1 then we attempt 10X 
    // tighter accuracy if we can get it. But we won't fail as long as
    // we manage to reach consAccuracy.
    const Real overshootFactor = opts.getOvershootFactor();
    const Real consAccuracyToTryFor = 
        std::max(overshootFactor*consAccuracy, SignificantReal);

    // Check whether we should stop if we see the solution diverging
    // which should not happen when we're in the neighborhood of a solution
    // on entry. This is always set while integrating, except during
    // initialization.
    const bool localOnly = opts.isOptionSet(ProjectOptions::LocalOnly);
    // We are permitted to use an out-of-date Jacobian for projection unless
    // this is set. TODO: always using full Newton at the moment.
    const bool forceFullNewton =
        opts.isOptionSet(ProjectOptions::ForceFullNewton);

    // Get problem dimensions.
    const SBInstanceCache& ic = getInstanceCache(s);
    const int nq     = getNQ(s);
    const int nfq    = ic.getTotalNumFreeQ();
    const int nu     = getNU(s);
    const bool hasPrescribedMotion = (nfq != nq);

    /* Using only enforce position constraint equations and free generalized
    coordinates q, solve:
           (Tp Pq Wq\) dq_WLS  = Tp perr
                           dq  = Wq\ dq_WLS
                            q -= dq
    until RMS(Tp perr) <= 0.1*accuracy.
    
    But Pq=P N\, Wq\ = N Wu\ N\ so Pq Wq\ = P Wu\ N\, since N\ N=I. So we can 
    rewrite the above:
        
         (Tp P Wu\ N\) dq_WLS  = Tp perr
                           dq  = N Wu\ N\ dq_WLS
                            q -= dq
    
    We define Pqwt = ~Pqw = ~(Tp P Wu\ N\) = ~N\ Wu\ ~P Tp (diagonal weights are
    symmetric). We only retain the nfq rows that correspond to free (non 
    prescribed) q's, and retain only the mfp columns corresponding to enforced 
    position constraint equations.
    
    This is a nonlinear least squares problem. Below is a full Newton iteration 
    since we recalculate the iteration matrix each time around the loop. TODO: 
    a solution could be found using the same iteration matrix, since we are 
    projecting from (presumably) not too far away. Q1: will it be the same 
    solution? Q2: if not, does it matter? */

    // These will be updated as we go.
    Real perrNormAchieved = perrNormOnEntry;
    Real quatNormAchieved = quatNormOnEntry;

    // We always use absolute scaling for q's, derived from the absolute
    // scaling of u's.
    const Vector& Wu    = getUWeights(s);         // 1/unit change (Wu)
    const Vector  WuInv = Wu.elementwiseInvert(); // Wu\

    Real lastChangeMadeWRMS = 0; // size of last change in weighted dq
    int nItsUsed = 0;


    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is sloppy; should depend on constraint tolerance
    // and rank should be saved and reused in velocity and acceleration
    // constraint methods (or should be calculated elsewhere and passed in).
    const Real conditioningTol = mfp         
      //* SignificantReal; -- too tight
        * SqrtEps;

    // Keep the starting q in case we have to restore it, which we'll do
    // if the attempts here make the constraint norm worse.
    const Vector saveQ = getQ(s);

    Matrix Pqwrt(nfq,mfp);
    Vector dfq_WLS(nfq), du(nu), dq(nq); // = Wq\ dq_WLS
    Vector udfq_WLS(hasPrescribedMotion ? nq : 0); // unpacked if needed
    udfq_WLS.setToZero(); // must initialize unwritten elements
    FactorQTZ Pqwr_qtz;
    Real prevPerrNormAchieved = perrNormAchieved; // watch for divergence
    bool diverged = false;
    const int MaxIterations  = 20;
    do {
        calcWeightedPqrTranspose(s, enforcedPosConsEqns,
                                 perrWeights, WuInv, Pqwrt);//nfq X mfp

        // This factorization acts like a pseudoinverse.
        Pqwr_qtz.factor<Real>(~Pqwrt, conditioningTol); 

        //printf("projectQ %d: m=%d condTol=%g rank=%d rcond=%g\n",
        //    nItsUsed, Pqwrt.ncol(), conditioningTol, Pqwr_qtz.getRank(),
        //    Pqwr_qtz.getRCondEstimate());

        Pqwr_qtz.solve(scaledPerrs, dfq_WLS); // this is weighted dq_WLS=Wq*dq
        lastChangeMadeWRMS = dfq_WLS.normRMS(); // change in weighted norm

        // Switch back to unweighted dq = Wq\ dq_WLS = N Wu\ N\ dq_WLS.
        if (hasPrescribedMotion) {
            unpackFreeQ(s, dfq_WLS, udfq_WLS); // zeroes in q_p slots
            multiplyByNInv(s,false,udfq_WLS,du);
        } else {
            multiplyByNInv(s,false,dfq_WLS,du);
        }
        // Here du = du_WLS = N\ dq_WLS.
        du.rowScaleInPlace(WuInv); // Now du = Wu\ du_WLS.
        multiplyByN(s,false,du,dq);    // dq = N du

        // This causes quaternions to become unnormalized, but it doesn't
        // matter because N is calculated from the unnormalized q's so
        // scales dq to match.
        updQ(s) -= dq; // this is unweighted dq
        results.setAnyChangeMade(true);

        // Now recalculate the position constraint errors at the new q.
        realizeSubsystemPosition(s); // pErrs changes here

        packActiveMultipliers(enforcedPosConsEqns,
                              pErrs.rowScale(perrWeights), // Tp * pErrs
                              scaledPerrs);
        perrNormAchieved = useNormInf ? scaledPerrs.normInf()
                                      : scaledPerrs.normRMS();
        ++nItsUsed;

        if (localOnly && nItsUsed >= 2 
            && perrNormAchieved > prevPerrNormAchieved) {
            // perr norm got worse; restore to end of previous iteration
            updQ(s) += dq;
            realizeSubsystemPosition(s); // pErrs changes here
            packActiveMultipliers(enforcedPosConsEqns,
                                  pErrs.rowScale(perrWeights),
                                  scaledPerrs);
            perrNormAchieved = useNormInf ? scaledPerrs.normInf()
                                          : scaledPerrs.normRMS();
            diverged = true;
            break; // diverging -- quit now to prevent a bad solution
        }

        prevPerrNormAchieved = perrNormAchieved;

    } while (   perrNormAchieved > consAccuracyToTryFor
             && nItsUsed < MaxIterations);

    results.setNumIterations(nItsUsed);

    //printf("        perrNormAchieved=%g in %d its\n",perrNormAchieved, nItsUsed);

    // Make sure we achieved at least the required constraint accuracy. If not 
    // we'll return with an error. If we see that the norm has been made worse
    // than it was on entry, we'll restore the state to what it was on entry. 
    // Otherwise we'll return with the improved-but-not-good-enough result.
    if (perrNormAchieved > consAccuracy) {
        if (perrNormAchieved >= perrNormOnEntry) { // made it worse
            updQ(s) = saveQ; // revert
            realizeSubsystemPosition(s);
            perrNormAchieved = perrNormOnEntry;
        }
     
        results.setNormOnExit(perrNormAchieved);

        if (diverged) {
            results.setExitStatus(ProjectResults::FailedToConverge);
            if (!dontThrow) {
                SimTK_ERRCHK_ALWAYS(!diverged,
                    "SimbodyMatterSubsystem::projectQ()",
                    "Attempt to project constraints locally diverged.");
            }
        } else {
            results.setExitStatus(ProjectResults::FailedToAchieveAccuracy);
            if (!dontThrow) {
                SimTK_ERRCHK3_ALWAYS(perrNormAchieved <= consAccuracy, 
                    "SimbodyMatterSubsystem::projectQ()",
                    "Failed to achieve required accuracy %g. Norm on entry "
                    " was %g; norm on exit %g. You might need a better"
                    " starting configuration, or if there are prescribed or "
                    " locked q's you might have to free some of them.", 
                    consAccuracy, perrNormOnEntry, perrNormAchieved);
            }
        }

        return 1;
    }

    // Position constraint errors were successfully driven to consAccuracy.

    /* Next, remove the corresponding error from the integrator's error 
    estimate.   
        (Tp Pq Wq\)_r dqr_WLS = (Tp Pq Wq\)_r (Wq qErrest)_r
                       dq_WLS = unpack(dqr_WLS) (with zero fill)
                          dq  = Wq\ dq_WLS
                              = N Wu\ N\ dq_WLS
                     qErrest -= dq
    No iteration is required. Note that we only include the enforced constraints
    that we used above, not all enabled constraint equations.
    
    We can simplify the RHS of the first equation above:
           (Tp Pq Wq\)_r (Wq qErrest)_r = Tp Pq unpack(qErrest_r)
    for which we have an O(n) operator to use for the matrix-vector product. 
    Proof: expand Pq, Wq\, and Wq and cancel Wu\ Wu and N\ N. */
    if (qErrest.size()) {
        // Work in Wq-norm
        Vector Tp_Pq_qErrest0, packed_Tp_Pq_qErrest0, bias_p;
        calcBiasForMultiplyByPq(s, bias_p);

        // Switch back to unweighted dq = Wq\ dq_WLS = N Wu\ N\ dq_WLS.
        if (hasPrescribedMotion) {
            Vector qErrest_0(qErrest);
            zeroKnownQ(s, qErrest_0); // zero out prescribed entries
            multiplyByPq(s, bias_p, qErrest_0, Tp_Pq_qErrest0); // (Pq*qErrest)_r
            Tp_Pq_qErrest0.rowScaleInPlace(perrWeights); // now Tp*(Pq*qErrest)_r
            packActiveMultipliers(enforcedPosConsEqns,
                                  Tp_Pq_qErrest0,
                                  packed_Tp_Pq_qErrest0);         
            Pqwr_qtz.solve(packed_Tp_Pq_qErrest0, dfq_WLS); // weighted
            unpackFreeQ(s, dfq_WLS, udfq_WLS); // zeroes in q_p slots
            multiplyByNInv(s,false,udfq_WLS,du);
        } else {
            multiplyByPq(s, bias_p, qErrest, Tp_Pq_qErrest0); // Pq*qErrest
            Tp_Pq_qErrest0.rowScaleInPlace(perrWeights); // now Tp*Pq*qErrest
            packActiveMultipliers(enforcedPosConsEqns,
                                  Tp_Pq_qErrest0,
                                  packed_Tp_Pq_qErrest0);         
            Pqwr_qtz.solve(packed_Tp_Pq_qErrest0, dfq_WLS); // weighted
            multiplyByNInv(s,false,dfq_WLS,du);
        }
        // Here du = du_WLS = N\ dq_WLS.
        du.rowScaleInPlace(WuInv); // now du = Wu\ du_WLS
        multiplyByN(s,false,du,dq);    // dq = N du
        qErrest -= dq; // unweighted
    }

    //cout << "!!!! perr TRMS achieved " << normAchievedTRMS << " in " 
    //     << nItsUsed << " iterations"  << endl;

    /* By design, normalization of quaternions can't have any effect on the 
    constraints we just fixed (because we normalize internally for
    calculations). So now we can simply normalize the quaternions. We can't 
    touch any q's that are prescribed, though, so it is possible that we'll fail
    to achieve the required tolerance if some quaterion is prescribed but not up
    to date. */
    if (mQuats) {
        const bool anyQuatChange = normalizeQuaternions(s,qErrest);
        if (anyQuatChange) results.setAnyChangeMade(true);
        quatNormAchieved = useNormInf ? quatErrs.normInf(&worstQuatErr)
                                      : quatErrs.normRMS(&worstQuatErr);
        if (quatNormAchieved > consAccuracy) {
            results.setNormOnExit(quatNormAchieved);
            results.setExitStatus(ProjectResults::FailedToAchieveAccuracy);
            if (!dontThrow) {
                SimTK_ERRCHK2_ALWAYS(quatNormAchieved <= consAccuracy, 
                        "SimbodyMatterSubsystem::projectQ()",
                        "Failed to normalize quaternions. Norm achieved=%g"
                        " but required norm=%g. Did you forget to call"
                        " prescribeQ()?", quatNormAchieved, consAccuracy);
            }
            return 1;
        }
    }
    results.setNormOnExit(std::max(perrNormAchieved, quatNormAchieved));
    results.setExitStatus(ProjectResults::Succeeded);
    return 0;
}
//................................. PROJECT Q ..................................



//==============================================================================
//                         NORMALIZE QUATERNIONS
//==============================================================================
/* Project quaternions onto their constraint manifold by normalizing them. Also 
removes any error component along the length of the quaternions if given a 
qErrest. Returns true if any change was made. */
bool SimbodyMatterSubsystemRep::normalizeQuaternions
   (State& s, Vector& qErrest) const
{
    SBStateDigest sbs(s, *this, Stage::Model);
    const SBInstanceCache& ic = getInstanceCache(s);

    Vector& q  = updQ(s); // invalidates q's. TODO: see below.

    bool anyChangeMade = false;
    for (int i=0; i<(int)rbNodeLevels.size(); ++i) 
        for (int j=0; j<(int)rbNodeLevels[i].size(); ++j) { 
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            const SBInstancePerMobodInfo& mobodInfo =
                ic.mobodInstanceInfo[node.getNodeNum()];
            if (mobodInfo.qMethod != Motion::Free)
                continue;

            if (node.enforceQuaternionConstraints(sbs,q,qErrest))
                anyChangeMade = true;
        }

    // This will recalculate the qnorms (all 1), qerrs (all 0). The only
    // other quaternion dependency is the N matrix (and NInv, NDot).
    // TODO: better if these updates could be made without invalidating
    // Position stage in general. I *think* N is always calculated on the fly.
    realizeSubsystemPosition(s);
    return anyChangeMade;
}
//.......................... NORMALIZE QUATERNIONS .............................



//==============================================================================
//                                 PROJECT U
//==============================================================================
int SimbodyMatterSubsystemRep::projectU
   (State&                  s, 
    Vector&                 uErrest,        // u error estimate or empty 
    const ProjectOptions&   opts, 
    ProjectResults&         results) const
{
    SimTK_STAGECHECK_GE(getStage(s), Stage::Velocity,
        "SimbodyMatterSubsystemRep::projectU()");

    results.clear();
    const Real consAccuracy = opts.getRequiredAccuracy();
    // Normally we'll use an RMS norm for the perrs.
    const bool useNormInf = opts.isOptionSet(ProjectOptions::UseInfinityNorm);
    // Force projection even if accuracy ok on entry.
    const bool forceOneIter = opts.isOptionSet(ProjectOptions::ForceProjection);
    // Normally we'll throw an exception with a helpful message. If this is
    // set we'll quietly return status instead.
    const bool dontThrow = opts.isOptionSet(ProjectOptions::DontThrow);

    // Here we deal with the nonholonomic (velocity) constraints and the 
    // derivatives of the holonomic constraints.
    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);

    // These arrays contain a unique MultiplierIndex for each perrdot and verr 
    // that is to be enforced by projection. Baumgarte-stabilized constraint
    // equations don't count since then position and velocity constraints are 
    // softly enforced through the acceleration error equation instead.
    const Array_<MultiplierIndex> enforcedPosConsEqns =
                                            findEnforcedPosConstraintEqns(s);
    const Array_<MultiplierIndex> enforcedVelConsEqns =
                                            findEnforcedVelConstraintEqns(s);

    const int mfp = (int)enforcedPosConsEqns.size();
    const int mfv = (int)enforcedVelConsEqns.size();
    const int mfpv = mfp + mfv;
    
    // This is a const view into the State; the contents it refers to will 
    // change though. These are all the velocity-level constraint errors,
    // mHolo from differentiating holonomic constraints and mNonholo directly
    // from the nonholonomic constraints.
    const Vector& pvErrs = getUErr(s); // mHolo+mNonholo of these
    const Vector& pverrWeights = getUErrWeights(s); // 1/unit err (Tpv)

    // Determine norm on entry.
    int worstPVerr;
    Vector scaledPVerrs(mfpv); // scaled and packed
    packActiveMultipliers2(enforcedPosConsEqns, enforcedVelConsEqns,
                           pvErrs.rowScale(pverrWeights),
                           scaledPVerrs);
    const Real pverrNormOnEntry = useNormInf ? scaledPVerrs.normInf(&worstPVerr)
                                             : scaledPVerrs.normRMS(&worstPVerr);
    // Map worstPVerr back to real index.
    if (worstPVerr >= 0) {// will be -1 if mfpv==0
        if (worstPVerr < mfp)
            worstPVerr = enforcedPosConsEqns[worstPVerr];  
        else 
            worstPVerr = enforcedVelConsEqns[worstPVerr - mfp];  
    }
    results.setNormOnEntrance(pverrNormOnEntry, worstPVerr);

    if (pverrNormOnEntry > opts.getProjectionLimit()) {
        results.setProjectionLimitExceeded(true);
        results.setExitStatus(ProjectResults::FailedToConverge);
        return 1;
    }
    
    //cout << "!!!! initially @" << s.getTime() << ", verr TRMS=" 
    //     << normAchievedTRMS << " consAcc=" << consAccuracy;
    //if (uErrest.size())
    //    cout << " uErrest WRMS=" << uErrest.rowScale(Wu).normRMS();
    //else cout << " NO U ERROR ESTIMATE";
    //cout << endl;
    

    // Return quickly if (a) constraint norm is zero (probably because there
    // aren't any), or (b) constraints are already satisfied and we're
    // being forced to go ahead anyway.
    if (    pverrNormOnEntry == 0
        || (pverrNormOnEntry <= consAccuracy && !forceOneIter)) {
        // numIterations==0.
        results.setAnyChangeMade(false);
        results.setNormOnExit(pverrNormOnEntry);
        results.setExitStatus(ProjectResults::Succeeded);
        return 0;
    }

    // We're going to have to project constraints. Get the remaining options.

    // This is the factor by which we try to achieve a tighter accuracy
    // than requested. E.g. if overshootFactor=0.1 then we attempt 10X 
    // tighter accuracy if we can get it. But we won't fail as long as
    // we manage to reach consAccuracy.
    const Real overshootFactor = opts.getOvershootFactor();
    const Real consAccuracyToTryFor = 
        std::max(overshootFactor*consAccuracy, SignificantReal);

    // Check whether we should stop if we see the solution diverging
    // which should not happen when we're in the neighborhood of a solution
    // on entry. This is always set while integrating, except during
    // initialization.
    const bool localOnly = opts.isOptionSet(ProjectOptions::LocalOnly);
    // We are permitted to use an out-of-date Jacobian for projection unless
    // this is set. TODO: always using modified Newton at the moment.
    const bool forceFullNewton =
        opts.isOptionSet(ProjectOptions::ForceFullNewton);

    // Get problem dimensions.
    const SBInstanceCache& ic = getInstanceCache(s);
    const int nq       = getNQ(s);
    const int nu       = getNU(s);
    const int nfu      = ic.getTotalNumFreeU();
    bool hasPrescribedMotion = (nfu != nu);

    // Solve 
    //   (Tpv [P;V] Wu^-1) du_WLS  = Tpv uerr
    //                         du  = Wu^-1*du_WLS
    //                          u -= du
    // Note that although this is a nonlinear least squares problem since uerr 
    // is a function of u, we do not need to refactor the matrix since it does 
    // not depend on u.
    // TODO: I don't think that's true -- V can depend on u (rarely). That
    // doesn't mean we need to refactor it, but then this is a modified Newton
    // iteration (rather than full) if we're not updating V when we could be.
    //
    // This is a nonlinear least squares problem, but we only need to factor 
    // once since only the RHS is dependent on u (TODO: see above).

    // This will be updated as we go.
    Real pverrNormAchieved = pverrNormOnEntry;


    // Calculate relative scaling for changes to u.
    const Vector& u = getU(s);
    const Vector& Wu = getUWeights(s); // 1/unit change (Wu)
    Vector uRelScale(nu);
    for (int i=0; i<nu; ++i) {
        const Real ui = std::abs(u[i]);
        const Real wi = Wu[i];
        uRelScale[i] = ui*wi > 1 ? ui : 1/wi; // max(unit error, u) (1/Eu)
    }

    Real lastChangeMadeWRMS = 0;
    int nItsUsed = 0;

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is much too tight; should depend on constraint tolerance
    // and should be consistent with the holonomic rank.
    const Real conditioningTol = (mHolo+mNonholo) 
        //* SignificantReal;
        * SqrtEps;

    // Keep the starting u in case we have to restore it, which we'll do
    // if the attempts here make the constraint norm worse.
    const Vector saveU = getU(s);

    Matrix PVwrt(nfu, mfpv);
    Vector dfu_WLS(nfu);
    Vector du(nu); // unpacked into here if necessary
    if (hasPrescribedMotion)
        du.setToZero(); // must initialize unwritten elements

    calcWeightedPVrTranspose(s, enforcedPosConsEqns, enforcedVelConsEqns,
                             pverrWeights, uRelScale, PVwrt);
    // PVwrt is now Eu^-1 (Pt Vt) Tpv

    // Calculate pseudoinverse (just once)
    FactorQTZ PVwr_qtz;
    PVwr_qtz.factor<Real>(~PVwrt, conditioningTol);

    //printf("projectU m=%d condTol=%g rank=%d rcond=%g\n",
    //    PVwrt.ncol(), conditioningTol, PVwr_qtz.getRank(),
    //    PVwr_qtz.getRCondEstimate());

    Real prevPVerrNormAchieved = pverrNormAchieved; // watch for divergence
    bool diverged = false;
    const int MaxIterations  = 7;
    do {
        PVwr_qtz.solve(scaledPVerrs, dfu_WLS);
        lastChangeMadeWRMS = dfu_WLS.normRMS(); // change in weighted norm

        // switch back to unweighted du=Eu^-1*du_WLS
        if (hasPrescribedMotion) {
            unpackFreeU(s, dfu_WLS, du);    // zeroes in u_p slots
            du.rowScaleInPlace(uRelScale); // du=Eu^-1*unpack(dfu_WLS)
        } else {
            du = dfu_WLS.rowScale(uRelScale); // unscale: du=Eu^-1*du_WLS
        }
        updU(s) -= du;
        results.setAnyChangeMade(true);

        // Recalculate the constraint errors for the new u's.
        realizeSubsystemVelocity(s);

        packActiveMultipliers2(enforcedPosConsEqns, enforcedVelConsEqns,
                               pvErrs.rowScale(pverrWeights), // Tpv * pvErrs
                               scaledPVerrs);
        pverrNormAchieved = useNormInf ? scaledPVerrs.normInf()
                                       : scaledPVerrs.normRMS();
        ++nItsUsed;

        if (localOnly && nItsUsed >= 2 
            && pverrNormAchieved > prevPVerrNormAchieved) {
            // Velocity norm worse -- restore to end of previous iteration.
            updU(s) += du;
            realizeSubsystemVelocity(s); // pvErrs changes here
            packActiveMultipliers2(enforcedPosConsEqns, enforcedVelConsEqns,
                                   pvErrs.rowScale(pverrWeights),
                                   scaledPVerrs);
            pverrNormAchieved = useNormInf ? scaledPVerrs.normInf()
                                           : scaledPVerrs.normRMS();
            diverged = true;
            break; // diverging -- quit now to prevent a bad solution
        }

        prevPVerrNormAchieved = pverrNormAchieved;

    } while (pverrNormAchieved > consAccuracyToTryFor
                && nItsUsed < MaxIterations);

    results.setNumIterations(nItsUsed);

    // Make sure we achieved at least the required constraint accuracy. If not 
    // we'll return with an error. If we see that the norm has been made worse
    // than it was on entry, we'll restore the state to what it was on entry. 
    // Otherwise we'll return with the improved-but-not-good-enough result.
    if (pverrNormAchieved > consAccuracy) {
        if (pverrNormAchieved >= pverrNormOnEntry) { // made it worse
            updU(s) = saveU; // revert
            realizeSubsystemVelocity(s);
            pverrNormAchieved = pverrNormOnEntry;
        }
     
        results.setNormOnExit(pverrNormAchieved);

        if (diverged) {
            results.setExitStatus(ProjectResults::FailedToConverge);
            if (!dontThrow) {
                SimTK_ERRCHK_ALWAYS(!diverged,
                    "SimbodyMatterSubsystem::projectU()",
                    "Attempt to project constraints locally diverged.");
            }
        } else {
            results.setExitStatus(ProjectResults::FailedToAchieveAccuracy);
            if (!dontThrow) {
                SimTK_ERRCHK3_ALWAYS(pverrNormAchieved <= consAccuracy, 
                    "SimbodyMatterSubsystem::projectU()",
                    "Failed to achieve required accuracy %g. Norm on entry "
                    " was %g; norm on exit %g. You might need a better"
                    " starting configuration, or if there are prescribed or "
                    " locked u's you might have to free some of them.", 
                    consAccuracy, pverrNormOnEntry, pverrNormAchieved);
            }
        }

        return 1;
    }

    // Velocity constraint errors were successfully driven to consAccuracy.

    /* Next, if we projected out the velocity constraint errors, remove the
    corresponding error from the integrator's error estimate.
    
      (Tpv [P;V] Wu^-1)_f dfu_WLS  = (Tpv [P;V] Wu^-1)_f (Wu uErrest)_f
                           du_WLS  = unpack(dfu_WLS) (with zero fill)
                            du  = Wu^-1 du_WLS
                       uErrest -= du
    No iteration is required. Note that we only include the enforced constraints
    that we used above, not all enabled constraint equations.
    
    We can simplify the RHS of the first equation above:
      (Tpv [P;V] Wu^-1)_f (Wu uErrest)_f = Tpv [P;V] unpack(uErrest_f)
    for which we have an O(n) operator to compute the matrix-vector 
    product. */

    if (uErrest.size()) {
        // Work in Wu-norm
        Vector Tpv_PV_uErrest0, packed_Tpv_PV_uErrest0, bias_pv;
        calcBiasForMultiplyByPVA(s,true,true,false,bias_pv); // just P,V
        if (hasPrescribedMotion) {
            Vector uErrest_0(uErrest);
            zeroKnownU(s, uErrest_0); // zero out prescribed entries
            multiplyByPVA(s,true,true,false,bias_pv,
                            uErrest_0,Tpv_PV_uErrest0);
            Tpv_PV_uErrest0.rowScaleInPlace(pverrWeights); // = Tpv*PV*uErrest_0
            packActiveMultipliers2(enforcedPosConsEqns,enforcedVelConsEqns,
                                   Tpv_PV_uErrest0,
                                   packed_Tpv_PV_uErrest0);
            PVwr_qtz.solve(packed_Tpv_PV_uErrest0, dfu_WLS);
            unpackFreeU(s, dfu_WLS, du); // still weighted
        } else {
            multiplyByPVA(s,true,true,false,bias_pv,uErrest,Tpv_PV_uErrest0);
            Tpv_PV_uErrest0.rowScaleInPlace(pverrWeights); // = Tpv PV uErrEst
            packActiveMultipliers2(enforcedPosConsEqns,enforcedVelConsEqns,
                                   Tpv_PV_uErrest0,
                                   packed_Tpv_PV_uErrest0);
            PVwr_qtz.solve(packed_Tpv_PV_uErrest0, du);
        }
        du.rowScaleInPlace(uRelScale); // now du=Eu^-1*unpack(dfu_WLS)
        uErrest -= du; // this is unweighted now
    }
   
    //cout << "!!!! verr achieved " << pverrNormAchieved << " in " 
    //     << nItsUsed << " iterations" << endl;
    //if (uErrest.size())
    //    cout << " uErrest WRMS=" << uErrest.rowScale(Wu).normRMS() << endl;

    results.setNormOnExit(pverrNormAchieved);
    results.setExitStatus(ProjectResults::Succeeded);
    return 0;
}
//................................ PROJECT U ...................................



//==============================================================================
//                     CALC TREE FORWARD DYNAMICS OPERATOR
//==============================================================================
// Given a State realized through Stage::Dynamics, and a complete set of applied 
// forces, calculate all acceleration results into the return arguments here. 
// This routine *does not* affect the State cache -- it is an operator. In 
// typical usage, the output arguments actually will be part of the state cache 
// to effect a response, but this method can also be used to effect an operator.
//
// Note that although acceleration constraint errors will be calculated, the 
// returned accelerations will not obey the constraints, unless the supplied 
// forces already account for constraints. The argument list allows for some 
// extra forces to be supplied, with the intent that these will be used to deal 
// with internal forces generated by constraints. Note that the extra forces 
// here are treated with opposite sign from the applied forces, as is 
// appropriate for constraint forces.
void SimbodyMatterSubsystemRep::calcTreeForwardDynamicsOperator(
    const State&                    s,
    const Vector&                   mobilityForces,
    const Vector_<Vec3>&            particleForces,
    const Vector_<SpatialVec>&      bodyForces,
    const Vector*                   extraMobilityForces,
    const Vector_<SpatialVec>*      extraBodyForces,
    SBTreeAccelerationCache&        tac,  // accels, prescribed forces go here
    Vector&                         udot, // in/out (in for prescribed udot)
    Vector&                         qdotdot,
    Vector&                         udotErr) const
{
    SBStateDigest sbs(s, *this, Stage::Acceleration);

    const SBModelCache&         mc  = sbs.getModelCache();
    const SBInstanceCache&      ic  = sbs.getInstanceCache();
    const SBTreePositionCache&  tpc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  tvc = sbs.getTreeVelocityCache();
    const SBDynamicsCache&      dc  = sbs.getDynamicsCache();

    // Ensure that output arguments have been allocated properly.
    tac.allocate(topologyCache, mc, ic);
    udot.resize(topologyCache.nDOFs);
    qdotdot.resize(topologyCache.maxNQs);

    Vector              totalMobilityForces;
    Vector_<SpatialVec> totalBodyForces;

    // inputs

    // First assume we'll use the input forces as-is.
    const Vector*              mobilityForcesToUse  = &mobilityForces;
    const Vector_<SpatialVec>* bodyForcesToUse      = &bodyForces;

    if (extraMobilityForces) {
        totalMobilityForces = mobilityForces - *extraMobilityForces; // note sign
        mobilityForcesToUse = &totalMobilityForces;
    }

    if (extraBodyForces) {
        totalBodyForces = bodyForces - *extraBodyForces;    // note sign
        bodyForcesToUse = &totalBodyForces;
    }

    // outputs
    Vector&              netHingeForces = tac.epsilon;
    Array_<SpatialVec,MobilizedBodyIndex>&
                         abForcesZ      = tac.z;
    Array_<SpatialVec,MobilizedBodyIndex>&
                         abForcesZPlus  = tac.zPlus;
    Vector_<SpatialVec>& A_GB           = tac.bodyAccelerationInGround;
    Vector&              tau            = tac.presMotionForces;

    // Calculate accelerations produced by these forces in three forms:
    // body accelerations A_GB, u-space generalized accelerations udot,
    // and q-space generalized accelerations qdotdot.
    calcTreeAccelerations
       (s, *mobilityForcesToUse, *bodyForcesToUse, dc.presUDotPool,
        netHingeForces, abForcesZ, abForcesZPlus,
        A_GB, udot, qdotdot, tau);

    // Feed the accelerations into the constraint error methods to determine
    // the acceleration constraint errors they generate.
    calcConstraintAccelerationErrors(s, A_GB, udot, qdotdot, udotErr);
}
//......................CALC TREE FORWARD DYNAMICS OPERATOR ....................



//==============================================================================
//                        REALIZE TREE FORWARD DYNAMICS
//==============================================================================
// This is the response version of the above operator; that is, it uses
// the operator but puts the results in the State cache. Note that this
// only makes sense if the force arguments also come from the State
// somewhere else in the System that includes this Subsystem.
void SimbodyMatterSubsystemRep::realizeTreeForwardDynamics(
    const State&               s,
    const Vector&              mobilityForces,
    const Vector_<Vec3>&       particleForces,
    const Vector_<SpatialVec>& bodyForces,
    const Vector*              extraMobilityForces,
    const Vector_<SpatialVec>* extraBodyForces) const
{
    // Output goes into State's global cache and our AccelerationCache.
    SBTreeAccelerationCache&        tac     = updTreeAccelerationCache(s);
    Vector&                         udot    = updUDot(s);
    Vector&                         qdotdot = updQDotDot(s);
    Vector&                         udotErr = updUDotErr(s);

    calcTreeForwardDynamicsOperator
       (s, mobilityForces, particleForces, bodyForces,
        extraMobilityForces, extraBodyForces,
        tac, udot, qdotdot, udotErr);

    // Since we're realizing, mark the resulting cache entry valid.
    markCacheValueRealized(s, topologyCache.treeAccelerationCacheIndex);
}
//....................... REALIZE TREE FORWARD DYNAMICS ........................



//==============================================================================
//                    CALC LOOP FORWARD DYNAMICS OPERATOR
//==============================================================================
// Given a State realized through Stage::Dynamics, and a complete set of applied 
// forces, calculate all acceleration results resulting from those forces AND 
// enforcement of the acceleration constraints. The results go into the return 
// arguments here. This routine *does not* affect the State cache (except for
// lazy realization of ABIs) -- it is an operator. In typical usage, the 
// output arguments actually will be part of the state cache to effect a 
// response, but this method can also be used to effect an operator.
void SimbodyMatterSubsystemRep::calcLoopForwardDynamicsOperator
   (const State& s, 
    const Vector&                   mobilityForces,
    const Vector_<Vec3>&            particleForces,
    const Vector_<SpatialVec>&      bodyForces,
    SBTreeAccelerationCache&        tac,
    SBConstrainedAccelerationCache& cac,
    Vector&                         udot,
    Vector&                         qdotdot,
    Vector&                         multipliers,
    Vector&                         udotErr) const
{
    assert(getStage(s) >= Stage::Acceleration-1);

    // Calculate acceleration results ignoring Constraints, except to have
    // them calculate the resulting constraint errors.
    calcTreeForwardDynamicsOperator
       (s, mobilityForces, particleForces, bodyForces,
        nullptr, nullptr, tac, udot, qdotdot, udotErr);

    // Next, determine how many acceleration-level constraint equations 
    // need to be obeyed.

    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = getNumAccelerationOnlyConstraintEquationsInUse(s);
    const int m        = mHolo+mNonholo+mAccOnly;
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    multipliers.resize(m);
    multipliers.setToZero();
    if (m==0 || nu==0)
        return;

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is probably too tight; should depend on constraint tolerance
    // and should be consistent with position and velocity projection ranks.
    // Tricky here because conditioning depends on mass matrix as well as
    // constraints.
    const Real conditioningTol = m * Eps34;

    // Calculate multipliers lambda as
    //     C lambda = aerr
    // where C=G M^-1 ~G is the constraint compliance matrix. 
    //
    // The method here calculates the mXm matrix G*M^-1*G^T as fast as 
    // I know how to do, O(m*n) with O(n) temporary memory, using a series
    // of O(n) operators. Then we'll factor it here in O(m^3) time. 
    //
    // TODO: Currently computing C for all enabled constraints and then select
    // only the active ones afterwards. Better just to compute the active
    // subset in the first place.
    Matrix GMInvGt(m,m);
    calcGMInvGt(s, GMInvGt); // All enabled constraints.

    Array_<MultiplierIndex> active = findActiveMultipliers(s);
    const int mActive = (int)active.size();
    if (mActive == 0) {
        SimTK_DEBUG2(
        "calcFwdDyn @%g: None of the %d enabled constraint equations was active.\n",
            s.getTime(), m);
        return; // all multipliers are zero; accelerations unconstrained
    }

    // At least one constraint is active.

    #ifndef NDEBUG
    printf("calcFwdDyn @%g: %d constraints active out of %d enabled\n", 
           s.getTime(), mActive, m);
    cout << "  Active: " << active << endl;
    #endif

    // Here we don't care about the particular constraint status (rolling,
    // sliding, impending slip); just whether it is active at all. If so we
    // generate C for the normal and rolling constraints.
    Matrix C(mActive,mActive);
    formActiveSubmatrix(active, GMInvGt, C);

    // TODO: Slip/impending slip are dealt with by replacing rows before
    // factoring.
    
    // specify 1/cond at which we declare rank deficiency
    FactorQTZ qtz(C, conditioningTol); 

    //printf("fwdDynamics: m=%d condTol=%g rank=%d rcond=%g\n",
    //    GMInvGt.nrow(), conditioningTol, qtz.getRank(),
    //    qtz.getRCondEstimate());

    Vector packedUdotErr(mActive), activeMultipliers(mActive);
    packActiveMultipliers(active, udotErr, packedUdotErr);
    qtz.solve(packedUdotErr, activeMultipliers);

    // Unpack active multipliers into the full size multiplier array (which was
    // initialized to zero above).
    unpackActiveMultipliers(active, activeMultipliers, multipliers);

    // We have the multipliers, now turn them into forces.

    Vector_<SpatialVec> bodyForcesInG;
    Vector              mobilityF;
    calcConstraintForcesFromMultipliers(s,multipliers,bodyForcesInG,mobilityF,
        cac.constrainedBodyForcesInG, cac.constraintMobilityForces);
    // Note that constraint forces have the opposite sign from applied forces
    // so must be subtracted to calculate the total forces.

    // Recalculate the accelerations applying the constraint forces in addition
    // to the applied forces that were passed in. The constraint errors 
    // calculated now should be within numerical noise of zero.
    calcTreeForwardDynamicsOperator
       (s, mobilityForces, particleForces, bodyForces,
        &mobilityF, &bodyForcesInG, tac, udot, qdotdot, udotErr);
}
//................... CALC LOOP FORWARD DYNAMICS OPERATOR ......................



//==============================================================================
//                       REALIZE LOOP FORWARD DYNAMICS
//==============================================================================
// Given the set of forces in the state, calculate accelerations resulting from
// those forces and enforcement of acceleration constraints.
void SimbodyMatterSubsystemRep::realizeLoopForwardDynamics(const State& s, 
    const Vector&               mobilityForces,
    const Vector_<Vec3>&        particleForces,
    const Vector_<SpatialVec>&  bodyForces) const 
{
    // Because we are realizing, we want to direct the output of the operator
    // back into the State cache.
    SBTreeAccelerationCache&        tac         = updTreeAccelerationCache(s);
    SBConstrainedAccelerationCache& cac         = updConstrainedAccelerationCache(s);
    Vector&                         udot        = updUDot(s);
    Vector&                         qdotdot     = updQDotDot(s);
    Vector&                         udotErr     = updUDotErr(s);
    Vector&                         multipliers = updMultipliers(s);

    calcLoopForwardDynamicsOperator
       (s, mobilityForces, particleForces, bodyForces,
        tac, cac, udot, qdotdot, multipliers, udotErr);

    // Since we're realizing, note that we're done with these cache entries.
    markCacheValueRealized(s, topologyCache.treeAccelerationCacheIndex);
    markCacheValueRealized(s, topologyCache.constrainedAccelerationCacheIndex);
}
//....................... REALIZE LOOP FORWARD DYNAMICS ........................



// =============================================================================
//                        CALC COMPOSITE BODY INERTIAS
// =============================================================================
// Given a State realized to Position stage, calculate the composite
// body inertias seen by each mobilizer. A composite body inertia is
// the inertia of the rigid body created by locking all joints outboard
// of a particular mobilized body. (Constraints have no effect on the result.)
//
void SimbodyMatterSubsystemRep::calcCompositeBodyInertias(const State& s,
    Array_<SpatialInertia,MobilizedBodyIndex>& R) const
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    R.resize(getNumBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcCompositeBodyInertiasInward(tpc,R);
}
//....................... CALC COMPOSITE BODY INERTIAS .........................



// =============================================================================
//                                  REALIZE Y
// =============================================================================
// Y is used for length constraints: sweep from base to tip. You can call this
// after Position stage but it may have to realize articulated bodies first.
void SimbodyMatterSubsystemRep::realizeY(const State& s) const {
    realizeArticulatedBodyInertias(s);

    const SBInstanceCache&                ic  = getInstanceCache(s);
    const SBTreePositionCache&            tpc = getTreePositionCache(s);
    const SBArticulatedBodyInertiaCache&  abc = getArticulatedBodyInertiaCache(s);
    SBDynamicsCache&                      dc  = updDynamicsCache(s);

    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->realizeYOutward(ic,tpc,abc,dc);
}
//.................................. REALIZE Y .................................



//==============================================================================
//                           CALC KINETIC ENERGY
//==============================================================================
Real SimbodyMatterSubsystemRep::calcKineticEnergy(const State& s) const {
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const SBTreeVelocityCache& tvc = getTreeVelocityCache(s);

    Real ke = 0;

    // Skip ground level 0!
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            ke += rbNodeLevels[i][j]->calcKineticEnergy(tpc,tvc);

    return ke;
}



//==============================================================================
//                          CALC TREE ACCELERATIONS
//==============================================================================
// Operator for open-loop forward dynamics.
// This Subsystem must have already realized VelocityKinematics so that 
// Coriolis terms are available, and articulated body inertias and articulated
// body velocities are realized here if necessary. All vectors must use 
// contiguous storage.
void SimbodyMatterSubsystemRep::calcTreeAccelerations(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    const Array_<Real>&        presUDots, // packed
    Vector&                    netHingeForces,
    Array_<SpatialVec,MobilizedBodyIndex>& allZ, 
    Array_<SpatialVec,MobilizedBodyIndex>& allZPlus, 
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot,    // in/out (in for prescribed udots)
    Vector&                    qdotdot,
    Vector&                    tau) const 
{
    // Note that realize(Acceleration) depends on getting here to fulfill the
    // promise of these cache entries' computed-by stage.
    realizeArticulatedBodyInertias(s); // might already be done
    realizeArticulatedBodyVelocity(s);

    const SBArticulatedBodyInertiaCache& abc = 
        getArticulatedBodyInertiaCache(s);
    const SBArticulatedBodyVelocityCache& abvc = 
        getArticulatedBodyVelocityCache(s);

    SBStateDigest sbs(s, *this, Stage::Acceleration);

    const SBInstanceCache&      ic  = sbs.getInstanceCache();
    const SBTreePositionCache&  tpc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  tvc = sbs.getTreeVelocityCache();
    const SBDynamicsCache&      dc  = sbs.getDynamicsCache();

    assert(mobilityForces.size() == getTotalDOF());
    assert(bodyForces.size() == getNumBodies());

    netHingeForces.resize(getTotalDOF());
    allZ.resize(getNumBodies());
    allZPlus.resize(getNumBodies());
    A_GB.resize(getNumBodies());
    udot.resize(getTotalDOF());
    qdotdot.resize(getTotalQAlloc());
    tau.resize(ic.getTotalNumPresForces());

    assert(mobilityForces.hasContiguousData());
    assert(bodyForces.hasContiguousData());
    assert(netHingeForces.hasContiguousData());
    assert(A_GB.hasContiguousData());
    assert(udot.hasContiguousData());
    assert(qdotdot.hasContiguousData());
    assert(tau.hasContiguousData());

    const Real*       mobilityForcePtr = mobilityForces.size() 
                                            ? &mobilityForces[0] : nullptr;
    const SpatialVec* bodyForcePtr     = bodyForces.size() 
                                            ? &bodyForces[0] : nullptr;
    Real*             hingeForcePtr    = netHingeForces.size() 
                                            ? &netHingeForces[0] : nullptr;
    SpatialVec*       aPtr             = A_GB.size()    ? &A_GB[0] : nullptr;
    Real*             udotPtr          = udot.size()    ? &udot[0] : nullptr;
    Real*             qdotdotPtr       = qdotdot.size() ? &qdotdot[0] : nullptr;
    Real*             tauPtr           = tau.size()     ? &tau[0] : nullptr;
    SpatialVec*       zPtr             = allZ.begin();    
    SpatialVec*       zPlusPtr         = allZPlus.begin(); 

    // If there are any prescribed udots, scatter them into the appropriate
    // udot entries now. We must also set known-zero udots to zero here.
    assert(presUDots.size() == ic.getTotalNumPresUDot());
    for (PresUDotPoolIndex i(0); i < presUDots.size(); ++i)
        udotPtr[ic.presUDot[i]] = presUDots[i];
    for (int i=0; i < (int)ic.zeroUDot.size(); ++i)
        udotPtr[ic.zeroUDot[i]] = 0;

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass1Inward(ic,tpc,abc,abvc,
                mobilityForcePtr, bodyForcePtr, udotPtr, zPtr, zPlusPtr,
                hingeForcePtr);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass2Outward(ic,tpc,abc,tvc,dc, 
                hingeForcePtr, aPtr, udotPtr, tauPtr);
            node.calcQDotDot(sbs, &udotPtr[node.getUIndex()], 
                             &qdotdotPtr[node.getQIndex()]);
        }
}
//......................... CALC TREE ACCELERATIONS ............................



//==============================================================================
//                            MULTIPLY BY M INV
//==============================================================================
// Calculate udot = M^-1 f. We also get spatial accelerations A_GB for 
// each body as a side effect.
// This Subsystem must already be realized through Position stage or at
// least have PositionKinematics already available; we'll 
// realize articulated body inertias here if necessary.
// All vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyByMInv(const State& s,
    const Vector&                                           f,
    Vector&                                                 MInvf) const 
{
    const SBInstanceCache&                  ic  = getInstanceCache(s);
    const SBTreePositionCache&              tpc = getTreePositionCache(s);

    realizeArticulatedBodyInertias(s); // (may already have been realized)
    const SBArticulatedBodyInertiaCache&    abc = getArticulatedBodyInertiaCache(s);

    const int nb = getNumBodies();
    const int nu = getNU(s);

    assert(f.size() == nu);

    MInvf.resize(nu);
    if (nu==0)
        return;

    assert(f.hasContiguousData());
    assert(MInvf.hasContiguousData());

    // Temporaries
    Array_<Real>        eps(nu);
    Array_<SpatialVec>  z(nb), zPlus(nb), A_GB(nb);

    // Point to raw data of input arguments.
    const Real* fPtr     = &f[0];       
    Real*       MInvfPtr = &MInvf[0];

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyByMInvPass1Inward(ic,tpc,abc,
                fPtr, z.begin(), zPlus.begin(), eps.begin());
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyByMInvPass2Outward(ic,tpc,abc, 
                eps.cbegin(), A_GB.begin(), MInvfPtr);
        }
}
//............................. CALC M INVERSE F ...............................



//==============================================================================
//                              MULTIPLY BY M
//==============================================================================
// Calculate f = M a.
// This Subsystem must already have been realized to Position stage.
// All vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyByM(const State&    s,
                                            const Vector&   a,
                                            Vector&         Ma) const 
{

    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const int nb = getNumBodies();
    const int nu = getNU(s);

    assert(a.size() == nu);
    Ma.resize(nu);

    if (nu == 0)
        return;

    assert(a.hasContiguousData());
    assert(Ma.hasContiguousData());

    // Temporaries
    Array_<SpatialVec>  fTmp(nb), A_GB(nb);

    // Point to raw data of input arguments.
    const Real* aPtr    = &a[0];       
    Real*       MaPtr   = &Ma[0];

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyByMPass1Outward(tpc, aPtr, A_GB.begin());
        }

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyByMPass2Inward(tpc,A_GB.cbegin(),fTmp.begin(),MaPtr);
        }
}



//==============================================================================
//                                  CALC M
//==============================================================================
// Calculate the mass matrix M in O(n^2) time. This Subsystem must already have
// been realized to Position stage.
// It is OK if M's data is not contiguous.
void SimbodyMatterSubsystemRep::calcM(const State& s, Matrix& M) const {
    const int nu = getTotalDOF();
    M.resize(nu,nu);
    if (nu==0) return;

    // This could be calculated much faster by doing it directly and calculating
    // only half of it. As a placeholder, however, we're doing this with 
    // repeated O(n) calls to multiplyByM() to get M one column at a time.

    // If M's columns are contiguous we can avoid copying.
    const bool isContiguous = M(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : nu);

    Vector v(nu); v.setToZero();
    for (int i=0; i < nu; ++i) {
        v[i] = 1;
        if (isContiguous) {
            multiplyByM(s, v, M(i));
        } else {
            multiplyByM(s, v, contig_col);
            M(i) = contig_col;
        }
        v[i] = 0;
    }
}



//==============================================================================
//                                CALC MInv
//==============================================================================
// Calculate the mass matrix inverse MInv(=M^-1) in O(n^2) time. This Subsystem
// must already have been realized to Position stage.
// It is OK if MInv's data is not contiguous.
void SimbodyMatterSubsystemRep::calcMInv(const State& s, Matrix& MInv) const {
    const int nu = getTotalDOF();
    MInv.resize(nu,nu);
    if (nu==0) return;

    // This could probably be calculated faster by doing it directly and
    // filling in only half. For now we're doing it with repeated calls to
    // the O(n) operator multiplyByMInv().

    // If M's columns are contiguous we can avoid copying.
    const bool isContiguous = MInv(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : nu);

    Vector f(nu); f.setToZero();
    for (int i=0; i < nu; ++i) {
        f[i] = 1;
        if (isContiguous) {
            multiplyByMInv(s, f, MInv(i));
        } else {
            multiplyByMInv(s, f, contig_col);
            MInv(i) = contig_col;
        }
        f[i] = 0;
    }
}



//==============================================================================
//                          CALC TREE RESIDUAL FORCES
//==============================================================================
// Operator for tree system inverse dynamics. 
// Note that this includes the effects of inertial forces.
// This calculates
//      f_resid = M(q) udot + f_inertial(q,u) - f_applied
// given udot and f_applied as arguments, with the rest from the state.
// No constraint forces are included unless the caller has included them in 
// f_applied.
//
// This Subsystem must already have been realized to Velocity stage in the
// supplied state. All vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcTreeResidualForces(const State& s,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    residualMobilityForces) const
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const SBTreeVelocityCache& tvc = getTreeVelocityCache(s);

    // We allow the input Vectors to be zero length, meaning they are to be
    // considered as though they were full length but all zero. For now we have
    // to make explicit Vectors of zero for them in that case; better would
    // be to have a special operator.
    const Vector*              pAppliedMobForces  = &appliedMobilityForces;
    const Vector_<SpatialVec>* pAppliedBodyForces = &appliedBodyForces;
    const Vector*              pKnownUdot         = &knownUdot;

    Vector              zeroPerMobility;
    Vector_<SpatialVec> zeroPerBody;
    if (appliedMobilityForces.size()==0 || knownUdot.size()==0) {
        zeroPerMobility.resize(getNumMobilities());
        zeroPerMobility = 0;
        if (appliedMobilityForces.size()==0) 
            pAppliedMobForces = &zeroPerMobility;
        if (knownUdot.size()==0)             
            pKnownUdot        = &zeroPerMobility;
    }
    if (appliedBodyForces.size()==0) {
        zeroPerBody.resize(getNumBodies());
        zeroPerBody = SpatialVec(Vec3(0),Vec3(0));
        pAppliedBodyForces = &zeroPerBody;
    }

    // At this point the three pointers point either to the original arguments
    // or to appropriate-sized arrays of zero. Any non-zero length original 
    // arguments should already have been verified by the caller (a method
    // in the SimTK API) to be the correct length.

    // Check input sizes.
    assert(pAppliedMobForces->size()  == getNumMobilities());
    assert(pAppliedBodyForces->size() == getNumBodies());
    assert(pKnownUdot->size()         == getNumMobilities());

    // Resize outputs if necessary.
    A_GB.resize(getNumBodies());
    residualMobilityForces.resize(getNumMobilities());

    assert(pAppliedMobForces->hasContiguousData());
    assert(pAppliedBodyForces->hasContiguousData());
    assert(pKnownUdot->hasContiguousData());
    assert(A_GB.hasContiguousData());
    assert(residualMobilityForces.hasContiguousData());


    // Allocate temporary.
    Vector_<SpatialVec> allFTmp(getNumBodies());

    // Make pointers to (contiguous) Vector data for fast access.
    const Real* knownUdotPtr = &(*pKnownUdot)[0];
    SpatialVec* aPtr = A_GB.size() ? &A_GB[0] : NULL;
    const Real* mobilityForcePtr = pAppliedMobForces->size() 
                                   ? &(*pAppliedMobForces)[0] : NULL;
    const SpatialVec* bodyForcePtr = pAppliedBodyForces->size() 
                                     ? &(*pAppliedBodyForces)[0] : NULL;
    Real *residualPtr = residualMobilityForces.size() 
                        ? &residualMobilityForces[0] : NULL;
    SpatialVec* tempPtr = allFTmp.size() ? &allFTmp[0] : NULL;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcBodyAccelerationsFromUdotOutward
               (tpc,tvc,knownUdotPtr,aPtr);
        }

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInverseDynamicsPass2Inward(
                tpc,tvc,aPtr,
                mobilityForcePtr,bodyForcePtr,
                tempPtr,residualPtr);
        }
}
//........................ CALC TREE RESIDUAL FORCES ...........................



//==============================================================================
//                               MULTIPLY BY N
//==============================================================================
// q=Nu or u=~Nq
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyByN
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Position
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next());

    assert(in.size() == (transpose?getTotalQAlloc():getTotalDOF()));
    out.resize(transpose?getTotalDOF():getTotalQAlloc());

    assert(in.hasContiguousData());
    assert(out.hasContiguousData());

    const Real*     inp  = in.size()  ? &in[0]  : 0;
    Real*           outp = out.size() ? &out[0] : 0;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& rbn = *rbNodeLevels[i][j];
            const int maxNQ = rbn.getMaxNQ();

            // Skip weld joints: no q's, no work to do here.
            if (maxNQ == 0)
                continue;

            // Find the right piece of the vectors to work with.
            const int qx = rbn.getQIndex();
            const int ux = rbn.getUIndex();
            const int inpx  = transpose ? qx : ux;
            const int outpx = transpose ? ux : qx;

            // TODO: kludge: for now q-like output may have an unused element 
            // because we always allocate the max space. Set the last element 
            // to zero in case it doesn't get written.
            if (!transpose) outp[outpx + maxNQ-1] = 0;

            rbn.multiplyByN(sbState, transpose, &inp[inpx], &outp[outpx]);
        }
}



//==============================================================================
//                              MULTIPLY BY NDOT
//==============================================================================
// q=NDot*u or u=~NDot*q
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyByNDot
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Velocity
    const SBStateDigest sbState(s, *this, Stage(Stage::Velocity).next());

    assert(in.size() == (transpose?getTotalQAlloc():getTotalDOF()));
    out.resize(transpose?getTotalDOF():getTotalQAlloc());

    assert(in.hasContiguousData());
    assert(out.hasContiguousData());

    const Real*     inp  = in.size()  ? &in[0]  : 0;
    Real*           outp = out.size() ? &out[0] : 0;

    // Skip ground; it doesn't have q's or u's!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& rbn = *rbNodeLevels[i][j];
            const int maxNQ = rbn.getMaxNQ();

            // Skip weld joints: no q's, no work to do here.
            if (maxNQ == 0)
                continue;

            // Find the right piece of the vectors to work with.
            const int qx = rbn.getQIndex();
            const int ux = rbn.getUIndex();
            const int inpx  = transpose ? qx : ux;
            const int outpx = transpose ? ux : qx;

            // TODO: kludge: for now q-like output may have an unused element 
            // because we always allocate the max space. Set the last element 
            // to zero in case it doesn't get written.
            if (!transpose) outp[outpx + maxNQ-1] = 0;

            rbn.multiplyByNDot(sbState, transpose, &inp[inpx], &outp[outpx]);
        }
}



//==============================================================================
//                              MULTIPLY BY NINV
//==============================================================================
// u= NInv * q or q = ~NInv * u
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyByNInv
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Position
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next());

    assert(in.size() == (transpose?getTotalDOF():getTotalQAlloc()));
    out.resize(transpose?getTotalQAlloc():getTotalDOF());

    assert(in.hasContiguousData());
    assert(out.hasContiguousData());

    const Real*     inp  = in.size()  ? &in[0]  : 0;
    Real*           outp = out.size() ? &out[0] : 0;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& rbn = *rbNodeLevels[i][j];
            const int maxNQ = rbn.getMaxNQ();

            // Skip weld joints: no q's, no work to do here.
            if (maxNQ == 0)
                continue;

            // Find the right piece of the vectors to work with.
            const int qx = rbn.getQIndex();
            const int ux = rbn.getUIndex();
            const int inpx  = transpose ? ux : qx;
            const int outpx = transpose ? qx : ux;

            // TODO: kludge: for now q-like output may have an unused element 
            // because we always allocate the max space. Set the last element 
            // to zero in case it doesn't get written.
            if (transpose) outp[outpx + maxNQ-1] = 0;

            rbn.multiplyByNInv(sbState, transpose, &inp[inpx], &outp[outpx]);
        }
}



// =============================================================================
//                        CALC MOBILIZER REACTION FORCES
// =============================================================================
// This method calculates mobilizer reaction forces using repeated application
// of the equation
//     F_reaction = PPlus*APlus + zPlus
// where P is an articulated body inertia, A is the spatial acceleration of
// that body, and z is the articulated body residual force. The "Plus" 
// indicates that the quantity is as seen on the *inboard* side of the
// mobilizer, although still measured about Bo and expressed in G.
// All of these quantities are already available at Stage::Acceleration except
// APlus= ~phi * A_GP, the parent body's acceleration shifted to the child.
//
// See Abhi Jain's 2011 book "Robot and Multibody Dynamics", Eq. 7.34 on
// page 128: reaction = P(A-a)+z = PPlus*APlus + zPlus. The first equation
// is not correct if there is prescribed motion (you'd have to remove H*udot
// also), but the "Plus" version works regardless.
//
// After calculating the reaction at the body frame origin Bo, we shift it to
// the mobilizer's outboard frame M and report it there, though expressed in G.
// Note that any generalized forces applied at mobilities end up included in
// the reaction forces.
//
// Cost is 114 flops/body plus lots of memory access to dredge up the 
// already-calculated goodies. If you don't need all the reactions, you can 
// calculate them one at a time as needed just as efficiently.
void SimbodyMatterSubsystemRep::calcMobilizerReactionForces
   (const State& s, Vector_<SpatialVec>& FM_G) const 
{
    const int nb = getNumBodies();
    // We're going to work with forces in Ground, applied at the body frame
    // of each body. Then at the end we'll shift to the M frame as promised
    // (though still expressed in Ground).
    FM_G.resize(nb);

    const Array_<ArticulatedInertia,MobilizedBodyIndex>& PPlus = 
                                            getArticulatedBodyInertiasPlus(s);
    const Array_<SpatialVec,MobilizedBodyIndex>& zPlus =
                                            getArticulatedBodyForcesPlus(s);

    for (MobodIndex mbx(0); mbx < nb; ++mbx) {
        const MobilizedBody& body = getMobilizedBody(mbx);
        const Transform& X_GB = body.getBodyTransform(s);
     
        SpatialVec FB_G = zPlus[mbx];
        if (mbx != GroundIndex) {
            const MobilizedBody& parent = body.getParentMobilizedBody();
            const Transform&  X_GP = parent.getBodyTransform(s);
            const SpatialVec& A_GP = parent.getBodyAcceleration(s);
            const Vec3& p_PB_G = X_GB.p() - X_GP.p(); // 3 flops
            SpatialVec APlus( A_GP[0],
                              A_GP[1] + A_GP[0] % p_PB_G ); // 12 flops
            FB_G += PPlus[mbx]*APlus; // 72 flops
        }
        // Shift to M
        const Vec3&      p_BM   = body.getOutboardFrame(s).p();
        const Vec3       p_BM_G = X_GB.R()*p_BM; // p_BM in G, 15 flops
        FM_G[mbx] = shiftForceBy(FB_G, p_BM_G);  // 12 flops
    }
}
//....................... CALC MOBILIZER REACTION FORCES .......................



// =============================================================================
//            CALC MOBILIZER REACTION FORCES USING FREEBODY METHOD
// =============================================================================
// This method provides an alternative way to calculate mobilizer reaction
// forces. It is about 3X slower than calcMobilizerReactionForces() so should
// not be used except for Simbody debugging and regression testing purposes.
void SimbodyMatterSubsystemRep::calcMobilizerReactionForcesUsingFreebodyMethod
   (const State& s, Vector_<SpatialVec>& FM_G) const 
{
    const int nb = getNumBodies();
    // We're going to work with forces in Ground, applied at the body frame
    // of each body. Then at the end we'll shift to the M frame as promised
    // (though still expressed in Ground).
    FM_G.resize(nb);
    
    // Find the body forces on every body from all sources *other* than 
    // mobilizer reaction forces; we accumulate them in otherForces_G.
    
    // First, get the applied body forces (at Bo).
    Vector_<SpatialVec> otherFB_G = 
        getMultibodySystem().getRigidBodyForces(s, Stage::Dynamics);

    // Plus body forces applied by constraints (watch the sign).
    Vector_<SpatialVec> constrainedBodyForces_G(getNumBodies());
    Vector constrainedMobilizerForces(s.getNU());
    calcConstraintForcesFromMultipliers(s, s.getMultipliers(), 
        constrainedBodyForces_G, constrainedMobilizerForces);
    otherFB_G -= constrainedBodyForces_G;

    // We'll account below for gyroscopic forces due to angular velocity.

    // Starting from the leaf nodes and working back toward ground, take the 
    // difference between the total force and the other forces we can account
    // for to find the reaction force, then apply that to the parent as a
    // known force.  
    for (int i = (int)rbNodeLevels.size()-1; i >= 0; --i)
        for (int j = 0; j < (int)rbNodeLevels[i].size(); ++j) {
            const MobilizedBodyIndex mbx = rbNodeLevels[i][j]->getNodeNum();
            const MobilizedBody& body   = getMobilizedBody(mbx);

            SpatialVec totalFB_G(Vec3(0)); // Force from freebody f=ma
            if (mbx != GroundIndex) {
                const SpatialVec&     A_GB = body.getBodyAcceleration(s);
                const SpatialInertia& MB_G = body.getBodySpatialInertiaInGround(s);
                totalFB_G = MB_G * A_GB; // f = ma (45 flops)
            }

            // Body B's reaction force, but applied at Bo
            const SpatialVec FB_G = totalFB_G 
                                    - (otherFB_G[mbx]-getGyroscopicForce(s, mbx));
            // Shift to M
            const Transform& X_GB   = body.getBodyTransform(s);
            const Vec3&      p_BM   = body.getOutboardFrame(s).p();
            const Vec3       p_BM_G = X_GB.R()*p_BM; // p_BM, but expressed in G
            FM_G[mbx] = shiftForceBy(FB_G, p_BM_G);
            if (i==0) continue; // no parent

            // Now apply reaction to parent as an "otherForce". Apply equal and 
            // opposite force & torque to parent at Po, by shifting the forces
            // at Bo.
            const MobilizedBody& parent = body.getParentMobilizedBody();
            const MobilizedBodyIndex px = parent.getMobilizedBodyIndex();
            const Transform& X_GP = parent.getBodyTransform(s);
            const Vec3 p_BP_G = X_GP.p() - X_GB.p();
            otherFB_G[px] -= shiftForceBy(FB_G, p_BP_G);
        }
}
//....................... CALC MOBILIZER REACTION FORCES .......................



//==============================================================================
//                             CALC MOTION ERRORS
//==============================================================================
// Return one error per known value, in order of mobilized bodies.
Vector SimbodyMatterSubsystemRep::
calcMotionErrors(const State& s, const Stage& stage) const
{
    const SBInstanceCache& ic = getInstanceCache(s);
    Vector errs;

    if (stage == Stage::Position) {
        const SBTimeCache& tc = getTimeCache(s);
        assert(tc.presQPool.size() == ic.getTotalNumPresQ());
        errs.resize(ic.getTotalNumPresQ()+ic.getTotalNumZeroQ());
        int nxt = 0;
        for (MobodIndex mbx(1); mbx < getNumBodies(); ++mbx) {
            const Mobod& mobod = getMobilizedBody(mbx);
            const int nq = mobod.getNumQ(s);
            const MobilizedBodyImpl& mbimpl = mobod.getImpl();
            const SBInstancePerMobodInfo& mbinfo = mbimpl.getMyInstanceInfo(s);
            if (mbinfo.qMethod == Motion::Zero) {
                errs(nxt, nq) = mobod.getQAsVector(s); //TODO not right for quats
                nxt += nq;
            } else if (mbinfo.qMethod == Motion::Prescribed) {
                Vector pres(nq, &tc.presQPool[mbinfo.firstPresQ], true);
                errs(nxt, nq) = mobod.getQAsVector(s) - pres;
                nxt += nq;
            }
        }
        return errs;
    }

    if (stage == Stage::Velocity) {
        const SBConstrainedPositionCache& cpc = getConstrainedPositionCache(s);
        assert(cpc.presUPool.size() == ic.getTotalNumPresU());
        errs.resize(ic.getTotalNumPresU()+ic.getTotalNumZeroU());
        int nxt = 0;
        for (MobodIndex mbx(1); mbx < getNumBodies(); ++mbx) {
            const Mobod& mobod = getMobilizedBody(mbx);
            const int nu = mobod.getNumU(s);
            const MobilizedBodyImpl& mbimpl = mobod.getImpl();
            const SBInstancePerMobodInfo& mbinfo = mbimpl.getMyInstanceInfo(s);
            if (mbinfo.uMethod == Motion::Zero) {
                errs(nxt, nu) = mobod.getUAsVector(s);
                nxt += nu;
            } else if (mbinfo.uMethod == Motion::Prescribed) {
                Vector pres(nu, &cpc.presUPool[mbinfo.firstPresU], true);
                errs(nxt, nu) = mobod.getUAsVector(s) - pres;
                nxt += nu;
            }
        }
        return errs;
    }

    if (stage == Stage::Acceleration) {
        const SBDynamicsCache& dc = getDynamicsCache(s);
        assert(dc.presUDotPool.size() == ic.getTotalNumPresUDot());
        errs.resize(ic.getTotalNumPresUDot()+ic.getTotalNumZeroUDot());
        int nxt = 0;
        for (MobodIndex mbx(1); mbx < getNumBodies(); ++mbx) {
            const Mobod& mobod = getMobilizedBody(mbx);
            const int nu = mobod.getNumU(s);
            const MobilizedBodyImpl& mbimpl = mobod.getImpl();
            const SBInstancePerMobodInfo& mbinfo = mbimpl.getMyInstanceInfo(s);
            if (mbinfo.udotMethod == Motion::Zero) {
                errs(nxt, nu) = mobod.getUDotAsVector(s);
                nxt += nu;
            } else if (mbinfo.udotMethod == Motion::Prescribed) {
                Vector pres(nu, &dc.presUDotPool[mbinfo.firstPresUDot], true);
                errs(nxt, nu) = mobod.getUDotAsVector(s) - pres;
                nxt += nu;
            }
        }
        return errs;
    }

    SimTK_APIARGCHECK1_ALWAYS
       (Stage::Position <= stage && stage <= Stage::Acceleration,
        "SimbodyMatterSubsystem", "calcMotionErrors()",
        "Requested stage must be Position, Velocity, or Acceleration but"
        " was %s.", stage.getName().c_str());

    /*NOTREACHED*/
    return errs;
}



//==============================================================================
//                    FIND MOTION FORCES, CALC MOTION POWER
//==============================================================================
void SimbodyMatterSubsystemRep::
findMotionForces(const State&   s,
                 Vector&        mobilityForces) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int nu = getNU(s);
    const int nu_p = ic.getTotalNumPresForces(); // num prescribed udots

    mobilityForces.resize(nu); 
    mobilityForces.setToZero(); // assume nothing is prescribed
    if (nu_p == 0)
        return; // nothing more to do

    const Vector& tau = getMotionMultipliers(s); // nu_p of these, contiguous
    assert(tau.size() == nu_p);

    const Real* taup = &tau[0];
    for (PresForcePoolIndex i(0); i < nu_p; ++i)
        mobilityForces[ic.presForce[i]] = taup[i];
}

Real SimbodyMatterSubsystemRep::
calcMotionPower(const State& s) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int nu_p = ic.getTotalNumPresForces(); // num prescribed udots
    if (nu_p == 0)
        return 0;

    const Vector& tau = getMotionMultipliers(s); // nu_p of these, contiguous
    assert(tau.size() == nu_p);
    assert(tau.hasContiguousData());

    const Vector& u = getU(s);
    assert(u.hasContiguousData());

    const Real* taup = &tau[0];
    const Real* up   = &u[0];
    Real power = 0;

    // Note that we're negating tau to get +power to mean adding energy
    // and -power to mean removing energy. That's required because tau's
    // appear on the LHS of the equations of motion so have the opposite
    // sign from applied forces.
    for (PresForcePoolIndex i(0); i < nu_p; ++i)
        power -= taup[i] * up[ic.presForce[i]]; 

    return power;
}
//................... FIND MOTION FORCES, CALC MOTION POWER ....................



//==============================================================================
//               FIND CONSTRAINT FORCES, CALC CONSTRAINT POWER
//==============================================================================
void SimbodyMatterSubsystemRep::
findConstraintForces(const State&           s, 
                     Vector_<SpatialVec>&   bodyForcesInG,
                     Vector&                mobilityForces) const 
{
    const SBInstanceCache& ic = getInstanceCache(s);
    const int nb = getNumBodies();
    const int nu = getNU(s);

    // Set all forces to zero here; we'll only update the ones that are
    // affected by some active constraint.
    bodyForcesInG.resize(getNumBodies()); bodyForcesInG.setToZero();
    mobilityForces.resize(getNU(s));      mobilityForces.setToZero();

    // Loop over all enabled constraints, get their forces, and
    // accumulate the results in the global problem return vectors.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const ConstraintImpl& crep = constraints[cx]->getImpl();
        const SBInstancePerConstraintInfo& 
                              cInfo = ic.getConstraintInstanceInfo(cx);

        // No heap allocation is being done here. These are views directly
        // into the proper segment of the longer array.
        const ArrayViewConst_<SpatialVec,ConstrainedBodyIndex> bodyF1_G = 
            crep.getConstrainedBodyForcesInGFromState(s);
        const ArrayViewConst_<Real,ConstrainedUIndex>          mobilityF1 = 
            crep.getConstrainedMobilityForcesFromState(s);

        const int ncb = bodyF1_G.size();
        const int ncu = mobilityF1.size();

        // Unpack constrained body forces and add them to the proper slots 
        // in the global body forces array. They are already expressed in
        // the Ground frame.
        for (ConstrainedBodyIndex cbx(0); cbx < ncb; ++cbx)
            bodyForcesInG[crep.getMobilizedBodyIndexOfConstrainedBody(cbx)] 
                += bodyF1_G[cbx];       // 6 flops

        // Unpack constrained mobility forces and add them into global array.
        for (ConstrainedUIndex cux(0); cux < ncu; ++cux) 
            mobilityForces[cInfo.getUIndexFromConstrainedU(cux)] 
                += mobilityF1[cux];     // 1 flop
    }
}

Real SimbodyMatterSubsystemRep::
calcConstraintPower(const State& s) const {
    Real power = 0;

    // Loop over all enabled constraints, get their forces, and
    // accumulate the results in the global problem return vectors.
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const Constraint& constraint = getConstraint(cx);
        power += constraint.calcPower(s); // sign already correct
    }

    return power;
}
//............... FIND CONSTRAINT FORCES, CALC CONSTRAINT POWER ................



//==============================================================================
//                                CALC QDOT
//==============================================================================
// Must be done with Position stage to calculate qdot = N*u.
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcQDot
   (const State& s, const Vector& u, Vector& qdot) const 
{
    SBStateDigest sbs(s, *this, Stage::Velocity);

    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    assert(u.hasContiguousData());
    assert(qdot.hasContiguousData());

    const Real* uPtr    = u.size() ? &u[0] : NULL;
    Real*       qdotPtr = qdot.size() ? &qdot[0] : NULL;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcQDot(sbs, &uPtr[node.getUIndex()], &qdotPtr[node.getQIndex()]);
        }
}
//............................... CALC QDOT ....................................



// =============================================================================
//                              CALC QDOTDOT
// =============================================================================
// Must be done with Velocity stage to calculate qdotdot = Ndot*u + N*udot.
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcQDotDot
   (const State& s, const Vector& udot, Vector& qdotdot) const 
{
    SBStateDigest sbs(s, *this, Stage::Dynamics);

    assert(udot.size() == getTotalDOF());
    qdotdot.resize(getTotalQAlloc());

    assert(udot.hasContiguousData());
    assert(qdotdot.hasContiguousData());

    const Real* udotPtr    = udot.size() ? &udot[0] : NULL;
    Real*       qdotdotPtr = qdotdot.size() ? &qdotdot[0] : NULL;

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcQDotDot(sbs, &udotPtr[node.getUIndex()], 
                             &qdotdotPtr[node.getQIndex()]);
        }
}
//............................. CALC QDOTDOT ...................................




// State must be in Stage::Position.
void SimbodyMatterSubsystemRep::
calcMobilizerQDotFromU(const State& s, MobilizedBodyIndex mb, int nu, const Real* u, 
                       int nq, Real* qdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Position);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQInUse(sbState.getModelVars()));

    n.calcQDot(sbState, u, qdot);
}

// State must be realized to Stage::Velocity, so that we can extract N(q),
// NDot(q,u), and u from it to calculate qdotdot=N(q)*udot + NDot(q,u)*u for 
// this mobilizer.
void SimbodyMatterSubsystemRep::
calcMobilizerQDotDotFromUDot(const State& s, MobilizedBodyIndex mb, 
                             int nu, const Real* udot, 
                             int nq, Real* qdotdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Velocity);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQInUse(sbState.getModelVars()));

    n.calcQDotDot(sbState, udot, qdotdot);
}

// State must be realized through Stage::Instance. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of q's is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized coordinates.
// Returns X_FM(q).
Transform SimbodyMatterSubsystemRep::
calcMobilizerTransformFromQ(const State& s, MobilizedBodyIndex mb, 
                            int nq, const Real* q) const {
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
calcMobilizerVelocityFromU(const State& s, MobilizedBodyIndex mb, 
                           int nu, const Real* u) const {
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
calcMobilizerAccelerationFromUDot(const State& s, MobilizedBodyIndex mb, 
                                  int nu, const Real* udot) const{
    const SBStateDigest sbState(s, *this, Stage::Velocity);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());

    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}

// These perform the same computations as above but then transform the results 
// so that they relate the child body's frame B to its parent body's frame P, 
// rather than the M and F frames which are attached to B and P respectively 
// but differ by a constant transform.
Transform SimbodyMatterSubsystemRep::
calcParentToChildTransformFromQ(const State& s, MobilizedBodyIndex mb, 
                                int nq, const Real* q) const {
    const Transform& X_PF = getMobilizerFrameOnParent(s,mb);
    const Transform& X_BM = getMobilizerFrame(s,mb);

    const Transform X_FM = calcMobilizerTransformFromQ(s,mb,nq,q);
    return X_PF * X_FM * ~X_BM; // X_PB
}

SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildVelocityFromU(const State& s, MobilizedBodyIndex mb, 
                               int nu, const Real* u) const {
    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}
SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildAccelerationFromUDot(const State& s, MobilizedBodyIndex mb, 
                                      int nu, const Real* udot) const {
    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}




// =============================================================================
//                         MULTIPLY BY SYSTEM JACOBIAN
// =============================================================================
// We have V_GB = J u where J=~Phi*~H is the kinematic Jacobian (partial 
// velocity matrix) that maps generalized speeds to spatial velocities. 
// This method performs the multiplication J*u in O(n) time (i.e., without
// actually forming J).
// The input and output vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyBySystemJacobian(const State& s,
    const Vector&              v,
    Vector_<SpatialVec>&       Jv) const 
{
    Jv.resize(getNumBodies());

    assert(v.size() == getNU(s));
    assert(v.hasContiguousData() && Jv.hasContiguousData());

    const SBTreePositionCache& tpc = getTreePositionCache(s);

    const Real* vPtr = v.size() ? &v[0] : NULL;
    SpatialVec* jvPtr = Jv.size() ? &Jv[0] : NULL;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyBySystemJacobian(tpc, vPtr, jvPtr);
        }
}
//......................... MULTIPLY BY SYSTEM JACOBIAN ........................



// =============================================================================
//                     MULTIPLY BY SYSTEM JACOBIAN TRANSPOSE
// =============================================================================
// The system Jacobian J (a.k.a. partial velocity matrix) is the kinematic 
// mapping between generalized speeds u and body spatial velocities V. Its
// transpose ~J maps body spatial forces to generalized forces. This method
// calculates in O(n) time the product of ~J and a "spatial force-like" 
// vector X.
// The input and output vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::multiplyBySystemJacobianTranspose
   (const State&                s, 
    const Vector_<SpatialVec>&  X,
    Vector&                     JtX) const
{
    assert(X.size() == getNumBodies());
    JtX.resize(getNU(s));

    assert(X.hasContiguousData() && JtX.hasContiguousData());

    const SBTreePositionCache& tpc = getTreePositionCache(s);

    Vector_<SpatialVec> zTemp(getNumBodies()); zTemp.setToZero();
    const SpatialVec* xPtr = X.size() ? &X[0] : NULL;
    Real* jtxPtr = JtX.size() ? &JtX[0] : NULL;
    SpatialVec* zPtr = zTemp.size() ? &zTemp[0] : NULL;

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.multiplyBySystemJacobianTranspose(tpc, zPtr, xPtr, jtxPtr);
        }
}
//................... MULTIPLY BY SYSTEM JACOBIAN TRANSPOSE ....................



// =============================================================================
//                     CALC TREE EQUIVALENT MOBILITY FORCES
// =============================================================================
// This routine does the same thing as the above but accounts for centrifugal
// forces induced by velocities. The equivalent joint forces returned include
// both the applied forces and the centrifugal ones. Constraints are ignored.
// Both vectors must use contiguous storage.
// TODO: is this useful for anything?
void SimbodyMatterSubsystemRep::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const
{
    const SBTreePositionCache&  tpc = getTreePositionCache(s);
    const SBTreeVelocityCache&  tvc = getTreeVelocityCache(s);

    assert(bodyForces.size() == getNumBodies());
    mobilityForces.resize(getTotalDOF());

    assert(bodyForces.hasContiguousData());
    assert(mobilityForces.hasContiguousData());

    Vector_<SpatialVec> allZ(getNumBodies());
    const SpatialVec* bodyForcePtr = bodyForces.size() ? &bodyForces[0] : NULL;
    Real* mobilityForcePtr = mobilityForces.size() ? &mobilityForces[0] : NULL;
    SpatialVec* zPtr = allZ.size() ? &allZ[0] : NULL;

    // Don't do ground's level since ground has no inboard joint.
    for (int i=rbNodeLevels.size()-1 ; i>0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcEquivalentJointForces(tpc,tvc,
                bodyForcePtr, zPtr,
                mobilityForcePtr);
        }
}
//.................... CALC TREE EQUIVALENT MOBILITY FORCES ....................



bool SimbodyMatterSubsystemRep::getShowDefaultGeometry() const 
{   return m_showDefaultGeometry; }
void SimbodyMatterSubsystemRep::setShowDefaultGeometry(bool show) 
{   m_showDefaultGeometry = show; }


bool SimbodyMatterSubsystemRep::getUseEulerAnglesByDefault() const 
{   return m_useEulerAnglesByDefault; }
void SimbodyMatterSubsystemRep::setUseEulerAnglesByDefault(bool useAngles) 
{   m_useEulerAnglesByDefault = useAngles; }

std::ostream& SimTK::operator<<(std::ostream& o, 
                                const SimbodyMatterSubsystemRep& tree) {
    o << "SimbodyMatterSubsystemRep has " << tree.getNumBodies() << " bodies (incl. G) in "
      << tree.rbNodeLevels.size() << " levels." << std::endl;
    o << "NodeNum->level,offset;stored nodeNum,level (stateOffset:dim)" << std::endl;
    for (MobilizedBodyIndex i(0); i < tree.getNumBodies(); ++i) {
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

    // This is *not* from the State; the copy in State is not used.
    const SBTopologyCache& topo = matter.getMatterTopologyCache();

    if (state.getSystemStage() >= Stage::Model) {
        q = &matter.getQ(state);
        u = &matter.getU(state);
    }
    if (g >= Stage::Model) {
        mv = &matter.getModelVars(state);
        mc = &matter.updModelCache(state);
        iv = &matter.getInstanceVars(state);
    }
    if (g >= Stage::Instance) {
        if (topo.instanceCacheIndex.isValid())
            ic = &matter.updInstanceCache(state);
        
        // All cache entries, for any stage, can be modified at instance stage 
        // or later.
        
        if (topo.timeCacheIndex.isValid())
            tc = &matter.updTimeCache(state);
        if (topo.treePositionCacheIndex.isValid())
            tpc = &matter.updTreePositionCache(state);
        if (topo.constrainedPositionCacheIndex.isValid())
            cpc = &matter.updConstrainedPositionCache(state);
        if (topo.treeVelocityCacheIndex.isValid())
            tvc = &matter.updTreeVelocityCache(state);
        if (topo.constrainedVelocityCacheIndex.isValid())
            cvc = &matter.updConstrainedVelocityCache(state);
        if (topo.dynamicsCacheIndex.isValid())
            dc = &matter.updDynamicsCache(state);
        if (topo.treeAccelerationCacheIndex.isValid())
            tac = &matter.updTreeAccelerationCache(state);
        if (topo.constrainedAccelerationCacheIndex.isValid())
            cac = &matter.updConstrainedAccelerationCache(state);

        qdot = &matter.updQDot(state);
    }
    if (g >= Stage::Time) {
        if (mc->timeVarsIndex.isValid())
            tv = &matter.getTimeVars(state);

        // These are unavailable until Instance stage is done.
        qErr = &matter.updQErr(state);
        uErr = &matter.updUErr(state);
    }
    if (g >= Stage::Position) {
        if (mc->qVarsIndex.isValid())
            pv = &matter.getPositionVars(state);
    }
    if (g >= Stage::Velocity) {
        if (mc->uVarsIndex.isValid())
            vv = &matter.getVelocityVars(state);
    }
    if (g >= Stage::Dynamics) {
        if (mc->dynamicsVarsIndex.isValid())
            dv = &matter.getDynamicsVars(state);

        // Prescribed accelerations are filled in at Dynamics
        // stage so we may need these now.
        udot = &matter.updUDot(state);
        qdotdot = &matter.updQDotDot(state);
    }
    if (g >= Stage::Acceleration) {
        if (mc->accelerationVarsIndex.isValid())
            av = &matter.getAccelerationVars(state);

        udotErr = &matter.updUDotErr(state);
    }

    stage = g;
}

