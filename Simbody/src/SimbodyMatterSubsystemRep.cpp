/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Portions derived from NIH IVM code written by                *
 *               Charles Schwieters                                           *
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
 */

#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "simbody/internal/common.h"

#include "SimbodyMatterSubsystemRep.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "MultibodySystemRep.h"
#include "MobilizedBodyImpl.h"
#include "ConstraintImpl.h"

#include <string>
#include <iostream>
using std::cout; using std::endl;

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

    showDefaultGeometry = true;
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
                                                       MobilizedBodyIndex(),    // no parent 
                                                       MobilizedBodyIndex(0));
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
SimbodyMatterSubsystemRep::getCoriolisAcceleration(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getCoriolisAcceleration(getTreeVelocityCache(s));
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
SimbodyMatterSubsystemRep::getCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getCentrifugalForces(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getTotalCentrifugalForces(getDynamicsCache(s));
}



//------------------------------------------------------------------------------
//                               REALIZE TOPOLOGY
//------------------------------------------------------------------------------

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
        nodeNum2NodeMap.push_back(RigidBodyNodeIndex(level, nodeIndexWithinLevel));

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
}

int SimbodyMatterSubsystemRep::realizeSubsystemTopologyImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "SimbodyMatterSubsystem::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this 
    // call*. Afterwards they are truly const so we don't declare them mutable,
    // but cheat here instead.
    SimbodyMatterSubsystemRep* mutableThis = const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!subsystemTopologyHasBeenRealized()) mutableThis->endConstruction(s); // no more bodies after this!

    // Fill in the local copy of the topologyCache from the information
    // calculated in endConstruction(). Also ask the State for some room to
    // put Modeling variables & cache and remember the indices in our 
    // construction cache.

    mutableThis->topologyCache.nBodies      = nodeNum2NodeMap.size();
    mutableThis->topologyCache.nConstraints = constraints.size();
    mutableThis->topologyCache.nParticles   = 0; // TODO
    mutableThis->topologyCache.nAncestorConstrainedBodies = 
        (int)nextAncestorConstrainedBodyPoolSlot;

    mutableThis->topologyCache.nDOFs        = DOFTotal;
    mutableThis->topologyCache.maxNQs       = maxNQTotal;
    mutableThis->topologyCache.sumSqDOFs    = SqDOFTotal;

    SBModelVars mvars;
    mvars.allocate(topologyCache);
    setDefaultModelValues(topologyCache, mvars);
    mutableThis->topologyCache.modelingVarsIndex  = 
        allocateDiscreteVariable(s,Stage::Model, new Value<SBModelVars>(mvars));

    mutableThis->topologyCache.modelingCacheIndex = 
        allocateCacheEntry(s,Stage::Model, new Value<SBModelCache>());

    SBInstanceVars iv;
    iv.allocate(topologyCache);
    setDefaultInstanceValues(mvars, iv);
    mutableThis->topologyCache.topoInstanceVarsIndex = 
        allocateDiscreteVariable(s, Stage::Instance, new Value<SBInstanceVars>(iv));

    mutableThis->topologyCache.valid = true;

    // Allocate a cache entry for the topologyCache, and save a copy there.
    mutableThis->topologyCacheIndex = 
        allocateCacheEntry(s,Stage::Topology, new Value<SBTopologyCache>(topologyCache));

    return 0;
}



//------------------------------------------------------------------------------
//                                 REALIZE MODEL
//------------------------------------------------------------------------------
// Here we lock in modeling choices as conveyed by the values of Model-stage
// state variables which now all have values. These choices determine the number 
// and types of state variables we're going to use to represent the changeable 
// properties of this matter subsystem. This is the last realization stage at
// which we are given a writable State. That means all the state variables 
// we'll ever need must be allocated here (although cache entries can be added
// later since they are mutable.)
//
// Variables we'll settle on here include:
//  - whether to use 4 quaternions or 3 Euler angles to represent orientation
//
// We allocate and fill in the Model-stage cache with information that can be
// calculated now that we have values for the Model-stage state variables.

int SimbodyMatterSubsystemRep::realizeSubsystemModelImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Topology, 
        "SimbodyMatterSubsystem::realizeModel()");

    SBStateDigest sbs(s, *this, Stage::Model);
    const SBModelVars& mv = sbs.getModelVars();

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
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            SBModelCache::PerMobilizedBodyModelInfo& 
                mbInfo = mc.updMobilizedBodyModelInfo(node.getNodeNum());

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

    // SBInstanceVars are allocated at topology stage; could also have instance
    // vars that aren't allocated until here but there aren't any right now.

    mc.instanceCacheIndex = 
        allocateCacheEntry(s, Stage::Instance, new Value<SBInstanceCache>());

    mc.timeVarsIndex = 
        allocateDiscreteVariable(s, Stage::Time, new Value<SBTimeVars>());
    mc.timeCacheIndex = 
        allocateCacheEntry(s, Stage::Time, new Value<SBTimeCache>());

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
    mc.qVarsIndex.invalidate(); // no position-stage vars other than q

    // Basic tree position kinematics can be calculated any time after Time
    // stage and should be filled in first during realizePosition() and then
    // marked valid so later computations during the same realization can
    // access these quantities.
    mc.treePositionCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Time,
                               new Value<SBTreePositionCache>());

    // Here is where later computations during realizePosition() go; these
    // will assume that the TreePositionCache is available. So you can 
    // calculate these prior to Position stage's completion but not until
    // the TreePositionCache has been marked valid.
    mc.constrainedPositionCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Time,
                               new Value<SBConstrainedPositionCache>());

    // Composite body inertias *can* be calculated any time after Position
    // stage but we want to put them off as long as possible since
    // they may never be needed. These will only be valid if they are 
    // explicitly realized at some point.
    mc.compositeBodyInertiaCacheIndex =
        allocateLazyCacheEntry(s, Stage::Position,
                               new Value<SBCompositeBodyInertiaCache>());

    // Articulated body inertias *can* be calculated any time after Position 
    // stage but we want to put them off until Dynamics stage if possible.
    mc.articulatedBodyInertiaCacheIndex =
        allocateCacheEntry(s, Stage::Position, Stage::Dynamics, 
                           new Value<SBArticulatedBodyInertiaCache>());

    // Velocity variables are just the generalized speeds u, which the State 
    // knows how to deal with. Zero is always a reasonable value for velocity,
    // so we'll initialize it here.
    Vector uInit(DOFTotal);
    setDefaultVelocityValues(mv, uInit);

    mc.uIndex = s.allocateU(getMySubsystemIndex(), uInit);
    mc.uVarsIndex.invalidate(); // no velocity-stage vars other than u

    // Basic tree velocity kinematics can be calculated any time after Position
    // stage and should be filled in first during realizeVelocity() and then
    // marked valid so later computations during the same realization can
    // access these quantities. Note that qdots are automatically allocated in 
    // the State's Velocity-stage cache.
    mc.treeVelocityCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Position,
                               new Value<SBTreeVelocityCache>());

    // Here is where later computations during realizeVelocity() go; these
    // will assume that the TreeVelocityCache is available. So you can 
    // calculate these prior to Velocity stage's completion but not until
    // the TreeVelocityCache has been marked valid.
    mc.constrainedVelocityCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Position,
                               new Value<SBConstrainedVelocityCache>());

    // no z's
    // Probably no dynamics-stage variables but we'll allocate anyway.
    SBDynamicsVars dvars;
    dvars.allocate(topologyCache);
    setDefaultDynamicsValues(mv, dvars);

    mc.dynamicsVarsIndex = 
        allocateDiscreteVariable(s, Stage::Dynamics, new Value<SBDynamicsVars>(dvars));
    mc.dynamicsCacheIndex = 
        allocateCacheEntry(s, Stage::Dynamics, new Value<SBDynamicsCache>());

    // No Acceleration variables that I know of. But we can go through the
    // charade here anyway.
    SBAccelerationVars rvars;
    rvars.allocate(topologyCache);
    setDefaultAccelerationValues(mv, rvars);
    mc.accelerationVarsIndex = 
        allocateDiscreteVariable(s, Stage::Acceleration, new Value<SBAccelerationVars>(rvars));

    // Tree acceleration kinematics can be calculated any time after Dynamics
    // stage and should be filled in first during realizeAcceleration() and then
    // marked valid so later computations during the same realization can
    // access these quantities.
    // Note that qdotdots, udots, zdots are automatically allocated by
    // the State when we advance the stage past Model.
    mc.treeAccelerationCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Dynamics,
                               new Value<SBTreeAccelerationCache>());

    // Here is where later computations during realizeAcceleration() go; these
    // will assume that the TreeAccelerationCache is available. So you can 
    // calculate these prior to Acceleration stage's completion but not until
    // the TreeAccelerationCache has been marked valid.
    mc.constrainedAccelerationCacheIndex = 
        allocateLazyCacheEntry(s, Stage::Dynamics,
                               new Value<SBConstrainedAccelerationCache>());

    return 0;
}



//------------------------------------------------------------------------------
//                               REALIZE INSTANCE
//------------------------------------------------------------------------------
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

int SimbodyMatterSubsystemRep::realizeSubsystemInstanceImpl(const State& s) const {
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
    for (int i=0; i<getNumBodies(); ++i)
        ic.totalMass += iv.bodyMassProperties[i].getMass();
    // TODO: central and principal inertias

    // Mobilizer geometry is now available in the InstanceVars.
    // TODO: reference configuration

    // Count position-, velocity-, and acceleration- prescribed motions generated by
    // Motion objects associated with MobilizedBodies and allocate pools to hold the
    // associated values in the state cache. When position is prescribed
    // (by specifying q(t)), the corresponding qdot and qdotdot are also prescribed
    // and we use them to set u and udot (via u=N^-1 qdot and udot = N^-1(qdotdot-NDot*u)).
    // Each prescribed udot will have a corresponding force calculated, and
    // other known udots (zero or discrete; anything but free) will also need force
    // slots although they don't get UDotPool slots.
    //
    // There is no built-in support in the State for these pools, so we allocate them
    // in Simbody's cache entries at the appropriate stages. The prescribed q pool
    // is in the TimeCache, prescribed u (dependent on the TreePositionCache) is 
    // written into the ConstrainedPositionCache, prescribed udots are written 
    // directly into the udot array in the State, and the prescribed forces tau 
    // are calculated at the same time as the unprescribed udots and are thus in 
    // the TreeAccelerationCache.
    //
    // NOTE: despite appearances here, each pool is in MobilizedBodyIndex order, 
    // meaning that the prescribed position, velocity, and acceleration, and the
    // other known udot entries, will be intermingled here rather than neatly 
    // lined up as I've drawn them.
    //
    //            --------------------
    //     QPool |       nPresQ       |      NOTE: not really ordered like this
    //            \------------------/ 
    //             \----------------/-------------
    //     UPool   |              nPresU          |
    //             |----------------|-------------|
    //             |------------------------------|--------------
    //      UDot   |                    nPresUDot                |
    //             |------------------------------|--------------|
    //             |---------------------------------------------|--------------
    // ForcePool   |                           nPresForces                      |
    //              ---------------------------------------------|--------------
    //
    // Note that there are no slots allocated for q's, u's, or udots that are 
    // known to be zero; only the explicitly prescribed ones get a slot. And 
    // udots known for any reason get a slot in the ForcePool.

    // Motion options for all mobilizers are now available in the InstanceVars.
    // We need to figure out what cache resources are required, including
    // slots for the constraint forces that implement prescribed motion.
    ic.presQ.clear();       ic.zeroQ.clear();       ic.freeQ.clear();
    ic.presU.clear();       ic.zeroU.clear();       ic.freeU.clear();
    ic.presUDot.clear();    ic.zeroUDot.clear();    ic.freeUDot.clear();
    ic.totalNPresForce = 0;
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx) {
        const MobilizedBody& mobod = getMobilizedBody(mbx);
        const SBModelCache::PerMobilizedBodyModelInfo& 
            modelInfo    = mc.getMobilizedBodyModelInfo(mbx);
        SBInstanceCache::PerMobodInstanceInfo&
            instanceInfo = ic.updMobodInstanceInfo(mbx);

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

        // Not Ground or a Weld.

        if (mobod.hasMotion()) {
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
        }

        // Assign udots to appropriate index vectors for convenient
        // manipulation. Prescribed udots also need pool allocation,
        // and any non-Free udot needs a prescribed force (tau) slot.
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
        }

        // Count mobilities that need a slot to hold the calculated force due
        // to a known udot, whether prescribed or known for some other reason.
        if (instanceInfo.udotMethod != Motion::Free) {
            instanceInfo.firstPresForce = 
                PresForcePoolIndex(ic.totalNPresForce);
            ic.totalNPresForce += nu;
        }
    }

    // Now instantiate the implementing RigidBodyNodes.
    SBStateDigest stateDigest(s, *this, Stage::Instance);
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeInstance(stateDigest); 


    
    // CONSTRAINT INSTANCE


    // Count position, velocity, and acceleration constraint equations
    // generated by each Constraint that has not been disabled. The State's 
    // QErr, UErr, UDotErr/Multiplier arrays are laid out like this:
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
    // Note that a Constraint with both holonomic and nonholonomic constraint 
    // equations will get two disjoint segments in UErr (and UDotErr).

    // Each Constraint's realizeInstance() method will add its contribution to
    // these.
    ic.totalNHolonomicConstraintEquationsInUse = ic.totalNNonholonomicConstraintEquationsInUse = 
        ic.totalNAccelerationOnlyConstraintEquationsInUse = 0;


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
    updDynamicsCache(s).allocate(topologyCache, mc, ic);
    updTreeAccelerationCache(s).allocate(topologyCache, mc, ic);
    updConstrainedAccelerationCache(s).allocate(topologyCache, mc, ic);
    
    return 0;
}

//------------------------------------------------------------------------------
//                                REALIZE TIME
//------------------------------------------------------------------------------
int SimbodyMatterSubsystemRep::realizeSubsystemTimeImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "SimbodyMatterSubsystem::realizeTime()");

    const SBStateDigest stateDigest(s, *this, Stage::Time);
    const SBModelCache&     mc = stateDigest.getModelCache();
    const SBInstanceCache&  ic = stateDigest.getInstanceCache();
    SBTimeCache&            tc = stateDigest.updTimeCache();

    // the multibody tree cannot have time dependence

    // Now realize the MobilizedBodies, which will realize prescribed positions
    // if there are any; those will go into the TimeCache.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeTime(stateDigest);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeTime(stateDigest);

    // We're done with the TimeCache now.
    markCacheValueRealized(s, mc.timeCacheIndex);
    return 0;
}



//------------------------------------------------------------------------------
//                             REALIZE POSITION
//------------------------------------------------------------------------------
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

    // Set up StateDigest for calculating position information.
    const SBStateDigest stateDigest(s, *this, Stage::Position);
    const SBModelCache&         mc   = stateDigest.getModelCache();
    const SBInstanceCache&      ic   = stateDigest.getInstanceCache();
    SBTreePositionCache&        tpc  = stateDigest.updTreePositionCache();

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
            .calcConstrainedBodyTransformInAncestor(stateDigest, tpc);

    // Now we're done with the TreePositionCache.
    markCacheValueRealized(s, mc.treePositionCacheIndex);

    // MobilizedBodies
    // This will include writing the prescribed velocities (u's) into
    // the PresUPool in the ConstrainedPositionCache.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizePosition(stateDigest);

    // Constraints
    // TODO: should include writing the qErr's
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizePosition(stateDigest);

    // Put position constraint equation errors in qErr
    Vector& qErr = stateDigest.updQErr();
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& pseg = cInfo.holoErrSegment;
        if (pseg.length)
            constraints[cx]->getImpl().realizePositionErrors
               (s, pseg.length, &qErr[pseg.offset]);
    }

    // Now we're done with the ConstrainedPositionCache.
    markCacheValueRealized(s, mc.constrainedPositionCacheIndex);
    return 0;
}



//------------------------------------------------------------------------------
//                      REALIZE COMPOSITE BODY INERTIAS
//------------------------------------------------------------------------------
void SimbodyMatterSubsystemRep::realizeCompositeBodyInertias(const State& state) const {
    const CacheEntryIndex cbx = getModelCache(state).compositeBodyInertiaCacheIndex;

    if (isCacheValueRealized(state, cbx))
        return; // already realized

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position, 
        "SimbodyMatterSubsystem::realizeCompositeBodyInertias()");

    SBCompositeBodyInertiaCache& cbc = 
        Value<SBCompositeBodyInertiaCache>::updDowncast(updCacheEntry(state, cbx));

    calcCompositeBodyInertias(state, cbc.compositeBodyInertia);
    markCacheValueRealized(state, cbx);
}



//------------------------------------------------------------------------------
//                     REALIZE ARTICULATED BODY INERTIAS
//------------------------------------------------------------------------------
void SimbodyMatterSubsystemRep::realizeArticulatedBodyInertias(const State& state) const {
    const CacheEntryIndex abx = getModelCache(state).articulatedBodyInertiaCacheIndex;

    if (isCacheValueRealized(state, abx))
        return; // already realized

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position, 
        "SimbodyMatterSubsystem::realizeArticulatedBodyInertias()");

    const SBInstanceCache&          ic  = getInstanceCache(state);
    const SBTreePositionCache&      tpc = getTreePositionCache(state);
    SBArticulatedBodyInertiaCache&  abc = updArticulatedBodyInertiaCache(state);

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; --i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeArticulatedBodyInertiasInward(ic,tpc,abc);

    markCacheValueRealized(state, abx);
}



//------------------------------------------------------------------------------
//                               REALIZE VELOCITY
//------------------------------------------------------------------------------
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
// In step (2) we calculate the matter subsystem's other position dependencies
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

    const SBStateDigest stateDigest(s, *this, Stage::Velocity);
    const SBModelCache&    mc   = stateDigest.getModelCache();
    const SBInstanceCache& ic   = stateDigest.getInstanceCache();
    SBTreeVelocityCache&   tvc  = stateDigest.updTreeVelocityCache();

    // realize tree velocity kinematics
    // This includes all local cross-mobilizer velocities (M in F, B in P)
    // and all global velocities relative to Ground (G).

    // Set generalized speeds: sweep from base to tips.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeVelocity(stateDigest); 

    // Ask the constraints to calculate ancestor-relative velocity kinematics 
    // (still goes in TreePositionCache).
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl()
            .calcConstrainedBodyVelocityInAncestor(stateDigest, tvc);

    // Now we're done with the TreeVelocityCache.
    markCacheValueRealized(s, mc.treeVelocityCacheIndex);

    // MobilizedBodies
    // probably not much to do here
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeVelocity(stateDigest);

    // Constraints
    // TODO: should include writing the uErr's
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeVelocity(stateDigest);

    // Put velocity constraint equation errors in uErr
    Vector& uErr = stateDigest.updUErr();
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);

        const Segment& holoseg    = cInfo.holoErrSegment; // for derivs of holo constraints
        const Segment& nonholoseg = cInfo.nonholoErrSegment; // includes holo+nonholo
        const int mHolo = holoseg.length, mNonholo = nonholoseg.length;
        if (mHolo)
            constraints[cx]->getImpl().realizePositionDotErrors
               (s, mHolo,    &uErr[holoseg.offset]);
        if (mNonholo)
            constraints[cx]->getImpl().realizeVelocityErrors
               (s, mNonholo, &uErr[ic.totalNHolonomicConstraintEquationsInUse 
                                   + nonholoseg.offset]);
    }

    // Now we're done with the ConstrainedVelocityCache.
    markCacheValueRealized(s, mc.constrainedVelocityCacheIndex);
    return 0;
}


//------------------------------------------------------------------------------
//                               REALIZE DYNAMICS
//------------------------------------------------------------------------------
// Prepare for dynamics by calculating position-dependent quantities like the 
// articulated body inertias P, and velocity-dependent quantities like the 
// Coriolis acceleration.

int SimbodyMatterSubsystemRep::realizeSubsystemDynamicsImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "SimbodyMatterSubsystem::realizeDynamics()");

    // tip-to-base calculation
    realizeArticulatedBodyInertias(s); // (may already have been realized)
    const SBArticulatedBodyInertiaCache& abc = getArticulatedBodyInertiaCache(s);

    SBStateDigest stateDigest(s, *this, Stage::Dynamics);

    // Get the Dynamics-stage cache; it was already allocated at Instance stage.
    SBDynamicsCache& dc = stateDigest.updDynamicsCache();

    // realize velocity-dependent articulated body quantities needed for dynamics
    // base-to-tip
    for (int i=0; i < (int)rbNodeLevels.size(); ++i)
        for (int j=0; j < (int)rbNodeLevels[i].size(); ++j)
            rbNodeLevels[i][j]->realizeDynamics(abc, stateDigest);

    // MobilizedBodies
    // This will include writing the prescribed accelerations into
    // the udot array in the State.
    for (MobilizedBodyIndex mbx(0); mbx < mobilizedBodies.size(); ++mbx)
        getMobilizedBody(mbx).getImpl().realizeDynamics(stateDigest);

    // realize Constraint dynamics

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeDynamics(stateDigest);

    return 0;
}



//------------------------------------------------------------------------------
//                             REALIZE ACCELERATION
//------------------------------------------------------------------------------
int SimbodyMatterSubsystemRep::realizeSubsystemAccelerationImpl(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "SimbodyMatterSubsystem::realizeAcceleration()");

    SBStateDigest stateDigest(s, *this, Stage::Acceleration);

    // Get the Acceleration-stage cache entries. They were all allocated by
    // the end of Instance stage.
    Vector&                         udot    = stateDigest.updUDot();
    Vector&                         qdotdot = stateDigest.updQDotDot();
    SBTreeAccelerationCache&        tac     = stateDigest.updTreeAccelerationCache();
    SBConstrainedAccelerationCache& cac     = stateDigest.updConstrainedAccelerationCache();

    // We ask our containing MultibodySystem for a reference to the cached forces
    // accumulated from all the force subsystems. We use these to compute accelerations,
    // with all results going into the AccelerationCache.
    const MultibodySystem& mbs = getMultibodySystem();  // owner of this subsystem
    realizeLoopForwardDynamics(s,
        mbs.getMobilityForces(s, Stage::Dynamics),
        mbs.getParticleForces(s, Stage::Dynamics),
        mbs.getRigidBodyForces(s, Stage::Dynamics));

    calcQDotDot(s, udot, qdotdot);
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeAcceleration(stateDigest); 

    return 0;
}



//------------------------------------------------------------------------------
//                               REALIZE REPORT
//------------------------------------------------------------------------------
int SimbodyMatterSubsystemRep::realizeSubsystemReportImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "SimbodyMatterSubsystem::realizeReport()");

    // realize MobilizedBody report
    SBStateDigest stateDigest(s, *this, Stage::Report);
    for (int i=0 ; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeReport(stateDigest);

    // realize Constraint report
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeReport(s);

    return 0;
}



//------------------------------------------------------------------------------
//                    CALC DECORATIVE GEOMETRY AND APPEND
//------------------------------------------------------------------------------
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
// DON'T USE THIS
void SimbodyMatterSubsystemRep::calcQUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const {
    bounds[0] = bounds[2] = bounds[4] = Infinity;
    bounds[1] = bounds[3] = bounds[5] = -Infinity;
    
    // Call this method recursively on the body's children, and build up a bounding box
    // for everything downstream of it.
    
    for (int i = 0; i < body.getNumChildren(); ++i) {
        Vec6 childBounds;
        calcQUnitWeightsRecursively(s, tempState, weights, childBounds, *body.getChild(i));
        bounds[0] = std::min(bounds[0], childBounds[0]);
        bounds[2] = std::min(bounds[2], childBounds[2]);
        bounds[4] = std::min(bounds[4], childBounds[4]);
        bounds[1] = std::max(bounds[1], childBounds[1]);
        bounds[3] = std::max(bounds[3], childBounds[3]);
        bounds[5] = std::max(bounds[5], childBounds[5]);
    }
    const SBTreePositionCache& pc = getTreePositionCache(s);
    const SBModelVars& mv = getModelVars(s);
    Vec3 origin = body.getX_GB(pc).p();
    bounds[0] = std::min(bounds[0], origin[0]);
    bounds[2] = std::min(bounds[2], origin[1]);
    bounds[4] = std::min(bounds[4], origin[2]);
    bounds[1] = std::max(bounds[1], origin[0]);
    bounds[3] = std::max(bounds[3], origin[1]);
    bounds[5] = std::max(bounds[5], origin[2]);
    Array_<Vec3> corners;
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
        const SBTreePositionCache& tempCache = getTreePositionCache(tempState);
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
// DON'T USE THIS
void SimbodyMatterSubsystemRep::calcUUnitWeightsRecursively(const State& s, State& tempState, Vector& weights, Vec6& bounds, const RigidBodyNode& body) const {
    bounds[0] = bounds[2] = bounds[4] = Infinity;
    bounds[1] = bounds[3] = bounds[5] = -Infinity;
    
    // Call this method recursively on the body's children, and build up a bounding box
    // for everything downstream of it.
    
    for (int i = 0; i < body.getNumChildren(); ++i) {
        Vec6 childBounds;
        calcUUnitWeightsRecursively(s, tempState, weights, childBounds, *body.getChild(i));
        bounds[0] = std::min(bounds[0], childBounds[0]);
        bounds[2] = std::min(bounds[2], childBounds[2]);
        bounds[4] = std::min(bounds[4], childBounds[4]);
        bounds[1] = std::max(bounds[1], childBounds[1]);
        bounds[3] = std::max(bounds[3], childBounds[3]);
        bounds[5] = std::max(bounds[5], childBounds[5]);
    }
    const SBTreePositionCache& pc = getTreePositionCache(s);
    const SBModelVars& mv = getModelVars(s);
    Vec3 origin = body.getX_GB(pc).p();
    bounds[0] = std::min(bounds[0], origin[0]);
    bounds[2] = std::min(bounds[2], origin[1]);
    bounds[4] = std::min(bounds[4], origin[2]);
    bounds[1] = std::max(bounds[1], origin[0]);
    bounds[3] = std::max(bounds[3], origin[1]);
    bounds[5] = std::max(bounds[5], origin[2]);
    Array_<Vec3> corners;
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
        const SBTreePositionCache& tempCache = getTreePositionCache(tempState);
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

// We are in the process of realizeTopology() when we need to make this call.
// We pass in the partially-completed Topology-stage cache, which must have all
// the dimensions properly filled in at this point.
void SimbodyMatterSubsystemRep::setDefaultModelValues(const SBTopologyCache& topologyCache, 
                                                      SBModelVars& modelVars) const 
{
    // Tree-level defaults
    modelVars.useEulerAngles = false;

    assert((int)modelVars.prescribed.size() == getNumMobilizedBodies());
    modelVars.prescribed[0] = true; // Ground
    for (MobilizedBodyIndex i(1); i < getNumMobilizedBodies(); ++i)
        modelVars.prescribed[i] = false;

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultModelValues(topologyCache, modelVars);

}

void SimbodyMatterSubsystemRep::setDefaultInstanceValues(const SBModelVars& mv, 
                                                         SBInstanceVars& instanceVars) const 
{
    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setNodeDefaultInstanceValues(mv, instanceVars);

    assert((int)instanceVars.disabled.size() == getNumConstraints());
    for (ConstraintIndex i(0); i < getNumConstraints(); ++i)
        instanceVars.disabled[i] = getConstraint(i).isDisabledByDefault();

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

void SimbodyMatterSubsystemRep::setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void SimbodyMatterSubsystemRep::setMobilizerIsPrescribed(State& s, MobilizedBodyIndex body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void SimbodyMatterSubsystemRep::setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disable) const {
    SBInstanceVars& instanceVars = updInstanceVars(s); // check/adjust stage
    instanceVars.disabled[constraint] = disable;   
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
    const SBInstanceVars& instanceVars = getInstanceVars(s); // check stage
    return instanceVars.disabled[constraint];
}

void SimbodyMatterSubsystemRep::
convertToEulerAngles(const State& inputState, State& outputState) const {
    outputState = inputState;
    if (!getUseEulerAngles(inputState)) {
        setUseEulerAngles(outputState, true);
        getMultibodySystem().realizeModel(outputState);
        const Vector& inputQ  = inputState.getQ();
        Vector&       outputQ = outputState.updQ();
        for (int i = 0; i < getNumMobilizedBodies(); ++i) {
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
        for (int i = 0; i < getNumMobilizedBodies(); ++i) {
            const RigidBodyNode& node = getRigidBodyNode(i);
            node.convertToQuaternions(inputQ, outputQ);
        }
    }
}


bool SimbodyMatterSubsystemRep::isUsingQuaternion(const State& s, MobilizedBodyIndex body) const {
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
    return mc.getMobilizedBodyModelInfo(body).quaternionPoolIndex;
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

void SimbodyMatterSubsystemRep::calcHolonomicConstraintMatrixPNInv(const State& s, Matrix& PNInv) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nq = getNQ(s);

    PNInv.resize(mp,nq);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
            cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into qErr and mHolo (mp)

        PNInv(holoSeg.offset, 0, holoSeg.length, nq) = 
            constraints[cx]->calcPositionConstraintMatrixPNInv(s);
    }
}
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixP(const State& s, Matrix& P) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mHolo = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    P.resize(mHolo,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        P(holoSeg.offset, 0, holoSeg.length, nu) = 
            constraints[cx]->calcPositionConstraintMatrixP(s);
    }
}
void SimbodyMatterSubsystemRep::calcHolonomicVelocityConstraintMatrixPt(const State& s, Matrix& Pt) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mHolo = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Pt.resize(nu,mHolo);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg = cInfo.holoErrSegment; // offset into uErr and mHolo (mp)

        // Fill in columns of Pt
        Pt(0, holoSeg.offset, nu, holoSeg.length) = 
            constraints[cx]->calcPositionConstraintMatrixPt(s);
    }
}
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixV(const State& s, Matrix& V) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    V.resize(mNonholo,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        V(nonholoSeg.offset, 0, nonholoSeg.length, nu) = 
            constraints[cx]->calcVelocityConstraintMatrixV(s);
    }
}
void SimbodyMatterSubsystemRep::calcNonholonomicConstraintMatrixVt(const State& s, Matrix& Vt) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int nu = getNU(s);

    Vt.resize(nu,mNonholo);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& nonholoSeg = cInfo.nonholoErrSegment; // after holo derivs, offset into uerr

        // Fill in columns of Vt
        Vt(0, nonholoSeg.offset, nu, nonholoSeg.length) = 
            constraints[cx]->calcVelocityConstraintMatrixVt(s);
    }
}
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixA (const State& s, Matrix& A) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    A.resize(mAccOnly,nu);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        A(accOnlySeg.offset, 0, accOnlySeg.length, nu) = 
            constraints[cx]->calcAccelerationConstraintMatrixA(s);
    }
}
void SimbodyMatterSubsystemRep::calcAccelerationOnlyConstraintMatrixAt(const State& s, Matrix& At) const {
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int nu = getNU(s); // also nudot

    At.resize(nu,mAccOnly);

    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& 
                       cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment; // after holo&nonholo derivs, offset into udoterr

        At(0, accOnlySeg.offset, nu, accOnlySeg.length) = 
            constraints[cx]->calcAccelerationConstraintMatrixAt(s);
    }
}



//------------------------------------------------------------------------------
//                  CALC CONSTRAINT FORCES FROM MULTIPLIERS
//------------------------------------------------------------------------------
// Must be realized to Stage::Position.
// Note that constraint forces have the opposite sign from applied forces,
// because we calculate multipliers from
//    M udot + ~G lambda = f_applied
// If you want to view the constraint generated forces as though they
// were applied forces, negate lambda before the call here.
void SimbodyMatterSubsystemRep::calcConstraintForcesFromMultipliers
  (const State& s, const Vector& lambda,
   Vector_<SpatialVec>& bodyForcesInG,
   Vector&              mobilityForces) const
{
    const SBInstanceCache& ic = getInstanceCache(s);
    const int mHolo    = ic.totalNHolonomicConstraintEquationsInUse;
    const int mNonholo = ic.totalNNonholonomicConstraintEquationsInUse;
    const int mAccOnly = ic.totalNAccelerationOnlyConstraintEquationsInUse;
    const int ma = mHolo+mNonholo+mAccOnly;

    assert(lambda.size() == ma);

    bodyForcesInG.resize(getNumBodies()); bodyForcesInG.setToZero();
    mobilityForces.resize(getNU(s));      mobilityForces.setToZero();

    Vector_<SpatialVec> bodyF1;          // per constraint
    Vector              mobilityF1;
    Vector              lambda1; // multipliers for 1 constraint
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;

        const SBInstanceCache::PerConstraintInstanceInfo& cInfo = ic.getConstraintInstanceInfo(cx);
        const Segment& holoSeg    = cInfo.holoErrSegment;
        const Segment& nonholoSeg = cInfo.nonholoErrSegment;
        const Segment& accOnlySeg = cInfo.accOnlyErrSegment;
        const int mh=holoSeg.length, mnh=nonholoSeg.length, mao=accOnlySeg.length;

        lambda1.resize(mh+mnh+mao);
        lambda1(0,mh)        = lambda(holoSeg.offset, mh);
        lambda1(mh, mnh)     = lambda(mHolo+nonholoSeg.offset, mnh);
        lambda1(mh+mnh, mao) = lambda(mHolo+mNonholo+accOnlySeg.offset, mao);

        const ConstraintImpl& crep = constraints[cx]->getImpl();
 
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
    matter.multiplyByNInv(s, false, qErrest, qhatErrest); // qhatErrest = N+ qErrest
    qhatErrest.rowScaleInPlace(uWeights);                 // qhatErrest = Wu N+ qErrest
    return qhatErrest.normRMS();
}


// -----------------------------------------------------------------------------
//                                PRESCRIBE
// -----------------------------------------------------------------------------
// This is a solver that sets continuous state variables q, or u (depending 
// on stage) to their prescribed values that will already have been computed. 
// Note that prescribed udot=udot(t,q,u) is not dealt with here because it does 
// not involve a state change.
bool SimbodyMatterSubsystemRep::prescribe(State& s, Stage g) const {
    const SBModelCache&    mc = getModelCache(s);
    const SBInstanceCache& ic = getInstanceCache(s);
    switch(g) {

    // Prescribe position.
    case Stage::PositionIndex: {
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
    } break;

    // Prescribe velocity.
    case Stage::VelocityIndex: {
        const int npu = ic.getTotalNumPresU();
        const int nzu = ic.getTotalNumZeroU();
        if (npu==0 && nzu==0) return false; // don't invalidate positions

        // copy prescribed u's from cache to state
        // set known-zero u's to zero
        const SBConstrainedPositionCache& cpc = getConstrainedPositionCache(s);
        Vector& u = updU(s); // this Subsystem's u's, now invalidated

        for (int i=0; i < npu; ++i)
            u[ic.presU[i]] = cpc.presUPool[i];

        for (int i=0; i < nzu; ++i)
            u[ic.zeroU[i]] = 0;
    } break;

    default:
        SimTK_ASSERT1_ALWAYS(!"bad stage",
            "SimbodyMatterSubsystemRep::prescribe(): bad stage argument %s.", 
            g.getName().c_str());
    }

    return true;
}
// .............................. PRESCRIBE ....................................



//------------------------------------------------------------------------------
//                       ENFORCE POSITION CONSTRAINTS
//------------------------------------------------------------------------------
void SimbodyMatterSubsystemRep::enforcePositionConstraints
   (State& s, Real consAccuracy, const Vector& yWeights,
    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{
    assert(getStage(s) >= Stage::Position-1);
    realizeSubsystemPosition(s);

    // First work only with the holonomic (position) constraints, which appear 
    // first in the QErr array. Don't work on the quaternion constraints in 
    // this first section.
    const int mHolo  = getNumHolonomicConstraintEquationsInUse(s);
    const int mQuats = getNumQuaternionsInUse(s);
    const int nq     = getNQ(s);
    const int nu     = getNU(s);

    const VectorView uWeights   = yWeights(nq,nu);
    const Vector     ooUWeights = uWeights.elementwiseInvert();
    //Real oow[] = {1., 1., 1., .05, .05, .05};
   // Vector     ooUWeights(nu, oow);
    const VectorView ooPTols  = ooTols(0,mHolo);

    VectorView qErrest = yErrest.size() ? yErrest(0,nq) : yErrest(0,0);

    // This is a const view into the State; the contents it refers to will 
    // change though.
    const VectorView pErrs = getQErr(s)(0,mHolo); // just leave off quaternions

    bool anyChange = false;

    // Check whether we should stop if we see the solution diverging
    // which should not happen when we're in the neighborhood of a solution
    // on entry. This is always set while integrating, except during
    // initialization.
    const bool localOnly = opts.isOptionSet(System::ProjectOptions::LocalOnly);

    // Solve 
    //   (Tp P Wu^-1) dqhat_WLS = T perr, q -= N*Wu^-1*dqhat_WLS
    // until perr(q)_TRMS <= 0.1*accuracy.
    //
    // This is a nonlinear least squares problem. This is a full Newton 
    // iteration since we recalculate the iteration matrix each time around the
    // loop. TODO: a solution could be found using the same iteration matrix, 
    // since we are projecting from (presumably) not too far away. Q1: will it
    // be the same solution? Q2: if not, does it matter?
    Vector scaledPerrs = pErrs.rowScale(ooPTols);
    Real normAchievedTRMS = scaledPerrs.normRMS();

    
    //cout << "!!!! initially @" << s.getTime() << ", perr TRMS=" 
    //     << normAchievedTRMS << " consAcc=" << consAccuracy;
    //if (qErrest.size())
    //    cout << " qErrest WRMS=" 
    //         << calcQErrestWeightedNorm(*this,s,qErrest,uWeights);
    //else cout << " NO Q ERROR ESTIMATE";
    //cout << endl;
    //cout << "!!!! PERR=" << pErrs << endl;
 

    Real lastChangeMadeWRMS = 0; // size of last change in weighted dq
    int nItsUsed = 0;

    // Set how far past the required tolerance we'll attempt to go. 
    // We only fail if we can't achieve consAccuracy, but while we're
    // solving we'll see if we can get consAccuracyToTryFor.
    const Real consAccuracyToTryFor = 
        std::max(0.1*consAccuracy, SignificantReal);

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is much too tight; should depend on constraint tolerance
    // and rank should be saved and reused in velocity and acceleration
    // constraint methods (or should be calculated elsewhere and passed in).
    const Real conditioningTol = mHolo         
      //* SignificantReal;
        * SqrtEps;
    //cout << "  Conditioning tolerance=" << conditioningTol << endl;

    if (normAchievedTRMS > consAccuracyToTryFor) {
        Vector saveQ = getQ(s);
        Matrix Pt(nu,mHolo);
        Vector dqhat_WLS(nu);
        Vector dq(nq); // = N Wu^-1 dqhat_WLS
        FactorQTZ P_qtz;
        Real prevNormAchievedTRMS = normAchievedTRMS; // watch for divergence
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

            P_qtz.factor<Real>(~Pt, conditioningTol); // this acts like a pseudoinverse

            //std::cout << "POSITION PROJECTION TOL=" << conditioningTol
            //          << " RANK=" << P_qtz.getRank() 
            //          << " RCOND=" << P_qtz.getRCondEstimate() << std::endl;

            P_qtz.solve(scaledPerrs, dqhat_WLS);
            lastChangeMadeWRMS = dqhat_WLS.normRMS(); // size of change in weighted norm
            //cout << "!!!! dqhat weighted=" << dqhat_WLS << endl;

            // switch back to unweighted dqhat=Wu^-1*dqhat_WLS
            dqhat_WLS.rowScaleInPlace(ooUWeights);
            //cout << "!!!! dqhat unweighted=" << dqhat_WLS << endl;

            multiplyByN(s, false, dqhat_WLS, dq); // N*(Wu^-1 dqhat_WLS)
            //cout << "!!!! dq unweighted=" << dq << endl;
            updQ(s) -= dq; // this is actually unweighted dq
            anyChange = true;

            // Now recalculate the position constraint errors at the new q.
            realizeSubsystemPosition(s);
            //cout << "!!!! PERR=" << pErrs << endl;
            scaledPerrs = pErrs.rowScale(ooPTols);
            normAchievedTRMS = scaledPerrs.normRMS();
            ++nItsUsed;

            //std::cout << "  !! iter " << nItsUsed << ": TRMS(perr)=" 
            //    << normAchievedTRMS << ", WRMS(dq)=" 
            //    << lastChangeMadeWRMS << endl;

            if (localOnly && nItsUsed >= 2 
                && normAchievedTRMS > prevNormAchievedTRMS) {
                //std::cout << "   POS NORM GOT WORSE -- STOP\n";
                // restore to end of previous iteration
                updQ(s) += dq;
                realizeSubsystemPosition(s);
                scaledPerrs = pErrs.rowScale(ooPTols);
                normAchievedTRMS = scaledPerrs.normRMS();
                break; // diverging -- quit now to prevent a bad solution
            }


            prevNormAchievedTRMS = normAchievedTRMS;

        } while (normAchievedTRMS > consAccuracyToTryFor
                 && nItsUsed < MaxIterations);

        // Make sure we achieved at least the required constraint accuracy.
        if (normAchievedTRMS > consAccuracy) {
            updQ(s) = saveQ;
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "Failed to converge in position projection");
        }


        // Next, if we projected out the position constraint errors, remove the
        // corresponding error from the integrator's error estimate.
        //   (Tp P Wu^-1) dqhat_WLS = (Tp P Wu^-1) Wu N^+ qErrest, 
        //                                  qErrest -= N Wu^-1 dqhat_WLS
        // No iteration is required.
        if (qErrest.size()) {
            // Work in qhat W-norm
            Vector qhatErrest(nu), dqErrest(nq);
            multiplyByNInv(s, false, qErrest, qhatErrest);
            qhatErrest.rowScaleInPlace(uWeights); // qbarErrest = Wu * N^+ * qErrest

            P_qtz.solve(~Pt*qhatErrest, dqhat_WLS);
            const Real normOfAdjustment_WRMS = dqhat_WLS.normRMS();

            dqhat_WLS.rowScaleInPlace(ooUWeights); // unscale the result
            multiplyByN(s, false, dqhat_WLS, dqErrest);
            qErrest -= dqErrest; // unweighted
            //cout << "  !! FIXUP: now WRMS(qerrest)=" 
            //     << calcQErrestWeightedNorm(*this,s,qErrest,uWeights) 
            //     << " using WRMS(dq)=" << normOfAdjustment_WRMS << endl;
        }
    }

    //cout << "!!!! perr TRMS achieved " << normAchievedTRMS << " in " 
    //     << nItsUsed << " iterations"  << endl;
    //cout << "!!!! ... PERR=" << pErrs << endl;

    // By design, normalization of quaternions can't have any effect on the 
    // length constraints we just fixed (because we normalize internally for 
    // calculations). So now we can simply normalize the quaternions.
    SBStateDigest sbs(s, *this, Stage::Model);

    if (mQuats) {
        //cout << "!!!! QUAT START: errs=" << getQErr(s)(mHolo, mQuats) 
        //     << " RMS(qErrest)=" << qErrest.normRMS() << endl;
        Vector& q  = updQ(s); //TODO: this invalidates q's already

        for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
            for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
                if (rbNodeLevels[i][j]->enforceQuaternionConstraints(sbs,q,qErrest))
                    anyChange = true;

        //TODO: shouldn't need this
        realizeSubsystemPosition(s);

        //cout << "!!!! ... QUAT END: errs=" << getQErr(s)(mHolo, mQuats) 
        //     << " RMS(qErrest)=" << qErrest.normRMS() << endl;
        //TODO: quaternion constraints shouldn't invalidate anything except
        // the qnorms, which will be all 1 now
    }

    if (anyChange)
        s.invalidateAll(Stage::Position);
}



//------------------------------------------------------------------------------
//                          ENFORCE VELOCITY CONSTRAINTS
//------------------------------------------------------------------------------
void SimbodyMatterSubsystemRep::enforceVelocityConstraints
   (State& s, Real consAccuracy, const Vector& yWeights,
    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{
    assert(getStage(s) >= Stage::Velocity-1);
    realizeSubsystemVelocity(s);

    // Here we deal with the nonholonomic (velocity) constraints and the 
    // derivatives of the holonomic constraints.
    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    const VectorView uWeights   = yWeights(nq,nu); // skip the first nq weights
    const Vector     ooUWeights = uWeights.elementwiseInvert();
    const VectorView ooPVTols   = ooTols(0,mHolo+mNonholo);
    //TODO: scale holo part by time scale

    VectorView uErrest = yErrest.size() ? yErrest(nq,nu) : yErrest(0,0);

    // This is a const view into the State; the contents it refers to will 
    // change though.
    const Vector& vErrs = getUErr(s); // all velocity constraint errors (mHolo+mNonholo)

    bool anyChange = false;

    // Solve 
    //   (Tpv [P;V] Wu^-1) du_WLS = Tpv uerr, u -= Wu^-1*du_WLS
    // Note that although this is a nonlinear least squares problem since uerr 
    // is a function of u, we do not need to refactor the matrix since it does 
    // not depend on u.
    // TODO: Tp P Wu^-1 should already have been calculated for position 
    // projection (at least if any position projection occurred)
    //
    // This is a nonlinear least squares problem, but we only need to factor 
    // once since only the RHS is dependent on u. 
    Vector scaledVerrs = vErrs.rowScale(ooPVTols);
    Real normAchievedTRMS = scaledVerrs.normRMS();
    
    //
    //cout << "!!!! initially @" << s.getTime() << ", verr TRMS=" 
    //     << normAchievedTRMS << " consAcc=" << consAccuracy;
    //if (uErrest.size())
    //    cout << " uErrest WRMS=" << uErrest.rowScale(uWeights).normRMS();
    //else cout << " NO U ERROR ESTIMATE";
    //cout << endl;
    //cout << "!!!! VERR=" << vErrs << endl;
    

    Real lastChangeMadeWRMS = 0;
    int nItsUsed = 0;

    // Check whether we should stop if we see the solution diverging
    // which should not happen when we're in the neighborhood of a solution
    // on entry.
    const bool localOnly = opts.isOptionSet(System::ProjectOptions::LocalOnly);

    // Set how far past the required tolerance we'll attempt to go. 
    // We only fail if we can't achieve consAccuracy, but while we're
    // solving we'll see if we can get consAccuracyToTryFor.
    const Real consAccuracyToTryFor = 
        std::max(0.1*consAccuracy, SignificantReal);

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is much too tight; should depend on constraint tolerance
    // and should be consistent with the holonomic rank.
    const Real conditioningTol = (mHolo+mNonholo) 
        //* SignificantReal;
        * SqrtEps;
    //cout << "  Conditioning tolerance=" << conditioningTol << endl;

    if (normAchievedTRMS > consAccuracyToTryFor) {
        const Vector saveU = getU(s);
        Matrix PVt(nu, mHolo+mNonholo);

        calcHolonomicVelocityConstraintMatrixPt(s, PVt(  0,  0,nu, mHolo));    // nu X mp
        calcNonholonomicConstraintMatrixVt     (s, PVt(0,mHolo,nu, mNonholo)); // nu X mv
        PVt.rowAndColScaleInPlace(ooUWeights, ooPVTols); // PVt is now Wu^-1 (Pt Vt) Tpv

        // Calculate pseudoinverse (just once)
        FactorQTZ PVqtz;
        PVqtz.factor<Real>(~PVt, conditioningTol);

        //std::cout << "VELOCITY PROJECTION TOL=" << conditioningTol
        //          << " RANK=" << PVqtz.getRank() 
        //          << " RCOND=" << PVqtz.getRCondEstimate() << std::endl;

        Vector du_WLS(nu);
        Real prevNormAchievedTRMS = normAchievedTRMS; // watch for divergence
        const int MaxIterations  = 7;
        do {
            PVqtz.solve(scaledVerrs, du_WLS);
            lastChangeMadeWRMS = du_WLS.normRMS(); // change size in weighted norm
            du_WLS.rowScaleInPlace(ooUWeights); // remove scaling: du=Wu^-1*du_WLS
            updU(s) -= du_WLS; // this is actually unweighted du
            anyChange = true;

            // Recalculate the constraint errors for the new u's.
            realizeSubsystemVelocity(s);
            scaledVerrs = vErrs.rowScale(ooPVTols);
            normAchievedTRMS=scaledVerrs.normRMS();
            ++nItsUsed;

            //cout << "  !! iter " << nItsUsed << ": TRMS(verr)=" 
            //     << normAchievedTRMS << ", WRMS(du)=" 
            //     << lastChangeMadeWRMS << endl;

            if (localOnly && nItsUsed >= 2 
                && normAchievedTRMS > prevNormAchievedTRMS) {
                //std::cout << "   VEL NORM GOT WORSE -- STOP\n";
                // restore to end of previous iteration
                updU(s) += du_WLS;
                realizeSubsystemVelocity(s);
                scaledVerrs = vErrs.rowScale(ooPVTols);
                normAchievedTRMS=scaledVerrs.normRMS();
                break; // diverging -- quit now to prevent a bad solution
            }

            prevNormAchievedTRMS = normAchievedTRMS;

        } while (normAchievedTRMS > consAccuracyToTryFor
                 && nItsUsed < MaxIterations);

        // Make sure we achieved at least the required constraint accuracy.
        if (normAchievedTRMS > consAccuracy) {
            updU(s) = saveU;
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "Failed to converge in velocity projection");
        }

        // Next, if we projected out the velocity constraint errors, remove the
        // corresponding error from the integrator's error estimate.
        //   (Tpv [P;V] Wu^-1) du_WLS = (Tpv [P;V] Wu^-1) Wu*uErrEst, 
        //                                  uErrEst = Wu^-1(Wu*uErrEst - du_WLS)
        // No iteration is required.
        if (uErrest.size()) {
            uErrest.rowScaleInPlace(uWeights); // uErrest = Wu*uErrEst
            PVqtz.solve(~PVt*uErrest, du_WLS);
            uErrest -= du_WLS; // still weighted
            //cout << "  !! U FIXUP: now WRMS(uerrest)=" << uErrest.normRMS() 
            //     << " using WRMS(du)=" << du_WLS.normRMS() << endl;
            uErrest.rowScaleInPlace(ooUWeights); // back to unscaled err. est.
        }
    }

    
    //cout << "!!!! verr achieved " << normAchievedTRMS << " in " 
    //     << nItsUsed << " iterations" << endl;
    //cout << "!!!! ... VERR=" << vErrs << endl;

    if (anyChange)
        s.invalidateAll(Stage::Velocity);
}



//------------------------------------------------------------------------------
//                        CALC TREE FORWARD DYNAMICS
//------------------------------------------------------------------------------
//
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
    SBTreeAccelerationCache&        tac,  // accelerations and prescribed forces go here
    Vector&                         udot, // in/out (in for prescribed udot)
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
    Vector_<SpatialVec>& A_GB           = tac.bodyAccelerationInGround;
    Vector&              tau            = tac.presMotionForces;

    calcTreeAccelerations(s, *mobilityForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot, tau);

    // Ask the constraints to calculate ancestor-relative acceleration 
    // kinematics (still goes in TreeAccelerationCache).
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl()
            .calcConstrainedBodyAccelerationInAncestor(sbs, tac);


    // TODO: we have to say we're done with the TreeAccelerationCache in the 
    // State but that is wrong if "tac" is not from the State. So this alleged 
    // "operator" currently can be used only for realization of the State; i.e.,
    // the passed-in TreeAccelerationCache better have come from inside this 
    // State! Need a different design for the constraint interface to support 
    // operators.
    markCacheValueRealized(s, mc.treeAccelerationCacheIndex);

    // Put acceleration constraint equation errors in udotErr
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx) {
        if (isConstraintDisabled(s,cx))
            continue;
        const SBInstanceCache::PerConstraintInstanceInfo& cInfo = ic.getConstraintInstanceInfo(cx);

        const Segment& holoseg    = cInfo.holoErrSegment;    // for 2nd derivatives of holonomic constraints
        const Segment& nonholoseg = cInfo.nonholoErrSegment; // for 1st derivatives of nonholonomic constraints
        const Segment& acconlyseg = cInfo.accOnlyErrSegment; // for acceleration-only constraints
        const int mHolo = holoseg.length, mNonholo = nonholoseg.length, mAccOnly = acconlyseg.length;
        if (mHolo)
            constraints[cx]->getImpl().realizePositionDotDotErrors(s, mHolo,
                &udotErr[holoseg.offset]);
        if (mNonholo)
            constraints[cx]->getImpl().realizeVelocityDotErrors(s, mNonholo, 
                &udotErr[  ic.totalNHolonomicConstraintEquationsInUse 
                         + nonholoseg.offset]);
        if (mAccOnly)
            constraints[cx]->getImpl().realizeAccelerationErrors(s, mAccOnly, 
                &udotErr[  ic.totalNHolonomicConstraintEquationsInUse
                         + ic.totalNNonholonomicConstraintEquationsInUse 
                         + acconlyseg.offset]);
    }
}

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
    Vector&                         udotErr = updUDotErr(s);

    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    extraMobilityForces, extraBodyForces,
                                    tac, udot, udotErr);
}



//------------------------------------------------------------------------------
//                          CALC LOOP FORWARD DYNAMICS
//------------------------------------------------------------------------------
// 
// Given a State realized through Stage::Dynamics, and a complete set of applied 
// forces, calculate all acceleration results resulting from those forces AND 
// enforcement of the acceleration constraints. The results go into the return 
// arguments here. This routine *does not* affect the State cache -- it is an 
// operator. In typical usage, the output arguments actually will be part of 
// the state cache to effect a response, but this method can also be used to 
// effect an operator.
void SimbodyMatterSubsystemRep::calcLoopForwardDynamicsOperator(const State& s, 
    const Vector&                   mobilityForces,
    const Vector_<Vec3>&            particleForces,
    const Vector_<SpatialVec>&      bodyForces,
    SBTreeAccelerationCache&        tac,
    Vector&                         udot,
    Vector&                         multipliers,
    Vector&                         udotErr) const
{
    assert(getStage(s) >= Stage::Acceleration-1);

    // Calculate acceleration results ignoring Constraints, except to have
    // them calculate the resulting constraint errors.
    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    0, 0, tac, udot, udotErr);

    // Next, determine how many acceleration-level constraint equations 
    // need to be obeyed.

    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = getNumAccelerationOnlyConstraintEquationsInUse(s);
    const int ma = mHolo+mNonholo+mAccOnly;
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    if (ma==0 || nu==0) {
        multipliers.resize(0);
        return;
    }

    multipliers.resize(ma);

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is probably too tight; should depend on constraint tolerance
    // and should be consistent with position and velocity projection ranks.
    // Tricky here because conditioning depends on mass matrix as well as
    // constraints.
    const Real conditioningTol = ma 
        //* SignificantReal;
        * SqrtEps*std::sqrt(SqrtEps); // Eps^(3/4)

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
    // if we've partitioned it into little subblocks, this is all very 
    // reasonable. One slip up and you'll toss in a factor of mn^2 or m^2n and
    // screw this up -- be careful! (see below)
    Matrix MInvGt(nu, ma);
    Vector_<SpatialVec> A_GB(getNumBodies()); // dummy
    for (int j=0; j<ma; ++j) // This is O(mn)
        calcMInverseF(s, Gt(j), A_GB, MInvGt(j));

    // TODO: Toldya! Check out this m^2n bit here ...
    Matrix GMInvGt = (~Gt)*MInvGt; // TODO: BAD!!! O(m^2n) -- Use m x G udot operators instead for O(mn)
    FactorQTZ qtz(GMInvGt, conditioningTol); // specify 1/cond at which we declare rank deficiency
    qtz.solve(udotErr, multipliers);

        //std::cout << "MULTIPLIER SOLVE TOL=" << conditioningTol
        //          << " RANK=" << qtz.getRank() 
        //          << " RCOND=" << qtz.getRCondEstimate() << std::endl;
        //std::cout << " multipliers=" << multipliers << std::endl;

    // We have the multipliers, now turn them into forces.

    Vector_<SpatialVec> bodyF;
    Vector mobilityF;
    calcConstraintForcesFromMultipliers(s,multipliers,bodyF,mobilityF);
    // Note that constraint forces have the opposite sign from applied forces
    // so must be subtracted to calculate the total forces.

    // Recalculate the accelerations applying the constraint forces in addition
    // to the applied forces that were passed in. The constraint errors calculated
    // now should be within noise of zero.

    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    &mobilityF, &bodyF, tac, udot, udotErr);
}


// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void SimbodyMatterSubsystemRep::realizeLoopForwardDynamics(const State& s, 
    const Vector&               mobilityForces,
    const Vector_<Vec3>&        particleForces,
    const Vector_<SpatialVec>&  bodyForces) const 
{
    // Because we are realizing, we want to direct the output of the operator
    // back into the State cache.
    SBTreeAccelerationCache&        tac         = updTreeAccelerationCache(s);
    Vector&                         udot        = updUDot(s);
    Vector&                         udotErr     = updUDotErr(s);
    Vector&                         multipliers = updMultipliers(s);

    calcLoopForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    tac, udot, multipliers, udotErr);
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
    Vector_<SpatialVec> bodyForces(getNumBodies());
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


// -----------------------------------------------------------------------------
//                        CALC COMPOSITE BODY INERTIAS
// -----------------------------------------------------------------------------
// Given a State realized to Position stage, calculate the composite
// body inertias seen by each mobilizer. A composite body inertia is
// the inertia of the rigid body created by locking all joints outboard
// of a particular mobilized body. (Constraints have no effect on the result.)
//
void SimbodyMatterSubsystemRep::calcCompositeBodyInertias(const State& s,
    Array_<SpatialInertia>& R) const
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    R.resize(getNumBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcCompositeBodyInertiasInward(tpc,R);
}
//....................... CALC COMPOSITE BODY INERTIAS .........................



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

// Process forces for subsequent use by realizeTreeAccel() below.
void SimbodyMatterSubsystemRep::realizeZ(const State& s, 
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces) const
{
    const SBStateDigest sbs(s, *this, Stage::Acceleration);
    const SBTreePositionCache&  tpc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  tvc = sbs.getTreeVelocityCache();
    const SBDynamicsCache&      dc  = sbs.getDynamicsCache();
    SBTreeAccelerationCache&    tac = sbs.updTreeAccelerationCache();

    const SBArticulatedBodyInertiaCache&    
            abc = getArticulatedBodyInertiaCache(s);

    // TODO: does this need to do level 0?
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.realizeZ(tpc,abc,tvc,dc,tac,mobilityForces,bodyForces);
        }
}

// Calc acceleration: sweep from base to tip. This uses the forces
// that were last supplied to realizeZ()above.
void SimbodyMatterSubsystemRep::realizeTreeAccel(const State& s) const {

    SBStateDigest sbs(s, *this, Stage::Acceleration);
    const SBTreePositionCache&  tpc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  tvc = sbs.getTreeVelocityCache();
    const SBDynamicsCache&      dc  = sbs.getDynamicsCache();
    SBTreeAccelerationCache&    tac = sbs.updTreeAccelerationCache();

    const SBArticulatedBodyInertiaCache&    
            abc     = getArticulatedBodyInertiaCache(s);
    Vector& udot    = updUDot(s);
    Vector& qdotdot = updQDotDot(s);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            rbNodeLevels[i][j]->realizeAccel(tpc,abc,tvc,dc,tac,udot);
            rbNodeLevels[i][j]->calcQDotDot(sbs, udot, qdotdot);
        }
}



//------------------------------------------------------------------------------
//                           CALC KINETIC ENERGY
//------------------------------------------------------------------------------
Real SimbodyMatterSubsystemRep::calcKineticEnergy(const State& s) const {
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const SBTreeVelocityCache& tvc = getTreeVelocityCache(s);

    Real ke = 0.;

    // Skip ground level 0!
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            ke += rbNodeLevels[i][j]->calcKineticEnergy(tpc,tvc);

    return ke;
}



//------------------------------------------------------------------------------
//                          CALC TREE ACCELERATIONS
//------------------------------------------------------------------------------
// Operator for open-loop forward dynamics.
// This Subsystem must have already been realized to Dynamics stage so that 
// dynamics quantities like articulated body inertias are available.
void SimbodyMatterSubsystemRep::calcTreeAccelerations(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot,    // in/out (in for prescribed udots)
    Vector&                    tau) const 
{
    const SBInstanceCache&               ic  = getInstanceCache(s);
    const SBTreePositionCache&           tpc = getTreePositionCache(s);
    const SBTreeVelocityCache&           tvc = getTreeVelocityCache(s);
    const SBDynamicsCache&               dc  = getDynamicsCache(s);
    const SBArticulatedBodyInertiaCache& abc = getArticulatedBodyInertiaCache(s);

    assert(mobilityForces.size() == getTotalDOF());
    assert(bodyForces.size() == getNumBodies());

    netHingeForces.resize(getTotalDOF());
    A_GB.resize(getNumBodies());
    udot.resize(getTotalDOF());
    tau.resize(ic.totalNPresForce);

    // Temporaries
    Vector_<SpatialVec> allZ(getNumBodies());
    Vector_<SpatialVec> allGepsilon(getNumBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass1Inward(ic,tpc,abc,dc,
                mobilityForces, bodyForces, udot, allZ, allGepsilon,
                netHingeForces);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass2Outward(ic,tpc,abc,tvc,dc, 
                                      netHingeForces, A_GB, udot, tau);
        }
}



//------------------------------------------------------------------------------
//                           CALC M INVERSE F
//------------------------------------------------------------------------------
// Calculate udot = M^-1 f. We also get spatial accelerations A_GB for 
// each body as a side effect.
// This Subsystem must already be realized through Dynamics stage.
void SimbodyMatterSubsystemRep::calcMInverseF(const State& s,
    const Vector&              f,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBInstanceCache&                  ic  = getInstanceCache(s);
    const SBTreePositionCache&              tpc = getTreePositionCache(s);
    const SBDynamicsCache&                  dc  = getDynamicsCache(s);
    const SBArticulatedBodyInertiaCache&    abc = getArticulatedBodyInertiaCache(s);

    assert(f.size() == getTotalDOF());

    A_GB.resize(getNumBodies());
    udot.resize(getTotalDOF());

    // Temporaries
    Vector              allEpsilon(getTotalDOF());
    Vector_<SpatialVec> allZ(getNumBodies());
    Vector_<SpatialVec> allGepsilon(getNumBodies());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass1Inward(ic,tpc,abc,dc,
                f, allZ, allGepsilon,
                allEpsilon);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass2Outward(ic,tpc,abc,dc, allEpsilon, A_GB, udot);
        }
}



//------------------------------------------------------------------------------
//                          CALC TREE RESIDUAL FORCES
//------------------------------------------------------------------------------
// Operator for tree system inverse dynamics. 
// Note that this includes the effects of inertial forces.
// This Subsystem must already have been realized to Velocity stage.
void SimbodyMatterSubsystemRep::calcTreeResidualForces(const State& s,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    residualMobilityForces) const
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);
    const SBTreeVelocityCache& tvc = getTreeVelocityCache(s);

    // We allow the input Vectors to be zero length. For now we have
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
        if (appliedMobilityForces.size()==0) pAppliedMobForces = &zeroPerMobility;
        if (knownUdot.size()==0)             pKnownUdot        = &zeroPerMobility;
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

    // Allocate temporary.
    Vector_<SpatialVec> allFTmp(getNumBodies());

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInverseDynamicsPass1Outward(tpc,tvc,*pKnownUdot,A_GB);
        }

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInverseDynamicsPass2Inward(
                tpc,tvc,A_GB,
                *pAppliedMobForces,*pAppliedBodyForces,
                allFTmp,residualMobilityForces);
        }
}



//------------------------------------------------------------------------------
//                                 CALC M V
//------------------------------------------------------------------------------
// Calculate x = M v. If the vector v is a generalized acceleration
// udot, then we also get spatial accelerations A_GB for 
// each body as a side effect.
// This Subsystem must already have been realized to Position stage.
void SimbodyMatterSubsystemRep::calcMV(const State& s,
    const Vector&              v,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    f) const 
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);

    assert(v.size() == getTotalDOF());

    A_GB.resize(getNumBodies());
    f.resize(getTotalDOF());

    // Temporary
    Vector_<SpatialVec> allFTmp(getNumBodies());

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMVPass1Outward(tpc, v, A_GB);
        }

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMVPass2Inward(tpc,A_GB,allFTmp,f);
        }
}

//------------------------------------------------------------------------------
//                                  CALC M
//------------------------------------------------------------------------------
// Calculate the mass matrix M in O(n^2) time. This Subsystem must already have
// been realized to Position stage.
void SimbodyMatterSubsystemRep::calcM(const State& s, Matrix& M) const {
    const int nu = getTotalDOF();
    M.resize(nu,nu);

    // This could be calculated much faster by doing it directly and calculating
    // only half of it. As a placeholder, however, we're doing this with 
    // repeated O(n) calls to calcMV() to get M one column at a time.
    Vector_<SpatialVec> A_GB(getNumBodies()); // unused dummy needed
    Vector v(nu); v = 0;
    for (int i=0; i < nu; ++i) {
        v[i] = 1;
        calcMV(s, v, A_GB, M(i));
        v[i] = 0;
    }
}

//------------------------------------------------------------------------------
//                                CALC MInv
//------------------------------------------------------------------------------
// Calculate the mass matrix inverse MInv(=M^-1) in O(n^2) time. This Subsystem
// must already have been realized to Position stage.
void SimbodyMatterSubsystemRep::calcMInv(const State& s, Matrix& MInv) const {
    const int nu = getTotalDOF();
    MInv.resize(nu,nu);

    // This could probably be calculated faster by doing it directly and
    // filling in only half. For now we're doing it with repeated calls to
    // the O(n) operator calcMInverseF().
    Vector_<SpatialVec> A_GB(getNumBodies()); // unused dummy needed
    Vector f(nu); f = 0;
    for (int i=0; i < nu; ++i) {
        f[i] = 1;
        calcMInverseF(s, f, A_GB, MInv(i));
        f[i] = 0;
    }
}

//------------------------------------------------------------------------------
//                               MULTIPLY BY N
//------------------------------------------------------------------------------
// q=Nu or u=~Nq
void SimbodyMatterSubsystemRep::multiplyByN
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Position
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next());

    assert(in.size() == (transpose?getTotalQAlloc():getTotalDOF()));
    out.resize(transpose?getTotalDOF():getTotalQAlloc());

    //TODO: this shouldn't be necessary
    assert(in.size()  < 2 || &in[1]  == &in[0] +1); // for now must be contiguous in memory
    assert(out.size() < 2 || &out[1] == &out[0]+1); 

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



//------------------------------------------------------------------------------
//                              MULTIPLY BY NDOT
//------------------------------------------------------------------------------
// q=NDot*u or u=~NDot*q
void SimbodyMatterSubsystemRep::multiplyByNDot
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Velocity
    const SBStateDigest sbState(s, *this, Stage(Stage::Velocity).next());

    assert(in.size() == (transpose?getTotalQAlloc():getTotalDOF()));
    out.resize(transpose?getTotalDOF():getTotalQAlloc());

    //TODO: this shouldn't be necessary
    assert(in.size()  < 2 || &in[1]  == &in[0] +1); // for now must be contiguous in memory
    assert(out.size() < 2 || &out[1] == &out[0]+1); 

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



//------------------------------------------------------------------------------
//                              MULTIPLY BY NINV
//------------------------------------------------------------------------------
// u= NInv * q or q = ~NInv * u
void SimbodyMatterSubsystemRep::multiplyByNInv
   (const State& s, bool transpose, const Vector& in, Vector& out) const
{
    // i.e., we must be *done* with Stage::Position
    const SBStateDigest sbState(s, *this, Stage(Stage::Position).next());

    assert(in.size() == (transpose?getTotalDOF():getTotalQAlloc()));
    out.resize(transpose?getTotalQAlloc():getTotalDOF());

    //TODO: this shouldn't be necessary
    assert(in.size()  < 2 || &in[1]  == &in[0] +1); // for now must be contiguous in memory
    assert(out.size() < 2 || &out[1] == &out[0]+1); 

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



// -----------------------------------------------------------------------------
//                        CALC MOBILIZER REACTION FORCES
// -----------------------------------------------------------------------------
void SimbodyMatterSubsystemRep::calcMobilizerReactionForces
   (const State& s, Vector_<SpatialVec>& forces) const 
{
    forces.resize(getNumBodies());
    
    // Find the total body force on every body from all sources *other* than 
    // mobilizer reaction forces.
    
    Vector_<SpatialVec> otherForces = getMultibodySystem().getRigidBodyForces(s, Stage::Dynamics);
    Vector_<SpatialVec> constrainedBodyForces(getNumBodies());
    Vector constrainedMobilizerForces(s.getNU());
    calcConstraintForcesFromMultipliers(s, s.getMultipliers(), constrainedBodyForces, constrainedMobilizerForces);
    otherForces -= constrainedBodyForces;
    for (MobilizedBodyIndex index(0); index < getNumBodies(); index++)
        otherForces[index] -= getGyroscopicForce(s, index);
    
    // Find the total force that was actually applied.
    
    Vector_<SpatialVec> totalForce(forces.size());
    for (MobilizedBodyIndex index(0); index < getNumBodies(); index++) {
        const MobilizedBody& body = getMobilizedBody(index);
        const MassProperties& mass = body.getBodyMassProperties(s);
        const SpatialVec& acceleration = body.getBodyAcceleration(s);
        if (mass.getMass() == Infinity)
            totalForce[index] = SpatialVec(Vec3(0), Vec3(0));
        else
            totalForce[index] = body.calcBodySpatialInertiaMatrixInGround(s)
                                        * acceleration;
    }

    // Starting from the leaf nodes and working back toward ground, take the 
    // difference to find the reaction forces, then apply them to the parent.
    
    for (int i = (int) rbNodeLevels.size()-1; i >= 0; i--)
        for (int j = 0; j < (int) rbNodeLevels[i].size(); j++) {
            MobilizedBodyIndex index = rbNodeLevels[i][j]->getNodeNum();
            forces[index] = totalForce[index]-otherForces[index];
            if (i > 0) {
                const MobilizedBody& body = getMobilizedBody(index);
                MobilizedBodyIndex parentIndex = rbNodeLevels[i][j]->getParent()->getNodeNum();
                const MobilizedBody& parent = getMobilizedBody(parentIndex);
                Vec3 parentPos = parent.findStationAtAnotherBodyStation(s, body, body.getOutboardFrame(s).p());
                parent.applyForceToBodyPoint(s, parentPos, -forces[index][1], otherForces);
                Vec3 offset = parent.getBodyTransform(s).R()*(body.getMobilizerTransform(s).R()*body.getOutboardFrame(s).p());
                otherForces[parentIndex][0] -= forces[index][0]-offset%forces[index][1];
            }
        }
    
    // Transform the force to be reported at the outboard joint location.
    
    for (MobilizedBodyIndex index(0); index < getNumBodies(); index++) {
        const MobilizedBody& body = getMobilizedBody(index);
        Vec3 localForce = ~body.getBodyTransform(s).R()*forces[index][1];
        forces[index][0] -= body.getBodyTransform(s).R()*(body.getOutboardFrame(s).p()%localForce);
    }
}
//....................... CALC MOBILIZER REACTION FORCES .......................



// -----------------------------------------------------------------------------
//                                CALC QDOT
// -----------------------------------------------------------------------------
// Must be done with Position stage to calculate qdot = N*u.
void SimbodyMatterSubsystemRep::calcQDot
   (const State& s, const Vector& u, Vector& qdot) const 
{
    SBStateDigest sbs(s, *this, Stage::Position.next());

    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(sbs, u, qdot);
}
//............................... CALC QDOT ....................................



// -----------------------------------------------------------------------------
//                              CALC QDOTDOT
// -----------------------------------------------------------------------------
// Must be done with Velocity stage to calculate qdotdot = Ndot*u + N*udot.
void SimbodyMatterSubsystemRep::calcQDotDot
   (const State& s, const Vector& udot, Vector& qdotdot) const 
{
    SBStateDigest sbs(s, *this, Stage::Velocity.next());

    assert(udot.size() == getTotalDOF());
    qdotdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDotDot(sbs, udot, qdotdot);
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

    n.calcLocalQDotFromLocalU(sbState, u, qdot);
}

// State must be realized to Stage::Velocity, so that we can extract N(q),
// NDot(q,u), and u from it to calculate qdotdot=N(q)*udot + NDot(q,u)*u for 
// this mobilizer.
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
    const SBTreePositionCache& tpc = getTreePositionCache(s);

    assert(v.size() == getTotalDOF());

    Jv.resize(getNumBodies());

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcSpatialKinematicsFromInternal(tpc,v, Jv);
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
    assert(X.size() == getNumBodies());

    const SBTreePositionCache& tpc = getTreePositionCache(s);

    Vector_<SpatialVec> zTemp(getNumBodies()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(tpc, zTemp, X, JX);
        }
}

// This routine does the same thing as the above but accounts for centrifugal
// forces induced by velocities. The equivalent joint forces returned include
// both the applied forces and the centrifugal ones. Constraints are ignored.
void SimbodyMatterSubsystemRep::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const
{
    const SBTreePositionCache&  tpc = getTreePositionCache(s);
    const SBDynamicsCache&      dc  = getDynamicsCache(s);

    assert(bodyForces.size() == getNumBodies());
    mobilityForces.resize(getTotalDOF());

    Vector_<SpatialVec> allZ(getNumBodies());

    // Don't do ground's level since ground has no inboard joint.
    for (int i=rbNodeLevels.size()-1 ; i>0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcEquivalentJointForces(tpc,dc,
                bodyForces, allZ,
                mobilityForces);
        }
}

bool SimbodyMatterSubsystemRep::getShowDefaultGeometry() const {
    return showDefaultGeometry;
}

void SimbodyMatterSubsystemRep::setShowDefaultGeometry(bool show) {
    showDefaultGeometry = show;
}

std::ostream& operator<<(std::ostream& o, const SimbodyMatterSubsystemRep& tree) {
    o << "SimbodyMatterSubsystemRep has " << tree.getNumBodies() << " bodies (incl. G) in "
      << tree.rbNodeLevels.size() << " levels." << std::endl;
    o << "NodeNum->level,offset;stored nodeNum,level (stateOffset:dim)" << std::endl;
    for (int i=0; i < tree.getNumBodies(); ++i) {
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
        if (mc->instanceCacheIndex.isValid())
            ic = &matter.updInstanceCache(state);
    }
    if (g >= Stage::Time) {
        if (mc->timeVarsIndex.isValid())
            tv = &matter.getTimeVars(state);
        if (mc->timeCacheIndex.isValid())
            tc = &matter.updTimeCache(state);
    }
    if (g >= Stage::Position) {
        if (mc->qVarsIndex.isValid())
            pv = &matter.getPositionVars(state);

        qErr = &matter.updQErr(state);
        if (mc->treePositionCacheIndex.isValid())
            tpc = &matter.updTreePositionCache(state);
        if (mc->constrainedPositionCacheIndex.isValid())
            cpc = &matter.updConstrainedPositionCache(state);
    }
    if (g >= Stage::Velocity) {
        if (mc->uVarsIndex.isValid())
            vv = &matter.getVelocityVars(state);

        qdot = &matter.updQDot(state);
        uErr = &matter.updUErr(state);
        if (mc->treeVelocityCacheIndex.isValid())
            tvc = &matter.updTreeVelocityCache(state);
        if (mc->constrainedVelocityCacheIndex.isValid())
            cvc = &matter.updConstrainedVelocityCache(state);
    }
    if (g >= Stage::Dynamics) {
        if (mc->dynamicsVarsIndex.isValid())
            dv = &matter.getDynamicsVars(state);
        if (mc->dynamicsCacheIndex.isValid())
            dc = &matter.updDynamicsCache(state);

        // Prescribed accelerations are filled in at Dynamics
        // stage so we may need these now.
        udot = &matter.updUDot(state);
        qdotdot = &matter.updQDotDot(state);
    }
    if (g >= Stage::Acceleration) {
        if (mc->accelerationVarsIndex.isValid())
            av = &matter.getAccelerationVars(state);

        udotErr = &matter.updUDotErr(state);
        if (mc->treeAccelerationCacheIndex.isValid())
            tac = &matter.updTreeAccelerationCache(state);
        if (mc->constrainedAccelerationCacheIndex.isValid())
            cac = &matter.updConstrainedAccelerationCache(state);
    }

    stage = g;
}

