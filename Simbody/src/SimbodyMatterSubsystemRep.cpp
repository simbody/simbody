/* -------------------------------------------------------------------------- *
 *                              SimTK Simbody(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-11 Stanford University and the Authors.        *
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
SimbodyMatterSubsystemRep::getMobilizerCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getMobilizerCentrifugalForces(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCentrifugalForces(const State& s, MobilizedBodyIndex body) const {
  return getRigidBodyNode(body).getTotalCentrifugalForces(getDynamicsCache(s));
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
    SimbodyMatterSubsystemRep* mutableThis = 
        const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!subsystemTopologyHasBeenRealized()) 
        mutableThis->endConstruction(s); // no more bodies after this!

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
            const RigidBodyNode& node  = *rbNodeLevels[i][j];
            SBModelPerMobodInfo& mbInfo = 
                mc.updMobilizedBodyModelInfo(node.getNodeNum());

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
        allocateDiscreteVariable(s, Stage::Dynamics, 
                                 new Value<SBDynamicsVars>(dvars));
    mc.dynamicsCacheIndex = 
        allocateCacheEntry(s, Stage::Dynamics, 
                           new Value<SBDynamicsCache>());

    // No Acceleration variables that I know of. But we can go through the
    // charade here anyway.
    SBAccelerationVars rvars;
    rvars.allocate(topologyCache);
    setDefaultAccelerationValues(mv, rvars);
    mc.accelerationVarsIndex = 
        allocateDiscreteVariable(s, Stage::Acceleration, 
                                 new Value<SBAccelerationVars>(rvars));

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
        const SBModelPerMobodInfo& 
            modelInfo    = mc.getMobilizedBodyModelInfo(mbx);
        SBInstancePerMobodInfo&
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
    markCacheValueRealized(s, mc.constrainedPositionCacheIndex);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizePosition(stateDigest);
    return 0;
}



//==============================================================================
//                      REALIZE COMPOSITE BODY INERTIAS
//==============================================================================
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



//==============================================================================
//                     REALIZE ARTICULATED BODY INERTIAS
//==============================================================================
void SimbodyMatterSubsystemRep::
realizeArticulatedBodyInertias(const State& state) const {
    const CacheEntryIndex abx = 
        getModelCache(state).articulatedBodyInertiaCacheIndex;

    if (isCacheValueRealized(state, abx))
        return; // already realized

    SimTK_STAGECHECK_GE_ALWAYS(getStage(state), Stage::Position, 
        "SimbodyMatterSubsystem::realizeArticulatedBodyInertias()");

    const SBInstanceCache&          ic  = getInstanceCache(state);
    const SBTreePositionCache&      tpc = getTreePositionCache(state);
    SBArticulatedBodyInertiaCache&  abc = updArticulatedBodyInertiaCache(state);

    // tip-to-base sweep
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; --i) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; ++j)
            rbNodeLevels[i][j]->realizeArticulatedBodyInertiasInward(ic,tpc,abc);

    markCacheValueRealized(state, abx);
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
    markCacheValueRealized(s, mc.constrainedVelocityCacheIndex);

    // Constraints
    for (ConstraintIndex cx(0); cx < constraints.size(); ++cx)
        getConstraint(cx).getImpl().realizeVelocity(stateDigest);
    return 0;
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

    // tip-to-base calculation
    realizeArticulatedBodyInertias(s); // (may already have been realized)
    const SBArticulatedBodyInertiaCache& abc = getArticulatedBodyInertiaCache(s);

    SBStateDigest stateDigest(s, *this, Stage::Dynamics);

    // Get the Dynamics-stage cache; it was already allocated at Instance stage.
    SBDynamicsCache& dc = stateDigest.updDynamicsCache();

    // Realize velocity-dependent articulated body quantities needed for 
    // dynamics: base-to-tip.
    for (int i=0; i < (int)rbNodeLevels.size(); ++i)
        for (int j=0; j < (int)rbNodeLevels[i].size(); ++j)
            rbNodeLevels[i][j]->realizeDynamics(abc, stateDigest);

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

    SBStateDigest stateDigest(s, *this, Stage::Acceleration);

    // Get the Acceleration-stage cache entries. They were all allocated by
    // the end of Instance stage.
    Vector&                         udot    = stateDigest.updUDot();
    Vector&                         qdotdot = stateDigest.updQDotDot();
    SBTreeAccelerationCache&        tac     = stateDigest.updTreeAccelerationCache();
    SBConstrainedAccelerationCache& cac     = stateDigest.updConstrainedAccelerationCache();

    // We ask our containing MultibodySystem for a reference to the cached 
    // forces accumulated from all the force subsystems. We use these to 
    // compute accelerations, with all results going into the AccelerationCache.
    const MultibodySystem& mbs = getMultibodySystem(); // owner of this subsystem
    realizeLoopForwardDynamics(s,
        mbs.getMobilityForces(s, Stage::Dynamics),
        mbs.getParticleForces(s, Stage::Dynamics),
        mbs.getRigidBodyForces(s, Stage::Dynamics));

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

// TODO: the weight for u_i should be something like the largest Dv_j/Du_i in 
// the system Jacobian, where v_j is body j's origin speed, with a lower limit
// given by length scale (e.g. 1 unit).
// Then the q_i weight should be obtained via dq = N*du.
int SimbodyMatterSubsystemRep::calcQUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNQ(s));
    weights = 1;
    return 0;
}
int SimbodyMatterSubsystemRep::calcUUnitWeightsImpl(const State& s, Vector& weights) const {
    weights.resize(getNU(s));
    weights = 1;
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

void SimbodyMatterSubsystemRep::
setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void SimbodyMatterSubsystemRep::
setMobilizerIsPrescribed(State& s, MobilizedBodyIndex body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void SimbodyMatterSubsystemRep::
setConstraintIsDisabled(State& s, ConstraintIndex constraint, bool disable) const {
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
// the whole multiplication can be done in O(m) time where m=mp+mv+ma is the 
// total number of active constraint equations.
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
//                         CALC WEIGHTED Pq TRANSPOSE
//==============================================================================
// This is a private method.
// We want constraint matrix Pq, but row scaled by 1/constraint tolerances
// and column scaled by 1/q weights (and we're actually computing the
// transpose). We call the result Pqw (or Pqwt for its transpose).
//
//   Pq = P N^-1
//   Pqw = Tp Pq     Wq^-1 
//       = Tp Pq (N Wu^-1 N^-1) 
//       = Tp P Wu^-1 N^-1
//   Pqwt = ~Pqw = ~N^-1 Wu^-1 ~P Tp  (weights are symmetric)
//       (= ~N^-1 ~Pw)
//
// We calculate one column at a time to avoid any matrix ops. We can do the
// column scaling of ~P by Tp for free, but the row scaling requires nq*mp flops.
// Pqwt must have contiguous-data columns.
void SimbodyMatterSubsystemRep::
calcWeightedPqTranspose( 
        const State&     s,
        const Vector&    Tp,   // 1/perr tols
        const Vector&    ooWu, // 1/u weights
        Matrix&          Pqwt) const // nq X mp
{
    const SBInstanceCache& ic = getInstanceCache(s);

    // Global problem dimensions.
    const int mp = ic.totalNHolonomicConstraintEquationsInUse;
    const int nu = getNU(s);
    const int nq = getNQ(s);

    assert(Tp.size() == mp);
    assert(ooWu.size() == nu);

    Pqwt.resize(nq, mp);
    if (mp==0 || nq==0)
        return;

    assert(Pqwt(0).hasContiguousData());

    Vector Ptcol(nu);
    Vector lambdap(mp, Real(0));

    for (int i=0; i < mp; ++i) {
        lambdap[i] = Tp[i]; // this gives column i of ~P scaled by Tp[i]
        multiplyByPVATranspose(s, true, false, false, lambdap, Ptcol);
        lambdap[i] = 0;
        Ptcol.rowScaleInPlace(ooWu); // now (Wu^-1 ~P Tp)_i
        // Calculate (~Pqw)(i) = ~(N^-1) * (~Pw)(i)
        multiplyByNInv(s, true/*transpose*/,Ptcol,Pqwt(i));
    }
}



//==============================================================================
//                      CALC BIAS FOR MULTIPLY BY PVA
//==============================================================================
// We have these constraint equations available:
// (1)  pverr = Pq*qdot - Pt            (Pq==P*N^-1, Pt==c(t,q))
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
//                              MULTIPLY BY Pq
//==============================================================================
// We have these mp constraint equations available:
// (1)  pverr(t,q;qdot) = Pq*qdot - Pt     (where Pq=P*N^-1, Pt=c(t,q))      
//                      = P *u    - Pt
// with P=P(t,q), N=N(q), and Pt=Pt(t,q). Individual constraint
// error equations are calculated in constant time, so the whole set can be 
// evaluated in O(nq+mp) time where mp is the total number of holonomic
// constraint equations. We expect to be given bias_p=-Pt as a 
// precalculated argument. (See calcBiasForMultiplyByPVA().)
//
// Given a q-like vector we can calculate
//      PqXqlike = Pq*qlike = P*N^-1*qlike = pverr(t,q;qlike) - bias_p.
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
// (1)  pverr = Pq*qdot - Pt       (Pq==P*N^-1, Pt==c(t,q))
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
// Calculate multipliers lambda as
//     (G M^-1 ~G) lambda = aerr
// Optimally, we would calculate this mXm matrix in O(m^2) time. I don't know 
// how to calculate it that fast, but using m calls to operator sequence:
//     Gt_j       = Gt* lambda_j        O(n)
//     MInvGt_j   = M^-1* Gt_j          O(n)
//     GMInvGt(j) = G* MInvGt_j         O(n)
// we can calculate it in O(mn) time. As long as m << n, and
// especially if m is a small constant independent of n, and even better
// if we've partitioned it into little subblocks, this is all very 
// reasonable. One slip up and you'll toss in a factor of mn^2 or m^2n and
// screw this up -- be careful!
//
// Note that we do not require contiguous storage for GMInvGt's columns, 
// although we'll take advantage of it if they are. If not, we'll work in
// a contiguous temp and then copy back. This is because we want to allow
// any matrix at the Simbody API level and we don't want to force the API
// method to have to allocate a whole new mXm matrix when all we need for
// a temporary here is an m-length temporary.
//
// Complexity is O(m^2 + m*n) = O(m*n).
//
// TODO: as long as the force transmission matrix for all constraints is G^T
// the resulting matrix is symmetric. But (a) I don't know how to take 
// advantage of that in forming the matrix, and (b) some constraints may
// result in the force transmission matrix != G (this occurs for example for
// some kinds of "working" constraints like sliding friction).
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

    // This dummy is needed for calcMInverseF().
    Vector_<SpatialVec> A_GB(getNumBodies());

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
        calcMInverseF(s, Gtcol, A_GB, MInvGtcol);
        if (columnsAreContiguous)
            multiplyByPVA(s, true, true, true, bias, MInvGtcol, GMInvGt(j));
        else {
            multiplyByPVA(s, true, true, true, bias, MInvGtcol, GMInvGt_j);
            GMInvGt(j) = GMInvGt_j;
        }
    }
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
//                                PRESCRIBE
// =============================================================================
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
        if (npu==0 && nzu==0) return false; // don't invalidate velocities

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
//............................... PRESCRIBE ....................................



//==============================================================================
//                       ENFORCE POSITION CONSTRAINTS
//==============================================================================
// A note on state variable weights:
// - q and u weights are not independent
// - we consider u weights Wu primary and want the weighted variables to be 
//   related like the unweighted ones: qdot=N*u so qdotw=N*uw, where
//   uw = Wu*u. So qdotw should be Wq*qdot=N*uw=N*Wu*u=N*Wu*N^-1*qdot ==>
//   Wq = N*Wu*N^-1. (and Wq^-1=N*Wu^-1*N^-1)
// - Wu is diagonal, but Wq is block diagonal. We have fast operators for
//   multiplying these matrices by columns, but not for producing Wq so we 
//   just create it operationally as we go.

// These statics are for debugging use only.
static Real calcQErrestWeightedNormU(const SimbodyMatterSubsystemRep& matter, 
    const State& s, const Vector& qErrest, const Vector& uWeights) {
    Vector qhatErrest(uWeights.size());
    matter.multiplyByNInv(s, false, qErrest, qhatErrest); // qhatErrest = N+ qErrest
    qhatErrest.rowScaleInPlace(uWeights);                 // qhatErrest = Wu N+ qErrest
    return qhatErrest.normRMS();
}
static Real calcQErrestWeightedNormQ(const SimbodyMatterSubsystemRep& matter, 
    const State& s, const Vector& qErrest, const Vector& qWeights) {
    Vector Wq_qErrest(qErrest.size());
    Wq_qErrest = qErrest.rowScale(qWeights); // Wq*qErrest
    return Wq_qErrest.normRMS();
}

void SimbodyMatterSubsystemRep::enforcePositionConstraints
   (State& s, Real consAccuracy, const Vector& yWeights,
    const Vector& ooTols, Vector& yErrest, System::ProjectOptions opts) const
{
    assert(getStage(s) >= Stage::Position-1);
    SBStateDigest sbs(s, *this, Stage::Model);

    realizeSubsystemPosition(s);

    // First work only with the holonomic (position) constraints, which appear 
    // first in the QErr array. Don't work on the quaternion constraints in 
    // this first section.
    const int mHolo  = getNumHolonomicConstraintEquationsInUse(s);
    const int mQuats = getNumQuaternionsInUse(s);
    const int nq     = getNQ(s);
    const int nu     = getNU(s);

    // Wq    = N * Wu    * N^-1
    // Wq^-1 = N * Wu^-1 * N^-1
    const VectorView uWeights   = yWeights(nq,nu);
    const Vector     ooUWeights = uWeights.elementwiseInvert(); //TODO: precalc
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
    //       (Tp Pq Wq^-1) dq_WLS  = Tp perr
    //                         dq  = Wq^-1 dq_WLS
    //                          q -= dq
    // until RMS(Tp perr) <= 0.1*accuracy.
    //
    // But Pq=P*N^-1, Wq^-1=N*Wu^-1*N^-1 so Pq Wq^-1=P*Wu^-1*N^-1 so we can
    // rewrite the above:
    //     
    //   (Tp P Wu^-1 N^-1) dq_WLS  = Tp perr
    //                         dq  = N Wu^-1 N^-1 dq_WLS
    //                          q -= dq
    //
    // We define Pqwt = ~Pqw = ~(Tp P Wu^-1 N^-1) = ~N^-1 Wu^-1 ~P Tp
    // because diagonal weights are symmetric.
    //
    // This is a nonlinear least squares problem. Below is a full Newton 
    // iteration since we recalculate the iteration matrix each time around the
    // loop. TODO: a solution could be found using the same iteration matrix, 
    // since we are projecting from (presumably) not too far away. Q1: will it
    // be the same solution? Q2: if not, does it matter?
    Vector scaledPerrs = pErrs.rowScale(ooPTols);
    Real normAchievedTRMS = scaledPerrs.normRMS();

    Real lastChangeMadeWRMS = 0; // size of last change in weighted dq
    int nItsUsed = 0;

    // Set how far past the required tolerance we'll attempt to go. 
    // We only fail if we can't achieve consAccuracy, but while we're
    // solving we'll see if we can get consAccuracyToTryFor.
    const Real consAccuracyToTryFor = 
        std::max(0.1*consAccuracy, SignificantReal);

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is sloppy; should depend on constraint tolerance
    // and rank should be saved and reused in velocity and acceleration
    // constraint methods (or should be calculated elsewhere and passed in).
    const Real conditioningTol = mHolo         
      //* SignificantReal; -- too tight
        * SqrtEps;

    if (normAchievedTRMS > consAccuracyToTryFor) {
        Vector saveQ = getQ(s);
        Matrix Pqwt(nq,mHolo);
        Vector dq_WLS(nq), du(nu), dq(nq); // = Wq^-1 dq_WLS
        FactorQTZ Pqw_qtz;
        Real prevNormAchievedTRMS = normAchievedTRMS; // watch for divergence
        const int MaxIterations  = 20;
        do {
            calcWeightedPqTranspose(s, ooPTols, ooUWeights, Pqwt); // nq X mp

            // This factorization acts like a pseudoinverse.
            Pqw_qtz.factor<Real>(~Pqwt, conditioningTol); 

            //std::cout << "POSITION PROJECTION TOL=" << conditioningTol
            //          << " RANK=" << Pq_qtz.getRank() 
            //          << " RCOND=" << Pq_qtz.getRCondEstimate() << std::endl;

            Pqw_qtz.solve(scaledPerrs, dq_WLS); // this is weighted dq_WLS=Wq*dq
            lastChangeMadeWRMS = dq_WLS.normRMS(); // change in weighted norm

            // switch back to unweighted dq=Wq^-1*dq_WLS
            // = N * Wu^-1 * N^-1 * dq_WLS
            multiplyByNInv(s,false,dq_WLS,du);
            du.rowScaleInPlace(ooUWeights); // in place to save memory
            multiplyByN(s,false,du,dq);

            // This causes quaternions to become unnormalized, but it doesn't
            // matter because N is calculated from the unnormalized q's so
            // scales dq to match.
            updQ(s) -= dq; // this is unweighted dq
            anyChange = true;

            // Now recalculate the position constraint errors at the new q.
            realizeSubsystemPosition(s); // pErrs changes here

            scaledPerrs = pErrs.rowScale(ooPTols); // Tp * pErrs
            normAchievedTRMS = scaledPerrs.normRMS();
            ++nItsUsed;

            if (localOnly && nItsUsed >= 2 
                && normAchievedTRMS > prevNormAchievedTRMS) {
                // perr norm got worse; restore to end of previous iteration
                updQ(s) += dq;
                realizeSubsystemPosition(s); // pErrs changes here
                scaledPerrs = pErrs.rowScale(ooPTols);
                normAchievedTRMS = scaledPerrs.normRMS();
                break; // diverging -- quit now to prevent a bad solution
            }

            prevNormAchievedTRMS = normAchievedTRMS;

        } while (normAchievedTRMS > consAccuracyToTryFor
                 && nItsUsed < MaxIterations);

        // Make sure we achieved at least the required constraint accuracy.
        if (normAchievedTRMS > consAccuracy) {
            updQ(s) = saveQ; // revert
            realizeSubsystemPosition(s);
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "Failed to converge in position projection");
        }


        // Next, if we projected out the position constraint errors, remove the
        // corresponding error from the integrator's error estimate.
        //
        //   (Tp Pq Wq^-1) dq_WLS  = (Tp Pq Wq^-1) Wq*qErrest
        //                     dq  = Wq^-1 dq_WLS
        //                         = N Wu^-1 N^-1 dq_WLS
        //                qErrest -= dq
        // No iteration is required.
        //
        // We can simplify the RHS of the first equation above:
        //        (Tp Pq Wq^-1) Wq qErrest = Tp Pq qErrest
        // for which we have an O(n) operator to use for the matrix-vector product.
        if (qErrest.size()) {
            // Work in Wq-norm
            Vector Tp_Pq_qErrest, bias_p;
            calcBiasForMultiplyByPq(s, bias_p);
            multiplyByPq(s, bias_p, qErrest, Tp_Pq_qErrest); // Pq*qErrest
            Tp_Pq_qErrest.rowScaleInPlace(ooPTols); // now Tp*Pq*qErrest

            Pqw_qtz.solve(Tp_Pq_qErrest, dq_WLS); // weighted
            const Real normOfAdjustment_WRMS = dq_WLS.normRMS();
            // Switch back to unweighted dq=Wq^-1*dq_WLS
            // = N * Wu^-1 * N^-1 * dq_WLS
            multiplyByNInv(s,false,dq_WLS,du);
            du.rowScaleInPlace(ooUWeights); // in place to save memory
            multiplyByN(s,false,du,dq);
            qErrest -= dq; // unweighted
        }
    }

    //cout << "!!!! perr TRMS achieved " << normAchievedTRMS << " in " 
    //     << nItsUsed << " iterations"  << endl;

    // By design, normalization of quaternions can't have any effect on the 
    // length constraints we just fixed (because we normalize internally for 
    // calculations). So now we can simply normalize the quaternions.
    if (mQuats) {
        Vector& q  = updQ(s); // invalidates q's. TODO: see below.

        for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
            for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
                if (rbNodeLevels[i][j]->enforceQuaternionConstraints(sbs,q,qErrest))
                    anyChange = true;

        // This will recalculate the qnorms (all 1), qerrs (all 0). The only
        // other quaternion dependency is the N matrix (and NInv, NDot).
        // TODO: better if these updates could be made without invalidating
        // Position stage in general. I *think* N is always calculated on they fly.
        realizeSubsystemPosition(s);
    }
}
//........................ ENFORCE POSITION CONSTRAINTS ........................



//==============================================================================
//                         ENFORCE VELOCITY CONSTRAINTS
//==============================================================================
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
    //   (Tpv [P;V] Wu^-1) du_WLS  = Tpv uerr
    //                         du  = Wu^-1*du_WLS
    //                          u -= du
    // Note that although this is a nonlinear least squares problem since uerr 
    // is a function of u, we do not need to refactor the matrix since it does 
    // not depend on u.
    // TODO: I don't think that's true -- V can depend on u (rarely). That
    // doesn't mean we need to refactor it, but then this is a modified Newton
    // iteration (rather than full) if we're not updating V when we could be.
    // TODO: Tp P Wu^-1 should already have been calculated for position 
    // projection (at least if any position projection occurred)
    //
    // This is a nonlinear least squares problem, but we only need to factor 
    // once since only the RHS is dependent on u (TODO: see above).
    Vector scaledVerrs = vErrs.rowScale(ooPVTols);
    Real normAchievedTRMS = scaledVerrs.normRMS();
    
    //cout << "!!!! initially @" << s.getTime() << ", verr TRMS=" 
    //     << normAchievedTRMS << " consAcc=" << consAccuracy;
    //if (uErrest.size())
    //    cout << " uErrest WRMS=" << uErrest.rowScale(uWeights).normRMS();
    //else cout << " NO U ERROR ESTIMATE";
    //cout << endl;
    

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

    if (normAchievedTRMS > consAccuracyToTryFor) {
        const Vector saveU = getU(s);
        Matrix PVt(nu, mHolo+mNonholo);

        calcPVATranspose(s, true, true, false, PVt); // just P,V
        PVt.rowAndColScaleInPlace(ooUWeights, ooPVTols); 
        // PVt is now Wu^-1 (Pt Vt) Tpv

        // Calculate pseudoinverse (just once)
        FactorQTZ PVqtz;
        PVqtz.factor<Real>(~PVt, conditioningTol);

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

            if (localOnly && nItsUsed >= 2 
                && normAchievedTRMS > prevNormAchievedTRMS) {
                // Velocity norm got worse -- restore to end of previous iteration.
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
            realizeSubsystemVelocity(s);
            SimTK_THROW1(Exception::NewtonRaphsonFailure, 
                         "Failed to converge in velocity projection");
        }

        // Next, if we projected out the velocity constraint errors, remove the
        // corresponding error from the integrator's error estimate.
        //   (Tpv [P;V] Wu^-1) du_WLS  = (Tpv [P;V] Wu^-1) Wu uErrest
        //                         du  = Wu^-1 du_WLS
        //                    uErrest -= du
        // No iteration is required.
        // We can simplify the RHS of the first equation above:
        //   (Tpv [P;V] Wu^-1) Wu uErrest = Tpv [P;V] uErrest
        // for which we have an O(n) operator to compute the matrix-vector product.

        if (uErrest.size()) {
            Vector PV_uErrest(mHolo+mNonholo);
            Vector bias_pv(mHolo+mNonholo);
            calcBiasForMultiplyByPVA(s,true,true,false,bias_pv); // just P,V
            multiplyByPVA(s,true,true,false,bias_pv,uErrest,PV_uErrest);
            PV_uErrest.rowScaleInPlace(ooPVTols); // PV_uErrest = Tpv [P;V] uErrEst
            PVqtz.solve(PV_uErrest, du_WLS);
            du_WLS.rowScaleInPlace(ooUWeights); // now du (=Wu^-1*du_WLS)
            uErrest -= du_WLS; // this is really unweighted du
        }
    }
   
    //cout << "!!!! verr achieved " << normAchievedTRMS << " in " 
    //     << nItsUsed << " iterations" << endl;
    //if (uErrest.size())
    //    cout << " uErrest WRMS=" << uErrest.rowScale(uWeights).normRMS() << endl;
}
//........................ ENFORCE VELOCITY CONSTRAINTS ........................



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
    Vector_<SpatialVec>& A_GB           = tac.bodyAccelerationInGround;
    Vector&              tau            = tac.presMotionForces;

    // Calculate accelerations produced by these forces in three forms:
    // body accelerations A_GB, u-space generalized acceleratiosn udot,
    // and q-space generalized accelerations qdotdot.
    calcTreeAccelerations(s, *mobilityForcesToUse, *bodyForcesToUse,
                          netHingeForces, abForcesZ, A_GB, udot, qdotdot, tau);

    // Feed the accelerations into the constraint error methods to determine
    // the acceleratin constraint errors they generate.
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
    const SBModelCache& mc  = getModelCache(s);

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
    markCacheValueRealized(s, mc.treeAccelerationCacheIndex);
}
//....................... REALIZE TREE FORWARD DYNAMICS ........................



//==============================================================================
//                    CALC LOOP FORWARD DYNAMICS OPERATOR
//==============================================================================
// Given a State realized through Stage::Dynamics, and a complete set of applied 
// forces, calculate all acceleration results resulting from those forces AND 
// enforcement of the acceleration constraints. The results go into the return 
// arguments here. This routine *does not* affect the State cache -- it is an 
// operator. In typical usage, the output arguments actually will be part of 
// the state cache to effect a response, but this method can also be used to 
// effect an operator.
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
        0, 0, tac, udot, qdotdot, udotErr);

    // Next, determine how many acceleration-level constraint equations 
    // need to be obeyed.

    const int mHolo    = getNumHolonomicConstraintEquationsInUse(s);
    const int mNonholo = getNumNonholonomicConstraintEquationsInUse(s);
    const int mAccOnly = getNumAccelerationOnlyConstraintEquationsInUse(s);
    const int m        = mHolo+mNonholo+mAccOnly;
    const int nq       = getNQ(s);
    const int nu       = getNU(s);

    multipliers.resize(m);
    if (m==0) return;
    if (nu==0) {multipliers.setToZero(); return;}

    // Conditioning tolerance. This determines when we'll drop a 
    // constraint. 
    // TODO: this is probably too tight; should depend on constraint tolerance
    // and should be consistent with position and velocity projection ranks.
    // Tricky here because conditioning depends on mass matrix as well as
    // constraints.
    const Real conditioningTol = m 
        //* SignificantReal;
        * SqrtEps*std::sqrt(SqrtEps); // Eps^(3/4)

    // Calculate multipliers lambda as
    //     (G M^-1 ~G) lambda = aerr
    // The method here calculates the mXm matrix G*M^-1*G^T as fast as 
    // I know how to do, O(m*n) with O(n) temporary memory, using a series
    // of O(n) operators. Then we'll factor it here in O(m^3) time. 
    Matrix GMInvGt(m,m);
    calcGMInvGt(s, GMInvGt);
    
    // specify 1/cond at which we declare rank deficiency
    FactorQTZ qtz(GMInvGt, conditioningTol); 
    qtz.solve(udotErr, multipliers);

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
    const SBModelCache& mc  = getModelCache(s);

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
    markCacheValueRealized(s, mc.treeAccelerationCacheIndex);
    markCacheValueRealized(s, mc.constrainedAccelerationCacheIndex);
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
    const Real* mobilityForcePtr = mobilityForces.size() ? &mobilityForces[0] : NULL;
    const SpatialVec* bodyForcePtr = bodyForces.size() ? &bodyForces[0] : NULL;

    // TODO: does this need to do level 0?
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.realizeZ(tpc,abc,tvc,dc,tac,mobilityForcePtr,bodyForcePtr);
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
    Real* udot    = updUDot(s).size() ? &updUDot(s)[0] : NULL;
    Real* qdotdot = updQDotDot(s).size() ? &updQDotDot(s)[0] : NULL;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.realizeAccel(tpc,abc,tvc,dc,tac,&udot[node.getUIndex()]);
            node.calcQDotDot(sbs, &udot[node.getUIndex()], &qdotdot[node.getQIndex()]);
        }
}



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
// This Subsystem must have already been realized to Dynamics stage so that 
// dynamics quantities like articulated body inertias are available.
// All vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcTreeAccelerations(const State& s,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Array_<SpatialVec,MobilizedBodyIndex>& abForcesZ, 
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot,    // in/out (in for prescribed udots)
    Vector&                    qdotdot,
    Vector&                    tau) const 
{
    SBStateDigest sbs(s, *this, Stage::Dynamics.next());
    const SBArticulatedBodyInertiaCache& abc = getArticulatedBodyInertiaCache(s);

    const SBInstanceCache&      ic  = sbs.getInstanceCache();
    const SBTreePositionCache&  tpc = sbs.getTreePositionCache();
    const SBTreeVelocityCache&  tvc = sbs.getTreeVelocityCache();
    const SBDynamicsCache&      dc  = sbs.getDynamicsCache();

    assert(mobilityForces.size() == getTotalDOF());
    assert(bodyForces.size() == getNumBodies());

    netHingeForces.resize(getTotalDOF());
    abForcesZ.resize(getNumBodies());
    A_GB.resize(getNumBodies());
    udot.resize(getTotalDOF());
    qdotdot.resize(getTotalQAlloc());
    tau.resize(ic.totalNPresForce);

    assert(mobilityForces.hasContiguousData());
    assert(bodyForces.hasContiguousData());
    assert(netHingeForces.hasContiguousData());
    assert(A_GB.hasContiguousData());
    assert(udot.hasContiguousData());
    assert(qdotdot.hasContiguousData());
    assert(tau.hasContiguousData());

    // Temporary
    Vector_<SpatialVec> allGepsilon(getNumBodies());

    const Real*       mobilityForcePtr = mobilityForces.size() 
                                            ? &mobilityForces[0] : NULL;
    const SpatialVec* bodyForcePtr     = bodyForces.size() 
                                            ? &bodyForces[0] : NULL;
    Real*             hingeForcePtr    = netHingeForces.size() 
                                            ? &netHingeForces[0] : NULL;
    SpatialVec*       aPtr             = A_GB.size()    ? &A_GB[0] : NULL;
    Real*             udotPtr          = udot.size()    ? &udot[0] : NULL;
    Real*             qdotdotPtr       = qdotdot.size() ? &qdotdot[0] : NULL;
    Real*             tauPtr           = tau.size()     ? &tau[0] : NULL;
    SpatialVec*       zPtr             = abForcesZ.begin();    
    SpatialVec*       gepsPtr          = allGepsilon.size() 
                                            ? &allGepsilon[0] : NULL;    

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass1Inward(ic,tpc,abc,dc,
                mobilityForcePtr, bodyForcePtr, udotPtr, zPtr, gepsPtr,
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
//                            CALC M INVERSE F
//==============================================================================
// Calculate udot = M^-1 f. We also get spatial accelerations A_GB for 
// each body as a side effect.
// This Subsystem must already be realized through Dynamics stage.
// All vectors must use contiguous storage.
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

    assert(f.hasContiguousData());
    assert(A_GB.hasContiguousData());
    assert(udot.hasContiguousData());

    // Temporaries
    Vector              allEpsilon(getTotalDOF());
    Vector_<SpatialVec> allZ(getNumBodies());
    Vector_<SpatialVec> allGepsilon(getNumBodies());
    const Real* fPtr = f.size() ? &f[0] : NULL;
    SpatialVec* aPtr = A_GB.size() ? &A_GB[0] : NULL;
    Real* udotPtr = udot.size() ? &udot[0] : NULL;
    SpatialVec* zPtr = allZ.size() ? &allZ[0] : NULL;
    SpatialVec* gepsPtr = allGepsilon.size() ? &allGepsilon[0] : NULL;
    Real* epsPtr = allEpsilon.size() ? &allEpsilon[0] : NULL;

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass1Inward(ic,tpc,abc,dc,
                fPtr, zPtr, gepsPtr, epsPtr);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass2Outward(ic,tpc,abc,dc, 
                epsPtr, aPtr, udotPtr);
        }
}
//............................. CALC M INVERSE F ...............................



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
//                                 CALC M V
//==============================================================================
// Calculate x = M v. If the vector v is a generalized acceleration
// udot, then we also get spatial accelerations A_GB for 
// each body as a side effect.
// This Subsystem must already have been realized to Position stage.
// All vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcMV(const State& s,
    const Vector&              v,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    f) const 
{
    const SBTreePositionCache& tpc = getTreePositionCache(s);

    assert(v.size() == getTotalDOF());

    A_GB.resize(getNumBodies());
    f.resize(getTotalDOF());

    assert(v.hasContiguousData());
    assert(A_GB.hasContiguousData());
    assert(f.hasContiguousData());

    const Real* vPtr = v.size() ? &v[0] : NULL;
    SpatialVec* aPtr = A_GB.size() ? &A_GB[0] : NULL;
    Real* fPtr = f.size() ? &f[0] : NULL;

    // Temporary
    Vector_<SpatialVec> allFTmp(getNumBodies());
    SpatialVec* tmpPtr = allFTmp.size() ? &allFTmp[0] : NULL;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMVPass1Outward(tpc, vPtr, aPtr);
        }

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMVPass2Inward(tpc,aPtr,tmpPtr,fPtr);
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
    // repeated O(n) calls to calcMV() to get M one column at a time.

    // If M's columns are contiguous we can avoid copying.
    const bool isContiguous = M(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : nu);

    Vector_<SpatialVec> A_GB(getNumBodies()); // unused dummy needed
    Vector v(nu); v.setToZero();
    for (int i=0; i < nu; ++i) {
        v[i] = 1;
        if (isContiguous) {
            calcMV(s, v, A_GB, M(i));
        } else {
            calcMV(s, v, A_GB, contig_col);
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
    // the O(n) operator calcMInverseF().

    // If M's columns are contiguous we can avoid copying.
    const bool isContiguous = MInv(0).hasContiguousData();
    Vector contig_col(isContiguous ? 0 : nu);

    Vector_<SpatialVec> A_GB(getNumBodies()); // unused dummy needed
    Vector f(nu); f = 0;
    for (int i=0; i < nu; ++i) {
        f[i] = 1;
        if (isContiguous) {
            calcMInverseF(s, f, A_GB, MInv(i));
        } else {
            calcMInverseF(s, f, A_GB, contig_col);
            MInv(i) = contig_col;
        }
        f[i] = 0;
    }
}



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
//     reaction = P*(A-a) + z
// where P is an articulated body inertia, A is the spatial acceleration of
// that body, a is the Coriolis acceleration, and z is the articulated body
// force. All of these quantities are already available at Stage::Acceleration.
// See Abhi Jain's 2011 book "Robot and Multibody Dynamics", Eq. 7.34 on
// page 128.
//
// After calculating the reaction at the body frame origin Bo, we shift it to
// the mobilizer's outboard frame M and report it there, though expressed in G.
// Note that any generalized forces applied at mobilities end up included in
// the reaction forces.
//
// Cost is 105 flops/body plus lots of memory access to dredge up the 
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

    const Array_<ArticulatedInertia,MobilizedBodyIndex>& P = 
                                            getArticulatedBodyInertias(s);
    const Array_<SpatialVec,MobilizedBodyIndex>& z =
                                            getArticulatedBodyForces(s);

    for (MobodIndex mbx(0); mbx < nb; ++mbx) {
        const MobilizedBody& body   = getMobilizedBody(mbx);
        const SpatialVec& A_GB = body.getBodyAcceleration(s);
        const SpatialVec& a    = getMobilizerCoriolisAcceleration(s,mbx);
        SpatialVec FB_G = z[mbx];
        if (mbx != GroundIndex) FB_G += P[mbx]*(A_GB-a); // 78 flops
        // Shift to M
        const Transform& X_GB   = body.getBodyTransform(s);
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



// =============================================================================
//                                CALC QDOT
// =============================================================================
// Must be done with Position stage to calculate qdot = N*u.
// Both vectors must use contiguous storage.
void SimbodyMatterSubsystemRep::calcQDot
   (const State& s, const Vector& u, Vector& qdot) const 
{
    SBStateDigest sbs(s, *this, Stage::Position.next());

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
    SBStateDigest sbs(s, *this, Stage::Velocity.next());

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
    const SBDynamicsCache&      dc  = getDynamicsCache(s);

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
            node.calcEquivalentJointForces(tpc,dc,
                bodyForcePtr, zPtr,
                mobilityForcePtr);
        }
}
//.................... CALC TREE EQUIVALENT MOBILITY FORCES ....................



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
        
        // All cache entries, for any stage, can be modified at instance stage 
        // or later.
        
        if (mc->timeCacheIndex.isValid())
            tc = &matter.updTimeCache(state);
        if (mc->treePositionCacheIndex.isValid())
            tpc = &matter.updTreePositionCache(state);
        if (mc->constrainedPositionCacheIndex.isValid())
            cpc = &matter.updConstrainedPositionCache(state);
        if (mc->treeVelocityCacheIndex.isValid())
            tvc = &matter.updTreeVelocityCache(state);
        if (mc->constrainedVelocityCacheIndex.isValid())
            cvc = &matter.updConstrainedVelocityCache(state);
        if (mc->dynamicsCacheIndex.isValid())
            dc = &matter.updDynamicsCache(state);
        if (mc->treeAccelerationCacheIndex.isValid())
            tac = &matter.updTreeAccelerationCache(state);
        if (mc->constrainedAccelerationCacheIndex.isValid())
            cac = &matter.updConstrainedAccelerationCache(state);
    }
    if (g >= Stage::Time) {
        if (mc->timeVarsIndex.isValid())
            tv = &matter.getTimeVars(state);
    }
    if (g >= Stage::Position) {
        if (mc->qVarsIndex.isValid())
            pv = &matter.getPositionVars(state);

        qErr = &matter.updQErr(state);
    }
    if (g >= Stage::Velocity) {
        if (mc->uVarsIndex.isValid())
            vv = &matter.getVelocityVars(state);

        qdot = &matter.updQDot(state);
        uErr = &matter.updUErr(state);
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

