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
#include "simbody/internal/common.h"

#include "SimbodyMatterSubsystemRep.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "ConstraintNode.h"
#include "LengthConstraints.h"
#include "MultibodySystemRep.h"
#include "MobilizedBodyRep.h"
#include "ConstraintRep.h"

#include <string>

SimbodyMatterSubsystemRep::SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep& src)
  : SimTK::Subsystem::Guts("SimbodyMatterSubsystemRep", "X.X.X")
{
    assert(!"SimbodyMatterSubsystemRep copy constructor ... TODO!");
}


void SimbodyMatterSubsystemRep::clearTopologyState() {
    // Constraints are independent from one another, so any deletion order
    // is fine. However, they depend on bodies and not vice versa so we'll
    // delete them first just to be neat.
    for (int i=0; i < (int)constraints.size(); ++i)
        delete constraints[i];
    constraints.clear();

    // These are the owner handles, so this deletes the MobilizedBodyReps also.
    // We'll delete from the terminal nodes inward just to be neat.
    for (int i=(int)mobilizedBodies.size()-1; i >= 0; --i)
        delete mobilizedBodies[i];
    mobilizedBodies.clear();
}

void SimbodyMatterSubsystemRep::clearTopologyCache() {
    invalidateSubsystemTopologyCache();

    // TODO: state indices really shouldn't be dealt out until Stage::Model.
    // At the moment they are part of the topology.
    nextUSlot=nextUSqSlot=nextQSlot        = 0;
    nextQErrSlot=nextUErrSlot=nextMultSlot = 0;
    DOFTotal=SqDOFTotal=maxNQTotal         = -1;
    topologyCache.clear();
    topologyCacheIndex = -1;

    delete lConstraints; lConstraints=0;

    for (int i=0; i<(int)distanceConstraints.size(); ++i)
        delete distanceConstraints[i];
    distanceConstraints.clear();

    for (int i=0; i<(int)pointInPlaneConstraints.size(); ++i)
        delete pointInPlaneConstraints[i];
    pointInPlaneConstraints.clear();

    // RigidBodyNodes themselves are owned by the MobilizedBodyReps and will
    // be deleted when the MobilizedBodyRep objects are.
    rbNodeLevels.clear();
    nodeNum2NodeMap.clear();
}

MobilizedBodyId SimbodyMatterSubsystemRep::adoptMobilizedBody
   (MobilizedBodyId parentId, MobilizedBody& child) 
{
    invalidateSubsystemTopologyCache();

    const MobilizedBodyId id((int)mobilizedBodies.size());
    assert(parentId < id);

    mobilizedBodies.push_back(new MobilizedBody()); // grow
    MobilizedBody& m = *mobilizedBodies.back(); // refer to the empty handle we just created

    child.disown(m); // transfer ownership to m

    // Now tell the MobilizedBody object its owning MatterSubsystem, id within
    // that Subsystem, and parent MobilizedBody object.
    m.updRep().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), parentId, id);
    return id;
}

ConstraintId SimbodyMatterSubsystemRep::adoptConstraint(Constraint& child) {
    invalidateSubsystemTopologyCache();

    const ConstraintId id((int)constraints.size());

    constraints.push_back(new Constraint()); // grow
    Constraint& c = *constraints.back(); // refer to the empty handle we just created

    child.disown(c); // transfer ownership to c

    // Now tell the Constraint object its owning MatterSubsystem and id within
    // that Subsystem.
    c.updRep().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), id);
    return id;
}


void SimbodyMatterSubsystemRep::createGroundBody() {
    assert(mobilizedBodies.empty());
    invalidateSubsystemTopologyCache();

    mobilizedBodies.push_back(new MobilizedBody::Ground());
    mobilizedBodies[0]->updRep().setMyMatterSubsystem(updMySimbodyMatterSubsystemHandle(), 
                                                      MobilizedBodyId(), 
                                                      MobilizedBodyId(0));
}

// Add a distance constraint and assign it to use a particular slot in the
// qErr, uErr, and multiplier arrays.
// Return the assigned distance constraint index for caller's use.
int SimbodyMatterSubsystemRep::addOneDistanceConstraintEquation(
    const RBStation& s1, const RBStation& s2, const Real& d,
    int qerrIndex, int uerrIndex, int multIndex)
{
    RBDistanceConstraint* dc = new RBDistanceConstraint(s1,s2,d);
    dc->setQErrIndex(qerrIndex);
    dc->setUErrIndex(uerrIndex);
    dc->setMultIndex(multIndex);
    dc->setDistanceConstraintNum(distanceConstraints.size());
    distanceConstraints.push_back(dc);
    return distanceConstraints.size()-1;
}


// Add a point-in-plane constraint and assign it to use a particular slot in the
// qErr, uErr, and multiplier arrays.
// Return the assigned point-in-plane constraint index for caller's use.
int SimbodyMatterSubsystemRep::addOnePointInPlaneEquation(
        const RBDirection& d, Real height, const RBStation& s,
        int qerrIndex, int uerrIndex, int multIndex)
{
    RBPointInPlaneConstraint* pipc = new RBPointInPlaneConstraint(d,height,s);
    pipc->setQErrIndex(qerrIndex);
    pipc->setUErrIndex(uerrIndex);
    pipc->setMultIndex(multIndex);
    pipc->setPointInPlaneConstraintNum(pointInPlaneConstraints.size());
    pointInPlaneConstraints.push_back(pipc);
    return pointInPlaneConstraints.size()-1;
}

MobilizedBodyId SimbodyMatterSubsystemRep::getParent(MobilizedBodyId body) const { 
    return MobilizedBodyId(getRigidBodyNode(body).getParent()->getNodeNum());
}

Array<MobilizedBodyId> SimbodyMatterSubsystemRep::getChildren(MobilizedBodyId body) const {
    const RigidBodyNode& node = getRigidBodyNode(body);
    Array<MobilizedBodyId> children;
    for (MobilizedBodyId i(0); i < node.getNChildren(); ++i)
        children += MobilizedBodyId(node.getChild(i)->getNodeNum());
    return children;
}

const MassProperties&
SimbodyMatterSubsystemRep::getDefaultBodyMassProperties(MobilizedBodyId body) const
  { return getRigidBodyNode(body).getMassProperties_OB_B(); }

const Transform&
SimbodyMatterSubsystemRep::getDefaultMobilizerFrame(MobilizedBodyId body) const
  { return getRigidBodyNode(body).getX_BM(); }

const Transform&
SimbodyMatterSubsystemRep::getDefaultMobilizerFrameOnParent(MobilizedBodyId body) const
  { return getRigidBodyNode(body).getX_PF(); }

const Transform&
SimbodyMatterSubsystemRep::getBodyTransform(const State& s, MobilizedBodyId body) const
  { return getRigidBodyNode(body).getX_GB(getPositionCache(s)); }

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyVelocity(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getV_GB(getVelocityCache(s));
}

const SpatialVec&
SimbodyMatterSubsystemRep::getCoriolisAcceleration(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCoriolisAcceleration(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getTotalCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getGyroscopicForce(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getGyroscopicForce(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getCentrifugalForces(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getCentrifugalForces(getDynamicsCache(s));
}

const SpatialMat&
SimbodyMatterSubsystemRep::getArticulatedBodyInertia(const State& s, MobilizedBodyId body) const {
  return getRigidBodyNode(body).getArticulatedBodyInertia(getDynamicsCache(s));
}

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes we're going to need later for state
// variables and cache entries. We allocate and initialize all the
// Modeling variables here.
void SimbodyMatterSubsystemRep::endConstruction() {
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
    nextUSlot = nextUSqSlot = nextQSlot = 0; // state allocation

    //Must do these in order from lowest number (ground) to highest. 
    for (int i=0; i<getNumMobilizedBodies(); ++i) {
        // Create the RigidBodyNode properly linked to its parent.
        const MobilizedBody::MobilizedBodyRep& mbr = getMobilizedBody(MobilizedBodyId(i)).getRep();
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
    // This creates a ConstraintNode owned by the Topology cache of each Constraint, and 
    // doles out the multipliers and constraint error slots for topological
    // constraints. 
    // TODO: currently our constraints are all holonomic, meaning
    // position-level, so that they occupy one slot in the qErr array, then
    // their time derivatives need one slot in uErr and their 2nd time 
    // derivatives need one acceleration-level multiplier. Later we will
    // have constraints which start at the velocity level (nonholonomic) so
    // they won't use up a qErr slot.
    // Also, quaternion normalization constraints exist only at the 
    // position level, however they are not topological since modeling
    // choices affect whether we use them. See realizeModel() below.
    nextQErrSlot = nextUErrSlot = nextMultSlot = 0; // state cache allocation
    for (int i=0; i<getNumConstraints(); ++i) {
        const Constraint::ConstraintRep& crep = getConstraint(ConstraintId(i)).getRep();
        crep.realizeTopology(nextQErrSlot,nextUErrSlot,nextMultSlot);
    }
    
    // Now create the computational data structure for the length constraint
    // equations which currently are used to implement all the Constraints.
    // This is owned by the SimbodyMatterSubsystemRep.
    lConstraints = new LengthConstraints(*this, 0);
    lConstraints->construct(distanceConstraints, pointInPlaneConstraints);
}

int SimbodyMatterSubsystemRep::realizeSubsystemTopologyImpl(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "SimbodyMatterSubsystem::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this call*.
    // Afterwards they are truly const so we don't declare them mutable, but cheat
    // here instead.
    SimbodyMatterSubsystemRep* mutableThis = const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!subsystemTopologyHasBeenRealized()) mutableThis->endConstruction(); // no more bodies after this!

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
    mutableThis->topologyCache.nPointInPlaneConstraints = pointInPlaneConstraints.size();

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

    // Get the Modeling-stage cache and make sure it has been allocated and initialized if needed.
    // It is OK to hold a reference here because the discrete variables (and cache entries) in
    // the State are stable, that is, they don't change even if more variables are added.
    SBModelCache& mc = updModelCache(s);
    mc.allocate(topologyCache);

    // Count quaternions, and assign a "quaternion index" to each body that
    // needs one. We can't do this until Model stage because there is a modeling
    // variable which decides whether ball joints get quaternions or Euler angles).
    mc.nQuaternionsInUse = 0;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            if (node.isUsingQuaternion(mv)) {
                mc.quaternionIndex[node.getNodeNum()] = mc.nQuaternionsInUse;
                mc.nQuaternionsInUse++;
            }
        }

    mc.firstQuaternionQErrSlot = nextQErrSlot; // begins after last topological constraint qerr

    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModel(mv,mc); 

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
    for (MobilizedBodyId i(1); i < (int)mobilizedBodies.size(); ++i) {
        const MobilizedBody::MobilizedBodyRep& mb = mobilizedBodies[i]->getRep();
        mb.copyOutDefaultQ(s, qInit);
    }

    mc.qIndex = s.allocateQ(getMySubsystemId(), qInit);
    mc.qVarsIndex = -1; // no config vars other than q
    mc.qCacheIndex = allocateCacheEntry(s, Stage::Position, new Value<SBPositionCache>());

    // We'll store the the physical constraint errors (which consist solely of distance
    // constraint equations at the moment), followed by the quaternion constraints.
    mc.qErrIndex = allocateQErr(s, nextQErrSlot + mc.nQuaternionsInUse);

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.

    Vector uInit(DOFTotal);
    setDefaultVelocityValues(mv, uInit);

    mc.uIndex = s.allocateU(getMySubsystemId(), uInit);
    mc.uVarsIndex = -1; // no velocity vars other than u
    mc.uCacheIndex = allocateCacheEntry(s, Stage::Velocity, new Value<SBVelocityCache>());
    // Note that qdots are automatically allocated in the Velocity stage cache.

    // Only physical constraints exist at the velocity and acceleration levels; 
    // the quaternion normalization constraints are gone.
    mc.uErrIndex    = allocateUErr   (s, nextUErrSlot);
    mc.udotErrIndex = allocateUDotErr(s, nextMultSlot);

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

    // Put position constraint equation errors in qErr
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(qErr,pc);
    for (int i=0; i < (int)pointInPlaneConstraints.size(); ++i)
        pointInPlaneConstraints[i]->calcPosInfo(qErr,pc);

    return 0;
}

// Set generalized speeds: sweep from base to tip.
// realizePosition() must have been called already.
int SimbodyMatterSubsystemRep::realizeSubsystemVelocityImpl(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "SimbodyMatterSubsystem::realizeVelocity()");

    const SBModelVars&     mv   = getModelVars(s);
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

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcVelInfo(pc,uErr,vc);
    for (int i=0; i < (int)pointInPlaneConstraints.size(); ++i)
        pointInPlaneConstraints[i]->calcVelInfo(pc,uErr,vc);

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
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    // Let the bodies and mobilizers have a chance to generate some geometry.
    for (MobilizedBodyId i(0); i<(int)mobilizedBodies.size(); ++i)
        mobilizedBodies[i]->getRep().calcDecorativeGeometryAndAppend(s,stage,geom);

    // Likewise for the constraints
    for (ConstraintId i(0); i<(int)constraints.size(); ++i)
        constraints[i]->getRep().calcDecorativeGeometryAndAppend(s,stage,geom);

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
        //geom.push_back(DecorativeSphere(0.02).setBodyId(GroundId).setTransform(com)
        //                .setColor(Green).setRepresentation(DecorativeGeometry::DrawPoints)
         //               .setResolution(1));
    }
    default: 
        assert(getStage(s) >= stage);
    }

    return 0;
}

int SimbodyMatterSubsystemRep::getQIndex(MobilizedBodyId body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getQIndex();
}
int SimbodyMatterSubsystemRep::getQAlloc(MobilizedBodyId body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getMaxNQ();
}
int SimbodyMatterSubsystemRep::getUIndex(MobilizedBodyId body) const {
    assert(subsystemTopologyHasBeenRealized());
    return getRigidBodyNode(body).getUIndex();
}
int SimbodyMatterSubsystemRep::getDOF(MobilizedBodyId body) const {
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
    modelVars.prescribed = false;
    modelVars.prescribed[0] = true; // ground
    //modelVars.enabled.assign(getNConstraints(), false);
    modelVars.enabled = false;

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
void SimbodyMatterSubsystemRep::setMobilizerIsPrescribed(State& s, MobilizedBodyId body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void SimbodyMatterSubsystemRep::setConstraintIsEnabled(State& s, ConstraintId constraint, bool enable) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.enabled[constraint] = enable;   
}

bool SimbodyMatterSubsystemRep::getUseEulerAngles(const State& s) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.useEulerAngles;
}
bool SimbodyMatterSubsystemRep::isMobilizerPrescribed(const State& s, MobilizedBodyId body) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.prescribed[body];
}
bool SimbodyMatterSubsystemRep::isConstraintEnabled(const State& s, ConstraintId constraint) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.enabled[constraint];
}

bool SimbodyMatterSubsystemRep::isUsingQuaternion(const State& s, MobilizedBodyId body) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    return n.isUsingQuaternion(getModelVars(s));
}

int SimbodyMatterSubsystemRep::getNQuaternionsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.nQuaternionsInUse;
}

int SimbodyMatterSubsystemRep::getQuaternionIndex(const State& s, MobilizedBodyId body) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.quaternionIndex[body];
}

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyAcceleration(const State& s, MobilizedBodyId body) const
  { return getRigidBodyNode(body).getA_GB(getAccelerationCache(s)); }

void SimbodyMatterSubsystemRep::enforcePositionConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const {
    const SBModelVars& mv = getModelVars(s);
    Vector&            q  = updQ(s); //TODO: this invalidates q's already

    bool anyChange = false;
 
    // First, fix the position constraints produced by defined length constraints.
    if (lConstraints->enforcePositionConstraints(s, requiredTol, desiredTol))
        anyChange = true;

    // By design, normalization of quaternions can't have any effect on the length
    // constraints we just fixed (because we normalize internally for calculations).
    // So now we can simply normalize the quaternions.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            if (rbNodeLevels[i][j]->enforceQuaternionConstraints(mv,q))
                anyChange = true;
    //TODO: quaternion constraints shouldn't invalidate anything except
    // the qnorms, which will be all 1 now

    if (anyChange)
        s.invalidateAll(Stage::Position);
}

void SimbodyMatterSubsystemRep::enforceVelocityConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const {
    assert(getStage(s) >= Stage::Velocity-1);

    // Currently there are no coordinate constraints for velocity.

    // Fix the velocity constraints produced by defined length constraints.
    const bool anyChange = lConstraints->enforceVelocityConstraints(s, requiredTol, desiredTol);

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
    
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(pc,vc,udotErr,ac);
    for (int i=0; i < (int)pointInPlaneConstraints.size(); ++i)
        pointInPlaneConstraints[i]->calcAccInfo(pc,vc,udotErr,ac);
}


// LOOP FORWARD DYNAMICS OPERATOR
// 
// Given a State realized through Stage::Dynamics, and a complete
// set of applied forces, calculate all acceleration results resulting
// from those forces AND enforcement of the accelerati constraints.
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

    Vector_<SpatialVec> cFrc(getNBodies()); 
    cFrc.setToZero();

    calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                    0, 0, ac, udot, udotErr);
    if (lConstraints->calcConstraintForces(s, udotErr, multipliers, ac)) {
        lConstraints->addInCorrectionForces(s, ac, cFrc);
        calcTreeForwardDynamicsOperator(s, mobilityForces, particleForces, bodyForces,
                                        0, &cFrc, ac, udot, udotErr);
    }
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
calcMobilizerQDotFromU(const State& s, MobilizedBodyId mb, int nu, const Real* u, 
                       int nq, Real* qdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Position);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQ(sbState.getModelVars()));

    n.calcLocalQDotFromLocalU(sbState, u, qdot);
}

// State must be realized to Stage::Velocity, so that we can extract Q(q), QDot(q,u), and u from it to calculate
// qdotdot=Q(q)*udot + QDot(q,u)*u for this mobilizer.
void SimbodyMatterSubsystemRep::
calcMobilizerQDotDotFromUDot(const State& s, MobilizedBodyId mb, int nu, const Real* udot, 
                                  int nq, Real* qdotdot) const
{
    const SBStateDigest sbState(s, *this, Stage::Velocity);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nu == n.getDOF());
    assert(nq == n.getNQ(sbState.getModelVars()));

    n.calcLocalQDotDotFromLocalUDot(sbState, udot, qdotdot);
}

// State must be realized through Stage::Instance. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of q's is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized coordinates.
// Returns X_FM(q).
Transform SimbodyMatterSubsystemRep::
calcMobilizerTransformFromQ(const State& s, MobilizedBodyId mb, int nq, const Real* q) const {
    const SBStateDigest sbState(s, *this, Stage::Instance);
    const RigidBodyNode& n  = getRigidBodyNode(mb);

    assert(nq == n.getNQ(sbState.getModelVars()));

    return n.calcMobilizerTransformFromQ(sbState, q);
}

// State must be realized through Stage::Position. Neither the State nor its
// cache are modified by this method, since it is an operator.
// The number of u's is passed in as a sanity check, to make sure the caller
// and the called mobilizer agree on the generalized speeds.
// Returns V_FM(q,u)=H_FM(q)*u, where the q dependency is extracted from the State via
// the hinge transition matrix H_FM(q).
SpatialVec SimbodyMatterSubsystemRep::
calcMobilizerVelocityFromU(const State& s, MobilizedBodyId mb, int nu, const Real* u) const {
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
calcMobilizerAccelerationFromUDot(const State& s, MobilizedBodyId mb, int nu, const Real* udot) const{
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
calcParentToChildTransformFromQ(const State& s, MobilizedBodyId mb, int nq, const Real* q) const {
    const Transform& X_PF = getMobilizerFrameOnParent(s,mb);
    const Transform& X_BM = getMobilizerFrame(s,mb);

    const Transform X_FM = calcMobilizerTransformFromQ(s,mb,nq,q);
    return X_PF * X_FM * ~X_BM; // X_PB
}

SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildVelocityFromU(const State& s, MobilizedBodyId mb, int nu, const Real* u) const {
    assert(!"not implemented yet");
    return SpatialVec(Vec3(0),Vec3(0));
}
SpatialVec SimbodyMatterSubsystemRep::
calcParentToChildAccelerationFromUDot(const State& s, MobilizedBodyId mb, int nu, const Real* udot) const {
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
    assert(g <= matter.getStage(state).next());
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
        q = &matter.getQ(state)[0];
        if (mc->qVarsIndex >= 0)
            pv = &matter.getPositionVars(state);
        if (mc->qCacheIndex >= 0)
            pc = &matter.updPositionCache(state);
    }
    if (g >= Stage::Velocity) {
        u = &matter.getU(state)[0];
        if (mc->uVarsIndex >= 0)
            vv = &matter.getVelocityVars(state);
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
        if (mc->accelerationCacheIndex >= 0)
            ac = &matter.updAccelerationCache(state);
    }

    stage = g;
}

