/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors: Derived from IVM code written by Charles Schwieters.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of SimbodyMatterSubsystemRep.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "ConstraintNode.h"
#include "LengthConstraints.h"
#include "MultibodySystemRep.h"

#include <string>

SimbodyMatterSubsystemRep::SimbodyMatterSubsystemRep(const SimbodyMatterSubsystemRep& src)
   : SimTK::MatterSubsystemRep("SimbodyMatterSubsystemRep", "X.X.X")
{
    assert(!"SimbodyMatterSubsystemRep copy constructor ... TODO!");
}


SimbodyMatterSubsystemRep::~SimbodyMatterSubsystemRep() {
    delete lConstraints; lConstraints=0;

    for (int i=0; i<(int)constraintNodes.size(); ++i)
        delete constraintNodes[i];
    constraintNodes.resize(0);

    for (int i=0; i<(int)distanceConstraints.size(); ++i)
        delete distanceConstraints[i];
    distanceConstraints.resize(0);

    for (int i=0; i<(int)rbNodeLevels.size(); ++i) {
        for (int j=0; j<(int)rbNodeLevels[i].size(); ++j) 
            delete rbNodeLevels[i][j];
        rbNodeLevels[i].resize(0);
    }
    rbNodeLevels.resize(0);
}

int SimbodyMatterSubsystemRep::addRigidBodyNode
    (RigidBodyNode&           parent,
     const MassProperties&    m,            // mass properties in body frame
     const Transform&         X_PMb,        // parent's frame for attaching this joint
     const Transform&         X_BM,         // inboard joint frame J in body frame
     Mobilizer::MobilizerType type,
     bool                     isReversed,   // child-to-parent orientation?
     int&                     nxtU,
     int&                     nxtUSq,
     int&                     nxtQ)
{
    RigidBodyNode* n = RigidBodyNode::create(m,X_PMb,X_BM,type,isReversed,nxtU,nxtUSq,nxtQ);
    const int level = parent.getLevel() + 1;
    n->setLevel(level);

    // Put node in tree at the right level
    if ((int)rbNodeLevels.size()<=level) rbNodeLevels.resize(level+1);
    const int nxt = rbNodeLevels[level].size();
    rbNodeLevels[level].push_back(n);

    // Assign a unique reference integer to this node, for use by caller
    const int nodeNum = nodeNum2NodeMap.size();
    nodeNum2NodeMap.push_back(RigidBodyNodeIndex(level,nxt));
    n->setNodeNum(nodeNum);

    // Link in to the tree topology (bidirectional).
    parent.addChild(n);
    n->setParent(&parent);

    return nodeNum;
}

// Add a new ground node. Must be first node added during construction.
void SimbodyMatterSubsystemRep::addGroundNode() {
    // Make sure this is the first body
    assert(nodeNum2NodeMap.size() == 0);
    assert(rbNodeLevels.size() == 0);

    RigidBodyNode* n = 
        RigidBodyNode::create(MassProperties(), Transform(), Transform(), 
                              Mobilizer::ThisIsGround,
                              false, nextUSlot, nextUSqSlot, nextQSlot);
    n->setLevel(0);

    // Put ground node in tree at level 0
    rbNodeLevels.resize(1);
    rbNodeLevels[0].push_back(n);
    nodeNum2NodeMap.push_back(RigidBodyNodeIndex(0,0));
    n->setNodeNum(0);
}

ConstraintId SimbodyMatterSubsystemRep::addConstantDistanceConstraint(
    const RigidBodyNode& parent, const Vec3& stationInP,
    const RigidBodyNode& child,  const Vec3& stationInC,
    const Real& distance)
{
    ConstraintNode* cn = new ConstantDistanceConstraintNode(parent,stationInP,child,stationInC,distance);
    return addConstraintNode(cn);
}

ConstraintId SimbodyMatterSubsystemRep::addCoincidentStationsConstraint(
    const RigidBodyNode& parent, const Vec3& stationInP,
    const RigidBodyNode& child,  const Vec3& stationInC)
{
    ConstraintNode* cn = new CoincidentStationsConstraintNode(parent,stationInP,child,stationInC);
    return addConstraintNode(cn);
}


ConstraintId SimbodyMatterSubsystemRep::addWeldConstraint(
    const RigidBodyNode& parent, const Transform& frameInP,
    const RigidBodyNode& child,  const Transform& frameInC)
{
    ConstraintNode* cn = new WeldConstraintNode(parent,frameInP,child,frameInC);
    return addConstraintNode(cn);
}

// Store an already-allocated abstract constraint in the RigidBody tree, assigning
// it a constraint number which is returned. The SimbodyMatterSubsystemRep takes over ownership
// of the ConstraintNode; don't use the pointer any more!
ConstraintId SimbodyMatterSubsystemRep::addConstraintNode(ConstraintNode*& cn) {
    cn->setConstraintNum(constraintNodes.size());
    constraintNodes.push_back(cn);
    cn = 0; // it's all mine now!
    return ConstraintId(constraintNodes.size()-1);
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

BodyId SimbodyMatterSubsystemRep::getParent(BodyId body) const { 
    return BodyId(getRigidBodyNode(body).getParent()->getNodeNum());
}

Array<BodyId> SimbodyMatterSubsystemRep::getChildren(BodyId body) const {
    const RigidBodyNode& node = getRigidBodyNode(body);
    Array<BodyId> children;
    for (BodyId i(0); i < node.getNChildren(); ++i)
        children += BodyId(node.getChild(i)->getNodeNum());
    return children;
}

const MassProperties&
SimbodyMatterSubsystemRep::getBodyMassProperties(const State&, BodyId body) const
  { return getRigidBodyNode(body).getMassProperties_OB_B(); }

const Transform&
SimbodyMatterSubsystemRep::getMobilizerFrame(const State&, BodyId body) const
  { return getRigidBodyNode(body).getX_BM(); }

const Transform&
SimbodyMatterSubsystemRep::getMobilizerFrameOnParent(const State&, BodyId body) const
  { return getRigidBodyNode(body).getX_PMb(); }

const Transform&
SimbodyMatterSubsystemRep::getBodyTransform(const State& s, BodyId body) const
  { return getRigidBodyNode(body).getX_GB(getPositionCache(s)); }

const SpatialVec&
SimbodyMatterSubsystemRep::getBodyVelocity(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getV_GB(getVelocityCache(s));
}

const SpatialVec&
SimbodyMatterSubsystemRep::getCoriolisAcceleration(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getTotalCoriolisAcceleration(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getTotalCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getGyroscopicForce(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getGyroscopicForce(getDynamicsCache(s));
}
const SpatialVec&
SimbodyMatterSubsystemRep::getCentrifugalForces(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getCentrifugalForces(getDynamicsCache(s));
}

const SpatialMat&
SimbodyMatterSubsystemRep::getArticulatedBodyInertia(const State& s, BodyId body) const {
  return getRigidBodyNode(body).getArticulatedBodyInertia(getDynamicsCache(s));
}

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes we're going to need later for state
// variables and cache entries. We allocate and initialize all the
// Modeling variables here.
void SimbodyMatterSubsystemRep::endConstruction() {
    if (built) return; // already done

    // Not built yet. Let's count topological things.

    DOFTotal = SqDOFTotal = maxNQTotal = 0;
    for (int i=0; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0; j<(int)rbNodeLevels[i].size(); ++j) {
            const int ndof = rbNodeLevels[i][j]->getDOF();
            DOFTotal += ndof; SqDOFTotal += ndof*ndof;
            maxNQTotal += rbNodeLevels[i][j]->getMaxNQ();
        }

    // Dole out the multipliers and constraint error slots for topological
    // constraints. TODO: currently our constraints are all holonomic, meaning
    // position-level, so that they occupy one slot in the qErr array, then
    // their time derivatives need one slot in uErr and their 2nd time 
    // derivatives need one acceleration-level multiplier. Later we will
    // have constraints which start at the velocity level (nonholonomic) so
    // they won't use up a qErr slot.
    // Also, quaternion normalization constraints exist only at the 
    // position level, however they are not topological since modeling
    // choices affect whether we use them. See realizeModel() below.
    for (int i=0; i<(int)constraintNodes.size(); ++i) {
        ConstraintNode& cons = *constraintNodes[i];
        cons.setQErrIndex(nextQErrSlot);
        cons.setUErrIndex(nextUErrSlot);
        cons.setMultIndex(nextMultSlot);
        const int nConsEqns = cons.getNConstraintEquations();
        nextQErrSlot += nConsEqns;
        nextUErrSlot += nConsEqns;
        nextMultSlot += nConsEqns;
        cons.finishConstruction(*this);
    }

    lConstraints = new LengthConstraints(*this, 0);
    lConstraints->construct(distanceConstraints);
    built = true;
}

void SimbodyMatterSubsystemRep::realizeTopology(State& s) const {
    // This is a long-winded way of saying that the Stage must be exactly Empty.
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Topology).prev(), 
        "SimbodyMatterSubsystemRep::realizeTopology()");
    SimTK_STAGECHECK_LT_ALWAYS(getStage(s), Stage(Stage::Topology), 
        "SimbodyMatterSubsystemRep::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this call*.
    // Afterwards they are truly const so we don't declare them mutable, but cheat
    // here instead.
    SimbodyMatterSubsystemRep* mutableThis = const_cast<SimbodyMatterSubsystemRep*>(this);

    if (!built) mutableThis->endConstruction(); // no more bodies after this!

    // Fill in the local copy of the topologyCache from the information
    // calculated in endConstruction(). Also ask the State for some room to
    // put Modeling variables & cache and remember the indices in our construction
    // cache.

    mutableThis->topologyCache.nBodies      = nodeNum2NodeMap.size();
    mutableThis->topologyCache.nConstraints = constraintNodes.size();
    mutableThis->topologyCache.nDOFs        = DOFTotal;
    mutableThis->topologyCache.maxNQs       = maxNQTotal;
    mutableThis->topologyCache.sumSqDOFs    = SqDOFTotal;
    mutableThis->topologyCache.nDistanceConstraints = distanceConstraints.size();

    SBModelVars mvars;
    mvars.allocate(topologyCache);
    setDefaultModelValues(topologyCache, mvars);
    mutableThis->topologyCache.modelingVarsIndex  = 
        s.allocateDiscreteVariable(getMySubsystemIndex(),Stage::Model, new Value<SBModelVars>(mvars));

    mutableThis->topologyCache.modelingCacheIndex = 
        s.allocateCacheEntry(getMySubsystemIndex(),Stage::Model, new Value<SBModelCache>());

    mutableThis->topologyCache.valid = true;

    // Allocate a cache entry for the topologyCache, and save a copy there.
    mutableThis->topologyCacheIndex = 
        s.allocateCacheEntry(getMySubsystemIndex(),Stage::Topology, new Value<SBTopologyCache>(topologyCache));
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
void SimbodyMatterSubsystemRep::realizeModel(State& s) const {
    // This is a long-winded way of saying that the Stage must be exactly Built.
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Model).prev(), 
        "SimbodyMatterSubsystemRep::realizeModel()");
    SimTK_STAGECHECK_LT_ALWAYS(getStage(s), Stage(Stage::Model), 
        "SimbodyMatterSubsystemRep::realizeModel()");

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
        s.allocateDiscreteVariable(getMySubsystemIndex(),Stage::Instance, 
                                   new Value<SBInstanceVars>(iv));
    mc.instanceCacheIndex = 
        s.allocateCacheEntry(getMySubsystemIndex(),Stage::Instance, 
                             new Value<SBInstanceCache>());

    // No time vars or cache
    mc.timeVarsIndex = -1;
    mc.timeCacheIndex = -1;

    // Position variables are just q's, which the State knows how to deal with. 

    Vector qInit(maxNQTotal);
    setDefaultPositionValues(mv, qInit);

    mc.qIndex = s.allocateQ(getMySubsystemIndex(), qInit);
    mc.qVarsIndex = -1; // no config vars other than q
    mc.qCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(),Stage::Position, 
        new Value<SBPositionCache>());

    // We'll store the the physical constraint errors (which consist solely of distance
    // constraint equations at the moment), followed by the quaternion constraints.
    mc.qErrIndex = s.allocateQErr(getMySubsystemIndex(), 
                                  nextQErrSlot + mc.nQuaternionsInUse);

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.

    Vector uInit(DOFTotal);
    setDefaultVelocityValues(mv, uInit);

    mc.uIndex = s.allocateU(getMySubsystemIndex(), uInit);
    mc.uVarsIndex = -1; // no velocity vars other than u
    mc.uCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(),Stage::Velocity, 
        new Value<SBVelocityCache>());
    // Note that qdots are automatically allocated in the Velocity stage cache.

    // Only physical constraints exist at the velocity and acceleration levels; 
    // the quaternion normalization constraints are gone.
    mc.uErrIndex    = s.allocateUErr(getMySubsystemIndex(),    nextUErrSlot);
    mc.udotErrIndex = s.allocateUDotErr(getMySubsystemIndex(), nextMultSlot);

    // no z's
    // We do have dynamic vars for now for forces & pres. accel. but those will
    // probably go away. TODO
    SBDynamicsVars dvars;
    dvars.allocate(topologyCache);
    setDefaultDynamicsValues(mv, dvars);
    mc.dynamicsVarsIndex = 
        s.allocateDiscreteVariable(getMySubsystemIndex(),Stage::Dynamics, 
                                   new Value<SBDynamicsVars>(dvars));
    mc.dynamicsCacheIndex  = 
        s.allocateCacheEntry(getMySubsystemIndex(),Stage::Dynamics, 
                             new Value<SBDynamicsCache>());

    // No reaction variables that I know of. But we can go through the
    // charade here anyway.
    SBAccelerationVars rvars;
    rvars.allocate(topologyCache);
    setDefaultAccelerationValues(mv, rvars);

    mc.accelerationVarsIndex = 
        s.allocateDiscreteVariable(getMySubsystemIndex(),Stage::Acceleration, 
                                   new Value<SBAccelerationVars>(rvars));
    mc.accelerationCacheIndex = 
        s.allocateCacheEntry(getMySubsystemIndex(),Stage::Acceleration, 
                             new Value<SBAccelerationCache>());

    // Note that qdots, qdotdots, udots, zdots are automatically allocated by
    // the State when we advance the stage past modeling.
}

// Here we lock in parameterization of the model, such as body masses.
void SimbodyMatterSubsystemRep::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "SimbodyMatterSubsystemRep::realizeInstance()");

    const SBModelVars&    mv = getModelVars(s);
    const SBInstanceVars& iv = getInstanceVars(s);

    // Get the Instance-stage cache and make sure it has been allocated and initialized if needed.
    SBInstanceCache& ic = updInstanceCache(s);
    ic.allocate(topologyCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeInstance(mv,iv,ic); 
}

void SimbodyMatterSubsystemRep::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "SimbodyMatterSubsystemRep::realizeTime()");

    // nothing yet
}

// Set generalized coordinates: sweep from base to tips.
void SimbodyMatterSubsystemRep::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "SimbodyMatterSubsystemRep::realizePosition()");

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

    // Constraint errors go in qErr after the quaternion constraints.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(qErr,pc); //TODO: qErr
}

// Set generalized speeds: sweep from base to tip.
// realizePosition() must have been called already.
void SimbodyMatterSubsystemRep::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "SimbodyMatterSubsystemRep::realizeVelocity()");

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
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.
// Then go ask around to collect up all the applied forces from any
// force subsystems.

void SimbodyMatterSubsystemRep::realizeDynamics(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "SimbodyMatterSubsystemRep::realizeDynamics()");
    const SBPositionCache& pc = getPositionCache(s);
    const SBVelocityCache& vc = getVelocityCache(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    SBDynamicsCache&            dc = updDynamicsCache(s);
    dc.allocate(topologyCache);

    calcArticulatedBodyInertias(s);
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcJointIndependentDynamicsVel(pc,vc,dc);

    // Now total up all the forces
    // TODO: this shouldn't be copying!!
    const MultibodySystem& mbs = getMultibodySystem();  // owner of this subsystem
    mbs.getRep().updKineticEnergy(s) += calcKineticEnergy(s);

    dc.appliedMobilityForces  = mbs.getMobilityForces(s);
    dc.appliedParticleForces  = mbs.getParticleForces(s);
    dc.appliedRigidBodyForces = mbs.getRigidBodyForces(s);
}

void SimbodyMatterSubsystemRep::realizeAcceleration(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "SimbodyMatterSubsystemRep::realizeAcceleration()");

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    Vector&              udot    = updUDot(s);
    Vector&              qdotdot = updQDotDot(s);
    SBAccelerationCache& ac      = updAccelerationCache(s);
    ac.allocate(topologyCache);

    calcLoopForwardDynamics(s);
    calcQDotDot(s, udot, qdotdot);
}


int SimbodyMatterSubsystemRep::getQIndex(BodyId body) const 
  { assert(built);return getRigidBodyNode(body).getQIndex();}
int SimbodyMatterSubsystemRep::getQAlloc(BodyId body) const 
  { assert(built);return getRigidBodyNode(body).getMaxNQ();}
int SimbodyMatterSubsystemRep::getUIndex(BodyId body) const
  { assert(built);return getRigidBodyNode(body).getUIndex();}
int SimbodyMatterSubsystemRep::getDOF   (BodyId body) const
  { assert(built);return getRigidBodyNode(body).getDOF();}

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
            rbNodeLevels[i][j]->setDefaultModelValues(topologyCache, modelVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultInstanceValues(const SBModelVars& mv, 
                                             SBInstanceVars& paramVars) const 
{
    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultInstanceValues(mv, paramVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultTimeValues(const SBModelVars& mv, 
                                         SBTimeVars& timeVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultTimeValues(mv, timeVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultPositionValues(mv, q);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultVelocityValues(const SBModelVars& mv, Vector& u) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultVelocityValues(mv, u);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultDynamicsValues(const SBModelVars& mv, 
                                             SBDynamicsVars& dynamicsVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultDynamicsValues(mv, dynamicsVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setDefaultAccelerationValues(const SBModelVars& mv, 
                                             SBAccelerationVars& reactionVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultAccelerationValues(mv, reactionVars);

    // TODO: constraint defaults
}

void SimbodyMatterSubsystemRep::setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void SimbodyMatterSubsystemRep::setMobilizerIsPrescribed(State& s, BodyId body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void SimbodyMatterSubsystemRep::setConstraintIsEnabled(State& s, int constraint, bool enable) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.enabled[constraint] = enable;   
}

bool SimbodyMatterSubsystemRep::getUseEulerAngles(const State& s) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.useEulerAngles;
}
bool SimbodyMatterSubsystemRep::isMobilizerPrescribed(const State& s, BodyId body) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.prescribed[body];
}
bool SimbodyMatterSubsystemRep::isConstraintEnabled(const State& s, int constraint) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.enabled[constraint];
}

bool SimbodyMatterSubsystemRep::isUsingQuaternion(const State& s, BodyId body) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    return n.isUsingQuaternion(getModelVars(s));
}

int SimbodyMatterSubsystemRep::getNQuaternionsInUse(const State& s) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.nQuaternionsInUse;
}

int SimbodyMatterSubsystemRep::getQuaternionIndex(const State& s, BodyId body) const {
    const SBModelCache& mc = getModelCache(s); // must be >=Model stage
    return mc.quaternionIndex[body];
}

const Transform& SimbodyMatterSubsystemRep::getMobilizerTransform(const State& s, BodyId body) const { 
    const RigidBodyNode& n = getRigidBodyNode(body);
    const SBPositionCache& cc = getPositionCache(s);
    return n.getX_MbM(cc);
}
const SpatialVec& SimbodyMatterSubsystemRep::getMobilizerVelocity(const State& s, BodyId body) const { 
    const RigidBodyNode& n  = getRigidBodyNode(body);
    const SBVelocityCache& mc = getVelocityCache(s);
    return n.getV_MbM(mc);
}
void SimbodyMatterSubsystemRep::setMobilizerTransform(State& s, BodyId body, const Transform& X_MbM) const { 
    const RigidBodyNode& n  = getRigidBodyNode(body);
    const SBModelVars&   mv = getModelVars(s);
    Vector& q = updQ(s);
    n.setMobilizerTransform(mv, X_MbM, q);
}
void SimbodyMatterSubsystemRep::setMobilizerVelocity(State& s, BodyId body, const SpatialVec& V_MbM) const { 
    const RigidBodyNode& n  = getRigidBodyNode(body);
    const SBModelVars&   mv = getModelVars(s);
    Vector& u = updU(s);
    n.setMobilizerVelocity(mv, V_MbM, u);
}

const Vector& 
SimbodyMatterSubsystemRep::getAppliedMobilityForces(const State& s) const {
    const SBDynamicsCache& dc = getDynamicsCache(s);
    return dc.appliedMobilityForces;
}
const Vector_<SpatialVec>& 
SimbodyMatterSubsystemRep::getAppliedBodyForces(const State& s) const {
    const SBDynamicsCache& dc = getDynamicsCache(s);
    return dc.appliedRigidBodyForces;
}


const SpatialVec&
SimbodyMatterSubsystemRep::getBodyAcceleration(const State& s, BodyId body) const
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

// Given a forces in the state, calculate accelerations ignoring
// constraints, and leave the results in the state. 
// Must have already called realizeDynamics().
// We also allow some extra forces to be supplied, with the intent
// that these will be used to deal with internal forces generated
// by constraints. 
void SimbodyMatterSubsystemRep::calcTreeForwardDynamics(
    const State&               s,
    const Vector*              extraMobilityForces,
    const Vector_<SpatialVec>* extraBodyForces) const
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBVelocityCache& mc = getVelocityCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);

    Vector              totalMobilityForces;
    Vector_<SpatialVec> totalBodyForces;

    // inputs
    const Vector&              mobForces  = dc.appliedMobilityForces;
    const Vector_<SpatialVec>& bodyForces = dc.appliedRigidBodyForces;

    const Vector*              mobForcesToUse  = &mobForces;
    const Vector_<SpatialVec>* bodyForcesToUse = &bodyForces;

    if (extraMobilityForces) {
        totalMobilityForces = mobForces + *extraMobilityForces;
        mobForcesToUse      = &totalMobilityForces;
    }

    if (extraBodyForces) {
        totalBodyForces = bodyForces + *extraBodyForces;
        bodyForcesToUse = &totalBodyForces;
    }

    // outputs
    Vector&              udotErr        = updUDotErr(s);
    SBAccelerationCache& rc             = updAccelerationCache(s);
    Vector&              netHingeForces = rc.netHingeForces;
    Vector_<SpatialVec>& A_GB           = rc.bodyAccelerationInGround;

    Vector&              udot           = updUDot(s);

    calcTreeAccelerations(s, *mobForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot);
    
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(cc,mc,udotErr,rc);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void SimbodyMatterSubsystemRep::calcLoopForwardDynamics(const State& s) const 
{
    assert(getStage(s) >= Stage::Acceleration-1);

    Vector_<SpatialVec> cFrc(getNBodies()); 
    cFrc.setToZero();

    calcTreeForwardDynamics(s, 0, 0);
    if (lConstraints->calcConstraintForces(s)) {
        lConstraints->addInCorrectionForces(s, cFrc);
        calcTreeForwardDynamics(s, 0, &cFrc);
    }
}

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
void SimbodyMatterSubsystemRep::calcZ(const State& s, 
                          const SpatialVecList& spatialForces) const
{
    const SBPositionCache& pc = getPositionCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);
    SBAccelerationCache&   ac = updAccelerationCache(s);


    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(pc,dc,spatialForces[node.getNodeNum()], ac);
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
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBPositionCache& pc = getPositionCache(s);
    const SBDynamicsCache& dc = getDynamicsCache(s);

    assert(jointForces.size() == getTotalDOF());
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
                jointForces, bodyForces, allZ, allGepsilon,
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

// If V is a spatial velocity, and you have a X=d(something)/dV (one per body)
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
            RigidBodyNode& node = *rbNodeLevels[i][j];
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

