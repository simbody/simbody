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
 * Implementation of RigidBodyTree.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "RigidBodyTree.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "ConstraintNode.h"
#include "LengthConstraints.h"
#include "MultibodySystemRep.h"

#include <string>

RigidBodyTree::RigidBodyTree(const RigidBodyTree& src)
   : SimTK::MatterSubsystemRep("RigidBodyTree", "X.X.X")
{
    assert(!"RigidBodyTree copy constructor ... TODO!");
}


RigidBodyTree::~RigidBodyTree() {
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

int RigidBodyTree::addRigidBodyNode
    (RigidBodyNode&           parent,
     const MassProperties&    m,            // mass properties in body frame
     const Transform&         X_PJb,        // parent's frame for attaching this joint
     const Transform&         X_BJ,         // inboard joint frame J in body frame
     Mobilizer::MobilizerType type,
     bool                     isReversed,   // child-to-parent orientation?
     int&                     nxtU,
     int&                     nxtUSq,
     int&                     nxtQ)
{
    RigidBodyNode* n = RigidBodyNode::create(m,X_PJb,X_BJ,type,isReversed,nxtU,nxtUSq,nxtQ);
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
void RigidBodyTree::addGroundNode() {
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

int RigidBodyTree::addConstantDistanceConstraint(
    const RigidBodyNode& parent, const Vec3& stationInP,
    const RigidBodyNode& child,  const Vec3& stationInC,
    const Real& distance)
{
    ConstraintNode* cn = new ConstantDistanceConstraintNode(parent,stationInP,child,stationInC,distance);
    return addConstraintNode(cn);
}

int RigidBodyTree::addCoincidentStationsConstraint(
    const RigidBodyNode& parent, const Vec3& stationInP,
    const RigidBodyNode& child,  const Vec3& stationInC)
{
    ConstraintNode* cn = new CoincidentStationsConstraintNode(parent,stationInP,child,stationInC);
    return addConstraintNode(cn);
}


int RigidBodyTree::addWeldConstraint(
    const RigidBodyNode& parent, const Transform& frameInP,
    const RigidBodyNode& child,  const Transform& frameInC)
{
    ConstraintNode* cn = new WeldConstraintNode(parent,frameInP,child,frameInC);
    return addConstraintNode(cn);
}

// Store an already-allocated abstract constraint in the RigidBody tree, assigning
// it a constraint number which is returned. The RigidBodyTree takes over ownership
// of the ConstraintNode; don't use the pointer any more!
int RigidBodyTree::addConstraintNode(ConstraintNode*& cn) {
    cn->setConstraintNum(constraintNodes.size());
    constraintNodes.push_back(cn);
    cn = 0; // it's all mine now!
    return constraintNodes.size()-1;
}

// Add a distance constraint and assign it to use a particular multiplier. 
// Return the assigned distance constraint index for caller's use.
int RigidBodyTree::addOneDistanceConstraintEquation(
    const RBStation& s1, const RBStation& s2, const Real& d,
    int multIndex)
{
    RBDistanceConstraint* dc = new RBDistanceConstraint(s1,s2,d);
    dc->setMultIndex(multIndex);
    dc->setDistanceConstraintNum(distanceConstraints.size());
    distanceConstraints.push_back(dc);
    return distanceConstraints.size()-1;
}

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes we're going to need later for state
// variables and cache entries. We allocate and initialize all the
// Modeling variables here.
void RigidBodyTree::endConstruction() {
    if (built) return; // already done

    // Not built yet. Let's count topological things.

    DOFTotal = SqDOFTotal = maxNQTotal = 0;
    for (int i=0; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0; j<(int)rbNodeLevels[i].size(); ++j) {
            const int ndof = rbNodeLevels[i][j]->getDOF();
            DOFTotal += ndof; SqDOFTotal += ndof*ndof;
            maxNQTotal += rbNodeLevels[i][j]->getMaxNQ();
        }

    int nxtMultIndex = 0; // dole out the multipliers
    for (int i=0; i<(int)constraintNodes.size(); ++i) {
        constraintNodes[i]->setMultIndex(nxtMultIndex); // must assign mults first
        nxtMultIndex += constraintNodes[i]->getNMult();
        constraintNodes[i]->finishConstruction(*this);
    }

    lConstraints = new LengthConstraints(*this, 0);
    lConstraints->construct(distanceConstraints);
    built = true;
}


int RigidBodyTree::getParent(int body) const { 
    return getRigidBodyNode(body).getParent()->getNodeNum();
}

Array<int> RigidBodyTree::getChildren(int body) const {
    const RigidBodyNode& node = getRigidBodyNode(body);
    Array<int> children;
    for (int i=0; i < node.getNChildren(); ++i)
        children += node.getChild(i)->getNodeNum();
    return children;
}

const Real&
RigidBodyTree::getBodyMass(const State&, int body) const
  { return getRigidBodyNode(body).getMass(); }

const Vec3&
RigidBodyTree::getBodyCenterOfMassStation(const State&, int body) const
  { return getRigidBodyNode(body).getCOM_B(); }


const Transform&
RigidBodyTree::getMobilizerFrame(const State&, int body) const
  { return getRigidBodyNode(body).getX_BJ(); }

const Transform&
RigidBodyTree::getMobilizerFrameOnParent(const State&, int body) const
  { return getRigidBodyNode(body).getX_PJb(); }

const Transform&
RigidBodyTree::getBodyPosition(const State& s, int body) const
  { return getRigidBodyNode(body).getX_GB(getPositionCache(s)); }

const SpatialVec&
RigidBodyTree::getBodyVelocity(const State& s, int body) const {
  return getRigidBodyNode(body).getV_GB(getMotionCache(s));
}

const SpatialVec&
RigidBodyTree::getCoriolisAcceleration(const State& s, int body) const {
  return getRigidBodyNode(body).getCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
RigidBodyTree::getTotalCoriolisAcceleration(const State& s, int body) const {
  return getRigidBodyNode(body).getTotalCoriolisAcceleration(getDynamicsCache(s));
}
const SpatialVec&
RigidBodyTree::getGyroscopicForce(const State& s, int body) const {
  return getRigidBodyNode(body).getGyroscopicForce(getDynamicsCache(s));
}
const SpatialVec&
RigidBodyTree::getCentrifugalForces(const State& s, int body) const {
  return getRigidBodyNode(body).getCentrifugalForces(getDynamicsCache(s));
}

const SpatialMat&
RigidBodyTree::getArticulatedBodyInertia(const State& s, int body) const {
  return getRigidBodyNode(body).getArticulatedBodyInertia(getDynamicsCache(s));
}

void RigidBodyTree::realizeTopology(State& s) const {
    // This is a long-winded way of saying that the Stage must be exactly Empty.
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Topology).prev(), 
        "RigidBodyTree::realizeTopology()");
    SimTK_STAGECHECK_LT_ALWAYS(getStage(s), Stage(Stage::Topology), 
        "RigidBodyTree::realizeTopology()");

    // Some of our 'const' values must be treated as mutable *just for this call*.
    // Afterwards they are truly const so we don't declare them mutable, but cheat
    // here instead.
    RigidBodyTree* mutableThis = const_cast<RigidBodyTree*>(this);

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
void RigidBodyTree::realizeModel(State& s) const {
    // This is a long-winded way of saying that the Stage must be exactly Built.
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Model).prev(), 
        "RigidBodyTree::realizeModel()");
    SimTK_STAGECHECK_LT_ALWAYS(getStage(s), Stage(Stage::Model), 
        "RigidBodyTree::realizeModel()");

    const SBModelVars&   mv = getModelVars(s);

    // Get the Modeling-stage cache and make sure it has been allocated and initialized if needed.
    // It is OK to hold a reference here because the discrete variables (and cache entries) in
    // the State are stable, that is, they don't change even if more variables are added.
    SBModelCache&        mc = updModelCache(s);
    mc.allocate(topologyCache);

    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModel(mv,mc); 

    // Now allocate all remaining variables and cache entries. We can properly initialize only
    // the next stage, the parameters. Everything else gets some kind of meaningless initial
    // values but those could change as we realize the higher stages.
    SBInstanceVars pv;
    pv.allocate(topologyCache);
    setDefaultInstanceValues(mv, pv);

    mc.instanceVarsIndex = 
        s.allocateDiscreteVariable(getMySubsystemIndex(),Stage::Instance, 
                                   new Value<SBInstanceVars>(pv));
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

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.

    Vector uInit(DOFTotal);
    setDefaultVelocityValues(mv, uInit);

    mc.uIndex = s.allocateU(getMySubsystemIndex(), uInit);
    mc.uVarsIndex = -1; // no velocity vars other than u
    mc.uCacheIndex = s.allocateCacheEntry(getMySubsystemIndex(),Stage::Velocity, 
        new Value<SBVelocityCache>());
    // Note that qdots are automatically allocated in the Velocity stage cache.

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
void RigidBodyTree::realizeInstance(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "RigidBodyTree::realizeInstance()");

    const SBModelVars&  mv = getModelVars(s);
    const SBInstanceVars& pv = getInstanceVars(s);

    // Get the Instance-stage cache and make sure it has been allocated and initialized if needed.
    SBInstanceCache& pc = updInstanceCache(s);
    pc.allocate(topologyCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeInstance(mv,pv,pc); 
}

void RigidBodyTree::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "RigidBodyTree::realizeTime()");

    // nothing yet
}

// Set generalized coordinates: sweep from base to tips.
void RigidBodyTree::realizePosition(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "RigidBodyTree::realizePosition()");

    const SBModelVars& mv = getModelVars(s);
    const Vector&         q  = getQ(s);

    // Get the Position-stage cache and make sure it has been allocated and initialized if needed.
    SBPositionCache& cc = updPositionCache(s);
    cc.allocate(topologyCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizePosition(mv,q,cc); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(cc);
}

// Set generalized speeds: sweep from base to tip.
// realizePosition() must have been called already.
void RigidBodyTree::realizeVelocity(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "RigidBodyTree::realizeVelocity()");

    const SBModelVars&       mv = getModelVars(s);
    const Vector&               q  = getQ(s);
    const SBPositionCache& cc = getPositionCache(s);
    const Vector&               u  = getU(s);

    // Get the Motion-stage cache and make sure it has been allocated and initialized if needed.
    SBVelocityCache&              mc   = updMotionCache(s);
    mc.allocate(topologyCache);
    Vector&                     qdot = updQDot(s);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeVelocity(mv,q,cc,u,mc,qdot); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcVelInfo(cc,mc);
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.
// Then go ask around to collect up all the applied forces from any
// force subsystems.

void RigidBodyTree::realizeDynamics(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "RigidBodyTree::realizeDynamics()");
    const SBPositionCache& cc = getPositionCache(s);
    const SBVelocityCache&        mc = getMotionCache(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    SBDynamicsCache&            dc = updDynamicsCache(s);
    dc.allocate(topologyCache);

    calcArticulatedBodyInertias(s);
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcJointIndependentDynamicsVel(cc,mc,dc);

    // Now total up all the forces
    // TODO: this shouldn't be copying!!
    const MultibodySystem& mbs = getMultibodySystem();  // owner of this subsystem
    mbs.getRep().updKineticEnergy(s) += calcKineticEnergy(s);

    dc.appliedMobilityForces  = mbs.getMobilityForces(s);
    dc.appliedParticleForces  = mbs.getParticleForces(s);
    dc.appliedRigidBodyForces = mbs.getRigidBodyForces(s);
}

void RigidBodyTree::realizeAcceleration(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "RigidBodyTree::realizeAcceleration()");

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    Vector&             udot    = updUDot(s);
    Vector&             qdotdot = updQDotDot(s);
    SBAccelerationCache&    rc      = updAccelerationCache(s);
    rc.allocate(topologyCache);

    calcLoopForwardDynamics(s);
    calcQDotDot(s, udot, qdotdot);
}


int RigidBodyTree::getQIndex(int body) const 
  { assert(built);return getRigidBodyNode(body).getQIndex();}
int RigidBodyTree::getQAlloc(int body) const 
  { assert(built);return getRigidBodyNode(body).getMaxNQ();}
int RigidBodyTree::getUIndex(int body) const
  { assert(built);return getRigidBodyNode(body).getUIndex();}
int RigidBodyTree::getDOF   (int body) const
  { assert(built);return getRigidBodyNode(body).getDOF();}

// We are in the process of realizingConstruction() when we need to make this call.
// We pass in the partially-completed Construction-stage cache, which must have all
// the dimensions properly filled in at this point.
void RigidBodyTree::setDefaultModelValues(const SBTopologyCache& topologyCache, 
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

void RigidBodyTree::setDefaultInstanceValues(const SBModelVars& mv, 
                                              SBInstanceVars& paramVars) const 
{
    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultInstanceValues(mv, paramVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultTimeValues(const SBModelVars& mv, 
                                         SBTimeVars& timeVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultTimeValues(mv, timeVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultPositionValues(mv, q);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultVelocityValues(const SBModelVars& mv, Vector& u) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultVelocityValues(mv, u);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultDynamicsValues(const SBModelVars& mv, 
                                             SBDynamicsVars& dynamicsVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultDynamicsValues(mv, dynamicsVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultAccelerationValues(const SBModelVars& mv, 
                                             SBAccelerationVars& reactionVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultAccelerationValues(mv, reactionVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setUseEulerAngles(State& s, bool useAngles) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void RigidBodyTree::setMobilizerIsPrescribed(State& s, int body, bool prescribe) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.prescribed[body] = prescribe;
}
void RigidBodyTree::setConstraintIsEnabled(State& s, int constraint, bool enable) const {
    SBModelVars& modelVars = updModelVars(s); // check/adjust stage
    modelVars.enabled[constraint] = enable;   
}

bool RigidBodyTree::getUseEulerAngles(const State& s) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.useEulerAngles;
}
bool RigidBodyTree::isMobilizerPrescribed(const State& s, int body) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.prescribed[body];
}
bool RigidBodyTree::isConstraintEnabled(const State& s, int constraint) const {
    const SBModelVars& modelVars = getModelVars(s); // check stage
    return modelVars.enabled[constraint];
}

const Real& RigidBodyTree::getMobilizerQ(const State& s, int body, int axis) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getNQ(getModelVars(s)));
    return getQ(s)[n.getQIndex()+axis];
}


const Real& RigidBodyTree::getMobilizerU(const State& s, int body, int axis) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    return getU(s)[n.getUIndex()+axis];
}

void RigidBodyTree::setMobilizerQ(State& s, int body, int axis, const Real& r) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getNQ(getModelVars(s)));
    updQ(s)[n.getQIndex()+axis] = r;
}


void RigidBodyTree::setMobilizerU(State& s, int body, int axis, const Real& r) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    updU(s)[n.getUIndex()+axis] = r;
}


const Transform& RigidBodyTree::getMobilizerPosition(const State& s, int body) const { 
    const RigidBodyNode& n = getRigidBodyNode(body);
    const SBPositionCache& cc = getPositionCache(s);
    return n.getX_JbJ(cc);
}
const SpatialVec& RigidBodyTree::getMobilizerVelocity(const State& s, int body) const { 
    const RigidBodyNode& n  = getRigidBodyNode(body);
    const SBVelocityCache& mc = getMotionCache(s);
    return n.getV_JbJ(mc);
}
void RigidBodyTree::setMobilizerPosition(State& s, int body, const Transform& X_JbJ) const { 
    const RigidBodyNode& n = getRigidBodyNode(body);
    const SBModelVars& mv = getModelVars(s);
    Vector& q = updQ(s);
    n.setMobilizerPosition(mv, X_JbJ, q);
}
void RigidBodyTree::setMobilizerVelocity(State& s, int body, const SpatialVec& V_JbJ) const { 
    const RigidBodyNode& n = getRigidBodyNode(body);
    const SBModelVars& mv = getModelVars(s);
    Vector& u = updU(s);
    n.setMobilizerVelocity(mv, V_JbJ, u);
}

void RigidBodyTree::addInStationForce(const State& s, int body, 
                                      const Vec3& stationInB, const Vec3& forceInG, 
                                      Vector_<SpatialVec>& rigidBodyForces) const
{
    assert(rigidBodyForces.size() == getNBodies());
    const SBPositionCache& cc = getPositionCache(s);
    const Rotation& R_GB = getRigidBodyNode(body).getX_GB(cc).R();
    rigidBodyForces[body] += SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void RigidBodyTree::addInBodyTorque(const State& s, int body, const Vec3& torqueInG, 
                                    Vector_<SpatialVec>& rigidBodyForces) const 
{
    assert(rigidBodyForces.size() == getNBodies());
    rigidBodyForces[body][0] += torqueInG; // no force
}

void RigidBodyTree::addInMobilityForce(const State& s, int body, int axis, const Real& r,
                                       Vector& mobilityForces) const 
{
    assert(mobilityForces.size() == getTotalDOF());
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    mobilityForces[n.getUIndex()+axis] += r;
}

const Vector& 
RigidBodyTree::getAppliedMobilityForces(const State& s) const {
    const SBDynamicsCache& dc = getDynamicsCache(s);
    return dc.appliedMobilityForces;
}
const Vector_<SpatialVec>& 
RigidBodyTree::getAppliedBodyForces(const State& s) const {
    const SBDynamicsCache& dc = getDynamicsCache(s);
    return dc.appliedRigidBodyForces;
}


const SpatialVec&
RigidBodyTree::getBodyAcceleration(const State& s, int body) const
  { return getRigidBodyNode(body).getA_GB(getAccelerationCache(s)); }

void RigidBodyTree::enforcePositionConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const {
    const SBModelVars& mv = getModelVars(s);
    Vector&               q  = updQ(s); //TODO: this invalidates q's already

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

void RigidBodyTree::enforceVelocityConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const {
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
void RigidBodyTree::calcTreeForwardDynamics(
    const State&               s,
    const Vector*              extraMobilityForces,
    const Vector_<SpatialVec>* extraBodyForces) const
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBVelocityCache&        mc = getMotionCache(s);
    const SBDynamicsCache&      dc = getDynamicsCache(s);

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
    SBAccelerationCache&     rc             = updAccelerationCache(s);
    Vector&              netHingeForces = rc.netHingeForces;
    Vector_<SpatialVec>& A_GB           = rc.bodyAccelerationInGround;

    Vector&              udot           = updUDot(s);

    calcTreeAccelerations(s, *mobForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot);
    
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(cc,mc,rc);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const State& s) const 
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
void RigidBodyTree::calcArticulatedBodyInertias(const State& s) const {
    const SBPositionCache& cc = getPositionCache(s);
    SBDynamicsCache&            dc = updDynamicsCache(s);

    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcArticulatedBodyInertiasInward(cc,dc);
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcZ(const State& s, 
                          const SpatialVecList& spatialForces) const
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBDynamicsCache&      dc = getDynamicsCache(s);
    SBAccelerationCache&            rc = updAccelerationCache(s);


    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(cc,dc,spatialForces[node.getNodeNum()], rc);
        }
}

// Y is used for length constraints: sweep from base to tip.
void RigidBodyTree::calcY(const State& s) const {
    const SBPositionCache& cc = getPositionCache(s);
    SBDynamicsCache&            dc = updDynamicsCache(s);

    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcYOutward(cc,dc);
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcTreeAccel(const State& s) const {
    const SBModelVars&       mv      = getModelVars(s);
    const Vector&               q       = getQ(s);
    const SBPositionCache& cc      = getPositionCache(s);
    const Vector&               u       = getU(s);
    const SBDynamicsCache&      dc      = getDynamicsCache(s);

    SBAccelerationCache&            rc      = updAccelerationCache(s);
    Vector&                     udot    = updUDot(s);
    Vector&                     qdotdot = updQDotDot(s);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel(mv,q,cc,u,dc,rc,udot,qdotdot);
}

void RigidBodyTree::fixVel0(State& s, Vector& vel) const {
    lConstraints->fixVel0(s, vel);
}

Real RigidBodyTree::calcKineticEnergy(const State& s) const {
    const SBPositionCache& cc = getPositionCache(s);
    const SBVelocityCache&        mc = getMotionCache(s);

    Real ke = 0.;

    // Skip ground level 0!
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            ke += rbNodeLevels[i][j]->calcKineticEnergy(cc,mc);

    return ke;
}

//
// Operator for open-loop dynamics.
//
void RigidBodyTree::calcTreeAccelerations(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBDynamicsCache&      dc = getDynamicsCache(s);

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
            node.calcUDotPass1Inward(cc,dc,
                jointForces, bodyForces, allZ, allGepsilon,
                netHingeForces);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass2Outward(cc,dc, netHingeForces, A_GB, udot);
        }
}

//
// Calculate udot = M^-1 f. We also get spatial accelerations A_GB for 
// each body as a side effect.
//
void RigidBodyTree::calcMInverseF(const State& s,
    const Vector&              f,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBDynamicsCache&      dc = getDynamicsCache(s);

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
            node.calcMInverseFPass1Inward(cc,dc,
                f, allZ, allGepsilon,
                allEpsilon);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcMInverseFPass2Outward(cc,dc, allEpsilon, A_GB, udot);
        }
}


// Must be in ConfigurationStage to calculate qdot = Q*u.
void RigidBodyTree::calcQDot(const State& s, const Vector& u, Vector& qdot) const {
    const SBModelVars&       mv = getModelVars(s);
    const Vector&               q  = getQ(s);
    const SBPositionCache& cc = getPositionCache(s);

    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(mv,q,cc, u, qdot);
}

// Must be in Stage::Velocity to calculate qdotdot = Qdot*u + Q*udot.
void RigidBodyTree::calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
    const SBModelVars&       mv = getModelVars(s);
    const Vector&               q  = getQ(s);
    const SBPositionCache& cc = getPositionCache(s);
    const Vector&               u  = getU(s);

    assert(udot.size() == getTotalDOF());
    qdotdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDotDot(mv,q,cc,u, udot, qdotdot);
}

// If V is a spatial velocity, and you have a X=d(something)/dV (one per body)
// this routine will return d(something)/du for internal generalized speeds u. If
// instead you have d(something)/dR where R is a spatial configuration, this routine
// returns d(something)/dq PROVIDED that dq/dt = u for all q's. That's not true for
// quaternions, so be careful how you use this routine.
// In Kane's terminology, we are calculating the product of a (generalized)
// partial velocity with some vector.
void RigidBodyTree::calcInternalGradientFromSpatial(const State& s, 
                                                    const Vector_<SpatialVec>& X,
                                                    Vector& JX) const
{
    assert(X.size() == getNBodies());

    const SBPositionCache& cc = getPositionCache(s);

    Vector_<SpatialVec> zTemp(getNBodies()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(cc, zTemp, X, JX);
        }
}

// This routine does the same thing as the above but accounts for centrifugal
// forces induced by velocities. The equivalent joint forces returned include
// both the applied forces and the centrifugal ones. Constraints are ignored.
void RigidBodyTree::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const
{
    const SBPositionCache& cc = getPositionCache(s);
    const SBDynamicsCache&      dc = getDynamicsCache(s);

    assert(bodyForces.size() == getNBodies());
    mobilityForces.resize(getTotalDOF());

    Vector_<SpatialVec> allZ(getNBodies());

    // Don't do ground's level since ground has no inboard joint.
    for (int i=rbNodeLevels.size()-1 ; i>0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcEquivalentJointForces(cc,dc,
                bodyForces, allZ,
                mobilityForces);
        }
}

// Pass in a set of internal forces in T; we'll modify them here.
void RigidBodyTree::calcConstraintCorrectedInternalForces(const State& s, Vector& T) {
    lConstraints->projectUVecOntoMotionConstraints(s, T);
}

std::ostream& operator<<(std::ostream& o, const RigidBodyTree& tree) {
    o << "RigidBodyTree has " << tree.getNBodies() << " bodies (incl. G) in "
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

