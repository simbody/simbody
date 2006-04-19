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

#include <string>


SBState::~SBState() {
    delete rep;
}

SBState::SBState(const SBState& src) : rep(0) {
    if (src.rep) {
        rep = new SBStateRep(*src.rep);
        rep->setMyHandle(*this);
    }
}

SBState& SBState::operator=(const SBState& src) {
    if (&src != this) {
        delete rep; rep = 0;
        if (src.rep) {
            rep = new SBStateRep(*src.rep);
            rep->setMyHandle(*this);
        }
    }
    return *this;
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
    (RigidBodyNode&          parent,
     const MassProperties&   m,            // mass properties in body frame
     const Transform&     X_PJb,        // parent's frame for attaching this joint
     const Transform&     X_BJ,         // inboard joint frame J in body frame
     JointSpecification::JointType        
                             type,
     bool                    isReversed,   // child-to-parent orientation?
     int&                    nxtU,
     int&                    nxtUSq,
     int&                    nxtQ)
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
                              JointSpecification::ThisIsGround,
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
void RigidBodyTree::realizeConstruction(State& s) {
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

    lConstraints = new LengthConstraints(*this, 1e-8,0); // TODO: get rid of these numbers
    lConstraints->construct(distanceConstraints);
    built = true;

    // Now allocate in the State the variables we need for the Modeling stage,
    // and set them to their defaults.

    modelingVarsIndex  = s.allocateDiscreteVariable(Stage::Modeled, Value<SBModelingVars>());
    modelingCacheIndex = s.allocateCacheEntry      (Stage::Modeled, Value<SBModelingCache>());

    setDefaultModelingValues(Value<SBModelingVars>::downcast(s.updDiscreteVariable(modelingVarsIndex)));
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
void RigidBodyTree::realizeModeling(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Modeled).prev(), 
        "RigidBodyTree::realizeModeling()");

    // Get the Modeling-stage cache and make sure it has been allocated and initialized if needed.
    SBModelingCache& mCache = updModelingCache(s);
    mCache.allocate();

    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModeling(s, mCache); 

    // Now allocate all remaining variables and cache entries. We can properly initialize only
    // the next stage, the parameters. Everything else gets some kind of meaningless initial
    // values but those could change as we realize the higher stages.
    SBParameterVars pvars;
    initializeParameters(pvars);
    mCache.parametersIndex = 
        s.allocateDiscreteVariable(Stage::Parametrized, Value<SBParameterVars>(pvars));
    mCache.parameterCacheIndex = s.allocateCacheEntry(Stage::Parametrized, Value<SBParameterCache>());

    // No time vars or cache

    // Position variables are just q's, which the State knows how to deal with. We don't know
    // what values are reasonable for q's so we'll set them to NaN here.
    Vector q(maxNQTotal); q.setToNaN();
    mCache.qIndex = s.allocateQRange(q);
    mCache.qCacheIndex = s.allocateCacheEntry(Stage::Configured, Value<SBConfigurationCache>());

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.
    Vector u(DOFTotal); u.setToZero();
    mCache.uIndex = s.allocateURange(u);
    mCache.uCacheIndex = s.allocateCacheEntry(Stage::Moving, Value<SBMotionCache>());
    // Note that qdots are automatically allocated in the Moving stage cache.

    // no z's or other dynamics vars
    mCache.dynamicCacheIndex  = s.allocateCacheEntry(Stage::Dynamics, Value<SBDynamicsCache>());

    // Reaction variables are forces and prescribed accelerations
    SBReactionVars rvars;
    initializeReactionVars(rvars);
    mCache.reactionIndex = 
        s.allocateDiscreteVariable(Stage::Reacting, Value<SBReactionVars>(rvars));
    mCache.reactionCacheIndex = s.allocateCacheEntry(Stage::Reacting, Value<SBReactionCache>());
    // Note that qdotdots, udots, zdots are automatically allocated in the Reacting stage cache.
}

// Here we lock in parameterization of the model, such as body masses.
void RigidBodyTree::realizeParameters(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Parametrized).prev(), 
        "RigidBodyTree::realizeParameters()");

    // We need the Modeling results to tell us where to find the parameters.
    const SBModelingCache& mCache = getModelingCache(s);

    const SBParameterVars& pVars =
        Value<SBParameterVars>::downcast(s.getDiscreteVariable(mCache.parametersIndex)).get();

    // Get the Parameter-stage cache and make sure it has been allocated and initialized if needed.
    SBParameterCache& pCache = 
        Value<SBParameterCache>::downcast(s.updCacheEntry(mCache.parameterCacheIndex)).upd();
    pCache.allocate();

    // Calculate whether we should apply gravity.
    pCache.applyGravity = !(pVars.gravity == Vec3(0.));

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeParameters(s, pCache); 
}

void RigidBodyTree::realizeTime(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Timed).prev(), 
        "RigidBodyTree::realizeTime()");

    // nothing yet
}

// Set generalized coordinates: sweep from base to tips.
void RigidBodyTree::realizeConfiguration(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Configured).prev(), 
        "RigidBodyTree::realizeConfiguration()");

    // We need the Modeling results to tell us where to find things.
    const SBModelingCache& mCache = getModelingCache(s);

    // TODO: does Vector& work here? Would have to cast from VectorView without copy.
    const VectorView q = s.getQ()(mCache.qIndex, maxNQTotal);

    // Get the Configured-stage cache and make sure it has been allocated and initialized if needed.
    SBConfigurationCache& qCache = 
        Value<SBConfigurationCache>::downcast(s.updCacheEntry(mCache.qCacheIndex)).upd();
    qCache.allocate();

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeConfiguration(s,q,qCache); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(s,q,qCache);
}

// Set generalized speeds: sweep from base to tip.
// realizeConfiguration() must have been called already.
void RigidBodyTree::realizeMotion(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Moving).prev(), 
        "RigidBodyTree::realizeMotion()");

    // We need the Modeling results to tell us where to find things.
    const SBModelingCache& mCache = getModelingCache(s);

    const VectorView u = s.getU()(mCache.uIndex, DOFTotal);

    // Get the Motion-stage cache and make sure it has been allocated and initialized if needed.
    VectorView qdot = s.updQDot()(mCache.qIndex, maxNQTotal);
    SBMotionCache& uCache = 
        Value<SBMotionCache>::downcast(s.updCacheEntry(mCache.uCacheIndex)).upd();
    uCache.allocate();

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeMotion(s,u,uCache); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcVelInfo(s,u,uCache);
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.

void RigidBodyTree::realizeDynamics(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Dynamics).prev(), 
        "RigidBodyTree::realizeDynamics()");

    // We need the Modeling results to tell us where to find things.
    const SBModelingCache& mCache = getModelingCache(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    SBDynamicsCache& dynamicsCache = 
        Value<SBDynamicsCache>::downcast(s.updCacheEntry(mCache.dynamicsCacheIndex)).upd();
    dynamicsCache.allocate();

    calcArticulatedBodyInertias(s, dynamicsCache);
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcJointIndependentDynamicsVel(s, dynamicsCache);
}

void RigidBodyTree::realizeReaction(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Reaction).prev(), 
        "RigidBodyTree::realizeReaction()");

    // We need the Modeling results to tell us where to find things.
    const SBModelingCache& mCache = getModelingCache(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    VectorView udot = s.updUDot()(mCache.uIndex, DOFTotal);
    VectorView qdotdot = s.updQDotDot()(mCache.qIndex, maxNQTotal);
    SBReactionCache& reactionCache = 
        Value<SBReactionCache>::downcast(s.updCacheEntry(mCache.reactionCacheIndex)).upd();
    reactionCache.allocate();

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

void RigidBodyTree::setDefaultModelingValues(const SBStateRep& s, 
                                             SBModelingVars& modelVars) const 
{
    assert(s.getStage(*this) >= Stage::Built);

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
            rbNodeLevels[i][j]->setDefaultModelingValues(s, modelVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultParameterValues(const SBStateRep& s, 
                                              SBParameterVars& paramVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults
    paramVars.gravity = 0.;

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultParameterValues(s, paramVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultTimeValues(const SBStateRep& s, 
                                         SBTimeVars& timeVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultTimeValues(s, timeVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultConfigurationValues(const SBStateRep& s, 
                                                  SBConfigurationVars& configVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultConfigurationValues(s, configVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultMotionValues(const SBStateRep& s, 
                                           SBMotionVars& motionVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultMotionValues(s, motionVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultDynamicsValues(const SBStateRep& s, 
                                             SBDynamicsVars& dynamicsVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultDynamicsValues(s, dynamicsVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultReactionValues(const SBStateRep& s, 
                                             SBReactionVars& reactionVars) const 
{
    assert(s.getStage(*this) >= Stage::Modeled);

    // Tree-level defaults
    reactionVars.appliedJointForces.setToZero();
    reactionVars.appliedBodyForces.setToZero();
    reactionVars.prescribedUdot.setToZero();

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultReactionValues(s, reactionVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setUseEulerAngles(SBStateRep& s, bool useAngles) const {
    SBModelingVars& modelVars = s.updModelingVars(*this); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void RigidBodyTree::setJointIsPrescribed(SBStateRep& s, int joint, bool prescribe) const {
    SBModelingVars& modelVars = s.updModelingVars(*this); // check/adjust stage
    modelVars.prescribed[joint] = prescribe;
}
void RigidBodyTree::setConstraintIsEnabled(SBStateRep& s, int constraint, bool enable) const {
    SBModelingVars& modelVars = s.updModelingVars(*this); // check/adjust stage
    modelVars.enabled[constraint] = enable;   
}

bool RigidBodyTree::getUseEulerAngles(const SBStateRep& s) const {
    const SBModelingVars& modelVars = s.getModelingVars(*this); // check stage
    return modelVars.useEulerAngles;
}
bool RigidBodyTree::isJointPrescribed(const SBStateRep& s, int joint) const {
    const SBModelingVars& modelVars = s.getModelingVars(*this); // check stage
    return modelVars.prescribed[joint];
}
bool RigidBodyTree::isConstraintEnabled(const SBStateRep& s, int constraint) const {
    const SBModelingVars& modelVars = s.getModelingVars(*this); // check stage
    return modelVars.enabled[constraint];
}

void RigidBodyTree::setJointQ(SBStateRep& s, int body, int axis, const Real& r) const {
    SBConfigurationVars& configVars = s.updConfigurationVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getNQ(s));
    configVars.q[n.getQIndex()+axis] = r;
}


void RigidBodyTree::setJointU(SBStateRep& s, int body, int axis, const Real& r) const {
    SBMotionVars& motionVars = s.updMotionVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    motionVars.u[n.getUIndex()+axis] = r;
}


void RigidBodyTree::setPrescribedUdot(SBStateRep& s, int body, int axis, const Real& r) const {
    SBReactionVars& reactionVars = s.updReactionVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    reactionVars.prescribedUdot[n.getUIndex()+axis] = r;
}

void RigidBodyTree::clearAppliedForces(SBStateRep& s) const {
    SBReactionVars& reactionVars = s.updReactionVars(*this); // check/adjust stage
    reactionVars.appliedJointForces.setToZero();
    reactionVars.appliedBodyForces.setToZero();
}

void RigidBodyTree::applyGravity(SBStateRep& s, const Vec3& g) const {
    assert(s.getStage(*this) >= Stage::Configured);

    for (int body=1; body<getNBodies(); ++body) {
        const RigidBodyNode& n = getRigidBodyNode(body);
        applyPointForce(s, body, n.getCOM_B(s), n.getMass(s)*g);
    }
}

void RigidBodyTree::applyPointForce(SBStateRep& s, int body, const Vec3& stationInB, 
                        const Vec3& forceInG) const
{
    SBReactionVars& reactionVars = s.updReactionVars(*this); // check/adjust stage

    const RotationMat& R_GB = getRigidBodyNode(body).getX_GB(s).R();
    reactionVars.appliedBodyForces[body] += 
        SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void RigidBodyTree::applyBodyTorque(SBStateRep& s, int body, const Vec3& torqueInG) const {
    SBReactionVars& reactionVars = s.updReactionVars(*this); // check/adjust stage

    reactionVars.appliedBodyForces[body] += SpatialVec(torqueInG, Vec3(0));
}

void RigidBodyTree::applyJointForce(SBStateRep& s, int body, int axis, const Real& r) const {
    SBReactionVars& reactionVars = s.updReactionVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    reactionVars.appliedJointForces[n.getUIndex()+axis] = r;
}

// Access to continuous state variables and their derivatives.
void RigidBodyTree::setQ(SBStateRep& s, const Vector& q) const {
    assert(q.size() == getTotalQAlloc());

    SBConfigurationVars& configVars = s.updConfigurationVars(*this);
    configVars.q = q;
}

void RigidBodyTree::setU(SBStateRep& s, const Vector& u) const {
    assert(u.size() == getTotalDOF());
    SBMotionVars& motionVars = s.updMotionVars(*this);
    motionVars.u = u;
}

Vector& RigidBodyTree::updQ(SBStateRep& s) const {
    SBConfigurationVars& configVars = s.updConfigurationVars(*this);
    return configVars.q;
}

Vector& RigidBodyTree::updU(SBStateRep& s) const {
    SBMotionVars& motionVars = s.updMotionVars(*this);
    return motionVars.u;
}

const Vector& RigidBodyTree::getQ(const SBStateRep& s) const {
    const SBConfigurationVars& configVars = s.getConfigurationVars(*this);
    return configVars.q;
}

const Vector& RigidBodyTree::getU(const SBStateRep& s) const {
    const SBMotionVars& motionVars = s.getMotionVars(*this);
    return motionVars.u;
}

const Vector& 
RigidBodyTree::getAppliedJointForces(const SBStateRep& s) const {
    const SBReactionVars& reactionVars = s.getReactionVars(*this);
    return reactionVars.appliedJointForces;
}
const Vector_<SpatialVec>& 
RigidBodyTree::getAppliedBodyForces(const SBStateRep& s) const {
    const SBReactionVars& reactionVars = s.getReactionVars(*this);
    return reactionVars.appliedBodyForces;
}

const Vector& RigidBodyTree::getQDot(const SBStateRep& s) const {
    const SBMotionCache& motionCache = s.getMotionCache(*this);
    return motionCache.qdot;
}
const Vector& RigidBodyTree::getUDot(const SBStateRep& s) const {
    const SBReactionCache& reactionCache = s.getReactionCache(*this);
    return reactionCache.udot;
}
const Vector& RigidBodyTree::getQDotDot(const SBStateRep& s) const {
    const SBReactionCache& reactionCache = s.getReactionCache(*this);
    return reactionCache.qdotdot;
}

void RigidBodyTree::enforceConfigurationConstraints(SBStateRep& s) const {
    assert(s.getStage(*this) >= Stage::Configured-1);

    // Fix coordinates first.
    bool anyChange = false;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            if (rbNodeLevels[i][j]->enforceQuaternionConstraints(s))
                anyChange = true;
 
    // Now fix the position constraints produced by defined length constraints.
    if (lConstraints->enforceConfigurationConstraints(s))
        anyChange = true;

    if (anyChange && s.getStage(*this) >= Stage::Configured)
        s.setStage(*this, Stage(Stage::Configured).prev());
}

void RigidBodyTree::enforceMotionConstraints(SBStateRep& s) const {
    assert(s.getStage(*this) >= Stage::Moving-1);

    // Currently there are no coordinate constraints for velocity.

    // Fix the velocity constraints produced by defined length constraints.
    const bool anyChange = lConstraints->enforceMotionConstraints(s);

    if (anyChange && s.getStage(*this) >= Stage::Moving)
        s.setStage(*this, Stage(Stage::Moving).prev());
}

// Given a forces in the state, calculate accelerations ignoring
// constraints, and leave the results in the state. 
// Must have already called realizeDynamics().
// We also allow some extra forces to be supplied, with the intent
// that these will be used to deal with internal forces generated
// by constraints. 
void RigidBodyTree::calcTreeForwardDynamics (const SBStateRep& s,
     const Vector*              extraJointForces,
     const Vector_<SpatialVec>* extraBodyForces) const
{
    assert(s.getStage(*this) >= Stage::Reacting-1);

    Vector              totalJointForces;
    Vector_<SpatialVec> totalBodyForces;

    // inputs
    const Vector& jointForces = s.reactionVars.appliedJointForces;
    const Vector_<SpatialVec>&
                  bodyForces  = s.reactionVars.appliedBodyForces;

    const Vector*              jointForcesToUse = &jointForces;
    const Vector_<SpatialVec>* bodyForcesToUse = &bodyForces;

    if (extraJointForces) {
        totalJointForces = jointForces + *extraJointForces;
        jointForcesToUse = &totalJointForces;
    }

    if (extraBodyForces) {
        totalBodyForces = bodyForces + *extraBodyForces;
        bodyForcesToUse = &totalBodyForces;
    }

    // outputs
    Vector& netHingeForces  = s.reactionCache.netHingeForces;
    Vector& udot            = s.reactionCache.udot;
    Vector_<SpatialVec>& 
            A_GB            = s.reactionCache.bodyAccelerationInGround;

    calcTreeAccelerations(s, *jointForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot);
    
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(s);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const SBStateRep& s) const 
{
    assert(s.getStage(*this) >= Stage::Reacting-1);

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
void RigidBodyTree::calcArticulatedBodyInertias(const SBStateRep& s) const {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcArticulatedBodyInertiasInward(s);
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcZ(const SBStateRep& s, 
                          const SpatialVecList& spatialForces) const
{
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(s, spatialForces[node.getNodeNum()]);
        }
}

// Y is used for length constraints: sweep from base to tip.
void RigidBodyTree::calcY(const SBStateRep& s) const {
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcYOutward(s);
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcTreeAccel(const SBStateRep& s) const {
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel(s);
}

void RigidBodyTree::fixVel0(SBStateRep& s, Vector& vel) const {
    lConstraints->fixVel0(s, vel);
}

Real RigidBodyTree::calcKineticEnergy(const SBStateRep& s) const {
    assert(s.getStage(*this) >= Stage::Moving);

    Real ke = 0.;

    // Skip ground level 0!
    for (int i=1 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            ke += rbNodeLevels[i][j]->calcKineticEnergy(s);

    return ke;
}

//
// Operator for open-loop dynamics.
//
void RigidBodyTree::calcTreeAccelerations(const SBStateRep& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    assert(s.getStage(*this) >= Stage::Dynamics);
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
            node.calcUDotPass1Inward(s,
                jointForces, bodyForces, allZ, allGepsilon,
                netHingeForces);
        }

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcUDotPass2Outward(s, netHingeForces, A_GB, udot);
        }
}


// Must be in ConfigurationStage to calculate qdot = Q*u.
void RigidBodyTree::calcQDot(const SBStateRep& s, const Vector& u, Vector& qdot) const {
    assert(s.getStage(*this) >= Stage::Configured);
    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(s, u, qdot);
}

// Must be in Stage::Moving to calculate qdotdot = Qdot*u + Q*udot.
void RigidBodyTree::calcQDotDot(const SBStateRep& s, const Vector& udot, Vector& qdotdot) const {
    assert(s.getStage(*this) >= Stage::Moving);
    assert(udot.size() == getTotalDOF());
    qdotdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDotDot(s, udot, qdotdot);
}

// If V is a spatial velocity, and you have a X=d(something)/dV (one per body)
// this routine will return d(something)/du for internal generalized speeds u. If
// instead you have d(something)/dR where R is a spatial configuration, this routine
// returns d(something)/dq PROVIDED that dq/dt = u for all q's. That's not true for
// quaternions, so be careful how you use this routine.
// In Kane's terminology, we are calculating the product of a (generalized)
// partial velocity with some vector.
void RigidBodyTree::calcInternalGradientFromSpatial(const SBStateRep& s, 
                                                    const Vector_<SpatialVec>& X,
                                                    Vector& JX) 
{
    assert(X.size() == getNBodies());
    assert(s.getStage(*this) >= Stage::Configured);

    Vector_<SpatialVec> zTemp(getNBodies()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(s, zTemp, X, JX);
        }
}

void RigidBodyTree::calcTreeEquivalentJointForces(const SBStateRep& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    jointForces)
{
    assert(s.getStage(*this) >= Stage::Dynamics);
    assert(bodyForces.size() == getNBodies());
    jointForces.resize(getTotalDOF());

    Vector_<SpatialVec> allZ(getNBodies());
    Vector_<SpatialVec> allGepsilon(getNBodies());

    // Don't do ground's level since ground has no inboard joint.
    for (int i=rbNodeLevels.size()-1 ; i>0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcEquivalentJointForces(s,
                bodyForces, allZ, allGepsilon,
                jointForces);
        }
}

// Pass in a set of internal forces in T; we'll modify them here.
void RigidBodyTree::calcConstraintCorrectedInternalForces(const SBStateRep& s, Vector& T) {
    lConstraints->fixGradient(s, T);
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

