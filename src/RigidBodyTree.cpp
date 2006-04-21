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

RigidBodyTree::RigidBodyTree(const RigidBodyTree& src) {
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
    n->setTreeAndParent(this, &parent);

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
    RBDistanceConstraint* dc = new RBDistanceConstraint(*this,s1,s2,d);
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

    lConstraints = new LengthConstraints(*this, 1e-8,0); // TODO: get rid of these numbers
    lConstraints->construct(distanceConstraints);
    built = true;
}

void RigidBodyTree::realizeConstruction(State& s) const {
    // This is a long-winded way of saying that the Stage must be exactly Allocated.
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Built).prev(), 
        "RigidBodyTree::realizeConstruction()");
    SimTK_STAGECHECK_LT_ALWAYS(s.getStage(), Stage(Stage::Built), 
        "RigidBodyTree::realizeConstruction()");

    // Some of our 'const' values must be treated as mutable *just for this call*.
    // Afterwards they are truly const so we don't declare them mutable, but cheat
    // here instead.
    RigidBodyTree* mutableThis = const_cast<RigidBodyTree*>(this);

    // Fill in the local copy of the constructionCache from the information
    // calculated in endConstruction(). Also ask the State for some room to
    // put Modeling variables & cache and remember the indices in our construction
    // cache.

    mutableThis->constructionCache.nBodies      = nodeNum2NodeMap.size();
    mutableThis->constructionCache.nConstraints = constraintNodes.size();
    mutableThis->constructionCache.nDOFs        = DOFTotal;
    mutableThis->constructionCache.maxNQs       = maxNQTotal;
    mutableThis->constructionCache.sumSqDOFs    = SqDOFTotal;
    mutableThis->constructionCache.nDistanceConstraints = distanceConstraints.size();

    SBModelingVars mvars;
    mvars.allocate(constructionCache);
    setDefaultModelingValues(constructionCache, mvars);
    mutableThis->constructionCache.modelingVarsIndex  = 
        s.allocateDiscreteVariable(Stage::Modeled, new Value<SBModelingVars>(mvars));

    mutableThis->constructionCache.modelingCacheIndex = 
        s.allocateCacheEntry(Stage::Modeled, new Value<SBModelingCache>());

    mutableThis->constructionCache.valid = true;

    // Allocate a cache entry for the constructionCache, and save a copy there.
    mutableThis->constructionCacheIndex = 
        s.allocateCacheEntry(Stage::Built, new Value<SBConstructionCache>(constructionCache));
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
void RigidBodyTree::realizeModeling(State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Modeled).prev(), 
        "RigidBodyTree::realizeModeling()");

    // Get the Modeling-stage cache and make sure it has been allocated and initialized if needed.
    // It is OK to hold a reference here because the discrete variables (and cache entries) in
    // the State are stable, that is, they don't change even if more variables are added.
    SBModelingCache& mCache = updModelingCache(s);
    mCache.allocate(constructionCache);

    // Give the bodies a chance to put something in the cache if they need to.
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModeling(s); 

    // Now allocate all remaining variables and cache entries. We can properly initialize only
    // the next stage, the parameters. Everything else gets some kind of meaningless initial
    // values but those could change as we realize the higher stages.
    SBParameterVars pvars;
    pvars.allocate(constructionCache);
    setDefaultParameterValues(s, pvars);

    mCache.parameterVarsIndex = 
        s.allocateDiscreteVariable(Stage::Parametrized, 
                                   new Value<SBParameterVars>(pvars));
    mCache.parameterCacheIndex = 
        s.allocateCacheEntry(Stage::Parametrized, new Value<SBParameterCache>());

    // No time vars or cache
    mCache.timeVarsIndex = -1;
    mCache.timeCacheIndex = -1;

    // Position variables are just q's, which the State knows how to deal with. We don't know
    // what values are reasonable for q's so we'll set them to NaN here.
    Vector q(maxNQTotal); q.setToNaN();
    setDefaultConfigurationValues(s, q);

    mCache.qIndex = s.allocateQ(q);
    mCache.qVarsIndex = -1; // no config vars other than q
    mCache.qCacheIndex = s.allocateCacheEntry(Stage::Configured, 
        new Value<SBConfigurationCache>());

    // Velocity variables are just the generalized speeds u, which the State knows how to deal
    // with. Zero is always a reasonable value for velocity, so we'll initialize it here.
    Vector u(DOFTotal); u.setToZero();
    setDefaultMotionValues(s, u);

    mCache.uIndex = s.allocateU(u);
    mCache.uVarsIndex = -1; // no velocity vars other than u
    mCache.uCacheIndex = s.allocateCacheEntry(Stage::Moving, 
        new Value<SBMotionCache>());
    // Note that qdots are automatically allocated in the Moving stage cache.

    // no z's or other dynamics vars
    mCache.dynamicsVarsIndex = -1;
    mCache.dynamicsCacheIndex  = s.allocateCacheEntry(Stage::Dynamics, 
        new Value<SBDynamicsCache>());

    // Reaction variables are forces and prescribed accelerations
    SBReactionVars rvars;
    rvars.allocate(constructionCache);
    setDefaultReactionValues(s, rvars);

    mCache.reactionVarsIndex = 
        s.allocateDiscreteVariable(Stage::Reacting, new Value<SBReactionVars>(rvars));
    mCache.reactionCacheIndex = 
        s.allocateCacheEntry(Stage::Reacting, new Value<SBReactionCache>());

    // Note that qdots, qdotdots, udots, zdots are automatically allocated by
    // the State when we advance the stage past modeling.
}

// Here we lock in parameterization of the model, such as body masses.
void RigidBodyTree::realizeParameters(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Parametrized).prev(), 
        "RigidBodyTree::realizeParameters()");

    const SBParameterVars& pVars = getParameterVars(s);

    // Get the Parameter-stage cache and make sure it has been allocated and initialized if needed.
    SBParameterCache& pCache = updParameterCache(s);
    pCache.allocate(constructionCache);

    // Calculate whether we should apply gravity.
    pCache.applyGravity = !(pVars.gravity == Vec3(0.));

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeParameters(s); 
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

    // Get the Configured-stage cache and make sure it has been allocated and initialized if needed.
    SBConfigurationCache& qCache = updConfigurationCache(s);
    qCache.allocate(constructionCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeConfiguration(s); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcPosInfo(s);
}

// Set generalized speeds: sweep from base to tip.
// realizeConfiguration() must have been called already.
void RigidBodyTree::realizeMotion(const State& s) const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Moving).prev(), 
        "RigidBodyTree::realizeMotion()");

    // Get the Motion-stage cache and make sure it has been allocated and initialized if needed.
    SBMotionCache& uCache = updMotionCache(s);
    uCache.allocate(constructionCache);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeMotion(s); 

    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcVelInfo(s);
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.

void RigidBodyTree::realizeDynamics(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Dynamics).prev(), 
        "RigidBodyTree::realizeDynamics()");

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    SBDynamicsCache& dynamicsCache = updDynamicsCache(s);
    dynamicsCache.allocate(constructionCache);

    calcArticulatedBodyInertias(s);
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcJointIndependentDynamicsVel(s);
}

void RigidBodyTree::realizeReaction(const State& s)  const {
    SimTK_STAGECHECK_GE_ALWAYS(s.getStage(), Stage(Stage::Reacting).prev(), 
        "RigidBodyTree::realizeReaction()");

    // We need the Modeling results to tell us where to find things.
    const SBModelingCache& mCache = getModelingCache(s);

    // Get the Dynamics-stage cache and make sure it has been allocated and initialized if needed.
    Vector& udot    = s.updUDot()(mCache.uIndex, DOFTotal);
    Vector& qdotdot = s.updQDotDot()(mCache.qIndex, maxNQTotal);
    SBReactionCache& reactionCache = updReactionCache(s);
    reactionCache.allocate(constructionCache);

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
void RigidBodyTree::setDefaultModelingValues(const SBConstructionCache& constructionCache, 
                                             SBModelingVars& modelVars) const 
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
            rbNodeLevels[i][j]->setDefaultModelingValues(constructionCache, modelVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultParameterValues(const State& s, 
                                              SBParameterVars& paramVars) const 
{
    // Tree-level defaults
    paramVars.gravity = 0.;

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultParameterValues(s, paramVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultTimeValues(const State& s, 
                                         SBTimeVars& timeVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultTimeValues(s, timeVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultConfigurationValues(const State& s, Vector& q) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultConfigurationValues(s, q);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultMotionValues(const State& s, Vector& u) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultMotionValues(s, u);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultDynamicsValues(const State& s, 
                                             SBDynamicsVars& dynamicsVars) const 
{
    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultDynamicsValues(s, dynamicsVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultReactionValues(const State& s, 
                                             SBReactionVars& reactionVars) const 
{
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

void RigidBodyTree::setUseEulerAngles(State& s, bool useAngles) const {
    SBModelingVars& modelVars = updModelingVars(s); // check/adjust stage
    modelVars.useEulerAngles = useAngles;
}
void RigidBodyTree::setJointIsPrescribed(State& s, int joint, bool prescribe) const {
    SBModelingVars& modelVars = updModelingVars(s); // check/adjust stage
    modelVars.prescribed[joint] = prescribe;
}
void RigidBodyTree::setConstraintIsEnabled(State& s, int constraint, bool enable) const {
    SBModelingVars& modelVars = updModelingVars(s); // check/adjust stage
    modelVars.enabled[constraint] = enable;   
}

bool RigidBodyTree::getUseEulerAngles(const State& s) const {
    const SBModelingVars& modelVars = getModelingVars(s); // check stage
    return modelVars.useEulerAngles;
}
bool RigidBodyTree::isJointPrescribed(const State& s, int joint) const {
    const SBModelingVars& modelVars = getModelingVars(s); // check stage
    return modelVars.prescribed[joint];
}
bool RigidBodyTree::isConstraintEnabled(const State& s, int constraint) const {
    const SBModelingVars& modelVars = getModelingVars(s); // check stage
    return modelVars.enabled[constraint];
}

const Real& RigidBodyTree::getJointQ(const State& s, int body, int axis) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getNQ(s));
    return getQ(s)[n.getQIndex()+axis];
}


const Real& RigidBodyTree::getJointU(const State& s, int body, int axis) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    return getU(s)[n.getUIndex()+axis];
}

void RigidBodyTree::setJointQ(State& s, int body, int axis, const Real& r) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getNQ(s));
    updQ(s)[n.getQIndex()+axis] = r;
}


void RigidBodyTree::setJointU(State& s, int body, int axis, const Real& r) const {
    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    updU(s)[n.getUIndex()+axis] = r;
}


void RigidBodyTree::setPrescribedUdot(State& s, int body, int axis, const Real& r) const {
    SBReactionVars& reactionVars = updReactionVars(s); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    reactionVars.prescribedUdot[n.getUIndex()+axis] = r;
}

void RigidBodyTree::clearAppliedForces(State& s) const {
    SBReactionVars& reactionVars = updReactionVars(s); // check/adjust stage
    reactionVars.appliedJointForces.setToZero();
    reactionVars.appliedBodyForces.setToZero();
}

void RigidBodyTree::applyGravity(State& s, const Vec3& g) const {
    assert(s.getStage() >= Stage::Configured);

    for (int body=1; body<getNBodies(); ++body) {
        const RigidBodyNode& n = getRigidBodyNode(body);
        applyPointForce(s, body, n.getCOM_B(s), n.getMass(s)*g);
    }
}

void RigidBodyTree::applyPointForce(State& s, int body, const Vec3& stationInB, 
                        const Vec3& forceInG) const
{
    SBReactionVars& reactionVars = updReactionVars(s); // check/adjust stage

    const RotationMat& R_GB = getRigidBodyNode(body).getX_GB(s).R();
    reactionVars.appliedBodyForces[body] += 
        SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void RigidBodyTree::applyBodyTorque(State& s, int body, const Vec3& torqueInG) const {
    SBReactionVars& reactionVars = updReactionVars(s); // check/adjust stage

    reactionVars.appliedBodyForces[body] += SpatialVec(torqueInG, Vec3(0));
}

void RigidBodyTree::applyJointForce(State& s, int body, int axis, const Real& r) const {
    SBReactionVars& reactionVars = updReactionVars(s); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    reactionVars.appliedJointForces[n.getUIndex()+axis] = r;
}

const Vector& 
RigidBodyTree::getAppliedJointForces(const State& s) const {
    const SBReactionVars& reactionVars = getReactionVars(s);
    return reactionVars.appliedJointForces;
}
const Vector_<SpatialVec>& 
RigidBodyTree::getAppliedBodyForces(const State& s) const {
    const SBReactionVars& reactionVars = getReactionVars(s);
    return reactionVars.appliedBodyForces;
}

const Transform&
RigidBodyTree::getBodyConfiguration(const State& s, int body) const
  { return getRigidBodyNode(body).getX_GB(s); }

const SpatialVec&
RigidBodyTree::getBodyVelocity(const State& s, int body) const
  { return getRigidBodyNode(body).getV_GB(s); }

const SpatialVec&
RigidBodyTree::getBodyAcceleration(const State& s, int body) const
  { return getRigidBodyNode(body).getA_GB(s); }

void RigidBodyTree::enforceConfigurationConstraints(State& s) const {
    assert(s.getStage() >= Stage::Configured-1);

    // Fix coordinates first.
    bool anyChange = false;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            if (rbNodeLevels[i][j]->enforceQuaternionConstraints(s))
                anyChange = true;
    //TODO: quaternion constraints shouldn't invalidate anything except
    // the qnorms, which will be all 1 now
 
    // Now fix the position constraints produced by defined length constraints.
    if (lConstraints->enforceConfigurationConstraints(s))
        anyChange = true;

    if (anyChange)
        s.invalidateStage(Stage::Configured);
}

void RigidBodyTree::enforceMotionConstraints(State& s) const {
    assert(s.getStage() >= Stage::Moving-1);

    // Currently there are no coordinate constraints for velocity.

    // Fix the velocity constraints produced by defined length constraints.
    const bool anyChange = lConstraints->enforceMotionConstraints(s);

    if (anyChange)
        s.invalidateStage(Stage::Moving);
}

// Given a forces in the state, calculate accelerations ignoring
// constraints, and leave the results in the state. 
// Must have already called realizeDynamics().
// We also allow some extra forces to be supplied, with the intent
// that these will be used to deal with internal forces generated
// by constraints. 
void RigidBodyTree::calcTreeForwardDynamics (const State& s,
     const Vector*              extraJointForces,
     const Vector_<SpatialVec>* extraBodyForces) const
{
    assert(s.getStage() >= Stage::Reacting-1);

    Vector              totalJointForces;
    Vector_<SpatialVec> totalBodyForces;

    // inputs
    const Vector& jointForces = getReactionVars(s).appliedJointForces;
    const Vector_<SpatialVec>&
                  bodyForces  = getReactionVars(s).appliedBodyForces;

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
    Vector& netHingeForces  = updReactionCache(s).netHingeForces;
    Vector& udot            = updUDot(s);
    Vector_<SpatialVec>& 
            A_GB            = updReactionCache(s).bodyAccelerationInGround;

    calcTreeAccelerations(s, *jointForcesToUse, *bodyForcesToUse,
                          netHingeForces, A_GB, udot);
    
    // Calculate constraint acceleration errors.
    for (int i=0; i < (int)distanceConstraints.size(); ++i)
        distanceConstraints[i]->calcAccInfo(s);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const State& s) const 
{
    assert(s.getStage() >= Stage::Reacting-1);

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
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcArticulatedBodyInertiasInward(s);
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcZ(const State& s, 
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
void RigidBodyTree::calcY(const State& s) const {
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcYOutward(s);
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcTreeAccel(const State& s) const {
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel(s);
}

void RigidBodyTree::fixVel0(State& s, Vector& vel) const {
    lConstraints->fixVel0(s, vel);
}

Real RigidBodyTree::calcKineticEnergy(const State& s) const {
    assert(s.getStage() >= Stage::Moving);

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
void RigidBodyTree::calcTreeAccelerations(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    netHingeForces,
    Vector_<SpatialVec>&       A_GB,
    Vector&                    udot) const 
{
    assert(s.getStage() >= Stage::Dynamics);
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
void RigidBodyTree::calcQDot(const State& s, const Vector& u, Vector& qdot) const {
    assert(s.getStage() >= Stage::Configured);
    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(s, u, qdot);
}

// Must be in Stage::Moving to calculate qdotdot = Qdot*u + Q*udot.
void RigidBodyTree::calcQDotDot(const State& s, const Vector& udot, Vector& qdotdot) const {
    assert(s.getStage() >= Stage::Moving);
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
void RigidBodyTree::calcInternalGradientFromSpatial(const State& s, 
                                                    const Vector_<SpatialVec>& X,
                                                    Vector& JX) 
{
    assert(X.size() == getNBodies());
    assert(s.getStage() >= Stage::Configured);

    Vector_<SpatialVec> zTemp(getNBodies()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(s, zTemp, X, JX);
        }
}

void RigidBodyTree::calcTreeEquivalentJointForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    jointForces)
{
    assert(s.getStage() >= Stage::Dynamics);
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
void RigidBodyTree::calcConstraintCorrectedInternalForces(const State& s, Vector& T) {
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

