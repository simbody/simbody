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

void RBStation::calcPosInfo(const SBStateRep& s, RBStationRuntime& rt) const {
    rt.station_G = getNode().getX_GB(s).R() * station_B;
    rt.pos_G     = getNode().getX_GB(s).T() + rt.station_G;
}

void RBStation::calcVelInfo(const SBStateRep& s, RBStationRuntime& rt) const {
    const Vec3& w_G = getNode().getSpatialAngVel(s);
    const Vec3& v_G = getNode().getSpatialLinVel(s);
    rt.stationVel_G = cross(w_G, rt.station_G);
    rt.vel_G = v_G + rt.stationVel_G;
}

void RBStation::calcAccInfo(const SBStateRep& s, RBStationRuntime& rt) const {
    const Vec3& w_G  = getNode().getSpatialAngVel(s);
    const Vec3& v_G  = getNode().getSpatialLinVel(s);
    const Vec3& aa_G = getNode().getSpatialAngAcc(s);
    const Vec3& a_G  = getNode().getSpatialLinAcc(s);
    rt.acc_G = a_G + cross(aa_G, rt.station_G)
                   + cross(w_G, rt.stationVel_G); // i.e., w X (wXr)
}

std::ostream& operator<<(std::ostream& o, const RBStation& s) {
    o << "station " << s.getStation() << " on node " << s.getNode().getNodeNum();
    return o;
}

void RBDistanceConstraint::calcPosInfo(const SBStateRep& s, RBDistanceConstraintRuntime& rt) const
{
    assert(isValid() && runtimeIndex >= 0);
    for (int i=0; i<=1; ++i) stations[i].calcPosInfo(s, rt.stationRuntimes[i]);

    rt.fromTip1ToTip2_G = rt.stationRuntimes[1].pos_G - rt.stationRuntimes[0].pos_G;
    const double separation = rt.fromTip1ToTip2_G.norm();
    rt.unitDirection_G = rt.fromTip1ToTip2_G / separation;
    rt.posErr = distance - separation;
}

void RBDistanceConstraint::calcVelInfo(const SBStateRep& s, RBDistanceConstraintRuntime& rt) const
{
    assert(isValid() && runtimeIndex >= 0);
    for (int i=0; i<=1; ++i) stations[i].calcVelInfo(s, rt.stationRuntimes[i]);

    rt.relVel_G = rt.stationRuntimes[1].vel_G - rt.stationRuntimes[0].vel_G;
    rt.velErr = ~rt.unitDirection_G * rt.relVel_G;
}

void RBDistanceConstraint::calcAccInfo(const SBStateRep& s, RBDistanceConstraintRuntime& rt) const
{
    assert(isValid() && runtimeIndex >= 0);
    for (int i=0; i<=1; ++i) stations[i].calcAccInfo(s, rt.stationRuntimes[i]);

//XXX this doesn't look right
    const Vec3 relAcc_G = rt.stationRuntimes[1].acc_G - rt.stationRuntimes[0].acc_G;
    rt.accErr = rt.relVel_G.normSqr() + (~relAcc_G * rt.fromTip1ToTip2_G);
}

RigidBodyTree::~RigidBodyTree() {
    delete lConstraints; lConstraints=0;

    for (int i=0; i<(int)constraintNodes.size(); ++i)
        delete constraintNodes[i];
    constraintNodes.resize(0);

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
     const TransformMat&     X_PJb,        // parent's frame for attaching this joint
     const TransformMat&     X_BJ,         // inboard joint frame J in body frame
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
        RigidBodyNode::create(MassProperties(), TransformMat(), TransformMat(), 
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
    const RigidBodyNode& parent, const TransformMat& frameInP,
    const RigidBodyNode& child,  const TransformMat& frameInC)
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

// Add a distance constraint and allocate slots to hold the runtime information for
// its stations. Return the assigned distance constraint index for caller's use.
int RigidBodyTree::addOneDistanceConstraintEquation(const RBStation& s1, const RBStation& s2, const double& d)
{
    RBDistanceConstraint dc(s1,s2,d);
    dc.setRuntimeIndex(dcRuntimeInfo.size());
    dcRuntimeInfo.push_back(RBDistanceConstraintRuntime());
    distanceConstraints.push_back(dc);
    return distanceConstraints.size()-1;
}

SBStage RigidBodyTree::getStage(const SBStateRep& s) const {
    return s.getStage(*this);
}

void RigidBodyTree::realize(const SBStateRep& s, SBStage stage) const {
    while (s.getStage(*this) < stage) {
        switch (s.getStage(*this)) {
        case BuiltStage:         realizeModeling(s);     break;
        case ModeledStage:       realizeParameters(s);   break;
        case ParametrizedStage:  realizeTime(s);         break;
        case TimedStage:         realizeConfiguration(s);break;
        case ConfiguredStage:    realizeMotion(s);       break;
        case MovingStage:        realizeDynamics(s);     break;
        case DynamicsStage:      realizeReaction(s);     break;
        default: assert(false);
        }
    }
}

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes for state variables. We construct
// the initialState here also.
void RigidBodyTree::realizeConstruction() {
    if (built) return; // already done

    // Not built yet. Let's count topological things.

    DOFTotal = SqDOFTotal = maxNQTotal = 0;
    for (int i=0; i<(int)rbNodeLevels.size() ; ++i) 
        for (int j=0; j<(int)rbNodeLevels[i].size(); ++j) {
            const int ndof = rbNodeLevels[i][j]->getDOF();
            DOFTotal += ndof; SqDOFTotal += ndof*ndof;
            maxNQTotal += rbNodeLevels[i][j]->getMaxNQ();
        }

    for (int i=0; i<(int)constraintNodes.size(); ++i)
        constraintNodes[i]->finishConstruction(*this);

    lConstraints = new LengthConstraints(*this, 1e-8,0); // TODO: get rid of these numbers
    lConstraints->construct(distanceConstraints, dcRuntimeInfo);
    built = true;

    // Now allocate in initialState the variable we need for the modeling stage,
    // and set them to their defaults.

    assert(initialState.rep == 0);
    initialState.rep = new SBStateRep();
    initialState.rep->setMyHandle(initialState);
    initialState.rep->setStage(*this, BuiltStage);
    initialState.rep->initializeModelingVars(*this);
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
void RigidBodyTree::realizeModeling(const SBStateRep& s) const {
    assert(s.getStage(*this) >= ModeledStage-1);
    if (s.getStage(*this) >= ModeledStage) return;

    s.allocateCacheIfNeeded(*this, ModeledStage);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModeling(s); 

    s.setStage(*this, ModeledStage);
    s.initializeAllVars(*this);
}

// Here we lock in parameterization of the model, such as body masses.
void RigidBodyTree::realizeParameters(const SBStateRep& s) const {
    assert(s.getStage(*this) >= ParametrizedStage-1);
    if (s.getStage(*this) >= ParametrizedStage) return;

    s.allocateCacheIfNeeded(*this, ParametrizedStage);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeParameters(s); 

    s.setStage(*this, ParametrizedStage);
}

void RigidBodyTree::realizeTime(const SBStateRep& s) const {
    assert(s.getStage(*this) >= TimedStage-1);
    if (s.getStage(*this) >= TimedStage) return;

    s.allocateCacheIfNeeded(*this, TimedStage);

    // nothing yet 

    s.setStage(*this, TimedStage);
}

// Set generalized coordinates: sweep from base to tips.
void RigidBodyTree::realizeConfiguration(const SBStateRep& s) const {
    assert(s.getStage(*this) >= ConfiguredStage-1);
    if (s.getStage(*this) >= ConfiguredStage) return;

    s.allocateCacheIfNeeded(*this, ConfiguredStage);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeConfiguration(s); 

    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcPosInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);

    s.setStage(*this, ConfiguredStage);
}

// Set generalized speeds: sweep from base to tip.
// realizeConfiguration() must have been called already.
void RigidBodyTree::realizeMotion(const SBStateRep& s) const {
    assert(s.getStage(*this) >= MovingStage-1);
    if (s.getStage(*this) >= MovingStage) return;

    s.allocateCacheIfNeeded(*this, MovingStage);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeMotion(s); 

    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcVelInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);

    s.setStage(*this, MovingStage);
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P, and velocity-dependent
// quantities like the Coriolis acceleration.

void RigidBodyTree::realizeDynamics(const SBStateRep& s)  const {
    assert(s.getStage(*this) >= DynamicsStage-1);
    if (s.getStage(*this) >= DynamicsStage) return;

    s.allocateCacheIfNeeded(*this, DynamicsStage);

    calcArticulatedBodyInertias(s);
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcJointIndependentDynamicsVel(s);

    s.setStage(*this, DynamicsStage);
}

void RigidBodyTree::realizeReaction(const SBStateRep& s)  const {
    assert(s.getStage(*this) >= ReactingStage-1);
    if (s.getStage(*this) >= ReactingStage) return;

    s.allocateCacheIfNeeded(*this, ReactingStage);

    calcLoopForwardDynamics(s);
    calcQDotDot(s, s.reactionCache.udot, s.reactionCache.qdotdot);

    s.setStage(*this, ReactingStage);
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
    assert(s.getStage(*this) >= BuiltStage);

    // Tree-level defaults
    modelVars.useEulerAngles = false;
    modelVars.prescribed.assign(getNBodies(), false);
    modelVars.prescribed[0] = true; // ground
    modelVars.enabled.assign(getNConstraints(), false);

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultModelingValues(s, modelVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultParameterValues(const SBStateRep& s, 
                                              SBParameterVars& paramVars) const 
{
    assert(s.getStage(*this) >= ModeledStage);

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultParameterValues(s, paramVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultTimeValues(const SBStateRep& s, 
                                         SBTimeVars& timeVars) const 
{
    assert(s.getStage(*this) >= ModeledStage);

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
    assert(s.getStage(*this) >= ModeledStage);

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
    assert(s.getStage(*this) >= ModeledStage);

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
    assert(s.getStage(*this) >= ModeledStage);

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
    assert(s.getStage(*this) >= ModeledStage);

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
    assert(s.getStage(*this) >= ConfiguredStage);

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
    assert(s.getStage(*this) >= ConfiguredStage-1);

    // Fix coordinates first.
    bool anyChange = false;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            if (rbNodeLevels[i][j]->enforceQuaternionConstraints(s))
                anyChange = true;
 
    // Now fix the position constraints produced by defined length constraints.
    if (lConstraints->enforceConfigurationConstraints(s))
        anyChange = true;

    if (anyChange && s.getStage(*this) >= ConfiguredStage)
        s.setStage(*this, SBStage(ConfiguredStage-1));
}

void RigidBodyTree::enforceMotionConstraints(SBStateRep& s) const {
    assert(s.getStage(*this) >= MovingStage-1);

    // Currently there are no coordinate constraints for velocity.

    // Fix the velocity constraints produced by defined length constraints.
    const bool anyChange = lConstraints->enforceMotionConstraints(s);

    if (anyChange && s.getStage(*this) >= MovingStage)
        s.setStage(*this, SBStage(MovingStage-1));
}

// Enforce loop constraints. TODO: OBSOLETE
void RigidBodyTree::enforceLengthConstraints(SBStateRep& s) const {
    Vector& pos = updQ(s);
    Vector& vel = updU(s);
    lConstraints->enforce(s,pos,vel); //FIX: previous constraints still obeyed? (CDS)
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
    assert(s.getStage(*this) >= ReactingStage-1);

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
    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcAccInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);
}

// Given the set of forces in the state, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const SBStateRep& s) const 
{
    assert(s.getStage(*this) >= ReactingStage-1);

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
    assert(s.getStage(*this) >= MovingStage);

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
    assert(s.getStage(*this) >= DynamicsStage);
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
    assert(s.getStage(*this) >= ConfiguredStage);
    assert(u.size() == getTotalDOF());
    qdot.resize(getTotalQAlloc());

    // Skip ground; it doesn't have qdots!
    for (int i=1; i<(int)rbNodeLevels.size(); i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcQDot(s, u, qdot);
}

// Must be in MovingStage to calculate qdotdot = Qdot*u + Q*udot.
void RigidBodyTree::calcQDotDot(const SBStateRep& s, const Vector& udot, Vector& qdotdot) const {
    assert(s.getStage(*this) >= MovingStage);
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
    assert(s.getStage(*this) >= ConfiguredStage);

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
    assert(s.getStage(*this) >= DynamicsStage);
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

