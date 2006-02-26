/**@file
 *
 * Implementation of RigidBodyTree.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "RigidBodyTree.h"
#include "SimbodyTreeState.h"
#include "RigidBodyNode.h"
#include "LengthConstraints.h"

#include <string>

SBState::~SBState() {
    delete rep;
}

SBState::SBState(const SBState& src) : rep(0) {
    if (src.rep) rep = new SBStateRep(*src.rep);
}

SBState& SBState::operator=(const SBState& src) {
    if (&src != this) {
        delete rep; rep = 0;
        if (src.rep) rep = new SBStateRep(*src.rep);
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
    delete lConstraints;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) {
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
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

// Add a distance constraint and allocate slots to hold the runtime information for
// its stations. Return the assigned distance constraint index for caller's use.
int RigidBodyTree::addDistanceConstraint(const RBStation& s1, const RBStation& s2, const double& d)
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
        case MovingStage:        realizeReaction(s);     break;
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
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const int ndof = rbNodeLevels[i][j]->getDOF();
            DOFTotal += ndof; SqDOFTotal += ndof*ndof;
            maxNQTotal += rbNodeLevels[i][j]->getMaxNQ();
        }

    lConstraints = new LengthConstraints(*this, 1e-6,0); // TODO: get rid of these numbers
    lConstraints->construct(distanceConstraints, dcRuntimeInfo);
    built = true;

    // Now allocate in initialState the variable we need for the modeling stage,
    // and set them to their defaults.

    assert(initialState.rep == 0);
    initialState.rep = new SBStateRep();
    initialState.rep->setStage(*this, BuiltStage);
    initialState.rep->allocateVarsIfNeeded(*this, ModeledStage);
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
    s.allocateVarsIfNeeded(*this, SBStage(ModeledStage + 1));
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
    s.allocateVarsIfNeeded(*this, SBStage(ParametrizedStage + 1));
}

void RigidBodyTree::realizeTime(const SBStateRep& s) const {
    assert(s.getStage(*this) >= TimedStage-1);
    if (s.getStage(*this) >= TimedStage) return;

    s.allocateCacheIfNeeded(*this, TimedStage);

    // nothing yet 

    s.setStage(*this, TimedStage);
    s.allocateVarsIfNeeded(*this, SBStage(TimedStage + 1));
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
    s.allocateVarsIfNeeded(*this, SBStage(ConfiguredStage + 1));
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
    s.allocateVarsIfNeeded(*this, SBStage(MovingStage + 1));
}

void RigidBodyTree::realizeReaction(const SBStateRep& s)  const {
    assert(s.getStage(*this) >= ReactingStage-1);
    if (s.getStage(*this) >= ReactingStage) return;

    s.allocateCacheIfNeeded(*this, ReactingStage);

    prepareForDynamics(s);
    calcLoopForwardDynamics(s, s.dynamicVars.appliedBodyForces);

    s.setStage(*this, ReactingStage);
    // no more variables to allocate
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
    assert(s.getStage(*this) >= ModeledStage-1);
    s.setStage(*this, SBStage(ModeledStage-1)); // back up if necessary

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
    assert(s.getStage(*this) >= ParametrizedStage-1);
    s.setStage(*this, SBStage(ParametrizedStage-1)); // back up if necessary

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
    assert(s.getStage(*this) >= TimedStage-1);
    s.setStage(*this, SBStage(TimedStage-1)); // back up if necessary

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
    assert(s.getStage(*this) >= ConfiguredStage-1);
    s.setStage(*this, SBStage(ConfiguredStage-1)); // back up if necessary

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
    assert(s.getStage(*this) >= MovingStage-1);
    s.setStage(*this, SBStage(MovingStage-1)); // back up if necessary

    // Tree-level defaults (none)

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultMotionValues(s, motionVars);

    // TODO: constraint defaults
}

void RigidBodyTree::setDefaultDynamicValues(const SBStateRep& s, 
                                            SBDynamicVars& dynamicVars) const 
{
    assert(s.getStage(*this) >= ReactingStage-1);
    s.setStage(*this, SBStage(ReactingStage-1)); // back up if necessary

    // Tree-level defaults
    dynamicVars.appliedJointForces.setToZero();
    dynamicVars.appliedBodyForces.setToZero();
    dynamicVars.prescribedUdot.setToZero();

    // Node/joint-level defaults
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->setDefaultDynamicValues(s, dynamicVars);

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
    SBDynamicVars& dynamicVars = s.updDynamicVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    dynamicVars.prescribedUdot[n.getUIndex()+axis] = r;
}

void RigidBodyTree::clearAppliedForces(SBStateRep& s) const {
    SBDynamicVars& dynamicVars = s.updDynamicVars(*this); // check/adjust stage
    dynamicVars.appliedJointForces.setToZero();
    dynamicVars.appliedBodyForces.setToZero();
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
    SBDynamicVars& dynamicVars = s.updDynamicVars(*this); // check/adjust stage

    const RotationMat& R_GB = getRigidBodyNode(body).getX_GB(s).R();
    dynamicVars.appliedBodyForces[body] += 
        SpatialVec((R_GB*stationInB) % forceInG, forceInG);
}

void RigidBodyTree::applyBodyTorque(SBStateRep& s, int body, const Vec3& torqueInG) const {
    SBDynamicVars& dynamicVars = s.updDynamicVars(*this); // check/adjust stage

    dynamicVars.appliedBodyForces[body] += SpatialVec(torqueInG, Vec3(0));
}

void RigidBodyTree::applyJointForce(SBStateRep& s, int body, int axis, const Real& r) const {
    SBDynamicVars& dynamicVars = s.updDynamicVars(*this); // check/adjust stage

    const RigidBodyNode& n = getRigidBodyNode(body);
    assert(0 <= axis && axis < n.getDOF());
    dynamicVars.appliedJointForces[n.getUIndex()+axis] = r;
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
    const SBDynamicVars& dynamicVars = s.getDynamicVars(*this);
    return dynamicVars.appliedJointForces;
}
const Vector_<SpatialVec>& 
RigidBodyTree::getAppliedBodyForces(const SBStateRep& s) const {
    const SBDynamicVars& dynamicVars = s.getDynamicVars(*this);
    return dynamicVars.appliedBodyForces;
}

const Vector& RigidBodyTree::getQDot(const SBStateRep& s) const {
    const SBMotionCache& motionCache = s.getMotionCache(*this);
    return motionCache.qdot;
}
const Vector& RigidBodyTree::getUDot(const SBStateRep& s) const {
    const SBDynamicCache& dynamicCache = s.getDynamicCache(*this);
    return dynamicCache.udot;
}
const Vector& RigidBodyTree::getQDotDot(const SBStateRep& s) const {
    const SBDynamicCache& dynamicCache = s.getDynamicCache(*this);
    return dynamicCache.qdotdot;
}

// Enforce coordinate constraints -- order doesn't matter.
void RigidBodyTree::enforceQuaternionConstraints(SBStateRep& s) const {
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->enforceQuaternionConstraints(s);
}

// Enforce loop constraints.
void RigidBodyTree::enforceLengthConstraints(SBStateRep& s) const {
    Vector& pos = updQ(s);
    Vector& vel = updU(s);
    lConstraints->enforce(s,pos,vel); //FIX: previous constraints still obeyed? (CDS)
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P.
void RigidBodyTree::prepareForDynamics(const SBStateRep& s) const {
    calcP(s);
}

// Given a set of spatial forces, calculate accelerations ignoring
// constraints. Must have already called prepareForDynamics().
// TODO: also applies stored internal forces (hinge torques) which
// will cause surprises if non-zero.
void RigidBodyTree::calcTreeForwardDynamics(const SBStateRep& s, 
                                            const SpatialVecList& spatialForces) const
{
    calcZ(s,spatialForces);
    calcTreeAccel(s);
    
    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcAccInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);
}

// Given a set of spatial forces, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const SBStateRep& s, 
                                            const SpatialVecList& spatialForces) const 
{
    SpatialVecList sFrc = spatialForces;
    calcTreeForwardDynamics(s, sFrc);
    if (lConstraints->calcConstraintForces(s)) {
        lConstraints->addInCorrectionForces(s, sFrc);
        calcTreeForwardDynamics(s, sFrc);
    }
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcP(const SBStateRep& s) const {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcP(s);
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
            rbNodeLevels[i][j]->calcY(s);
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

    Vector_<SpatialVec> zTemp(getTotalDOF()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(s, zTemp, X, JX);
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

