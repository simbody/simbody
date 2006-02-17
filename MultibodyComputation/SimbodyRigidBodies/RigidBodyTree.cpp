/**@file
 *
 * Implementation of RigidBodyTree.
 * Note: there must be no mention of atoms anywhere in this code.
 */

#include "RigidBodyTree.h"
#include "RigidBodyNode.h"
#include "LengthConstraints.h"

#include <string>

SBState::~SBState() {
    delete vars;
    delete cache;
}

SBState::SBState(const SBState& src) : vars(0), cache(0) {
    if (src.vars)  vars  = new SimbodyTreeVariables(*src.vars);
    if (src.cache) cache = new SimbodyTreeResults(*src.cache);
}

SBState& SBState::operator=(const SBState& src) {
    if (&src != this) {
        delete vars; delete cache; vars = 0; cache = 0;
        if (src.vars)  vars  = new SimbodyTreeVariables(*src.vars);
        if (src.cache) cache = new SimbodyTreeResults(*src.cache);
    }
    return *this;
}

void SBState::allocate() {
    delete vars; delete cache; // just in case
    vars  = new SimbodyTreeVariables();
    cache = new SimbodyTreeResults();   // stage is still UninitializedStage;
}

SBStage SBState::getStage() const { 
    if (!cache) return UninitializedStage;
    return cache->stage;
}


void RBStation::calcPosInfo(const SBState& s, RBStationRuntime& rt) const {
    rt.station_G = getNode().getR_GB(s) * station_B;
    rt.pos_G     = getNode().getOB_G(s) + rt.station_G;
}

void RBStation::calcVelInfo(const SBState& s, RBStationRuntime& rt) const {
    const Vec3& w_G = getNode().getSpatialAngVel(s);
    const Vec3& v_G = getNode().getSpatialLinVel(s);
    rt.stationVel_G = cross(w_G, rt.station_G);
    rt.vel_G = v_G + rt.stationVel_G;
}

void RBStation::calcAccInfo(const SBState& s, RBStationRuntime& rt) const {
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

void RBDistanceConstraint::calcPosInfo(const SBState& s, RBDistanceConstraintRuntime& rt) const
{
    assert(isValid() && runtimeIndex >= 0);
    for (int i=0; i<=1; ++i) stations[i].calcPosInfo(s, rt.stationRuntimes[i]);

    rt.fromTip1ToTip2_G = rt.stationRuntimes[1].pos_G - rt.stationRuntimes[0].pos_G;
    const double separation = rt.fromTip1ToTip2_G.norm();
    rt.unitDirection_G = rt.fromTip1ToTip2_G / separation;
    rt.posErr = distance - separation;
}

void RBDistanceConstraint::calcVelInfo(const SBState& s, RBDistanceConstraintRuntime& rt) const
{
    assert(isValid() && runtimeIndex >= 0);
    for (int i=0; i<=1; ++i) stations[i].calcVelInfo(s, rt.stationRuntimes[i]);

    rt.relVel_G = rt.stationRuntimes[1].vel_G - rt.stationRuntimes[0].vel_G;
    rt.velErr = ~rt.unitDirection_G * rt.relVel_G;
}

void RBDistanceConstraint::calcAccInfo(const SBState& s, RBDistanceConstraintRuntime& rt) const
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

// Add a new node, taking over the heap space.
int RigidBodyTree::addRigidBodyNode(RigidBodyNode&      parent,
                                    const TransformMat& referenceConfig, // body frame in parent
                                    RigidBodyNode*&     nodep)
{
    RigidBodyNode* n = nodep; nodep=0;  // take ownership
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

    n->refOrigin_P = referenceConfig.T(); // ignore rotation for now, it's always identity
    //n->X_GB = TransformMat(X_GB.getRotation(), 
    //                           X_GB.getTranslation() + n->refOrigin_P);
    //n->COM_G = child->X_GB.getTranslation() + n->COMstation_G;

    return nodeNum;
}

// Add a new ground node. Must be first node added during construction.
void RigidBodyTree::addGroundNode() {
    // Make sure this is the first body
    assert(nodeNum2NodeMap.size() == 0);
    assert(rbNodeLevels.size() == 0);

    RigidBodyNode* n = 
        RigidBodyNode::create(MassProperties(), TransformMat(), Joint::ThisIsGround,
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

// Here we lock in the topological structure of the multibody system,
// and compute allocation sizes for state variables.
void RigidBodyTree::realizeConstruction(const double& ctol, int verbose) {
    DOFTotal = SqDOFTotal = maxNQTotal = 0;
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            const int ndof = rbNodeLevels[i][j]->getDOF();
            DOFTotal += ndof; SqDOFTotal += ndof*ndof;
            maxNQTotal += rbNodeLevels[i][j]->getMaxNQ();
        }

    lConstraints = new LengthConstraints(*this, ctol, verbose);
    lConstraints->construct(distanceConstraints, dcRuntimeInfo);

    defaultState.allocate();
    defaultState.vars->allocateModelingVars(getNBodies(), getNConstraints());
    defaultState.cache->stage = BuiltStage;
    setDefaultModelingValues(defaultState);
    realizeModeling(defaultState);
    defaultState.vars->allocateAllVars(getTotalDOF(), getTotalQAlloc(), getNConstraints());

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getDefaultParameters(defaultState);
    realizeParameters(defaultState);

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getDefaultConfiguration(defaultState); 
    realizeConfiguration(defaultState);

    defaultState.vars->setVelocitiesToZero();
    realizeMotion(defaultState);

    defaultState.vars->clearForces();

    // just to be tidy -- shouldn't do anything
    defaultState.vars->setPrescribedAccelerationsToZero();
    realizeReaction(defaultState);
}

// Here we lock in modeling choices like whether to use quaternions or Euler
// angles; what joints are prescribed, etc.
void RigidBodyTree::realizeModeling(const SBState& s) const {
    assert(s.getStage() >= BuiltStage);
    if (s.getStage() >= ModeledStage) return;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeModeling(s); 

    s.cache->stage = ModeledStage;
}

// Here we lock in parameterization of the model, such as body masses.
void RigidBodyTree::realizeParameters(const SBState& s) const {
    assert(s.getStage() >= ModeledStage);
    if (s.getStage() >= ParametrizedStage) return;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeParameters(s); 
    s.cache->stage = ParametrizedStage;
}

// Set generalized coordinates: sweep from base to tips.
void RigidBodyTree::realizeConfiguration(const SBState& s)  {
    assert(s.getStage() >= ParametrizedStage);
    if (s.getStage() >= ConfiguredStage) return;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeConfiguration(s); 

    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcPosInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);

    s.cache->stage = ConfiguredStage;
}

// Set generalized speeds: sweep from base to tip.
// realizeConfiguration() must have been called already.
void RigidBodyTree::realizeMotion(const SBState& s)  {
    assert(s.getStage() >= ConfiguredStage);
    if (s.getStage() >= MovingStage) return;

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->realizeVelocity(s); 

    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcVelInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);

    s.cache->stage = MovingStage;
}

void RigidBodyTree::realizeReaction(const SBState& s)  const {
    assert(s.getStage() >= MovingStage);
    if (s.getStage() >= ReactingStage) return;

    assert(false); // TODO

    s.cache->stage = ReactingStage;
}

// Access to continuous state variables and their derivatives.

const Vector& RigidBodyTree::getQ(const SBState& s) const {
    assert(s.getStage() >= ConfiguredStage);
    return s.vars->q;
}
Vector& RigidBodyTree::updQ(SBState& s) const {
    assert(s.getStage() >= ParametrizedStage);
    s.cache->stage = ParametrizedStage; // back up if necessary
    return s.vars->q;
}
const Vector& RigidBodyTree::getU(const SBState& s) const {
    assert(s.getStage() >= MovingStage);
    return s.vars->u;
}
Vector& RigidBodyTree::updU(SBState& s) const {
    assert(s.getStage() >= ConfiguredStage);
    s.cache->stage = ConfiguredStage; // back up if necessary
    return s.vars->u;
}
const Vector& RigidBodyTree::getQdot(const SBState& s) const {
    assert(s.getStage() >= MovingStage);
    return s.cache->qdot;
}
const Vector& RigidBodyTree::getUdot(const SBState& s) const {
    assert(s.getStage() >= ReactingStage);
    return s.cache->udot;
}
const Vector& RigidBodyTree::getQdotDot(const SBState& s) const {
    assert(s.getStage() >= ReactingStage);
    return s.cache->qdotdot;
}

// Enforce coordinate constraints -- order doesn't matter.
void RigidBodyTree::enforceQuaternionConstraints(SBState& s) {
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) 
            rbNodeLevels[i][j]->enforceQuaternionConstraints(s);
}

// Enforce loop constraints.
void RigidBodyTree::enforceLengthConstraints(SBState& s) {
    Vector& pos = updQ(s);
    Vector& vel = updU(s);
    lConstraints->enforce(pos,vel); //FIX: previous constraints still obeyed? (CDS)
}


// Prepare for dynamics by calculating position-dependent quantities
// like the articulated body inertias P.
void RigidBodyTree::prepareForDynamics(const SBState& s) {
    calcP(s);
}

// Given a set of spatial forces, calculate accelerations ignoring
// constraints. Must have already called prepareForDynamics().
// TODO: also applies stored internal forces (hinge torques) which
// will cause surprises if non-zero.
void RigidBodyTree::calcTreeForwardDynamics(const SBState& s, const SpatialVecList& spatialForces) {
    calcZ(s,spatialForces);
    calcTreeAccel(s);
    
    for (size_t i=0; i < distanceConstraints.size(); ++i)
        distanceConstraints[i].calcAccInfo(s,
            dcRuntimeInfo[distanceConstraints[i].getRuntimeIndex()]);
}

// Given a set of spatial forces, calculate acclerations resulting from
// those forces and enforcement of acceleration constraints.
void RigidBodyTree::calcLoopForwardDynamics(const SBState& s, const SpatialVecList& spatialForces) {
    SpatialVecList sFrc = spatialForces;
    calcTreeForwardDynamics(s, sFrc);
    if (lConstraints->calcConstraintForces()) {
        lConstraints->addInCorrectionForces(sFrc);
        calcTreeForwardDynamics(s, sFrc);
    }
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcP(const SBState& s) {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcP(s);
}

// should be:
//   foreach tip {
//     traverse back to node which has more than one child hinge.
//   }
void RigidBodyTree::calcZ(const SBState& s, const SpatialVecList& spatialForces) {
    // level 0 for atoms whose position is fixed
    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcZ(s, spatialForces[node.getNodeNum()]);
        }
}

// Y is used for length constraints: sweep from base to tip.
void RigidBodyTree::calcY(const SBState& s) {
    for (int i=0; i < (int)rbNodeLevels.size(); i++)
        for (int j=0; j < (int)rbNodeLevels[i].size(); j++)
            rbNodeLevels[i][j]->calcY(s);
}

// Calc acceleration: sweep from base to tip.
void RigidBodyTree::calcTreeAccel(const SBState& s) {
    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->calcAccel(s);
}

void RigidBodyTree::fixVel0(Vector& vel) {
    lConstraints->fixVel0(vel);
}

// If V is a spatial velocity, and you have a X=d(something)/dV (one per body)
// this routine will return d(something)/du for internal generalized speeds u. If
// instead you have d(something)/dR where R is a spatial configuration, this routine
// returns d(something)/dq PROVIDED that dq/dt = u for all q's. That's not true for
// quaternions, so be careful how you use this routine.
// In Kane's terminology, we are calculating the product of a (generalized)
// partial velocity with some vector.
void RigidBodyTree::calcInternalGradientFromSpatial(const SBState& s, 
                                                    const Vector_<SpatialVec>& X,
                                                    Vector& JX) 
{
    assert(X.size() == getNBodies());
    assert(s.getStage() >= ConfiguredStage);

    Vector_<SpatialVec> zTemp(getTotalDOF()); zTemp.setToZero();
    JX.resize(getTotalDOF());

    for (int i=rbNodeLevels.size()-1 ; i>=0 ; i--)
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++) {
            RigidBodyNode& node = *rbNodeLevels[i][j];
            node.calcInternalGradientFromSpatial(s, zTemp, X, JX);
        }
}

// Pass in a set of internal forces in T; we'll modify them here.
void RigidBodyTree::calcConstraintCorrectedInternalForces(const SBState& s, Vector& T) {
    lConstraints->fixGradient(T);
}

// Get default parameter values for each node and write them into
// the supplied State. We expect the State to have been realized
// through ModeledStage.
void RigidBodyTree::getDefaultParameters(SBState& s) const {
    assert(s.getStage() >= ModeledStage);
    s.cache->stage = ModeledStage; // back up if necessary

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getDefaultConfiguration(s); 
}

// Get default position states for each node and write them into
// the supplied State. We expect the State to have been realized
// through ParametrizedStage since parameters might affect the
// configuration defaults.
void RigidBodyTree::getDefaultConfiguration(SBState& s) const {
    assert(s.getStage() >= ParametrizedStage);
    s.cache->stage = ParametrizedStage; // back up if necessary

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getDefaultConfiguration(s); 
}

// Get default velocity states for each node and write them into
// the supplied State. We expect the State to have been realized
// through ConfiguredStage. TODO: this is weird. Shouldn't the
// default velocities just be zero? Or perhaps some kind of 
// prescribed motion initialization should take place here?
void RigidBodyTree::getDefaultVelocity(SBState& s) const {
    assert(s.getStage() >= ConfiguredStage);
    s.cache->stage = ConfiguredStage; // back up if necessary

    for (int i=0 ; i<(int)rbNodeLevels.size() ; i++) 
        for (int j=0 ; j<(int)rbNodeLevels[i].size() ; j++)
            rbNodeLevels[i][j]->getDefaultVelocity(s); 
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

