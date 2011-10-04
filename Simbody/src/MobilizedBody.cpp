/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-8 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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
 * Private implementation of MobilizedBody, and its included subclasses which
 * represent the built-in mobilizer types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"

#include "MobilizedBodyImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {

    ////////////////////
    // MOBILIZED BODY //
    ////////////////////

MobilizedBody::MobilizedBody() {
}

MobilizedBody::MobilizedBody(MobilizedBodyImpl* r) : HandleBase(r) {
}

MobilizedBody& MobilizedBody::addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
    updImpl().addOutboardDecoration(X_MD,g);
    return *this;
}
MobilizedBody& MobilizedBody::addInboardDecoration (const Transform& X_MbD, const DecorativeGeometry& g) {
    updImpl().addInboardDecoration(X_MbD,g);
    return *this;
}

const SimbodyMatterSubsystem& MobilizedBody::getMatterSubsystem() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMatterSubsystem() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}


MobilizedBodyIndex MobilizedBody::getMobilizedBodyIndex() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMobilizedBodyIndex() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().getMyMobilizedBodyIndex();
}

const MobilizedBody& MobilizedBody::getParentMobilizedBody() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getParentMobilizedBody() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().getMyMatterSubsystemRep().getMobilizedBody(getImpl().getMyParentMobilizedBodyIndex());
}

const MobilizedBody& MobilizedBody::getBaseMobilizedBody() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getBaseMobilizedBody() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().getMyMatterSubsystemRep().getMobilizedBody(getImpl().getMyBaseBodyMobilizedBodyIndex());
}

bool MobilizedBody::isInSubsystem() const {
    return !isEmptyHandle() && getImpl().isInSubsystem();
}

bool MobilizedBody::isInSameSubsystem(const MobilizedBody& otherBody) const {
    return isInSubsystem() && otherBody.isInSubsystem()
           && getMatterSubsystem().isSameSubsystem(otherBody.getMatterSubsystem());
}

bool MobilizedBody::isSameMobilizedBody(const MobilizedBody& otherBody) const {
    return !isEmptyHandle() && isSameHandle(otherBody);
}

bool MobilizedBody::isGround() const {
    return isInSubsystem() && isSameMobilizedBody(getMatterSubsystem().getGround());
}

int MobilizedBody::getLevelInMultibodyTree() const {
    return getImpl().getMyLevel();
}

SimbodyMatterSubsystem& MobilizedBody::updMatterSubsystem() {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "updMatterSubsystem() called on a MobilizedBody that is not part of a subsystem.");
    return updImpl().updMyMatterSubsystemRep().updMySimbodyMatterSubsystemHandle();
}

const Body& MobilizedBody::getBody() const {
    return getImpl().theBody;
}
Body& MobilizedBody::updBody() {
    getImpl().invalidateTopologyCache();
    return updImpl().theBody;
}
MobilizedBody& MobilizedBody::setBody(const Body& b) {
    updBody() = b;
    return *this;
}
MobilizedBody& MobilizedBody::setDefaultInboardFrame (const Transform& X_PMb) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultInboardFrame = X_PMb;
    return *this;
}
MobilizedBody& MobilizedBody::setDefaultOutboardFrame(const Transform& X_BM) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultOutboardFrame = X_BM;
    return *this;
}

const Transform& MobilizedBody::getDefaultInboardFrame() const {
    return getImpl().defaultInboardFrame;
}
const Transform& MobilizedBody::getDefaultOutboardFrame() const {
    return getImpl().defaultOutboardFrame;
}

// Access to State

const MassProperties& MobilizedBody::getBodyMassProperties(const State& s) const {
    return getImpl().getBodyMassProperties(s);
}

const SpatialInertia& MobilizedBody::
getBodySpatialInertiaInGround(const State& s) const {
    return getImpl().getBodySpatialInertiaInGround(s);
}

const Transform& MobilizedBody::getInboardFrame (const State& s) const {
    return getImpl().getInboardFrame(s);
}
const Transform& MobilizedBody::getOutboardFrame(const State& s) const {
    return getImpl().getOutboardFrame(s);
}

void MobilizedBody::setInboardFrame (State& s, const Transform& X_PMb) const {
    getImpl().setInboardFrame(s, X_PMb);
}
void MobilizedBody::setOutboardFrame(State& s, const Transform& X_BM) const {
    getImpl().setOutboardFrame(s, X_BM);
}

const Transform& MobilizedBody::getBodyTransform(const State& s) const {
    return getImpl().getBodyTransform(s);
}
const SpatialVec& MobilizedBody::getBodyVelocity(const State& s) const {
    return getImpl().getBodyVelocity(s);
}
const SpatialVec& MobilizedBody::getBodyAcceleration(const State& s) const {
    return getImpl().getBodyAcceleration(s);
}

const Transform& MobilizedBody::getMobilizerTransform(const State& s) const {
    return getImpl().getMobilizerTransform(s);
}
const SpatialVec& MobilizedBody::getMobilizerVelocity(const State& s) const {
    return getImpl().getMobilizerVelocity(s);
}

void MobilizedBody::setQToFitTransform(State& s, const Transform& X_MbM) const { 
    getImpl().setQToFitTransform(s,X_MbM); 
}
void MobilizedBody::setQToFitRotation(State& s, const Rotation& R_MbM) const { 
    getImpl().setQToFitRotation(s,R_MbM); 
}
void MobilizedBody::setQToFitTranslation(State& s, const Vec3& p_MbM) const { 
    getImpl().setQToFitTranslation(s,p_MbM);
}
void MobilizedBody::setUToFitVelocity(State& s, const SpatialVec& V_MbM) const { 
    getImpl().setUToFitVelocity(s,V_MbM);
}
void MobilizedBody::setUToFitAngularVelocity(State& s, const Vec3& w_MbM) const { 
    getImpl().setUToFitAngularVelocity(s,w_MbM);
}
void MobilizedBody::setUToFitLinearVelocity(State& s, const Vec3& v_MbM) const { 
    getImpl().setUToFitLinearVelocity(s,v_MbM);
}

SpatialVec MobilizedBody::getHCol(const State& s, MobilizerUIndex ux) const {
    SimTK_INDEXCHECK(ux, getNumU(s), "MobilizedBody::getHCol()");
    return getImpl().getHCol(s,ux);
}

SpatialVec MobilizedBody::getH_FMCol(const State& s, MobilizerUIndex ux) const {
    SimTK_INDEXCHECK(ux, getNumU(s), "MobilizedBody::getH_FMCol()");
    return getImpl().getH_FMCol(s,ux);
}

int MobilizedBody::getNumQ(const State& s) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s, qStart, nq);
    return nq;
}

QIndex MobilizedBody::getFirstQIndex(const State& s) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s, qStart, nq);
    return qStart;
}

int MobilizedBody::getNumU(const State& s) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s, uStart, nu);
    return nu;
}

UIndex MobilizedBody::getFirstUIndex(const State& s) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s, uStart, nu);
    return uStart;
}

Motion::Method MobilizedBody::getQMotionMethod(const State& s) const 
{   return getImpl().getQMotionMethod(s); }
Motion::Method MobilizedBody::getUMotionMethod(const State& s) const
{   return getImpl().getUMotionMethod(s); }
Motion::Method MobilizedBody::getUDotMotionMethod(const State& s) const
{   return getImpl().getUDotMotionMethod(s); }

Real  MobilizedBody::getOneFromQPartition(const State& s, int which, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s, qStart, nq);
    assert(0 <= which && which < nq);
    return qlike[qStart+which];
}
Real& MobilizedBody::updOneFromQPartition(const State& s, int which, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s, qStart, nq);
    assert(0 <= which && which < nq);
    return qlike[qStart+which];
}

Real  MobilizedBody::getOneFromUPartition(const State& s, int which, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s, uStart, nu);
    assert(0 <= which && which < nu);
    return ulike[uStart+which];
}
Real& MobilizedBody::updOneFromUPartition(const State& s, int which, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s, uStart, nu);
    assert(0 <= which && which < nu);
    return ulike[uStart+which];
}

void MobilizedBody::
convertQForceToUForce(  const State&                        state,
                        const Array_<Real,MobilizerQIndex>& fq,
                        Array_<Real,MobilizerUIndex>&       fu) const
{
    const MobilizedBodyImpl& impl = getImpl();
    const SimbodyMatterSubsystemRep& matter = impl.getMyMatterSubsystemRep();

    UIndex uStart; int nu;
    QIndex qStart; int nq;
    matter.findMobilizerUs(state, impl.getMyMobilizedBodyIndex(), uStart, nu);
    matter.findMobilizerQs(state, impl.getMyMobilizedBodyIndex(), qStart, nq);
    assert(fq.size() == nq);

    fu.resize(nu);
    const RigidBodyNode& node = impl.getMyRigidBodyNode();
    const SBStateDigest digest(state, matter, Stage::Velocity);
    node.multiplyByN(digest, true, fq.cbegin(), fu.begin());
}


void MobilizedBody::applyBodyForce(const State& s, const SpatialVec& spatialForceInG, 
                                   Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNumBodies());
    bodyForces[getMobilizedBodyIndex()] += spatialForceInG;
}

void MobilizedBody::applyBodyTorque(const State& s, const Vec3& torqueInG, 
                     Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNumBodies());
    bodyForces[getMobilizedBodyIndex()][0] += torqueInG; // don't change force
}

void MobilizedBody::applyForceToBodyPoint(const State& s, const Vec3& pointInB, const Vec3& forceInG,
                           Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNumBodies());
    const Rotation& R_GB = getBodyTransform(s).R();
    bodyForces[getMobilizedBodyIndex()] += SpatialVec((R_GB*pointInB) % forceInG, forceInG);
}

Real MobilizedBody::getOneQ(const State& s, int which) const {
    return getOneFromQPartition(s,which,getImpl().getMyMatterSubsystemRep().getQ(s));
}

void MobilizedBody::setOneQ(State& s, int which, Real value) const {
    updOneFromQPartition(s,which,getImpl().getMyMatterSubsystemRep().updQ(s)) = value;
}

Real MobilizedBody::getOneU(const State& s, int which) const {
    return getOneFromUPartition(s,which,getImpl().getMyMatterSubsystemRep().getU(s));
}
void MobilizedBody::setOneU(State& s, int which, Real value) const {
    updOneFromUPartition(s,which,getImpl().getMyMatterSubsystemRep().updU(s)) = value;
}

Real MobilizedBody::getOneQDot(const State& s, int which) const {
    return getOneFromQPartition(s,which,getImpl().getMyMatterSubsystemRep().getQDot(s));
}
Real MobilizedBody::getOneUDot(const State& s, int which) const {
    return getOneFromUPartition(s,which,getImpl().getMyMatterSubsystemRep().getUDot(s));
}
Real MobilizedBody::getOneQDotDot(const State& s, int which) const {
    return getOneFromQPartition(s,which,getImpl().getMyMatterSubsystemRep().getQDotDot(s));
}
Real MobilizedBody::getOneTau(const State& s, MobilizerUIndex which) const {
    const MobilizedBodyImpl&         mbimpl = MobilizedBody::getImpl();
    const MobilizedBodyIndex         mbx    = mbimpl.getMyMobilizedBodyIndex();
    const SimbodyMatterSubsystemRep& matter = mbimpl.getMyMatterSubsystemRep();
    const SBModelCache&              mc     = matter.getModelCache(s);
    const SBModelPerMobodInfo&
        mobodModelInfo = mc.getMobodModelInfo(mbx);
    const int nu = mobodModelInfo.nUInUse;

    SimTK_INDEXCHECK(which, nu, "MobilizedBody::getOneTau()");

    const SBInstanceCache& ic = matter.getInstanceCache(s);
    const SBInstancePerMobodInfo& 
        mobodInstanceInfo = ic.getMobodInstanceInfo(mbx);
    if (mobodInstanceInfo.udotMethod == Motion::Free)
        return 0; // not prescribed

    const SBTreeAccelerationCache& tac = matter.getTreeAccelerationCache(s);
    return tac.presMotionForces[mobodInstanceInfo.firstPresForce + which];
}

Vector MobilizedBody::getTauAsVector(const State& s) const {
    const MobilizedBodyImpl&         mbimpl = MobilizedBody::getImpl();
    const MobilizedBodyIndex         mbx    = mbimpl.getMyMobilizedBodyIndex();
    const SimbodyMatterSubsystemRep& matter = mbimpl.getMyMatterSubsystemRep();
    const SBModelCache&              mc     = matter.getModelCache(s);
    const SBModelPerMobodInfo&
        mobodModelInfo = mc.getMobodModelInfo(mbx);
    const int nu = mobodModelInfo.nUInUse;

    const SBInstanceCache& ic = matter.getInstanceCache(s);
    const SBInstancePerMobodInfo& 
        mobodInstanceInfo = ic.getMobodInstanceInfo(mbx);
    if (mobodInstanceInfo.udotMethod == Motion::Free)
        return Vector(nu, Real(0)); // not prescribed

    const SBTreeAccelerationCache& tac = matter.getTreeAccelerationCache(s);
    return tac.presMotionForces(mobodInstanceInfo.firstPresForce, nu);
}

Vector MobilizedBody::getQAsVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQ(s)(qStart,nq);
}

void MobilizedBody::setQFromVector(State& s, const Vector& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s, qStart, nq);
    assert(q.size() == nq);
    mbr.getMyMatterSubsystemRep().updQ(s)(qStart,nq) = q;
}

Vector MobilizedBody::getUAsVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s, uStart, nu);
    return mbr.getMyMatterSubsystemRep().getU(s)(uStart,nu);
}

void MobilizedBody::setUFromVector(State& s, const Vector& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s, uStart, nu);
    assert(u.size() == nu);
    mbr.getMyMatterSubsystemRep().updU(s)(uStart,nu) = u;
}

Vector MobilizedBody::getQDotAsVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQDot(s)(qStart,nq);
}

Vector MobilizedBody::getUDotAsVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s, uStart, nu);
    return mbr.getMyMatterSubsystemRep().getUDot(s)(uStart,nu);
}

Vector MobilizedBody::getQDotDotAsVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)(qStart,nq);
}

MobilizedBody& MobilizedBody::cloneForNewParent(MobilizedBody& parent) const {
    MobilizedBody copyBody;
    copyBody.copyAssign(*this);
    copyBody.updImpl().myMatterSubsystemRep = 0;
    copyBody.updImpl().myRBnode = 0;
    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(), copyBody);
    return parent.updMatterSubsystem().updMobilizedBody(copyBody.getMobilizedBodyIndex());
}

void MobilizedBody::adoptMotion(Motion& ownerHandle) {
   updImpl().adoptMotion(ownerHandle);
}

void MobilizedBody::clearMotion() {
   updImpl().clearMotion();
}

bool MobilizedBody::hasMotion() const {
    return getImpl().hasMotion();
}

const Motion& MobilizedBody::getMotion() const {
    return getImpl().getMotion();
}

    /////////////////////////
    // MOBILIZED BODY IMPL //
    /////////////////////////

void MobilizedBodyImpl::findMobilizerQs(const State& s, QIndex& qStart, int& nq) const {
    getMyMatterSubsystemRep()
        .findMobilizerQs(s, myMobilizedBodyIndex, qStart, nq);
}
void MobilizedBodyImpl::findMobilizerUs(const State& s, UIndex& uStart, int& nu) const {
    getMyMatterSubsystemRep()
        .findMobilizerUs(s, myMobilizedBodyIndex, uStart, nu);
}

Motion::Method MobilizedBodyImpl::getQMotionMethod(const State& s) const 
{   return getMyInstanceInfo(s).qMethod; }
Motion::Method MobilizedBodyImpl::getUMotionMethod(const State& s) const 
{   return getMyInstanceInfo(s).uMethod; }
Motion::Method MobilizedBodyImpl::getUDotMotionMethod(const State& s) const 
{   return getMyInstanceInfo(s).udotMethod; }

void MobilizedBodyImpl::copyOutDefaultQ(const State& s, Vector& qDefault) const {
    SimTK_STAGECHECK_GE_ALWAYS(getMyMatterSubsystemRep().getStage(s), Stage::Topology,
        "MobilizedBody::copyOutDefaultQ()");
    QIndex qStart; int nq;
    findMobilizerQs(s, qStart, nq);
    if (nq)
        copyOutDefaultQImpl(nq, &qDefault[qStart]);
}

    // TODO: currently we delegate these requests to the RigidBodyNodes. 
    // Probably most of this functionality should be handled directly
    // by the MobilizedBody objects.

void MobilizedBodyImpl::setQToFitTransform(State& s, const Transform& X_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitTransform(digest, X_MbM, q);
}
void MobilizedBodyImpl::setQToFitRotation(State& s, const Rotation& R_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitRotation(digest, R_MbM, q);
}
void MobilizedBodyImpl::setQToFitTranslation(State& s, const Vec3& p_MbM) const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitTranslation(digest, p_MbM, q);
}

void MobilizedBodyImpl::setUToFitVelocity(State& s, const SpatialVec& V_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitVelocity(digest, q, V_MbM, u);
}
void MobilizedBodyImpl::setUToFitAngularVelocity(State& s, const Vec3& w_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitAngularVelocity(digest, q, w_MbM, u);
}
void MobilizedBodyImpl::setUToFitLinearVelocity(State& s, const Vec3& v_MbM)  const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitLinearVelocity(digest, q, v_MbM, u);
}

    // REALIZE TOPOLOGY
const RigidBodyNode& MobilizedBodyImpl::realizeTopology
   (State& s, UIndex& nxtU, USquaredIndex& nxtUSq, QIndex& nxtQ) const
{
    delete myRBnode;
    myRBnode = createRigidBodyNode(nxtU,nxtUSq,nxtQ);

    assert(myMobilizedBodyIndex.isValid());
    assert(myParentIndex.isValid() || myMobilizedBodyIndex == GroundIndex);

    if (myParentIndex.isValid()) {
        // not ground
        const MobilizedBodyImpl& parent = 
            myMatterSubsystemRep->getMobilizedBody(myParentIndex).getImpl();
        assert(myLevel == parent.myRBnode->getLevel() + 1);
        parent.myRBnode->addChild(myRBnode);
        myRBnode->setParent(parent.myRBnode);
    }

    myRBnode->setLevel(myLevel);
    myRBnode->setNodeNum(myMobilizedBodyIndex);

    // Realize Motion topology.
    if (hasMotion())
        getMotion().getImpl().realizeTopology(s);

    // Realize MobilizedBody-specific topology.
    realizeTopologyVirtual(s);
    return *myRBnode;
}

    // REALIZE MODEL
void MobilizedBodyImpl::realizeModel(State& s) const {
    if (hasMotion())
        getMotion().getImpl().realizeModel(s);
    realizeModelVirtual(s);
}

    // REALIZE INSTANCE
void MobilizedBodyImpl::realizeInstance(const SBStateDigest& sbs) const {
    if (hasMotion())
        getMotion().getImpl().realizeInstance(sbs.getState());
    realizeInstanceVirtual(sbs.getState());
}

    // REALIZE TIME
void MobilizedBodyImpl::realizeTime(const SBStateDigest& sbs) const {
    const MobilizedBodyIndex      mbx       = getMyMobilizedBodyIndex();
    const SBInstanceCache&        ic        = sbs.getInstanceCache();
    const SBInstancePerMobodInfo& instInfo  = ic.getMobodInstanceInfo(mbx);

    // Note that we only need to deal with explicitly prescribed motion;
    // if the mobilizer is prescribed to zero it will be dealt with elsewhere.
    if (instInfo.qMethod==Motion::Prescribed) {
        const SBModelCache& mc = sbs.getModelCache();
        const SBModelPerMobodInfo& 
            modelInfo = mc.getMobodModelInfo(mbx);
        const int            nq  = modelInfo.nQInUse;
        const PresQPoolIndex pqx = instInfo.firstPresQ;
        SBTimeCache& tc = sbs.updTimeCache();
        const MotionImpl& motion = getMotion().getImpl();

        motion.calcPrescribedPosition(sbs.getState(), nq, 
                                      &tc.presQPool[pqx]);
    }

    realizeTimeVirtual(sbs.getState());
}

    // REALIZE POSITION
void MobilizedBodyImpl::realizePosition(const SBStateDigest& sbs) const {
    const MobilizedBodyIndex    mbx = getMyMobilizedBodyIndex();
    const SBInstanceCache&      ic  = sbs.getInstanceCache();
    const SBInstancePerMobodInfo& 
        instInfo = ic.getMobodInstanceInfo(mbx);

    // Note that we only need to deal with explicitly prescribed motion;
    // if the mobilizer is prescribed to zero it will be dealt with elsewhere.
    // Prescribed u may be due to nonholonomic prescribed u, or to derivative
    // of holonomic prescribed q.
    if (instInfo.uMethod==Motion::Prescribed) {
        const SBModelCache& mc = sbs.getModelCache();
        const SBModelPerMobodInfo& 
            modelInfo = mc.getMobodModelInfo(mbx);
        const int            nu  = modelInfo.nUInUse;
        const PresUPoolIndex pux = instInfo.firstPresU;
        SBConstrainedPositionCache& cpc = sbs.updConstrainedPositionCache();
        const MotionImpl& motion = getMotion().getImpl();

        if (instInfo.qMethod==Motion::Prescribed) {
            // Holonomic
            const int nq = modelInfo.nUInUse;
            const RigidBodyNode& rbn = getMyRigidBodyNode();

            if (rbn.isQDotAlwaysTheSameAsU()) {
                assert(nq==nu);
                motion.calcPrescribedPositionDot(sbs.getState(), nu,
                                                 &cpc.presUPool[pux]);
            } else {
                Real qdot[8]; // we won't use all of these -- max is 7
                motion.calcPrescribedPositionDot(sbs.getState(), nq, qdot);
                // u = N^-1 qdot
                rbn.multiplyByNInv(sbs, false, qdot, &cpc.presUPool[pux]);
            }
        } else { // Non-holonomic
            motion.calcPrescribedVelocity(sbs.getState(), nu,
                                          &cpc.presUPool[pux]);
        }
    }

    realizePositionVirtual(sbs.getState());
}

    // REALIZE VELOCITY
void MobilizedBodyImpl::realizeVelocity(const SBStateDigest& sbs) const {
    // We don't deal with prescribed accelerations until Dynamics stage.
    realizeVelocityVirtual(sbs.getState());
}

    // REALIZE DYNAMICS
void MobilizedBodyImpl::realizeDynamics(const SBStateDigest& sbs) const {
    const MobilizedBodyIndex      mbx      = getMyMobilizedBodyIndex();
    const SBInstanceCache&        ic       = sbs.getInstanceCache();
    const SBInstancePerMobodInfo& instInfo = ic.getMobodInstanceInfo(mbx);

    // Note that we only need to deal with explicitly prescribed motion;
    // if the mobilizer is prescribed to zero it will be dealt with elsewhere.
    // Prescribed udot may be due to acceleration-only prescribed udot,
    // derivative of nonholonomic prescribed u, or to 2nd derivative
    // of holonomic prescribed q.
    if (instInfo.udotMethod==Motion::Prescribed) {
        const SBModelCache&         mc = sbs.getModelCache();
        const SBModelPerMobodInfo&  modelInfo = mc.getMobodModelInfo(mbx);
        const int                   nu = modelInfo.nUInUse;
        const UIndex                ux = modelInfo.firstUIndex;
        const PresUDotPoolIndex     pudx = instInfo.firstPresUDot;
        SBDynamicsCache& dc = sbs.updDynamicsCache();
        const MotionImpl& motion = getMotion().getImpl();
        Real* presUDotp = &dc.presUDotPool[pudx];

        if (instInfo.qMethod==Motion::Prescribed) {
            // Holonomic
            const int nq = modelInfo.nQInUse;
            const RigidBodyNode& rbn = getMyRigidBodyNode();

            if (rbn.isQDotAlwaysTheSameAsU()) {
                assert(nq==nu);
                motion.calcPrescribedPositionDotDot(sbs.getState(), nu,
                                                    &dc.presUDotPool[pudx]);
            } else {
                Real ndotU[8]; // remainder term NDot*u (nq of these; max is 7)
                const Vector& u = sbs.getU();
                rbn.multiplyByNDot(sbs, false, &u[ux], ndotU);

                Real qdotdot[8]; // nq of these -- max is 7
                motion.calcPrescribedPositionDotDot(sbs.getState(), nq, qdotdot);

                for (int i=0; i < nq; ++i)
                    qdotdot[i] -= ndotU[i]; // velocity correction

                // udot = N^-1 (qdotdot - NDot*u)
                rbn.multiplyByNInv(sbs, false, qdotdot, presUDotp);
            }
        } else if (instInfo.uMethod==Motion::Prescribed) { 
            // Non-holonomic
            motion.calcPrescribedVelocityDot(sbs.getState(), nu, presUDotp);
        } else { 
            // Acceleration-only
            motion.calcPrescribedAcceleration(sbs.getState(), nu, presUDotp);
        }
    }

    realizeDynamicsVirtual(sbs.getState());
}

    // REALIZE ACCELERATION
void MobilizedBodyImpl::realizeAcceleration(const SBStateDigest& sbs) const {
    realizeAccelerationVirtual(sbs.getState());
}

    // REALIZE REPORT
void MobilizedBodyImpl::realizeReport(const SBStateDigest& sbs) const {
    realizeReportVirtual(sbs.getState());
}

    /////////////////////////
    // MOBILIZED BODY::PIN //
    /////////////////////////

MobilizedBody::Pin::Pin(Direction d) : MobilizedBody(new PinImpl(d)) {
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new PinImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Transform& inbFrame,
                        const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new PinImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}


Real MobilizedBody::Pin::getDefaultQ() const {
    return getImpl().defaultQ;
}

MobilizedBody::Pin& MobilizedBody::Pin::setDefaultQ(Real q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

Real MobilizedBody::Pin::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Pin::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Pin::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Pin::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Pin::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Pin::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Pin::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Pin::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Pin::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Pin::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Pin::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Pin, MobilizedBody::PinImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::SLIDER //
    ////////////////////////////

MobilizedBody::Slider::Slider(Direction d) : MobilizedBody(new SliderImpl(d)) {
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new SliderImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Transform& inbFrame,
                              const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new SliderImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

Real MobilizedBody::Slider::getDefaultQ() const {
    return getImpl().defaultQ;
}

MobilizedBody::Slider& MobilizedBody::Slider::setDefaultQ(Real q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

Real MobilizedBody::Slider::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Slider::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Slider::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Slider::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Slider::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Slider::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Slider::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Slider::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Slider::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Slider::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Slider::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Slider, MobilizedBody::SliderImpl, MobilizedBody);

    ///////////////////////////////
    // MOBILIZED BODY::UNIVERSAL //
    ///////////////////////////////

MobilizedBody::Universal::Universal(Direction d) : MobilizedBody(new UniversalImpl(d)) {
}


MobilizedBody::Universal::Universal(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new UniversalImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Universal::Universal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new UniversalImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Universal, MobilizedBody::UniversalImpl, MobilizedBody);

    //////////////////////////////
    // MOBILIZED BODY::CYLINDER //
    //////////////////////////////

MobilizedBody::Cylinder::Cylinder(Direction d) : MobilizedBody(new CylinderImpl(d)) {
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new CylinderImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new CylinderImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Cylinder, MobilizedBody::CylinderImpl, MobilizedBody);

    //////////////////////////////////
    // MOBILIZED BODY::BEND STRETCH //
    //////////////////////////////////

MobilizedBody::BendStretch::BendStretch(Direction d) : MobilizedBody(new BendStretchImpl(d)) {
}


MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new BendStretchImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                                        const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new BendStretchImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::BendStretch, MobilizedBody::BendStretchImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::PLANAR //
    ////////////////////////////

MobilizedBody::Planar::Planar(Direction d) : MobilizedBody(new PlanarImpl(d)) {
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new PlanarImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Transform& inbFrame,
                              const Body& body, const Transform& outbFrame,
                              Direction d) 
:   MobilizedBody(new PlanarImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

const Vec3& MobilizedBody::Planar::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::Planar& MobilizedBody::Planar::setDefaultQ(const Vec3& q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

const Vec3& MobilizedBody::Planar::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Planar::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::Planar::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::Planar::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Planar::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Planar::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Planar::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::Planar::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Planar::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::Planar::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Planar::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Planar, MobilizedBody::PlanarImpl, MobilizedBody);


    //////////////////////////////////////
    // MOBILIZED BODY::SPHERICAL COORDS //
    //////////////////////////////////////

MobilizedBody::SphericalCoords::SphericalCoords(Direction d) 
:   MobilizedBody(new SphericalCoordsImpl(d)) {}

MobilizedBody::SphericalCoords::SphericalCoords
   (MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new SphericalCoordsImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::SphericalCoords::SphericalCoords
   (MobilizedBody&  parent, const Transform& inbFrame,
    const Body&     body,   const Transform& outbFrame,
    Direction       d) 
:   MobilizedBody(new SphericalCoordsImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::SphericalCoords::SphericalCoords
   (MobilizedBody&  parent,         const Transform&    inbFrame,
    const Body&     body,           const Transform&    outbFrame,
    Real            azimuthOffset,  bool                azimuthNegated,
    Real            zenithOffset,   bool                zenithNegated,
    CoordinateAxis  radialAxis,     bool                radialNegated,
    Direction       d) 
:   MobilizedBody(new SphericalCoordsImpl(azimuthOffset, azimuthNegated,
                                          zenithOffset,  zenithNegated,
                                          radialAxis,    radialNegated,
                                          d))
{
    SimTK_APIARGCHECK_ALWAYS(radialAxis != YAxis, "MobilizedBody::SphericalCoords", "ctor",
        "SphericalCoords translation (radial) axis must be X or Z; Y is not allowed.");

    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::SphericalCoords& 
MobilizedBody::SphericalCoords::setRadialAxis(CoordinateAxis axis) {
    SimTK_APIARGCHECK_ALWAYS(axis != YAxis, "MobilizedBody::SphericalCoords", "setRadialAxis",
        "SphericalCoords translation (radial) axis must be X or Z; Y is not allowed.");
    getImpl().invalidateTopologyCache();
    updImpl().axisT = axis;
    return *this;
}

MobilizedBody::SphericalCoords& 
MobilizedBody::SphericalCoords::setNegateAzimuth(bool shouldNegate) {
    getImpl().invalidateTopologyCache();
    updImpl().negAz = shouldNegate;
    return *this;
}
MobilizedBody::SphericalCoords& 
MobilizedBody::SphericalCoords::setNegateZenith(bool shouldNegate) {
    getImpl().invalidateTopologyCache();
    updImpl().negZe = shouldNegate;
    return *this;
}
MobilizedBody::SphericalCoords& 
MobilizedBody::SphericalCoords::setNegateRadial(bool shouldNegate) {
    getImpl().invalidateTopologyCache();
    updImpl().negT = shouldNegate;
    return *this;
}

CoordinateAxis MobilizedBody::SphericalCoords::getRadialAxis()    const {return getImpl().axisT;}
bool           MobilizedBody::SphericalCoords::isAzimuthNegated() const {return getImpl().negAz;}
bool           MobilizedBody::SphericalCoords::isZenithNegated()  const {return getImpl().negZe;}
bool           MobilizedBody::SphericalCoords::isRadialNegated()  const {return getImpl().negT;}

const Vec3& MobilizedBody::SphericalCoords::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::SphericalCoords& MobilizedBody::SphericalCoords::setDefaultQ(const Vec3& q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

const Vec3& MobilizedBody::SphericalCoords::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::SphericalCoords::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::SphericalCoords::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::SphericalCoords::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::SphericalCoords::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::SphericalCoords::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::SphericalCoords::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::SphericalCoords::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::SphericalCoords::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::SphericalCoords::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::SphericalCoords::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::SphericalCoords, MobilizedBody::SphericalCoordsImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::GIMBAL //
    ////////////////////////////

MobilizedBody::Gimbal::Gimbal(Direction d) : MobilizedBody(new GimbalImpl(d)) {
}


MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new GimbalImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new GimbalImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Gimbal& MobilizedBody::Gimbal::setDefaultRadius(Real r) {
    getImpl().invalidateTopologyCache();
    updImpl().setDefaultRadius(r);
    return *this;
}

Real MobilizedBody::Gimbal::getDefaultRadius() const {
    return getImpl().getDefaultRadius();
}

const Vec3& MobilizedBody::Gimbal::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::Gimbal& MobilizedBody::Gimbal::setDefaultQ(const Vec3& q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

const Vec3& MobilizedBody::Gimbal::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Gimbal::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::Gimbal::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::Gimbal::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Gimbal::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Gimbal::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Gimbal::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::Gimbal::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Gimbal::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::Gimbal::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Gimbal::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Gimbal, MobilizedBody::GimbalImpl, MobilizedBody);

    // GimbalImpl

void MobilizedBody::GimbalImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // Call the superclass implementation to get standard default geometry.
    
    MobilizedBodyImpl::calcDecorativeGeometryAndAppendImpl(s, stage, geom);
    
    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the parent and child mobilizer frame
    // placement on the body, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Transform& X_PMb = getInboardFrame(s);
        const Transform& X_BM  = getOutboardFrame(s);

        // On the inboard body, draw a solid sphere and a wireframe one attached to it for
        // easier visualization of its rotation. These are at about 90% of the radius.
        geom.push_back(DecorativeSphere(0.92*getDefaultRadius())
                                            .setColor(Gray)
                                            .setRepresentation(DecorativeGeometry::DrawSurface)
                                            .setOpacity(0.5)
                                            .setResolution(0.75)
                                            .setBodyId(getMyParentMobilizedBodyIndex())
                                            .setTransform(X_PMb));
        geom.push_back(DecorativeSphere(0.90*getDefaultRadius())
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMyParentMobilizedBodyIndex())
            .setTransform(X_PMb));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(getMyMobilizedBodyIndex())
                                            .setTransform(X_BM));
    }
}

    ///////////////////////////////////////////////////
    // MOBILIZED BODY::BALL (ORIENTATION, SPHERICAL) //
    ///////////////////////////////////////////////////

MobilizedBody::Ball::Ball(Direction d) : MobilizedBody(new BallImpl(d)) {
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new BallImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new BallImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ball& MobilizedBody::Ball::setDefaultRadius(Real r) {
    getImpl().invalidateTopologyCache();
    updImpl().setDefaultRadius(r);
    return *this;
}

Real MobilizedBody::Ball::getDefaultRadius() const {
    return getImpl().getDefaultRadius();
}

const Quaternion& MobilizedBody::Ball::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::Ball& MobilizedBody::Ball::setDefaultQ(const Quaternion& quat) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = quat;
    return *this;
}

const Vec4& MobilizedBody::Ball::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Ball::setQ(State& s, const Vec4& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    Vec4::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec4& MobilizedBody::Ball::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec4& MobilizedBody::Ball::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}

const Vec3& MobilizedBody::Ball::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Ball::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Ball::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec4& MobilizedBody::Ball::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Ball::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec4& MobilizedBody::Ball::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Ball::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ball, MobilizedBody::BallImpl, MobilizedBody);

    // BallImpl

void MobilizedBody::BallImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // Call the superclass implementation to get standard default geometry.
    
    MobilizedBodyImpl::calcDecorativeGeometryAndAppendImpl(s, stage, geom);

    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the parent and child mobilizer frame
    // placement on the body, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Transform& X_PMb = getInboardFrame(s);
        const Transform& X_BM  = getOutboardFrame(s);

        // On the inboard body, draw a solid sphere and a wireframe one attached to it for
        // easier visualization of its rotation. These are at about 90% of the radius.
        geom.push_back(DecorativeSphere(0.92*getDefaultRadius())
                                            .setColor(Gray)
                                            .setRepresentation(DecorativeGeometry::DrawSurface)
                                            .setOpacity(0.5)
                                            .setResolution(0.75)
                                            .setBodyId(getMyParentMobilizedBodyIndex())
                                            .setTransform(X_PMb));
        geom.push_back(DecorativeSphere(0.90*getDefaultRadius())
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMyParentMobilizedBodyIndex())
            .setTransform(X_PMb));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(getMyMobilizedBodyIndex())
                                            .setTransform(X_BM));
    }
}

    ///////////////////////////////
    // MOBILIZED BODY::ELLIPSOID //
    ///////////////////////////////

MobilizedBody::Ellipsoid::Ellipsoid(Direction d) : MobilizedBody(new EllipsoidImpl(d)) {
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new EllipsoidImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body,      const Transform& outbFrame, Direction d)
:   MobilizedBody(new EllipsoidImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body,      const Transform& outbFrame,
                                    const Vec3& radii, Direction d)
:   MobilizedBody(new EllipsoidImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);
    setDefaultRadii(radii);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ellipsoid& MobilizedBody::Ellipsoid::setDefaultRadii(const Vec3& r) {
    updImpl().setDefaultRadii(r);
    return *this;
}

const Vec3& MobilizedBody::Ellipsoid::getDefaultRadii() const {
    return getImpl().getDefaultRadii();
}

const Quaternion& MobilizedBody::Ellipsoid::getDefaultQ() const {
    return getImpl().defaultQ;
}
Quaternion& MobilizedBody::Ellipsoid::updDefaultQ() {
    return updImpl().defaultQ;
}

const Vec4& MobilizedBody::Ellipsoid::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Ellipsoid::setQ(State& s, const Vec4& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    Vec4::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec4& MobilizedBody::Ellipsoid::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec4& MobilizedBody::Ellipsoid::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}

const Vec3& MobilizedBody::Ellipsoid::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Ellipsoid::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Ellipsoid::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec4& MobilizedBody::Ellipsoid::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Ellipsoid::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec4& MobilizedBody::Ellipsoid::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Ellipsoid::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

void MobilizedBody::EllipsoidImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
{
    // Call the superclass implementation to get standard default geometry.
    
    MobilizedBodyImpl::calcDecorativeGeometryAndAppendImpl(s, stage, geom);

    // We can't generate the ellipsoid until we know the radius, and we can't place either
    // piece of geometry on the bodies until we know the parent and child mobilizer frame
    // placements, which might not be until Instance stage.
    if (stage == Stage::Instance && getMyMatterSubsystemRep().getShowDefaultGeometry()) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Transform& X_PMb = getInboardFrame(s);
        const Transform& X_BM  = getOutboardFrame(s);

        //TODO: this should come from the State.
        const Vec3 radii = getDefaultRadii();

        // Put an ellipsoid on the parent, and some wires to make it easier to track.
        geom.push_back(DecorativeEllipsoid(radii)
            .setColor(Gray)
            .setRepresentation(DecorativeGeometry::DrawSurface)
            .setOpacity(0.5)
            .setResolution(1.25)
            .setBodyId(getMyParentMobilizedBodyIndex())
            .setTransform(X_PMb));
        geom.push_back(DecorativeEllipsoid(radii*.99)
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMyParentMobilizedBodyIndex())
            .setTransform(X_PMb));

        // Calculate the follower plate dimensions from the ellipsoid dimensions.
        const Real minr = std::min(radii[0],std::min(radii[1],radii[2]));
        const Real hw = minr/3;  // half width of follower plate in x
        const Real hh = minr/30; // half height of follower plate

        // Still on the inboard body draw, an orange mesh ellipsoid "trapping" 
        // the follower plate.
        /*
        geom.push_back(DecorativeEllipsoid(radii + 9*Vec3(hh))
            .setColor(Orange)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            //.setOpacity(.2)
            .setResolution(0.75)
            .setLineThickness(1)
            .setBodyId(getMyParentMobilizedBodyIndex())
            .setTransform(X_PMb));
        */

        // raise up so bottom is on xy plane
        const Transform X_BFollower(X_BM.R(), X_BM.p() + Vec3(0,0,hh));
        geom.push_back(DecorativeBrick(Vec3(hw,2*hw/3.,hh))
            .setColor(Orange)
            .setBodyId(getMyMobilizedBodyIndex())
            .setTransform(X_BFollower));
    }
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ellipsoid, MobilizedBody::EllipsoidImpl, MobilizedBody);

    /////////////////////////////////
    // MOBILIZED BODY::TRANSLATION //
    /////////////////////////////////

MobilizedBody::Translation::Translation(Direction d) : MobilizedBody(new TranslationImpl(d)) {
}


MobilizedBody::Translation::Translation(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new TranslationImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Translation::Translation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new TranslationImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

const Vec3& MobilizedBody::Translation::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::Translation& MobilizedBody::Translation::setDefaultQ(const Vec3& q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

const Vec3& MobilizedBody::Translation::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Translation::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::Translation::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::Translation::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Translation::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Translation::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Translation::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::Translation::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Translation::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::Translation::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Translation::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Translation, MobilizedBody::TranslationImpl, MobilizedBody);

    //////////////////////////
    // MOBILIZED BODY::FREE //
    //////////////////////////

MobilizedBody::Free::Free(Direction d) : MobilizedBody(new FreeImpl(d)) {
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new FreeImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new FreeImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Free& MobilizedBody::Free::setDefaultTranslation(const Vec3& p_FM) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQTranslation = p_FM;
    return *this;
}


MobilizedBody::Free& MobilizedBody::Free::setDefaultQuaternion(const Quaternion& R_FM) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQOrientation = R_FM;
    return *this;
}

MobilizedBody::Free& MobilizedBody::Free::setDefaultRotation(const Rotation& R_FM) {
    setDefaultQuaternion(R_FM.convertRotationToQuaternion());
    return *this;
}

MobilizedBody::Free& MobilizedBody::Free::setDefaultTransform(const Transform& X_FM) {
    setDefaultTranslation(X_FM.p());
    setDefaultQuaternion(X_FM.R().convertRotationToQuaternion());
    return *this;
}

const Vec3& MobilizedBody::Free::getDefaultTranslation() const {
    return getImpl().defaultQTranslation;
}

const Quaternion& MobilizedBody::Free::getDefaultQuaternion() const {
    return getImpl().defaultQOrientation;
}

const Vec7& MobilizedBody::Free::getDefaultQ() const {
    // assuming struct is packed so (Orientation,Translation) are contiguous
    return Vec7::getAs((const Real*)&getImpl().defaultQOrientation);
}
MobilizedBody::Free& MobilizedBody::Free::setDefaultQ(const Vec7& q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQOrientation = Quaternion(q.getSubVec<4>(0));
    updImpl().defaultQTranslation = q.getSubVec<3>(4);
    return *this;
}

const Vec7& MobilizedBody::Free::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Free::setQ(State& s, const Vec7& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    Vec7::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec7& MobilizedBody::Free::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec7& MobilizedBody::Free::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec6& MobilizedBody::Free::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Free::setU(State& s, const Vec6& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    Vec6::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec6& MobilizedBody::Free::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec7& MobilizedBody::Free::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 7);
    return Vec7::getAs(&qlike[qStart]);
}

const Vec6& MobilizedBody::Free::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&ulike[uStart]);
}

Vec7& MobilizedBody::Free::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 7);
    return Vec7::updAs(&qlike[qStart]);
}

Vec6& MobilizedBody::Free::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Free, MobilizedBody::FreeImpl, MobilizedBody);

    //////////////////////////////////////
    // MOBILIZED BODY::LINE ORIENTATION //
    //////////////////////////////////////

MobilizedBody::LineOrientation::LineOrientation(Direction d) : MobilizedBody(new LineOrientationImpl(d)) {
}


MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new LineOrientationImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new LineOrientationImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

const Quaternion& MobilizedBody::LineOrientation::getDefaultQ() const {
    return getImpl().defaultQ;
}
MobilizedBody::LineOrientation& MobilizedBody::LineOrientation::setDefaultQ(const Quaternion& quat) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = quat;
    return *this;
}

const Vec4& MobilizedBody::LineOrientation::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::LineOrientation::setQ(State& s, const Vec4& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    Vec4::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec4& MobilizedBody::LineOrientation::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec4& MobilizedBody::LineOrientation::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}

const Vec2& MobilizedBody::LineOrientation::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec2::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::LineOrientation::setU(State& s, const Vec2& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec2::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec2& MobilizedBody::LineOrientation::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec2::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec4& MobilizedBody::LineOrientation::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&qlike[qStart]);
}

const Vec2& MobilizedBody::LineOrientation::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec2::getAs(&ulike[uStart]);
}

Vec4& MobilizedBody::LineOrientation::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::updAs(&qlike[qStart]);
}

Vec2& MobilizedBody::LineOrientation::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec2::updAs(&ulike[uStart]);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::LineOrientation, MobilizedBody::LineOrientationImpl, MobilizedBody);

    ///////////////////////////////
    // MOBILIZED BODY::FREE LINE //
    ///////////////////////////////

MobilizedBody::FreeLine::FreeLine(Direction d) : MobilizedBody(new FreeLineImpl(d)) {
}


MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Body& body, Direction d) 
:   MobilizedBody(new FreeLineImpl(d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame, Direction d) 
:   MobilizedBody(new FreeLineImpl(d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}
SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::FreeLine, MobilizedBody::FreeLineImpl, MobilizedBody);

    //////////////////////////
    // MOBILIZED BODY::WELD //
    //////////////////////////

MobilizedBody::Weld::Weld() : MobilizedBody(new WeldImpl(MobilizedBody::Forward)) {
}


MobilizedBody::Weld::Weld(MobilizedBody& parent, const Body& body) 
:   MobilizedBody(new WeldImpl(MobilizedBody::Forward)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Weld::Weld(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) 
:   MobilizedBody(new WeldImpl(MobilizedBody::Forward)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Weld, MobilizedBody::WeldImpl, MobilizedBody);

    ////////////////////////////////
    // (IM)MOBILIZED BODY::GROUND //
    ////////////////////////////////

MobilizedBody::Ground::Ground() : MobilizedBody(new GroundImpl()) {
    setBody(Body::Ground());
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ground, MobilizedBody::GroundImpl, MobilizedBody);

    ///////////////////////////
    // MOBILIZED BODY::SCREW //
    ///////////////////////////

MobilizedBody::Screw::Screw(Real pitch, Direction d) : MobilizedBody(new ScrewImpl(pitch,d)) {
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Body& body, Real pitch, Direction d) 
:   MobilizedBody(new ScrewImpl(pitch,d)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Transform& inbFrame,
                            const Body& body, const Transform& outbFrame,
                            Real pitch, Direction d) 
:   MobilizedBody(new ScrewImpl(pitch,d)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Screw& MobilizedBody::Screw::setDefaultPitch(Real pitch) {
    updImpl().setDefaultPitch(pitch);
    return *this;
}

Real MobilizedBody::Screw::getDefaultPitch() const {
    return getImpl().getDefaultPitch();
}

Real MobilizedBody::Screw::getDefaultQ() const {
    return getImpl().defaultQ;
}

MobilizedBody::Screw& MobilizedBody::Screw::setDefaultQ(Real q) {
    getImpl().invalidateTopologyCache();
    updImpl().defaultQ = q;
    return *this;
}

Real MobilizedBody::Screw::getQ(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Screw::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Screw::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Screw::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    QIndex qStart; int nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Screw::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Screw::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Screw::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    UIndex uStart; int nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Screw::getMyPartQ(const State& s, const Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Screw::getMyPartU(const State& s, const Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Screw::updMyPartQ(const State& s, Vector& qlike) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Screw::updMyPartU(const State& s, Vector& ulike) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Screw, MobilizedBody::ScrewImpl, MobilizedBody);

////////////////////////////
// MOBILIZED BODY::CUSTOM //
////////////////////////////

// We are given an Implementation object which is already holding a CustomImpl
// object for us. We'll first take away ownership of the CustomImpl, then
// make the CustomImpl take over ownership of the Implementation object.
MobilizedBody::Custom::Custom
   (MobilizedBody& parent, MobilizedBody::Custom::Implementation* implementation, 
    const Body& body, Direction d)
:   MobilizedBody(implementation ? implementation->updImpl().removeOwnershipOfCustomImpl() : 0)
{
    SimTK_ASSERT_ALWAYS(implementation,
        "MobilizedBody::Custom::Custom(): Implementation pointer was NULL.");
    setBody(body);

    updImpl().setDirection(d);
    
    // Now store the Implementation pointer in our CustomImpl. The Implementation
    // object retains its original pointer to the CustomImpl object so it can
    // operate as a proxy for the CustomImpl. However the Custom handle now owns the
    // CustomImpl and the CustomImpl owns the Implementation.
    updImpl().takeOwnershipOfImplementation(implementation);
    
    updImpl().updMyMatterSubsystemRep().adoptMobilizedBody(parent.getMobilizedBodyIndex(), *this);
}

MobilizedBody::Custom::Custom
   (MobilizedBody& parent, MobilizedBody::Custom::Implementation* implementation, 
    const Transform& inbFrame, const Body& body, const Transform& outbFrame,
    Direction d)
:   MobilizedBody(implementation ? implementation->updImpl().removeOwnershipOfCustomImpl() : 0)
{
    SimTK_ASSERT_ALWAYS(implementation,
        "MobilizedBody::Custom::Custom(): Implementation pointer was NULL.");
    setBody(body);
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);

    updImpl().setDirection(d);
    
    // Now store the Implementation pointer in our CustomImpl. The Implementation
    // object retains its original pointer to the CustomImpl object so it can
    // operate as a proxy for the CustomImpl. However the Custom handle now owns the
    // CustomImpl and the CustomImpl owns the Implementation.
    updImpl().takeOwnershipOfImplementation(implementation);
    
    updImpl().updMyMatterSubsystemRep().adoptMobilizedBody(parent.getMobilizedBodyIndex(), *this);
}

const MobilizedBody::Custom::Implementation& MobilizedBody::Custom::getImplementation() const {
    return getImpl().getImplementation();
}

MobilizedBody::Custom::Implementation& MobilizedBody::Custom::updImplementation() {
    return updImpl().updImplementation();
}

SimTK_INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Custom, MobilizedBody::CustomImpl, MobilizedBody);

// MobilizedBody::CustomImpl

// The Implementation object should already contain a pointer to this CustomImpl object.
void MobilizedBody::CustomImpl::takeOwnershipOfImplementation(Custom::Implementation* userImpl) {
    assert(!implementation); // you can only do this once!
    assert(userImpl);
    const Custom::ImplementationImpl& impImpl = userImpl->getImpl();
    assert(&impImpl.getCustomImpl() == this && !impImpl.isOwnerOfCustomImpl());
    implementation = userImpl;
}  

////////////////////////////////////////////
// MOBILIZED BODY::CUSTOM::IMPLEMENTATION //
////////////////////////////////////////////

// We create the initial CustomImpl as though it the mobilizer will be in the forward direction.
// However, that may be changed later when this implementation is put to use.
MobilizedBody::Custom::Implementation::Implementation(SimbodyMatterSubsystem& matter, int nu, int nq, int nAngles) 
: PIMPLHandle<Implementation,ImplementationImpl>(new ImplementationImpl(new CustomImpl(MobilizedBody::Forward), nu, nq, nAngles)) {
    assert(nu > 0);
    assert(nq > 0);
    assert(nAngles >= 0);
    assert(nu <= 6);
    assert(nq <= 7);
    assert(nAngles <= 4);
    assert(nu <= nq);
    assert(nAngles <= nq);
    // We don't know the MobilizedBodyIndex yet since this hasn't been adopted by the MatterSubsystem.
    updImpl().updCustomImpl().setMyMatterSubsystem(matter, MobilizedBodyIndex(0), MobilizedBodyIndex(0));
}

MobilizedBody::Custom::Implementation::~Implementation() {
}

void MobilizedBody::Custom::Implementation::invalidateTopologyCache() const {
    getImpl().getCustomImpl().invalidateTopologyCache();
}

bool MobilizedBody::Custom::Implementation::getUseEulerAngles(const State& state) const {
    return getImpl().getCustomImpl().getMyMatterSubsystemRep().getUseEulerAngles(state);
}

Vector MobilizedBody::Custom::Implementation::getQ(const State& state) const {
    const SBModelVars& mv = getImpl().getCustomImpl().getMyMatterSubsystemRep().getModelVars(state);
    const RigidBodyNode& body = getImpl().getCustomImpl().getMyRigidBodyNode();
    const int indexBase = body.getQIndex();
    Vector q(body.getNQInUse(mv));
    for (int i = 0; i < q.size(); ++i) {
        int index = indexBase+i;
        q[i] = state.getQ()[index];
    }
    return q;
}

Vector MobilizedBody::Custom::Implementation::getU(const State& state) const {
    const SBModelVars& mv = getImpl().getCustomImpl().getMyMatterSubsystemRep().getModelVars(state);
    const RigidBodyNode& body = getImpl().getCustomImpl().getMyRigidBodyNode();
    const int indexBase = body.getUIndex();
    Vector u(body.getNUInUse(mv));
    for (int i = 0; i < u.size(); ++i) {
        int index = indexBase+i;
        u[i] = state.getU()[index];
    }
    return u;
}

Vector MobilizedBody::Custom::Implementation::getQDot(const State& state) const {
    const SBModelVars& mv = getImpl().getCustomImpl().getMyMatterSubsystemRep().getModelVars(state);
    const RigidBodyNode& body = getImpl().getCustomImpl().getMyRigidBodyNode();
    const int indexBase = body.getQIndex();
    Vector qdot(body.getNQInUse(mv));
    for (int i = 0; i < qdot.size(); ++i) {
        int index = indexBase+i;
        qdot[i] = state.getQDot()[index];
    }
    return qdot;
}

Vector MobilizedBody::Custom::Implementation::getUDot(const State& state) const {
    const SBModelVars& mv = getImpl().getCustomImpl().getMyMatterSubsystemRep().getModelVars(state);
    const RigidBodyNode& body = getImpl().getCustomImpl().getMyRigidBodyNode();
    const int indexBase = body.getUIndex();
    Vector udot(body.getNUInUse(mv));
    for (int i = 0; i < udot.size(); ++i) {
        int index = indexBase+i;
        udot[i] = state.getUDot()[index];
    }
    return udot;
}

Vector MobilizedBody::Custom::Implementation::getQDotDot(const State& state) const {
    const SBModelVars& mv = getImpl().getCustomImpl().getMyMatterSubsystemRep().getModelVars(state);
    const RigidBodyNode& body = getImpl().getCustomImpl().getMyRigidBodyNode();
    const int indexBase = body.getQIndex();
    Vector qdotdot(body.getNQInUse(mv));
    for (int i = 0; i < qdotdot.size(); ++i) {
        int index = indexBase+i;
        qdotdot[i] = state.getQDotDot()[index];
    }
    return qdotdot;
}

// Careful: must return the transform X_F0M0 using the "as defined" frames, rather than
// X_FM which might be reversed due to mobilizer reversal.
Transform MobilizedBody::Custom::Implementation::getMobilizerTransform(const State& s) const {
    // Use "upd" instead of "get" here because the custom mobilizer definition needs access to
    // this local information during realizePosition() so get would throw a stage violation.
    const SBTreePositionCache& pc = getImpl().getCustomImpl().getMyMatterSubsystemRep().updTreePositionCache(s);
    const RigidBodyNode& node = getImpl().getCustomImpl().getMyRigidBodyNode();
    return node.findX_F0M0(pc);
}

// Careful: must return the velocity V_F0M0 using the "as defined" frames, rather than
// V_FM which might be reversed due to mobilizer reversal.
SpatialVec MobilizedBody::Custom::Implementation::getMobilizerVelocity(const State& s) const {
    const SBTreePositionCache& pc = getImpl().getCustomImpl().getMyMatterSubsystemRep().getTreePositionCache(s);
    // Use "upd" instead of "get" here because the custom mobilizer definition needs access to
    // this local information during realizeVelocity() so get would throw a stage violation.
    const SBTreeVelocityCache& vc = getImpl().getCustomImpl().getMyMatterSubsystemRep().updTreeVelocityCache(s);
    const RigidBodyNode& node = getImpl().getCustomImpl().getMyRigidBodyNode();
    return node.findV_F0M0(pc,vc);
}

void MobilizedBody::Custom::Implementation::multiplyByN(const State& s, bool transposeMatrix, 
                       int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNumAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = in[i];
}

void MobilizedBody::Custom::Implementation::multiplyByNInv(const State& s, bool transposeMatrix, 
                                       int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNumAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = in[i];
}

void MobilizedBody::Custom::Implementation::multiplyByNDot(const State& s, bool transposeMatrix, 
                                         int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNumAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = 0;
}

void MobilizedBody::Custom::Implementation::setQToFitTransform(const State& state, const Transform& X_FM, int nq, Real* q) const {
    class OptimizerFunction : public OptimizerSystem {
    public:
        OptimizerFunction(const MobilizedBody::Custom::Implementation& impl, const State& state, int nq, const Transform& X_FM) :
                OptimizerSystem(nq), impl(impl), state(state), X_FM(X_FM) {
        }
        int objectiveFunc(const Vector& parameters, const bool new_parameters, Real& f) const {
            Transform transform = impl.calcMobilizerTransformFromQ(state, parameters.size(), &parameters[0]);
            f = (transform.p()-X_FM.p()).norm();
            f += std::abs((~transform.R()*X_FM.R()).convertRotationToAngleAxis()[0]);
            return 0;
        }
    private:
        const MobilizedBody::Custom::Implementation& impl;
        const State& state;
        const Transform& X_FM;
    };
    OptimizerFunction function(*this, state, nq, X_FM);
    Optimizer opt(function);
    opt.useNumericalJacobian(true);
    opt.useNumericalGradient(true);
    opt.setLimitedMemoryHistory(100);
    Vector qvec(nq);
    
    // Pick initiial values which are 1) deterministic and 2) unlikely to correspond to a local
    // maximum or inflection point, which could cause the optimizer to fail.
    
    for (int i = 0; i < nq; i++)
        qvec[i] = i+0.12354;
    opt.optimize(qvec);
    for (int i = 0; i < nq; i++)
        q[i] = qvec[i];
}

void MobilizedBody::Custom::Implementation::setUToFitVelocity(const State& state, const SpatialVec& V_FM, int nu, Real* u) const {
    class OptimizerFunction : public OptimizerSystem {
    public:
        OptimizerFunction(const MobilizedBody::Custom::Implementation& impl, const State& state, int nu, const SpatialVec& V_FM) :
                OptimizerSystem(nu), impl(impl), state(state), V_FM(V_FM) {
        }
        int objectiveFunc(const Vector& parameters, const bool new_parameters, Real& f) const {
            SpatialVec v = impl.multiplyByHMatrix(state, parameters.size(), &parameters[0]);
            f = (v[0]-V_FM[0]).norm();
            f += (v[1]-V_FM[1]).norm();
            return 0;
        }
    private:
        const MobilizedBody::Custom::Implementation& impl;
        const State& state;
        const SpatialVec& V_FM;
    };
    OptimizerFunction function(*this, state, nu, V_FM);
    Optimizer opt(function);
    opt.useNumericalJacobian(true);
    opt.useNumericalGradient(true);
    opt.setLimitedMemoryHistory(100);
    Vector uvec(nu);
    
    // Pick initiial values which are 1) deterministic and 2) unlikely to correspond to a local
    // maximum or inflection point, which could cause the optimizer to fail.
    
    for (int i = 0; i < nu; i++)
        uvec[i] = i+0.12354;
    opt.optimize(uvec);
    for (int i = 0; i < nu; i++)
        u[i] = uvec[i];
}

// Constructors without user-specified axes for function-based mobilized body
MobilizedBody::FunctionBased::FunctionBased
   (MobilizedBody& parent, const Body& body, 
    int nmobilities, const Array_<const Function*>& functions,
    const Array_<Array_<int> >& coordIndices,
    Direction direction)
:   Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices), body, direction) 
{
}

MobilizedBody::FunctionBased::FunctionBased
   (MobilizedBody& parent, const Transform& inbFrame, 
    const Body& body, const Transform& outbFrame, 
    int nmobilities, const Array_<const Function*>& functions,
    const Array_<Array_<int> >& coordIndices,
    Direction direction)
:   Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices), body, direction) 
{
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
}

// Constructors that allow user-specified axes for function-based mobilized body
MobilizedBody::FunctionBased::FunctionBased
   (MobilizedBody& parent, const Body& body, 
    int nmobilities, const Array_<const Function*>& functions,
    const Array_<Array_<int> >& coordIndices, const Array_<Vec3>& axes,
    Direction direction)
:   Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices, axes), body, direction) {
}

MobilizedBody::FunctionBased::FunctionBased
   (MobilizedBody& parent, const Transform& inbFrame, 
    const Body& body, const Transform& outbFrame, 
    int nmobilities, const Array_<const Function*>& functions,
    const Array_<Array_<int> >& coordIndices, const Array_<Vec3>& axes,
    Direction direction)
:   Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices, axes), body, direction) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
}

} // namespace SimTK

