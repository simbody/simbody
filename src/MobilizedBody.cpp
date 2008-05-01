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
void MobilizedBody::setQToFitTranslation(State& s, const Vec3& T_MbM) const { 
    getImpl().setQToFitTranslation(s,T_MbM);
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

int MobilizedBody::getNumQ(const State& s) const {
    QIndex qStart; int nq; getImpl().findMobilizerQs(s, qStart, nq);
    return nq;
}

int MobilizedBody::getNumU(const State& s) const {
    UIndex uStart; int nu; getImpl().findMobilizerUs(s, uStart, nu);
    return nu;
}

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


void MobilizedBody::applyBodyForce(const State& s, const SpatialVec& spatialForceInG, 
                                   Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
    bodyForces[getMobilizedBodyIndex()] += spatialForceInG;
}

void MobilizedBody::applyBodyTorque(const State& s, const Vec3& torqueInG, 
                     Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
    bodyForces[getMobilizedBodyIndex()][0] += torqueInG; // don't change force
}

void MobilizedBody::applyForceToBodyPoint(const State& s, const Vec3& pointInB, const Vec3& forceInG,
                           Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
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
    copyBody = *this;
    copyBody.updImpl().myMatterSubsystemRep = 0;
    copyBody.updImpl().myRBnode = 0;
    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(), copyBody);
    return parent.updMatterSubsystem().updMobilizedBody(copyBody.getMobilizedBodyIndex());
}

    ////////////////////////
    // MOBILIZED BODY REP //
    ////////////////////////

void MobilizedBodyImpl::findMobilizerQs(const State& s, QIndex& qStart, int& nq) const {
    getMyMatterSubsystemRep()
        .findMobilizerQs(s, myMobilizedBodyIndex, qStart, nq);
}
void MobilizedBodyImpl::findMobilizerUs(const State& s, UIndex& uStart, int& nu) const {
    getMyMatterSubsystemRep()
        .findMobilizerUs(s, myMobilizedBodyIndex, uStart, nu);
}

void MobilizedBodyImpl::copyOutDefaultQ(const State& s, Vector& qDefault) const {
    SimTK_STAGECHECK_GE_ALWAYS(getMyMatterSubsystemRep().getStage(s), Stage::Topology,
        "MobilizedBody::copyOutDefaultQ()");
    QIndex qStart; int nq;
    findMobilizerQs(s, qStart, nq);
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
void MobilizedBodyImpl::setQToFitTranslation(State& s, const Vec3& T_MbM) const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBStateDigest digest(s, matterRep, Stage::Instance);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitTranslation(digest, T_MbM, q);
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
    realizeTopologyImpl(s);
    return *myRBnode;
}

    /////////////////////////
    // MOBILIZED BODY::PIN //
    /////////////////////////

MobilizedBody::Pin::Pin() : MobilizedBody(new PinImpl()) {
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Body& body) : MobilizedBody(new PinImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Transform& inbFrame,
                        const Body& body, const Transform& outbFrame) : MobilizedBody(new PinImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Pin, MobilizedBody::PinImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::SLIDER //
    ////////////////////////////

MobilizedBody::Slider::Slider() : MobilizedBody(new SliderImpl()) {
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Body& body) : MobilizedBody(new SliderImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Transform& inbFrame,
                        const Body& body, const Transform& outbFrame) : MobilizedBody(new SliderImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Slider, MobilizedBody::SliderImpl, MobilizedBody);

    ///////////////////////////////
    // MOBILIZED BODY::UNIVERSAL //
    ///////////////////////////////

MobilizedBody::Universal::Universal() : MobilizedBody(new UniversalImpl()) {
}


MobilizedBody::Universal::Universal(MobilizedBody& parent, const Body& body) : MobilizedBody(new UniversalImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Universal::Universal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame) : MobilizedBody(new UniversalImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Universal, MobilizedBody::UniversalImpl, MobilizedBody);

    //////////////////////////////
    // MOBILIZED BODY::CYLINDER //
    //////////////////////////////

MobilizedBody::Cylinder::Cylinder() : MobilizedBody(new CylinderImpl()) {
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Body& body) : MobilizedBody(new CylinderImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : MobilizedBody(new CylinderImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Cylinder, MobilizedBody::CylinderImpl, MobilizedBody);

    //////////////////////////////////
    // MOBILIZED BODY::BEND STRETCH //
    //////////////////////////////////

MobilizedBody::BendStretch::BendStretch() : MobilizedBody(new BendStretchImpl()) {
}


MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Body& body) : MobilizedBody(new BendStretchImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                                        const Body& body, const Transform& outbFrame) : MobilizedBody(new BendStretchImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::BendStretch, MobilizedBody::BendStretchImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::PLANAR //
    ////////////////////////////

MobilizedBody::Planar::Planar() : MobilizedBody(new PlanarImpl()) {
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Body& body) : MobilizedBody(new PlanarImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : MobilizedBody(new PlanarImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Planar, MobilizedBody::PlanarImpl, MobilizedBody);

    ////////////////////////////
    // MOBILIZED BODY::GIMBAL //
    ////////////////////////////

MobilizedBody::Gimbal::Gimbal() : MobilizedBody(new GimbalImpl()) {
}


MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Body& body) : MobilizedBody(new GimbalImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame) : MobilizedBody(new GimbalImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Gimbal, MobilizedBody::GimbalImpl, MobilizedBody);

    // GimbalImpl

void MobilizedBody::GimbalImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
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

MobilizedBody::Ball::Ball() : MobilizedBody(new BallImpl()) {
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Body& body) : MobilizedBody(new BallImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : MobilizedBody(new BallImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ball, MobilizedBody::BallImpl, MobilizedBody);

    // BallImpl

void MobilizedBody::BallImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
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

MobilizedBody::Ellipsoid::Ellipsoid() : MobilizedBody(new EllipsoidImpl()) {
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Body& body) : MobilizedBody(new EllipsoidImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : MobilizedBody(new EllipsoidImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

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

void MobilizedBody::EllipsoidImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
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
        const Transform X_BFollower(X_BM.R(), X_BM.T() + Vec3(0,0,hh));
        geom.push_back(DecorativeBrick(Vec3(hw,2*hw/3.,hh))
            .setColor(Orange)
            .setBodyId(getMyMobilizedBodyIndex())
            .setTransform(X_BFollower));
    }
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ellipsoid, MobilizedBody::EllipsoidImpl, MobilizedBody);

    /////////////////////////////////
    // MOBILIZED BODY::TRANSLATION //
    /////////////////////////////////

MobilizedBody::Translation::Translation() : MobilizedBody(new TranslationImpl()) {
}


MobilizedBody::Translation::Translation(MobilizedBody& parent, const Body& body) : MobilizedBody(new TranslationImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Translation::Translation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : MobilizedBody(new TranslationImpl()) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Translation, MobilizedBody::TranslationImpl, MobilizedBody);

    //////////////////////////
    // MOBILIZED BODY::FREE //
    //////////////////////////

MobilizedBody::Free::Free() : MobilizedBody(new FreeImpl()) {
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Body& body) : MobilizedBody(new FreeImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : MobilizedBody(new FreeImpl()) {
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
    setDefaultTranslation(X_FM.T());
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Free, MobilizedBody::FreeImpl, MobilizedBody);

    //////////////////////////////////////
    // MOBILIZED BODY::LINE ORIENTATION //
    //////////////////////////////////////

MobilizedBody::LineOrientation::LineOrientation() : MobilizedBody(new LineOrientationImpl()) {
}


MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Body& body) : MobilizedBody(new LineOrientationImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : MobilizedBody(new LineOrientationImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::LineOrientation, MobilizedBody::LineOrientationImpl, MobilizedBody);

    ///////////////////////////////
    // MOBILIZED BODY::FREE LINE //
    ///////////////////////////////

MobilizedBody::FreeLine::FreeLine() : MobilizedBody(new FreeLineImpl()) {
}


MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Body& body) : MobilizedBody(new FreeLineImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : MobilizedBody(new FreeLineImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}
INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::FreeLine, MobilizedBody::FreeLineImpl, MobilizedBody);

    //////////////////////////
    // MOBILIZED BODY::WELD //
    //////////////////////////

MobilizedBody::Weld::Weld() : MobilizedBody(new WeldImpl()) {
}


MobilizedBody::Weld::Weld(MobilizedBody& parent, const Body& body) : MobilizedBody(new WeldImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Weld::Weld(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : MobilizedBody(new WeldImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Weld, MobilizedBody::WeldImpl, MobilizedBody);

    ////////////////////////////////
    // (IM)MOBILIZED BODY::GROUND //
    ////////////////////////////////

MobilizedBody::Ground::Ground() : MobilizedBody(new GroundImpl()) {
    setBody(Body::Ground());
}

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Ground, MobilizedBody::GroundImpl, MobilizedBody);

    ///////////////////////////
    // MOBILIZED BODY::SCREW //
    ///////////////////////////

MobilizedBody::Screw::Screw(Real pitch) : MobilizedBody(new ScrewImpl(pitch)) {
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Body& body, Real pitch) : MobilizedBody(new ScrewImpl(pitch)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyIndex(),
                                                   *this);
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Transform& inbFrame,
                            const Body& body, const Transform& outbFrame,
                            Real pitch) : MobilizedBody(new ScrewImpl(pitch)) {
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Screw, MobilizedBody::ScrewImpl, MobilizedBody);

////////////////////////////
// MOBILIZED BODY::CUSTOM //
////////////////////////////

// We are given an Implementation object which is already holding a CustomImpl
// object for us. We'll first take away ownership of the CustomImpl, then
// make the CustomImpl take over ownership of the Implementation object.
MobilizedBody::Custom::Custom(MobilizedBody& parent, MobilizedBody::Custom::Implementation* implementation, const Body& body)
    : MobilizedBody(implementation ? implementation->updImpl().removeOwnershipOfCustomImpl() : 0)
{
    SimTK_ASSERT_ALWAYS(implementation,
        "MobilizedBody::Custom::Custom(): Implementation pointer was NULL.");
    setBody(body);
    
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

INSERT_DERIVED_HANDLE_DEFINITIONS(MobilizedBody::Custom, MobilizedBody::CustomImpl, MobilizedBody);

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

MobilizedBody::Custom::Implementation::Implementation(SimbodyMatterSubsystem& matter, int nu, int nq, int nAngles) 
: PIMPLHandle<Implementation,ImplementationImpl>(new ImplementationImpl(new CustomImpl(), nu, nq, nAngles)) {
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

const Transform&  MobilizedBody::Custom::Implementation::getMobilizerTransform(const State& s) const {
    const SBPositionCache& pc = getImpl().getCustomImpl().getMyMatterSubsystemRep().getPositionCache(s, s.getSystemStage() < Stage::Position);
    const RigidBodyNode& node = getImpl().getCustomImpl().getMyRigidBodyNode();
    return node.getX_FM(pc);
}

const SpatialVec& MobilizedBody::Custom::Implementation::getMobilizerVelocity(const State& s) const {
    return getImpl().getCustomImpl().getMobilizerVelocity(s);
}

void MobilizedBody::Custom::Implementation::multiplyByQMatrix(const State& s, bool transposeMatrix, 
                       int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = in[i];
}

void MobilizedBody::Custom::Implementation::multiplyByQInverse(const State& s, bool transposeMatrix, 
                                       int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = in[i];
}

void MobilizedBody::Custom::Implementation::multiplyByQDotMatrix(const State& s, bool transposeMatrix, 
                                         int nIn, const Real* in, int nOut, Real* out) const {
    // Default implementation
    assert((nIn==0 || in) && (nOut==0 || out));
    int nu = getImpl().getNU(), nq = getImpl().getNQ(), nAngles = getImpl().getNAngles();
    assert(nq==nu && nIn==nu && nOut==nu && nAngles < 4);
    for (int i=0; i<nu; ++i) out[i] = 0;
}

MobilizedBody::FunctionBased::FunctionBased(MobilizedBody& parent, const Body& body, int nmobilities, const std::vector<Function<1>*>& functions, const std::vector<std::vector<int> >& coordIndices)
        : Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices), body) {
}

MobilizedBody::FunctionBased::FunctionBased(MobilizedBody& parent, const Transform& inbFrame, const Body& body, const Transform& outbFrame, int nmobilities, const std::vector<Function<1>*>& functions, const std::vector<std::vector<int> >& coordIndices)
        : Custom(parent, new FunctionBasedImpl(parent.updMatterSubsystem(), nmobilities, functions, coordIndices), body) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
}

} // namespace SimTK

