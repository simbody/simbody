/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

// This suppresses the 'extern template' instantiations in MobilizedBody.h so that
// we can instantiate them for real here.
#define SimTK_DEFINING_MOBILIZED_BODY

#include "SimTKcommon.h"
#include "SimTKcommon/internal/PrivateImplementation_Defs.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"

#include "MobilizedBodyImpl.h"
#include "SimbodyMatterSubsystemRep.h"

namespace SimTK {

template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl>;
template class PIMPLImplementation<MobilizedBody, MobilizedBodyImpl>;
template class PIMPLDerivedHandle<MobilizedBody::Ball, MobilizedBody::BallImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::BendStretch, MobilizedBody::BendStretchImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Custom, MobilizedBody::CustomImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Cylinder, MobilizedBody::CylinderImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Ellipsoid, MobilizedBody::EllipsoidImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Free, MobilizedBody::FreeImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::FreeLine, MobilizedBody::FreeLineImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Gimbal, MobilizedBody::GimbalImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Ground, MobilizedBody::GroundImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::LineOrientation, MobilizedBody::LineOrientationImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Pin, MobilizedBody::PinImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Planar, MobilizedBody::PlanarImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Screw, MobilizedBody::ScrewImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Slider, MobilizedBody::SliderImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Translation, MobilizedBody::TranslationImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Universal, MobilizedBody::UniversalImpl, MobilizedBody>;
template class PIMPLDerivedHandle<MobilizedBody::Weld, MobilizedBody::WeldImpl, MobilizedBody>;

    ////////////////////
    // MOBILIZED BODY //
    ////////////////////

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


MobilizedBodyId MobilizedBody::getMobilizedBodyId() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMobilizedBodyId() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().myMobilizedBodyId;
}

const MobilizedBody& MobilizedBody::getParentMobilizedBody() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getParentMobilizedBody() called on a MobilizedBody that is not part of a subsystem.");
    return getImpl().getMyMatterSubsystemRep().getMobilizedBody(getImpl().myParentId);
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
    return getImpl().getMyRigidBodyNode().getLevel();
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
    getImpl().setQToFitTranslation(s,T_MbM,false); // allow rotation
}
void MobilizedBody::setQToFitTranslationOnly(State& s, const Vec3& T_MbM) const { 
    getImpl().setQToFitTranslation(s,T_MbM,true);  // prevent rotation
}

void MobilizedBody::setUToFitVelocity(State& s, const SpatialVec& V_MbM) const { 
    getImpl().setUToFitVelocity(s,V_MbM);
}
void MobilizedBody::setUToFitAngularVelocity(State& s, const Vec3& w_MbM) const { 
    getImpl().setUToFitAngularVelocity(s,w_MbM);
}
void MobilizedBody::setUToFitLinearVelocity(State& s, const Vec3& v_MbM) const { 
    getImpl().setUToFitLinearVelocity(s,v_MbM,false); // allow angular velocity change
}
void MobilizedBody::setUToFitLinearVelocityOnly(State& s, const Vec3& v_MbM) const { 
    getImpl().setUToFitLinearVelocity(s,v_MbM,true);  // prevent angular velocity change
}

int MobilizedBody::getNumQ(const State& s) const {
    int qStart, nq; getImpl().findMobilizerQs(s, qStart, nq);
    return nq;
}

int MobilizedBody::getNumU(const State& s) const {
    int uStart, nu; getImpl().findMobilizerUs(s, uStart, nu);
    return nu;
}

Real  MobilizedBody::getOneFromQPartition(const State& s, int which, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s, qStart, nq);
    assert(0 <= which && which < nq);
    return qlike[qStart+which];
}
Real& MobilizedBody::updOneFromQPartition(const State& s, int which, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s, qStart, nq);
    assert(0 <= which && which < nq);
    return qlike[qStart+which];
}

Real  MobilizedBody::getOneFromUPartition(const State& s, int which, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s, uStart, nu);
    assert(0 <= which && which < nu);
    return ulike[uStart+which];
}
Real& MobilizedBody::updOneFromUPartition(const State& s, int which, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s, uStart, nu);
    assert(0 <= which && which < nu);
    return ulike[uStart+which];
}


void MobilizedBody::applyBodyForce(const State& s, const SpatialVec& spatialForceInG, 
                                   Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
    bodyForces[getMobilizedBodyId()] += spatialForceInG;
}

void MobilizedBody::applyBodyTorque(const State& s, const Vec3& torqueInG, 
                     Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
    bodyForces[getMobilizedBodyId()][0] += torqueInG; // don't change force
}

void MobilizedBody::applyForceToBodyPoint(const State& s, const Vec3& pointInB, const Vec3& forceInG,
                           Vector_<SpatialVec>& bodyForces) const 
{
    assert(bodyForces.size() == getMatterSubsystem().getNBodies());
    const Rotation& R_GB = getBodyTransform(s).R();
    bodyForces[getMobilizedBodyId()] += SpatialVec((R_GB*pointInB) % forceInG, forceInG);
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

Vector MobilizedBody::getQVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQ(s)(qStart,nq);
}

void MobilizedBody::setQVector(State& s, const Vector& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s, qStart, nq);
    assert(q.size() == nq);
    mbr.getMyMatterSubsystemRep().updQ(s)(qStart,nq) = q;
}

Vector MobilizedBody::getUVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s, uStart, nu);
    return mbr.getMyMatterSubsystemRep().getU(s)(uStart,nu);
}

void MobilizedBody::setUVector(State& s, const Vector& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s, uStart, nu);
    assert(u.size() == nu);
    mbr.getMyMatterSubsystemRep().updU(s)(uStart,nu) = u;
}

Vector MobilizedBody::getQDotVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQDot(s)(qStart,nq);
}

Vector MobilizedBody::getUDotVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s, uStart, nu);
    return mbr.getMyMatterSubsystemRep().getUDot(s)(uStart,nu);
}

Vector MobilizedBody::getQDotDotVector(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s, qStart, nq);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)(qStart,nq);
}

MobilizedBody* MobilizedBody::cloneForNewParent(MobilizedBody& parent) const {
    MobilizedBody* copyBody = new MobilizedBody();
    *copyBody = *this;
    copyBody->updImpl().myMatterSubsystemRep = 0;
    copyBody->updImpl().myRBnode = 0;
    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(), *copyBody);
    return copyBody;
}

    ////////////////////////
    // MOBILIZED BODY REP //
    ////////////////////////

void MobilizedBodyImpl::findMobilizerQs(const State& s, int& qStart, int& nq) const {
    getMyMatterSubsystemRep()
        .findMobilizerQs(s, myMobilizedBodyId, qStart, nq);
}
void MobilizedBodyImpl::findMobilizerUs(const State& s, int& uStart, int& nu) const {
    getMyMatterSubsystemRep()
        .findMobilizerUs(s, myMobilizedBodyId, uStart, nu);
}

void MobilizedBodyImpl::copyOutDefaultQ(const State& s, Vector& qDefault) const {
    SimTK_STAGECHECK_GE_ALWAYS(getMyMatterSubsystemRep().getStage(s), Stage::Topology,
        "MobilizedBody::copyOutDefaultQ()");
    int qStart, nq;
    findMobilizerQs(s, qStart, nq);
    copyOutDefaultQImpl(nq, &qDefault[qStart]);
}

    // TODO: currently we delegate these requests to the RigidBodyNodes. 
    // Probably most of this functionality should be handled directly
    // by the MobilizedBody objects.

void MobilizedBodyImpl::setQToFitTransform(State& s, const Transform& X_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitTransform(mv, X_MbM, q);
}
void MobilizedBodyImpl::setQToFitRotation(State& s, const Rotation& R_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitRotation(mv, R_MbM, q);
}
void MobilizedBodyImpl::setQToFitTranslation(State& s, const Vec3& T_MbM, 
                             bool dontChangeOrientation) const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    Vector& q = matterRep.updQ(s);
    return getMyRigidBodyNode().setQToFitTranslation(mv, T_MbM, q, dontChangeOrientation);
}

void MobilizedBodyImpl::setUToFitVelocity(State& s, const SpatialVec& V_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitVelocity(mv, q, V_MbM, u);
}
void MobilizedBodyImpl::setUToFitAngularVelocity(State& s, const Vec3& w_MbM) const {
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitAngularVelocity(mv, q, w_MbM, u);
}
void MobilizedBodyImpl::setUToFitLinearVelocity(State& s, const Vec3& v_MbM,
                                bool dontChangeAngularVelocity)  const
{
    const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
    const SBModelVars& mv = matterRep.getModelVars(s);
    const Vector& q = matterRep.updQ(s);
    Vector&       u = matterRep.updU(s);
    return getMyRigidBodyNode().setUToFitLinearVelocity(mv, q, v_MbM, u, dontChangeAngularVelocity);
}


const RigidBodyNode& MobilizedBodyImpl::realizeTopology
   (int& nxtU, int& nxtUSq, int& nxtQ) const
{
    delete myRBnode;
    myRBnode = createRigidBodyNode(nxtU,nxtUSq,nxtQ);

    int level;
    if (!myParentId.isValid()) {
        // this is ground
        assert(myMobilizedBodyId == 0);
        level = 0;
    } else {
        // not ground
        const MobilizedBodyImpl& parent = 
            myMatterSubsystemRep->getMobilizedBody(myParentId).getImpl();
        level = parent.myRBnode->getLevel() + 1;
        parent.myRBnode->addChild(myRBnode);
        myRBnode->setParent(parent.myRBnode);
    }

    myRBnode->setLevel(level);
    myRBnode->setNodeNum(myMobilizedBodyId);
    return *myRBnode;
}

    /////////////////////////
    // MOBILIZED BODY::PIN //
    /////////////////////////

MobilizedBody::Pin::Pin() : PIMPLDerivedHandleBase(new PinImpl()) {
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new PinImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Pin::Pin(MobilizedBody& parent, const Transform& inbFrame,
                        const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new PinImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Pin::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Pin::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Pin::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Pin::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Pin::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Pin::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Pin::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Pin::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Pin::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Pin::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

    ////////////////////////////
    // MOBILIZED BODY::SLIDER //
    ////////////////////////////

MobilizedBody::Slider::Slider() : PIMPLDerivedHandleBase(new SliderImpl()) {
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new SliderImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Slider::Slider(MobilizedBody& parent, const Transform& inbFrame,
                        const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new SliderImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Slider::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Slider::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Slider::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Slider::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Slider::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Slider::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Slider::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Slider::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Slider::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Slider::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

    ///////////////////////////////
    // MOBILIZED BODY::UNIVERSAL //
    ///////////////////////////////

MobilizedBody::Universal::Universal() : PIMPLDerivedHandleBase(new UniversalImpl()) {
}


MobilizedBody::Universal::Universal(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new UniversalImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Universal::Universal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new UniversalImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    //////////////////////////////
    // MOBILIZED BODY::CYLINDER //
    //////////////////////////////

MobilizedBody::Cylinder::Cylinder() : PIMPLDerivedHandleBase(new CylinderImpl()) {
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new CylinderImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Cylinder::Cylinder(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new CylinderImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    //////////////////////////////////
    // MOBILIZED BODY::BEND STRETCH //
    //////////////////////////////////

MobilizedBody::BendStretch::BendStretch() : PIMPLDerivedHandleBase(new BendStretchImpl()) {
}


MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new BendStretchImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::BendStretch::BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                                        const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new BendStretchImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    ////////////////////////////
    // MOBILIZED BODY::PLANAR //
    ////////////////////////////

MobilizedBody::Planar::Planar() : PIMPLDerivedHandleBase(new PlanarImpl()) {
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new PlanarImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Planar::Planar(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new PlanarImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Planar::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::Planar::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::Planar::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Planar::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Planar::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Planar::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::Planar::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Planar::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::Planar::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Planar::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

    ////////////////////////////
    // MOBILIZED BODY::GIMBAL //
    ////////////////////////////

MobilizedBody::Gimbal::Gimbal() : PIMPLDerivedHandleBase(new GimbalImpl()) {
}


MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new GimbalImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Gimbal::Gimbal(MobilizedBody& parent, const Transform& inbFrame,
                                    const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new GimbalImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    ///////////////////////////////////////////////////
    // MOBILIZED BODY::BALL (ORIENTATION, SPHERICAL) //
    ///////////////////////////////////////////////////

MobilizedBody::Ball::Ball() : PIMPLDerivedHandleBase(new BallImpl()) {
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new BallImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Ball::Ball(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new BallImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Ball::setQ(State& s, const Vec4& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    Vec4::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec4& MobilizedBody::Ball::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec4& MobilizedBody::Ball::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Ball::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Ball::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Ball::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec4& MobilizedBody::Ball::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Ball::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec4& MobilizedBody::Ball::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 4);
    return Vec4::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Ball::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

    // BallImpl

void MobilizedBody::BallImpl::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the parent and child mobilizer frame
    // placement on the body, which might not be until Instance stage.
    if (stage == Stage::Instance) {
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
                                            .setBodyId(getMyParentMobilizedBodyId())
                                            .setTransform(X_PMb));
        geom.push_back(DecorativeSphere(0.90*getDefaultRadius())
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMyParentMobilizedBodyId())
            .setTransform(X_PMb));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(getMyMobilizedBodyId())
                                            .setTransform(X_BM));
    }
}

    ///////////////////////////////
    // MOBILIZED BODY::ELLIPSOID //
    ///////////////////////////////

MobilizedBody::Ellipsoid::Ellipsoid() : PIMPLDerivedHandleBase(new EllipsoidImpl()) {
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new EllipsoidImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Ellipsoid::Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new EllipsoidImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    // We can't generate the ellipsoid until we know the radius, and we can't place either
    // piece of geometry on the bodies until we know the parent and child mobilizer frame
    // placements, which might not be until Instance stage.
    if (stage == Stage::Instance) {
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
            .setBodyId(getMyParentMobilizedBodyId())
            .setTransform(X_PMb));
        geom.push_back(DecorativeEllipsoid(radii*.99)
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(getMyParentMobilizedBodyId())
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
            .setBodyId(getMyParentMobilizedBodyId())
            .setTransform(X_PMb));
        */

        // raise up so bottom is on xy plane
        const Transform X_BFollower(X_BM.R(), X_BM.T() + Vec3(0,0,hh));
        geom.push_back(DecorativeBrick(Vec3(hw,2*hw/3.,hh))
            .setColor(Orange)
            .setBodyId(getMyMobilizedBodyId())
            .setTransform(X_BFollower));
    }
}

    /////////////////////////////////
    // MOBILIZED BODY::TRANSLATION //
    /////////////////////////////////

MobilizedBody::Translation::Translation() : PIMPLDerivedHandleBase(new TranslationImpl()) {
}


MobilizedBody::Translation::Translation(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new TranslationImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Translation::Translation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new TranslationImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Translation::setQ(State& s, const Vec3& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec3& MobilizedBody::Translation::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec3& MobilizedBody::Translation::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec3& MobilizedBody::Translation::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Translation::setU(State& s, const Vec3& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    Vec3::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec3& MobilizedBody::Translation::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec3& MobilizedBody::Translation::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::getAs(&qlike[qStart]);
}

const Vec3& MobilizedBody::Translation::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::getAs(&ulike[uStart]);
}

Vec3& MobilizedBody::Translation::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 3);
    return Vec3::updAs(&qlike[qStart]);
}

Vec3& MobilizedBody::Translation::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 3);
    return Vec3::updAs(&ulike[uStart]);
}

    //////////////////////////
    // MOBILIZED BODY::FREE //
    //////////////////////////

MobilizedBody::Free::Free() : PIMPLDerivedHandleBase(new FreeImpl()) {
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new FreeImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Free::Free(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new FreeImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQ(s)[qStart]);
}
void MobilizedBody::Free::setQ(State& s, const Vec7& q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    Vec7::updAs(&mbr.getMyMatterSubsystemRep().updQ(s)[qStart]) = q;
}
const Vec7& MobilizedBody::Free::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQDot(s)[qStart]);
}
const Vec7& MobilizedBody::Free::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq==7);
    return Vec7::getAs(&mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart]);
}


const Vec6& MobilizedBody::Free::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&mbr.getMyMatterSubsystemRep().getU(s)[uStart]);
}
void MobilizedBody::Free::setU(State& s, const Vec6& u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    Vec6::updAs(&mbr.getMyMatterSubsystemRep().updU(s)[uStart]) = u;
}
const Vec6& MobilizedBody::Free::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&mbr.getMyMatterSubsystemRep().getUDot(s)[uStart]);
}

const Vec7& MobilizedBody::Free::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 7);
    return Vec7::getAs(&qlike[qStart]);
}

const Vec6& MobilizedBody::Free::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::getAs(&ulike[uStart]);
}

Vec7& MobilizedBody::Free::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 7);
    return Vec7::updAs(&qlike[qStart]);
}

Vec6& MobilizedBody::Free::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 6);
    return Vec6::updAs(&ulike[uStart]);
}

    //////////////////////////////////////
    // MOBILIZED BODY::LINE ORIENTATION //
    //////////////////////////////////////

MobilizedBody::LineOrientation::LineOrientation() : PIMPLDerivedHandleBase(new LineOrientationImpl()) {
}


MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new LineOrientationImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::LineOrientation::LineOrientation(MobilizedBody& parent, const Transform& inbFrame,
                          const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new LineOrientationImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    ///////////////////////////////
    // MOBILIZED BODY::FREE LINE //
    ///////////////////////////////

MobilizedBody::FreeLine::FreeLine() : PIMPLDerivedHandleBase(new FreeLineImpl()) {
}


MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new FreeLineImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::FreeLine::FreeLine(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new FreeLineImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    //////////////////////////
    // MOBILIZED BODY::WELD //
    //////////////////////////

MobilizedBody::Weld::Weld() : PIMPLDerivedHandleBase(new WeldImpl()) {
}


MobilizedBody::Weld::Weld(MobilizedBody& parent, const Body& body) : PIMPLDerivedHandleBase(new WeldImpl()) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Weld::Weld(MobilizedBody& parent, const Transform& inbFrame,
                                  const Body& body, const Transform& outbFrame) : PIMPLDerivedHandleBase(new WeldImpl()) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

    ////////////////////////////////
    // (IM)MOBILIZED BODY::GROUND //
    ////////////////////////////////

MobilizedBody::Ground::Ground() : PIMPLDerivedHandleBase(new GroundImpl()) {
    setBody(Body::Ground());
}

    ///////////////////////////
    // MOBILIZED BODY::SCREW //
    ///////////////////////////

MobilizedBody::Screw::Screw(Real pitch) : PIMPLDerivedHandleBase(new ScrewImpl(pitch)) {
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Body& body, Real pitch) : PIMPLDerivedHandleBase(new ScrewImpl(pitch)) {
    // inb & outb frames are just the parent body's frame and new body's frame
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
                                                   *this);
}

MobilizedBody::Screw::Screw(MobilizedBody& parent, const Transform& inbFrame,
                            const Body& body, const Transform& outbFrame,
                            Real pitch) : PIMPLDerivedHandleBase(new ScrewImpl(pitch)) {
    setDefaultInboardFrame(inbFrame);
    setDefaultOutboardFrame(outbFrame);
    setBody(body);

    parent.updMatterSubsystem().adoptMobilizedBody(parent.getMobilizedBodyId(),
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
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQ(s)[qStart];
}
void MobilizedBody::Screw::setQ(State& s, Real q) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    mbr.getMyMatterSubsystemRep().updQ(s)[qStart] = q;
}
Real MobilizedBody::Screw::getQDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDot(s)[qStart];
}
Real MobilizedBody::Screw::getQDotDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int qStart, nq; mbr.findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return mbr.getMyMatterSubsystemRep().getQDotDot(s)[qStart];
}


Real MobilizedBody::Screw::getU(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getU(s)[uStart];
}
void MobilizedBody::Screw::setU(State& s, Real u) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    mbr.getMyMatterSubsystemRep().updU(s)[uStart] = u;
}
Real MobilizedBody::Screw::getUDot(const State& s) const {
    const MobilizedBodyImpl& mbr = MobilizedBody::getImpl();
    int uStart, nu; mbr.findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return mbr.getMyMatterSubsystemRep().getUDot(s)[uStart];
}

Real MobilizedBody::Screw::getMyPartQ(const State& s, const Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real MobilizedBody::Screw::getMyPartU(const State& s, const Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

Real& MobilizedBody::Screw::updMyPartQ(const State& s, Vector& qlike) const {
    int qStart, nq; getImpl().findMobilizerQs(s,qStart,nq); assert(nq == 1);
    return qlike[qStart];
}

Real& MobilizedBody::Screw::updMyPartU(const State& s, Vector& ulike) const {
    int uStart, nu; getImpl().findMobilizerUs(s,uStart,nu); assert(nu == 1);
    return ulike[uStart];
}

    ////////////////////////////
    // MOBILIZED BODY::CUSTOM //
    ////////////////////////////

MobilizedBody::Custom::Custom(int nMobilities, int nCoordinates) : PIMPLDerivedHandleBase(new CustomImpl(nMobilities, nCoordinates)) {
}

} // namespace SimTK

