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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintRep.h"
#include "SimbodyMatterSubsystemRep.h"
#include "ConstraintNode.h"

namespace SimTK {


    ////////////////
    // CONSTRAINT //
    ////////////////

bool Constraint::isEmptyHandle() const {return rep==0;}
bool Constraint::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}

void Constraint::disown(Constraint& newOwnerHandle) {
    SimTK_ASSERT_ALWAYS(rep && rep->myHandle==this,
        "disown() not allowed for an empty or non-owner Constraint handle.");
    SimTK_ASSERT_ALWAYS(!newOwnerHandle.rep,
        "disown() can only transfer ownership to an empty Constraint handle.");

    newOwnerHandle.setRep(*rep);
    rep->setMyHandle(newOwnerHandle);
}

Constraint::~Constraint() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

// Make this Constraint a non-owner handle referring to the same
// object as the source.
Constraint::Constraint(Constraint& src) : rep(src.rep) {
}

// Make this empty or non-owner handle refer to the same object
// as the source. This is illegal if the current handle is an
// owner.
Constraint& Constraint::operator=(Constraint& src) {
    if (&src != this) {
        SimTK_ASSERT_ALWAYS(!(rep && rep->myHandle==this),
            "You can't reassign the owner handle of a Constraint.");
        rep = src.rep;
    }
    return *this;
}

const SimbodyMatterSubsystem& Constraint::getMatterSubsystem() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return getRep().getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

ConstraintId Constraint::getConstraintId() const {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "getConstraintId() called on a Constraint that is not part of a subsystem.");
    return rep->myConstraintId;
}

SimbodyMatterSubsystem& Constraint::updMatterSubsystem() {
    SimTK_ASSERT_ALWAYS(isInSubsystem(),
        "updMatterSubsystem() called on a Constraint that is not part of a subsystem.");
    return updRep().updMyMatterSubsystemRep().updMySimbodyMatterSubsystemHandle();
}

bool Constraint::isInSubsystem() const {
    return rep && rep->isInSubsystem();
}

bool Constraint::isInSameSubsystem(const MobilizedBody& body) const {
    return rep && rep->isInSameSubsystem(body);
}

int Constraint::getNumConstrainedBodies() const {
    assert(getRep().subsystemTopologyHasBeenRealized());
    return (int)getRep().myConstrainedBodies.size();
}

const MobilizedBody& Constraint::getConstrainedBody(ConstrainedBodyId b) const {
    assert(getRep().subsystemTopologyHasBeenRealized());
    return getMatterSubsystem().getMobilizedBody(getRep().myConstrainedBodies[b]);
}

const SimbodyMatterSubsystem::Subtree& Constraint::getSubtree() const {
    assert(getRep().subsystemTopologyHasBeenRealized());
    return getRep().mySubtree;
}

    /////////////////////
    // CONSTRAINT::ROD //
    /////////////////////

Constraint::Rod::Rod(MobilizedBody& body1, MobilizedBody& body2, Real defaultRodLength)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    rep = new RodRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();
    updRep().defaultRodLength = defaultRodLength;

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);
}

Constraint::Rod::Rod(MobilizedBody& body1, const Vec3& point1,
                     MobilizedBody& body2, const Vec3& point2, Real defaultRodLength)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Rod(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Rod(): both bodies to be connected must be in the same MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(defaultRodLength > 0,
        "Constraint::Rod(): Rod length must always be greater than zero");

    rep = new RodRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();
    updRep().defaultPoint1 = point1;
    updRep().defaultPoint2 = point2;
    updRep().defaultRodLength = defaultRodLength;

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().B1 = updRep().addConstrainedBody(body1);
    updRep().B2 = updRep().addConstrainedBody(body2);
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody1(const Vec3& p1) {
    updRep().defaultPoint1 = p1;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultPointOnBody2(const Vec3& p2) {
    updRep().defaultPoint2 = p2;
    return *this;
}

Constraint::Rod& Constraint::Rod::setDefaultRodLength(Real length) {
    updRep().defaultRodLength = length;
    return *this;
}


MobilizedBodyId Constraint::Rod::getBody1Id() const {
    return getRep().body1;
}
MobilizedBodyId Constraint::Rod::getBody2Id() const {
    return getRep().body2;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody1() const {
    return getRep().defaultPoint1;
}
const Vec3& Constraint::Rod::getDefaultPointOnBody2() const {
    return getRep().defaultPoint2;
}
Real Constraint::Rod::getDefaultRodLength() const {
    return getRep().defaultRodLength;
}



    // Rod bookkeeping //

bool Constraint::Rod::isInstanceOf(const Constraint& s) {
    return RodRep::isA(s.getRep());
}
const Constraint::Rod& Constraint::Rod::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Rod&>(s);
}
Constraint::Rod& Constraint::Rod::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Rod&>(s);
}
const Constraint::Rod::RodRep& Constraint::Rod::getRep() const {
    return dynamic_cast<const RodRep&>(*rep);
}
Constraint::Rod::RodRep& Constraint::Rod::updRep() {
    return dynamic_cast<RodRep&>(*rep);
}


    ////////////////////////////////
    // CONSTRAINT::POINT IN PLANE //
    ////////////////////////////////

Constraint::PointInPlane::PointInPlane
   (MobilizedBody& planeBody,    const UnitVec3& defPlaneNormal, Real defPlaneHeight,
    MobilizedBody& followerBody, const Vec3&     defFollowerPoint)
{
    SimTK_ASSERT_ALWAYS(planeBody.isInSubsystem() && followerBody.isInSubsystem(),
        "Constraint::PointInPlane(): both bodies must already be in a SimbodyMatterSubsystem.");
    SimTK_ASSERT_ALWAYS(planeBody.isInSameSubsystem(followerBody),
        "Constraint::PointInPlane(): both bodies to be connected must be in the same SimbodyMatterSubsystem.");

    rep = new PointInPlaneRep(); rep->setMyHandle(*this);

    updRep().planeBody    = planeBody.getMobilizedBodyId();
    updRep().followerBody = followerBody.getMobilizedBodyId();
    updRep().defaultPlaneNormal   = defPlaneNormal;
    updRep().defaultPlaneHeight   = defPlaneHeight;
    updRep().defaultFollowerPoint = defFollowerPoint;

    planeBody.updMatterSubsystem().adoptConstraint(*this);

    updRep().addConstrainedBody(planeBody);
    updRep().addConstrainedBody(followerBody);
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneNormal(const UnitVec3& n) {
    getRep().invalidateTopologyCache();
    updRep().defaultPlaneNormal = n;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultPlaneHeight(Real h) {
    getRep().invalidateTopologyCache();
    updRep().defaultPlaneHeight = h;
    return *this;
}

Constraint::PointInPlane& Constraint::PointInPlane::setDefaultFollowerPoint(const Vec3& p) {
    getRep().invalidateTopologyCache();
    updRep().defaultFollowerPoint = p;
    return *this;
}

MobilizedBodyId Constraint::PointInPlane::getPlaneBodyId() const {
    return getRep().planeBody;
}
MobilizedBodyId Constraint::PointInPlane::getFollowerBodyId() const {
    return getRep().followerBody;
}
const UnitVec3& Constraint::PointInPlane::getDefaultPlaneNormal() const {
    return getRep().defaultPlaneNormal;
}
Real Constraint::PointInPlane::getDefaultPlaneHeight() const {
    return getRep().defaultPlaneHeight;
}
const Vec3& Constraint::PointInPlane::getDefaultFollowerPoint() const {
    return getRep().defaultFollowerPoint;
}

Constraint::PointInPlane& Constraint::PointInPlane::setPlaneDisplayHalfWidth(Real h) {
    updRep().setPlaneDisplayHalfWidth(h);
    return *this;
}
Constraint::PointInPlane& Constraint::PointInPlane::setPointDisplayRadius(Real r) {
    updRep().setPointDisplayRadius(r);
    return *this;
}

Real Constraint::PointInPlane::getPlaneDisplayHalfWidth() const {
    return getRep().getPlaneDisplayHalfWidth();
}

Real Constraint::PointInPlane::getPointDisplayRadius() const {
    return getRep().getPointDisplayRadius();
}

    // PointInPlane bookkeeping //

bool Constraint::PointInPlane::isInstanceOf(const Constraint& s) {
    return PointInPlaneRep::isA(s.getRep());
}
const Constraint::PointInPlane& Constraint::PointInPlane::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const PointInPlane&>(s);
}
Constraint::PointInPlane& Constraint::PointInPlane::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<PointInPlane&>(s);
}
const Constraint::PointInPlane::PointInPlaneRep& Constraint::PointInPlane::getRep() const {
    return dynamic_cast<const PointInPlaneRep&>(*rep);
}

Constraint::PointInPlane::PointInPlaneRep& Constraint::PointInPlane::updRep() {
    return dynamic_cast<PointInPlaneRep&>(*rep);
}

    // PointInPlaneRep

void Constraint::PointInPlane::PointInPlaneRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    // We can't generate the artwork until we know the normal, height, and follower
    // point location, which might not be until Instance stage.
    if (stage == Stage::Instance) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        // TODO: should be instance-stage data from State rather than topological data
        // This makes z axis point along plane normal
        const Transform X_B1(Rotation(defaultPlaneNormal,ZAxis), defaultPlaneHeight*defaultPlaneNormal);
        const Transform X_B2(Rotation(), defaultFollowerPoint);

        if (planeHalfWidth > 0 && pointRadius > 0) {
            // On the inboard body, draw a solid sphere and a wireframe one attached to it for
            // easier visualization of its rotation. These are at about 90% of the radius.
            geom.push_back(DecorativeBrick(Vec3(planeHalfWidth,planeHalfWidth,pointRadius/2))
                                                .setColor(Gray)
                                                .setRepresentation(DecorativeGeometry::DrawSurface)
                                                .setOpacity(0.3)
                                                .setBodyId(planeBody)
                                                .setTransform(X_B1));
            geom.push_back(DecorativeBrick(Vec3(planeHalfWidth,planeHalfWidth,pointRadius/2))
                                                .setColor(Black)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setBodyId(planeBody)
                                                .setTransform(X_B1));

            // On the follower body draw an orange mesh sphere at the ball radius.
            geom.push_back(DecorativeSphere(pointRadius)
                                                .setColor(Orange)
                                                .setRepresentation(DecorativeGeometry::DrawWireframe)
                                                .setResolution(0.5)
                                                .setBodyId(followerBody)
                                                .setTransform(X_B2));
        }
    }
}



    //////////////////////
    // CONSTRAINT::BALL //
    //////////////////////

Constraint::Ball::Ball(MobilizedBody& body1, MobilizedBody& body2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new BallRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().addConstrainedBody(body1);
    updRep().addConstrainedBody(body2);
}

Constraint::Ball::Ball(MobilizedBody& body1, const Vec3& point1,
                       MobilizedBody& body2, const Vec3& point2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Ball(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Ball(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new BallRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();
    updRep().defaultPoint1 = point1;
    updRep().defaultPoint2 = point2;

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().addConstrainedBody(body1);
    updRep().addConstrainedBody(body2);
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody1(const Vec3& p1) {
    getRep().invalidateTopologyCache();
    updRep().defaultPoint1 = p1;
    return *this;
}

Constraint::Ball& Constraint::Ball::setDefaultPointOnBody2(const Vec3& p2) {
    getRep().invalidateTopologyCache();
    updRep().defaultPoint2 = p2;
    return *this;
}

MobilizedBodyId Constraint::Ball::getBody1Id() const {
    return getRep().body1;
}
MobilizedBodyId Constraint::Ball::getBody2Id() const {
    return getRep().body2;
}
const Vec3& Constraint::Ball::getDefaultPointOnBody1() const {
    return getRep().defaultPoint1;
}
const Vec3& Constraint::Ball::getDefaultPointOnBody2() const {
    return getRep().defaultPoint2;
}

Constraint::Ball& Constraint::Ball::setDefaultRadius(Real r) {
    getRep().invalidateTopologyCache();
    updRep().setDefaultRadius(r);
    return *this;
}

Real Constraint::Ball::getDefaultRadius() const {
    return getRep().getDefaultRadius();
}


    // Ball bookkeeping //

bool Constraint::Ball::isInstanceOf(const Constraint& s) {
    return BallRep::isA(s.getRep());
}
const Constraint::Ball& Constraint::Ball::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Ball&>(s);
}
Constraint::Ball& Constraint::Ball::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Ball&>(s);
}
const Constraint::Ball::BallRep& Constraint::Ball::getRep() const {
    return dynamic_cast<const BallRep&>(*rep);
}

Constraint::Ball::BallRep& Constraint::Ball::updRep() {
    return dynamic_cast<BallRep&>(*rep);
}

    // BallRep

void Constraint::Ball::BallRep::calcDecorativeGeometryAndAppendImpl
   (const State& s, Stage stage, Array<DecorativeGeometry>& geom) const
{
    // We can't generate the ball until we know the radius, and we can't place
    // the geometry on the body until we know the body1 and body2 point
    // placements on the bodies, which might not be until Instance stage.
    if (stage == Stage::Instance && getDefaultRadius() > 0) {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const Transform X_B1(Rotation(), defaultPoint1); // should be point from State
        const Transform X_B2(Rotation(), defaultPoint2);

        // On the inboard body, draw a solid sphere and a wireframe one attached to it for
        // easier visualization of its rotation. These are at about 90% of the radius.
        geom.push_back(DecorativeSphere(0.92*getDefaultRadius())
                                            .setColor(Gray)
                                            .setRepresentation(DecorativeGeometry::DrawSurface)
                                            .setOpacity(0.5)
                                            .setResolution(0.75)
                                            .setBodyId(body1)
                                            .setTransform(X_B1));
        geom.push_back(DecorativeSphere(0.90*getDefaultRadius())
            .setColor(White)
            .setRepresentation(DecorativeGeometry::DrawWireframe)
            .setResolution(0.75)
            .setLineThickness(3)
            .setOpacity(0.1)
            .setBodyId(body1)
            .setTransform(X_B1));

        // On the outboard body draw an orange mesh sphere at the ball radius.
        geom.push_back(DecorativeSphere(getDefaultRadius())
                                            .setColor(Orange)
                                            .setRepresentation(DecorativeGeometry::DrawWireframe)
                                            .setOpacity(0.5)
                                            .setResolution(0.5)
                                            .setBodyId(body2)
                                            .setTransform(X_B2));
    }
}

    //////////////////////
    // CONSTRAINT::WELD //
    //////////////////////

Constraint::Weld::Weld(MobilizedBody& body1, MobilizedBody& body2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new WeldRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().addConstrainedBody(body1);
    updRep().addConstrainedBody(body2);
}

Constraint::Weld::Weld(MobilizedBody& body1, const Transform& frame1,
                       MobilizedBody& body2, const Transform& frame2)
{
    SimTK_ASSERT_ALWAYS(body1.isInSubsystem() && body2.isInSubsystem(),
        "Constraint::Weld(): both bodies must already be in a MatterSubsystem.");
    SimTK_ASSERT_ALWAYS(body1.isInSameSubsystem(body2),
        "Constraint::Weld(): both bodies to be connected must be in the same MatterSubsystem.");

    rep = new WeldRep(); rep->setMyHandle(*this);

    updRep().body1 = body1.getMobilizedBodyId();
    updRep().body2 = body2.getMobilizedBodyId();
    updRep().defaultFrame1 = frame1;
    updRep().defaultFrame2 = frame2;

    body1.updMatterSubsystem().adoptConstraint(*this);

    updRep().addConstrainedBody(body1);
    updRep().addConstrainedBody(body2);
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody1(const Transform& f1) {
    updRep().defaultFrame1 = f1;
    return *this;
}

Constraint::Weld& Constraint::Weld::setDefaultFrameOnBody2(const Transform& f2) {
    updRep().defaultFrame2 = f2;
    return *this;
}

MobilizedBodyId Constraint::Weld::getBody1Id() const {
    return getRep().body1;
}
MobilizedBodyId Constraint::Weld::getBody2Id() const {
    return getRep().body2;
}
const Transform& Constraint::Weld::getDefaultFrameOnBody1() const {
    return getRep().defaultFrame1;
}
const Transform& Constraint::Weld::getDefaultFrameOnBody2() const {
    return getRep().defaultFrame2;
}


    // Weld bookkeeping //

bool Constraint::Weld::isInstanceOf(const Constraint& s) {
    return WeldRep::isA(s.getRep());
}
const Constraint::Weld& Constraint::Weld::downcast(const Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const Weld&>(s);
}
Constraint::Weld& Constraint::Weld::updDowncast(Constraint& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<Weld&>(s);
}
const Constraint::Weld::WeldRep& Constraint::Weld::getRep() const {
    return dynamic_cast<const WeldRep&>(*rep);
}
Constraint::Weld::WeldRep& Constraint::Weld::updRep() {
    return dynamic_cast<WeldRep&>(*rep);
}

    ////////////////////
    // CONSTRAINT REP //
    ////////////////////

/*virtual*/ Constraint::ConstraintRep::~ConstraintRep() {
    delete myConstraintNode; myConstraintNode=0;
}

void Constraint::ConstraintRep::realizeTopology(int& nxtQErr, int& nxtUErr, int& nxtMult) const
{
    // Calculate the relevant Subtree.
    mySubtree.clear();
    mySubtree.setSimbodyMatterSubsystem(getMyMatterSubsystem());
    for (ConstrainedBodyId b(0); b < (int)myConstrainedBodies.size(); ++b)
        mySubtree.addTerminalBody(myConstrainedBodies[b]);
    mySubtree.realizeTopology();

    // Create a constraint node for dealing with computational issues.
    delete myConstraintNode;
    myConstraintNode = createConstraintNode();
    myConstraintNode->setConstraintNum(myConstraintId);
    myConstraintNode->setQErrIndex(nxtQErr);
    myConstraintNode->setUErrIndex(nxtUErr);
    myConstraintNode->setMultIndex(nxtMult);
    const int nConsEqns = myConstraintNode->getNConstraintEquations();
    nxtQErr += nConsEqns;
    nxtUErr += nConsEqns;
    nxtMult += nConsEqns;
    myConstraintNode->finishConstruction(*myMatterSubsystemRep);
}

void Constraint::ConstraintRep::invalidateTopologyCache() const {
    delete myConstraintNode; myConstraintNode=0;
    if (myMatterSubsystemRep)
        myMatterSubsystemRep->invalidateSubsystemTopologyCache();
}

bool Constraint::ConstraintRep::subsystemTopologyHasBeenRealized() const {
    return myMatterSubsystemRep && myMatterSubsystemRep->subsystemTopologyHasBeenRealized();
}

const ConstraintNode& Constraint::ConstraintRep::getMyConstraintNode() const {
    SimTK_ASSERT(myConstraintNode && myMatterSubsystemRep && 
                 myMatterSubsystemRep->subsystemTopologyHasBeenRealized(),
      "An operation on a Constraint was illegal because realizeTopology() has "
      "not been performed on the containing Subsystem since the last topological change."
    );
    return *myConstraintNode;
}

void Constraint::ConstraintRep::setMyMatterSubsystem
   (SimbodyMatterSubsystem& matter, ConstraintId id)
{
    assert(!isInSubsystem());
    myMatterSubsystemRep = &matter.updRep();
    myConstraintId = id;
}

const SimbodyMatterSubsystem& 
Constraint::ConstraintRep::getMyMatterSubsystem() const {
    return getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
}

// These are measured from and expressed in the ancestor (A) frame.
//TODO: should precalculate in State, return reference
Transform Constraint::ConstraintRep::getBodyTransform(const State& s, ConstrainedBodyId B) const { // X_AB
    const Transform& X_GB = getMyMatterSubsystemRep().getBodyTransform(s, myConstrainedBodies[B]);
    const Transform& X_GA = getMyMatterSubsystemRep().getBodyTransform(s, mySubtree.getAncestorBody());
    return ~X_GA*X_GB;
}

SpatialVec Constraint::ConstraintRep::getBodyVelocity(const State& s, ConstrainedBodyId B) const { // V_AB
    const Transform&  X_GB = getMyMatterSubsystemRep().getBodyTransform(s, myConstrainedBodies[B]);
    const Transform&  X_GA = getMyMatterSubsystemRep().getBodyTransform(s, mySubtree.getAncestorBody());
    const SpatialVec& V_GB = getMyMatterSubsystemRep().getBodyVelocity(s, myConstrainedBodies[B]);
    const SpatialVec& V_GA = getMyMatterSubsystemRep().getBodyVelocity(s, mySubtree.getAncestorBody());
    const Vec3 p_AB_G     = X_GB.T() - X_GA.T();
    const Vec3 p_AB_G_dot = V_GB[1]  - V_GA[1];        // time deriv of p taken in G

    const Vec3 w_AB_G = V_GB[0] - V_GA[0];             // relative angular velocity
    const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // time deriv of p in A, exp in G
    return ~X_GA.R() * SpatialVec(w_AB_G, v_AB_G);     // re-express in A
}

SpatialVec Constraint::ConstraintRep::getBodyAcceleration(const State& s, ConstrainedBodyId B) const { // A_AB
    const Transform&  X_GB = getMyMatterSubsystemRep().getBodyTransform(s, myConstrainedBodies[B]);
    const Transform&  X_GA = getMyMatterSubsystemRep().getBodyTransform(s, mySubtree.getAncestorBody());
    const SpatialVec& V_GB = getMyMatterSubsystemRep().getBodyVelocity(s, myConstrainedBodies[B]);
    const SpatialVec& V_GA = getMyMatterSubsystemRep().getBodyVelocity(s, mySubtree.getAncestorBody());
    const SpatialVec& A_GB = getMyMatterSubsystemRep().getBodyAcceleration(s, myConstrainedBodies[B]);
    const SpatialVec& A_GA = getMyMatterSubsystemRep().getBodyAcceleration(s, mySubtree.getAncestorBody());
    const Vec3 p_AB_G        = X_GB.T() - X_GA.T();
    const Vec3 p_AB_G_dot    = V_GB[1]  - V_GA[1];     // taken in G
    const Vec3 p_AB_G_dotdot = A_GB[1]  - A_GA[1];     // taken in G

    const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G; // taken in A, exp. in G
    const Vec3 b_AB_G = A_GB[0] - A_GA[0];             // relative angular acceleration
    const Vec3 a_AB_G = p_AB_G_dotdot - (A_GA[0] % p_AB_G + V_GA[0] % p_AB_G_dot); // taken in A, exp. in G
    return ~X_GA.R() * SpatialVec(b_AB_G, a_AB_G);     // re-express in A
}

// Default implementations for ConstraintRep virtuals throw "unimplemented"
// exceptions. These shouldn't be called unless the concrete constraint has
// given a non-zero value for mp, mv, and/or ma which is a promise to 
// implement the associated methods.

    // These must be defined if there are any positin (holonomic) constraints defined.

void Constraint::ConstraintRep::
calcPositionErrorsVirtual(const State&, int mp,  Real* perr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcPositionErrors");
}

void Constraint::ConstraintRep::
calcPositionDotErrorsVirtual(const State&, int mp,  Real* pverr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcPositionDotErrors");
}

void Constraint::ConstraintRep::
calcPositionDotDotErrorsVirtual(const State&, int mp,  Real* paerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcPositionDotDotErrors");
}


void Constraint::ConstraintRep::
applyPositionConstraintForcesVirtual
   (const State&, int mp, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyPositionConstraintForces");
}

    // These must be defined if there are any velocity (nonholonomic) constraints defined.

void Constraint::ConstraintRep::
calcVelocityErrorsVirtual(const State&, int mv,  Real* verr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcVelocityErrors");
}


void Constraint::ConstraintRep::
calcVelocityDotErrorsVirtual(const State&, int mv,  Real* vaerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcVelocityDotErrors");
}


void Constraint::ConstraintRep::
applyVelocityConstraintForcesVirtual
   (const State&, int mv, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyVelocityConstraintForces");
}



// These must be defined if there are any acceleration-only constraints defined.
void Constraint::ConstraintRep::
calcAccelerationErrorsVirtual(const State&, int ma,  Real* aerr) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "calcAccelerationErrors");
}

void Constraint::ConstraintRep::
applyAccelerationConstraintForcesVirtual
   (const State&, int ma, const Real* multipliers,
    Vector_<SpatialVec>& bodyForces,
    Vector&              mobilityForces) const
{
    SimTK_THROW2(Exception::UnimplementedVirtualMethod,
        "Constraint::ConstraintRep", "applyAccelerationConstraintForces");
}


} // namespace SimTK

