/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * Implementation of the IVM Simbody interface.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/IVMSimbodyInterface.h"

#include "IVMSimbodyInterfaceRep.h"

#include "IVMRigidBodyTree.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>

static CDSVec3 toCDSVec3(const Vec3& v) {return CDSVec3(v[0],v[1],v[2]);}
static Vec3    toVec3(const CDSVec3& v) {return Vec3(v[0],v[1],v[2]);}
static CDSVec6 toCDSVec6(const SpatialVec& v) { 
    CDSVec6 cv;
    for (int i=0; i<3; ++i) cv[i] = v[0][i], cv[i+3] = v[1][i];
    return cv;
}
static SpatialVec toSpatialVector(const CDSVec6& v) 
  { return SpatialVec(Vec3(v[0],v[1],v[2]),Vec3(v[3],v[4],v[5])); }

static CDSVecVec6 toCDSVecVec6(const Vector_<SpatialVec>& a) {
    CDSVecVec6 vv(a.size());
    for (int i=0; i < (int)a.size(); ++i)
        vv(i) = toCDSVec6(a[i]);
    return vv;
}

static CDSMat33 toCDSMat33(const Mat33& m) {
    return CDSMat33(m(0,0), m(0,1), m(0,2),
                    m(1,0), m(1,1), m(1,2),
                    m(2,0), m(2,1), m(2,2));
}
static RVec toRVec(const Vector& v) {
    RVec r(v.size());   // a 1-based vector
    for (int i=0; i < v.size(); ++i)
        r[i+1] = v[i];
    return r;
}
static Vector toVector(const RVec& r) {
    return Vector(r.size(), r.pointer());
}
static Mat33 toMat33(const CDSMat33& m) {
    return Mat33(Row3(m(0,0), m(0,1), m(0,2)),
                 Row3(m(1,0), m(1,1), m(1,2)),
                 Row3(m(2,0), m(2,1), m(2,2)));
}
static IVMInertia toIVMInertia(const InertiaMat& i) {
    return IVMInertia(toCDSMat33(i.toMat33()));
}

static RotationMat toMatRotation(const CDSMat33& m) {
    const Mat33 m33 = toMat33(m);
    return reinterpret_cast<const RotationMat&>(m33);
}

static IVMMassProperties toIVMMassProperties(const Real& m, const Vec3& c, const InertiaMat& i) {
    return IVMMassProperties(m, toCDSVec3(c), toIVMInertia(i));
}

static IVMFrame toIVMFrame(const TransformMat& f) {
    return IVMFrame(toCDSMat33(f.getRotation().asMat33()), toCDSVec3(f.getTranslation()));
}

static TransformMat toFrame(const IVMFrame& f) {
    return TransformMat(toMatRotation(f.getRot_RF()), toVec3(f.getLoc_RF()));
}
    // IVM SIMBODY INTERFACE //

IVMSimbodyInterface::IVMSimbodyInterface(const Multibody& m, bool oldStyle) : rep(0) {
    rep = oldStyle ? (IVMSimbodyInterfaceRep*)new OldIVMSimbodyInterfaceRep(m)
                   : (IVMSimbodyInterfaceRep*)new NewIVMSimbodyInterfaceRep(m);
    rep->setMyHandle(*this);
}

int IVMSimbodyInterface::getNBodies()     const {return rep->getNBodies();}
int IVMSimbodyInterface::getNParameters() const {return rep->getNBodies();}
int IVMSimbodyInterface::getNQ()          const {return rep->getNQ();}
int IVMSimbodyInterface::getNU()          const {return rep->getNU();}

void IVMSimbodyInterface::realizeParameters(const State& s) const {
    rep->realizeParameters(s);
}

void IVMSimbodyInterface::realizeConfiguration(const State& s) const {
    rep->realizeConfiguration(s);
}

void IVMSimbodyInterface::realizeMotion(const State& s) const {
    rep->realizeMotion(s);
}

void IVMSimbodyInterface::realizeReaction(const State& s) const {
    rep->realizeReaction(s);
}

void IVMSimbodyInterface::enforceConfigurationConstraints(State&) const {
    assert(false); //TODO
}
void IVMSimbodyInterface::enforceMotionConstraints(State&) const {
    assert(false); //TODO
}

const Vector& 
IVMSimbodyInterface::getQDot(const State& s) const {
    return rep->getQDot(s);
}

Vector 
IVMSimbodyInterface::calcUDot(const State& s, 
                              const Vector_<SpatialVec>& bodyForces,
                              const Vector& hingeForces) const
{
    return rep->calcUDot(s,bodyForces,hingeForces);
}



// TODO: this should return a reference to a cache entry but that is awkward at the moment
// since we and IVM don't agree on the body frame. IVM uses the computed reference frame;
// Simbody uses a frame chosen by the user. Note that SD/FAST works like IVM in this regard.
TransformMat
IVMSimbodyInterface::getBodyConfiguration(const State& s, const Body& b) const
{
    return rep->getBodyConfiguration(s,b);
}

SpatialVec // TODO: this should return a reference to a cache entry
IVMSimbodyInterface::getBodyVelocity(const State& s, const Body& b) const
{
    return rep->getBodyVelocity(s,b);
}

SpatialVec // TODO: this should return a reference to a cache entry
IVMSimbodyInterface::getBodyAcceleration(const State& s, const Body& b) const
{
    return rep->getBodyAcceleration(s,b);
}

State IVMSimbodyInterface::getDefaultState() const {
    return rep->getDefaultState();
}

// Accumulate the forces at the IVM reference frame origin (the inboard joint location).
// The passed in point is measured from the Simbody origin OB and expressed in B. We 
// want it measured from the IVM origin OR and expressed in G. And we want the force in G.
void IVMSimbodyInterface::applyPointForce(const State& s, const Body& b,
                                          const Vec3& pt_B, const Vec3& frc_B, 
                                          Vector_<SpatialVec>& bodyForces) const
{
    const RBTreeMap& info  = rep->getBodyInfo(b);
    const int        index = info.getRBIndex();
    const TransformMat&     F_BR  = info.getRefFrameInBody();
    const TransformMat      F_GB  = getBodyConfiguration(s,b);
    const TransformMat      F_GR  = F_GB.compose(F_BR);

    const Vec3 pt_R  = F_BR.shiftBaseStationToFrame(pt_B); // vector from OR to pt, expressed in R
    const Vec3 pt_RG = F_GR.xformFrameVecToBase(pt_R);     // re-express in G (but still measured from OR)
    const Vec3 frc_G = F_GB.xformFrameVecToBase(frc_B);    // re-express in G

    bodyForces[index][0] += pt_RG % frc_G; // shift to OR: introduces moment v X f
    bodyForces[index][1] += frc_G;         // frc unchanged
}

void IVMSimbodyInterface::applyBodyTorque(const State& s, const Body& b, const Vec3& trq_B, 
                                          Vector_<SpatialVec>& bodyForces) const
{
    const RBTreeMap& info    = rep->getBodyInfo(b);
    const int        rbIndex = info.getRBIndex();

    const TransformMat F_GB  = getBodyConfiguration(s,b);
    const Vec3  trq_G = F_GB.xformFrameVecToBase(trq_B);    // re-express in G
    bodyForces[rbIndex][0] += trq_G;
}

// Given g in the ground frame, apply a force mg to the center of mass of each body.
void IVMSimbodyInterface::applyGravity(const State& s, const Vec3& g, 
                                       Vector_<SpatialVec>& bodyForces) const
{
    // skip ground
    for (int i=1; i < getNBodies(); ++i) {
        const RBTreeMap& info    = rep->getBodyInfoByIndex(i);
        const Body&      body    = info.getBody();
        const TransformMat&     F_BR    = info.getRefFrameInBody();
        const Vec3&      com_R   = info.getCOMInRef();
        const Real&      mass    = info.getMass();

        const TransformMat      F_GB    = getBodyConfiguration(s,body);
        const TransformMat      F_GR    = F_GB.compose(F_BR);
        const int        rbIndex = info.getRBIndex();

        const Vec3 com_RG = F_GR.xformFrameVecToBase(com_R); // OR to COM, in G
        const Vec3 frc_G  = mass*g;

        bodyForces[i][0] += com_RG % frc_G; // shift to OR: introduces moment v X f
        bodyForces[i][1] += frc_G;          // frc unchanged
    }
}

void IVMSimbodyInterface::applyHingeForce(const State&, const Joint& j, int axis, const Real& frc, 
                     Vector& hingeForces) const
{
    assert(false); //TODO
}

    // IVM SIMBODY INTERFACE REP //

IVMSimbodyInterfaceRep::IVMSimbodyInterfaceRep(const Multibody& m) 
  : handle(0), mbs(m)
{
    // Make sure the subsystem has been realized so we can get numerical
    // values out of it.
    mbs.realize(Stage::Startup);

    // First find a tree within the multibody system.

    std::vector<const Joint*> joints;
    for (int i=0; i<mbs.getNSubsystems(); ++i)
        if (Joint::isInstanceOf(mbs[i]))
            joints.push_back(&Joint::downcast(mbs[i]));


    size_t nxt=0;
    mbs2tree.push_back(RBTreeMap(&Body::downcast(mbs["Ground"]),TransformMat(),TransformMat(),
                       NTraits<Real>::getInfinity(), Vec3(0), InertiaMat(), // mass, com, inertia
                       0,0,0)); // no parent, no inboard joint, level 0

    while (nxt < mbs2tree.size()) {
        const Body& parentBody = mbs2tree[nxt].getBody();

        // Find all the joints whose reference bodies are this parent.
        for (size_t i=0; i<joints.size(); ++i) {
            if (!Body::getPlacementBody(joints[i]->getReferenceFrame())
                                        .isSameSubsystem(parentBody))
                continue;

            const Body& childBody = Body::getPlacementBody(joints[i]->getMovingFrame());

            // Calculate the reference frame, which is located at origin of J
            // and aligned with its parent's reference frame.
            const TransformMat& fBJ  = joints[i]->getMovingFrame().getValue();       // on B
            const TransformMat& fPJi = joints[i]->getReferenceFrame().getValue();    // on P
            const TransformMat fBR(fBJ.getRotation()*~fPJi.getRotation(), 
                                   fBJ.getTranslation());
            const TransformMat fRJ(fPJi.getRotation(), Vec3(0));

            const Real&       mass      = childBody.getMass().getValue();
            const Vec3&       com_B     = childBody.getMassCenter().getValue();
            const InertiaMat& iner_OB_B = childBody.getInertia().getValue();

            const Vec3 com_R = fBR.shiftBaseStationToFrame(com_B);
            const InertiaMat iner_CB_B = iner_OB_B.shiftToCOM(com_B,mass);
            const InertiaMat iner_CB_R = iner_CB_B.changeAxes(fBR.getRotation());
            const InertiaMat iner_OR_R = iner_CB_R.shiftFromCOM(-com_R,mass);

            mbs2tree.push_back(RBTreeMap(&childBody,
                                         fBR, fRJ, mass, com_R, iner_OR_R,
                                         nxt, // index of parent
                                         joints[i],
                                         mbs2tree[nxt].getLevel() + 1));
        }
        ++nxt;
    }
}
   
/*static*/ IVMJointType
IVMSimbodyInterfaceRep::mapToIVMJointType(Joint::JointType jt) {
    switch (jt) {
    case Joint::UnknownJointType:   return IVMUnknownJointType;
    case Joint::ThisIsGround:       return IVMThisIsGround;
    case Joint::Weld:               return IVMWeldJoint;
    case Joint::Torsion:            return IVMTorsionJoint;  // aka PinJoint
    case Joint::Sliding:            return IVMSlidingJoint;
    case Joint::Universal:          return IVMUJoint;
    case Joint::Cylinder:           return IVMCylinderJoint;
    case Joint::Planar:             return IVMPlanarJoint;
    case Joint::Gimbal:             return IVMGimbalJoint;
    case Joint::Orientation:        return IVMOrientationJoint; // aka BallJoint
    case Joint::Cartesian:          return IVMCartesianJoint;
    case Joint::FreeLine:           return IVMFreeLineJoint;
    case Joint::Free:               return IVMFreeJoint;
    default: assert(false);
    }
    //NOTREACHED
    return IVMUnknownJointType;
}

    // OLD IVM SIMBODY INTERFACE REP

State OldIVMSimbodyInterfaceRep::getDefaultState() const {
    const IVMRigidBodyTree& t = getRigidBodyTree();
    State s(getNQ(),getNU());
    RVec pos(getNQ()), vel(getNQ());
    t.getPos(pos); t.getVel(vel);
    s.updQ() = toVector(pos);
    s.updU() = toVector(vel);
    return s;
}

void OldIVMSimbodyInterfaceRep::realizeConfiguration(const State& s) const {
    IVMRigidBodyTree& t = const_cast<IVMRigidBodyTree&>(tree);
    t.setPos(toRVec(s.getQ()));
}
void OldIVMSimbodyInterfaceRep::realizeMotion(const State& s) const {
    IVMRigidBodyTree& t = const_cast<IVMRigidBodyTree&>(tree);
    t.setVel(toRVec(s.getU()));
}

Vector
OldIVMSimbodyInterfaceRep::calcUDot(const State& s, 
                        const Vector_<SpatialVec>& bodyForces,
                        const Vector& hingeForces) const
{
    IVMRigidBodyTree& t = const_cast<IVMRigidBodyTree&>(tree);
    t.prepareForDynamics();
    t.calcLoopForwardDynamics(toCDSVecVec6(bodyForces));
    RVec a(s.getU().size());
    t.getAcc(a);
    return toVector(a);
}

TransformMat
OldIVMSimbodyInterfaceRep::getBodyConfiguration(const State& s, const Body& body) const {
    // Get from IVM the body reference frame R in ground.
    const RBTreeMap& info = getBodyInfo(body);
    const TransformMat& F_BR = info.getRefFrameInBody();
    const IVMRigidBodyNode& n = getRigidBodyTree().getRigidBodyNode(info.getRBIndex());

    const TransformMat F_GR(toMatRotation(n.getR_GB()), toVec3(n.getOB_G()));
    const TransformMat F_RB(F_BR.invert());
    return F_GR.compose(F_RB);
}

void OldIVMSimbodyInterfaceRep::buildTree() {
    cout << "**** OLD IVM TREE ****" << endl;
    int nextStateOffset = 1; // Because RVecs are 1-based
    for (size_t i=0; i<mbs2tree.size(); ++i) {
        RBTreeMap& childEntry = mbs2tree[i];

        // Deal with ground.
        cout << childEntry.getLevel() << ": " << childEntry.getBody().getFullName() << endl;
        if (!childEntry.getLevel()) {
            IVMRigidBodyNode* rb = IVMRigidBodyNode::create(IVMMassProperties(), IVMFrame(), IVMThisIsGround,
                                      false, false, nextStateOffset);
            const int rbIndex = tree.addGroundNode(rb);
            childEntry.setRBIndex(rbIndex);
            continue;
        }

        // Not ground -- we must have a parent.
        const RBTreeMap& parentEntry = mbs2tree[childEntry.getParentIndex()];

        cout << Joint::getJointTypeName(childEntry.getJoint().getJointType())
             << " Joint: "  << childEntry.getJoint().getFullName() 
             << " Parent: " << parentEntry.getBody().getFullName() << endl;

        cout << " mass=" << childEntry.getMass() 
             << " com_R=" << childEntry.getCOMInRef()
             << " iner_OR_R=" << childEntry.getInertiaAboutRef() 
             << endl;

        cout << "frame_BJ=" << childEntry.getJoint().getMovingFrame().getValue() << endl;
        cout << "frame_BR=" << childEntry.getRefFrameInBody()   << endl;
        cout << "frame_RJ=" << childEntry.getJointFrameInRef()  << endl;

        cout << "JointType=" << mbs2tree[i].getJoint().getJointType()
             << "  IVMJointType=" << mapToIVMJointType(mbs2tree[i].getJoint().getJointType())
             << endl;

        const int save = nextStateOffset;
        IVMRigidBodyNode* rb = IVMRigidBodyNode::create(
            toIVMMassProperties(childEntry.getMass(), 
                                childEntry.getCOMInRef(), 
                                childEntry.getInertiaAboutRef()),
            toIVMFrame(childEntry.getJointFrameInRef()),
            mapToIVMJointType(childEntry.getJoint().getJointType()),
            false,
            false,
            nextStateOffset);
        cout << "CREATED: states " << save << "-" << nextStateOffset-1 << endl;


        IVMRigidBodyNode& parent = tree.updRigidBodyNode(parentEntry.getRBIndex());
        const int rbIndex = tree.addRigidBodyNode(parent,IVMFrame()/*XXX*/,rb);
        childEntry.setRBIndex(rbIndex);
    }

    tree.finishConstruction(1e-6, 0);

    std::cout << "*** IVMRigidBodyTree:" << std::endl << tree << std::endl;
}

    // NEW IVM SIMBODY INTERFACE REP

State NewIVMSimbodyInterfaceRep::getDefaultState() const {
    const RigidBodyTree& t = getRigidBodyTree();
    State s(getNQ(),getNU());
    t.getPos(s.updQ()); t.getVel(s.updU());
    return s;
}


Vector
NewIVMSimbodyInterfaceRep::calcUDot(const State& s, 
                        const Vector_<SpatialVec>& bodyForces,
                        const Vector& hingeForces) const
{
    RigidBodyTree& t = const_cast<RigidBodyTree&>(tree);
    t.prepareForDynamics();
    t.calcLoopForwardDynamics(bodyForces);
    Vector a(s.getU().size());
    tree.getAcc(a);
    return a;
}

TransformMat
NewIVMSimbodyInterfaceRep::getBodyConfiguration(const State& s, const Body& body) const {
    // Get from IVM the body reference frame R in ground.
    const RBTreeMap& info = getBodyInfo(body);
    const TransformMat& F_BR = info.getRefFrameInBody();
    const RigidBodyNode& n = getRigidBodyTree().getRigidBodyNode(info.getRBIndex());

    const TransformMat F_GR(n.getR_GB(), n.getOB_G());
    const TransformMat F_RB(F_BR.invert());
    return F_GR.compose(F_RB);
}


void NewIVMSimbodyInterfaceRep::buildTree() {
    cout << "**** NEW RB TREE ****" << endl;
    int nextStateOffset = 0; // Because Vectors are 0-based
    for (size_t i=0; i<mbs2tree.size(); ++i) {
        RBTreeMap& childEntry = mbs2tree[i];

        // Deal with ground.
        cout << childEntry.getLevel() << ": " << childEntry.getBody().getFullName() << endl;
        if (!childEntry.getLevel()) {
            RigidBodyNode* rb = RigidBodyNode::create(MassProperties(), TransformMat(), Joint::ThisIsGround,
                                      false, false, nextStateOffset);
            const int rbIndex = tree.addGroundNode(rb);
            childEntry.setRBIndex(rbIndex);
            continue;
        }

        // Not ground -- we must have a parent.
        const RBTreeMap& parentEntry = mbs2tree[childEntry.getParentIndex()];

        cout << Joint::getJointTypeName(childEntry.getJoint().getJointType())
             << " Joint: "  << childEntry.getJoint().getFullName() 
             << " Parent: " << parentEntry.getBody().getFullName() << endl;

        cout << " mass=" << childEntry.getMass() 
             << " com_R=" << childEntry.getCOMInRef()
             << " iner_OR_R=" << childEntry.getInertiaAboutRef() 
             << endl;

        cout << "frame_BJ=" << childEntry.getJoint().getMovingFrame().getValue() << endl;
        cout << "frame_BR=" << childEntry.getRefFrameInBody()   << endl;
        cout << "frame_RJ=" << childEntry.getJointFrameInRef()  << endl;

        cout << "JointType=" << mbs2tree[i].getJoint().getJointType()
             << endl;

        const int save = nextStateOffset;
        RigidBodyNode* rb = RigidBodyNode::create(
            MassProperties(childEntry.getMass(), 
                           childEntry.getCOMInRef(), 
                           childEntry.getInertiaAboutRef()),
            childEntry.getJointFrameInRef(),
            childEntry.getJoint().getJointType(),
            false,
            false,
            nextStateOffset);
        cout << "CREATED: states " << save << "-" << nextStateOffset-1 << endl;

        RigidBodyNode& parent = tree.updRigidBodyNode(parentEntry.getRBIndex());
        const int rbIndex = tree.addRigidBodyNode(parent,TransformMat()/*XXX*/,rb);
        childEntry.setRBIndex(rbIndex);
    }

    tree.finishConstruction(1e-6, 0);
    std::cout << "*** RigidBodyTree:" << std::endl << tree << std::endl;
}
