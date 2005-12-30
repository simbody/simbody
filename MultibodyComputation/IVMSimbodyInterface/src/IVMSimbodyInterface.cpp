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

#include "RigidBodyTree.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>

static CDSVec3 toCDSVec3(const Vec3& v) {return CDSVec3(v[0],v[1],v[2]);}
static Vec3    toVec3(const CDSVec3& v) {return Vec3(v[0],v[1],v[2]);}

static CDSMat33 toCDSMat33(const Mat33& m) {
    return CDSMat33(m(0,0), m(0,1), m(0,2),
                    m(1,0), m(1,1), m(1,2),
                    m(2,0), m(2,1), m(2,2));
}
static Mat33 toMat33(const CDSMat33& m) {
    return Mat33(Row3(m(0,0), m(0,1), m(0,2)),
                 Row3(m(1,0), m(1,1), m(1,2)),
                 Row3(m(2,0), m(2,1), m(2,2)));
}
static RBInertia toRBInertia(const MatInertia& i) {
    return RBInertia(toCDSMat33(i.toMat33()));
}
static MatInertia toMatInertia(const RBInertia& i) {
    return MatInertia(toMat33(i));
}
static MatRotation toMatRotation(const CDSMat33& m) {
    return reinterpret_cast<const MatRotation&>(toMat33(m));
}

static RBMassProperties toRBMassProperties(const Real& m, const Vec3& c, const MatInertia& i) {
    return RBMassProperties(m, toCDSVec3(c), toRBInertia(i));
}

static RBFrame toRBFrame(const Frame& f) {
    return RBFrame(toCDSMat33(f.getAxes().asMat33()), toCDSVec3(f.getOrigin()));
}
static Frame toFrame(const RBFrame& f) {
    return Frame(toMatRotation(f.getRot_RF()), toVec3(f.getLoc_RF()));
}
    // IVM SIMBODY INTERFACE //

IVMSimbodyInterface::IVMSimbodyInterface(const Multibody& m) : rep(0) {
    rep = new IVMSimbodyInterfaceRep(m);
    rep->setMyHandle(*this);
}

int IVMSimbodyInterface::getNBodies()     const {return rep->getRigidBodyTree().getNBodies();}
int IVMSimbodyInterface::getNParameters() const {return 0;}
int IVMSimbodyInterface::getNQ()          const {return rep->getRigidBodyTree().getDim();}
int IVMSimbodyInterface::getNU()          const {return rep->getRigidBodyTree().getDOF();}


Frame // TODO: this should return a reference to a cache entry
IVMSimbodyInterface::getBodyConfiguration(const State&, const Body& b) const
{
    const RBTreeMap& info = rep->getBodyInfo(b);
    const Frame& BR = info.getRefFrameInBody();
    const RigidBodyNode& n = rep->getRigidBodyTree().getRigidBodyNode(info.getRBIndex());
    const MatRotation GR = toMatRotation(n.getR_GB()); // B is really R
    const Vec3        OR_G = toVec3(n.getOB_G());
    return Frame(GR * ~BR.getAxes(), OR_G - BR.getOrigin());
}

SpatialVector // TODO: this should return a reference to a cache entry
IVMSimbodyInterface::getBodyVelocity(const State&, const Body& b) const
{
    return SpatialVector();
}

void IVMSimbodyInterface::applyPointForce(const State&, const Body& b, const Vec3& pt, const Vec3& frc, 
                     Array<SpatialVector>& bodyForces) const
{
}

void IVMSimbodyInterface::applyBodyTorque(const State&, const Body& b, const Vec3& trq, 
                     Array<SpatialVector>& bodyForces) const
{
}

void IVMSimbodyInterface::applyGravity(const State&, const Vec3& g, Array<SpatialVector>& bodyForces) const
{
}

void IVMSimbodyInterface::applyHingeForce(const State&, const Joint& j, int axis, const Real& frc, 
                     Vector& hingeForces) const
{
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
    mbs2tree.push_back(RBTreeMap(&Body::downcast(mbs["Ground"]),Frame(),Frame(),0,0,0));
    while (nxt < mbs2tree.size()) {
        for (size_t i=0; i<joints.size(); ++i) {
            if (!Body::getPlacementBody(joints[i]->getReferenceFrame())
                     .isSameSubsystem(mbs2tree[nxt].getBody())) continue;


            // Calculate the reference frame, which is located at origin of J
            // and aligned with its parent's reference frame.
            const Frame& BJ = joints[i]->getMovingFrame().getValue();        // on B
            const Frame& PJi = joints[i]->getReferenceFrame().getValue();    // on P
            const Frame BR(BJ.getAxes()*~PJi.getAxes(), BJ.getOrigin());

            mbs2tree.push_back(RBTreeMap(&Body::getPlacementBody(joints[i]->getMovingFrame()),
                                         BR, Frame(PJi.getAxes(), Vec3(0)), // RJ
                                         nxt, // index of parent
                                         joints[i],
                                         mbs2tree[nxt].getLevel() + 1));
        }
        ++nxt;
    }



    cout << "**** TREE ****" << endl;
    int nextStateOffset = 0;
    for (size_t i=0; i<mbs2tree.size(); ++i) {
        cout << mbs2tree[i].getLevel() << ": " << mbs2tree[i].getBody().getFullName() << endl;
        if (!mbs2tree[i].getLevel()) {
            RigidBodyNode* rb = RigidBodyNode::create(RBMassProperties(), RBFrame(), RBThisIsGround,
                                      false, false, nextStateOffset);
            const int rbIndex = tree.addGroundNode(rb);
            mbs2tree[i].setRBIndex(rbIndex);
            continue;
        }
        RBTreeMap& parentEntry = mbs2tree[mbs2tree[i].getParentIndex()];
        cout << Joint::getJointTypeName(mbs2tree[i].getJoint().getJointType())
             << " Joint: "  << mbs2tree[i].getJoint().getFullName() 
             << " Parent: " << parentEntry.getBody().getFullName() << endl;

        const Real& mass = mbs2tree[i].getBody().getMass().getValue();
        const Vec3& com  = mbs2tree[i].getBody().getMassCenter().getValue();
        const MatInertia& iner = mbs2tree[i].getBody().getInertia().getValue();
        cout << "mass=" << mass << " com=" << com << " iner=" << iner << endl;

        const Frame& fBJ = mbs2tree[i].getJoint().getMovingFrame().getValue();
        cout << "frame_BJ=" << fBJ << endl;
        cout << "frame_BR=" << mbs2tree[i].getRefFrameInBody() << endl;
        cout << "frame_RJ=" << mbs2tree[i].getJointFrameInRef() << endl;

        cout << "JointType=" << mbs2tree[i].getJoint().getJointType()
             << "  RBJointType=" << mapToRBJointType(mbs2tree[i].getJoint().getJointType())
             << endl;

        const int save = nextStateOffset;
        RigidBodyNode* rb = RigidBodyNode::create(
            toRBMassProperties(mass,com,iner),
            toRBFrame(mbs2tree[i].getJointFrameInRef()),
            mapToRBJointType(mbs2tree[i].getJoint().getJointType()),
            false,
            false,
            nextStateOffset);
        cout << "CREATED: states " << save << "-" << nextStateOffset-1 << endl;


        RigidBodyNode& parent = tree.updRigidBodyNode(parentEntry.getRBIndex());
        const int rbIndex = tree.addRigidBodyNode(parent,RBFrame()/*XXX*/,rb);
        mbs2tree[i].setRBIndex(rbIndex);
    }

    tree.finishConstruction(1e-6, 0);

    std::cout << "*** RigidBodyTree:" << std::endl << tree << std::endl;
}
   
/*static*/ RBJointType
IVMSimbodyInterfaceRep::mapToRBJointType(Joint::JointType jt) {
    switch (jt) {
    case Joint::UnknownJointType:   return RBUnknownJointType;
    case Joint::ThisIsGround:       return RBThisIsGround;
    case Joint::Weld:               return RBWeldJoint;
    case Joint::Torsion:            return RBTorsionJoint;  // aka PinJoint
    case Joint::Sliding:            return RBSlidingJoint;
    case Joint::Universal:          return RBUJoint;
    case Joint::Cylinder:           return RBCylinderJoint;
    case Joint::Planar:             return RBPlanarJoint;
    case Joint::Gimbal:             return RBGimbalJoint;
    case Joint::Orientation:        return RBOrientationJoint; // aka BallJoint
    case Joint::Cartesian:          return RBCartesianJoint;
    case Joint::FreeLine:           return RBFreeLineJoint;
    case Joint::Free:               return RBFreeJoint;
    default: assert(false);
    }
    //NOTREACHED
    return RBUnknownJointType;
}

