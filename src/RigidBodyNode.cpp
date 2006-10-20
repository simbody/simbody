/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors: Derived from IVM code written by Charles Schwieters.
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
 * This file contains all the multibody mechanics code that involves a single body and
 * its inboard joint, that is, one node in the multibody tree.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

#include "RigidBodyTree.h"
#include "RigidBodyNode.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

//////////////////////////////////////////////
// Implementation of RigidBodyNode methods. //
//////////////////////////////////////////////

void RigidBodyNode::addChild(RigidBodyNode* child) {
    children.push_back( child );
}

//
// Calc posCM, mass, Mk
//      phi, inertia
// Should be calc'd from base to tip.
// We depend on transforms X_PB and X_GB being available.
void RigidBodyNode::calcJointIndependentKinematicsPos(
    SBPositionCache&   cc) const
{
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 T_PB_G = getX_GP(cc).R() * getX_PB(cc).T();

    // The Phi matrix conveniently performs child-to-parent shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // shifting for velocities.
    updPhi(cc) = PhiMatrix(T_PB_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    updInertia_OB_G(cc) = getInertia_OB_B().changeAxes(~getX_GB(cc).R());
    updCB_G(cc)         = getX_GB(cc).R()*getCOM_B();

    updCOM_G(cc) = getX_GB(cc).T() + getCB_G(cc);

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    // Note: we need to calculate this now so that we'll be able to calculate
    // kinetic energy without going past the Motion stage.
    const Mat33 offDiag = getMass()*crossMat(getCB_G(cc));
    updMk(cc) = SpatialMat( getInertia_OB_G(cc).toMat33() ,     offDiag ,
                                   -offDiag             , getMass()*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (sVel). This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBPositionCache& cc,
    SBVelocityCache&       mc) const
{
    updV_GB(mc) = ~getPhi(cc)*parent->getV_GB(mc) + getV_PB_G(mc);
}

Real RigidBodyNode::calcKineticEnergy(
    const SBPositionCache& cc,
    const SBVelocityCache& mc) const 
{
    const Real ret = dot(getV_GB(mc) , getMk(cc)*getV_GB(mc));
    return 0.5*ret;
}

// Calculate velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and coriolis acceleration.
// This routine expects that all spatial velocities & spatial inertias are
// already available.
// Must be called base to tip.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(
    const SBPositionCache& cc,
    const SBVelocityCache& mc,
    SBDynamicsCache&       dc) const
{
    if (nodeNum == 0) { // ground, just in case
        updGyroscopicForce(dc)           = SpatialVec(Vec3(0), Vec3(0));
        updCoriolisAcceleration(dc)      = SpatialVec(Vec3(0), Vec3(0));
        updTotalCoriolisAcceleration(dc) = SpatialVec(Vec3(0), Vec3(0));
        updCentrifugalForces(dc)         = SpatialVec(Vec3(0), Vec3(0));
        updTotalCentrifugalForces(dc)    = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    const Vec3& omega = getV_GB(mc)[0];  // spatial angular velocity
    const Vec3& vel   = getV_GB(mc)[1];  // spatial linear velocity

    updGyroscopicForce(dc) = 
        SpatialVec(    omega % (getInertia_OB_G(cc)*omega),     // gyroscopic moment
                    getMass()*(omega % (omega % getCB_G(cc)))); // gyroscopic force

    // Parent velocity.
    const Vec3& pOmega = parent->getV_GB(mc)[0];
    const Vec3& pVel   = parent->getV_GB(mc)[1];

    // Calc a: coriolis acceleration.
    // Sherm TODO: See Schwieters & Clore Eq. [16], which uses *this*
    // body's omega (w_k) while the code below uses the *parent* pOmega
    // (w_k-1). Compare with Jain, Vaidehi, & Rodriguez 1991, Eq. 4.4,
    // esp. the following paragraph saying that the cross product
    // will be the same whether we use omega or pOmega in the 2nd term,
    // because they can only differ along H which is constant between P & B.
    // (caution: JV&R number backwards so the parent is w_k+1 there).
    // Is that also true in the first term? I.e., can we use pOmega, omega,
    // or both below and get the same answers? Anyway the code below is
    // consistent with JV&R (and works) but not obviously consistent with S&C
    // (which doesn't).

    // (sherm 060728) My guess is that I introduced the need for the cross-joint
    // term here when I separated the body frame from the joint frame.
    // Dan Rosenthal's description in abandoned patent app 10/053,348
    // paragraph [0097] is easier to understand than S&C or JV&R, although
    // it is also missing the cross-joint term.

    // Note: (vel-pVel) is the total relative velocity; V_PB_G is that 
    // portion due just to this joint's u's.

    // (sherm 060912) Caution: using T_JbB_G here is assuming that all the joints do their
    // rotations in the Jb frame *followed* by translations in the J frame.
    // I *think* you could switch the joint to work the other way and change to T_JB_G
    // in the joints and here, but I haven't tried it. TODO: this should be written
    // in terms of the H matrix somehow so that this computation and the joint definition
    // do not have to be synchronized this way.

    const Vec3 T_PB_G = getX_GB(cc).T() - getX_GP(cc).T();
    const Vec3 T_JbB_G = T_PB_G - getX_GP(cc).R()*getX_PJb().T();
    //const Vec3 T_JB_G  = -getX_GB(cc).R()*getX_BJ().T();
    const Vec3 w_PB_G = getV_PB_G(mc)[0];
    const Vec3 v_PB_G = getV_PB_G(mc)[1]; // includes w_PB_G % T_JbB_G term

    const SpatialVec wwr = SpatialVec(Vec3(0), 
                  pOmega % (pOmega % T_PB_G) - w_PB_G % (w_PB_G % T_JbB_G));
    const SpatialVec wv  = SpatialVec( pOmega % w_PB_G, 
                                       2*(omega % v_PB_G)); 
    
    updCoriolisAcceleration(dc) = wwr + wv;

    updTotalCoriolisAcceleration(dc) =
        ~getPhi(cc) * parent->getTotalCoriolisAcceleration(dc)
        + getCoriolisAcceleration(dc); // just calculated above

    updCentrifugalForces(dc) =
        getP(dc) * getCoriolisAcceleration(dc) + getGyroscopicForce(dc);

    updTotalCentrifugalForces(dc) = 
        getP(dc) * getTotalCoriolisAcceleration(dc) + getGyroscopicForce(dc);

}

void RigidBodyNode::nodeDump(std::ostream& o) const {
    o << "NODE DUMP level=" << level << " type=" << type() << std::endl;
    nodeSpecDump(o);
    o << "END OF NODE type=" << type() << std::endl;
}

std::ostream& operator<<(std::ostream& o, const RigidBodyNode& node) {
    node.nodeDump(o);
    return o;
}


////////////////////////////////////////////////
// Define classes derived from RigidBodyNode. //
////////////////////////////////////////////////

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class RBGroundBody : public RigidBodyNode {
public:
    RBGroundBody() // TODO: should set mass properties to infinity
      : RigidBodyNode(MassProperties(),Transform(),Transform()) {}
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }
    /*virtual*/int getDOF() const { return 0; }
    /*virtual*/int getMaxNQ() const { return 0; }
    /*virtual*/int getNQ(const SBModelVars&) const { return 0; }

    /*virtual*/void calcZ(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const SpatialVec& spatialForce,
        SBAccelerationCache&               ) const {} 

    /*virtual*/void calcYOutward(
        const SBPositionCache& cc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/void calcAccel(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        const SBDynamicsCache& dc,
        SBAccelerationCache&   rc,
        Vector&                udot,
        Vector&                qdotdot) const {}

    /*virtual*/void realizeModel(
        const SBModelVars& mv,
        SBModelCache&      mc) const {}

    /*virtual*/void realizeInstance(
        const SBModelVars&    mv,
        const SBInstanceVars& pv,
        SBInstanceCache&      pc) const {}

    /*virtual*/void realizePosition(
        const SBModelVars&  mv,
        const Vector&       q,
        SBPositionCache&    cc) const {}

    /*virtual*/void realizeVelocity(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        SBVelocityCache&       mc,
        Vector&                qdot) const {}

    /*virtual*/void setVelFromSVel(
        const SBPositionCache& cc, 
        const SBVelocityCache& mc,
        const SpatialVec&      sVel, 
        Vector&                u) const {}

    /*virtual*/bool enforceQuaternionConstraints(
        const SBModelVars& mv,
        Vector&             q) const {return false;}
    
    /*virtual*/void calcArticulatedBodyInertiasInward(
        const SBPositionCache& cc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/ void calcInternalGradientFromSpatial(
        const SBPositionCache&      cc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const { }

    /*virtual*/ void calcEquivalentJointForces(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector&                    jointForces) const 
    { 
        allZ[0] = bodyForces[0];
    }

    /*virtual*/void calcUDotPass1Inward(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const
    {
        allZ[0] = -bodyForces[0]; // TODO sign is weird
        allGepsilon[0] = SpatialVec(Vec3(0), Vec3(0));
    } 
    /*virtual*/void calcUDotPass2Outward(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&              epsilonTmp,
        Vector_<SpatialVec>&       allA_GB,
        Vector&                    allUDot) const
    {
        allA_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    }

    /*virtual*/void calcMInverseFPass1Inward(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&              f,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const
    {
        allZ[0] = SpatialVec(Vec3(0), Vec3(0));
        allGepsilon[0] = SpatialVec(Vec3(0), Vec3(0));
    } 

    /*virtual*/void calcMInverseFPass2Outward(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const
    {
        allA_GB[0] = SpatialVec(Vec3(0), Vec3(0));
    }

    /*virtual*/void setDefaultModelValues(const SBTopologyCache&, 
                                          SBModelVars& v) const
    {
        v.prescribed[0] = true; // ground's motion is prescribed to zero
    }

    // /*virtual*/ const SpatialRow& getHRow(int i) const;
};

template<int dof>
class RigidBodyNodeSpec : public RigidBodyNode {
public:
    RigidBodyNodeSpec(const MassProperties& mProps_B,
                      const Transform&      X_PJb,
                      const Transform&      X_BJ,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNode(mProps_B, X_PJb, X_BJ)
    {
        // don't call any virtual methods in here!
        uIndex   = nextUSlot;
        uSqIndex = nextUSqSlot;
        qIndex   = nextQSlot;
    }

    void updateSlots(int& nextUSlot, int& nextUSqSlot, int& nextQSlot) {
        // OK to call virtual method here.
        nextUSlot   += getDOF();
        nextUSqSlot += getDOF()*getDOF();
        nextQSlot   += getMaxNQ();
    }

    virtual ~RigidBodyNodeSpec() {}

    // This is the type of the joint transition matrix H.
    typedef Mat<dof,2,Row3,1,2> HType;


    // The following routines calculate joint-specific position kinematic
    // quantities. They may assume that *all* position kinematics (not just
    // joint-specific) has been done for the parent, and that the position
    // state variables q are available. Each routine may also assume that the
    // previous routines have been called, in the order below.
    // The routines are structured as operators -- they use the State but
    // do not change anything in it, including the cache. Instead they are
    // passed argument to write their results into. In practice, these
    // arguments will typically be in the State cache (see below).

    /// This mandatory routine performs expensive floating point operations sin,cos,sqrt
    /// in one place so we don't end up repeating them. sin&cos are used only
    /// for joints which have angular coordinates, and qnorm is only for joints
    /// which are using quaternions. Other joints can provide a null routine.
    /// Each of the passed-in Vectors is a "q-like" object, that is, allocated
    /// to the bodies in a manner parallel to the q state variable.
    virtual void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const=0;

    /// This mandatory routine calculates the across-joint transform X_JbJ generated
    /// by the current q values. This may depend on sines & cosines or normalized
    /// quaternions already being available in the State cache.
    virtual void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const=0;

    /// This routine is NOT joint specific, but cannot be called until the across-joint
    /// transform X_JbJ has been calculated and is available in the State cache.
    void calcBodyTransforms(
        const SBPositionCache& cc, 
        Transform&             X_PB, 
        Transform&             X_GB) const 
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_JbJ = getX_JbJ(cc); // just calculated
        const Transform& X_GP  = getX_GP(cc);  // already calculated

        X_PB = X_PJb * X_JbJ * ~X_BJ; // TODO: precalculate X_JB
        X_GB = X_GP * X_PB;
    }

    /// This mandatory routine calcluates the joint transition matrix H, giving the
    /// change of *spatial* velocity induced by the generalized speeds u for this
    /// joint. It may depend on X_PB and X_GB having been calculated already.
    virtual void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const=0;

    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (u) are available. The
    /// quanitites that must be computed are:
    ///   V_PB_G  relative velocity of B in P, expr. in G
    /// The code is the same for all joints, although parametrized by dof.
    void calcJointKinematicsVel(
        const SBPositionCache& cc,
        const Vector&          u,
        SBVelocityCache&       mc) const 
    {
        updV_PB_G(mc) = ~getH(cc) * fromU(u);
    }

    // These next two routines are options, but if you supply one you
    // must supply the other. (That is, ball-containing joints provide
    // both of these routines.)
    virtual void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        Vector&                qdot) const
    {
        toQ(qdot) = fromU(u);        // default is qdot=u
    }

    virtual void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const
    {
        toQ(qdotdot) = fromU(udot);  // default is qdotdot=udot
    }

    void realizeModel(
        const SBModelVars& mv,
        SBModelCache&      mc) const 
    {
    }

    void realizeInstance(
        const SBModelVars&    mv,
        const SBInstanceVars& pv,
        SBInstanceCache&      pc) const
    {
    }

    // Set a new configuration and calculate the consequent kinematics.
    // Must call base-to-tip.
    void realizePosition(
        const SBModelVars& mv,
        const Vector&      q,
        SBPositionCache&   cc) const 
    {
        calcJointSinCosQNorm(mv, q, cc.sq, cc.cq, cc.qnorm);
        calcAcrossJointTransform (mv, q, updX_JbJ(cc));
        calcBodyTransforms       (cc, updX_PB(cc), updX_GB(cc));
        calcJointTransitionMatrix(cc, updH(cc));
        calcJointIndependentKinematicsPos(cc);
    }

    // Set new velocities for the current configuration, and calculate
    // all the velocity-dependent terms. Must call base-to-tip.
    void realizeVelocity(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        SBVelocityCache&       mc,
        Vector&                qdot) const 
    {
        calcQDot(mv,q,cc,u,qdot);
        calcJointKinematicsVel(cc,u,mc);
        calcJointIndependentKinematicsVel(cc,mc);
    }

    // This is a dynamics-stage calculation and must be called tip-to-base (inward).
    void calcArticulatedBodyInertiasInward(
        const SBPositionCache& cc,
        SBDynamicsCache&       dc) const;

    // calcJointIndependentDynamicsVel() must be called after ArticulatedBodyInertias.

    // This dynamics-stage calculation is needed for handling constraints. It
    // must be called base-to-tip (outward);
    void calcYOutward(
        const SBPositionCache& cc,
        SBDynamicsCache&       dc) const;

    // These routines give each node a chance to set appropriate defaults in a piece
    // of the state corresponding to a particular stage. Default implementations here
    // assume non-ball joint; override if necessary.
    virtual void setDefaultModelValues (const SBTopologyCache&, SBModelVars&)  const {}
    virtual void setDefaultInstanceValues(const SBModelVars&, SBInstanceVars&) const {}
    virtual void setDefaultTimeValues     (const SBModelVars&, SBTimeVars&)      const {}

    virtual void setDefaultPositionValues(const SBModelVars& s, Vector& q) const 
    {
        toQ(q) = 0.;
    }
    virtual void setDefaultVelocityValues(const SBModelVars&, Vector& u) const 
    {
        toU(u) = 0.;
    }
    virtual void setDefaultDynamicsValues(const SBModelVars&, SBDynamicsVars&) const {}
    virtual void setDefaultAccelerationValues(const SBModelVars&, 
                                          SBDynamicsVars& v) const {}

    // setQ and setU extract this node's values from the supplied
    // q-sized or u-sized array and put them in the corresponding
    // locations in the output variable. Joints which need quaternions should
    // override setQ to copy the extra q.
    virtual void setQ(
        const SBModelVars& mv, 
        const Vector&      qIn, 
        Vector&            q) const
    {
        toQ(q) = fromQ(qIn);
    }

    virtual void setU(
        const SBModelVars& mv, 
        const Vector&      uIn, 
        Vector&            u) const
    {
        toU(u) = fromU(uIn);
    }

    int          getDOF()            const { return dof; }
    virtual int  getMaxNQ()          const { return dof; } // maxNQ can be larger than dof
    virtual int  getNQ(const SBModelVars&) const { return dof; } // DOF <= NQ <= maxNQ

    virtual void setVelFromSVel(
        const SBPositionCache& cc, 
        const SBVelocityCache& mc,
        const SpatialVec&      sVel, 
        Vector&               u) const;

    // Return true if any change is made to the output variable.
    virtual bool enforceQuaternionConstraints(
        const SBModelVars& mv,
        Vector&            q) const 
    {
        return false;
    }

    const SpatialRow& getHRow(const SBPositionCache& cc, int i) const {
        return getH(cc)[i];
    }

    // Access to body-oriented state and cache entries is the same for all nodes,
    // and joint oriented access is almost the same but parametrized by dof. There is a special
    // case for quaternions because they use an extra state variable, and although we don't
    // have to we make special scalar routines available for 1-dof joints. Note that all State access
    // routines are inline, not virtual, so the cost is just an indirection and an index.

    // General joint-dependent select-my-goodies-from-the-pool routines.
    const Vec<dof>&     fromQ  (const Vector& q)   const {return Vec<dof>::getAs(&q[qIndex]);}
    Vec<dof>&           toQ    (      Vector& q)   const {return Vec<dof>::updAs(&q[qIndex]);}
    const Vec<dof>&     fromU  (const Vector& u)   const {return Vec<dof>::getAs(&u[uIndex]);}
    Vec<dof>&           toU    (      Vector& u)   const {return Vec<dof>::updAs(&u[uIndex]);}
    const Mat<dof,dof>& fromUSq(const Vector& uSq) const {return Mat<dof,dof>::getAs(&uSq[uSqIndex]);}
    Mat<dof,dof>&       toUSq  (      Vector& uSq) const {return Mat<dof,dof>::updAs(&uSq[uSqIndex]);}

    // Same, but specialized for the common case where dof=1 so everything is scalar.
    const Real& from1Q  (const Vector& q)   const {return q[qIndex];}
    Real&       to1Q    (      Vector& q)   const {return q[qIndex];}
    const Real& from1U  (const Vector& u)   const {return u[uIndex];}
    Real&       to1U    (      Vector& u)   const {return u[uIndex];}
    const Real& from1USq(const Vector& uSq) const {return uSq[uSqIndex];}
    Real&       to1USq  (      Vector& uSq) const {return uSq[uSqIndex];}

    // Same, specialized for quaternions. We're assuming that the quaternion comes first in the coordinates.
    const Vec4& fromQuat(const Vector& q) const {return Vec4::getAs(&q[qIndex]);}
    Vec4&       toQuat  (      Vector& q) const {return Vec4::updAs(&q[qIndex]);}

    // Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
    const Vec3& fromQVec3(const Vector& q, int offs) const {return Vec3::getAs(&q[qIndex+offs]);}
    Vec3&       toQVec3  (      Vector& q, int offs) const {return Vec3::updAs(&q[qIndex+offs]);}
    const Vec3& fromUVec3(const Vector& u, int offs) const {return Vec3::getAs(&u[uIndex+offs]);}
    Vec3&       toUVec3  (      Vector& u, int offs) const {return Vec3::updAs(&u[uIndex+offs]);}

    // Applications of the above extraction routines to particular interesting items in the State. Note
    // that you can't use these for quaternions since they extract "dof" items.

    // Applied forces from cache.
    const Vec<dof>&   getAppliedJointForce(const SBDynamicsCache& dc) const 
        {return fromU(dc.appliedMobilityForces);}
    //TODO
    const Vec<dof>&   getPrescribedUdot   (const SBDynamicsVars& dv) const 
        {return fromU(dv.prescribedUdot);}

    // Special case state access for 1-dof joints
    const Real& get1AppliedJointForce(const SBDynamicsCache& dc) const {return from1U(dc.appliedMobilityForces);}
    const Real& get1PrescribedUdot   (const SBDynamicsVars& dv) const {return from1U(dv.prescribedUdot);}

    // Cache entries (cache is mutable in a const State)

        // Configuration

    // TODO: should store as H or else always reference Ht
    const Mat<dof,2,Row3,1,2>& getH(const SBPositionCache& cc) const
      { return ~Mat<2,dof,Vec3>::getAs(&cc.storageForHt(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH(SBPositionCache& cc) const
      { return ~Mat<2,dof,Vec3>::updAs(&cc.storageForHt(0,uIndex)); }

    // These are sines and cosines of angular qs. The rest of the slots are garbage.
    const Vec<dof>&   getSinQ (const SBPositionCache& cc) const {return fromQ (cc.sq);}
    Vec<dof>&         updSinQ (SBPositionCache&       cc) const {return toQ   (cc.sq);}
    const Real&       get1SinQ(const SBPositionCache& cc) const {return from1Q(cc.sq);}
    Real&             upd1SinQ(SBPositionCache&       cc) const {return to1Q  (cc.sq);}

    const Vec<dof>&   getCosQ (const SBPositionCache& cc) const {return fromQ (cc.cq);}
    Vec<dof>&         updCosQ (SBPositionCache&       cc) const {return toQ   (cc.cq);}
    const Real&       get1CosQ(const SBPositionCache& cc) const {return from1Q(cc.cq);}
    Real&             upd1CosQ(SBPositionCache&       cc) const {return to1Q  (cc.cq);}

    // These are normalized quaternions in slots for balls. Everything else is garbage.
    const Vec4&       getQNorm(const SBPositionCache& cc) const {return fromQuat(cc.qnorm);}
    Vec4&             updQNorm(SBPositionCache&       cc) const {return toQuat  (cc.qnorm);}

        // Motion

        // Dynamics
    const Mat<dof,dof>& getD(const SBDynamicsCache& dc) const {return fromUSq(dc.storageForD);}
    Mat<dof,dof>&       updD(SBDynamicsCache&       dc) const {return toUSq  (dc.storageForD);}

    const Mat<dof,dof>& getDI(const SBDynamicsCache& dc) const {return fromUSq(dc.storageForDI);}
    Mat<dof,dof>&       updDI(SBDynamicsCache&       dc) const {return toUSq  (dc.storageForDI);}

    const Mat<2,dof,Vec3>& getG(const SBDynamicsCache& dc) const
      { return Mat<2,dof,Vec3>::getAs(&dc.storageForG(0,uIndex)); }
    Mat<2,dof,Vec3>&       updG(SBDynamicsCache&       dc) const
      { return Mat<2,dof,Vec3>::updAs(&dc.storageForG(0,uIndex)); }

        // Reaction

    const Vec<dof>&   getNetHingeForce (const SBAccelerationCache& rc) const {return fromU (rc.netHingeForces);}
    Vec<dof>&         updNetHingeForce (SBAccelerationCache&       rc) const {return toU   (rc.netHingeForces);}
    const Real&       get1NetHingeForce(const SBAccelerationCache& rc) const {return from1U(rc.netHingeForces);}
    Real&             upd1NetHingeForce(SBAccelerationCache&       rc) const {return to1U  (rc.netHingeForces);}


    const Vec<dof>&   getNu (const SBAccelerationCache& rc) const {return fromU (rc.nu);}
    Vec<dof>&         updNu (SBAccelerationCache&       rc) const {return toU   (rc.nu);}
    const Real&       get1Nu(const SBAccelerationCache& rc) const {return from1U(rc.nu);}
    Real&             upd1Nu(SBAccelerationCache&       rc) const {return to1U  (rc.nu);}

    const Vec<dof>&   getEpsilon (const SBAccelerationCache& rc) const {return fromU (rc.epsilon);}
    Vec<dof>&         updEpsilon (SBAccelerationCache&       rc) const {return toU   (rc.epsilon);}
    const Real&       get1Epsilon(const SBAccelerationCache& rc) const {return from1U(rc.epsilon);}
    Real&             upd1Epsilon(SBAccelerationCache&       rc) const {return to1U  (rc.epsilon);}

    void calcZ(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const SpatialVec& spatialForce,
        SBAccelerationCache&               ) const;

    void calcAccel(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        const SBDynamicsCache& dc,
        SBAccelerationCache&   rc,
        Vector&                udot,
        Vector&                qdotdot) const;

    void calcInternalGradientFromSpatial(
        const SBPositionCache&      cc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const;

    void calcEquivalentJointForces(
        const SBPositionCache&      cc,
        const SBDynamicsCache&      dc,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector&                     jointForces) const;

    void calcUDotPass1Inward(
        const SBPositionCache&      cc,
        const SBDynamicsCache&      dc,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const;

    void calcUDotPass2Outward(
        const SBPositionCache&      cc,
        const SBDynamicsCache&      dc,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const;

    void calcMInverseFPass1Inward(
        const SBPositionCache&      cc,
        const SBDynamicsCache&      dc,
        const Vector&               f,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const;

    void calcMInverseFPass2Outward(
        const SBPositionCache&      cc,
        const SBDynamicsCache&      dc,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const;

    /*
    void nodeSpecDump(std::ostream& o, const State& s) const {
        o << "uIndex=" << uIndex << " mass=" << getMass() 
            << " COM_G=" << getCOM_G(s) << std::endl;
        o << "inertia_OB_G=" << getInertia_OB_G(s) << std::endl;
        o << "H=" << getH(s) << std::endl;
        o << "SVel=" << getV_GB(s) << std::endl;
        o << "a=" << getCoriolisAcceleration(s) << std::endl;
        o << "b=" << getGyroscopicForce(s) << std::endl;
        o << "Th  =" << getQ(s) << std::endl;
        o << "dTh =" << getU(s) << std::endl;
        o << "ddTh=" << getUDot(s) << std::endl;
        o << "SAcc=" << getA_GB(s) << std::endl;
    }
    */
};

//////////////////////////////////////////
// Derived classes for each joint type. //
//////////////////////////////////////////


/**
 * Translate (Cartesian) joint. This provides three degrees of
 * translational freedom which is suitable (e.g.) for connecting a
 * free atom to ground. The Cartesian directions are the axes of
 * the parent body's Jb frame, with J=Jb when all 3 coords are 0,
 * and the orientation of J in Jb is 0 forever.
 */
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const MassProperties& mProps_B,
                    const Transform&      X_PJb,
                    const Transform&      X_BJ,
                    int&                  nextUSlot,
                    int&                  nextUSqSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setMobilizerPosition(const SBModelVars&, const Transform& X_JbJ,
                                   Vector& q) const 
    {
        toQ(q) = X_JbJ.T();
    }
    void setMobilizerVelocity(const SBModelVars&, const SpatialVec& V_JbJ,
                              Vector& u) const
    {
        toU(u) = V_JbJ[1];
    }

    // This is required but does nothing here since there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_JbJ = Transform(Rotation(), fromQ(q));
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const 
    {
        const Transform& X_PJb   = getX_PJb();      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP    = getX_GP(cc); // parent orientation in ground

        // Note that H is spatial. The current spatial directions for our qs are
        // the axes of the Jb frame expressed in Ground.
        const Rotation R_GJb = X_GP.R()*X_PJb.R();
        H[0] = SpatialRow( Row3(0), ~R_GJb.x() );
        H[1] = SpatialRow( Row3(0), ~R_GJb.y() );
        H[2] = SpatialRow( Row3(0), ~R_GJb.z() );
    }
};



/**
 * Sliding joint (1 dof translation). The translation is along the x
 * axis of the parent body's Jb frame, with J=Jb when the coordinate
 * is zero and the orientation of J in Jb frozen at 0 forever.
 */
class RBNodeSlider : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "slider"; }

    RBNodeSlider(const MassProperties& mProps_B,
                 const Transform&      X_PJb,
                 const Transform&      X_BJ,
                 int&                  nextUSlot,
                 int&                  nextUSqSlot,
                 int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }
        // Implementations of virtual methods.

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const { }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        // Translation vector q is expressed in Jb (and J since they have same orientation).
        // A sliding joint can't change orientation, and only translates along x. 
        X_JbJ = Transform(Rotation(), Vec3(from1Q(q),0,0));
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_PJb   = getX_PJb();      // fixed config of Jb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP    = getX_GP(cc); // parent configuration in ground

        // Note that H is spatial. The current spatial directions for our q is
        // the x axis of the Jb frame expressed in Ground.
        const Vec3 x_GJb = X_GP.R()*X_PJb.x();
        H[0] = SpatialRow( Row3(0), ~x_GJb );
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis, the z axis of the parent's Jb frame, which is 
 * aligned forever with the z axis of the body's J frame. In addition, the
 * origin points of J and Jb are identical forever.
 */
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const MassProperties& mProps_B,
                  const Transform&      X_PJb,
                  const Transform&      X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        const Real& angle = from1Q(q); // angular coordinate
        to1Q(sine)    = std::sin(angle);
        to1Q(cosine)  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        const Real& theta  = from1Q(q);    // angular coordinate

        // We're only updating the orientation here because a torsion joint
        // can't translate (it is defined as a rotation about the z axis).
        X_JbJ.updR().setToRotationAboutZ(theta);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // Calc H matrix in space-fixed coords.
        // This works because the joint z axis is the same in J & Jb
        // since that's what we rotate around.
        const Vec3 z_G = X_GP.R() * X_PJb.z();
        H[0] = SpatialRow( ~z_G, ~(z_G % T_JB_G) );
    }
};

/**
 * This is a "cylinder" joint, meaning one degree of rotational freedom
 * about a particular axis, and one degree of translational freedom
 * along the same axis. For molecules you can think of this as a combination
 * of torsion and bond stretch. The axis used is the z axis of the parent's
 * Jb frame, which is aligned forever with the z axis of the body's J frame.
 * In addition, the origin points of J and Jb are separated only along the
 * z axis; i.e., they have the same x & y coords in the Jb frame. The two
 * generalized coordinates are the rotation and the translation, in that order.
 */
class RBNodeCylinder : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "cylinder"; }

    RBNodeCylinder(const MassProperties& mProps_B,
                   const Transform&      X_PJb,
                   const Transform&      X_BJ,
                   int&                  nextUSlot,
                   int&                  nextUSqSlot,
                   int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        const Real& angle = fromQ(q)[0];
        toQ(sine)[0]    = std::sin(angle);
        toQ(cosine)[0]  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        const Vec2& coords  = fromQ(q);

        X_JbJ.updR().setToRotationAboutZ(coords[0]);
        X_JbJ.updT() = Vec3(0,0,coords[1]);
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // Calc H matrix in space-fixed coords.
        // This works because the joint z axis is the same in J & Jb
        // since that's what we rotate around.
        const Vec3 z_GJb = X_GP.R()*X_PJb.z();
        H[0] = SpatialRow( ~z_GJb, ~(z_GJb % T_JB_G) );
        H[1] = SpatialRow( Row3(0), ~z_GJb );
    }
};


/**
 * This is a "bend-stretch" joint, meaning one degree of rotational freedom
 * about a particular axis, and one degree of translational freedom
 * along a perpendicular axis. The z axis of the parent's Jb frame is 
 * used for rotation (and that is always aligned with the J frame z axis).
 * The x axis of the *J* frame is used for translation; that is, first
 * we rotate around z, which moves J's x with respect to Jb's x. Then
 * we slide along the rotated x axis. The two
 * generalized coordinates are the rotation and the translation, in that order.
 */
class RBNodeBendStretch : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "bendstretch"; }

    RBNodeBendStretch(const MassProperties& mProps_B,
                      const Transform&      X_PJb,
                      const Transform&      X_BJ,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        const Real& angle = fromQ(q)[0];
        toQ(sine)[0]    = std::sin(angle);
        toQ(cosine)[0]  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        const Vec2& coords  = fromQ(q);    // angular coordinate

        X_JbJ.updR().setToRotationAboutZ(coords[0]);
        X_JbJ.updT() = X_JbJ.R()*Vec3(coords[1],0,0); // because translation is in J frame
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Transform X_GJb = X_GP * X_PJb;
        const Transform X_GJ  = X_GB * X_BJ;

        const Vec3 T_JbB_G = X_GB.T() - X_GJb.T();
        //const Vec3 T_JB_G = -X_GB.R()*X_BJ.T();

        // Rotate around Jb's z axis, *then* translate along J's new x axis.
        H[0] = SpatialRow( ~X_GJb.z(), ~(X_GJb.z() % T_JbB_G) );
        H[1] = SpatialRow( Row3(0), ~X_GJ.x());
    }
};


/**
 * U-joint like joint type which allows rotation about the two axes
 * perpendicular to zDir. This is appropriate for diatoms and for allowing 
 * torsion+bond angle bending.
 */
class RBNodeRotate2 : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "rotate2"; }

    RBNodeRotate2(const MassProperties& mProps_B,
                  const Transform&      X_PJb,
                  const Transform&      X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        const Vec2& a = fromQ(q); // angular coordinates
        toQ(sine)   = Vec2(std::sin(a[0]), std::sin(a[1]));
        toQ(cosine) = Vec2(std::cos(a[0]), std::cos(a[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        const Vec2& angles  = fromQ(q); // angular coordinates

        // We're only updating the orientation here because a U-joint
        // can't translate.
        X_JbJ.updR().setToSpaceFixed12(angles);
        X_JbJ.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The coordinates are defined in the space-fixed (that is, Jb) frame, so
        // the orientation of Jb in ground gives the instantaneous spatial 
        // meaning of the coordinates.
        const Rotation R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
    }
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis.
 * The coordinate definitions are a combination of a rotate2 joint and a
 * Cartesian joint. The first 2 are rotational, the next 3 are translations.
 * However, the rotations don't affect the translations.
 */
class RBNodeTranslateRotate2 : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "diatom"; }

    RBNodeTranslateRotate2(const MassProperties& mProps_B,
                           const Transform&      X_PJb,
                           const Transform&      X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    // TODO: partial implementation; just translation
    void setMobilizerPosition(const SBModelVars&, const Transform& X_JbJ,
                              Vector& q) const 
    {
        toQ(q).updSubVec<3>(2) = X_JbJ.T();
    }
    // TODO: partial implementation; just translation
    void setMobilizerVelocity(const SBModelVars&, const SpatialVec& V_JbJ,
                              Vector& u) const
    {
        toU(u).updSubVec<3>(2) = V_JbJ[1];
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        const Vec2& a = fromQ(q).getSubVec<2>(0); // angular coordinates
        toQ(sine).updSubVec<2>(0)   = Vec2(std::sin(a[0]), std::sin(a[1]));
        toQ(cosine).updSubVec<2>(0) = Vec2(std::cos(a[0]), std::cos(a[1]));
        // no quaternions
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const 
    {
        const Vec<5>& coords = fromQ(q);     // joint coordinates
        const Vec<2>& angles = coords.getSubVec<2>(0);

        //X_JbJ.updR().setToSpaceFixed12(coords.getSubVec<2>(0));
        X_JbJ.updR() = Rotation::aboutXThenNewY(angles[0], angles[1]);
        //X_JbJ.updT() = X_JbJ.R()*coords.getSubVec<3>(2); // because translation is in J
        X_JbJ.updT() = coords.getSubVec<3>(2); // because translation is in Jb
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Transform X_GJb = X_GP * X_PJb;
        const Transform X_GJ  = X_GB * X_BJ;

        const Vec3 T_JB_G  = X_GB.T() - X_GJ.T();
        const Vec3 T_JbB_G = X_GB.T() - X_GJb.T();

        // The rotational coordinates are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        // *Then* we translate along the new J frame axes.

        //H[0] = SpatialRow(~X_GJb.x(), ~(X_GJb.x() % T_JbB_G));
        //H[1] = SpatialRow(~X_GJb.y(), ~(X_GJb.y() % T_JbB_G));
        H[0] = SpatialRow(~X_GJ.x(), ~(X_GJ.x() % T_JbB_G));
        H[1] = SpatialRow(~X_GJ.y(), ~(X_GJ.y() % T_JbB_G));
        H[2] = SpatialRow(  Row3(0) ,     ~X_GJb.x());
        H[3] = SpatialRow(  Row3(0) ,     ~X_GJb.y());
        H[4] = SpatialRow(  Row3(0) ,     ~X_GJb.z());    
    }
};

/// Ball joint. This provides three degrees of rotational freedom, i.e.,
/// unrestricted orientation of the body's J frame in the parent's Jb frame.
/// The 3 u's are the cross-joint angular velocity vector of J in Jb, and
/// udot's are the angular velocity time derivative. The q's, however, are
/// either 3 Euler angles in a 3-2-1 body-fixed sequence, or 4 quaternions.
/// In that case we calculate either 3 or 4 qdots from the u's.
class RBNodeRotate3 : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "rotate3"; }

    RBNodeRotate3(const MassProperties& mProps_B,
                  const Transform&      X_PJb,
                  const Transform&      X_BJ,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerPosition(const SBModelVars& mv, const Transform& X_JbJ,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv)) {
            //TODO
        } else {
            toQuat(q) = X_JbJ.R().convertToQuaternion().asVec4();
        }
    }
    void setMobilizerVelocity(const SBModelVars&, const SpatialVec& V_JbJ,
                              Vector& u) const
    {
            toU(u) = V_JbJ[0]; // relative angular velocity always used as generalized speeds
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q); // angular coordinates
            toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            toQuat(qnorm) = quat / quat.norm();
        }
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_JbJ) const
    {
        X_JbJ.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(mv))
            X_JbJ.updR().setToBodyFixed123(fromQ(q));
        else {
            // TODO: should use qnorm pool
            X_JbJ.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Vec3 T_JB_G = -X_GB.R()*X_BJ.T(); // vec from OJ to OB, expr. in G

        // The rotational coordinates are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const Rotation R_GJb = X_GP.R() * X_PJb.R();
        H[0] = SpatialRow(~R_GJb.x(), ~(R_GJb.x() % T_JB_G));
        H[1] = SpatialRow(~R_GJb.y(), ~(R_GJb.y() % T_JB_G));
        H[2] = SpatialRow(~R_GJb.z(), ~(R_GJb.z() % T_JB_G));
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const Vec3& w_JbJ = fromU(u); // angular velocity of J in Jb 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_JbJ = getX_JbJ(cc).R();
            toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(q),
                                        ~R_JbJ*w_JbJ); // need w in *body*, not parent
        } else
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(q),w_JbJ);
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3& w_JbJ     = fromU(u); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_JbJ_dot = fromU(udot);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_JbJ = getX_JbJ(cc).R();
            toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                                  (fromQ(q), ~R_JbJ*w_JbJ, ~R_JbJ*w_JbJ_dot);
        } else
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(q),w_JbJ,w_JbJ_dot);
    }

    void setQ(
        const SBModelVars& mv, 
        const Vector&      qIn, 
        Vector&            q) const 
    {
        if (getUseEulerAngles(mv))
            toQ(q) = fromQ(qIn);
        else
            toQuat(q) = fromQuat(qIn);
    }

    int getMaxNQ()              const {return 4;}
    int getNQ(const SBModelVars& mv) const {
        return getUseEulerAngles(mv) ? 3 : 4;
    } 

    void setDefaultPositionValues(
        const SBModelVars& mv,
        Vector&            q) const 
    {
        if (getUseEulerAngles(mv)) {
            //TODO: kludge
            toQuat(q) = Vec4(0); // clear unused element
            toQ(q) = 0.;
        }
        else toQuat(q) = Vec4(1.,0.,0.,0.);
    }

    bool enforceQuaternionConstraints(
        const SBModelVars& mv,
        Vector&            q) const 
    {
        if (getUseEulerAngles(mv)) 
            return false;   // no change
        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();
        return true;
    }

    void getInternalForce(const SBAccelerationCache&, Vector&) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
        /* OLD BALL CODE:
        Vector& f = s.cache->netHingeForces;
        //dependency: calcR_PB must be called first
        assert( useEuler );

        const Vec<3,Vec2>& scq = getSinCosQ(s);
        const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
        const Real sTheta = scq[1][0], cTheta = scq[1][1];
        const Real sPsi   = scq[2][0], cPsi   = scq[2][1];

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[uIndex]) = eTorque;
        */
    }
};

/// Free joint. This provides six degrees of freedom, three rotational and
/// three translational. The rotation is like the ball joint above; the
/// translation is like the Cartesian joint above.
class RBNodeTranslateRotate3 : public RigidBodyNodeSpec<6> {
public:
    virtual const char* type() { return "full"; }

    RBNodeTranslateRotate3(const MassProperties& mProps_B,
                           const Transform&      X_PJb,
                           const Transform&      X_BJ,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,X_PJb,X_BJ,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerPosition(const SBModelVars& mv, const Transform& X_JbJ,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv)) {
            //TODO orientation
            toQVec3(q,3) = X_JbJ.T(); // translation
        } else {
            toQuat(q) = X_JbJ.R().convertToQuaternion().asVec4();
            toQVec3(q,4) = X_JbJ.T();
        }
    }
    void setMobilizerVelocity(const SBModelVars&, const SpatialVec& V_JbJ,
                              Vector& u) const
    {
        toUVec3(u,0) = V_JbJ[0]; // relative angular velocity always used as generalized speeds
        toUVec3(u,3) = V_JbJ[1];
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars& mv, 
        const Vector&      q, 
        Vector&            sine, 
        Vector&            cosine, 
        Vector&            qnorm) const
    {
        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            toQuat(qnorm) = quat / quat.norm();
        }
    }

    // Calculate X_JbJ.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform& X_JbJ) const 
    {
        if (getUseEulerAngles(mv)) {
            X_JbJ.updR().setToBodyFixed123(fromQVec3(q,0));
            X_JbJ.updT() = X_JbJ.R()*fromQVec3(q,3);
        } else {
            X_JbJ.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
            X_JbJ.updT() = X_JbJ.R()*fromQVec3(q,4);  // because translation is in J frame
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BJ  = getX_BJ();  // fixed
        const Transform& X_PJb = getX_PJb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Transform X_GJb = X_GP * X_PJb;
        const Transform X_GJ  = X_GB * X_BJ;

        const Vec3 T_JbB_G = X_GB.T() - X_GJb.T();

        // The rotational speeds are defined in the space-fixed 
        // (that is, Jb) frame, so the orientation of Jb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        // *Then* we translate along the new J axes.

        H[0] = SpatialRow(~X_GJb.x(), ~(X_GJb.x() % T_JbB_G));
        H[1] = SpatialRow(~X_GJb.y(), ~(X_GJb.y() % T_JbB_G));
        H[2] = SpatialRow(~X_GJb.z(), ~(X_GJb.z() % T_JbB_G));
        H[3] = SpatialRow(  Row3(0) ,     ~X_GJ.x());
        H[4] = SpatialRow(  Row3(0) ,     ~X_GJ.y());
        H[5] = SpatialRow(  Row3(0) ,     ~X_GJ.z());
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        Vector&                qdot) const
    {
        const Vec3& w_JbJ = fromUVec3(u,0); // Angular velocity
        const Vec3& v_JbJ = fromUVec3(u,3); // Linear velocity
        if (getUseEulerAngles(mv)) {
            const Rotation& R_JbJ = getX_JbJ(cc).R();
            const Vec3& theta = fromQVec3(q,0); // Euler angles
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
                                            ~R_JbJ*w_JbJ); // need w in *body*, not parent
            toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdot,3) = v_JbJ;
        } else {
            const Vec4& quat = fromQuat(q);
            toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,w_JbJ);
            toQVec3(qdot,4) = v_JbJ;
        }
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3& w_JbJ     = fromUVec3(u,0); // angular velocity of J in Jb
        const Vec3& v_JbJ     = fromUVec3(u,3); // linear velocity
        const Vec3& w_JbJ_dot = fromUVec3(udot,0);
        const Vec3& v_JbJ_dot = fromUVec3(udot,3);
        if (getUseEulerAngles(mv)) {
            const Rotation& R_JbJ = getX_JbJ(cc).R();
            const Vec3& theta  = fromQVec3(q,0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, ~R_JbJ*w_JbJ, ~R_JbJ*w_JbJ_dot);
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_JbJ_dot;
            //cout << "   w_JbJ_dot=" << w_JbJ_dot << "  v_JbJ_dot=" << v_JbJ_dot << endl;
            //cout << "   qdotdot=" << fromQVec3(qdotdot,0) << endl;
        } else {
            const Vec4& quat  = fromQuat(q);
            toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                             (quat,w_JbJ,w_JbJ_dot);
            toQVec3(qdotdot,4) = v_JbJ_dot;
        }
    }

    void setQ(Vector& q, const SBModelVars& mv, const Vector& qIn) const {
        if (getUseEulerAngles(mv))
            toQ(q) = fromQ(qIn);
        else {
            toQuat(q)    = fromQuat(qIn);
            toQVec3(q,4) = fromQVec3(qIn,4);
        }
    }

    int getMaxNQ()                   const {return 7;}
    int getNQ(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 

    void setDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
    {
        if (getUseEulerAngles(mv)) {
            toQVec3(q,4) = Vec3(0); // TODO: kludge, clear unused element
            toQ(q) = 0.;
        } else {
            toQuat(q) = Vec4(1.,0.,0.,0.);
            toQVec3(q,4) = 0.;
        }
    }

    bool enforceQuaternionConstraints(const SBModelVars& mv, Vector& q) const {
        if (getUseEulerAngles(mv)) 
            return false; // no change
        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();
        return true;
    }

    void getInternalForce(const SBAccelerationCache&, Vector&) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
        /* OLD BALL CODE:
        Vector& f = s.cache->netHingeForces;
        //dependency: calcR_PB must be called first
        assert( useEuler );

        const Vec<3,Vec2>& scq = getSinCosQ(s);
        const Real sPhi   = scq[0][0], cPhi   = scq[0][1];
        const Real sTheta = scq[1][0], cTheta = scq[1][1];
        const Real sPsi   = scq[2][0], cPsi   = scq[2][1];

        Vec3 torque = forceInternal;
        const Mat33 M( 0.          , 0.          , 1.    ,
                      -sPhi        , cPhi        , 0.    ,
                       cPhi*cTheta , sPhi*cTheta ,-sTheta );
        Vec3 eTorque = RigidBodyNode::DEG2RAD * M * torque;

        Vec3::updAs(&v[uIndex]) = eTorque;
        */
    }
};

////////////////////////////////////////////////////
// RigidBodyNode factory based on mobilizer type. //
////////////////////////////////////////////////////

/*static*/ RigidBodyNode*
RigidBodyNode::create(
    const MassProperties&    m,            // mass properties in body frame
    const Transform&         X_PJb,        // parent's attachment frame for this joint
    const Transform&         X_BJ,         // inboard joint frame J in body frame
    Mobilizer::MobilizerType type,
    bool                     isReversed,   // child-to-parent orientation?
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot)  
{
    assert(!isReversed);

    switch(type) {
    case Mobilizer::ThisIsGround:
        return new RBGroundBody();
    case Mobilizer::Torsion:
        return new RBNodeTorsion(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Universal:        
        return new RBNodeRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Orientation:
        return new RBNodeRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Cartesian:
        return new RBNodeTranslate(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::FreeLine:
        return new RBNodeTranslateRotate2(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Free:
        return new RBNodeTranslateRotate3(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Sliding:
        return new RBNodeSlider(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Cylinder:
        return new RBNodeCylinder(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::BendStretch:
        return new RBNodeBendStretch(m,X_PJb,X_BJ,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Planar:
    case Mobilizer::Gimbal:
    case Mobilizer::Weld:

    default: 
        assert(false);
    };

    return 0;
}

/////////////////////////////////////////////////////////////
// Implementation of RigidBodyNodeSpec base class methods. //
/////////////////////////////////////////////////////////////


//
// to be called from base to tip.
//
template<int dof> void
RigidBodyNodeSpec<dof>::setVelFromSVel(
    const SBPositionCache& cc, 
    const SBVelocityCache& mc,
    const SpatialVec&      sVel, 
    Vector&                u) const 
{
    toU(u) = getH(cc) * (sVel - (~getPhi(cc) * parent->getV_GB(mc)));
}

//
// Given only position-related quantities from the State 
//      Mk  (this body's spatial inertia matrix)
//      Phi (composite body child-to-parent shift matrix)
//      H   (joint transition matrix)
// we calculate dynamic quantities 
//      P   (articulated body inertia)
//      D   (factored mass matrix LDL' diagonal part D=H*P*~H)
//      DI  (inverse of D)
//      G   (P * ~H * DI)
//   tauBar (I-G*H, a temporary not reused elsewhere)
//      Psi (Phi*(I-G*H), articulated body child-to-parent shift matrix)
// and put them in the state cache.
// This must be called tip-to-base (inward).
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcArticulatedBodyInertiasInward(
    const SBPositionCache& cc,
    SBDynamicsCache&       dc) const 
{
    updP(dc) = getMk(cc);
    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialMat& tauBarChild = children[i]->getTauBar(dc);
        const SpatialMat& PChild      = children[i]->getP(dc);
        const PhiMatrix&  phiChild    = children[i]->getPhi(cc);

        // TODO: this is around 450 flops but could be cut in half by
        // exploiting symmetry.
        updP(dc) += phiChild * (tauBarChild * PChild) * ~phiChild;
    }

    const Mat<2,dof,Vec3> PHt = getP(dc) * ~getH(cc);
    updD(dc)  = getH(cc) * PHt;
    // this will throw an exception if the matrix is ill conditioned
    updDI(dc) = getD(dc).invert();
    updG(dc)  = PHt * getDI(dc);

    // TODO: change sign on tau to make it GH-I instead, which only requires
    // subtractions on the diagonal rather than negating all the off-diag stuff.
    // That would save 30 flops here (I know, not much).
    updTauBar(dc)  = 1.; // identity matrix
    updTauBar(dc) -= getG(dc) * getH(cc);
    updPsi(dc)     = getPhi(cc) * getTauBar(dc);
}



// To be called base to tip.
// sherm 060723: As best I can tell this is calculating the inverse of
// the "operational space inertia" at the body frame origin for each body.
// See Equation 20 in Rodriguez,Jain, & Kreutz-Delgado: A spatial operator algebra 
// for manipulator modeling and control. Intl. J. Robotics Research 
// 10(4):371-381 (1991).
template<int dof> void
RigidBodyNodeSpec<dof>::calcYOutward(
    const SBPositionCache& cc,
    SBDynamicsCache&       dc) const 
{
    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?
    updY(dc) = (~getH(cc) * getDI(dc) * getH(cc)) 
                + (~getPsi(dc) * parent->getY(dc) * getPsi(dc));
}

//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(
    const SBPositionCache& cc,
    const SBDynamicsCache& dc,
    const SpatialVec&      spatialForce,
    SBAccelerationCache&   rc) const 
{
    SpatialVec& z = updZ(rc);
    z = getCentrifugalForces(dc) - spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(rc);
        const PhiMatrix&  phiChild  = children[i]->getPhi(cc);
        const SpatialVec& GepsChild = children[i]->getGepsilon(rc);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(rc)  = getAppliedJointForce(dc) - getH(cc)*z; // TODO: pass in hinge forces
    updNu(rc)       = getDI(dc) * getEpsilon(rc);
    updGepsilon(rc) = getG(dc)  * getEpsilon(rc);
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to calcZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel(
    const SBModelVars&     mv,
    const Vector&          allQ,
    const SBPositionCache& cc,
    const Vector&          allU,
    const SBDynamicsCache& dc,
    SBAccelerationCache&   rc,
    Vector&                allUdot,
    Vector&                allQdotdot) const 
{
    Vec<dof>&        udot   = toU(allUdot);
    const SpatialVec alphap = ~getPhi(cc) * parent->getA_GB(rc); // ground A_GB is 0

    udot        = getNu(rc) - (~getG(dc)*alphap);
    updA_GB(rc) = alphap + ~getH(cc)*udot + getCoriolisAcceleration(dc);  

    calcQDotDot(mv, allQ, cc, allU, allUdot, allQdotdot);  
}

 
//
// To be called from tip to base.
// Temps do not need to be initialized.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcUDotPass1Inward(
    const SBPositionCache&      cc,
    const SBDynamicsCache&      dc,
    const Vector&               jointForces,
    const Vector_<SpatialVec>&  bodyForces,
    Vector_<SpatialVec>&        allZ,
    Vector_<SpatialVec>&        allGepsilon,
    Vector&                     allEpsilon) const 
{
    const Vec<dof>&   myJointForce = fromU(jointForces);
    const SpatialVec& myBodyForce  = fromB(bodyForces);
    SpatialVec&       z            = toB(allZ);
    SpatialVec&       Geps         = toB(allGepsilon);
    Vec<dof>&         eps          = toU(allEpsilon);

    z = getCentrifugalForces(dc) - myBodyForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(cc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - getH(cc)*z;
    Geps = getG(dc)  * eps;
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcUDotPass2Outward(
    const SBPositionCache& cc,
    const SBDynamicsCache& dc,
    const Vector&          allEpsilon,
    Vector_<SpatialVec>&   allA_GB,
    Vector&                allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = parent->getNodeNum()== 0 
        ? SpatialVec(Vec3(0), Vec3(0))
        : ~getPhi(cc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + ~getH(cc)*udot + getCoriolisAcceleration(dc);  
}

 
//
// To be called from tip to base.
// Temps do not need to be initialized.
//
// This calculates udot = M^-1 f in two O(N) passes. Note that
// we are ignoring velocities; if there are any velocity-dependent
// forces they should already be in f.
//
// (sherm 060727) TODO: surely this can be tightened up?
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcMInverseFPass1Inward(
    const SBPositionCache& cc,
    const SBDynamicsCache& dc,
    const Vector&          f,
    Vector_<SpatialVec>&   allZ,
    Vector_<SpatialVec>&   allGepsilon,
    Vector&                allEpsilon) const 
{
    const Vec<dof>&   myJointForce = fromU(f);
    SpatialVec&       z            = toB(allZ);
    SpatialVec&       Geps         = toB(allGepsilon);
    Vec<dof>&         eps          = toU(allEpsilon);

    z = SpatialVec(Vec3(0), Vec3(0));

    for (int i=0 ; i<(int)children.size() ; i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(cc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - getH(cc)*z;
    Geps = getG(dc)  * eps;
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcMInverseFPass2Outward(
    const SBPositionCache& cc,
    const SBDynamicsCache& dc,
    const Vector&          allEpsilon,
    Vector_<SpatialVec>&   allA_GB,
    Vector&                allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = parent->getNodeNum()== 0 
        ? SpatialVec(Vec3(0), Vec3(0))
        : ~getPhi(cc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + ~getH(cc)*udot;  
}

//
// Calculate product of partial velocities J and a gradient vector on each of the
// outboard bodies. Requires that Phi and H are available, so this
// should only be called in Stage::Position or higher. This does not change the cache at all.
// NOTE (sherm 060214): I reworked this from the original. This one no longer incorporates
// applied hinge gradients if there are any; just add those in at the end if you want them.
//
// (sherm 060727) In spatial operators, this calculates H*Phi*F where F are the spatial forces
// applied to each body. See Schwieters Eq. 41.
//
// Call tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcInternalGradientFromSpatial(
    const SBPositionCache&      cc,
    Vector_<SpatialVec>&        zTmp,
    const Vector_<SpatialVec>&  X, 
    Vector&                     JX) const
{
    const SpatialVec& in  = X[getNodeNum()];
    Vec<dof>&         out = Vec<dof>::updAs(&JX[getUIndex()]);
    SpatialVec&       z   = zTmp[getNodeNum()];

    z = in;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild = children[i]->getPhi(cc);

        z += phiChild * zChild;
    }

    out = getH(cc) * z; 
}

//
// To be called from tip to base.
// Temps do not need to be initialized.
// (sherm 060727) In spatial operators, this calculates H*Phi*(F-(Pa+b))
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcEquivalentJointForces(
    const SBPositionCache&     cc,
    const SBDynamicsCache&     dc,
    const Vector_<SpatialVec>& bodyForces,
    Vector_<SpatialVec>&       allZ,
    Vector&                    jointForces) const 
{
    const SpatialVec& myBodyForce  = fromB(bodyForces);
    SpatialVec&       z            = toB(allZ);
    Vec<dof>&         eps          = toU(jointForces);

    // Centrifugal forces are PA+b where P is articulated body inertia,
    // A is total coriolis acceleration, and b is gyroscopic force.
    z = myBodyForce - getTotalCentrifugalForces(dc);

    for (int i=0 ; i<(int)children.size() ; i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(cc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];

        z += phiChild * zChild; 
    }

    eps  = getH(cc) * z;
}

