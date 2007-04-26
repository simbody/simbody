/* Portions copyright (c) 2005-6 Stanford University and Michael Sherman.
 * Contributors: Derived from public domain IVM code written by
 *               Charles Schwieters at NIH.
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
 * its mobilizer (inboard joint), that is, one node in the multibody tree.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

#include "SimbodyMatterSubsystemRep.h"
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

    // The Phi matrix conveniently performs child-to-parent (inward) shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // (outward) shifting for velocities.
    updPhi(cc) = PhiMatrix(T_PB_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    updInertia_OB_G(cc) = getInertia_OB_B().reexpress(~getX_GB(cc).R());
    updCB_G(cc)         = getX_GB(cc).R()*getCOM_B();

    updCOM_G(cc) = getX_GB(cc).T() + getCB_G(cc);

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    // Note: we need to calculate this now so that we'll be able to calculate
    // kinetic energy without going past the Velocity stage.
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
    // because they can only differ along H *which is constant* between P & B
    // (I don't think that is always true for us.)
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

    // (sherm 060912) Caution: using T_MbB_G here is assuming that all the mobilizers do their
    // rotations in the Mb frame *followed* by translations in the M frame.
    // I *think* you could switch the mobilizer to work the other way and change to T_MB_G
    // in the mobilizers and here, but I haven't tried it. TODO: this should be written
    // in terms of the H matrix somehow so that this computation and the joint definition
    // do not have to be synchronized this way.

    // TODO: this is all very fishy and needs to be rederived rigorously.
    // (sherm 070425) My current wild theory is that this code is trying
    // to compute the Hdot*u term of the coriolis acceleration "a"
    // that arises from differentiation of
    // the velocity term H*u, but is leaving out a term of the derivative
    // that arises when H is not constant in the parent (Mb) frame.

    const Vec3 T_PB_G = getX_GB(cc).T() - getX_GP(cc).T();
    const Vec3 T_MbB_G = T_PB_G - getX_GP(cc).R()*getX_PMb().T();
    const Vec3 w_PB_G = getV_PB_G(mc)[0];
    const Vec3 v_PB_G = getV_PB_G(mc)[1]; // includes w_PB_G % T_MbB_G term

    const SpatialVec wwr = SpatialVec(Vec3(0), 
                  pOmega % (pOmega % T_PB_G) - w_PB_G % (w_PB_G % T_MbB_G));
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
    /*virtual*/int  getDOF()   const {return 0;}
    /*virtual*/int  getMaxNQ() const {return 0;}
    /*virtual*/int  getNQ(const SBModelVars&) const {return 0;}
    /*virtual*/bool isUsingQuaternion(const SBModelVars&) const {return false;}


    /*virtual*/void calcZ(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const SpatialVec& spatialForce,
        SBAccelerationCache&               ) const {} 

    /*virtual*/void calcYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/void calcAccel(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBDynamicsCache& dc,
        SBAccelerationCache&   ac,
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
        const SBModelCache& mc,
        const Vector&       q,
        Vector&             qErr,
        SBPositionCache&    pc) const {}

    /*virtual*/void realizeVelocity(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        SBVelocityCache&       vc,
        Vector&                qdot) const {}

    /*virtual*/ void setMobilizerTransform
       (const SBModelVars&, const Transform& X_MbM, Vector& q) const {}
    /*virtual*/ void setMobilizerRotation
       (const SBModelVars&, const Rotation& R_MbM, Vector& q) const {}
    /*virtual*/ void setMobilizerTranslation
       (const SBModelVars&, const Vec3& T_MbM, Vector& q,
        bool dontChangeOrientation) const {}

    /*virtual*/ void setMobilizerVelocity
       (const SBModelVars&, const Vector& q, const SpatialVec& V_MbM, Vector& u) const {}
    /*virtual*/ void setMobilizerAngularVelocity
       (const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {}
    /*virtual*/ void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector& q, const Vec3& v_MbM, Vector& u,
        bool dontChangeAngularVelocity) const {}


    /*virtual*/void setVelFromSVel(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
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
                      const Transform&      X_PMb,
                      const Transform&      X_BM,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNode(mProps_B, X_PMb, X_BM)
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


    // Provide default implementations for setMobilizerTransform() and setMobilizerVelocity() 
    // which are implemented using the rotational and translational quantity routines. These assume
    // that the rotational and translational coordinates are independent, with rotation handled
    // first and then left alone. If a mobilizer type needs to deal with rotation and
    // translation simultaneously, it should provide a specific implementation for these two routines.
    // *Each* mobilizer must implement setMobilizer{Rotation,Translation,AngularVelocity,LinearVelocity};
    // there are no defaults.

    virtual void setMobilizerTransform(const SBModelVars& mv, const Transform& X_MbM, Vector& q) const {
        setMobilizerRotation   (mv,X_MbM.R(),q);
        setMobilizerTranslation(mv,X_MbM.T(),q,true); // don't fiddle with the rotation
    }

    virtual void setMobilizerVelocity(const SBModelVars& mv, const Vector& q, const SpatialVec& V_MbM, Vector& u) const {
        setMobilizerAngularVelocity(mv,q,V_MbM[0],u);
        setMobilizerLinearVelocity (mv,q,V_MbM[1],u,true); // don't fiddle with the angular velocity
    }

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
    /// for mobilizers which have angular coordinates, and qErr and qnorm are only for
    /// mobilizers using quaternions. Other mobilizers can provide a null routine.
    /// Each of the passed-in Vectors is a "q-like" object, that is, allocated
    /// to the bodies in a manner parallel to the q state variable, except that qErr
    /// has just one slot per quaternion and must be accessed using the node's 
    /// quaternionIndex which is in the Model cache.
    virtual void calcJointSinCosQNorm(
        const SBModelVars&  mv, 
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const=0;

    /// This mandatory routine calculates the across-joint transform X_MbM generated
    /// by the current q values. This may depend on sines & cosines or normalized
    /// quaternions already being available in the State cache.
    virtual void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const=0;

    /// This routine is NOT joint specific, but cannot be called until the across-joint
    /// transform X_MbM has been calculated and is available in the State cache.
    void calcBodyTransforms(
        const SBPositionCache& cc, 
        Transform&             X_PB, 
        Transform&             X_GB) const 
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_MbM = getX_MbM(cc); // just calculated
        const Transform& X_GP  = getX_GP(cc);  // already calculated

        X_PB = X_PMb * X_MbM * ~X_BM; // TODO: precalculate X_MB
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
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q,
        Vector&             qErr,
        SBPositionCache&    pc) const 
    {
        calcJointSinCosQNorm(mv, mc, q, pc.sq, pc.cq, qErr, pc.qnorm);

        calcAcrossJointTransform (mv, q, updX_MbM(pc));
        calcBodyTransforms       (pc, updX_PB(pc), updX_GB(pc));
        calcJointTransitionMatrix(pc, updH(pc));
        calcJointIndependentKinematicsPos(pc);
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

    // copyQ and copyU extract this node's values from the supplied
    // q-sized or u-sized array and put them in the corresponding
    // locations in the output variable. Joints which need quaternions should
    // override copyQ to copy the extra q.
    virtual void copyQ(
        const SBModelVars& mv, 
        const Vector&      qIn, 
        Vector&            q) const
    {
        toQ(q) = fromQ(qIn);
    }

    virtual void copyU(
        const SBModelVars& mv, 
        const Vector&      uIn, 
        Vector&            u) const
    {
        toU(u) = fromU(uIn);
    }

    int          getDOF()            const {return dof;}
    virtual int  getMaxNQ()          const {return dof;} // maxNQ can be larger than dof
    virtual int  getNQ(const SBModelVars&) const { return dof; } // DOF <= NQ <= maxNQ
    virtual bool isUsingQuaternion(const SBModelVars&) const {return false;}

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
    const Vec<dof>&   getPrescribedUdot   (const SBDynamicsCache& dc) const 
        {return fromU(dc.prescribedUdot);}

    // Special case state access for 1-dof joints
    const Real& get1AppliedJointForce(const SBDynamicsCache& dc) const {return from1U(dc.appliedMobilityForces);}
    const Real& get1PrescribedUdot   (const SBDynamicsCache& dc) const {return from1U(dc.prescribedUdot);}

    // Cache entries (cache is mutable in a const State)

        // Position

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

        // Velocity

        // Dynamics
    const Mat<dof,dof>& getD(const SBDynamicsCache& dc) const {return fromUSq(dc.storageForD);}
    Mat<dof,dof>&       updD(SBDynamicsCache&       dc) const {return toUSq  (dc.storageForD);}

    const Mat<dof,dof>& getDI(const SBDynamicsCache& dc) const {return fromUSq(dc.storageForDI);}
    Mat<dof,dof>&       updDI(SBDynamicsCache&       dc) const {return toUSq  (dc.storageForDI);}

    const Mat<2,dof,Vec3>& getG(const SBDynamicsCache& dc) const
      { return Mat<2,dof,Vec3>::getAs(&dc.storageForG(0,uIndex)); }
    Mat<2,dof,Vec3>&       updG(SBDynamicsCache&       dc) const
      { return Mat<2,dof,Vec3>::updAs(&dc.storageForG(0,uIndex)); }

        // Acceleration

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
 * the parent body's Mb frame, with M=Mb when all 3 coords are 0,
 * and the orientation of M in Mb is 0 (identity) forever.
 */
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const MassProperties& mProps_B,
                    const Transform&      X_PMb,
                    const Transform&      X_BM,
                    int&                  nextUSlot,
                    int&                  nextUSqSlot,
                    int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // the only rotation this mobilizer can represent is identity
    }
    void setMobilizerTranslation(const SBModelVars&, const Vec3&  T_MbM, Vector& q, bool only) const {
        // here's what this joint is really good at!
        toQ(q) = T_MbM;
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // The only angular velocity this can represent is zero.
    }
    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // linear velocity is in a Cartesian joint's sweet spot
        toU(u) = v_MbM;
    }

    // This is required but does nothing here since there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const { }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        // Translation vector q is expressed in Mb (and M since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_MbM = Transform(Rotation(), fromQ(q));
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const 
    {
        const Transform& X_PMb = getX_PMb();      // fixed config of Mb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP = getX_GP(cc); // parent orientation in ground

        // Note that H is spatial. The current spatial directions for our qs are
        // the axes of the Mb frame expressed in Ground.
        const Rotation R_GMb = X_GP.R()*X_PMb.R();
        H[0] = SpatialRow( Row3(0), ~R_GMb.x() );
        H[1] = SpatialRow( Row3(0), ~R_GMb.y() );
        H[2] = SpatialRow( Row3(0), ~R_GMb.z() );
    }
};



/**
 * Sliding joint (1 dof translation). The translation is along the x
 * axis of the parent body's Mb frame, with M=Mb when the coordinate
 * is zero and the orientation of M in Mb frozen at 0 forever.
 */
class RBNodeSlider : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "slider"; }

    RBNodeSlider(const MassProperties& mProps_B,
                 const Transform&      X_PMb,
                 const Transform&      X_BM,
                 int&                  nextUSlot,
                 int&                  nextUSqSlot,
                 int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }
        // Implementations of virtual methods.

    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation a slider can represent is identity.
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // We can only represent the x coordinate with this joint.
        to1Q(q) = T_MbM[0];
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // The only angular velocity a slider can represent is zero.
    }

    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // We can only represent a velocity along x with this joint.
        to1U(u) = v_MbM[0];
    }

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const { }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        // Translation vector q is expressed in Mb (and M since they have same orientation).
        // A sliding joint can't change orientation, and only translates along x. 
        X_MbM = Transform(Rotation(), Vec3(from1Q(q),0,0));
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_PMb   = getX_PMb();      // fixed config of Mb in P

        // Calculated already since we're going base to tip.
        const Transform& X_GP    = getX_GP(cc); // parent configuration in ground

        // Note that H is spatial. The current spatial directions for our q is
        // the x axis of the Mb frame expressed in Ground.
        const Vec3 x_GMb = X_GP.R()*X_PMb.x();
        H[0] = SpatialRow( Row3(0), ~x_GMb );
    }
};

/**
 * This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
 * about a particular axis, the z axis of the parent's Mb frame, which is 
 * aligned forever with the z axis of the body's M frame. In addition, the
 * origin points of M and Mb are identical forever.
 */
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const MassProperties& mProps_B,
                  const Transform&      X_PMb,
                  const Transform&      X_BM,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our pin joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        to1Q(q) = angles123[2];
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_MbM[2]; // project angular velocity onto z axis
    }

    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Real& angle = from1Q(q); // angular coordinate
        to1Q(sine)    = std::sin(angle);
        to1Q(cosine)  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        const Real& theta  = from1Q(q);    // angular coordinate

        // We're only updating the orientation here because a torsion joint
        // can't translate (it is defined as a rotation about the z axis).
        X_MbM.updR().setToRotationAboutZ(theta);
        X_MbM.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        // Vec from OM to OB, expr. in G (note minus sign). This is the
        // amount by which the B origin is off the Mz axis.
        const Vec3 T_MB_G = -X_GB.R()*X_BM.T(); 

        // Calc H matrix in space-fixed coords.
        // This works because the joint z axis is the same in M & Mb
        // since that's what we rotate around.
        const Vec3 z_G = X_GP.R() * X_PMb.z();
        H[0] = SpatialRow( ~z_G, ~(z_G % T_MB_G) );
    }
};

/**
 * This is a "cylinder" joint, meaning one degree of rotational freedom
 * about a particular axis, and one degree of translational freedom
 * along the same axis. For molecules you can think of this as a combination
 * of torsion and bond stretch. The axis used is the z axis of the parent's
 * Mb frame, which is aligned forever with the z axis of the body's M frame.
 * In addition, the origin points of M and Mb are separated only along the
 * z axis; i.e., they have the same x & y coords in the Mb frame. The two
 * generalized coordinates are the rotation and the translation, in that order.
 */
class RBNodeCylinder : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "cylinder"; }

    RBNodeCylinder(const MassProperties& mProps_B,
                   const Transform&      X_PMb,
                   const Transform&      X_BM,
                   int&                  nextUSlot,
                   int&                  nextUSqSlot,
                   int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }


    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our cylinder joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        toQ(q)[0] = angles123[2];
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // Because the M and Mb origins must lie along their shared z axis, there is no way to
        // create a translation by rotating around z. So the only translation we can represent
        // is that component which is along z.
        toQ(q)[1] = T_MbM[2];
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_MbM[2];
    }

    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // Because the M and Mb origins must lie along their shared z axis, there is no way to
        // create a linear velocity by rotating around z. So the only linear velocity we can represent
        // is that component which is along z.
        toU(u)[1] = v_MbM[2];
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Real& angle = fromQ(q)[0];
        toQ(sine)[0]    = std::sin(angle);
        toQ(cosine)[0]  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        const Vec2& coords  = fromQ(q);

        X_MbM.updR().setToRotationAboutZ(coords[0]);
        X_MbM.updT() = Vec3(0,0,coords[1]);
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        // Vec from OM to OB, expr. in G (note minus sign). This is the
        // amount by which the B origin is off the Mz axis.
        const Vec3 T_MB_G = -X_GB.R()*X_BM.T();

        // Calc H matrix in space-fixed coords.
        // This works because the joint z axis is the same in M & Mb
        // since that's what we rotate around.
        const Vec3 z_GMb = X_GP.R()*X_PMb.z();
        H[0] = SpatialRow( ~z_GMb, ~(z_GMb % T_MB_G) );
        H[1] = SpatialRow( Row3(0), ~z_GMb );
    }
};


/**
 * This is a "bend-stretch" joint, meaning one degree of rotational freedom
 * about a particular axis, and one degree of translational freedom
 * along a perpendicular axis. The z axis of the parent's Mb frame is 
 * used for rotation (and that is always aligned with the M frame z axis).
 * The x axis of the *M* frame is used for translation; that is, first
 * we rotate around z, which moves M's x with respect to Mb's x. Then
 * we slide along the rotated x axis. The two
 * generalized coordinates are the rotation and the translation, in that order.
 */
class RBNodeBendStretch : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "bendstretch"; }

    RBNodeBendStretch(const MassProperties& mProps_B,
                      const Transform&      X_PMb,
                      const Transform&      X_BM,
                      int&                  nextUSlot,
                      int&                  nextUSqSlot,
                      int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }


    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our bend-stretch joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        toQ(q)[0] = angles123[2];
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // We can represent any translation that puts the M origin in the x-y plane of Mb,
        // by a suitable rotation around z followed by translation along x.
        const Vec2 r = T_MbM.getSubVec<2>(0); // (rx, ry)

        // If we're not allowed to change rotation then we can only move along Mx.
        if (only) {
            const Real angle = fromQ(q)[0];
            const Vec2 Mx(std::cos(angle), std::sin(angle)); // a unit vector
            toQ(q)[1] = dot(r,Mx);
            return;
        }

        const Real d = r.norm();

        // If there is no translation worth mentioning, we'll leave the rotational
        // coordinate alone, otherwise rotate so M's x axis is aligned with r.
        if (d >= 4*NTraits<Real>::Eps) {
            const Real angle = std::atan2(r[1],r[0]);
            toQ(q)[0] = angle;
            toQ(q)[1] = d;
        } else
            toQ(q)[1] = 0;
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_MbM[2];
    }

    // If the translational coordinate is zero, we can only represent a linear velocity 
    // of OM in Mb which is along M's current x axis direction. Otherwise, we can 
    // represent any velocity in the x-y plane by introducing angular velocity about z.
    // We can never represent a linear velocity along z.
    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector& q, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // Decompose the requested v into "along Mx" and "along My" components.
        const Rotation R_MbM = Rotation::aboutZ(fromQ(q)[0]); // =[ Mx My Mz ] in Mb
        const Vec3 v_MbM_M = ~R_MbM*v_MbM; // re-express in M frame

        toU(u)[1] = v_MbM_M[0]; // velocity along Mx we can represent directly

        if (only) {
            // We can't do anything about My velocity if we're not allowed to change
            // angular velocity, so we're done.
            return;
        }

        const Real x = fromQ(q)[1]; // translation along Mx (signed)
        if (std::abs(x) < NTraits<Real>::Eps_78) {
            // No translation worth mentioning; we can only do x velocity, which we just set above.
            return;
        }

        // significant translation
        toU(u)[0] = v_MbM_M[1] / x; // set angular velocity about z to produce vy
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Real& angle = fromQ(q)[0];
        toQ(sine)[0]    = std::sin(angle);
        toQ(cosine)[0]  = std::cos(angle);
        // no quaternions
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        const Vec2& coords  = fromQ(q);    // angular coordinate

        X_MbM.updR().setToRotationAboutZ(coords[0]);
        X_MbM.updT() = X_MbM.R()*Vec3(coords[1],0,0); // because translation is in M frame
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        const Transform X_GMb = X_GP * X_PMb;
        const Transform X_GM  = X_GB * X_BM;

        const Vec3 T_MbB_G = X_GB.T() - X_GMb.T();
        //const Vec3 T_MB_G = -X_GB.R()*X_BM.T();

        // Rotate around Mb's z axis, *then* translate along M's new x axis.
        H[0] = SpatialRow( ~X_GMb.z(), ~(X_GMb.z() % T_MbB_G) );
        H[1] = SpatialRow( Row3(0), ~X_GM.x());
    }
};


/**
 * Somewhat odd U-joint-like mobilizer which allows rotation about the two axes
 * perpendicular to the outboard body's Mz axis (that is, about Mx and My), but
 * never allows rotation about Mz. This is appropriate for diatoms and for allowing 
 * torsion+bond angle bending. The generalized coordinates are a 1-2
 * *space* (Mb) fixed Euler angle sequence (same as 1-2-3 with the 3rd rotation zero).
 * The generalized speeds are the x and y coordinates of w_MbM_M, that is, 
 * the angular velocity of M in Mb, expressed in M (the Mz component of 
 * the angular velocity is always zero).
 * The qdots are just u's because (apparently!) the angular velocity x and y
 * components in M are also the derivatives of the Mb-fixed rotation angles.
 * TODO: WHY DOES THIS WORK???
 */
class RBNodeRotate2 : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "rotate2"; }

    RBNodeRotate2(const MassProperties& mProps_B,
                  const Transform&      X_PMb,
                  const Transform&      X_BM,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotations this joint can handle are about Mx and My.
        // TODO: isn't there a better way to come up with "the rotation around x&y that
        // best approximates a rotation R"? Here we're just hoping that the supplied
        // rotation matrix can be decomposed into (x,y) rotations.
        const Vec2 angles12 = R_MbM.convertToSpaceFixed12();
        toQ(q) = angles12;
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    // We can only express angular velocity that can be produced with our generalized
    // speeds which are Mx and My rotations rates.
    void setMobilizerAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {
        const Rotation R_MbM = Rotation::aboutXThenOldY(fromQ(q)[0], fromQ(q)[1]); // space fixed 1-2 sequence
        const Vec3     w_MbM_M = ~R_MbM*w_MbM;
        toU(u) = w_MbM.getSubVec<2>(0); // project angular velocity onto xy axes
    }

    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Vec2& a = fromQ(q); // angular coordinates
        toQ(sine)   = Vec2(std::sin(a[0]), std::sin(a[1]));
        toQ(cosine) = Vec2(std::cos(a[0]), std::cos(a[1]));
        // no quaternions
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        // We're only updating the orientation here because a U-joint can't translate.
        X_MbM.updR() = Rotation::aboutXThenOldY(fromQ(q)[0], fromQ(q)[1]); // space fixed 1-2 sequence
        X_MbM.updT() = 0.;
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Rotation&  R_BM  = getX_BM().R();     // fixed
        const Vec3       T_BM  = getX_BM().T();     //   "
        const Transform& X_PMb = getX_PMb();        // fixed
        const Transform& X_GP  = getX_GP(cc);       // calculated earlier
        const Transform& X_GB  = getX_GB(cc);       // just calculated
        const Rotation&  R_GB  = X_GB.R();

        const Rotation   R_GM  = R_GB * R_BM;

        // Vec from OM (and OMb) to OB, expr. in G (note minus sign). Cross
        // this with a rotation axis to get the linear velocity induced by
        // the off-axis positioning of OB.
        const Vec3 r_MB_G = -(R_GB*T_BM);
        const Vec3 r_MP_G = -X_GP.R()*X_PMb.T();

        // The rotational speeds ux,uy are the Euler angle derivatives for
        // the space-fixed 1-2 sequence represented by qx,qy. These are also
        // (apparently!) the x,y measure numbers of w_MbM_M so the generalized speeds
        // translate directly into velocities.

        //TODO: this doesn't work when M != B, although it works fine when
        //Mb != P.
        H[0] = SpatialRow(~R_GM.x(), ~(R_GM.x() % r_MB_G)); // == r_MbB_G
        H[1] = SpatialRow(~R_GM.y(), ~(R_GM.y() % r_MB_G));
    }
};

/**
 * The "diatom" joint is the equivalent of a free joint for a body with no inertia in
 * one direction, such as one composed of just two atoms. It allows unrestricted
 * translation but rotation only about directions perpendicular to the body's
 * inertialess axis, which we define to be the M frame's z axis Mz.
 *
 * So rotation is allowed about the two axes
 * perpendicular to the outboard body's Mz axis (that is, about Mx and My), but
 * never allows rotation about Mz. Thus the generalized coordinates are:
 *   - a 1-2 body fixed Euler angle sequence (same as 1-2-3 with the 3rd rotation zero)
 *   - then 3 translations, comprising the vector r_OMb_OM, but expressed in M.
 * The generalized speeds are:
 *   - the 1-2 body fixed Euler angle derivatives
 *   - then the 3 vector velocity of OM in Mb, expressed in M
 * This means that all the qdots are just u's, but that the rotational u's are
 * not angular velocity components. The translational qdot's are just the
 * translational u's (but remember that these are in the M frame).
 *
 */
class RBNodeTranslateRotate2 : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "diatom"; }

    RBNodeTranslateRotate2(const MassProperties& mProps_B,
                           const Transform&      X_PMb,
                           const Transform&      X_BM,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setMobilizerRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotations this joint can handle are about Mx and My.
        // TODO: isn't there a better way to come up with "the rotation around x&y that
        // best approximates a rotation R"?
        //const Vec2 angles12 = R_MbM.convertToBodyFixed12();
        const Vec2 angles12 = R_MbM.convertToSpaceFixed12();
        toQ(q).updSubVec<2>(0) = angles12;
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // This joint allows general translation, so we never have to modify any rotations here.
        // However, the translational generalized coordinates are in the M frame here, while the
        // user gives us the translation in the Mb frame, so we have to reexpress.
        //const Rotation R_MbM = Rotation::aboutXThenNewY(fromQ(q)[0], fromQ(q)[1]); // body 1-2
        const Rotation R_MbM = Rotation::aboutXThenOldY(fromQ(q)[0], fromQ(q)[1]); // space 1-2
        toQVec3(q,2) = ~R_MbM*T_MbM; // skip the 2 Euler angles
    }

    // We can only express angular velocity with no component along Mz. The user provides the
    // angular velocity in the Mb frame, so we have to reexpress in M and then project on the
    // (Mx,My) plane. TODO
    void setMobilizerAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {
        //const Rotation R_MbM = Rotation::aboutXThenNewY(fromQ(q)[0], fromQ(q)[1]); // body fixed 1-2 sequence
        const Rotation R_MbM = Rotation::aboutXThenOldY(fromQ(q)[0], fromQ(q)[1]); // space fixed 1-2 sequence
        const Vec3     w_MbM_M = ~R_MbM*w_MbM;
        toU(u).updSubVec<2>(0) = w_MbM.getSubVec<2>(0); // project angular velocity onto xy plane
    }

    // TODO -- angular velocity must be calculated from Euler angle derivatives
    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector& q, const Vec3& v_MbM, Vector& u, bool only) const
    {
        //const Rotation R_MbM   = Rotation::aboutXThenNewY(fromQ(q)[0], fromQ(q)[1]); // body fixed 1-2 sequence
        const Rotation R_MbM   = Rotation::aboutXThenOldY(fromQ(q)[0], fromQ(q)[1]); // space fixed 1-2 sequence
        Vec3           r_MbM_M = fromQVec3(q,2); // skip 2 euler angles to get to translation
        Vec3           v_MbM_M = ~R_MbM*v_MbM;

        // Now we have the cross-joint transform, except the translation is expressed in M.

        const Vec3 w_MbM_M(fromU(u)[0], fromU(u)[1], 0);
        const Vec3 v_wXr_MbM_M = w_MbM_M % r_MbM_M; // linear velocity of OM induced by angular velocity, exp in M

        // Calculate the residual velocity we have to add, then express in outboard frame
        toUVec3(u,3) = v_MbM_M - v_wXr_MbM_M;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Vec2& a = fromQ(q).getSubVec<2>(0); // angular coordinates
        toQ(sine).updSubVec<2>(0)   = Vec2(std::sin(a[0]), std::sin(a[1]));
        toQ(cosine).updSubVec<2>(0) = Vec2(std::cos(a[0]), std::cos(a[1]));
        // no quaternions
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const 
    {
        const Vec<2>& angles  = fromQ(q).getSubVec<2>(0);
        const Vec<3>& r_MbM_M = fromQ(q).getSubVec<3>(2); // skip 2 euler angles to get to translation

        //X_MbM.updR() = Rotation::aboutXThenNewY(angles[0], angles[1]); // body fixed 1-2
        X_MbM.updR() = Rotation::aboutXThenOldY(angles[0], angles[1]); // space fixed 1-2
        //X_MbM.updT() = X_MbM.R()*r_MbM_M; // because translation is in M
        X_MbM.updT() = r_MbM_M; // because translation is in Mb
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();   // fixed
        const Transform& X_PMb = getX_PMb();  // fixed
        const Transform& X_GP  = getX_GP(cc); // calculated earlier
        const Transform& X_GB  = getX_GB(cc); // just calculated
        const Transform& X_MbM = getX_MbM(cc);

        const Transform X_GMb = X_GP * X_PMb;
        const Transform X_GM  = X_GB * X_BM;
        const Rotation& R_GM  = X_GM.R();
        const Rotation& R_GMb = X_GMb.R();
        const Rotation& R_MbM = X_MbM.R();

        const Vec3 r_MbB_G = X_GB.T() - X_GMb.T();
        const Vec3 r_MB_G  = -(X_GB.R()*X_BM.T());

        // The rotational speeds ux,uy are the Euler angle derivatives for
        // the body-fixed 1-2 sequence represented by qx,qy. The angular
        // velocity w=w_MbM_M induced by these speeds is
        //  [ wx ]   [ cqy 0 ]   [ ux ]
        //  [ wy ] = [  0  1 ] * [ uy ]
        //  [ wz ]   [ sqy 0 ]
        // The first column cqy,0,sqy is the first row of R_MbM.
        //
        // Translation occurs about the new M frame axes.
        const Vec3 x_M = /*~(R_MbM[0])*/Vec3(1,0,0);
        const Vec3 y_M = Vec3(0,1,0);
        const Vec3 x = R_GM*x_M; // axes in ground
        const Vec3 y = R_GM*y_M;

        H[0] = SpatialRow(~R_GM.x(), ~(R_GM.x() % r_MbB_G));
        H[1] = SpatialRow(~R_GM.y(), ~(R_GM.y() % r_MbB_G));
        //H[0] = SpatialRow(~x, ~(x % r_MbB_G));
        //H[1] = SpatialRow(~y, ~(y % r_MbB_G));
        H[2] = SpatialRow(  Row3(0) ,     ~R_GM.x());
        H[3] = SpatialRow(  Row3(0) ,     ~R_GM.y());
        H[4] = SpatialRow(  Row3(0) ,     ~R_GM.z());    
    }

};

/// Ball joint. This provides three degrees of rotational freedom,  i.e.,
/// unrestricted orientation of the body's M frame in the parent's Mb frame.
/// The generalized coordinates are:
///   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
/// and generalized speeds are:
///   * angular velocity w_MbM as a vector expressed in the Mb frame.
/// Thus rotational qdots have to be derived from the generalized speeds to
/// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeRotate3 : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "rotate3"; }

    RBNodeRotate3(const MassProperties& mProps_B,
                  const Transform&      X_PMb,
                  const Transform&      X_BM,
                  int&                  nextUSlot,
                  int&                  nextUSqSlot,
                  int&                  nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQ(q)    = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    void setMobilizerTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setMobilizerAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM,
                                     Vector& u) const
    {
            toU(u) = w_MbM[0]; // relative angular velocity always used as generalized speeds
    }

    void setMobilizerLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q); // angular coordinates
            toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(mc.quaternionIndex[nodeNum] >= 0);
            qErr[mc.firstQuaternionQErrSlot+mc.quaternionIndex[nodeNum]] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const
    {
        X_MbM.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(mv))
            X_MbM.updR().setToBodyFixed123(fromQ(q));
        else {
            // TODO: should use qnorm pool
            X_MbM.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_GP  = getX_GP(cc);  // calculated earlier
        const Transform& X_GB  = getX_GB(cc);  // just calculated

        // Vec from OM (and OMb) to OB, expr. in G (note minus sign). Cross
        // this with a rotation axis to get the linear velocity induced by
        // the off-axis positioning of OB.
        const Vec3 T_MB_G = -X_GB.R()*X_BM.T();

        // The generalized speeds are defined in the space-fixed 
        // (that is, Mb) frame, so the orientation of Mb in ground gives
        // the instantaneous spatial meaning of those coordinates. 
        const Rotation R_GMb = X_GP.R() * X_PMb.R();
        H[0] = SpatialRow(~R_GMb.x(), ~(R_GMb.x() % T_MB_G));
        H[1] = SpatialRow(~R_GMb.y(), ~(R_GMb.y() % T_MB_G));
        H[2] = SpatialRow(~R_GMb.z(), ~(R_GMb.z() % T_MB_G));
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const Vec3& w_MbM = fromU(u); // angular velocity of M in Mb 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_MbM = getX_MbM(cc).R();
            toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(q),
                                        ~R_MbM*w_MbM); // need w in *body*, not parent
        } else
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(q),w_MbM);
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3& w_MbM     = fromU(u); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_MbM_dot = fromU(udot);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_MbM = getX_MbM(cc).R();
            toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                                  (fromQ(q), ~R_MbM*w_MbM, ~R_MbM*w_MbM_dot);
        } else
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(q),w_MbM,w_MbM_dot);
    }

    void copyQ(
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
    bool isUsingQuaternion(const SBModelVars& mv) const {
        return !getUseEulerAngles(mv);
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
/// TODO: to get this to work I had to make the translations be in the outboard
/// frame (M, not Mb). So currently the generalized coordinates are:
///   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
///   * translation from OMb to OM as a 3-vector in the outboard body mobilizer (M) frame
/// and generalized speeds are:
///   * angular velocity w_MbM as a vector expressed in the Mb frame
///   * linear velocity of the M origin in Mb (v_MbM), expressed in M
/// Thus translational qdots are just generalized speeds, but rotational
/// qdots have to be derived from the generalized speeds to be turned into
/// either 4 quaternion derivatives or 3 Euler angle derivatives.
///   
class RBNodeTranslateRotate3 : public RigidBodyNodeSpec<6> {
public:
    virtual const char* type() { return "full"; }

    RBNodeTranslateRotate3(const MassProperties& mProps_B,
                           const Transform&      X_PMb,
                           const Transform&      X_BM,
                           int&                  nextUSlot,
                           int&                  nextUSqSlot,
                           int&                  nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setMobilizerRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQVec3(q,0) = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    // The user gives us the translation vector from OMb to OM as a vector expressed in Mb,
    // but the generalized coordinates here constitute a vector in M. So we have to re-express
    // the translation using the current orientation. However, with a free joint 
    // we never have to *change* orientation coordinates in order to achieve a translation.
    // Note: a quaternion from a state is not necessarily normalized so can't be used
    // direction as though it were a set of Euler parameters; it must be normalized first.
    void setMobilizerTranslation(const SBModelVars& mv, const Vec3& T_MbM, Vector& q, bool only) const {
        Rotation R_MbM;
        if (getUseEulerAngles(mv)) {
            R_MbM.setToBodyFixed123(fromQVec3(q,0));
            toQVec3(q,3) = ~R_MbM*T_MbM; // skip the 3 Euler angles
        } else {
            R_MbM.setToQuaternion(Quaternion(fromQuat(q))); // must normalize q before use!
            toQVec3(q,4) = ~R_MbM*T_MbM; // skip the 4 quaternions
        }
    }

    // Our 3 rotational generalized speeds are just the angular velocity vector of M in Mb,
    // expressed in Mb, which is exactly what the user provides here.
    void setMobilizerAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM,
                                     Vector& u) const
    {
            toUVec3(u,0) = w_MbM; // relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in Mb,
    // but expressed in M. The user gives us that same vector, but expressed in Mb so we
    // have to reexpress it using the orientation currently contained in the passed-in qs.
    // Then we have to subtract off the linear velocity that is already present because
    // of the w X r term from the current angular velocity of M in Mb and the offset from
    // the Mb origin to M's origin.
    void setMobilizerLinearVelocity
       (const SBModelVars& mv, const Vector& q, const Vec3& v_MbM, Vector& u, bool only) const
    {
        Rotation R_MbM;
        Vec3     r_MbM; // translation vector from OMb to OM, in Mb
        if (getUseEulerAngles(mv)) {
            R_MbM.setToBodyFixed123(fromQVec3(q,0));
            r_MbM = R_MbM*fromQVec3(q,3); // skip 3 euler angles to get to translation
        } else {
            R_MbM.setToQuaternion(Quaternion(fromQuat(q))); // normalizing
            r_MbM = R_MbM*fromQVec3(q,4); // skip 4 quaternions to get to translation
        }

        // Now we have the cross-joint transform X_MbM=(R_MbM,r_MbM)

        const Vec3 w_MbM = fromUVec3(u,0);
        const Vec3 v_wXr_MbM = w_MbM % r_MbM; // linear velocity of OM induced by angular velocity

        // Calculate the residual velocity we have to add, then express in outboard frame
        toUVec3(u,3) = ~R_MbM * (v_MbM - v_wXr_MbM);
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(mc.quaternionIndex[nodeNum] >= 0);
            qErr[mc.firstQuaternionQErrSlot+mc.quaternionIndex[nodeNum]] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_MbM.
    void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform& X_MbM) const 
    {
        if (getUseEulerAngles(mv)) {
            X_MbM.updR().setToBodyFixed123(fromQVec3(q,0));
            X_MbM.updT() = X_MbM.R()*fromQVec3(q,3);
        } else {
            X_MbM.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
            X_MbM.updT() = X_MbM.R()*fromQVec3(q,4);  // because translation is in M frame
        }
    }

    // Calculate H.
    void calcJointTransitionMatrix(
        const SBPositionCache& cc, 
        HType&                 H) const
    {
        const Transform& X_BM  = getX_BM();   // fixed
        const Transform& X_PMb = getX_PMb();  // fixed
        const Transform& X_GP  = getX_GP(cc); // calculated earlier
        const Transform& X_GB  = getX_GB(cc); // just calculated

        const Transform X_GMb = X_GP * X_PMb;
        const Transform X_GM  = X_GB * X_BM;

        const Vec3 T_MbB_G = X_GB.T() - X_GMb.T();

        // The rotational speeds (angular velocity of M in Mb)
        // are defined in the space-fixed (that is, Mb) frame,
        // so the orientation of Mb in ground gives
        // the instantaneous spatial meaning of those speeds. 
        // *Then* we translate along the new M axes.

        H[0] = SpatialRow(~X_GMb.x(), ~(X_GMb.x() % T_MbB_G));
        H[1] = SpatialRow(~X_GMb.y(), ~(X_GMb.y() % T_MbB_G));
        H[2] = SpatialRow(~X_GMb.z(), ~(X_GMb.z() % T_MbB_G));
        H[3] = SpatialRow(  Row3(0) ,     ~X_GM.x());
        H[4] = SpatialRow(  Row3(0) ,     ~X_GM.y());
        H[5] = SpatialRow(  Row3(0) ,     ~X_GM.z());
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& cc,
        const Vector&          u,
        Vector&                qdot) const
    {
        const Vec3& w_MbM = fromUVec3(u,0); // Angular velocity
        const Vec3& v_MbM = fromUVec3(u,3); // Linear velocity
        if (getUseEulerAngles(mv)) {
            const Rotation& R_MbM = getX_MbM(cc).R();
            const Vec3& theta = fromQVec3(q,0); // Euler angles
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
                                            ~R_MbM*w_MbM); // need w in *body*, not parent
            toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdot,3) = v_MbM;
        } else {
            const Vec4& quat = fromQuat(q);
            toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,w_MbM);
            toQVec3(qdot,4) = v_MbM;
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
        const Vec3& w_MbM     = fromUVec3(u,0); // angular velocity of M in Mb
        const Vec3& v_MbM     = fromUVec3(u,3); // linear velocity
        const Vec3& w_MbM_dot = fromUVec3(udot,0);
        const Vec3& v_MbM_dot = fromUVec3(udot,3);
        if (getUseEulerAngles(mv)) {
            const Rotation& R_MbM = getX_MbM(cc).R();
            const Vec3& theta  = fromQVec3(q,0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, ~R_MbM*w_MbM, ~R_MbM*w_MbM_dot);
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_MbM_dot;
            //cout << "   w_MbM_dot=" << w_MbM_dot << "  v_MbM_dot=" << v_MbM_dot << endl;
            //cout << "   qdotdot=" << fromQVec3(qdotdot,0) << endl;
        } else {
            const Vec4& quat  = fromQuat(q);
            toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                             (quat,w_MbM,w_MbM_dot);
            toQVec3(qdotdot,4) = v_MbM_dot;
        }
    }

    void copyQ(Vector& q, const SBModelVars& mv, const Vector& qIn) const {
        if (getUseEulerAngles(mv))
            toQ(q) = fromQ(qIn);
        else {
            toQuat(q)    = fromQuat(qIn);
            toQVec3(q,4) = fromQVec3(qIn,4);
        }
    }

    int  getMaxNQ()                   const {return 7;}
    int  getNQ(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 
    bool isUsingQuaternion(const SBModelVars& mv) const {
        return !getUseEulerAngles(mv);
    }

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
    const Transform&         X_PMb,        // parent's attachment frame for this joint
    const Transform&         X_BM,         // inboard joint frame J in body frame
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
        return new RBNodeTorsion(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Universal:        
        return new RBNodeRotate2(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Orientation:
        return new RBNodeRotate3(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Cartesian:
        return new RBNodeTranslate(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::FreeLine:
        return new RBNodeTranslateRotate2(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Free:
        return new RBNodeTranslateRotate3(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Sliding:
        return new RBNodeSlider(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::Cylinder:
        return new RBNodeCylinder(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
    case Mobilizer::BendStretch:
        return new RBNodeBendStretch(m,X_PMb,X_BM,nxtUSlot,nxtUSqSlot,nxtQSlot);
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

