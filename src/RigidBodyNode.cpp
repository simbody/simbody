/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
#include "MobilizedBodyRep.h"

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
    SBPositionCache&   pc) const
{
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 T_PB_G = getX_GP(pc).R() * getX_PB(pc).T();

    // The Phi matrix conveniently performs child-to-parent (inward) shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // (outward) shifting for velocities.
    updPhi(pc) = PhiMatrix(T_PB_G);

    // Calculate spatial mass properties. That means we need to transform
    // the local mass moments into the Ground frame and reconstruct the
    // spatial inertia matrix Mk.

    updInertia_OB_G(pc) = getInertia_OB_B().reexpress(~getX_GB(pc).R());
    updCB_G(pc)         = getX_GB(pc).R()*getCOM_B();

    updCOM_G(pc) = getX_GB(pc).T() + getCB_G(pc);

    // Calc Mk: the spatial inertia matrix about the body origin.
    // Note that this is symmetric; offDiag is *skew* symmetric so
    // that transpose(offDiag) = -offDiag.
    // Note: we need to calculate this now so that we'll be able to calculate
    // kinetic energy without going past the Velocity stage.
    const Mat33 offDiag = getMass()*crossMat(getCB_G(pc));
    updMk(pc) = SpatialMat( getInertia_OB_G(pc).toMat33() ,     offDiag ,
                                   -offDiag             , getMass()*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (V_GB). This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBPositionCache& pc,
    SBVelocityCache&       mc) const
{
    updV_GB(mc) = ~getPhi(pc)*parent->getV_GB(mc) + getV_PB_G(mc);
}

Real RigidBodyNode::calcKineticEnergy(
    const SBPositionCache& pc,
    const SBVelocityCache& mc) const 
{
    const Real ret = dot(getV_GB(mc) , getMk(pc)*getV_GB(mc));
    return 0.5*ret;
}

// Calculate velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and coriolis acceleration.
// This routine expects that all spatial velocities & spatial inertias are
// already available.
// Must be called base to tip.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(
    const SBPositionCache& pc,
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

    const Vec3& w_GB = getV_GB(mc)[0];  // spatial angular velocity
    const Vec3& v_GB = getV_GB(mc)[1];  // spatial linear velocity (of B origin in G)

    updGyroscopicForce(dc) = 
        SpatialVec(    w_GB % (getInertia_OB_G(pc)*w_GB),     // gyroscopic moment
                    getMass()*(w_GB % (w_GB % getCB_G(pc)))); // gyroscopic force

    // Parent velocity.
    const Vec3& w_GP = parent->getV_GB(mc)[0];
    const Vec3& v_GP = parent->getV_GB(mc)[1];

    // Calc a: coriolis acceleration.
    // The coriolis acceleration "a" is a 
    // "remainder" term in the spatial acceleration, depending only on velocities,
    // but involving time derivatives of the Phi and H matrices. Specifically,
    //   a = ~PhiDot * V_GP + ~HDot * u
    // As correctly calculated in Schwieters' paper, Eq [16], the first term above
    // simplifies to SpatialVec( 0, w_GP % (v_GB-v_GP) ). However, Schwieters' second
    // term in [16] is correct only if H is constant in P, in which case the derivative
    // just accounts for the rotation of P in G. In general H is not constant in P,
    // so we don't try to calculate the derivative here but assume that ~HDot*u has
    // already been calculated for us and stored in VD_PB_G. (That is,
    // V_PB_G = ~H*u, VD_PB_G = ~HDot*u.)

    updCoriolisAcceleration(dc) =
        SpatialVec( Vec3(0), w_GP % (v_GB-v_GP) ) + getVD_PB_G(dc);

    updTotalCoriolisAcceleration(dc) =
        ~getPhi(pc) * parent->getTotalCoriolisAcceleration(dc)
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
        const SBInstanceVars& iv,
        SBInstanceCache&      ic) const {}

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

    /*virtual*/ void realizeDynamics(
        const SBModelVars&     mv,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBVelocityCache& vc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/ void setQToFitTransform
       (const SBModelVars&, const Transform& X_MbM, Vector& q) const {}
    /*virtual*/ void setQToFitRotation
       (const SBModelVars&, const Rotation& R_MbM, Vector& q) const {}
    /*virtual*/ void setQToFitTranslation
       (const SBModelVars&, const Vec3& T_MbM, Vector& q,
        bool dontChangeOrientation) const {}

    /*virtual*/ void setUToFitVelocity
       (const SBModelVars&, const Vector& q, const SpatialVec& V_MbM, Vector& u) const {}
    /*virtual*/ void setUToFitAngularVelocity
       (const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {}
    /*virtual*/ void setUToFitLinearVelocity
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
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/ void calcInternalGradientFromSpatial(
        const SBPositionCache&      pc, 
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

    /*virtual*/void setMobilizerDefaultModelValues(const SBTopologyCache&, 
                                          SBModelVars& v) const
    {
        v.prescribed[0] = true; // ground's motion is prescribed to zero
    }

    // /*virtual*/ const SpatialRow& getHRow(int i) const;
};

// This still-abstract class is a skeleton implementation of a built-in mobilizer, with the
// number of mobilities (dofs, u's) specified as a template parameter. That way all
// the code that is common except for the dimensionality of the mobilizer can be
// written once, and the compiler generates specific implementations for each of the
// six possible dimensionalities (1-6 mobilities). Each implementation works only on
// fixed size Vec<> and Mat<> types, so can use very high speed inline operators.
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

    // This is the type of the joint transition matrix H, which we're
    // viewing as the transpose of the "more important" ~H, whose type
    // is consequently nicer: Mat<2,dof,Vec3>.
    typedef Mat<dof,2,Row3,1,2> HType;


    // Provide default implementations for setQToFitTransform() and setQToFitVelocity() 
    // which are implemented using the rotational and translational quantity routines. These assume
    // that the rotational and translational coordinates are independent, with rotation handled
    // first and then left alone. If a mobilizer type needs to deal with rotation and
    // translation simultaneously, it should provide a specific implementation for these two routines.
    // *Each* mobilizer must implement setQToFit{Rotation,Translation} and 
    // setUToFit{AngularVelocity,LinearVelocity}; there are no defaults.

    virtual void setQToFitTransform(const SBModelVars& mv, const Transform& X_MbM, Vector& q) const {
        setQToFitRotation   (mv,X_MbM.R(),q);
        setQToFitTranslation(mv,X_MbM.T(),q,true); // don't fiddle with the rotation
    }

    virtual void setUToFitVelocity(const SBModelVars& mv, const Vector& q, const SpatialVec& V_MbM, Vector& u) const {
        setUToFitAngularVelocity(mv,q,V_MbM[0],u);
        setUToFitLinearVelocity (mv,q,V_MbM[1],u,true); // don't fiddle with the angular velocity
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
    /// This method constitutes the *definition* of the generalized coordinates for
    /// a particular joint.
    virtual void calcAcrossJointTransform(
        const SBModelVars& mv,
        const Vector&      q,
        Transform&         X_MbM) const=0;


    /// This mandatory routine calculates the joint transition matrix H_MbM, giving
    /// the change of velocity induced by the generalized speeds u for this 
    /// mobilizer, expressed in the mobilizer's inboard frame Mb (attached to
    /// this body's parent). 
    /// This method constitutes the *definition* of the generalized speeds for
    /// a particular joint.
    /// This routine can depend on X_MbM having already
    /// been calculated and available in the PositionCache but must NOT depend
    /// on any quantities involving Ground or other bodies.
    virtual void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const=0;

    /// This mandatory routine calculates the time derivative taken in Mb of
    /// the matrix H_MbM (call it H_MbM_Dot). This is zero if the generalized
    /// speeds are all defined in terms of the Mb frame, which is often the case.
    /// This routine can depend on X_MbM and H_MbM being available already in
    /// the state PositionCache, and V_MbM being already in the state VelocityCache.
    /// However, it must NOT depend on any quantities involving Ground or other bodies.
    virtual void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const = 0;

    /// This routine is NOT joint specific, but cannot be called until the across-joint
    /// transform X_MbM has been calculated and is available in the State cache.
    void calcBodyTransforms(
        const SBPositionCache& pc, 
        Transform&             X_PB, 
        Transform&             X_GB) const 
    {
        const Transform& X_BM  = getX_BM();  // fixed
        const Transform& X_PMb = getX_PMb(); // fixed
        const Transform& X_MbM = getX_MbM(pc); // just calculated
        const Transform& X_GP  = getX_GP(pc);  // already calculated

        X_PB = X_PMb * X_MbM * ~X_BM; // TODO: precalculate X_MB
        X_GB = X_GP * X_PB;
    }

    // Same for all mobilizers.
    void calcParentToChildVelocityJacobianInGround(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType& H_PB_G) const;

    // Same for all mobilizers.
    void calcParentToChildVelocityJacobianInGroundDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        const SBDynamicsCache& dc, 
        HType& H_PB_G_Dot) const;


    /// Calculate joint-specific kinematic quantities dependent on
    /// on velocities. This routine may assume that *all* position 
    /// kinematics (not just joint-specific) has been done for this node,
    /// that all velocity kinematics has been done for the parent, and
    /// that the velocity state variables (u) are available. The
    /// quanitites that must be computed are:
    ///   V_MbM   relative velocity of B's M frame in P's Mb frame, 
    ///             expressed in Mb (note: this is also V_PM_Mb since
    ///             Mb is fixed on P).
    ///   V_PB_G  relative velocity of B in P, expr. in G
    /// The code is the same for all joints, although parametrized by ndof.
    void calcJointKinematicsVel(
        const SBPositionCache& pc,
        const Vector&          u,
        SBVelocityCache&       vc) const 
    {
        updV_MbM(vc)  = ~getH_MbM(pc) * fromU(u);
        updV_PB_G(vc) = ~getH(pc) * fromU(u);
    }

    /// Calculate joint-specific dynamics quantities dependent on velocities.
    /// This method may assume that *all* position & velocity kinematics
    /// (not just joint-specific) has been done for this node, and that
    /// dynamics has been done for the parent.
    void calcJointDynamics(
        const SBPositionCache& pc,
        const Vector&          u,
        const SBVelocityCache& vc, 
        SBDynamicsCache&       dc) const 
    {
        updVD_PB_G(dc) = ~getHDot(dc) * fromU(u);
    }

    // These next two routines are options, but if you supply one you
    // must supply the other. (That is, ball-containing joints provide
    // both of these routines.)
    virtual void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        Vector&                qdot) const
    {
        toQ(qdot) = fromU(u);        // default is qdot=u
    }

    virtual void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
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
        const SBInstanceVars& iv,
        SBInstanceCache&      ic) const
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

        calcAcrossJointVelocityJacobian(mv, pc, updH_MbM(pc));
        calcParentToChildVelocityJacobianInGround(mv,pc, updH(pc));

        //calcJointTransitionMatrix(pc, updH(pc)); // TODO: soon to be obsolete

        calcJointIndependentKinematicsPos(pc);
    }

    // Set new velocities for the current configuration, and calculate
    // all the velocity-dependent terms. Must call base-to-tip.
    void realizeVelocity(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        SBVelocityCache&       vc,
        Vector&                qdot) const 
    {
        calcQDot(mv,q,pc,u,qdot);
        calcJointKinematicsVel(pc,u,vc);
        calcJointIndependentKinematicsVel(pc,vc);
    }

    void realizeDynamics(
        const SBModelVars&     mv,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBVelocityCache& vc,
        SBDynamicsCache&       dc) const
    {
        // Mobilizer-specific.
        calcAcrossJointVelocityJacobianDot          (mv,pc,vc, updH_MbM_Dot(dc));
        calcParentToChildVelocityJacobianInGroundDot(mv,pc,vc, dc, updHDot(dc));

        calcJointDynamics(pc,u,vc,dc);

        // Mobilizer independent.
        calcJointIndependentDynamicsVel(pc,vc,dc);
    }

    // This is a dynamics-stage calculation and must be called tip-to-base (inward).
    void calcArticulatedBodyInertiasInward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const;

    // calcJointIndependentDynamicsVel() must be called after ArticulatedBodyInertias.

    // This dynamics-stage calculation is needed for handling constraints. It
    // must be called base-to-tip (outward);
    void calcYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const;

    // These routines give each node a chance to set appropriate defaults in a piece
    // of the state corresponding to a particular stage. Default implementations here
    // assume non-ball joint; override if necessary.
    virtual void setMobilizerDefaultModelValues   (const SBTopologyCache&, SBModelVars&)  const {}
    virtual void setMobilizerDefaultInstanceValues(const SBModelVars&, SBInstanceVars&) const {}
    virtual void setMobilizerDefaultTimeValues    (const SBModelVars&, SBTimeVars&)    const {}

    virtual void setMobilizerDefaultPositionValues(const SBModelVars& s, Vector& q) const 
    {
        toQ(q) = 0.;
    }
    virtual void setMobilizerDefaultVelocityValues(const SBModelVars&, Vector& u) const 
    {
        toU(u) = 0.;
    }
    virtual void setMobilizerDefaultDynamicsValues(const SBModelVars&, SBDynamicsVars&) const {}
    virtual void setMobilizerDefaultAccelerationValues(const SBModelVars&, 
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
        const SBPositionCache& pc, 
        const SBVelocityCache& mc,
        const SpatialVec&      sVel, 
        Vector&               u) const;

    // Return true if any change is made to the q vector.
    virtual bool enforceQuaternionConstraints(
        const SBModelVars& mv,
        Vector&            q) const 
    {
        return false;
    }

    const SpatialRow& getHRow(const SBPositionCache& pc, int i) const {
        return getH(pc)[i];
    }

    // Access to body-oriented state and cache entries is the same for all nodes,
    // and joint oriented access is almost the same but parametrized by dof. There is a special
    // case for quaternions because they use an extra state variable, and although we don't
    // have to we make special scalar routines available for 1-dof joints. Note that all State access
    // routines are inline, not virtual, so the cost is just an indirection and an index.
    //
    // TODO: these inner-loop methods probably shouldn't be indexing a Vector, which requires
    // several indirections. Instead, the top-level caller should find the Real* data contained in the
    // Vector and then pass that to the RigidBodyNode routines which will call these ones.

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

    // Applied forces from acceleration variables.
    const Vec<dof>&   getAppliedJointForce(const SBAccelerationVars& av) const 
        {return fromU(av.appliedMobilityForces);}
    //TODO
    const Vec<dof>&   getPrescribedUdot   (const SBAccelerationVars& av) const 
        {return fromU(av.prescribedUdot);}

    // Special case state access for 1-dof joints
    const Real& get1AppliedJointForce(const SBAccelerationVars& av) const {return from1U(av.appliedMobilityForces);}
    const Real& get1PrescribedUdot   (const SBAccelerationVars& av) const {return from1U(av.prescribedUdot);}

    // Cache entries (cache is mutable in a const State)

        // Position


    // TODO: should store as H_MbM or else always reference Ht_MbM
    const Mat<dof,2,Row3,1,2>& getH_MbM(const SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::getAs(&pc.storageForHtMbM(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH_MbM(SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::updAs(&pc.storageForHtMbM(0,uIndex)); }

    // "H" here should really be H_PB_G, that is, cross joint transition
    // matrix relating parent and body frames, but expressed in Ground.
    // TODO: should store as H or else always reference Ht
    const Mat<dof,2,Row3,1,2>& getH(const SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::getAs(&pc.storageForHt(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH(SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::updAs(&pc.storageForHt(0,uIndex)); }

    // These are sines and cosines of angular qs. The rest of the slots are garbage.
    const Vec<dof>&   getSinQ (const SBPositionCache& pc) const {return fromQ (pc.sq);}
    Vec<dof>&         updSinQ (SBPositionCache&       pc) const {return toQ   (pc.sq);}
    const Real&       get1SinQ(const SBPositionCache& pc) const {return from1Q(pc.sq);}
    Real&             upd1SinQ(SBPositionCache&       pc) const {return to1Q  (pc.sq);}

    const Vec<dof>&   getCosQ (const SBPositionCache& pc) const {return fromQ (pc.cq);}
    Vec<dof>&         updCosQ (SBPositionCache&       pc) const {return toQ   (pc.cq);}
    const Real&       get1CosQ(const SBPositionCache& pc) const {return from1Q(pc.cq);}
    Real&             upd1CosQ(SBPositionCache&       pc) const {return to1Q  (pc.cq);}

    // These are normalized quaternions in slots for balls. Everything else is garbage.
    const Vec4&       getQNorm(const SBPositionCache& pc) const {return fromQuat(pc.qnorm);}
    Vec4&             updQNorm(SBPositionCache&       pc) const {return toQuat  (pc.qnorm);}

        // Velocity

        // Dynamics

    const Mat<dof,2,Row3,1,2>& getH_MbM_Dot(const SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::getAs(&dc.storageForHtMbMDot(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH_MbM_Dot(SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::updAs(&dc.storageForHtMbMDot(0,uIndex)); }

    const Mat<dof,2,Row3,1,2>& getHDot(const SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::getAs(&dc.storageForHtDot(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updHDot(SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::updAs(&dc.storageForHtDot(0,uIndex)); }

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
        const SBAccelerationVars&,
        SBAccelerationCache&               ) const;

    void calcAccel(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        const SBDynamicsCache& dc,
        SBAccelerationCache&   rc,
        Vector&                udot,
        Vector&                qdotdot) const;

    void calcInternalGradientFromSpatial(
        const SBPositionCache&      pc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const;

    void calcEquivalentJointForces(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector&                     jointForces) const;

    void calcUDotPass1Inward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const;

    void calcUDotPass2Outward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const;

    void calcMInverseFPass1Inward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               f,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const;

    void calcMInverseFPass2Outward(
        const SBPositionCache&      pc,
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


    // TRANSLATION (CARTESIAN) //

// Translate (Cartesian) joint. This provides three degrees of
// translational freedom which is suitable (e.g.) for connecting a
// free atom to ground. The Cartesian directions are the axes of
// the parent body's Mb frame, with M=Mb when all 3 coords are 0,
// and the orientation of M in Mb is 0 (identity) forever.
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

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // the only rotation this mobilizer can represent is identity
    }
    void setQToFitTranslation(const SBModelVars&, const Vec3&  T_MbM, Vector& q, bool only) const {
        // here's what this joint is really good at!
        toQ(q) = T_MbM;
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // The only angular velocity this can represent is zero.
    }
    void setUToFitLinearVelocity
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

    // Generalized speeds together are the velocity of M's origin in the Mb frame,
    // expressed in Mb. So individually they produce velocity along Mb's x,y,z
    // axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0), Row3(1,0,0) );
        H_MbM[1] = SpatialRow( Row3(0), Row3(0,1,0) );
        H_MbM[2] = SpatialRow( Row3(0), Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

};



    // SLIDING (PRISMATIC) //

// Sliding joint (1 dof translation). The translation is along the x
// axis of the parent body's Mb frame, with M=Mb when the coordinate
// is zero and the orientation of M in Mb frozen at 0 forever.
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

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation a slider can represent is identity.
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // We can only represent the x coordinate with this joint.
        to1Q(q) = T_MbM[0];
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // The only angular velocity a slider can represent is zero.
    }

    void setUToFitLinearVelocity
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

    // The generalized speed is the velocity of M's origin in the Mb frame,
    // along Mb's x axis, expressed in Mb.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0), Row3(1,0,0) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // PIN (TORSION) //

// This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
// about a particular axis, the z axis of the parent's Mb frame, which is 
// aligned forever with the z axis of the body's M frame. In addition, the
// origin points of M and Mb are identical forever.
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

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our pin joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        to1Q(q) = angles123[2];
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_MbM[2]; // project angular velocity onto z axis
    }

    void setUToFitLinearVelocity
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


    // The generalized speed is the angular velocity of M in the Mb frame,
    // about Mb's z axis, expressed in Mb. (This axis is also constant in M.)
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0,0,1), Row3(0) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative in Mb is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};


    // SCREW //

// This is a one-dof "screw" joint, meaning one degree of rotational freedom
// about a particular axis, coupled to translation along that same axis.
// Here we use the common z axis of the Mb and M frames, which remains
// aligned forever. 
// For the generalized coordinate q, we use the rotation angle. For the
// generalized speed u we use the rotation rate, which is also the
// angular velocity of M in Mb (about the z axis). We compute the
// translational position as pitch*q, and the translation rate as pitch*u.
class RBNodeScrew : public RigidBodyNodeSpec<1> {
    Real pitch;
public:
    virtual const char* type() { return "screw"; }

    RBNodeScrew(const MassProperties& mProps_B,
                const Transform&      X_PMb,
                const Transform&      X_BM,
                Real                  p,  // the pitch
                int&                  nextUSlot,
                int&                  nextUSqSlot,
                int&                  nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot),
        pitch(p)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our screw joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        to1Q(q) = angles123[2];
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        to1Q(q) = T_MbM[2]/pitch;
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_MbM[2]; // project angular velocity onto z axis
    }

    void setUToFitLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        to1U(u) = v_MbM[2]/pitch;
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

        X_MbM.updR().setToRotationAboutZ(theta);
        X_MbM.updT() = Vec3(0,0,theta*pitch);
    }


    // The generalized speed is the angular velocity of M in the Mb frame,
    // about Mb's z axis, expressed in Mb. (This axis is also constant in M.)
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0,0,1), Row3(0,0,pitch) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative in Mb is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // CYLINDER //

// This is a "cylinder" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along the same axis. For molecules you can think of this as a combination
// of torsion and bond stretch. The axis used is the z axis of the parent's
// Mb frame, which is aligned forever with the z axis of the body's M frame.
// In addition, the origin points of M and Mb are separated only along the
// z axis; i.e., they have the same x & y coords in the Mb frame. The two
// generalized coordinates are the rotation and the translation, in that order.
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


    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our cylinder joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        toQ(q)[0] = angles123[2];
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // Because the M and Mb origins must lie along their shared z axis, there is no way to
        // create a translation by rotating around z. So the only translation we can represent
        // is that component which is along z.
        toQ(q)[1] = T_MbM[2];
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_MbM[2];
    }

    void setUToFitLinearVelocity
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


    // The generalized speeds are (1) the angular velocity of M in the Mb frame,
    // about Mb's z axis, expressed in Mb, and (2) the velocity of M's origin
    // in Mb, along Mb's z axis. (The z axis is also constant in M for this joint.)
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0,0,1), Row3(0)     );
        H_MbM[1] = SpatialRow( Row3(0),     Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative in Mb is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
    }

};


    // BEND-STRETCH //

// This is a "bend-stretch" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along a perpendicular axis. The z axis of the parent's Mb frame is 
// used for rotation (and that is always aligned with the M frame z axis).
// The x axis of the *M* frame is used for translation; that is, first
// we rotate around z, which moves M's x with respect to Mb's x. Then
// we slide along the rotated x axis. The two
// generalized coordinates are the rotation and the translation, in that order.
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


    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our bend-stretch joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        toQ(q)[0] = angles123[2];
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
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

    void setUToFitAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_MbM[2];
    }

    // If the translational coordinate is zero, we can only represent a linear velocity 
    // of OM in Mb which is along M's current x axis direction. Otherwise, we can 
    // represent any velocity in the x-y plane by introducing angular velocity about z.
    // We can never represent a linear velocity along z.
    void setUToFitLinearVelocity
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

    // The generalized speeds for this bend-stretch joint are (1) the angular
    // velocity of M in the Mb frame, about Mb's z axis, expressed in Mb, and
    // (2) the (linear) velocity of M's origin in Mb, along *M's* current x axis
    // (that is, after rotation about z). (The z axis is also constant in M for this joint.)
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb

        const Vec3&     T_MbM = getX_MbM(pc).T();
        H_MbM[0] = SpatialRow( Row3(0,0,1), ~(Vec3(0,0,1) % T_MbM)   );
        H_MbM[1] = SpatialRow( Row3(0),              ~Mx_Mb          );
    }

    // Since the the Jacobian above is not constant in Mb,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb

        const Vec3&     w_MbM = getV_MbM(vc)[0]; // angular velocity of M in Mb
        const Vec3&     v_MbM = getV_MbM(vc)[1]; // linear velocity of OM in Mb

        H_MbM_Dot[0] = SpatialRow( Row3(0), ~(Vec3(0,0,1) % v_MbM) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), ~(w_MbM % Mx_Mb) );
    }

};

    // UNIVERSAL (U-JOINT, HOOKE'S JOINT) //

// This is a Universal Joint (U-Joint), also known as a Hooke's joint.
// This is identical to the joint that connects pieces of a driveshaft
// together under a car. Physically, you can think of this as
// a parent body P, hinged to a massless cross-shaped coupler, which is then
// hinged to the child body B. The massless coupler doesn't actually appear in the
// model. Instead, we use a body-fixed 1-2 Euler rotation sequence for
// orientation, which has the same effect: starting with frames B and P
// aligned (when q0=q1=0), rotate frame B about the Px(=Bx) axis by q0; then, 
// rotate frame B further about the new By(!=Py) by q1. For generalized
// speeds u we use the Euler angle derivatives qdot, which are *not* the
// same as angular velocity components because u0 is a rotation rate 
// around Px(!=Bx any more) while u1 is a rotation rate about By.
//
// To summarize,
//    q's: a two-angle body-fixed rotation sequence about x, then new y
//    u's: time derivatives of the q's
//
// Note that the U-Joint degrees of freedom relating the parent's Mb frame
// to the child's M frame are about x and y, with the "long" axis of the
// driveshaft along z.
class RBNodeUJoint : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "ujoint"; }

    RBNodeUJoint(const MassProperties& mProps_B,
                 const Transform&      X_PMb,
                 const Transform&      X_BM,
                 int&                  nextUSlot,
                 int&                  nextUSqSlot,
                 int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotations this joint can handle are about Mx and My.
        // TODO: isn't there a better way to come up with "the rotation around x&y that
        // best approximates a rotation R"? Here we're just hoping that the supplied
        // rotation matrix can be decomposed into (x,y) rotations.
        const Vec2 angles12 = R_MbM.convertToBodyFixed12();
        toQ(q) = angles12;
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    // We can only express angular velocity that can be produced with our generalized
    // speeds which are Mbx and My rotations rates. So we'll take the supplied angular velocity
    // expressed in Mb, project it on Mbx and use that as the first generalized speed. Then
    // take whatever angular velocity is unaccounted for, express it in M, and project onto
    // My and use that as the second generalized speed.
    void setUToFitAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM, Vector& u) const {
        const Rotation R_MbM = Rotation::aboutXThenNewY(fromQ(q)[0], fromQ(q)[1]); // body fixed 1-2 sequence
        const Vec3     wyz_MbM_M = ~R_MbM*Vec3(0,w_MbM[1],w_MbM[2]);
        toU(u) = Vec2(w_MbM[0], wyz_MbM_M[1]);
    }

    void setUToFitLinearVelocity
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
        X_MbM.updR() = Rotation::aboutXThenNewY(fromQ(q)[0], fromQ(q)[1]); // body fixed 1-2 sequence
        X_MbM.updT() = 0.;
    }

    // The generalized speeds for this 2-dof rotational joint are the time derivatlves of
    // the body-fixed 1-2 rotation sequence defining the orientation. That is, the first speed
    // is just a rotation rate about Mbx. The second is a rotation rate about the current My, so
    // we have to transform it into Mb to make H_MbM uniformly expressed in Mb.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();

        H_MbM[0] = SpatialRow(  Row3(1,0,0) , Row3(0) );
        H_MbM[1] = SpatialRow( ~(R_MbM.y()) , Row3(0) );
    }

    // Since the second row of the Jacobian H_MbM above is not constant in Mb,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     w_MbM = getV_MbM(vc)[0]; // angular velocity of M in Mb

        H_MbM_Dot[0] = SpatialRow(        Row3(0)      , Row3(0) );
        H_MbM_Dot[1] = SpatialRow( ~(w_MbM % R_MbM.y()), Row3(0) );
    }

};



    // PLANAR //

// This provides free motion (translation and rotation) in a plane. We use
// the 2d coordinate system formed by the x,y axes of Mb as the translations,
// and the common z axis of Mb and M as the rotational axis. The generalized
// coordinates are theta,x,y interpreted as rotation around z and translation
// along the (space fixed) Mbx and Mby axes.
class RBNodePlanar : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "planar"; }

    RBNodePlanar(const MassProperties& mProps_B,
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

    void setQToFitRotation(const SBModelVars&, const Rotation& R_MbM, Vector& q) const {
        // The only rotation our planar joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_MbM.convertToBodyFixed123();
        toQ(q)[0] = angles123[2];
    }
    void setQToFitTranslation(const SBModelVars&, const Vec3&  T_MbM, Vector& q, bool only) const {
        // Ignore translation in the z direction.
        toQ(q)[1] = T_MbM[0]; // x
        toQ(q)[2] = T_MbM[1]; // y
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM, Vector& u) const {
        // We can represent the z angular velocity exactly, but nothing else.
        toU(u)[0] = w_MbM[2];
    }
    void setUToFitLinearVelocity
       (const SBModelVars&, const Vector&, const Vec3& v_MbM, Vector& u, bool only) const
    {
        // Ignore translational velocity in the z direction.
        toU(u)[1] = v_MbM[0]; // x
        toU(u)[2] = v_MbM[1]; // y
    }

    // This is required but does nothing here since there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Real& angle = fromQ(q)[0]; // angular coordinate
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
        // Rotational q is about common z axis, translational q's along Mbx and Mby.
        X_MbM = Transform(Rotation::aboutZ(fromQ(q)[0]), 
                          Vec3(fromQ(q)[1], fromQ(q)[2], 0));
    }

    // The rotational generalized speed is about the common z axis; translations
    // are along Mbx and Mby so all axes are constant in Mb.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(0,0,1),   Row3(0) );
        H_MbM[1] = SpatialRow(   Row3(0),   Row3(1,0,0) );
        H_MbM[2] = SpatialRow(   Row3(0),   Row3(0,1,0) );
    }

    // Since the Jacobian above is constant in Mb, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // ORIENTATION (BALL) //

// Ball joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's Mb frame.
// The generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_MbM as a vector expressed in the Mb frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
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

    void setQToFitRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQ(q)    = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBModelVars&, const Vector&, const Vec3& w_MbM,
                                     Vector& u) const
    {
            toU(u) = w_MbM[0]; // relative angular velocity always used as generalized speeds
    }

    void setUToFitLinearVelocity
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

    // Generalized speeds are the angular velocity expressed in Mb, so they
    // cause rotations around Mb x,y,z axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(1,0,0), Row3(0) );
        H_MbM[1] = SpatialRow( Row3(0,1,0), Row3(0) );
        H_MbM[2] = SpatialRow( Row3(0,0,1), Row3(0) );
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const Vec3& w_MbM = fromU(u); // angular velocity of M in Mb 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_MbM = getX_MbM(pc).R();
            toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(q),
                                        ~R_MbM*w_MbM); // need w in *body*, not parent
        } else
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(q),w_MbM);
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3& w_MbM     = fromU(u); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_MbM_dot = fromU(udot);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_MbM = getX_MbM(pc).R();
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

    void setMobilizerDefaultPositionValues(
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


    // FREE //

// Free joint. This provides six degrees of freedom, three rotational and
// three translational. The rotation is like the ball joint above; the
// translation is like the Cartesian joint above.
// TODO: to get this to work I had to make the translations be in the outboard
// frame (M, not Mb). So currently the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//   * translation from OMb to OM as a 3-vector in the outboard body mobilizer (M) frame
// and generalized speeds are:
//   * angular velocity w_MbM as a vector expressed in the Mb frame
//   * linear velocity of the M origin in Mb (v_MbM), expressed in M
// Thus translational qdots are just generalized speeds, but rotational
// qdots have to be derived from the generalized speeds to be turned into
// either 4 quaternion derivatives or 3 Euler angle derivatives.
//   
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

    void setQToFitRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQVec3(q,0) = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    // The user gives us the translation vector from OMb to OM as a vector expressed in Mb, which
    // is what we use as translational generalized coordinates. Also, with a free joint 
    // we never have to change orientation coordinates in order to achieve a translation.
    void setQToFitTranslation(const SBModelVars& mv, const Vec3& T_MbM, Vector& q, bool only) const {
        if (getUseEulerAngles(mv))
            toQVec3(q,3) = T_MbM; // skip the 3 Euler angles
        else
            toQVec3(q,4) = T_MbM; // skip the 4 quaternions
    }

    // Our 3 rotational generalized speeds are just the angular velocity vector of M in Mb,
    // expressed in Mb, which is exactly what the user provides here.
    void setUToFitAngularVelocity(const SBModelVars&, const Vector& q, const Vec3& w_MbM,
                                     Vector& u) const
    {
        toUVec3(u,0) = w_MbM; // relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in Mb,
    // expressed in Mb, which is just what the user gives us.
    void setUToFitLinearVelocity
       (const SBModelVars& mv, const Vector& q, const Vec3& v_MbM, Vector& u, bool only) const
    {
        toUVec3(u,3) = v_MbM;
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
            X_MbM.updT() = fromQVec3(q,3); // translation is in Mb already
        } else {
            X_MbM.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
            X_MbM.updT() = fromQVec3(q,4); // translation is in Mb already
        }
    }


    // The generalized speeds for this 6-dof ("free") joint are 
    //   (1) the angular velocity of M in the Mb frame, expressed in Mb, and
    //   (2) the (linear) velocity of M's origin in Mb, expressed in Mb.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        H_MbM[0] = SpatialRow( Row3(1,0,0),   Row3(0)   );  // rotations
        H_MbM[1] = SpatialRow( Row3(0,1,0),   Row3(0)   );
        H_MbM[2] = SpatialRow( Row3(0,0,1),   Row3(0)   );

        H_MbM[3] = SpatialRow(   Row3(0),   Row3(1,0,0) );  // translations
        H_MbM[4] = SpatialRow(   Row3(0),   Row3(0,1,0) );
        H_MbM[5] = SpatialRow(   Row3(0),   Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in Mb, its derivative in Mb is 0.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        H_MbM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[2] = SpatialRow( Row3(0), Row3(0) );

        H_MbM_Dot[3] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[4] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[5] = SpatialRow( Row3(0), Row3(0) );
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        Vector&                qdot) const
    {
        const Vec3& w_MbM = fromUVec3(u,0); // Angular velocity in Mb
        const Vec3& v_MbM = fromUVec3(u,3); // Linear velocity in Mb
        if (getUseEulerAngles(mv)) {
            const Rotation& R_MbM = getX_MbM(pc).R();
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
        const SBPositionCache& pc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3& w_MbM     = fromUVec3(u,0); // angular velocity of M in Mb
        const Vec3& v_MbM     = fromUVec3(u,3); // linear velocity of M in Mb, expressed in Mb
        const Vec3& w_MbM_dot = fromUVec3(udot,0);
        const Vec3& v_MbM_dot = fromUVec3(udot,3);
        if (getUseEulerAngles(mv)) {
            const Rotation& R_MbM = getX_MbM(pc).R();
            const Vec3& theta  = fromQVec3(q,0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, ~R_MbM*w_MbM, ~R_MbM*w_MbM_dot);
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_MbM_dot;
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

    void setMobilizerDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
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

    // LINE ORIENTATION //

// LineOrientation joint. Like a Ball joint, this provides full rotational
// freedom, but for a degenerate body which is thin (inertialess) along its
// own z axis. These arise in molecular modeling for linear molecules formed
// by pairs of atoms, or by multiple atoms in a linear arrangement like
// carbon dioxide (CO2) whose structure is O=C=O in a straight line. We are
// assuming that there is no meaning to a rotation about the linear axis,
// so free orientation requires just *two* degrees of freedom, not *three*
// as is required for general rigid bodies. And in fact we can get away with
// just two generalized speeds. But so far, no one has been able to come up
// with a way to manage with only two generalized coordinates, so this joint
// has the same q's as a regular Orientation (Ball) joint: either a quaternion
// for unconditional stability, or a three-angle (body fixed 1-2-3)
// Euler sequence which will be dynamically singular when the middle (y) axis
// is 90 degrees. Use the Euler sequence only for small motions or for kinematics
// problems (and note that only the first two are meaningful).
//
// To summarize, the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * the x,y components of the angular velocity w_MbM_M, that is, the angular
//     velocity of M in Mb expressed in M (where we want wz=0).
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeLineOrientation : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "lineOrientation"; }

    RBNodeLineOrientation(const MassProperties& mProps_B,
                          const Transform&      X_PMb,
                          const Transform&      X_BM,
                          int&                  nextUSlot,
                          int&                  nextUSqSlot,
                          int&                  nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQVec3(q,0)    = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    void setQToFitTranslation(const SBModelVars&, const Vec3& T_MbM, Vector& q, bool only) const {
        // M and Mb frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBModelVars& mv, const Vector& q, const Vec3& w_MbM,
                                     Vector& u) const
    {
        Rotation R_MbM;
        if (getUseEulerAngles(mv))
            R_MbM.setToBodyFixed123(fromQVec3(q,0));
        else {
            // TODO: should use qnorm pool
            R_MbM.setToQuaternion(Quaternion(fromQuat(q))); // normalize
        }
        const Vec3 w_MbM_M = ~R_MbM*w_MbM;
        toU(u) = Vec2(w_MbM_M[0], w_MbM_M[1]); // (x,y) of relative angular velocity always used as generalized speeds
    }

    void setUToFitLinearVelocity
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
            const Vec3& a = fromQVec3(q,0); // angular coordinates
            toQVec3(sine,0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQVec3(cosine,0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
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
            X_MbM.updR().setToBodyFixed123(fromQVec3(q,0));
        else {
            // TODO: should use qnorm pool
            X_MbM.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
        }
    }

    // The generalized speeds for this 2-dof rotational joint are the x and y
    // components of the angular velocity of M in the Mb frame, expressed in the *M*
    // frame.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb
        const Vec3&     My_Mb = R_MbM.y(); // M's y axis, expressed in Mb

        H_MbM[0] = SpatialRow( ~Mx_Mb, Row3(0) );
        H_MbM[1] = SpatialRow( ~My_Mb, Row3(0) );
    }

    // Since the Jacobian above is not constant in Mb,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb
        const Vec3&     My_Mb = R_MbM.y(); // M's y axis, expressed in Mb

        const Vec3&     w_MbM = getV_MbM(vc)[0]; // angular velocity of M in Mb

        H_MbM_Dot[0] = SpatialRow( ~(w_MbM % Mx_Mb), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( ~(w_MbM % My_Mb), Row3(0) );
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const Vec3 w_MbM_M = fromU(u).append1(0); // angular velocity of M in Mb, exp in M (with wz=0) 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(fromQVec3(q,0),
                                        w_MbM_M); // need w in *body*, not parent
        } else {
            const Rotation& R_MbM = getX_MbM(pc).R();
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(q),
                                        R_MbM*w_MbM_M); // need w in *parent* frame
        }
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3 w_MbM_M     = fromU(u).append1(0); // angular velocity of M in Mb, exp in M (with wz=0)
        const Vec3 w_MbM_M_dot = fromU(udot).append1(0);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                       (fromQVec3(q,0), w_MbM_M, w_MbM_M_dot); // body frame
        } else {
            const Rotation& R_MbM = getX_MbM(pc).R();
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(q),R_MbM*w_MbM_M,R_MbM*w_MbM_M_dot); // parent frame
        }
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

    void setMobilizerDefaultPositionValues(
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
};


    // FREE LINE //

// FreeLine joint. Like a Free joint, this provides full rotational and
// translational freedom, but for a degenerate body which is thin (inertialess)
// along its own z axis. These arise in molecular modeling for linear molecules formed
// by pairs of atoms, or by multiple atoms in a linear arrangement like
// carbon dioxide (CO2) whose structure is O=C=O in a straight line. We are
// assuming that there is no meaning to a rotation about the linear axis,
// so free orientation requires just *two* degrees of freedom, not *three*
// as is required for general rigid bodies. And in fact we can get away with
// just two rotational generalized speeds so this joint provides only 5 mobilities.
// But so far, no one has been able to come up
// with a way to manage with only two rotational generalized *coordinates*, so this joint
// has the same q's as a regular Free joint: either a quaternion
// for unconditional stability, or a three-angle (body fixed 1-2-3)
// Euler sequence which will be dynamically singular when the middle (y) axis
// is 90 degrees. Use the Euler sequence only for small motions or for kinematics
// problems (and note that only the first two are meaningful). Translations here
// are treated exactly as for a Free joint (or for a Cartesian joint for that matter).
//
// To summarize, the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//   * 3 components of the translation vector T_MbM (that is, vector from origin
//     of Mb to origin of M, expressed in Mb)
// and generalized speeds are:
//   * the x,y components of the angular velocity w_MbM_M, that is, the angular
//     velocity of M in Mb expressed in *M* (where we want wz=0).
//   * 3 components of the linear velocity of origin of M in Mb, expressed in Mb.
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeFreeLine : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "full"; }

    RBNodeFreeLine(const MassProperties& mProps_B,
                   const Transform&      X_PMb,
                   const Transform&      X_BM,
                   int&                  nextUSlot,
                   int&                  nextUSqSlot,
                   int&                  nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PMb,X_BM,nextUSlot,nextUSqSlot,nextQSlot)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBModelVars& mv, const Rotation& R_MbM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(mv))
            toQVec3(q,0) = R_MbM.convertToBodyFixed123();
        else
            toQuat(q) = R_MbM.convertToQuaternion().asVec4();
    }

    // The user gives us the translation vector from OMb to OM as a vector expressed in Mb.
    // With a free joint we never have to *change* orientation coordinates in order to achieve a translation.
    // Note: a quaternion from a state is not necessarily normalized so can't be used
    // direction as though it were a set of Euler parameters; it must be normalized first.
    void setQToFitTranslation(const SBModelVars& mv, const Vec3& T_MbM, Vector& q, bool only) const {
        if (getUseEulerAngles(mv))
            toQVec3(q,3) = T_MbM; // skip the 3 Euler angles
        else
            toQVec3(q,4) = T_MbM; // skip the 4 quaternions
    }

    // Our 2 rotational generalized speeds are just the (x,y) components of the
    // angular velocity vector of M in Mb, expressed in M.
    void setUToFitAngularVelocity(const SBModelVars& mv, const Vector& q, const Vec3& w_MbM,
                                     Vector& u) const
    {
        Rotation R_MbM;
        if (getUseEulerAngles(mv))
            R_MbM.setToBodyFixed123(fromQVec3(q,0));
        else {
            // TODO: should use qnorm pool
            R_MbM.setToQuaternion(Quaternion(fromQuat(q))); // normalize
        }
        const Vec3 w_MbM_M = ~R_MbM*w_MbM;
        toU(u).updSubVec<2>(0) = Vec2(w_MbM_M[0], w_MbM_M[1]); // (x,y) of relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in Mb,
    // expressed in Mb. The user gives us that same vector.
    void setUToFitLinearVelocity
       (const SBModelVars& mv, const Vector& q, const Vec3& v_MbM, Vector& u, bool only) const
    {
        toUVec3(u,2) = v_MbM;
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
            X_MbM.updT() = fromQVec3(q,3); // translation is in Mb
        } else {
            X_MbM.updR().setToQuaternion(Quaternion(fromQuat(q))); // normalize
            X_MbM.updT() = fromQVec3(q,4);  // translation is in Mb
        }
    }


    // The generalized speeds for this 5-dof ("free line") joint are 
    //   (1) the (x,y) components of angular velocity of M in the Mb frame, expressed in M, and
    //   (2) the (linear) velocity of M's origin in Mb, expressed in Mb.
    void calcAcrossJointVelocityJacobian(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        HType&                 H_MbM) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb
        const Vec3&     My_Mb = R_MbM.y(); // M's y axis, expressed in Mb

        H_MbM[0] = SpatialRow( ~Mx_Mb, Row3(0) );        // x,y angular velocity in M, re-expressed im Mb
        H_MbM[1] = SpatialRow( ~My_Mb, Row3(0) );

        H_MbM[2] = SpatialRow( Row3(0), Row3(1,0,0) );   // translations in Mb
        H_MbM[3] = SpatialRow( Row3(0), Row3(0,1,0) );
        H_MbM[4] = SpatialRow( Row3(0), Row3(0,0,1) );
    }

    // Since the first two rows of the Jacobian above are not constant in Mb,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBModelVars&     mv,
        const SBPositionCache& pc, 
        const SBVelocityCache& vc, 
        HType&                 H_MbM_Dot) const
    {
        const Rotation& R_MbM = getX_MbM(pc).R();
        const Vec3&     Mx_Mb = R_MbM.x(); // M's x axis, expressed in Mb
        const Vec3&     My_Mb = R_MbM.y(); // M's y axis, expressed in Mb

        const Vec3&     w_MbM = getV_MbM(vc)[0]; // angular velocity of M in Mb

        H_MbM_Dot[0] = SpatialRow( ~(w_MbM % Mx_Mb), Row3(0) );
        H_MbM_Dot[1] = SpatialRow( ~(w_MbM % My_Mb), Row3(0) );

        // For translation in Mb.
        H_MbM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[3] = SpatialRow( Row3(0), Row3(0) );
        H_MbM_Dot[4] = SpatialRow( Row3(0), Row3(0) );
    }

    void calcQDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u,
        Vector&                qdot) const
    {
        const Vec3  w_MbM_M = Vec3(fromU(u)[0], fromU(u)[1], 0); // Angular velocity in M
        const Vec3& v_MbM   = fromUVec3(u,2);                    // Linear velocity in Mb

        if (getUseEulerAngles(mv)) {
            const Vec3& theta = fromQVec3(q,0); // Euler angles
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
                                            w_MbM_M); // need w in *body*, not parent
            toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdot,3) = v_MbM;
        } else {
            const Rotation& R_MbM = getX_MbM(pc).R();
            const Vec4& quat = fromQuat(q);
            toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,
                                            R_MbM*w_MbM_M); // need w in *parent* frame here
            toQVec3(qdot,4) = v_MbM;
        }
    }
 
    void calcQDotDot(
        const SBModelVars&     mv,
        const Vector&          q,
        const SBPositionCache& pc,
        const Vector&          u, 
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const Vec3  w_MbM_M     = Vec3(fromU(u)[0], fromU(u)[1], 0); // Angular velocity of M in Mb, exp. in M
        const Vec3& v_MbM       = fromUVec3(u,2); // linear velocity of M in Mb, expressed in M
        const Vec3  w_MbM_M_dot = Vec3(fromU(udot)[0], fromU(udot)[1], 0);
        const Vec3& v_MbM_dot   = fromUVec3(udot,2);

        if (getUseEulerAngles(mv)) {
            const Vec3& theta  = fromQVec3(q,0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, w_MbM_M, w_MbM_M_dot); // needed in body frame here
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_MbM_dot;
        } else {
            const Rotation& R_MbM = getX_MbM(pc).R();
            const Vec4& quat  = fromQuat(q);
            toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                             (quat,R_MbM*w_MbM_M,R_MbM*w_MbM_M_dot); // needed in parent frame
            toQVec3(qdotdot,4) = v_MbM_dot;
        }
    }

    void copyQ(Vector& q, const SBModelVars& mv, const Vector& qIn) const {
        if (getUseEulerAngles(mv)) {
            toQVec3(q,0) = fromQVec3(qIn,0); // euler angles
            toQVec3(q,3) = fromQVec3(qIn,3); // translations
        } else {
            toQuat(q)    = fromQuat(qIn);    // quaternion
            toQVec3(q,4) = fromQVec3(qIn,4); // translations
        }
    }

    int  getMaxNQ()                   const {return 7;}
    int  getNQ(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 
    bool isUsingQuaternion(const SBModelVars& mv) const {
        return !getUseEulerAngles(mv);
    }

    void setMobilizerDefaultPositionValues(const SBModelVars& mv, Vector& q) const 
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
};


/////////////////////////////////////////////////////////////
// Implementation of RigidBodyNodeSpec base class methods. //
/////////////////////////////////////////////////////////////


// Same for all mobilizers.
template<int dof> void
RigidBodyNodeSpec<dof>::calcParentToChildVelocityJacobianInGround(
    const SBModelVars&     mv,
    const SBPositionCache& pc, 
    HType& H_PB_G) const
{
    const HType& H_MbM = getH_MbM(pc);

    // want r_MB_Mb, that is, the vector from OM to OB, expressed in Mb
    const Vec3&     r_MB    = getX_MB().T();    // fixed
    const Rotation& R_MbM   = getX_MbM(pc).R(); // just calculated
    const Vec3      r_MB_Mb = R_MbM*r_MB;

    // Reminder: round brackets () applied to a matrix select columns.
    HType H_MB;
    H_MB(0) = Row3(0); // fills column with zero
    H_MB(1) = H_MbM(0) * crossMat(r_MB_Mb);

    // Now we want R_GMb so we can reexpress the cross-joint velocity V_MbB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PMb = getX_PMb().R();      // fixed config of Mb in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GMb = R_GP * R_PMb;

    H_PB_G = (H_MbM + H_MB) * ~R_GMb;
}

// Same for all mobilizers.
template<int dof> void
RigidBodyNodeSpec<dof>::calcParentToChildVelocityJacobianInGroundDot(
    const SBModelVars&     mv,
    const SBPositionCache& pc, 
    const SBVelocityCache& vc, 
    const SBDynamicsCache& dc,
    HType& H_PB_G_Dot) const
{
    const HType& H_MbM = getH_MbM(pc);
    const HType& H_MbM_Dot = getH_MbM_Dot(dc);

    HType H_MB, H_MB_Dot;

    // want r_MB_Mb, that is, the vector from OM to OB, expressed in Mb
    const Vec3&     r_MB    = getX_MB().T();    // fixed
    const Rotation& R_MbM   = getX_MbM(pc).R(); // just calculated
    const Vec3      r_MB_Mb = R_MbM*r_MB;

    const Vec3& w_MbM = getV_MbM(vc)[0]; // local angular velocity

    // Reminder: round brackets () applied to a matrix select columns.
    H_MB(0) = Row3(0); // fills column with zero
    H_MB(1) = H_MbM(0) * crossMat(r_MB_Mb);

    H_MB_Dot(0) = Row3(0);
    H_MB_Dot(1) =   H_MbM_Dot(0) * crossMat(r_MB_Mb) 
                  + H_MbM(0)     * crossMat(w_MbM % r_MB_Mb);

    // Now we want R_GMb so we can reexpress the cross-joint velocity V_MbB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PMb = getX_PMb().R();      // fixed config of Mb in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GMb = R_GP * R_PMb;

    const Vec3& w_GMb = getV_GP(vc)[0]; // Mb and P have same angular velocity

    // Note: time derivative of R_GMb is crossMat(w_GMb)*R_GMb, so derivative
    // of ~R_GMb is ~(crossMat(w_GMb)*R_GMb) = - ~R_GMb*crossMat(w_GMb) since
    // crossMat's are skew-symmetric.
    //      H_PB_G = (H_MbM + H_MB) * ~R_GMb (see above method)
    const HType& H_PB_G = getH(pc);
    H_PB_G_Dot =  (H_MbM_Dot + H_MB_Dot) * ~R_GMb
                 - H_PB_G * crossMat(w_GMb);
}

//
// to be called from base to tip.
//
template<int dof> void
RigidBodyNodeSpec<dof>::setVelFromSVel(
    const SBPositionCache& pc, 
    const SBVelocityCache& mc,
    const SpatialVec&      sVel, 
    Vector&                u) const 
{
    toU(u) = getH(pc) * (sVel - (~getPhi(pc) * parent->getV_GB(mc)));
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
    const SBPositionCache& pc,
    SBDynamicsCache&       dc) const 
{
    updP(dc) = getMk(pc);
    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialMat& tauBarChild = children[i]->getTauBar(dc);
        const SpatialMat& PChild      = children[i]->getP(dc);
        const PhiMatrix&  phiChild    = children[i]->getPhi(pc);

        // TODO: this is around 450 flops but could be cut in half by
        // exploiting symmetry.
        updP(dc) += phiChild * (tauBarChild * PChild) * ~phiChild;
    }

    const Mat<2,dof,Vec3> PHt = getP(dc) * ~getH(pc);
    updD(dc)  = getH(pc) * PHt;
    // this will throw an exception if the matrix is ill conditioned
    updDI(dc) = getD(dc).invert();
    updG(dc)  = PHt * getDI(dc);

    // TODO: change sign on tau to make it GH-I instead, which only requires
    // subtractions on the diagonal rather than negating all the off-diag stuff.
    // That would save 30 flops here (I know, not much).
    updTauBar(dc)  = 1.; // identity matrix
    updTauBar(dc) -= getG(dc) * getH(pc);
    updPsi(dc)     = getPhi(pc) * getTauBar(dc);
}



// To be called base to tip.
// sherm 060723: As best I can tell this is calculating the inverse of
// the "operational space inertia" at the body frame origin for each body.
// See Equation 20 in Rodriguez,Jain, & Kreutz-Delgado: A spatial operator algebra 
// for manipulator modeling and control. Intl. J. Robotics Research 
// 10(4):371-381 (1991).
template<int dof> void
RigidBodyNodeSpec<dof>::calcYOutward(
    const SBPositionCache& pc,
    SBDynamicsCache&       dc) const 
{
    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?
    updY(dc) = (~getH(pc) * getDI(dc) * getH(pc)) 
                + (~getPsi(dc) * parent->getY(dc) * getPsi(dc));
}

//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(
    const SBPositionCache&    pc,
    const SBDynamicsCache&    dc,
    const SpatialVec&         spatialForce,
    const SBAccelerationVars& av,
    SBAccelerationCache&      ac) const 
{
    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - spatialForce;

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(ac);
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& GepsChild = children[i]->getGepsilon(ac);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(ac)  = getAppliedJointForce(av) - getH(pc)*z; // TODO: pass in hinge forces
    updNu(ac)       = getDI(dc) * getEpsilon(ac);
    updGepsilon(ac) = getG(dc)  * getEpsilon(ac);
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
    const SBPositionCache& pc,
    const Vector&          allU,
    const SBDynamicsCache& dc,
    SBAccelerationCache&   ac,
    Vector&                allUdot,
    Vector&                allQdotdot) const 
{
    Vec<dof>&        udot   = toU(allUdot);
    const SpatialVec alphap = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

    udot        = getNu(ac) - (~getG(dc)*alphap);
    updA_GB(ac) = alphap + ~getH(pc)*udot + getCoriolisAcceleration(dc);  

    calcQDotDot(mv, allQ, pc, allU, allUdot, allQdotdot);  
}

 
//
// To be called from tip to base.
// Temps do not need to be initialized.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcUDotPass1Inward(
    const SBPositionCache&      pc,
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
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - getH(pc)*z;
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
    const SBPositionCache& pc,
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
        : ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + ~getH(pc)*udot + getCoriolisAcceleration(dc);  
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
    const SBPositionCache& pc,
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
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - getH(pc)*z;
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
    const SBPositionCache& pc,
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
        : ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + ~getH(pc)*udot;  
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
    const SBPositionCache&      pc,
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
        const PhiMatrix&  phiChild = children[i]->getPhi(pc);

        z += phiChild * zChild;
    }

    out = getH(pc) * z; 
}

//
// To be called from tip to base.
// Temps do not need to be initialized.
// (sherm 060727) In spatial operators, this calculates H*Phi*(F-(Pa+b))
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcEquivalentJointForces(
    const SBPositionCache&     pc,
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
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];

        z += phiChild * zChild; 
    }

    eps  = getH(pc) * z;
}


// The Ground node is special because it doesn't need a mobilizer.
/*static*/ RigidBodyNode*
RigidBodyNode::createGroundNode() {
    return new RBGroundBody();
}


    //////////////////////////////////////////////////////////////////////
    // Implementation of MobilizedBodyRep createRigidBodyNode() methods //
    //////////////////////////////////////////////////////////////////////

RigidBodyNode* MobilizedBody::Pin::PinRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeTorsion(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Slider::SliderRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeSlider(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Universal::UniversalRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeUJoint(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Cylinder::CylinderRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeCylinder(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::BendStretch::BendStretchRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeBendStretch(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Planar::PlanarRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodePlanar(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Gimbal::GimbalRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    assert(!"Gimbal MobilizedBody not implemented yet"); return 0;
    // return new RBNodeGimbal(getDefaultRigidBodyMassProperties(),
    //     getDefaultInboardFrame(),getDefaultOutboardFrame(),
    //     nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Ball::BallRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeRotate3(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Translation::TranslationRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeTranslate(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Free::FreeRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeTranslateRotate3(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::LineOrientation::LineOrientationRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeLineOrientation(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::FreeLine::FreeLineRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeFreeLine(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nxtUSlot,nxtUSqSlot,nxtQSlot);
}


RigidBodyNode* MobilizedBody::Screw::ScrewRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBNodeScrew(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        getDefaultPitch(),nxtUSlot,nxtUSqSlot,nxtQSlot);
}

RigidBodyNode* MobilizedBody::Weld::WeldRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    assert(!"Weld MobilizedBody not implemented yet"); return 0;
    // return new RBNodeWeld(
    //     getDefaultRigidBodyMassProperties(),
    //     getDefaultInboardFrame(),getDefaultOutboardFrame(),
    //     nxtUSlot,nxtUSqSlot,nxtQSlot);
}


RigidBodyNode* MobilizedBody::Ground::GroundRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    return new RBGroundBody();
}

RigidBodyNode* MobilizedBody::Custom::CustomRep::createRigidBodyNode(
    int&                     nxtUSlot,
    int&                     nxtUSqSlot,
    int&                     nxtQSlot) const
{
    assert(!"Custom MobilizedBody not implemented yet"); return 0;
    // return new RBNodeCustom(
    //     getDefaultRigidBodyMassProperties(),
    //     getDefaultInboardFrame(),getDefaultOutboardFrame(),
    //     nxtUSlot,nxtUSqSlot,nxtQSlot);
}



