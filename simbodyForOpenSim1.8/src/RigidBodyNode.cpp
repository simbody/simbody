/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *    Charles Schwieters (NIH): wrote the public domain IVM code from which   *
 *                              this was derived.                             *
 *    Ajay Seth: wrote the Ellipsoid joint.                                   *
 *    Peter Eastman: wrote the Custom mobilizer and Euler Angle<->Quaternion  *
 *                   conversion                                               *
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
 * This file contains all the multibody mechanics code that involves a single body and
 * its mobilizer (inboard joint), that is, one node in the multibody tree. These methods
 * constitute the inner loops of the multibody calculations, and much suffering is
 * undergone here to make them run fast. In particular most calculations are templatized
 * by the number of mobilities, so that compile-time sizes are known for everything.
 *
 * Most methods here expect to be called in a particular order during traversal of the
 * tree -- either base to tip or tip to base.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "MobilizedBodyImpl.h"

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


////////////////////////////////////////////////
// Define classes derived from RigidBodyNode. //
////////////////////////////////////////////////

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class RBGroundBody : public RigidBodyNode {
public:
    RBGroundBody(const MassProperties& mProps_B, const Transform& X_PF, const Transform& X_BM) : 
        RigidBodyNode(mProps_B, X_PF, X_BM, QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed) 
    {
        uIndex   = UIndex(0);
        uSqIndex = USquaredIndex(0);
        qIndex   = QIndex(0);
    }
    ~RBGroundBody() {}

    /*virtual*/const char* type() const { return "ground"; }
    /*virtual*/int  getDOF()   const {return 0;}
    /*virtual*/int  getMaxNQ() const {return 0;}
    /*virtual*/int  getNUInUse(const SBModelVars&) const {return 0;}
    /*virtual*/int  getNQInUse(const SBModelVars&) const {return 0;}
    /*virtual*/bool isUsingQuaternion(const SBStateDigest&, MobilizerQIndex& ix) const {
        ix.invalidate();
        return false;
    }
    /*virtual*/bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& ix, int& nAngles) const {
        ix.invalidate(); nAngles = 0;
        return false;
    }
    /*virtual*/bool enforceQuaternionConstraints(
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const {return false;}

    /*virtual*/void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {}
    /*virtual*/void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {}

    /*virtual*/void setMobilizerDefaultModelValues(const SBTopologyCache&, 
                                          SBModelVars& v) const
    {
        v.prescribed[0] = true; // ground's motion is prescribed to zero
    }

    /*virtual*/ void setQToFitTransform
       (const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {}
    /*virtual*/ void setQToFitRotation
       (const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {}
    /*virtual*/ void setQToFitTranslation
       (const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {}

    /*virtual*/ void setUToFitVelocity
       (const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {}
    /*virtual*/ void setUToFitAngularVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {}
    /*virtual*/ void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const {}


    /*virtual*/void realizeModel(SBStateDigest& sbs) const {}

    /*virtual*/void realizeInstance(SBStateDigest& sbs) const {}

    /*virtual*/void realizeTime(SBStateDigest& sbs) const {}

    /*virtual*/void realizePosition(SBStateDigest& sbs) const {}

    /*virtual*/void realizeVelocity(SBStateDigest& sbs) const {}

    /*virtual*/ void realizeDynamics(SBStateDigest& sbs) const {}

    /*virtual*/ void realizeAcceleration(SBStateDigest& sbs) const {}

    /*virtual*/ void realizeReport(SBStateDigest& sbs) const {}

    /*virtual*/void calcArticulatedBodyInertiasInward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const {}


    /*virtual* void calcZ(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        SBAccelerationCache&) const {} */


    /*virtual*/void calcZ(
        const SBStateDigest&,
        const SBDynamicsCache&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces) const {} 

    /*virtual*/void calcYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const {}

    /*virtual*/void calcAccel(
            const SBStateDigest&   sbs,
            Vector&                udot,
            Vector&                qdotdot) const {}

    /*virtual*/ void calcSpatialKinematicsFromInternal(
        const SBPositionCache&      pc,
        const Vector&               v,
        Vector_<SpatialVec>&        Jv) const    
    {
        Jv[0] = SpatialVec(Vec3(0), Vec3(0));
    }

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

    /*virtual*/void setVelFromSVel(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        const SpatialVec&      sVel, 
        Vector&                u) const {}
    
    /*virtual*/ void multiplyByQBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                  bool matrixOnRight, 
                                  const Real* in, Real* out) const {}
    /*virtual*/ void multiplyByQInvBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                     bool matrixOnRight,
                                     const Real* in, Real* out) const {}
    /*virtual*/ void multiplyByQDotBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q, const Real* u,
                                     bool matrixOnRight,
                                     const Real* in, Real* out) const {}


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
                      const Transform&      X_PF,
                      const Transform&      X_BM,
                      UIndex&               nextUSlot,
                      USquaredIndex&        nextUSqSlot,
                      QIndex&               nextQSlot,
                      QDotHandling          qdotHandling,
                      QuaternionUse         quaternionUse)
      : RigidBodyNode(mProps_B, X_PF, X_BM, qdotHandling, quaternionUse)
    {
        // don't call any virtual methods in here!
        uIndex   = nextUSlot;
        uSqIndex = nextUSqSlot;
        qIndex   = nextQSlot;
    }

    void updateSlots(UIndex& nextUSlot, USquaredIndex& nextUSqSlot, QIndex& nextQSlot) {
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

    virtual void setQToFitTransform(const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {
        setQToFitRotation   (sbs,X_FM.R(),q);
        setQToFitTranslation(sbs,X_FM.T(),q);
    }

    virtual void setUToFitVelocity(const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {
        setUToFitAngularVelocity(sbs,q,V_FM[0],u);
        setUToFitLinearVelocity (sbs,q,V_FM[1],u);
    }

    // The following routines calculate joint-specific position kinematic
    // quantities. They may assume that *all* position kinematics (not just
    // joint-specific) has been done for the parent, and that the position
    // state variables q are available. Each routine may also assume that the
    // previous routines have been called, in the order below.
    // The routines are structured as operators -- they use the State but
    // do not change anything in it, including the cache. Instead they are
    // passed arguments to write their results into. In practice, these
    // arguments will typically be in the State cache (see below).

    // This mandatory routine performs expensive floating point operations sin,cos,sqrt
    // in one place so we don't end up repeating them. sin&cos are used only
    // for mobilizers which have angular coordinates, and qErr and qnorm are only for
    // mobilizers using quaternions. Other mobilizers can provide a null routine.
    // Each of the passed-in Vectors is a "q-like" object, that is, allocated
    // to the bodies in a manner parallel to the q state variable, except that qErr
    // has just one slot per quaternion and must be accessed using the node's 
    // quaternionIndex which is in the Model cache.
    virtual void calcJointSinCosQNorm(
        const SBModelVars&  mv, 
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const=0;

    // This mandatory routine calculates the across-joint transform X_FM generated
    // by the current q values. This may depend on sines & cosines or normalized
    // quaternions already being available in the State cache.
    // This method constitutes the *definition* of the generalized coordinates for
    // a particular joint.
    virtual void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const=0;


    // This mandatory routine calculates the joint transition matrix H_FM, giving
    // the change of velocity induced by the generalized speeds u for this 
    // mobilizer, expressed in the mobilizer's inboard "fixed" frame F (attached to
    // this body's parent). 
    // This method constitutes the *definition* of the generalized speeds for
    // a particular joint.
    // This routine can depend on X_FM having already
    // been calculated and available in the PositionCache but must NOT depend
    // on any quantities involving Ground or other bodies.
    virtual void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const=0;

    // This mandatory routine calculates the time derivative taken in F of
    // the matrix H_FM (call it H_FM_Dot). This is zero if the generalized
    // speeds are all defined in terms of the F frame, which is often the case.
    // This routine can depend on X_FM and H_FM being available already in
    // the state PositionCache, and V_FM being already in the state VelocityCache.
    // However, it must NOT depend on any quantities involving Ground or other bodies.
    virtual void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const = 0;

    // This routine is NOT joint specific, but cannot be called until the across-joint
    // transform X_FM has been calculated and is available in the State cache.
    void calcBodyTransforms(
        const SBPositionCache& pc, 
        Transform&             X_PB, 
        Transform&             X_GB) const 
    {
        const Transform& X_BM = getX_BM();   // fixed
        const Transform& X_PF = getX_PF();   // fixed
        const Transform& X_FM = getX_FM(pc); // just calculated
        const Transform& X_GP = getX_GP(pc); // already calculated

        X_PB = X_PF * X_FM * ~X_BM; // TODO: precalculate X_MB
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


    // Calculate joint-specific kinematic quantities dependent on
    // velocities. This routine may assume that *all* position 
    // kinematics (not just joint-specific) has been done for this node,
    // that all velocity kinematics has been done for the parent, and
    // that the velocity state variables (u) are available. The
    // quanitites that must be computed are:
    //   V_FM   relative velocity of B's M frame in P's F frame, 
    //             expressed in F (note: this is also V_PM_F since
    //             F is fixed on P).
    //   V_PB_G  relative velocity of B in P, expr. in G
    // The code is the same for all joints, although parametrized by ndof.
    void calcJointKinematicsVel(
        const SBPositionCache& pc,
        const Vector&          u,
        SBVelocityCache&       vc) const 
    {
        updV_FM(vc)  = ~getH_FM(pc) * fromU(u);
        updV_PB_G(vc) = ~getH(pc) * fromU(u);
    }

    // Calculate joint-specific dynamics quantities dependent on velocities.
    // This method may assume that *all* position & velocity kinematics
    // (not just joint-specific) has been done for this node, and that
    // dynamics has been done for the parent.
    void calcJointDynamics(
        const SBPositionCache& pc,
        const Vector&          u,
        const SBVelocityCache& vc, 
        SBDynamicsCache&       dc) const 
    {
        updVD_PB_G(dc) = ~getHDot(dc) * fromU(u);
    }

    // These next two routines are optional, but if you supply one you
    // must supply the other. (That is, ball-containing joints provide
    // both of these routines.)
    virtual void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u,
        Vector&                qdot) const
    {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        toQ(qdot) = fromU(u);        // default is qdot=u
    }

    virtual void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const
    {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        toQ(qdotdot) = fromU(udot);  // default is qdotdot=udot
    }

    void realizeModel(SBStateDigest& sbs) const 
    {
    }

    void realizeInstance(SBStateDigest& sbs) const
    {
    }

    void realizeTime(SBStateDigest& sbs) const
    {
    }

    // Set a new configuration and calculate the consequent kinematics.
    // Must call base-to-tip.
    void realizePosition(SBStateDigest& sbs) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBModelCache& mc = sbs.getModelCache();
        const SBInstanceCache& ic = sbs.getInstanceCache();
        SBPositionCache& pc = sbs.updPositionCache();
        calcJointSinCosQNorm(mv, mc, ic, sbs.getQ(), pc.sq, pc.cq, sbs.updQErr(), pc.qnorm);

        calcAcrossJointTransform (sbs, sbs.getQ(), updX_FM(pc));
        calcBodyTransforms       (pc, updX_PB(pc), updX_GB(pc));

        calcAcrossJointVelocityJacobian(sbs, updH_FM(pc));
        calcParentToChildVelocityJacobianInGround(mv,pc, updH(pc));

        calcJointIndependentKinematicsPos(pc);
    }

    // Set new velocities for the current configuration, and calculate
    // all the velocity-dependent terms. Must call base-to-tip.
    void realizeVelocity(SBStateDigest& sbs) const 
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        SBVelocityCache& vc = sbs.updVelocityCache();
        calcQDot(sbs, sbs.getU(), sbs.updQDot());
        calcJointKinematicsVel(pc,sbs.getU(),vc);
        calcJointIndependentKinematicsVel(pc,vc);
    }

    void realizeDynamics(SBStateDigest& sbs) const
    {
        // Mobilizer-specific.
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        SBDynamicsCache& dc = sbs.updDynamicsCache();
        calcAcrossJointVelocityJacobianDot          (sbs, updH_FM_Dot(dc));
        calcParentToChildVelocityJacobianInGroundDot(mv,pc,vc, dc, updHDot(dc));

        calcJointDynamics(pc,sbs.getU(),vc,dc);

        // Mobilizer independent.
        calcJointIndependentDynamicsVel(pc,vc,dc);
    }

    void realizeAcceleration(SBStateDigest& sbs) const
    {
    }

    void realizeReport(SBStateDigest& sbs) const
    {
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
        assert(quaternionUse == QuaternionIsNeverUsed);
        toQ(q) = fromQ(qIn);
    }

    // Not virtual -- the number of u's is always the template argument dof.
    void copyU(
        const SBModelVars& mv, 
        const Vector&      uIn, 
        Vector&            u) const
    {
        toU(u) = fromU(uIn);
    }

    int          getDOF()            const {return dof;}
    virtual int  getMaxNQ()          const {
        assert(quaternionUse == QuaternionIsNeverUsed);
        return dof; // maxNQ can be larger than dof if there's a quaternion
    }
    virtual int  getNQInUse(const SBModelVars&) const {
        assert(quaternionUse == QuaternionIsNeverUsed); // method must be overridden otherwise
        return dof; // DOF <= NQ <= maxNQ
    }
    virtual int  getNUInUse(const SBModelVars&) const {
        // Currently NU is always just the Mobilizer's DOFs (the template argument)
        // Later we may want to offer modeling options to lock joints, or perhaps
        // break them.
        return dof;
    }
    virtual bool isUsingQuaternion(const SBStateDigest&, MobilizerQIndex& startOfQuaternion) const {
        assert(quaternionUse == QuaternionIsNeverUsed); // method must be overridden otherwise
        startOfQuaternion.invalidate();
        return false;
    }

    // Most mobilizers do use angles, so we're not going to provide a default implementation
    // of the pure virtual isUsingAngles() method here.

    // State digest should be at Stage::Position.
    virtual void calcLocalQDotFromLocalU(const SBStateDigest&, const Real* u, Real* qdot) const {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        Vec<dof>::updAs(qdot) = Vec<dof>::getAs(u); // default says qdot=u
    }

    // State digest should be at Stage::Velocity.
    virtual void calcLocalQDotDotFromLocalUDot(const SBStateDigest&, const Real* udot, Real* qdotdot) const {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        Vec<dof>::updAs(qdotdot) = Vec<dof>::getAs(udot); // default says qdotdot=udot
    }

    // State digest should be at Stage::Position for calculating Q (the matrix which maps
    // generalized speeds u to coordinate derivatives qdot, qdot=Qu). The default implementation
    // assumes nq==nu and the nqXnu block of Q corresponding to this mobilizer is identity. Then 
    // either operation (regardless of side) just copies nu numbers from in to out.
    //
    // THIS MUST BE OVERRIDDEN by any mobilizer for which nq != nu, or qdot != u.
    virtual void multiplyByQBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                  bool matrixOnRight,  const Real* in, Real* out) const
    {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        Vec<dof>::updAs(out) = Vec<dof>::getAs(in);
    }

    virtual void multiplyByQInvBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                     bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        Vec<dof>::updAs(out) = Vec<dof>::getAs(in);
    }

    // State digest should be at Stage::Velocity for calculating QDot (the matrix which is used in
    // mapping generalized accelerations udot to coordinate 2nd derivatives qdotdot, 
    // qdotdot=Q udot + QDot u. The default implementation assumes nq==nu and the nqXnu block
    // of Q corresponding to this mobilizer is identity, so QDot is an nuXnu block of zeroes.
    // Then either operation (regardless of side) just copies nu zeros to out.
    //
    // THIS MUST BE OVERRIDDEN by any mobilizer for which nq != nu, or qdot != u.
    /* TODO
    virtual void multiplyByQDotBlock(const SBStateDigest&, bool useEulerAnglesIfPossible, 
                                     const Real* q, const Real* u,
                                     bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(qdotHandling == QDotIsAlwaysTheSameAsU);
        Vec<dof>::updAs(out) = 0;
    }
    */

    // No default implementations here for:
    // calcMobilizerTransformFromQ
    // calcMobilizerVelocityFromU
    // calcMobilizerAccelerationFromUDot
    // calcParentToChildTransformFromQ
    // calcParentToChildVelocityFromU
    // calcParentToChildAccelerationFromUDot


    virtual void setVelFromSVel(
        const SBPositionCache& pc, 
        const SBVelocityCache& mc,
        const SpatialVec&      sVel, 
        Vector&               u) const;

    // Return true if any change is made to the q vector.
    virtual bool enforceQuaternionConstraints(
        const SBStateDigest& sbs,
        Vector&            q,
        Vector&            qErrest) const 
    {
        assert(quaternionUse == QuaternionIsNeverUsed);
        return false;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        // The default implementation just copies Q.  Subclasses may override this.
        assert(quaternionUse == QuaternionIsNeverUsed);
        toQ(outputQ) = fromQ(inputQ);
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        // The default implementation just copies Q.  Subclasses may override this.
        assert(quaternionUse == QuaternionIsNeverUsed);
        toQ(outputQ) = fromQ(inputQ);
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
    const Vec<dof>&     fromQ  (const Real* q)   const {return Vec<dof>::getAs(&q[qIndex]);}
    Vec<dof>&           toQ    (      Real* q)   const {return Vec<dof>::updAs(&q[qIndex]);}
    const Vec<dof>&     fromU  (const Real* u)   const {return Vec<dof>::getAs(&u[uIndex]);}
    Vec<dof>&           toU    (      Real* u)   const {return Vec<dof>::updAs(&u[uIndex]);}
    const Mat<dof,dof>& fromUSq(const Real* uSq) const {return Mat<dof,dof>::getAs(&uSq[uSqIndex]);}
    Mat<dof,dof>&       toUSq  (      Real* uSq) const {return Mat<dof,dof>::updAs(&uSq[uSqIndex]);}

    // Same, but specialized for the common case where dof=1 so everything is scalar.
    const Real& from1Q  (const Real* q)   const {return q[qIndex];}
    Real&       to1Q    (      Real* q)   const {return q[qIndex];}
    const Real& from1U  (const Real* u)   const {return u[uIndex];}
    Real&       to1U    (      Real* u)   const {return u[uIndex];}
    const Real& from1USq(const Real* uSq) const {return uSq[uSqIndex];}
    Real&       to1USq  (      Real* uSq) const {return uSq[uSqIndex];}

    // Same, specialized for quaternions. We're assuming that the quaternion is stored
    // in the *first* four coordinates, if there are more than four altogether.
    const Vec4& fromQuat(const Real* q) const {return Vec4::getAs(&q[qIndex]);}
    Vec4&       toQuat  (      Real* q) const {return Vec4::updAs(&q[qIndex]);}

    // Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
    const Vec3& fromQVec3(const Real* q, int offs) const {return Vec3::getAs(&q[qIndex+offs]);}
    Vec3&       toQVec3  (      Real* q, int offs) const {return Vec3::updAs(&q[qIndex+offs]);}
    const Vec3& fromUVec3(const Real* u, int offs) const {return Vec3::getAs(&u[uIndex+offs]);}
    Vec3&       toUVec3  (      Real* u, int offs) const {return Vec3::updAs(&u[uIndex+offs]);}

    // These provide an identical interface for when q,u are given as Vectors rather than Reals.

    const Vec<dof>&     fromQ  (const Vector& q)   const {return fromQ(&q[0]);} // convert to array of Real
    Vec<dof>&           toQ    (      Vector& q)   const {return toQ  (&q[0]);}
    const Vec<dof>&     fromU  (const Vector& u)   const {return fromU(&u[0]);}
    Vec<dof>&           toU    (      Vector& u)   const {return toU  (&u[0]);}
    const Mat<dof,dof>& fromUSq(const Vector& uSq) const {return fromUSq(&uSq[0]);}
    Mat<dof,dof>&       toUSq  (      Vector& uSq) const {return toUSq  (&uSq[0]);}

    const Real& from1Q  (const Vector& q)   const {return from1Q  (&q[0]);}
    Real&       to1Q    (      Vector& q)   const {return to1Q    (&q[0]);}
    const Real& from1U  (const Vector& u)   const {return from1U  (&u[0]);}
    Real&       to1U    (      Vector& u)   const {return to1U    (&u[0]);}
    const Real& from1USq(const Vector& uSq) const {return from1USq(&uSq[0]);}
    Real&       to1USq  (      Vector& uSq) const {return to1USq  (&uSq[0]);}

    // Same, specialized for quaternions. We're assuming that the quaternion comes first in the coordinates.
    const Vec4& fromQuat(const Vector& q) const {return fromQuat(&q[0]);}
    Vec4&       toQuat  (      Vector& q) const {return toQuat  (&q[0]);}

    // Extract a Vec3 from a Q-like or U-like object, beginning at an offset from the qIndex or uIndex.
    const Vec3& fromQVec3(const Vector& q, int offs) const {return fromQVec3(&q[0], offs);}
    Vec3&       toQVec3  (      Vector& q, int offs) const {return toQVec3  (&q[0], offs);}
    const Vec3& fromUVec3(const Vector& u, int offs) const {return fromUVec3(&u[0], offs);}
    Vec3&       toUVec3  (      Vector& u, int offs) const {return toUVec3  (&u[0], offs);}

    // Applications of the above extraction routines to particular interesting items in the State. Note
    // that you can't use these for quaternions since they extract "dof" items.

    // Cache entries (cache is mutable in a const State)

        // Position


    // TODO: should store as H_FM or else always reference Ht_FM
    const Mat<dof,2,Row3,1,2>& getH_FM(const SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::getAs(&pc.storageForHtFM(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH_FM(SBPositionCache& pc) const
      { return ~Mat<2,dof,Vec3>::updAs(&pc.storageForHtFM(0,uIndex)); }

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

    const Mat<dof,2,Row3,1,2>& getH_FM_Dot(const SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::getAs(&dc.storageForHtFMDot(0,uIndex)); }
    Mat<dof,2,Row3,1,2>&       updH_FM_Dot(SBDynamicsCache& dc) const
      { return ~Mat<2,dof,Vec3>::updAs(&dc.storageForHtFMDot(0,uIndex)); }

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

    /*void calcZ(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces,
        SBAccelerationCache& ) const;*/
    void calcZ(
        const SBStateDigest&,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces) const;

    void calcAccel(
        const SBStateDigest&   sbs,
        Vector&                udot,
        Vector&                qdotdot) const;

    void calcSpatialKinematicsFromInternal(
        const SBPositionCache&      pc,
        const Vector&               v,
        Vector_<SpatialVec>&        Jv) const;

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

};

    //////////////////////////////////////////
    // Derived classes for each joint type. //
    //////////////////////////////////////////


    // TRANSLATION (CARTESIAN) //

// Translate (Cartesian) joint. This provides three degrees of
// translational freedom which is suitable (e.g.) for connecting a
// free atom to ground. The Cartesian directions are the axes of
// the parent body's F frame, with M=F when all 3 coords are 0,
// and the orientation of M in F is 0 (identity) forever.
class RBNodeTranslate : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "translate"; }

    RBNodeTranslate(const MassProperties& mProps_B,
                    const Transform&      X_PF,
                    const Transform&      X_BM,
                    UIndex&               nextUSlot,
                    USquaredIndex&        nextUSqSlot,
                    QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // the only rotation this mobilizer can represent is identity
    }
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3&  T_FM, Vector& q) const {
        // here's what this joint is really good at!
        toQ(q) = T_FM;
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // The only angular velocity this can represent is zero.
    }
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // linear velocity is in a Cartesian joint's sweet spot
        toU(u) = v_FM;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        startOfAngles.invalidate(); nAngles=0; // no angles for a Cartesian mobilizer
        return false;
    }

    // This is required but does nothing here since there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const { }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        // Translation vector q is expressed in F (and M since they have same orientation).
        // A Cartesian joint can't change orientation. 
        X_FM = Transform(Rotation(), fromQ(q));
    }

    // Generalized speeds together are the velocity of M's origin in the F frame,
    // expressed in F. So individually they produce velocity along F's x,y,z
    // axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0), Row3(1,0,0) );
        H_FM[1] = SpatialRow( Row3(0), Row3(0,1,0) );
        H_FM[2] = SpatialRow( Row3(0), Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in F, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

};



    // SLIDING (PRISMATIC) //

// Sliding joint (1 dof translation). The translation is along the x
// axis of the parent body's F frame, with M=F when the coordinate
// is zero and the orientation of M in F frozen at 0 forever.
class RBNodeSlider : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "slider"; }

    RBNodeSlider(const MassProperties& mProps_B,
                 const Transform&      X_PF,
                 const Transform&      X_BM,
                 UIndex&               nextUSlot,
                 USquaredIndex&        nextUSqSlot,
                 QIndex&               nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }
        // Implementations of virtual methods.

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation a slider can represent is identity.
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // We can only represent the x coordinate with this joint.
        to1Q(q) = T_FM[0];
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // The only angular velocity a slider can represent is zero.
    }

    void setUToFitLinearVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // We can only represent a velocity along x with this joint.
        to1U(u) = v_FM[0];
    }


    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        startOfAngles.invalidate(); nAngles=0; // no angles for a Slider
        return false;
    }

    // This is required but does nothing here since we there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const { }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&      q,
        Transform&         X_FM) const
    {
        // Translation vector q is expressed in F (and M since they have same orientation).
        // A sliding joint can't change orientation, and only translates along x. 
        X_FM = Transform(Rotation(), Vec3(from1Q(q),0,0));
    }

    // The generalized speed is the velocity of M's origin in the F frame,
    // along F's x axis, expressed in F.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0), Row3(1,0,0) );
    }

    // Since the Jacobian above is constant in F, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // PIN (TORSION) //

// This is a "pin" or "torsion" joint, meaning one degree of rotational freedom
// about a particular axis, the z axis of the parent's F frame, which is 
// aligned forever with the z axis of the body's M frame. In addition, the
// origin points of M and F are identical forever.
class RBNodeTorsion : public RigidBodyNodeSpec<1> {
public:
    virtual const char* type() { return "torsion"; }

    RBNodeTorsion(const MassProperties& mProps_B,
                  const Transform&      X_PF,
                  const Transform&      X_BM,
                  UIndex&               nextUSlot,
                  USquaredIndex&        nextUSqSlot,
                  QIndex&               nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our pin joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        to1Q(q) = angles123[2];
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_FM[2]; // project angular velocity onto z axis
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }


    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        startOfAngles = MobilizerQIndex(0); nAngles=1; // torsion mobilizer
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Real& theta  = from1Q(q);    // angular coordinate

        // We're only updating the orientation here because a torsion joint
        // can't translate (it is defined as a rotation about the z axis).
        X_FM.updR().setRotationFromAngleAboutZ(theta);
        X_FM.updT() = 0.;
    }


    // The generalized speed is the angular velocity of M in the F frame,
    // about F's z axis, expressed in F. (This axis is also constant in M.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0,0,1), Row3(0) );
    }

    // Since the Jacobian above is constant in F, its time derivative in F is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};


    // SCREW //

// This is a one-dof "screw" joint, meaning one degree of rotational freedom
// about a particular axis, coupled to translation along that same axis.
// Here we use the common z axis of the F and M frames, which remains
// aligned forever. 
// For the generalized coordinate q, we use the rotation angle. For the
// generalized speed u we use the rotation rate, which is also the
// angular velocity of M in F (about the z axis). We compute the
// translational position as pitch*q, and the translation rate as pitch*u.
class RBNodeScrew : public RigidBodyNodeSpec<1> {
    Real pitch;
public:
    virtual const char* type() { return "screw"; }

    RBNodeScrew(const MassProperties& mProps_B,
                const Transform&      X_PF,
                const Transform&      X_BM,
                Real                  p,  // the pitch
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
      : RigidBodyNodeSpec<1>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed),
        pitch(p)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our screw joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        to1Q(q) = angles123[2];
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        to1Q(q) = T_FM[2]/pitch;
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        to1U(u) = w_FM[2]; // project angular velocity onto z axis
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        to1U(u) = v_FM[2]/pitch;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // We're currently using an angle as the generalized coordinate for the screw joint
        // but could just as easily have used translation or some non-physical coordinate. It
        // might make sense to offer a Model stage option to set the coordinate meaning.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Real& theta  = from1Q(q);    // angular coordinate

        X_FM.updR().setRotationFromAngleAboutZ(theta);
        X_FM.updT() = Vec3(0,0,theta*pitch);
    }


    // The generalized speed is the angular velocity of M in the F frame,
    // about F's z axis, expressed in F. (This axis is also constant in M.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0,0,1), Row3(0,0,pitch) );
    }

    // Since the Jacobian above is constant in F, its time derivative in F is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // CYLINDER //

// This is a "cylinder" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along the same axis. For molecules you can think of this as a combination
// of torsion and bond stretch. The axis used is the z axis of the parent's
// F frame, which is aligned forever with the z axis of the body's M frame.
// In addition, the origin points of M and F are separated only along the
// z axis; i.e., they have the same x & y coords in the F frame. The two
// generalized coordinates are the rotation and the translation, in that order.
class RBNodeCylinder : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "cylinder"; }

    RBNodeCylinder(const MassProperties& mProps_B,
                   const Transform&      X_PF,
                   const Transform&      X_BM,
                   UIndex&               nextUSlot,
                   USquaredIndex&        nextUSqSlot,
                   QIndex&               nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }


    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our cylinder joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        toQ(q)[0] = angles123[2];
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // Because the M and F origins must lie along their shared z axis, there is no way to
        // create a translation by rotating around z. So the only translation we can represent
        // is that component which is along z.
        toQ(q)[1] = T_FM[2];
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_FM[2];
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // Because the M and F origins must lie along their shared z axis, there is no way to
        // create a linear velocity by rotating around z. So the only linear velocity we can represent
        // is that component which is along z.
        toU(u)[1] = v_FM[2];
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Cylinder joint has one angular coordinate, which comes first.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Vec2& coords  = fromQ(q);

        X_FM.updR().setRotationFromAngleAboutZ(coords[0]);
        X_FM.updT() = Vec3(0,0,coords[1]);
    }


    // The generalized speeds are (1) the angular velocity of M in the F frame,
    // about F's z axis, expressed in F, and (2) the velocity of M's origin
    // in F, along F's z axis. (The z axis is also constant in M for this joint.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0,0,1), Row3(0)     );
        H_FM[1] = SpatialRow( Row3(0),     Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in F, its time derivative in F is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
    }

};


    // BEND-STRETCH //

// This is a "bend-stretch" joint, meaning one degree of rotational freedom
// about a particular axis, and one degree of translational freedom
// along a perpendicular axis. The z axis of the parent's F frame is 
// used for rotation (and that is always aligned with the M frame z axis).
// The x axis of the *M* frame is used for translation; that is, first
// we rotate around z, which moves M's x with respect to F's x. Then
// we slide along the rotated x axis. The two
// generalized coordinates are the rotation and the translation, in that order.
class RBNodeBendStretch : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "bendstretch"; }

    RBNodeBendStretch(const MassProperties& mProps_B,
                      const Transform&      X_PF,
                      const Transform&      X_BM,
                      UIndex&               nextUSlot,
                      USquaredIndex&        nextUSqSlot,
                      QIndex&               nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }


    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our bend-stretch joint can handle is about z.
        // TODO: this code is bad -- see comments for Torsion joint above.
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        toQ(q)[0] = angles123[2];
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // We can represent any translation that puts the M origin in the x-y plane of F,
        // by a suitable rotation around z followed by translation along x.
        const Vec2 r = T_FM.getSubVec<2>(0); // (rx, ry)

        // If we're not allowed to change rotation then we can only move along Mx.
//        if (only) {
//            const Real angle = fromQ(q)[0];
//            const Vec2 Mx(std::cos(angle), std::sin(angle)); // a unit vector
//            toQ(q)[1] = dot(r,Mx);
//            return;
//        }

        const Real d = r.norm();

        // If there is no translation worth mentioning, we'll leave the rotational
        // coordinate alone, otherwise rotate so M's x axis is aligned with r.
        if (d >= 4*Eps) {
            const Real angle = std::atan2(r[1],r[0]);
            toQ(q)[0] = angle;
            toQ(q)[1] = d;
        } else
            toQ(q)[1] = 0;
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {
        // We can only represent an angular velocity along z with this joint.
        toU(u)[0] = w_FM[2];
    }

    // If the translational coordinate is zero, we can only represent a linear velocity 
    // of OM in F which is along M's current x axis direction. Otherwise, we can 
    // represent any velocity in the x-y plane by introducing angular velocity about z.
    // We can never represent a linear velocity along z.
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
    {
        // Decompose the requested v into "along Mx" and "along My" components.
        const Rotation R_FM = Rotation( fromQ(q)[0], ZAxis ); // =[ Mx My Mz ] in F
        const Vec3 v_FM_M = ~R_FM*v_FM; // re-express in M frame

        toU(u)[1] = v_FM_M[0]; // velocity along Mx we can represent directly

//        if (only) {
//            // We can't do anything about My velocity if we're not allowed to change
//            // angular velocity, so we're done.
//            return;
//        }

        const Real x = fromQ(q)[1]; // translation along Mx (signed)
        if (std::abs(x) < SignificantReal) {
            // No translation worth mentioning; we can only do x velocity, which we just set above.
            return;
        }

        // significant translation
        toU(u)[0] = v_FM_M[1] / x; // set angular velocity about z to produce vy
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Bend-stretch joint has one angular coordinate, which comes first.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const Vec2& coords  = fromQ(q);    // angular coordinate

        X_FM.updR().setRotationFromAngleAboutZ(coords[0]);
        X_FM.updT() = X_FM.R()*Vec3(coords[1],0,0); // because translation is in M frame
    }

    // The generalized speeds for this bend-stretch joint are (1) the angular
    // velocity of M in the F frame, about F's z axis, expressed in F, and
    // (2) the (linear) velocity of M's origin in F, along *M's* current x axis
    // (that is, after rotation about z). (The z axis is also constant in M for this joint.)
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F

        const Vec3&     T_FM = getX_FM(pc).T();
        H_FM[0] = SpatialRow( Row3(0,0,1), ~(Vec3(0,0,1) % T_FM)   );
        H_FM[1] = SpatialRow( Row3(0),              ~Mx_F          );
    }

    // Since the the Jacobian above is not constant in F,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F

        const Vec3&     w_FM = getV_FM(vc)[0]; // angular velocity of M in F
        const Vec3&     v_FM = getV_FM(vc)[1]; // linear velocity of OM in F

        H_FM_Dot[0] = SpatialRow( Row3(0), ~(Vec3(0,0,1) % v_FM) );
        H_FM_Dot[1] = SpatialRow( Row3(0), ~(w_FM % Mx_F) );
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
// Note that the U-Joint degrees of freedom relating the parent's F frame
// to the child's M frame are about x and y, with the "long" axis of the
// driveshaft along z.
class RBNodeUJoint : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "ujoint"; }

    RBNodeUJoint(const MassProperties& mProps_B,
                 const Transform&      X_PF,
                 const Transform&      X_BM,
                 UIndex&               nextUSlot,
                 USquaredIndex&        nextUSqSlot,
                 QIndex&               nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotations this joint can handle are about Mx and My.
        // TODO: isn't there a better way to come up with "the rotation around x&y that
        // best approximates a rotation R"? Here we're just hoping that the supplied
        // rotation matrix can be decomposed into (x,y) rotations.
        const Vec2 angles12 = R_FM.convertRotationToBodyFixedXY();
        toQ(q) = angles12;
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    // We can only express angular velocity that can be produced with our generalized
    // speeds which are Fx and My rotations rates. So we'll take the supplied angular velocity
    // expressed in F, project it on Fx and use that as the first generalized speed. Then
    // take whatever angular velocity is unaccounted for, express it in M, and project onto
    // My and use that as the second generalized speed.
    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {
        const Rotation R_FM = Rotation( BodyRotationSequence, fromQ(q)[0], XAxis, fromQ(q)[1], YAxis );  // body fixed 1-2 sequence
        const Vec3     wyz_FM_M = ~R_FM*Vec3(0,w_FM[1],w_FM[2]);
        toU(u) = Vec2(w_FM[0], wyz_FM_M[1]);
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // U-joint has two angular coordinates.
        startOfAngles = MobilizerQIndex(0); nAngles=2; 
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        // We're only updating the orientation here because a U-joint can't translate.
        X_FM.updR() = Rotation( BodyRotationSequence, fromQ(q)[0], XAxis, fromQ(q)[1], YAxis );  // body fixed 1-2 sequence
        X_FM.updT() = 0.;
    }

    // The generalized speeds for this 2-dof rotational joint are the time derivatlves of
    // the body-fixed 1-2 rotation sequence defining the orientation. That is, the first speed
    // is just a rotation rate about Fx. The second is a rotation rate about the current My, so
    // we have to transform it into F to make H_FM uniformly expressed in F.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache();
        const Rotation& R_FM = getX_FM(pc).R();

        H_FM[0] = SpatialRow(  Row3(1,0,0) , Row3(0) );
        H_FM[1] = SpatialRow( ~(R_FM.y()) , Row3(0) );
    }

    // Since the second row of the Jacobian H_FM above is not constant in F,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     w_FM = getV_FM(vc)[0]; // angular velocity of M in F

        H_FM_Dot[0] = SpatialRow(        Row3(0)      , Row3(0) );
        H_FM_Dot[1] = SpatialRow( ~(w_FM % R_FM.y()), Row3(0) );
    }

};



    // PLANAR //

// This provides free motion (translation and rotation) in a plane. We use
// the 2d coordinate system formed by the x,y axes of F as the translations,
// and the common z axis of F and M as the rotational axis. The generalized
// coordinates are theta,x,y interpreted as rotation around z and translation
// along the (space fixed) Fx and Fy axes.
class RBNodePlanar : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "planar"; }

    RBNodePlanar(const MassProperties& mProps_B,
                 const Transform&      X_PF,
                 const Transform&      X_BM,
                 UIndex&               nextUSlot,
                 USquaredIndex&        nextUSqSlot,
                 QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

        // Implementations of virtual methods.

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        // The only rotation our planar joint can handle is about z.
        // TODO: should use 321 to deal with singular configuration (angle2==pi/2) better;
        // in that case 1 and 3 are aligned and the conversion routine allocates all the
        // rotation to whichever comes first.
        // TODO: isn't there a better way to come up with "the rotation around z that
        // best approximates a rotation R"?
        const Vec3 angles123 = R_FM.convertRotationToBodyFixedXYZ();
        toQ(q)[0] = angles123[2];
    }
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3&  T_FM, Vector& q) const {
        // Ignore translation in the z direction.
        toQ(q)[1] = T_FM[0]; // x
        toQ(q)[2] = T_FM[1]; // y
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM, Vector& u) const {
        // We can represent the z angular velocity exactly, but nothing else.
        toU(u)[0] = w_FM[2];
    }
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // Ignore translational velocity in the z direction.
        toU(u)[1] = v_FM[0]; // x
        toU(u)[2] = v_FM[1]; // y
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Planar joint has one angular coordinate, which comes first.
        startOfAngles = MobilizerQIndex(0); nAngles=1; 
        return true;
    }

    // This is required but does nothing here since there are no rotations for this joint.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
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

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& mv,
        const Vector&        q,
        Transform&           X_FM) const
    {
        // Rotational q is about common z axis, translational q's along Fx and Fy.
        X_FM = Transform(Rotation( fromQ(q)[0], ZAxis ), 
                          Vec3(fromQ(q)[1], fromQ(q)[2], 0));
    }

    // The rotational generalized speed is about the common z axis; translations
    // are along Fx and Fy so all axes are constant in F.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(0,0,1),   Row3(0) );
        H_FM[1] = SpatialRow(   Row3(0),   Row3(1,0,0) );
        H_FM[2] = SpatialRow(   Row3(0),   Row3(0,1,0) );
    }

    // Since the Jacobian above is constant in F, its time derivative is zero.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

};

    // GIMBAL //

// Gimbal joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's F frame.
// The generalized coordinates are:
//   * 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into 3 Euler angle derivatives.
//
// NOTE: This joint has a singularity when the middle angle is near 90 degrees.
// In most cases you should use a Ball joint instead, which by default uses
// a quaternion as its generalized coordinates to avoid this singularity.

class RBNodeGimbal : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "gimbal"; }

    RBNodeGimbal( const MassProperties& mProps_B,
                  const Transform&      X_PF,
                  const Transform&      X_BM,
                  UIndex&               nextUSlot,
                  USquaredIndex&        nextUSqSlot,
                  QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionIsNeverUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const 
    {
        toQ(q) = R_FM.convertRotationToBodyFixedXYZ();
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM,
                                  Vector& u) const
    {
        toU(u) = w_FM; // relative angular velocity always used as generalized speeds
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Gimbal joint has three angular coordinates.
        startOfAngles = MobilizerQIndex(0); nAngles=3; 
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const Vec3& a = fromQ(q); // angular coordinates
        toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
        toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
        // no quaternions
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        X_FM.updT() = 0.; // This joint can't translate.
        X_FM.updR().setRotationToBodyFixedXYZ( fromQ(q) );
    }

    // Generalized speeds are the angular velocity expressed in F, so they
    // cause rotations around F x,y,z axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(1,0,0), Row3(0) );
        H_FM[1] = SpatialRow( Row3(0,1,0), Row3(0) );
        H_FM[2] = SpatialRow( Row3(0,0,1), Row3(0) );
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM = fromU(u); // angular velocity of M in F 
        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(sbs.getQ()),
                                    ~R_FM*w_FM); // need w in *body*, not parent
    }

    void calcLocalQDotFromLocalU(const SBStateDigest& sbs, const Real* u, Real* qdot) const {
        assert(sbs.getStage() >= Stage::Position);
        assert(u && qdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();

        const Vec3&            w_FM = Vec3::getAs(u);

        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(allQ),
                                         ~R_FM*w_FM); // need w in *body*, not parent
    }

    // Compute out_q = Q * in_u
    //   or    out_u = in_q * Q
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        // TODO: it's annoying that this Q block is only available in the Body (M) frame,
        // because this mobilizer uses angular velocity in the Parent (F) frame
        // as generalized speeds. So we have to do an expensive conversion here. It
        // would be just as easy to compute this matrix in the Parent frame in 
        // the first place.
        const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
        const Mat33    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
        if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * Q;
        else               Vec3::updAs(out) = Q * Vec3::getAs(in);
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        // TODO: see above regarding the need for this R_FM kludge
        const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
        const Mat33    QInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
        if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * QInv;
        else               Vec3::updAs(out) = QInv * Vec3::getAs(in);
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM     = fromU(sbs.getU()); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_FM_dot = fromU(udot);

        const Rotation& R_FM = getX_FM(pc).R();
        toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                              (fromQ(sbs.getQ()), ~R_FM*w_FM, ~R_FM*w_FM_dot);
    }

    void calcLocalQDotDotFromLocalUDot(const SBStateDigest& sbs, const Real* udot, Real* qdotdot) const {
        assert(sbs.getStage() >= Stage::Velocity);
        assert(udot && qdotdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();
        const Vector&          allU = sbs.getU();

        const Vec3& w_FM     = fromU(allU);
        const Vec3& w_FM_dot = Vec3::getAs(udot);

        const Rotation& R_FM = getX_FM(pc).R();
        Vec3::updAs(qdotdot) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                    (fromQ(allQ), ~R_FM*w_FM, ~R_FM*w_FM_dot);
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

    // ORIENTATION (BALL) //

// Ball joint. This provides three degrees of rotational freedom,  i.e.,
// unrestricted orientation of the body's M frame in the parent's F frame.
// The generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame.
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeBall : public RigidBodyNodeSpec<3> {
public:
    virtual const char* type() { return "ball"; }

    RBNodeBall(const MassProperties& mProps_B,
               const Transform&      X_PF,
               const Transform&      X_BM,
               UIndex&               nextUSlot,
               USquaredIndex&        nextUSqSlot,
               QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionMayBeUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQ(q)    = R_FM.convertRotationToBodyFixedXYZ();
        else
            toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector&, const Vec3& w_FM,
                                     Vector& u) const
    {
            toU(u) = w_FM; // relative angular velocity always used as generalized speeds
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Ball joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q); // angular coordinates
            toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const SBModelVars& mv = sbs.getModelVars();
        X_FM.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(mv))
            X_FM.updR().setRotationToBodyFixedXYZ( fromQ(q) );
        else {
            // TODO: should use qnorm pool
            X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
        }
    }

    // Generalized speeds are the angular velocity expressed in F, so they
    // cause rotations around F x,y,z axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(1,0,0), Row3(0) );
        H_FM[1] = SpatialRow( Row3(0,1,0), Row3(0) );
        H_FM[2] = SpatialRow( Row3(0,0,1), Row3(0) );
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM = fromU(u); // angular velocity of M in F 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(sbs.getQ()),
                                        ~R_FM*w_FM); // need w in *body*, not parent
        } else
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),w_FM);
    }

    void calcLocalQDotFromLocalU(const SBStateDigest& sbs, const Real* u, Real* qdot) const {
        assert(sbs.getStage() >= Stage::Position);
        assert(u && qdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();

        const Vec3&            w_FM = Vec3::getAs(u);

        if (getUseEulerAngles(mv)) {
            Vec4::updAs(qdot) = 0; // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            Vec3::updAs(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(allQ),
                                             ~R_FM*w_FM); // need w in *body*, not parent
        } else
            Vec4::updAs(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(allQ),w_FM);
    }

    // CAUTION: we do not zero the unused 4th element of q for Euler angles; it
    // is up to the caller to do that if it is necessary.
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: it's annoying that this Q block is only available in the Body (M) frame,
            // because this mobilizer uses angular velocity in the Parent (F) frame
            // as generalized speeds. So we have to do an expensive conversion here. It
            // would be just as easy to compute this matrix in the Parent frame in 
            // the first place.
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * Q;
            else               Vec3::updAs(out) = Q * Vec3::getAs(in);
        } else {
            // Quaternion
            const Mat43 Q = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) Row3::updAs(out) = Row4::getAs(in) * Q;
            else               Vec4::updAs(out) = Q * Vec3::getAs(in);
        }
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: see above regarding the need for this R_FM kludge
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    QInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * QInv;
            else               Vec3::updAs(out) = QInv * Vec3::getAs(in);
        } else {
            
            // Quaternion
            const Mat34 QInv = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) Row4::updAs(out) = Row3::getAs(in) * QInv;
            else               Vec3::updAs(out) = QInv * Vec4::getAs(in);
        }
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM     = fromU(sbs.getU()); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_FM_dot = fromU(udot);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                                  (fromQ(sbs.getQ()), ~R_FM*w_FM, ~R_FM*w_FM_dot);
        } else
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(sbs.getQ()),w_FM,w_FM_dot);
    }

    void calcLocalQDotDotFromLocalUDot(const SBStateDigest& sbs, const Real* udot, Real* qdotdot) const {
        assert(sbs.getStage() >= Stage::Velocity);
        assert(udot && qdotdot);

        const SBModelVars&     mv   = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const Vector&          allQ = sbs.getQ();
        const Vector&          allU = sbs.getU();

        const Vec3& w_FM     = fromU(allU);
        const Vec3& w_FM_dot = Vec3::getAs(udot);

        if (getUseEulerAngles(mv)) {
            Vec4::updAs(qdotdot) = 0; // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            Vec3::updAs(qdotdot) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                        (fromQ(allQ), ~R_FM*w_FM, ~R_FM*w_FM_dot);
        } else
            Vec4::updAs(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                        (fromQuat(allQ),w_FM,w_FM_dot);
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
    int getNQInUse(const SBModelVars& mv) const {
        return getUseEulerAngles(mv) ? 3 : 4;
    } 
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
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
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const 
    {
        if (getUseEulerAngles(sbs.getModelVars())) 
            return false;   // no change

        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();

        if (qErrest.size()) {
            Vec4& qerr = toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }

        return true;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        toQuat(outputQ) = Vec4(0); // clear unused element
        toQ(outputQ) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(fromQ(inputQ));
        toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
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



    // ELLIPSOID //

// ELLIPSOID mobilizer. This provides three degrees of rotational freedom, i.e.,
// unrestricted orientation, of the body's M frame in the parent's F frame, 
// along with coordinated translation that keeps the M frame origin on
// the surface of an ellipsoid fixed in F and centered on the F origin.
// The surface point is chosen for
// a given orientation of M in F as the unique point on the ellipsoid surface
// where the surface normal is aligned with Mz. That is, the z axis of M is 
// assumed to be normal to the ellipsoid at all times, and the translation is
// chosen to make that true.
//
// Unlike most joints, the reference configuration (i.e., X_FM when q=0 is
// not the identity transform. Instead, although the frames are aligned the
// M frame origin is offset from F along their shared +z axis, so that it lies
// on the ellipsoid surface at the point (0,0,rz) where rz is the z-radius
// (semiaxis) of the ellipsoid.
//
// The generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//     In Euler angles, axis 3 is just the spin of the outboard body about
//     its Mz axis, which is always normal to the ellipse. For small 1,2 angles
//     you can think of angle 1 (about x) as latitude, and angle 2 (about y)
//     as longitude when viewing the ellipsoid by looking down the z axis.
//     (That would be true for large angles too if we were using space fixed
//     angles, that is, angles defined about F rather than M axes.)
//
// Generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame.
//
// Thus rotational qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
//
// This mobilizer was written by Ajay Seth and hacked somewhat by Sherm.

//
class RBNodeEllipsoid : public RigidBodyNodeSpec<3> {
	Vec3 semi; // semi axis dimensions in x,y,z resp.
public:
    virtual const char* type() { return "ellipsoid"; }

    RBNodeEllipsoid(const MassProperties& mProps_B,
                  const Transform&      X_PF,
                  const Transform&      X_BM,
                  const Vec3&           radii, // x,y,z
                  UIndex&               nextUSlot,
                  USquaredIndex&        nextUSqSlot,
                  QIndex&               nextQSlot)
      : RigidBodyNodeSpec<3>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionMayBeUsed),
        semi(radii)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                           Vector& q) const 
    {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQ(q)    = R_FM.convertRotationToBodyFixedXYZ();
        else
            toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
    }

    // We can't hope to represent arbitrary translations with a joint that has only
    // rotational coordinates! However, since F is at the center of the ellipsoid and
    // M on its surface, we can at least obtain a translation in the *direction* of
    // the requested translation. The magnitude must of course be set to end up with
    // the M origin right on the surface of the ellipsoid, and Mz will be the normal
    // at that point.
    //
    // Expressed as an x-y-z body fixed Euler sequence, the z rotation is just the
    // spin around the Mz (surface normal) and could be anything, so we'll just leave
    // it at its current value. The x and y rotations act like polar coordinates
    // to get the M origin point on the direction indicated by the requested translation.
    //
    // If the requested translation is near zero we can't do anything since we can't find
    // a direction to align with. And of course we can't do anything if "only" is true
    // here -- that means we aren't allowed to touch the rotations, and for this
    // joint that's all there is.
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        if (T_FM.norm() < Eps) return;

        const UnitVec3 e(T_FM); // direction from F origin towards desired M origin
        const Real latitude  = std::atan2(-e[1],e[2]); // project onto F's yz plane
        const Real longitude = std::atan2( e[0],e[2]); // project onto F's xz plane

		// Calculate the current value of the spin coordinate (3rd Euler angle).
        Real spin;
		if (getUseEulerAngles(sbs.getModelVars())){
            spin = fromQ(q)[2];
		} else {
			const Rotation R_FM_now(Quaternion(fromQuat(q)));
            spin = R_FM_now.convertRotationToBodyFixedXYZ()[2];
        }

        // Calculate the desired rotation, which is a space-fixed 1-2 sequence for
        // latitude and longitude, followed by a body fixed rotation for spin.
        const Rotation R_FM = Rotation( SpaceRotationSequence, latitude, XAxis, longitude, YAxis ) * Rotation(spin,ZAxis);

        if (getUseEulerAngles(sbs.getModelVars())) {
            const Vec3 q123 = R_FM.convertRotationToBodyFixedXYZ();
            toQ(q) = q123;
        } else {
            const Quaternion quat = R_FM.convertRotationToQuaternion();
			toQuat(q) = quat.asVec4();
        }
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
                                  Vector& u) const
    {
            toU(u) = w_FM; // relative angular velocity always used as generalized speeds
    }

    // We can't do general linear velocity with this rotation-only mobilizer, but we can
    // express any velocity which is tangent to the ellipsoid surface. So we'll find the
    // current surface normal (Mz) and ignore any component of the requested velocity
    // which is along that direction. (The resulting vz won't be zero, though, but it
    // is completely determined by vx,vy.)
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
    {
        Transform X_FM;
        calcAcrossJointTransform(sbs,q,X_FM);

        const Vec3 v_FM_M    = ~X_FM.R()*v_FM; // we can only do vx and vy in this frame
        const Vec3 r_FM_M    = ~X_FM.R()*X_FM.T(); 
        const Vec3 wnow_FM_M = ~X_FM.R()*fromU(u); // preserve z component

        // Now vx can only result from angular velocity about y, vy from x.
        // TODO: THIS IS ONLY RIGHT FOR A SPHERE!!!
        const Real wx = -v_FM_M[1]/r_FM_M[2];
        const Real wy = v_FM_M[0]/r_FM_M[2];
        const Vec3 w_FM_M(wx, wy, wnow_FM_M[2]);
        const Vec3 w_FM = X_FM.R()*w_FM_M;

        toU(u) = w_FM;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Ellipsoid joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q); // angular coordinates
            toQ(sine)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();

            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());

            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
		// Calcuate the rotation R_FM first.
        const SBModelVars& mv = sbs.getModelVars();
		if (getUseEulerAngles(mv)){
			const Vec3& a = fromQ(q);		
			X_FM.updR().setRotationToBodyFixedXYZ(a);
		} else {
            // TODO: should use qnorm pool
            //       Conversion to Quaternion here involves expensive normalization
            //       because state variables q can never be assumed normalized.
			Quaternion quat = Quaternion(fromQuat(q));
            X_FM.updR().setRotationFromQuaternion(quat);
        }

		const Vec3& n = X_FM.R().z();
		X_FM.updT() = Vec3(semi[0]*n[0], semi[1]*n[1], semi[2]*n[2]);
	}

    // Generalized speeds are the angular velocity expressed in F, so they
    // cause rotations around F x,y,z axes respectively.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache();
		const Vec3& n = getX_FM(pc).z();

		H_FM[0] = SpatialRow( Row3(1,0,0), Row3(      0,      -n[2]*semi[1], n[1]*semi[2]) );
		H_FM[1] = SpatialRow( Row3(0,1,0), Row3( n[2]*semi[0],       0,     -n[0]*semi[2]) );
		H_FM[2] = SpatialRow( Row3(0,0,1), Row3(-n[1]*semi[0], n[0]*semi[1],       0     ) );
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
		const Vec3& n     = getX_FM(pc).z(); // ellipsoid normal
		const Vec3& w_FM = getV_FM(vc)[0];  // angular velocity of M in F

		const Vec3 ndot = w_FM % n; // time derivative of n in F

		H_FM_Dot[0] = SpatialRow( Row3(0), Row3(      0,         -ndot[2]*semi[1], ndot[1]*semi[2]) );
		H_FM_Dot[1] = SpatialRow( Row3(0), Row3( ndot[2]*semi[0],       0,        -ndot[0]*semi[2]) );
		H_FM_Dot[2] = SpatialRow( Row3(0), Row3(-ndot[1]*semi[0], ndot[0]*semi[1],       0        ) );
    }

    // CAUTION: we do not zero the unused 4th element of q for Euler angles; it
    // is up to the caller to do that if it is necessary.
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: it's annoying that this Q block is only available in the Body (M) frame,
            // because this mobilizer uses angular velocity in the Parent (F) frame
            // as generalized speeds. So we have to do an expensive conversion here. It
            // would be just as easy to compute this matrix in the Parent frame in 
            // the first place.
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * Q;
            else               Vec3::updAs(out) = Q * Vec3::getAs(in);
        } else {
            // Quaternion
            const Mat43 Q = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) Row3::updAs(out) = Row4::getAs(in) * Q;
            else               Vec4::updAs(out) = Q * Vec3::getAs(in);
        }
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: see above regarding the need for this R_FM kludge
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    QInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * QInv;
            else               Vec3::updAs(out) = QInv * Vec3::getAs(in);
        } else {
            
            // Quaternion
            const Mat34 QInv = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) Row4::updAs(out) = Row3::getAs(in) * QInv;
            else               Vec3::updAs(out) = QInv * Vec4::getAs(in);
        }
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM = fromU(u); // angular velocity of M in F 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            toQ(qdot) = Rotation::convertAngVelToBodyFixed123Dot(fromQ(sbs.getQ()),
                                        ~R_FM*w_FM); // need w in *body*, not parent
        } else
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),w_FM);
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM     = fromU(sbs.getU()); // angular velocity of J in Jb, expr in Jb
        const Vec3& w_FM_dot = fromU(udot);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot) = Vec4(0); // TODO: kludge, clear unused element
            const Rotation& R_FM = getX_FM(pc).R();
            toQ(qdotdot)    = Rotation::convertAngVelDotToBodyFixed123DotDot
                                  (fromQ(sbs.getQ()), ~R_FM*w_FM, ~R_FM*w_FM_dot);
        } else
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(sbs.getQ()),w_FM,w_FM_dot);
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

    int getMaxNQ() const {return 4;}
    int getNQInUse(const SBModelVars& mv) const {
        return getUseEulerAngles(mv) ? 3 : 4;
    } 
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
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
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const 
    {
        if (getUseEulerAngles(sbs.getModelVars())) 
            return false;   // no change

        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();

        if (qErrest.size()) {
            Vec4& qerr = toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }

        return true;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        toQuat(outputQ) = Vec4(0); // clear unused element
        toQ(outputQ) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(fromQ(inputQ));
        toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
    }

    void getInternalForce(const SBAccelerationCache&, Vector&) const {
        assert(false); // TODO: decompose cross-joint torque into 123 gimbal torques
    }
};

    // FREE //

// Free joint. This provides six degrees of freedom, three rotational and
// three translational. The rotation is like the ball joint above; the
// translation is like the Cartesian joint above.
// TODO: to get this to work I had to make the translations be in the outboard
// frame (M, not F). So currently the generalized coordinates are:
//   * 4 quaternions or 3 1-2-3 body fixed Euler angles (that is, fixed in M)
//   * translation from OF to OM as a 3-vector in the outboard body mobilizer (M) frame
// and generalized speeds are:
//   * angular velocity w_FM as a vector expressed in the F frame
//   * linear velocity of the M origin in F (v_FM), expressed in M
// Thus translational qdots are just generalized speeds, but rotational
// qdots have to be derived from the generalized speeds to be turned into
// either 4 quaternion derivatives or 3 Euler angle derivatives.
//   
class RBNodeFree : public RigidBodyNodeSpec<6> {
public:
    virtual const char* type() { return "free"; }

    RBNodeFree(const MassProperties& mProps_B,
               const Transform&      X_PF,
               const Transform&      X_BM,
               UIndex&               nextUSlot,
               USquaredIndex&        nextUSqSlot,
               QIndex&               nextQSlot)
      : RigidBodyNodeSpec<6>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionMayBeUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,0) = R_FM.convertRotationToBodyFixedXYZ();
        else
            toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
    }

    // The user gives us the translation vector from OF to OM as a vector expressed in F, which
    // is what we use as translational generalized coordinates. Also, with a free joint 
    // we never have to change orientation coordinates in order to achieve a translation.
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,3) = T_FM; // skip the 3 Euler angles
        else
            toQVec3(q,4) = T_FM; // skip the 4 quaternions
    }

    // Our 3 rotational generalized speeds are just the angular velocity vector of M in F,
    // expressed in F, which is exactly what the user provides here.
    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
                                     Vector& u) const
    {
        toUVec3(u,0) = w_FM; // relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in F,
    // expressed in F, which is just what the user gives us.
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
    {
        toUVec3(u,3) = v_FM;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // Free joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        if (getUseEulerAngles(mv)) {
            X_FM.updR().setRotationToBodyFixedXYZ( fromQVec3(q,0) );
            X_FM.updT() = fromQVec3(q,3); // translation is in F already
        } else {
            X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
            X_FM.updT() = fromQVec3(q,4); // translation is in F already
        }
    }


    // The generalized speeds for this 6-dof ("free") joint are 
    //   (1) the angular velocity of M in the F frame, expressed in F, and
    //   (2) the (linear) velocity of M's origin in F, expressed in F.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        H_FM[0] = SpatialRow( Row3(1,0,0),   Row3(0)   );  // rotations
        H_FM[1] = SpatialRow( Row3(0,1,0),   Row3(0)   );
        H_FM[2] = SpatialRow( Row3(0,0,1),   Row3(0)   );

        H_FM[3] = SpatialRow(   Row3(0),   Row3(1,0,0) );  // translations
        H_FM[4] = SpatialRow(   Row3(0),   Row3(0,1,0) );
        H_FM[5] = SpatialRow(   Row3(0),   Row3(0,0,1) );
    }

    // Since the Jacobian above is constant in F, its derivative in F is 0.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        H_FM_Dot[0] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[1] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );

        H_FM_Dot[3] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[4] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[5] = SpatialRow( Row3(0), Row3(0) );
    }

    // CAUTION: we do not zero the unused 4th element of q for Euler angles; it
    // is up to the caller to do that if it is necessary.
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: it's annoying that this Q block is only available in the Body (M) frame,
            // because this mobilizer uses angular velocity in the Parent (F) frame
            // as generalized speeds. So we have to do an expensive conversion here. It
            // would be just as easy to compute this matrix in the Parent frame in 
            // the first place.
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q)) * ~R_FM;
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * Q;
            else               Vec3::updAs(out) = Q * Vec3::getAs(in);
            // translational part of Q block is identity
            Vec3::updAs(out+3) = Vec3::getAs(in+3);
        } else {
            // Quaternion
            const Mat43 Q = Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) {
                Row3::updAs(out)   = Row4::getAs(in) * Q;
                Row3::updAs(out+3) = Row3::getAs(in+4); // translational part of Q block is identity
            } else { // matrix on left
                Vec4::updAs(out) = Q * Vec3::getAs(in);
                Vec3::updAs(out+4) = Vec3::getAs(in+3); // translational part of Q block is identity
            }
        }
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            // TODO: see above regarding the need for this R_FM kludge
            const Rotation R_FM(BodyRotationSequence, q[0], XAxis, q[1], YAxis, q[2], ZAxis);
            const Mat33    QInv(R_FM*Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q)));
            if (matrixOnRight) Row3::updAs(out) = Row3::getAs(in) * QInv;
            else               Vec3::updAs(out) = QInv * Vec3::getAs(in);
            // translational part of QInv block is identity
            Vec3::updAs(out+3) = Vec3::getAs(in+3);
        } else {           
            // Quaternion
            const Mat34 QInv = Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q));
            if (matrixOnRight) {
                Row4::updAs(out) = Row3::getAs(in) * QInv;
                Row3::updAs(out+4) = Row3::getAs(in+3); // translational part of QInv block is identity
            } else { // matrix on left
                Vec3::updAs(out) = QInv * Vec4::getAs(in);
                Vec3::updAs(out+3) = Vec3::getAs(in+4); // translational part of QInv block is identity
            }
        }
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u,
        Vector&                qdot) const
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM = fromUVec3(u,0); // Angular velocity in F
        const Vec3& v_FM = fromUVec3(u,3); // Linear velocity in F
        if (getUseEulerAngles(mv)) {
            const Rotation& R_FM = getX_FM(pc).R();
            const Vec3& theta = fromQVec3(sbs.getQ(),0); // Euler angles
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
                                            ~R_FM*w_FM); // need w in *body*, not parent
            toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdot,3) = v_FM;
        } else {
            const Vec4& quat = fromQuat(sbs.getQ());
            toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,w_FM);
            toQVec3(qdot,4) = v_FM;
        }
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3& w_FM     = fromUVec3(sbs.getU(),0); // angular velocity of M in F
        const Vec3& v_FM     = fromUVec3(sbs.getU(),3); // linear velocity of M in F, expressed in F
        const Vec3& w_FM_dot = fromUVec3(udot,0);
        const Vec3& v_FM_dot = fromUVec3(udot,3);
        if (getUseEulerAngles(mv)) {
            const Rotation& R_FM = getX_FM(pc).R();
            const Vec3& theta  = fromQVec3(sbs.getQ(),0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, ~R_FM*w_FM, ~R_FM*w_FM_dot);
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_FM_dot;
        } else {
            const Vec4& quat  = fromQuat(sbs.getQ());
            toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                             (quat,w_FM,w_FM_dot);
            toQVec3(qdotdot,4) = v_FM_dot;
        }
    }

    void copyQ(const SBModelVars& mv, const Vector& qIn, Vector& q) const {
        if (getUseEulerAngles(mv))
            toQ(q) = fromQ(qIn);
        else {
            toQuat(q)    = fromQuat(qIn);
            toQVec3(q,4) = fromQVec3(qIn,4);
        }
    }

    int  getMaxNQ()                   const {return 7;}
    int  getNQInUse(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
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

    bool enforceQuaternionConstraints(
        const SBStateDigest& sbs, 
        Vector&             q,
        Vector&             qErrest) const
    {
        if (getUseEulerAngles(sbs.getModelVars())) 
            return false; // no change

        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();

        if (qErrest.size()) {
            Vec4& qerr = toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }

        return true;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = Vec3(0); // clear unused element
        toQVec3(outputQ, 3) = fromQVec3(inputQ, 4);
        toQVec3(outputQ, 0) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = fromQVec3(inputQ, 3);
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(fromQVec3(inputQ, 0));
        toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
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
//   * the x,y components of the angular velocity w_FM_M, that is, the angular
//     velocity of M in F expressed in M (where we want wz=0).
//     NOTE: THAT IS A DIFFERENT FRAME THAN IS USED FOR BALL AND GIMBAL
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeLineOrientation : public RigidBodyNodeSpec<2> {
public:
    virtual const char* type() { return "lineOrientation"; }

    RBNodeLineOrientation(const MassProperties& mProps_B,
                          const Transform&      X_PF,
                          const Transform&      X_BM,
                          UIndex&               nextUSlot,
                          USquaredIndex&        nextUSqSlot,
                          QIndex&               nextQSlot)
      : RigidBodyNodeSpec<2>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionMayBeUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,0)    = R_FM.convertRotationToBodyFixedXYZ();
        else
            toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
    }

    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a translation by rotating. So the only translation we can represent is 0.
    }

    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
                                     Vector& u) const
    {
        Rotation R_FM;
        if (getUseEulerAngles(sbs.getModelVars()))
            R_FM.setRotationToBodyFixedXYZ( fromQVec3(q,0) );
        else {
            // TODO: should use qnorm pool
            R_FM.setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
        }
        const Vec3 w_FM_M = ~R_FM*w_FM;
        toU(u) = Vec2(w_FM_M[0], w_FM_M[1]); // (x,y) of relative angular velocity always used as generalized speeds
    }

    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector&, const Vec3& v_FM, Vector& u) const
    {
        // M and F frame origins are always coincident for this mobilizer so there is no
        // way to create a linear velocity by rotating. So the only linear velocity
        // we can represent is 0.
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // LineOrientation joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQVec3(q,0); // angular coordinates
            toQVec3(sine,0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQVec3(cosine,0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const
    {
        const SBModelVars& mv = sbs.getModelVars();
        X_FM.updT() = 0.; // This joint can't translate.
        if (getUseEulerAngles(mv))
            X_FM.updR().setRotationToBodyFixedXYZ( fromQVec3(q,0) );
        else {
            // TODO: should use qnorm pool
            X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
        }
    }

    // The generalized speeds for this 2-dof rotational joint are the x and y
    // components of the angular velocity of M in the F frame, expressed in the *M*
    // frame.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
        const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

        H_FM[0] = SpatialRow( ~Mx_F, Row3(0) );
        H_FM[1] = SpatialRow( ~My_F, Row3(0) );
    }

    // Since the Jacobian above is not constant in F,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
        const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

        const Vec3&     w_FM = getV_FM(vc)[0]; // angular velocity of M in F

        H_FM_Dot[0] = SpatialRow( ~(w_FM % Mx_F), Row3(0) );
        H_FM_Dot[1] = SpatialRow( ~(w_FM % My_F), Row3(0) );
    }

    // CAUTION: we do not zero the unused 4th element of q for Euler angles; it
    // is up to the caller to do that if it is necessary.
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            const Mat32    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                    .getSubMat<3,2>(0,0); // drop 3rd column
            if (matrixOnRight) Row2::updAs(out) = Row3::getAs(in) * Q;
            else               Vec3::updAs(out) = Q * Vec2::getAs(in);
        } else {
            // Quaternion: Q block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat42 Q = (Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q))*R_FM)
                                .getSubMat<4,2>(0,0); // drop 3rd column
            if (matrixOnRight) Row2::updAs(out) = Row4::getAs(in) * Q;
            else               Vec4::updAs(out) = Q * Vec2::getAs(in);
        }
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            const Mat23    QInv = Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                        .getSubMat<2,3>(0,0); // drop 3rd row
            if (matrixOnRight) Row3::updAs(out) = Row2::getAs(in) * QInv;
            else               Vec2::updAs(out) = QInv * Vec3::getAs(in);
        } else {
            // Quaternion: QInv block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat24 QInv = (~R_FM*Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q)))
                                    .getSubMat<2,4>(0,0);   // drop 3rd row
            if (matrixOnRight) Row4::updAs(out) = Row2::getAs(in) * QInv;
            else               Vec2::updAs(out) = QInv * Vec4::getAs(in);
        }
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u, 
        Vector&                qdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3 w_FM_M = fromU(u).append1(0); // angular velocity of M in F, exp in M (with wz=0) 
        if (getUseEulerAngles(mv)) {
            toQuat(qdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(fromQVec3(sbs.getQ(),0),
                                        w_FM_M); // need w in *body*, not parent
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            toQuat(qdot) = Rotation::convertAngVelToQuaternionDot(fromQuat(sbs.getQ()),
                                        R_FM*w_FM_M); // need w in *parent* frame
        }
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3 w_FM_M     = fromU(sbs.getU()).append1(0); // angular velocity of M in F, exp in M (with wz=0)
        const Vec3 w_FM_M_dot = fromU(udot).append1(0);

        if (getUseEulerAngles(mv)) {
            toQuat(qdotdot)    = Vec4(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                       (fromQVec3(sbs.getQ(),0), w_FM_M, w_FM_M_dot); // body frame
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            toQuat(qdotdot) = Rotation::convertAngVelDotToQuaternionDotDot
                                  (fromQuat(sbs.getQ()),R_FM*w_FM_M,R_FM*w_FM_M_dot); // parent frame
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
    int getNQInUse(const SBModelVars& mv) const {
        return getUseEulerAngles(mv) ? 3 : 4;
    } 
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
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
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const 
    {
        if (getUseEulerAngles(sbs.getModelVars())) 
            return false;   // no change

        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();

        if (qErrest.size()) {
            Vec4& qerr = toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }

        return true;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = Vec3(0); // clear unused element
        toQVec3(outputQ, 2) = fromQVec3(inputQ, 3);
        toQVec3(outputQ, 0) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = fromQVec3(inputQ, 3);
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(fromQVec3(inputQ, 0));
        toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
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
//   * 3 components of the translation vector T_FM (that is, vector from origin
//     of F to origin of M, expressed in F)
// and generalized speeds are:
//   * the x,y components of the angular velocity w_FM_M, that is, the angular
//     velocity of M in F expressed in *M* (where we want wz=0).
//   * 3 components of the linear velocity of origin of M in F, expressed in F.
//     NOTE: THAT IS NOT THE SAME FRAME AS FOR A FREE JOINT
// Thus the qdots have to be derived from the generalized speeds to
// be turned into either 4 quaternion derivatives or 3 Euler angle derivatives.
class RBNodeFreeLine : public RigidBodyNodeSpec<5> {
public:
    virtual const char* type() { return "full"; }

    RBNodeFreeLine(const MassProperties& mProps_B,
                   const Transform&      X_PF,
                   const Transform&      X_BM,
                   UIndex&               nextUSlot,
                   USquaredIndex&        nextUSqSlot,
                   QIndex&               nextQSlot)
      : RigidBodyNodeSpec<5>(mProps_B,X_PF,X_BM,nextUSlot,nextUSqSlot,nextQSlot,
                             QDotMayDifferFromU, QuaternionMayBeUsed)
    {
        updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }

    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM,
                              Vector& q) const 
    {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,0) = R_FM.convertRotationToBodyFixedXYZ();
        else
            toQuat(q) = R_FM.convertRotationToQuaternion().asVec4();
    }

    // The user gives us the translation vector from OF to OM as a vector expressed in F.
    // With a free joint we never have to *change* orientation coordinates in order to achieve a translation.
    // Note: a quaternion from a state is not necessarily normalized so can't be used
    // direction as though it were a set of Euler parameters; it must be normalized first.
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        if (getUseEulerAngles(sbs.getModelVars()))
            toQVec3(q,3) = T_FM; // skip the 3 Euler angles
        else
            toQVec3(q,4) = T_FM; // skip the 4 quaternions
    }

    // Our 2 rotational generalized speeds are just the (x,y) components of the
    // angular velocity vector of M in F, expressed in M.
    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM,
                                     Vector& u) const
    {
        Rotation R_FM;
        if (getUseEulerAngles(sbs.getModelVars()))
            R_FM.setRotationToBodyFixedXYZ( fromQVec3(q,0) );
        else {
            // TODO: should use qnorm pool
            R_FM.setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
        }
        const Vec3 w_FM_M = ~R_FM*w_FM;
        toU(u).updSubVec<2>(0) = Vec2(w_FM_M[0], w_FM_M[1]); // (x,y) of relative angular velocity always used as generalized speeds
    }

    // Our 3 translational generalized speeds are the linear velocity of M's origin in F,
    // expressed in F. The user gives us that same vector.
    void setUToFitLinearVelocity
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const
    {
        toUVec3(u,2) = v_FM;
    }

    // This is required for all mobilizers.
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& nAngles) const {
        // FreeLine joint has three angular coordinates when Euler angles are being used, 
        // none when quaternions are being used.
        if (!getUseEulerAngles(sbs.getModelVars())) {startOfAngles.invalidate(); nAngles=0; return false;} 
        startOfAngles = MobilizerQIndex(0);
        nAngles = 3;
        return true;
    }

    // Precalculate sines and cosines.
    void calcJointSinCosQNorm(
        const SBModelVars&  mv,
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const
    {
        const SBModelCache::PerMobilizedBodyModelInfo& bInfo = mc.getMobilizedBodyModelInfo(nodeNum);

        if (getUseEulerAngles(mv)) {
            const Vec3& a = fromQ(q).getSubVec<3>(0); // angular coordinates
            toQ(sine).updSubVec<3>(0)   = Vec3(std::sin(a[0]), std::sin(a[1]), std::sin(a[2]));
            toQ(cosine).updSubVec<3>(0) = Vec3(std::cos(a[0]), std::cos(a[1]), std::cos(a[2]));
            // no quaternions
        } else {
            // no angles
            const Vec4& quat = fromQuat(q); // unnormalized quaternion from state
            const Real  quatLen = quat.norm();
            assert(bInfo.hasQuaternionInUse && bInfo.quaternionPoolIndex.isValid());
            qErr[ic.firstQuaternionQErrSlot+bInfo.quaternionPoolIndex] = quatLen - Real(1);
            toQuat(qnorm) = quat / quatLen;
        }
    }

    // Calculate X_FM.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        if (getUseEulerAngles(mv)) {
            X_FM.updR().setRotationToBodyFixedXYZ( fromQVec3(q,0) );
            X_FM.updT() = fromQVec3(q,3); // translation is in F
        } else {
            X_FM.updR().setRotationFromQuaternion( Quaternion(fromQuat(q)) ); // normalize
            X_FM.updT() = fromQVec3(q,4);  // translation is in F
        }
    }


    // The generalized speeds for this 5-dof ("free line") joint are 
    //   (1) the (x,y) components of angular velocity of M in the F frame, expressed in M, and
    //   (2) the (linear) velocity of M's origin in F, expressed in F.
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&               H_FM) const
    {
        const SBPositionCache& pc = sbs.updPositionCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
        const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

        H_FM[0] = SpatialRow( ~Mx_F, Row3(0) );        // x,y angular velocity in M, re-expressed im F
        H_FM[1] = SpatialRow( ~My_F, Row3(0) );

        H_FM[2] = SpatialRow( Row3(0), Row3(1,0,0) );   // translations in F
        H_FM[3] = SpatialRow( Row3(0), Row3(0,1,0) );
        H_FM[4] = SpatialRow( Row3(0), Row3(0,0,1) );
    }

    // Since the first two rows of the Jacobian above are not constant in F,
    // its time derivative is non zero. Here we use the fact that for
    // a vector r_B_A fixed in a moving frame B but expressed in another frame A,
    // its time derivative in A is the angular velocity of B in A crossed with
    // the vector, i.e., d_A/dt r_B_A = w_AB % r_B_A.
    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&               H_FM_Dot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        const Rotation& R_FM = getX_FM(pc).R();
        const Vec3&     Mx_F = R_FM.x(); // M's x axis, expressed in F
        const Vec3&     My_F = R_FM.y(); // M's y axis, expressed in F

        const Vec3&     w_FM = getV_FM(vc)[0]; // angular velocity of M in F

        H_FM_Dot[0] = SpatialRow( ~(w_FM % Mx_F), Row3(0) );
        H_FM_Dot[1] = SpatialRow( ~(w_FM % My_F), Row3(0) );

        // For translation in F.
        H_FM_Dot[2] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[3] = SpatialRow( Row3(0), Row3(0) );
        H_FM_Dot[4] = SpatialRow( Row3(0), Row3(0) );
    }

    // CAUTION: we do not zero the unused 4th element of q for Euler angles; it
    // is up to the caller to do that if it is necessary.
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                          bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Model);
        assert(q && in && out);

        if (useEulerAnglesIfPossible) {
            const Mat32    Q = Rotation::calcQBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                    .getSubMat<3,2>(0,0); // drop 3rd column
            if (matrixOnRight) {
                Row2::updAs(out)   = Row3::getAs(in) * Q;
                Row3::updAs(out+2) = Row3::getAs(in+3);// translational part of Q block is identity
            } else {
                Vec3::updAs(out)   = Q * Vec2::getAs(in);        
                Vec3::updAs(out+3) = Vec3::getAs(in+2);// translational part of Q block is identity
            }

        } else {
            // Quaternion: Q block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat42 Q = (Rotation::calcUnnormalizedQBlockForQuaternion(Vec4::getAs(q))*R_FM)
                                .getSubMat<4,2>(0,0); // drop 3rd column
            if (matrixOnRight) {
                Row2::updAs(out)   = Row4::getAs(in) * Q;
                Row3::updAs(out+2) = Row3::getAs(in+4); // translational part of Q block is identity
            } else { // matrix on left
                Vec4::updAs(out)   = Q * Vec2::getAs(in);
                Vec3::updAs(out+4) = Vec3::getAs(in+2); // translational part of Q block is identity
            }
        }
    }

    // Compute out_u = inv(Q) * in_q
    //   or    out_q = in_u * inv(Q)
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q,
                             bool matrixOnRight, const Real* in, Real* out) const
    {
        assert(sbs.getStage() >= Stage::Position);
        assert(in && out);

        if (useEulerAnglesIfPossible) {
            const Mat23    QInv = Rotation::calcQInvBlockForBodyXYZInBodyFrame(Vec3::getAs(q))
                                    .getSubMat<2,3>(0,0);   // drop 3rd row
            if (matrixOnRight) {
                Row3::updAs(out)   = Row2::getAs(in) * QInv;
                Row3::updAs(out+3) = Row3::getAs(in+2); // translational part of QInv block is identity
            } else {
                Vec2::updAs(out)   = QInv * Vec3::getAs(in);
                Vec3::updAs(out+2) = Vec3::getAs(in+3); // translational part of QInv block is identity
            }
        } else {           
            // Quaternion: QInv block is only available expecting angular velocity in the
            // parent frame F, but we have it in M for this joint.
            const Rotation R_FM(Quaternion(Vec4::getAs(q)));
            const Mat24 QInv = (~R_FM*Rotation::calcUnnormalizedQInvBlockForQuaternion(Vec4::getAs(q)))
                                    .getSubMat<2,4>(0,0);   // drop 3rd row
            if (matrixOnRight) {
                Row4::updAs(out)   = Row2::getAs(in) * QInv;
                Row3::updAs(out+4) = Row3::getAs(in+2); // translational part of QInv block is identity
            } else { // matrix on left
                Vec2::updAs(out)   = QInv * Vec4::getAs(in);
                Vec3::updAs(out+2) = Vec3::getAs(in+4); // translational part of QInv block is identity
            }
        }
    }

    void calcQDot(
        const SBStateDigest&   sbs,
        const Vector&          u,
        Vector&                qdot) const
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3  w_FM_M = Vec3(fromU(u)[0], fromU(u)[1], 0); // Angular velocity in M
        const Vec3& v_FM   = fromUVec3(u,2);                    // Linear velocity in F

        if (getUseEulerAngles(mv)) {
            const Vec3& theta = fromQVec3(sbs.getQ(),0); // Euler angles
            toQVec3(qdot,0) = Rotation::convertAngVelToBodyFixed123Dot(theta,
                                            w_FM_M); // need w in *body*, not parent
            toQVec3(qdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdot,3) = v_FM;
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            const Vec4& quat = fromQuat(sbs.getQ());
            toQuat (qdot)   = Rotation::convertAngVelToQuaternionDot(quat,
                                            R_FM*w_FM_M); // need w in *parent* frame here
            toQVec3(qdot,4) = v_FM;
        }
    }
 
    void calcQDotDot(
        const SBStateDigest&   sbs,
        const Vector&          udot, 
        Vector&                qdotdot) const 
    {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc = sbs.getPositionCache();
        const Vec3  w_FM_M     = Vec3(fromU(sbs.getU())[0], fromU(sbs.getU())[1], 0); // Angular velocity of M in F, exp. in M
        const Vec3& v_FM       = fromUVec3(sbs.getU(),2); // linear velocity of M in F, expressed in M
        const Vec3  w_FM_M_dot = Vec3(fromU(udot)[0], fromU(udot)[1], 0);
        const Vec3& v_FM_dot   = fromUVec3(udot,2);

        if (getUseEulerAngles(mv)) {
            const Vec3& theta  = fromQVec3(sbs.getQ(),0); // Euler angles
            toQVec3(qdotdot,0) = Rotation::convertAngVelDotToBodyFixed123DotDot
                                             (theta, w_FM_M, w_FM_M_dot); // needed in body frame here
            toQVec3(qdotdot,4) = Vec3(0); // TODO: kludge, clear unused element
            toQVec3(qdotdot,3) = v_FM_dot;
        } else {
            const Rotation& R_FM = getX_FM(pc).R();
            const Vec4& quat  = fromQuat(sbs.getQ());
            toQuat(qdotdot)   = Rotation::convertAngVelDotToQuaternionDotDot
                                             (quat,R_FM*w_FM_M,R_FM*w_FM_M_dot); // needed in parent frame
            toQVec3(qdotdot,4) = v_FM_dot;
        }
    }

    void copyQ(const SBModelVars& mv, const Vector& qIn, Vector& q) const {
        if (getUseEulerAngles(mv)) {
            toQVec3(q,0) = fromQVec3(qIn,0); // euler angles
            toQVec3(q,3) = fromQVec3(qIn,3); // translations
        } else {
            toQuat(q)    = fromQuat(qIn);    // quaternion
            toQVec3(q,4) = fromQVec3(qIn,4); // translations
        }
    }

    int  getMaxNQ()                   const {return 7;}
    int  getNQInUse(const SBModelVars& mv) const {return getUseEulerAngles(mv) ? 6 : 7;} 
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (getUseEulerAngles(sbs.getModelVars())) {startOfQuaternion.invalidate(); return false;}
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
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

    bool enforceQuaternionConstraints(
        const SBStateDigest& sbs, 
        Vector&             q,
        Vector&             qErrest) const 
    {
        if (getUseEulerAngles(sbs.getModelVars())) 
            return false; // no change

        Vec4& quat = toQuat(q);
        quat = quat / quat.norm();

        if (qErrest.size()) {
            Vec4& qerr = toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }

        return true;
    }

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = Vec3(0); // clear unused element
        toQVec3(outputQ, 3) = fromQVec3(inputQ, 4);
        toQVec3(outputQ, 0) = Rotation(Quaternion(fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
    }
    
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        toQVec3(outputQ, 4) = fromQVec3(inputQ, 3);
        Rotation rot;
        rot.setRotationToBodyFixedXYZ(fromQVec3(inputQ, 0));
        toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
    }
};


    // WELD //

// This is a "joint" with no degrees of freedom, that simply forces
// the two reference frames to be identical.
class RBNodeWeld : public RBGroundBody {
public:
    const char* type() { return "weld"; }

    RBNodeWeld(const MassProperties& mProps_B, const Transform& X_PF, const Transform& X_BM) : RBGroundBody(mProps_B, X_PF, X_BM) {
    }

    void realizePosition(SBStateDigest& sbs) const {
        SBPositionCache& pc = sbs.updPositionCache();

        const Transform& X_BM = getX_BM();   // fixed
        const Transform& X_PF = getX_PF();   // fixed
        const Transform& X_GP = getX_GP(pc); // already calculated

        updX_FM(pc).setToZero();
        updX_PB(pc) = X_PF * ~X_BM; // TODO: precalculate X_MB
        updX_GB(pc) = X_GP * getX_PB(pc);
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
    
    void realizeVelocity(SBStateDigest& sbs) const {
        const SBPositionCache& pc = sbs.getPositionCache();
        SBVelocityCache& vc = sbs.updVelocityCache();

        updV_FM(vc) = 0;
        updV_PB_G(vc) = 0;
        calcJointIndependentKinematicsVel(pc,vc);
    }

    void realizeDynamics(SBStateDigest& sbs) const {
        // Mobilizer-specific.
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        SBDynamicsCache& dc = sbs.updDynamicsCache();
        
        updVD_PB_G(dc) = 0;

        // Mobilizer independent.
        calcJointIndependentDynamicsVel(pc,vc,dc);
    }

    void calcArticulatedBodyInertiasInward(const SBPositionCache& pc, SBDynamicsCache& dc) const {
        updP(dc) = getMk(pc);
        updTauBar(dc)  = 1.; // identity matrix
        updPsi(dc)     = getPhi(pc) * getTauBar(dc);
    }
    
    void calcQDotDot(const SBStateDigest& sbs, const Vector& udot, Vector& qdotdot) const {
    }
    
    void calcUDotPass1Inward(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector_<SpatialVec>&        allGepsilon,
        Vector&                     allEpsilon) const 
    {
        const SpatialVec& myBodyForce  = fromB(bodyForces);
        SpatialVec&       z            = toB(allZ);
        SpatialVec&       Geps         = toB(allGepsilon);

        z = getCentrifugalForces(dc) - myBodyForce;

        for (int i=0 ; i<(int)children.size() ; i++) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

            z += phiChild * (zChild + GepsChild);
        }

        Geps = 0;
    }

    void calcUDotPass2Outward(
        const SBPositionCache& pc,
        const SBDynamicsCache& dc,
        const Vector&          allEpsilon,
        Vector_<SpatialVec>&   allA_GB,
        Vector&                allUDot) const
    {
        SpatialVec&     A_GB = toB(allA_GB);

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = parent->getNodeNum()== 0 
            ? SpatialVec(Vec3(0), Vec3(0))
            : ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP + getCoriolisAcceleration(dc);  
    }
    
    void calcMInverseFPass1Inward(
        const SBPositionCache& pc,
        const SBDynamicsCache& dc,
        const Vector&          f,
        Vector_<SpatialVec>&   allZ,
        Vector_<SpatialVec>&   allGepsilon,
        Vector&                allEpsilon) const 
    {
        SpatialVec&       z            = toB(allZ);
        SpatialVec&       Geps         = toB(allGepsilon);

        z = SpatialVec(Vec3(0), Vec3(0));

        for (int i=0 ; i<(int)children.size() ; i++) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

            z += phiChild * (zChild + GepsChild);
        }

        Geps = 0;
    }

    void calcMInverseFPass2Outward(
        const SBPositionCache& pc,
        const SBDynamicsCache& dc,
        const Vector&          allEpsilon,
        Vector_<SpatialVec>&   allA_GB,
        Vector&                allUDot) const
    {
        SpatialVec&     A_GB = toB(allA_GB);

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = parent->getNodeNum()== 0 
            ? SpatialVec(Vec3(0), Vec3(0))
            : ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP;
    }
};


/////////////////////////////////////////
// RigidBodyNode for custom mobilizers //
/////////////////////////////////////////

template <int nu>
class RBNodeCustom : public RigidBodyNodeSpec<nu> {
public:
    RBNodeCustom(const MobilizedBody::Custom::Implementation& impl, const MassProperties& mProps_B,
            const Transform& X_PF, const Transform& X_BM,UIndex& nextUSlot, USquaredIndex& nextUSqSlot, QIndex& nextQSlot) : 
            RigidBodyNodeSpec<nu>(mProps_B, X_PF, X_BM, nextUSlot, nextUSqSlot, nextQSlot, RigidBodyNode::QDotMayDifferFromU,
            impl.getImpl().getNAngles() == 4 ? RigidBodyNode::QuaternionMayBeUsed : RigidBodyNode::QuaternionIsNeverUsed),
            impl(impl), nq(impl.getImpl().getNQ()), nAngles(impl.getImpl().getNAngles()) {
        this->updateSlots(nextUSlot,nextUSqSlot,nextQSlot);
    }
    const char* type() const {
        return "custom";
    }
    int  getMaxNQ() const {
        return nq;
    }
    int getNQInUse(const SBModelVars& mv) const {
        return (nAngles == 4 && this->getUseEulerAngles(mv) ? nq-1 : nq);
    }
    virtual int getNUInUse(const SBModelVars& mv) const {
        return nu;
    }
    bool isUsingQuaternion(const SBStateDigest& sbs, MobilizerQIndex& startOfQuaternion) const {
        if (nAngles < 4 || this->getUseEulerAngles(sbs.getModelVars())) {
            startOfQuaternion.invalidate();
            return false;
        }
        startOfQuaternion = MobilizerQIndex(0); // quaternion comes first
        return true;
    }
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& startOfAngles, int& numAngles) const {
        if (nAngles == 0 || (nAngles == 4 && !this->getUseEulerAngles(sbs.getModelVars()))) {
            startOfAngles.invalidate();
            numAngles = 0;
            return false;
        }
        startOfAngles = MobilizerQIndex(0);
        numAngles = std::min(nAngles, 3); 
        return true;
    }
    void copyQ(const SBModelVars& mv, const Vector& qIn, Vector& q) const {
        const int n = getNQInUse(mv);
        for (int i = 0; i < n; ++i)
            q[i] = qIn[i];
    }
    void calcLocalQDotFromLocalU(const SBStateDigest& sbs, const Real* u, Real* qdot) const {
        impl.multiplyByQMatrix(sbs.getState(), false, nu, u, getNQInUse(sbs.getModelVars()), qdot);
    }
    void calcLocalQDotDotFromLocalUDot(const SBStateDigest& sbs, const Real* udot, Real* qdotdot) const {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const int nqInUse = getNQInUse(sbs.getModelVars());
        const Real* u = &sbs.getU()[this->getUIndex()];
        impl.multiplyByQMatrix(sbs.getState(), false, nu, udot, nqInUse, qdotdot);
        Real temp[7];
        impl.multiplyByQDotMatrix(sbs.getState(), false, nu, u, nqInUse, temp);
        for (int i = 0; i < nqInUse; ++i)
            qdotdot[i] += temp[i];
    }
    void multiplyByQBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q, bool matrixOnRight, 
                                  const Real* in, Real* out) const {
        const SBModelVars& mv = sbs.getModelVars();
        int nIn, nOut;
        if (matrixOnRight) {
            nIn = getNQInUse(mv);
            nOut = getNUInUse(mv);
        }
        else {
            nIn = getNUInUse(mv);
            nOut = getNQInUse(mv);
        }
        impl.multiplyByQMatrix(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
    }
    void multiplyByQInvBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q, bool matrixOnRight,
                                     const Real* in, Real* out) const {
        const SBModelVars& mv = sbs.getModelVars();
        int nIn, nOut;
        if (matrixOnRight) {
            nIn = getNUInUse(mv);
            nOut = getNQInUse(mv);
        }
        else {
            nIn = getNQInUse(mv);
            nOut = getNUInUse(mv);
        }
        impl.multiplyByQInverse(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
    }
    void multiplyByQDotBlock(const SBStateDigest& sbs, bool useEulerAnglesIfPossible, const Real* q, const Real* u,
                                     bool matrixOnRight, const Real* in, Real* out) const {
        const SBModelVars& mv = sbs.getModelVars();
        int nIn, nOut;
        if (matrixOnRight) {
            nIn = getNQInUse(mv);
            nOut = getNUInUse(mv);
        }
        else {
            nIn = getNUInUse(mv);
            nOut = getNQInUse(mv);
        }
        impl.multiplyByQDotMatrix(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
    }

    void calcQDot(const SBStateDigest& sbs, const Vector& u, Vector& qdot) const {
        const int nqInUse = getNQInUse(sbs.getModelVars());
        const int qindex = this->getQIndex();
        impl.multiplyByQMatrix(sbs.getState(), false, nu, &u[this->getUIndex()], nqInUse, &qdot[qindex]);
        for (int i = nqInUse; i < nq; ++i)
            qdot[qindex+i] = 0.0;
    }

    void calcQDotDot(const SBStateDigest& sbs, const Vector& udot, Vector& qdotdot) const {
        const SBModelVars& mv = sbs.getModelVars();
        const SBPositionCache& pc   = sbs.getPositionCache();
        const int nqInUse = getNQInUse(sbs.getModelVars());
        const int qindex = this->getQIndex();
        const Real* u = &sbs.getU()[this->getUIndex()];
        impl.multiplyByQMatrix(sbs.getState(), false, nu, &udot[this->getUIndex()], nqInUse, &qdotdot[qindex]);
        Real temp[7];
        impl.multiplyByQDotMatrix(sbs.getState(), false, nu, u, nqInUse, temp);
        for (int i = 0; i < nqInUse; ++i)
            qdotdot[qindex+i] += temp[i];
        for (int i = nqInUse; i < nq; ++i)
            qdotdot[qindex+i] = 0.0;
    }
    bool enforceQuaternionConstraints(const SBStateDigest& sbs, Vector& q, Vector& qErrest) const {
        if (nAngles != 4 || this->getUseEulerAngles(sbs.getModelVars())) 
            return false;
        Vec4& quat = this->toQuat(q);
        quat = quat / quat.norm();
        if (qErrest.size()) {
            Vec4& qerr = this->toQuat(qErrest);
            qerr -= dot(qerr,quat) * quat;
        }
        return true;
    }
    
    // Convert from quaternion to Euler angle representations.
    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
        int indexBase = this->getQIndex();
        if (nAngles != 4) {
            for (int i = 0; i < nq; ++i)
                outputQ[indexBase+i] = inputQ[indexBase+i];
        }
        else {
            this->toQVec3(outputQ, 0) = Rotation(Quaternion(this->fromQuat(inputQ))).convertRotationToBodyFixedXYZ();
            for (int i = 3; i < nq-1; ++i)
                outputQ[indexBase+i] = inputQ[indexBase+i+1];
            outputQ[indexBase+nq-1] = 0.0;
        }
    }
    // Convert from Euler angle to quaternion representations.
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
        int indexBase = this->getQIndex();
        if (nAngles != 4) {
            for (int i = 0; i < nq; ++i)
                outputQ[indexBase+i] = inputQ[indexBase+i];
        }
        else {
            Rotation rot;
            rot.setRotationToBodyFixedXYZ(Vec3(inputQ[indexBase], inputQ[indexBase+1], inputQ[indexBase+2]));
            this->toQuat(outputQ) = rot.convertRotationToQuaternion().asVec4();
            for (int i = 4; i < nq; ++i)
                outputQ[indexBase+i] = inputQ[indexBase+i-1];
        }
    };

    void setQToFitTransform(const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {
        impl.setQToFitTransform(sbs.getState(), X_FM, this->getNQInUse(sbs.getModelVars()), &q[this->getQIndex()]);
    }
    void setQToFitRotation(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        setQToFitTransform(sbs, Transform(R_FM), q);
    }
    void setQToFitTranslation(const SBStateDigest& sbs, const Vec3& T_FM, Vector& q) const {
        setQToFitTransform(sbs, Transform(T_FM), q);
    }

    void setUToFitVelocity(const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {
        impl.setUToFitVelocity(sbs.getState(), V_FM, nu, &u[this->getUIndex()]);
    }
    void setUToFitAngularVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {
        setUToFitVelocity(sbs, q, SpatialVec(w_FM, Vec3(0)), u);
    }
    void setUToFitLinearVelocity(const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const {
        setUToFitVelocity(sbs, q, SpatialVec(Vec3(0), v_FM), u);
    }

        // VIRTUAL METHODS FOR SINGLE-NODE OPERATOR CONTRIBUTIONS //

    void realizeModel(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeModel(sbs);
        impl.realizeModel(sbs.updState());
    }

    void realizeInstance(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeInstance(sbs);
        impl.realizeInstance(sbs.getState());
    }

    void realizeTime(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeTime(sbs);
        impl.realizeTime(sbs.getState());
    }

    void realizePosition(SBStateDigest& sbs) const {
        impl.realizePosition(sbs.getState());
        RigidBodyNodeSpec<nu>::realizePosition(sbs);
    }

    void realizeVelocity(SBStateDigest& sbs) const {
        impl.realizeVelocity(sbs.getState());
        RigidBodyNodeSpec<nu>::realizeVelocity(sbs);
    }

    void realizeDynamics(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeDynamics(sbs);
        impl.realizeDynamics(sbs.getState());
    }

    void realizeAcceleration(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeAcceleration(sbs);
        impl.realizeAcceleration(sbs.getState());
    }

    void realizeReport(SBStateDigest& sbs) const {
        RigidBodyNodeSpec<nu>::realizeReport(sbs);
        impl.realizeReport(sbs.getState());
    }

    void getInternalForce(const SBAccelerationCache& ac, Vector& tau) const {
        assert(false);
    }

    void calcJointSinCosQNorm(
        const SBModelVars&  mv, 
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const {
        
    }
    
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_FM) const {
        int nq = getNQInUse(sbs.getModelVars());
        if (nAngles == 4 && !this->getUseEulerAngles(sbs.getModelVars())) {
            Vec<nu+1> localQ = Vec<nu+1>::getAs(&q[this->getQIndex()]);
            Vec4::updAs(&localQ[0]) = Vec4::getAs(&localQ[0]).normalize(); // Normalize the quaternion
            X_FM = impl.calcMobilizerTransformFromQ(sbs.getState(), nq, &(localQ[0]));
        }
        else
            X_FM = impl.calcMobilizerTransformFromQ(sbs.getState(), nq, &(q[this->getQIndex()]));
    }
    
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        Mat<nu,2,Row3,1,2>&  H_FM) const {
        for (int i = 0; i < nu; ++i) {
            Vec<nu> u(0);
            u[i] = 1;
            H_FM.row(i) = ~impl.multiplyByHMatrix(sbs.getState(), nu, &u[0]);
        }
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        Mat<nu,2,Row3,1,2>&  H_FM_Dot) const {
        for (int i = 0; i < nu; ++i) {
            Vec<nu> u(0);
            u[i] = 1;
            H_FM_Dot.row(i) = ~impl.multiplyByHDotMatrix(sbs.getState(), nu, &u[0]);
        }
    }

private:
    const MobilizedBody::Custom::Implementation& impl;
    const int nq, nAngles;
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
    const HType& H_FM = getH_FM(pc);

    // want r_MB_F, that is, the vector from OM to OB, expressed in F
    const Vec3&     r_MB   = getX_MB().T();    // fixed
    const Rotation& R_FM   = getX_FM(pc).R(); // just calculated
    const Vec3      r_MB_F = R_FM*r_MB;

    // Reminder: round brackets () applied to a matrix select columns.
    HType H_MB;
    H_MB(0) = Row3(0); // fills column with zero
    H_MB(1) = H_FM(0) * crossMat(r_MB_F);

    // Now we want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();      // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GF = R_GP * R_PF;

    H_PB_G = (H_FM + H_MB) * ~R_GF;
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
    const HType& H_FM     = getH_FM(pc);
    const HType& H_FM_Dot = getH_FM_Dot(dc);

    HType H_MB, H_MB_Dot;

    // want r_MB_F, that is, the vector from OM to OB, expressed in F
    const Vec3&     r_MB   = getX_MB().T();    // fixed
    const Rotation& R_FM   = getX_FM(pc).R(); // just calculated
    const Vec3      r_MB_F = R_FM*r_MB;

    const Vec3& w_FM = getV_FM(vc)[0]; // local angular velocity

    // Reminder: round brackets () applied to a matrix select columns.
    H_MB(0) = Row3(0); // fills column with zero
    H_MB(1) = H_FM(0) * crossMat(r_MB_F);

    H_MB_Dot(0) = Row3(0);
    H_MB_Dot(1) =   H_FM_Dot(0) * crossMat(r_MB_F) 
                  + H_FM(0)     * crossMat(w_FM % r_MB_F);

    // Now we want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();      // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GF = R_GP * R_PF;

    const Vec3& w_GF = getV_GP(vc)[0]; // F and P have same angular velocity

    // Note: time derivative of R_GF is crossMat(w_GF)*R_GF, so derivative
    // of ~R_GF is ~(crossMat(w_GF)*R_GF) = - ~R_GF*crossMat(w_GF) since
    // crossMat's are skew-symmetric.
    //      H_PB_G = (H_FM + H_MB) * ~R_GF (see above method)
    const HType& H_PB_G = getH(pc);
    H_PB_G_Dot =  (H_FM_Dot + H_MB_Dot) * ~R_GF
                 - H_PB_G * crossMat(w_GF);
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

    // TODO: change sign on taubar to make it GH-I instead, which only requires
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
    const SBStateDigest&       sbs,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces) const 
{
    const SBPositionCache& pc = sbs.getPositionCache();
    const SBDynamicsCache& dc = sbs.getDynamicsCache();
    SBAccelerationCache&   ac = sbs.updAccelerationCache();

    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - fromB(bodyForces);

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(ac);
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& GepsChild = children[i]->getGepsilon(ac);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(ac)  = fromU(mobilityForces) - getH(pc)*z;
    updNu(ac)       = getDI(dc) * getEpsilon(ac);
    updGepsilon(ac) = getG(dc)  * getEpsilon(ac);
}
/*
template<int dof> void
RigidBodyNodeSpec<dof>::calcZ(
    const SBPositionCache&     pc,
    const SBDynamicsCache&     dc,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces,
    SBAccelerationCache&       ac) const 
{
    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - fromB(bodyForces);

    for (int i=0 ; i<(int)children.size() ; i++) {
        const SpatialVec& zChild    = children[i]->getZ(ac);
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& GepsChild = children[i]->getGepsilon(ac);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(ac)  = fromU(mobilityForces) - getH(pc)*z;
    updNu(ac)       = getDI(dc) * getEpsilon(ac);
    updGepsilon(ac) = getG(dc)  * getEpsilon(ac);
}
*/
//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to calcZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcAccel(
    const SBStateDigest&   sbs,
    Vector&                allUdot,
    Vector&                allQdotdot) const 
{
    const SBPositionCache& pc = sbs.getPositionCache();
    const SBDynamicsCache& dc = sbs.getDynamicsCache();
    SBAccelerationCache&   ac = sbs.updAccelerationCache();
    Vec<dof>&        udot   = toU(allUdot);
    const SpatialVec alphap = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

    udot        = getNu(ac) - (~getG(dc)*alphap);
    updA_GB(ac) = alphap + ~getH(pc)*udot + getCoriolisAcceleration(dc);  

    calcQDotDot(sbs, allUdot, allQdotdot);  
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
// Calculate product of kinematic Jacobian J=~Phi*~H and a mobility-space vector. Requires 
// that Phi and H are available, so this should only be called in Stage::Position or higher.
// This does not change the cache at all.
//
// Note that if the vector v==u (generalized speeds) then the result is V_GB=J*u, the
// body spatial velocities generated by those speeds.
//
// Call base to tip (outward).
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcSpatialKinematicsFromInternal(
    const SBPositionCache&      pc,
    const Vector&               v,
    Vector_<SpatialVec>&        Jv) const
{
    const Vec<dof>& in  = fromU(v);
    SpatialVec&     out = toB(Jv);

    // Shift parent's result outward (ground result is 0).
    const SpatialVec outP = parent->getNodeNum()== 0 
        ? SpatialVec(Vec3(0), Vec3(0))
        : ~getPhi(pc) * parent->fromB(Jv);

    out = outP + ~getH(pc)*in;  
}

//
// Calculate product of kinematic Jacobian transpose ~J=H*Phi and a gradient vector on each of the
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
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);

        z += phiChild * zChild; 
    }

    eps  = getH(pc) * z;
}


// The Ground node is special because it doesn't need a mobilizer.
/*static*/ RigidBodyNode*
RigidBodyNode::createGroundNode() {
    return new RBGroundBody(MassProperties(Infinity, Vec3(0), Infinity*Inertia(1)), Transform(), Transform());
}


    ///////////////////////////////////////////////////////////////////////
    // Implementation of MobilizedBodyImpl createRigidBodyNode() methods //
    ///////////////////////////////////////////////////////////////////////

RigidBodyNode* MobilizedBody::PinImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeTorsion(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::SliderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeSlider(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::UniversalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeUJoint(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::CylinderImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeCylinder(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::BendStretchImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeBendStretch(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::PlanarImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodePlanar(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::GimbalImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeGimbal(getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::BallImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeBall(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::EllipsoidImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeEllipsoid(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        getDefaultRadii(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::TranslationImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeTranslate(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::FreeImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeFree(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::LineOrientationImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeLineOrientation(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::FreeLineImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeFreeLine(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot,nextUSqSlot,nextQSlot);
}


RigidBodyNode* MobilizedBody::ScrewImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeScrew(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        getDefaultPitch(),nextUSlot,nextUSqSlot,nextQSlot);
}

RigidBodyNode* MobilizedBody::WeldImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeWeld(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame());
}


RigidBodyNode* MobilizedBody::GroundImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return RigidBodyNode::createGroundNode();
}

RigidBodyNode* MobilizedBody::CustomImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    switch (getImplementation().getImpl().getNU()) {
    case 1:
        return new RBNodeCustom<1>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    case 2:
        return new RBNodeCustom<2>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    case 3:
        return new RBNodeCustom<3>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    case 4:
        return new RBNodeCustom<4>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    case 5:
        return new RBNodeCustom<5>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    case 6:
        return new RBNodeCustom<6>(getImplementation(), getDefaultRigidBodyMassProperties(),
            getDefaultInboardFrame(), getDefaultOutboardFrame(), nextUSlot, nextUSqSlot, nextQSlot);
    default:
        assert(!"Illegal number of degrees of freedom for custom MobilizedBody");
        return 0;
    }
}



