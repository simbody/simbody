/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *    Charles Schwieters (NIH): wrote the public domain IVM code from which   *
 *                              this was derived.                             *
 *    Peter Eastman wrote the Weld joint.                                     *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**@file
 * This file contains the implementations for RigidBodyNodes which have
 * no degrees of freedom -- Ground and Weld mobilizers. These cannot be
 * derived from the usual RigidBodyNodeSpec<dof> class because dof==0
 * is problematic there. Also, these can have very efficient
 * implementations here since they know they have no dofs.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "MobilizedBodyImpl.h"


////////////////////////////////////////////////
// Define classes derived from RigidBodyNode. //
////////////////////////////////////////////////

/**
 * This still-abstract class is the common base for any MobilizedBody which
 * has no mobilities. Currently that is only the unique Ground body and
 * MobilizedBody::Weld, but it is conceivable that others could crop up.
 *
 * The base class overrides all the virtual methods which deal in q's and u's
 * so have common, trivial implementations if there are no q's and no u's.
 */
class ImmobileRigidBodyNode : public RigidBodyNode {
public:
    ImmobileRigidBodyNode(const MassProperties& mProps_B, const Transform& X_PF, const Transform& X_BM,
        const UIndex& uIx, const USquaredIndex& usqIx, const QIndex& qIx)
    :   RigidBodyNode(mProps_B, X_PF, X_BM, QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed)
    {
        uIndex   = uIx;
        uSqIndex = usqIx;
        qIndex   = qIx;
    }
    ~ImmobileRigidBodyNode() {}

    int  getDOF()   const {return 0;}
    int  getMaxNQ() const {return 0;}
    int  getNUInUse(const SBModelVars&) const {return 0;}
    int  getNQInUse(const SBModelVars&) const {return 0;}
    bool isUsingQuaternion(const SBStateDigest&, MobilizerQIndex& ix) const
    {   ix.invalidate(); return false; }
    int getPosPoolSize(const SBStateDigest&) const {return 0;}
    int getVelPoolSize(const SBStateDigest&) const {return 0;}

    int calcQPoolSize(const SBModelVars&) const {return 0;}

    void performQPrecalculations(const SBStateDigest& sbs,
                                 const Real* q, int nq,
                                 Real* qCache,  int nQCache,
                                 Real* qErr,    int nQErr) const
    {   assert(nq==0 && nQCache==0 && nQErr==0); }


    // An immobile mobilizer holds the mobilized body's M frame coincident
    // with the parent body's F frame forever.
    void calcX_FM(const SBStateDigest& sbs,
                  const Real* q,      int nq,
                  const Real* qCache, int nQCache,
                  Transform&  X_FM) const
    {   assert(nq==0 && nQCache==0);
        X_FM.setToZero(); }

    void setQToFitTransformImpl
       (const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {}
    void setQToFitRotationImpl
       (const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {}
    void setQToFitTranslationImpl
       (const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {}

    void setUToFitVelocityImpl
       (const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {}
    void setUToFitAngularVelocityImpl
       (const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {}
    void setUToFitLinearVelocityImpl
       (const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const {}


    void multiplyByN(const SBStateDigest&, bool matrixOnRight,
                     const Real* in, Real* out) const {}
    void multiplyByNInv(const SBStateDigest&, bool matrixOnRight,
                        const Real* in, Real* out) const {}
    void multiplyByNDot(const SBStateDigest&, bool matrixOnRight,
                        const Real* in, Real* out) const {}

    void calcQDot(const SBStateDigest&, const Real* udot, Real* qdotdot) const {}
    void calcQDotDot(const SBStateDigest&, const Real* udot, Real* qdotdot) const {}

    bool enforceQuaternionConstraints(
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const {return false;}

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {}
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {}

    void setVelFromSVel(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const SpatialVec&           sVel,
        Vector&                     u) const {}
};

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class RBGroundBody : public ImmobileRigidBodyNode {
public:
    RBGroundBody()
    :   ImmobileRigidBodyNode(MassProperties(Infinity, Vec3(0), Inertia(Infinity)),
                              Transform(), Transform(),
                              UIndex(0), USquaredIndex(0), QIndex(0)) {}

    const char* type() const { return "ground"; }

    // TODO: should ground set the various cache entries here?
    void realizeModel   (SBStateDigest&) const {}
    void realizeInstance(const SBStateDigest& sbs) const {
        // Initialize cache entries that will never be changed at later stages.

        SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();
        SBDynamicsCache& dc = sbs.updDynamicsCache();
        SBTreeAccelerationCache& ac = sbs.updTreeAccelerationCache();
        updY(dc) = SpatialMat(Mat33(0));
        updA_GB(ac) = 0;
    }
    void realizePosition(const SBStateDigest&) const {}
    void realizeVelocity(const SBStateDigest&) const {}
    void realizeDynamics(const SBArticulatedBodyInertiaCache&, const SBStateDigest&) const {}
    // There is no realizeAcceleration().
    void realizeReport  (const SBStateDigest&) const {}

    // Ground's "composite" body inertia is still the infinite mass
    // and inertia it started with; no need to look at the children.
    // This overrides the base class default implementation.
    void calcCompositeBodyInertiasInward(
        const SBTreePositionCache&                  pc,
        Array_<SpatialInertia,MobilizedBodyIndex>&  R) const
    {   R[GroundIndex] = SpatialInertia(Infinity, Vec3(0), UnitInertia(1)); }

    // Ground's "articulated" body inertia is still the infinite mass and
    // inertia it started with; no need to look at the children.
    void realizeArticulatedBodyInertiasInward(
        const SBInstanceCache&,
        const SBTreePositionCache&,
        SBArticulatedBodyInertiaCache& abc) const
    {   ArticulatedInertia& P = updP(abc);
        P = ArticulatedInertia(SymMat33(Infinity), Mat33(Infinity),
                                       SymMat33(0));
        updPPlus(abc) = P;
    }

    void realizeYOutward(
        const SBInstanceCache&,
        const SBTreePositionCache&,
        const SBArticulatedBodyInertiaCache&,
        SBDynamicsCache&                        dc) const
    {
    }


    // Treat Ground as though welded to the universe at the ground
    // origin. The reaction there collects the effects of all the
    // base bodies and of any forces applied directly to Ground.
    void calcUDotPass1Inward(
        const SBInstanceCache&     ic,
        const SBTreePositionCache& pc,
        const SBArticulatedBodyInertiaCache&,
        const SBDynamicsCache&,
        const Real*                jointForces,
        const SpatialVec*          bodyForces,
        const Real*                allUDot,
        SpatialVec*                allZ,
        SpatialVec*                allZPlus,
        Real*                      allEpsilon) const
    {
        const SpatialVec& F            = bodyForces[0];
        SpatialVec&       z            = allZ[0];
        SpatialVec&       zPlus        = allZPlus[0];

        z = -F;

        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild   = children[i]->getPhi(pc);
            const SpatialVec& zPlusChild = allZPlus[children[i]->getNodeNum()];
            z += phiChild * zPlusChild; // 18 flops
        }

        zPlus = z;
    }

    void calcUDotPass2Outward(
        const SBInstanceCache&,
        const SBTreePositionCache&,
        const SBArticulatedBodyInertiaCache&,
        const SBTreeVelocityCache&,
        const SBDynamicsCache&,
        const Real*                epsilonTmp,
        SpatialVec*                allA_GB,
        Real*                      allUDot,
        Real*                      allTau) const
    {
        allA_GB[0] = 0;
    }

    // Ground doesn't contribute to M^-1*f. Inward pass does nothing since
    // Ground can't be the child of any body.
    void multiplyByMInvPass1Inward(
        const SBInstanceCache&     ic,
        const SBTreePositionCache& pc,
        const SBArticulatedBodyInertiaCache&,
        const Real*                f,
        SpatialVec*                allZ,
        SpatialVec*                allZPlus,
        Real*                      allEpsilon) const override
    {
    }

    // Outward pass must make sure A_GB[0] is zero so it can be propagated
    // outwards properly.
    void multiplyByMInvPass2Outward(
        const SBInstanceCache&,
        const SBTreePositionCache&,
        const SBArticulatedBodyInertiaCache&,
        const Real*                 epsilonTmp,
        SpatialVec*                 allA_GB,
        Real*                       allUDot) const override
    {
        allA_GB[0] = 0;
    }

    // Also serves as pass 1 for inverse dynamics.
    void calcBodyAccelerationsFromUdotOutward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const
    {
        allA_GB[0] = 0;
    }

    // Here Ground is the last body processed. Although it has no mobility forces
    // we can still collect up all the forces from the base bodies to Ground
    // in case anyone cares.
    void calcInverseDynamicsPass2Inward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const SpatialVec*           allA_GB,
        const Real*                 jointForces,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allF,
        Real*                       allTau) const
    {
        allF[0] = -bodyForces[0];

        // Add in forces on base bodies, shifted to Ground.
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
            allF[0] += phiChild * FChild;
        }

        // no taus
    }

    void multiplyByMPass1Outward(
        const SBTreePositionCache&  pc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const
    {
        allA_GB[0] = 0;
    }

    void multiplyByMPass2Inward(
        const SBTreePositionCache&  pc,
        const SpatialVec*           allA_GB,
        SpatialVec*                 allF,
        Real*                       allTau) const
    {
        allF[0] = 0;

        // Add in forces on base bodies, shifted to Ground.
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
            allF[0] += phiChild * FChild;
        }

        // no taus
    }


    void multiplyBySystemJacobian(
        const SBTreePositionCache&  pc,
        const Real*                 v,
        SpatialVec*                 Jv) const
    {
        Jv[0] = SpatialVec(Vec3(0));
    }

    void multiplyBySystemJacobianTranspose(
        const SBTreePositionCache&  pc,
        SpatialVec*                 zTmp,
        const SpatialVec*           X,
        Real*                       JtX) const
    {
        zTmp[0] = X[0];
        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild = children[i]->getPhi(pc);
            zTmp[0] += phiChild * zChild;
        }
        // No generalized speeds so no contribution to JtX.
    }

    void calcEquivalentJointForces(
        const SBTreePositionCache&  pc,
        const SBDynamicsCache&,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allZ,
        Real*                       jointForces) const
    {
        allZ[0] = bodyForces[0];
        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            allZ[0] += phiChild * zChild;
        }
    }
};

    // WELD //

// This is a "joint" with no degrees of freedom, that simply forces
// the two reference frames to be identical. A Weld node always has a parent
// but has no q's and no u's.
class RBNodeWeld : public ImmobileRigidBodyNode {
public:
    RBNodeWeld(const MassProperties& mProps_B, const Transform& X_PF, const Transform& X_BM,
        const UIndex& uIx, const USquaredIndex& usqIx, const QIndex& qIx)
    :   ImmobileRigidBodyNode(mProps_B, X_PF, X_BM, uIx, usqIx, qIx) {}

    const char* type() { return "weld"; }

    void realizeModel(SBStateDigest& sbs) const {}
    void realizeInstance(const SBStateDigest& sbs) const {
        // Initialize cache entries that will never be changed at later stages.

        SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();
        SBTreeAccelerationCache& ac = sbs.updTreeAccelerationCache();
        updV_FM(vc) = 0;
        updV_PB_G(vc) = 0;
        updVD_PB_G(vc) = 0;
    }

    void realizePosition(const SBStateDigest& sbs) const {
        SBTreePositionCache& pc = sbs.updTreePositionCache();

        const Transform& X_MB = getX_MB();   // fixed
        const Transform& X_PF = getX_PF();   // fixed
        const Transform& X_GP = getX_GP(pc); // already calculated

        updX_FM(pc).setToZero();
        updX_PB(pc) = X_PF * X_MB;
        updX_GB(pc) = X_GP * getX_PB(pc);
        const Vec3 p_PB_G = getX_GP(pc).R() * getX_PB(pc).p();

        // The Phi matrix conveniently performs child-to-parent (inward) shifting
        // on spatial quantities (forces); its transpose does parent-to-child
        // (outward) shifting for velocities.
        updPhi(pc) = PhiMatrix(p_PB_G);

        // Calculate spatial mass properties. That means we need to transform
        // the local mass moments into the Ground frame and reconstruct the
        // spatial inertia matrix Mk.

        const Rotation& R_GB = getX_GB(pc).R();
        const Vec3&     p_GB = getX_GB(pc).p();

        // reexpress inertia in ground (57 flops)
        const UnitInertia G_Bo_G  = getUnitInertia_OB_B().reexpress(~R_GB);
        const Vec3        p_BBc_G = R_GB*getCOM_B(); // 15 flops

        updCOM_G(pc) = p_GB + p_BBc_G; // 3 flops

        // Calc Mk: the spatial inertia matrix about the body origin.
        // Note: we need to calculate this now so that we'll be able to calculate
        // kinetic energy without going past the Velocity stage.
        updMk_G(pc) = SpatialInertia(getMass(), p_BBc_G, G_Bo_G);
    }

    void realizeVelocity(const SBStateDigest& sbs) const {
        const SBTreePositionCache& pc = sbs.getTreePositionCache();
        SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();
        calcJointIndependentKinematicsVel(pc,vc);
    }

    void realizeDynamics(const SBArticulatedBodyInertiaCache&   abc,
                         const SBStateDigest&                   sbs) const {
        // Mobilizer-specific.
        const SBTreePositionCache&  pc = sbs.getTreePositionCache();
        const SBTreeVelocityCache&  vc = sbs.getTreeVelocityCache();
        SBDynamicsCache&            dc = sbs.updDynamicsCache();

        // Mobilizer independent.
        calcJointIndependentDynamicsVel(pc,abc,vc,dc);
    }

    // There is no realizeAcceleration().

    void realizeReport(const SBStateDigest& sbs) const {}

    // Weld uses base class implementation of calcCompositeBodyInertiasInward() since
    // that is independent of mobilities.

    void realizeArticulatedBodyInertiasInward
       (const SBInstanceCache&          ic,
        const SBTreePositionCache&      pc,
        SBArticulatedBodyInertiaCache&  abc) const
    {
        ArticulatedInertia& P = updP(abc);
        P = ArticulatedInertia(getMk_G(pc));
        for (unsigned i=0 ; i<children.size() ; i++) {
            const PhiMatrix&          phiChild   = children[i]->getPhi(pc);
            const ArticulatedInertia& PPlusChild = children[i]->getPPlus(abc);

            P += PPlusChild.shift(phiChild.l());
        }
        updPPlus(abc) = P;
    }


    void realizeYOutward(
        const SBInstanceCache&,
        const SBTreePositionCache&              pc,
        const SBArticulatedBodyInertiaCache&    abc,
        SBDynamicsCache&                        dc) const
    {
        // This psi actually has the wrong sign, but it doesn't matter since we multiply
        // by it twice.

        SpatialMat psi = getPhi(pc).toSpatialMat();
        updY(dc) = ~psi * parent->getY(dc) * psi;
    }


    void calcUDotPass1Inward(
        const SBInstanceCache&      ic,
        const SBTreePositionCache&  pc,
        const SBArticulatedBodyInertiaCache&,
        const SBDynamicsCache&      dc,
        const Real*                 jointForces,
        const SpatialVec*           bodyForces,
        const Real*                 allUDot,
        SpatialVec*                 allZ,
        SpatialVec*                 allZPlus,
        Real*                       allEpsilon) const
    {
        const SpatialVec& myBodyForce  = bodyForces[nodeNum];
        SpatialVec&       z            = allZ[nodeNum];
        SpatialVec&       zPlus        = allZPlus[nodeNum];

        z = getMobilizerCentrifugalForces(dc) - myBodyForce;

        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild   = children[i]->getPhi(pc);
            const SpatialVec& zPlusChild = allZPlus[children[i]->getNodeNum()];

            z += phiChild * zPlusChild;                 // 18 flops
        }

        zPlus = z;
    }

    void calcUDotPass2Outward(
        const SBInstanceCache&,
        const SBTreePositionCache&  pc,
        const SBArticulatedBodyInertiaCache&,
        const SBTreeVelocityCache&  vc,
        const SBDynamicsCache&      dc,
        const Real*                 allEpsilon,
        SpatialVec*                 allA_GB,
        Real*                       allUDot,
        Real*                       allTau) const
    {
        SpatialVec& A_GB = allA_GB[nodeNum];

        const PhiMatrix&    phi = getPhi(pc);
        const SpatialVec&   a   = getMobilizerCoriolisAcceleration(vc);

        // Shift parent's acceleration outward (Ground==0). 12 flops
        const SpatialVec& A_GP  = allA_GB[parent->getNodeNum()];
        const SpatialVec  APlus = ~phi * A_GP;

        A_GB = APlus + a;  // no udot for weld
    }

    // A weld doesn't have udots but we still have to calculate z, zPlus,
    // for use by the parent of this body.
    void multiplyByMInvPass1Inward(
        const SBInstanceCache&      ic,
        const SBTreePositionCache&  pc,
        const SBArticulatedBodyInertiaCache&,
        const Real*                 f,
        SpatialVec*                 allZ,
        SpatialVec*                 allZPlus,
        Real*                       allEpsilon) const override
    {
        SpatialVec& z       = allZ[nodeNum];
        SpatialVec& zPlus   = allZPlus[nodeNum];

        z = 0;

        for (unsigned i=0; i<children.size(); i++) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zPlusChild = allZPlus[children[i]->getNodeNum()];
            z += phiChild * zPlusChild; // 18 flops
        }

        zPlus = z;
    }

    // Must set A_GB properly for propagation to children.
    void multiplyByMInvPass2Outward(
        const SBInstanceCache&,
        const SBTreePositionCache&  pc,
        const SBArticulatedBodyInertiaCache&,
        const Real*                 allEpsilon,
        SpatialVec*                 allA_GB,
        Real*                       allUDot) const override
    {
        SpatialVec&      A_GB = allA_GB[nodeNum];
        const PhiMatrix& phi  = getPhi(pc);

        // Shift parent's acceleration outward (Ground==0). 12 flops
        const SpatialVec& A_GP  = allA_GB[parent->getNodeNum()];
        const SpatialVec  APlus = ~phi * A_GP;

        A_GB = APlus;
    }

    // Also serves as pass 1 for inverse dynamics.
    void calcBodyAccelerationsFromUdotOutward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const
    {
        SpatialVec& A_GB = allA_GB[nodeNum];

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP + getMobilizerCoriolisAcceleration(vc); // no udot for weld
    }

    void calcInverseDynamicsPass2Inward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const SpatialVec*           allA_GB,
        const Real*                 jointForces,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allF,
        Real*                       allTau) const
    {
        const SpatialVec& myBodyForce   = bodyForces[nodeNum];
        const SpatialVec& A_GB          = allA_GB[nodeNum];
        SpatialVec&       F             = allF[nodeNum];

        // Start with rigid body force from desired body acceleration and
        // gyroscopic forces due to angular velocity, minus external forces
        // applied directly to this body.
        F = getMk_G(pc)*A_GB + getGyroscopicForce(vc) - myBodyForce;

        // Add in forces on children, shifted to this body.
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
            F += phiChild * FChild;
        }

        // no taus.
    }

    void multiplyByMPass1Outward(
        const SBTreePositionCache&  pc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const
    {
        SpatialVec& A_GB = allA_GB[nodeNum];

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP;
    }

    void multiplyByMPass2Inward(
        const SBTreePositionCache&  pc,
        const SpatialVec*           allA_GB,
        SpatialVec*                 allF,   // temp
        Real*                       allTau) const
    {
        const SpatialVec& A_GB  = allA_GB[nodeNum];
        SpatialVec&       F     = allF[nodeNum];

        F = getMk_G(pc)*A_GB;

        for (int i=0 ; i<(int)children.size() ; i++) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
            F += phiChild * FChild;
        }
    }

    void multiplyBySystemJacobian(
        const SBTreePositionCache&  pc,
        const Real*                 v,
        SpatialVec*                 Jv) const
    {
        SpatialVec& out = Jv[nodeNum];

        // Shift parent's result outward (ground result is 0).
        const SpatialVec outP = ~getPhi(pc) * Jv[parent->getNodeNum()];

        out = outP;
    }

    void multiplyBySystemJacobianTranspose(
        const SBTreePositionCache&  pc,
        SpatialVec*                 zTmp,
        const SpatialVec*           X,
        Real*                       JtX) const
    {
        const SpatialVec& in  = X[getNodeNum()];
        SpatialVec&       z   = zTmp[getNodeNum()];

        z = in;

        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild = children[i]->getPhi(pc);
            z += phiChild * zChild;
        }
        // No generalized speeds so no contribution to JtX.
    }

    void calcEquivalentJointForces(
        const SBTreePositionCache&  pc,
        const SBDynamicsCache&      dc,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allZ,
        Real*                       jointForces) const
    {
        const SpatialVec& myBodyForce  = bodyForces[nodeNum];
        SpatialVec&       z            = allZ[nodeNum];

        // Centrifugal forces are PA+b where P is articulated body inertia,
        // A is total coriolis acceleration, and b is gyroscopic force.
        z = myBodyForce - getTotalCentrifugalForces(dc);

        for (int i=0 ; i<(int)children.size() ; i++) {
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            z += phiChild * zChild;
        }
    }
};


// The Ground node is special because it doesn't need a mobilizer.
/*static*/ RigidBodyNode*
RigidBodyNode::createGroundNode() {
    return new RBGroundBody();
}

RigidBodyNode* MobilizedBody::WeldImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return new RBNodeWeld(
        getDefaultRigidBodyMassProperties(),
        getDefaultInboardFrame(),getDefaultOutboardFrame(),
        nextUSlot, nextUSqSlot, nextQSlot);
}

RigidBodyNode* MobilizedBody::GroundImpl::createRigidBodyNode(
    UIndex&        nextUSlot,
    USquaredIndex& nextUSqSlot,
    QIndex&        nextQSlot) const
{
    return RigidBodyNode::createGroundNode();
}

