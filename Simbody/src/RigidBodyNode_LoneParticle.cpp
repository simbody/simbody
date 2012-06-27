/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "MobilizedBodyImpl.h"
#include "RigidBodyNodeSpec_Translation.h"

/**
 * This is a specialized class used for MobilizedBody::Translation objects that satisfy
 * all of the following requirements:
 * 
 * <ul>
 * <li>The body has no children.</li>
 * <li>The body's parent is ground.</li>
 * <li>The inboard and outboard transforms are both identities.</li>
 * <li>The body is not reversed.</li>
 * </ul>
 * 
 * These assumptions allow lots of routines to be implemented in simpler, faster ways.
 */
class RBNodeLoneParticle : public RigidBodyNode {
public:
RBNodeLoneParticle(const MassProperties& mProps_B,
                UIndex&               nextUSlot,
                USquaredIndex&        nextUSqSlot,
                QIndex&               nextQSlot)
  : RigidBodyNode(mProps_B, Vec3(0), Vec3(0), QDotIsAlwaysTheSameAsU, QuaternionIsNeverUsed, false) {
    uIndex = nextUSlot;
    uSqIndex = nextUSqSlot;
    qIndex = nextQSlot;
    nextUSlot += 3;
    nextUSqSlot += 9;
    nextQSlot += 3;
}

const char* type() const {return "loneparticle";}
int  getDOF() const {return 3;}
int  getMaxNQ() const {return 3;}

int getNQInUse(const SBModelVars&) const {return 3;}
int getNUInUse(const SBModelVars&) const {return 3;}

bool isUsingQuaternion(const SBStateDigest&, MobilizerQIndex& startOfQuaternion) const {
    return false;
}

bool isUsingAngles(const SBStateDigest&, MobilizerQIndex& startOfAngles, int& nAngles) const {
    return false;
}

int calcQPoolSize(const SBModelVars&) const {
    return 0;
}

void performQPrecalculations(const SBStateDigest& sbs,
                             const Real* q, int nq,
                             Real* qCache,  int nQCache,
                             Real* qErr,    int nQErr) const {
    assert(q && nq==3 && nQCache==0 && nQErr==0);
}

void calcX_FM(const SBStateDigest& sbs,
              const Real* q,      int nq,
              const Real* qCache, int nQCache,
              Transform&  X_FM) const {
    assert(q && nq==3 && nQCache==0);
    X_FM = Transform(Rotation(), Vec3::getAs(&q[0]));
}

void calcQDot(const SBStateDigest&, const Real* u, Real* qdot) const {
    Vec3::updAs(qdot) = Vec3::getAs(u);
}
void calcQDotDot(const SBStateDigest&, const Real* udot, Real* qdotdot) const {
    Vec3::updAs(qdotdot) = Vec3::getAs(udot);
}

void multiplyByN(const SBStateDigest&, bool matrixOnRight, const Real* in, Real* out) const {
    Vec3::updAs(out) = Vec3::getAs(in);
}
void multiplyByNInv(const SBStateDigest&, bool matrixOnRight, const Real* in, Real* out) const {
    Vec3::updAs(out) = Vec3::getAs(in);
}
void multiplyByNDot(const SBStateDigest&, bool matrixOnRight, const Real* in, Real* out) const {
    Vec3::updAs(out) = 0;
}

bool enforceQuaternionConstraints(
    const SBStateDigest& sbs,
    Vector&            q,
    Vector&            qErrest) const {
    return false;
}

void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {
    Vec3::updAs(&outputQ[qIndex]) = Vec3::getAs(&inputQ[qIndex]);
}

void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {
    Vec3::updAs(&outputQ[qIndex]) = Vec3::getAs(&inputQ[qIndex]);
}

void setMobilizerDefaultModelValues
   (const SBTopologyCache&, SBModelVars&)        const {}

void setMobilizerDefaultInstanceValues    
   (const SBModelVars&,     SBInstanceVars&)     const {}
void setMobilizerDefaultTimeValues        
   (const SBModelVars&,     SBTimeVars&)         const {}
void setMobilizerDefaultPositionValues    
   (const SBModelVars&,     Vector& q)           const {q = 0;}
void setMobilizerDefaultVelocityValues    
   (const SBModelVars&,     Vector& u)           const {u = 0;}
void setMobilizerDefaultDynamicsValues    
   (const SBModelVars&,     SBDynamicsVars&)     const {}
void setMobilizerDefaultAccelerationValues
   (const SBModelVars&,     SBAccelerationVars&) const {}

void realizeModel(SBStateDigest& sbs) const {
}

void realizeInstance(const SBStateDigest& sbs) const {
    // Initialize cache entries that will never be changed at later stages.
    
    SBTreePositionCache& pc = sbs.updTreePositionCache();
    SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();
    SBDynamicsCache& dc = sbs.updDynamicsCache();
    SBTreeAccelerationCache& ac = sbs.updTreeAccelerationCache();
    Transform& X_FM = toB(pc.bodyJointInParentJointFrame);
    X_FM.updR().setRotationToIdentityMatrix();
    updV_FM(vc)[0] = Vec3(0);
    updV_PB_G(vc)[0] = Vec3(0);
    updVD_PB_G(vc)[0] = Vec3(0);
    updV_GB(vc)[0] = Vec3(0);
    updGyroscopicForce(vc) = SpatialVec(Vec3(0), Vec3(0));
    updMobilizerCoriolisAcceleration(vc) = SpatialVec(Vec3(0), Vec3(0));
    updTotalCoriolisAcceleration(vc) = SpatialVec(Vec3(0), Vec3(0));
    updMobilizerCentrifugalForces(dc) = SpatialVec(Vec3(0), Vec3(0));
    updTotalCentrifugalForces(dc) = SpatialVec(Vec3(0), Vec3(0));
    updY(dc) = SpatialMat(Mat33(0), Mat33(0), Mat33(0), Mat33(1.0/getMass()));
    updA_GB(ac)[0] = Vec3(0);
}

void realizePosition(const SBStateDigest& sbs) const {
    SBTreePositionCache& pc = sbs.updTreePositionCache();
    Transform& X_FM = toB(pc.bodyJointInParentJointFrame);
    const Vec3& q = Vec3::getAs(&sbs.getQ()[qIndex]);
    X_FM.updP() = q;
    updX_PB(pc) = X_FM;
    updX_GB(pc) = X_FM;
    updPhi(pc) = PhiMatrix(q);
    updCOM_G(pc) = q + getCOM_B(); // 3 flops
    updMk_G(pc) = SpatialInertia(getMass(), getCOM_B(), getUnitInertia_OB_B());
}

void realizeVelocity(const SBStateDigest& sbs) const {
    SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();
    const Vector& allU = sbs.getU();
    const Vec3& u = Vec3::getAs(&allU[uIndex]);

    Vec3::updAs(&sbs.updQDot()[qIndex]) = Vec3::getAs(&allU[uIndex]);
    updV_FM(vc)[1] = u;
    updV_PB_G(vc)[1] = u;
    updV_GB(vc)[1] = u;
}

void realizeDynamics(const SBArticulatedBodyInertiaCache& abc, const SBStateDigest& sbs) const {
}

// There is no realizeAcceleration().

void realizeReport(const SBStateDigest& sbs) const {
}

void realizeArticulatedBodyInertiasInward(
        const SBInstanceCache&          ic,
        const SBTreePositionCache&      pc,
        SBArticulatedBodyInertiaCache&  abc) const {
    ArticulatedInertia& P     = updP(abc);
    ArticulatedInertia& PPlus = updPPlus(abc);

    PPlus = P = ArticulatedInertia(getMk_G(pc));
}

void realizeYOutward(
            const SBInstanceCache&                ic,
            const SBTreePositionCache&            pc,
            const SBArticulatedBodyInertiaCache&  abc,
            SBDynamicsCache&                      dc) const {
}

void calcCompositeBodyInertiasInward
   (const SBTreePositionCache& pc, 
    Array_<SpatialInertia,MobilizedBodyIndex>& R) const {
    toB(R) = getMk_G(pc);
}

void multiplyBySystemJacobian(
        const SBTreePositionCache&  pc,
        const Real*                 v,
        SpatialVec*                 Jv) const {
    const Vec3& in = Vec3::getAs(&v[uIndex]);
    SpatialVec& out = Jv[nodeNum];
    out = ~getPhi(pc) * Jv[parent->getNodeNum()];
    out[1] += in;
}

void multiplyBySystemJacobianTranspose(
        const SBTreePositionCache&  pc, 
        SpatialVec*                 zTmp,
        const SpatialVec*           X, 
        Real*                       JtX) const {
    const SpatialVec& in = X[getNodeNum()];
    Vec3& out = Vec3::updAs(&JtX[getUIndex()]);
    SpatialVec& z = zTmp[getNodeNum()];
    z = in;
    out = z[1]; 
}

void calcEquivalentJointForces(
        const SBTreePositionCache&  pc,
        const SBDynamicsCache&      dc,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allZ,
        Real*                       jointForces) const {
    const SpatialVec& myBodyForce = bodyForces[nodeNum];
    SpatialVec& z = allZ[nodeNum];
    Vec3& eps = Vec3::updAs(&jointForces[uIndex]);
    z = myBodyForce;
    eps = z[1];
}

void calcUDotPass1Inward(
        const SBInstanceCache&                  ic,
        const SBTreePositionCache&              pc,
        const SBArticulatedBodyInertiaCache&    abc,
        const SBDynamicsCache&                  dc,
        const Real*                             jointForces,
        const SpatialVec*                       bodyForces,
        const Real*                             allUDot,
        SpatialVec*                             allZ,
        SpatialVec*                             allZPlus,
        Real*                                   allEpsilon) const 
{
    const Vec3&         f       = Vec3::getAs(&jointForces[uIndex]);
    const SpatialVec&   F       = bodyForces[nodeNum];
    SpatialVec&         z       = allZ[nodeNum];
    SpatialVec&         zPlus   = allZPlus[nodeNum];
    Vec3&               eps     = Vec3::updAs(&allEpsilon[uIndex]);

    const bool isPrescribed = isUDotKnown(ic);

    z = -F;

    if (isPrescribed && !isUDotKnownToBeZero(ic)) {
        const Vec3& udot = Vec3::getAs(&allUDot[uIndex]);
        z[1] += getMass()*udot; // == P*H*udot
    }

    // Lone particle has no children.

    eps = f - z[1];

    zPlus = z;
    if (!isPrescribed)
        zPlus[1] += eps;
}

void calcUDotPass2Outward(
        const SBInstanceCache&                  ic,
        const SBTreePositionCache&              pc,
        const SBArticulatedBodyInertiaCache&    abc,
        const SBTreeVelocityCache&              vc,
        const SBDynamicsCache&                  dc,
        const Real*                             allEpsilon,
        SpatialVec*                             allA_GB,
        Real*                                   allUDot,
        Real*                                   allTau) const {
    const Vec3& eps = Vec3::getAs(&allEpsilon[uIndex]);
    SpatialVec& A_GB = allA_GB[nodeNum];
    Vec3& udot = Vec3::updAs(&allUDot[uIndex]);

    const bool isPrescribed = isUDotKnown(ic);

    if (isPrescribed) {
        const PresForcePoolIndex tauIx = 
            ic.getMobodInstanceInfo(nodeNum).firstPresForce;
        assert(tauIx.isValid());
        Vec3& tau = Vec3::updAs(&allTau[tauIx]);
        tau = eps; // our sign convention
    } else 
        udot = eps/getMass();

    A_GB = SpatialVec(Vec3(0), udot);
}

// Note that we're not setting z temporaries here; you can't count on that as
// a side effect the M^-1*f kernel.
void multiplyByMInvPass1Inward(
        const SBInstanceCache&                  ic,
        const SBTreePositionCache&              pc,
        const SBArticulatedBodyInertiaCache&    abc,
        const SBDynamicsCache&                  dc,
        const Real*                             jointForces,
        SpatialVec*                             allZ,
        SpatialVec*                             allZPlus,
        Real*                                   allEpsilon) const 
{
    if (isUDotKnown(ic)) // prescribed
        return;

    // We promised not to look at f if it is part of f_p (prescribed).
    const Vec3& f     = Vec3::getAs(&jointForces[uIndex]);
    Vec3&       eps   = Vec3::updAs(&allEpsilon[uIndex]);

    eps = f;
}

// Outward pass must set A_GB properly so it can be propagated
// to children.
void multiplyByMInvPass2Outward(
        const SBInstanceCache&                  ic,
        const SBTreePositionCache&              pc,
        const SBArticulatedBodyInertiaCache&    abc,
        const SBDynamicsCache&                  dc,
        const Real*                             allEpsilon,
        SpatialVec*                             allA_GB,
        Real*                                   allUDot) const 
{
    const bool isPrescribed = isUDotKnown(ic);
    const Vec3& eps = Vec3::getAs(&allEpsilon[uIndex]);
    SpatialVec& A_GB = allA_GB[nodeNum];
    Vec3&       udot = Vec3::updAs(&allUDot[uIndex]);

    if (isPrescribed) {
        udot = 0;
        A_GB = SpatialVec(Vec3(0), Vec3(0));
    } else {
        udot = eps/getMass();
        A_GB = SpatialVec(Vec3(0), udot);
    }
}

// Also serves as pass 1 for inverse dynamics.
void calcBodyAccelerationsFromUdotOutward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const {
    const Vec3& udot = Vec3::getAs(&allUDot[uIndex]);
    SpatialVec& A_GB = allA_GB[nodeNum];
    A_GB = SpatialVec(Vec3(0), udot);
}
void calcInverseDynamicsPass2Inward(
        const SBTreePositionCache&  pc,
        const SBTreeVelocityCache&  vc,
        const SpatialVec*           allA_GB,
        const Real*                 jointForces,
        const SpatialVec*           bodyForces,
        SpatialVec*                 allF,
        Real*                       allTau) const {
    const Vec3& myJointForce = Vec3::getAs(&jointForces[uIndex]);
    const SpatialVec& myBodyForce = bodyForces[nodeNum];
    const SpatialVec& A_GB = allA_GB[nodeNum];
    SpatialVec& F = allF[nodeNum];
    Vec3& tau = Vec3::updAs(&allTau[uIndex]);
    F = getMk_G(pc)*A_GB - myBodyForce;
    tau = F[1] - myJointForce;
}

void multiplyByMPass1Outward(
        const SBTreePositionCache&  pc,
        const Real*                 allUDot,
        SpatialVec*                 allA_GB) const {
    const Vec3& udot = Vec3::getAs(&allUDot[uIndex]);
    SpatialVec& A_GB = allA_GB[nodeNum];
    A_GB = SpatialVec(Vec3(0), udot);
}
void multiplyByMPass2Inward(
        const SBTreePositionCache&  pc,
        const SpatialVec*           allA_GB,
        SpatialVec*                 allF,
        Real*                       allTau) const {
    const SpatialVec& A_GB = allA_GB[nodeNum];
    SpatialVec& F = allF[nodeNum];
    Vec3& tau = Vec3::updAs(&allTau[uIndex]);
    F = getMk_G(pc)*A_GB; 
    tau = F[1];
}

const SpatialVec& getHCol(const SBTreePositionCache& pc, int j) const {
    Mat<2,3,Vec3> H = Mat<2,3,Vec3>::getAs(&pc.storageForH[2*uIndex]);
    SpatialVec& col = H(j);
    col = SpatialVec(Vec3(0), Vec3(0));
    col[1][j] = 1;
    return col;
}

const SpatialVec& getH_FMCol(const SBTreePositionCache& pc, int j) const {
    Mat<2,3,Vec3> H = Mat<2,3,Vec3>::getAs(&pc.storageForH_FM[2*uIndex]);
    SpatialVec& col = H(j);
    col = SpatialVec(Vec3(0), Vec3(0));
    col[1][j] = 1;
    return col;
}

void setQToFitTransformImpl(const SBStateDigest&, const Transform& X_F0M0, Vector& q) const {
    Vec3::updAs(&q[qIndex]) = X_F0M0.p();
}
void setQToFitRotationImpl(const SBStateDigest&, const Rotation& R_F0M0, Vector& q) const {
}
void setQToFitTranslationImpl(const SBStateDigest&, const Vec3& p_F0M0, Vector& q) const {
    Vec3::updAs(&q[qIndex]) = p_F0M0;
}

void setUToFitVelocityImpl(const SBStateDigest&, const Vector& q, const SpatialVec& V_F0M0, Vector& u) const {
    Vec3::updAs(&u[uIndex]) = V_F0M0[1];
}
void setUToFitAngularVelocityImpl(const SBStateDigest&, const Vector& q, const Vec3& w_F0M0, Vector& u) const {    
}
void setUToFitLinearVelocityImpl(const SBStateDigest&, const Vector& q, const Vec3& v_F0M0, Vector& u) const {
    Vec3::updAs(&u[uIndex]) = v_F0M0;
}

};

RigidBodyNode* MobilizedBody::TranslationImpl::createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const {
    if (!hasChildren && getMyParentMobilizedBodyIndex() == 0 && !isReversed() &&
            getDefaultInboardFrame().p() == 0 && getDefaultInboardFrame().R() == Mat33(1) &&
            getDefaultOutboardFrame().p() == 0 && getDefaultOutboardFrame().R() == Mat33(1)) {
        // This satisfies all the requirements to use RBNodeLoneParticle.
        
        return new RBNodeLoneParticle(getDefaultRigidBodyMassProperties(), nextUSlot,nextUSqSlot,nextQSlot);
    }
    
    // Use RBNodeTranslate for the general case.
    
    bool noX_MB = (getDefaultOutboardFrame().p() == 0 && getDefaultOutboardFrame().R() == Mat33(1));
    bool noR_PF = (getDefaultInboardFrame().R() == Mat33(1));
    if (noX_MB) {
        if (noR_PF)
            return new RBNodeTranslate<true, true> (
                getDefaultRigidBodyMassProperties(),
                getDefaultInboardFrame(),getDefaultOutboardFrame(),
                isReversed(),
                nextUSlot,nextUSqSlot,nextQSlot);
        else
            return new RBNodeTranslate<true, false> (
                getDefaultRigidBodyMassProperties(),
                getDefaultInboardFrame(),getDefaultOutboardFrame(),
                isReversed(),
                nextUSlot,nextUSqSlot,nextQSlot);
    }
    else {
        if (noR_PF)
            return new RBNodeTranslate<false, true> (
                getDefaultRigidBodyMassProperties(),
                getDefaultInboardFrame(),getDefaultOutboardFrame(),
                isReversed(),
                nextUSlot,nextUSqSlot,nextQSlot);
        else
            return new RBNodeTranslate<false, false> (
                getDefaultRigidBodyMassProperties(),
                getDefaultInboardFrame(),getDefaultOutboardFrame(),
                isReversed(),
                nextUSlot,nextUSqSlot,nextQSlot);
    }
}
