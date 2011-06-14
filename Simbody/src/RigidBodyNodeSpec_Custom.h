#ifndef SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CUSTOM_H_
#define SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CUSTOM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-9 Stanford University and the Authors.         *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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
 * Define the RigidBodyNode that implements Custom mobilizers.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"
#include "MobilizedBodyImpl.h" // need Custom::ImplementationImpl


/**
 * RigidBodyNodeSpec for Custom Mobilizers. This is still templatized
 * by the number of u's (mobilities) in the user-defined Mobilizer.
 */
template <int nu, bool noX_MB, bool noR_PF>
class RBNodeCustom : public RigidBodyNodeSpec<nu, false, noX_MB, noR_PF> {
    typedef typename RigidBodyNodeSpec<nu, false, noX_MB, noR_PF>::HType HType;
public:
    RBNodeCustom(const MobilizedBody::Custom::Implementation& impl,
                 const MassProperties&  mProps_B, 
                 const Transform&       X_PF, 
                 const Transform&       X_BM,
                 bool                   isReversed,
                 UIndex&                nextUSlot, 
                 USquaredIndex&         nextUSqSlot, 
                 QIndex&                nextQSlot)
    :   RigidBodyNodeSpec<nu, false, noX_MB, noR_PF>(mProps_B, X_PF, X_BM, nextUSlot, nextUSqSlot, nextQSlot, 
                              RigidBodyNode::QDotMayDifferFromU,
                              impl.getImpl().getNumAngles() == 4 ? RigidBodyNode::QuaternionMayBeUsed 
                                                                 : RigidBodyNode::QuaternionIsNeverUsed,
                              isReversed),
        impl(impl), nq(impl.getImpl().getNQ()), nAngles(impl.getImpl().getNumAngles()) 
    {
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
    void calcQDot(const SBStateDigest& sbs, 
                                 const Real* u, Real* qdot) const 
    {
        const int nqInUse = getNQInUse(sbs.getModelVars());
        impl.multiplyByN(sbs.getState(), false, nu, u, getNQInUse(sbs.getModelVars()), qdot);
        for (int i = nqInUse; i < nq; ++i)
            qdot[i] = 0.0;
    }
    void calcQDotDot(const SBStateDigest& sbs, 
                                       const Real* udot, Real* qdotdot) const 
    {
        const SBModelVars&          mv = sbs.getModelVars();
        const SBTreePositionCache&  pc = sbs.getTreePositionCache();
        const int nqInUse = getNQInUse(sbs.getModelVars());
        const Real* u = &sbs.getU()[this->getUIndex()];
        impl.multiplyByN(sbs.getState(), false, nu, udot, nqInUse, qdotdot);
        Real temp[7];
        impl.multiplyByNDot(sbs.getState(), false, nu, u, nqInUse, temp);
        for (int i = 0; i < nqInUse; ++i)
            qdotdot[i] += temp[i];
        for (int i = nqInUse; i < nq; ++i)
            qdotdot[i] = 0.0;
    }
    void multiplyByN(const SBStateDigest& sbs, bool matrixOnRight, 
                     const Real* in, Real* out) const 
    {
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
        impl.multiplyByN(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
    }
    void multiplyByNInv(const SBStateDigest& sbs, bool matrixOnRight,
                        const Real* in, Real* out) const 
    {
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
        impl.multiplyByNInv(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
    }
    void multiplyByNDot(const SBStateDigest& sbs, bool matrixOnRight, 
                        const Real* in, Real* out) const 
    {
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
        impl.multiplyByNDot(sbs.getState(), matrixOnRight, nIn, in, nOut, out);
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

    void setQToFitTransformImpl(const SBStateDigest& sbs, const Transform& X_FM, Vector& q) const {
        impl.setQToFitTransform(sbs.getState(), X_FM, this->getNQInUse(sbs.getModelVars()), &q[this->getQIndex()]);
    }
    void setQToFitRotationImpl(const SBStateDigest& sbs, const Rotation& R_FM, Vector& q) const {
        setQToFitTransformImpl(sbs, Transform(R_FM), q);
    }
    void setQToFitTranslationImpl(const SBStateDigest& sbs, const Vec3& p_FM, Vector& q) const {
        setQToFitTransformImpl(sbs, Transform(p_FM), q);
    }

    void setUToFitVelocityImpl(const SBStateDigest& sbs, const Vector& q, const SpatialVec& V_FM, Vector& u) const {
        impl.setUToFitVelocity(sbs.getState(), V_FM, nu, &u[this->getUIndex()]);
    }
    void setUToFitAngularVelocityImpl(const SBStateDigest& sbs, const Vector& q, const Vec3& w_FM, Vector& u) const {
        setUToFitVelocityImpl(sbs, q, SpatialVec(w_FM, Vec3(0)), u);
    }
    void setUToFitLinearVelocityImpl(const SBStateDigest& sbs, const Vector& q, const Vec3& v_FM, Vector& u) const {
        setUToFitVelocityImpl(sbs, q, SpatialVec(Vec3(0), v_FM), u);
    }

        // VIRTUAL METHODS FOR SINGLE-NODE OPERATOR CONTRIBUTIONS //

    // We're not going to attempt to cache any q precalculations.
    int calcQPoolSize(const SBModelVars&) const {return 0;}

    void performQPrecalculations(const SBStateDigest& sbs,
                                 const Real* q, int nq,
                                 Real* qCache,  int nQCache,
                                 Real* qErr,    int nQErr) const
    {
        assert(nq==getNQInUse(sbs.getModelVars()) && nQCache==0); 
        if (nAngles == 4 && !this->getUseEulerAngles(sbs.getModelVars())) {
            // Need to calculate qerr
            assert(nQErr==1);
            const Real quatLen = Vec4::getAs(&q[0]).norm();
            qErr[0] = quatLen - Real(1);    // normalization error
        }
    }
    
    void calcX_FM(const SBStateDigest& sbs,
                  const Real* q,      int nq,
                  const Real* qCache, int nQCache,
                  Transform&  X_F0M0) const
    {
        assert(nq==getNQInUse(sbs.getModelVars()) && nQCache==0); 
        // Note: quaternion will be unnormalized.
        X_F0M0 = impl.calcMobilizerTransformFromQ(sbs.getState(), nq, q);
    }
    
    void calcAcrossJointVelocityJacobian(
        const SBStateDigest& sbs,
        HType&  H_F0M0) const {
        for (int i = 0; i < nu; ++i) {
            Vec<nu> u(0);
            u[i] = 1;
            H_F0M0(i) = impl.multiplyByHMatrix(sbs.getState(), nu, &u[0]);
        }
    }

    void calcAcrossJointVelocityJacobianDot(
        const SBStateDigest& sbs,
        HType&  HDot_F0M0) const {
        for (int i = 0; i < nu; ++i) {
            Vec<nu> u(0);
            u[i] = 1;
            HDot_F0M0(i) = impl.multiplyByHDotMatrix(sbs.getState(), nu, &u[0]);
        }
    }

private:
    const MobilizedBody::Custom::Implementation& impl;
    const int nq, nAngles;
};


#endif // SimTK_SIMBODY_RIGID_BODY_NODE_SPEC_CUSTOM_H_

