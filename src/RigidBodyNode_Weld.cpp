/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *    Charles Schwieters (NIH): wrote the public domain IVM code from which   *
 *                              this was derived.                             *
 *    Peter Eastman wrote the Weld joint.                                     *
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
    bool isUsingAngles(const SBStateDigest& sbs, MobilizerQIndex& ix, int& nAngles) const
    {   ix.invalidate(); nAngles = 0; return false; }
    void copyQ(const SBStateDigest&, const Vector&, Vector&) const {}
    void copyU(const SBStateDigest&, const Vector&, Vector&) const {}

    void calcJointSinCosQNorm(
        const SBModelVars&  mv, 
        const SBModelCache& mc,
        const SBInstanceCache& ic,
        const Vector&       q, 
        Vector&             sine, 
        Vector&             cosine, 
        Vector&             qErr,
        Vector&             qnorm) const {}


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

    
    void multiplyByN(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                  bool matrixOnRight, 
                                  const Real* in, Real* out) const {}
    void multiplyByNInv(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q,
                                     bool matrixOnRight,
                                     const Real* in, Real* out) const {}
    void multiplyByNDot(const SBStateDigest&, bool useEulerAnglesIfPossible, const Real* q, const Real* u,
                                     bool matrixOnRight,
                                     const Real* in, Real* out) const {}

    void calcQDot(const SBStateDigest&,const Vector&,Vector&) const {}
    void calcQDotDot(const SBStateDigest&, const Vector&, Vector&) const {}

    bool enforceQuaternionConstraints(
        const SBStateDigest& sbs,
        Vector&             q,
        Vector&             qErrest) const {return false;}

    void convertToEulerAngles(const Vector& inputQ, Vector& outputQ) const {}
    void convertToQuaternions(const Vector& inputQ, Vector& outputQ) const {}

    void setVelFromSVel(
        const SBPositionCache& pc, 
        const SBVelocityCache& vc,
        const SpatialVec&      sVel, 
        Vector&                u) const {}
};

/**
 * This is the distinguished body representing the immobile ground frame. Other bodies may
 * be fixed to this one, but only this is the actual Ground.
 */
class RBGroundBody : public ImmobileRigidBodyNode {
public:
    RBGroundBody()
    :   ImmobileRigidBodyNode(MassProperties(Infinity, Vec3(0), Infinity*Inertia(1)), 
                              Transform(), Transform(),
                              UIndex(0), USquaredIndex(0), QIndex(0)) {}

    const char* type() const { return "ground"; }

    // Imagine that Ground is welded to the universe at its origin. Its "inboard"
    // joint is a Weld which keeps the Ground M frame aligned forever with the
    // universe's F frame.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_F0M0) const 
    {   X_F0M0.setToZero(); }

    // Ground's motion is prescribed at zero.
    void setMobilizerDefaultModelValues(const SBTopologyCache&, SBModelVars& mv) const
    {   mv.prescribed[0] = true; }

    // TODO: should ground set the various cache entries here?
    void realizeModel(SBStateDigest& sbs) const {}
    void realizeInstance(SBStateDigest& sbs) const {}
    void realizeTime(SBStateDigest& sbs) const {}
    void realizePosition(SBStateDigest& sbs) const {}
    void realizeVelocity(SBStateDigest& sbs) const {}
    void realizeDynamics(SBStateDigest& sbs) const {}
    void realizeAcceleration(SBStateDigest& sbs) const {}
    void realizeReport(SBStateDigest& sbs) const {}

    // Ground's "composite" body inertia is still the infinite mass
    // and inertia it started with; no need to look at the children.
    // This overrides the base class default implementation.
    void calcCompositeBodyInertiasInward(
        const SBPositionCache& pc,
        Vector_<SpatialMat>&   R) const
    {   R[0] = SpatialMat(Mat33(Infinity)); }

    // Ground's "articulated" body inertia is still the infinite mass and
    // inertia it started with; no need to look at the children.
    void realizeArticulatedBodyInertiasInward(
        const SBPositionCache&,
        SBDynamicsCache&       dc) const 
    {   updP(dc) = SpatialMat(Mat33(Infinity)); }

    void realizeYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const 
    {   updY(dc) = SpatialMat(Mat33(0)); }

    void realizeZ(
        const SBStateDigest&       sbs,
        const Vector&              mobilityForces,
        const Vector_<SpatialVec>& bodyForces) const 
    {   
        const SBPositionCache& pc = sbs.getPositionCache();
        SBAccelerationCache&   ac = sbs.updAccelerationCache(); 
        updZ(ac) = -bodyForces[0];
        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild    = children[i]->getZ(ac);
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& GepsChild = children[i]->getGepsilon(ac);
            updZ(ac) += phiChild * (zChild + GepsChild);
        }
        updGepsilon(ac) = SpatialVec(Vec3(0));
    }
    void realizeAccel(
            const SBStateDigest&   sbs,
            Vector&                udot,
            Vector&                qdotdot) const 
    {   SBAccelerationCache& ac = sbs.updAccelerationCache(); 
        updA_GB(ac) = 0; }

    void calcUDotPass1Inward(
        const SBPositionCache&     pc,
        const SBDynamicsCache&,
        const Vector&              jointForces,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const
    {
        allZ[0] = -bodyForces[0];
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];
            allZ[0] += phiChild * (zChild + GepsChild);
        }
        allGepsilon[0] = 0;
    } 
    void calcUDotPass2Outward(
        const SBPositionCache&,
        const SBVelocityCache&,
        const SBDynamicsCache&,
        const Vector&              epsilonTmp,
        Vector_<SpatialVec>&       allA_GB,
        Vector&                    allUDot) const
    {
        allA_GB[0] = 0;
    }

    void calcMInverseFPass1Inward(
        const SBPositionCache&     pc,
        const SBDynamicsCache&,
        const Vector&              f,
        Vector_<SpatialVec>&       allZ,
        Vector_<SpatialVec>&       allGepsilon,
        Vector&                    allEpsilon) const
    {
        allZ[0] = 0;
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];
            allZ[0] += phiChild * (zChild + GepsChild);
        }
        allGepsilon[0] = 0;
    } 

    void calcMInverseFPass2Outward(
        const SBPositionCache&,
        const SBDynamicsCache&,
        const Vector&               epsilonTmp,
        Vector_<SpatialVec>&        allA_GB,
        Vector&                     allUDot) const
    {
        allA_GB[0] = 0;
    }

    void calcInverseDynamicsPass1Outward(
        const SBPositionCache& pc,
        const SBVelocityCache& vc,
        const Vector&          allUDot,
        Vector_<SpatialVec>&   allA_GB) const 
    {
        allA_GB[0] = 0;
    }

    // Here Ground is the last body processed. Although it has no mobility forces
    // we can still collect up all the forces from the base bodies to Ground
    // in case anyone cares.
    void calcInverseDynamicsPass2Inward(
        const SBPositionCache&      pc,
        const SBVelocityCache&      vc,
        const Vector_<SpatialVec>&  allA_GB,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allF,
        Vector&                     allTau) const
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

	void calcMVPass1Outward(
		const SBPositionCache& pc,
		const Vector&          allUDot,
		Vector_<SpatialVec>&   allA_GB) const
    {
        allA_GB[0] = 0;
    }

	void calcMVPass2Inward(
		const SBPositionCache& pc,
		const Vector_<SpatialVec>& allA_GB,
		Vector_<SpatialVec>&       allF,
		Vector&                    allTau) const
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


    void calcSpatialKinematicsFromInternal(
        const SBPositionCache&      pc,
        const Vector&               v,
        Vector_<SpatialVec>&        Jv) const    
    {
        Jv[0] = SpatialVec(Vec3(0));
    }

    void calcInternalGradientFromSpatial(
        const SBPositionCache&      pc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const 
    {
        zTmp[0] = X[0];
        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild = children[i]->getPhi(pc);
            zTmp[0] += phiChild * zChild;
        }
    }

    void calcEquivalentJointForces(
        const SBPositionCache& pc,
        const SBDynamicsCache&,
        const Vector_<SpatialVec>& bodyForces,
        Vector_<SpatialVec>&       allZ,
        Vector&                    jointForces) const 
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

    // A Weld holds the mobilized body's M frame coincident with the
    // parent body's F frame forever.
    void calcAcrossJointTransform(
        const SBStateDigest& sbs,
        const Vector&        q,
        Transform&           X_F0M0) const
    {   X_F0M0.setToZero(); }

    // If you want to think of the Weld mobilizer as being always prescribed,
    // that's fine.
    void setMobilizerDefaultModelValues(const SBTopologyCache&, SBModelVars& mv) const
    {   mv.prescribed[getNodeNum()] = true; }

    void realizeModel(SBStateDigest& sbs) const {}
    void realizeInstance(SBStateDigest& sbs) const {}
    void realizeTime(SBStateDigest& sbs) const {}

    void realizePosition(SBStateDigest& sbs) const {
        SBPositionCache& pc = sbs.updPositionCache();

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

        updInertia_OB_G(pc) = getInertia_OB_B().reexpress(~getX_GB(pc).R());
        updCB_G(pc)         = getX_GB(pc).R()*getCOM_B();
        updCOM_G(pc) = getX_GB(pc).p() + getCB_G(pc);

        // Calc Mk: the spatial inertia matrix about the body origin.
        // Note that this is symmetric; offDiag is *skew* symmetric so
        // that transpose(offDiag) = -offDiag.
        // Note: we need to calculate this now so that we'll be able to calculate
        // kinetic energy without going past the Velocity stage.
        
        const Mat33 offDiag = crossMat(getMass()*getCB_G(pc));
        updMk(pc) = SpatialMat( getInertia_OB_G(pc).toMat33() ,     offDiag ,
                                       -offDiag             ,   Mat33(getMass()) );
    }
    
    void realizeVelocity(SBStateDigest& sbs) const {
        const SBPositionCache& pc = sbs.getPositionCache();
        SBVelocityCache& vc = sbs.updVelocityCache();

        updV_FM(vc) = 0;
        updV_PB_G(vc) = 0;
        updVD_PB_G(vc) = 0;
        calcJointIndependentKinematicsVel(pc,vc);
    }

    void realizeDynamics(SBStateDigest& sbs) const {
        // Mobilizer-specific.
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        SBDynamicsCache& dc = sbs.updDynamicsCache();
        
        // Mobilizer independent.
        calcJointIndependentDynamicsVel(pc,vc,dc);
    }


    void realizeAcceleration(SBStateDigest& sbs) const {}
    void realizeReport(SBStateDigest& sbs) const {}

    // Weld uses base class implementation of calcCompositeBodyInertiasInward() since
    // that is independent of mobilities.

    void realizeArticulatedBodyInertiasInward(const SBPositionCache& pc, SBDynamicsCache& dc) const {
        updP(dc) = getMk(pc);
		for (unsigned i=0 ; i<children.size() ; i++) {
            const PhiMatrix&  phiChild    = children[i]->getPhi(pc);
            const SpatialMat& tauBarChild = children[i]->getTauBar(dc);
            const SpatialMat& PChild      = children[i]->getP(dc);
            const SpatialMat& psiChild    = children[i]->getPsi(dc);

			// TODO: too slow -- can get a 50% speedup by exploiting
            // symmetry; see the RigidBodyNodeSpec<dof> 
            // implementation for more info.
            // (Subtracting here because our Psi has reverse sign convention
            // from Jain's.)
            updP(dc) -= psiChild * (PChild * ~phiChild);
		}

        // Note our backwards sign convention for TauBar and Psi (from Jain's).
        updTauBar(dc) = -1; // -identity
        // TODO: wasting 33 flops negating; just re-create from -p.
        updPsi(dc) = -getPhi(pc).toSpatialMat();
    }


    void realizeYOutward(
        const SBPositionCache& pc,
        SBDynamicsCache&       dc) const
    {
        updY(dc) = ~getPsi(dc) * parent->getY(dc) * getPsi(dc);
    }

    void realizeZ(
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
        updGepsilon(ac) = SpatialVec(Vec3(0));
    }


    void realizeAccel(
            const SBStateDigest&   sbs,
            Vector&                udot,
            Vector&                qdotdot) const
    {
        const SBPositionCache& pc = sbs.getPositionCache();
        const SBVelocityCache& vc = sbs.getVelocityCache();
        SBAccelerationCache&   ac = sbs.updAccelerationCache();
        const SpatialVec alphap = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

        updA_GB(ac) = alphap + getCoriolisAcceleration(vc);  
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

        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];
            z += phiChild * (zChild + GepsChild);
        }

        Geps = 0;
    }

    void calcUDotPass2Outward(
        const SBPositionCache& pc,
        const SBVelocityCache& vc,
        const SBDynamicsCache& dc,
        const Vector&          allEpsilon,
        Vector_<SpatialVec>&   allA_GB,
        Vector&                allUDot) const
    {
        SpatialVec& A_GB = toB(allA_GB);

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP + getCoriolisAcceleration(vc);  
    }
    
    void calcMInverseFPass1Inward(
        const SBPositionCache& pc,
        const SBDynamicsCache& dc,
        const Vector&          f,
        Vector_<SpatialVec>&   allZ,
        Vector_<SpatialVec>&   allGepsilon,
        Vector&                allEpsilon) const 
    {
        SpatialVec& z       = toB(allZ);
        SpatialVec& Geps    = toB(allGepsilon);

        z = 0;
        for (unsigned i=0; i<children.size(); ++i) {
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
        SpatialVec& A_GB = toB(allA_GB);

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP;
    }

    void calcInverseDynamicsPass1Outward(
        const SBPositionCache& pc,
        const SBVelocityCache& vc,
        const Vector&          allUDot,
        Vector_<SpatialVec>&   allA_GB) const 
    {
        SpatialVec& A_GB = toB(allA_GB);

        // Shift parent's A_GB outward. (Ground A_GB is zero.)
        const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

        A_GB = A_GP + getCoriolisAcceleration(vc); 
    }

    void calcInverseDynamicsPass2Inward(
        const SBPositionCache&      pc,
        const SBVelocityCache&      vc,
        const Vector_<SpatialVec>&  allA_GB,
        const Vector&               jointForces,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allF,
        Vector&                     allTau) const
    {
        const SpatialVec& myBodyForce   = fromB(bodyForces);
        const SpatialVec& A_GB          = fromB(allA_GB);
        SpatialVec&       F		        = toB(allF);

        // Start with rigid body force from desired body acceleration and
        // gyroscopic forces due to angular velocity, minus external forces
        // applied directly to this body.
	    F = getMk(pc)*A_GB + getGyroscopicForce(vc) - myBodyForce;

        // Add in forces on children, shifted to this body.
        for (unsigned i=0; i<children.size(); ++i) {
            const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
            const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
            F += phiChild * FChild;
        }

        // no taus.
    }

	void calcMVPass1Outward(
		const SBPositionCache& pc,
		const Vector&          allUDot,
		Vector_<SpatialVec>&   allA_GB) const
	{
		SpatialVec& A_GB = toB(allA_GB);

		// Shift parent's A_GB outward. (Ground A_GB is zero.)
		const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

		A_GB = A_GP;  
	}

	void calcMVPass2Inward(
		const SBPositionCache& pc,
		const Vector_<SpatialVec>& allA_GB,
		Vector_<SpatialVec>&       allF,	// temp
		Vector&                    allTau) const 
	{
		const SpatialVec& A_GB  = fromB(allA_GB);
		SpatialVec&       F		= toB(allF);

		F = getMk(pc)*A_GB;

		for (int i=0 ; i<(int)children.size() ; i++) {
			const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
			const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
			F += phiChild * FChild;
		}
	}

    void calcSpatialKinematicsFromInternal(
        const SBPositionCache&      pc,
        const Vector&               v,
        Vector_<SpatialVec>&        Jv) const    
    {
        SpatialVec& out = toB(Jv);

        // Shift parent's result outward (ground result is 0).
        const SpatialVec outP = ~getPhi(pc) * parent->fromB(Jv);

        out = outP;  
    }

    void calcInternalGradientFromSpatial(
        const SBPositionCache&      pc, 
        Vector_<SpatialVec>&        zTmp,
        const Vector_<SpatialVec>&  X, 
        Vector&                     JX) const
    {
        const SpatialVec& in  = X[getNodeNum()];
        SpatialVec&       z   = zTmp[getNodeNum()];

        z = in;

        for (unsigned i=0; i<children.size(); ++i) {
            const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
            const PhiMatrix&  phiChild = children[i]->getPhi(pc);
            z += phiChild * zChild;
        }
    }

    void calcEquivalentJointForces(
        const SBPositionCache&      pc,
        const SBDynamicsCache&      dc,
        const Vector_<SpatialVec>&  bodyForces,
        Vector_<SpatialVec>&        allZ,
        Vector&                     jointForces) const 
    {
        const SpatialVec& myBodyForce  = fromB(bodyForces);
        SpatialVec&       z            = toB(allZ);

        // Centrifugal forces are Pa+b where P is articulated body inertia,
        // a is total coriolis acceleration, and b is gyroscopic force.
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

