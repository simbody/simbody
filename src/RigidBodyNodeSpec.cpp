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
 * This file contains implementations of the base class methods for the
 * templatized class RigidBodyNodeSpec<dof>, and instantiations of the 
 * class for all possible values of dof (1-6).
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

#include "RigidBodyNodeSpec_Pin.h"
#include "RigidBodyNodeSpec_Slider.h"
#include "RigidBodyNodeSpec_Ball.h"
#include "RigidBodyNodeSpec_Free.h"
#include "RigidBodyNodeSpec_Custom.h"


// Same for all mobilizers.
// CAUTION: our H matrix definition is transposed from Jain and Schwieters.
template<int dof> void
RigidBodyNodeSpec<dof>::calcParentToChildVelocityJacobianInGround(
    const SBModelVars&     mv,
    const SBPositionCache& pc, 
    HType& H_PB_G) const
{
    const HType& H_FM = getH_FM(pc);

    // want r_MB_F, that is, the vector from OM to OB, expressed in F
    const Vec3&     r_MB   = getX_MB().p();     // fixed
    const Rotation& R_FM   = getX_FM(pc).R();   // just calculated
    const Vec3      r_MB_F = R_FM*r_MB;         // 15 flops

    HType H_MB;
    H_MB[0] = Vec3(0); // fills top row with zero
    H_MB[1] = ~crossMat(r_MB_F) * H_FM[0];      // 15*dof + 3 flops

    // Now we want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();      // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GF = R_GP * R_PF;     // 45 flops

    H_PB_G =  R_GF * (H_FM + H_MB);         // 36*dof flops
}

// Same for all mobilizers.
// CAUTION: our H matrix definition is transposed from Jain and Schwieters.
template<int dof> void
RigidBodyNodeSpec<dof>::calcParentToChildVelocityJacobianInGroundDot(
    const SBModelVars&     mv,
    const SBPositionCache& pc, 
    const SBVelocityCache& vc,
    HType& HDot_PB_G) const
{
    const HType& H_FM     = getH_FM(pc);
    const HType& HDot_FM = getHDot_FM(vc);

    HType H_MB, HDot_MB;

    // want r_MB_F, that is, the vector from OM to OB, expressed in F
    const Vec3&     r_MB   = getX_MB().p();     // fixed
    const Rotation& R_FM   = getX_FM(pc).R();   // just calculated
    const Vec3      r_MB_F = R_FM*r_MB;         // 15 flops

    const Vec3& w_FM = getV_FM(vc)[0]; // local angular velocity

    H_MB[0] = Vec3(0); // fills top row with zero
    H_MB[1] = ~crossMat(r_MB_F) * H_FM[0];      // 15*dof + 3 flops

    HDot_MB[0] = Vec3(0);
    HDot_MB[1] =   ~crossMat(r_MB_F)        * HDot_FM[0] // 30*dof + 18 flops
                 + ~crossMat(w_FM % r_MB_F) * H_FM[0];

    // Now we want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();       // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R();     // parent orientation in ground
    const Rotation  R_GF = R_GP * R_PF;         // 45 flops

    const Vec3& w_GF = getV_GP(vc)[0]; // F and P have same angular velocity

    // Note: time derivative of R_GF is crossMat(w_GF)*R_GF.
    //      H = H_PB_G =  R_GF * (H_FM + H_MB) (see above method)
    const HType& H_PB_G = getH(pc);
    HDot_PB_G = R_GF * (HDot_FM + HDot_MB) // 66*dof + 3 flops
                 + crossMat(w_GF) * H_PB_G;
}

// This is the default implementation for turning H_MF into H_FM. 
// A mobilizer can override this to do it faster.
// From the Simbody theory manual,
//      H_FM_w = - R_FM*H_MF_w
//      H_FM_v = -(R_FM*H_MF_v + p_FM x H_FM_w)
//             
template<int dof> void
RigidBodyNodeSpec<dof>::calcReverseMobilizerH_FM(
    const SBStateDigest& sbs,
    HType&               H_FM) const
{
    // Must use "upd" here rather than "get" because this is
    // called during realize(Position).
    const SBPositionCache& pc = sbs.updPositionCache();

    HType H_MF;
    calcAcrossJointVelocityJacobian(sbs, H_MF);

    const Rotation& R_FM = getX_FM(pc).R();
    const Vec3&     p_FM = getX_FM(pc).p();

    // Make cross product matrix (3 flops). Saves a few flops (3*dof) to
    // transpose this and get the negative of the cross product matrix.
    const Mat33     np_FM_x  = ~crossMat(p_FM);

    H_FM[0] = -(R_FM * H_MF[0]);                // 18*dof flops
    H_FM[1] = np_FM_x * H_FM[0] - R_FM*H_MF[1]; // 33*dof flops
}

// This is the default implementation for turning HDot_MF into HDot_FM. 
// A mobilizer can override this to do it faster.
// We depend on H_FM having already been calculated.
//
// From the Simbody theory manual,
//      HDot_FM_w = -R_FM * HDot_MF_w + w_FM_x H_FM_w
//                  
//      HDot_FM_v = -R_FM * HDot_MF_v + w_FM_x H_FM_v 
//                  - (p_FM_x HDot_FM_w
//                     + (v_FM_x - w_FM_x p_FM_x)H_FM_w)
//
// where "a_x" indicates the cross product matrix of vector a.
//  
template<int dof> void
RigidBodyNodeSpec<dof>::calcReverseMobilizerHDot_FM(
    const SBStateDigest& sbs,
    HType&               HDot_FM) const
{
    const SBPositionCache& pc = sbs.getPositionCache();
    // Must use "upd" here rather than "get" because this is
    // called during realize(Velocity).
    const SBVelocityCache& vc = sbs.updVelocityCache();

    HType HDot_MF;
    calcAcrossJointVelocityJacobianDot(sbs, HDot_MF);

    const Rotation& R_FM    = getX_FM(pc).R();
    const Vec3&     p_FM    = getX_FM(pc).p();
    const HType&    H_FM    = getH_FM(pc);

    const Vec3&     w_FM    = getV_FM(vc)[0];
    const Vec3&     v_FM    = getV_FM(vc)[1];
    
    // Make cross product matrices.
    const Mat33     p_FM_x  = crossMat(p_FM);   // 3 flops
    const Mat33     w_FM_x  = crossMat(w_FM);   // 3 flops
    const Mat33     v_FM_x  = crossMat(v_FM);   // 3 flops
    const Mat33     vwp     = v_FM_x - w_FM_x*p_FM_x;   // 54 flops

    // Initialize both rows with the first two terms above.
    HDot_FM = w_FM_x*H_FM - R_FM*HDot_MF;               // 66*dof flops

    // Add in the additional terms in the second row.
    HDot_FM[1] -= p_FM_x * HDot_FM[0] + vwp * H_FM[0];  // 36*dof flops
}


//
// Given only position-related quantities from the State 
//      Mk  (this body's spatial inertia matrix)
//      Phi (composite body child-to-parent shift matrix)
//      H   (joint transition matrix; sense is transposed from Jain and Schwieters)
// we calculate dynamic quantities 
//      P   (articulated body inertia)
//      D   (factored mass matrix LDL' diagonal part D=~H*P*H)
//      DI  (inverse of D)
//      G   (P * H * DI)
//   tauBar (I-G*~H, a temporary not reused elsewhere)
//      Psi (Phi*(I-G*~H), articulated body child-to-parent shift matrix)
// and put them in the state cache.
// This must be called tip-to-base (inward).
//
template<int dof> void
RigidBodyNodeSpec<dof>::realizeArticulatedBodyInertiasInward(
    const SBPositionCache& pc,
    SBDynamicsCache&       dc) const 
{
    SpatialMat& P = updP(dc);
    // Start with the spatial inertia of the current body.
    P = getMk(pc);

    // For each child, take its articulated body inertia P and 
    // remove the portion of that inertia that can't be felt from 
    // the parent because of the joint mobilities. That is, we
    // calculate Pnew = P - P H DI ~H P, where I believe the
    // second term is the projection of P into the mobility space.
    // Then we shift Pnew from child to parent: Pparent += Phi*Pnew*~Phi.
    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild    = children[i]->getPhi(pc);
        const SpatialMat& PChild      = children[i]->getP(dc);
        const SpatialMat& psiChild    = children[i]->getPsi(dc);

        // TODO: this is around 550 flops, 396 due to the 6x6
        // multiply with Psi, about 100 for the P*Phi shift, and 36
        // for the +=. The 6x6 multiply yields a symmetric result
        // so can be cut almost in half:
        // Psi is -Phi*(1 - P H DI ~H) so the multiply by P*Phi gives
        // the negative of the projected result we want,
        // Phi*(P - P H DI ~H P)*Phi. Subtracting here fixes the sign.
        P -= psiChild * (PChild * ~phiChild);
    }

    const HType&  H  = getH(pc);
    HType&        G  = updG(dc);
    Mat<dof,dof>& D  = updD(dc);
    Mat<dof,dof>& DI = updDI(dc);

    const Mat<2,dof,Vec3> PH = P * H;       // 66*dof   flops
    D  = ~H * PH;                           // 11*dof^2 flops (symmetric result)

    // this will throw an exception if the matrix is ill conditioned
    DI = D.invert();                        // ~dof^3 flops (symmetric)
    G  = PH * DI;                           // 11*dof^2 flops

    // Note: Jain has TauBar=I-G*~H but we're negating that to G*~H-I because we
    // can save 30 flops by just subtracting 1 from the diagonals rather than having
    // to negate all the off-diagonals. Then Psi ends up with the wrong sign here
    // also, which gets handled in the shifting loop above by subtracting instead
    // of adding. So the P's we calculate all have the right sign in the end.
    SpatialMat& tauBar = updTauBar(dc);
    tauBar = G * ~H;                        // 11*dof^2 flops
    tauBar(0,0) -= 1; // subtract identity matrix (only touches diagonals -- 3 flops)
    tauBar(1,1) -= 1; //    "    (3 flops)
    updPsi(dc)   = getPhi(pc) * tauBar;     // ~100 flops
}

// To be called base to tip.
// sherm 060723: As best I can tell this is calculating the inverse of
// the "operational space inertia" at the body frame origin for each body.
// See Equation 20 in Rodriguez,Jain, & Kreutz-Delgado: A spatial operator algebra 
// for manipulator modeling and control. Intl. J. Robotics Research 
// 10(4):371-381 (1991).
template<int dof> void
RigidBodyNodeSpec<dof>::realizeYOutward(
    const SBPositionCache& pc,
    SBDynamicsCache&       dc) const 
{
    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?

    // Psi here has the opposite sign from Jain's, but we're multiplying twice
    // by it here so it doesn't matter.
    updY(dc) = (getH(pc) * getDI(dc) * ~getH(pc)) 
                + (~getPsi(dc) * parent->getY(dc) * getPsi(dc));
}

//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::realizeZ(
    const SBStateDigest&       sbs,
    const Vector&              mobilityForces,
    const Vector_<SpatialVec>& bodyForces) const 
{
    const SBPositionCache& pc = sbs.getPositionCache();
    const SBDynamicsCache& dc = sbs.getDynamicsCache();
    SBAccelerationCache&   ac = sbs.updAccelerationCache();

    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - fromB(bodyForces);

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = children[i]->getZ(ac);
        const SpatialVec& GepsChild = children[i]->getGepsilon(ac);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(ac)  = fromU(mobilityForces) - ~getH(pc)*z;
    updNu(ac)       = getDI(dc) * getEpsilon(ac);
    updGepsilon(ac) = getG(dc)  * getEpsilon(ac);
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to realizeZ (as embodied in 'nu').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::realizeAccel(
    const SBStateDigest&   sbs,
    Vector&                allUdot,
    Vector&                allQdotdot) const 
{
    const SBPositionCache& pc = sbs.getPositionCache();
    const SBVelocityCache& vc = sbs.getVelocityCache();
    const SBDynamicsCache& dc = sbs.getDynamicsCache();
    SBAccelerationCache&   ac = sbs.updAccelerationCache();

    Vec<dof>&        udot   = toU(allUdot);
    const SpatialVec A_GP = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

    udot        = getNu(ac) - (~getG(dc)*A_GP);
    updA_GB(ac) = A_GP + getH(pc)*udot + getCoriolisAcceleration(vc);  

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

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - ~getH(pc)*z;
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
    const SBVelocityCache& vc,
    const SBDynamicsCache& dc,
    const Vector&          allEpsilon,
    Vector_<SpatialVec>&   allA_GB,
    Vector&                allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + getH(pc)*udot + getCoriolisAcceleration(vc);  
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

    z = 0;

    for (unsigned i=0; i<children.size(); i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];

        z += phiChild * (zChild + GepsChild);
    }

    eps  = myJointForce - ~getH(pc)*z;
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
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    udot = getDI(dc) * eps - (~getG(dc)*A_GP);
    A_GB = A_GP + getH(pc)*udot;  
}

// The next two methods calculate 
//      f = M*udot + C(u) - f_applied 
// in O(N) time, where C(u) are the velocity-dependent coriolis, centrifugal and gyroscopic forces,
// f_applied are the applied body and joint forces, and udot is given.
//
// Pass 1 is base to tip and calculates the body accelerations arising
// from udot and the coriolis accelerations.
template<int dof> void 
RigidBodyNodeSpec<dof>::calcInverseDynamicsPass1Outward(
    const SBPositionCache& pc,
    const SBVelocityCache& vc,
    const Vector&          allUDot,
    Vector_<SpatialVec>&   allA_GB) const
{
    const Vec<dof>& udot = fromU(allUDot);
    SpatialVec&     A_GB = toB(allA_GB);

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    A_GB = A_GP + getH(pc)*udot + getCoriolisAcceleration(vc); 
}

// Pass 2 is tip to base. It takes body accelerations from pass 1 (including coriolis
// accelerations as well as hinge accelerations), and the applied forces,
// and calculates the additional hinge forces that would be necessary to produce 
// the observed accelerations.
template<int dof> void
RigidBodyNodeSpec<dof>::calcInverseDynamicsPass2Inward(
    const SBPositionCache&      pc,
    const SBVelocityCache&      vc,
    const Vector_<SpatialVec>&  allA_GB,
    const Vector&               jointForces,
    const Vector_<SpatialVec>&  bodyForces,
    Vector_<SpatialVec>&        allF,	// temp
    Vector&                     allTau) const 
{
    const Vec<dof>&   myJointForce  = fromU(jointForces);
    const SpatialVec& myBodyForce   = fromB(bodyForces);
    const SpatialVec& A_GB          = fromB(allA_GB);
    SpatialVec&       F		        = toB(allF);
    Vec<dof>&         tau	        = toU(allTau);

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

    // Project body forces into hinge space and subtract any hinge forces already
    // being applied to get the remaining hinge forces needed.
    tau = ~getH(pc)*F - myJointForce;
}

// The next two methods calculate x=M*v (or f=M*a) in O(N) time, given v.
// The first one is called base to tip.
template<int dof> void 
RigidBodyNodeSpec<dof>::calcMVPass1Outward(
    const SBPositionCache& pc,
    const Vector&          allUDot,
    Vector_<SpatialVec>&   allA_GB) const
{
    const Vec<dof>& udot = fromU(allUDot);
    SpatialVec&     A_GB = toB(allA_GB);

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    A_GB = A_GP + getH(pc)*udot;  
}

// Call tip to base after calling calcMVPass1Outward() for each
// rigid body.
template<int dof> void
RigidBodyNodeSpec<dof>::calcMVPass2Inward(
    const SBPositionCache& pc,
    const Vector_<SpatialVec>& allA_GB,
    Vector_<SpatialVec>&       allF,	// temp
    Vector&                    allTau) const 
{
    const SpatialVec& A_GB  = fromB(allA_GB);
    SpatialVec&       F		= toB(allF);
    Vec<dof>&         tau	= toU(allTau);

    F = getMk(pc)*A_GB;

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
        F += phiChild * FChild;
    }

    tau = ~getH(pc)*F;
}

//
// Calculate product of kinematic Jacobian J=~Phi*H and a mobility-space vector. Requires 
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
    const SpatialVec outP = ~getPhi(pc) * parent->fromB(Jv);

    out = outP + getH(pc)*in;  
}

//
// Calculate product of kinematic Jacobian transpose ~J=~H*Phi and a gradient vector on each of the
// outboard bodies. Requires that Phi and H are available, so this
// should only be called in Stage::Position or higher. This does not change the cache at all.
// NOTE (sherm 060214): I reworked this from the original. This one no longer incorporates
// applied hinge gradients if there are any; just add those in at the end if you want them.
//
// (sherm 060727) In spatial operators, this calculates ~H*Phi*F where F are the spatial forces
// applied to each body. See Schwieters Eq. 41. (But with sense of H transposed.)
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

    for (unsigned i=0; i<children.size(); ++i) {
        const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild = children[i]->getPhi(pc);

        z += phiChild * zChild;
    }

    out = ~getH(pc) * z; 
}

//
// To be called from tip to base.
// Temps do not need to be initialized.
// (sherm 060727) In spatial operators, this calculates ~H*Phi*(F-(Pa+b))
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

    for (unsigned i=0; i<children.size(); ++i) {
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);

        z += phiChild * zChild; 
    }

    eps  = ~getH(pc) * z;
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
    toU(u) = ~getH(pc) * (sVel - (~getPhi(pc) * parent->getV_GB(mc)));
}

	////////////////////
	// INSTANTIATIONS //
	////////////////////

template class RigidBodyNodeSpec<1>;
template class RigidBodyNodeSpec<2>;
template class RigidBodyNodeSpec<3>;
template class RigidBodyNodeSpec<4>;
template class RigidBodyNodeSpec<5>;
template class RigidBodyNodeSpec<6>;
