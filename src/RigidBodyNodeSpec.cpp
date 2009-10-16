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
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    HType&                      H_PB_G) const
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
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    const SBTreeVelocityCache&  vc,
    HType&                      HDot_PB_G) const
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
    const SBTreePositionCache& pc = sbs.updTreePositionCache();

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
    const SBTreePositionCache& pc = sbs.getTreePositionCache();
    // Must use "upd" here rather than "get" because this is
    // called during realize(Velocity).
    const SBTreeVelocityCache& vc = sbs.updTreeVelocityCache();

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
//   tauBar (G*~H - I, a temporary not reused elsewhere)
//      Psi (Phi*(G*~H - I), articulated body child-to-parent shift matrix)
// and put them in the state cache.
// This must be called tip-to-base (inward).
//
template<int dof> void
RigidBodyNodeSpec<dof>::realizeArticulatedBodyInertiasInward(
    const SBInstanceCache&          ic,
    const SBTreePositionCache&      pc,
    SBArticulatedBodyInertiaCache&  abc) const 
{
    SpatialMat& P = updP(abc);
    // Start with the spatial inertia of the current body.
    P = getMk(pc);

    // For each child, take its articulated body inertia P and 
    // remove the portion of that inertia that can't be felt from 
    // the parent because of the joint mobilities. That is, we
    // calculate Pnew = P - P H DI ~H P, where I believe the
    // second term is the projection of P into the mobility space.
    // Then we shift Pnew from child to parent: Pparent += Phi*Pnew*~Phi.
    // TODO: can this be optimized for the common case where the
    // child is a terminal body and hence its P is an ordinary
    // rigid body inertia?
    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild    = children[i]->getPhi(pc);
        const SpatialMat& PChild      = children[i]->getP(abc);
        const SpatialMat& psiChild    = children[i]->getPsi(abc);

        // If the accelerations are prescribed, this is just a
        // rigid body (composite) inertia shift, which is much
        // cheaper than an articulated body shift.
        // Can be done in 93 flops (72 for the shift and 21 to add it in).
        if (children[i]->isUDotKnown(ic)) {
            P += phiChild * (PChild * ~phiChild); // ~250 flops; TODO: symmetric result
            continue;
        }

        // TODO: this is around 550 flops, 396 due to the 6x6
        // multiply with Psi, about 100 for the P*Phi shift, and 36
        // for the -=. The 6x6 multiply yields a symmetric result
        // so can be cut almost in half:
        // Psi is -Phi*(1 - P H DI ~H) so the multiply by P*~Phi gives
        // the negative of the projected result we want,
        // Phi*(P - P H DI ~H P)*~Phi. Subtracting here fixes the sign.
        P -= psiChild * (PChild * ~phiChild);
    }

    // If this is a prescribed mobilizer, leave G, D, DI, tauBar, and psi
    // untouched -- they should have been set to NaN at Instance stage.
    if (isUDotKnown(ic))
        return;

    const HType&  H  = getH(pc);
    HType&        G  = updG(abc);
    Mat<dof,dof>& D  = updD(abc);
    Mat<dof,dof>& DI = updDI(abc);

    const Mat<2,dof,Vec3> PH = P * H;       // 66*dof   flops
    D  = ~H * PH;                           // 11*dof^2 flops (symmetric result)

    // this will throw an exception if the matrix is ill conditioned
    DI = D.invert();                        // ~dof^3 flops (symmetric)
    G  = PH * DI;                           // 12*dof^2-6*dof flops

    // Note: Jain has TauBar=I-G*~H but we're negating that to G*~H-I because we
    // can save 30 flops by just subtracting 1 from the diagonals rather than having
    // to negate all the off-diagonals. Then Psi ends up with the wrong sign here
    // also, which gets handled in the shifting loop above by subtracting instead
    // of adding. So the P's we calculate all have the right sign in the end.
    SpatialMat& tauBar = updTauBar(abc);
    tauBar = G * ~H;                        // 11*dof^2 flops
    tauBar(0,0) -= 1; // subtract identity matrix (only touches diagonals -- 3 flops)
    tauBar(1,1) -= 1; //    "    (3 flops)
    updPsi(abc)  = getPhi(pc) * tauBar;     // ~100 flops
}

// To be called base to tip.
// sherm 060723: As best I can tell this is calculating the inverse of
// the "operational space inertia" at the body frame origin for each body.
// See Equation 20 in Rodriguez,Jain, & Kreutz-Delgado: A spatial operator algebra 
// for manipulator modeling and control. Intl. J. Robotics Research 
// 10(4):371-381 (1991).
template<int dof> void
RigidBodyNodeSpec<dof>::realizeYOutward
   (const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    SBDynamicsCache&                        dc) const
{
    if (isUDotKnown(ic)) {
        //TODO: (sherm 090810) is this right?
        assert(false);
        //updY(dc) = (~getPhi(pc) * parent->getY(dc)) * getPhi(pc); // rigid shift
        return;
    }

    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?


    // Psi here has the opposite sign from Jain's, but we're multiplying twice
    // by it here so it doesn't matter.
    updY(dc) = (getH(pc) * getDI(abc) * ~getH(pc)) 
                + (~getPsi(abc) * parent->getY(dc) * getPsi(abc));
}

//
// To be called from tip to base.
//
template<int dof> void
RigidBodyNodeSpec<dof>::realizeZ(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    const Vector&                           mobilityForces,
    const Vector_<SpatialVec>&              bodyForces) const 
{
    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - fromB(bodyForces);

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = children[i]->getZ(ac);
        const SpatialVec& GepsChild = children[i]->getGepsilon(ac);

        z += phiChild * (zChild + GepsChild);
    }

    updEpsilon(ac)  = fromU(mobilityForces) - ~getH(pc)*z;
    updGepsilon(ac) = getG(abc)  * getEpsilon(ac);
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were fed to realizeZ (as embodied in 'epsilon').
// (Base to tip)
//
template<int dof> void 
RigidBodyNodeSpec<dof>::realizeAccel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    Vector&                                 allUdot) const 
{
    Vec<dof>&        udot   = toU(allUdot);
    const SpatialVec A_GP = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

    udot        = getDI(abc) * getEpsilon(ac) - (~getG(abc)*A_GP);
    updA_GB(ac) = A_GP + getH(pc)*udot + getCoriolisAcceleration(vc);  
}

 
//
// To be called from tip to base.
// Temps do not need to be initialized.
//
template<int dof> void
RigidBodyNodeSpec<dof>::calcUDotPass1Inward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Vector&                           jointForces,
    const Vector_<SpatialVec>&              bodyForces,
    const Vector&                           allUDot,
    Vector_<SpatialVec>&                    allZ,
    Vector_<SpatialVec>&                    allGepsilon,
    Vector&                                 allEpsilon) const 
{
    const Vec<dof>&   myJointForce = fromU(jointForces);
    const SpatialVec& myBodyForce  = fromB(bodyForces);
    SpatialVec&       z            = toB(allZ);
    SpatialVec&       Geps         = toB(allGepsilon);
    Vec<dof>&         eps          = toU(allEpsilon);

    // Pa+b - F (TODO: include P*udot_p here also?)
    z = getCentrifugalForces(dc) - myBodyForce; // 6 flops

    if (isUDotKnown(ic) && !isUDotKnownToBeZero(ic)) {
        const Vec<dof>& udot = fromU(allUDot);
        z += getP(abc)*(getH(pc)*udot); // add in P*udot_p (66+12*dof flops)
    }

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];

        if (children[i]->isUDotKnown(ic))
            z += phiChild * zChild;
        else {
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];
            z += phiChild * (zChild + GepsChild);
        }
    }

    eps  = myJointForce - ~getH(pc)*z; // 12*dof flops

    if (!isUDotKnown(ic))
        Geps = getG(abc) * eps; // 12*dof-6 flops
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcUDotPass2Outward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    const Vector&                           allEpsilon,
    Vector_<SpatialVec>&                    allA_GB,
    Vector&                                 allUDot,
    Vector&                                 allTau) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot);    // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    if (isUDotKnown(ic)) {
        Vec<dof>& tau = updTau(ic,allTau);  // pull out this node's tau
        tau = ~getH(pc)*(getP(abc)*A_GP) - eps; // 66 + 12*dof flops
    } else
        udot = getDI(abc) * eps - (~getG(abc)*A_GP); // 2*dof^2 + 11*dof

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
    const SBTreePositionCache&           pc,
    const SBArticulatedBodyInertiaCache& abc,
    const SBDynamicsCache&               dc,
    const Vector&                        f,
    Vector_<SpatialVec>&                 allZ,
    Vector_<SpatialVec>&                 allGepsilon,
    Vector&                              allEpsilon) const 
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

        z += phiChild * (zChild + GepsChild);   // 24 flops
    }

    eps  = myJointForce - ~getH(pc)*z;          // 12*dof flops
    Geps = getG(abc) * eps;                     // 12*dof-6 flops
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof> void 
RigidBodyNodeSpec<dof>::calcMInverseFPass2Outward(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Vector&                           allEpsilon,
    Vector_<SpatialVec>&                    allA_GB,
    Vector&                                 allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = toB(allA_GB);
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()]; // 12 flops

    udot = getDI(abc) * eps - (~getG(abc)*A_GP);    // 2dof^2 + 11 dof flops
    A_GB = A_GP + getH(pc)*udot;                    // 12 dof flops
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
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const Vector&               allUDot,
    Vector_<SpatialVec>&        allA_GB) const
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
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
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
// Cost is: pass1 12*dof+12 flops
//          pass2 11*dof+84 flops (TODO: could be 24 flops cheaper; see below).
//          total 23*dof + 96     (TODO: 23*dof + 72)
// TODO: Possibly this could be much faster if composite body inertias R were
// precalculated along with D = ~H*R*H. Then I *think* this could be a single
// pass with tau = D*udot. Whether this is worth it would depend on how much
// this gets re-used after R and D are calculated.
template<int dof> void 
RigidBodyNodeSpec<dof>::calcMVPass1Outward(
    const SBTreePositionCache&  pc,
    const Vector&               allUDot,
    Vector_<SpatialVec>&        allA_GB) const
{
    const Vec<dof>& udot = fromU(allUDot);
    SpatialVec&     A_GB = toB(allA_GB);

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()]; // 12 flops

    A_GB = A_GP + getH(pc)*udot;    // 12*dof flops
}

// Call tip to base after calling calcMVPass1Outward() for each
// rigid body.
template<int dof> void
RigidBodyNodeSpec<dof>::calcMVPass2Inward(
    const SBTreePositionCache&  pc,
    const Vector_<SpatialVec>&  allA_GB,
    Vector_<SpatialVec>&        allF,	// temp
    Vector&                     allTau) const 
{
    const SpatialVec& A_GB  = fromB(allA_GB);
    SpatialVec&       F		= toB(allF);
    Vec<dof>&         tau	= toU(allTau);

    F = getMk(pc)*A_GB; // 66 flops TODO: 42 if you take advantage of Mk's structure

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
        F += phiChild * FChild; // 18 flops
    }

    tau = ~getH(pc)*F;          // 11*dof flops
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
    const SBTreePositionCache&  pc,
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
    const SBTreePositionCache&  pc,
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
    const SBTreePositionCache&  pc,
    const SBDynamicsCache&      dc,
    const Vector_<SpatialVec>&  bodyForces,
    Vector_<SpatialVec>&        allZ,
    Vector&                     jointForces) const 
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
    const SBTreePositionCache&  pc, 
    const SBTreeVelocityCache&  mc,
    const SpatialVec&           sVel, 
    Vector&                     u) const 
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
