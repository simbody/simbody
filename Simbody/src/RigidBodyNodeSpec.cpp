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
 * templatized class RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>, and instantiations of the 
 * class for all possible values of the arguments.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"

#include "RigidBodyNodeSpec_Pin.h"
#include "RigidBodyNodeSpec_Slider.h"
#include "RigidBodyNodeSpec_Ball.h"
#include "RigidBodyNodeSpec_Free.h"
#include "RigidBodyNodeSpec_Custom.h"

//------------------------------------------------------------------------------
//                              CALC H_PB_G
//------------------------------------------------------------------------------
// Same for all mobilizers.
// CAUTION: our H matrix definition is transposed from Jain and Schwieters.
// Cost: 60 + 45*dof flops
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcParentToChildVelocityJacobianInGround(
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    HType&                      H_PB_G) const
{
    const HType& H_FM = getH_FM(pc);

    // We want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();      // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R(); // parent orientation in ground
    const Rotation  R_GF = (noR_PF ? R_GP : R_GP * R_PF);     // 45 flops

    if (noX_MB || noR_FM)
        H_PB_G = R_GF * H_FM;       // 3*dof flops
    else {
        // want r_MB_F, that is, the vector from Mo to Bo, expressed in F
        const Vec3&     r_MB   = getX_MB().p();     // fixed
        const Rotation& R_FM   = getX_FM(pc).R();   // just calculated
        const Vec3      r_MB_F = (noR_FM ? r_MB : R_FM*r_MB);         // 15 flops
        HType H_MB_F;
        H_MB_F[0] =  Vec3(0); // fills top row with zero
        H_MB_F[1] = -r_MB_F % H_FM[0]; // 9*dof (negation not actually done)
        H_PB_G = R_GF * (H_FM + H_MB_F); // 36*dof flops
    }
}

//------------------------------------------------------------------------------
//                            CALC H_PB_G_DOT
//------------------------------------------------------------------------------
// Same for all mobilizers.
// CAUTION: our H matrix definition is transposed from Jain and Schwieters. 
// Cost is 69 + 65*dof flops
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcParentToChildVelocityJacobianInGroundDot(
    const SBModelVars&          mv,
    const SBTreePositionCache&  pc, 
    const SBTreeVelocityCache&  vc,
    HType&                      HDot_PB_G) const
{
    const HType& H_FM    = getH_FM(pc);
    const HType& HDot_FM = getHDot_FM(vc);

    // We want R_GF so we can reexpress the cross-joint velocity V_FB (==V_PB)
    // in the ground frame, to get V_PB_G.

    const Rotation& R_PF = getX_PF().R();       // fixed config of F in P

    // Calculated already since we're going base to tip.
    const Rotation& R_GP = getX_GP(pc).R();     // parent orientation in ground
    const Rotation  R_GF = (noR_PF ? R_GP : R_GP * R_PF); // 45 flops (TODO: again??)

    const Vec3& w_GF = getV_GP(vc)[0]; // F and P have same angular velocity

    // Note: time derivative of R_GF is crossMat(w_GF)*R_GF.
    //      H = H_PB_G =  R_GF * (H_FM + H_MB_F) (see above method)
    const HType& H_PB_G = getH(pc);
    if (noX_MB || noR_FM)
        HDot_PB_G = R_GF * HDot_FM // 48*dof
                  + HType(w_GF % H_PB_G[0], 
                          w_GF % H_PB_G[1]);
    else {
        // want r_MB_F, that is, the vector from OM to OB, expressed in F 
        const Vec3&     r_MB   = getX_MB().p();     // fixed
        const Rotation& R_FM   = getX_FM(pc).R();   // just calculated
        const Vec3      r_MB_F = (noR_FM ? r_MB : R_FM*r_MB); // 15 flops

        const Vec3& w_FM = getV_FM(vc)[0]; // local angular velocity

        HType HDot_MB_F;
        HDot_MB_F[0] = Vec3(0);
        HDot_MB_F[1] =          -r_MB_F  % HDot_FM[0] // 21*dof + 9 flops
                       - (w_FM % r_MB_F) % H_FM[0];


        HDot_PB_G =   R_GF * (HDot_FM + HDot_MB_F) // 54*dof
                    + HType(w_GF % H_PB_G[0], 
                            w_GF % H_PB_G[1]);
    }
}

//------------------------------------------------------------------------------
//                       CALC REVERSE MOBILIZER H_FM
//------------------------------------------------------------------------------
// This is the default implementation for turning H_MF into H_FM. 
// A mobilizer can override this to do it faster.
// From the Simbody theory manual,
//      H_FM_w = - R_FM*H_MF_w
//      H_FM_v = -(R_FM*H_MF_v + p_FM x H_FM_w)
//             
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcReverseMobilizerH_FM(
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

    if (noR_FM) {
        H_FM[0] = -H_MF[0];
        H_FM[1] = np_FM_x * H_FM[0] - H_MF[1]; // 15*dof flops
    }
    else {
        H_FM[0] = -(R_FM * H_MF[0]);                // 18*dof flops
        H_FM[1] = np_FM_x * H_FM[0] - R_FM*H_MF[1]; // 33*dof flops
    }
}

//------------------------------------------------------------------------------
//                     CALC REVERSE MOBILIZER HDOT_FM
//------------------------------------------------------------------------------
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
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcReverseMobilizerHDot_FM(
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
    HDot_FM = w_FM_x*H_FM - (noR_FM ? HDot_MF : R_FM*HDot_MF); // 66*dof flops

    // Add in the additional terms in the second row.
    HDot_FM[1] -= p_FM_x * HDot_FM[0] + vwp * H_FM[0];  // 36*dof flops
}


//------------------------------------------------------------------------------
//                     REALIZE ARTICULATED BODY INERTIAS
//------------------------------------------------------------------------------
// Given only position-related quantities from the State 
//      Mk  (this body's spatial inertia matrix)
//      Phi (composite body child-to-parent shift matrix)
//      H   (joint transition matrix; sense is transposed from Jain & Schwieters)
// we calculate dynamic quantities 
//      P   (articulated body inertia)
//      D   (factored mass matrix LDL' diagonal part D=~H*P*H)
//      DI  (inverse of D)
//      G   (P * H * DI)
//    PPlus (P - P * H * DI * ~H * P)
// and put them in the state cache.
// This must be called tip-to-base (inward).
//
// Cost is 93 flops per child plus
//   n^3 + 23*n^2 + 115*n + 12
//   e.g. pin=143, ball=591 (197/dof), free=1746 (291/dof)
// Note that per-child cost is paid just once for each non-base body in
// the whole tree.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::realizeArticulatedBodyInertiasInward(
    const SBInstanceCache&          ic,
    const SBTreePositionCache&      pc,
    SBArticulatedBodyInertiaCache&  abc) const 
{
    ArticulatedInertia& P = updP(abc);

    // Start with the spatial inertia of the current body (in Ground frame).
    P = ArticulatedInertia(getMk_G(pc)); // 12 flops

    // For each child, we previously took its articulated body inertia P and 
    // removed the portion of that inertia that can't be felt from  the parent
    // because of the joint mobilities. That is, we calculated 
    // P+ = P - P H DI ~H P, where the second term is the projection of P into
    // the mobility space of the child's inboard mobilizer. Note that if the
    // child's mobilizer is prescribed, then the entire inertia will be felt
    // by the parent so P+ = P in that case. Now we're going to shift P+
    // from child to parent: Pparent += Phi*P+*~Phi.
    // TODO: can this be optimized for the common case where the
    // child is a terminal body and hence its P is an ordinary
    // spatial inertia? (Spatial inertia shift is 37 flops vs. 72 here; not
    // really much help.)
    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&          phiChild   = children[i]->getPhi(pc);
        const ArticulatedInertia& PChild     = children[i]->getP(abc);
        const ArticulatedInertia& PPlusChild = children[i]->getPPlus(abc);

        // Apply the articulated body shift.
        // This takes 93 flops (72 for the shift and 21 to add it in).
        if (children[i]->isUDotKnown(ic))
            P += PChild.shift(phiChild.l());
        else
            P += PPlusChild.shift(phiChild.l());
    }

    // If this is a prescribed mobilizer, leave G, D, DI, and PPlus
    // untouched -- they should have been set to NaN at Instance stage.
    if (isUDotKnown(ic))
        return;

    const HType&  H  = getH(pc);
    HType&        G  = updG(abc);
    Mat<dof,dof>& D  = updD(abc);
    Mat<dof,dof>& DI = updDI(abc);

    const HType PH = P.toSpatialMat() * H;  // 66*dof   flops
    D  = ~H * PH;                           // 11*dof^2 flops (symmetric result)

    // this will throw an exception if the matrix is ill conditioned
    DI = D.invert();                        // ~dof^3 flops (symmetric)
    G  = PH * DI;                           // 12*dof^2-6*dof flops

    // Want P+ = P - P H DI ~H P = P - G*~PH. 
    // We can do this in about 55*dof flops.
    ArticulatedInertia& PPlus = updPPlus(abc);
    // These require 9 dot products of length dof. The symmetric ones could
    // be done with 6 dot products instead for a small savings but this gives
    // us a chance to symmetrize and hopefully clean up some numerical errors.
    // The full price for all three is 54*dof-27 flops.
    Mat33 massMoment = G.row(0)*~PH.row(1); // (full)            9*(2*dof-1) flops
    Mat33 mass       = G.row(1)*~PH.row(1); // symmetric result  9*(2*dof-1) flops
    Mat33 inertia    = G.row(0)*~PH.row(0); // symmetric result  9*(2*dof-1) flops
    // These must be symmetrized due to numerical errors for 12 more flops. 
    SymMat33 symMass( mass(0,0), 
                     (mass(1,0)+mass(0,1))/2,  mass(1,1), 
                     (mass(2,0)+mass(0,2))/2, (mass(2,1)+mass(1,2))/2, mass(2,2));
    SymMat33 symInertia( 
         inertia(0,0), 
        (inertia(1,0)+inertia(0,1))/2,  inertia(1,1), 
        (inertia(2,0)+inertia(0,2))/2, (inertia(2,1)+inertia(1,2))/2, inertia(2,2));
    PPlus = P - ArticulatedInertia(symMass, massMoment, symInertia); // 21 flops
}

//------------------------------------------------------------------------------
//                                  REALIZE Y
//------------------------------------------------------------------------------
// To be called base to tip.
// sherm 060723: As best I can tell this is calculating the inverse of
// the "operational space inertia" at the body frame origin for each body.
// See Equation 20 in Rodriguez,Jain, & Kreutz-Delgado: A spatial operator algebra 
// for manipulator modeling and control. Intl. J. Robotics Research 
// 10(4):371-381 (1991).
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::realizeYOutward
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

    // Compute psi. Jain has TauBar=I-G*~H but we're negating that to G*~H-I because we
    // can save 30 flops by just subtracting 1 from the diagonals rather than having
    // to negate all the off-diagonals. Then Psi ends up with the wrong sign here
    // also, which doesn't matter because we multiply by it twice.

    SpatialMat tauBar = getG(abc)*~getH(pc);// 11*dof^2 flops
    tauBar(0,0) -= 1; // subtract identity matrix (only touches diagonals -- 3 flops)
    tauBar(1,1) -= 1; //    "    (3 flops)
    SpatialMat psi = getPhi(pc)*tauBar; // ~100 flops

    // TODO: this is very expensive (~1000 flops?) Could cut be at least half
    // by exploiting symmetry. Also, does Psi have special structure?
    // And does this need to be computed for every body or only those
    // which are loop "base" bodies or some such?


    // Psi here has the opposite sign from Jain's, but we're multiplying twice
    // by it here so it doesn't matter.
    updY(dc) = (getH(pc) * getDI(abc) * ~getH(pc)) 
                + (~psi * parent->getY(dc) * psi);
}

//------------------------------------------------------------------------------
//                  REALIZE Z and REALIZE ACCEL (obsolete)
//------------------------------------------------------------------------------
// To be called from tip to base.
//
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::realizeZ(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    const Real*                             mobilityForces,
    const SpatialVec*                       bodyForces) const 
{
    SpatialVec& z = updZ(ac);
    z = getCentrifugalForces(dc) - bodyForces[nodeNum];

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
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void 
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::realizeAccel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    SBTreeAccelerationCache&                ac,
    Real*                                   udot) const 
{
    const SpatialVec A_GP = ~getPhi(pc) * parent->getA_GB(ac); // ground A_GB is 0

    Vec<dof>::updAs(udot) = getDI(abc) * getEpsilon(ac) - (~getG(abc)*A_GP);
    updA_GB(ac) = A_GP + getH(pc)*Vec<dof>::getAs(udot) + getCoriolisAcceleration(vc);  
}

//------------------------------------------------------------------------------
//                                   CALC UDOT
//------------------------------------------------------------------------------
// To be called from tip to base.
// Temps do not need to be initialized.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcUDotPass1Inward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Real*                             jointForces,
    const SpatialVec*                       bodyForces,
    const Real*                             allUDot,
    SpatialVec*                             allZ,
    SpatialVec*                             allGepsilon,
    Real*                                   allEpsilon) const 
{
    const Vec<dof>&   myJointForce = fromU(jointForces);
    const SpatialVec& myBodyForce  = bodyForces[nodeNum];
    SpatialVec&       z            = allZ[nodeNum];
    SpatialVec&       Geps         = allGepsilon[nodeNum];
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
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void 
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcUDotPass2Outward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    const SBDynamicsCache&                  dc,
    const Real*                             allEpsilon,
    SpatialVec*                             allA_GB,
    Real*                                   allUDot,
    Real*                                   allTau) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = allA_GB[nodeNum];
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

 
//------------------------------------------------------------------------------
//                              CALC M INVERSE F
//------------------------------------------------------------------------------
// To be called from tip to base.
// Temps do not need to be initialized.
//
// This calculates udot = M^-1 f in two O(N) passes. Note that we are ignoring 
// velocities; if there are any velocity-dependent forces they should already be
// in f.
//
// When there is prescribed motion, it is already reflected in the articulated
// body inertias, which were formed with rigid body shifts across the prescribed
// mobilizers. In that case you can think of the udots and generalized forces f
// as partitioned into two sets each: udot={udot_p, udot_r}, f={f_p, f_r}. We
// can discuss them as though the partitions were contiguous although they are
// not and the code deals with that properly. So now you can view the system as
//     [ Mpp ~Mrp ] [udot_p]   [tau_p]   [f_p]
//     [ Mrp  Mrr ] [udot_r] + [  0  ] = [f_r]
// where udot_p is given and the unknowns are udot_r (the free accelerations)
// and tau_p (the unknown forces that enforce the prescribed motion). This 
// produces two equations when multiplied out:
//     (1) Mrr udot_r = f_r - Mrp udot_p
//     (2)     tau_p  = f_p - Mpp udot_p - ~Mrp udot_r 
//
// Now we can define what the present method does: it assumes udot_p is zero,
// and then calculates udot_r = Mrr^-1 f_r. f_p is ignored and won't be 
// examined; udot_p is ignored and won't be written.
//
// Cost per body is 
//      30 + 47*ndof_r + 2*ndof_r^2
// where ndof_r is the number of u's for a free mobilizer; 0 for a prescribed 
// one.
//
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcMInverseFPass1Inward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Real*                             f,
    SpatialVec*                             allZ,
    SpatialVec*                             allGepsilon,
    Real*                                   allEpsilon) const
{
    const Vec<dof>&   myJointForce = fromU(f);
    SpatialVec&       z            = allZ[nodeNum];
    SpatialVec&       Geps         = allGepsilon[nodeNum];
    Vec<dof>&         eps          = toU(allEpsilon);

    z = 0;

    for (unsigned i=0; i<children.size(); i++) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& zChild    = allZ[children[i]->getNodeNum()];

        if (children[i]->isUDotKnown(ic))
            z += phiChild * zChild;                 // 18 flops
        else {
            const SpatialVec& GepsChild = allGepsilon[children[i]->getNodeNum()];
            z += phiChild * (zChild + GepsChild);   // 24 flops
        }
    }

    if (!isUDotKnown(ic)) {
        eps  = myJointForce - ~getH(pc)*z;          // 12*dof flops
        Geps = getG(abc) * eps;                     // 12*dof-6 flops
    }
}

//
// Calculate acceleration in internal coordinates, based on the last set
// of forces that were reduced into epsilon (e.g., see above).
// Base to tip: temp allA_GB does not need to be initialized before
// beginning the iteration.
//
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void 
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcMInverseFPass2Outward(
    const SBInstanceCache&                  ic,
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBDynamicsCache&                  dc,
    const Real*                             allEpsilon,
    SpatialVec*                             allA_GB,
    Real*                                   allUDot) const
{
    const Vec<dof>& eps  = fromU(allEpsilon);
    SpatialVec&     A_GB = allA_GB[nodeNum];
    Vec<dof>&       udot = toU(allUDot); // pull out this node's udot

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    A_GB = ~getPhi(pc) * allA_GB[parent->getNodeNum()]; // 12 flops

    if (!isUDotKnown(ic)) {
        udot = getDI(abc) * eps - (~getG(abc)*A_GB);    // 2dof^2 + 11 dof flops
        A_GB += getH(pc)*udot;                          // 12 dof flops
    }
}

//------------------------------------------------------------------------------
//                            CALC INVERSE DYNAMICS
//------------------------------------------------------------------------------
// The next two methods calculate 
//      f = M*udot + C(u) - f_applied 
// in O(N) time, where C(u) are the velocity-dependent coriolis, centrifugal and 
// gyroscopic forces, f_applied are the applied body and joint forces, and udot 
// is given.
//      pass1: 12*dof + 18 flops
//      pass2: 12*dof + 75 flops
//      total: 24*dof + 93 flops
// Pass 1 is base to tip and calculates the body accelerations arising
// from udot and the coriolis accelerations.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void 
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcInverseDynamicsPass1Outward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const
{
    const Vec<dof>& udot = fromU(allUDot);
    SpatialVec&     A_GB = allA_GB[nodeNum];

    // Shift parent's A_GB outward. (Ground A_GB is zero.) 12 flops.
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()];

    A_GB = A_GP + getH(pc)*udot + getCoriolisAcceleration(vc); // 12*dof+6 flops.
}

// Pass 2 is tip to base. It takes body accelerations from pass 1 (including coriolis
// accelerations as well as hinge accelerations), and the applied forces,
// and calculates the additional hinge forces that would be necessary to produce 
// the observed accelerations.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcInverseDynamicsPass2Inward(
    const SBTreePositionCache&  pc,
    const SBTreeVelocityCache&  vc,
    const SpatialVec*           allA_GB,
    const Real*                 jointForces,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allF,   // temp
    Real*                       allTau) const 
{
    const Vec<dof>&   myJointForce  = fromU(jointForces);
    const SpatialVec& myBodyForce   = bodyForces[nodeNum];
    const SpatialVec& A_GB          = allA_GB[nodeNum];
    SpatialVec&       F             = allF[nodeNum];
    Vec<dof>&         tau           = toU(allTau);

    // Start with rigid body force from desired body acceleration and
    // gyroscopic forces due to angular velocity, minus external forces
    // applied directly to this body.
    F = getMk_G(pc)*A_GB + getGyroscopicForce(vc) - myBodyForce; // 57 flops

    // Add in forces on children, shifted to this body.
    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
        F += phiChild * FChild;         // 18 flops
    }

    // Project body forces into hinge space and subtract any hinge forces already
    // being applied to get the remaining hinge forces needed.
    tau = ~getH(pc)*F - myJointForce;   // 12*dof flops
}

//------------------------------------------------------------------------------
//                                 CALC M V
//------------------------------------------------------------------------------
// The next two methods calculate x=M*v (or f=M*a) in O(N) time, given v.
// The first one is called base to tip.
// Cost is: pass1 12*dof+12 flops
//          pass2 11*dof+63 flops
//          total 23*dof + 75
// TODO: Possibly this could be much faster if composite body inertias R were
// precalculated along with D = ~H*R*H. Then I *think* this could be a single
// pass with tau = D*udot. Whether this is worth it would depend on how much
// this gets re-used after R and D are calculated.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void 
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcMVPass1Outward(
    const SBTreePositionCache&  pc,
    const Real*                 allUDot,
    SpatialVec*                 allA_GB) const
{
    const Vec<dof>& udot = fromU(allUDot);
    SpatialVec&     A_GB = allA_GB[nodeNum];

    // Shift parent's A_GB outward. (Ground A_GB is zero.)
    const SpatialVec A_GP = ~getPhi(pc) * allA_GB[parent->getNodeNum()]; // 12 flops

    A_GB = A_GP + getH(pc)*udot;    // 12*dof flops
}

// Call tip to base after calling calcMVPass1Outward() for each
// rigid body.
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcMVPass2Inward(
    const SBTreePositionCache&  pc,
    const SpatialVec*           allA_GB,
    SpatialVec*                 allF,   // temp
    Real*                       allTau) const 
{
    const SpatialVec& A_GB  = allA_GB[nodeNum];
    SpatialVec&       F     = allF[nodeNum];
    Vec<dof>&         tau   = toU(allTau);

    // 45 flops because Mk has a nice structure
    F = getMk_G(pc)*A_GB; 

    for (unsigned i=0; i<children.size(); ++i) {
        const PhiMatrix&  phiChild  = children[i]->getPhi(pc);
        const SpatialVec& FChild    = allF[children[i]->getNodeNum()];
        F += phiChild * FChild; // 18 flops
    }

    tau = ~getH(pc)*F;          // 11*dof flops
}

//------------------------------------------------------------------------------
//                       MULTIPLY BY SYSTEM JACOBIAN
//------------------------------------------------------------------------------
// Calculate product of kinematic Jacobian J=~Phi*H and a mobility-space vector. Requires 
// that Phi and H are available, so this should only be called in Stage::Position or higher.
// This does not change the cache at all.
//
// Note that if the vector v==u (generalized speeds) then the result is V_GB=J*u, the
// body spatial velocities generated by those speeds.
//
// Call base to tip (outward).
//
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::
multiplyBySystemJacobian(
    const SBTreePositionCache&  pc,
    const Real*                 v,
    SpatialVec*                 Jv) const
{
    const Vec<dof>& in  = fromU(v);
    SpatialVec&     out = Jv[nodeNum];

    // Shift parent's result outward (ground result is 0).
    const SpatialVec outP = ~getPhi(pc) * Jv[parent->getNodeNum()]; // 12 flops

    out = outP + getH(pc)*in;  // 12*dof flops
}

//------------------------------------------------------------------------------
//                   MULTIPLY BY SYSTEM JACOBIAN TRANSPOSE
//------------------------------------------------------------------------------
// Calculate product of kinematic Jacobian transpose ~J=~H*Phi and a spatial
// forces vector on each of the outboard bodies. Requires that Phi and H are 
// available, so this should only be called in Stage::Position or higher. This 
// does not  change the cache at all.
// NOTE (sherm 060214): I reworked this from the original. This one no longer 
// incorporates applied hinge gradients if there are any; just add those in at 
// the end if you want them.
//
// (sherm 060727) In spatial operators, this calculates ~H*Phi*F where F are the 
// spatial forces applied to each body. See Schwieters Eq. 41. (But with sense 
// of H transposed.)
//
// Call tip to base.
//
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::
multiplyBySystemJacobianTranspose(
    const SBTreePositionCache&  pc,
    SpatialVec*                 zTmp,
    const SpatialVec*           X, 
    Real*                       JtX) const
{
    const SpatialVec& in  = X[getNodeNum()];
    Vec<dof>&         out = Vec<dof>::updAs(&JtX[getUIndex()]);
    SpatialVec&       z   = zTmp[getNodeNum()];

    z = in;

    for (unsigned i=0; i<children.size(); ++i) {
        const SpatialVec& zChild   = zTmp[children[i]->getNodeNum()];
        const PhiMatrix&  phiChild = children[i]->getPhi(pc);

        z += phiChild * zChild; // 18 flops
    }

    out = ~getH(pc) * z; // 11*dof flops
}

//------------------------------------------------------------------------------
//                       CALC EQUIVALENT JOINT FORCES
//------------------------------------------------------------------------------
// To be called from tip to base.
// Temps do not need to be initialized.
// (sherm 060727) In spatial operators, this calculates ~H*Phi*(F-(Pa+b))
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::calcEquivalentJointForces(
    const SBTreePositionCache&  pc,
    const SBDynamicsCache&      dc,
    const SpatialVec*           bodyForces,
    SpatialVec*                 allZ,
    Real*                       jointForces) const 
{
    const SpatialVec& myBodyForce  = bodyForces[nodeNum];
    SpatialVec&       z            = allZ[nodeNum];
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
template<int dof, bool noR_FM, bool noX_MB, bool noR_PF> void
RigidBodyNodeSpec<dof, noR_FM, noX_MB, noR_PF>::setVelFromSVel(
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

#define INSTANTIATE(dof, noR_FM) \
template class RigidBodyNodeSpec<dof, noR_FM, false, false>; \
template class RigidBodyNodeSpec<dof, noR_FM, false, true>; \
template class RigidBodyNodeSpec<dof, noR_FM, true, false>; \
template class RigidBodyNodeSpec<dof, noR_FM, true, true>;

INSTANTIATE(1, false)
INSTANTIATE(2, false)
INSTANTIATE(3, false)
INSTANTIATE(4, false)
INSTANTIATE(5, false)
INSTANTIATE(6, false)
INSTANTIATE(1, true)
INSTANTIATE(3, true)
