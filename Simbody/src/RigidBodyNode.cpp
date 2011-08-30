/* -------------------------------------------------------------------------- *
 *                              SimTK Simbody(tm)                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-11 Stanford University and the Authors.        *
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
 * This file contains implementations for non-inline methods of the
 * base RigidBodyNode class. The methods here are those for which we do
 * not need to know the number of mobilities associated with a node in order
 * to compute efficiently.
 */

#include "SimbodyMatterSubsystemRep.h"
#include "RigidBodyNode.h"
#include "RigidBodyNodeSpec.h"


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
    SBTreePositionCache&    pc) const
{
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 p_PB_G = getX_GP(pc).R() * getX_PB(pc).p(); // 15 flops

    // The Phi matrix conveniently performs child-to-parent (inward) shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // (outward) shifting for velocities and accelerations.
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

// Calculate velocity-related quantities: spatial velocity (V_GB), 
// and coriolis acceleration. This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G and
// velocity-dependent acceleration remainder term VD_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBTreePositionCache& pc,
    SBTreeVelocityCache&       vc) const
{
    if (nodeNum == 0) { // ground, just in case
        updV_GB(vc)                          = SpatialVec(Vec3(0), Vec3(0));
        updGyroscopicForce(vc)               = SpatialVec(Vec3(0), Vec3(0));
        updMobilizerCoriolisAcceleration(vc) = SpatialVec(Vec3(0), Vec3(0));
        updTotalCoriolisAcceleration(vc)     = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    const SpatialVec& V_GP   = parent->getV_GB(vc); // parent P's velocity in G
    const SpatialVec& V_PB_G = getV_PB_G(vc); // child B's vel in P, exp. in G
    const PhiMatrixTranspose PhiT = ~getPhi(pc); // shift outwards

    // calc spatial velocity of B's origin Bo in G (angular,linear)
    const SpatialVec V_GB = PhiT*V_GP + V_PB_G; // 18 flops
    updV_GB(vc) = V_GB;

    const Vec3& w_GB = V_GB[0]; // for convenience
    const Vec3& v_GB = V_GB[1];

    // Calculate gyroscopic moment and force (48 flops). Although this is 
    // really a dynamic quantity (requires spatial inertia and is itself
    // a force), we calculate this here rather than in stage Dynamics because 
    // it is needed in inverse dynamics but does not require articulated body 
    // inertias to be calculated. This could be deferred until needed but
    // probably isn't worth the trouble.
    updGyroscopicForce(vc) = getMass() *                   // 6 flops
        SpatialVec( w_GB % (getUnitInertia_OB_G(pc)*w_GB), // moment (24 flops)
                   (w_GB % (w_GB % getCB_G(pc))));         // force  (18 flops)


    // Parent velocity.
    const Vec3& w_GP = V_GP[0]; // for convenience
    const Vec3& v_GP = V_GP[1];

    // Calculate this mobilizer's incremental contribution to coriolis 
    // acceleration, and this body's total coriolis acceleration (it parent's 
    // coriolis  acceleration plus the incremental contribution).
    //
    // We just calculated 
    // (1)  V_GB = J*u = ~Phi * V_GP + H*u 
    // above. Eventually we will  want to compute 
    // (2)  A_GB = d/dt V_GB 
    // (3)       = J*udot + Jdot*u
    // (4)       = (~Phi*A_GP + H*udot) + (~Phidot*V_GP+Hdot*u).
    // (Don't be tempted to match the "J" terms in (3) with the two terms in (4)
    // because A_GP already includes coriolis terms up to the parent.)
    // That second term in (4) is just velocity dependent so we can calculate 
    // it here. That is what we're calling the "incremental contribution to 
    // coriolis acceleration" of mobilizer B.
    // CAUTION: our definition of H is transposed from Jain's and Schwieters'.
    //
    // So the incremental contribution to the coriolis acceleration is
    //   Amob = ~PhiDot * V_GP + HDot * u
    // As correctly calculated in Schwieters' paper, Eq [16], the first term 
    // above simplifies to SpatialVec( 0, w_GP % (v_GB-v_GP) ). However, 
    // Schwieters' second term in [16] is correct only if H is constant in P, 
    // in which case the derivative just accounts for the rotation of P in G. 
    // In general H is not constant in P, so we don't try to calculate the 
    // derivative here but assume that HDot*u has already been calculated for 
    // us and stored in VD_PB_G. (That is, V_PB_G = H*u, VD_PB_G = HDot*u.)
    //
    // Note: despite all the ground-relative velocities here, this is just
    // the contribution of the cross-joint velocity, but reexpressed in G.
    const SpatialVec& VD_PB_G = getVD_PB_G(vc);
    const SpatialVec  Amob(VD_PB_G[0], 
                           VD_PB_G[1] + w_GP % (v_GB-v_GP)); // 15 flops

    updMobilizerCoriolisAcceleration(vc) = Amob;

    // Finally, the total coriolis acceleration (normally just called "coriolis
    // acceleration"!) of body B is the total coriolis acceleration of its 
    // parent shifted outward, plus B's local contribution that we just 
    // calculated.
    updTotalCoriolisAcceleration(vc) =
        PhiT * parent->getTotalCoriolisAcceleration(vc) + Amob; // 18 flops
}

Real RigidBodyNode::calcKineticEnergy(
    const SBTreePositionCache& pc,
    const SBTreeVelocityCache& vc) const 
{
    const Real ret = dot(getV_GB(vc) , getMk_G(pc)*getV_GB(vc));
    return ret/2;
}



//==============================================================================
//                     CALC JOINT INDEPENDENT DYNAMICS VEL
//==============================================================================
// Calculate mass and velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and centrifugal 
// forces due to coriolis acceleration.
// This routine expects that spatial velocities, spatial inertias, and
// articulated body inertias are already available.
// Must be called base to tip.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    SBDynamicsCache&                        dc) const
{
    if (nodeNum == 0) { // ground, just in case
        updMobilizerCentrifugalForces(dc)    = SpatialVec(Vec3(0), Vec3(0));
        updTotalCentrifugalForces(dc)       = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    // 72 flops
    updMobilizerCentrifugalForces(dc) =
        getP(abc) * getMobilizerCoriolisAcceleration(vc) + getGyroscopicForce(vc);

    // 72 flops
    updTotalCentrifugalForces(dc) = 
        getP(abc) * getTotalCoriolisAcceleration(vc) + getGyroscopicForce(vc);

}


//
// Given only position-related quantities from the State 
//      Mk_G  (this body's spatial inertia matrix, exp. in Ground)
//      Phi   (composite body child-to-parent shift matrix)
// calculate the inverse dynamics quantity
//      R     (composite body inertia)
// This must be called tip-to-base (inward).
//
// Note that this method does not depend on the mobilities of
// the joint so can be implemented here rather than in RigidBodyNodeSpec.
// Ground should override this implementation though.
//
void
RigidBodyNode::calcCompositeBodyInertiasInward(
    const SBTreePositionCache&  pc,
    Array_<SpatialInertia>& allR) const
{
    SpatialInertia& R = toB(allR);
    R = getMk_G(pc);
    for (unsigned i=0; i<children.size(); ++i) {
        const SpatialInertia& RChild  = children[i]->fromB(allR);
        const PhiMatrix&  phiChild    = children[i]->getPhi(pc);
        R += RChild.shift(-phiChild.l()); // ~80 flops
    }
}
