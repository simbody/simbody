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
 * Contributors: Derived from IVM code written by Charles Schwieters          *
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



//==============================================================================
//                    CALC JOINT INDEPENDENT KINEMATICS POS
//==============================================================================
// Calc posCM, mass, Mk
//      phi, inertia
// Should be calc'd from base to tip.
// We depend on transforms X_PB and X_GB being available.
// Cost is 90 flops.
void RigidBodyNode::calcJointIndependentKinematicsPos(
    SBTreePositionCache&    pc) const
{
    assert(nodeNum != 0); // Don't call this for Ground.

    // Re-express parent-to-child shift vector (Bo-Po) into the ground frame.
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



//==============================================================================
//                   CALC JOINT INDEPENDENT KINEMATICS VEL
//==============================================================================
// Calculate velocity-related quantities: spatial velocity (V_GB),
// and coriolis acceleration. This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G and
// velocity-dependent acceleration remainder term VD_PB_G.
// Cost is about 100 flops.
void
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBTreePositionCache& pc,
    SBTreeVelocityCache&       vc) const
{
    assert(nodeNum != 0); // Don't call this for Ground.

    const SpatialVec& V_GP   = parent->getV_GB(vc); // parent P's velocity in G
    const SpatialVec& V_PB_G = getV_PB_G(vc); // child B's vel in P, exp. in G
    const PhiMatrixTranspose PhiT = ~getPhi(pc); // shift outwards

    // calc spatial velocity of B's origin Bo in G (angular,linear)
    const SpatialVec V_GB = PhiT*V_GP + V_PB_G; // 18 flops
    updV_GB(vc) = V_GB;

    const Vec3& w_GB = V_GB[0]; // for convenience
    const Vec3& v_GB = V_GB[1];

    // Calculate gyroscopic moment and force b (48 flops). Although this is
    // really a dynamic quantity (requires spatial inertia and is itself
    // a force), we calculate this here rather than in stage Dynamics because
    // it is needed in inverse dynamics but does not require articulated body
    // inertias to be calculated. This could be deferred until needed but
    // probably isn't worth the trouble.
    const SpatialVec b  = getMass() *                       // 6 flops
        SpatialVec(w_GB % (getUnitInertia_OB_G(pc)*w_GB),   // moment (24 flops)
                   w_GB % (w_GB % getCB_G(pc)));            // force  (18 flops)
    updGyroscopicForce(vc) = b;

    // Parent velocity.
    const Vec3& w_GP = V_GP[0]; // for convenience
    const Vec3& v_GP = V_GP[1];

    // Calculate this mobilizer's incremental contribution to coriolis
    // acceleration a, and this body's total coriolis acceleration A (it
    // parent's coriolis acceleration plus the incremental contribution).
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
    //   A = ~PhiDot * V_GP + HDot * u
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
    const SpatialVec  A(VD_PB_G[0],
                        VD_PB_G[1] + w_GP % (v_GB-v_GP)); // 15 flops
    updMobilizerCoriolisAcceleration(vc) = A;

    // Next, the total coriolis acceleration a (normally just called "coriolis
    // acceleration"!) of body B is the total coriolis acceleration of its
    // parent shifted outward, plus B's local contribution A that we just
    // calculated (18 flops).
    const SpatialVec a = PhiT * parent->getTotalCoriolisAcceleration(vc) + A;
    updTotalCoriolisAcceleration(vc) = a;

    // Finally, calculate the total of the rotational velocity-dependent forces
    // acting on this body (45 flops).
    updTotalCentrifugalForces(vc) =  getMk_G(pc) * a + b;
}



//==============================================================================
//                          CALC KINETIC ENERGY
//==============================================================================
// Cost is 57 flops.
Real RigidBodyNode::calcKineticEnergy(
    const SBTreePositionCache& pc,
    const SBTreeVelocityCache& vc) const
{
    const Real ret = dot(getV_GB(vc) , getMk_G(pc)*getV_GB(vc));
    return ret/2;
}



//==============================================================================
//                 REALIZE ARTICULATED BODY VELOCITY CACHE
//==============================================================================
// Calculate velocity-related quantities that also depend on articulated body
// inertias. This routine expects that coriolis accelerations, gyroscopic
// forces, and articulated body inertias are already available.
// As written, the calling order doesn't matter, but Ground's entries must
// have been precalculated so don't call this on the Ground body.
// Cost is 72 flops.
void
RigidBodyNode::realizeArticulatedBodyVelocityCache
   (const SBTreePositionCache&              pc,
    const SBTreeVelocityCache&              vc,
    const SBArticulatedBodyInertiaCache&    abc,
    SBArticulatedBodyVelocityCache&         abvc) const
{
    assert(nodeNum != 0); // Don't call this for Ground.

    // 72 flops (P*a + b)
    updArticulatedBodyCentrifugalForces(abvc) =
        getP(abc)*getMobilizerCoriolisAcceleration(vc) + getGyroscopicForce(vc);
}



//==============================================================================
//                      CALC COMPOSITE BODY INERTIAS
//==============================================================================
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
// Cost is about 80 flops.
void
RigidBodyNode::calcCompositeBodyInertiasInward(
    const SBTreePositionCache&  pc,
    Array_<SpatialInertia,MobilizedBodyIndex>& allR) const
{
    SpatialInertia& R = toB(allR);
    R = getMk_G(pc);
    for (unsigned i=0; i<children.size(); ++i) {
        const SpatialInertia& RChild  = children[i]->fromB(allR);
        const PhiMatrix&  phiChild    = children[i]->getPhi(pc);
        R += RChild.shift(-phiChild.l()); // ~80 flops
    }
}
