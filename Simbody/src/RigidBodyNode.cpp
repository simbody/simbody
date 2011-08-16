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
// gyroscopic forces, coriolis acceleration. This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G and
// velocity-dependent acceleration remainder term VD_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBTreePositionCache& pc,
    SBTreeVelocityCache&       vc) const
{
    if (nodeNum == 0) { // ground, just in case
        updV_GB(vc)                      = SpatialVec(Vec3(0), Vec3(0));
        updGyroscopicForce(vc)           = SpatialVec(Vec3(0), Vec3(0));
        updCoriolisAcceleration(vc)      = SpatialVec(Vec3(0), Vec3(0));
        updTotalCoriolisAcceleration(vc) = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    // 18 flops
    updV_GB(vc) = ~getPhi(pc)*parent->getV_GB(vc) + getV_PB_G(vc);

    const Vec3& w_GB = getV_GB(vc)[0];  // spatial angular velocity
    const Vec3& v_GB = getV_GB(vc)[1];  // spatial linear velocity (of B origin in G)

    updGyroscopicForce(vc) = 
        getMass()*SpatialVec(w_GB % (getUnitInertia_OB_G(pc)*w_GB), // gyroscopic moment (24 flops)
                             (w_GB % (w_GB % getCB_G(pc)))); // gyroscopic force  (21 flops)

    // Parent velocity.
    const Vec3& w_GP = parent->getV_GB(vc)[0];
    const Vec3& v_GP = parent->getV_GB(vc)[1];

    // Calc a: coriolis acceleration.
    // The coriolis acceleration "a" is a 
    // "remainder" term in the spatial acceleration, depending only on velocities,
    // but involving time derivatives of the Phi and H matrices. 
    // CAUTION: our definition of H is transposed from Jain's and Schwieters'.
    //
    // Specifically,
    //   a = ~PhiDot * V_GP + HDot * u
    // As correctly calculated in Schwieters' paper, Eq [16], the first term above
    // simplifies to SpatialVec( 0, w_GP % (v_GB-v_GP) ). However, Schwieters' second
    // term in [16] is correct only if H is constant in P, in which case the derivative
    // just accounts for the rotation of P in G. In general H is not constant in P,
    // so we don't try to calculate the derivative here but assume that HDot*u has
    // already been calculated for us and stored in VD_PB_G. (That is,
    // V_PB_G = H*u, VD_PB_G = HDot*u.)

    updCoriolisAcceleration(vc) =
        SpatialVec( Vec3(0), w_GP % (v_GB-v_GP) ) + getVD_PB_G(vc); // 18 flops

    // 18 flops
    updTotalCoriolisAcceleration(vc) =
        ~getPhi(pc) * parent->getTotalCoriolisAcceleration(vc)
        + getCoriolisAcceleration(vc); // just calculated above
}

Real RigidBodyNode::calcKineticEnergy(
    const SBTreePositionCache& pc,
    const SBTreeVelocityCache& vc) const 
{
    const Real ret = dot(getV_GB(vc) , getMk_G(pc)*getV_GB(vc));
    return 0.5*ret;
}

// Calculate velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and coriolis acceleration.
// This routine expects that all spatial velocities & spatial inertias are
// already available.
// Must be called base to tip.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(
    const SBTreePositionCache&              pc,
    const SBArticulatedBodyInertiaCache&    abc,
    const SBTreeVelocityCache&              vc,
    SBDynamicsCache&                        dc) const
{
    if (nodeNum == 0) { // ground, just in case
        updCentrifugalForces(dc)         = SpatialVec(Vec3(0), Vec3(0));
        updTotalCentrifugalForces(dc)    = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    // 72 flops
    updCentrifugalForces(dc) =
        getP(abc) * getCoriolisAcceleration(vc) + getGyroscopicForce(vc);

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
