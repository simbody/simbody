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
    SBPositionCache&   pc) const
{
    // Re-express parent-to-child shift vector (OB-OP) into the ground frame.
    const Vec3 p_PB_G = getX_GP(pc).R() * getX_PB(pc).p();

    // The Phi matrix conveniently performs child-to-parent (inward) shifting
    // on spatial quantities (forces); its transpose does parent-to-child
    // (outward) shifting for velocities and accelerations.
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
    const Mat33 offDiag = getMass()*crossMat(getCB_G(pc));
    updMk(pc) = SpatialMat( getInertia_OB_G(pc).toMat33() ,     offDiag ,
                                   -offDiag             , getMass()*Mat33(1) );
}

// Calculate velocity-related quantities: spatial velocity (V_GB). This must be
// called base to tip: depends on parent's spatial velocity, and
// the just-calculated cross-joint spatial velocity V_PB_G.
void 
RigidBodyNode::calcJointIndependentKinematicsVel(
    const SBPositionCache& pc,
    SBVelocityCache&       mc) const
{
    updV_GB(mc) = ~getPhi(pc)*parent->getV_GB(mc) + getV_PB_G(mc);
}

Real RigidBodyNode::calcKineticEnergy(
    const SBPositionCache& pc,
    const SBVelocityCache& mc) const 
{
    const Real ret = dot(getV_GB(mc) , getMk(pc)*getV_GB(mc));
    return 0.5*ret;
}

// Calculate velocity-related quantities that are needed for building
// our dynamics operators, namely the gyroscopic force and coriolis acceleration.
// This routine expects that all spatial velocities & spatial inertias are
// already available.
// Must be called base to tip.
void 
RigidBodyNode::calcJointIndependentDynamicsVel(
    const SBPositionCache& pc,
    const SBVelocityCache& mc,
    SBDynamicsCache&       dc) const
{
    if (nodeNum == 0) { // ground, just in case
        updGyroscopicForce(dc)           = SpatialVec(Vec3(0), Vec3(0));
        updCoriolisAcceleration(dc)      = SpatialVec(Vec3(0), Vec3(0));
        updTotalCoriolisAcceleration(dc) = SpatialVec(Vec3(0), Vec3(0));
        updCentrifugalForces(dc)         = SpatialVec(Vec3(0), Vec3(0));
        updTotalCentrifugalForces(dc)    = SpatialVec(Vec3(0), Vec3(0));
        return;
    }

    const Vec3& w_GB = getV_GB(mc)[0];  // spatial angular velocity
    const Vec3& v_GB = getV_GB(mc)[1];  // spatial linear velocity (of B origin in G)

    updGyroscopicForce(dc) = 
        SpatialVec(    w_GB % (getInertia_OB_G(pc)*w_GB),     // gyroscopic moment
                    getMass()*(w_GB % (w_GB % getCB_G(pc)))); // gyroscopic force

    // Parent velocity.
    const Vec3& w_GP = parent->getV_GB(mc)[0];
    const Vec3& v_GP = parent->getV_GB(mc)[1];

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

    updCoriolisAcceleration(dc) =
        SpatialVec( Vec3(0), w_GP % (v_GB-v_GP) ) + getVD_PB_G(dc);

    updTotalCoriolisAcceleration(dc) =
        ~getPhi(pc) * parent->getTotalCoriolisAcceleration(dc)
        + getCoriolisAcceleration(dc); // just calculated above

    updCentrifugalForces(dc) =
        getP(dc) * getCoriolisAcceleration(dc) + getGyroscopicForce(dc);

    updTotalCentrifugalForces(dc) = 
        getP(dc) * getTotalCoriolisAcceleration(dc) + getGyroscopicForce(dc);

}


