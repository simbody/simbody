/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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
 * This file contains the code used to build the various constraint types.
 */
#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintRep.h"
#include "SimbodyMatterSubsystemRep.h"
#include "ConstraintNode.h"

#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setprecision;

///////////////////////////////////////////////
// Implementation of ConstraintNode methods. //
///////////////////////////////////////////////

/*
 * TODO: How to specify a constraint equation:
 *
 * Required info:
 *
 *    Dependencies:
 *      - constraint level: position, velocity, acceleration
 *      - has time dependence?
 *    A list of bodies
 *      - this is just those bodies to which constraint
 *        forces & torques are *directly* applied to 
 *        enforce the constraint (can be spatial force
 *        or mobility force)
 *
 * All constraints:
 *    Acceleration error (given udot)
 *    Constraint forces (given lambda)
 * Position or velocity level:
 *    Velocity error (given u)
 * Position level only:
 *    Position error (given q)
 */


    //////////////////////////////
    // CONSTRAINT::ROD::RODREP //
    /////////////////////////////

ConstraintNode* Constraint::Rod::RodRep::createConstraintNode() const
{
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B1);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(B2);
    return new ConstantDistanceConstraintNode(
        sbdyrep.getRigidBodyNode(mobilizedBody1), defaultPoint1,
        sbdyrep.getRigidBodyNode(mobilizedBody2), defaultPoint2,
        defaultRodLength);
}



    ////////////////////////////////////////////////////
    // CONSTRAINT::POINT IN PLANE::POINT IN PLANE REP //
    ////////////////////////////////////////////////////

ConstraintNode* Constraint::PointInPlane::PointInPlaneRep::createConstraintNode() const
{
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex planeMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(planeBody);
    const MobilizedBodyIndex followerMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(followerBody);

    return new PointInPlaneConstraintNode(
        sbdyrep.getRigidBodyNode(planeMobilizedBody), defaultPlaneNormal, defaultPlaneHeight,
        sbdyrep.getRigidBodyNode(followerMobilizedBody), defaultFollowerPoint);
}


    //////////////////////////////////////////////////
    // CONSTRAINT::POINT ON LINE::POINT ON LINE REP //
    //////////////////////////////////////////////////

ConstraintNode* Constraint::PointOnLine::PointOnLineRep::createConstraintNode() const
{
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex lineMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(lineBody);
    const MobilizedBodyIndex followerMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(followerBody);

    return new PointOnLineConstraintNode(
        sbdyrep.getRigidBodyNode(lineMobilizedBody), defaultLineDirection, defaultPointOnLine,
        sbdyrep.getRigidBodyNode(followerMobilizedBody), defaultFollowerPoint);
}

    ////////////////////////////////////////////////////
    // CONSTRAINT::CONSTANT ANGLE::CONSTANT ANGLE REP //
    ////////////////////////////////////////////////////

ConstraintNode* Constraint::ConstantAngle::ConstantAngleRep::createConstraintNode() const {
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex baseMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(B);
    const MobilizedBodyIndex followerMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(F);
    return new ConstantAngleConstraintNode(
        sbdyrep.getRigidBodyNode(baseMobilizedBody), defaultAxisB,
        sbdyrep.getRigidBodyNode(followerMobilizedBody), defaultAxisF, defaultAngle);
}

    ////////////////////////////////////////////////////////////////
    // CONSTRAINT::CONSTANT ORIENTATION::CONSTANT ORIENTATION REP //
    ////////////////////////////////////////////////////////////////

ConstraintNode* Constraint::ConstantOrientation::ConstantOrientationRep::createConstraintNode() const {
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex baseMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(B);
    const MobilizedBodyIndex followerMobilizedBody = getMobilizedBodyIndexOfConstrainedBody(F);
    return new ConstantOrientationConstraintNode(
        sbdyrep.getRigidBodyNode(baseMobilizedBody), defaultRB,
        sbdyrep.getRigidBodyNode(followerMobilizedBody), defaultRF);
}

    ///////////////////////////////
    // CONSTRAINT::BALL::BALLREP //
    ///////////////////////////////

ConstraintNode* Constraint::Ball::BallRep::createConstraintNode() const {
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B1);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(B2);
    return new CoincidentStationsConstraintNode(
        sbdyrep.getRigidBodyNode(mobilizedBody1), defaultPoint1,
        sbdyrep.getRigidBodyNode(mobilizedBody2), defaultPoint2);
}


    ///////////////////////////////
    // CONSTRAINT::WELD::WELDREP //
    ///////////////////////////////

ConstraintNode* Constraint::Weld::WeldRep::createConstraintNode() const {
    assert(isInSubsystem());
    const SimbodyMatterSubsystemRep& sbdyrep = getMyMatterSubsystemRep();
    const MobilizedBodyIndex mobilizedBody1 = getMobilizedBodyIndexOfConstrainedBody(B);
    const MobilizedBodyIndex mobilizedBody2 = getMobilizedBodyIndexOfConstrainedBody(F);
    return new WeldConstraintNode(
        sbdyrep.getRigidBodyNode(mobilizedBody1), defaultFrameB,
        sbdyrep.getRigidBodyNode(mobilizedBody2), defaultFrameF);
}

