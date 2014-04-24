#ifndef SimTK_SIMBODY_CONSTRAINT_POINT_IN_PLANE_H_
#define SimTK_SIMBODY_CONSTRAINT_POINT_IN_PLANE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/** @file
Declares the Constraint::PointInPlane class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                      CONSTRAINT::  POINT IN PLANE
//==============================================================================

/** One constraint equation. This constraint enforces that a point fixed to
one body (the "follower body") must travel in a plane fixed on another body
(the "plane body"). The constraint is enforced by an internal (non-working)
scalar force acting at the spatial location of the follower point, directed 
along the plane normal, and equal and opposite on the two bodies.

The assembly condition is the same as the run-time constraint: the point
has to be moved into the plane.
**/
class SimTK_SIMBODY_EXPORT Constraint::PointInPlane : public Constraint  {
public:
// no default constructor
PointInPlane(MobilizedBody& planeBody_B, 
             const UnitVec3& defaultPlaneNormal_B, Real defaultHeight,
             MobilizedBody& followerBody_F, const Vec3& defaultFollowerPoint_F);
    
/** Default constructor creates an empty handle. **/
PointInPlane() {}

// These affect only generated decorative geometry for visualization;
// the plane is really infinite in extent with zero depth and the
// point is really of zero radius.
PointInPlane& setPlaneDisplayHalfWidth(Real);
PointInPlane& setPointDisplayRadius(Real);
Real getPlaneDisplayHalfWidth() const;
Real getPointDisplayRadius() const;

// Defaults for Instance variables.
PointInPlane& setDefaultPlaneNormal(const UnitVec3&);
PointInPlane& setDefaultPlaneHeight(Real);
PointInPlane& setDefaultFollowerPoint(const Vec3&);

// Stage::Topology
MobilizedBodyIndex getPlaneMobilizedBodyIndex() const;
MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

const UnitVec3& getDefaultPlaneNormal() const;
Real            getDefaultPlaneHeight() const;
const Vec3&     getDefaultFollowerPoint() const;

// Stage::Instance
const UnitVec3& getPlaneNormal(const State&) const;
Real            getPlaneHeight(const State&) const;
const Vec3&     getFollowerPoint(const State&) const;

// Stage::Position, Velocity
Real getPositionError(const State&) const;
Real getVelocityError(const State&) const;

// Stage::Acceleration
Real getAccelerationError(const State&) const;
Real getMultiplier(const State&) const;
Real getForceOnFollowerPoint(const State&) const; // in normal direction

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (PointInPlane, PointInPlaneImpl, Constraint);
/** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_POINT_IN_PLANE_H_



