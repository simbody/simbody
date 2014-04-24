#ifndef SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_
#define SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_

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

/**@file
Private implementation of Constraint::Rod. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Constraint.h"

#include "ConstraintImpl.h"

class SimbodyMatterSubsystemRep;

namespace SimTK {

class SimbodyMatterSubsystem;
class SimbodyMatterSubtree;
class MobilizedBody;

//==============================================================================
//                                  ROD IMPL
//==============================================================================
// TODO: should use distance rather than distance^2 to improve scaling.
class Constraint::RodImpl : public ConstraintImpl {
public:
RodImpl() 
    : ConstraintImpl(1,0,0), defaultPoint1(0), defaultPoint2(0), defaultRodLength(1),
    pointRadius(-1) // this means "use default point radius"
{ 
    // Rod constructor sets all the data members here directly
}
RodImpl* clone() const { return new RodImpl(*this); }

// Draw some end points and a rubber band line.
void calcDecorativeGeometryAndAppendVirtual
    (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const;

void setPointDisplayRadius(Real r) {
    // r == 0 means don't display point, r < 0 means use default which is some fraction of rod length
    invalidateTopologyCache();
    pointRadius= r > 0 ? r : 0;
}
Real getPointDisplayRadius() const {return pointRadius;}

// Implementation of virtuals required for holonomic constraints.

// TODO: this is badly scaled; consider using length instead of length^2
// perr = (p^2 - d^2)/2
void calcPositionErrorsVirtual      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   X_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const
{
    assert(X_AB.size()==2 && constrainedQ.size()==0 && perr.size()==1);
    const Vec3 p1 = findStationLocation(X_AB, B1, defaultPoint1); // meas from & expr in ancestor
    const Vec3 p2 = findStationLocation(X_AB, B2, defaultPoint2);
    const Vec3 p = p2 - p1;
    //TODO: save p in state

    perr[0] = (dot(p, p) - square(defaultRodLength)) / 2;
}

// pverr = d/dt perr = pdot*p = v*p, where v=v2-v1 is relative velocity
void calcPositionDotErrorsVirtual      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const 
{
    assert(V_AB.size()==2 && constrainedQDot.size()==0 && pverr.size()==1);
    //TODO: should be able to get p from State
    const Vec3 p1 = findStationLocationFromState(s, B1, defaultPoint1); // meas from & expr in ancestor
    const Vec3 p2 = findStationLocationFromState(s, B2, defaultPoint2);
    const Vec3 p = p2 - p1;

    const Vec3 v1 = findStationVelocity(s, V_AB, B1, defaultPoint1); // meas & expr in ancestor
    const Vec3 v2 = findStationVelocity(s, V_AB, B2, defaultPoint2);
    const Vec3 v = v2 - v1;
    pverr[0] = dot(v, p);
}

// paerr = d/dt verr = vdot*p + v*pdot =a*p+v*v, where a=a2-a1 is relative 
// acceleration
void calcPositionDotDotErrorsVirtual      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const
{
    assert(A_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size()==1);
    //TODO: should be able to get p and v from State
    const Vec3 p1 = findStationLocationFromState(s, B1, defaultPoint1); // meas from & expr in ancestor
    const Vec3 p2 = findStationLocationFromState(s, B2, defaultPoint2);
    const Vec3 p = p2 - p1;
    const Vec3 v1 = findStationVelocityFromState(s, B1, defaultPoint1); // meas & expr in ancestor
    const Vec3 v2 = findStationVelocityFromState(s, B2, defaultPoint2);
    const Vec3 v = v2 - v1;

    const Vec3 a1 = findStationAcceleration(s, A_AB, B1, defaultPoint1); // meas & expr in ancestor
    const Vec3 a2 = findStationAcceleration(s, A_AB, B2, defaultPoint2);
    const Vec3 a = a2 - a1;

    paerr[0] = dot(a, p) + dot(v, v);
}

// Write this routine by inspection of the pdot routine, looking for terms 
// involving velocity. On point2 we see v2*p, on point1 we see -v1*p, so forces
// are m*p and -m*p, respectively.
void addInPositionConstraintForcesVirtual
   (const State&                                    s,      // Stage::Position
    const Array_<Real>&                             multipliers, // mp of these
    Array_<SpatialVec,ConstrainedBodyIndex>&        bodyForcesInA,
    Array_<Real,      ConstrainedQIndex>&           qForces) 
    const
{
    assert(multipliers.size()==1 && bodyForcesInA.size()==2 
           && qForces.size()==0);
    const Real lambda = multipliers[0];
    //TODO: should be able to get p from State
    const Vec3 p1 = findStationLocationFromState(s, B1, defaultPoint1); // meas from & expr in ancestor
    const Vec3 p2 = findStationLocationFromState(s, B2, defaultPoint2);
    const Vec3 p = p2 - p1;

    const Vec3 f2 = lambda * p;

    // The forces on either point have the same line of action because they are
    // aligned with the vector between the points. Applying the forces to any 
    // point along the line would have the same effect (e.g., same point in 
    // space on both bodies) so this is the same as an equal and opposite force
    // applied to the same point and this constraint will do no work even if 
    // the position or velocity constraints are not satisfied.
    addInStationForce(s, B2, defaultPoint2,  f2, bodyForcesInA);
    addInStationForce(s, B1, defaultPoint1, -f2, bodyForcesInA);
}

SimTK_DOWNCAST(RodImpl, ConstraintImpl);
//------------------------------------------------------------------------------
                                    private:
friend class Constraint::Rod;

ConstrainedBodyIndex    B1; // must be 0 and 1
ConstrainedBodyIndex    B2;
Vec3                    defaultPoint1; // on body 1, exp. in B1 frame
Vec3                    defaultPoint2; // on body 2, exp. in B2 frame
Real                    defaultRodLength;

// This is just for visualization
Real                    pointRadius;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_ROD_IMPL_H_



