#ifndef SimTK_SIMBODY_CONSTRAINT_WELD_H_
#define SimTK_SIMBODY_CONSTRAINT_WELD_H_

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
Declares the Constraint::Weld class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                  WELD (COINCIDENT FRAMES) CONSTRAINT
//==============================================================================

/** Six constraint equations. This constraint enforces coincidence between
a frame on one body and a frame on another body. This is a combination
of a ConstantOrientation constraint and a Ball constraint. The first three
equations correspond to the perpendicularity constraints associated with
the orientation constraint, the last three equations are the
coincident point conditions.

The constraint is enforced by an internal (non-working) force applied at the
spatial location of the frame origin on body 2, on material points of each body
that are coincident with that spatial location. Note that this is somewhat
asymmetric when the Weld is not properly assembled -- it acts as though the
contact occurs at the origin of the frame on body 2, *not* at the origin of the
frame on body 1. The orientation constraints on the other hand are symmetric,
they are three "constant angle" constraints enforcing perpendicularity between
body2's x,y,z axes with body1's y,z,x axes respectively, via an internal
(non-working) torque vector applied equal and opposite on each body.

TODO: Although the frame origins can be brought together by the Ball constraint,
the perpendicularity conditions can be satisfied with antiparallel axes in
addition to the parallel ones we want. Therefore the assembly conditions must
include additional (redundant) constraints requiring parallel axes.
**/
class SimTK_SIMBODY_EXPORT Constraint::Weld : public Constraint {
public:
    // no default constructor

/// Make the body frame of one body coincident with the body frame
/// of the other body.
Weld(MobilizedBody& body1, MobilizedBody& body2);

/// Make a particular frame attached to one body coincident with
/// a particular frame attached to the other body. The frames are
/// specified by giving the transform X_BF which expresses the
/// position and orientation of frame F relative to the body frame B.
Weld(MobilizedBody& body1, const Transform& frame1,
        MobilizedBody& body2, const Transform& frame2);

/** Default constructor creates an empty handle. **/
Weld() {}

    // Control over generated decorative geometry.

/// This is used only for visualization. Set r <= 0 to disable
/// default frame drawing. Default axis length is r=1. This is a
/// topology-stage variable, not changeable later.
Weld& setAxisDisplayLength(Real r);

/// Report the length being used for display of the frames being
/// connected by this Weld. If this returns 0 then no geometry is
/// being generated for the frames.
Real getAxisDisplayLength() const;

    // Defaults for Instance variables.

/// Explicitly set the default value for the frame on body1 which
/// is to be made coincident with a frame on body2. Note that this is
/// topology-stage value so requires non-const access to the Constraint.
Weld& setDefaultFrameOnBody1(const Transform&);

/// Retrieve the default transform for the frame on body 1.
const Transform& getDefaultFrameOnBody1() const;

/// Explicitly set the default value for the frame on body2 which
/// is to be made coincident with a frame on body1. Note that this is
/// topology-stage value so requires non-const access to the Constraint.
Weld& setDefaultFrameOnBody2(const Transform&);

/// Retrieve the default transform for the frame on body 2.
const Transform& getDefaultFrameOnBody2() const;


    // Stage::Topology

/// Report the MobilizedBodyIndex of body 1 for this Weld constraint.
MobilizedBodyIndex getBody1MobilizedBodyIndex() const;

/// Report the MobilizedBodyIndex of body 2 for this Weld constraint.
MobilizedBodyIndex getBody2MobilizedBodyIndex() const;


    // Stage::Instance
const Transform& getFrameOnBody1(const State&) const;
const Transform& getFrameOnBody2(const State&) const;

    // Stage::Position, Velocity, Acceleration
Vec6 getPositionErrors(const State&) const;
Vec6 getVelocityErrors(const State&) const;

    // Stage::Acceleration
Vec6 getAccelerationErrors(const State&) const;
Vec6 getMultipliers(const State&) const;

    // Forces are reported expressed in the body frame of the indicated body.
const SpatialVec& getWeldReactionOnBody1(const State&) const;
const SpatialVec& getWeldReactionOnBody2(const State&) const;

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Weld, WeldImpl, Constraint);
/** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_WELD_H_



