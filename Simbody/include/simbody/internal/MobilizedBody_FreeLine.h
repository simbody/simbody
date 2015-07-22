#ifndef SimTK_SIMBODY_MOBILIZED_BODY_FREELINE_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_FREELINE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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
Declares the MobilizedBody::FreeLine class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** Five mobilities, representing unrestricted motion for a body which is
inertialess along its own z axis. The rotational generalized coordinates are the
same as for the LineOrientation mobilizer. The translational coordinates are
the same as in a Free mobilizer, or a Cartesian (Translation) mobilizer.

LineOrientation and FreeLine are special "ball" and "free" mobilizers designed
to allow arbitrary orientations for "linear" bodies, such as a CO2 molecule
consisting only of point masses arranged along a straight line. Such bodies have
no inertia about the line and cause singularities in the equations of motion if
attached to Ball (Spherical) or Free mobilizers. Instead, use the
LineOrientation and LineFree mobilizers, making sure that the inertialess
direction is along the outboard body's z axis (that is, Mz). These mobilizers
introduce only two rotational mobilities (generalized speeds u), being incapable
of representing non-zero angular velocity of M in F about Mz. The generalized
speeds are in fact the wx and wy components of w_FM_M, that is, the x and y
components of the angular velocity of M in F <em>expressed in M</em>. However,
at least three generalized coordinates (q's) are required to represent the
orientation. By default we use four quaternions for unconditional stability.
Alternatively, you can request a 1-2-3 body fixed Euler angle sequence (that is,
about x, then new y, then new z) which will suffer a singularity when the y
rotation is 90 degrees since that aligns the first rotation axis (x) with the
last (z) which is the inertialess direction.

@see MobilizedBody::LineOrientation, MobilizedBody::Free **/
class SimTK_SIMBODY_EXPORT MobilizedBody::FreeLine : public MobilizedBody {
public:
    /** Default constructor provides an empty handle that can be assigned to
    reference any %MobilizedBody::Ball. **/
    FreeLine() {}

    /** Create a %FreeLine mobilizer between an existing parent (inboard) body P
    and a new child (outboard) body B created by copying the given \a bodyInfo
    into a privately-owned Body within the constructed %MobilizedBody object.
    Specify the mobilizer frames F fixed to parent P and M fixed to child B.
    @see MobilizedBody for a diagram and explanation of terminology. **/
    FreeLine(MobilizedBody& parent, const Transform& X_PF,
             const Body& bodyInfo,  const Transform& X_BM, Direction=Forward);

    /** Abbreviated constructor you can use if the mobilizer frames are
    coincident with the parent and child body frames. **/
    FreeLine(MobilizedBody& parent, const Body& bodyInfo, Direction=Forward);

    FreeLine& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    FreeLine& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    FreeLine& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    FreeLine& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    FreeLine& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(FreeLine, FreeLineImpl,
                                             MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_FREELINE_H_



