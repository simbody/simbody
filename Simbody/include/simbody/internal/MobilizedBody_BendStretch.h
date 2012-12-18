#ifndef SimTK_SIMBODY_MOBILIZED_BODY_BENDSTRETCH_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_BENDSTRETCH_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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
Declares the MobilizedBody::BendStretch class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Two mobilities: The z axis of the parent's F frame is 
/// used for rotation (and that is always aligned with the M frame z axis).
/// The x axis of the *M* (outboard) frame is then used for translation;
/// that is, first we rotate around z, which moves M's x with respect to F's x. Then
/// we slide along the rotated x axis. The two generalized coordinates are the
/// rotation and the translation, in that order.
/// This can also be viewed a a 2D polar coordinate mobilizer since the coordinates
/// are (theta, r) about perpendicular axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::BendStretch : public MobilizedBody {
public:
    explicit BendStretch(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    BendStretch(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame, Direction=Forward);

    BendStretch& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    BendStretch& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    BendStretch& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    BendStretch& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    BendStretch& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    /** @cond **/ // Don't let doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(BendStretch, BendStretchImpl, 
                                             MobilizedBody);
    /** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_BENDSTRETCH_H_



