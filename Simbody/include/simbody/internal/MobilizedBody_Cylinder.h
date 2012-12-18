#ifndef SimTK_SIMBODY_MOBILIZED_BODY_CYLINDER_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_CYLINDER_H_

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
 * Contributors: Paul Mitiguy, Peter Eastman                                  *
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
Declares the MobilizedBody::Cylinder class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Two mobilities -- rotation and translation along the common z axis
/// of the inboard and outboard mobilizer frames.
class SimTK_SIMBODY_EXPORT MobilizedBody::Cylinder : public MobilizedBody {
public:
    explicit Cylinder(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Cylinder(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Cylinder(MobilizedBody& parent, const Transform& inbFrame,
             const Body&,           const Transform& outbFrame, Direction=Forward);

    Cylinder& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Cylinder& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Cylinder& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Cylinder& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Cylinder& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    /** @cond **/ // Don't let doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Cylinder, CylinderImpl, 
                                             MobilizedBody);
    /** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_CYLINDER_H_



