#ifndef SimTK_SIMBODY_MOBILIZED_BODY_WELD_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_WELD_H_

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
Declares the MobilizedBody::Weld class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Zero mobilities. This degenerate "mobilizer" serves only to weld together
/// the M frame of a body to the F frame on its parent.
class SimTK_SIMBODY_EXPORT MobilizedBody::Weld : public MobilizedBody {
public:
    /// Note: there is no "reverse" weld, because "reverse" refers to
    /// how the q's and u's are defined and there are none.
    Weld();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Weld(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Weld(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Weld& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Weld& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Weld& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Weld& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Weld& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Weld, WeldImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_WELD_H_



