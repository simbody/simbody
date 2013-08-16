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
 * Portions copyright (c) 2007-13 Stanford University and the Authors.        *
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

/** Zero mobilities. This degenerate "mobilizer" serves only to weld together
the M frame of a body to the F frame on its parent. It has no generalized
coordinates or speeds. Note that there is no "reverse" weld, because "reverse" 
for a mobilizer refers to how the q's and u's are defined and there are none.

You can use this (im)mobilizer to create a composite rigid body from simpler
rigid bodies. Although the effect is as though the bodies were combined, they
are still tracked separately so per-body information remains available. You can
also get the reaction force at the weld in the usual manner. **/
class SimTK_SIMBODY_EXPORT MobilizedBody::Weld : public MobilizedBody {
public:
    /** Default constructor provides an empty handle that can be assigned to
    reference any %MobilizedBody::Weld. **/
    Weld() {};

    /** Create a %Weld mobilizer between an existing parent (inboard) body P 
    and a new child (outboard) body B created by copying the given \a bodyInfo 
    into a privately-owned Body within the constructed %MobilizedBody object. 
    Specify the mobilizer frames F fixed to parent P and M fixed to child B. 
    @see MobilizedBody for a diagram and explanation of terminology. **/
    Weld(MobilizedBody& parent, const Transform& X_PF,
         const Body& bodyInfo,  const Transform& X_BM);

    /** Abbreviated constructor you can use if the mobilizer frames are 
    coincident with the parent and child body frames. **/
    Weld(MobilizedBody& parent, const Body& bodyInfo);

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



