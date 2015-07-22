#ifndef SimTK_SIMBODY_MOBILIZED_BODY_SCREW_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_SCREW_H_

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
Declares the MobilizedBody::Screw class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** One mobility -- coordinated rotation and translation along the
common z axis of the inboard and outboard mobilizer frames. A
"pitch" is specified relating the two. The generalized coordinate
q is the rotation angle in radians, the translation is always
pitch*q. **/
class SimTK_SIMBODY_EXPORT MobilizedBody::Screw : public MobilizedBody {
public:
    /** Default constructor provides an empty handle that can be assigned to
    reference any %MobilizedBody::Screw. **/
    Screw() {}

    /** Create a %Screw mobilizer between an existing parent (inboard) body P
    and a new child (outboard) body B created by copying the given \a bodyInfo
    into a privately-owned Body within the constructed %MobilizedBody object.
    Specify the mobilizer frames F fixed to parent P and M fixed to child B.
    @see MobilizedBody for a diagram and explanation of terminology. **/
    Screw(MobilizedBody& parent, const Transform& X_PF,
          const Body& bodyInfo,  const Transform& X_BM,
          Real pitch, Direction=Forward);

    /** Abbreviated constructor you can use if the mobilizer frames are
    coincident with the parent and child body frames. **/
    Screw(MobilizedBody& parent, const Body& bodyInfo, Real pitch,
          Direction=Forward);


    Screw& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Screw& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Screw& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Screw& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Screw& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    Screw& setDefaultPitch(Real pitch);
    Real   getDefaultPitch() const;

    Screw& setDefaultQ(Real);
    Real   getDefaultQ() const;

    Real getQ(const State&) const;
    Real getQDot(const State&) const;
    Real getQDotDot(const State&) const;
    Real getU(const State&) const;
    Real getUDot(const State&) const;

    void setQ(State&, Real) const;
    void setU(State&, Real) const;

    Real getMyPartQ(const State&, const Vector& qlike) const;
    Real getMyPartU(const State&, const Vector& ulike) const;

    Real& updMyPartQ(const State&, Vector& qlike) const;
    Real& updMyPartU(const State&, Vector& ulike) const;

    /** @cond **/ // Don't let doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Screw, ScrewImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_SCREW_H_



