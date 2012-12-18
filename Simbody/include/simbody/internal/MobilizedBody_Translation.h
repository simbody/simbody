#ifndef SimTK_SIMBODY_MOBILIZED_BODY_TRANSLATION_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_TRANSLATION_H_

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
Declares the MobilizedBody::Translation class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Three translational mobilities. The generalized coordinates are
/// x,y,z translations along the parent (inboard) F frame axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Translation : public MobilizedBody {
public:
    explicit Translation(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Translation(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Translation(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame, Direction=Forward);

    Translation& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Translation& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Translation& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Translation& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Translation& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for setting this mobilizer's generalized coordinates,
    // which together constitute the vector from the F frame's origin to the M
    // frame's origin, expressed in F.

    // Set the topological default values for the initial q's.
    Translation& setDefaultTranslation(const Vec3& p_FM) {
        return setDefaultQ(p_FM);
    }

    // Get the topological default values for the initial q's.
    const Vec3& getDefaultTranslation() const {
        return getDefaultQ();
    }

    // Set the current value of q's in the given State. Note that this is
    // the *cross-mobilizer* translation, not location in the Ground frame.
    void setMobilizerTranslation(State& s, const Vec3& p_FM) const {
        setQ(s,p_FM);
    }

    // Get the current value of the q's for this mobilizer from the given State.
    const Vec3& getMobilizerTranslation(const State& s) const {
        return getQ(s);
    }


    // Set the current value of u's in the given State. Note that this is
    // the *cross-mobilizer* velocity v_FM, not velocity in the Ground frame.
    void setMobilizerVelocity(State& s, const Vec3& v_FM) const {
        setU(s,v_FM);
    }

    // Get the current value of the u's for this mobilizer from the given State.
    const Vec3& getMobilizerVelocity(const State& s) const {
        return getU(s);
    }

    // Get the value of the udot's for this mobilizer from the given State.
    const Vec3& getMobilizerAcceleration(const State& s) const {
        return getUDot(s);
    }

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Translation& setDefaultQ(const Vec3& q);

    const Vec3& getQ(const State&) const;
    const Vec3& getQDot(const State&) const;
    const Vec3& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec3&) const;
    void setU(State&, const Vec3&) const;

    const Vec3& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec3& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Translation, TranslationImpl, 
                                             MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_TRANSLATION_H_



