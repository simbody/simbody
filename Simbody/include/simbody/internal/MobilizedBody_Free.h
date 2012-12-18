#ifndef SimTK_SIMBODY_MOBILIZED_BODY_FREE_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_FREE_H_

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
Declares the MobilizedBody::Free class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Unrestricted motion for a rigid body (six mobilities). Orientation
/// is modeled the same as for the Orientation mobilizer, that is, using
/// quaternions to avoid singularities. A modeling option exists to 
/// have the joint modeled with a 1-2-3 body fixed Euler sequence like
/// a Gimbal mobilizer. Translational generalized coordinates are
/// x,y,z translations along the F (inboard) axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Free : public MobilizedBody {
public:
    explicit Free(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Free(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Free(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame,
         Direction=Forward);

    Free& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Free& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Free& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Free& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Free& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Leaves rotation unchanged.
    Free& setDefaultTranslation(const Vec3&);

    // Leaves translation unchanged. The internal representation is a quaternion
    // so we guarantee that the stored value is numerically identical to the
    // supplied one.
    Free& setDefaultQuaternion(const Quaternion&);

    // Leaves translation unchanged. The Rotation matrix will be converted to
    // a quaternion for storage.
    Free& setDefaultRotation(const Rotation&);
    // Sets both translation and rotation. The Rotation part of the Transform 
    // will be converted to a quaternion for storage.
    Free& setDefaultTransform(const Transform&);

    // These return references to the stored default values.
    const Vec3& getDefaultTranslation() const;
    const Quaternion& getDefaultQuaternion() const;

    // These next two are derived from the stored values.
    Rotation getDefaultRotation() const {
        return Rotation(getDefaultQuaternion());
    }
    Transform getDefaultTransform() const {
        return Transform(Rotation(getDefaultQuaternion()), getDefaultTranslation());
    }

    // Generic default state Topology methods.

    // Returns (Vec4,Vec3) where the Vec4 is a normalized quaternion.
    const Vec7& getDefaultQ() const;

    // Interprets the supplied q as (Vec4,Vec3) where the Vec4 is a possibly
    // unnormalized quaternion. The quaternion will be normalized before it is
    // stored here, so you may not get back exactly the value supplied here if
    // you call getDefaultQ().
    Free& setDefaultQ(const Vec7& q);

    // Note that there is no guarantee that the quaternion part of the returned Q is normalized.
    const Vec7& getQ(const State&) const;
    const Vec7& getQDot(const State&) const;
    const Vec7& getQDotDot(const State&) const;

    const Vec6& getU(const State&) const;
    const Vec6& getUDot(const State&) const;

    // The Q's in the state are set exactly as supplied without normalization.
    void setQ(State&, const Vec7&) const;
    void setU(State&, const Vec6&) const;

    const Vec7& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec6& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec7& updMyPartQ(const State&, Vector& qlike) const;
    Vec6& updMyPartU(const State&, Vector& ulike) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Free, FreeImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_FREE_H_



