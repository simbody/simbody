#ifndef SimTK_SIMBODY_MOBILIZED_BODY_ELLIPSOID_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_ELLIPSOID_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
 * Authors: Ajay Seth, Michael Sherman                                        *
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
Declares the MobilizedBody::Ellipsoid class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Three mobilities -- coordinated rotation and translation along the
/// surface of an ellipsoid fixed to the parent (inboard) body.
/// The generalized coordinates are the same as for a Ball (Orientation)
/// joint, that is, a quaternion or 1-2-3 Euler sequence.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ellipsoid : public MobilizedBody {
public:
    /// The ellipsoid is placed on the mobilizer's inboard frame F, with
    /// half-axis dimensions along F's x,y,z respectively.
    explicit Ellipsoid(Direction=Forward); // not very useful until radii are set, but has some defaults

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ellipsoid(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame,
              Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames, and give the radii at
    /// the same time.
    Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame,
              const Vec3& radii,
              Direction=Forward);

    Ellipsoid& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Ellipsoid& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Ellipsoid& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Ellipsoid& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Ellipsoid& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for the generalized coordinate.
    Ellipsoid& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToQuaternion());
    }
    Rotation getDefaultRotation() const {return Rotation(getDefaultQ());}

    Ellipsoid& setDefaultRadii(const Vec3& r);
    const Vec3& getDefaultRadii() const;

    // Generic default state Topology methods.
    const Quaternion& getDefaultQ() const;
    Quaternion& updDefaultQ();
    Ellipsoid& setDefaultQ(const Quaternion& q) {updDefaultQ()=q; return *this;}

    const Vec4& getQ(const State&) const;
    const Vec4& getQDot(const State&) const;
    const Vec4& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec4&) const;
    void setU(State&, const Vec3&) const;

    const Vec4& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec4& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ellipsoid, EllipsoidImpl, 
                                             MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_ELLIPSOID_H_



