#ifndef SimTK_SIMBODY_MOBILIZED_BODY_PLANAR_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_PLANAR_H_

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
Declares the MobilizedBody::Planar class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/// Three mobilities -- z rotation and x,y translation. The generalized
/// coordinates are rotation about the shared z axis of the F and M
/// frame, translation along the F frame's x axis, and translation along
/// its y axis, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::Planar : public MobilizedBody {
public:
    explicit Planar(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Planar(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Planar(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame,
           Direction=Forward);

    Planar& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Planar& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Planar& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Planar& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Planar& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Friendly, mobilizer-specific access to coordinates and speeds.
    Planar& setDefaultAngle(Real a) {
        Vec3 q = getDefaultQ(); q[0] = a; setDefaultQ(q);
        return *this;
    }
    Planar& setDefaultTranslation(const Vec2& r) {
        Vec3 q = getDefaultQ(); q.updSubVec<2>(1) = r; setDefaultQ(q);
        return *this;
    }

    Real getDefaultAngle() const {return getDefaultQ()[0];}
    const Vec2& getDefaultTranslation() const {return getDefaultQ().getSubVec<2>(1);}

    void setAngle      (State& s, Real        a) {setOneQ(s,0,a);}
    void setTranslation(State& s, const Vec2& r) {setOneQ(s,1,r[0]); setOneQ(s,2,r[1]);}

    Real getAngle(const State& s) const {return getQ(s)[0];}
    const Vec2& getTranslation(const State& s) const {return getQ(s).getSubVec<2>(1);}

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Planar& setDefaultQ(const Vec3& q);

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

    /** @cond **/ // Don't let doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Planar, PlanarImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_PLANAR_H_



