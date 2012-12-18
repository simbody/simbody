#ifndef SimTK_SIMBODY_MOBILIZED_BODY_SPHERICALCOORDS_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_SPHERICALCOORDS_H_

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
 * Contributors: Christopher Bruns                                            *
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
Declares the MobilizedBody::SphericalCoords class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** Three mobilities -- body fixed 3-2 (z-y) rotation followed by translation
along body z or body x. Interpreted as spherical coordinates the first rotation
is the azimuth angle, the second is the zenith, and the translation is
the radius. We permit a simple mapping from generalized coordinates to
(azimuth, zenith, radius):
<pre>
    azimuth = s0*q0 + az0   (about Fz==Mz)
    zenith  = s1*q1 + ze0   (about My)
    radius  = s2*q2         (along Mz or Mx; Mz is default)
</pre>
where s0,s1,s2 are signs (1 or -1) and az0 and ze0 are offset angles. The
F and M frames are coincident when azimuth==zenith==radius==0. But note
that with non-zero offsets the F and M frames will not be aligned in
the reference configuration where q0==q1==q2==0. The F and M origins
will always be coincident when q2==0, however.

With this you can define a "geographical" coordinate system where Mx is the 
Greenwich line, a is latitude and z longitude (with north positive):
<pre>
     v  = Mx
     s0 =  1, az0 = 0
     s1 = -1, ze0 = 0
     s2 =  1
</pre>
If you want the translation direction to be in Mz (the default) but would like 
q1=0 to mean the equatorial position rather than the (possibly singular) north
pole which should be -90, define
<pre>
     v  = Mz
     s0 = 1, az0 = 0
     s1 = 1, ze0 = Pi/2
     s2 = 1
</pre>
One common convention for atomic (torsion,bend,stretch) uses the default 
spherical coordinate system but the final stretch is along the -z direction. 
For that, take all defaults but set s2=-1.

This mobilizer can be used to give unrestricted 3-d motion to inertialess 
particles (as with a Cartesian mobilizer but parameterized torsion,bend,stretch
instead of x,y,z) but in this case you must watch for two possible 
singularities: (1) radius==0, and (2) zenith==n*Pi (or equivalently 
q1==n*Pi-s1*ze0). If your operating range steers clear of those singularities, 
you're fine. **/
class SimTK_SIMBODY_EXPORT MobilizedBody::SphericalCoords : public MobilizedBody {
public:
    explicit SphericalCoords(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    SphericalCoords(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames. This gives you a pure spherical
    /// coordinate system in which q0=azimuth about Fz(==Mz), q1=zenith about My, 
    /// and q2=radius along Mz.
    SphericalCoords(MobilizedBody& parent, const Transform& inbFrame,
                    const Body&,           const Transform& outbFrame,
                    Direction=Forward);

    /// Use this constructor to specify the general case described above.
    SphericalCoords(MobilizedBody& parent,      const Transform& inbFrame,
                    const Body&,                const Transform& outbFrame,
                    Real azimuthOffset,         bool azimuthNegated,
                    Real zenithOffset,          bool zenithNegated,
                    CoordinateAxis radialAxis,  bool radialNegated,
                    Direction=Forward);

    SphericalCoords& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    SphericalCoords& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    SphericalCoords& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    SphericalCoords& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    SphericalCoords& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Friendly, mobilizer-specific access to coordinates and speeds.
    SphericalCoords& setDefaultAngles(const Vec2& a) {
        Vec3 q = getDefaultQ(); q.updSubVec<2>(0) = a; setDefaultQ(q);
        return *this;
    }
    SphericalCoords& setDefaultRadius(Real r) {
        Vec3 q = getDefaultQ(); q[2] = r; setDefaultQ(q);
        return *this;
    }
    SphericalCoords& setRadialAxis(CoordinateAxis);
    SphericalCoords& setNegateAzimuth(bool);
    SphericalCoords& setNegateZenith(bool);
    SphericalCoords& setNegateRadial(bool);

    const Vec2&    getDefaultAngles()      const {return getDefaultQ().getSubVec<2>(0);}
    Real           getDefaultTranslation() const {return getDefaultQ()[2];}
    
    CoordinateAxis getRadialAxis()    const;
    bool           isAzimuthNegated() const;
    bool           isZenithNegated()  const;
    bool           isRadialNegated()  const;

    void setAngles(State& s, const Vec2& a) {setOneQ(s,0,a[0]); setOneQ(s,1,a[1]);}
    void setRadius(State& s, Real        r) {setOneQ(s,2,r);}

    const Vec2& getAngles(const State& s) const {return getQ(s).getSubVec<2>(0);}
    Real        getRadius(const State& s) const {return getQ(s)[2];}

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    SphericalCoords& setDefaultQ(const Vec3& q);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(SphericalCoords, 
                                             SphericalCoordsImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_SPHERICALCOORDS_H_



