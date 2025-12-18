#ifndef SimTK_SIMBODY_MOBILIZED_BODY_CANTILEVER_FREE_BEAM_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_CANTILEVER_FREE_BEAM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-25 Stanford University and the Authors.        *
 * Authors: Nicholas Bianco                                                   *
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
Declares the MobilizedBody::CantileverFreeBeam class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** Although this mobilizer has only three mobilities (all rotational), the
mobilized frame M translates relative to the fixed frame F as a specified
function of this mobilizers rotational generalized coordinates q₀, q₁, q₂. The
specified function is based on the deflection shape of a cantilever-free beam,
whose free end is subjected to a transverse point load. Although this shape is
inspired by a cantilever-free beam, this mobilizer does not actually model a
beam, in that there is no elastic restoring forces/torques when q₀ or q₁ or
q₂ ≠ 0. In the undeformed state (when q₀ = q₁ = q₂ = 0), the position from Fo
(frame F's origin) to Mo (frame M's origin) is L Fz.

\section cf_beam_eqs Cantilever-free beam equations

To define this mobilizer, we first need to derive the equations for the
deflection of a cantilever-free beam under a transverse point load. The general
differential equation for finite elastic deformations of a beam subject to a
bending moment, M, is given by EI d²y⁄dx² = M (equation 8.1-5 in [1]),
where E is the modulus of elasticity, I is the area moment of
inertia, x is the position along the axis of the undeflected beam, and
y is the vertical deflection.

A cantilever-free beam (with length L) is rigidly fixed at one end (i.e., y = 0
and dy⁄dx = 0 at x = 0) and free at the other end (i.e., d²y⁄dx² = 0 and
d³y⁄dx³ = 0 at x = L). When subjected to a transverse point load, P, at the end
of the beam, the bending moment is given by M = P(L-x). Substituting this
loading condition into the differential equation and applying the boundary
conditions, we can solve for the vertical deflection at the end of the beam:

y(L) = PL³ ⁄ (3EI)

Similarly, we can derive the angle of deflection, θ, at the end of the beam:

θ(L) = PL² ⁄ (2EI)

\subsection beam_deflection Beam deflection equation

Using the two previously derived expressions, we can define a kinematic
relationship between the deflection and angle of the beam's endpoint:

y(L) = (2⁄3) θ L

Note that this equation is only a function of the deflection angle and the
beam length, and not E, I, or P, so we can use it to define
the kinematic coupling between the X-Y translations and the X-Y rotations of the
mobilizer defined below (see \ref mobilizer).

\subsection beam_displacement Beam displacement equation

The shortening of the beam can be estimated by equation 8.1-14 from "Roark's
Formulas for Stress and Strain" [1]:

ΔL = −0.5 ∫₀ᴸ (dy⁄dx)² dx

After integrating, applying the boundary conditions, and substituting for
θ, we get the following relationship relating the shortening of the
beam to the angle of deflection:

ΔL(θ) = −(4⁄15) θ² L

This equation will drive the kinematic coupling between the Z translation and
the X-Y rotations of the mobilizer defined below (see \ref mobilizer).

\section mobilizer Mobilizer definition

The generalized coordinates q are the same as for a Gimbal mobilizer, that is,
an X-Y-Z body-fixed Euler sequence. The three generalized speeds u for
this mobilizer are also the same as for a Gimbal mobilizer: the time derivatives
of the generalized coordinates, that is, uᵢ = q̇ᵢ (i.e., uᵢ = qdotᵢ for i = 0,
1, 2).

The first two rotations, q₀ and q₁, induce translations according to the beam
deflection formula for a cantilever beam under a transverse point load applied
at the end of the beam (see \ref beam_deflection). The change in the beam's
endpoint position in the Fz direction is governed by apparent shortening of a
beam due to bending (see \ref beam_displacement). Therefore, the components of
the position vector from frame F's origin to frame M's origin, p₀, p₁, and p₂,
are given by:

p₀ = 2⁄3 q₁ L

p₁ = −2⁄3 q₀ L

p₂ = L − 4⁄15 (q₀² + q₁²) L

The third rotation, q₂, induces a rotation about the Mz axis, which is always
tangent to the beam at the beam's endpoint.

Note that while the endpoint shortens in the Fz direction, the total length of
the beam actually lengthens slightly (about 3% for a rotation of π/4 radians
about X or Y), since the beam formula for an end-loaded cantilever beam assumes
finite deflections and is derived from a simplified form of the beam curvature
term in the Euler-Bernoulli beam equation. However, since the primary utility of
this mobilizer is to provide a lightweight way for modeling flexible structures
(e.g., the bending of the spinal column in a human or animal skeleton), and not
necessarily to accurately model large beam deflections, the beam lengthening can
be accounted for with appropriate modeling adjustments.

While this mobilizer provides arbitrary orientation, the Euler angle
derivatives are singular when q₁ (the middle rotation) is near ±π/2 radians.
That means you should not attempt to do any dynamics in that configuration.

\par References:

[1] Young, W. C., and R. G. Budynas. (2002). Roark's Formulas for Stress and
    Strain (7th Edition). McGraw-Hill.

**/
class SimTK_SIMBODY_EXPORT MobilizedBody::CantileverFreeBeam :
        public MobilizedBody {
public:
    /** Default constructor provides an empty handle that can be assigned to
    reference any %MobilizedBody::CantileverFreeBeam. **/
    CantileverFreeBeam() {}

    /** Create a %CantileverFreeBeam mobilizer between an existing parent
    (inboard) body P and a new child (outboard) body B created by copying the
    given \a bodyInfo into a privately-owned Body within the constructed
    %MobilizedBody object. Specify the mobilizer frames F fixed to parent P and
    M fixed to child B. The origin of the beam is placed on the mobilizer's
    inboard frame F, and the endpoint is placed on the mobilizer's outboard
    frame M. In the default configuration, the Fz and Mz and the parent frame
    origin is located at (0,0,L) in F, where L is the length of the beam.
    @see MobilizedBody for a diagram and explanation of terminology. **/
    CantileverFreeBeam(MobilizedBody& parent, const Transform& X_PF,
        const Body& bodyInfo, const Transform& X_BM, const Real& length,
        Direction=Forward);

    /** Abbreviated constructor you can use if the mobilizer frames are
    coincident with the parent and child body frames. **/
    CantileverFreeBeam(MobilizedBody& parent,  const Body& bodyInfo,
        const Real& length, Direction=Forward);

    /** This constructor assumes you'll set the beam length later; meanwhile it
    uses a default length. **/
    CantileverFreeBeam(MobilizedBody& parent, const Transform& X_PF,
        const Body& bodyInfo, const Transform& X_BM, Direction=Forward);

    /** This constructor assumes you'll set the beam length later; meanwhile it
    uses a default length. The parent and child body frames are used as the
    mobilizer frames. **/
    CantileverFreeBeam(MobilizedBody& parent, const Body& bodyInfo,
        Direction=Forward);

    /** Modify the default length of the beam. This is usually set on
    construction. **/
    CantileverFreeBeam& setDefaultLength(const Real& length);

    /** Get the default length of the beam as specified during construction or
    via setDefaultLength(). **/
    const Real& getDefaultLength() const;

    /** Provide a default orientation for this mobilizer if you don't want to
    start with the identity rotation (that is, alignment of the F and M
    frames). This is the orientation the mobilizer will have in the default
    state for the containing System. The supplied Rotation will be converted
    to a body fixed x-y-z Euler sequence and used as the three generalized
    coordinates q. **/
    CantileverFreeBeam& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToBodyFixedXYZ());
    }

    /** Get the value for the default orientation of this mobilizer; unless
    changed by setDefaultRotation() it will be the identity rotation. **/
    Rotation getDefaultRotation() const {
        Rotation R_FM;
        R_FM.setRotationToBodyFixedXYZ(getDefaultQ());
        return R_FM;
    }

    CantileverFreeBeam& addBodyDecoration(const Transform& X_BD,
            const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }

    CantileverFreeBeam& addOutboardDecoration(const Transform& X_MD,
            const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }

    CantileverFreeBeam& addInboardDecoration (const Transform& X_FD,
            const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    CantileverFreeBeam& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    CantileverFreeBeam& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Vec3& updDefaultQ();
    CantileverFreeBeam& setDefaultQ(const Vec3& q) {
        updDefaultQ() = q;
        return *this;
    }

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(CantileverFreeBeam,
        CantileverFreeBeamImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_CANTILEVER_FREE_BEAM_H_
