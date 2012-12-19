#ifndef SimTK_SIMBODY_MOBILIZED_BODY_GIMBAL_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_GIMBAL_H_

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
Declares the MobilizedBody::Gimbal class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** Three mobilities -- unrestricted orientation modeled as a 1-2-3 body-fixed 
Euler angle sequence, though with a singularity when the middle angle is +/- 90 
degrees. Use MobilizedBody::Ball for a similar joint that is free of 
singularities.

This mobilizer does not provide any translation, so the parent's F frame and
child's M frame origins will always be coincident. (MobilizedBody::Bushing is
defined like Gimbal but adds translation.) When q0=q1=q2=0, the F and M 
frames are aligned. Then the generalized coordinates q should be interpreted
as a series of rotation angles in radians: q0 is a rotation about the x axis, 
then q1 is a rotation about the now-rotated y axis, and then q2 is a rotation 
about the now twice-rotated z axis. The generalized speeds u for the %Gimbal
mobilizer are the time derivatives of
the generalized coordinates, that is, u=qdot. Note that this is different than
the choice of generalized speeds used for MobilizedBody::Ball, where the
angular velocity vector is used as the three generalized speeds. (Euler angle 
derivatives are \e not the same as angular velocity, except when q=0.)

While this mobilizer provides arbitrary orientation, the Euler angle 
derivatives are singular when q1 (the middle rotation) is near +/- Pi/2. That
means you should not attempt to do any dynamics in that configuration; if you
can't be sure the motion will remain away from that region, you should use
a MobilizedBody::Ball instead, which uses quaternions to ensure 
singularity-free rotation.

Here is another way to describe the meaning of a %Gimbal mobilizer:
A %Gimbal between frames F and M is exactly equivalent to a series 
of three pin joints using two massless intermediate bodies B1 and B2. The
first pin joint rotates about the common Fx and B1x axes. The second rotates
about the common B1y and B2y axes. The last one rotates about the common B2z
and Mz axes. Other than performance (which is somewhat better if you use the
%Gimbal), the pins-and-intermediate bodies system generates exactly the same
generalized coordinates and generalized speeds as the %Gimbal. 
 
@see MobilizedBody::Ball, MobilizedBody::Gimbal **/
class SimTK_SIMBODY_EXPORT MobilizedBody::Gimbal : public MobilizedBody {
public:
    /** Create a %Gimbal mobilizer between an existing parent (inboard) body P 
    and a new child (outboard) body B created by copying the given \a bodyInfo 
    into a privately-owned Body within the constructed %MobilizedBody::%Gimbal 
    object. Specify the mobilizer frames F fixed to parent P and M fixed to 
    child B. **/
    Gimbal(MobilizedBody& parent, const Transform& X_PF,
           const Body& bodyInfo,  const Transform& X_BM, Direction=Forward);

    /** Abbreviated constructor you can use if the mobilizer frames are 
    coincident with the parent and child body frames. **/
    Gimbal(MobilizedBody& parent, const Body& bodyInfo, Direction=Forward);

    /** Change the default inboard ("fixed") frame F on the parent body of this
    mobilizer. The supplied frame is given as a Transform from the parent
    body's frame P and replaces the F frame specified in the constructor. This
    is a topological change, meaning that you must call realizeTopology() again 
    if you call this method. **/
    Gimbal& setDefaultInboardFrame(const Transform& X_PF) 
    {   (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this; }

    /** Change the default outboard ("moving") frame M on this mobilized body
    (that is, on this mobilizer's child body). The supplied frame is given as
    a Transform from the child body's frame B and replaces the M frame specified
    in the constructor. This is a topological change, meaning that you must 
    call realizeTopology() again if you use this method. **/
    Gimbal& setDefaultOutboardFrame(const Transform& X_BM) 
    {   (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this; }

    /** Override the default orientation for this mobilizer. Normally the
    default is q=0. Here you can provide a Rotation giving the default 
    orientation of the mobilizer's M frame in the parent body's F frame. That
    orientation will be converted to a body fixed 1-2-3 Euler sequence and
    used to set the default values for the generalized coordinates q that will
    appear in a State that has just been created by a realizeTopology() 
    call. This is a topological change, meaning that you must call 
    realizeTopology() again if you use this method. **/
    Gimbal& setDefaultRotation(const Rotation& R_FM) 
    {   return setDefaultQ(R_FM.convertRotationToBodyFixedXYZ()); }

    /** Return the default orientation for this mobilizer as a Rotation 
    matrix. This is the initial orientation of the mobilizer's M frame in
    the parent body's F frame that you will see in a State that has just been
    created via a realizeTopology() call. This is a topological change, 
    meaning that you must call realizeTopology() again if you use this
    method. **/ 
    Rotation getDefaultRotation() const {
        const Vec3& q = getDefaultQ();
        return Rotation(BodyRotationSequence,
            q[0], XAxis, q[1], YAxis, q[2], ZAxis);
    }

    /** Add DecorativeGeometry to the privately-owned Body contained in this
    %MobilizedBody, positioned relative to the body's frame. **/
    Gimbal& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this; }

    /** Add DecorativeGeometry to the privately-owned Body contained in this
    %MobilizedBody, positioned relative to the mobilizer's M frame on that 
    body. **/
    Gimbal& addOutboardDecoration(const Transform& X_MD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this; }

    /** Add DecorativeGeometry to the parent (inboard) body of this mobilizer, 
    positioned relative to the mobilizer's F frame on the parent body. **/
    Gimbal& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this; }


    /** This affects the radius of the ball that will be drawn as default
    visualization geometry for a %Gimbal mobilizer. It has no effect whatsoever
    on simulation or other computations. **/
    Gimbal& setDefaultRadius(Real r);
    /** Get the current value of the default visualization ball radius. **/
    Real getDefaultRadius() const;

    /** Obtain the default value for the generalized coordinates q that has
    been set for this mobilizer. The value is returned as a 3-vector but note
    that the returned quantity is \e not a vector; it is a body fixed 1-2-3
    Euler angle sequence in radians. The default is q=0 unless overridden. 
    @see getDefaultRotation() **/
    const Vec3& getDefaultQ() const; // X,Y,Z body-fixed Euler angles
    /** Set the default value for the generalized coordinates q for this 
    mobilizer. The value is given as a 3-vector but note
    that the supplied quantity is \e not a vector; it is a body fixed 1-2-3
    Euler angle sequence in radians. The default is q=0 unless overridden.  
    @see getDefaultRotation() **/
    Gimbal& setDefaultQ(const Vec3& q);

    /** Obtain the current value for this mobilizer's generalized coordinates 
    q from the given State, which must be realized through Stage::Model (these
    are state variables). **/ 
    const Vec3& getQ(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized coordinate 
    time derivatives qdot from the given State, which must be realized through
    Stage::Velocity. These are Euler angle time derivatives and, for the %Gimbal
    mobilizer, have the same numerical values as the generalized speeds u. **/ 
    const Vec3& getQDot(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized coordinate 
    second time derivatives qdotdot from the given State, which must be 
    realized through Stage::Acceleration. These are Euler angle second time
    derivatives and, for the %Gimbal mobilizer, have the same numerical value
    as the generalized accelerations udot. **/ 
    const Vec3& getQDotDot(const State& state) const;

    /** Obtain the current value for this mobilizer's generalized speeds u
    from the given State, which must be realized through Stage::Model (these are
    state variables). For the %Gimbal joint, these are the time derivatives of  
    q, that is, qdot=u. **/ 
    const Vec3& getU(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized accelerations
    udot (d/dt u) from the given State, which must be realized through
    Stage::Acceleration. For the %Gimbal joint, these are the time derivatives 
    of qdot, that is, qdotdot=udot. **/ 
    const Vec3& getUDot(const State& state) const;

    /** Set new values for this mobilizer's generalized coordinates q in the
    given \a state. This invalidates Stage::Position and above. The new
    value is given as a 3-vector, but is interpreted as a body fixed 1-2-3
    Euler angle sequence in radians. **/
    void setQ(State& state, const Vec3& q) const;
    /** Set new values for this mobilizer's generalized speeds u in the
    given \a state. This invalidates Stage::Velocity and above. The new
    value is given as a 3-vector, but is interpreted as the time derivatives
    of a body fixed 1-2-3 Euler angle sequence in radians/time. **/
    void setU(State& state, const Vec3& u) const;

    /** @name               Advanced/Obscure
    Most users won't use these methods. **/
    /**@{**/
    /** Create a disembodied %Gimbal mobilizer that is not part of any System. **/
    explicit Gimbal(Direction=Forward);

    /** Given a vector in the System's generalized coordinate basis, extract
    the three q's belonging to this mobilizer. **/
    const Vec3& getMyPartQ(const State& state, const Vector& qlike) const;
    /** Given a vector in the System's generalized speed basis, extract
    the three u's belonging to this mobilizer. **/
    const Vec3& getMyPartU(const State& state, const Vector& ulike) const;
   
    /** Given a writable vector in the System's generalized coordinate basis, 
    extract a writable reference to the three q's in that vector that belong to
    this mobilizer. **/
    Vec3& updMyPartQ(const State& state, Vector& qlike) const;
    /** Given a writable vector in the System's generalized speed basis, 
    extract a writable reference to the three u's in that vector that belong to
    this mobilizer. **/
    Vec3& updMyPartU(const State& state, Vector& ulike) const;
    /**@}**/

    /** @cond **/ // Hide from Doxygen.
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Gimbal, GimbalImpl, MobilizedBody);
    /** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_GIMBAL_H_



