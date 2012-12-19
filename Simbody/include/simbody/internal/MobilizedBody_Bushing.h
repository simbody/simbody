#ifndef SimTK_SIMBODY_MOBILIZED_BODY_BUSHING_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_BUSHING_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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
Declares the MobilizedBody::Bushing class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** Six mobilities -- arbitrary relative motion modeled as x-y-z translation
followed by an x-y-z body-fixed Euler angle sequence, though with a 
singularity when the middle angle is +/- 90 degrees. Use MobilizedBody::Free 
for a similar six-dof joint that is free of singularities.

The six generalized coordinates q are defined q={qx,qy,qz,px,py,pz} where the
first three coordinates are the rotation angles and the last three are the
Cartesian translations. Note that we follow the usual Simbody convention of
having the rotational coordinates precede the translational ones; that does
not imply any particular order of operations. In fact, as a sequence this
mobilizer should be thought of as first translating the M frame by a vector
p_FM=[px,py,pz] (expressed in F), then reorienting M about its new origin 
location using the rotation angles. In the reference configuration where q=0, 
the F and M frames are aligned and their origins coincident (Mo=Fo). The
convention used here exactly matches the convention used by the force element 
Force::LinearBushing which infers an identical set of coordinates to express the
relative pose of its two connected frames (also called F and M).

The first three generalized coordinates q should be interpreted as a series of 
rotation angles in radians: q[0]=qx is a rotation about the Mx (=Fx) axis, then 
q[1]=qy is a rotation about the now-rotated My axis, and then q[2]=qz is a 
rotation about the now twice-rotated Mz axis. These rotations do not affect the 
location of M's origin point Mo. The generalized speeds u for the %Bushing
mobilizer are just the time derivatives of the generalized coordinates, that is, 
u=qdot={qxdot,qydot,qzdot,vx,vy,vz}, where qxdot=d/dt qx, etc. and 
v_FM=[vx,vy,vz]=d/dt p_FM is the velocity of Mo measured and expressed in frame
F. Note that this is different than the choice of generalized speeds used for
MobilizedBody::Free, where the angular velocity vector w_FM is used as the three
rotational generalized speeds instead of Euler angle derivatives. (Euler angle 
derivatives are \e not the same as angular velocity, except when q=0.)

While this mobilizer provides arbitrary orientation, the Euler angle derivatives 
are singular when q[1] (qy, the middle rotation) is near +/- Pi/2. That means 
you should not attempt to do any dynamics in that configuration; if you can't 
be sure the motion will remain away from that region, you should use a 
MobilizedBody::Free instead, which uses quaternions to ensure singularity-free 
rotation.

Here is another way to describe the meaning of a %Bushing mobilizer: except
for the ordering of the generalized coordinates, a %Bushing between frames F 
on parent body P and frame M on child body B is exactly equivalent to a 
translation (Cartesian) mobilizer between F and a massless body B0, followed by 
a series of three pin joints using two massless intermediate bodies B1 and B2. 
The first pin joint rotates about the common B0x(=Fx) and B1x axes. The second 
rotates about the common B1y and B2y axes. The last one rotates about the common
B2z and Mz axes. Other than performance (which is somewhat better if you use the
%Bushing since there are fewer bodies), the translation+pins+intermediate 
bodies system generates exactly the same mobility as the %Bushing, but with 
coordinates q={px,py,pz,qx,qy,qz} (flipped from the %Bushing coordinate 
ordering).
 
@see Force::LinearBushing, MobilizedBody::Free, MobilizedBody::Gimbal **/
class SimTK_SIMBODY_EXPORT MobilizedBody::Bushing : public MobilizedBody {
public:
    /** Create a %Bushing mobilizer between an existing parent (inboard) body P 
    and a new child (outboard) body B created by copying the given \a bodyInfo 
    into a privately-owned Body within the constructed %MobilizedBody::%Bushing 
    object. Specify the mobilizer frames F fixed to parent P and M fixed to 
    child B. **/
    Bushing(MobilizedBody& parent, const Transform& X_PF,
           const Body& bodyInfo,  const Transform& X_BM, Direction=Forward);

    /** Abbreviated constructor you can use if the mobilizer frames are 
    coincident with the parent and child body frames. **/
    Bushing(MobilizedBody& parent, const Body& bodyInfo, Direction=Forward);

    /** Change the default inboard ("fixed") frame F on the parent body of this
    mobilizer. The supplied frame is given as a Transform from the parent
    body's frame P and replaces the F frame specified in the constructor. This
    is a topological change, meaning that you must call realizeTopology() again 
    if you call this method. **/
    Bushing& setDefaultInboardFrame(const Transform& X_PF) 
    {   (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this; }

    /** Change the default outboard ("moving") frame M on this mobilized body
    (that is, on this mobilizer's child body). The supplied frame is given as
    a Transform from the child body's frame B and replaces the M frame specified
    in the constructor. This is a topological change, meaning that you must 
    call realizeTopology() again if you use this method. **/
    Bushing& setDefaultOutboardFrame(const Transform& X_BM) 
    {   (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this; }

    /** Override the default pose for this mobilizer. Normally the
    default is q=0, meaning that the F and M frames are aligned. Here you can 
    provide a Transform giving the default pose of the mobilizer's M frame in 
    the parent body's F frame. The given orientation will be converted to a 
    body fixed x-y-z Euler sequence; the given position vector will be used
    directly. Those six values will be used to set the default values for the 
    generalized coordinates q that will appear in a State that has just been 
    created by a realizeTopology() call. This is a topological change, meaning 
    that you must call realizeTopology() again if you use this method. **/
    Bushing& setDefaultTransform(const Transform& X_FM) 
    {   Vec6 q;
        q.updSubVec<3>(0) = X_FM.R().convertRotationToBodyFixedXYZ();
        q.updSubVec<3>(3) = X_FM.p();      
        return setDefaultQ(q); }

    /** Return the default pose for this mobilizer as a Transform. See 
    setDefaultTransform() for more information. **/ 
    Transform getDefaultTransform() const {
        const Vec6& q = getDefaultQ();
        return Transform(Rotation(BodyRotationSequence,
                                  q[0], XAxis, q[1], YAxis, q[2], ZAxis),
                         q.getSubVec<3>(3));
    }

    /** Add DecorativeGeometry to the privately-owned Body contained in this
    %MobilizedBody, positioned relative to the body's frame. **/
    Bushing& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this; }

    /** Add DecorativeGeometry to the privately-owned Body contained in this
    %MobilizedBody, positioned relative to the mobilizer's M frame on that 
    body. **/
    Bushing& addOutboardDecoration(const Transform& X_MD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this; }

    /** Add DecorativeGeometry to the parent (inboard) body of this mobilizer, 
    positioned relative to the mobilizer's F frame on the parent body. **/
    Bushing& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) 
    {   (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this; }

    /** Obtain the default value for the generalized coordinates q that has
    been set for this mobilizer. The value is returned as a Vec6 but note
    that the rotational part (first three elements) of the returned quantity 
    is \e not a vector; it is a body fixed x-y-z Euler angle sequence in 
    radians. The last three elements are a position vector p_FM. The default 
    is q=0 unless overridden. 
    @see getDefaultTransform() **/
    const Vec6& getDefaultQ() const;
    /** Set the default value for the generalized coordinates q for this 
    mobilizer. The value is given as a Vec6 but note
    that the supplied quantity is \e not a pair of vectors; it is a body fixed x-y-z
    Euler angle sequence in radians, followed by a position vector in the
    last three elements. The default is q=0 unless overridden.  
    @see setDefaultTransform() **/
    Bushing& setDefaultQ(const Vec6& q);

    /** Obtain the current value for this mobilizer's generalized coordinates 
    q from the given State, which must be realized through Stage::Model (these
    are state variables). **/ 
    const Vec6& getQ(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized coordinate 
    time derivatives qdot from the given State, which must be realized through
    Stage::Velocity. These are the three Euler angle time derivatives followed
    by the three measure numbers of the velocity vector v_FM giving the 
    velocity of M's origin Mo taken in the F frame. For the %Bushing
    mobilizer, qdot has the same numerical values as the generalized speeds u,
    but unlike u, qdot is not a state variable. **/ 
    const Vec6& getQDot(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized coordinate 
    second time derivatives qdotdot from the given State, which must be 
    realized through Stage::Acceleration. These are Euler angle second time
    derivatives followed by the three measure numbers of the vector a_FM
    giving the acceleration of M's origin Mo taken in the F frame. For the 
    %Bushing mobilizer, qdotdot=udot. **/ 
    const Vec6& getQDotDot(const State& state) const;

    /** Obtain the current value for this mobilizer's generalized speeds u
    from the given State, which must be realized through Stage::Model (these are
    state variables). For the %Bushing joint, these are the time derivatives of  
    q, that is, u=qdot. **/ 
    const Vec6& getU(const State& state) const;
    /** Obtain the current value for this mobilizer's generalized accelerations
    udot (=d/dt u) from the given State, which must be realized through
    Stage::Acceleration. For the %Bushing joint, these are the time derivatives 
    of qdot, that is, udot=qdotdot. **/ 
    const Vec6& getUDot(const State& state) const;

    /** Set new values for this mobilizer's generalized coordinates q in the
    given \a state. This invalidates Stage::Position and above. The new
    value is given as a Vec6, but is interpreted as a body fixed x-y-z
    Euler angle sequence in radians, followed by the three measure numbers of
    the position vector p_FM giving the location of the M frame origin Mo 
    measured from and expressed in the F frame. 
    @see MobilizedBody::setQToFitTransform() **/
    void setQ(State& state, const Vec6& q) const;
    /** Set new values for this mobilizer's generalized speeds u in the given 
    \a state. This invalidates Stage::Velocity and above. The new value is given
    as a Vec6, but is interpreted as the time derivatives of a body fixed x-y-z 
    Euler angle sequence in radians/time (\e not an angular velocity!), followed
    by the three measure numbers of the velocity vector v_FM giving the velocity
    of the M frame origin Mo measured and expressed in the F frame. For the
    %Bushing mobilizer, these are just the time derivatives of the generalized
    coordinates q, that is, u=qdot. 
    @see MobilizedBody::setUToFitVelocity() **/
    void setU(State& state, const Vec6& u) const;

    /** @name               Advanced/Obscure
    Most users won't use these methods. **/
    /**@{**/
    /** Create a disembodied %Bushing mobilizer that is not part of any System. **/
    explicit Bushing(Direction=Forward);

    /** Given a vector in the System's generalized coordinate basis, extract
    the six q's belonging to this mobilizer. **/
    const Vec6& getMyPartQ(const State& state, const Vector& qlike) const;
    /** Given a vector in the System's generalized speed basis, extract
    the six u's belonging to this mobilizer. **/
    const Vec6& getMyPartU(const State& state, const Vector& ulike) const;
   
    /** Given a writable vector in the System's generalized coordinate basis, 
    extract a writable reference to the six q's in that vector that belong to
    this mobilizer. **/
    Vec6& updMyPartQ(const State& state, Vector& qlike) const;
    /** Given a writable vector in the System's generalized speed basis, 
    extract a writable reference to the six u's in that vector that belong to
    this mobilizer. **/
    Vec6& updMyPartU(const State& state, Vector& ulike) const;
    /**@}**/

    /** @cond **/ // Hide from Doxygen.
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Bushing, BushingImpl, MobilizedBody);
    /** @endcond **/
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_BUSHING_H_



