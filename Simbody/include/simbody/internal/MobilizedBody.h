#ifndef SimTK_SIMBODY_MOBILIZED_BODY_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_H_

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
This defines the MobilizedBody class, which associates a new body (the 
"child", "outboard", or "successor" body) with a mobilizer and a reference 
frame on an existing body (the "parent", "inboard", or "predecessor" body)
that is already part of a SimbodyMatterSubsystem.

MobilizedBody is an abstract base class handle, with concrete classes 
defined for each kind of mobilizer. There are a set of built-in mobilizers
and a generic "Custom" mobilizer (an actual abstract base class) from
which advanced users may derive their own mobilizers.
**/

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/Motion.h"

#include <cassert>

namespace SimTK {

class SimbodyMatterSubsystem;
class Motion;
class MobilizedBody;
class MobilizedBodyImpl;

// We only want the template instantiation to occur once. This symbol is 
// defined in the SimTK core compilation unit that instantiates the mobilized 
// body class but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
    extern template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl, true>;
#endif

/** %Mobod is the approved abbreviation for MobilizedBody.\ Feel free to
use it if you get tired of typing or seeing the full name. 
@relates SimTK::MobilizedBody **/
typedef MobilizedBody Mobod;



//==============================================================================
//                             MOBILIZED BODY
//==============================================================================
/** A %MobilizedBody is Simbody's fundamental body-and-joint object used to
parameterize a system's motion by constructing a multibody tree containing
each body and its unique "mobilizer" (internal coordinate joint). A 
%MobilizedBody connects a new body (the "child", "outboard", or "successor" 
body) with a mobilizer and a reference frame to an existing %MobilizedBody (the 
"parent", "inboard", or "predecessor" body) that is already part of a 
SimbodyMatterSubsystem. In the topology of the multibody tree, the parent 
mobilized body is always closer to Ground ("inboard") than is the child.

This is the base class for all %MobilizedBody classes, which include a body
and a particular kind of mobilizer. Each built-in %MobilizedBody type is a local 
subclass within %MobilizedBody, so the built-ins have names like 
MobilizedBody::Pin. All concrete %MobilizedBody classes, including the 
built-ins, are derived from MobilizedBody. There is a MobilizedBody::Custom 
class available for defining your own mobilizer types.

Normally, a %MobilizedBody's motion is unknown and Simbody calculates it using
forward dynamics where the motion results from external or internal forces.
However, a %MobilizedBody may be locked or controlled by an associated Motion 
object which defines how it is to move; those are inverse dynamics conditions
in which the motion is known but the generalized forces required to produce that
motion are to be determined. In general a system may contain a mix of forward-
and inverse-dynamics mobilizers, and mobilizers may be switched between forward
and inverse dynamics during simulations.

The %MobilizedBody class supports a large number of methods useful for
obtaining and modifying body-specific and mobilizer-specific information and
for performing useful computations. These are organized in three sets:
   - %State Access
   - Basic Operators
   - High Level Operators

<em>%State Access</em> methods simply extract already-calculated data from the 
State or %State Cache, or set state values. They involve no additional 
computation, have names beginning with "get" and "upd" (update) and return 
references to the requested quantities rather than calculated values. We 
divide these into routines which deal with bodies and routines which deal 
with mobilizers and mobilities.

<em>Basic Operators</em> use State Access methods to compute basic quantities
which cannot be precomputed, such as the velocity of an arbitrary point, 
using an inline combination of basic floating point operations which can be 
reliably determined at compile time. These have names beginning with "find" 
or a more specific verb, as a reminder that they do not require a great deal 
of computation. 

<em>High Level Operators</em> combine responses and basic operators with 
run-time tests to calculate more complex quantities, with more complicated 
implementations that can exploit special cases at run time. These begin with 
"calc" (calculate) as a reminder that they may involve substantial run time 
computation.

There is also a set of methods used for construction, and miscellaneous 
utilities.

<h3>Mobilizer Terminology and Notation</h3>

Refer to the figure below for the terminology we use when discussing mobilizers 
and mobilized bodies. 

@image html MobilizerTerminology.png "Terminology and notation for mobilized bodies"

The figure shows the coordinate frames used in describing the mobility of 
%MobilizedBody B with respect to its inboard parent body P. Everything blue 
in the figure is associated with B. The origin point O of each frame is labeled.
A new %MobilizedBody with body frame B is added to the multibody tree by 
choosing a parent body P that is already present in the tree. There are two 
frames associated with the mobilizer: the "fixed" frame F that is attached to 
the parent, and the "moving" frame M that is attached to the new body B. Frame 
F is specified by giving its transform X_PF relative to the P frame. Frame M is 
specified by giving its transform X_BM relative to the B frame. At run time
the transform X_FM between the two mobilizer frames represents translation and
rotation of the mobilizer. That motion is parameterized via generalized 
coordinates q and generalized speeds u, the specific meaning of which is a
unique property of each type of mobilizer.

In the API below, we'll refer to the current ("this") MobilizedBody as "body B". 
It is the "object" or "main" body with which we are concerned. Often there will 
be another body mentioned in the argument list as a target for some conversion. 
That "another" body will be called "body A". The Ground body is abbreviated "G".

We use Fo to mean "the origin of frame F", Bc is "the mass center of body B". 
R_AF is the rotation matrix giving frame F's orientation in frame A, such that a 
vector v expressed in F is reexpressed in A by v_A = R_AF*v_F. X_AF is the 
spatial transform giving frame F's origin location and orientation in frame A, 
such that a point P whose location is measured from F's origin Fo and expressed 
in F by position vector p_FP (or more explicitly p_FoP) is remeasured from frame
A's origin Ao and reexpressed in A via p_AP = X_AF*p_FP, where p_AP==p_AoP. 

<h3>Theory</h3>
For the mathematical and computational theory behind Simbody's mobilizers, see
  - the <a href="https://simtk.org/docman/view.php/47/1536/Seth-2010-ShermanEastmanDelp-MinimalJointFormulationForBiomechanisms-NonlinearDyn-v62-p291.pdf">
    paper</a> Seth, A.; Sherman, M.A.; Eastman, P.; Delp, S.L. "Minimal formulation 
    of joint motion for biomechanisms" Nonlinear Dynamics 62:291-303 (2010), or
  - the Simbody Theory Manual.
**/

class SimTK_SIMBODY_EXPORT MobilizedBody 
:   public PIMPLHandle<MobilizedBody, MobilizedBodyImpl, true> {
public:

/** Constructors can take an argument of this type to indicate that the 
mobilizer is being defined in the reverse direction, meaning from the outboard
(child) body to the inboard (parent) body. That means that the mobilizer 
coordinates and speeds will be defined as though the tree had been built in the 
opposite direction. It does not actually affect which body is the inboard one
and which is the outboard; it just affects the definitions of the generalized
coordinates q and speeds u that parameterize the motion of the outboard body
with respect to the inboard one. This is a topological setting and can't be 
changed dynamically. **/
enum Direction {
    Forward = 0, ///< Use default definitions for q and u (inboard to outboard).
    Reverse = 1  ///< Define q and u in the reverse order (outboard to inboard).
};

//------------------------------------------------------------------------------
/**@name                Mobilizer locking and unlocking

Every %MobilizedBody object supports locking and unlocking of the mobilizer it
contains. You can lock the mobilizer's position, velocity, or acceleration. In 
all cases the generalized accelerations udot of a locked mobilizer are 
prescribed. If you lock just the accelerations, then velocity and position 
remain free. If you lock velocity, then the generalized speeds u are prescribed 
to specified values, accelerations udot are prescribed to zero, and positions 
will remain free. If you lock position (the default locking level) then the 
generalized coordinates q are prescribed to specified values and the speeds u 
and accelerations udot are prescribed to zero. Prescribed values may be obtained
from the current state, or set explicitly.

You can also specify that a mobilizer is locked by default (at acceleration,
velocity, or position level). In that case the prescribed value is recorded
when the state is realized to Stage::Model, using values taken from the
state at that time.

If this mobilizer is driven by a Motion object, locking overrides that while
the lock is active; when unlocked the Motion object resumes control. **/
/**@{**/
/** Lock the contained mobilizer's position or velocity, or lock the 
acceleration to zero, depending on the \p level parameter. The prescribed 
position or velocity is taken from \p state and recorded. If the mobilizer was 
already locked, the new lock specification overrides it. The \p state must 
already have been realized to at least Stage::Model and will be lowered to 
Stage::Model on return since a lock is an Instance-stage change. **/ 
void lock(State& state, Motion::Level level=Motion::Position) const;

/** Lock this mobilizer's q, u, or udot to the given scalar \p value, depending 
on \p level. This is only permitted for mobilizers that have just a single
q, u, or udot, such as a Pin or Slider mobilizer. **/
void lockAt(State& state, Real value, 
            Motion::Level level=Motion::Position) const;

/** Lock this mobilizer's q, u, or udot to the given Vector \p value, depending 
on \p level. The Vector must be the expected length, either nq=getNumQ() for
the default Motion::Position level, or nu=getNumU() for Motion::Velocity or
Motion::Acceleration levels. 
@see getNumQ(), getNumU() **/
void lockAt(State& state, const Vector& value, 
            Motion::Level level=Motion::Position) const;

/** Lock this mobilizer's q, u, or udot to the given Vec\<N\> \p value, 
depending on the \p level. The size N of the given Vec must match the expected 
length, either nq=getNumQ() for the default Motion::Position level, or 
nu=getNumU() for Motion::Velocity or Motion::Acceleration levels.  
@see getNumQ(), getNumU() **/
template <int N> SimTK_SIMBODY_EXPORT // instantiated in library
void lockAt(State& state, const Vec<N>& value,
            Motion::Level level=Motion::Position) const;

/** Unlock this mobilizer, returning it to its normal behavior which may be
free motion or may be controlled by a Motion object. If the mobilizer is
already unlocked nothing happens. **/
void unlock(State& state) const;

/** Check whether this mobilizer is currently locked in the given \p state. **/
bool isLocked(const State& state) const 
{   return getLockLevel(state)!=Motion::NoLevel; }

/** Returns the lock level if the mobilizer is locked in the given \p state, 
otherwise Motion::NoLevel. **/
Motion::Level getLockLevel(const State& state) const;

/** Return the q, u, or udot value at which this mobilizer is locked, depending
on the lock level, as a Vector of the appropriate length. If the mobilizer is
not currently locked, a zero-length Vector is returned. **/
Vector getLockValueAsVector(const State& state) const;

/** Change whether this mobilizer is initially locked. On construction the
mobilizer is unlocked by default; you can change that here. To respecify that 
the motion is unlocked, use Motion::NoLevel as the \p level. This is a 
topological change; you'll have to call realizeTopology() again and get a new 
State if you call this method. If the default lock is at the position or 
velocity level, the required q or u value is recorded at the time a \p state is
realized to Stage::Model. A default Motion::Acceleration lock always prescribes
the accelerations udot to zero. **/
MobilizedBody& lockByDefault(Motion::Level level=Motion::Position);

/** Check whether this mobilizer is to be locked in the default state. **/
bool isLockedByDefault() const 
{   return getLockByDefaultLevel()!=Motion::NoLevel; }

/** Returns the level at which the mobilizer is locked by default, if it is
locked by default, otherwise Motion::NoLevel. **/
Motion::Level getLockByDefaultLevel() const;
/**@}**/

    //////////////////////////
    // STATE ACCESS METHODS //
    //////////////////////////

//------------------------------------------------------------------------------
/** @name                    State Access - Bodies
These methods extract already-computed information from the State or State 
cache, or set values in the State. **/
/**@{**/

/** Extract from the state cache the already-calculated spatial configuration 
X_GB of body B's body frame, measured with respect to the Ground frame and 
expressed in the Ground frame. That is, we return the location of the body 
frame's origin, and the orientation of its x, y, and z axes, as the Transform 
X_GB. This notation is intended to convey unambiguously the sense of this 
transform, which is as follows: if you have a station (body fixed point) S on 
body B, represented by position vector p_BS (a.k.a. p_BoS) from the origin Bo 
of B to the point S and expressed in the B frame, then p_GS=X_GB*p_BS where 
p_GS (== p_GoS) is the position vector from the Ground origin Go to the point in
space currently coincident with S and expressed in the Ground frame. The inverse 
transformation is obtained using the "~" operator where ~X_GB=X_BG, so that 
p_BS = ~X_GB*p_GS. This response is available at Position stage. **/
const Transform& getBodyTransform(const State& state) const; // X_GB

/** Extract from the state cache the already-calculated spatial orientation R_GB 
of body B's body frame x, y, and z axes expressed in the Ground frame, as the 
Rotation matrix R_GB. The sense of this rotation matrix is such that if you have
a vector v fixed on body B, represented by the vector v_B expressed in the B 
frame, then v_G=R_GB*v_B where v_G is the same vector but re-expressed in the 
Ground frame. The inverse transformation is obtained using the "~" operator 
where ~R_GB=R_BG, so that v_B = ~R_GB*v_G. This response is available at 
Position stage. **/
const Rotation& getBodyRotation(const State& state) const {
    return getBodyTransform(state).R();
}
/** Extract from the state cache the already-calculated spatial location of body
B's body frame origin Bo, measured from the Ground origin Go and expressed in 
the Ground frame, as the position vector p_GB (== p_GoBo). This response is 
available at Position stage. **/
const Vec3& getBodyOriginLocation(const State& state) const {
    return getBodyTransform(state).p();
}

/** At stage Position or higher, return the cross-mobilizer transform X_FM, the 
body's inboard mobilizer frame M measured and expressed in the parent body's 
corresponding outboard frame F. **/
const Transform& getMobilizerTransform(const State& state) const;   // X_FM

/** Extract from the state cache the already-calculated spatial velocity V_GB of
this body's reference frame B, measured with respect to the Ground frame and 
expressed in the Ground frame. That is, we return the linear velocity v_GB of 
the body frame's origin in G, and the body's angular velocity w_GB as the 
spatial velocity vector V_GB = {w_GB, v_GB}. This response is available at 
Velocity stage. **/
const SpatialVec& getBodyVelocity(const State& state) const;        // V_GB

/** Extract from the state cache the already-calculated inertial angular
velocity vector w_GB of this body B, measured with respect to the Ground frame
and expressed in the Ground frame. This response is available at Velocity 
stage. **/
const Vec3& getBodyAngularVelocity(const State& state) const {      // w_GB
    return getBodyVelocity(state)[0]; 
}
/** Extract from the state cache the already-calculated inertial linear velocity
vector v_GB (more explicitly, v_GBo) of this body B's origin point Bo, measured 
with respect to the Ground frame and expressed in the Ground frame. This 
response is available at Velocity stage. **/
const Vec3& getBodyOriginVelocity(const State& state) const {       // v_GB
    return getBodyVelocity(state)[1];
}

/** At stage Velocity or higher, return the cross-mobilizer velocity V_FM, the 
relative velocity of this body's "moving" mobilizer frame M in the parent body's 
corresponding "fixed" frame F, measured and expressed in F. Note that this isn't
the usual spatial velocity since it isn't expressed in G. **/
const SpatialVec& getMobilizerVelocity(const State& state) const;   // V_FM

/** Extract from the state cache the already-calculated spatial acceleration 
A_GB of this body's reference frame B, measured with respect to the Ground frame
and expressed in the Ground frame. That is, we return the linear acceleration 
a_GB of the body frame's origin in G, and the body's angular acceleration b_GB 
as the spatial acceleration vector A_GB = {b_GB, a_GB}. This response is 
available at Acceleration stage. **/
const SpatialVec& getBodyAcceleration(const State& state) const;    // A_GB

/** Extract from the state cache the already-calculated inertial angular
acceleration vector b_GB of this body B, measured with respect to the Ground 
frame and expressed in the Ground frame. This response is available at 
Acceleration stage. **/
const Vec3& getBodyAngularAcceleration(const State& state) const {  // b_GB
    return getBodyAcceleration(state)[0]; 
}

/** Extract from the state cache the already-calculated inertial linear
acceleration vector a_GB (more explicitly, a_GBo) of this body B's origin point 
Bo, measured with respect to the Ground frame and expressed in the Ground frame. 
This response is available at Acceleration stage. **/
const Vec3& getBodyOriginAcceleration(const State& state) const {   // a_GB
    return getBodyAcceleration(state)[1];
}

/** TODO: Not implemented yet -- any volunteers? At stage Acceleration, return
the cross-mobilizer acceleration A_FM, the relative acceleration of body B's 
"moving" mobilizer frame M in the parent body's corresponding "fixed" frame F, 
measured and expressed in F. Note that this isn't the usual spatial acceleration 
since it isn't expressed in G. **/
const SpatialVec& getMobilizerAcceleration(const State& state) const { // A_FM
    SimTK_ASSERT_ALWAYS(!"unimplemented method", 
"MobilizedBody::getMobilizerAcceleration() not yet implemented -- volunteers?");
    return *(new SpatialVec());
}

/** Return a reference to this body's mass properties in the State cache. The 
State must have been realized to Stage::Instance or higher. **/
const MassProperties& getBodyMassProperties(const State& state) const;

/** Return a reference to the already-calculated SpatialInertia of this  body, 
taken about the body's origin (\e not its mass center), and expressed in the 
Ground frame. The State must have been realized to Stage::Position or higher. **/
const SpatialInertia& getBodySpatialInertiaInGround(const State& state) const;

/** Return the mass of this body. The State must have been realized to 
Stage::Instance. **/
Real getBodyMass(const State& state) const {
    return getBodyMassProperties(state).getMass();
}

/** Return this body's center of mass station (i.e., the vector fixed in the 
body, going from body origin to body mass center, expressed in the body frame.)
The State must have been realized to Stage::Instance or higher. **/
const Vec3& getBodyMassCenterStation(const State& state) const {
    return getBodyMassProperties(state).getMassCenter();
}

/** Return a reference to this body's unit inertia matrix in the State cache, 
taken about the body origin and expressed in the body frame. The State must have
been realized to Stage::Instance or higher. **/
const UnitInertia& getBodyUnitInertiaAboutBodyOrigin(const State& state) const {
    return getBodyMassProperties(state).getUnitInertia();
}

/** Return a reference to this mobilizer's frame F fixed on the parent body P, 
as the fixed Transform from P's body frame to the frame F fixed to P. If this 
frame is changeable, the result comes from the State cache, otherwise it is from
the MobilizedBody object itself. The State must have been realized to 
Stage::Instance or higher. **/
const Transform& getInboardFrame (const State& state) const;    // X_PF
/** Return a reference to this MobilizedBody's mobilizer frame M, as the fixed 
Transform from this body B's frame to the frame M fixed on B. If this frame is 
changeable, the result comes from the State cache, otherwise it is from the 
MobilizedBody object itself. The State must have been realized to 
Stage::Instance or higher. **/
const Transform& getOutboardFrame(const State& state) const;    // X_BM

/** TODO: not implemented yet. Set the location and orientation of the inboard 
(parent) mobilizer frame F, fixed to this mobilizer's parent body P.
@see setDefaultInboardFrame() **/
void setInboardFrame (State& state, const Transform& X_PF) const;
/** TODO: not implemented yet. Set the location and orientation of the outboard 
mobilizer frame M, fixed to this body B.
@see setDefaultOutboardFrame() **/
void setOutboardFrame(State& state, const Transform& X_BM) const;

// End of State Access - Bodies
/**@}**/


//------------------------------------------------------------------------------
/** @name   State Access - Mobilizer generalized coordinates q and speeds u
These methods extract q- or u-related information from the State or State cache, 
or set q or u values in the State. **/
/**@{**/
/** Return the number of generalized coordinates q currently in use by this 
mobilizer. State must have been realized to Stage::Model. **/
int getNumQ(const State& state) const;
/** Return the number of generalized speeds u currently in use by this 
mobilizer. State must have been realized to Stage::Model. **/
int getNumU(const State& state) const;

/** Return the global QIndex of the first q for this mobilizer; all the q's
range from getFirstQIndex() to QIndex(getFirstQIndex()+getNumQ()-1). **/
QIndex getFirstQIndex(const State& state) const;

/** Return the global UIndex of the first u for this mobilizer; all the u's
range from getFirstUIndex() to UIndex(getFirstUIndex()+getNumU()-1). **/
UIndex getFirstUIndex(const State& state) const;

/** Determine how generalized coordinate q values are being determined.
@param[in] state    Must be realized to Instance stage. **/
Motion::Method getQMotionMethod(const State& state) const;
/** Determine how generalized speed u values are being determined.
@param[in] state    Must be realized to Instance stage. **/
Motion::Method getUMotionMethod(const State& state) const;
/** Determine how generalized acceleration udot values are being determined.
@param[in] state    Must be realized to Instance stage. **/
Motion::Method getUDotMotionMethod(const State& state) const;

/** Return one of the generalized coordinates q from this mobilizer's partition 
of the matter subsystem's full q vector in the State. The particular coordinate 
is selected using the \p which parameter, numbering from zero to 
getNumQ()-1. **/
Real getOneQ(const State& state, int which) const;

/** Return one of the generalized speeds u from this mobilizer's partition of 
the matter subsystem's full u vector in the State. The particular coordinate is 
selected using the \p which parameter, numbering from zero to getNumU()-1. **/
Real getOneU(const State& state, int which) const;

/** Return as a Vector of length getNumQ() all the generalized coordinates q
currently in use by this mobilizer, from this mobilizer's partion in the matter 
subsystem's full q vector in the State. **/
Vector getQAsVector(const State& state) const;
/** Return as a Vector of length getNumU() all the generalized speeds u 
currently in use by this mobilizer, from this mobilizer's partion in the matter 
subsystem's full u vector in the State. **/
Vector getUAsVector(const State& state) const;

/** Return one of the generalized coordinate derivatives qdot from this 
mobilizer's partition of the matter subsystem's full qdot vector in the State 
cache. The particular coordinate is selected using the \p which parameter, 
numbering from zero to getNumQ()-1. **/
Real getOneQDot   (const State& state, int which) const;
/** Return as a Vector of length getNumQ() all the generalized coordinate 
derivatives qdot currently in use by this mobilizer, from this mobilizer's 
partition in the matter subsystem's full qdot vector in the State cache. **/
Vector getQDotAsVector(const State& state) const;

/** Return one of the generalized accelerations udot from this mobilizer's 
partition of the matter subsystem's full udot vector in the State cache. The 
particular coordinate is selected using the \p which parameter, numbering from 
zero to getNumU()-1. **/
Real getOneUDot(const State& state, int which) const;
/** Return one of the generalized coordinate second derivatives qdotdot from 
this mobilizer's partition of the matter subsystem's full qdotdot vector in the 
State cache. The particular coordinate is selected using the \p which parameter,
numbering from zero to getNumQ()-1. **/
Real getOneQDotDot(const State& state, int which) const;
/** Return as a Vector of length getNumU() all the generalized accelerations 
udot currently in use by this mobilizer, from this mobilizer's partion in the 
matter subsystem's full udot vector in the State cache. **/
Vector getUDotAsVector(const State& state) const;
/** Return as a Vector of length getNumQ() all the generalized coordinate 
second derivatives qdotdot currently in use by this mobilizer, from this 
mobilizer's partion in the matter subsystem's full qdotdot vector in the 
State cache. **/
Vector getQDotDotAsVector(const State& state) const;

/** Return the generalized forces tau resulting from prescribed (known) 
acceleration, corresponding to each of this mobilizer's mobilities, as a Vector 
of length nu=getNumU().

If this mobilizer has prescribed accelerations udot due to an active lock or
Motion object, the set of generalized forces tau that must be added in order to 
produce those accelerations is calculated at Acceleration stage. There is one 
scalar tau per mobility and they can be returned individually or as a Vector. 
The return value is zero if this mobilizer is currently free. **/
Vector getTauAsVector(const State& state) const;

/** Return one of the tau forces resulting from prescribed (known) acceleration,
corresponding to one of this mobilizer's mobilities as selected here using the 
\p which parameter, numbered from zero to getNumU()-1.
@see getTauAsVector() for more information **/
Real getOneTau(const State& state, MobilizerUIndex which) const;

/** Set one of the generalized coordinates q to value \p v, in this mobilizer's 
partition of the matter subsystem's full q vector in the State. The particular 
coordinate is selected using the \p which parameter, numbering from zero to 
getNumQ()-1. **/
void setOneQ(State& state, int which, Real v) const;
/** Set one of the generalized speeds u to value \p v, in this mobilizer's 
partition of the matter subsystem's full u vector in the State. The particular 
coordinate is selected using the \p which parameter, numbering from zero to 
getNumU()-1. **/
void setOneU(State& state, int which, Real v) const;

/** Set all of the generalized coordinates q to value \p v (a Vector of length 
getNumQ()), in this mobilizer's partition of the matter subsystem's full q 
vector in the State. **/
void setQFromVector(State& state, const Vector& v) const;
/** Set all of the generalized speeds u to value \p v (a Vector of length 
getNumU()), in this mobilizer's partition of the matter subsystem's full u 
vector in the State. **/
void setUFromVector(State& state, const Vector& v) const;

/** Adjust this mobilizer's q's to best approximate the supplied Transform
which requests a particular relative orientation and translation between the 
F "fixed" frame and M "moving" frame connected by this mobilizer.

This set of methods sets the generalized coordinates, or speeds (state
variables) for just the mobilizer associated with this MobilizedBody (ignoring 
all other mobilizers and constraints), without requiring knowledge of the 
meanings of the individual state variables. The idea here is to provide a 
physically-meaningful quantity relating the mobilizer's inboard and outboard 
frames, and then ask the mobilizer to set its state variables to reproduce that 
quantity to the extent it can.

These methods can be called in Stage::Model, however the routines may consult 
the current values of the state variables in some cases, so you must make sure 
they have been set to reasonable, or at least innocuous values (zero will work). 
In no circumstance will any of these methods look at any state variables that 
belong to another mobilizer; they are limited to working locally with just the
current mobilizer.

Routines which specify only translation (linear velocity) may use rotational 
coordinates to help satisfy the translation requirement. An alternate "Only" 
method is available to forbid modification of purely rotational coordinates in 
that case. When a mobilizer uses state variables which have combined rotational 
and translational character (e.g. a screw joint) consult the documentation for 
the derived MobilizedBody class to find out how that mobilizer responds to these
routines.

There is no guarantee that the desired physical quantity will be achieved by 
these routines; you can check on return if you're worried. Individual mobilizers 
make specific promises about what they will do; consult the documentation. These
routines do not throw exceptions even for absurd requests like specifying a
rotation for a sliding mobilizer. Nothing happens if there are no mobilities 
here, i.e. Ground or a Weld mobilizer. **/
void setQToFitTransform(State& state, const Transform& X_FM) const;

/** Adjust this mobilizer's q's to best approximate the supplied Rotation matrix
which requests a particular relative orientation between the "fixed" frame F and
"moving" frame M connected by this mobilizer.
@see setQToFitTransform() **/
void setQToFitRotation(State& state, const Rotation& R_FM) const;

/** Adjust this mobilizer's q's to best approximate the supplied position vector
which requests a particular offset between the origins of the F "fixed" frame
and M "moving" frame connected by this mobilizer, with \e any q's (rotational
or translational) being modified if doing so helps satisfy the request.
@see setQToFitTransform() **/
void setQToFitTranslation(State& state, const Vec3& p_FM) const;

/** Adjust this mobilizer's u's (generalized speeds) to best approximate the 
supplied spatial velocity \p V_FM which requests the relative angular and linear 
velocity between the "fixed" and "moving" frames connected by this mobilizer. 
Routines which affect generalized speeds u depend on the generalized coordinates
q already having been set; they never change these coordinates.
@see setQToFitTransform() **/
void setUToFitVelocity(State& state, const SpatialVec& V_FM) const;

/** Adjust this mobilizer's u's (generalized speeds) to best approximate the 
supplied angular velocity \p w_FM which requests a particular relative angular
between the "fixed" and "moving" frames connected by this mobilizer.
@see setQToFitTransform(), setUToFitVelocity() **/
void setUToFitAngularVelocity(State& state, const Vec3& w_FM) const;

/** Adjust <em>any</em> of this mobilizer's u's (generalized speeds) to best 
approximate the supplied linear velocity \p v_FM which requests a particular 
velocity for the "moving" frame M origin in the "fixed" frame F on the parent 
where these are the frames connected by this mobilizer.
@see setQToFitTransform(), setUToFitVelocity() **/
void setUToFitLinearVelocity(State& state, const Vec3& v_FM) const;

/** Expert use only: obtain a column of the hinge matrix H corresponding to one 
of this mobilizer's mobilities (actually a column of H_PB_G; what Jain calls H* 
and Schwieters calls H^T). This is the matrix that maps generalized speeds u to 
the cross-body relative spatial velocity V_PB_G via V_PB_G=H*u. Note that 
although H relates child body B to parent body B, it is expressed in the ground 
frame G so the resulting cross-body velocity of B in P is also expressed in G. 
The supplied state must have been realized through Position stage because H 
varies with this mobilizer's generalized coordinates q.
@see getH_FMCol() **/
SpatialVec getHCol(const State& state, MobilizerUIndex ux) const;

/** Expert use only: obtain a column of the mobilizer-local hinge matrix H_FM 
which maps generalized speeds u to cross-mobilizer spatial velocity V_FM via 
V_FM=H_FM*u. Note that H and V here are expressed in the parent body's (inboard)
frame F. The supplied state must have been realized through Position stage 
because H varies with this mobilizer's generalized coordinates q.
@see getHCol() **/
SpatialVec getH_FMCol(const State& state, MobilizerUIndex ux) const;

// End of State Access Methods.
/**@}**/ 

    /////////////////////
    // BASIC OPERATORS //
    /////////////////////

//------------------------------------------------------------------------------
/** @name                        Basic Operators

These methods use state variables and Response methods to compute basic
quantities which cannot be precomputed, but which can be implemented with an 
inline combination of basic floating point operations which can be reliably 
determined at compile time. The method names and descriptions use the following 
terms:
 - %Body or ThisBody: the %Body B associated with the current %MobilizedBody. 
   ThisBody is implied when no other Body is mentioned.
 - %Ground: the "MobilizedBody" G representing the %Ground reference frame which 
   never moves.
 - AnotherBody: the %Body A being referenced, which in general is neither 
   ThisBody nor %Ground.
 - Station: a point S fixed on ThisBody B, located by a position vector p_BS (or
   more explicitly, p_BoS) from the B-frame origin Bo to the point S, expressed 
   in the B-frame coordinate system.
 - %Vector: a vector v fixed on ThisBody B, given by a vector v_B expressed in 
   the B-frame coordinate system.
 - Direction: a unit vector u fixed on ThisBody B, given by a unit vector u_B 
   expressed in the B-frame coordinate system.
 - Frame: an origin and coordinate axes F fixed on ThisBody B, given by a 
   transform X_BF that locates F's origin (a Station) in B and expresses each of
   F's axes (Directions) in B.
 - Origin: the Station located at (0,0,0) in ThisBody frame B, that is, body B's
   origin point.
 - MassCenter: the Station on ThisBody B which is the center of mass for B.
 - GroundPoint, GroundVector: a Point P or Vector v on the %Ground "Body" G. 
   These are measured and expressed in the %Ground frame, as p_GP or v_G.
 - AnotherBodyStation, AnotherBodyVector, etc.: a Station S or %Vector v on 
   AnotherBody A. These are measured and expressed in the A frame, as p_AS 
   or v_A.
 - Mobilizer frame M: the mobilizer's outboard "moving" frame, fixed to 
   ThisBody B.
 - Mobilizer frame F: the mobilizer's inboard "fixed" frame, fixed to the parent 
   body P. **/

/**@{**/ 

/** Return X_AB, the spatial transform giving this body B's frame in body A's 
frame. Cost is 63 flops. If you know that one of the bodies is Ground, use the 
0-cost response getBodyTransform() instead of this operators. This operator is 
available in Position stage.
@see getBodyTransform() **/
Transform findBodyTransformInAnotherBody(const State& state, 
                                         const MobilizedBody& inBodyA) const
{
    const Transform& X_GA = inBodyA.getBodyTransform(state);
    const Transform& X_GB = this->getBodyTransform(state);

    return ~X_GA*X_GB; // X_AB=X_AG*X_GB
}

/** Return R_AB, the rotation matrix giving this body B's axes in body A's 
frame. Cost is 45 flops. If you know that one of the bodies is Ground, use the 
0-cost response getBodyRotation() instead of this operators. This operator is 
available in Position stage.
@see getBodyRotation() **/
Rotation findBodyRotationInAnotherBody(const State& state, 
                                       const MobilizedBody& inBodyA) const
{
    const Rotation& R_GA = inBodyA.getBodyRotation(state);
    const Rotation& R_GB = this->getBodyRotation(state);

    return ~R_GA*R_GB; // R_AB=R_AG*R_GB
}

/** Return the station on another body A (that is, a point measured and 
expressed in A) that is  currently coincident in space with the origin Bo of 
this body B. Cost is 18 flops. This operator is available at Position stage. 
Note: "findBodyOriginLocationInGround" doesn't exist because it would be the 
same as the no-cost response method getBodyOriginLocation().
@see getBodyOriginLocation() **/
Vec3 findBodyOriginLocationInAnotherBody
   (const State& state, const MobilizedBody& toBodyA) const {   
    return toBodyA.findStationAtGroundPoint(state,
                                            getBodyOriginLocation(state)); 
}

/** Return the angular and linear velocity of body B's frame in body A's frame, 
expressed in body A, and arranged as a SpatialVec. Cost is 51 flops. If you know
inBodyA is Ground, don't use this routine; use the response method 
getBodyVelocity() which is free. This operator is available in Velocity stage.
@see getBodyVelocity() **/
SpatialVec findBodyVelocityInAnotherBody
   (const State& state, const MobilizedBody& inBodyA) const
{
    const SpatialVec& V_GB   = this->getBodyVelocity(state);
    const SpatialVec& V_GA   = inBodyA.getBodyVelocity(state);
    // angular velocity of B in A, exp in G
    const Vec3        w_AB_G = V_GB[0]-V_GA[0];             //  3 flops 

    // Angular vel. was easy, but for linear vel. we have to add in an wXr term.

    const Transform&  X_GB       = getBodyTransform(state);
    const Transform&  X_GA       = inBodyA.getBodyTransform(state);
    // vector from Ao to Bo, exp in G
    const Vec3        p_AB_G     = X_GB.p() - X_GA.p();     //  3 flops
    // d/dt p taken in G
    const Vec3        p_AB_G_dot = V_GB[1]  - V_GA[1];      //  3 flops

    // d/dt p taken in A, exp in G 
    const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G;      // 12 flops

    // We're done, but answer is expressed in Ground. Reexpress in A and return.
    return ~X_GA.R()*SpatialVec(w_AB_G, v_AB_G);            // 30 flops
}

/** Return the angular velocity w_AB of body B's frame in body A's frame, 
expressed in body A. Cost is 18 flops. If you know inBodyA is Ground, don't use 
this routine; use the response method getBodyAngularVelocity() which is free. 
This operator is available in Velocity stage.
@see getBodyAngularVelocity() **/
Vec3 findBodyAngularVelocityInAnotherBody(const State& state,
                                          const MobilizedBody& inBodyA) const 
{
    const Vec3& w_GB   = getBodyAngularVelocity(state);
    const Vec3& w_GA   = inBodyA.getBodyAngularVelocity(state);
    // angular velocity of B in A, exp in G
    const Vec3  w_AB_G = w_GB-w_GA;                                 //  3 flops

    // Now reexpress in A.
    return inBodyA.expressGroundVectorInBodyFrame(state, w_AB_G);   // 15 flops
}

/** Return the velocity of body B's origin point in body A's frame, expressed in
body A. Cost is 51 flops. If you know inBodyA is Ground, don't use this routine; 
use the response method getBodyOriginVelocity() which is free. This operator is 
available in Velocity stage.
@see getBodyOriginVelocity() **/
Vec3 findBodyOriginVelocityInAnotherBody(const State& state,
                                         const MobilizedBody& inBodyA) const
{
    // Doesn't save much to special case this one.
    return findBodyVelocityInAnotherBody(state,inBodyA)[1];
}

/** Return the angular and linear acceleration of body B's frame in body A's 
frame, expressed in body A, and arranged as a SpatialVec. Cost is 105 flops. If 
you know that inBodyA is Ground, don't use this operator; instead, use the 
response method getBodyAcceleration() which is free. This operator is available 
in Acceleration stage.
@see getBodyAcceleration() **/
SpatialVec findBodyAccelerationInAnotherBody(const State& state,
                                             const MobilizedBody& inBodyA) const
{
    const Transform&  X_GA = inBodyA.getBodyTransform(state);
    const SpatialVec& V_GA = inBodyA.getBodyVelocity(state);
    const SpatialVec& A_GA = inBodyA.getBodyAcceleration(state);
    const Transform&  X_GB = this->getBodyTransform(state);
    const SpatialVec& V_GB = this->getBodyVelocity(state);
    const SpatialVec& A_GB = this->getBodyAcceleration(state);

    return findRelativeAcceleration(X_GA, V_GA, A_GA,
                                    X_GB, V_GB, A_GB);
}

/** Return the angular acceleration of body B's frame in body A's frame, 
expressed in body A. Cost is 33 flops. If you know \p inBodyA is Ground, don't 
use this operator; instead use the response method getBodyAngularAccleration() 
which is free. This operator is available at AccelerationStage.
@see getBodyAngularAcceleration() **/
Vec3 findBodyAngularAccelerationInAnotherBody
   (const State& state, const MobilizedBody& inBodyA) const
{
    const Rotation& R_GA = inBodyA.getBodyRotation(state);
    const Vec3&     w_GA = inBodyA.getBodyAngularVelocity(state);
    const Vec3&     w_GB = this->getBodyAngularVelocity(state);
    const Vec3&     b_GA = inBodyA.getBodyAngularAcceleration(state);
    const Vec3&     b_GB = this->getBodyAngularAcceleration(state);

    // relative ang. vel. of B in A, exp. in G
    const Vec3 w_AB_G     = w_GB - w_GA; // 3 flops
    const Vec3 w_AB_G_dot = b_GB - b_GA; // d/dt of w_AB_G taken in G; 3 flops

    // We have the derivative in G; change it to derivative in A by adding 
    // in contribution caused by motion of G in A, that is w_AG X w_AB_G. 
    // (Note that w_AG=-w_GA.)
    const Vec3 b_AB_G = w_AB_G_dot - w_GA % w_AB_G; // ang. accel. of B in A
                                                    // 12 flops

    return ~R_GA * b_AB_G; // taken in A, expressed in A; 15 flops
}

/** Return the acceleration of body B's origin point in body A's frame, 
expressed in body A. Cost is 105 flops. If you know that inBodyA is Ground, 
don't use this operator; instead, use the response method 
getBodyOriginAcceleration() which is free. This operator is available at 
Acceleration stage.
@see getBodyOriginAcceleration() **/
Vec3 findBodyOriginAccelerationInAnotherBody(const State& state, 
                                             const MobilizedBody& inBodyA) const
{
    // Not much to be saved by trying to optimize this since the linear part
    // is the most expensive to calculate.
    return findBodyAccelerationInAnotherBody(state,inBodyA)[1];
}

/** Return the spatial reaction force (moment and force) applied by the 
mobilizer to body B at the location of the mobilizer frame M (fixed to body B, 
but not necessarily at the body frame origin), expressed in Ground. This 
operator is available at Acceleration stage. Cost is about 120 flops.
@see findMobilizerReactionOnParentAtFInGround()
@see findMobilizerReactionOnBodyAtOriginInGround()
@see SimTK::SimbodyMatterSubsystem::calcMobilizerReactionForces() **/
SpatialVec findMobilizerReactionOnBodyAtMInGround(const State& state) const;

/** Return the spatial reaction force (moment and force) applied by the 
mobilizer to body B but shifted to the B frame origin, and expressed in Ground. 
This operator is available at Acceleration stage. Cost is about 90 flops.
@see findMobilizerReactionOnParentAtOriginInGround()
@see findMobilizerReactionOnBodyAtMInGround()
@see SimTK::SimbodyMatterSubsystem::calcMobilizerReactionForces() **/
SpatialVec findMobilizerReactionOnBodyAtOriginInGround
   (const State& state) const;

/** Return the spatial reaction force (moment and force) applied by the 
mobilizer to the parent (inboard) body P at the location of the inboard "fixed" 
mobilizer frame F (fixed to body P, but not necessarily at the P frame origin), 
expressed in Ground. This operator is available at Acceleration stage. Cost is 
about 140 flops.
@see findMobilizerReactionOnBodyAtMInGround()
@see findMobilizerReactionOnParentAtOriginInGround()
@see SimTK::SimbodyMatterSubsystem::calcMobilizerReactionForces() **/
SpatialVec findMobilizerReactionOnParentAtFInGround(const State& state) const;

/** Return the spatial reaction force (moment and force) applied by the 
mobilizer to the parent (inboard) body P at the location of the P frame origin, 
and expressed in Ground. This operator is available at Acceleration stage. Cost 
is about 110 flops.
@see findMobilizerReactionOnBodyAtOriginInGround()
@see findMobilizerReactionOnParentAtFInGround()
@see SimTK::SimbodyMatterSubsystem::calcMobilizerReactionForces() **/
SpatialVec findMobilizerReactionOnParentAtOriginInGround
   (const State& state) const;

/** Return the Cartesian (ground) location that is currently coincident with a 
station (point) S fixed on body B. That is, we return 
locationInG=X_GB*stationOnB which means the result is measured from the Ground 
origin and expressed in Ground. In more precise notation, we're calculating 
p_GS = X_GB * p_BS for position vectors p_GS and p_BS. Cost is 18 flops. This 
operator is available at Position stage. **/
Vec3 findStationLocationInGround
   (const State& state, const Vec3& stationOnB) const {
    return getBodyTransform(state) * stationOnB;
}


/** Given a station S on this body B, return the location on another body A 
which is at the same location in space. That is, we return 
locationOnA=X_AB*locationOnB, which means the result is measured from the body A 
origin and expressed in body A. In more precise notation, we're calculating 
p_AS = X_AB * p_BS, which we actually calculate as p_AS = X_AG*(X_GB*p_BS). Cost
is 36 flops. Note: if you know that one of the bodies is Ground, use one of the 
routines above which is specialized for Ground to avoid half the work. This 
operator is available at Position stage or higher. **/
Vec3 findStationLocationInAnotherBody
   (const State& state, const Vec3& stationOnB, const MobilizedBody& toBodyA) 
    const 
{
    return toBodyA.findStationAtGroundPoint(state, 
                                findStationLocationInGround(state,stationOnB));
}

/** Given a station fixed on body B, return its inertial (Cartesian) velocity, 
that is, its velocity relative to the Ground frame, expressed in the Ground 
frame. Cost is 27 flops. If you know the station is the body origin (0,0,0) 
don't use this routine; use the response getBodyOriginVelocity() which is free. 
This operator is available at Velocity stage.
@see getBodyOriginVelocity() **/
Vec3 findStationVelocityInGround
   (const State& state, const Vec3& stationOnB) const {
    const Vec3& w = getBodyAngularVelocity(state); // in G
    const Vec3& v = getBodyOriginVelocity(state);  // in G
    const Vec3  r = expressVectorInGroundFrame(state,stationOnB);   // 15 flops
    return v + w % r;                                               // 12 flops
}


/** Return the velocity of a station S fixed on body B, in body A's frame, 
expressed in body A. Cost is 93 flops. If you know \p inBodyA is Ground, don't 
use this operator; instead use the operator findStationVelocityInGround() which 
is much cheaper. This operator is available in Velocity stage.
@see findStationVelocityInGround() **/
Vec3 findStationVelocityInAnotherBody(const State&         state, 
                                      const Vec3&          stationOnBodyB,//p_BS
                                      const MobilizedBody& inBodyA) const
{
    const SpatialVec V_AB = 
        findBodyVelocityInAnotherBody(state, inBodyA); //51 flops
    // Bo->S rexpressed in A but not shifted to Ao
    const Vec3 p_BS_A = 
        expressVectorInAnotherBodyFrame(state, stationOnBodyB, inBodyA); //30
    return V_AB[1] + (V_AB[0] % p_BS_A);                              //12 flops
}

      
/** Given a station fixed on body B, return its inertial (Cartesian) 
acceleration, that is, its acceleration relative to the ground frame, expressed 
in the ground frame. Cost is 48 flops. If you know the station is the body 
origin (0,0,0) don't use this routine; use the response 
getBodyOriginAcceleration() which is free. This operator is available at 
Acceleration stage.
@see getBodyOriginAcceleration() **/
Vec3 findStationAccelerationInGround
   (const State& state, const Vec3& stationOnB) const {
    const Vec3& w = getBodyAngularVelocity(state);     // in G
    const Vec3& b = getBodyAngularAcceleration(state); // in G
    const Vec3& a = getBodyOriginAcceleration(state);  // in G

    const Vec3  r = expressVectorInGroundFrame(state,stationOnB); // 15 flops
    return a + b % r + w % (w % r);                               // 33 flops
}

/** Return the acceleration of a station S fixed on body B, in another body A's 
frame, expressed in body A. Cost is 186 flops. If you know that \p inBodyA is 
Ground, don't use this operator; instead, use the operator 
findStationAccelerationInGround() which is much cheaper. This operator is 
available in Acceleration stage. **/
Vec3 findStationAccelerationInAnotherBody(const State&         state,
                                          const Vec3&          stationOnBodyB, 
                                          const MobilizedBody& inBodyA) const 
{
    const Vec3       w_AB = 
        findBodyAngularVelocityInAnotherBody(state,inBodyA); // 18 flops
    const SpatialVec A_AB = 
        findBodyAccelerationInAnotherBody(state,inBodyA);    // 105 flops
    // Bo->S rexpressed in A but not shifted to Ao
    const Vec3       p_BS_A = 
        expressVectorInAnotherBodyFrame(state, stationOnBodyB, inBodyA); // 30

    return A_AB[1] + (A_AB[0] % p_BS_A) + w_AB % (w_AB % p_BS_A);    // 33 flops
}

/** It is cheaper to calculate a station's ground location and velocity together 
than to do them separately. Here we can return them both in 30 flops, vs. 45 to 
do them in two calls. This operator is available at Velocity stage. **/
void findStationLocationAndVelocityInGround
   (const State& state, const Vec3& locationOnB,
    Vec3& locationOnGround, Vec3& velocityInGround) const
{
    const Vec3& p_GB   = getBodyOriginLocation(state);
    const Vec3  p_BS_G = expressVectorInGroundFrame(state,locationOnB);//15flops
    locationOnGround = p_GB + p_BS_G;                                  // 3flops

    const Vec3& w_GB = getBodyAngularVelocity(state);
    const Vec3& v_GB = getBodyOriginVelocity(state);
    velocityInGround = v_GB + w_GB % p_BS_G;                         // 12 flops
}


/** It is cheaper to calculate a station's ground location, velocity, and 
acceleration together than to do them separately. Here we can return them all in
54 flops, vs. 93 to do them in three calls. This operator is available at 
Acceleration stage. **/
void findStationLocationVelocityAndAccelerationInGround
   (const State& state, const Vec3& locationOnB,
    Vec3& locationOnGround, Vec3& velocityInGround, Vec3& accelerationInGround) 
    const
{
    const Rotation&  R_GB = getBodyRotation(state);
    const Vec3&      p_GB = getBodyOriginLocation(state);

    // re-express station vector p_BS in G
    const Vec3 r = R_GB*locationOnB;  // 15 flops
    locationOnGround  = p_GB + r;     //  3 flops

    const Vec3& w = getBodyAngularVelocity(state);      // in G
    const Vec3& v = getBodyOriginVelocity(state);       // in G
    const Vec3& b = getBodyAngularAcceleration(state);  // in G
    const Vec3& a = getBodyOriginAcceleration(state);   // in G

    const Vec3 wXr = w % r; // "whipping" velocity w X r due to ang vel; 9 flops
    velocityInGround     = v + wXr;                 // v + w X r (3 flops)
    accelerationInGround = a + b % r + w % wXr;     // 24 flops
}

/** Return the Cartesian (ground) location of this body B's mass center. Cost is
18 flops. This operator is available at Position stage. **/
Vec3 findMassCenterLocationInGround(const State& state) const {
    return findStationLocationInGround(state,getBodyMassCenterStation(state));
}

/** Return the point of another body A that is currently coincident in space 
with the mass center CB of this body B. Cost is 36 flops. This operator is 
available at Position stage. **/
Vec3 findMassCenterLocationInAnotherBody
   (const State& state, const MobilizedBody& toBodyA) const {
    return findStationLocationInAnotherBody(state,
                                    getBodyMassCenterStation(state),toBodyA);
}

/** Return the station (point) S of this body B that is coincident with the
given Ground location. That is we return locationOnB = X_BG*locationInG, which 
means the result is measured from the body origin Bo and expressed in the body 
frame. In more precise notation, we're calculating p_BS = X_BG*p_GS. Cost is 18 
flops. This operator is available at 
Position stage. **/
Vec3 findStationAtGroundPoint
   (const State& state, const Vec3& locationInG) const {
    return ~getBodyTransform(state) * locationInG;
}

/** Return the station (point) on this body B that is coincident with the given 
station on another body A. That is we return stationOnB = X_BA * stationOnA, 
which means the result is measured from the body origin Bo and expressed in the 
body frame. Cost is 36 flops. This operator is available at Position stage.
@see findStationLocationInAnotherBody() **/
Vec3 findStationAtAnotherBodyStation
   (const State& state, const MobilizedBody& fromBodyA, 
    const Vec3& stationOnA) const {
    return fromBodyA.findStationLocationInAnotherBody(state, stationOnA, *this);
}

/** Return the station S of this body that is currently coincident in space with
the origin Ao of another body A. Cost is 18 flops. This operator is available at
Position stage. **/
Vec3 findStationAtAnotherBodyOrigin
   (const State& state, const MobilizedBody& fromBodyA) const {
    return findStationAtGroundPoint(state,
                                    fromBodyA.getBodyOriginLocation(state));
}

/** Return the station S of this body that is currently coincident in space with 
the mass center Ac of another body A. Cost is 36 flops. This operator is 
available at Position stage. **/
Vec3 findStationAtAnotherBodyMassCenter
   (const State& state, const MobilizedBody& fromBodyA) const {
    return fromBodyA.findStationLocationInAnotherBody(state,
                                        getBodyMassCenterStation(state),*this);
}

/** Return the current Ground-frame pose (position and orientation) of a frame 
F that is fixed to body B. That is, we return X_GF=X_GB*X_BF. Cost is 63 flops. 
This operator is available at Position stage. **/
Transform findFrameTransformInGround
    (const State& state, const Transform& frameOnB) const {
    return getBodyTransform(state) * frameOnB;
}

/** Return the current Ground-frame spatial velocity V_GF (that is, angular and 
linear velocity) of a frame F that is fixed to body B. The angular velocity of F
is the same as that of B, but the linear velocity is the velocity of F's origin 
Fo rather than B's origin Bo. This operator is available at Velocity stage. Cost
is 27 flops. **/
SpatialVec findFrameVelocityInGround
   (const State& state, const Transform& frameOnB) const {
    return SpatialVec(getBodyAngularVelocity(state),
                      findStationVelocityInGround(state,frameOnB.p()));
}

/** Return the current Ground-frame spatial acceleration A_GF (that is, angular 
and linear acceleration) of a frame F that is fixed to body B. The angular 
acceleration of F is the same as that of B, but the linear acceleration is the 
acceleration of F's origin Fo rather than B's origin Bo. This operator is 
available at Acceleration stage. Cost is 48 flops. **/
SpatialVec findFrameAccelerationInGround
    (const State& state, const Transform& frameOnB) const {
    return SpatialVec(getBodyAngularAcceleration(state),
                      findStationAccelerationInGround(state,frameOnB.p()));
}

/** Re-express a vector expressed in this body B's frame into the same vector in 
G, by applying only a rotation. That is, we return vectorInG = R_GB * vectorInB. 
Cost is 15 flops. This operator is available at Position stage. **/
Vec3 expressVectorInGroundFrame
   (const State& state, const Vec3& vectorInB) const {
    return getBodyRotation(state)*vectorInB;
}

/** Re-express a vector expressed in Ground into the same vector expressed in 
this body B, by applying only rotation. That is, we return 
vectorInB = R_BG * vectorInG. Cost is 15 flops. This operator is available at 
Position stage. **/
Vec3 expressGroundVectorInBodyFrame
   (const State& state, const Vec3& vectorInG) const {
    return ~getBodyRotation(state)*vectorInG;
}

/** Re-express a vector expressed in this body B into the same vector expressed 
in body A, by applying only a rotation. That is, we return 
vectorInA = R_AB * vectorInB. Cost is 30 flops. This operator is available at 
Position stage. Note: if you know one of the bodies is Ground, call one of the 
specialized methods above to save 15 flops. **/
Vec3 expressVectorInAnotherBodyFrame
   (const State& state, const Vec3& vectorInB, 
    const MobilizedBody& inBodyA) const
{
    return inBodyA.expressGroundVectorInBodyFrame(state, 
                                expressVectorInGroundFrame(state,vectorInB));
}

/** Re-express this body B's mass properties in Ground by applying only a 
rotation, not a shift of reference point. The mass properties remain measured in
the body B frame, taken about the body B origin Bo, but are reexpressed in 
Ground. **/
MassProperties expressMassPropertiesInGroundFrame(const State& state) const {
    const MassProperties& M_Bo_B = getBodyMassProperties(state);
    const Rotation&       R_GB   = getBodyRotation(state);
    return M_Bo_B.reexpress(~R_GB);
}

/** Re-express this body B's mass properties in another body A's frame by 
applying only a rotation, not a shift of reference point. The mass properties 
remain measured in the body B frame, taken about the body B origin Bo, but are 
reexpressed in A. **/
MassProperties expressMassPropertiesInAnotherBodyFrame
   (const State& state, const MobilizedBody& inBodyA) const {
    const MassProperties& M_Bo_B = getBodyMassProperties(state);
    const Rotation        R_AB   = findBodyRotationInAnotherBody(state,inBodyA);
    return M_Bo_B.reexpress(~R_AB);
}

// End of Basic Operators.
/**@}**/ 

//------------------------------------------------------------------------------
/** @name                    High-Level Operators
High level operators combine State Access and Basic Operators with run-time 
tests to calculate more complex %MobilizedBody-specific quantities, with more 
complicated implementations that can exploit special cases at run time. **/

/**@{**/
/** Return the mass properties of body B, measured from and about the B 
origin Bo, but expressed in Ground and then returned as a Spatial Inertia 
Matrix. The mass properties are arranged in the SpatialMat like this:
<pre>
        M=[      I_Bo      crossMat(m*Bc) ]
        [ ~crossMat(m*Bc)   diag(m)     ]
</pre>
where I_Bo is the inertia taken about the B frame origin Bo, and Bc is the 
vector p_BoBc from B's origin to its mass center.
    
The Spatial Inertia Matrix for Ground has infinite mass and inertia, with 
the cross terms set to zero. That is, it looks like a 6x6 diagonal matrix 
with Infinity on the diagonals.
    
@par Required stage
\c Stage::Position, unless \p objectBodyB == \c GroundIndex **/
SpatialMat calcBodySpatialInertiaMatrixInGround(const State& state) const
{
    if (isGround())
        return SpatialMat(Mat33(Infinity)); // sets diagonals to Inf

    const MassProperties& mp   = getBodyMassProperties(state);
    const Rotation&       R_GB = getBodyRotation(state);
        // re-express in Ground without shifting, convert to spatial mat.
    return mp.reexpress(~R_GB).toSpatialMat();
}

/** Return the central inertia for body B, that is, the inertia taken about
body B's mass center Bc, and expressed in B.
///
@par Required stage
  \c Stage::Instance **/
Inertia calcBodyCentralInertia(const State& state, 
                                MobilizedBodyIndex objectBodyB) const
{
    return getBodyMassProperties(state).calcCentralInertia();
}

/** Return the inertia of this body B, taken about an arbitrary point PA of 
body A, and expressed in body A.
TODO: this needs testing! **/
Inertia calcBodyInertiaAboutAnotherBodyStation
   (const State& state, const MobilizedBody& inBodyA, 
    const Vec3& aboutLocationOnBodyA) const
{
    // get B's mass props MB, measured about Bo, exp. in B
    const MassProperties& MB_Bo_B = getBodyMassProperties(state);

    // Calculate the vector from the body B origin (current "about" point) to 
    // the new "about" point PA, expressed in B.
    const Vec3 p_Bo_PA = 
        findStationAtAnotherBodyStation(state, inBodyA, aboutLocationOnBodyA);

    // Now shift the "about" point for body B's inertia IB to PA, but still 
    // expressed in B.
    const Inertia IB_PA_B = MB_Bo_B.calcShiftedInertia(p_Bo_PA);
        
    // Finally reexpress the inertia in the A frame.
    const Rotation R_BA    = 
        inBodyA.findBodyRotationInAnotherBody(state, *this);
    const Inertia  IB_PA_A = IB_PA_B.reexpress(R_BA);
    return IB_PA_A;
}


/** Calculate body B's momentum (angular, linear) measured and expressed in 
Ground, but taken about the body origin Bo. **/
SpatialVec calcBodyMomentumAboutBodyOriginInGround(const State& state) {
    const MassProperties M_Bo_G = expressMassPropertiesInGroundFrame(state);
    const SpatialVec&    V_GB   = getBodyVelocity(state);
    return M_Bo_G.toSpatialMat() * V_GB;
}

/** Calculate body B's momentum (angular, linear) measured and expressed in 
Ground, but taken about the body mass center Bc. **/
SpatialVec calcBodyMomentumAboutBodyMassCenterInGround
   (const State& state) const {
    const MassProperties& M_Bo_B = getBodyMassProperties(state);
    const Rotation&       R_GB   = getBodyRotation(state);

    // Given a central inertia matrix I, angular velocity w, and mass center 
    // velocity v, the central angular momentum is Iw and linear momentum is mv.
    const Inertia I_Bc_B = M_Bo_B.calcCentralInertia();
    const Inertia I_Bc_G = I_Bc_B.reexpress(~R_GB);
    const Real  mb    = M_Bo_B.getMass();
    const Vec3& w_GB  = getBodyAngularVelocity(state);
    Vec3        v_GBc = 
        findStationVelocityInGround(state, M_Bo_B.getMassCenter());

    return SpatialVec( I_Bc_G*w_GB, mb*v_GBc );
}

/** Calculate the distance from a station PB on body B to a station PA on 
body A. We are given the location vectors (stations) p_Bo_PB and p_Ao_PA, 
expressed in their respective frames. We return |p_PB_PA|. **/
Real calcStationToStationDistance(const State&         state,
                                  const Vec3&          locationOnBodyB,
                                  const MobilizedBody& bodyA,
                                  const Vec3&          locationOnBodyA) const
{
    if (isSameMobilizedBody(bodyA))
        return (locationOnBodyA-locationOnBodyB).norm();

    const Vec3 r_Go_PB = 
        this->findStationLocationInGround(state,locationOnBodyB);
    const Vec3 r_Go_PA = 
        bodyA.findStationLocationInGround(state,locationOnBodyA);
    return (r_Go_PA - r_Go_PB).norm();
}

/** Calculate the time rate of change of distance from a fixed point PB on 
body B to a fixed point PA on body A. We are given the location vectors 
p_Bo_PB and p_Ao_PA, expressed in their respective frames. We return 
d/dt |p_BoAo|, under the assumption that the time derivatives of the two given 
vectors in their own frames is zero. **/
Real calcStationToStationDistanceTimeDerivative
   (const State&         state,
    const Vec3&          locationOnBodyB,
    const MobilizedBody& bodyA,
    const Vec3&          locationOnBodyA) const
{
    if (isSameMobilizedBody(bodyA))
        return 0;

    Vec3 rB, rA, vB, vA;
    this->findStationLocationAndVelocityInGround(state,locationOnBodyB,rB,vB);
    bodyA.findStationLocationAndVelocityInGround(state,locationOnBodyA,rA,vA);
    const Vec3 r = rA-rB, v = vA-vB;
    const Real d = r.norm();

    // When the points are coincident, the rate of change of distance is just 
    // their relative speed. Otherwise, it is the speed along the direction of 
    // separation. 
    if (d==0) return v.norm();
    else return dot(v, r/d);
}


/** Calculate the second time derivative of distance from a fixed point PB on 
body B to a fixed point PA on body A. We are given the position vectors 
(stations) p_Bo_PB and p_Ao_PA, expressed in their respective frames. We return
d^2/dt^2 |p_PB_PA|, under the assumption that the time derivatives of the two 
given vectors in their own frames is zero. **/
Real calcStationToStationDistance2ndTimeDerivative
   (const State&         state,
    const Vec3&          locationOnBodyB,
    const MobilizedBody& bodyA,
    const Vec3&          locationOnBodyA) const
{
    if (isSameMobilizedBody(bodyA))
        return 0;

    Vec3 rB, rA, vB, vA, aB, aA;
    this->findStationLocationVelocityAndAccelerationInGround
                                            (state,locationOnBodyB,rB,vB,aB);
    bodyA.findStationLocationVelocityAndAccelerationInGround
                                            (state,locationOnBodyA,rA,vA,aA);

    const Vec3 r = rA-rB, v = vA-vB, a = aA-aB;
    const Real d = r.norm();
        
    // This method is the time derivative of 
    // calcFixedPointToPointDistanceTimeDerivative(), so it must follow the same
    // two cases. That is, when the points are coincident the change in 
    // separation rate is the time derivative of the speed |v|, otherwise it is 
    // the time derivative of the speed along the separation vector.

    if (d==0) {
        // Return d/dt |v|. This has two cases: if |v| is zero, the rate of 
        // change of speed is just the points' relative acceleration magnitude. 
        // Otherwise, it is the acceleration in the direction of the current 
        // relative velocity vector.
        const Real s = v.norm(); // speed
        if (s==0) return a.norm();
        else return dot(a, v/s);
    }

    // Points are separated.
    const Vec3 u = r/d; // u is separation direction (a unit vector from B to A) 
    const Vec3 vp = v - dot(v,u)*u; // velocity perp. to separation direction
    return dot(a,u) + dot(vp,v)/d;
}


/** TODO: not implemented yet -- any volunteers? Return the velocity of a 
point P moving on body B, in body A's frame, expressed in body A. **/
Vec3 calcBodyMovingPointVelocityInBody
   (const State&         state,
    const Vec3&          locationOnBodyB, 
    const Vec3&          velocityOnBodyB,
    const MobilizedBody& inBodyA) const
{
    SimTK_ASSERT_ALWAYS(!"unimplemented method", 
                        "MobilizedBody::calcBodyMovingPointVelocityInBody()"
                        " is not yet implemented -- any volunteers?");
    return Vec3::getNaN();
}


/** TODO: not implemented yet -- any volunteers? Return the velocity of a point
P moving (and possibly accelerating) on body B, in body A's frame, expressed in
body A. **/
Vec3 calcBodyMovingPointAccelerationInBody
   (const State&         state, 
    const Vec3&          locationOnBodyB, 
    const Vec3&          velocityOnBodyB, 
    const Vec3&          accelerationOnBodyB,
    const MobilizedBody& inBodyA) const
{
    SimTK_ASSERT_ALWAYS(!"unimplemented method", 
                        "MobilizedBody::calcBodyMovingPointAccelerationInBody()"
                        " is not yet implemented -- any volunteers?");
    return Vec3::getNaN();
}

/** TODO: not implemented yet -- any volunteers? Calculate the time rate of 
change of distance from a moving point PB on body B to a moving point PA on body
A. We are given the location vectors p_Bo_PB and p_Ao_PA, and the velocities of
PB in B and PA in A, all expressed in their respective frames. We return 
d/dt |p_BoAo|, taking into account the (given) time derivatives of the locations
in their local frames, as well as the relative velocities of the bodies. **/
Real calcMovingPointToPointDistanceTimeDerivative
   (const State&         state,
    const Vec3&          locationOnBodyB,
    const Vec3&          velocityOnBodyB,
    const MobilizedBody& bodyA,
    const Vec3&          locationOnBodyA,
    const Vec3&          velocityOnBodyA) const
{
    SimTK_ASSERT_ALWAYS(!"unimplemented method", 
                "MobilizedBody::calcMovingPointToPointDistanceTimeDerivative()"
                " is not yet implemented -- any volunteers?");
    return NaN;
}

/** TODO: not implemented yet -- any volunteers? Calculate the second time 
derivative of distance from a moving point PB on body B to a moving point PA on 
body A. We are given the location vectors p_Bo_PB and p_Ao_PA, and the 
velocities and accelerations of PB in B and PA in A, all expressed in their 
respective frames. We return d^2/dt^2 |p_AoBo|, taking into account the time 
derivatives of the locations in their local frames, as well as the relative 
velocities and accelerations of the bodies. **/
Real calcMovingPointToPointDistance2ndTimeDerivative
   (const State&         state,
    const Vec3&          locationOnBodyB,
    const Vec3&          velocityOnBodyB,
    const Vec3&          accelerationOnBodyB,
    const MobilizedBody& bodyA,
    const Vec3&          locationOnBodyA,
    const Vec3&          velocityOnBodyA,
    const Vec3&          accelerationOnBodyA) const
{
    SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::calcMovingPointToPointDistance2ndTimeDerivative()"
            " is not yet implemented -- any volunteers?");
    return NaN;
}


// End of High Level Operators.
/**@}**/


        //////////////////////////
        // CONSTRUCTION METHODS //
        //////////////////////////

/** The default constructor provides an empty %MobilizedBody handle that
can be assigned to reference any type of %MobilizedBody. **/
MobilizedBody() {}

/** Internal use only. **/
explicit MobilizedBody(MobilizedBodyImpl* r);

//------------------------------------------------------------------------------
/** @name                Construction and Misc Methods
These methods are the base class services which are used while building a 
concrete %MobilizedBody, or to query a %MobilizedBody to find out how it was 
built. These are unlikely to be used by end users of MobilizedBodies. **/
/**@{**/

/** Return a const reference to the Body contained within this %MobilizedBody. 
This refers to an internal copy of the Body that is owned by the 
%MobilizedBody. **/
const Body& getBody() const;

/** Return a writable reference to the Body contained within this 
%MobilizedBody. This refers to an internal copy of the Body that is owned by the
%MobilizedBody. Calling this method invalidates the %MobilizedBody's topology, 
so the containing System's realizeTopology() method must be called again. **/
Body& updBody();

/** Replace the Body contained within this %MobilizedBody with a new one. 
Calling this method invalidates the %MobilizedBody's topology, so the containing
System's realizeTopology() method must be called again. A reference to this 
%MobilizedBody is returned so that this can be chained like an assignment 
operator. **/
MobilizedBody& setBody(const Body&);

/** Add decorative geometry specified relative to the new (outboard) body's 
reference frame B. Note that the body itself may already have had some 
decorative geometry on it when it was first put into this MobilizedBody; this 
just adds more and the returned index is larger. Use the underlying Body 
object's accessors to find this decorative geometry again. The given 
\p geometry object is \e copied here; we do not keep a reference to the 
supplied object. **/
int addBodyDecoration(const Transform&          X_BD, 
                      const DecorativeGeometry& geometry) {
    return updBody().addDecoration(X_BD, geometry);
}
/** Convenience method for use when the \p geometry is supplied in the body 
frame. This is the same as addBodyDecoration(Transform(),geometry). **/
int addBodyDecoration(const DecorativeGeometry& geometry) {
    return updBody().addDecoration(geometry);
}

/** Add decorative geometry specified relative to the outboard mobilizer frame M
attached to body B. If body B already has decorative geometry on it, this just 
adds some more, but kept in a separate list from the body decorations and 
inboard decorations. Returns a unique index that can be used to identify this 
outboard decoration later (numbered starting from zero for outboard decorations 
only). **/
int addOutboardDecoration(const Transform&          X_MD, 
                          const DecorativeGeometry& geometry);
/** Return the count of decorations added with addOutboardDecoration(). **/
int getNumOutboardDecorations() const;
/** Return a const reference to the i'th outboard decoration. **/
const DecorativeGeometry& getOutboardDecoration(int i) const;
/** Return a writable reference to the i'th outboard decoration. **/
DecorativeGeometry& updOutboardDecoration(int i);

/** Add decorative geometry specified relative to the inboard mobilizer frame F 
attached to the parent body P. If body P already has decorative geometry on it, 
this just adds some more, but kept in a separate list from the body decorations 
and outboard decorations. Returns a unique index that can be used to identify 
this inboard decoration later (numbered starting from zero for inboard 
decorations only). **/
int addInboardDecoration(const Transform&          X_FD, 
                         const DecorativeGeometry& geometry);
/** Return the count of decorations added with addInboardDecoration(). **/
int getNumInboardDecorations() const;
/** Return a const reference to the i'th inboard decoration. **/
const DecorativeGeometry& getInboardDecoration(int i) const;
/** Return a writable reference to the i'th inboard decoration. **/
DecorativeGeometry& updInboardDecoration(int i);

/** If the contained Body can have its mass properties set to the supplied value 
\p m its mass properties are changed, otherwise the method fails. Calling this 
method invalidates the MobilizedBody's topology, so the containing matter 
subsystem's realizeTopology() method must be called again. A reference to this 
%MobilizedBody is returned so that this can be chained like an assignment
operator. **/
MobilizedBody& setDefaultMassProperties(const MassProperties& m) {
    updBody().setDefaultRigidBodyMassProperties(m); // might not be allowed
    return *this;
}

/** Return the mass properties of the Body stored within this MobilizedBody. **/
const MassProperties& getDefaultMassProperties() const {
     // every body type can do this
    return getBody().getDefaultRigidBodyMassProperties();
}

/** Provide a unique Motion object for this %MobilizedBody. The %MobilizedBody 
takes over ownership of the Motion object and is responsible for cleaning up 
its heap space when the time comes. This is a Topology-changing operation and
consequently requires write access to the %MobilizedBody which will propagate
to invalidate the containing Subsystem and System's topology. There can only
be one Motion object per mobilizer; this method will throw an exception if
there is already one here. **/
void adoptMotion(Motion& ownerHandle);

/** If there is a Motion object associated with this MobilizedBody it is 
removed; otherwise, nothing happens. If a Motion is deleted, the containing 
System's topology is invalidated. **/
void clearMotion();

/** Check whether this %MobilizedBody has an associated Motion object. This does
not tell you whether the Motion object is currently enabled or in use; just
whether it is available. **/
bool hasMotion() const;

/** If there is a Motion object associated with this MobilizedBody, this returns
a const reference to it. Otherwise it will throw an exception. You can check 
first using hasMotion(). Note that there is no provision to obtain a writable
reference to the contained Motion object; if you want to change it clear the
existing object instead and replace it with a new one.
@see hasMotion() **/
const Motion& getMotion() const;

/** Change this mobilizer's frame F on the parent body P. Calling this method
invalidates the %MobilizedBody's topology, so the containing matter subsystem's 
realizeTopology() method must be called again. A reference to this 
%MobilizedBody is returned so that this can be chained like an assignment
operator. **/
MobilizedBody& setDefaultInboardFrame (const Transform& X_PF);
/** Change this mobilizer's frame M fixed on this (the outboard) body B. Calling
this method invalidates the %MobilizedBody's topology, so the containing
matter subsystem's realizeTopology() method must be called again. A reference
to this %MobilizedBody is returned so that this can be chained like an 
assignment operator. **/
MobilizedBody& setDefaultOutboardFrame(const Transform& X_BM);

/** Return a reference to this mobilizer's default for the frame F fixed on the 
parent (inboard) body P, as the fixed Transform from P's body frame to the frame
F fixed to P. This default Transform is stored with the %MobilizedBody object, 
not the State. **/
const Transform& getDefaultInboardFrame()  const; // X_PF
/** Return a reference to this %MobilizedBody's default for mobilizer frame M, 
as the fixed Transform from this body B's frame to the frame M fixed on B. This 
default Transform is stored with the %MobilizedBody object, not the State. **/
const Transform& getDefaultOutboardFrame() const; // X_BM

/** This is an implicit conversion from %MobilizedBody to MobilizedBodyIndex 
when needed. This will fail unless this %MobilizedBody is owned by some 
SimbodyMatterSubsystem. We guarantee that the MobilizedBodyIndex of a mobilized 
body is numerically larger than the MobilizedBodyIndex of its parent. **/
operator MobilizedBodyIndex() const {return getMobilizedBodyIndex();}

/** Return the MobilizedBodyIndex of this MobilizedBody within the owning 
SimbodyMatterSubsystem. This will fail unless this MobilizedBody is owned by 
some SimbodyMatterSubsystem. We guarantee that the MobilizedBodyIndex of a 
mobilized body is numerically larger than the MobilizedBodyIndex of its 
parent. **/
MobilizedBodyIndex getMobilizedBodyIndex()  const;

/** Return a reference to the MobilizedBody serving as the parent body of the 
current MobilizedBody. This call will fail if the current MobilizedBody is 
Ground, since Ground has no parent. **/
const MobilizedBody& getParentMobilizedBody() const;

/** Return a reference to this %MobilizedBody's oldest ancestor other than 
Ground, or return Ground if this %MobilizedBody is Ground. That is, we return 
the "base" %MobilizedBody for this %MobilizedBody, meaning the one which 
connects this branch of the multibody tree directly to Ground. **/
const MobilizedBody& getBaseMobilizedBody()   const;

/** Obtain a reference to the SimbodyMatterSubsystem which contains this 
%MobilizedBody. This will fail unless this MobilizedBody is owned by some 
SimbodyMatterSubsystem. **/
const SimbodyMatterSubsystem& getMatterSubsystem() const;
/** Obtain a writable reference to the SimbodyMatterSubsystem which contains 
this %MobilizedBody. This will fail unless this %MobilizedBody is owned by some 
SimbodyMatterSubsystem. **/
SimbodyMatterSubsystem& updMatterSubsystem();

/** Determine whether the current MobilizedBody object is owned by a matter 
subsystem. **/
bool isInSubsystem() const;

/** Determine whether a given MobilizedBody \p mobod is in the same matter 
subsystem as the current body. If the bodies are not in a subsystem, this 
routine will return \c false. **/
bool isInSameSubsystem(const MobilizedBody& mobod) const;

/** Determine whether a given %MobilizedBody \p mobod is the same %MobilizedBody 
as this one. For this to be true the handles must not be empty, and the 
implementation objects must be <em>the same object</em> not separate objects 
with identical contents. **/
bool isSameMobilizedBody(const MobilizedBody& mobod) const;

/** Determine whether this body is Ground, meaning that it is actually 
body 0 of some matter subsytem, not just that its body type is Ground. **/
bool isGround() const;

/** Return this body's level in the tree of bodies, starting with ground at 0, 
bodies directly connected to ground at 1, bodies directly connected to those at 
2, etc. This is callable after realizeTopology(). This is the graph distance of 
the body from Ground. **/
int getLevelInMultibodyTree() const;
    
/** Create a new MobilizedBody which is identical to this one, except that it 
has a different parent (and consequently might belong to a different 
MultibodySystem). **/
MobilizedBody& cloneForNewParent(MobilizedBody& parent) const;


    // Utility operators //

/** This utility selects one of the q's (generalized coordinates) associated 
with this mobilizer from a supplied "q-like" Vector, meaning a Vector which is 
the same length as the Vector of q's for the containing matter subsystem. **/
Real getOneFromQPartition
   (const State& state, int which, const Vector& qlike) const;

/** This utility returns a writable reference to one of the q's (generalized 
coordinates) associated with this mobilizer from a supplied "q-like" Vector, 
meaning a Vector which is the same length as the Vector of q's for the 
containing matter subsystem. **/
Real& updOneFromQPartition
   (const State& state, int which, Vector& qlike) const;

/** This utility selects one of the u's (generalized speeds) associated with 
this mobilizer from a supplied "u-like" Vector, meaning a Vector which is 
the same length as the Vector of u's for the containing matter subsystem. **/
Real getOneFromUPartition
   (const State& state, int which, const Vector& ulike) const;

/** This utility returns a writable reference to one of the u's (generalized 
speeds) associated with this mobilizer from a supplied "u-like" Vector, meaning 
a Vector which is the same length as the Vector of u's for the containing matter
subsystem. **/
Real& updOneFromUPartition(const State& state, int which, Vector& ulike) const;

/** This utility adds in the supplied generalized force \p f (a scalar) to
the appropriate slot of the supplied \p mobilityForces Vector, which is a 
"u-like" Vector. Note that we are <em>adding</em> this not <em>setting</em> it 
so it important that \p mobilityForces be initialized to zero before making a 
set of calls to applyOneMobilityForce(). **/
void applyOneMobilityForce(const State& state, int which, Real f, 
                           Vector& mobilityForces) const
{
    updOneFromUPartition(state,which,mobilityForces) += f;
}

/** Given a generalized force in the q-space of this mobilizer, convert it to 
the equivalent generalized mobility force (u-space force). This uses the 
kinematic coupling matrix N that appears in equation (1) qdot=N*u. Here we 
compute (2) fu = ~N*fq (that's N transpose, not inverse).

Simbody deals with generalized forces in mobility (u) space, but sometimes these
are more convenient to generate in generalized coordinate (q) space. In that 
case this utility method is useful to perform the conversion from q space to u 
space that is necessary for communicating the force to Simbody.

@param[in]      state
    A State already realized through Position stage, from which this mobilizer's 
    kinematic coupling matrix N(q) is obtained.
@param[in]      fq
    This is a generalized force in the space of the generalized coordinates q 
    rather than the generalized speeds u. The length of \p fq must be nq, the 
    number of q's currently being used by this mobilizer in the given \p state. 
    (This can depend on a Model-stage state variable.)
@param[out]     fu
    This is the generalized force in mobility space (the space of the 
    generalized speeds u) that is equivalent to \p fq. \p fu will be resized if 
    necessary to length nu, the number of u's being used by this mobilizer.

<h3>Theory</h3>
The physical quantity power (force times velocity) must not change as a result 
of a change of coordinates. Hence we must have ~fq*qdot==~fu*u which follows 
from equations (1) and (2): multiply (1) by ~fq to get 
<pre>
    ~fq*qdot= ~fq*N*u
            = ~(~N*fq)*u
            = ~fu*u           from equation (2).
</pre>
For any mobilizer where qdot==u this simply copies the input to the output. 
Otherwise a multiplication by ~N is done, but that is very fast since N has 
already been computed. Cost depends on type of mobilizer but is unlikely to 
exceed 25 flops. **/
void convertQForceToUForce(const State&                        state,
                           const Array_<Real,MobilizerQIndex>& fq,
                           Array_<Real,MobilizerUIndex>&       fu) const;

/** This utility adds in the supplied spatial force \p spatialForceInG 
(consisting of a torque vector, and a force vector to be applied at the current 
body's origin) to the appropriate slot of the supplied \p bodyForcesInG Vector. 
Note that we are <em>adding</em> this not <em>setting</em> it so it is important 
that \p bodyForcesInG be initialized to zero before making a set of calls to 
applyBodyForce(). **/
void applyBodyForce(const State& state, const SpatialVec& spatialForceInG, 
                    Vector_<SpatialVec>& bodyForcesInG) const;

/** This utility adds in the supplied pure torque \p torqueInG to the 
appropriate slot of the supplied \p bodyForcesInG Vector. Note that we are 
<em>adding</em> this not <em>setting</em> it so it is important that 
\p bodyForcesInG be initialized to zero before making a set of calls to 
applyBodyTorque(). **/
void applyBodyTorque(const State& state, const Vec3& torqueInG, 
                     Vector_<SpatialVec>& bodyForcesInG) const;

/** This utility adds in the supplied force \p forceInG applied at a point 
\p pointInB to the appropriate slot of the supplied \p bodyForcesInG Vector. 
Notes: 
 - we are <em>adding</em> this not <em>setting</em> it so it is important that
   \p bodyForcesInG be initialized to zero before making a set of calls to 
   applyForceToBodyPoint().
 - \p pointInB represents a fixed station of B and is provided by giving the
   vector from body B's origin to the point, expressed in the B frame, while
   the applied force (and resulting body forces and torques) are expressed
   in the ground frame. **/
void applyForceToBodyPoint
   (const State& state, const Vec3& pointInB, const Vec3& forceInG,
    Vector_<SpatialVec>& bodyForcesInG) const;

// End of Construction and Misc Methods.
/**@}**/

    /////////////////////////////////////
    // BUILT IN MOBILIZER DECLARATIONS //
    /////////////////////////////////////

//------------------------------------------------------------------------------
// These are the built-in MobilizedBody types. Each of these has a known 
// number of coordinates and speeds (at least a default number) so
// can define routines which return and accept specific-size arguments, e.g.
// Real (for 1-dof mobilizer) and Vec5 (for 5-dof mobilizer). Here is the
// conventional interface that each built-in should provide. The base type
// provides similar routines but using variable-sized or "one at a time"
// arguments. (Vec<1> here will actually be a Real; assume the built-in
// MobilizedBody class is "BuiltIn")
//
//    BuiltIn&       setDefaultQ(const Vec<nq>&);
//    const Vec<nq>& getDefaultQ() const;
//
//    const Vec<nq>& getQ[Dot[Dot]](const State&) const;
//    const Vec<nu>& getU[Dot](const State&) const;
//
//    void setQ(State&, const Vec<nq>&) const;
//    void setU(State&, const Vec<nu>&) const;
//
//    const Vec<nq>& getMyPartQ(const State&, const Vector& qlike) const;
//    const Vec<nu>& getMyPartU(const State&, const Vector& ulike) const;
//   
//    Vec<nq>& updMyPartQ(const State&, Vector& qlike) const;
//    Vec<nu>& updMyPartU(const State&, Vector& ulike) const;
//
// Each built in mobilized body type is declared in its own header
// file using naming convention MobilizedBody_Pin.h, for example. All the
// built-in headers are collected in MobilizedBody_BuiltIns.h; you should
// include new ones there also.


class Pin;              
typedef Pin Torsion;  ///< Synonym for Pin mobilizer.
typedef Pin Revolute; ///< Synonym for Pin mobilizer.

class Universal;
class Cylinder;
class Weld;

class Slider;           
typedef Slider Prismatic; ///< Synonym for Slider mobilizer.

class Translation;      
typedef Translation Cartesian;       ///< Synonym for Translation mobilizer.
typedef Translation CartesianCoords; ///< Synonym for Translation mobilizer.

class BendStretch;      
typedef BendStretch PolarCoords; ///< Synonym for BendStretch mobilizer.

class SphericalCoords;
class LineOrientation;

class Planar;
class Gimbal;
class Bushing;

class Ball;             
typedef Ball Orientation;   ///< Synonym for Ball mobilizer.
typedef Ball Spherical;     ///< Synonym for Ball mobilizer.

class Free;
class FreeLine;
class Screw;
class Ellipsoid;
class Custom;
class Ground;
class FunctionBased;
    
// Internal use only.
class PinImpl;
class SliderImpl;
class UniversalImpl;
class CylinderImpl;
class BendStretchImpl;
class PlanarImpl;
class GimbalImpl;
class BushingImpl;
class BallImpl;
class TranslationImpl;
class SphericalCoordsImpl;
class FreeImpl;
class LineOrientationImpl;
class FreeLineImpl;
class WeldImpl;
class ScrewImpl;
class EllipsoidImpl;
class CustomImpl;
class GroundImpl;
class FunctionBasedImpl;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_H_



