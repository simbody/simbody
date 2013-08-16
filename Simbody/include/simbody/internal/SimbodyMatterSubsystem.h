#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MobilizedBody.h"

#include <cassert>
#include <vector>
#include <iostream>

class SimbodyMatterSubsystemRep;

namespace SimTK {

class MobilizedBody;
class MultibodySystem;
class Constraint;

/** This subsystem contains the bodies ("matter") in the multibody system,
the mobilizers (joints) that define the generalized coordinates used to 
represent the motion of those bodies, and constraints that must be satisfied
by the values of those coordinates.

There are many methods in the API for this class. For whole-system information
and calculations, the methods here are the right ones to use. For information
associated with individual objects contained in the subsystem, such as 
MobilizedBody and Constraint objects, it is generally easier to obtain the
information through the contained objects' APIs instead.

This class is is a "handle" containing only an opaque reference to the 
underlying implementation class.

<h3>Theory discussion</h3>
The bodies, mobilizers, and constraints are represented mathematically with
the following set of equations:
<pre>
                     qdot = N u                 Kinematic differential eqns.
                     zdot = zdot(t,q,u,z)       Auxiliary states

         M udot + ~G mult = f(t,q,u,z)          Equations of motion
         G udot           = b(t,q,u) 

                 where

          [P]    [bp]
        G=[V]  b=[bv]  f = T + ~J*(F-C)
          [A]    [ba]

          pdotdot = P udot - bp(t,q,u) = 0      Acceleration constraints
             vdot = V udot - bv(t,q,u) = 0
    a(t,q,u,udot) = A udot - ba(t,q,u) = 0
          
                   pdot = P u - c(t,q) = 0      Velocity constraints
                              v(t,q,u) = 0

                                p(t,q) = 0      Position constraints
                                  n(q) = 0      Normalization constraints
</pre>
where M(q) is the mass matrix, G(t,q,u) the acceleration constraint matrix, 
C(q,u) the coriolis and gyroscopic forces, T is user-applied joint mobility
forces, F is user-applied body forces and torques and gravity. J(q) is the
%System Jacobian (partial velocity matrix) whose transpose ~J maps spatial
forces to joint mobility forces. p(t,q) are the holonomic (position) 
constraints, v(t,q,u) the non-holonomic (velocity) constraints, and 
a(t,q,u,udot) the acceleration-only constraints, which must be linear in udot, 
with A(t,q,u) the coefficient matrix for a(). pdot, pdotdot are obtained by 
differentiation of p(), vdot by differentiation of v(). P(t,q)=Dpdot/Du 
(yes, that's u, not q -- we can get Pq=Dp/Dq when we need it as Pq=P*N^-1) and
V(t,q,u)=Dv/Du. (We use capital "D" to indicate partial derivative.) n(q) is 
the set of quaternion normalization constraints, which exist only at the 
position level and are uncoupled from everything else.

We calculate the constraint multipliers like this:
<pre>
          G M^-1 ~G mult = G udot0 - b
          where    udot0 = M^-1 f
</pre>
using the pseudo inverse of G M^-1 ~G to give a least squares solution for
mult: mult = pinv(G M^-1 ~G)(G M^-1 f - b). Then the real udot is
udot = udot0 - udotC, with udotC = M^-1 ~G mult. Note: M^-1* is an
O(n) operator that provides the desired result; it <em>does not</em> require
forming or factoring M.

NOTE: only the following constraint matrices have to be formed and factored:
<pre>
   [G M^-1 ~G]   to calculate multipliers

   [P N^-1]      for projection onto position manifold (a.k.a. Pq)

   [ P ]         for projection onto velocity manifold
   [ V ]  
</pre>

When working in a weighted norm with weights W on the state variables and
weights T (1/tolerance) on the constraint errors, the matrices we need are
actually [Tp Pq Wq^-1], [Tpv [P;V] Wu^-1], etc. with T and W diagonal
weighting matrices. These can then be used to find least squares solutions
in the weighted norms.

In many cases these matrices consist of decoupled blocks which can
be solved independently. (TODO: take advantage of that whenever possible
to solve a set of smaller systems rather than one large one.) Also, in the
majority of biosimulation applications we are likely to have only holonomic
(position) constraints, so there is no V or A and G=P is the whole story.
**/
class SimTK_SIMBODY_EXPORT SimbodyMatterSubsystem : public Subsystem {
public:

//==============================================================================
/** @name      Construction, Destruction, Topological information

Methods in this section are used in the extended construction phase for
a %SimbodyMatterSubsystem which we call defining the "topology" of the
multibody system. This includes adding mobilized bodies and constraints. 
Topological information is always state-independent since it is kept in
the %SimbodyMatterSubsystem object directly. The construction phase ends
when realizeTopology() is called on the containing System. **/
/**@{**/

/** Create a matter subsystem containing only the Ground body (mobilized 
body 0), and add the subsystem to the indicated MultibodySystem. The 
MultibodySystem takes over ownership of the subsystem, which is not 
copied. The MultibodySystem and this subsystem handle both refer to the
same subsystem after this call. **/
explicit SimbodyMatterSubsystem(MultibodySystem&);
/** Create an orphan matter subsystem containing only the Ground body 
(mobilized body 0); normally use the other constructor to place the 
subsystem in a MultibodySystem. **/
SimbodyMatterSubsystem();
/** The destructor destroys the subsystem implementation object only if this
handle is the last reference. Normally, there is a MultibodySystem that holds
a reference to the subsystem implementation, so this destruction will do 
nothing. **/
~SimbodyMatterSubsystem() {} // invokes ~Subsystem()

/** Given a MobilizedBodyIndex, return a read-only (const) reference to the 
corresponding MobilizedBody within this matter subsystem. This method will 
fail if the index is invalid or out of range. MobilizedBodyIndex(0) selects
the Ground mobilized body. **/
const MobilizedBody& getMobilizedBody(MobilizedBodyIndex) const;

/** Given a MobilizedBodyIndex, return a writable reference to the 
corresponding MobilizedBody within this matter subsystem. This method will 
fail if the index is invalid or out of range. MobilizedBodyIndex(0) selects
the Ground mobilized body. **/
MobilizedBody& updMobilizedBody(MobilizedBodyIndex);

/** Return a read-only (const) reference to the Ground MobilizedBody
within this matter subsystem. **/
const MobilizedBody::Ground& getGround() const;
/** Return a writable reference to the Ground MobilizedBody within this
matter subsystem; you need a writable reference if you're adding a
mobilized body that is directly connected to Ground. **/
MobilizedBody::Ground& updGround();
/** This is a synonym for updGround() that makes for nicer-looking examples.
Note: topology is not marked invalid upon returning a writable reference
here; that will be done only if a non-const method of the returned 
MobilizedBody is called. That means it is OK to use Ground() to satisfy 
a const argument; it won't have an "invalidate topology" side effect.
@see updGround() **/
MobilizedBody::Ground& Ground() {return updGround();}


/** Given a ConstraintIndex, return a read-only (const) reference to the 
corresponding Constraint within this matter subsystem. This method will 
fail if the index is invalid or out of range. **/
const Constraint& getConstraint(ConstraintIndex) const;

/** Given a ConstraintIndex, return a writable reference to the corresponding 
Constraint within this matter subsystem. This method will 
fail if the index is invalid or out of range. **/
Constraint& updConstraint(ConstraintIndex);

/** Normally the matter subsystem will attempt to generate some decorative
geometry as a sketch of the defined multibody system; you can disable that
with this method. **/
void setShowDefaultGeometry(bool show);
/** Get whether this matter subsystem is set to generate default decorative 
geometry that can be used to visualize this multibody system. **/
bool getShowDefaultGeometry() const;

/** The number of bodies includes all mobilized bodies \e including Ground,
which is the 0th mobilized body. (Note: if special particle handling were
implmemented, the count here would \e not include particles.) Bodies and their
inboard mobilizers have the same index since they are grouped together as a 
MobilizedBody. MobilizedBody numbering (using unique integer type
MobilizedBodyIndex) starts with Ground at MobilizedBodyIndex(0) with a regular
labeling such that children have higher indices than their parents. Ground
does not have a mobilizer (or I suppose you could think of its mobilizer as 
the Weld joint that attaches it to the universe), but otherwise 
every mobilized body has a unique body and mobilizer. **/
int getNumBodies() const;

/** This is the total number of defined constraints, each of which may
generate more than one constraint equation. This is the number of Constraint
objects that were defined; in a given State some of these may be disabled. **/
int getNumConstraints() const;

/** The sum of all the mobilizer degrees of freedom. This is also the length
of the state variable vector u and the mobility forces array. **/
int getNumMobilities() const;

/** The sum of all the q vector allocations for each joint. There may be
some that are not in use for particular modeling options. **/
int getTotalQAlloc() const;

/** This is the sum of all the allocations for constraint multipliers, one per
acceleration constraint equation. There may be some that are not in use due
to the corresonding Constraint elements being disabled in a given State. **/
int getTotalMultAlloc() const;


/** Attach new matter by attaching it to the indicated parent body (not
normally called by users -- see MobilizedBody). The mobilizer and mass 
properties are provided by \a child. A new MobilizedBodyIndex is assigned for
the child; it is guaranteed to be numerically larger than the 
MobilizedBodyIndex of the parent. We take over ownership of \a child's 
implementation object from the given MobilizedBody handle, leaving that handle
as a reference to the implementation object now owned by the matter subsystem. 
It is an error if the given MobilizedBody handle wasn't the owner of the 
implementation object to which it refers.
@note
This method is usually called by concrete MobilizedBody constructors;
it does not normally need to be called by end users. **/
MobilizedBodyIndex   adoptMobilizedBody(MobilizedBodyIndex  parent, 
                                        MobilizedBody&      child);

/** Add a new Constraint object to the matter subsystem (not normally called
by users -- see Constraint). The details of the Constraint are opaque here. 
A new ConstraintIndex is assigned. We take  over ownership of the 
implementation object from the given Constraint handle, leaving that handle as
a reference to the implementation object now owned by the matter subsystem. It
is an error if the given Constraint handle wasn't the owner of the 
implementation object to which it refers.
@note
This method is usually called by concrete Constraint constructors; it does not
normally need to be called by end users. **/
ConstraintIndex   adoptConstraint(Constraint&);

/** Copy constructor is not very useful. **/
SimbodyMatterSubsystem(const SimbodyMatterSubsystem& ss) : Subsystem(ss) {}
/** Copy assignment is not very useful. **/
SimbodyMatterSubsystem& operator=(const SimbodyMatterSubsystem& ss) 
{   Subsystem::operator=(ss); return *this; }


//==============================================================================
/** @name                  Set/get modeling options

Methods in this section involve setting and getting various modeling options
that may be selected. This includes whether to use quaternions or Euler angles
to represent rotations, and enabling/disabling constraints. **/

/**@{**/
/** For all mobilizers offering unrestricted orientation, decide what
method we should use to model their orientations. Choices are: 
quaternions (best for dynamics), or rotation angles (1-2-3 Euler 
sequence, good for optimization). Changing this flag invalidates Model
stage and above in the supplied \a state, leaving it realized only through
Topology stage, so you must call realizeModel() on the containing 
MultibodySystem prior to using this \a state in further calculations. **/
void setUseEulerAngles(State& state, bool useEulerAngles) const;

/** Return the current setting of the "use Euler angles" model variable as
set in the supplied \a state. **/
bool getUseEulerAngles  (const State& state) const;

/** Return the number of quaternions in use by the mobilizers of this system, 
given the current setting of the "use Euler angles" flag in the supplied
\a state, and the types of mobilizers in the multibody tree. 
@see isUsingQuaternion(), getQuaternionPoolIndex() **/
int  getNumQuaternionsInUse(const State& state) const;

/** Check whether a given mobilizer is currently using quaternions, based
on the type of mobilizer and the setting of the "use Euler angles" flag in
the supplied \a state. 
@see getNumQuaternionsInUse(), getQuaternionPoolIndex() **/
bool isUsingQuaternion(const State& state, MobilizedBodyIndex mobodIx) const;
/** If the given mobilizer is currently using a quaternion to represent
orientation, return the QuaternionPoolIndex (a small integer) assigned to that
quaternion. This is used, for example, to find which normalization constraint
error is associated with which quaternion. 
@see isUsingQuaternion(), getNumQuaternionsInUse() **/
QuaternionPoolIndex getQuaternionPoolIndex(const State& state, 
                                           MobilizedBodyIndex mobodIx) const;

/** Disable or enable the Constraint whose ConstraintIndex is supplied within
the supplied \a state. Whether a Constraint is disabled is an Instance-stage
state variable so enabling or disabling invalidates Instance stage and higher
in the given \a state, leaving the \a state realized no higher than Model
stage.
@see isConstraintDisabled() **/
void setConstraintIsDisabled(State&          state, 
                             ConstraintIndex constraintIx, 
                             bool            shouldDisableConstraint) const;

/** Determine whether a particular Constraint is currently disabled in the
given \a state. 
@see setConstraintIsDisabled() **/
bool isConstraintDisabled(const State&, ConstraintIndex constraint) const;

/** Given a State which may be modeled using quaternions, copy it to another
State which represents the same configuration using Euler angles instead. If
the \a inputState already uses Euler angles, the output will just be a
duplicate. All continuous and discrete State variables will be copied to the
\a outputState but they will not necessarily have been realized to the same
level as the \a inputState. **/
void convertToEulerAngles(const State& inputState, State& outputState) const;

/** Given a State which may be modeled using Euler angles, copy it to another
State which represents the same configuration using quaternions instead. If
the \a inputState already uses quaternions, the output will just be a
duplicate. All continuous and discrete State variables will be copied to the
\a outputState but they will not necessarily have been realized to the same
level as the \a inputState. **/
void convertToQuaternions(const State& inputState, State& outputState) const; 

/**@}**/


//==============================================================================
/** @name               Calculate whole-system properties

These methods perform calculations that yield properties of the system as
a whole. These are \e operators, meaning that they make use of the supplied
State but do not modify the State. They simply calculate a result and return
it to you without storing it internally. Each method requires that the
State has already been realized to at least a particular stage which is 
documented with the method. **/

/**@{**/
/** Calculate the total system mass.
@par Required stage
  \c Stage::Instance **/
Real calcSystemMass(const State& s) const;

/** Return the position vector p_GC of the system mass center C, measured from 
the Ground origin, and expressed in Ground. 
@par Required stage
  \c Stage::Position **/
Vec3 calcSystemMassCenterLocationInGround(const State& s) const;

/** Return total system mass, mass center location measured from the Ground 
origin, and system inertia taken about the Ground origin, expressed in Ground.
@par Required stage
  \c Stage::Position **/
MassProperties calcSystemMassPropertiesInGround(const State& s) const;

/** Return the system inertia matrix taken about the system center of mass,
expressed in Ground.
@par Required stage
  \c Stage::Position **/
Inertia calcSystemCentralInertiaInGround(const State& s) const;

/** Return the velocity v_GC = d/dt p_GC of the system mass center C in the
Ground frame G, measured from Ground origin and expressed in G.
@par Required stage
  \c Stage::Velocity **/
Vec3 calcSystemMassCenterVelocityInGround(const State& s) const;

/** Return the acceleration a_GC = d/dt p_GC of the system mass center C in the
Ground frame G, measured from Ground origin and expressed in G.
@par Required stage
  \c Stage::Acceleration **/
Vec3 calcSystemMassCenterAccelerationInGround(const State& s) const;

/** Return the momentum of the system as a whole (angular, linear) measured
in the Ground frame, taken about the Ground origin and expressed in Ground.
(The linear component is independent of the "about" point.)
@see calcSystemCentralMomentum()
@par Required stage
  \c Stage::Velocity **/
SpatialVec calcSystemMomentumAboutGroundOrigin(const State& s) const;

/** Return the momentum of the system as a whole (angular, linear) measured
in the Ground frame, taken about the current system center of mass
location C and expressed in Ground. (The linear component is independent of the
"about" point.)
@see calcSystemMomentumAboutGroundOrigin()
@par Required stage
  \c Stage::Velocity **/
SpatialVec calcSystemCentralMomentum(const State& s) const;

/** Calculate the total kinetic energy of all the mobilized bodies in this
matter subsystem, given the configuration and velocities in \a state.
@par Required stage
  \c Stage::Velocity **/
Real calcKineticEnergy(const State& state) const;
/**@}**/

//==============================================================================
/** @name         System and Task Space Kinematic Jacobian Operators 

The system kinematic Jacobian maps between mobility space (generalized speeds 
and generalized forces) and Cartesian body space (mobilized body frame spatial 
velocities and spatial forces). A task space Jacobian maps between mobility 
space and a specified set of task points or frames fixed to a subset of the 
bodies, and generally located away from the body frame. A task space Jacobian
J can be used to construct various task space matrices such as the task space 
compliance matrix J M^-1 ~J or its inverse, the task space (or operational
space) inertia matrix. 

The system Jacobian J(q) maps n generalized speeds u to spatial velocities V of 
each of the nb mobilized bodies (including Ground), measured at the body frame 
origin relative to Ground, and expressed in the Ground frame. The transpose ~J 
of this matrix maps nb spatial forces to n generalized forces, where the spatial
forces are applied at the body frame origin and expressed in Ground. Similarly, 
task space Jacobians map from n generalized speeds to nt task frame spatial 
velocities (expressed in Ground), and transposed task space Jacobians map 
between task frame spatial forces (or impulses), expressed in Ground, and 
generalized forces (or generalized impulses).

Simbody provides fast O(n) methods ("operators") that can form matrix-vector
products like J*u or ~J*F without forming J. The "bias" term Jdot*u (also known
as the Coriolis acceleration) is also available; this arises when working at
the acceleration level because d/dt J*u = J*udot+Jdot*u (where dot means time
derivative). The computational cost of these operators is O(n+nt) so it is 
\e much more efficient to work with a group of tasks simultaneously than to 
process one at a time, which would have complexity O(n*nt). Alternatively, we 
provide methods that will return all or part of J explicitly; in general it 
is \e much more efficient computationally to work with the O(n) matrix-vector 
multiply operators rather than to form explicit matrices and then perform O(n^2) 
matrix-vector products. Performance estimates are given with each method so that
you can determine which methods to use. If you can, you should use the O(n) 
methods -- it is a good habit to get into when using an O(n) multibody code like
Simbody!

Note that the Jacobian is associated with an expressed-in frame for the
velocity or force vector and a designated station (point) on each body. We 
always use the Ground frame for Jacobians. For the system Jacobian, the body 
origin is always the designated station; for task Jacobians different stations
may be specified. We provide three different sets of methods for working with
    - the full %System Jacobian: J, nb X n 6-vectors (or 6*nb X n scalars)
    - the Station Jacobian for a set of nt task stations (points): JS, nt rows 
      of n 3-vectors (or a 3*nt X n Matrix of scalars)
    - the Frame Jacobian for a set of nt task frames fixed to a body: JF, nt 
      rows of n 6-vectors (or a 6*nt X n Matrix of scalars)

The rotational part of a Jacobian is the same for any frame fixed to the same 
body. So for Frame Jacobians you need specify only a station on the body (the 
frame's origin point). That means if you want a 3*nt X n Orientation Jacobian, 
you can obtain it from alternate rows of a Frame Jacobian. Using the above 
terminology, the complete %System Jacobian is a Frame Jacobian for which the 
task frames are the body frames, with each MobilizedBody appearing only once 
and in order of MobilizedBodyIndex (starting with Ground). 

It is acceptable for the same body to appear more than once in a list of tasks;
these are likely to conflict but that can be dealt with elsewhere. **/

/**@{**/

/** Calculate the product of the %System kinematic Jacobian J (also known as the 
partial velocity matrix) and a mobility-space vector u in O(n) time. If the 
vector u is a set of generalized speeds, then this produces the body spatial 
velocities that result from those generalized speeds. That is, the result is 
V_GB = J*u where V_GB[i] is the spatial velocity of the i'th body's body frame 
origin (in Ground) that results from the given set of generalized speeds. 

@param[in]      state
    A State compatible with this System that has already been realized to
    Stage::Position.
@param[in]      u
    A mobility-space Vector, such as a set of generalized speeds. The length
    and order must match the mobilities of this system (that is n, the number
    of generalized speeds u, \e not nq, the number of generalized 
    coordinates q).
@param[out]     Ju
    This is the product V=J*u as described above. Each element is a spatial
    vector, one per mobilized body, to be indexed by MobilizedBodyIndex.
    If the input vector is a set of generalized speeds u, then the results
    are nb spatial velocities V_GBi (that is, a pair of vectors w_GBi and v_GBi 
    giving angular and linear velocity). Note that Ground is body 0 so the 0th 
    element V_GB0=V_GG=Ju[0] is always zero on return.

The kinematic Jacobian (partial velocity matrix) J is defined as follows:
<pre>
      partial(V)                                 T                        T
  J = ----------, V = [V_GB0 V_GB1 ... V_GB nb-1] ,  u = [u0 u1 ... u n-1]
      partial(u)
</pre>
Thus the element J(i,j)=partial(V_GBi)/partial(uj) (each element of J is a
spatial vector). The transpose of this matrix maps spatial forces to 
generalized forces; see multiplyBySystemJacobianTranspose().

Note that we're using "monogram" notation for the spatial velocities, where
<pre>
            G Bi
    V_GBi =  V
</pre>
the spatial velocity of body i's body frame Bi (at its origin), measured and
expressed in the Ground frame G.

<h3>Performance discussion</h3>
This is a very fast operator, costing about 12*(nb+n) flops, where nb is the
number of bodies and n the number of mobilities (degrees of freedom) u. In 
contrast, even if you have already calculated the entire nbXnX6 matrix J, the 
multiplication J*u would cost 12*nb*n flops. As an example, for a 20 body 
system with a free flying base and 19 pin joints (25 dofs altogether), this 
method takes 12*(20+25)=540 flops while the explicit matrix-vector multiply 
would take 12*20*25=6000 flops. So this method is already >10X faster for 
that small system; for larger systems the difference grows rapidly.

@see multiplyBySystemJacobianTranspose(), calcSystemJacobian() **/
void multiplyBySystemJacobian( const State&         state,
                               const Vector&        u,
                               Vector_<SpatialVec>& Ju) const;

/** Calculate the acceleration bias term for the %System Jacobian, that is, the
part of the acceleration that is due only to velocities. This term is also
known as the Coriolis acceleration, and it is returned here as a spatial
acceleration of each body frame in Ground.

@param[in]      state
    A State that has already been realized through Velocity stage.
@param[out]     JDotu
    The product JDot*u where JDot = d/dt J, and u is the vector of generalized
    speeds taken from \a state. This is a Vector of nb SpatialVec elements.

<h3>Theory</h3>
The spatial velocity V_GBi of each body i can be obtained from the generalized
speeds u by V = {V_GBi} = J*u. Taking the time derivative in G gives
<pre>
    A = d/dt V = {A_GBi} = J*udot + JDot*u
</pre>
where JDot=JDot(q,u). This method returns JDot*u, which depends only on 
configuration q and speeds u. Note that the same u is used to calculate JDot, 
which is linear in u, so this term is quadratic in u.

<h3>Implementation</h3>
This method simply extracts the total Coriolis acceleration for each body that
is already available in the \a state cache so there is no computation done
here.
@see getTotalCoriolisAcceleration()
**/
void calcBiasForSystemJacobian(const State&         state,
                               Vector_<SpatialVec>& JDotu) const;


/** Alternate signature that returns the bias as a 6*nb-vector of scalars 
rather than as an nb-vector of 2x3 spatial vectors. See the other signature for
documentation. **/
void calcBiasForSystemJacobian(const State&         state,
                               Vector&              JDotu) const;

/** Calculate the product of the transposed kinematic Jacobian ~J (==J^T) and
a vector F_G of spatial force-like elements, one per body, in O(n) time to 
produce a generalized force-like result f=~J*F. If F_G is actually a set of
spatial forces applied at the body frame origin of each body, and expressed
in the Ground frame, then the result is the equivalent set of generalized
forces f that would produce the same accelerations as F_G.

@param[in]      state
    A State compatible with this System that has already been realized to
    Stage::Position.
@param[in]      F_G
    This is a vector of SpatialVec elements, one per mobilized body and in
    order of MobilizedBodyIndex (with the 0th entry a force on Ground; hence
    ignored). Each SpatialVec is a spatial force-like pair of 3-vectors 
    (torque,force) with the force applied at the body origin and the vectors
    expressed in Ground.
@param[out]     f
    This is the product f=~J*F_G as described above. This result is in the
    generalized force space, that is, it has one scalar entry for each of the
    n system mobilities (velocity degrees of freedom). Resized if necessary.

The kinematic Jacobian (partial velocity matrix) J is defined as follows:
<pre>
      partial(V)                                 T                        T
  J = ----------, V = [V_GB0 V_GB1 ... V_GB nb-1] ,  u = [u0 u1 ... u n-1]
      partial(u)
</pre>
Thus the element J(i,j)=partial(V_GBi)/partial(uj) (each element of J is a
spatial vector). J maps generalized speeds to spatial velocities (see
multiplyBySystemJacobian()); its transpose ~J maps spatial forces 
to generalized forces.

Note that we're using "monogram" notation for the spatial velocities, where
<pre>
            G Bi
    V_GBi =  V
</pre>
the spatial velocity of body i's body frame Bi (at its origin), measured and
expressed in the Ground frame G.

<h3>Performance discussion</h3>
This is a very fast operator, costing about 18*nb+11*n flops, where nb is the
number of bodies and n the number of mobilities (degrees of freedom) u. In 
contrast, even if you have already calculated the entire 6*nbXnu matrix J, the
multiplication ~J*F would cost 12*nb*n flops. As an example, for a 20 body 
system with a free flying base and 19 pin joints (25 dofs altogether), this 
method takes 18*20+11*25=635 flops while the explicit matrix-vector multiply 
would take 12*20*25=6000 flops. So this method is already >9X faster for 
that small system; for larger systems the difference grows rapidly. 

@see multiplyBySystemJacobian(), calcSystemJacobian() **/
void multiplyBySystemJacobianTranspose( const State&                state,
                                        const Vector_<SpatialVec>&  F_G,
                                        Vector&                     f) const;


/** Explicitly calculate and return the nb x nu whole-system kinematic 
Jacobian J_G, with each element a 2x3 spatial vector (SpatialVec). This matrix 
maps generalized speeds to the spatial velocities of all the bodies, which 
will be at the body origins, measured and expressed 
in Ground. That is, if you have a set of n generalized speeds u, you can 
find the spatial velocities of all nb bodies as V_G = J_G*u. The transpose of 
this matrix maps a set of spatial forces F_G, applied at the body frame 
origins and expressed in Ground, to the equivalent set of n generalized 
forces f: f = ~J_G*F_G. 

@note The 0th row of the returned Jacobian is always zero since it represents
the spatial velocity of Ground.

<h3>Performance discussion</h3>
Before using this method, consider whether you really need to form this
very large matrix which necessarily will take O(n^2) space and time; it will 
almost always be \e much faster to use the multiplyBySystemJacobian() method 
that directly calculate the matrix-vector product in O(n) time without explictly 
forming the matrix. Here are the details:

As currently implemented, forming the full Jacobian J costs about
12*n*(nb+n) flops. Assuming nb ~= n, this is about 24*n^2 flops. Then
if you want to form a product J*u explicitly, the matrix-vector multiply will 
cost about 12*n^2 flops each time you do it. In contrast the J*u product is 
calculated using multiplyBySystemJacobian() in about 24*n flops. Even for
very small systems it is cheaper to make repeated calls to 
multiplyBySystemJacobian() than to form J explicitly and multiply by it.
See the Performance section for multiplyBySystemJacobian() for more
comparisons.

@see multiplyBySystemJacobian(), multiplyBySystemJacobianTranspose()
@see calcSystemJacobian() alternate signature using scalar elements **/
void calcSystemJacobian(const State&            state,
                        Matrix_<SpatialVec>&    J_G) const; // nb X nu

/** Alternate signature that returns a system Jacobian as a 6*nb X n Matrix 
of scalars rather than as an nb X n matrix of 2x3 spatial vectors. See
the other signature for documentation and important performance 
considerations. **/
void calcSystemJacobian(const State&            state,
                        Matrix&                 J_G) const; // 6 nb X nu


/** Calculate the Cartesian ground-frame velocities of a set of task stations 
(points fixed on bodies) that results from a particular set of generalized 
speeds u. The result is the station velocities measured and expressed in Ground.

@param[in]      state
    A State that has already been realized through Position stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the stations of 
    interest are fixed.
@param[in]      stationPInB
    The array of nt station points P of interest (one per task), each
    corresponding to one of the bodies B from \a onBodyB, given as vectors 
    from each body B's origin Bo to its station P, expressed in frame B.
@param[in]      u
    A mobility-space Vector, such as a set of generalized speeds. The length
    and order must match the mobilities of this system (that is n, the number
    of generalized speeds u, \e not nq, the number of generalized 
    coordinates q).
@param[out]     JSu
    The resulting product JS*u, where JS is the station task Jacobian. Resized
    to nt if needed.

<h3>Performance discussion</h3>
It is almost always better to use this method than to form an explicit 3*nt X n 
station task Jacobian explicitly and then multiply by it. If you have only one 
or two tasks, so that the matrix is only 3xn or 6xn, and then perform many 
multiplies with that matrix, it might be slightly cheaper to form it. For 
example, it is about 4X cheaper to use this method than to form a one-task 
Station Jacobian JS explicitly and use it once. However, because this would be 
such a skinny matrix (3 X n) explicit multiplication is cheap so if you will 
re-use this same Jacobian repeatedly before recalculating (at least 6 times) 
then it may be worth calculating and saving it. Here are the details:

A call to this method costs 27*nt + 12*(nb+n) flops. If you assume that 
nb ~= n >> 1, you could say this is about 27*nt + 24*n flops. In
contrast, assuming you already have the 3*nt X n station Jacobian JS available,
you can compute the JS*u product in about 6*nt*n flops, 3X faster for one task,
about even for three tasks, and slower for more than three tasks.
However forming JS costs about 40*nt+90*n flops (see calcStationJacobian()).
So to form a one-task Jacobian and use it once is 4X more expensive (96*n vs
24*n), but if you use it more than 5 times it is cheaper to do it
explicitly. Forming a one-task JS and using it 100 times costs about 690*n 
flops while calling this method 100 times would cost about 2400*n flops.

@see multiplyByStationJacobianTranspose(), calcStationJacobian() **/
void multiplyByStationJacobian(const State&                      state,
                               const Array_<MobilizedBodyIndex>& onBodyB,
                               const Array_<Vec3>&               stationPInB,
                               const Vector&                     u,
                               Vector_<Vec3>&                    JSu) const;

/** Alternate signature for when you just have a single station task. 
@returns JS*u, where JS is the station task Jacobian. **/
Vec3 multiplyByStationJacobian(const State&         state,
                               MobilizedBodyIndex   onBodyB,
                               const Vec3&          stationPInB,
                               const Vector&        u) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&stationPInB, &stationPInB+1);
    Vector_<Vec3>                       JSu(1);
    multiplyByStationJacobian(state, bodies, stations, u, JSu);
    return JSu[0];
}


/** Calculate the generalized forces resulting from a single force applied
to a set of nt station tasks (points fixed to bodies) P. The applied forces 
f_GP should be 3-vectors expressed in Ground. This is considerably faster than 
forming the Jacobian explicitly and then performing the matrix-vector multiply.

<h3>Performance discussion</h3>
Cost is about 30*nt + 18*nb + 11*n. Assuming nb ~= n, this is roughly
30*(n+nt). In contrast, forming the complete 3*nt X n matrix would cost about
90*(n+nt/2), and subsequent explicit matrix-vector multiplies would cost
about 6*nt*n each.

@see multiplyByStationJacobian(), calcStationJacobian() **/
void multiplyByStationJacobianTranspose
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 stationPInB,
    const Vector_<Vec3>&                f_GP,
    Vector&                             f) const;

/** Alternate signature for when you just have a single station task. **/
void multiplyByStationJacobianTranspose
   (const State&                        state,
    MobilizedBodyIndex                  onBodyB,
    const Vec3&                         stationPInB,
    const Vec3&                         f_GP,
    Vector&                             f) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&stationPInB, &stationPInB+1);
    Vector_<Vec3>                       forces(1, f_GP);
    multiplyByStationJacobianTranspose(state, bodies, stations, forces, f);
}

/** Explicitly calculate and return the 3*nt x n kinematic Jacobian JS for a 
set of nt station tasks P (a station is a point fixed on a particular mobilized 
body). This matrix maps generalized speeds to the Cartesian velocity of each
station, measured and expressed in Ground. That is, if you have a set of n 
generalized speeds u, you can find the Cartesian velocities of stations P as 
v_GP = JS*u, where v_GP is a 3*nt column vector. The transpose of this 
matrix maps a 3*nt vector of forces f_GP (expressed in Ground and applied 
to P) to the equivalent set of n generalized forces f: f = ~JS*f_GP.

@note It is almost always far more efficient to use multiplyByStationJacobian()
or multiplyByStationJacobianTranspose() to form matrix-vector products rather 
than to use this method to form the Jacobian explicitly. See the performance 
discussions there.

Overloaded signatures of this method are available to allow you to obtain the
Jacobian either as an nt X n Matrix with Vec3 elements, or as 3*nt X n Matrix
with scalar elements.

@param[in]      state
    A State that has already been realized through Position stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the stations of 
    interest are fixed.
@param[in]      stationPInB
    The array of nt station points P of interest (one per task), each
    corresponding to one of the bodies B from \a onBodyB, given as vectors 
    from each body B's origin Bo to its station P, expressed in frame B.
@param[out]     JS
    The resulting nt X n station task Jacobian. Resized if necessary.

<h3>Performance discussion</h3>
The cost of a call to this method is about 42*nt + 54*nb + 33*n flops. If we 
assume that nb ~= n >> 1, this is roughly 90*(n+nt/2) flops. Then once the 
Station Jacobian JS has been formed, each JS*u matrix-vector product costs
6*nt*n flops to form. When nt is small enough (say one or two tasks), and you
plan to re-use it a lot, this can be computationally efficient; but for single
use or more than a few tasks you can do much better with 
multiplyByStationJacobian() or multiplyByStationJacobianTranspose().

@see multiplyByStationJacobian(), multiplyByStationJacobianTranspose() **/
void calcStationJacobian(const State&                        state,
                         const Array_<MobilizedBodyIndex>&   onBodyB,
                         const Array_<Vec3>&                 stationPInB,
                         Matrix_<Vec3>&                      JS) const;

/** Alternate signature for when you just have a single station task. **/
void calcStationJacobian(const State&       state,
                         MobilizedBodyIndex onBodyB,
                         const Vec3&        stationPInB,
                         RowVector_<Vec3>&  JS) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&stationPInB, &stationPInB+1);
    calcStationJacobian(state, bodies, stations, JS);
}

/** Alternate signature that returns a station Jacobian as a 3*nt x n Matrix 
rather than as a Matrix of Vec3 elements. See the other signature for 
documentation and important performance considerations. **/
void calcStationJacobian(const State&                        state,
                         const Array_<MobilizedBodyIndex>&   onBodyB,
                         const Array_<Vec3>&                 stationPInB,
                         Matrix&                             JS) const;

/** Alternate signature for when you just have a single station task. **/
void calcStationJacobian(const State&       state,
                         MobilizedBodyIndex onBodyB,
                         const Vec3&        stationPInB,
                         Matrix&            JS) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&stationPInB, &stationPInB+1);
    calcStationJacobian(state, bodies, stations, JS);
}


/** Calculate the acceleration bias term for a station Jacobian, that is, the
part of the station's acceleration that is due only to velocities. This term 
is also known as the Coriolis acceleration, and it is returned here as a linear
acceleration of the station in Ground.

@param[in]      state
    A State that has already been realized through Velocity stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the stations of 
    interest are fixed.
@param[in]      stationPInB
    The array of nt station points P of interest (one per task), each
    corresponding to one of the bodies B from \a onBodyB, given as vectors 
    from each body B's origin Bo to its station P, expressed in frame B.
@param[out]     JSDotu
    The resulting product JSDot*u, where JSDot is the time derivative of JS,
    the station task Jacobian. Resized to nt if needed.

<h3>Theory</h3>
The velocity v_GP of a station point P in the Ground frame G can be obtained 
from the generalized speeds u using the station Jacobian for P, as <pre>
    v_GP = JS_P*u
</pre> Taking the time derivative in G gives <pre>
    a_GP = JS_P*udot + JSDot_P*u
</pre>
This method returns JSDot_P*u, which depends only on configuration and 
velocities. We allow for a set of task points P so that all their bias terms
can be calculated in a single sweep of the multibody tree. Note that u is taken
from the \a state and that the same u shown above is also used to calculate 
JSDot_P, which is linear in u, so the bias term is quadratic in u.

<h3>Implementation</h3>
This method just obtains body B's total Coriolis acceleration already available
in the \a state cache and shifts it to station point P. Cost is 48*nt flops.
@see getTotalCoriolisAcceleration(), shiftAccelerationBy()
**/
void calcBiasForStationJacobian(const State&                      state,
                                const Array_<MobilizedBodyIndex>& onBodyB,
                                const Array_<Vec3>&               stationPInB,
                                Vector_<Vec3>&                    JSDotu) const;

/** Alternate signature that returns the bias as a 3*nt-vector of scalars 
rather than as an nt-vector of Vec3s. See the other signature for
documentation. **/
void calcBiasForStationJacobian(const State&                      state,
                                const Array_<MobilizedBodyIndex>& onBodyB,
                                const Array_<Vec3>&               stationPInB,
                                Vector&                           JSDotu) const;

/** Alternate signature for when you just have a single station task. 
@returns JSDot*u, where JSDot is the station Jacobian time derivative. **/
Vec3 calcBiasForStationJacobian(const State&         state,
                                MobilizedBodyIndex   onBodyB,
                                const Vec3&          stationPInB) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&stationPInB, &stationPInB+1);
    Vector_<Vec3>                       JSDotu(1);
    calcBiasForStationJacobian(state, bodies, stations, JSDotu);
    return JSDotu[0];
}


/** Calculate the spatial velocities of a set of nt task frames A={Ai} fixed to 
nt bodies B={Bi}, that result from a particular set of n generalized speeds u.

The result is each task frame's angular and linear velocity measured and 
expressed in Ground. Using this method is considerably faster than forming the 
6*nt X n Frame Jacobian explicitly and then performing the matrix-vector 
multiply. See the performance analysis below for details.

There is a simplified signature of this method available if you have only a
single frame task.

@param[in]      state
    A State that has already been realized through Position stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the task frames of 
    interest are fixed. These may be in any order and the same body may appear
    more than once if there are multiple task frames on it.
@param[in]      originAoInB
    An array of nt frame origin points Ao for the task frames interest (one 
    per task), each corresponding to one of the bodies B from \a onBodyB, given
    as vectors from each body B's origin Bo to its task frame origin Ao, 
    expressed in frame B.
@param[in]      u
    A mobility-space Vector, such as a set of generalized speeds. The length
    and order must match the mobilities of this system (that is n, the number
    of generalized speeds u, \e not nq, the number of generalized 
    coordinates q).
@param[out]     JFu
    The resulting product JF*u, where JF is the frame task Jacobian. Resized
    if needed to a Vector of nt SpatialVec entries.

@note All frames A fixed to a given body B have the same angular velocity so 
we do not actually need to know the task frames' orientations here, just the
location on B of their origin points Ao. If you have a Transform X_BA giving 
the pose of frame A in the body frame B, you can extract the position vector 
for the origin point Ao using X_BA.p() and pass that as the \a originAoInB 
parameter here.

<h3>Performance discussion</h3>
A call to this method costs 27*nt + 12*(nb+n) flops. If you assume that 
nb ~= n >> 1, you could say this is about 25*(nt+n) flops. In contrast, assuming 
you already have the 6*nt X n Frame Jacobian JF available, you can compute the
JF*u product in about 12*nt*n flops. If you have just one task (nt==1) this
explicit multiplication is about twice as fast; at two tasks it is about even
and for more than two it is more expensive. However forming JF costs about 
180*(n+nt/4) flops (see calcFrameJacobian()). So to form a one-task Jacobian 
and use it once is almost 8X more expensive (192*n vs 25*n), but if you use it 
more than 16 times it is (marginally) cheaper to do it explicitly (for one
task). For example, forming a one-task JF and using it 100 times costs 1392*n 
flops while calling this method 100 times would cost about 2500*n flops.

Conclusion: in almost all practical cases you are better off using this operator
rather than forming JF, even if you have only a single frame task and certainly 
if you have more than two tasks.

@see multiplyByFrameJacobianTranspose(), calcFrameJacobian() **/
void multiplyByFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 originAoInB,
    const Vector&                       u,
    Vector_<SpatialVec>&                JFu) const;

/** Simplified signature for when you just have a single frame task; see the
main signature for documentation.
@returns JF*u, where JF is the single frame task Jacobian. **/
SpatialVec multiplyByFrameJacobian( const State&         state,
                                    MobilizedBodyIndex   onBodyB,
                                    const Vec3&          originAoInB,
                                    const Vector&        u) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&originAoInB, &originAoInB+1);
    Vector_<SpatialVec>                 JFu(1);
    multiplyByFrameJacobian(state, bodies, stations, u, JFu);
    return JFu[0];
}


/** Calculate the n generalized forces f resulting from a set of spatial forces 
(torque,force pairs) F applied at nt task frames Ai fixed to nt bodies Bi. The 
applied forces are spatial vectors (pairs of 3-vectors) expressed in Ground. Use
of this O(n) method is considerably faster than forming the 6*nt X n Jacobian 
explicitly and then performing an O(n^2) matrix-vector multiply.

@param[in]      state
    A State that has already been realized through Position stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the task frames of 
    interest are fixed. These may be in any order and the same body may appear
    more than once if there are multiple task frames on it.
@param[in]      originAoInB
    An array of nt frame origin points Ao for the task frames interest (one 
    per task), each corresponding to one of the bodies B from \a onBodyB, given
    as vectors from each body B's origin Bo to its task frame origin Ao, 
    expressed in frame B.
@param[in]      F_GAo
    A Vector of nt spatial forces, each applied one of the task frames. These
    are expressed in Ground.
@param[out]     f
    The Vector of n generalized forces that results from applying the forces
    \a F_GAo to the task frames. Resized if necessary.

<h3>Performance discussion</h3>
A call to this method costs 33*nt + 18*nb + 11*n flops. If you assume that 
nb ~= n >> 1, you could say this is about 30*(n+nt) flops. In contrast, assuming 
you already have the 6*nt X n Frame Jacobian JF available, you can compute the
~JF*F product in about 12*nt*n flops. For one or two tasks that would be faster
than applying the operator. However forming JF costs about 180*(n+nt/4) flops 
(see calcFrameJacobian()). So to form even a one-task Frame Jacobian and use 
it once is about 6X more expensive than using the operator (192*n vs 30*n), 
but if you use it more than 10 times it is (marginally) cheaper to do it 
explicitly. For example, forming a one-task JF and using it 100 times costs 
around 1392*n flops while calling this method 100 times would cost about 
3000*n flops.

Conclusion: in almost all practical cases you are better off using this operator
rather than forming JF, even if you have only a single frame task and certainly 
if you have more than two tasks.

@see multiplyByFrameJacobian(), calcFrameJacobian() **/
void multiplyByFrameJacobianTranspose
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 originAoInB,
    const Vector_<SpatialVec>&          F_GAo,
    Vector&                             f) const;

/** Simplified signature for when you just have a single frame task. See the
other signature for documentation. **/
void multiplyByFrameJacobianTranspose(  const State&        state,
                                        MobilizedBodyIndex  onBodyB,
                                        const Vec3&         originAoInB,
                                        const SpatialVec&   F_GAo,
                                        Vector&             f) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&originAoInB, &originAoInB+1);
    const Vector_<SpatialVec>           forces(1, F_GAo);
    multiplyByFrameJacobianTranspose(state, bodies, stations, forces, f);
}



/** Explicitly calculate and return the 6*nt x n frame task Jacobian JF for a 
set of nt frame tasks A={Ai} fixed to nt bodies B={Bi}. This matrix maps 
generalized speeds to the Cartesian spatial velocity (angular and linear
velocity) of each frame, measured and expressed in Ground. That is, if you have
a set of n generalized speeds u, you can find the Cartesian spatial velocities 
of task frames A as V_GA = JF*u, where V_GA is a 6*nt column vector. The 
transpose of this matrix maps a 6*nt vector of spatial forces F_GA (expressed 
in Ground and applied to the origins of frames A) to the equivalent set of n 
generalized forces f: f = ~JF*F_GA.

@note It is almost always far more efficient to use multiplyByFrameJacobian() or
multiplyByFrameJacobianTranspose() to form matrix-vector products rather than to
use this method to form the Jacobian explicitly. See the performance discussion 
there.

Overloaded signatures of this method are available to allow you to obtain the
Jacobian either as an nt X n Matrix with SpatialVec elements, or as 6*nt X n 
Matrix with scalar elements.

@param[in]      state
    A State that has already been realized through Position stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the task frames of 
    interest are fixed. These may be in any order and the same body may appear
    more than once if there are multiple task frames on it.
@param[in]      originAoInB
    An array of nt frame origin points Ao for the task frames of interest (one 
    per task), each corresponding to one of the bodies B from \a onBodyB, given
    as vectors from each body B's origin Bo to its task frame origin Ao, 
    expressed in frame B.
@param[out]     JF
    The resulting nt X n frame task Jacobian, with each element a SpatialVec. 
    Resized if necessary.

<h3>Performance discussion</h3>
The cost of a call to this method is about 42*nt + 108*nb + 66*n flops. If we 
assume that nb ~= n >> 1, this is roughly 180*(n+nt/4) flops. Then once the 
Frame Jacobian JF has been formed, each JF*u matrix-vector product costs about
12*nt*n flops to form. When nt is small enough (say one or two tasks), and you
plan to re-use it a lot, this can be computationally efficient; but for single
use or more than a few tasks you can do much better with 
multiplyByFrameJacobian() or multiplyByFrameJacobianTranspose().

@see multiplyByFrameJacobian(), multiplyByFrameJacobianTranspose() **/
void calcFrameJacobian(const State&                         state,
                       const Array_<MobilizedBodyIndex>&    onBodyB,
                       const Array_<Vec3>&                  originAoInB,
                       Matrix_<SpatialVec>&                 JF) const;

/** Simplified signature for when you just have a single frame task. See the
other signature for documentation. **/
void calcFrameJacobian(const State&             state,
                       MobilizedBodyIndex       onBodyB,
                       const Vec3&              originAoInB,
                       RowVector_<SpatialVec>&  JF) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&originAoInB, &originAoInB+1);
    calcFrameJacobian(state, bodies, stations, JF);
}

/** Alternate signature that returns a frame Jacobian as a 6*nt X n Matrix 
rather than as an nt X n Matrix of SpatialVecs. See the other signature for
documentation and important performance considerations.**/
void calcFrameJacobian(const State&                         state,
                       const Array_<MobilizedBodyIndex>&    onBodyB,
                       const Array_<Vec3>&                  originAoInB,
                       Matrix&                              JF) const;

/** Simplified signature for when you just have a single frame task. See the
other signature for documentation. **/
void calcFrameJacobian(const State&             state,
                       MobilizedBodyIndex       onBodyB,
                       const Vec3&              originAoInB,
                       Matrix&                  JF) const 
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&originAoInB, &originAoInB+1);
    calcFrameJacobian(state, bodies, stations, JF);
}

/** Calculate the acceleration bias term for a task frame Jacobian, that is, the
parts of the frames' accelerations that are due only to velocities. This term 
is also known as the Coriolis acceleration, and it is returned here as spatial
accelerations of the frames in Ground.

There is a simplified signature of this method available if you have only a
single frame task.

@param[in]      state
    A State that has already been realized through Velocity stage.
@param[in]      onBodyB
    An array of nt mobilized bodies (one per task) to which the task frames of 
    interest are fixed. These may be in any order and the same body may appear
    more than once if there are multiple task frames on it.
@param[in]      originAoInB
    An array of nt frame origin points Ao for the task frames interest (one 
    per task), each corresponding to one of the bodies B from \a onBodyB, given
    as vectors from each body B's origin Bo to its task frame origin Ao, 
    expressed in frame B.
@param[out]     JFDotu
    The result JFDot*u, where JF is the task frame Jacobian and JFDot its
    time derivative, and u is the set of generalized speeds taken from the
    the supplied \a state.

<h3>Theory</h3>
The spatial velocity V_GA of frame A can be obtained from the generalized
speeds u using the frame Jacobian for A, as V_GA = JF*u. Taking the time 
derivative in G gives
<pre>
    A_GA = JF*udot + JFDot*u
</pre>
This method returns JFDot*u, which depends only on configuration and 
velocities. Note that the same u is used to calculate JFDot, which is linear
in u, so the term JFDot*u is quadratic in u.

<h3>Implementation</h3>
This method just obtains body B's total Coriolis acceleration already available
in the \a state cache and shifts it to the A frame's origin Ao, for each of the
nt task frames. Cost is 48*nt flops.

@see getTotalCoriolisAcceleration(), shiftAccelerationBy()
**/
void calcBiasForFrameJacobian
   (const State&                        state,
    const Array_<MobilizedBodyIndex>&   onBodyB,
    const Array_<Vec3>&                 originAoInB,
    Vector_<SpatialVec>&                JFDotu) const;

/** Alternate signature that returns the bias as a 6*nt-vector of scalars 
rather than as an nt-vector of SpatialVec elements. See the other signature for
documentation. **/
void calcBiasForFrameJacobian(const State&                      state,
                              const Array_<MobilizedBodyIndex>& onBodyB,
                              const Array_<Vec3>&               originAoInB,
                              Vector&                           JFDotu) const;

/** Simplified signature for when you just have a single frame task. 
@returns JFDot*u, where JFDot is the frame task Jacobian time derivative and
u the generalized speeds taken from \a state. **/
SpatialVec calcBiasForFrameJacobian(const State&         state,
                                    MobilizedBodyIndex   onBodyB,
                                    const Vec3&          originAoInB) const
{
    ArrayViewConst_<MobilizedBodyIndex> bodies(&onBodyB, &onBodyB+1);
    ArrayViewConst_<Vec3>               stations(&originAoInB, &originAoInB+1);
    Vector_<SpatialVec>                 JFDotu(1);
    calcBiasForFrameJacobian(state, bodies, stations, JFDotu);
    return JFDotu[0];
}

/**@}**/

//==============================================================================
/** @name               System matrix manipulation
The documentation for the SimbodyMatterSubsystem describes the system equations
in matrix notion, although internal computations are generally matrix-free.
The operators in this section provide the ability to perform fast operations
that can be described in terms of those matrices (e.g., multiply by the mass
matrix) but are actually done using O(n), matrix-free algorithms. There are
also routines here for obtaining the matrices explicitly, although working with
explicit matrices should be avoided whenever performance is an issue.

The mass matrix M and constraint matrix G are the most significant. G=[P;V;A]
is composed of submatrices P for position (holonomic), V for velocity 
(nonholonomic), and A for acceleration-only constraints. These matrices are
sometimes needed separately. Also, these matrices are all in mobility space
(generalized speeds u). When qdot != u, the matrix N in the equation
qdot = N*u becomes important and operators for working with it efficiently
are also provided here. In that case, the position constraint matrix 
in generalized coordinate q space, Pq, can also be accessed. (In terms of
the other matrices, Pq=P*N^-1.) **/
/**@{**/

/** This operator calculates in O(n) time the product M*v where M is the 
system mass matrix and v is a supplied mobility-space vector (that is, it has
one entry for each of the n mobilities). If v is a set of mobility accelerations
(generalized accelerations udot), then the result is a generalized force 
(f=M*udot). Only the supplied vector is used, and M depends only on position 
states, so the result here is not affected by velocities in the State.
Constraints and prescribed motions are ignored.

The current implementation requires about 120*n flops and does not require 
realization of composite-body or articulated-body inertias. 
@par Required stage
  \c Stage::Position **/
void multiplyByM(const State& state, const Vector& a, Vector& Ma) const;

/** This operator calculates in O(n) time the product M^-1*v where M is the 
system mass matrix and v is a supplied vector with one entry per u-space
mobility. If v is a set of generalized forces f, the result is a generalized 
acceleration (udot=M^-1*f). Only the supplied vector is used, and M depends 
only on position states, so the result here is not affected by velocities in 
\a state. In particular, you'll have to obtain your own inertial forces and 
put them in f if you want them included.

@param[in]      state
    This is a State that has been realized through Position stage, from which
    the current system configuration and articulated body inertias are 
    obtained. If necessary, the articulated body inertias will be realized in
    the state the first time this is called. They will then be retained in the
    \a state cache for speed.
@param[in]      v
    This is a generalized-force like vector in mobility space (u-space). If 
    there is any prescribed motion specified using Motion objects or mobilizer
    locking (see below), then only the entries of v corresponding to 
    non-prescribed mobilities are examined by this method; the prescribed ones 
    are not referenced at all. 
@param[out]     MinvV
    This is the result M^-1*v. If there is any prescribed motion specified
    using Motion objects or mobilizer locks (see below), then only the 
    non-prescribed entries in MinvV are calculated; the prescribed ones are set 
    to zero.

<h3>Behavior with prescribed motion</h3>
If you prescribe the motion of one or more mobilizers using Motion objects or
mobilizer locking, the behavior of this method is altered. (This does \e not 
apply if you use Constraint objects to specify the motion.) With prescribed
motion enabled, this method works only with the free (non-prescribed) 
mobilities. Only the entries in \a v corresponding to free mobilities are 
examined, and only the entries in the result \a MinvV corresponding to free 
mobilities are calculated; the others are set to zero.

<h3>Theory</h3>
View the unconstrained, prescribed zero-velocity equations of motion 
M udot + tau = f as partitioned into "free" and "prescribed" variables like 
this:
<pre>
    [M_ff ~M_fp] [udot_f]   [ 0 ]   [f_f]
    [          ] [      ] + [   ] = [   ]
    [M_fp  M_pp] [udot_p]   [tau]   [f_p]
</pre>
The free and prescribed variables have been grouped here for clarity but
in general they are interspersed among the columns and rows of M.

Given that decomposition, this method returns
<pre>
    [udot_f]   [udot_f]   [M_ff^-1  0  ][f_f]
    [      ] = [      ] = [            ][   ]
    [udot_p]   [  0   ]   [   0     0  ][f_p]
</pre>
When there is no prescribed motion M_ff is the entire mass matrix, and the 
result is udot_f=udot=M^-1*f. When there is prescribed motion, M_ff is a 
submatrix of M, and the result is the nf elements of udot_f, with udot_p=0.

<h3>Implementation</h3>
This is a stripped-down version of forward dynamics. It requires the hybrid
free/prescribed articulated body inertias to have been realized and will 
initiate that calculation if necessary the first time it is called for a given 
configuration q. The M^-1*f calculation requires two sweeps of the multibody 
tree, an inward sweep to accumulate forces, followed by an outward sweep to 
propagate accelerations.

<h3>Performance</h3>
If the supplied State does not already contain realized values for the
articulated body inertias, then they will be realized when this operator is 
first called for a new set of positions. Calculating articulated body inertias
is O(n) but relatively expensive. Once the appropriate articulated body 
inertias are available, repeated calls to this operator are very fast, with 
worst case around 80*n flops when all mobilizers have 1 dof. If you want to
force realization of the articulated body inertias, call the method
realizeArticulatedBodyInertias().

@par Required stage
  \c Stage::Position 

@see multiplyByM(), calcMInv(), realizeArticulatedBodyInertias() **/ 
void multiplyByMInv(const State& state,
    const Vector&               v,
    Vector&                     MinvV) const;

/** This operator explicitly calculates the n X n mass matrix M. Note that this
is inherently an O(n^2) operation since the mass matrix has n^2 elements 
(although only n(n+1)/2 are unique due to symmetry). <em>DO NOT USE THIS CALL 
DURING NORMAL DYNAMICS</em>. To do so would change an O(n) operation into an 
O(n^2) one. Instead, see if you can accomplish what you need with O(n) operators
like multiplyByM() which calculates the matrix-vector product M*v in O(n)
without explicitly forming M. Also, don't invert this matrix numerically to get
M^-1. Instead, call the method calcMInv() which can produce M^-1 directly.
@see multiplyByM(), calcMInv() **/
void calcM(const State&, Matrix& M) const;

/** This operator explicitly calculates the inverse of the part of the system
mobility-space mass matrix corresponding to free (non-prescribed)
mobilities. The returned matrix is always n X n, but rows and columns 
corresponding to prescribed mobilities are zero. This is an O(n^2) operation, 
which is of course within a constant factor of optimal for returning a matrix 
with n^2 elements explicitly. (There are actually only n(n+1)/2 unique elements
since the matrix is symmetric.) <em>DO NOT USE THIS CALL DURING NORMAL 
DYNAMICS</em>. To do so would change an O(n) operation into an O(n^2) one. 
Instead, see if you can accomplish what you need with O(n) operators like 
multiplyByMInv() which calculates the matrix-vector product M^-1*v in O(n) 
without explicitly forming M or M^-1. If you need M explicitly, you can get it
with the calcM() method.
@see multiplyByMInv(), calcM() **/
void calcMInv(const State&, Matrix& MInv) const;

/** This operator calculates in O(m*n) time the m X m "projected inverse mass 
matrix" or "constraint compliance matrix" W=G*M^-1*~G, where G (mXn) is the 
acceleration-level constraint Jacobian mapped to generalized coordinates,
and M (nXn) is the unconstrained system mass matrix. In case there is prescribed
motion specified with Motion objects or mobilizer locking, M^-1 here is really
M_ff^-1, that is, it is restricted to the free (non-prescribed) mobilities, but 
scattered into a full n X n matrix (conceptually). See multiplyByMInv() and 
calcMInv() for more information.

W is the projection of the inverse mass matrix into the constraint coordinate
space (that is, the vector space of the multipliers lambda). It can be used to 
solve for the constraint forces that will eliminate a given constraint 
acceleration error:
<pre>
    (1)     W * lambda = aerr
    (2)     aerr = G*udot - b(t,q,u)
</pre>
where udot is an unconstrained generalized acceleration. Note that you can
view equation (1) as a dynamic system in a reduced set of m generalized
coordinates, with the caveat that W may be singular.

In general W is singular and does not uniquely determine lambda. Simbody 
normally calculates a least squares solution for lambda so that loads are 
distributed among redundant constraints. 

@note If you just need to multiply W by a vector or matrix, you do not need
to form W explicitly. Instead you can use the method described in the 
Implementation section to produce a W*v product in the O(n) time it takes to 
compute a single column of W.

<h3>Implementation</h3>
We are able to form W without forming G or M^-1 and without performing any 
matrix-matrix multiplies. Instead, W is calculated using m applications of 
O(n) operators:
    - multiplyByGTranspose() by a unit vector to form a column of ~G
    - multiplyByMInv() to form a column of M^-1 ~G
    - multiplyByG() to form a column of W
    
Even if G and M^-1 were already available, computing W by matrix multiplication
would cost O(m^2*n + m*n^2) time and O(m*n) intermediate storage. Here we do 
it in O(m*n) time with O(n) intermediate storage, which is a \e lot better.
     
@see multiplyByG(), calcG(), multiplyByGTranspose(), calcGTranspose()
@see multiplyByMInv(), calcMInv() **/
void calcProjectedMInv(const State&   s,
                       Matrix&        GMInvGt) const;

/** Given a set of desired constraint-space speed changes, calculate the
corresponding constraint-space impulses that would cause those changes. Here we 
are solving the equation
<pre>
    W * impulse = deltaV
</pre>
for \e impulse, where W=G*M^-1*~G is the "projected inverse mass matrix" as 
described for calcProjectedMInv(). In general W is singular due to constraint
redundancies, so the solution for \e impulse is not unique. Simbody handles 
redundant constraints by finding least squares solutions, and this operator 
method duplicates the method Simbody uses for determining the rank and 
performing the factorization of W. 

@param[in]      state
    The State whose generalized coordinates and speeds define the matrix W.
    Must already be realized to Dynamics stage.
@param[in]      deltaV
    The set of desired velocity changes to be produced by the impulse, in 
    constraint space. These will consist of observed velocity constraint 
    violations (-verr) and constraint violations that would be generated by
    impulsive applied forces (-G*M^-1*f).
@param[out]     impulse
    The set of constraint multiplier-space impulses that will produce the
    desired velocity changes without violating the constraints.

To convert these constraint-space impulses into updates to the mobility-space
generalized speeds u, use code like this:
@code
    const SimbodyMatterSubsystem& matter=...;
    Vector deltaV=...;  // constraint space speed change desired; length m
    Vector impulse;     // constraint space impulses; length m
    solveForConstraintImpulses(state, deltaV, impulse);
    Vector f;           // mobility space impulses; length n
    Vector du;          // change to generalized speeds u; length n
    matter.multiplyByGTranspose(s,impulse,f);
    matter.multiplyByMInv(s,f,du);
    state.updU() += du; // update generalized speeds
@endcode 

Note that the length of the constraint-space vectors is m=mp+mv+ma, the total 
number of acceleration-level constraints including the second time derivatives
of the position (holonomic) constraints, the first time derivatives of the 
velocity (nonholonomic) constraints, and the acceleration-only constraints. 
@see calcProjectedMInv(), multiplyByGTranspose(), multiplyByMInv() **/
void solveForConstraintImpulses(const State&     state,
                                const Vector&    deltaV,
                                Vector&          impulse) const;


/** Returns Gulike = G*ulike, the product of the mXn acceleration 
constraint Jacobian G and a "u-like" (mobility space) vector of length n. 
m is the number of active acceleration-level constraint equations, n is the 
number of mobilities. This is an O(m+n) operation.

If you are going to call this method repeatedly at the same time, positions and
velocities, you should precalculate the bias term once and supply it to the
alternate signature of this method. See the Implementation section for more 
information.

@pre \a state realized to Velocity stage
@par Implementation
This is accomplished by treating the input vector \a ulike as though it were
a set of generalized accelerations (for nonholonomic and acceleration-only
constraints) or generalized speeds (for holonomic constraints). These are 
mapped to body accelerations (or velocities) in O(n) time. See
calcBodyAccelerationFromUDot() for more information (converting from 
generalized speeds to velocities is just multiplying by the System Jacobian).
The method calcBiasForMultiplyByG() is used to determine the state-dependent
term of the constraint error equations. Then a second call is made to
evaluate the bias term aerr(t,q,u;0)=-b(t,q,u). We then calculate 
Gulike = aerr(t,q,u;ulike)-aerr(t,q,u;0) in O(m) time.
@see calcBiasForMultiplyByG() **/
void multiplyByG(const State&  state,
                 const Vector& ulike,
                 Vector&       Gulike) const {
    Vector bias;
    calcBiasForMultiplyByG(state, bias);
    multiplyByG(state, ulike, bias, Gulike);
}


/** Multiply Gulike=G*ulike using the supplied precalculated bias vector to 
improve performance (approximately 2X) over the other signature. 
@see calcBiasForMultiplyByG() **/
void multiplyByG(const State&  state,
                 const Vector& ulike,
                 const Vector& bias,
                 Vector&       Gulike) const;

/** Calculate the bias vector needed for the higher-performance signature of
the multiplyByG() method above. 

@param[in]      state
    Provides time t, positions q, and speeds u; must be realized through
    Velocity stage so that all body spatial velocities are known.
@param[out]     bias
    This is the bias vector for use in repeated calls to multiplyByG(). It
    will be resized if necessary to length m=mp+mv+ma, the total number of 
    active acceleration-level constraint equations. 

@pre \a state realized to Velocity stage
@par Implementation
This method uses either velocity- or acceleration- level constraint error
functions with zero input to determine the bias term for use in 
multiplyByG(). Body quantities and generalized quantities are supplied to each 
of the m active constraints' (constant time) error methods to calculate
<pre>
   pverr(t,q,u;ulike)=G*ulike - c(t,q)    (holonomic) 
or aerr(t,q,u;ulike)=G*ulike - b(t,q,u)   (nonholonomic or acceleration-only)
</pre>
with ulike=0, giving the bias term in O(m) time. 

If you want the acceleration-level bias terms b for all the constraints, even
if they are holonomic, use calcBiasForAccelerationConstraints(). **/
void calcBiasForMultiplyByG(const State& state,
                            Vector&      bias) const;

/** This O(m*n) operator explicitly calculates the m X n acceleration-level 
constraint Jacobian G which appears in the system equations of 
motion. Consider using the multiplyByG() method instead of this one, 
which forms the matrix-vector product G*v in O(m+n) time without explicitly 
forming G.

@par Implementation
This method generates G columnwise using repeated calls to multiplyByG(), 
which makes use of the constraint error methods to perform a G*v product
in O(m+n) time. To within numerical error, for non-working constraints
this should be identical to the transpose of the matrix returned by calcGt() 
which uses the constraint force methods instead. 
@see multiplyByG(), calcGt(), calcPq() **/
void calcG(const State& state, Matrix& G) const;


/** Calculate the acceleration constraint bias vector, that is, the terms in
the acceleration constraints that are independent of the accelerations.

@param[in]      state
    Provides time t, positions q, and speeds u; must be realized through
    Velocity stage so that all body spatial velocities are known.
@param[out]     bias
    This is the bias vector for all the acceleration constraint equations
    together. It will be resized if necessary to length m=mp+mv+ma, the total 
    number of active acceleration-level constraint equations. 

@pre \a state realized to Velocity stage
@par Implementation
We have constant-time constraint acceleration error methods 
<pre>   
paerr(t,q,u;udot)=P*udot - b_p(t,q,u) 
vaerr(t,q,u;udot)=V*udot - b_v(t,q,u) 
 aerr(t,q,u;udot)=A*udot - b_a(t,q,u)   
</pre>
that together define the acceleration constraint equation G*udot-b=0
where G=[P;V;A] and b=[b_p b_v b_a]. There is one of these error functions 
for each %Constraint, with paerr() the twice-differentiated position (holonomic)
constraints, vaerr() the once-differentiated velocity (nonholonomic)
constraints, and aerr() the acceleration-only constraints. This method 
sets \c udot = 0 and invokes each of those methods to obtain 
bias = -[b_p b_v b_a].

<h3>Performance note</h3>
The actual acceleration constraint functions require both udot and body 
accelerations for the constrained bodies; even with udot==0 body accelerations 
may have a non-zero velocity-dependent component (the coriolis accelerations). 
Those are already available in the state, but only as accelerations in Ground. 
For constraints that have a non-Ground Ancestor, we have to convert the 
accelerations to A at a cost of 105 flops/constrained body. **/
void calcBiasForAccelerationConstraints(const State& state,
                                        Vector&      bias) const;


/** Returns f = ~G*lambda, the product of the n X m transpose of the 
acceleration constraint Jacobian G (=[P;V;A]) and a multiplier-like vector 
\a lambda of length m, returning a generalized-force like quantity \a f of
length n. m=mp+mv+ma is the total number of active constraint equations, 
n (==nu) is the number of mobilities (generalized speeds u). If lambda is a set
of constraint multipliers, then f=~G*lambda is the set of forces generated by
the constraints, mapped into generalized forces. This is an O(m+n) operation.

Because the velocity (non-holonomic) or acceleration-only constraint Jacobians
V and A can have velocity dependence, the \a state supplied here must generally
be realized through Velocity stage. If the system has only position (holonomic)
constraints then the \a state need be realized only through Position stage.

@param[in]      state
    A State that has been realized through Velocity stage (or Position stage
    if the system has only position constraints). Time, configuration,
    and velocities if needed are taken from \a state.
@param[in]      lambda
    A multiplier-like vector to be multiplied by ~G. Its length must be the
    same as the total number of active constraint equations m.
@param[out]     f
    This is the generalized force-like output. It will be resized if necessary
    to length equal to the number of mobilities (generalized speeds) n (==nu). 

@par Implementation
This is accomplished by treating the input vector \a lambda as though it were
a set of Lagrange multipliers, then calling each of the active Constraints' 
(constant time) force generation methods, providing the appropriate subset of 
the multipliers each time. That gives body forces F0 and mobility forces f0 in 
O(m) time. We then use the equivalent of multiplyBySystemJacobianTranspose() 
to convert the returned body spatial forces to generalized forces in O(n) 
time, and finally return the generalized force-like result f = ~J*F0 + f0. 
@see multiplyByG(), multiplyBySystemJacobianTranspose() **/
void multiplyByGTranspose(const State&  state,
                          const Vector& lambda,
                          Vector&       f) const;
    
/** This O(nm) operator explicitly calculates the n X m transpose of the 
acceleration-level constraint Jacobian G = [P;V;A] which appears in the system 
equations of motion. This method generates ~G columnwise use the constraint 
force generating methods which map constraint multipliers to constraint forces.
To within numerical error, this should be identical to the transpose of
the matrix returned by calcG() which uses a different method. Consider using 
the multiplyByGTranspose() method instead of this one, which forms the 
matrix-vector product ~G*v in O(n) time without explicitly forming ~G.
@par Required stage
  \c Stage::Velocity
@see calcG(), multiplyByGTranspose() **/
void calcGTranspose(const State&, Matrix& Gt) const;


/** Calculate in O(n) time the product Pq*qlike where Pq is the mp X nq 
position (holonomic) constraint Jacobian and \a qlike is a "q-like" 
(generalized coordinate space) vector of length nq. Here mp is the number of 
active position-level constraint equations in this system.

If you are going to call this method repeatedly at the same time t and 
configuration q and want maximum efficiency, you can gain a factor of almost
2X by precalculating a bias term once using calcBiasForMultiplyByPq() and 
supplying it to the alternate signature of this method. See the Theory
section below for an explanation of the bias term.

@pre \a state realized to Position stage

<h3>Theory</h3>
Simbody's position (holonomic) constraints are defined by the constraint 
error equation 
<pre>
    (1)    perr(t;q) = p(t,q)
</pre>
where we try to maintain perr=0 at all times. We also have available time 
derivatives of equation (1); the first time derivative is relevant here:
<pre>
    (2)    pverr(t,q;qdot) = dperr/dt = Pq * qdot + Pt
</pre>
where Pq=Dperr/Dq and Pt=Dperr/Dt (capital "D" means partial derivative).
Pt=Pt(t,q) is called the "bias" term. (Note that because u=N^-1*qdot we also
have Pq=P*N^-1, where P=Dpverr/Du is the very useful mobility-space holonomic 
constraint Jacobian.) Eq. (2) can be used to perform efficient multiplication 
by Pq, since it can be used to calculate Pq*qlike+Pt, and a second evaluation 
at qlike=0 can be used to calculate the unwanted bias term for removal: 
<pre>
    (3)    Pq*qlike = pverr(t,q;qlike) - pverr(t,q;0)  
</pre>
Despite appearances, eq. (2) calculates its result in constant time per
constraint equation, for a total cost that is O(n) or more strictly O(mp+nq).
The matrix Pq is never actually formed; instead the matrix-vector product
is calculated directly.

<h3>Implementation</h3>
We treat the input vector \a qlike as though it were a set of generalized 
coordinate derivatives qdot. These are mapped to body velocities V in O(n) 
time, using V=Jq*qdot, where Jq is the coordinate space system Jacobian 
(partial velocity matrix), with Jq=J*N^-1. Then the body velocities and qdots 
are supplied to each of the mp active position constraints' (constant time) 
velocity error methods to get pverr(t,q;qlike)=Pq*qlike-Pt in O(n) time. A 
second call is made to evaluate the bias term pverr(t,q;0)=-Pt. We then 
calculate the result \a PqXqlike = pverr(t,q;qlike)-pverr(t,q;0) in O(n) time
using equation (3). 

@see calcBiasForMultiplyByPq() **/
void multiplyByPq(const State&  state,
                  const Vector& qlike,
                  Vector&       PqXqlike) const {
    Vector biasp;
    calcBiasForMultiplyByPq(state, biasp);
    multiplyByPq(state, qlike, biasp, PqXqlike);
}


/** Multiply Pq*qlike using the supplied precalculated bias vector to 
improve performance (approximately 2X) over the other signature. 
@see calcBiasForMultiplyByPq() **/
void multiplyByPq(const State&  state,
                 const Vector&  qlike,
                 const Vector&  biasp,
                 Vector&        PqXqlike) const;

/** Calculate the bias vector needed for the higher-performance signature of
the multiplyByPq() method above. 

@param[in]      state
    Provides time t, and positions q; must be realized through
    Position stage so that all body spatial poses are known.
@param[out]     biasp
    This is the bias vector for use in repeated calls to multiplyByPq(). It
    will be resized if necessary to length mp, the total number of 
    active position-level (holonomic) constraint equations. 

@pre \a state realized to Position stage

See multiplyByPq() for theory and implementation; this method is just 
performing the qlike=0 case described there for calculating the bias term Pt.
**/
void calcBiasForMultiplyByPq(const State& state,
                             Vector&      biasp) const;

/** This O(m*n) operator explicitly calculates the mp X nq position-level 
(holonomic) constraint Jacobian Pq (=P*N^-1), the partial derivative of the
position error equations with respect to q. Consider using the multiplyByPq() 
method instead of this one, which forms the matrix-vector product Pq*v in 
O(m+n) time without explicitly forming Pq.

Note that quaternion normalization constraints are \e not included in mp; we
do not consider those holonomic constraints.

@pre \a state realized to Position stage

@param[in]      state
    A State realized through Position stage so that time and the pose 
    (configuration) of each body is known.
@param[out]     Pq
    The position constraint Jacobian Dperr/Dq. This will be resized to
    mp X nq if necessary.

@par Implementation
This method generates Pq columnwise using repeated calls to multiplyByPq(), 
which makes use of the position constraint velocity-level error methods to 
perrform a Pq*v product in O(m+n) time. See multiplyByPq() for a more 
detailed explanation. If Pq's columns are in contiguous memory we'll work
in place, otherwise columns are generated into a contiguous temporary and
then copied into Pq.

@see multiplyByPq() **/
void calcPq(const State& state, Matrix& Pq) const;


/** Returns f = ~Pq*lambdap, the product of the n X mp transpose of the 
position (holonomic) constraint Jacobian Pq (=P*N^-1) and a multiplier-like 
vector \a lambdap of length mp, returning a generalized-force like quantity 
\a f of length n. mp is the number of active position constraint equations, 
n (==nu) is the number of mobilities (generalized speeds u). If lambdap is a set
of mp constraint multipliers, then f=~G*lambdap is the set of forces generated 
by the position constraints, mapped into generalized forces. This is an 
O(mp+n) operation.

A holonomic constraint Jacobian cannot have a velocity dependence, so the
\a state need be realized only to Position stage here.

@param[in]      state
    A State that has been realized through Position stage. Time and
    configuration are taken from \a state.
@param[in]      lambdap
    A multiplier-like vector to be multiplied by ~Pq. Its length must be the
    same as the number of active position constraint equations mp.
@param[out]     f
    This is the generalized force-like output. It will be resized if necessary
    to length equal to the number of mobilities (generalized speeds) n (==nu). 

@par Implementation
This is accomplished by treating the input vector \a lambdap as though it were
a set of Lagrange multipliers, then calling each of the active holonomic
Constraints' (constant time) force generation methods, providing the 
appropriate subset of the multipliers each time. That gives body forces F0 and 
mobility forces f0 in O(mp) time. We then use the equivalent of 
multiplyBySystemJacobianTranspose() to convert the returned body spatial 
forces to generalized forces in O(n) time, and finally return the generalized 
force-like result f = ~J*F0 + f0. 
@see multiplyByPq(), multiplyBySystemJacobianTranspose() **/
void multiplyByPqTranspose(const State&  state,
                           const Vector& lambdap,
                           Vector&       f) const;

/** This O(m*n) operator explicitly calculates the nq X mp transpose of 
the position-level (holonomic) constraint Jacobian Pq (=P*N^-1), the partial 
derivative of the position error equations with respect to q. Consider using 
the multiplyByPqTranspose() method instead of this one, which forms the 
matrix-vector product ~Pq*v in O(m+n) time without explicitly forming ~Pq.

Note that quaternion normalization constraints are \e not included in mp; we
do not consider those holonomic constraints.

@pre \a state realized to Position stage

@param[in]      state
    A State realized through Position stage so that time and the pose 
    (configuration) of each body is known.
@param[out]     Pqt
    The transposed position constraint Jacobian ~Pq=(Dperr/Dq)^T. This will be
    resized to nq X mp if necessary.

@par Implementation
This method generates \a Pqt columnwise using repeated calls to 
multiplyByPqTranspose(), which makes use of the position constraint force
generating methods to perform a ~Pq*v product in O(m+n) time. See 
multiplyByPqTranspose() for a more detailed explanation. If Pqt's columns 
are in contiguous memory we'll work in place, otherwise columns are generated 
into a contiguous temporary and then copied into Pqt.

@see multiplyByPqTranspose() **/
void calcPqTranspose(const State& state, Matrix& Pqt) const;

/** Returns the mp X nu matrix P which is the Jacobian of the first time
derivative of the holonomic (position) constraint errors with respect to the 
generalized speeds u; that is, P = partial( dperr/dt )/partial(u). Here mp is 
the number of holonomic constraint equations (not including quaternion 
normalization constraints) and nu is the total number of generalized speeds as 
found in the supplied State. P is resized if necessary; an error will be thrown
if the Matrix is not the right size and not resizeable.

@pre \a state is realized to Position stage
@par Complexity:
Calculates the m X n matrix in O(m*n) time, which is good if you really need
this matrix. However, in many cases what is really needed is the product
of this matrix with a vector which can be done in O(n) time; consider whether
you really need the whole matrix explicitly.
@par Required stage
  \c Stage::Position **/
void calcP(const State& state, Matrix& P) const;

/** Returns the nu X mp matrix ~P - see calcP() for a description. 
@par Required stage
  \c Stage::Position **/
void calcPt(const State& state, Matrix& Pt) const;


/** Calculate out_q = N(q)*in_u (like qdot=N*u) or out_u = ~N*in_q. Note that 
one of "in" and "out" is always "q-like" while the other is "u-like", but which
is which changes if the matrix is transposed. Note that the transposed 
operation here is the same as multiplying by N on the right, with the Vectors 
viewed as RowVectors instead. This is an O(n) operator since N is block 
diagonal.
@par Required stage
  \c Stage::Position **/
void multiplyByN(const State& s, bool transpose, 
                 const Vector& in, Vector& out) const;

/** Calculate out_u = NInv(q)*in_q (like u=NInv*qdot) or out_q = ~NInv*in_u. 
Note that one of "in" and "out" is always "q-like" while the other is "u-like",
but which is which changes if the matrix is transposed. Note that the 
transposed operation here is the same as multiplying by NInv on the right, 
with the Vectors viewed as RowVectors instead. This is an O(N) operator since 
NInv is block diagonal. The configuration q is taken from the supplied state.
@par Required stage
  \c Stage::Position **/
void multiplyByNInv(const State& s, bool transpose, 
                    const Vector& in, Vector& out) const;

/** Calculate out_q = NDot(q,u)*in_u or out_u = ~NDot(q,u)*in_q. This is used,
for example, as part of the conversion between udot and qdotdot. Note that one
of "in" and "out" is always "q-like" while the other is "u-like", but which is
which changes if the matrix is transposed. Note that the transposed operation 
here is the same as multiplying by NDot on the right, with the Vectors viewed
as RowVectors instead. This is an O(N) operator since NDot is block diagonal.
Configuration q and generalized speeds u are taken from the supplied state.
@par Required stage
  \c Stage::Velocity **/
void multiplyByNDot(const State& s, bool transpose, 
                    const Vector& in, Vector& out) const;

/**@}**/


//==============================================================================
/** @name                  Miscellaneous Operators

Operators make use of the State but do not write their results back
into the State, not even into the State cache. **/
/**@{**/

/** This is the primary forward dynamics operator. It takes a state which
has been realized to the Dynamics stage, a complete set of forces to apply,
and returns the accelerations that result. Only the forces supplied here,
and those calculated internally from prescribed motion, constraints, and
centrifugal effects, affect the results. Acceleration constraints are 
always satisfied on return as long as the constraints are consistent. 
If the position and velocity constraints aren't already satisified in the 
State, results are harder to interpret physically, but they will still be 
calculated and the acceleration constraints will still be satisfied. No 
attempt will be made to satisfy position and velocity constraints, or to 
set prescribed positions and velocities, nor even to check whether these 
are satisfied; position and velocity constraint and prescribed positions 
and velocities are simply irrelevant here.

Given applied forces f_applied, this operator solves this set of equations:
<pre>
     M udot + tau + ~G lambda + f_inertial = f_applied       (1)
                                  G udot   = b               (2)
                                    udot_p = udot_p(t,q,u,z) (3)
</pre>
where udot={udot_f,udot_p}, tau={0,tau_p}. The unknowns are: the free 
generalized accelerations udot_f, the constraint multipliers lambda, and the 
prescribed motion generalized forces tau_p. A subset udot_p of udot may have 
been prescribed as a known function of state via Motion objects or locks 
associated with the mobilized bodies. On return all the entries in udot will 
have been set to their calculated or prescribed values, and body spatial 
accelerations A_GB (that is, measured and expressed in Ground) are also 
returned. Lambda and tau_p are necessarily calculated but are not returned here.

f_applied is the set of generalized (mobility) forces equivalent to the 
\a appliedMobilityForces and \a appliedBodyForces arguments supplied here. 
That is,
<pre>
    f_applied = appliedMobilityForces + ~J * appliedBodyForces
</pre> 
where J is the system Jacobian mapping between spatial and generalized
coordinates. Typically these forces will have been calculated as a function of 
state so we will have f_applied(t,q,u,z).

M(t,q), G(t,q,u), and b(t,q,u) are defined by the mobilized bodies and 
constraints present in the system. f_inertial(q,u) includes the 
velocity-dependent gyroscopic and coriolis forces due to rigid body 
rotations and is extracted internally from the already-realized state. 

Note that this method does not allow you to specify your own prescribed udots; 
those are calculated from the mobilizers' state-dependent Motion specifications
(or are zero due to mobilizer locks) that are already part of the system.

This is an O(n*m + m^3) operator where n is the number of generalized speeds
and m the number of constraint equations (mobilities with prescribed motion are
counted in n, not m).

@par Required stage
  \c Stage::Dynamics **/ 
void calcAcceleration
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    Vector&                    udot,    // output only; no prescribed motions
    Vector_<SpatialVec>&       A_GB) const;

/** This operator is similar to calcAcceleration() but ignores the effects of
acceleration constraints although it obeys prescribed accelerations. The 
supplied forces, prescribed motion forces, and velocity-induced centrifugal 
and gyroscopic effects are properly accounted for, but any forces that would 
have resulted from enforcing the contraints are not present. This operator 
solves the equations
<pre>
            M udot + tau + f_inertial = f_applied           (1)
                               udot_p = udot_p(t,q,u,z)     (2)
</pre>
where udot={udot_f,udot_p}, tau={0,tau_p}. The unknowns are the free 
generalized accelerations udot_f and the prescribed motion generalized forces
tau_p. f_inertial contains the velocity-dependent gyroscopic and coriolis
forces due to rigid body rotations. No constraint forces are included.

On return all the entries in udot will have been set to their calculated or 
prescribed values, and body spatial accelerations A_GB (that is, measured and 
expressed in Ground) are also returned. tau_p is not returned.

This is an O(n) operator.

@par Required stage
  \c Stage::Dynamics **/ 
void calcAccelerationIgnoringConstraints
   (const State&                state,
    const Vector&               appliedMobilityForces,
    const Vector_<SpatialVec>&  appliedBodyForces,
    Vector&                     udot,    
    Vector_<SpatialVec>&        A_GB) const;



/** This is the inverse dynamics operator for the tree system; if there are
any constraints or prescribed motion they are ignored. This method solves
<pre>
     f_residual = M udot + f_inertial - f_applied
</pre>
for f_residual in O(n) time, meaning that the mass matrix M is never formed. 
Inverse dynamics is considerably faster than forward dynamics, even though 
both are O(n) in Simbody.

In the above equation we solve for the residual forces \c f_residual given
desired accelerations and (optionally) a set of applied forces. Here 
\c f_applied is the mobility-space equivalent of all the applied forces
(including mobility and body forces), \c f_inertial is the mobility-space
equivalent of the velocity-dependent inertial forces due to rigid 
body rotations (coriolis and gyroscopic forces), and \c udot is the 
given set of values for the desired generalized accelerations. The returned 
\c f_residual is the additional generalized force (that is, mobility 
force) that would have to be applied at each mobility to give the desired
\c udot. The inertial forces depend on the velocities \c u already realized 
in the State. Otherwise, only the explicitly-supplied forces affect the 
results of this operator; any forces that may be present elsewhere in 
the system are ignored.

@param[in] state
     A State valid for the containing System, already realized to
     Stage::Velocity.
@param[in] appliedMobilityForces
     One scalar generalized force applied per mobility. Can be zero
     length if there are no mobility forces; otherwise must have exactly 
     one entry per mobility in the matter subsystem.
@param[in] appliedBodyForces
     One spatial force for each body. A spatial force is a force applied
     to the body origin and a torque on the body, each expressed in the 
     Ground frame. Gravity, if present, is specified here as a body force.
     The supplied Vector must be either zero length (interpreted as all-zero)
     or have exactly one entry per body in the matter subsystem, starting with
     Ground as body zero.
@param[in] knownUdot
     These are the desired generalized accelerations, one per mobility. 
     If this is zero length it will be treated as all-zero; otherwise 
     it must have exactly one entry per mobility in the matter subsystem.
@param[out] residualMobilityForces
     These are the residual generalized forces which, if added to the applied
     forces, would produce the given \a knownUdot in forward dynamics (assuming
     the system is unconstrained). This will be resized if necessary to have 
     length nu; that is, one scalar entry per mobility. You can view this as a 
     measure of how much the given \a knownUdot fails to satisfy the equations 
     of motion.

@par Required stage
  \c Stage::Velocity 

@see calcResidualForce(), multiplyByM()
@see calcAcceleration(), calcAccelerationIgnoringConstraints() **/
void calcResidualForceIgnoringConstraints
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    Vector&                    residualMobilityForces) const;


/** This is the inverse dynamics operator for when you know both the 
accelerations and Lagrange multipliers for a constrained system. Prescribed
motion is ignored. Using position and velocity from the given state, a set of 
applied forces, and a known set of generalized accelerations udot and 
constraint multipliers lambda, it calculates the additional generalized forces
that would be required to satisfy Newton's 2nd law, f=Ma. That is, this 
operator returns
<pre>
    f_residual = M udot + ~G lambda + f_inertial - f_applied
</pre>
where f_applied is the mobility-space equivalent to all the applied forces 
(including mobility and body forces), f_inertial is the mobility-space 
equivalent of the velocity-dependent inertial forces due to rigid body 
rotations (coriolis and gyroscopic forces), and the udots and lambdas are given
values of the generalized accelerations and constraint multipliers, resp.

Note that there is no requirement that the given udots satisfy the constraint
equations; we simply solve the above equation for \c f_residual.

The inertial forces depend on the velocities \c u already realized in the State.
Otherwise, only the explicitly-supplied forces affect the results of this 
operator; any forces that may be defined elsewhere in the system are ignored 
here.

@param[in] state
     A State valid for the containing System, already realized to
     \c Stage::Velocity.
@param[in] appliedMobilityForces
     One scalar generalized force applied per mobility. Can be zero
     length if there are no mobility forces; otherwise must have exactly 
     one entry per mobility in the matter subsystem.
@param[in] appliedBodyForces
     One spatial force for each body. A spatial force is a force applied
     to the body origin and a torque on the body, each expressed in the 
     Ground frame. Gravity, if present, is specified here as a body force.
     The supplied Vector must be either zero length (interpreted as all-zero)
     or have exactly one entry per body in the matter subsystem, starting with
     Ground as body zero.
@param[in] knownUdot
     These are the specified generalized accelerations, one per mobility so
     the length should be nu. If this is zero length it will be treated as 
     all-zero of length nu; otherwise it must have exactly one entry per 
     mobility in the matter subsystem.
@param[in] knownLambda
     These are the specified Lagrange multipliers, one per constraint
     equation. If this is zero length it will be treated as all-zero; otherwise 
     it must have exactly m entries, where m=mp+mv+ma is the total number of
     position, velocity, and acceleration-only constraints. There are no
     entries here corresponding to quaternion constraints, which do not 
     generate forces.
@param[out] residualMobilityForces
     These are the residual generalized forces which, if added to the applied
     forces along with the constraint forces ~G*lambda, would produce the 
     given \a knownUdot in unconstrained forward dynamics. This will be resized
     if necessary to have length nu; that is, one scalar entry per mobility.
     You can view this as a measure of how much the given udot and lambda fail
     to satisfy the equations of motion.

@par Required stage
  \c Stage::Velocity 

@see calcResidualForceIgnoringConstraints() **/
void calcResidualForce
   (const State&               state,
    const Vector&              appliedMobilityForces,
    const Vector_<SpatialVec>& appliedBodyForces,
    const Vector&              knownUdot,
    const Vector&              knownLambda,
    Vector&                    residualMobilityForces) const;


/** This operator calculates the composite body inertias R given a State 
realized to Position stage. Composite body inertias are the spatial mass 
properties of the rigid body formed by a particular body and all bodies 
outboard of that body as if all the outboard mobilizers were welded in their 
current orientations. 

This is a very fast O(n) operator.

@par Required stage
  \c Stage::Position **/
void calcCompositeBodyInertias(const State&    state,
    Array_<SpatialInertia,MobilizedBodyIndex>& R) const;



/** Given a complete set of n generalized accelerations udot, this kinematic 
operator calculates in O(n) time the resulting body accelerations, including 
velocity-dependent terms taken from the supplied \a state.

@pre \a state must already be realized to Velocity stage
@param[in] state
    The State from which position- and velocity- related terms are taken; 
    must already have been realized to Velocity stage.
@param[in] knownUDot
    A complete set of generalized accelerations. Must have the same length 
    as the number of mobilities nu, or if length zero the udots will be taken 
    as all zero in which case only velocity-dependent (Coriolis) accelerations 
    will be returned in \a A_GB.
@param[out] A_GB
    Spatial accelerations of all the body frames measured and expressed in
    the Ground frame, resulting from supplied generalized accelerations 
    \a knownUDot and velocity-dependent acceleration terms taken from 
    \a state. This will be resized if necessary to the number of bodies 
    <em>including</em> Ground so that the returned array may be indexed by 
    MobilizedBodyIndex with A_GB[0]==0 always. The angular acceleration
    vector for MobilizedBody i is A_GB[i][0]; linear acceleration of the
    body's origin is A_GB[i][1].

@par Theory
The generalized speeds u and spatial velocities V are related by the system
Jacobian J as V=J*u. Thus the spatial accelerations A=Vdot=J*udot+Jdot*u.

@par Implementation
The Coriolis accelerations Jdot*u are already available in a State realized
to Velocity stage. The J*udot term is equivalent to an application of 
multiplyBySystemJacobian() to the \a knownUdot vector. The current 
implementation uses 12*nu + 18*nb flops to produce nb body accelerations.

@par Required stage
  \c Stage::Velocity 
  
@see multiplyBySystemJacobian(), getTotalCoriolisAcceleration() **/
void calcBodyAccelerationFromUDot(const State&         state,
                                  const Vector&        knownUDot,
                                  Vector_<SpatialVec>& A_GB) const;


/** Treating all Constraints together, given a comprehensive set of m 
Lagrange multipliers \e lambda, generate the complete set of body spatial forces
and mobility (generalized) forces applied by all the Constraints.

Spatial forces are applied at each body's origin and the moment and force
vectors therein are expressed in the Ground frame. Watch the 
sign -- normally constraint forces have opposite sign from applied forces, 
because our equations of motion are 
    <pre>   M udot + ~G lambda = f_applied  </pre>
If you want to take Simbody-calculated multipliers and use them to generate 
forces that look like applied forces, negate the multipliers in the argument
passed to this call.

State must be realized to Stage::Velocity to call this operator (although 
typically the multipliers are obtained by realizing to Stage::Acceleration).

This is an O(m) operator. In particular it does \e not involve forming or
multiplying by the constraint force matrix ~G. Instead, one constant-time call
is made to each %Constraint's calcConstraintForce methods.
    
@par Required stage
  \c Stage::Velocity **/
void calcConstraintForcesFromMultipliers
  (const State&         state, 
   const Vector&        multipliers,
   Vector_<SpatialVec>& bodyForcesInG,
   Vector&              mobilityForces) const;

/** Calculate the mobilizer reaction force generated at each MobilizedBody,
as felt at the mobilizer's outboard frame M, and expressed in Ground.

@param[in] state        
    A State compatible with this System that has already been realized 
    to Stage::Acceleration.
@param[out] forcesAtMInG    
    A Vector of spatial force vectors, indexed by MobilizedBodyIndex 
    (beginning with 0 for Ground), giving the reaction moment and force
    applied by each body's unique inboard mobilizer to that body. The
    force is returned as though it were applied at the origin of the 
    body's mobilizer frame M. The returned force is expressed in the
    Ground frame. Applied mobility (generalized) forces are \e included in the
    returned reaction forces.

A simple way to think of the reaction force is to think of cutting the 
mobilizer, then imagine the force required to make the system move in 
the same manner as when the mobilizer was present. This is what the 
reaction forces accomplish. With that definition, mobility forces (that is,
generalized forces as opposed to body forces) are \e included in the reactions.
Some conventions do not include the mobility forces in the definition of a 
reaction force. We chose to include them since this preserves Newton's 
3rd law of equal and opposite reactions between bodies. Ours is the same 
convention as used in SD/FAST.

@note You can think of the Ground body being welded to the universe at the
Ground origin. The reactions reported for Ground are the ones that would 
occur in that Weld mobilizer if it were really present. That is, it includes
the effects of all the base bodies on Ground.

<h3>How to find the reaction felt by the parent body</h3>

A mobilizer connects a frame F fixed on the parent (inboard) body P to a
frame M fixed on the child (outboard) body B. It exerts equal and opposite 
reaction forces on the two bodies, at a given location in space. This method 
reports the force on the child body, as though it were applied at the origin 
Mo of frame M, and expressed in the Ground frame. The force on the parent body
<em>at Mo</em> is just the negative of the returned value. However, it is 
more likely that you would want it as felt <em>at Fo</em>, the origin of the
F frame on the parent. Here is one way to calculate that from the returned
quantities:
@code
    matter.calcMobilizerReactionForces(state,forcesAtMInG); // This method.
    const int nb = matter.getNumBodies();
    Vector_<SpatialVec> forcesAtFInG(nb); // to hold the result
    forcesAtFInG[0] = -forcesAtMInG[0]; // Ground is "welded" at origin
    for (MobilizedBodyIndex i(1); i < nb; ++i) {
        const MobilizedBody& body   = matter.getMobilizedBody(i);
        const MobilizedBody& parent = body.getParentMobilizedBody();
        // Want to shift reaction by p_MF, the vector from M to F across the
        // mobilizer, and negate. Can get p_FM; must reexpress in G.
        const Vec3& p_FM = body.getMobilizerTransform(state).p();
        const Rotation& R_PF = body.getInboardFrame(state).R(); // In parent.
        const Rotation& R_GP = parent.getBodyTransform(state).R();
        Rotation R_GF   =   R_GP*R_PF;  // F frame orientation in Ground.
        Vec3     p_MF_G = -(R_GF*p_FM); // Re-express and negate shift vector. 
        forcesAtFInG[i] = -shiftForceBy(forcesAtMInG[i], p_MF_G);
    }
@endcode

<h3>Implementation</h3>
This method combines already-calculated quantities to calculate the reactions.
See Abhi Jain's 2011 book "Robot and Multibody Dynamics", Eq. 7.34 page 128: 
<pre>   F_reaction = PPlus*APlus + zPlus  </pre>
where P is the articulated body inertia, A is the spatial acceleration,
a the Coriolis acceleration and z the articulated body forces, and "Plus"
indicates that we evaluate these on the inboard (parent) side of the mobilizer
rather than on the body's side. (The alternative P(A-a)+z given there does not 
work for prescribed mobilizers unless you replace "a" with "a_underscore" from
equation 16.14.) After calculating F_reaction at the body frame
origin Bo, we shift it to M for reporting.

<h3>Performance</h3>
The cost of the above calculation is 114 flops/body. The code presented
above for converting from M to F costs an additional 81 flops/body if you
use it.
    
@par Required stage
  \c Stage::Acceleration 
 
@see SimTK::MobilizedBody::findMobilizerReactionOnBodyAtMInGround()
@see calcMobilizerReactionForcesUsingFreebodyMethod() **/
void calcMobilizerReactionForces
   (const State&         state, 
    Vector_<SpatialVec>& forcesAtMInG) const;

/** Return a reference to the prescribed motion multipliers tau that have 
already been calculated in the given \a state, which must have been realized 
through Acceleration stage. The result contains entries only for prescribed 
mobilities; if you want these unpacked into u-space mobility forces, use
findMotionForces() instead. A mobilizer may follow prescribed motion either
because of a Motion object or a call to MobilizedBody::lock(). **/
const Vector& getMotionMultipliers(const State& state) const;

/** Calculate the degree to which the supplied \a state does not satisfy the
prescribed motion requirements at a particular Stage. For Position and Velocity
stage, a call to the prescribe() solver using the same stage will eliminate
the error. Accelerations should have been calculated to satisfy all prescribed
accelerations, so the returned value should be zero always. The returned 
Vector has one element per known (prescribed) q, known u, or known udot. 

The \a state must be realized to Time stage to check Position errors,
Position stage to check Velocity errors, and Acceleration stage to check
Acceleration errors. 

Errors are calculated actualValue - prescribedValue so a positive error
indicates that the value in \a state is too large. **/
Vector calcMotionErrors(const State& state, const Stage& stage) const;

/** Find the generalized mobility space forces produced by all the Motion 
objects active in this system. These are the same values as returned by 
getMotionMultipliers() but unpacked into u-space slots, with zeroes 
corresponding to any "free" mobilities, that is, those whose motion is not 
prescribed. **/
void findMotionForces
   (const State&         state,
    Vector&              mobilityForces) const;

/** Return a reference to the constraint multipliers lambda that have already 
been calculated in the given \a state, which must have been realized through 
Acceleration stage. Constraint multipliers are not directly interpretable as 
forces; if you want the actual forces use findConstraintForces() instead. If
you want to know individual Constraint contributions to these forces, ask the 
Constraint objects rather than this SimbodyMatterSubsystem object. **/
const Vector& getConstraintMultipliers(const State& state) const;

/** Find the forces produced by all the active Constraint objects in this 
system. Constraints produce both body spatial forces and generalized 
mobility-space forces. The supplied \a state must have been realized through 
Acceleration stage. **/
void findConstraintForces
  (const State&         state, 
   Vector_<SpatialVec>& bodyForcesInG,
   Vector&              mobilityForces) const;

/** Calculate the power being generated or dissipated by all the Motion objects
currently active in this system. The sign is chosen so that a positive value for
power means the Motion is adding energy to the system; negative means it is
removing energy. The \a state must already have been realized through 
Acceleration stage so that the prescribed motion forces are available.

@param[in]      state
    A State realized through Acceleration stage from which we obtain the
    prescribed motion forces and the velocities needed to calculate power.

<h3>Implementation</h3>
We calculate power=-dot(tau, u) where tau is the set of mobility reaction 
forces generated by Motion objects and mobilizer locks (tau[i]==0 if mobility
i is free), and u is the set of all generalized speeds. 
@see calcConstraintPower() **/
Real calcMotionPower(const State& state) const;

/** Return the power begin generated or dissipated by all the Constraint
objects currently active in this system. The sign is chosen so that a positive 
value for power means the Constraints (taken together) are adding energy to the
system; negative means they are removing energy. The \a state must already have
been realized through Acceleration stage so that the constraint forces are 
available.

Note that if you want to know the power output of an individual Constraint,
you should call that Constraint's calcPower() method; here they are all summed
together.

@param[in]      state
    A State realized through Acceleration stage from which we obtain the
    constraint forces and the velocities needed to calculate power.
@return
    The signed sum over all the Constraint objects of the power being generated
    or dissipated by each Constraint. A positive value means that together the
    constraints are adding energy to the system; negative means they are 
    removing energy.

<h3>Implementation</h3>
We calculate power=-(dot(F,V)+dot(f,u)) where F is the set of body spatial
reaction forces produced by the Constraints, V is the body spatial velocities, 
f is the set of mobility reaction forces produced by the Constraints, and u is 
the set of generalized speeds. 
@see calcMotionPower() **/
Real calcConstraintPower(const State& state) const;

/** Accounts for applied forces and inertial forces produced by non-zero 
velocities in the State. Returns a set of mobility forces which replace both 
the applied bodyForces and the inertial forces.
@par Required stage
  \c Stage::Dynamics **/
void calcTreeEquivalentMobilityForces(const State&, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobilityForces) const;


/** Calculate qdot = N(q)*u in O(n) time (very fast). Note that q is taken
from the supplied state while u is an argument to this operator method.
@par Required stage
  \c Stage::Position **/
void calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const;

/** Calculate qdotdot = N(q)*udot + Ndot(q,u)*u in O(n) time (very fast). Note
that q and u are taken from the supplied state while udot is an argument to
this operator method.
@par Required stage
  \c Stage::Velocity **/
void calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const;

/** Add in to the given body forces vector a force applied to a station (fixed
point) S on a body B. The new force is added into the existing spatial force 
slot for the body. Note that this does not actually apply any forces to the
multibody system! This is just a "helper" utility that makes it easier to 
fill in a body forces array. This has no effect on the system unless you
later supply the body forces array for use.

Provide the station in the body frame, force in the Ground frame. 
@par Required stage
  \c Stage::Position **/
void addInStationForce(const State&         state, 
                       MobilizedBodyIndex   bodyB, 
                       const Vec3&          stationOnB, 
                       const Vec3&          forceInG, 
                       Vector_<SpatialVec>& bodyForcesInG) const;

/** Add in to the given body forces vector a torque applied to a body B. The 
new torque is added into the existing spatial force slot for the body. Note 
that this does not actually apply any torques to the multibody system! This is 
just a "helper" utility that makes it easier to fill in a body forces array. 
This has no effect on the system unless you later supply the body forces array 
for use. Provide the torque vector in the Ground frame. **/
void addInBodyTorque(const State&           state, 
                     MobilizedBodyIndex     mobodIx, 
                     const Vec3&            torqueInG, 
                     Vector_<SpatialVec>&   bodyForcesInG) const;

/** Add in to the given mobility forces vector a scalar generalized force, that
is a force or torque applied to a mobilizer generalized speed. Note 
that this does not actually apply any forces to the multibody system! This is 
just a "helper" utility that makes it easier to fill in a mobility forces array. 
This has no effect on the system unless you later supply the mobility forces
array for use. The meaning of a generalized force f is determined by the
definition of the corresponding generalized speed u, so that f*u has physical
units of power. **/
void addInMobilityForce(const State&        state,
                        MobilizedBodyIndex  mobodIx, 
                        MobilizerUIndex     which, 
                        Real                f,
                        Vector&             mobilityForces) const;
/**@}**/



//==============================================================================
/** @name              Realization and response methods

Realization methods request that some calculation be performed ("realized") if
it has not already been done since the last change to one of the state 
variables on which the result depends, with the result being placed in the
state cache. Methods beginning with "get" are called \e responses and are used
to extract pre-calculated information that has been realized into the cache.

Realization is normally initiated at the System level. However, there are some
"lazy" calculations in the %SimbodyMatterSubsystem whose computations are
delayed until needed; you can cause those calculations to be performed 
explicitly here if you want. **/
/**@{**/

    // POSITION STAGE realizations //

/** This method checks whether composite body inertias have already been 
computed since the last change to a Position stage state variable (q) and if so 
returns immediately at little cost; otherwise, it initiates computation of 
composite body inertias for all of the mobilized bodies. These are not otherwise
computed unless specifically requested.
@par Required stage
  \c Stage::Position **/
void realizeCompositeBodyInertias(const State&) const;

/** This method checks whether articulated body inertias have already been 
computed since the last change to a Position stage state variable (q) and if so 
returns immediately at little cost; otherwise, it initiates the relatively 
expensive computation of articulated body inertias for all of the mobilized 
bodies. These are not otherwise computed until they are needed at Dynamics 
stage.
@par Required stage
  \c Stage::Position **/
void realizeArticulatedBodyInertias(const State&) const;


    // INSTANCE STAGE responses //

const Array_<QIndex>& getFreeQIndex(const State& state) const;
const Array_<UIndex>& getFreeUIndex(const State& state) const;
const Array_<UIndex>& getFreeUDotIndex(const State& state) const;
const Array_<UIndex>& getKnownUDotIndex(const State& state) const;
void packFreeQ
   (const State& s, const Vector& allQ, Vector& packedFreeQ) const;
void unpackFreeQ
   (const State& s, const Vector& packedFreeQ, Vector& unpackedFreeQ) const;
void packFreeU
   (const State& s, const Vector& allU, Vector& packedFreeU) const;
void unpackFreeU
   (const State& s, const Vector& packedFreeU, Vector& unpackedFreeU) const;


    // POSITION STAGE responses //

/** Return the composite body inertia for a particular mobilized body. You can 
call this any time after the State has been realized to Position stage, however
it will first trigger realization of all the composite body inertias if they 
have not already been calculated. Ground is mobilized body zero; its composite
body inertia has infinite mass and principle moments of inertia, and zero 
center of mass.
@par Required stage
  \c Stage::Position
@see realizeCompositeBodyInertias() **/
const SpatialInertia& 
getCompositeBodyInertia(const State& state, MobilizedBodyIndex mbx) const;

/** Return the articulated body inertia for a particular mobilized body. You
can call this any time after the State has been realized to Position stage, 
however it will first trigger expensive realization of all the articulated body
inertias if they have not already been calculated. Ground is mobilized body 
zero; its articulated body inertia is the same as its composite body inertia --
an ordinary Spatial Inertia but with infinite mass and principle moments of
inertia, and zero center of mass.
@par Required stage
  \c Stage::Position
@see realizeArticulatedBodyInertias() **/
const ArticulatedInertia& 
getArticulatedBodyInertia(const State& state, MobilizedBodyIndex mbx) const;


    // VELOCITY STAGE responses //

/** This is the angular velocity-dependent force on the body due to rotational 
inertia.
@par Required stage
  \c Stage::Velocity **/
const SpatialVec& 
getGyroscopicForce(const State& state, MobilizedBodyIndex mbx) const;

/** This is the cross-mobilizer incremental contribution to coriolis (angular 
velocity dependent) acceleration; not too useful, see 
getTotalCoriolisAcceleration() instead.
@par Required stage
  \c Stage::Velocity **/
const SpatialVec& 
getMobilizerCoriolisAcceleration(const State&       state, 
                                 MobilizedBodyIndex mbx) const;

/** This is the total coriolis acceleration including the effect of the parent's
angular velocity as well as the joint's. This is Jdot*u where J is the system
kinematic Jacobian and u is the current set of generalized speeds in the 
supplied state. It is thus the remainder term in calculation of body 
accelerations from generalized accelerations udot: since V=J*u, 
A=J*udot + Jdot*u.
@par Required stage
  \c Stage::Velocity **/
const SpatialVec& 
getTotalCoriolisAcceleration(const State& state, MobilizedBodyIndex mbx) const;


    // DYNAMICS STAGE responses //

/** This is the angular velocity-dependent force accounting for gyroscopic 
forces plus coriolis forces due only to the cross-mobilizer velocity; this 
ignores the parent's velocity and is not too useful -- see 
getTotalCentrifugalForces() instead.
@par Required stage
  \c Stage::Dynamics **/
const SpatialVec& 
getMobilizerCentrifugalForces(const State& state, MobilizedBodyIndex mbx) const;

/** This is the total angular velocity-dependent force acting on this body, 
including forces due to coriolis acceleration and forces due to rotational
inertia.
@par Required stage
  \c Stage::Dynamics **/
const SpatialVec& 
getTotalCentrifugalForces(const State& state, MobilizedBodyIndex mbx) const;
/**@}**/



//==============================================================================
/** @name              Testing and debugging utilities

Methods in this section provide alternate ways of calculating quantities for
which we provide more efficient methods above. You should use the better 
methods normally, but these can be very useful for regression testing and
Simbody development because the answers are obtained differently. Numerical
results should agree with the faster methods to within numerical precision. **/
/**@{**/

/** This is a slower alternative to calcMobilizerReactionForces(), for use in
regression testing and Simbody development. This method builds a freebody
"diagram" for each body in turn to determine the unknown reaction force at
its inboard mobilizer.

The given \a state must have been realized through Acceleration stage.

<h3>Implementation</h3>
We use a tip-to-base recursion in which we assemble all applied body forces
from force elements, constraints, and gyroscopic effects and compare that with
the apparent rigid body force determined from F=M*A where M is a body's 
spatial inertia (in G) and A its already-calculated spatial acceleration. The 
difference is the missing force applied by the body's mobilizer, i.e. the 
reaction force. That is shifted to the M frame and reported. Then the equal
and opposite reaction is applied to the parent body and included in its 
collection of applied forces, which can be used to determine its reaction force
and so on.

This is is about 3X slower than the method used by 
calcMobilizerReactionForces(). 
@see calcMobilizerReactionForces() **/
void calcMobilizerReactionForcesUsingFreebodyMethod
   (const State&         state, 
    Vector_<SpatialVec>& forcesAtMInG) const;
/**@}**/


//==============================================================================
/** @name Proposed particle API

(NOT IMPLEMENTED YET) These methods are a proposed API for explicit handling
of particles. Currently a particle should be implemented as point mass with a 
Cartesian (translation) mobilizer to Ground instead. The idea here would be to
special-case particles to make them faster; there would be no additional 
functionality. **/

/**@{**/

/// TODO: total number of particles.
int getNumParticles() const;

// The generalized coordinates for a particle are always the three measure numbers
// (x,y,z) of the particle's Ground-relative Cartesian location vector. The generalized
// speeds are always the three corresponding measure numbers of the particle's
// Ground-relative Cartesian velocity. The generalized applied forces are
// always the three measure numbers of a Ground-relative force vector.
const Vector_<Vec3>& getAllParticleLocations    (const State&) const;
const Vector_<Vec3>& getAllParticleVelocities   (const State&) const;

const Vec3& getParticleLocation(const State& s, ParticleIndex p) const {
    return getAllParticleLocations(s)[p];
}
const Vec3& getParticleVelocity(const State& s, ParticleIndex p) const {
    return getAllParticleVelocities(s)[p];
}

Vector& updAllParticleMasses(State& s) const;

void setAllParticleMasses(State& s, const Vector& masses) const {
    updAllParticleMasses(s) = masses;
}

// Note that particle generalized coordinates, speeds, and applied forces
// are defined to be the particle Cartesian locations, velocities, and
// applied force vectors, so can be set directly at Stage::Model or higher.

// These are the only routines that must be provided by the concrete MatterSubsystem.
Vector_<Vec3>& updAllParticleLocations(State&)     const;
Vector_<Vec3>& updAllParticleVelocities(State&)    const;

// The following inline routines are provided by the generic MatterSubsystem 
// class for convenience.

Vec3& updParticleLocation(State& s, ParticleIndex p) const {
    return updAllParticleLocations(s)[p];
}
Vec3& updParticleVelocity(State& s, ParticleIndex p) const {
    return updAllParticleVelocities(s)[p];
}

void setParticleLocation(State& s, ParticleIndex p, const Vec3& r) const {
    updAllParticleLocations(s)[p] = r;
}
void setParticleVelocity(State& s, ParticleIndex p, const Vec3& v) const {
    updAllParticleVelocities(s)[p] = v;
}

void setAllParticleLocations(State& s, const Vector_<Vec3>& r) const {
    updAllParticleLocations(s) = r;
}
void setAllParticleVelocities(State& s, const Vector_<Vec3>& v) const {
    updAllParticleVelocities(s) = v;
}

const Vector& getAllParticleMasses(const State&) const;

const Vector_<Vec3>& getAllParticleAccelerations(const State&) const;

const Vec3& getParticleAcceleration(const State& s, ParticleIndex p) const {
    return getAllParticleAccelerations(s)[p];
}
/**@}**/

//==============================================================================
/** @name                      Obsolete methods

These methods are deprecated because there is a better way now to do what they
used to do. This may involve just a name change, calling signature, or something
more substantial; see the documentation for the individual obsolete methods.
**/

/**@{**/
private:
/** Obsolete synonym for multiplyBySystemJacobian().
@see multiplyBySystemJacobian() **/
void calcSpatialKinematicsFromInternal(const State&         state,
                                       const Vector&        u,
                                       Vector_<SpatialVec>& Ju) const
{   multiplyBySystemJacobian(state,u,Ju); }

/** Obsolete synonym for multiplyBySystemJacobianTranspose().
@see multiplyBySystemJacobianTranspose() **/
void calcInternalGradientFromSpatial(const State&               state,
                                     const Vector_<SpatialVec>& F_G,
                                     Vector&                    f) const
{   multiplyBySystemJacobianTranspose(state, F_G, f); }

/** Obsolete synonym for multiplyByM().
@see multiplyByM() **/
void calcMV(const State& state, const Vector& v, Vector& MV) const
{   multiplyByM(state,v,MV); }

/** Obsolete synonym for multiplyByMInv().
@see multiplyByMInv() **/
void calcMInverseV(const State& state,
    const Vector&               v,
    Vector&                     MinvV) const
{   multiplyByMInv(state,v,MinvV); }

/** Obsolete slow version of calcPq.
@see calcPq() **/
void calcPNInv(const State& state, Matrix& PNInv) const;

/** Obsolete slow version of calcGTranspose().
@see calcGTranspose() **/
void calcGt(const State&, Matrix& Gt) const;


/**@}**/


//==============================================================================
//     Bookkeeping methods and internal types -- hide from Doxygen
/** @cond **/
public:
class Subtree; // used for working with a connected subgraph of the MobilizedBody tree
class SubtreeResults;


SimTK_PIMPL_DOWNCAST(SimbodyMatterSubsystem, Subsystem);
const SimbodyMatterSubsystemRep& getRep() const;
SimbodyMatterSubsystemRep&       updRep();
/** @endcond **/

private:
};

/** Dump some debug information about the given subsystem to the given
output stream; this is \e not for serialization.
@relates SimbodyMatterSubsystem **/
SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem&);


} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_H_
