#ifndef SimTK_SIMBODY_CONSTRAINT_LINE_ON_LINE_CONTACT_H_
#define SimTK_SIMBODY_CONSTRAINT_LINE_ON_LINE_CONTACT_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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
Declares the Constraint::LineOnLineContact class. **/

#include "simbody/internal/Constraint.h"

namespace SimTK {

//==============================================================================
//                         LINE ON LINE CONTACT
//==============================================================================
/** This constraint represents a \e bilateral connection between an edge
on one body and a non-parallel edge on another.

On construction you may choose whether the connection enforces non-slipping
(which we'll call "rolling" here for consistency with other constraints);
otherwise the lines will slip against each other. There is always one position 
(holonomic) constraint equation: the lines containing the edges must touch at 
some point. That leaves the edges with five unconstrained degrees of freedom 
(dofs) to take on any relative orientation and to be touching at any point on 
their lengths. Optionally, there are also two velocity (nonholonomic) constraint 
equations that prevent relative slip between the edges. In that case there can
be five position dofs but only three unconstrained velocity dofs.

Note that this is a bilateral, unconditional connection and will push or pull 
as necessary to keep the lines in contact. If rolling is being enforced then 
whatever tangential forces are necessary to keep the lines from slipping 
will be generated, regardless of the normal force. These constraints can form 
the basis for unilateral edge-edge contact with Coulomb friction, but additional
conditions must be added to determine when they are active.

There are two mobilized bodies involved, we'll call them F and B. Either or both
bodies can move or either can be Ground; however, we will orient the signs in 
our calculations so that we can think of F as "fixed"; the contact normal
points towards the "outside" of F and towards the "inside" of B. That means
that if we translate B relative to F in the direction of the normal, the
signed distance between them increases (separation increases or penetration
decreases). That allows us to say, for example, that a positive position error 
means separation, that a negative normal velocity error means impact, and that
a positive force (negative multiplier) means compression.

Each edge E is defined by a center point P, edge direction d, and exterior 
"space" direction s used to define the sign convention. The s direction 
points outward from the solid whose two faces meet to form E, along the
midplane between those two faces. We also ask for the edge half-lengths
which can be used for visualization and for detecting separation caused by
slipping of the end of an edge, but don't actually affect the constraint here
which works on the lines containing the edges.

F has an edge Ef attached to it with center point Pf, direction df, half length
hf, and outward normal sf. B has edge Eb=(Pb,db,hb,sb). The vector between the 
defining center points is p_PfPb=Pb-Pf, oriented from Ef's center point to Bf's 
center point. The contact normal n will be a unit vector aligned with 
w=df X db, with n=sign*w/||w||. Sign is chosen by examining the dot product
of w with sf and sb so that translating Eb along +n would increase the signed
distance between the edges.

To locate the contact point Co, we want to find the closest points Qf and Qb
on each line, with Qf=Pf + tf*df and Qb=Pb + tb*db for some parameters tf and 
tb. So we are looking for the solution of this equation:
<pre>
    (1) Pf + tf*df + r*n = Pb + tb*db
</pre>
where r is the separation distance and n is the contact normal defined above. 
These are three linear equations in three unknowns, saying that we can get from 
the closest point Qf on Ef to the closest point Qb on Eb by moving up the 
contact normal direction from Qf to Qb by a signed distance r (that is, if r<0
then we are moving \e down the contact normal). Since n is perpendicular to both
df and db, we can dot both sides of Eqn. (1) with n to get
<pre>
    (2) (Pf + r*n) . n = Pb . n
==>     r*(n.n) = (Pb - Pf) . n
==> (3) r = (Pb - Pf) . n           [since n is a unit vector]
</pre>
Similarly, we can dot with n X db and n X df to solve for tf and tb, resp.
<pre>
    (4) (Pf + tf*df) . (n X db) = Pb . (n X db)
==>     tf * (df . (n X db)) = (Pb - Pf) . (n X db)

             (Pb-Pf).(n X db)
==> (5) tf = ----------------       [can use wXdb here instead of nXdb]
              df . (n X db)
</pre> and <pre>
    (6) Pf . (n X df) = (Pb + tb*db) . (n X df)
==>     tb * (db . (n X df)) = (Pf - Pb) . (n X df)

               (Pb-Pf).(n X df)
==> (7) tb = - ----------------     [can use wXdf here instead of nXdf]
                db . (n X df)
</pre>

Then the position constraint we want to enforce is r==0, requiring the lines 
to touch at their closest points. If the constraint were perfectly enforced,
then the contact point Co would be at the same location as the two
closest points. Since it won't be enforced perfectly, we'll put the contact
point at the midpoint between Qf and Qb, so Co=(Qf+Qb)/2. The exact position
of Co does not matter for the position constraint, because the normal force
is applied along the line including Qf and Qb. However, it will matter (a 
little) for the friction constraints.

For any given pose, we will define a contact frame C with origin Co whose z
axis Cz=n, Cx=df, Cy=n X df. The instantaneous contact frame in Ground is 
available as the transform X_GC.

The contact constraints here are enforced by a normal multiplier acting along
Cz, and optionally two tangential multipliers acting along Cx and Cy
respectively. Together these three multipliers can be interpreted as a force 
expressed in frame C, and applied equal and opposite to bodies F and B at
their material points coincident with contact point Co.

The assembly condition is the same as the position constraint: the two lines
must meet at some point. No attempt is made to force the contact point to be
within the two edges; we only make the lines containing the edges touch
somewhere. There is no assembly
condition for the tangential constraints since they do not restrict the
allowable pose during assembly. **/
class SimTK_SIMBODY_EXPORT Constraint::LineOnLineContact
:   public Constraint  {
public:
    // no default constructor

//------------------------------------------------------------------------------
/** @name                       Construction
Methods in this section refer both to constructors, and to methods that can
be used to set or change contruction (Topology-stage) parameters; these
specify the values assigned by default to the corresponding state variables. 
Note:
  - Changing one of these default parameters invalidates the containing 
    System's topology, meaning that realizeTopology() will have to be called 
    and a new State obtained before subsequent use. 
  - The set...() methods return a reference to "this" %LineOnLineContact
    element (in the manner of an assignment operator) so they can be chained in
    a single expression. 
 
The default parameters can be overridden in any given State, and modified
without affecting Topology; see the "Runtime Changes" section below. **/
/*@{*/

/** Construct a line-on-line constraint as described in the 
Constraint::LineOnLineContact class documentation.

@param      mobod_F   
    The first MobilizedBody object to which a contacting edge is attached.
    We'll call it F, the "fixed" body just to orient the contact; actually
    either or both bodies can be moving.
@param      defaultEdgeFrameF
    This Transform X_FEf defines the location and direction of the edge Ef 
    on \a mobod_F, measured from and expressed in F. The frame origin is the 
    center point Pf of the edge. Its x axis is direction df aligned with the 
    edge; y is unused; z is Ef's outward ("space") direction sf pointing away 
    from the polygonal solid for which Ef is an edge, midway between the two 
    faces whose intersection defines the edge. This parameter provides the edge 
    frame that will be present in a default State; you can modify it later in 
    any particular State.
@param      defaultHalfLengthF
    This is the half-length hf of edge Ef. The line segment representing the
    edge thus runs from Pf-hf*df to Pf+hf*df. This parameter provides the
    half-length that will be present in a default State; you can modify it later
    in any particular State.
@param      mobod_B   
    The second MobilizedBody object to which a contacting edge is attached.
@param      defaultEdgeFrameB
    This Transform X_BEb defines the location and direction of the edge Eb 
    on \a mobod_B, measured from and expressed in B. The frame origin is the 
    center point Pb of the edge. Its x axis is direction db aligned with the 
    edge; y is unused; z is Eb's outward ("space") direction sb pointing away 
    from the polygonal solid for which Eb is an edge, midway between the two 
    faces whose intersection defines the edge. This parameter provides the edge 
    frame that will be present in a default State; you can modify it later in 
    any particular State.
@param      defaultHalfLengthB
    This is the half-length hb of edge Eb. The line segment representing the
    edge thus runs from Pb-hb*db to Pb+hb*db. This parameter provides the
    half-length that will be present in a default State; you can modify it later
    in any particular State.
@param      enforceRolling
    Whether to generate tangential forces to prevent the lines from slipping
    against one another. Otherwise only a normal force is generated and the 
    lines are free to slip.
**/
LineOnLineContact(  MobilizedBody&      mobod_F, 
                    const Transform&    defaultEdgeFrameF, 
                    Real                defaultHalfLengthF, 
                    MobilizedBody&      mobod_B, 
                    const Transform&    defaultEdgeFrameB, 
                    Real                defaultHalfLengthB, 
                    bool                enforceRolling);
    
/** Default constructor creates an empty handle that can be used to
reference any %LineOnLineContact %Constraint.**/
LineOnLineContact() {}

/** Return a reference to the first MobilizedBody to which a line is
attached. This refers to the \a mobod_F that was given in the 
constructor and cannot be changed after construction. **/
const MobilizedBody& getMobilizedBodyF() const;

/** Return a reference to the second MobilizedBody to which a line is
attached. This refers to the \a mobod_B that was given in the 
constructor and cannot be changed after construction. **/
const MobilizedBody& getMobilizedBodyB() const;

/** Report whether this %Constraint was constructed to generate rolling
constraints (otherwise it is frictionless). This cannot be changed after
construction. **/
bool isEnforcingRolling() const;

/** Replace the default frame of the edge attached to the first body,
\a mobod_F, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setEdgeFrameF() **/
LineOnLineContact& 
setDefaultEdgeFrameF(const Transform& defaultEdgeFrameF);
/** Replace the default half-length for the edge attached to the first body,
\a mobod_F, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setHalfLengthF() **/
LineOnLineContact& 
setDefaultHalfLengthF(Real defaultHalfLengthF);
/** Replace the default frame of the edge attached to the second body,
\a mobod_B, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setEdgeFrameB() **/
LineOnLineContact& 
setDefaultEdgeFrameB(const Transform& defaultEdgeFrameB);
/** Replace the default half-length for the edge attached to the second body,
\a mobod_B, that was supplied on construction. This is a topological change;
you'll have to call realizeTopology() again if you call this method. 
@see setHalfLengthB() **/
LineOnLineContact& 
setDefaultHalfLengthB(Real defaultHalfLengthF);

/** Return the default frame of the edge attached to the first
body, \a mobod_F, as set during construction or by the most recent call to 
setDefaultEdgeFrameF(). 
@see getEdgeFrameF() **/
const Transform& getDefaultEdgeFrameF() const;
/** Return the default half-length for the edge attached to the first body, 
\a mobod_F, as set during construction or by the most recent call to 
setDefaultHalfLengthF(). 
@see getHalfLengthF() **/
Real getDefaultHalfLengthF() const;

/** Return the default frame of the edge attached to the second
body, \a mobod_B, as set during construction or by the most recent call to 
setDefaultEdgeFrameB(). 
@see getEdgeFrameB() **/
const Transform& getDefaultEdgeFrameB() const;
/** Return the default half-length for the edge attached to the second body, 
\a mobod_B, as set during construction or by the most recent call to 
setDefaultHalfLengthB(). 
@see getHalfLengthB() **/
Real getDefaultHalfLengthB() const;
/*@}............................ Construction ................................*/

//------------------------------------------------------------------------------
/** @name                       Visualization
If you are allowing the SimbodyMatterSubsystem to generate default geometry,
this constraint element will attempt to represent the lines that
it connects. Methods here give you some limited control over the generated
default geometry; if you need more control you should disable the generation
of default geometry and make your own. **/
/*@{*/
// no controls yet
/*@}........................... Visualization ................................*/


//------------------------------------------------------------------------------
/** @name                      Runtime Changes
These refer to Position-stage discrete state variables that determine the line
parameters to be used to calculate constraint forces from a given State object.
If these are not set explicitly, the parameters are set to those provided in the
constructor or via the correponding setDefault...() methods. 
Note:
  - Changing one of these parameters invalidates the given State's 
    Stage::Position, meaning that the State's stage will be no higher than
    Stage::Time after the parameter change. That ensures that position 
    constraint errors will be recalculated before they are used.
  - The set...() methods here return a const reference to "this" 
    %LineOnLineContact element (in the manner of an assignment operator, 
    except read-only) so they can be chained in a single expression. 
    
You can also modify and examine the default parameters; see the "Construction"
section above. **/
/*@{*/

/** Modify the frame of the edge on the first body, \a mobod_F, in this 
\a state by providing a new Transform X_FEf measured from and expressed in the 
F frame. The origin is the location of the edge center point Pf; the x axis
is the edge direction df; z is the outward direction sf; y is unused. This 
overrides the \a defaultEdgeFrameF in the given \a state, whose Stage::Position 
is invalidated. **/
const LineOnLineContact& 
setEdgeFrameF(State& state, const Transform& edgeFrameF) const;

/** Modify the half-length hf of the edge on the first body, \a mobod_F, in this
\a state. This overrides the \a defaultHalfLengthF in the given \a state, whose 
Stage::Position is invalidated. **/
const LineOnLineContact& 
setHalfLengthF(State& state, Real halfLengthF) const;

/** Modify the frame of the edge on the second body, \a mobod_B, in this 
\a state by providing a new Transform X_BEb measured from and expressed in the 
B frame. The origin is the location of the edge center point Pb; the x axis
is the edge direction db; z is the outward direction sb; y is unused. This 
overrides the \a defaultEdgeFrameB in the given \a state, whose Stage::Position 
is invalidated. **/
const LineOnLineContact& 
setEdgeFrameB(State& state, const Transform& edgeFrameB) const;

/** Modify the half-length hb of the edge on the second body, \a mobod_B, in 
this \a state. This overrides the \a defaultHalfLengthB in the given \a state, 
whose Stage::Position is invalidated. **/
const LineOnLineContact& 
setHalfLengthB(State& state, Real halfLengthB) const;

/** Return the frame of the edge Ef on the first body, \a mobod_F, as currently
set in the given \a state. The value is returned as a Transform X_FEf measured
from and expressed in the F frame. **/
const Transform& getEdgeFrameF(const State& state) const;
/** Return the half-length of the edge Ef on the first body, 
\a mobod_F, as currently set in the given \a state. The returned value is
the scalar hf. **/
Real getHalfLengthF(const State& state) const;

/** Return the frame of the edge Eb on the second body, \a mobod_B, as currently
set in the given \a state. The value is returned as a Transform X_BEb measured
from and expressed in the B frame. **/
const Transform& getEdgeFrameB(const State& state) const;
/** Return the half-length of the edge Eb on the second body, 
\a mobod_B, as currently set in the given \a state. The returned value is
the scalar hb. **/
Real getHalfLengthB(const State& state) const;
/*@}.......................... Runtime Changes ...............................*/

//------------------------------------------------------------------------------
/** @name                      Computations
Methods here provide access to values already calculated by this constraint
element, and provide operators you can call to calculate related values. **/
/*@{*/
/** The returned position error can be viewed as the signed distance between
the lines. It is positive when the lines are separated and measures their
closest approach distance. It is negative when the spheres are interpenetrating
and measures the penetration depth. The given \a state must have already been 
realized through Stage::Position. **/
Real getPositionError(const State& state) const;

/** The returned velocity error vector has the time derivative of the quantity
returned by getPositionError() in its z coordinate, and violation of the
rolling constraints in its x and y coordinates. If rolling is not being 
enforced then the x and y components are returned zero; they will not contain
the slip velocity in that case since any slip velocity is acceptable. Note
that the returned vector is expressed in the instantaneous contact frame C,
considered as fixed on body F. That is, this is the velocity of the contact
point on body B's sphere in the F frame, expressed in C. The given \a state must
have already been realized through Stage::Velocity. **/
Vec3 getVelocityErrors(const State& state) const;

/** This vector is the time derivative of the value returned by
getVelocityError(). Note that this is different than the acceleration of
the point of sphere B at the contact point because the contact point moves
with respect to that sphere. The given \a state must have already been realized
through Stage::Acceleration. **/
Vec3 getAccelerationErrors(const State& state) const;

/** These are the Lagrange multipliers required to enforce the constraint
equations generated here. For this %Constraint it has units of force, but 
recall that the sign convention for multipliers is the opposite of that for 
applied forces. Thus the returned value is the negative of the force being
applied to sphere B at the contact point, expressed in the contact frame C.
The x,y coordinates are the tangential force used to enforce rolling (or
zero if rolling is not being enforced), and the z coordinate is the force 
needed to enforce contact. Since this is an unconditional, bilateral 
constraint the multipliers may have any sign and magnitude. The given \a state 
must already be realized to Stage::Acceleration. **/
Vec3 getMultipliers(const State& state) const;

/** Return the force vector currently being applied by this constraint to the
point of body B that is coincident in space with the contact point Co. An equal
and opposite force is applied to body F at the same location. This is zero if 
the constraint is not currently enabled. Tangential forces are generated only 
when rolling is being enforced, but since this result is in the Ground frame all
three vector measure numbers may be non-zero regardless. The given \a state must
already be realized to Stage::Acceleration. **/
Vec3 findForceOnBodyBInG(const State& state) const;

/** Return the instantaneous contact frame C in the Ground frame. The actual
contact point is the origin Co of the returned frame, which is placed at a 
point midway along the line segment connecting the closest points of the two 
lines. This does not imply that the normal constraint is satisifed; point Co 
could be far from the edges. The z direction of this frame is the
contact normal; it is perpendicular to the plane formed by directions df and db,
aligned so that it points towards the outside of surface F if the objects are 
separated. The Cx direction is df, and Cy=Cz X Cx to make a right-handed frame.
This method calculates a valid value even if this constraint is currently 
disabled. The given \a state must already be realized to Stage::Position. **/
Transform findContactFrameInG(const State& state) const;

/** Calculate the closest points on each of the two lines, measured and 
expressed in Ground. When the constraint is perfectly satisfied, these two
points will be in the same location. When the lines are parallel these points
are not unique; the returned point on each line will be midway between the
line origin point and the projection of the other line's origin point on this
line, and \a linesAreParallel will be returned \c true.  This calculates a valid
value even if this constraint is currently disabled. The given \a state must be
realized to Stage::Position. **/
void findClosestPointsInG(const State& state, Vec3& Qf, Vec3& Qb,
                          bool& linesAreParallel) const;

/** Calculate the separation distance or penetration depth of the two edges.
It is positive when the closest point Qb on line Lb lies outside the surface
containing line Lf, as indicated by the adjacent-face normals provided for Lf.
It is negative when point Qb lies below both of line Lf's adjacent faces in
which case it measures the penetration depth. This calculates a valid value 
even if this constraint is currently disabled. The given \a state must be 
realized to Stage::Position. **/
Real findSeparation(const State& state) const;
/*@}............................ Computations ................................*/

/** @cond **/ // hide from Doxygen
SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
    (LineOnLineContact, LineOnLineContactImpl, Constraint);
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_LINE_ON_LINE_CONTACT_H_



