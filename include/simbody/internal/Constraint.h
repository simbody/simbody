#ifndef SimTK_SIMBODY_CONSTRAINT_H_
#define SimTK_SIMBODY_CONSTRAINT_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * This defines the Constraint class, which is used to specify limitations
 * on the mobility of the MobilizedBodies in a MatterSubsystem.
 *
 * Constraint is a PIMPL-style abstract base class, with concrete classes defined
 * for each kind of constraint. There are a set of built-in constraints
 * and a generic "Custom" constraint (an abstract base class) from
 * which advanced users may derive their own constraints.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <cassert>

namespace SimTK {

class MobilizedBody;


    ///////////////////////////
    // CONSTRAINT BASE CLASS //
    ///////////////////////////

/**
 * This is the base class for all Constraint classes, which is just a handle for the underlying
 * hidden implementation. Each built-in Constraint type is a local subclass within
 * Constraint, and is also derived from Constraint.
 *
 * TODO: should derive from PIMPLHandle
 */
class SimTK_SIMBODY_EXPORT Constraint {
public:
    Constraint() : rep(0) { }
    Constraint(Constraint&); // shallow copy
    Constraint& operator=(Constraint&); // shallow assignment
    ~Constraint();

    // These will fail unless this Constraint is owned by a MatterSubsystem.
    ConstraintIndex               getConstraintIndex()      const;
    const SimbodyMatterSubsystem& getMatterSubsystem()      const;
    SimbodyMatterSubsystem&       updMatterSubsystem();

    bool isInSubsystem() const;
    bool isInSameSubsystem(const MobilizedBody&) const;

        // TOPOLOGY STAGE (i.e., post-construction) //

    /// Return the number of bodies *directly* restricted by this
    /// constraint. Included are any bodies to which the Constraint may
    /// apply a body force (i.e., torque or point force). The Ancestor body is not
    /// included unless it was specified as a ConstrainedBody. This is the length
    /// of the bodyForces array for this Constraint.
    int getNumConstrainedBodies() const;
    /// Return a reference to the actual MobilizedBodies included in the count
    /// above. 0 <= index < getNumConstrainedBodies().
    const MobilizedBody& getMobilizedBodyFromConstrainedBody(ConstrainedBodyIndex) const;
    /// Return a reference to the actual MobilizedBody which is serving as
    /// the Ancestor body for the constrained bodies in this Constraint. This
    /// will fail if there are no constrained bodies
    /// (i.e., if getNumConstrainedBodies()==0).
    const MobilizedBody& getAncestorMobilizedBody() const;

    /// Return the number of mobilizers *directly* restricted by this
    /// constraint. Included are any mobilizers to which the Constraint may
    /// apply any mobility force. Like bodies, mobilizers are referenced 
    /// using the MobilizedBody containing them.
    int getNumConstrainedMobilizers() const;
    /// Return a reference to the actual MobilizedBodies included in the count
    /// of constrained mobilizers above. 0 <= index < getNumConstrainedMobilizers().
    const MobilizedBody& getMobilizedBodyFromConstrainedMobilizer(ConstrainedMobilizerIndex) const;

    const SimbodyMatterSubsystem::Subtree& getSubtree() const;

        // MODEL STAGE //

    /// Return the number of constrainable generalized coordinates q
    /// associated with a particular constrained mobilizer. This is just
    /// the number of generalized coordinates for that mobilizer.
    int getNumConstrainedQ(const State&, ConstrainedMobilizerIndex) const;
    /// Return the number of constrainable mobilities u associated with a particular
    /// constrained mobilizer. This is just the number of generalized speeds for that
    /// mobilizer.
    int getNumConstrainedU(const State&, ConstrainedMobilizerIndex) const;

    /// Return the index into the constrained mobilities u array corresponding to
    /// a particular mobility of the indicated ConstrainedMobilizer. Don't confuse
    /// this with the set of *participating* mobilities which also includes all 
    /// mobilities on each branch between the ancestor and a constrained body.
    /// The *constrained* mobilities are just those belonging to the mobilizers
    /// which are directly constrained.
    ConstrainedUIndex getConstrainedUIndex
        (const State&, ConstrainedMobilizerIndex, MobilizerUIndex which) const;
    /// Return the index into the constrained coordinates q array corresponding to
    /// a particular coordinate of the indicated ConstrainedMobilizer. Don't confuse
    /// this with the set of *participating* coordinates which also includes all 
    /// coordinates on each branch between the ancestor and a constrained body.
    /// The *constrained* coordinates are just those belonging to the mobilizers
    /// which are directly constrained.
    ConstrainedQIndex getConstrainedQIndex
        (const State&, ConstrainedMobilizerIndex, MobilizerQIndex which) const;

    /// Return the sum of the number of coordinates q associated with each of
    /// the constrained mobilizers.
    int getNumConstrainedQ(const State&) const;

    /// Return the sum of the number of mobilities u associated with each of
    /// the constrained mobilizers. These are the only mobilities to which 
    /// the constraint may directly apply a force, so this is the dimension
    /// of the mobilityForces array.
    int getNumConstrainedU(const State&) const;
    
    /// Find out how many holonomic (position), nonholonomic (velocity),
    /// and acceleration-only constraint equations are generated by this Constraint.
    void getNumConstraintEquations(const State&, int& mp, int& mv, int& ma) const;

        // INSTANCE STAGE //
    // nothing in base class currently

        // POSITION STAGE //
    Vector getPositionError(const State&) const;	// mp of these
	Vector calcPositionErrorFromQ(const State&, const Vector& q) const;

	// Matrix P = partial(perr_dot)/partial(u). (just the holonomic constraints)
	Matrix calcPositionConstraintMatrixP(const State&) const; // mp X nu
	Matrix calcPositionConstraintMatrixPt(const State&) const; // nu X mp

	// Matrix PQInv = partial(perr)/partial(q) = P*Q^-1
	Matrix calcPositionConstraintMatrixPQInverse(const State&) const; // mp X nq

    // This operator calculates this constraint's body and mobility forces
    // given the complete set of multipliers lambda. We expect that lambda
    // has been packed to include multipliers associated with the
    // second derivatives of the position
    // (holonomic) constraints, the first derivatives of the velocity
    // (nonholonomic) constraints, and the acceleration only constraints, in
    // that order.
    // The state must be realized already to Stage::Position. Returned body
    // forces correspond only to the *constrained bodies* and the mobility
    // forces correspond only to the *constrained mobilities*; they must 
    // be unpacked by the caller into the actual system mobilized bodies and 
    // actual system mobilities.
    // Note that the body forces are in the ancestor body frame A, not necessarily
    // the Ground frame G.
    void calcConstraintForcesFromMultipliers(const State&,
        const Vector&        lambda,                // mp+mv+ma of these
        Vector_<SpatialVec>& bodyForcesInA,         // numConstrainedBodies
        Vector&              mobilityForces) const; // numConstrainedU

        // VELOCITY STAGE //
    Vector getVelocityError(const State&) const;	// mp+mv of these
	Vector calcVelocityErrorFromU(const State&,     // mp+mv of these
                                  const Vector& u) const;   // numParticipatingU u's

	// Matrix V = partial(verr)/partial(u) for just the non-holonomic constraints.
	Matrix calcVelocityConstraintMatrixV(const State&) const;  // mv X nu
	Matrix calcVelocityConstraintMatrixVt(const State&) const; // nu X mv

        // DYNAMICS STAGE //
    // nothing in base class currently

        // ACCELERATION STAGE //
    Vector getAccelerationError(const State&) const;	// mp+mv+ma of these
	Vector calcAccelerationErrorFromUDot(const State&,  // mp+mv+ma of these
                                         const Vector& udot) const // numParticipatingU udot's
    {
		assert(!"calcAccelerationErrorFromUDot: Not implemented yet");
		return Vector();
	}

    Vector getMultipliers(const State&) const;			// mp+mv+ma of these   

	// Matrix A = partial(aerr)/partial(udot) for just the acceleration-only constraints.
	Matrix calcAccelerationConstraintMatrixA(const State&) const;  // ma X nu
	Matrix calcAccelerationConstraintMatrixAt(const State&) const; // nu X ma

    // These are the built-in Constraint types. Types on the same line are
    // synonymous.
    class Rod;  typedef Rod  ConstantDistance;
    class Ball; typedef Ball CoincidentPoints;
    class Weld; typedef Weld CoincidentFrames;
    class PointInPlane;  // translations perpendicular to plane normal only
    class PointOnLine;   // translations along a line only
    class ConstantAngle; // prevent rotation about common normal of two vectors
    class ConstantOrientation; // allows any translation but no rotation
    class NoSlip1D; // same velocity at a point along a direction
    class ConstantSpeed; // prescribe generalized speed value
    class Custom;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;


    // Internal use only

    // The current handle is the owner of the rep. After this call
    // the supplied handle is the owner and this one is just a reference.
    void disown(Constraint&);
    class ConstraintRep; // local subclass
    explicit Constraint(class ConstraintRep* r) : rep(r) { }
    bool                 hasRep() const {return rep!=0;}
    const ConstraintRep& getRep() const {assert(rep); return *rep;}
    ConstraintRep&       updRep() const {assert(rep); return *rep;}
	void setRep(ConstraintRep& r) {assert(!rep); rep = &r;}
protected:
    class ConstraintRep* rep;
};

    ////////////////////////////////////////
    // ROD (CONSTANT DISTANCE) CONSTRAINT //
    ////////////////////////////////////////

/**
 *  One constraint equation. This constraint enforces a constant distance between
 *  a point on one body and a point on another body. This is like connecting them
 *  by a rigid, massless rod with ball joints at either end. The constraint is
 *  enforced by a force acting along the rod with opposite signs at either end.
 *  When positive, this represents tension in the rod pulling the points together;
 *  when negative it represents compression keeping the points separated.
 * 
 *  Caution: you can't use this to enforce a distance of zero between two points.
 *  That takes three constraints because there is no restriction on the force direction.
 *  For a distance of zero (i.e., you want the points to be coincident) use a Ball
 *  constraint, a.k.a. CoincidentPoints constraint.
 */
class SimTK_SIMBODY_EXPORT Constraint::Rod : public Constraint {
public:
    // no default constructor
    Rod(MobilizedBody& body1, MobilizedBody& body2,
        Real defaultLength=1);
    Rod(MobilizedBody& body1, const Vec3& defaultPoint1,
        MobilizedBody& body2, const Vec3& defaultPoint2,
        Real defaultLength=1);

    // Defaults for Instance variables.
    Rod& setDefaultPointOnBody1(const Vec3&);
    Rod& setDefaultPointOnBody2(const Vec3&);
    Rod& setDefaultRodLength(Real);

    // Stage::Topology
    MobilizedBodyIndex getBody1MobilizedBodyIndex() const;
    MobilizedBodyIndex getBody2MobilizedBodyIndex() const;
    const Vec3& getDefaultPointOnBody1() const;
    const Vec3& getDefaultPointOnBody2() const;
    Real getDefaultRodLength() const;

    // Stage::Instance
    const Vec3& getPointOnBody1(const State&) const;
    const Vec3& getPointOnBody2(const State&) const;
    Real        getRodLength   (const State&) const;

    // Stage::Position, Velocity, Acceleration
    Real getPositionError(const State&) const;
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getRodTension(const State&) const; // negative means compression

    class RodRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Rod, Constraint);
private:
    class RodRep& updRep();
    const RodRep& getRep() const;
};

    ///////////////////////////////
    // POINT IN PLANE CONSTRAINT //
    ///////////////////////////////

/**
 *  One constraint equation. This constraint enforces that a point fixed to
 *  one body (the "follower body") must travel in a plane fixed on another body (the
 *  "plane body"). The constraint is enforced by an internal (non-working)
 *  scalar force acting at the spatial location of the follower point, directed along
 *  the plane normal, and equal and opposite on the two bodies.
 * 
 *  The assembly condition is the same as the run-time constraint: the point
 *  has to be moved into the plane.
 */
class SimTK_SIMBODY_EXPORT Constraint::PointInPlane : public Constraint {
public:
    // no default constructor
    PointInPlane(MobilizedBody& planeBody_B, const UnitVec3& defaultPlaneNormal_B, Real defaultHeight,
                 MobilizedBody& followerBody_F, const Vec3& defaultFollowerPoint_F);

    // These affect only generated decorative geometry for visualization;
    // the plane is really infinite in extent with zero depth and the
    // point is really of zero radius.
    PointInPlane& setPlaneDisplayHalfWidth(Real);
    PointInPlane& setPointDisplayRadius(Real);
    Real getPlaneDisplayHalfWidth() const;
    Real getPointDisplayRadius() const;

    // Defaults for Instance variables.
    PointInPlane& setDefaultPlaneNormal(const UnitVec3&);
    PointInPlane& setDefaultPlaneHeight(Real);
    PointInPlane& setDefaultFollowerPoint(const Vec3&);

    // Stage::Topology
    MobilizedBodyIndex getPlaneMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const UnitVec3& getDefaultPlaneNormal() const;
    Real            getDefaultPlaneHeight() const;
    const Vec3&     getDefaultFollowerPoint() const;

    // Stage::Instance
    const UnitVec3& getPlaneNormal(const State&) const;
    Real            getPlaneHeight(const State&) const;
    const Vec3&     getFollowerPoint(const State&) const;

    // Stage::Position, Velocity
    Real getPositionError(const State&) const;
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getForceOnFollowerPoint(const State&) const; // in normal direction

    class PointInPlaneRep; // local subclass

    SimTK_PIMPL_DOWNCAST(PointInPlane, Constraint);
private:
    class PointInPlaneRep& updRep();
    const PointInPlaneRep& getRep() const;
};

    //////////////////////////////
    // POINT ON LINE CONSTRAINT //
    //////////////////////////////

/**
 *  Two constraint equations. This constraint enforces that a point fixed to
 *  one body (the "follower body") must travel along a line fixed on another body (the
 *  "line body"). The constraint is enforced by an internal (non-working)
 *  scalar force acting at the spatial location of the follower point, directed in the
 *  plane for which the line is a normal, and equal and opposite on the two bodies.
 * 
 *  The assembly condition is the same as the run-time constraint: the point
 *  has to be moved onto the line.
 */
class SimTK_SIMBODY_EXPORT Constraint::PointOnLine : public Constraint {
public:
    // no default constructor
    PointOnLine(MobilizedBody& lineBody_B, const UnitVec3& defaultLineDirection_B, const Vec3& defaultPointOnLine_B,
                MobilizedBody& followerBody_F, const Vec3& defaultFollowerPoint_F);

    // These affect only generated decorative geometry for visualization;
    // the line is really infinite in extent and the
    // point is really of zero radius.
    PointOnLine& setLineDisplayHalfLength(Real);
    PointOnLine& setPointDisplayRadius(Real);
    Real getLineDisplayHalfLength() const;
    Real getPointDisplayRadius() const;

    // Defaults for Instance variables.
    PointOnLine& setDefaultLineDirection(const UnitVec3&);
    PointOnLine& setDefaultPointOnLine(const Vec3&);
    PointOnLine& setDefaultFollowerPoint(const Vec3&);

    // Stage::Topology
    MobilizedBodyIndex getLineMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const UnitVec3& getDefaultLineDirection() const;
    const Vec3&     getDefaultPointOnLine() const;
    const Vec3&     getDefaultFollowerPoint() const;

    // Stage::Instance
    const UnitVec3& getLineDirection(const State&) const;
    const Vec3&     getPointOnLine(const State&) const;
    const Vec3&     getFollowerPoint(const State&) const;

    // Stage::Position, Velocity
    const Vec2& getPositionError(const State&) const;
    const Vec2& getVelocityError(const State&) const;

    // Stage::Acceleration
    const Vec2& getAccelerationErrors(const State&) const;
    const Vec2& getMultipliers(const State&) const;
    const Vec2& getForceOnFollowerPoint(const State&) const; // in normal direction

    class PointOnLineRep; // local subclass

    SimTK_PIMPL_DOWNCAST(PointOnLine, Constraint);
private:
    class PointOnLineRep& updRep();
    const PointOnLineRep& getRep() const;
};

    ///////////////////////////////
    // CONSTANT ANGLE CONSTRAINT //
    ///////////////////////////////

/**
 *  One constraint equation. This constraint enforces that a vector fixed to
 *  one body (the "base body") must maintain a fixed angle with respect to
 *  a vector fixed on the other body (the "follower body"). That is, we have a
 *  single constraint equation that prohibits rotation about the mutual perpendicular
 *  to the two vectors.
 * 
 *  This constraint is enforced by an internal scalar torque applied equal and
 *  opposite on each body, about the mutual perpendicular to the two vectors.
 * 
 *  The assembly condition is the same as the run-time constraint: the 
 *  bodies must be rotated until the vectors have the right angle between them.
 */
class SimTK_SIMBODY_EXPORT Constraint::ConstantAngle : public Constraint {
public:
    // no default constructor
    ConstantAngle(MobilizedBody& baseBody_B,     const UnitVec3& defaultAxis_B,
                  MobilizedBody& followerBody_F, const UnitVec3& defaultAxis_F, 
                  Real angle = Pi/2);

    // These affect only generated decorative geometry for visualization.
    ConstantAngle& setAxisDisplayLength(Real);
    ConstantAngle& setAxisDisplayWidth(Real);
    Real getAxisDisplayLength() const;
    Real getAxisDisplayWidth() const;

    // Defaults for Instance variables.
    ConstantAngle& setDefaultBaseAxis(const UnitVec3&);
    ConstantAngle& setDefaultFollowerAxis(const UnitVec3&);
    ConstantAngle& setDefaultAngle(Real);

    // Stage::Topology
    MobilizedBodyIndex getBaseMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const UnitVec3& getDefaultBaseAxis() const;
    const UnitVec3& getDefaultFollowerAxis() const;
    Real getDefaultAngle() const;

    // Stage::Instance
    const UnitVec3& getBaseAxis(const State&) const;
    const UnitVec3& getFollowerAxis(const State&) const;
    Real getAngle(const State&) const;

    // Stage::Position, Velocity
    Real getPositionError(const State&) const;
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getTorqueOnFollowerBody(const State&) const; // about f X b

    class ConstantAngleRep; // local subclass

    SimTK_PIMPL_DOWNCAST(ConstantAngle, Constraint);
private:
    class ConstantAngleRep& updRep();
    const ConstantAngleRep& getRep() const;
};

    /////////////////////////////////////////
    // BALL (COINCIDENT POINTS) CONSTRAINT //
    /////////////////////////////////////////

/**
 *  Three constraint equations. This constraint enforces coincident location between
 *  a point on one body and a point on another body.
 * 
 *  The constraint is enforced by an internal (non-working) force applied at the
 *  spatial location of the point on body 2, on material points of each body that
 *  are coincident with that spatial location. Note that this is somewhat asymmetric
 *  when the ball is not properly assembled -- it acts as though the contact occurs
 *  at the point on body 2, *not* at the point on body 1.
 * 
 *  The assembly condition is the same as the runtime constraint -- the two points
 *  can be brought together by driving the perr to zero.
 */
class SimTK_SIMBODY_EXPORT Constraint::Ball : public Constraint {
public:
    // no default constructor
    Ball(MobilizedBody& body1, MobilizedBody& body2);
    Ball(MobilizedBody& body1, const Vec3& defaultPoint1,
         MobilizedBody& body2, const Vec3& defaultPoint2);

    // Defaults for Instance variables.
    Ball& setDefaultPointOnBody1(const Vec3&);
    Ball& setDefaultPointOnBody2(const Vec3&);

    // This is used only for visualization.
    Ball& setDefaultRadius(Real r);
    Real getDefaultRadius() const;

    // Stage::Topology
    MobilizedBodyIndex getBody1MobilizedBodyIndex() const;
    MobilizedBodyIndex getBody2MobilizedBodyIndex() const;
    const Vec3& getDefaultPointOnBody1() const;
    const Vec3& getDefaultPointOnBody2() const;

    // Stage::Instance
    const Vec3& getPointOnBody1(const State&) const;
    const Vec3& getPointOnBody2(const State&) const;

    // Stage::Position, Velocity, Acceleration
    const Vec3& getPositionErrors(const State&) const;
    const Vec3& getVelocityErrors(const State&) const;

    // Stage::Acceleration
    const Vec3& getAccelerationErrors(const State&) const;
    const Vec3& getMultipliers(const State&) const;

    // Forces are reported expressed in the body frame of the indicated body.
    const Vec3& getBallReactionForceOnBody1(const State&) const;
    const Vec3& getBallReactionForceOnBody2(const State&) const;

    class BallRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Ball, Constraint);
private:
    class BallRep& updRep();
    const BallRep& getRep() const;
};

    /////////////////////////////////////
    // CONSTANT ORIENTATION CONSTRAINT //
    /////////////////////////////////////

/**
 *  Three constraint equations. This constraint enforces that a reference frame
 *  fixed to one body (the "follower body") must have the same orientation as another
 *  reference frame fixed on another body (the "base body"). That is, we have three
 *  constraint equations that collectively prohibit any relative rotation between
 *  the base and follower. The run time equations we use are just three "constant angle"
 *  constraints enforcing perpendicularity between follower's x,y,z axes with the base
 *  y,z,x axes respectively.
 * 
 *  This constraint is enforced by an internal (non-working) torque vector applied equal and
 *  opposite on each body.
 * 
 *  TODO: The assembly condition is not the same as the run-time constraint, because the
 *  perpendicularity conditions can be satisfied with antiparallel axes. For assembly
 *  we must have additional (redundant) constraints requiring parallel axes.
 */
class SimTK_SIMBODY_EXPORT Constraint::ConstantOrientation : public Constraint {
public:
    // no default constructor
    ConstantOrientation(MobilizedBody& baseBody_B,     const Rotation& defaultRB,
                        MobilizedBody& followerBody_F, const Rotation& defaultRF); 

    //TODO: default visualization geometry?

    // Defaults for Instance variables.
    ConstantOrientation& setDefaultBaseRotation(const Rotation&);
    ConstantOrientation& setDefaultFollowerRotation(const Rotation&);

    // Stage::Topology
    MobilizedBodyIndex getBaseMobilizedBodyIndex() const;
    MobilizedBodyIndex getFollowerMobilizedBodyIndex() const;

    const Rotation& getDefaultBaseRotation() const;
    const Rotation& getDefaultFollowerRotation() const;

    // Stage::Instance
    const Rotation& getBaseRotation(const State&) const;
    const Rotation& getFollowerRotation(const State&) const;

    // Stage::Position, Velocity
    Vec3 getPositionError(const State&) const;
    Vec3 getVelocityError(const State&) const;

    // Stage::Acceleration
    Vec3 getAccelerationError(const State&) const;
    Vec3 getMultiplier(const State&) const;
    Vec3 getTorqueOnFollowerBody(const State&) const;

    class ConstantOrientationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(ConstantOrientation, Constraint);
private:
    class ConstantOrientationRep& updRep();
    const ConstantOrientationRep& getRep() const;
};

    /////////////////////////////////////////
    // WELD (COINCIDENT FRAMES) CONSTRAINT //
    /////////////////////////////////////////

/**
 *  Six constraint equations. This constraint enforces coincidence between
 *  a frame on one body and a frame on another body. This is a combination
 *  of a ConstantOrientation constraint and a Ball constraint. The first three
 *  equations correspond to the perpendicularity constraints associated with
 *  the orientation constraint, the last three equations are the 
 *  coincident point conditions.
 * 
 *  The constraint is enforced by an internal (non-working) force applied at the
 *  spatial location of the frame origin on body 2, on material points of each body that
 *  are coincident with that spatial location. Note that this is somewhat asymmetric
 *  when the Weld is not properly assembled -- it acts as though the contact occurs
 *  at the origin of the frame on body 2, *not* at the origin of the frame on body 1.
 *  The orientation constrains on the other hand are symmetric, they are three
 *  "constant angle" constraints enforcing perpendicularity between body2's
 *  x,y,z axes with body1's y,z,x axes respectively, via an internal (non-working)
 *  torque vector applied equal and opposite on each body.
 * 
 *  TODO: Although the frame origins can be brought together by the Ball constraint, the
 *  perpendicularity conditions can be satisfied with antiparallel axes in addition
 *  to the parallel ones we want. Therefore the assembly conditions must include
 *  additional (redundant) constraints requiring parallel axes.
 */
class SimTK_SIMBODY_EXPORT Constraint::Weld : public Constraint {
public:
        // no default constructor

    /// Make the body frame of one body coincident with the body frame
    /// of the other body.
    Weld(MobilizedBody& body1, MobilizedBody& body2);

    /// Make a particular frame attached to one body coincident with
    /// a particular frame attached to the other body. The frames are
    /// specified by giving the transform X_BF which expresses the
    /// position and orientation of frame F relative to the body frame B.
    Weld(MobilizedBody& body1, const Transform& frame1,
         MobilizedBody& body2, const Transform& frame2);

        // Control over generated decorative geometry.

    /// This is used only for visualization. Set r <= 0 to disable
    /// default frame drawing. Default axis length is r=1. This is a
    /// topology-stage variable, not changeable later.
    Weld& setAxisDisplayLength(Real r);

    /// Report the length being used for display of the frames being
    /// connected by this Weld. If this returns 0 then no geometry is
    /// being generated for the frames.
    Real getAxisDisplayLength() const;

        // Defaults for Instance variables.

    /// Explicitly set the default value for the frame on body1 which
    /// is to be made coincident with a frame on body2. Note that this is
    /// topology-stage value so requires non-const access to the Constraint.
    Weld& setDefaultFrameOnBody1(const Transform&);

    /// Retrieve the default transform for the frame on body 1.
    const Transform& getDefaultFrameOnBody1() const;

    /// Explicitly set the default value for the frame on body2 which
    /// is to be made coincident with a frame on body1. Note that this is
    /// topology-stage value so requires non-const access to the Constraint.
    Weld& setDefaultFrameOnBody2(const Transform&);

    /// Retrieve the default transform for the frame on body 2.
    const Transform& getDefaultFrameOnBody2() const;


        // Stage::Topology

    /// Report the MobilizedBodyIndex of body 1 for this Weld constraint.
    MobilizedBodyIndex getBody1MobilizedBodyIndex() const;

    /// Report the MobilizedBodyIndex of body 2 for this Weld constraint.
    MobilizedBodyIndex getBody2MobilizedBodyIndex() const;


        // Stage::Instance
    const Transform& getFrameOnBody1(const State&) const;
    const Transform& getFrameOnBody2(const State&) const;

        // Stage::Position, Velocity, Acceleration
    const Vec6& getPositionErrors(const State&) const;
    const Vec6& getVelocityErrors(const State&) const;

        // Stage::Acceleration
    const Vec6& getAccelerationErrors(const State&) const;
    const Vec6& getMultipliers(const State&) const;

        // Forces are reported expressed in the body frame of the indicated body.
    const SpatialVec& getWeldReactionOnBody1(const State&) const;
    const SpatialVec& getWeldReactionOnBody2(const State&) const;

    class WeldRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Weld, Constraint);
private:
    class WeldRep& updRep();
    const WeldRep& getRep() const;
};

    ///////////////////////////
    // NO SLIP 1D CONSTRAINT //
    ///////////////////////////

/**
 * One non-holonomic constraint equation. There is a contact point P and a no-slip 
 * direction n fixed in a case body C. There are two moving bodies B0 and B1. The 
 * material point of B0 and the material point of B1 which are each coincident 
 * with the contact point P must have identical velocities in C, along the direction n.
 * This can be used to implement simple rolling contact between disks, such as occurs
 * in gear trains.
 * 
 * The assembly condition is the same as the run-time constraint: the velocities must
 * be made to match.
 */
class SimTK_SIMBODY_EXPORT Constraint::NoSlip1D : public Constraint {
public:
    // no default constructor
    NoSlip1D(MobilizedBody& caseBodyC, const Vec3& P_C, const UnitVec3& n_C,
             MobilizedBody& movingBody0, MobilizedBody& movingBody1);

    // These affect only generated decorative geometry for visualization;
    // the plane is really infinite in extent with zero depth and the
    // point is really of zero radius.
    NoSlip1D& setDirectionDisplayLength(Real);
    NoSlip1D& setPointDisplayRadius(Real);
    Real getDirectionDisplayLength() const;
    Real getPointDisplayRadius() const;

    // Defaults for Instance variables.
    NoSlip1D& setDefaultDirection(const UnitVec3&);
    NoSlip1D& setDefaultContactPoint(const Vec3&);

    // Stage::Topology
    MobilizedBodyIndex getCaseMobilizedBodyIndex() const;
    MobilizedBodyIndex getMovingBodyMobilizedBodyIndex(int which) const;

    const UnitVec3& getDefaultDirection() const;
    const Vec3&     getDefaultContactPoint() const;

    // Stage::Instance
    const UnitVec3& getDirection(const State&) const;
    const Vec3&     getContactPoint(const State&) const;

    // Stage::Position, Velocity
        // no position error
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getForceAtContactPoint(const State&) const; // in normal direction, no body 2

    class NoSlip1DRep; // local subclass

    SimTK_PIMPL_DOWNCAST(NoSlip1D, Constraint);
private:
    class NoSlip1DRep& updRep();
    const NoSlip1DRep& getRep() const;
};

    ////////////////////
    // CONSTANT SPEED //
    ////////////////////

/**
 * One non-holonomic constraint equation. Some mobility u is required to be at a
 * particular value s.
 * 
 * The assembly condition is the same as the run-time constraint: u must be set to s.
 */
class SimTK_SIMBODY_EXPORT Constraint::ConstantSpeed : public Constraint {
public:
    // no default constructor
   ConstantSpeed(MobilizedBody& mobilizer, MobilizerUIndex, Real speed);
   ConstantSpeed(MobilizedBody& mobilizer, Real speed); // only if 1 dof mobilizer

    // Stage::Topology
    MobilizedBodyIndex getMobilizedBodyIndex() const;
    MobilizerUIndex    getWhichU() const;
    Real               getDefaultSpeed() const;

    // Stage::Position, Velocity
        // no position error
    Real getVelocityError(const State&) const;

    // Stage::Acceleration
    Real getAccelerationError(const State&) const;
    Real getMultiplier(const State&) const;
    Real getGeneralizedForce(const State&) const;

    class ConstantSpeedRep; // local subclass

    SimTK_PIMPL_DOWNCAST(ConstantSpeed, Constraint);
private:
    class ConstantSpeedRep& updRep();
    const ConstantSpeedRep& getRep() const;
};

// TODO: this is just a sketch of a Custom Constraint base class.
class SimTK_SIMBODY_EXPORT Constraint::Custom : public Constraint {
public:
    // no default constructor (?)

    // These are the topological defaults for the number of holonomic, nonholonomic,
    // and acceleration only constraint equations to be generated. The actual number
    // can be changed prior to realizeModel().
    Custom(int mp, int mv, int ma);

        // Topological information//

    void setDefaultNumConstraintEquations(int mp, int mv, int ma);
    void getDefaultNumConstraintEquations(int& mp, int& mv, int& ma) const;

    // Start numbering from 0 for each Constraint. The supplied MobilizedBody must be
    // in the Matter Subsystem of which this Constraint is a part.
    ConstrainedBodyIndex addConstrainedBody(const MobilizedBody&);

    // Alternatively, declare this as a global constraint. (Constant energy or temperature
    // might be an example.)
    void setAllBodiesAreConstrained(bool);
    
    // getNumConstrainedBodies() and getConstrainedBody() are in the base class.

        // Model stage information //

    // Turn off this constraint altogether, but don't forget mp,mv, and ma.
    void setEnableConstraint(State&, bool) const;
    // isConstraintEnabled() is in base class.

    // Set Model-stage state variables to reflect the actual number of equations to 
    // be generated.
    void changeNumConstraintEquations(State&, int mp, int mv, int ma) const;

    // getNumConstraintEquations(), getNumMobilities(), and getParticipatingMobilities()
    // are in the base class.

        // Constraint implementation //

    // These must be defined if there are any position (holonomic) constraints defined.
    virtual void calcPositionErrors      (const State&, int mp,  Real* perr) const=0;
    virtual void calcPositionDotErrors   (const State&, int mp,  Real* pverr) const=0;
    virtual void calcPositionDotDotErrors(const State&, int mp,  Real* paerr) const=0;
    virtual void applyPositionConstraintForces(const State&,
                                               int mp, const Real*  multipliers,
                                               Vector_<SpatialVec>& bodyForces,
                                               Vector&              mobilityForces) const=0;

    // These must be defined if there are any velocity (nonholonomic) constraints defined.
    virtual void calcVelocityErrors   (const State&, int mv,  Real* verr) const=0;
    virtual void calcVelocityDotErrors(const State&, int mv,  Real* vaerr) const=0;
    virtual void applyVelocityConstraintForces(const State&,
                                               int mv, const Real*  multipliers,
                                               Vector_<SpatialVec>& bodyForces,
                                               Vector&              mobilityForces) const=0;

    // These must be defined if there are any acceleration-only constraints defined.
    virtual void calcAccelerationErrors(const State&, int ma,  Real* aerr) const=0;
    virtual void applyAccelerationConstraintForces(const State&,
                                                   int ma, const Real*  multipliers,
                                                   Vector_<SpatialVec>& bodyForces,
                                                   Vector&              mobilityForces) const=0;


    //TODO
    // - Derivatives with respect to u.
    // - How to handle nonholonomic constraints?
    // - How to get topological information: what bodies, what mobilities?

    class CustomRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Custom, Constraint);
private:
    class CustomRep& updRep();
    const CustomRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CONSTRAINT_H_



