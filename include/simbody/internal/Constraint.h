#ifndef SimTK_SIMBODY_CONSTRAINT_H_
#define SimTK_SIMBODY_CONSTRAINT_H_

/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

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

#include <cassert>

namespace SimTK {

class SimbodyMatterSubsystem;
class MobilizedBody;

/**
 * This is the base class for all Constraint classes, which is just a handle for the underlying
 * hidden implementation. Each built-in Constraint type is a local subclass within
 * Constraint, and is also derived from Constraint.
 */
class SimTK_SIMBODY_EXPORT Constraint {
public:
    Constraint() : rep(0) { }
    Constraint(Constraint&); // shallow copy
    Constraint& operator=(Constraint&); // shallow assignment
    ~Constraint();

    // These will fail unless this Constraint is owned by a MatterSubsystem.
    const SimbodyMatterSubsystem& getMatterSubsystem()      const;
    SimbodyMatterSubsystem&       updMatterSubsystem();

    ConstraintId           getConstraintId()      const;

    bool isInSubsystem() const;
    bool isInSameSubsystem(const MobilizedBody&) const;

    // Topology stage (i.e., construction).
    // nothing in base class currently

    // Model stage
    int getNumConstraintEquations(const State&) const;

    // Instance stage
    // Calling these reduces stage to Stage::Model.
    // nothing in base class currently

    // Position stage
    Vector getPositionErr(const State&) const;

    // Velocity stage
    Vector getVelocityErr(const State&) const;

    // Acceleration stage
    Vector getAccelerationError(const State&) const;
    Vector getMultipliers(const State&) const;

    // These are the built-in Constraint types. Types on the same line are
    // synonymous.
    class Rod;  typedef Rod  ConstantDistance;
    class Ball; typedef Ball CoincidentPoints;
    class Weld; typedef Weld CoincidentFrames;
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

/// One constraint equation. This constraint enforces a constant distance between
/// a point on one body and a point on another body. This is like connecting them
/// by a rigid, massless rod with ball joints at either end.
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
    MobilizedBodyId getBody1Id() const;
    MobilizedBodyId getBody2Id() const;
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
    Real getRodTension(const State&) const;

    class RodRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Rod, Constraint);
private:
    class RodRep& updRep();
    const RodRep& getRep() const;
};

/// Three constraint equations. This constraint enforces coincident location between
/// a point on one body and a point on another body.
class SimTK_SIMBODY_EXPORT Constraint::Ball : public Constraint {
public:
    // no default constructor
    Ball(MobilizedBody& body1, MobilizedBody& body2);
    Ball(MobilizedBody& body1, const Vec3& defaultPoint1,
         MobilizedBody& body2, const Vec3& defaultPoint2);

    // Defaults for Instance variables.
    Ball& setDefaultPointOnBody1(const Vec3&);
    Ball& setDefaultPointOnBody2(const Vec3&);

    // Stage::Topology
    MobilizedBodyId getBody1Id() const;
    MobilizedBodyId getBody2Id() const;
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

/// Six constraint equations. This constraint enforces coincidence between
/// a frame on one body and a frame on another body.
class SimTK_SIMBODY_EXPORT Constraint::Weld : public Constraint {
public:
    // no default constructor
    Weld(MobilizedBody& body1, MobilizedBody& body2);
    Weld(MobilizedBody& body1, const Transform& frame1,
         MobilizedBody& body2, const Transform& frame2);

    // Defaults for Instance variables.
    Weld& setDefaultFrameOnBody1(const Transform&);
    Weld& setDefaultFrameOnBody2(const Transform&);

    // Stage::Topology
    MobilizedBodyId getBody1Id() const;
    MobilizedBodyId getBody2Id() const;
    const Transform& getDefaultFrameOnBody1() const;
    const Transform& getDefaultFrameOnBody2() const;

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

/// Six constraint equations. This constraint enforces coincidence between
/// a frame on one body and a frame on another body.
class SimTK_SIMBODY_EXPORT Constraint::Custom : public Constraint {
public:
    // no default constructor (?)
    explicit Custom(int nEquations);

    virtual void calcPositionErrs(const State&, Vector& errs) const=0;
    virtual void calcVelocityErrs(const State&, Vector& errs) const=0;
    virtual void calcAccelerationErrs(const State&, Vector& errs) const=0;
    virtual void applyConstraintForces(const State&,
                                       const Vector&        multipliers,
                                       Vector_<SpatialVec>& bodyForces,
                                       Vector_<Vec3>&       particleForces,
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



