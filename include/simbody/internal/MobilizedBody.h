#ifndef SimTK_SIMBODY_MOBILIZED_BODY_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_H_

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
 * This defines the MobilizedBody class, which associates some mass
 * (the "outboard" body) with a mobilizer and a reference frame (the
 * parent or "inboard" body), already present in a MatterSubsystem.
 *
 * MobilizedBody is an abstract base class, with concrete classes defined
 * for each kind of mobilizer. There are a set of built-in mobilizers
 * and a generic "Custom" mobilizer (an abstract base class) from
 * which advanced users may derive their own mobilizers.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

class SimbodyMatterSubsystem;
class Body;

/**
 * This is the base class for all MobilizedBody classes, which is just a handle for the underlying
 * hidden implementation. Each built-in MobilizedBody type is a local subclass within
 * MobilizedBody, and is also derived from MobilizedBody.
 */
class SimTK_SIMBODY_EXPORT MobilizedBody {
public:
    MobilizedBody() : rep(0) { }
    MobilizedBody(MobilizedBody&); // shallow copy
    MobilizedBody& operator=(MobilizedBody&); // shallow assignment
    ~MobilizedBody();

    // Implicit conversion to MobilizedBodyId when needed.
    operator MobilizedBodyId() const {return getMobilizedBodyId();}


    // These will fail unless this MobilizedBody is owned by a MatterSubsystem.
    const SimbodyMatterSubsystem& getMatterSubsystem()      const;
    SimbodyMatterSubsystem&       updMatterSubsystem();
    MobilizedBodyId        getMobilizedBodyId()      const;
    const MobilizedBody&   getInboardMobilizedBody() const;

    bool isInSubsystem() const;
    bool isInSameSubsystem(const MobilizedBody&) const;

    // Topology stage (i.e., construction).
    // Calling these means you are (re)constructing the system and will have to do
    // realizeTopology() and extract a new State before doing any analysis.
    const Body&    getBody() const;
    Body&          updBody();
    MobilizedBody& setBody(const Body&);

    MobilizedBody& setDefaultInboardFrame (const Transform& X_PMb);
    MobilizedBody& setDefaultOutboardFrame(const Transform& X_BM);
    const Transform& getDefaultInboardFrame() const;
    const Transform& getDefaultOutboardFrame() const;

    // Model stage
    int getNumQ(const State&) const;
    int getNumU(const State&) const;

    // Instance stage
    // Calling these reduces stage to Stage::Model.
    void setInboardFrame (State&, const Transform& X_PMb) const;
    void setOutboardFrame(State&, const Transform& X_BM ) const;

    const Transform& getInboardFrame (const State&) const;
    const Transform& getOutboardFrame(const State&) const;

    // Position stage
    const Transform& getBodyTransform(const State&) const;
    const Transform& getMobilizerTransform(const State& s) const;
    Vector getQ(const State&) const;
    Vector updQ(State&) const; // reduce to Stage::Time
    void setQ(State&, const Vector&) const;

    // Velocity stage
    const SpatialVec& getBodyVelocity(const State& s) const;
    const SpatialVec& getMobilizerVelocity(const State& s) const;
    Vector getU(const State&) const;
    Vector updU(State&) const; // reduce to Stage::Position
    void setU(State&, const Vector&) const;

    // Acceleration stage
    const SpatialVec& getBodyAppliedForces(const State&) const;
    SpatialVec&       updBodyAppliedForces(State&) const; // reduce to Stage::Dynamics
    void applyBodyForce(State& s, const SpatialVec& f) const {
        updBodyAppliedForces(s) += f;
    }

    const SpatialVec& getBodyAcceleration(const State& s) const;
    Vector getUDot(const State&) const;

    // These are the built-in MobilizedBody types. Types on the same line are
    // synonymous.
    class Pin;         typedef Pin    Torsion;
    class Slider;      typedef Slider Prismatic;
    class Universal;
    class Cylinder;
    class BendStretch;
    class Planar;
    class Gimbal;
    class Ball; typedef Ball Orientation, Spherical;
    class Translation; typedef Translation Cartesian;
    class Free;
    class LineOrientation;
    class FreeLine;
    class Weld;
    class Screw;
    class Ellipsoid;
    class Custom;
    class Ground;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;


    // Internal use only

    // The current handle is the owner of the rep. After this call
    // the supplied handle is the owner and this one is just a reference.
    void disown(MobilizedBody&);
    class MobilizedBodyRep; // local subclass
    explicit MobilizedBody(class MobilizedBodyRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const MobilizedBodyRep& getRep() const {assert(rep); return *rep;}
    MobilizedBodyRep&       updRep() const {assert(rep); return *rep;}
	void setRep(MobilizedBodyRep& r) {assert(!rep); rep = &r;}
protected:
    class MobilizedBodyRep* rep;
};

/// One mobility -- rotation about the common z axis of the inboard
/// and outboard mobilizer frames.
/// Synonym: Torsion
class SimTK_SIMBODY_EXPORT MobilizedBody::Pin : public MobilizedBody {
public:
    Pin();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Pin(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Pin(MobilizedBody& parent, const Transform& inbFrame,
        const Body&,           const Transform& outbFrame);

    Pin& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Pin& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    // This is just a nicer name for the generalized coordinate.
    Pin& setDefaultAngle(Real angleInRadians) {return setDefaultQ(angleInRadians);}
    Real getDefaultAngle() const {return getDefaultQ();}

    // Friendly, mobilizer-specific access to coordinates and speeds.
    // TODO

    void setAngle(State& s, Real angleInRadians) {setQ(s, angleInRadians);}
    Real getAngle(const State& s) {return getQ(s);}

    void applyPinTorque(State& s, Real t) const {
        updMobilizerForces(s) += t;
    }

    // Generic default state Topology methods.
    Pin&  setDefaultQ(Real q) {updDefaultQ()=q; return *this;}
    Real  getDefaultQ() const;
    Real& updDefaultQ();

    // Generic state access routines.
    Real getQ(const State&) const;
    Real getU(const State&) const;
    Real& updQ(State&) const;
    Real& updU(State&) const;
    void setQ(State& s, Real q) const {updQ(s)=q;}
    void setU(State& s, Real u) const {updU(s)=u;}

    // Acceleration stage state variables
    Real getMobilizerForces(const State&) const;
    Real& updMobilizerForces(State&) const;

    class PinRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Pin, MobilizedBody);
private:
    PinRep&       updRep();
    const PinRep& getRep() const;
};

/// One mobility -- translation along the common x axis of the
/// inboard and outboard mobilizer frames.
/// Synonym: Prismatic
class SimTK_SIMBODY_EXPORT MobilizedBody::Slider : public MobilizedBody {
public:
    Slider();
    class SliderRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Slider, MobilizedBody);
private:
    class SliderRep& updRep();
    const SliderRep& getRep() const;
};

/// One mobility -- coordinated rotation and translation along the
/// common z axis of the inboard and outboard mobilizer frames. A
/// "pitch" is specified relating the two. The generalized coordinate
/// q is the rotation angle in radians, the translation is always
/// pitch*q.
class SimTK_SIMBODY_EXPORT MobilizedBody::Screw : public MobilizedBody {
public:
    Screw(Real pitch);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Screw(MobilizedBody& parent, const Body&, Real pitch);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Screw(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame,
         Real pitch);

    Screw& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Screw& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    Screw& setDefaultPitch(Real pitch);

    class ScrewRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Screw, MobilizedBody);
private:
    class ScrewRep& updRep();
    const ScrewRep& getRep() const;
};

/// Two mobilities -- rotation about the x axis, followed by a rotation
/// about the new y axis. This mobilizer is badly behaved when the
/// second rotation is near 90 degrees.
class SimTK_SIMBODY_EXPORT MobilizedBody::Universal : public MobilizedBody {
public:
    Universal();
    class UniversalRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Universal, MobilizedBody);
private:
    class UniversalRep& updRep();
    const UniversalRep& getRep() const;
};

/// Two mobilities -- rotation and translation along the common z axis
/// of the inboard and outboard mobilizer frames.
class SimTK_SIMBODY_EXPORT MobilizedBody::Cylinder : public MobilizedBody {
public:
    Cylinder();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Cylinder(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Cylinder(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Cylinder& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Cylinder& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    class CylinderRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Cylinder, MobilizedBody);
private:
    class CylinderRep& updRep();
    const CylinderRep& getRep() const;
};

/// Two mobilities: The z axis of the parent's Mb frame is 
/// used for rotation (and that is always aligned with the M frame z axis).
/// The x axis of the *M* (outboard) frame is then used for translation;
/// that is, first we rotate around z, which moves M's x with respect to Mb's x. Then
/// we slide along the rotated x axis. The two generalized coordinates are the
/// rotation and the translation, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::BendStretch : public MobilizedBody {
public:
    BendStretch();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    BendStretch(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame);

    BendStretch& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    BendStretch& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    class BendStretchRep; // local subclass

    SimTK_PIMPL_DOWNCAST(BendStretch, MobilizedBody);
private:
    class BendStretchRep& updRep();
    const BendStretchRep& getRep() const;
};

/// Three mobilities -- z rotation and x,y translation. The generalized
/// coordinates are rotation about the shared z axis of the Mb and M
/// frame, translation along the Mb frame's x axis, and translation along
/// its y axis, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::Planar : public MobilizedBody {
public:
    Planar();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Planar(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Planar(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame);

    Planar& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Planar& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    // Friendly, mobilizer-specific access to coordinates and speeds.
    Planar& setDefaultAngle(Real a) {updDefaultQ()[0] = a; return *this;}
    Planar& setDefaultTranslation(const Vec2& r) {
        updDefaultQ().updSubVec<2>(1) = r; 
        return *this;
    }

    void setAngle      (State& s, Real        a) {updQ(s)[0]              = a;}
    void setTranslation(State& s, const Vec2& r) {updQ(s).updSubVec<2>(1) = r;}

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Vec3& updDefaultQ();
    Planar& setDefaultQ(const Vec3& v) {updDefaultQ()=v; return *this;}

    // Generic state access routines.
    const Vec3& getQ(const State&) const;
    const Vec3& getU(const State&) const;
    Vec3& updQ(State&) const;
    Vec3& updU(State&) const;
    void setQ(State& s, const Vec3& q) const {updQ(s)=q;}
    void setU(State& s, const Vec3& u) const {updU(s)=u;}

    const Vec3& getQDot(const State&) const;
    const Vec3& getUDot(const State&) const;
    const Vec3& getQDotDot(const State&) const;

    // Acceleration stage state variables
    const Vec3& getMobilizerForces(const State&) const;
    Vec3& updMobilizerForces(State&) const;

    class PlanarRep; // local subclass
    SimTK_PIMPL_DOWNCAST(Planar, MobilizedBody);
private:
    class PlanarRep& updRep();
    const PlanarRep& getRep() const;
};

/// Three mobilities -- unrestricted orientation modeled as a 1-2-3
/// body-fixed Euler angle sequence. This is singular when the middle
/// angle is 90 degrees.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT MobilizedBody::Gimbal : public MobilizedBody {
public:
    Gimbal();
    class GimbalRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Gimbal, MobilizedBody);
private:
    class GimbalRep& updRep();
    const GimbalRep& getRep() const;
};

/// Three mobilities -- unrestricted orientation modeled with a
/// quaternion which is never singular. A modeling option allows the
/// joint to use a 1-2-3 Euler sequence (identical to a Gimbal) 
/// instead.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ball : public MobilizedBody {
public:
    Ball();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ball(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ball(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Ball& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Ball& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }
    class BallRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Ball, MobilizedBody);
private:
    class BallRep& updRep();
    const BallRep& getRep() const;
};

/// Three mobilities -- coordinated rotation and translation along the
/// surface of an ellipsoid fixed to the parent (inboard) body.
/// The generalized coordinates are the same as for a Ball (Orientation)
/// joint, that is, a quaternion or 1-2-3 Euler sequence.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ellipsoid : public MobilizedBody {
public:
    // The ellipsoid is placed on the mobilizer's inboard frame Mb, with
    // half-axis dimensions along Mb's x,y,z respectively.
    Ellipsoid(const Vec3& radii);
    class EllipsoidRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Ellipsoid, MobilizedBody);
private:
    class EllipsoidRep& updRep();
    const EllipsoidRep& getRep() const;
};

/// Three translational mobilities. The generalized coordinates are
/// x,y,z translations along the parent (inboard) Mb frame axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Translation : public MobilizedBody {
public:
    Translation();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Translation(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Translation(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Translation& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Translation& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    class TranslationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Translation, MobilizedBody);
private:
    class TranslationRep& updRep();
    const TranslationRep& getRep() const;
};

/// Unrestricted motion for a rigid body (six mobilities). Orientation
/// is modeled the same as for the Orientation mobilizer, that is, using
/// quaternions to avoid singularities. A modeling option exists to 
/// have the joint modeled with a 1-2-3 body fixed Euler sequence like
/// a Gimbal mobilizer. Translational generalized coordinates are
/// x,y,z translations along the Mb (inboard) axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Free : public MobilizedBody {
public:
    Free();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Free(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Free(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Free& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    Free& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    class FreeRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Free, MobilizedBody);
private:
    class FreeRep& updRep();
    const FreeRep& getRep() const;
};


// These are special "ball" and "free" joints designed to allow arbitrary orientations
// for "linear" bodies, such as a CO2 molecule consisting only of point masses arranged
// along a straight line. Such bodies have no inertia about the line and cause singularities
// in the equations of motion if attached to Orientation or Free mobilizers. Instead, use the
// LineOrientation and LineFree moblizers, making sure that the inertialess direction is
// along the outboard body's z axis (that is, Mz). These mobilizers introduce only two
// mobilities (generalized speeds u), being incapable of representing non-zero angular
// velocity of M in Mb about Mz. The generalized speeds are in fact the wx and wy 
// components of w_MbM_M, that is, the x and y components of the angular velocity of M
// in Mb *expressed in M*. However, at least three generalized coordinates (q's)
// are required to represent the orientation. By default we use four quaternions for
// unconditional stability. Alternatively, you can request a 1-2-3 body fixed 
// Euler angle sequence (that is, about x, then new y, then new z) which will
// suffer a singularity when the y rotation is 90 degrees since that aligns the
// first rotation axis (x) with the last (z) which is the inertialess direction.

/// Two mobilities, representing unrestricted orientation for a body which is
/// inertialess along its own z axis. The generalized coordinates are the same
/// as for the general Orientation (Ball) mobilizer, but there are only
/// two generalized speeds. These are the x,y components of the angular velocity
/// of frame M in Mb, but expressed in the *M* (outboard frame).
class SimTK_SIMBODY_EXPORT MobilizedBody::LineOrientation : public MobilizedBody {
public:
    LineOrientation();
    class LineOrientationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(LineOrientation, MobilizedBody);
private:
    class LineOrientationRep& updRep();
    const LineOrientationRep& getRep() const;
};

/// Five mobilities, representing unrestricted motion for a body which is
/// inertialess along its own z axis. The rotational generalized coordinates are the same
/// as for the LineOrientation mobilizer. The translational coordinates are
/// the same as in a Free mobilizer, or a Cartesian (Translation) mobilizer.
class SimTK_SIMBODY_EXPORT MobilizedBody::FreeLine : public MobilizedBody {
public:
    FreeLine();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    FreeLine(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    FreeLine(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    FreeLine& setDefaultInboardFrame(const Transform& X_PMb) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PMb);
        return *this;
    }

    FreeLine& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM);
        return *this;
    }

    class FreeLineRep; // local subclass

    SimTK_PIMPL_DOWNCAST(FreeLine, MobilizedBody);
private:
    class FreeLineRep& updRep();
    const FreeLineRep& getRep() const;
};

/// Zero mobilities. This degenerate "mobilizer" serves only to weld together
/// the M frame of a body to the Mb frame on its parent.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT MobilizedBody::Weld : public MobilizedBody {
public:
    Weld();
    class WeldRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Weld, MobilizedBody);
private:
    class WeldRep& updRep();
    const WeldRep& getRep() const;
};


/// This is a special type of "mobilized" body used as a placeholder for Ground
/// in the 0th slot for a MatterSubsystem's mobilized bodies.
/// The body type will also be Ground.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ground : public MobilizedBody {
public:
    Ground();
    class GroundRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Ground, MobilizedBody);
private:
    class GroundRep& updRep();
    const GroundRep& getRep() const;
};

/// 1-6 mobilities. TODO: this will be an abstract class with virtual methods
/// that a user's derived class can implement to define a custom mobilizer.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT MobilizedBody::Custom : public MobilizedBody {
public:
    Custom(int nMobilities, int nCoordinates);
    class CustomRep; // local subclass

    // Get calculations through Stage::Instance from State.
    virtual void calcTransform(const State&, const Vector& q, 
                               Transform& X_MbM) const = 0;
    //TODO: should H be a nuX2 Matrix_<Vec3> instead? or Vector_<SpatialVec>?
    //      or nuX6 Matrix?
    virtual void calcTransitionMatrix(const State& s, 
                               Vector_<SpatialRow>& H_MbM) const = 0;
    virtual void calcTransitionMatrixTimeDerivative(const State& s,  
                               Vector_<SpatialRow>& H_MbM_Dot) const = 0;

    // get q and calculations through Stage::Position from State if needed
    virtual void calcQDot(const State&, const Vector& u, Vector& qdot) const {
        qdot = u; //TODO: only if sizes match
    }
    // get q,u and calculations through Stage::Dynamics from State if needed
    virtual void calcQDotDot(const State&, const Vector& udot, Vector& qdotdot) const {
        qdotdot = udot; //TODO: only if sizes match
    }

    SimTK_PIMPL_DOWNCAST(Custom, MobilizedBody);
private:
    class CustomRep& updRep();
    const CustomRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_H_



