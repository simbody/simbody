#ifndef SimTK_SIMBODY_MOBILIZER_H_
#define SimTK_SIMBODY_MOBILIZER_H_

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
 * This defines the Mobilizer base class, the built-in Mobilizers, and
 * the abstract UserMobilizer class.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

/**
 * This is the base class for all Mobilizers, which is just a handle for the underlying
 * hidden implementation. Each built-in Mobilizer type is a local subclass within
 * Mobilizer, and is also derived from Mobilizer.
 */
class SimTK_SIMBODY_EXPORT Mobilizer {
public:
    Mobilizer() : rep(0) { }
    ~Mobilizer();
    Mobilizer(const Mobilizer&);
    Mobilizer& operator=(const Mobilizer&);

    // These are the built-in Mobilizer types. Types on the same line are
    // synonymous.
    class Pin;         typedef Pin    Torsion;
    class Slider;      typedef Slider Prismatic;
    class Universal;
    class Cylinder;
    class BendStretch;
    class Planar;
    class Gimbal;
    class Orientation; typedef Orientation Ball, Spherical;
    class Translation; typedef Translation Cartesian;
    class Free;
    class LineOrientation;
    class FreeLine;
    class Weld;
    class Screw;
    class User;

    // Is this handle the owner of this rep? This is true if the
    // handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    class MobilizerRep; // local subclass
    explicit Mobilizer(class MobilizerRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const MobilizerRep& getRep() const {assert(rep); return *rep;}
    MobilizerRep&       updRep() const {assert(rep); return *rep;}
	void setRep(MobilizerRep& r) {assert(!rep); rep = &r;}
protected:
    class MobilizerRep* rep;
};

/// One mobility -- rotation about the common z axis of the inboard
/// and outboard mobilizer frames.
/// Synonym: Torsion
class SimTK_SIMBODY_EXPORT Mobilizer::Pin : public Mobilizer {
public:
    Pin();
    class PinRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Pin, Mobilizer);
private:
    class PinRep& updRep();
    const PinRep& getRep() const;
};

/// One mobility -- translation along the common x axis of the
/// inboard and outboard mobilizer frames.
/// Synonym: Prismatic
class SimTK_SIMBODY_EXPORT Mobilizer::Slider : public Mobilizer {
public:
    Slider();
    class SliderRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Slider, Mobilizer);
private:
    class SliderRep& updRep();
    const SliderRep& getRep() const;
};

/// One mobility -- coordinated rotation and translation along the
/// common z axis of the inboard and outboard mobilizer frames. A
/// "pitch" is specified relating the two. The generalized coordinate
/// q is the rotation angle in radians, the translation is always
/// pitch*q.
class SimTK_SIMBODY_EXPORT Mobilizer::Screw : public Mobilizer {
public:
    Screw(Real pitch);
    class ScrewRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Screw, Mobilizer);
private:
    class ScrewRep& updRep();
    const ScrewRep& getRep() const;
};

/// Two mobilities -- rotation about the x axis, followed by a rotation
/// about the new y axis. This mobilizer is badly behaved when the
/// second rotation is near 90 degrees.
class SimTK_SIMBODY_EXPORT Mobilizer::Universal : public Mobilizer {
public:
    Universal();
    class UniversalRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Universal, Mobilizer);
private:
    class UniversalRep& updRep();
    const UniversalRep& getRep() const;
};

/// Two mobilities -- rotation and translation along the common z axis
/// of the inboard and outboard mobilizer frames.
class SimTK_SIMBODY_EXPORT Mobilizer::Cylinder : public Mobilizer {
public:
    Cylinder();
    class CylinderRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Cylinder, Mobilizer);
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
class SimTK_SIMBODY_EXPORT Mobilizer::BendStretch : public Mobilizer {
public:
    BendStretch();
    class BendStretchRep; // local subclass

    SimTK_PIMPL_DOWNCAST(BendStretch, Mobilizer);
private:
    class BendStretchRep& updRep();
    const BendStretchRep& getRep() const;
};

/// Three mobilities -- x,y translation and a z rotation. The generalized
/// coordinates are rotation about the shared z axis of the Mb and M
/// frame, translation along the Mb frame's x axis, and translation along
/// its y axis, in that order.
class SimTK_SIMBODY_EXPORT Mobilizer::Planar : public Mobilizer {
public:
    Planar();
    class PlanarRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Planar, Mobilizer);
private:
    class PlanarRep& updRep();
    const PlanarRep& getRep() const;
};

/// Three mobilities -- unrestricted orientation modeled as a 1-2-3
/// body-fixed Euler angle sequence. This is singular when the middle
/// angle is 90 degrees.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT Mobilizer::Gimbal : public Mobilizer {
public:
    Gimbal();
    class GimbalRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Gimbal, Mobilizer);
private:
    class GimbalRep& updRep();
    const GimbalRep& getRep() const;
};

/// Three mobilities -- unrestricted orientation modeled with a
/// quaternion which is never singular. A modeling option allows the
/// joint to use a 1-2-3 Euler sequence (identical to a Gimbal) 
/// instead.
class SimTK_SIMBODY_EXPORT Mobilizer::Orientation : public Mobilizer {
public:
    Orientation();
    class OrientationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Orientation, Mobilizer);
private:
    class OrientationRep& updRep();
    const OrientationRep& getRep() const;
};

/// Three translational mobilities. The generalized coordinates are
/// x,y,z translations along the parent (inboard) Mb frame axes.
class SimTK_SIMBODY_EXPORT Mobilizer::Translation : public Mobilizer {
public:
    Translation();
    class TranslationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Translation, Mobilizer);
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
class SimTK_SIMBODY_EXPORT Mobilizer::Free : public Mobilizer {
public:
    Free();
    class FreeRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Free, Mobilizer);
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
class SimTK_SIMBODY_EXPORT Mobilizer::LineOrientation : public Mobilizer {
public:
    LineOrientation();
    class LineOrientationRep; // local subclass

    SimTK_PIMPL_DOWNCAST(LineOrientation, Mobilizer);
private:
    class LineOrientationRep& updRep();
    const LineOrientationRep& getRep() const;
};

/// Five mobilities, representing unrestricted motion for a body which is
/// inertialess along its own z axis. The rotational generalized coordinates are the same
/// as for the LineOrientation mobilizer. The translational coordinates are
/// the same as in a Free mobilizer, or a Cartesian (Translation) mobilizer.
class SimTK_SIMBODY_EXPORT Mobilizer::FreeLine : public Mobilizer {
public:
    FreeLine();
    class FreeLineRep; // local subclass

    SimTK_PIMPL_DOWNCAST(FreeLine, Mobilizer);
private:
    class FreeLineRep& updRep();
    const FreeLineRep& getRep() const;
};

/// Zero mobilities. This degenerate "mobilizer" serves only to weld together
/// the M frame of a body to the Mb frame on its parent.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT Mobilizer::Weld : public Mobilizer {
public:
    Weld();
    class WeldRep; // local subclass

    SimTK_PIMPL_DOWNCAST(Weld, Mobilizer);
private:
    class WeldRep& updRep();
    const WeldRep& getRep() const;
};

/// 1-6 mobilities. TODO: this will be an abstract class with virtual methods
/// that a user's derived class can implement to define a custom mobilizer.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT Mobilizer::User : public Mobilizer {
public:
    User(int nMobilities, int nCoordinates);
    class UserRep; // local subclass

    virtual void calcTransform(const State& s, const Vector& q, 
                               Transform& X_MbM) const = 0;
    virtual void calcTransitionMatrix(const State& s, 
                               Vector_<SpatialRow>& H_MbM) const = 0;
    virtual void calcTransitionMatrixTimeDerivative(const State& s,  
                               Vector_<SpatialRow>& H_MbM_Dot) const = 0;

    SimTK_PIMPL_DOWNCAST(User, Mobilizer);
private:
    class UserRep& updRep();
    const UserRep& getRep() const;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZER_H_



