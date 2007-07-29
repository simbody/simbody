#ifndef SimTK_SIMMATRIX_ORIENTATION_H_
#define SimTK_SIMMATRIX_ORIENTATION_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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
 *
 * These are numerical utility classes for dealing with the relative orientations
 * of geometric objects. These build on the basic arithmetic classes for small
 * vectors and matrices.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/Constants.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"

#include <iostream>

// Some handy conversion constants. Use the templatized inline conversion
// routines instead whenever possible.
// These are defined so that you can multiply by them. For example, if you have
// an angle in radians ar and want the same angle in degrees, write ad = ar*SimTK_RTD.
// Note that the expressions here will always be evaluated at compile time, yielding
// long double results which you can cast to smaller sizes if you want.

#define SimTK_RTD (180/SimTK_PI)
#define SimTK_DTR (SimTK_PI/180)

namespace SimTK {

class CoordinateAxis {
public:
    class X; class Y; class Z;
    operator int() const {return axis;}
protected:
    class XType{}; class YType{}; class ZType{};

    CoordinateAxis(const XType&) : axis(0) { }
    CoordinateAxis(const YType&) : axis(1) { }
    CoordinateAxis(const ZType&) : axis(2) { }
private:
    int axis;
};

class CoordinateAxis::X : public CoordinateAxis {
  public: X() : CoordinateAxis(XType()) { }
};
class CoordinateAxis::Y : public CoordinateAxis {
  public: Y() : CoordinateAxis(YType()) { }
};
class CoordinateAxis::Z : public CoordinateAxis {
  public: Z() : CoordinateAxis(ZType()) { }
};

// Predefine constants XAxis, YAxis, ZAxis. These implicitly
// convert to integers 0, 1, 2 respectively.
static const CoordinateAxis::X XAxis;
static const CoordinateAxis::Y YAxis;
static const CoordinateAxis::Z ZAxis;

// Templatized conversion routines. These should be used only for "precisions", i.e.,
// float, double, long double.
template <class P> P inline static convertRadiansToDegrees(const P& rad) {
    return rad*P(SimTK_RTD);
}
template <class P> P inline static convertDegreesToRadians(const P& deg) {
    return deg*P(SimTK_DTR);
}

template <int S> class UnitVec;
template <int S> class UnitRow;
class Rotation;
class InverseRotation;
class Transform;
class InverseTransform;

typedef UnitVec<1> UnitVec3;

/**
 * This class is a Vec3 plus an ironclad guarantee either that:
 *      - the length is one (to within a very small tolerance), or
 *      - all components are NaN.
 * Thus it is a pure direction.
 */
template <int S>
class UnitVec : public Vec<3,Real,S> {
public:
    typedef Vec<3,Real,S> BaseVec;
    typedef UnitRow<S>    TransposeType;

    UnitVec() : BaseVec(NaN) { }

    // Copy constructor.
    UnitVec(const UnitVec& u) 
      : BaseVec(static_cast<const BaseVec&>(u)) { }

    // Automatic conversion from UnitVec with different stride.
    template <int S2> UnitVec(const UnitVec<S2>& u)
        : BaseVec(static_cast<const typename UnitVec<S2>::BaseVec&>(u)) { }

    // Explicit conversion from Vec to UnitVec, requiring expensive normalization.
    explicit UnitVec(const BaseVec& v) : BaseVec(v/v.norm()) { }
    template <int S2> explicit UnitVec(const Vec<3,Real,S2>& v)
        : BaseVec(v/v.norm()) { }

    UnitVec(const Real& x, const Real& y, const Real& z) : BaseVec(x,y,z) {
        static_cast<BaseVec&>(*this) /= BaseVec::norm();
    }

    // Create a unit axis vector 100 010 001
    explicit UnitVec(int axis) : BaseVec(0) {
        assert(0 <= axis && axis <= 2);
        BaseVec::operator[](axis) = 1.;
    }

    UnitVec& operator=(const UnitVec& u) {
        BaseVec::operator=(static_cast<const BaseVec&>(u)); 
        return *this;
    }
    template <int S2> UnitVec& operator=(const UnitVec<S2>& u) {
        BaseVec::operator=(static_cast<const typename UnitVec<S2>::BaseVec&>(u));
        return *this;
    }

    const BaseVec& asVec3() const {return static_cast<const BaseVec&>(*this);}

    // Override Vec3 methods which preserve length. These return the 
    // packed UnitVec regardless of our stride.
    UnitVec3 negate()    const {return UnitVec3(-asVec3(),true);}
    UnitVec3 operator-() const {return negate();}

    const TransposeType& operator~() const {
        return *reinterpret_cast<const TransposeType*>(this);
    }
    TransposeType& operator~() {
        return *reinterpret_cast<TransposeType*>(this);
    }

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitVec!
    const Real& operator[](int i) const {return BaseVec::operator[](i);}
    const Real& operator()(int i) const {return BaseVec::operator()(i);}

    // Return a vector whose measure numbers are the absolute values
    // of the ones here. This will still have unit length but will be
    // a reflection of this unit vector into the first octant (+x,+y,+z).
    // Note that we are returning the packed form of UnitVec regardless
    // of our stride here.
    UnitVec3 abs() const {
        return UnitVec3(asVec3().abs(),true);
    }

    // Return a unit vector perpendicular to this one (arbitrary).
    inline UnitVec3 perp() const;

    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitVec(const BaseVec& v, bool) : BaseVec(v) { }
    template <int S2> UnitVec(const Vec<3,Real,S2>& v, bool) : BaseVec(v) { }
};

template <int S>
inline UnitVec3 UnitVec<S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitVec3 u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Vec3 result which is then normalized.
    return UnitVec3(*this % UnitVec3(minAxis));
}

/**
 * This type is used for the transpose of UnitVec, and as the returned row
 * type of a Rotation. Don't construct these directly.
 */
template <int S>
class UnitRow : public Row<3,Real,S> {
public:
    typedef Row<3,Real,S> BaseRow;
    typedef UnitVec<S>    TransposeType;

    UnitRow() : BaseRow(NaN) { }

    // Copy constructor.
    UnitRow(const UnitRow& u) 
      : BaseRow(static_cast<const BaseRow&>(u)) { }

    // Automatic conversion from UnitRow with different stride.
    template <int S2> UnitRow(const UnitRow<S2>& u)
        : BaseRow(static_cast<const typename UnitRow<S2>::BaseRow&>(u)) { }

    // Copy assignment.
    UnitRow& operator=(const UnitRow& u) {
        BaseRow::operator=(static_cast<const BaseRow&>(u)); 
        return *this;
    }
    // Assignment from UnitRow with different stride.
    template <int S2> UnitRow& operator=(const UnitRow<S2>& u) {
        BaseRow::operator=(static_cast<const typename UnitRow<S2>::BaseRow&>(u));
        return *this;
    }

    // Explicit conversion from Row to UnitRow, requiring expensive normalization.
    explicit UnitRow(const BaseRow& v) : BaseRow(v/v.norm()) { }
    template <int S2> explicit UnitRow(const Row<3,Real,S2>& v)
        : BaseRow(v/v.norm()) { }

    UnitRow(const Real& x, const Real& y, const Real& z) : BaseRow(x,y,z) {
        static_cast<BaseRow&>(*this) /= BaseRow::norm();
    }

    // Create a unit axis vector 100 010 001
    explicit UnitRow(int axis) : BaseRow(0) {
        assert(0 <= axis && axis <= 2);
        BaseRow::operator[](axis) = 1.;
    }

    const BaseRow& asRow3() const {return static_cast<const BaseRow&>(*this);}

    // Override Row3 methods which preserve length. These return the 
    // packed UnitRow regardless of our stride.
    UnitRow<1> negate()    const {return UnitRow<1>(-asRow3(),true);}
    UnitRow<1> operator-() const {return negate();}

    const TransposeType& operator~() const {
        return *reinterpret_cast<const TransposeType*>(this);
    }
    TransposeType& operator~() {
        return *reinterpret_cast<TransposeType*>(this);
    }

    // We have to define these here so that the non-const ones won't be
    // inherited. We don't trust anyone to write on one element of a UnitRow!
    const Real& operator[](int i) const {return BaseRow::operator[](i);}
    const Real& operator()(int i) const {return BaseRow::operator()(i);}

    // Return a vector whose measure numbers are the absolute values
    // of the ones here. This will still have unit length but will be
    // a reflection of this unit vector into the first octant (+x,+y,+z).
    // Note that we are returning the packed form of UnitVec regardless
    // of our stride here.
    UnitRow<1> abs() const {
        return UnitRow<1>(asRow3().abs(),true);
    }

    // Return a unit row vector perpendicular to this one (arbitrary).
    inline UnitRow<1> perp() const;

    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitRow(const BaseRow& v, bool) : BaseRow(v) { }
    template <int S2> UnitRow(const Row<3,Real,S2>& v, bool) : BaseRow(v) { }
};

template <int S>
inline UnitRow<1> UnitRow<S>::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitRow<1> u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Row3 result which is then normalized.
    return UnitRow<1>(*this % UnitRow<1>(minAxis));
}

/**
 * A Quaternion is a Vec4 with the following behavior:
 *   - its length is always 1 (or else it is all NaN)
 *   - it is equivalent to an angle/axis rotation for
 *     angle a, axis unit vector v, like this:
 *        q=[ cos(a/2) sin(a/2)*v ]
 * We consider a quaternion to be in "canonical form" when its
 * first element is nonnegative. That corresponds to rotation
 * angles in the range -180 < a <= 180 degrees. We don't require
 * quaternions to be in canonical form; continuity during integration
 * requires them to range more widely. However, when we're creating
 * them from scratch and have a choice, we'll return them in 
 * canonical form.
 *
 * The (angle,axis) form is handled here also. When these are in
 * their canonical form, they have -180 < angle <= 180 and |axis|=1.
 * However, (angle,axis) is meaningful for any value of angle and
 * for any axis where |axis| > 0.
 */
class Quaternion : public Vec4 {
public:
    typedef Vec4 BaseVec;

    /// Default constructor produces a 0-rotation (identity) quaternion
    Quaternion() : Vec4(1,0,0,0) { }

    /// Initialize this quaternion from an unnormalized quaternion
    /// stored in a Vec4. The argument is *not* interpreted as an
    /// (angle,axis) -- it is just a slightly mangled quaternion.
    /// If the passed-in vector is *exactly* zero, we will assume
    /// it indicates a "zero rotation" and set the Quaternion to
    /// [1 0 0 0]. If the length is 0 < len < eps (eps being machine
    /// tolerance) we consider that an error condition and set the
    /// Quaternion to NaN. Otherwise we normalize it and return.
    /// The constructed quaternion is NOT put in canonical form -- it is
    /// as close to the original as possible.
    /// Because of the normalization, this costs about 40 flops.
    explicit Quaternion(const Vec4& v) {
        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real len = v.norm();
        if      (len == 0)  setToZero();
        else if (len < eps) setToNaN();
        else BaseVec::operator=(v/len);
    }

    /// Initialize this quaternion from a rotation matrix. The result
    /// will be in canonical form. The cost is about 60 flops.
    SimTK_SimTKCOMMON_EXPORT explicit Quaternion(const Rotation&);

    /// Copy constructor copies the source as-is; it does not 
    /// convert to canonical form, or normalize, or anything else.
    /// Zero cost.
    Quaternion(const Quaternion& q) : BaseVec(q) { }

    /// Copy assignment copies the source as-is; it does not 
    /// convert to canonical form, or normalize, or anything else.
    /// Zero cost.
    Quaternion& operator=(const Quaternion& q) {
        BaseVec::operator=(q.asVec4());
        return *this;
    }
    
    /// By zero here we mean "zero rotation", i.e., an identity rotation
    /// represented as [1 0 0 0]. This is in canonical form; [-1 0 0 0] would
    /// mean the same thing.
    void setToZero() {
        BaseVec::operator=(Vec4(1,0,0,0));
    }

    /// This is the only exception to the "must be normalized" rule for
    /// quaternions -- all elements are set to NaN. Note that unlike 
    /// naked Vec4's, Quaternions do not start out NaN even in Debug mode.
    /// The default constructor sets the Quaternion to "zero rotation" instead.
    void setToNaN() {
        BaseVec::setToNaN();
    }

    /// Resulting 4-vector is [ a vx vy vz ] with (a,v) in canonical form.
    /// That is, -180 < a <= 180 and |v|=1. The cost of this operation is
    /// roughly one atan2, one sqrt, and one divide, say about 100 flops.
    SimTK_SimTKCOMMON_EXPORT Vec4 convertToAngleAxis() const;

    /// Assign the current quaternion to the rotation represented by the 
    /// passed-in (angle,axis) form. The resulting quaternion will be in
    /// canonical form regardless of the condition of the (angle,axis) input.
    /// The "axis" will be normalized here unless it has zero length on
    /// entry, in which case the quaternion will be all-NaN.
    /// Cost is a normalization, a sin and a cos, or about 120 flops.
    SimTK_SimTKCOMMON_EXPORT void setToAngleAxis(const Vec4& av);

    /// Assign the current quaternion to the rotation represented by the 
    /// passed-in (angle,unitVector) form. The resulting quaternion will be in
    /// canonical form regardless of the condition of the (angle,axis) input.
    /// This can't fail for any angle since we know we have a good axis.
    /// Cost is one sin, one cos or roughly 80 flops.
    SimTK_SimTKCOMMON_EXPORT void setToAngleAxis(const Real& a, const UnitVec3& v);

    /// Upcast this Quaternion to its parent class, a Vec4. This is inline
    /// and should generate no code. You can do the same thing with static_cast
    /// if you prefer. Zero cost.
    const BaseVec& asVec4() const {
        return *static_cast<const BaseVec*>(this);
    }

    /// Don't use this unless you are *sure* this is already normalized! This
    /// is much faster than the normal Quaternion(Vec4) constructor which
    /// expects the Vec4 to need cleaning up. The second argument here is just
    /// to allow you to force a call to the fast constructor; it is otherwise
    /// ignored. By convention however, you should call this with the second
    /// argument set to "true". Zero cost.
    Quaternion(const Vec4& v, bool) : Vec4(v) { }
};

/**
 * This class is a Mat33 plus an ironclad guarantee that the matrix represents
 * a pure rotation. Having this as a separate class allows us both to ensure
 * that we have a legitimate rotation matrix, and to take advantage of
 * knowledge about the properties of these matrices. For example, multiplication
 * by a rotation matrix preserves the length of a vector so unit vectors
 * are still unit vectors afterwards and don't need to be normalized again.
 * 
 * A rotation is an orthogonal matrix whose columns and rows
 * are directions (that is, unit vectors) which are mutually orthogonal. 
 * Furthermore, if the columns (or rows) are labeled x,y,z it always holds
 * that z = x X y (rather than -(x X y)) ensuring that this is a rotation
 * and not a reflection.
 *
 * A rotation matrix is formed from a Cartesian coordinate frame
 * simply by using the x,y,z axes of the coordinate frame as
 * the columns of the rotation matrix. That is, if you have 
 * a frame F with its axes column vectors x,y,z expressed in frame G, you can write
 * this as a matrix  rot_GF = [ x y z ] which when applied to
 * a vector with measure numbers expressed in F yields that same
 * vector re-expressed in G: v_G = rot_GF * v_F. Because a rotation
 * is orthogonal its transpose is its inverse, which is also a rotation.
 * So we write rot_FG = ~rot_GF (where "~" is SimTK for "transpose"), and
 * this matrix can be used to rotate in the other direction: 
 *      v_F = rot_FG * v_G = ~rot_GF * v_G.
 */
class Rotation : public Mat33 {
public:
    typedef Mat33 BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    Rotation() : BaseMat(1.) { }    // default is identity
    // default copy constructor, copy assignment and destructor

    // Like copy constructor and copy assign but for inverse rotation;
    // constructor is also implicit conversion from InverseRotation to Rotation.
    inline Rotation(const InverseRotation&);
    inline Rotation& operator=(const InverseRotation&);

    // This constructor takes a Mat33 and orthogonalizes it into a (hopefully
    // nearby) rotation matrix.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const Mat33&);

    // These static methods are like constructors with friendlier names.

    static Rotation zero() {return Rotation();}
    static Rotation NaN()  {Rotation r;r.setToNaN();return r;}

    // One-angle rotations.
    static Rotation aboutX(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( 1 ,  0 ,  0 ,
                         0 ,  c , -s ,
                         0 ,  s ,  c );
    }
    static Rotation aboutY(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( c ,  0 ,  s ,
                         0 ,  1 ,  0 ,
                        -s ,  0 ,  c );
    }
    static Rotation aboutZ(const Real& angInRad) {
        const Real s = std::sin(angInRad), c = std::cos(angInRad);
        return Rotation( c , -s ,  0 ,
                         s ,  c ,  0 ,
                         0 ,  0 ,  1 );
    }

    SimTK_SimTKCOMMON_EXPORT static Rotation aboutAxis
       (const Real& angInRad, const UnitVec3& axis);

    static Rotation aboutAxis(const Real& angInRad, const Vec3& axis)
      { return aboutAxis(angInRad, UnitVec3(axis)); }

    // Two-angle space-fixed rotations.
    static Rotation aboutXThenOldY(const Real& xInRad, const Real& yInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        return Rotation( cy   ,  sy*sx ,  sy*cx ,
                          0   ,   cx   ,  -sx   ,
                        -sy   ,  cy*sx ,  cy*cx );
    }
    static Rotation aboutYThenOldX(const Real& yInRad, const Real& xInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        return Rotation(  cy   ,    0   ,   sy   ,
                         sy*sx ,   cx   , -cy*sx ,
                        -sy*cx ,   sx   ,  cy*cx );
    }
    static Rotation aboutXThenOldZ(const Real& xInRad, const Real& zInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation(  cz   , -sz*cx ,  sz*sx ,
                          sz   ,  cz*cx , -cz*sx ,
                           0   ,   sx   ,   cx   );
    }
    static Rotation aboutZThenOldX(const Real& zInRad, const Real& xInRad) {
        const Real sx = std::sin(xInRad), cx = std::cos(xInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation(  cz   ,  -sz   ,    0   ,
                         sz*cx ,  cz*cx ,  -sx   ,
                         sz*sx ,  cz*sx ,   cx   );
    }
    static Rotation aboutYThenOldZ(const Real& yInRad, const Real& zInRad) {
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation( cz*cy ,  -sz   ,  cz*sy ,
                         sz*cy ,   cz   ,  sz*sy ,
                         -sy   ,    0   ,   cy   );
    }
    static Rotation aboutZThenOldY(const Real& zInRad, const Real& yInRad) {
        const Real sy = std::sin(yInRad), cy = std::cos(yInRad);
        const Real sz = std::sin(zInRad), cz = std::cos(zInRad);
        return Rotation( cz*cy , -sz*cy ,   sy   ,
                          sz   ,   cz   ,    0   ,
                        -cz*sy ,  sz*sy ,   cy   );
    }

    // Two-angle body fixed rotations. These are the same as the
    // reversed space-fixed ones.
    static Rotation aboutXThenNewY(const Real& xInRad, const Real& yInRad)
      { return aboutYThenOldX(yInRad, xInRad); }
    static Rotation aboutYThenNewX(const Real& yInRad, const Real& xInRad)
      { return aboutXThenOldY(xInRad, yInRad); }
    static Rotation aboutXThenNewZ(const Real& xInRad, const Real& zInRad)
      { return aboutZThenOldX(zInRad, xInRad); }
    static Rotation aboutZThenNewX(const Real& zInRad, const Real& xInRad)
      { return aboutXThenOldZ(xInRad, zInRad); }
    static Rotation aboutYThenNewZ(const Real& yInRad, const Real& zInRad)
      { return aboutZThenOldY(zInRad, yInRad); }
    static Rotation aboutZThenNewY(const Real& zInRad, const Real& yInRad)
      { return aboutYThenOldZ(yInRad, zInRad); }

    /// Create a Rotation matrix by specifying only its z axis. 
    /// The resulting x and y axes will be appropriately perpendicular
    /// but are otherwise arbitrary. This will work for any stride
    /// UnitVec because there is always an implicit conversion
    /// available to the packed form used as the argument.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const UnitVec3& z);

    /// Create a rotation matrix from a quaternion.
    SimTK_SimTKCOMMON_EXPORT explicit Rotation(const Quaternion&);

    /// Set this Rotation to represent a rotation of +q0 about
    /// the body frame's Z axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW X axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
    SimTK_SimTKCOMMON_EXPORT void setToBodyFixed321(const Vec3&);

    /// Set this Rotation to represent a rotation of +q0 about
    /// the body frame's X axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW Z axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 1-2-3.
    SimTK_SimTKCOMMON_EXPORT void setToBodyFixed123(const Vec3&);

    /// Set this Rotation to represent the same rotation as
    /// the passed-in quaternion.
    SimTK_SimTKCOMMON_EXPORT void setToQuaternion(const Quaternion&);

    /// Convert this Rotation matrix to the equivalent quaternion.
    SimTK_SimTKCOMMON_EXPORT Quaternion convertToQuaternion() const;

    /// Convert this Rotation matrix to the equivalent 1-2-3 body fixed
    /// Euler angle sequence.
    SimTK_SimTKCOMMON_EXPORT Vec3 convertToBodyFixed123() const;

    /// Convert this Rotation matrix to the equivalent 1-2 body fixed
    /// Euler angle sequence. The result is only meaningful if the
    /// Rotation matrix is one that can be produced by such a sequence.
    SimTK_SimTKCOMMON_EXPORT Vec2 convertToBodyFixed12() const;

    /// Convert this Rotation matrix to the equivalent 1-2 space fixed
    /// Euler angle sequence. The result is only meaningful if the
    /// Rotation matrix is one that can be produced by such a sequence.
    SimTK_SimTKCOMMON_EXPORT Vec2 convertToSpaceFixed12() const;

    /// Convert this Rotation matrix to an equivalent (angle,axis)
    /// representation: (a vx vy vz), with v a unit vector, a in radians.
    SimTK_SimTKCOMMON_EXPORT Vec4 convertToAngleAxis()  const;

    /// Return true if the passed-in Rotation is identical to the current one
    /// to within a cone whose half-angle is supplied (in radians). That half-angle
    /// is called the "pointing error", and we require 0 < pointing error < pi since
    /// 0 would require perfect precision and pi (180 degrees) would return true
    /// for any pair of Rotations.
    SimTK_SimTKCOMMON_EXPORT bool isSameRotationToWithinAngle(const Rotation&, Real okPointingError) const;

    /// Return true if the passed-in Rotation is identical to the current one
    /// to within what can be expected at the precision being used.
    SimTK_SimTKCOMMON_EXPORT bool isSameRotationToMachinePrecision(const Rotation&) const;


    Rotation(const Rotation& R) : BaseMat(R) { }
    Rotation& operator=(const Rotation& R) {
        BaseMat::operator=(R.asMat33()); return *this;
    }

    /// By zero we mean "zero rotation", i.e., an identity matrix.
    Rotation& setToZero() {
        BaseMat::operator=(Real(1)); return *this;
    }

    Rotation& setToNaN() {
        BaseMat::setToNaN(); return *this;
    }

    /// Set this Rotation to represent a rotation of +q radians
    /// around the base frame's 0,0,1 axis.
    void setToRotationAboutZ(const Real& q) {
        const Real sq = std::sin(q), cq = std::cos(q);
        BaseMat::operator=(BaseMat( cq , -sq , 0. ,
                                    sq ,  cq , 0. ,
                                    0. ,  0. , 1. ));
    }

    /// Set this Rotation to represent a rotation of +q0 about
    /// the base frame's X axis, followed by a rotation of +q1 about
    /// the base frame's (unchanged) Y axis.
    void setToSpaceFixed12(const Vec2& q) {
        const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
        const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
        BaseMat::operator=(BaseMat( cq1 , sq1*sq0 , sq1*cq0,
                                     0. ,   cq0   ,   -sq0 ,
                                   -sq1 , cq1*sq0 , cq1*cq0));
    }


    /// Given Euler angles forming a body-fixed 3-2-1 sequence, and the relative
    /// angular velocity vector of B in the parent frame, *BUT EXPRESSED IN
    /// THE BODY FRAME*, return the Euler angle
    /// derivatives. You are dead if q[1] gets near 90 degrees!
    /// See Kane's Spacecraft Dynamics, page 428, body-three: 3-2-1.
    static Vec3 convertAngVelToBodyFixed321Dot(const Vec3& q, const Vec3& w_PB_B) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1 = Real(1)/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E( 0. ,   s2oc1  ,  c2oc1  ,
                       0. ,     c2   ,   -s2   ,
                       1. , s1*s2oc1 , s1*c2oc1 );
        return E * w_PB_B;
    }


    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed321DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv(  -s1  ,  0. , 1. ,
                          c1*s2 ,  c2 , 0. ,
                          c1*c2 , -s2 , 0. );
        return Einv*qd;
    }
    
    /// Given Euler angles forming a body-fixed 1-2-3 sequence, and the relative
    /// angular velocity vector of B in the parent frame,  *BUT EXPRESSED IN
    /// THE BODY FRAME*, return the Euler angle
    /// derivatives. You are dead if q[1] gets near 90 degrees!
    /// See Kane's Spacecraft Dynamics, page 427, body-three: 1-2-3.
    static Vec3 convertAngVelToBodyFixed123Dot(const Vec3& q, const Vec3& w) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E(    c2oc1  , -s2oc1  , 0.,
                            s2   ,    c2   , 0.,
                       -s1*c2oc1 , s1*s2oc1, 1. );
        return E*w;
    }

    /// Inverse of the above routine. Returned angular velocity is B in P,
    /// expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed123DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 Einv( c1*c2 ,  s2 , 0. ,
                         -c1*s2 ,  c2 , 0. ,
                           s1   ,  0. , 1. );
        return Einv*qd;
    }

    // TODO: sherm: is this right? Warning: everything is measured in the
    // *PARENT* frame, but has to be expressed in the *BODY* frame.
    static Vec3 convertAngVelDotToBodyFixed321DotDot
        (const Vec3& q, const Vec3& w_PB_B, const Vec3& wdot)
    {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E( 0. ,   s2oc1  ,  c2oc1  ,
                       0. ,     c2   ,   -s2   ,
                       1. , s1*s2oc1 , s1*c2oc1 );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot( 0. ,       a           ,         b         ,
                          0. ,   -qdot[2]*s2     ,    -qdot[2]*c2    ,
                          0. , s1*a + qdot[1]*s2 , s1*b + qdot[1]*c2 );

        return E*wdot + Edot*w_PB_B;
    }

    // TODO: sherm: is this right? Warning: everything is measured in the
    // *PARENT* frame, but has to be expressed in the *BODY* frame.
    static Vec3 convertAngVelDotToBodyFixed123DotDot
        (const Vec3& q, const Vec3& w_PB_B, const Vec3& wdot)
    {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1  = 1/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E(    c2oc1  , -s2oc1  , 0.,
                            s2   ,    c2   , 0.,
                       -s1*c2oc1 , s1*s2oc1, 1. );
        const Vec3 qdot = E * w_PB_B;

        const Real t =  qdot[1]*qdot[2]*s1*ooc1;
        const Real a =  t*c2oc1; // d/dt s2oc1
        const Real b = -t*s2oc1; // d/dt c2oc1

        const Mat33 Edot(       b           ,        -a         , 0.,
                             qdot[2]*c2     ,    -qdot[2]*s2    , 0.,
                         -s1*b - qdot[1]*c2 , s1*a + qdot[1]*s2 , 0. );

        return E*wdot + Edot*w_PB_B;
    }


    /// Given a possibly unnormalized quaternion (0th element is the scalar) and the
    /// relative angular velocity vector of B in its parent, expressed 
    /// in the *PARENT*, return the quaternion derivatives. This is never singular.
    static Vec4 convertAngVelToQuaternionDot(const Vec4& q, const Vec3& w_PB_P) {
        const Mat43 E(-q[1],-q[2],-q[3],
                       q[0], q[3],-q[2],    // TODO: signs???
                      -q[3], q[0], q[1],
                       q[2],-q[1], q[0]);
        return 0.5*(E*w_PB_P);
    }

    /// Everything is measured and expressed in the parent.
    static Vec4 convertAngVelDotToQuaternionDotDot
        (const Vec4& q, const Vec3& w_PB_P, const Vec3& wdot)
    {
        const Mat43 E(-q[1],-q[2],-q[3],
                       q[0], q[3],-q[2],    // TODO: signs???
                      -q[3], q[0], q[1],
                       q[2],-q[1], q[0]);
        const Vec4  qdot = 0.5*E*w_PB_P;
        const Mat43 Edot(-qdot[1],-qdot[2],-qdot[3],
                          qdot[0], qdot[3],-qdot[2],
                         -qdot[3], qdot[0], qdot[1],
                          qdot[2],-qdot[1], qdot[0]);

        return 0.5*(Edot*w_PB_P + E*wdot);
    }

    /// Inverse of the above routine. Returned AngVel is expressed in
    /// the *PARENT* frame: w_PB_P.
    static Vec3 convertQuaternionDotToAngVel(const Vec4& q, const Vec4& qd) {
        const Mat34 Et(-q[1], q[0],-q[3], q[2],  // TODO: signs???
                       -q[2], q[3], q[0],-q[1],
                       -q[3],-q[2], q[1], q[0]);
        return 2.*(Et*qd);
    }

    const InverseRotation& invert() const {
        return *reinterpret_cast<const InverseRotation*>(this);
    }
    InverseRotation& updInvert() {
        return *reinterpret_cast<InverseRotation*>(this);
    }

    // Override the Mat33 versions of transpose.
    const InverseRotation& transpose() const {return invert();}
    InverseRotation&       updTranspose() {return updInvert();}

    // Note that this does not have unit stride.
    const RowType& row(int i) const {
        return reinterpret_cast<const RowType&>(asMat33()[i]);
    }
    const ColType& col(int j) const {
        return reinterpret_cast<const ColType&>(asMat33()(j));
    }
    const ColType& x() const {return col(0);}
    const ColType& y() const {return col(1);}
    const ColType& z() const {return col(2);}

    inline Rotation& operator*=(const Rotation& R);
    inline Rotation& operator/=(const Rotation& R);
    inline Rotation& operator*=(const InverseRotation&);
    inline Rotation& operator/=(const InverseRotation&);

    const InverseRotation& operator~() const {return invert();}
    InverseRotation&       operator~()       {return updInvert();}

    const RowType& operator[](int i) const {return row(i);}
    const ColType& operator()(int j) const {return col(j);}

    const BaseMat& asMat33() const {
        return *static_cast<const BaseMat*>(this);
    }

    /// Less efficient version of asMat33() since it copies, but you don't
    /// have to know the internal layout.
    BaseMat toMat33() const {
        return asMat33();
    }

    static Rotation trustMe(const Mat33& m) {return Rotation(m, true);}

private:
    // We're trusting that m is a rotation.
    explicit Rotation(const Mat33& m, bool) : Mat33(m) { }

    // This is only for the most trustworthy of callers, that is, methods
    // of the Rotation class. There are a lot of ways for this NOT to
    // be a legitimate rotation matrix -- be careful!!
    // Note that these are supplied in rows.
    Rotation( const Real& xx, const Real& xy, const Real& xz,
              const Real& yx, const Real& yy, const Real& yz,
              const Real& zx, const Real& zy, const Real& zz )
      : BaseMat(xx,xy,xz, yx,yy,yz, zx,zy,zz)
    {
    }
};

class InverseRotation : public Mat33::TransposeType {
public:
    typedef Mat33::TransposeType BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    // Don't construct one these; they should only occur as expression intermediates.
    // But if you must ...
    InverseRotation() : BaseMat(1) { }

    InverseRotation(const InverseRotation& R) : BaseMat(R) { }
    InverseRotation& operator=(const InverseRotation& R) {
        BaseMat::operator=(R.asMat33()); return *this;
    }

    const Rotation& invert() const {
        return *reinterpret_cast<const Rotation*>(this);
    }
    Rotation& updInvert() {
        return *reinterpret_cast<Rotation*>(this);
    }

    // Override the Mat33 versions of transpose.
    const Rotation& transpose() const {return invert();}
    Rotation&       updTranspose() {return updInvert();}

    // Note that this does not have unit stride.
    const RowType& row(int i) const {
        return reinterpret_cast<const RowType&>(asMat33()[i]);
    }
    const ColType& col(int j) const {
        return reinterpret_cast<const ColType&>(asMat33()(j));
    }
    const ColType& x() const {return col(0);}
    const ColType& y() const {return col(1);}
    const ColType& z() const {return col(2);}

    const Rotation& operator~() const {return invert();}
    Rotation&       operator~()       {return updInvert();}

    const RowType& operator[](int i) const {return row(i);}
    const ColType& operator()(int j) const {return col(j);}

    const BaseMat& asMat33() const {
        return *static_cast<const BaseMat*>(this);
    }

    /// Less efficient version of asMat33() since it copies, but you don't
    /// have to know the internal layout.
    BaseMat toMat33() const {
        return asMat33();
    }
};


SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const Rotation& m);

template <int S> inline UnitVec<1>
operator*(const Rotation& R, const UnitVec<S>& v) {
    return UnitVec<1>(R.asMat33()*v.asVec3(), true);
}
template <int S> inline UnitRow<1>
operator*(const UnitRow<S>& r, const Rotation& R) {
    return UnitRow<1>(r.asRow3(), R.asMat33(), true);
}

template <int S> inline UnitVec<1>
operator*(const InverseRotation& R, const UnitVec<S>& v) {
    return UnitVec<1>(R.asMat33()*v.asVec3(), true);
}
template <int S> inline UnitRow<1>
operator*(const UnitRow<S>& r, const InverseRotation& R) {
    return UnitRow<1>(r.asRow3(), R.asMat33(), true);
}

inline Rotation::Rotation(const InverseRotation& R)
  : Mat33(R.asMat33()) { 
}
inline Rotation& Rotation::operator=(const InverseRotation& R) {
    static_cast<BaseMat&>(*this) = R.asMat33();
    return *this;
}

inline Rotation& Rotation::operator*=(const Rotation& R) {
    static_cast<BaseMat&>(*this) *= R.asMat33();
    return *this;
}
inline Rotation& Rotation::operator/=(const Rotation& R) {
    static_cast<BaseMat&>(*this) *= (~R).asMat33();
    return *this;
}
inline Rotation& Rotation::operator*=(const InverseRotation& R) {
    static_cast<BaseMat&>(*this) *= R.asMat33();
    return *this;
}
inline Rotation& Rotation::operator/=(const InverseRotation& R) {
    static_cast<BaseMat&>(*this) *= (~R).asMat33();
    return *this;
}

inline Rotation operator*(const Rotation&        R1, const Rotation&        R2) {return Rotation(R1)*=R2;}
inline Rotation operator*(const Rotation&        R1, const InverseRotation& R2) {return Rotation(R1)*=R2;}
inline Rotation operator*(const InverseRotation& R1, const Rotation&        R2) {return Rotation(R1)*=R2;}
inline Rotation operator*(const InverseRotation& R1, const InverseRotation& R2) {return Rotation(R1)*=R2;}

inline Rotation operator/(const Rotation&        R1, const Rotation&        R2) {return Rotation(R1)/=R2;}
inline Rotation operator/(const Rotation&        R1, const InverseRotation& R2) {return Rotation(R1)/=R2;}
inline Rotation operator/(const InverseRotation& R1, const Rotation&        R2) {return Rotation(R1)/=R2;}
inline Rotation operator/(const InverseRotation& R1, const InverseRotation& R2) {return Rotation(R1)/=R2;}

/**
 * This class represents the rotate-and-shift transform which gives the 
 * location and orientation of a new frame F in a base (reference) frame
 * B. A frame is an orthogonal, right-handed set of three axes, and an
 * origin point. A transform X from frame B to F consists of 3 perpendicular
 * unit vectors defining F's axes as viewed from B (that is, as expressed in 
 * the basis formed by B's axes), and a vector from B's origin point OB to F's
 * origin point OF. Note that the meaning of "B" comes from the context in
 * which the transform is used. We use the phrase "frame F is in frame B" to
 * describe the above relationship, that is, "in" means both measured from
 * and expressed in. 
 *
 * The axis vectors constitute a Rotation. They are ordered 1-2-3 or x-y-z
 * as you prefer, with z = x X y, making a right-handed set. These axes are arranged
 * as columns of a 3x3 rotation matrix R_BF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame B and F. (The
 * columns of R_BF are F's coordinate axes, expressed in B.) For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in B is given by vB = R_BF*vF. F's origin point OF is 
 * stored as the translation vector T_BF=(OF-OB) and expressed in B.
 *
 * Transform is designed to behave as much as possible like the computer
 * graphics 4x4 transform X which would be arranged like this:
 *
 * @verbatim
 *
 *         [       |   ]
 *     X = [   R   | T ]    R is a 3x3 orthogonal rotation matrix
 *         [.......|...]    T os a 3x1 translation vector
 *         [ 0 0 0   1 ]
 *
 * @endverbatim
 *
 * These can be composed directly by matrix multiplication, but more 
 * importantly they have a particularly simple inverse:
 *
 * @verbatim
 *
 *    -1   [       |    ]
 *   X   = [  ~R   | T* ]   ~R is R transpose, T* = ~R(-T).
 *         [.......|....]
 *         [ 0 0 0   1  ] 
 *
 * @endverbatim
 *
 * This inverse is so simple that we compute it simply by defining another
 * type, InverseTransform, which is identical to Transform in memory but
 * behaves as though it contains the inverse. That way we invert just by
 * changing point of view (recasting) rather than computing.
 *
 * This is a "POD" (plain old data) class with a well-defined memory
 * layout on which a client of this class may depend: There are 
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,T.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,T1,T2,T3. It is expressly allowed
 * to reinterpret Transform objects in any appropriate manner that depends
 * on this memory layout.
 */
class Transform {
public:
    /// Default constructor gives an identity transform.
    Transform() : R_BF(),  T_BF(0) { }

    /// Combine a rotation and a translation into a transform.
    Transform(const Rotation& R, const Vec3& T) : R_BF(R), T_BF(T) { }

    /// Construct or default-convert a rotation into a transform
    /// containing that rotation and zero translation.
    Transform(const Rotation& R) : R_BF(R), T_BF(0) { }

    /// Construct or default-convert a translation (expressed as
    /// a Vec3) into a transform with that translation and a zero
    /// rotation.
    Transform(const Vec3& T) : R_BF(),  T_BF(T) { }

    // default copy, assignment, destructor

    /// Assignment from InverseTransform. This means that the 
    /// transform we're assigning to must end up with the same @em meaning
    /// as the inverse transform X has, so we'll need to end up with:
    ///   @li  T == X.T()
    ///   @li  R == X.R()
    ///
    /// Cost: one frame conversion and a negation, 18 flops.
    // (Definition is below after InverseTransform is declared.)
    inline Transform& operator=(const InverseTransform& X);

    /// Assign a new value to this transform, explicitly providing
    /// the rotation and translation separately. We return a reference
    /// to the now-modified transform as though this were an assignment
    /// operator.
    Transform& set(const Rotation& R, const Vec3& T) {
        T_BF=T; R_BF=R; return *this;
    }

    /// By zero we mean "zero transform", i.e., an identity rotation
    /// and zero translation. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform& setToZero() {
        R_BF.setToZero(); T_BF = 0.; return *this;
    }

    /// This fills both the rotation and translation with NaNs. Note: this is
    /// @em not the same as a default-constructed transform, which is a
    /// legitimate identity transform instead. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform& setToNaN() {
        R_BF.setToNaN(); T_BF.setToNaN(); return *this;
    }

    /// Return a read-only inverse of the current Transform, simply by casting it to
    /// the InverseTransform type. Zero cost.
    const InverseTransform& invert() const {
        return *reinterpret_cast<const InverseTransform*>(this);
    }

    /// Return a writable (lvalue) inverse of the current transform, simply by casting it to
    /// the InverseTransform type. That is, this is an lvalue. Zero cost.
    InverseTransform& updInvert() {
        return *reinterpret_cast<InverseTransform*>(this);
    }

    /// Overload transpose operator to mean inversion. @see invert
    const InverseTransform& operator~() const {return invert();}

    /// Overload transpose operator to mean inversion. @see updInvert
    InverseTransform&       operator~()       {return updInvert();}

    /// Compose the current transform (X_BF) with the given one. That is,
    /// return X_BY=X_BF*X_FY. Cost is 63 flops.
    Transform compose(const Transform& X_FY) const {
        return Transform(R_BF * X_FY.R(),
                         T_BF + R_BF * X_FY.T());
    }

    /// Compose the current transform (X_BF) with one that is supplied
    /// as an InverseTransform (typically as a result of applying
    /// the "~" operator to a transform). That is, return 
    /// X_BY=X_BF*X_FY, but now X_FY is represented as ~X_YF. Cost
    /// is an extra 18 flops to calculate X_FY.T(), total 81 flops.
    // (Definition is below after InverseTransform is declared.)
    inline Transform compose(const InverseTransform& X_FY) const;

    /// %Transform a vector expressed in our "F" frame to our "B" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec3 xformFrameVecToBase(const Vec3& vF) const {return R_BF*vF;}

    /// %Transform a vector expressed in our "B" frame to our "F" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec3 xformBaseVecToFrame(const Vec3& vB) const {return ~R_BF*vB;}

    /// %Transform a point (station) measured from and expressed in
    /// our "F" frame to that same point but measured from and
    /// expressed in our "B" frame. Cost is 18 flops.
    Vec3 shiftFrameStationToBase(const Vec3& sF) const {
        return T_BF + xformFrameVecToBase(sF);
    }

    /// %Transform a point (station) measured from and expressed in
    /// our "B" frame to that same point but measured from and
    /// expressed in our "F" frame. Cost is 18 flops.
    Vec3 shiftBaseStationToFrame(const Vec3& sB) const {
        return xformBaseVecToFrame(sB - T_BF);
    }

    /// Return a read-only reference to the contained rotation R_BF.
    const Rotation& R()    const { return R_BF; }

    /// Return a writable (lvalue) reference to the contained rotation R_BF.
    Rotation& updR()       { return R_BF; }

    /// Return a read-only reference to the x direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType& x() const {return R().x();}
    /// Return a read-only reference to the y direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType& y() const {return R().y();}
    /// Return a read-only reference to the z direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType& z() const {return R().z();}

    /// Return a read-only reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    const InverseRotation& RInv() const {return ~R_BF;}

    /// Return a writable (lvalue) reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    InverseRotation& updRInv() {return ~R_BF;}

    /// Return a read-only reference to our translation vector T_BF.
    const Vec3& T() const {return T_BF;}

    /// Return a writable (lvalue) reference to our translation vector T_BF.
    /// Caution: if you write through this reference you update the
    /// transform.
    Vec3& updT() {return T_BF;}

    /// Assign a new value to our translation vector. We expect the
    /// supplied vector @p T to be expressed in our B frame. A reference
    /// to the now-modified transform is returned as though this were
    /// an assignment operator.
    Transform& setT(const Vec3& T) {T_BF=T; return *this;}

    /// Calculate the inverse of the translation vector in this transform.
    /// The returned vector will be the negative of the original and will
    /// be expressed in the F frame rather than our B frame. Cost is 18 flops.
    Vec3 TInv() const {return -(~R_BF*T_BF);}

    /// Assign a value to the @em inverse of our translation vector.
    /// That is, we're given a vector in F which we invert and reexpress
    /// in B to store it in T, so that we get the original argument back if
    /// we ask for the inverse of T. Sorry, can't update TInv as an lvalue, but here we
    /// want -(~R_BF*T_BF)=T_FB => T_BF=-(R_BF*T_FB) so we can calculate
    /// it in 18 flops. A reference to the now-modified transform is returned
    /// as though this were an assignment operator.
    Transform& setTInv(const Vec3& T_FB) {
        T_BF = -(R_BF*T_FB); return *this;
    }

    /// Recast this transform as a read-only 3x4 matrix. This is just a
    /// reinterpretation of the data; the first three columns are the
    /// columns of the rotation and the last column is the translation.
    /// Zero cost.
    const Mat34& asMat34() const {
        return Mat34::getAs(reinterpret_cast<const Real*>(this));
    }

    /// Less efficient version of asMat34(); copies into return variable.
    Mat34 toMat34() const {
        return asMat34();
    }

    /// Return the equivalent 4x4 transformation matrix.
    Mat44 toMat44() const {
        Mat44 tmp;
        tmp.updSubMat<3,4>(0,0) = asMat34();
        tmp[3]                  = Row4(0,0,0,1);
        return tmp;
    }
private:
    //TODO: these might not pack correctly; should use an array of 12 Reals.
    Rotation R_BF;   // rotation matrix that expresses F's axes in R
    Vec3     T_BF;   // location of F's origin measured from B's origin, expressed in B 
};

/// If we multiply a transform by a 3-vector, we treat it as though it had
/// a 4th element "1" appended, that is, it is treated as a @em station rather
/// than a @em vector.
inline Vec3 operator*(const Transform& X_BF, const Vec3& s_F) {
    return X_BF.shiftFrameStationToBase(s_F);
}

/// If we multiply a transform by an augmented 4-vector, we use the 4th element to 
/// decide how to treat it. The 4th element must be 0 or 1. If 0 it is
/// treated as a vector only and the translation is ignored. If 1 it
/// is treated as a station and rotated & shifted.
inline Vec4 operator*(const Transform& X_BF, const Vec4& a_F) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec3& v_F = Vec3::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec4 out;
    if (a_F[3] == 0) {
        Vec3::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);
        out[3] = 0;
    } else {
        Vec3::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);
        out[3] = 1;
    }
    return out;
}

/**
 * %Transform from frame B to frame F, but with the internal representation inverted.
 * That is, we store R*,T* here but the transform this represents is
 * @verbatim
 *
 *  B F    [       |   ]
 *   X   = [   R   | T ]   where R=~(R*), T = - ~(R*)(T*).
 *         [.......|...]
 *         [ 0 0 0   1 ] 
 *
 * @endverbatim
 */
class InverseTransform {
public:
    InverseTransform() : R_FB(), T_FB(0) { }
    // default copy, assignment, destructor

    // Implicit conversion to Transform
    operator Transform() const {
        return Transform(R(), T());
    }

    // Assignment from Transform. This means that the inverse
    // transform we're assigning to must end up with the same meaning
    // as the inverse transform X has, so we'll need:
    //          T* == X.TInv()
    //          R* == X.RInv()
    // Cost: one frame conversion and a negation for TInv, 18 flops.
    InverseTransform& operator=(const Transform& X) {
        // Be careful to do this in the right order in case X and this
        // are the same object, i.e. ~X = X which is weird but has
        // the same meaning as X = ~X, i.e. invert X in place.
        T_FB = X.TInv(); // This might change X.T ...
        R_FB = X.RInv(); // ... but this doesn't depend on X.T.
        return *this;
    }

    // Inverting one of these just recasts it back to a Transform.
    const Transform& invert() const {
        return *reinterpret_cast<const Transform*>(this);
    }
    Transform& updInvert() {
        return *reinterpret_cast<Transform*>(this);
    }

    // Overload transpose to mean inversion.
    const Transform& operator~() const {return invert();}
    Transform&       operator~()       {return updInvert();}

    // Return X_BY=X_BF*X_FY, where X_BF (this) is represented here as ~X_FB. This
    // costs exactly the same as a composition of two Transforms (63 flops).
    Transform compose(const Transform& X_FY) const {
        return Transform(~R_FB * X_FY.R(),
                         ~R_FB *(X_FY.T() - T_FB));
    }
    // Return X_BY=X_BF*X_FY, but now both xforms are represented by their inverses.
    // This costs one extra vector transformation and a negation (18 flops) more
    // than a composition of two Transforms, for a total of 81 flops.
    Transform compose(const InverseTransform& X_FY) const {
        return Transform( ~R_FB * X_FY.R(),
                          ~R_FB *(X_FY.T() - T_FB));
    }

    // Forward and inverse vector transformations cost the same here as
    // for a Transform (or for that matter, a Rotation): 15 flops.
    Vec3 xformFrameVecToBase(const Vec3& vF) const {return ~R_FB*vF;}
    Vec3 xformBaseVecToFrame(const Vec3& vB) const {return  R_FB*vB;}

    // Forward and inverse station shift & transform cost the same here
    // as for a Transform: 18 flops.
    Vec3 shiftFrameStationToBase(const Vec3& sF) const {
        return ~R_FB*(sF-T_FB);
    }
    Vec3 shiftBaseStationToFrame(const Vec3& sB) const {
        return R_FB*sB + T_FB;
    }
    
    const InverseRotation& R()       const {return ~R_FB;}
    InverseRotation&       updR()          {return ~R_FB;}

    const InverseRotation::ColType& x() const {return R().x();}
    const InverseRotation::ColType& y() const {return R().y();}
    const InverseRotation::ColType& z() const {return R().z();}

    const Rotation&        RInv()    const {return R_FB;}
    Rotation&              updRInv()       {return R_FB;}

    // Costs 18 flops to look at the real translation vector.
    Vec3 T() const {return -(~R_FB*T_FB);}
    // no updT lvalue

    // Sorry, can't update translation as an lvalue, but here we
    // want -(R_BF*T_FB)=T_BF => T_FB=-(R_FB*T_BF). Cost: 18 flops.
    void setT(const Vec3& T_BF) {
        T_FB = -(R_FB*T_BF);
    }

    // Inverse translation is free.
    const Vec3& TInv() const           {return T_FB;}
    void        setTInv(const Vec3& T) {T_FB=T;}

    /// For compatibility with Transform, but we don't provide an "as"
    /// method here since the internal storage layout is somewhat odd.
    Mat34 toMat34() const {
        return Transform(*this).asMat34();
    }

    /// Return the equivalent 4x4 transformation matrix.
    Mat44 toMat44() const {
        return Transform(*this).toMat44();
    }

private:
    // DATA LAYOUT MUST BE IDENTICAL TO Transform !!
    // TODO: redo packing here when it is done for Transform.
    Rotation R_FB; // transpose of our rotation matrix, R_BF
    Vec3     T_FB; // our translation is -(R_BF*T_FB)=-(~R_FB*T_FB)
};

/// If we multiply a transform by a 3-vector, we treat it as though it had
/// a 4th element "1" appended, that is, it is treated as a *station* rather
/// than a *vector*.
inline Vec3 operator*(const InverseTransform& X_BF, const Vec3& s_F) {
    return X_BF.shiftFrameStationToBase(s_F);
}

/// If we multiply a transform by an augmented 4-vector, we use the 4th element to 
/// decide how to treat it. The 4th element must be 0 or 1. If 0 it is
/// treated as a vector only and the translation is ignored. If 1 it
/// is treated as a station and rotated & shifted.
inline Vec4 operator*(const InverseTransform& X_BF, const Vec4& a_F) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec3& v_F = Vec3::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec4 out;
    if (a_F[3] == 0) {
        Vec3::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);
        out[3] = 0;
    } else {
        Vec3::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);
        out[3] = 1;
    }
    return out;
}

// These Transform definitions had to wait for InverseTransform to be declared.

inline Transform& 
Transform::operator=(const InverseTransform& X) {
    // Be careful to do this in the right order in case X and this
    // are the same object, i.e. we're doing X = ~X, inverting X in place.
    T_BF = X.T(); // This might change X.T ...
    R_BF = X.R(); // ... but this doesn't depend on X.T.
    return *this;
}

inline Transform 
Transform::compose(const InverseTransform& X_FY) const {
    return Transform(R_BF * X_FY.R(),
                        T_BF + R_BF * X_FY.T());
}

inline Transform
operator*(const Transform& X1, const Transform& X2) {
    return X1.compose(X2);
}
inline Transform
operator*(const Transform& X1, const InverseTransform& X2) {
    return X1.compose(X2);
}
inline Transform
operator*(const InverseTransform& X1, const Transform& X2) {
    return X1.compose(X2);
}
inline Transform
operator*(const InverseTransform& X1, const InverseTransform& X2) {
    return X1.compose(X2);
}
inline bool
operator==(const Transform& X1, const Transform& X2) {
    return X1.R()==X2.R() && X1.T()==X2.T();
}
inline bool
operator==(const InverseTransform& X1, const InverseTransform& X2) {
    return X1.R()==X2.R() && X1.T()==X2.T();
}
inline bool
operator==(const Transform& X1, const InverseTransform& X2) {
    return X1.R()==X2.R() && X1.T()==X2.T();
}
inline bool
operator==(const InverseTransform& X1, const Transform& X2) {
    return X1.R()==X2.R() && X1.T()==X2.T();
}

SimTK_SimTKCOMMON_EXPORT std::ostream& 
operator<<(std::ostream& o, const Transform&);

} // namespace SimTK


#endif // SimTK_SIMMATRIX_ORIENTATION_H_
