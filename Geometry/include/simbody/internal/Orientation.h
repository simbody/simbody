#ifndef SIMTK_SIMBODY_ORIENTATION_H_
#define SIMTK_SIMBODY_ORIENTATION_H_

/** @file
 *
 * These are numerical utility classes for dealing with the relative orientations
 * of geometric objects. These build on the basic arithmetic classes for small
 * vectors and matrices.
 */

#include "simtk/SimTK.h"
#include "simmatrix/SmallMatrix.h"

#include <iostream>

namespace simtk {

template <int S> class UnitVec;
template <int S> class UnitRow;
class RotationMat;
class InverseRotationMat;
class TransformMat;
class InverseTransformMat;

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

    UnitVec() : BaseVec(NTraits<Real>::getNaN()) { }

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
        static_cast<BaseVec&>(*this) /= norm();
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
        BaseVec::operator=(static_cast<const UnitVec<S2>::BaseVec&>(u));
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
 * type of a RotationMat. Don't construct these directly.
 */
template <int S>
class UnitRow : public Row<3,Real,S> {
public:
    typedef Row<3,Real,S> BaseRow;
    typedef UnitVec<S>    TransposeType;

    UnitRow() : BaseRow(NTraits<Real>::getNaN()) { }

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
        BaseRow::operator=(static_cast<const UnitRow<S2>::BaseRow&>(u));
        return *this;
    }

    // Explicit conversion from Row to UnitRow, requiring expensive normalization.
    explicit UnitRow(const BaseRow& v) : BaseRow(v/v.norm()) { }
    template <int S2> explicit UnitRow(const Row<3,Real,S2>& v)
        : BaseRow(v/v.norm()) { }

    UnitRow(const Real& x, const Real& y, const Real& z) : BaseRow(x,y,z) {
        static_cast<BaseRow&>(*this) /= norm();
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
class RotationMat : public Mat33 {
public:
    typedef Mat33 BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    RotationMat() : BaseMat(1.) { }    // default is identity

    /// Create a Rotation matrix by specifying only its z axis. 
    /// The resulting x and y axes will be appropriately perpendicular
    /// but are otherwise arbitrary. This will work for any stride
    /// UnitVec because there is always an implicit conversion
    /// available to the packed form used as the argument.
    explicit RotationMat(const UnitVec3& z);

    RotationMat(const RotationMat& R) : BaseMat(R) { }
    RotationMat& operator=(const RotationMat& R) {
        BaseMat::operator=(R.asMat33()); return *this;
    }

    /// By zero we mean "zero rotation", i.e., an identity matrix.
    void setToZero() {
        *static_cast<Mat33*>(this) = 1.;
    }

    void setToNaN() {
        static_cast<Mat33*>(this)->setToNaN();
    }

    /// Set this RotationMat to represent a rotation of +q radians
    /// around the base frame's 0,0,1 axis.
    void setToRotationAboutZ(const Real& q) {
        const Real sq = std::sin(q), cq = std::cos(q);
        *static_cast<Mat33*>(this) = Mat33( cq , -sq , 0. ,
                                            sq ,  cq , 0. ,
                                            0. ,  0. , 1. );
    }

    /// Set this RotationMat to represent a rotation of +q0 about
    /// the base frame's X axis, followed by a rotation of +q1 about
    /// the base frame's (unchanged) Y axis.
    void setToSpaceFixed12(const Vec2& q) {
        const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
        const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
        *static_cast<Mat33*>(this) = Mat33( cq1 , sq1*sq0 , sq1*cq0,
                                             0. ,   cq0   ,   -sq0 ,
                                           -sq1 , cq1*sq0 , cq1*cq0);
    }

    /// Set this RotationMat to represent a rotation of +q0 about
    /// the body frame's Z axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW X axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 3-2-1.
    void setToBodyFixed321(const Vec3& q) {
        const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
        const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
        const Real sq2 = std::sin(q[2]), cq2 = std::cos(q[2]);
        *static_cast<Mat33*>(this) = 
            Mat33( cq0*cq1 , cq0*sq1*sq2-sq0*cq2 , cq0*sq1*cq2+sq0*sq2,
                   sq0*cq1 , sq0*sq1*sq2+cq0*cq2 , sq0*sq1*cq2-cq0*sq2,
                    -sq1   ,       cq1*sq2       ,      cq1*cq2        );
    }

    /// Set this RotationMat to represent a rotation of +q0 about
    /// the body frame's X axis, followed by a rotation of +q1 about
    /// the body frame's NEW Y axis, followed by a rotation of +q3
    /// about the body frame's NEW Z axis.
    /// See Kane, Spacecraft Dynamics, pg. 423, body-three: 1-2-3.
    void setToBodyFixed123(const Vec3& q) {
        const Real sq0 = std::sin(q[0]), cq0 = std::cos(q[0]);
        const Real sq1 = std::sin(q[1]), cq1 = std::cos(q[1]);
        const Real sq2 = std::sin(q[2]), cq2 = std::cos(q[2]);
        *static_cast<Mat33*>(this) = 
            Mat33(      cq1*cq2        ,       -cq1*sq2       ,  sq1    ,
                   sq0*sq1*cq2+cq0*sq2 , -sq0*sq1*sq2+cq0*cq2 , -sq0*cq1,
                   -cq0*sq1*cq2+sq0*sq2 , cq0*sq1*sq2+sq0*cq2 ,  cq0*cq1 );
    }

    /// Set this RotationMat to represent the same rotation as
    /// the passed-in quaternion. The 0th element is the quaternion
    /// scalar. The quaternion is normalized before use to ensure
    /// a non-distorting rotation matrix.
    void setToQuaternion(const Vec4& qin) {
        const Vec4 q = qin / qin.norm();
        const Real q00=q[0]*q[0], q11=q[1]*q[1], q22=q[2]*q[2], q33=q[3]*q[3];
        const Real q01=q[0]*q[1], q02=q[0]*q[2], q03=q[0]*q[3];
        const Real q12=q[1]*q[2], q13=q[1]*q[3], q23=q[2]*q[3];

        *static_cast<Mat33*>(this) = 
            Mat33(q00+q11-q22-q33,   2.*(q12-q03)  ,   2.*(q13+q02),
                    2.*(q12+q03)  , q00-q11+q22-q33,   2.*(q23-q01),
                    2.*(q13-q02)  ,   2.*(q23+q01)  , q00-q11-q22+q33);
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


    // Inverse of the above routine. Returned angular velocity is B in P,
    // expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed321DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1 = Real(1)/c1;
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
        const Real ooc1 = Real(1)/c1;
        const Real s2oc1 = s2*ooc1, c2oc1 = c2*ooc1;

        const Mat33 E(    c2oc1  , -s2oc1  , 0.,
                            s2   ,    c2   , 0.,
                       -s1*c2oc1 , s1*s2oc1, 1. );
        return E*w;
    }

    // Inverse of the above routine. Returned angular velocity is B in P,
    // expressed in *B*: w_PB_B.
    static Vec3 convertBodyFixed123DotToAngVel(const Vec3& q, const Vec3& qd) {
        const Real s1 = std::sin(q[1]), c1 = std::cos(q[1]);
        const Real s2 = std::sin(q[2]), c2 = std::cos(q[2]);
        const Real ooc1 = Real(1)/c1;
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
        const Real ooc1 = Real(1)/c1;
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
        const Real ooc1 = Real(1)/c1;
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

    const InverseRotationMat& invert() const {
        return *reinterpret_cast<const InverseRotationMat*>(this);
    }
    InverseRotationMat& updInvert() {
        return *reinterpret_cast<InverseRotationMat*>(this);
    }

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

    const InverseRotationMat& operator~() const {return invert();}
    InverseRotationMat&       operator~()       {return updInvert();}

    const RowType& operator[](int i) const {return row(i);}
    const ColType& operator()(int j) const {return col(j);}

    const BaseMat& asMat33() const {
        return *static_cast<const BaseMat*>(this);
    }

    static RotationMat trustMe(const Mat33& m) {return RotationMat(m);}

private:
    // We're trusting that m is a rotation.
    explicit RotationMat(const BaseMat& m) : BaseMat(m) { }
    template <int CS, int RS>
    explicit RotationMat(const Mat<3,3,Real,CS,RS>& m) : BaseMat(m) { }

    friend RotationMat operator*(const RotationMat&,const RotationMat&);
};

class InverseRotationMat : public Mat33::TransposeType {
public:
    typedef Mat33::TransposeType BaseMat;
    typedef UnitVec<BaseMat::RowSpacing> ColType;
    typedef UnitRow<BaseMat::ColSpacing> RowType;

    // Don't construct one these; they should only occur as expression intermediates.
    // But if you must ...
    InverseRotationMat() : BaseMat(1) { }

    InverseRotationMat(const InverseRotationMat& R) : BaseMat(R) { }
    InverseRotationMat& operator=(const InverseRotationMat& R) {
        BaseMat::operator=(R.asMat33()); return *this;
    }

    // Implicit conversion to RotationMat.
    operator RotationMat() const {
        return RotationMat::trustMe(asMat33());
    }

    const RotationMat& invert() const {
        return *reinterpret_cast<const RotationMat*>(this);
    }
    RotationMat& updInvert() {
        return *reinterpret_cast<RotationMat*>(this);
    }

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

    const RotationMat& operator~() const {return invert();}
    RotationMat&       operator~()       {return updInvert();}

    const RowType& operator[](int i) const {return row(i);}
    const ColType& operator()(int j) const {return col(j);}

    const BaseMat& asMat33() const {
        return *static_cast<const BaseMat*>(this);
    }
};


std::ostream& operator<<(std::ostream& o, const RotationMat& m);

template <int S> inline UnitVec<1>
operator*(const RotationMat& R, const UnitVec<S>& v) {
    return UnitVec<1>(R.asMat33()*v.asVec3(), true);
}
template <int S> inline UnitRow<1>
operator*(const UnitRow<S>& r, const RotationMat& R) {
    return UnitRow<1>(r.asRow3(), R.asMat33(), true);
}

template <int S> inline UnitVec<1>
operator*(const InverseRotationMat& R, const UnitVec<S>& v) {
    return UnitVec<1>(R.asMat33()*v.asVec3(), true);
}
template <int S> inline UnitRow<1>
operator*(const UnitRow<S>& r, const InverseRotationMat& R) {
    return UnitRow<1>(r.asRow3(), R.asMat33(), true);
}

inline RotationMat
operator*(const RotationMat& R1, const RotationMat& R2) {
    return RotationMat::trustMe(R1.asMat33()*R2.asMat33());
}

inline RotationMat
operator*(const RotationMat& R1, const InverseRotationMat& R2) {
    return RotationMat::trustMe(R1.asMat33()*R2.asMat33());
}

inline RotationMat
operator*(const InverseRotationMat& R1, const RotationMat& R2) {
    return RotationMat::trustMe(R1.asMat33()*R2.asMat33());
}

inline RotationMat
operator*(const InverseRotationMat& R1, const InverseRotationMat& R2) {
    return RotationMat::trustMe(R1.asMat33()*R2.asMat33());
}

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
 * The axis vectors constitute a RotationMat. They are ordered 1-2-3 or x-y-z
 & as you prefer, with z = x X y, making a right-handed set. These axes are arranged
 * as columns of a 3x3 rotation matrix R_BF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame B and F. (The
 * columns of R_BF are F's coordinate axes, expressed in B.) For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in B is given by vB = R_BF*vF. F's origin point OF is 
 * stored as the translation vector T_BF=(OF-OB) and expressed in B.
 *
 * TransformMat is designed to behave as much as possible like the computer
 * graphics 4x4 transform X which would be arranged like this:
 *
 *         [       |   ]
 *     X = [   R   | T ]    R is a 3x3 orthogonal rotation matrix
 *         [.......|...]    T os a 3x1 translation vector
 *         [ 0 0 0   1 ]
 *
 * These can be composed directly by matrix multiplication, but more 
 * importantly they have a particularly simple inverse:
 *
 *    -1   [       |    ]
 *   X   = [  ~R   | T* ]   ~R is R transpose, T* = ~R(-T).
 *         [.......|....]
 *         [ 0 0 0   1  ] 
 *
 * This inverse is so simple that we compute it simply by defining another
 * type, InverseTransformMat, which is identical to TransformMat in memory but
 * behaves as though it contains the inverse. That way we invert just by
 * changing point of view (recasting) rather than computing.
 *
 * This is a "POD" (plain old data) class with a well-defined memory
 * layout on which a client of this class may depend: There are 
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,T.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,T1,T2,T3. It is expressly allowed
 * to reinterpret TransformMat objects in any appropriate manner that depends
 * on this memory layout.
 */
class TransformMat {
public:
    TransformMat()                                    : R_BF(),  T_BF(0) { }
    TransformMat(const RotationMat& R, const Vec3& T) : R_BF(R), T_BF(T) { }
    explicit TransformMat(const RotationMat& R)       : R_BF(R), T_BF(0) { }
    explicit TransformMat(const Vec3& T)              : R_BF(),  T_BF(T) { }
    // default copy, assignment, destructor

    // Assignment from InverseTransformMat. This means that the 
    // transform we're assigning to must end up with the same meaning
    // as the inverse transform X has, so we'll need:
    //          T == X.T()
    //          R == X.R()
    // Cost: one frame conversion and a negation, 18 flops.
    // Definition is below after InverseTransformMat is declared.
    inline TransformMat& operator=(const InverseTransformMat& X);

    void set(const RotationMat& R, const Vec3& T) {T_BF=T; R_BF=R;}

    /// By zero we mean "zero transform", i.e., an identity rotation and zero translation.
    void setToZero() {
        R_BF.setToZero(); T_BF = 0.;
    }

    void setToNaN() {
        R_BF.setToNaN(); T_BF.setToNaN();
    }

    // Inverting one of these just casts it to a TransformInverseMat.
    const InverseTransformMat& invert() const {
        return *reinterpret_cast<const InverseTransformMat*>(this);
    }
    InverseTransformMat& updInvert() {
        return *reinterpret_cast<InverseTransformMat*>(this);
    }

    // Overload transpose to mean inversion.
    const InverseTransformMat& operator~() const {return invert();}
    InverseTransformMat&       operator~()       {return updInvert();}

    // Return X_BY=X_BF*X_FY. Cost is 63 flops.
    TransformMat compose(const TransformMat& X_FY) const {
        return TransformMat(R_BF * X_FY.R(),
                            T_BF + R_BF * X_FY.T());
    }

    // Return X_BY=X_BF*X_FY, but now X_FY is represented as ~X_YF. Cost
    // is an extra 18 flops to calculate X_FY.T(), total 81 flops.
    // Definition is below after InverseTransformMat is declared.
    inline TransformMat compose(const InverseTransformMat& X_FY) const;

    // Costs 15 flops to transform vectors in either direction.
    Vec3 xformFrameVecToBase(const Vec3& vF) const {return R_BF*vF;}
    Vec3 xformBaseVecToFrame(const Vec3& vR) const {return ~R_BF*vR;}

    // Costs 18 flops to transform & shift stations either direction.
    Vec3 shiftFrameStationToBase(const Vec3& sF) const {
        return T_BF + xformFrameVecToBase(sF);
    }
    Vec3 shiftBaseStationToFrame(const Vec3& sB) const {
        return xformBaseVecToFrame(sB - T_BF);
    }

    const RotationMat& R()    const { return R_BF; }
    RotationMat&       updR()       { return R_BF; }

    const RotationMat::ColType& x() const {return R().x();}
    const RotationMat::ColType& y() const {return R().y();}
    const RotationMat::ColType& z() const {return R().z();}

    const InverseRotationMat& RInv()    const { return ~R_BF; }
    InverseRotationMat&       updRInv()       { return ~R_BF; }

    const Vec3&  T()    const        {return T_BF;}
    Vec3&        updT()              {return T_BF;}
    void         setT(const Vec3& T) {T_BF=T;}

    // Costs 18 flops to calculate the inverse translation.
    Vec3 TInv() const { return -(~R_BF*T_BF); }

    // Sorry, can't update TInv as an lvalue, but here we
    // want -(~R_BF*T_BF)=T_FB => T_BF=-(R_BF*T_FB). Cost: 18 flops.
    void setTInv(const Vec3& T_FB) {
        T_BF = -(R_BF*T_FB);
    }

    const Mat34& asMat34() const {return Mat34::getAs(reinterpret_cast<const Real*>(this));}
private:
    RotationMat R_BF;   // rotation matrix that expresses F's axes in R
    Vec3        T_BF;   // location of F's origin measured from B's origin, expressed in B 
};

/*
 * Transform from frame B to frame F, but with the internal representation inverted.
 * That is, we store R*,T* here but the transform this represents is
 *
 *  B F    [       |   ]
 *   X   = [   R   | T ]   where R=~(R*), T = - ~(R*)(T*).
 *         [.......|...]
 *         [ 0 0 0   1 ] 
 */
class InverseTransformMat {
public:
    InverseTransformMat() : R_FB(), T_FB(0) { }
    // default copy, assignment, destructor

    // Implicit conversion to TransformMat
    operator TransformMat() const {
        return TransformMat(R(), T());
    }

    // Assignment from TransformMat. This means that the inverse
    // transform we're assigning to must end up with the same meaning
    // as the inverse transform X has, so we'll need:
    //          T* == X.TInv()
    //          R* == X.RInv()
    // Cost: one frame conversion and a negation for TInv, 18 flops.
    InverseTransformMat& operator=(const TransformMat& X) {
        // Be careful to do this in the right order in case X and this
        // are the same object, i.e. ~X = X which is weird but has
        // the same meaning as X = ~X, i.e. invert X in place.
        T_FB = X.TInv(); // This might change X.T ...
        R_FB = X.RInv(); // ... but this doesn't depend on X.T.
    }

    // Inverting one of these just recasts it back to a TransformMat.
    const TransformMat& invert() const {
        return *reinterpret_cast<const TransformMat*>(this);
    }
    TransformMat& updInvert() {
        return *reinterpret_cast<TransformMat*>(this);
    }

    // Overload transpose to mean inversion.
    const TransformMat& operator~() const {return invert();}
    TransformMat&       operator~()       {return updInvert();}

    // Return X_BY=X_BF*X_FY, where X_BF (this) is represented here as ~X_FB. This
    // costs exactly the same as a composition of two TransformMats (63 flops).
    TransformMat compose(const TransformMat& X_FY) const {
        return TransformMat(~R_FB * X_FY.R(),
                            ~R_FB *(X_FY.T() - T_FB));
    }
    // Return X_BY=X_BF*X_FY, but now both xforms are represented by their inverses.
    // This costs one extra vector transformation and a negation (18 flops) more
    // than a composition of two TransformMats, for a total of 81 flops.
    TransformMat compose(const InverseTransformMat& X_FY) const {
        return TransformMat( ~R_FB * X_FY.R(),
                             ~R_FB *(X_FY.T() - T_FB));
    }

    // Forward and inverse vector transformations cost the same here as
    // for a TransformMat (or for that matter, a RotationMat): 15 flops.
    Vec3 xformFrameVecToBase(const Vec3& vF) const {return ~R_FB*vF;}
    Vec3 xformBaseVecToFrame(const Vec3& vB) const {return  R_FB*vB;}

    // Forward and inverse station shift & transform cost the same here
    // as for a TransformMat: 18 flops.
    Vec3 shiftFrameStationToBase(const Vec3& sF) const {
        return ~R_FB*(sF-T_FB);
    }
    Vec3 shiftBaseStationToFrame(const Vec3& sB) const {
        return R_FB*sB + T_FB;
    }
    
    const InverseRotationMat& R()       const {return ~R_FB;}
    InverseRotationMat&       updR()          {return ~R_FB;}

    const InverseRotationMat::ColType& x() const {return R().x();}
    const InverseRotationMat::ColType& y() const {return R().y();}
    const InverseRotationMat::ColType& z() const {return R().z();}

    const RotationMat&        RInv()    const {return R_FB;}
    RotationMat&              updRInv()       {return R_FB;}

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

private:
    // DATA LAYOUT MUST BE IDENTICAL TO TransformMat !!
    RotationMat R_FB; // transpose of our rotation matrix, R_BF
    Vec3        T_FB; // our translation is -(R_BF*T_FB)=-(~R_FB*T_FB)
};

// These had to wait for InverseTransformMat to be declared.

inline TransformMat& 
TransformMat::operator=(const InverseTransformMat& X) {
    // Be careful to do this in the right order in case X and this
    // are the same object, i.e. we're doing X = ~X, inverting X in place.
    T_BF = X.T(); // This might change X.T ...
    R_BF = X.R(); // ... but this doesn't depend on X.T.
}

inline TransformMat 
TransformMat::compose(const InverseTransformMat& X_FY) const {
    return TransformMat(R_BF * X_FY.R(),
                        T_BF + R_BF * X_FY.T());
}

inline TransformMat
operator*(const TransformMat& X1, const TransformMat& X2) {
    return X1.compose(X2);
}
inline TransformMat
operator*(const TransformMat& X1, const InverseTransformMat& X2) {
    return X1.compose(X2);
}
inline TransformMat
operator*(const InverseTransformMat& X1, const TransformMat& X2) {
    return X1.compose(X2);
}
inline TransformMat
operator*(const InverseTransformMat& X1, const InverseTransformMat& X2) {
    return X1.compose(X2);
}

std::ostream& operator<<(std::ostream& o, const TransformMat&);

} // namespace simtk

#endif // SIMTK_SIMBODY_ORIENTATION_H_
