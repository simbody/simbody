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

class UnitVec3;
class UnitRow3;
class MatRotation;
class Frame;

/**
 * This class is a Vec3 plus an ironclad guarantee either that:
 *      - the length is one (to within a very small tolerance), or
 *      - all components are NaN.
 * Thus it is a pure direction.
 */
class UnitVec3 : public Vec3 {
public:
    UnitVec3() : Vec3(NTraits<Real>::getNaN()) { }
    UnitVec3(const UnitVec3& u) : Vec3(static_cast<const Vec3&>(u)) { }
    explicit UnitVec3(const Vec3& v) : Vec3(v/v.norm()) { }

    UnitVec3(const Real& x, const Real& y, const Real& z) : Vec3(x,y,z) {
        static_cast<Vec3&>(*this) /= norm();
    }

    // Create a unit axis vector 100 010 001
    explicit UnitVec3(int axis) : Vec3(0.) {
        assert(0 <= axis && axis <= 2);
        (*this)[axis] = 1.;
    }

    UnitVec3& operator=(const UnitVec3& u) {
        Vec3::operator=(static_cast<const Vec3&>(u)); 
        return *this;
    }

    const Vec3& asVec3() const {return static_cast<const Vec3&>(*this);}

    // Override Vec3 methods which preserve length.
    UnitVec3 negate()    const {return UnitVec3(-asVec3(),true);}
    UnitVec3 operator-() const {return negate();}

    // We have to define these here since we had to override the writable
    // versions to make them private.
    const Real& operator[](int i) const {return Vec3::operator[](i);}
    const Real& operator()(int i) const {return Vec3::operator()(i);}

    // Return a vector whose measure numbers are the absolute values
    // of the ones here. This will still have unit length but will be
    // a reflection of this unit vector into the first octant (+x,+y,+z).
    UnitVec3 abs() const {
        return UnitVec3(asVec3().abs(),true);
    }

    // Return a unit vector perpendicular to this one (arbitrary).
    inline UnitVec3 perp() const;

private:
    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitVec3(const Vec3& v, bool) : Vec3(v) { }

    // These must be overridden here to "privatize" them. Only the elite
    // few can be trusted to update a single measure number without
    // blowing the unit vector guarantee.
    Real& operator[](int i) {return Vec3::operator[](i);}
    Real& operator()(int i) {return Vec3::operator()(i);}

    friend UnitVec3 operator*(const MatRotation&, const UnitVec3&);
};
std::ostream& operator<<(std::ostream& o, const UnitVec3& v);

// Scalar multiply and divide don't preserve 'unitness'
inline Vec3  operator*(const UnitVec3& v, const Real& r) {return v.asVec3()*r;}
inline Vec3  operator*(const Real& r, const UnitVec3& v) {return v.asVec3()*r;}
inline Vec3  operator/(const UnitVec3& v, const Real& r) {return v.asVec3()/r;}

inline Real  operator*(const Row3&     r, const UnitVec3& u) {return r*u.asVec3();}
inline Mat33 operator*(const UnitVec3& u, const Row3&     r) {return u.asVec3()*r;}
inline Vec3  operator%(const UnitVec3& u, const UnitVec3& v) {return u.asVec3()%v.asVec3();}
inline Vec3  operator%(const Vec3&     v, const UnitVec3& u) {return v%u.asVec3();}
inline Vec3  operator%(const UnitVec3& u, const Vec3&     v) {return u.asVec3()%v;}
inline Row3  operator%(const Row3&     r, const UnitVec3& u) {return r%u.asVec3();}
inline Row3  operator%(const UnitVec3& u, const Row3&     r) {return u.asVec3()%r;}

inline UnitVec3 UnitVec3::perp() const {
    // Choose the coordinate axis which makes the largest angle
    // with this vector, that is, has the "least u" along it.
    const UnitVec3 u(abs());    // reflect to first octant
    const int minAxis = u[0] <= u[1] ? (u[0] <= u[2] ? 0 : 2)
                                     : (u[1] <= u[2] ? 1 : 2);
    // Cross returns a Vec3 result (see operator above), which is then normalized.
    return UnitVec3(*this % UnitVec3(minAxis));
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
class MatRotation {
    Mat33 rot_GF; // xyz axes of frame F expressed in frame G
public:
    MatRotation() : rot_GF(1.) { }    // default is identity

    /// Create a Rotation matrix by specifying only its z axis. 
    /// The resulting x and y axes will be appropriately perpendicular
    /// but are otherwise arbitrary.
    explicit MatRotation(const UnitVec3& z);

    const UnitVec3& getAxis(int i) const
      { return reinterpret_cast<const UnitVec3&>(rot_GF(i)); }

    // TODO: with much agony involving templates this could be made free.
    MatRotation operator~() const {return MatRotation(~rot_GF);}

    const UnitVec3& operator()(int i) const {
        return reinterpret_cast<const UnitVec3&>(rot_GF(i));
    }

    const Mat33& asMat33() const {return rot_GF;}

private:
    // We're trusting that m is a rotation.
    explicit MatRotation(const Mat33& m) : rot_GF(m) { }
    friend MatRotation operator*(const MatRotation&,const MatRotation&);
};
std::ostream& operator<<(std::ostream& o, const MatRotation& m);

inline MatRotation operator*(const MatRotation& l, const MatRotation& r) {
    return MatRotation(l.asMat33()*r.asMat33());
}
inline Mat33 operator*(const MatRotation& l, const Mat33& r) {
    return l.asMat33()*r;
}
inline Mat33 operator*(const Mat33& l, const MatRotation& r) {
    return l*r.asMat33();
}
inline UnitVec3 operator*(const MatRotation& R, const UnitVec3& v) {
    return UnitVec3(R.asMat33()*v.asVec3(), true);
}
inline Vec3 operator*(const MatRotation& R, const Vec3& v) {
    return R.asMat33()*v;
}
inline Row3 operator*(const Row3& r, const MatRotation& R) {
    return r*R.asMat33();
}

/**
 * This class represents an orthogonal, right-handed coordinate frame F, 
 * measured from and expressed in a base (reference) coordinate frame B.
 * F consists of 3 perpendicular unit vectors defining its axes XF as
 * viewed from B (that is, as expressed in B's axes XB), and a vector
 * from B's origin point OB to F's origin point OF. Note that
 * the meaning of "B" comes from the context in which frame F is used.
 * We use the phrase "frame F is in frame B" to describe the above relationship,
 * that is, "in" means both measured from and expressed in. 
 *
 * The axis vectors are ordered 1-2-3 or x-y-z as you prefer, with
 * z = x X y, making a right-handed set. These axes are arranged as
 * columns of a 3x3 matrix X_BF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame B and F. (The
 * columns of X_BF are F's coordinate axes, expressed in B.) For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in B is given by vB = X_BF*vF. F's origin point OF is 
 * stored as the vector OF_B=(OF-OB) and expressed in B.
 *
 * This is a "POD" (plain old data) class with a well-defined memory
 * layout on which a client of this class may depend: There are 
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,O.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,O1,O2,O3. It is expressly allowed
 * to reinterpret Frame objects in any appropriate manner that depends
 * on this memory layout.
 */
class Frame {
public:
    Frame() : X_BF(), OF_B(0.) { }
    Frame(const MatRotation& axesInB, const Vec3& orgInB) : X_BF(axesInB), OF_B(orgInB) { }
    explicit Frame(const MatRotation& axesInB) : X_BF(axesInB), OF_B(0.) { }
    explicit Frame(const Vec3& orgInB)         : X_BF(), OF_B(orgInB) { }
    // default copy, assignment, destructor

    void setFrame(const MatRotation& axesInB, const Vec3& orgInB) 
      { X_BF=axesInB; OF_B=orgInB; }

    // If this is frame F measured and expressed in B, return frame B
    // measured and expressed in F. This is the inverse of F in that
    // it maps from F to B rather than from B to F.
    Frame invert() const {
        const MatRotation rot_FR = ~X_BF;
        return Frame(rot_FR, rot_FR*(-OF_B));
    }
    // return frame_RX
    Frame compose(const Frame& frame_FX) const {
        const MatRotation rot_RX = X_BF * frame_FX.getAxes();
        return Frame(rot_RX, OF_B + X_BF * frame_FX.getOrigin());
    }
    Vec3 xformFrameVecToBase(const Vec3& vF) const {return X_BF*vF;}
    Vec3 xformBaseVecToFrame(const Vec3& vR) const {return ~X_BF*vR;}
    Vec3 shiftFrameStationToBase(const Vec3& sF) const {
        return OF_B + xformFrameVecToBase(sF);
    }
    Vec3 shiftBaseStationToFrame(const Vec3& sB) const {
        return xformBaseVecToFrame(sB - OF_B);
    }

    const MatRotation& getAxes() const { return X_BF; }
    MatRotation&       updAxes()       { return X_BF; }

    const Vec3&  getOrigin() const { return OF_B; }
    Vec3&        updOrigin()       { return OF_B; }

private:
    MatRotation X_BF;   // rotation matrix that expresses F's axes in R
    Vec3        OF_B;   // location of F's origin measured from B's origin, expressed in B 
};
std::ostream& operator<<(std::ostream& o, const Frame&);

} // namespace simtk

#endif // SIMTK_SIMBODY_ORIENTATION_H_
