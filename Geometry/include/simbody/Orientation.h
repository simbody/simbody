#ifndef ORIENTATION_H_
#define ORIENTATION_H_

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
class MatRotation;
class Frame;

/**
 * This class is a Vec3 plus an ironclad guarantee either that:
 *      - the length is one (to within a very small tolerance), or
 *      - all components are NaN.
 * Thus it is a pure direction.
 */
class UnitVec3 {
    Vec3 dir;   // a unit vector
public:
    UnitVec3() : dir(NTraits<Real>::getNaN()) { }
    explicit UnitVec3(const Vec3& v) : dir(v/v.norm()) { }

    const Vec3& asVec3() const {return dir;}

    // Implicit conversion to read-only Vec3 when appropriate.
    // Note that there is no such conversion to non-const Vec3.
    operator const Vec3&() const {return dir;}

    UnitVec3 negate()    const {return UnitVec3(-dir,true);}
    UnitVec3 operator-() const {return negate();}

    const Real& operator[](int i) const {return dir[i];}
private:
    // This constructor is only for our friends whom we trust to
    // give us an already-normalized vector.
    UnitVec3(const Vec3& v, bool) : dir(v) { }
    friend UnitVec3 operator*(const MatRotation&, const UnitVec3&);
};
std::ostream& operator<<(std::ostream& o, const UnitVec3& v);

// Scalar multiply and divide don't preserve 'unitness'
Vec3 operator*(const UnitVec3& v, const Real& r) {return v.asVec3()*r;}
Vec3 operator*(const Real& r, const UnitVec3& v) {return v.asVec3()*r;}
Vec3 operator/(const UnitVec3& v, const Real& r) {return v.asVec3()/r;}

/**
 * This class is a Mat33 plus an ironclad guarantee that the matrix represents
 * a pure rotation. A rotation is an orthogonal matrix whose columns and rows
 * are directions (that is, unit vectors) which are mutually
 * orthogonal. Furthermore, if the columns (or rows) are
 * labeled x,y,z it always holds that z = x X y ensuring that
 * this is a rotation and not a reflection.
 *
 * A rotation matrix is formed from a Cartesian coordinate frame
 * simply by using the x,y,z axes of the coordinate frame as
 * the columns of the rotation matrix. That is, if you have 
 * a frame F with its axes expressed in frame G, you can write
 * this as a matrix  R_GF = [ x y z ] which when applied to
 * a vector with measure numbers expressed in F yields that same
 * vector re-expressed in G: v_G = R_GF * v_F. Because a rotation
 * is orthogonal its transpose is its inverse, which is
 * also a rotation. So we write R_FG = transpose(R_GF), and this
 * matrix can be used to rotate in the other direction: 
 * v_F = R_FG * v_G = ~R_GF*v_G.
 */
class MatRotation {
    Mat33 R_GF; // xyz axes of frame F expressed in frame G
public:
    MatRotation() : R_GF(1.) { }

    /// Create a Rotation matrix by specifying only its z axis. 
    /// The resulting x and y axes will be appropriately perpendicular
    /// but are otherwise arbitrary.
    explicit MatRotation(const UnitVec3& z);

    const UnitVec3& getAxis(int i)
      { return reinterpret_cast<const UnitVec3&>(R_GF(i)); }

    // TODO: with much agony involving templates this could be made free.
    MatRotation operator~() const {return MatRotation(~R_GF);}

    const UnitVec3& operator()(int i) const {
        return reinterpret_cast<const UnitVec3&>(R_GF(i));
    }

    const Mat33& asMat33() const {return R_GF;}

    // Implicit conversion to read-only Vec3 when appropriate.
    // Note that there is no such conversion to non-const Vec3.
    operator const Mat33&() const {return R_GF;}
private:
    // We're trusting that m is a rotation.
    explicit MatRotation(const Mat33& m) : R_GF(m) { }
    friend MatRotation operator*(const MatRotation&,const MatRotation&);
};
std::ostream& operator<<(std::ostream& o, const MatRotation& m);

MatRotation operator*(const MatRotation& l, const MatRotation& r) {
    return MatRotation(l.asMat33()*r.asMat33());
}
UnitVec3 operator*(const MatRotation& R, const UnitVec3& v) {
    return UnitVec3(R.asMat33()*v.asVec3(), true);
}
Vec3 operator*(const MatRotation& R, const Vec3& v) {
    return R.asMat33()*v;
}
/**
 * This class represents an orthogonal, right-handed coordinate frame F, 
 * measured from and expressed in a reference coordinate frame R. F consists of
 * 3 perpendicular unit vectors defining its axes as viewed from R, 
 * and a vector from R's origin point OR to F's origin point OF. Note that
 * the meaning of "R" comes from the context in which frame F is used.
 * We use the phrase "frame F is in frame R" to describe the above relationship,
 * that is, "in" means both measured in and expressed in. 
 *
 * The axis vectors are ordered 1-2-3 or x-y-z as you prefer, with
 * z = x X y, making a right-handed set. These axes are arranged as
 * columns of a 3x3 matrix Rot_RF = [ x y z ] which is a direction cosine
 * (rotation) matrix useful for conversions between frame R and F. For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in R is given by vR = Rot_RF*vF. The origin point Loc_RF is 
 * stored as the vector Loc_RF=(OF-OR) and expressed in R.
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
    Frame() : Rot_RF(), Loc_RF(0.) { }
    Frame(const MatRotation& axesInR, const Vec3& originInR)
        : Rot_RF(axesInR), Loc_RF(originInR) { }
    explicit Frame(const MatRotation& axesInR)
        : Rot_RF(axesInR), Loc_RF(0.) { }
    explicit Frame(const Vec3& originInR) 
        : Rot_RF(), Loc_RF(originInR) { }
    // default copy, assignment, destructor

    void setFrame(const MatRotation& axesInR, const Vec3& originInR) 
     { Rot_RF=axesInR; Loc_RF=originInR; }

    // Transform various items measured and expressed in F to those same
    // items measured and expressed in F's reference frame R.
    Vec3  xformVector2Ref  (const Vec3& vF)      const { return Rot_RF*vF; }
    Vec3  xformStation2Ref (const Vec3& sF)      const { return Loc_RF + xformVector2Ref(sF); }
    MatRotation xformRotation2Ref(const MatRotation& Rot_FX) const 
      { return Rot_RF*Rot_FX; }
    Frame xformFrame2Ref(const Frame& fF) const 
      { return Frame(xformRotation2Ref(fF.Rot_RF), xformStation2Ref(fF.Loc_RF)); }

    const MatRotation& getAxes() const { return Rot_RF; }
    MatRotation&       updAxes()       { return Rot_RF; }

    const Vec3&  getOrigin() const { return Loc_RF; }
    Vec3&        updOrigin()       { return Loc_RF; }

    // Computation-free conversions
    const Real*  getFrameAsArray(const Frame& f) const {return reinterpret_cast<const Real*>(&f);}
    Real*        updFrameAsArray(Frame& f)             {return reinterpret_cast<Real*>(&f);}
    const Frame& getArrayAsFrame(const Real* r)  const {return *reinterpret_cast<const Frame*>(r);}
    Frame&       updArrayAsFrame(Real* r)              {return *reinterpret_cast<Frame*>(r);}

private:
    MatRotation Rot_RF;   // rotation matrix that expresses F's axes in R
    Vec3        Loc_RF;   // location of F's origin measured from R's origin, expressed in R 
};

} // namespace simtk

#endif /* ORIENTATION_H_ */
