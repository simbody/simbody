//-----------------------------------------------------------------------------
// File:     Transform.h
// Class:    Transform and InverseTransform 
// Parent:   None:  Data contains Rotation (for orientation) and Vec3 (translation)
// Purpose:  Transform (orientation and translation) relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------
#ifndef SimTK_TRANSFORM_H 
#define SimTK_TRANSFORM_H 

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
 * Contributors: Paul Mitiguy                                                 *
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

//-----------------------------------------------------------------------------
#include "SimTKcommon/SmallMatrix.h"
#include "SimTKcommon/internal/BigMatrix.h"
#include "SimTKcommon/internal/UnitVec.h"
#include "SimTKcommon/internal/Quaternion.h"
#include "SimTKcommon/internal/Rotation.h"
//-----------------------------------------------------------------------------
#include <iosfwd>  // Forward declaration of iostream
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
namespace SimTK {

//-----------------------------------------------------------------------------
// Forward declarations
class InverseTransform;


//-----------------------------------------------------------------------------
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
 * stored as the translation vector p_BF=(OF-OB) and expressed in B.
 *
 * Transform is designed to behave as much as possible like the computer
 * graphics 4x4 transform X which would be arranged like this:
 *
 * @verbatim
 *
 *         [       |   ]
 *     X = [   R   | p ]    R is a 3x3 orthogonal rotation matrix
 *         [.......|...]    p os a 3x1 translation vector
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
 *   X   = [  ~R   | p* ]   ~R is R transpose, p* = ~R(-p).
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
 * exactly 4 consecutive, packed 3-vectors in the order x,y,z,p.
 * That is, this class is equivalent to an array of 12 Reals with 
 * the order x1,x2,x3,y1,y2,y3,z1,z2,z3,p1,p2,p3. It is expressly allowed
 * to reinterpret Transform objects in any appropriate manner that depends
 * on this memory layout.
 */
//-----------------------------------------------------------------------------
class Transform {
public:
    /// Default constructor gives an identity transform.
    Transform() : R_BF(),  p_BF(0) { }

    /// Combine a rotation and a translation into a transform.
    Transform( const Rotation& R, const Vec3& p ) : R_BF(R), p_BF(p) { }

    /// Construct or default-convert a rotation into a transform
    /// containing that rotation and zero translation.
    Transform( const Rotation& R ) : R_BF(R), p_BF(0) { }

    /// Construct or default-convert a translation (expressed as a Vec3)
    /// into a transform with that translation and a zero rotation.
    Transform( const Vec3& p ) : R_BF(),  p_BF(p) { }

    // default copy, assignment, destructor

    /// Assignment from InverseTransform. This means that the 
    /// transform we're assigning to must end up with the same @em meaning
    /// as the inverse transform X has, so we'll need to end up with:
    ///   @li  p == X.p()
    ///   @li  R == X.R()
    ///
    /// Cost: one frame conversion and a negation, 18 flops.
    // (Definition is below after InverseTransform is declared.)
    inline Transform&  operator=( const InverseTransform& X );

    /// Assign a new value to this transform, explicitly providing
    /// the rotation and translation separately. We return a reference to
    /// the now-modified transform as though this were an assignment operator.
    Transform&  set( const Rotation& R, const Vec3& p ) { p_BF=p; R_BF=R; return *this; }

    /// By zero we mean "zero transform", i.e., an identity rotation
    /// and zero translation. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform&  setToZero()  { R_BF.setRotationToIdentityMatrix();  p_BF = 0.;  return *this; }

    /// This fills both the rotation and translation with NaNs. Note: this is
    /// @em not the same as a default-constructed transform, which is a
    /// legitimate identity transform instead. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform&  setToNaN()  { R_BF.setRotationToNaN();  p_BF.setToNaN();  return *this; }

    /// Return a read-only inverse of the current Transform, simply by casting it to
    /// the InverseTransform type. Zero cost.
    const InverseTransform&  invert() const  { return *reinterpret_cast<const InverseTransform*>(this); }

    /// Return a writable (lvalue) inverse of the current transform, simply by casting it to
    /// the InverseTransform type. That is, this is an lvalue. Zero cost.
    InverseTransform&  updInvert()  { return *reinterpret_cast<InverseTransform*>(this); }

    /// Overload transpose operator to mean inversion. @see invert
    const InverseTransform&  operator~() const  {return invert();}

    /// Overload transpose operator to mean inversion. @see updInvert
    InverseTransform&        operator~()        {return updInvert();}

    /// Compose the current transform (X_BF) with the given one. That is,
    /// return X_BY=X_BF*X_FY. Cost is 63 flops.
    Transform compose(const Transform& X_FY) const {
        return Transform( R_BF * X_FY.R(),  p_BF + R_BF * X_FY.p() );
    }

    /// Compose the current transform (X_BF) with one that is supplied
    /// as an InverseTransform (typically as a result of applying
    /// the "~" operator to a transform). That is, return 
    /// X_BY=X_BF*X_FY, but now X_FY is represented as ~X_YF. Cost
    /// is an extra 18 flops to calculate X_FY.p(), total 81 flops.
    // (Definition is below after InverseTransform is declared.)
    inline Transform  compose( const InverseTransform& X_FY ) const;

    /// %Transform a vector expressed in our "F" frame to our "B" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec3  xformFrameVecToBase( const Vec3& vF ) const {return R_BF*vF;}

    /// %Transform a vector expressed in our "B" frame to our "F" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec3  xformBaseVecToFrame( const Vec3& vB ) const  { return ~R_BF*vB; }

    /// %Transform a point (station) measured from and expressed in
    /// our "F" frame to that same point but measured from and
    /// expressed in our "B" frame. Cost is 18 flops.
    Vec3  shiftFrameStationToBase( const Vec3& sF ) const { return p_BF + xformFrameVecToBase(sF); }

    /// %Transform a point (station) measured from and expressed in
    /// our "B" frame to that same point but measured from and
    /// expressed in our "F" frame. Cost is 18 flops.
    Vec3  shiftBaseStationToFrame( const Vec3& sB ) const { return xformBaseVecToFrame(sB - p_BF); }

    /// Return a read-only reference to the contained rotation R_BF.
    const Rotation&  R() const  { return R_BF; }

    /// Return a writable (lvalue) reference to the contained rotation R_BF.
    Rotation&  updR()           { return R_BF; }

    /// Return a read-only reference to the x direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType&  x() const  { return R().x(); }
    /// Return a read-only reference to the y direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType&  y() const  { return R().y(); }
    /// Return a read-only reference to the z direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const Rotation::ColType&  z() const  { return R().z(); }

    /// Return a read-only reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    const InverseRotation&  RInv() const  { return ~R_BF; }

    /// Return a writable (lvalue) reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    InverseRotation&  updRInv()  { return ~R_BF; }

    /// Return a read-only reference to our translation vector p_BF.
    const Vec3&  p() const  { return p_BF; }

    /// Return a writable (lvalue) reference to our translation vector p_BF.
    /// Caution: if you write through this reference you update the transform.
    Vec3&  updP()  { return p_BF; }

    /// Assign a new value to our translation vector. We expect the
    /// supplied vector @p p to be expressed in our B frame. A reference
    /// to the now-modified transform is returned as though this were
    /// an assignment operator.
    Transform&  setP( const Vec3& p )  { p_BF=p; return *this; }

    /// Calculate the inverse of the translation vector in this transform.
    /// The returned vector will be the negative of the original and will
    /// be expressed in the F frame rather than our B frame. Cost is 18 flops.
    Vec3  pInv() const  { return -(~R_BF*p_BF); }

    /// Assign a value to the @em inverse of our translation vector.
    /// That is, we're given a vector in F which we invert and reexpress
    /// in B to store it in p, so that we get the original argument back if
    /// we ask for the inverse of p. Sorry, can't update pInv as an lvalue, but here we
    /// want -(~R_BF*p_BF)=p_FB => p_BF=-(R_BF*p_FB) so we can calculate
    /// it in 18 flops. A reference to the now-modified transform is returned
    /// as though this were an assignment operator.
    Transform&  setPInv( const Vec3& p_FB )  { p_BF = -(R_BF*p_FB); return *this; }

    /// Recast this transform as a read-only 3x4 matrix. This is a zero-cost
    /// reinterpretation of the data; the first three columns are the
    /// columns of the rotation and the last column is the translation.
    const Mat34&  asMat34() const  { return Mat34::getAs(reinterpret_cast<const Real*>(this)); }

    /// Less efficient version of asMat34(); copies into return variable.
    Mat34  toMat34() const  { return asMat34(); }

    /// Return the equivalent 4x4 transformation matrix.
    Mat44 toMat44() const {
        Mat44 tmp;
        tmp.updSubMat<3,4>(0,0) = asMat34();
        tmp[3]                  = Row4(0,0,0,1);
        return tmp;
    }

    // OBSOLETE -- alternate name for p
    const Vec3& T() const {return p();}
    Vec3&  updT()  {return updP();}

private:
    //TODO: these might not pack correctly; should use an array of 12 Reals.
    Rotation R_BF;   // rotation matrix that expresses F's axes in R
    Vec3     p_BF;   // location of F's origin measured from B's origin, expressed in B 
};


//-----------------------------------------------------------------------------
/**
 * %Transform from frame B to frame F, but with the internal representation inverted.
 * That is, we store R*,p* here but the transform this represents is
 * @verbatim
 *
 *  B F    [       |   ]
 *   X   = [   R   | p ]   where R=~(R*), p = - ~(R*)(p*).
 *         [.......|...]
 *         [ 0 0 0   1 ] 
 *
 * @endverbatim
 */
//-----------------------------------------------------------------------------
class InverseTransform {
public:
    InverseTransform() : R_FB(), p_FB(0) { }
    // default copy, assignment, destructor

    // Implicit conversion to Transform
    operator Transform() const  { return Transform( R(), p() ); }

    // Assignment from Transform. This means that the inverse
    // transform we're assigning to must end up with the same meaning
    // as the inverse transform X has, so we'll need:
    //          p* == X.pInv()
    //          R* == X.RInv()
    // Cost: one frame conversion and a negation for pInv, 18 flops.
    InverseTransform&  operator=( const Transform& X ) {
        // Be careful to do this in the right order in case X and this
        // are the same object, i.e. ~X = X which is weird but has
        // the same meaning as X = ~X, i.e. invert X in place.
        p_FB = X.pInv(); // This might change X.p ...
        R_FB = X.RInv(); // ... but this doesn't depend on X.p.
        return *this;
    }

    // Inverting one of these just recasts it back to a Transform.
    const Transform&  invert() const  { return *reinterpret_cast<const Transform*>(this); }
    Transform&  updInvert()           { return *reinterpret_cast<Transform*>(this); }

    // Overload transpose to mean inversion.
    const Transform&  operator~() const  { return invert(); }
    Transform&        operator~()        { return updInvert(); }

    // Return X_BY=X_BF*X_FY, where X_BF (this) is represented here as ~X_FB. This
    // costs exactly the same as a composition of two Transforms (63 flops).
    Transform  compose(const Transform& X_FY) const {
        return Transform( ~R_FB * X_FY.R(),  ~R_FB *(X_FY.p() - p_FB) );
    }
    // Return X_BY=X_BF*X_FY, but now both xforms are represented by their inverses.
    // This costs one extra vector transformation and a negation (18 flops) more
    // than a composition of two Transforms, for a total of 81 flops.
    Transform  compose(const InverseTransform& X_FY) const { 
        return Transform(  ~R_FB * X_FY.R(),  ~R_FB *(X_FY.p() - p_FB)  ); 
    }

    // Forward and inverse vector transformations cost the same here as
    // for a Transform (or for that matter, a Rotation): 15 flops.
    Vec3  xformFrameVecToBase(const Vec3& vF) const {return ~R_FB*vF;}
    Vec3  xformBaseVecToFrame(const Vec3& vB) const {return  R_FB*vB;}

    // Forward and inverse station shift & transform cost the same here as for a Transform: 18 flops.
    Vec3  shiftFrameStationToBase(const Vec3& sF) const  { return ~R_FB*(sF-p_FB); }
    Vec3  shiftBaseStationToFrame(const Vec3& sB) const  { return R_FB*sB + p_FB; }
    
    const InverseRotation&  R() const  {return ~R_FB;}
    InverseRotation&        updR()     {return ~R_FB;}

    const InverseRotation::ColType&  x() const  {return R().x();}
    const InverseRotation::ColType&  y() const  {return R().y();}
    const InverseRotation::ColType&  z() const  {return R().z();}

    const Rotation&  RInv() const  {return R_FB;}
    Rotation&        updRInv()     {return R_FB;}

    // Costs 18 flops to look at the real translation vector.
    Vec3  p() const  { return -(~R_FB*p_FB); }


    // no updP lvalue

    // Sorry, can't update translation as an lvalue, but here we
    // want -(R_BF*p_FB)=p_BF => p_FB=-(R_FB*p_BF). Cost: 18 flops.
    void  setP( const Vec3& p_BF )  { p_FB = -(R_FB*p_BF); }

    // Inverse translation is free.
    const Vec3&  pInv() const              { return p_FB; }
    void         setPInv( const Vec3& p )  { p_FB = p; }

    /// For compatibility with Transform, but we don't provide an "as"
    /// method here since the internal storage layout is somewhat odd.
    Mat34  toMat34() const  { return Transform(*this).asMat34(); }

    /// Return the equivalent 4x4 transformation matrix.
    Mat44  toMat44() const  { return Transform(*this).toMat44(); }

    // OBSOLETE -- alternate name for p.
    Vec3 T() const {return p();}

private:
    // DATA LAYOUT MUST BE IDENTICAL TO Transform !!
    // TODO: redo packing here when it is done for Transform.
    Rotation R_FB; // transpose of our rotation matrix, R_BF
    Vec3     p_FB; // our translation is -(R_BF*p_FB)=-(~R_FB*p_FB)
};


/// If we multiply a transform by a 3-vector, we treat it as though it had a 4th element "1" appended,
/// that is, it is treated as a *station* rather than a *vector*.
inline Vec3  operator*( const Transform& X_BF,        const Vec3& s_F )  { return X_BF.shiftFrameStationToBase(s_F); }
inline Vec3  operator*( const InverseTransform& X_BF, const Vec3& s_F )  { return X_BF.shiftFrameStationToBase(s_F); }

/// If we multiply a transform by an augmented 4-vector, we use the 4th element to decide how to treat it.
/// The 4th element must be 0 or 1. If 0 it is treated as a vector only and the translation is ignored. 
/// If 1 it is treated as a station and rotated & shifted.
//-----------------------------------------------------------------------
inline Vec4  operator*( const Transform& X_BF, const Vec4& a_F ) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec3& v_F = Vec3::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec4 out;
    if( a_F[3] == 0 ) { Vec3::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);      out[3] = 0; } 
    else              { Vec3::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);  out[3] = 1; }
    return out;
}

//-----------------------------------------------------------------------
inline Vec4  operator*( const InverseTransform& X_BF, const Vec4& a_F ) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec3& v_F = Vec3::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec4 out;
    if( a_F[3] == 0 ) { Vec3::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);      out[3] = 0; } 
    else              { Vec3::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);  out[3] = 1; }
    return out;
}

/// Multiplying a matrix or vector by a Transform applies it to each element individually.
//@{
template <class E>
inline Vector_<E> operator*(const Transform& X, const VectorBase<E>& v) {
    Vector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class E>
inline Vector_<E> operator*(const VectorBase<E>& v, const Transform& X) {
    Vector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class E>
inline RowVector_<E> operator*(const Transform& X, const RowVectorBase<E>& v) {
    RowVector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class E>
inline RowVector_<E> operator*(const RowVectorBase<E>& v, const Transform& X) {
    RowVector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class E>
inline Matrix_<E> operator*(const Transform& X, const MatrixBase<E>& v) {
    Matrix_<E> result(v.nrow(), v.ncol());
    for (int i = 0; i < v.nrow(); ++i)
        for (int j = 0; j < v.ncol(); ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <class E>
inline Matrix_<E> operator*(const MatrixBase<E>& v, const Transform& X) {
    Matrix_<E> result(v.nrow(), v.ncol());
    for (int i = 0; i < v.nrow(); ++i)
        for (int j = 0; j < v.ncol(); ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <int N, class E>
inline Vec<N,E> operator*(const Transform& X, const Vec<N,E>& v) {
    Vec<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <int N, class E>
inline Vec<N,E> operator*(const Vec<N,E>& v, const Transform& X) {
    Vec<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <int N, class E>
inline Row<N,E> operator*(const Transform& X, const Row<N,E>& v) {
    Row<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <int N, class E>
inline Row<N,E> operator*(const Row<N,E>& v, const Transform& X) {
    Row<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <int M, int N, class E>
inline Mat<M,N,E> operator*(const Transform& X, const Mat<M,N,E>& v) {
    Mat<M,N,E> result;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <int M, int N, class E>
inline Mat<M,N,E> operator*(const Mat<M,N,E>& v, const Transform& X) {
    Mat<M,N,E> result;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            result(i, j) = X*v(i, j);
    return result;
}

//@}

// These Transform definitions had to wait for InverseTransform to be declared.
inline Transform&  Transform::operator=( const InverseTransform& X ) {
    // Be careful to do this in the right order in case X and this
    // are the same object, i.e. we're doing X = ~X, inverting X in place.
    p_BF = X.p(); // This might change X.p ...
    R_BF = X.R(); // ... but this doesn't depend on X.p.
    return *this;
}

inline Transform  Transform::compose( const InverseTransform& X_FY ) const {
    return Transform( R_BF * X_FY.R(), p_BF + R_BF * X_FY.p() );
}

inline Transform  operator*( const Transform& X1,        const Transform& X2 )         { return X1.compose(X2); }
inline Transform  operator*( const Transform& X1,        const InverseTransform& X2 )  { return X1.compose(X2); }
inline Transform  operator*( const InverseTransform& X1, const Transform& X2 )         { return X1.compose(X2); }
inline Transform  operator*( const InverseTransform& X1, const InverseTransform& X2 )  { return X1.compose(X2); }

inline bool  operator==( const Transform& X1,        const Transform& X2 )         { return X1.R()==X2.R() && X1.p()==X2.p(); }
inline bool  operator==( const InverseTransform& X1, const InverseTransform& X2 )  { return X1.R()==X2.R() && X1.p()==X2.p(); }
inline bool  operator==( const Transform& X1,        const InverseTransform& X2 )  { return X1.R()==X2.R() && X1.p()==X2.p(); }
inline bool  operator==( const InverseTransform& X1, const Transform& X2 )         { return X1.R()==X2.R() && X1.p()==X2.p(); }

SimTK_SimTKCOMMON_EXPORT std::ostream&  operator<<( std::ostream& o, const Transform& );


//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_TRANSFORM_H_
//--------------------------------------------------------------------------

