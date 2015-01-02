#ifndef SimTK_SimTKCOMMON_TRANSFORM_H 
#define SimTK_SimTKCOMMON_TRANSFORM_H 

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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
// Forward declarations (everything is templatized by precision).
template <class P> class Transform_;
template <class P> class InverseTransform_;

typedef Transform_<Real>    Transform;
typedef Transform_<float>   fTransform;
typedef Transform_<double>  dTransform;


//-----------------------------------------------------------------------------
/**
 * This class represents the rotate-and-shift transform which gives the 
 * location and orientation of a new frame F in a base (reference) frame
 * B. A frame is an orthogonal, right-handed set of three axes, and an
 * origin point. A transform X from frame B to F consists of 3 perpendicular
 * unit vectors defining F's axes as viewed from B (that is, as expressed in 
 * the basis formed by B's axes), and a vector from B's origin point Bo to F's
 * origin point Fo. Note that the meaning of "B" comes from the context in
 * which the transform is used. We use the phrase "frame F is in frame B" to
 * describe the above relationship, that is, "in" means both measured from
 * and expressed in. 
 *
 * The axis vectors constitute a Rotation_. They are ordered 1-2-3 or x-y-z
 * as you prefer, with z = x X y, making a right-handed set. These axes are 
 * arranged as columns of a 3x3 rotation matrix R_BF = [ x y z ] which is a 
 * direction cosine (rotation) matrix useful for conversions between frame 
 * B and F. (The columns of R_BF are F's coordinate axes, expressed in B.) For
 * example, given a vector vF expressed in the F frame, that same vector
 * re-expressed in B is given by vB = R_BF*vF. F's origin point OF is 
 * stored as the translation vector p_BF=(Fo-Bo) and expressed in B.
 *
 * Transform is designed to behave as much as possible like the computer
 * graphics 4x4 affine transform X which would be arranged like this:
 * <pre>
 *
 *         [       |   ]
 *     X = [   R   | p ]    R is a 3x3 orthogonal rotation matrix
 *         [.......|...]    p is a 3x1 translation vector
 *         [ 0 0 0   1 ]
 * </pre>
 *
 * These can be composed directly by matrix multiplication, but more 
 * importantly they have a particularly simple inverse:
 * <pre>
 *
 *    -1   [       |    ]
 *   X   = [  ~R   | p* ]   ~R is R transpose, p* = ~R(-p).
 *         [.......|....]
 *         [ 0 0 0   1  ] 
 * </pre>
 *
 * This inverse is so simple that we compute it simply by defining another
 * type, InverseTransform_, which is identical to %Transform_ in memory but
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
template <class P>
class Transform_ {
public:
    /// Default constructor gives an identity transform.
    Transform_() : R_BF(),  p_BF(0) { }

    /// Combine a rotation and a translation into a transform.
    Transform_( const Rotation_<P>& R, const Vec<3,P>& p ) : R_BF(R), p_BF(p) { }

    /// Construct or default-convert a rotation into a transform
    /// containing that rotation and zero translation.
    Transform_( const Rotation_<P>& R ) : R_BF(R), p_BF(0) { }

    /// Construct or default-convert a translation (expressed as a Vec3)
    /// into a transform with that translation and a zero rotation.
    Transform_( const Vec<3,P>& p ) : R_BF(),  p_BF(p) { }

    // default copy, assignment, destructor

    /// Assignment from InverseTransform. This means that the 
    /// transform we're assigning to must end up with the same @em meaning
    /// as the inverse transform X has, so we'll need to end up with:
    ///   @li  p == X.p()
    ///   @li  R == X.R()
    ///
    /// Cost: one frame conversion and a negation, 18 flops.
    // (Definition is below after InverseTransform is declared.)
    inline Transform_&  operator=( const InverseTransform_<P>& X );

    /// Add an offset to the position vector in this %Transform. Cost 
    /// is 3 flops.
    template <int S>
    Transform_& operator+=(const Vec<3,P,S>& offset_B)
    {   p_BF += offset_B; return *this; }

    /// Subtract an offset from the position vector in this %Transform. Cost
    /// is 3 flops.
    template <int S>
    Transform_& operator-=(const Vec<3,P,S>& offset_B)
    {   p_BF -= offset_B; return *this; }

    /// Assign a new value to this transform, explicitly providing
    /// the rotation and translation separately. We return a reference to
    /// the now-modified transform as though this were an assignment operator.
    Transform_&  set( const Rotation_<P>& R, const Vec<3,P>& p ) { p_BF=p; R_BF=R; return *this; }

    /// By zero we mean "zero transform", i.e., an identity rotation
    /// and zero translation. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform_&  setToZero()  { R_BF.setRotationToIdentityMatrix();  p_BF = P(0);  return *this; }

    /// This fills both the rotation and translation with NaNs. Note: this is
    /// @em not the same as a default-constructed transform, which is a
    /// legitimate identity transform instead. We return a reference to the now-modified
    /// transform as though this were an assignment operator. 
    Transform_&  setToNaN()  { R_BF.setRotationToNaN();  p_BF.setToNaN();  return *this; }

    /// Return a read-only inverse of the current Transform_<P>, simply by casting it to
    /// the InverseTransform_<P> type. Zero cost.
    const InverseTransform_<P>&  invert() const  { return *reinterpret_cast<const InverseTransform_<P>*>(this); }

    /// Return a writable (lvalue) inverse of the current transform, simply by casting it to
    /// the InverseTransform_<P> type. That is, this is an lvalue. Zero cost.
    InverseTransform_<P>&  updInvert()  { return *reinterpret_cast<InverseTransform_<P>*>(this); }

    /// Overload transpose operator to mean inversion. @see invert
    const InverseTransform_<P>&  operator~() const  {return invert();}

    /// Overload transpose operator to mean inversion. @see updInvert
    InverseTransform_<P>&        operator~()        {return updInvert();}

    /// Compose the current transform (X_BF) with the given one. That is,
    /// return X_BY=X_BF*X_FY. Cost is 63 flops.
    Transform_ compose(const Transform_& X_FY) const {
        return Transform_( R_BF * X_FY.R(),  p_BF + R_BF * X_FY.p() );
    }

    /// Compose the current transform (X_BF) with one that is supplied
    /// as an InverseTransform_ (typically as a result of applying
    /// the "~" operator to a transform). That is, return 
    /// X_BY=X_BF*X_FY, but now X_FY is represented as ~X_YF. Cost
    /// is an extra 18 flops to calculate X_FY.p(), total 81 flops.
    // (Definition is below after InverseTransform_ is declared.)
    inline Transform_  compose( const InverseTransform_<P>& X_FY ) const;

    /// %Transform a vector expressed in our "F" frame to our "B" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec<3,P>  xformFrameVecToBase( const Vec<3,P>& vF ) const {return R_BF*vF;}

    /// %Transform a vector expressed in our "B" frame to our "F" frame.
    /// Note that this involves rotation only; it is independent of
    /// the translation stored in this transform. Cost is 15 flops.
    Vec<3,P>  xformBaseVecToFrame( const Vec<3,P>& vB ) const  { return ~R_BF*vB; }

    /// %Transform a point (station) measured from and expressed in
    /// our "F" frame to that same point but measured from and
    /// expressed in our "B" frame. Cost is 18 flops.
    Vec<3,P>  shiftFrameStationToBase( const Vec<3,P>& sF ) const 
    {   return p_BF + xformFrameVecToBase(sF); }

    /// %Transform a point (station) measured from and expressed in
    /// our "B" frame to that same point but measured from and
    /// expressed in our "F" frame. Cost is 18 flops.
    Vec<3,P>  shiftBaseStationToFrame( const Vec<3,P>& sB ) const 
    {   return xformBaseVecToFrame(sB - p_BF); }

    /// Return a read-only reference to the contained rotation R_BF.
    const Rotation_<P>&  R() const  { return R_BF; }

    /// Return a writable (lvalue) reference to the contained rotation R_BF.
    Rotation_<P>&  updR()           { return R_BF; }

    /// Return a read-only reference to the x direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const typename Rotation_<P>::ColType&  x() const  { return R().x(); }
    /// Return a read-only reference to the y direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const typename Rotation_<P>::ColType&  y() const  { return R().y(); }
    /// Return a read-only reference to the z direction (unit vector)
    /// of the F frame, expressed in the B frame.
    const typename Rotation_<P>::ColType&  z() const  { return R().z(); }

    /// Return a read-only reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    const InverseRotation_<P>&  RInv() const  { return ~R_BF; }

    /// Return a writable (lvalue) reference to the inverse (transpose) of
    /// our contained rotation, that is R_FB.
    InverseRotation_<P>&  updRInv()  { return ~R_BF; }

    /// Return a read-only reference to our translation vector p_BF.
    const Vec<3,P>&  p() const  { return p_BF; }

    /// Return a writable (lvalue) reference to our translation vector p_BF.
    /// Caution: if you write through this reference you update the transform.
    Vec<3,P>&  updP()  { return p_BF; }

    /// Assign a new value to our translation vector. We expect the
    /// supplied vector @p p to be expressed in our B frame. A reference
    /// to the now-modified transform is returned as though this were
    /// an assignment operator.
    Transform_<P>&  setP( const Vec<3,P>& p )  { p_BF=p; return *this; }

    /// Calculate the inverse of the translation vector in this transform.
    /// The returned vector will be the negative of the original and will
    /// be expressed in the F frame rather than our B frame. Cost is 18 flops.
    Vec<3,P>  pInv() const  { return -(~R_BF*p_BF); }

    /// Assign a value to the @em inverse of our translation vector.
    /// That is, we're given a vector in F which we invert and reexpress
    /// in B to store it in p, so that we get the original argument back if
    /// we ask for the inverse of p. Sorry, can't update pInv as an lvalue, but here we
    /// want -(~R_BF*p_BF)=p_FB => p_BF=-(R_BF*p_FB) so we can calculate
    /// it in 18 flops. A reference to the now-modified transform is returned
    /// as though this were an assignment operator.
    Transform_<P>&  setPInv( const Vec<3,P>& p_FB )  { p_BF = -(R_BF*p_FB); return *this; }

    /// Recast this transform as a read-only 3x4 matrix. This is a zero-cost
    /// reinterpretation of the data; the first three columns are the
    /// columns of the rotation and the last column is the translation.
    const Mat<3,4,P>&  asMat34() const  { return Mat<3,4,P>::getAs(reinterpret_cast<const P*>(this)); }

    /// Less efficient version of asMat34(); copies into return variable.
    Mat<3,4,P>  toMat34() const  { return asMat34(); }

    /// Return the equivalent 4x4 transformation matrix.
    Mat<4,4,P> toMat44() const {
        Mat<4,4,P> tmp;
        tmp.template updSubMat<3,4>(0,0) = asMat34();
        tmp[3]                  = Row<4,P>(0,0,0,1);
        return tmp;
    }

    // OBSOLETE -- alternate name for p
    const Vec<3,P>& T() const {return p();}
    Vec<3,P>&  updT()  {return updP();}

private:
    //TODO: these might not pack correctly; should use an array of 12 Reals.
    Rotation_<P> R_BF;   // rotation matrix that expresses F's axes in R
    Vec<3,P>     p_BF;   // location of F's origin measured from B's origin, expressed in B 
};


//-----------------------------------------------------------------------------
/**
 * %Transform from frame B to frame F, but with the internal representation 
 * inverted. That is, we store R*,p* here but the transform this represents is
 * <pre>
 *
 *             B F    [       |   ]
 *      X_BF =  X   = [   R   | p ]   where R=~(R*), p = - ~(R*)(p*).
 *                    [.......|...]
 *                    [ 0 0 0   1 ] 
 * </pre>
 */
//-----------------------------------------------------------------------------
template <class P>
class InverseTransform_ {
public:
    /// Default constructor produces an identity transform.
    InverseTransform_() : R_FB(), p_FB(0) { }

    // default copy, assignment, destructor

    /// Implicit conversion from %InverseTransform_ to Transform_.
    operator Transform_<P>() const  { return Transform_<P>( R(), p() ); }

    // Assignment from Transform_. This means that the inverse
    // transform we're assigning to must end up with the same meaning
    // as the inverse transform X has, so we'll need:
    //          p* == X.pInv()
    //          R* == X.RInv()
    // Cost: one frame conversion and a negation for pInv, 18 flops.
    InverseTransform_&  operator=( const Transform_<P>& X ) {
        // Be careful to do this in the right order in case X and this
        // are the same object, i.e. ~X = X which is weird but has
        // the same meaning as X = ~X, i.e. invert X in place.
        p_FB = X.pInv(); // This might change X.p ...
        R_FB = X.RInv(); // ... but this doesn't depend on X.p.
        return *this;
    }

    // Inverting one of these just recasts it back to a Transform_<P>.
    const Transform_<P>&  invert() const  { return *reinterpret_cast<const Transform_<P>*>(this); }
    Transform_<P>&  updInvert()           { return *reinterpret_cast<Transform_<P>*>(this); }

    // Overload transpose to mean inversion.
    const Transform_<P>&  operator~() const  { return invert(); }
    Transform_<P>&        operator~()        { return updInvert(); }

    // Return X_BY=X_BF*X_FY, where X_BF (this) is represented here as ~X_FB. This
    // costs exactly the same as a composition of two Transforms (63 flops).
    Transform_<P>  compose(const Transform_<P>& X_FY) const {
        return Transform_<P>( ~R_FB * X_FY.R(),  ~R_FB *(X_FY.p() - p_FB) );
    }
    // Return X_BY=X_BF*X_FY, but now both xforms are represented by their inverses.
    // This costs one extra vector transformation and a negation (18 flops) more
    // than a composition of two Transforms, for a total of 81 flops.
    Transform_<P>  compose(const InverseTransform_<P>& X_FY) const { 
        return Transform_<P>(  ~R_FB * X_FY.R(),  ~R_FB *(X_FY.p() - p_FB)  ); 
    }

    // Forward and inverse vector transformations cost the same here as
    // for a Transform_<P> (or for that matter, a Rotation_<P>): 15 flops.
    Vec<3,P>  xformFrameVecToBase(const Vec<3,P>& vF) const {return ~R_FB*vF;}
    Vec<3,P>  xformBaseVecToFrame(const Vec<3,P>& vB) const {return  R_FB*vB;}

    // Forward and inverse station shift & transform cost the same here as for a Transform_<P>: 18 flops.
    Vec<3,P>  shiftFrameStationToBase(const Vec<3,P>& sF) const  { return ~R_FB*(sF-p_FB); }
    Vec<3,P>  shiftBaseStationToFrame(const Vec<3,P>& sB) const  { return R_FB*sB + p_FB; }
    
    const InverseRotation_<P>&  R() const  {return ~R_FB;}
    InverseRotation_<P>&        updR()     {return ~R_FB;}

    const typename InverseRotation_<P>::ColType&  x() const  {return R().x();}
    const typename InverseRotation_<P>::ColType&  y() const  {return R().y();}
    const typename InverseRotation_<P>::ColType&  z() const  {return R().z();}

    const Rotation_<P>&  RInv() const  {return R_FB;}
    Rotation_<P>&        updRInv()     {return R_FB;}

    /// Calculate the actual translation vector at a cost of 18 flops.
    /// It is better if you can just work with the %InverseTransform
    /// directly since then you'll never have to pay this cost.
    Vec<3,P>  p() const  { return -(~R_FB*p_FB); }


    // no updP lvalue

    // Sorry, can't update translation as an lvalue, but here we
    // want -(R_BF*p_FB)=p_BF => p_FB=-(R_FB*p_BF). Cost: 18 flops.
    void  setP( const Vec<3,P>& p_BF )  { p_FB = -(R_FB*p_BF); }

    // Inverse translation is free.
    const Vec<3,P>&  pInv() const              { return p_FB; }
    void         setPInv( const Vec<3,P>& p )  { p_FB = p; }

    /// For compatibility with Transform_<P>, but we don't provide an "as"
    /// method here since the internal storage layout is somewhat odd.
    Mat<3,4,P>  toMat34() const  { return Transform_<P>(*this).asMat34(); }

    /// Return the equivalent 4x4 transformation matrix.
    Mat<4,4,P>  toMat44() const  { return Transform_<P>(*this).toMat44(); }

    // OBSOLETE -- alternate name for p.
    Vec<3,P> T() const {return p();}

private:
    // DATA LAYOUT MUST BE IDENTICAL TO Transform_<P> !!
    // TODO: redo packing here when it is done for Transform_<P>.
    Rotation_<P> R_FB; // transpose of our rotation matrix, R_BF
    Vec<3,P>     p_FB; // our translation is -(R_BF*p_FB)=-(~R_FB*p_FB)
};

/// If we multiply a transform or inverse transform by a 3-vector, we treat 
/// the vector as though it had a 4th element "1" appended, that is, it is 
/// treated as a \e station rather than a \e vector. This way we use both
/// the rotational and translational components of the transform.
/// @relates Transform_
template <class P, int S> inline Vec<3,P>  
operator*(const Transform_<P>& X_BF,        const Vec<3,P,S>& s_F)  
{   return X_BF.shiftFrameStationToBase(s_F); }
template <class P, int S> inline Vec<3,P>  
operator*(const InverseTransform_<P>& X_BF, const Vec<3,P,S>& s_F)  
{   return X_BF.shiftFrameStationToBase(s_F); }
template <class P, int S> inline Vec<3,P>  
operator*(const Transform_<P>& X_BF,        const Vec<3,negator<P>,S>& s_F)  
{   return X_BF*Vec<3,P>(s_F); }
template <class P, int S> inline Vec<3,P>  
operator*(const InverseTransform_<P>& X_BF, const Vec<3,negator<P>,S>& s_F)  
{   return X_BF*Vec<3,P>(s_F); }

/// Adding a 3-vector to a Transform produces a new shifted transform.
/// @relates Transform_
template <class P, int S> inline Transform_<P>
operator+(const Transform_<P>& X_BF, const Vec<3,P,S>& offset_B)
{   return Transform_<P>(X_BF) += offset_B; }
/// Adding a 3-vector to a Transform produces a new shifted transform.
/// @relates Transform_
template <class P, int S> inline Transform_<P>
operator+(const Vec<3,P,S>& offset_B, const Transform_<P>& X_BF)
{   return Transform_<P>(X_BF) += offset_B; }

/// Subtracting a 3-vector from a Transform produces a new shifted transform.
/// @relates Transform_
template <class P, int S> inline Transform_<P>
operator-(const Transform_<P>& X_BF, const Vec<3,P,S>& offset_B)
{   return Transform_<P>(X_BF) -= offset_B; }

//-----------------------------------------------------------------------------
/// If we multiply a transform or inverse transform by an augmented 4-vector, 
/// we use the 4th element to decide how to treat it. The 4th element must be 
/// 0 or 1. If 0 it is treated as a vector only and the translation is ignored. 
/// If 1 it is treated as a station and rotated & shifted.
/// @relates Transform_
template <class P, int S> inline Vec<4,P> 
operator*(const Transform_<P>& X_BF, const Vec<4,P,S>& a_F) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec<3,P,S>& v_F = Vec<3,P,S>::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec<4,P> out;
    if( a_F[3] == 0 ) { Vec<3,P>::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);      out[3] = 0; } 
    else              { Vec<3,P>::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);  out[3] = 1; }
    return out;
}

template <class P, int S> inline Vec<4,P> 
operator*(const InverseTransform_<P>& X_BF, const Vec<4,P,S>& a_F ) {
    assert(a_F[3]==0 || a_F[3]==1);
    const Vec<3,P,S>& v_F = Vec<3,P,S>::getAs(&a_F[0]); // recast the 1st 3 elements as Vec3

    Vec<4,P> out;
    if( a_F[3] == 0 ) { Vec<3,P>::updAs(&out[0]) = X_BF.xformFrameVecToBase(v_F);      out[3] = 0; } 
    else              { Vec<3,P>::updAs(&out[0]) = X_BF.shiftFrameStationToBase(v_F);  out[3] = 1; }
    return out;
}
template <class P, int S> inline Vec<4,P>  
operator*(const Transform_<P>& X_BF,        const Vec<4,negator<P>,S>& s_F)  {return X_BF*Vec<4,P>(s_F);}
template <class P, int S> inline Vec<4,P>  
operator*(const InverseTransform_<P>& X_BF, const Vec<4,negator<P>,S>& s_F)  {return X_BF*Vec<4,P>(s_F);}
//-----------------------------------------------------------------------------

/// Multiplying a matrix or vector by a Transform_<P> applies it to each element 
/// individually.
/// @relates Transform_
template <class P, class E> inline Vector_<E> 
operator*(const Transform_<P>& X, const VectorBase<E>& v) {
    Vector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, class E> inline Vector_<E> 
operator*(const VectorBase<E>& v, const Transform_<P>& X) {
    Vector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, class E> inline RowVector_<E> 
operator*(const Transform_<P>& X, const RowVectorBase<E>& v) {
    RowVector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, class E> inline RowVector_<E> 
operator*(const RowVectorBase<E>& v, const Transform_<P>& X) {
    RowVector_<E> result(v.size());
    for (int i = 0; i < v.size(); ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, class E> inline Matrix_<E> 
operator*(const Transform_<P>& X, const MatrixBase<E>& v) {
    Matrix_<E> result(v.nrow(), v.ncol());
    for (int i = 0; i < v.nrow(); ++i)
        for (int j = 0; j < v.ncol(); ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <class P, class E> inline Matrix_<E> 
operator*(const MatrixBase<E>& v, const Transform_<P>& X) {
    Matrix_<E> result(v.nrow(), v.ncol());
    for (int i = 0; i < v.nrow(); ++i)
        for (int j = 0; j < v.ncol(); ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <class P, int N, class E, int S> inline Vec<N,E> 
operator*(const Transform_<P>& X, const Vec<N,E,S>& v) {
    Vec<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, int N, class E, int S> inline Vec<N,E> 
operator*(const Vec<N,E,S>& v, const Transform_<P>& X) {
    Vec<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, int N, class E, int S> inline Row<N,E> 
operator*(const Transform_<P>& X, const Row<N,E,S>& v) {
    Row<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, int N, class E, int S> inline Row<N,E> 
operator*(const Row<N,E,S>& v, const Transform_<P>& X) {
    Row<N,E> result;
    for (int i = 0; i < N; ++i)
        result[i] = X*v[i];
    return result;
}
template <class P, int M, int N, class E, int CS, int RS> inline Mat<M,N,E> 
operator*(const Transform_<P>& X, const Mat<M,N,E,CS,RS>& v) {
    Mat<M,N,E> result;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            result(i, j) = X*v(i, j);
    return result;
}
template <class P, int M, int N, class E, int CS, int RS> inline Mat<M,N,E> 
operator*(const Mat<M,N,E,CS,RS>& v, const Transform_<P>& X) {
    Mat<M,N,E> result;
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            result(i, j) = X*v(i, j);
    return result;
}

// These Transform definitions had to wait for InverseTransform to be declared.

template <class P> inline Transform_<P>&
Transform_<P>::operator=( const InverseTransform_<P>& X ) {
    // Be careful to do this in the right order in case X and this
    // are the same object, i.e. we're doing X = ~X, inverting X in place.
    p_BF = X.p(); // This might change X.p ...
    R_BF = X.R(); // ... but this doesn't depend on X.p.
    return *this;
}

template <class P> inline Transform_<P>
Transform_<P>::compose( const InverseTransform_<P>& X_FY ) const {
    return Transform_<P>( R_BF * X_FY.R(), p_BF + R_BF * X_FY.p() );
}

/// Composition of transforms. Operators are provided for all the combinations
/// of transform and inverse transform.
/// @relates Transform_
template <class P> inline Transform_<P>
operator*(const Transform_<P>& X1,        const Transform_<P>& X2)         {return X1.compose(X2);}
template <class P> inline Transform_<P>
operator*(const Transform_<P>& X1,        const InverseTransform_<P>& X2)  {return X1.compose(X2);}
template <class P> inline Transform_<P>
operator*(const InverseTransform_<P>& X1, const Transform_<P>& X2)         {return X1.compose(X2);}
template <class P> inline Transform_<P>
operator*(const InverseTransform_<P>& X1, const InverseTransform_<P>& X2)  {return X1.compose(X2);}

/// Comparison operators return true only if the two transforms are bit
/// identical; that's not too useful.
/// @relates Transform_
template <class P> inline bool
operator==(const Transform_<P>& X1,        const Transform_<P>& X2)         {return X1.R()==X2.R() && X1.p()==X2.p();}
template <class P> inline bool
operator==(const InverseTransform_<P>& X1, const InverseTransform_<P>& X2)  {return X1.R()==X2.R() && X1.p()==X2.p();}
template <class P> inline bool
operator==(const Transform_<P>& X1,        const InverseTransform_<P>& X2)  {return X1.R()==X2.R() && X1.p()==X2.p();}
template <class P> inline bool
operator==(const InverseTransform_<P>& X1, const Transform_<P>& X2)         {return X1.R()==X2.R() && X1.p()==X2.p();}

/// Generate formatted output of a Transform to an output stream.
/// @relates Transform_
template <class P> SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream&, const Transform_<P>&);
/// Generate formatted output of an InverseTransform to an output stream.
/// @relates InverseTransform_
template <class P> SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream&, const InverseTransform_<P>&);



//------------------------------------------------------------------------------
}  // End of namespace SimTK

//--------------------------------------------------------------------------
#endif // SimTK_SimTKCOMMON_TRANSFORM_H
//--------------------------------------------------------------------------

