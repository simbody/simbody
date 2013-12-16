#ifndef SimTK_SIMMATRIX_VECTORBASE_H_
#define SimTK_SIMMATRIX_VECTORBASE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/** @file
Define the SimTK::VectorBase class that is part of Simbody's BigMatrix 
toolset. **/

namespace SimTK {

//==============================================================================
//                                VECTOR BASE
//==============================================================================
/** @brief This is a dataless rehash of the MatrixBase class to specialize it 
for Vectors.

This mostly entails overriding a few of the methods. Note that all the 
MatrixBase operations remain available if you \c static_cast this up to a 
MatrixBase. **/
template <class ELT> class VectorBase : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                             Base;
    typedef typename Base::ScalarNormSq                 ScalarNormSq;
    typedef typename Base::EAbs                         EAbs;
    typedef typename CNT<ELT>::Scalar                   Scalar;
    typedef typename CNT<ELT>::Number                   Number;
    typedef typename CNT<ELT>::StdNumber                StdNumber;
    typedef VectorBase<ELT>                             T;
    typedef VectorBase<typename CNT<ELT>::TAbs>         TAbs;
    typedef VectorBase<typename CNT<ELT>::TNeg>         TNeg;
    typedef RowVectorView_<typename CNT<ELT>::THerm>    THerm;
public:  
    //  ------------------------------------------------------------------------
    /// @name       VectorBase "owner" construction
    ///
    /// These constructors create new VectorBase objects which own their
    /// own data and are (at least by default) resizable. The resulting matrices
    /// are m X 1 with the number of columns locked at 1. If there is any data
    /// allocated but not explicitly initialized, that data will be uninitialized
    /// garbage in Release builds but will be initialized to NaN (at a performance
    /// cost) in Debug builds.
    /// @{

    /// Default constructor makes a 0x1 matrix locked at 1 column; you can
    /// provide an initial allocation if you want.
    explicit VectorBase(int m=0) : Base(MatrixCommitment::Vector(), m, 1) {}

    /// Copy constructor is a deep copy (not appropriate for views!). That
    /// means it creates a new, densely packed vector whose elements are
    /// initialized from the source object.
    VectorBase(const VectorBase& source) : Base(source) {}

    /// Implicit conversion from compatible vector with negated elements.
    VectorBase(const TNeg& source) : Base(source) {}

    /// Construct an owner vector of length m, with each element initialized to
    /// the given value.
    VectorBase(int m, const ELT& initialValue)
    :   Base(MatrixCommitment::Vector(),m,1,initialValue) {}  

    /// Construct an owner vector of length m, with the elements initialized sequentially
    /// from a C++ array of elements which is assumed to be of length m. Note that we
    /// are expecting C++ packing; don't use this to initialize one Simmatrix vector
    /// from another because Simmatrix may pack its elements more densely than C++.
    VectorBase(int m, const ELT* cppInitialValues)
    :   Base(MatrixCommitment::Vector(),m,1,cppInitialValues) {}
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       VectorBase construction from pre-existing data
    ///
    /// Construct a non-resizeable, VectorBase view of externally supplied data. Note that
    /// stride should be interpreted as "the number of scalars between elements" and
    /// for composite elements may have a different value if the source is a C++ array 
    /// of elements vs. a Simmatrix packed data array. We provide constructors for
    /// both read-only and writable external data.
    /// @{

    /// Construct a read-only view of existing data.
    VectorBase(int m, int stride, const Scalar* s)
    :   Base(MatrixCommitment::Vector(m), MatrixCharacter::Vector(m),stride,s) { }
    /// Construct a writable view into existing data.
    VectorBase(int m, int stride, Scalar* s)
    :   Base(MatrixCommitment::Vector(m), MatrixCharacter::Vector(m),stride,s) { }
    /// @}
        
    //  ------------------------------------------------------------------------
    /// @name       VectorBase construction from an existing Helper.
    ///
    /// Create a new VectorBase from an existing helper. Both shallow (view) and deep 
    /// copies are possible. For shallow copies, there is a constructor providing a read-only
    /// view of the original data and one providing a writable view into the original data.
    /// @{

    /// Construct a writable view into the source data.
    VectorBase(MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::Vector(), h,s) { }
    /// Construct a read-only view of the source data.
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::Vector(), h,s) { }
    /// Construct a new owner vector initialized with the data from the source.
    VectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    
    :   Base(MatrixCommitment::Vector(), h,d) { }
    /// @}

    // This gives the resulting vector type when (v[i] op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef VectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

    /// Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    /// must match. If owner, we reallocate and copy rhs.
    VectorBase& operator=(const VectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor


    VectorBase& operator*=(const StdNumber& t)  { Base::operator*=(t); return *this; }
    VectorBase& operator/=(const StdNumber& t)  { Base::operator/=(t); return *this; }
    VectorBase& operator+=(const VectorBase& r) { Base::operator+=(r); return *this; }
    VectorBase& operator-=(const VectorBase& r) { Base::operator-=(r); return *this; }  


    template <class EE> VectorBase& operator=(const VectorBase<EE>& b) 
      { Base::operator=(b);  return *this; } 
    template <class EE> VectorBase& operator+=(const VectorBase<EE>& b) 
      { Base::operator+=(b); return *this; } 
    template <class EE> VectorBase& operator-=(const VectorBase<EE>& b) 
      { Base::operator-=(b); return *this; } 


    /// Fill current allocation with copies of element. Note that this is not the 
    /// same behavior as assignment for Matrices, where only the diagonal is set (and
    /// everything else is set to zero.)
    VectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; }  

    /// There's only one column here so it's a bit weird to use rowScale rather than
    /// elementwiseMultiply, but there's nothing really wrong with it. Using colScale
    /// would be really wacky since it is the same as a scalar multiply. We won't support
    /// colScale here except through inheritance where it will not be much use.
    template <class EE> VectorBase& rowScaleInPlace(const VectorBase<EE>& v)
	{ Base::template rowScaleInPlace<EE>(v); return *this; }
    template <class EE> inline void rowScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	{ Base::rowScale(v,out); }
    template <class EE> inline typename EltResult<EE>::Mul rowScale(const VectorBase<EE>& v) const
	{ typename EltResult<EE>::Mul out(nrow()); Base::rowScale(v,out); return out; }

    /** Return the root-mean-square (RMS) norm of a Vector of scalars, with 
    optional return of the index of the element of largest absolute value. 
    The RMS norm of a Vector v of length n is rms=sqrt(~v*v/n). If n==0 we
    define the RMS norm to be zero but return the element index as -1. **/
    typename CNT<ScalarNormSq>::TSqrt 
    normRMS(int* worstOne=0) const {
        if (!CNT<ELT>::IsScalar)
            SimTK_THROW1(Exception::Cant, 
                "Vector::normRMS() only defined for scalar elements.");
        const int n = nrow();
        if (n == 0) {
            if (worstOne) *worstOne = -1;
            return typename CNT<ScalarNormSq>::TSqrt(0);
        }

        ScalarNormSq sumsq = 0;
        if (worstOne) {
            *worstOne = 0;
            ScalarNormSq maxsq = 0; 
            for (int i=0; i<n; ++i) {
                const ScalarNormSq v2 = square((*this)[i]);
                if (v2 > maxsq) maxsq=v2, *worstOne=i;
                sumsq += v2;
            }
        } else { // don't track the worst element
            for (int i=0; i<n; ++i) {
                const ScalarNormSq v2 = square((*this)[i]);
                sumsq += v2;
            }
        }

        return CNT<ScalarNormSq>::sqrt(sumsq/n);
    }

    /** Return the weighted root-mean-square (WRMS) norm of a Vector of 
    scalars, with optional return of the index of the weighted element of 
    largest absolute value. The WRMS norm of a Vector v of length n with
    weights w is wrms=sqrt(sum_i((w_i*v_i)^2))/n). If n==0 we
    define the WRMS norm to be zero but return the element index as -1. **/
    template <class EE>
    typename CNT<ScalarNormSq>::TSqrt 
    weightedNormRMS(const VectorBase<EE>& w, int* worstOne=0) const {
        if (!CNT<ELT>::IsScalar || !CNT<EE>::IsScalar)
            SimTK_THROW1(Exception::Cant, 
            "Vector::weightedNormRMS() only defined for scalar elements"
            " and weights.");
        const int n = nrow();
        assert(w.nrow()==n);
        if (n == 0) {
            if (worstOne) *worstOne = -1;
            return typename CNT<ScalarNormSq>::TSqrt(0);
        }

        ScalarNormSq sumsq = 0;
        if (worstOne) {
            *worstOne = 0;
            ScalarNormSq maxsq = 0; 
            for (int i=0; i<n; ++i) {
                const ScalarNormSq wv2 = square(w[i]*(*this)[i]);
                if (wv2 > maxsq) maxsq=wv2, *worstOne=i;
                sumsq += wv2;
            }
        } else { // don't track the worst element
            for (int i=0; i<n; ++i) {
                const ScalarNormSq wv2 = square(w[i]*(*this)[i]);
                sumsq += wv2;
            }
        }

        return CNT<ScalarNormSq>::sqrt(sumsq/n);
    }

    /** Return the infinity norm (max absolute value) of a Vector of scalars, 
    with optional return of the index of the element of largest absolute value. 
    The Inf norm of a Vector v is inf=max_i(|v_i|). If n==0 we
    define the Inf norm to be zero but return the element index as -1. **/
    EAbs normInf(int* worstOne=0) const {
        if (!CNT<ELT>::IsScalar)
            SimTK_THROW1(Exception::Cant, 
                "Vector::normInf() only defined for scalar elements.");
        const int n = nrow();
        if (n == 0) {
            if (worstOne) *worstOne = -1;
            return EAbs(0);
        }

        EAbs maxabs = 0;
        if (worstOne) {
            *worstOne = 0;
            for (int i=0; i<n; ++i) {
                const EAbs a = std::abs((*this)[i]);
                if (a > maxabs) maxabs=a, *worstOne=i;
            }
        } else { // don't track the worst element
            for (int i=0; i<n; ++i) {
                const EAbs a = std::abs((*this)[i]);
                if (a > maxabs) maxabs=a;
            }
        }

        return maxabs;
    }

    /** Return the weighted infinity norm (max absolute value) WInf of a Vector
    of scalars, with optional return of the index of the weighted element of 
    largest absolute value. The WInf norm of a Vector v of length n with
    weights w is winf=max_i(|w_i*v_i|). If n==0 we
    define the WInf norm to be zero but return the element index as -1. **/
    template <class EE>
    EAbs weightedNormInf(const VectorBase<EE>& w, int* worstOne=0) const {
        if (!CNT<ELT>::IsScalar || !CNT<EE>::IsScalar)
            SimTK_THROW1(Exception::Cant, 
            "Vector::weightedNormInf() only defined for scalar elements"
            " and weights.");
        const int n = nrow();
        assert(w.nrow()==n);
        if (n == 0) {
            if (worstOne) *worstOne = -1;
            return EAbs(0);
        }

        EAbs maxabs = 0;
        if (worstOne) {
            *worstOne = 0;
            for (int i=0; i<n; ++i) {
                const EAbs wv = std::abs(w[i]*(*this)[i]);
                if (wv > maxabs) maxabs=wv, *worstOne=i;
            }
        } else { // don't track the worst element
            for (int i=0; i<n; ++i) {
                const EAbs wv = std::abs(w[i]*(*this)[i]);
                if (wv > maxabs) maxabs=wv;
            }
        }

        return maxabs;
    }

    /// Set this[i] = this[i]^-1.
    VectorBase& elementwiseInvertInPlace() {
        Base::elementwiseInvertInPlace();
        return *this;
    }

    /// Set supplied out[i] = this[i]^-1
    void elementwiseInvert(VectorBase<typename CNT<ELT>::TInvert>& out) const {
        Base::elementwiseInvert(out);
    }

    /// Return out[i]=this[i]^-1 as function return.
    VectorBase<typename CNT<ELT>::TInvert> elementwiseInvert() const {
        VectorBase<typename CNT<ELT>::TInvert> out(nrow());
        Base::elementwiseInvert(out);
        return out;
    }

    // elementwise multiply
    template <class EE> VectorBase& elementwiseMultiplyInPlace(const VectorBase<EE>& r)
	{ Base::template elementwiseMultiplyInPlace<EE>(r); return *this; }
    template <class EE> inline void elementwiseMultiply(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
	{ Base::template elementwiseMultiply<EE>(v,out); }
    template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const VectorBase<EE>& v) const
	{ typename EltResult<EE>::Mul out(nrow()); Base::template elementwiseMultiply<EE>(v,out); return out; }

    // elementwise multiply from left
    template <class EE> VectorBase& elementwiseMultiplyFromLeftInPlace(const VectorBase<EE>& r)
	{ Base::template elementwiseMultiplyFromLeftInPlace<EE>(r); return *this; }
    template <class EE> inline void 
    elementwiseMultiplyFromLeft(
        const VectorBase<EE>& v, 
        typename VectorBase<EE>::template EltResult<ELT>::Mul& out) const
	{ 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out);
    }
	template <class EE> inline typename VectorBase<EE>::template EltResult<ELT>::Mul 
    elementwiseMultiplyFromLeft(const VectorBase<EE>& v) const
	{ 
        typename VectorBase<EE>::template EltResult<ELT>::Mul out(nrow()); 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out); 
        return out;
    }

    // elementwise divide
    template <class EE> VectorBase& elementwiseDivideInPlace(const VectorBase<EE>& r)
	{ Base::template elementwiseDivideInPlace<EE>(r); return *this; }
    template <class EE> inline void elementwiseDivide(const VectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
	{ Base::template elementwiseDivide<EE>(v,out); }
    template <class EE> inline typename EltResult<EE>::Dvd elementwiseDivide(const VectorBase<EE>& v) const
	{ typename EltResult<EE>::Dvd out(nrow()); Base::template elementwiseDivide<EE>(v,out); return out; }

    // elementwise divide from left
    template <class EE> VectorBase& elementwiseDivideFromLeftInPlace(const VectorBase<EE>& r)
	{ Base::template elementwiseDivideFromLeftInPlace<EE>(r); return *this; }
    template <class EE> inline void 
    elementwiseDivideFromLeft(
        const VectorBase<EE>& v, 
        typename VectorBase<EE>::template EltResult<ELT>::Dvd& out) const
	{ 
        Base::template elementwiseDivideFromLeft<EE>(v,out);
    }
	template <class EE> inline typename VectorBase<EE>::template EltResult<ELT>::Dvd 
    elementwiseDivideFromLeft(const VectorBase<EE>& v) const
	{ 
        typename VectorBase<EE>::template EltResult<ELT>::Dvd out(nrow()); 
        Base::template elementwiseDivideFromLeft<EE>(v,out); 
        return out;
    }

    // Implicit conversions are allowed to Vector or Matrix, but not to RowVector.   
    operator const Vector_<ELT>&()     const { return *reinterpret_cast<const Vector_<ELT>*>(this); }
    operator       Vector_<ELT>&()           { return *reinterpret_cast<      Vector_<ELT>*>(this); }
    operator const VectorView_<ELT>&() const { return *reinterpret_cast<const VectorView_<ELT>*>(this); }
    operator       VectorView_<ELT>&()       { return *reinterpret_cast<      VectorView_<ELT>*>(this); }
    
    operator const Matrix_<ELT>&()     const { return *reinterpret_cast<const Matrix_<ELT>*>(this); }
    operator       Matrix_<ELT>&()           { return *reinterpret_cast<      Matrix_<ELT>*>(this); } 
    operator const MatrixView_<ELT>&() const { return *reinterpret_cast<const MatrixView_<ELT>*>(this); }
    operator       MatrixView_<ELT>&()       { return *reinterpret_cast<      MatrixView_<ELT>*>(this); } 


    // size() for Vectors is Base::nelt() but returns int instead of ptrdiff_t.
    int size() const { 
	assert(Base::nelt() <= (ptrdiff_t)std::numeric_limits<int>::max()); 
	assert(Base::ncol()==1);
	return (int)Base::nelt();
    }
    int       nrow() const {assert(Base::ncol()==1); return Base::nrow();}
    int       ncol() const {assert(Base::ncol()==1); return Base::ncol();}
    ptrdiff_t nelt() const {assert(Base::ncol()==1); return Base::nelt();}

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {TAbs result; Base::abs(result); return result;}
    
    // Override MatrixBase indexing operators          
    const ELT& operator[](int i) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(i));}
    ELT&       operator[](int i)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(i));}
    const ELT& operator()(int i) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(i));}
    ELT&       operator()(int i)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(i));}
         
    // Block (contiguous subvector) view creation      
    VectorView_<ELT> operator()(int i, int m) const {return Base::operator()(i,0,m,1).getAsVectorView();}
    VectorView_<ELT> operator()(int i, int m)       {return Base::operator()(i,0,m,1).updAsVectorView();}

    // Indexed view creation (arbitrary subvector). Indices must be 
    // monotonically increasing.
    VectorView_<ELT> index(const Array_<int>& indices) const {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), 
                               Base::getHelper(), indices);
        return VectorView_<ELT>(h);
    }
    VectorView_<ELT> updIndex(const Array_<int>& indices) {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), 
                               Base::updHelper(), indices);
        return VectorView_<ELT>(h);
    }

    VectorView_<ELT> operator()(const Array_<int>& indices) const {return index(indices);}
    VectorView_<ELT> operator()(const Array_<int>& indices)       {return updIndex(indices);}
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsRowVectorView();}
    THerm updTranspose()    {return Base::updTranspose().updAsRowVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const VectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    VectorBase& resize(int m)     {Base::resize(m,1); return *this;}
    VectorBase& resizeKeep(int m) {Base::resizeKeep(m,1); return *this;}

    //TODO: this is not re-locking the number of columns at 1.
    void clear() {Base::clear(); Base::resize(0,1);}

    ELT sum() const {ELT s; Base::getHelper().sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        
    VectorIterator<ELT, VectorBase<ELT> > begin() {
        return VectorIterator<ELT, VectorBase<ELT> >(*this, 0);
    }
    VectorIterator<ELT, VectorBase<ELT> > end() {
        return VectorIterator<ELT, VectorBase<ELT> >(*this, size());
    }

protected:
    // Create a VectorBase handle using a given helper rep. 
    explicit VectorBase(MatrixHelperRep<Scalar>* hrep) : Base(hrep) {}

private:
    // NO DATA MEMBERS ALLOWED
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_VECTORBASE_H_
