#ifndef SimTK_SIMMATRIX_ROWVECTORBASE_H_
#define SimTK_SIMMATRIX_ROWVECTORBASE_H_

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
Define the SimTK::RowVectorBase class that is part of Simbody's BigMatrix 
toolset. **/

namespace SimTK {

//==============================================================================
//                              ROW VECTOR BASE
//==============================================================================
/** @brief This is a dataless rehash of the MatrixBase class to specialize it 
for RowVectors.

This mostly entails overriding a few of the methods. Note that all the 
MatrixBase operations remain available if you \c static_cast this up to a 
MatrixBase. **/
template <class ELT> class RowVectorBase : public MatrixBase<ELT> {
    typedef MatrixBase<ELT>                             Base;
    typedef typename CNT<ELT>::Scalar                   Scalar;
    typedef typename CNT<ELT>::Number                   Number;
    typedef typename CNT<ELT>::StdNumber                StdNumber;
    typedef RowVectorBase<ELT>                          T;
    typedef RowVectorBase<typename CNT<ELT>::TAbs>      TAbs;
    typedef RowVectorBase<typename CNT<ELT>::TNeg>      TNeg;
    typedef VectorView_<typename CNT<ELT>::THerm>       THerm;
public: 
    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase "owner" construction
    ///
    /// These constructors create new RowVectorBase objects which own their
    /// own data and are (at least by default) resizable. The resulting matrices
    /// are 1 x n with the number of rows locked at 1. If there is any data
    /// allocated but not explicitly initialized, that data will be uninitialized
    /// garbage in Release builds but will be initialized to NaN (at a performance
    /// cost) in Debug builds.
    /// @{

    /// Default constructor makes a 1x0 matrix locked at 1 row; you can
    /// provide an initial allocation if you want.
    explicit RowVectorBase(int n=0) : Base(MatrixCommitment::RowVector(), 1, n) {}
    
    /// Copy constructor is a deep copy (not appropriate for views!). That
    /// means it creates a new, densely packed vector whose elements are
    /// initialized from the source object.    
    RowVectorBase(const RowVectorBase& source) : Base(source) {}

    /// Implicit conversion from compatible row vector with negated elements.
    RowVectorBase(const TNeg& source) : Base(source) {}

    /// Construct an owner row vector of length n, with each element initialized to
    /// the given value.
    RowVectorBase(int n, const ELT& initialValue)
    :   Base(MatrixCommitment::RowVector(),1,n,initialValue) {}  

    /// Construct an owner vector of length n, with the elements initialized sequentially
    /// from a C++ array of elements which is assumed to be of length n. Note that we
    /// are expecting C++ packing; don't use this to initialize one Simmatrix vector
    /// from another because Simmatrix may pack its elements more densely than C++.
    RowVectorBase(int n, const ELT* cppInitialValues)
    :   Base(MatrixCommitment::RowVector(),1,n,cppInitialValues) {}
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase construction from pre-existing data
    ///
    /// Construct a non-resizeable, RowVectorBase view of externally supplied data. Note that
    /// stride should be interpreted as "the number of scalars between elements" and
    /// for composite elements may have a different value if the source is a C++ array 
    /// of elements vs. a Simmatrix packed data array. We provide constructors for
    /// both read-only and writable external data.
    /// @{

    /// Construct a read-only view of existing data.
    RowVectorBase(int n, int stride, const Scalar* s)
    :   Base(MatrixCommitment::RowVector(n), MatrixCharacter::RowVector(n),stride,s) { }
    /// Construct a writable view into existing data.
    RowVectorBase(int n, int stride, Scalar* s)
    :   Base(MatrixCommitment::RowVector(n), MatrixCharacter::RowVector(n),stride,s) { }
    /// @}

    //  ------------------------------------------------------------------------
    /// @name       RowVectorBase construction from an existing Helper.
    ///
    /// Create a new RowVectorBase from an existing helper. Both shallow (view) and deep 
    /// copies are possible. For shallow copies, there is a constructor providing a read-only
    /// view of the original data and one providing a writable view into the original data.
    /// @{

    /// Construct a writable view into the source data.
    RowVectorBase(MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::RowVector(), h,s) { }
    /// Construct a read-only view of the source data.
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::ShallowCopy& s) 
    :   Base(MatrixCommitment::RowVector(), h,s) { }
    /// Construct a new owner vector initialized with the data from the source.
    RowVectorBase(const MatrixHelper<Scalar>& h, const typename MatrixHelper<Scalar>::DeepCopy& d)    
    :   Base(MatrixCommitment::RowVector(), h,d) { }
    /// @}

    // This gives the resulting rowvector type when (r(i) op P) is applied to each element.
    // It will have element types which are the regular composite result of ELT op P.
    template <class P> struct EltResult { 
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Mul> Mul;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Dvd> Dvd;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Add> Add;
        typedef RowVectorBase<typename CNT<ELT>::template Result<P>::Sub> Sub;
    };

    /// Copy assignment is deep copy but behavior depends on type of lhs: if view, rhs
    /// must match. If owner, we reallocate and copy rhs.
    RowVectorBase& operator=(const RowVectorBase& b) {
        Base::operator=(b); return *this;
    }

    // default destructor

    RowVectorBase& operator*=(const StdNumber& t)     {Base::operator*=(t); return *this;}
    RowVectorBase& operator/=(const StdNumber& t)     {Base::operator/=(t); return *this;}
    RowVectorBase& operator+=(const RowVectorBase& r) {Base::operator+=(r); return *this;}
    RowVectorBase& operator-=(const RowVectorBase& r) {Base::operator-=(r); return *this;}  

    template <class EE> RowVectorBase& operator=(const RowVectorBase<EE>& b) 
      { Base::operator=(b);  return *this; } 
    template <class EE> RowVectorBase& operator+=(const RowVectorBase<EE>& b) 
      { Base::operator+=(b); return *this; } 
    template <class EE> RowVectorBase& operator-=(const RowVectorBase<EE>& b) 
      { Base::operator-=(b); return *this; } 

    // default destructor
 
    /// Fill current allocation with copies of element. Note that this is not the 
    /// same behavior as assignment for Matrices, where only the diagonal is set (and
    /// everything else is set to zero.)
    RowVectorBase& operator=(const ELT& t) { Base::setTo(t); return *this; } 

    /// There's only one row here so it's a bit wierd to use colScale rather than
    /// elementwiseMultiply, but there's nothing really wrong with it. Using rowScale
    /// would be really wacky since it is the same as a scalar multiply. We won't support
    /// rowScale here except through inheritance where it will not be much use.
    template <class EE> RowVectorBase& colScaleInPlace(const VectorBase<EE>& v)
    { Base::template colScaleInPlace<EE>(v); return *this; }
    template <class EE> inline void colScale(const VectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
    { return Base::template colScale<EE>(v,out); }
    template <class EE> inline typename EltResult<EE>::Mul colScale(const VectorBase<EE>& v) const
    { typename EltResult<EE>::Mul out(ncol()); Base::template colScale<EE>(v,out); return out; }


    // elementwise multiply
    template <class EE> RowVectorBase& elementwiseMultiplyInPlace(const RowVectorBase<EE>& r)
    { Base::template elementwiseMultiplyInPlace<EE>(r); return *this; }
    template <class EE> inline void elementwiseMultiply(const RowVectorBase<EE>& v, typename EltResult<EE>::Mul& out) const
    { Base::template elementwiseMultiply<EE>(v,out); }
    template <class EE> inline typename EltResult<EE>::Mul elementwiseMultiply(const RowVectorBase<EE>& v) const
    { typename EltResult<EE>::Mul out(nrow()); Base::template elementwiseMultiply<EE>(v,out); return out; }

    // elementwise multiply from left
    template <class EE> RowVectorBase& elementwiseMultiplyFromLeftInPlace(const RowVectorBase<EE>& r)
    { Base::template elementwiseMultiplyFromLeftInPlace<EE>(r); return *this; }
    template <class EE> inline void 
    elementwiseMultiplyFromLeft(
        const RowVectorBase<EE>& v, 
        typename RowVectorBase<EE>::template EltResult<ELT>::Mul& out) const
    { 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out);
    }
    template <class EE> inline 
    typename RowVectorBase<EE>::template EltResult<ELT>::Mul 
    elementwiseMultiplyFromLeft(const RowVectorBase<EE>& v) const {
        typename RowVectorBase<EE>::template EltResult<ELT>::Mul out(nrow()); 
        Base::template elementwiseMultiplyFromLeft<EE>(v,out); 
        return out;
    }

    // elementwise divide
    template <class EE> RowVectorBase& elementwiseDivideInPlace(const RowVectorBase<EE>& r)
    { Base::template elementwiseDivideInPlace<EE>(r); return *this; }
    template <class EE> inline void elementwiseDivide(const RowVectorBase<EE>& v, typename EltResult<EE>::Dvd& out) const
    { Base::template elementwiseDivide<EE>(v,out); }
    template <class EE> inline typename EltResult<EE>::Dvd elementwiseDivide(const RowVectorBase<EE>& v) const
    { typename EltResult<EE>::Dvd out(nrow()); Base::template elementwiseDivide<EE>(v,out); return out; }

    // elementwise divide from left
    template <class EE> RowVectorBase& elementwiseDivideFromLeftInPlace(const RowVectorBase<EE>& r)
    { Base::template elementwiseDivideFromLeftInPlace<EE>(r); return *this; }
    template <class EE> inline void 
    elementwiseDivideFromLeft
       (const RowVectorBase<EE>& v, 
        typename RowVectorBase<EE>::template EltResult<ELT>::Dvd& out) const { 
        Base::template elementwiseDivideFromLeft<EE>(v,out);
    }
    template <class EE> inline 
    typename RowVectorBase<EE>::template EltResult<ELT>::Dvd 
    elementwiseDivideFromLeft(const RowVectorBase<EE>& v) const    { 
        typename RowVectorBase<EE>::template EltResult<ELT>::Dvd out(nrow()); 
        Base::template elementwiseDivideFromLeft<EE>(v,out); 
        return out;
    }

    // Implicit conversions are allowed to RowVector or Matrix, but not to Vector.   
    operator const RowVector_<ELT>&()     const {return *reinterpret_cast<const RowVector_<ELT>*>(this);}
    operator       RowVector_<ELT>&()           {return *reinterpret_cast<      RowVector_<ELT>*>(this);}
    operator const RowVectorView_<ELT>&() const {return *reinterpret_cast<const RowVectorView_<ELT>*>(this);}
    operator       RowVectorView_<ELT>&()       {return *reinterpret_cast<      RowVectorView_<ELT>*>(this);}
    
    operator const Matrix_<ELT>&()     const {return *reinterpret_cast<const Matrix_<ELT>*>(this);}
    operator       Matrix_<ELT>&()           {return *reinterpret_cast<      Matrix_<ELT>*>(this);} 
    operator const MatrixView_<ELT>&() const {return *reinterpret_cast<const MatrixView_<ELT>*>(this);}
    operator       MatrixView_<ELT>&()       {return *reinterpret_cast<      MatrixView_<ELT>*>(this);} 
    

    // size() for RowVectors is Base::nelt() but returns int instead of ptrdiff_t.
    int size() const { 
        assert(Base::nelt() <= (ptrdiff_t)std::numeric_limits<int>::max()); 
        assert(Base::nrow()==1);
        return (int)Base::nelt();
    }
    int       nrow() const {assert(Base::nrow()==1); return Base::nrow();}
    int       ncol() const {assert(Base::nrow()==1); return Base::ncol();}
    ptrdiff_t nelt() const {assert(Base::nrow()==1); return Base::nelt();}

    // Override MatrixBase operators to return the right shape
    TAbs abs() const {
        TAbs result; Base::abs(result); return result;
    }

    // Override MatrixBase indexing operators          
    const ELT& operator[](int j) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(j));}
    ELT&       operator[](int j)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(j));}
    const ELT& operator()(int j) const {return *reinterpret_cast<const ELT*>(Base::getHelper().getElt(j));}
    ELT&       operator()(int j)       {return *reinterpret_cast<ELT*>      (Base::updHelper().updElt(j));}
         
    // Block (contiguous subvector) creation      
    RowVectorView_<ELT> operator()(int j, int n) const {return Base::operator()(0,j,1,n).getAsRowVectorView();}
    RowVectorView_<ELT> operator()(int j, int n)       {return Base::operator()(0,j,1,n).updAsRowVectorView();}

    // Indexed view creation (arbitrary subvector). Indices must be monotonically increasing.
    RowVectorView_<ELT> index(const Array_<int>& indices) const {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::getHelper(), indices);
        return RowVectorView_<ELT>(h);
    }
    RowVectorView_<ELT> updIndex(const Array_<int>& indices) {
        MatrixHelper<Scalar> h(Base::getHelper().getCharacterCommitment(), Base::updHelper(), indices);
        return RowVectorView_<ELT>(h);
    }

    RowVectorView_<ELT> operator()(const Array_<int>& indices) const {return index(indices);}
    RowVectorView_<ELT> operator()(const Array_<int>& indices)       {return updIndex(indices);}
 
    // Hermitian transpose.
    THerm transpose() const {return Base::transpose().getAsVectorView();}
    THerm updTranspose()    {return Base::updTranspose().updAsVectorView();}

    THerm operator~() const {return transpose();}
    THerm operator~()       {return updTranspose();}

    const RowVectorBase& operator+() const {return *this; }

    // Negation

    const TNeg& negate()    const {return *reinterpret_cast<const TNeg*>(this); }
    TNeg&       updNegate()       {return *reinterpret_cast<TNeg*>(this); }

    const TNeg& operator-() const {return negate();}
    TNeg&       operator-()       {return updNegate();}

    RowVectorBase& resize(int n)     {Base::resize(1,n); return *this;}
    RowVectorBase& resizeKeep(int n) {Base::resizeKeep(1,n); return *this;}

    //TODO: this is not re-locking the number of rows at 1.
    void clear() {Base::clear(); Base::resize(1,0);}

    ELT sum() const {ELT s; Base::getHelper().sum(reinterpret_cast<Scalar*>(&s)); return s; } // add all the elements        
    VectorIterator<ELT, RowVectorBase<ELT> > begin() {
        return VectorIterator<ELT, RowVectorBase<ELT> >(*this, 0);
    }
    VectorIterator<ELT, RowVectorBase<ELT> > end() {
        return VectorIterator<ELT, RowVectorBase<ELT> >(*this, size());
    }

protected:
    // Create a RowVectorBase handle using a given helper rep. 
    explicit RowVectorBase(MatrixHelperRep<Scalar>* hrep) : Base(hrep) {}

private:
    // NO DATA MEMBERS ALLOWED
};

} //namespace SimTK

#endif // SimTK_SIMMATRIX_ROWVECTORBASE_H_
