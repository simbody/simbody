#ifndef SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_
#define SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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

/**@file
 * Declaration of class Vec<NROWS, ELEMENT_TYPE, STRIDE>.
 */

#include "SimTKcommon/internal/common.h"

namespace SimTK {


// The following functions are used internally by Vec.

// Hide from Doxygen.
/** @cond **/
namespace Impl {

// For those wimpy compilers that don't unroll short, constant-limit loops, 
// Peter Eastman added these recursive template implementations of 
// elementwise add, subtract, and copy. Sherm added multiply and divide.

template <class E1, int S1, class E2, int S2> void
conformingAdd(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, 
              Vec<1,typename CNT<E1>::template Result<E2>::Add>& result) {
    result[0] = r1[0] + r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingAdd(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, 
              Vec<N,typename CNT<E1>::template Result<E2>::Add>& result) {
    conformingAdd(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), 
                  reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), 
                  reinterpret_cast<Vec<N-1,typename CNT<E1>::
                              template Result<E2>::Add>&>(result));
    result[N-1] = r1[N-1] + r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
conformingSubtract(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, 
                   Vec<1,typename CNT<E1>::template Result<E2>::Sub>& result) {
    result[0] = r1[0] - r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
conformingSubtract(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, 
                   Vec<N,typename CNT<E1>::template Result<E2>::Sub>& result) {
    conformingSubtract(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), 
                       reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), 
                       reinterpret_cast<Vec<N-1,typename CNT<E1>::
                                   template Result<E2>::Sub>&>(result));
    result[N-1] = r1[N-1] - r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
elementwiseMultiply(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, 
              Vec<1,typename CNT<E1>::template Result<E2>::Mul>& result) {
    result[0] = r1[0] * r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
elementwiseMultiply(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, 
              Vec<N,typename CNT<E1>::template Result<E2>::Mul>& result) {
    elementwiseMultiply(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), 
                        reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), 
                        reinterpret_cast<Vec<N-1,typename CNT<E1>::
                                    template Result<E2>::Mul>&>(result));
    result[N-1] = r1[N-1] * r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
elementwiseDivide(const Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2, 
              Vec<1,typename CNT<E1>::template Result<E2>::Dvd>& result) {
    result[0] = r1[0] / r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
elementwiseDivide(const Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2, 
              Vec<N,typename CNT<E1>::template Result<E2>::Dvd>& result) {
    elementwiseDivide(reinterpret_cast<const Vec<N-1,E1,S1>&>(r1), 
                      reinterpret_cast<const Vec<N-1,E2,S2>&>(r2), 
                      reinterpret_cast<Vec<N-1,typename CNT<E1>::
                                  template Result<E2>::Dvd>&>(result));
    result[N-1] = r1[N-1] / r2[N-1];
}

template <class E1, int S1, class E2, int S2> void
copy(Vec<1,E1,S1>& r1, const Vec<1,E2,S2>& r2) {
    r1[0] = r2[0];
}
template <int N, class E1, int S1, class E2, int S2> void
copy(Vec<N,E1,S1>& r1, const Vec<N,E2,S2>& r2) {
    copy(reinterpret_cast<Vec<N-1,E1,S1>&>(r1), 
         reinterpret_cast<const Vec<N-1,E2,S2>&>(r2));
    r1[N-1] = r2[N-1];
}

}
/** @endcond **/

/** This is a fixed-length column vector designed for no-overhead inline 
computation.

@ingroup MatVecUtilities

@tparam     M       The number of rows in the vector.
@tparam     ELT     The element type. Must be a composite numerical type (CNT).
                    The default is ELT=Real.
@tparam     STRIDE  The spacing from one element to the next in memory, as an 
                    integer number of elements of type ELT. The default is 
                    STRIDE=1.

<b>Usage</b>

The %Vec and @ref SimTK::Vector_ "Vector" classes are commonly used to represent 
tuples of Real values, and have methods like %norm() to calculate the vector 
2-norm. Use %Vec for a small vector whose length is known at compile time; 
otherwise, use @ref SimTK::Vector_ "Vector". To collect elements of the same 
type that do not constitute a tuple, it is more appropriate to use the Array_ 
container. Some common %Vec use cases are provided below.

<b>Construction</b>

The 3-tuple <tt>~[0,0,0]</tt> can be created in the following ways:
\code
Vec<3,Real>(0,0,0);
Vec3(0,0,0);
Vec3(0);
\endcode
Note that the default element type is Real, and that Vec3 is a typedef for
%Vec<3>; analogous typedefs exist for vectors of up to 9 elements.

<b>Manipulation</b>

Standard arithmetic operators can be used, as well as methods like %sum() and
%normalize(). Here are some usage examples, each of which returns
<tt>~[1,2,3]</tt>:
\code
Vec9(0,0,0,0,0,1,2,3,0).getSubVec<3>(5);
Vec4(1,2,4,3).drop1(2);
Vec2(2,3).insert1(0,1);
Vec2(1,2).append1(3);
\endcode

<b>Conversion</b>

It may be necessary to convert between a %Vec and a Vector (to interface with
%FactorQTZ, for instance). In the example below, we print a Vec3 created from a
3-element Vector:
\code
Vector myVector(3);
for (int i=0; i<myVector.size(); ++i) myVector[i]=i+1;
std::cout << Vec3(&myVector[0]) << std::endl;
\endcode
Converting from a Vec3 to a Vector is also straightforward:
\code
Vec3 myVec3(1,2,3);
std::cout << Vector(myVec3) << std::endl;
\endcode

@see Vector_ for handling of large or variable-sized vectors.
@see Array_, Mat, Matrix_
**/
template <int M, class ELT, int STRIDE>
class Vec {
public:
    /** @name Advanced 
    These are obscure members of %Vec that are used for template metaprogramming
    and can be ignored by most users. **/
    /**@{**/
    /** Element type of this %Vec. **/
    typedef ELT                                 E;
    /** Negated version of this %Vec's element type; ENeg==negator< E >. **/
    typedef typename CNT<E>::TNeg               ENeg;
    /** Element type, stripped of negator<> if it has one. **/
    typedef typename CNT<E>::TWithoutNegator    EWithoutNegator;
    /** Type showing just the real part of an element of this %Vec if elements
    are complex; otherwise just the element type. **/
    typedef typename CNT<E>::TReal              EReal;
    /** Type showing the imaginary part of an element of this %Vec as real,
    if elements are complex; otherwise a type that can hold a zero of the 
    element type. **/
    typedef typename CNT<E>::TImag              EImag;
    /** Type that elements would have if complex, if E is currently real;
    otherwise just the element type E. **/
    typedef typename CNT<E>::TComplex           EComplex;
    /** Type of the Hermitian transpose of an element of this %Vec. **/
    typedef typename CNT<E>::THerm              EHerm;
    /** Type of a \e positional transpose of an element of this %Vec. **/
    typedef typename CNT<E>::TPosTrans          EPosTrans;
    /** Type of the expression ~E*E (default vector and matrix square;
    symmetric). **/
    typedef typename CNT<E>::TSqHermT           ESqHermT;
    /** Type of the expression E*~E ("row square"; symmetric). **/
    typedef typename CNT<E>::TSqTHerm           ESqTHerm;
    /** Type required to hold the result of sqrt(E). **/
    typedef typename CNT<E>::TSqrt              ESqrt;
    /** Type required to hold the result of abs(E). **/
    typedef typename CNT<E>::TAbs               EAbs;
    /** Return type of standardize(E) method; a packed type that can hold the
    value of an element after eliminating negator and conjugate types. **/
    typedef typename CNT<E>::TStandard          EStandard;
    /** Packed type that can hold the value returned from invert(E), the
    inverse type of an element. **/
    typedef typename CNT<E>::TInvert            EInvert;
    /** Packed type that can hold the value returned from normalize(E). **/
    typedef typename CNT<E>::TNormalize         ENormalize;

    typedef typename CNT<E>::Scalar             EScalar;
    typedef typename CNT<E>::ULessScalar        EULessScalar;
    typedef typename CNT<E>::Number             ENumber;
    typedef typename CNT<E>::StdNumber          EStdNumber;
    typedef typename CNT<E>::Precision          EPrecision;
    typedef typename CNT<E>::ScalarNormSq       EScalarNormSq;

    /** These compile-time constants are required of every Composite
    Numerical Type (CNT). **/
    #ifndef SWIG
    enum {
        NRows               = M,
        NCols               = 1,
        NPackedElements     = M,
        NActualElements     = M * STRIDE,   // includes trailing gap
        NActualScalars      = CNT<E>::NActualScalars * NActualElements,
        RowSpacing          = STRIDE,
        ColSpacing          = NActualElements,
        ImagOffset          = NTraits<ENumber>::ImagOffset,
        RealStrideFactor    = 1, // composite types don't change size when
                                 // cast from complex to real or imaginary
        ArgDepth            = ((int)CNT<E>::ArgDepth < (int)MAX_RESOLVED_DEPTH 
                                ? CNT<E>::ArgDepth + 1 
                                : MAX_RESOLVED_DEPTH),
        IsScalar            = 0,
        IsULessScalar       = 0,
        IsNumber            = 0,
        IsStdNumber         = 0,
        IsPrecision         = 0,
        SignInterpretation  = CNT<E>::SignInterpretation
    };
    #endif

    // These are reinterpretations of the current data, so have the
    // same packing (stride).

    /** The type of this Vec. **/
    typedef Vec<M,E,STRIDE>                 T;
    /** Type this Vec would have if its elements were interpreted as
    negated. **/
    typedef Vec<M,ENeg,STRIDE>              TNeg;
    /** Type of this Vec with negator removed from its element type, if
    the element is negated. **/
    typedef Vec<M,EWithoutNegator,STRIDE>   TWithoutNegator;
    /** Type of this Vec cast to show only the real part of its element;
    this might affect the stride. **/
    typedef Vec<M,EReal,STRIDE*CNT<E>::RealStrideFactor>         
                                            TReal;
    /** Type of this Vec cast to show only the imaginary part of its element;
    this might affect the stride. **/
    typedef Vec<M,EImag,STRIDE*CNT<E>::RealStrideFactor>         
                                            TImag;
    typedef Vec<M,EComplex,STRIDE>          TComplex;
    /** Type of this Vec after casting to its Hermitian transpose; that is,
    the Vec turns into a Row and each element turns into \e its Hermitian
    transpose. **/
    typedef Row<M,EHerm,STRIDE>             THerm;
    /** Type of this Vec after casting to its positional transpose; that is,
    the Vec turns into a Row but the element type remains unchanged. **/
    typedef Row<M,E,STRIDE>                 TPosTrans;
    /** Element type of this Vec. **/
    typedef E                               TElement;
    /** Type of a row of this CNT object (for a Vec, just its element type). **/
    typedef E                               TRow;
    /** Type of a column of this CNT object (for a Vec, the whole thing). **/
    typedef Vec                             TCol;

    // These are the results of calculations, so are returned in new, packed
    // memory. Be sure to refer to element types here which are also packed.
    typedef Vec<M,ESqrt,1>                  TSqrt;      // Note stride
    typedef Vec<M,EAbs,1>                   TAbs;       // Note stride
    typedef Vec<M,EStandard,1>              TStandard;
    typedef Row<M,EInvert,1>                TInvert;
    typedef Vec<M,ENormalize,1>             TNormalize;

    typedef ESqHermT                        TSqHermT;   // result of self dot product
    typedef SymMat<M,ESqTHerm>              TSqTHerm;   // result of self outer product

    // These recurse right down to the underlying scalar type no matter how
    // deep the elements are.
    typedef EScalar                         Scalar;
    typedef EULessScalar                    ULessScalar;
    typedef ENumber                         Number;
    typedef EStdNumber                      StdNumber;
    typedef EPrecision                      Precision;
    typedef EScalarNormSq                   ScalarNormSq;
    /**@}**/

    /** The number of elements in this Vec (note that stride does not 
    affect this number.) **/
    static int size() { return M; }
    /** The number of rows in a Vec is the number of elements. **/
    static int nrow() { return M; }
    /** The number of columns in a Vec is always 1. **/
    static int ncol() { return 1; }


    /** Scalar norm square is sum( conjugate squares of all underlying scalars ), 
    where conjugate square of scalar s is conj(s)*s. **/ 
    ScalarNormSq scalarNormSqr() const { 
        ScalarNormSq sum(0);
        for(int i=0;i<M;++i) sum += CNT<E>::scalarNormSqr(d[i*STRIDE]);
        return sum;
    }

    /** Elementwise square root; that is, the return value has the same
    length as this Vec but with each element replaced by whatever it thinks
    its square root is. The element type may have changed and the stride
    of the return Vec is always 1. **/
    TSqrt sqrt() const {
        TSqrt vsqrt;
        for(int i=0;i<M;++i) vsqrt[i] = CNT<E>::sqrt(d[i*STRIDE]);
        return vsqrt;
    }

    /** Elementwise absolute value; that is, the return value has the same
    dimension as this Vec but with each element replaced by whatever it thinks
    its absolute value is. The element type may have changed and the stride
    of the return Vec is always 1. **/
    TAbs abs() const {
        TAbs vabs;
        for(int i=0;i<M;++i) vabs[i] = CNT<E>::abs(d[i*STRIDE]);
        return vabs;
    }

    /** Return a copy of this Vec but with the underlying scalar type
    converted (if necessary) to one of the C++ standard real or complex
    floating point types. This may require floating point negations to
    occur to get read of negator or conjugate types. **/
    TStandard standardize() const {
        TStandard vstd;
        for(int i=0;i<M;++i) vstd[i] = CNT<E>::standardize(d[i*STRIDE]);
        return vstd;
    }

    /** Sum just adds up all the elements into a single return element that
    is the same type as this Vec's elements except standardized to use one
    of the C++ built-in real or complex types as its underlying scalars. **/
    EStandard sum() const {
        E sum(0);
        for (int i=0;i<M;++i) sum += d[i*STRIDE];
        return CNT<E>::standardize(sum);
    }


    // This gives the resulting vector type when (v[i] op P) is applied to 
    // each element of v. It is a vector of length M, stride 1, and element 
    // types which are the regular composite result of E op P. Typically P is 
    // a scalar type but it doesn't have to be.
    template <class P> struct EltResult { 
        typedef Vec<M, typename CNT<E>::template Result<P>::Mul, 1> Mul;
        typedef Vec<M, typename CNT<E>::template Result<P>::Dvd, 1> Dvd;
        typedef Vec<M, typename CNT<E>::template Result<P>::Add, 1> Add;
        typedef Vec<M, typename CNT<E>::template Result<P>::Sub, 1> Sub;
    };

    // This is the composite result for v op P where P is some kind of 
    // appropriately shaped non-scalar type.
    template <class P> struct Result { 
        typedef MulCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOp;
        typedef typename MulOp::Type Mul;

        typedef MulCNTsNonConforming<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> MulOpNonConforming;
        typedef typename MulOpNonConforming::Type MulNon;

        typedef DvdCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> DvdOp;
        typedef typename DvdOp::Type Dvd;

        typedef AddCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> AddOp;
        typedef typename AddOp::Type Add;

        typedef SubCNTs<M,1,ArgDepth,Vec,ColSpacing,RowSpacing,
            CNT<P>::NRows, CNT<P>::NCols, CNT<P>::ArgDepth,
            P, CNT<P>::ColSpacing, CNT<P>::RowSpacing> SubOp;
        typedef typename SubOp::Type Sub;
    };

    /** Shape-preserving element substitution (always packed). That is, if
    T1==%Vec\<M,E,S> is the current %Vec type, then type
    T2==T1::%Substitute\<P>::%Type is the type of a shape-compatible,
    packed %Vec whose elements are of type P rather than E, that is,
    T2==%Vec\<M,P,1>. **/
    template <class P> struct Substitute {
        typedef Vec<M,P> Type;
    };

    /** Default construction initializes %Vec's elements to NaN when debugging 
    but leaves them uninitialized garbage otherwise, so declarations have zero
    cost in Release builds. **/
    Vec(){ 
    #ifndef NDEBUG
        setToNaN();
    #endif
    }

    // It's important not to use the default copy constructor or copy
    // assignment because the compiler doesn't understand that we may
    // have noncontiguous storage and will try to copy the whole array.

    /** Copy constructor copies the logically-included elements from the
    source %Vec; gaps due to stride are not accessed in either source or
    destination. **/
    Vec(const Vec& src) {
        Impl::copy(*this, src);
    }
    /** Copy assignment operator copies the logically-included elements from 
    the source %Vec; gaps due to stride are not accessed in either source or
    destination. OK if source and destination are the same vector; results
    are unpredictable if they otherwise overlap with elements in common. **/
    Vec& operator=(const Vec& src) {    
        Impl::copy(*this, src);
        return *this;
    }

    /** This is an implicit conversion from a %Vec of the same length
    and element type but with a different stride. **/
    template <int SS> Vec(const Vec<M,E,SS>& src) {
        Impl::copy(*this, src);
    }

    /** This is an implicit conversion from a %Vec of the same length
    and \e negated element type (possibly with a different stride). **/
    template <int SS> Vec(const Vec<M,ENeg,SS>& src) {
        Impl::copy(*this, src);
    }

    /** Construct a Vec from a Vec of the same length, with any stride. Works 
    as long as the element types are assignment compatible. **/
    template <class EE, int SS> explicit Vec(const Vec<M,EE,SS>& src) {
        Impl::copy(*this, src);
    }

    /** Construction from a single value of this %Vec's element type assigns
    that value to each element. **/
    explicit Vec(const E& e) {for (int i=0;i<M;++i) d[i*STRIDE]=e;}

    /** Construction from a single value of this %Vec's negated element type 
    assigns that value to each element, requiring floating point negation
    to be performed once to compute the type-E representation of the 
    type negator<E> value provided. **/
    explicit Vec(const ENeg& ne) {
        const E e = ne; // requires floating point negation
        for (int i=0;i<M;++i) d[i*STRIDE]=e;
    }

    /** Given an int value, turn it into a suitable floating point number,
    convert that to element type E and then feed that to the above 
    single-element constructor. 
    @see Vec::Vec(const E&). **/
    explicit Vec(int i) {new (this) Vec(E(Precision(i)));}

    // A bevy of constructors for Vecs up to length 9.

    /** Construct a Vec<2,E> from two elements of type E, etc. **/
    Vec(const E& e0,const E& e1)
      { assert(M==2);(*this)[0]=e0;(*this)[1]=e1; }
    Vec(const E& e0,const E& e1,const E& e2)
      { assert(M==3);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3)
      { assert(M==4);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;(*this)[3]=e3; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4)
      { assert(M==5);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5)
      { assert(M==6);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6)
      { assert(M==7);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6, const E& e7)
      { assert(M==8);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7; }
    Vec(const E& e0,const E& e1,const E& e2,const E& e3,const E& e4,const E& e5, const E& e6, const E& e7, const E& e8)
      { assert(M==9);(*this)[0]=e0;(*this)[1]=e1;(*this)[2]=e2;
        (*this)[3]=e3;(*this)[4]=e4;(*this)[5]=e5;(*this)[6]=e6;(*this)[7]=e7;(*this)[8]=e8; }

    /** Construction from a pointer to elements of any type EE assumes we're 
    pointing at a C++ array of EE's of the right length, and that EE is
    assignment compatible with this %Vec's element type E. The supplied
    pointer cannot be null. **/
    template <class EE> explicit Vec(const EE* p)
    {   assert(p); for(int i=0;i<M;++i) d[i*STRIDE]=p[i]; }

    /** Assignment to a pointer to elements of any type EE assumes we're 
    pointing at a C++ array of EE's of the right length, and that EE is
    assignment compatible with this %Vec's element type E. The supplied
    pointer cannot be null. **/
    template <class EE> Vec& operator=(const EE* p)
    {   assert(p); for(int i=0;i<M;++i) d[i*STRIDE]=p[i]; return *this; }

    /** Assignment to a conforming %Vec, of any element type and stride,
    provided that the element types are assignment-compatible. **/
    template <class EE, int SS> Vec& operator=(const Vec<M,EE,SS>& vv) 
    {   Impl::copy(*this, vv); return *this; }

    /** Add in a conforming %Vec, of any element type and stride,
    provided that the element types are addition-compatible. **/
    template <class EE, int SS> Vec& operator+=(const Vec<M,EE,SS>& r)
    {   for(int i=0;i<M;++i) d[i*STRIDE] += r[i]; return *this; }
    /** Add in a conforming %Vec, of any negated element type and stride,
    provided that the element types are addition-compatible. The negation
    is removed at zero cost by subtracting rather than adding. **/
    template <class EE, int SS> Vec& operator+=(const Vec<M,negator<EE>,SS>& r)
    {   for(int i=0;i<M;++i) d[i*STRIDE] -= -(r[i]); return *this; }

    /** Subtract off a conforming %Vec, of any element type and stride,
    provided that the element types are addition-compatible. **/
    template <class EE, int SS> Vec& operator-=(const Vec<M,EE,SS>& r)
    {   for(int i=0;i<M;++i) d[i*STRIDE] -= r[i]; return *this; }
    /** Subtract off  a conforming %Vec, of any negated element type and stride,
    provided that the element types are addition-compatible. The negation
    is removed at zero cost by adding rather than subtracting. **/
    template <class EE, int SS> Vec& operator-=(const Vec<M,negator<EE>,SS>& r)
    {   for(int i=0;i<M;++i) d[i*STRIDE] += -(r[i]); return *this; }

    // Conforming binary ops with 'this' on left, producing new packed result.
    // Cases: v=v+v, v=v-v, m=v*r

    /** Vector addition -- use operator+ instead. **/
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Add>
    conformingAdd(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Add> result;
        Impl::conformingAdd(*this, r, result);
        return result;
    }
    /** Vector subtraction -- use operator- instead. **/
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Sub>
    conformingSubtract(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Sub> result;
        Impl::conformingSubtract(*this, r, result);
        return result;
    }

    /** Same as outer product (m = col*row) -- use operator* or outer() 
    instead. **/
    template <class EE, int SS> Mat<M,M,typename CNT<E>::template Result<EE>::Mul>
    conformingMultiply(const Row<M,EE,SS>& r) const {
        Mat<M,M,typename CNT<E>::template Result<EE>::Mul> result;
        for (int j=0;j<M;++j) result(j) = scalarMultiply(r(j));
        return result;
    }

    /** Elementwise multiply (Matlab " .* " operator). **/
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Mul>
    elementwiseMultiply(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Mul> result;
        Impl::elementwiseMultiply(*this, r, result);
        return result;
    }
    /** Elementwise divide (Matlab " ./ " operator). **/
    template <class EE, int SS> Vec<M,typename CNT<E>::template Result<EE>::Dvd>
    elementwiseDivide(const Vec<M,EE,SS>& r) const {
        Vec<M,typename CNT<E>::template Result<EE>::Dvd> result;
        Impl::elementwiseDivide(*this, r, result);
        return result;
    }

    /** Select an element of this %Vec and return a const reference to it.
    This is range-checked in Debug builds but has zero overhead in Release
    builds. **/
    const E& operator[](int i) const 
    {   assert(0 <= i && i < M); return d[i*STRIDE]; }
    /** Same as const operator[] above. **/
    const E& operator()(int i) const {return (*this)[i];}

    /** Select an element of this %Vec and return a writable reference 
    to it. This is range-checked in Debug builds but has zero overhead in 
    Release builds. **/
    E& operator[](int i) {assert(0 <= i && i < M); return d[i*STRIDE];}
    /** Same as non-const operator[] above. **/
    E& operator()(int i) {return (*this)[i];}

    ScalarNormSq normSqr() const { return scalarNormSqr(); }
    typename CNT<ScalarNormSq>::TSqrt 
        norm() const { return CNT<ScalarNormSq>::sqrt(scalarNormSqr()); }

    /** If the elements of this Vec are scalars, the result is what you get by
    dividing each element by the norm() calculated above. If the elements are
    \e not scalars, then the elements are *separately* normalized. That means
    you will get a different answer from Vec<2,Vec3>::normalize() than you
    would from a Vec<6>::normalize() containing the same scalars.
    
    Normalize returns a vector of the same dimension but in new, packed storage
    and with a return type that does not include negator<> even if the original
    Vec<> does, because we can eliminate the negation here almost for free.
    But we can't standardize (change conjugate to complex) for free, so we'll retain
    conjugates if there are any. **/
    TNormalize normalize() const {
        if (CNT<E>::IsScalar) {
            return castAwayNegatorIfAny() / (SignInterpretation*norm());
        } else {
            TNormalize elementwiseNormalized;
            for (int i=0; i<M; ++i) 
                elementwiseNormalized[i] = CNT<E>::normalize((*this)[i]);
            return elementwiseNormalized;
        }
    }

    /** This method is not supported for %Vec objects. **/
    TInvert invert() const {assert(false); return TInvert();} // TODO default inversion

    /** Unary plus does nothing. **/
    const Vec&   operator+() const { return *this; }
    /** Unary minus recasts this %Vec to a type that has the opposite 
    interpretation of the sign but is otherwise identical, so no computation
    or copying is performed here. **/
    const TNeg&  operator-() const { return negate(); }
    /** Recast to negated type and return a writable reference; writing to
    this will cause the negated result to be placed in the original %Vec. **/
    TNeg&        operator-()       { return updNegate(); }
    /** The Hermitian transpose operator recasts this %Vec to a type that
    specifies the opposite storage order (row vs.\ column) then returns a 
    reference, so no computation or copying is performed here. **/
    const THerm& operator~() const { return transpose(); }
    /** Recast to Hermitian transposed type and return a writable reference;
    the effect is that writing to elements of the result affects the transposed
    element of the original %Vec. **/
    THerm&       operator~()       { return updTranspose(); }

    /** Non-operator version of unary negation; just a recast. **/
    const TNeg&  negate() const { return *reinterpret_cast<const TNeg*>(this); }
    /** Non-operator version of unary negation; recasts and returns a 
    writable reference. **/
    TNeg&        updNegate()    { return *reinterpret_cast<      TNeg*>(this); }

    /** Non-operator version of Hermitian transpose; just a recast. **/
    const THerm& transpose()    const { return *reinterpret_cast<const THerm*>(this); }
    /** Non-operator version of Hermitian transpose; recasts and returns a 
    writable reference. **/
    THerm&       updTranspose()       { return *reinterpret_cast<      THerm*>(this); }

    /** Positional transpose turns this %Vec into a Row but does not transpose
    the individual elements. That is, a Vec<2,Vec3> becomes a Row<2,Vec3>,
    rather than a Row<2,Row3> as would happen with ordinary transpose(). This
    is just a recast; no copying or computation is performed here. **/
    const TPosTrans& positionalTranspose() const
        { return *reinterpret_cast<const TPosTrans*>(this); }
    /** Positional transpose returning a writable reference. **/
    TPosTrans&       updPositionalTranspose()
        { return *reinterpret_cast<TPosTrans*>(this); }

    /** Return a reference to the real portion of this %Vec if it has complex
    elements; otherwise the type doesn't change. This is just a recast; no
    copying or computation is done here. The result may have a different 
    stride than the original since the imaginary parts must be skipped. **/
    const TReal& real() const { return *reinterpret_cast<const TReal*>(this); }
    /** Recast to show only the real portion of this %Vec and return a writable
    reference. **/
    TReal&       real()       { return *reinterpret_cast<      TReal*>(this); }

    // Had to contort these next two routines to get them through VC++ 7.net

    /** Return a reference to the imaginary portion of this %Vec if it has 
    complex elements; otherwise the type doesn't change. This is just a recast; 
    no copying or computation is done here. The result may have a different 
    stride than the original since the real parts must be skipped. **/
    const TImag& imag()    const { 
        const int offs = ImagOffset;
        const EImag* p = reinterpret_cast<const EImag*>(this);
        return *reinterpret_cast<const TImag*>(p+offs);
    }
    /** Recast to show only the imaginary portion of this %Vec and return a 
    writable reference. **/
    TImag& imag() { 
        const int offs = ImagOffset;
        EImag* p = reinterpret_cast<EImag*>(this);
        return *reinterpret_cast<TImag*>(p+offs);
    }

    /** Recast to remove negators from this %Vec's type if present; this is
    handy for simplifying operations where we know the sign can be ignored
    such as squaring. **/
    const TWithoutNegator& castAwayNegatorIfAny() const 
    {   return *reinterpret_cast<const TWithoutNegator*>(this); }
    /** Recast to remove negators from this %Vec's type if present and return
    a writable reference. **/
    TWithoutNegator&       updCastAwayNegatorIfAny()    
    {   return *reinterpret_cast<TWithoutNegator*>(this); }

    // These are elementwise binary operators, (this op ee) by default but 
    // (ee op this) if 'FromLeft' appears in the name. The result is a packed 
    // Vec<M> but the element type may change. These are mostly used to 
    // implement global operators. We call these "scalar" operators but 
    // actually the "scalar" can be a composite type.

    //TODO: consider converting 'e' to Standard Numbers as precalculation and 
    // changing return type appropriately.
    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Mul>
    scalarMultiply(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Mul> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] * e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Mul>
    scalarMultiplyFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Mul> result;
        for (int i=0; i<M; ++i) result[i] = e * (*this)[i];
        return result;
    }

    // TODO: should precalculate and store 1/e, while converting to Standard 
    // Numbers. Note that return type should change appropriately.
    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Dvd>
    scalarDivide(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Dvd> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] / e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Dvd>
    scalarDivideFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Dvd> result;
        for (int i=0; i<M; ++i) result[i] = e / (*this)[i];
        return result;
    }

    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Add>
    scalarAdd(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Add> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] + e;
        return result;
    }
    // Add is commutative, so no 'FromLeft'.

    template <class EE> Vec<M, typename CNT<E>::template Result<EE>::Sub>
    scalarSubtract(const EE& e) const {
        Vec<M, typename CNT<E>::template Result<EE>::Sub> result;
        for (int i=0; i<M; ++i) result[i] = (*this)[i] - e;
        return result;
    }
    template <class EE> Vec<M, typename CNT<EE>::template Result<E>::Sub>
    scalarSubtractFromLeft(const EE& e) const {
        Vec<M, typename CNT<EE>::template Result<E>::Sub> result;
        for (int i=0; i<M; ++i) result[i] = e - (*this)[i];
        return result;
    }

    // Generic assignments for any element type not listed explicitly, including scalars.
    // These are done repeatedly for each element and only work if the operation can
    // be performed leaving the original element type.
    template <class EE> Vec& operator =(const EE& e) {return scalarEq(e);}
    template <class EE> Vec& operator+=(const EE& e) {return scalarPlusEq(e);}
    template <class EE> Vec& operator-=(const EE& e) {return scalarMinusEq(e);}
    template <class EE> Vec& operator*=(const EE& e) {return scalarTimesEq(e);}
    template <class EE> Vec& operator/=(const EE& e) {return scalarDivideEq(e);}

    // Generalized element assignment & computed assignment methods. These will work
    // for any assignment-compatible element, not just scalars.
    template <class EE> Vec& scalarEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee; return *this; }
    template <class EE> Vec& scalarPlusEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] += ee; return *this; }
    template <class EE> Vec& scalarMinusEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] -= ee; return *this; }
    template <class EE> Vec& scalarMinusEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee - d[i*STRIDE]; return *this; }
    template <class EE> Vec& scalarTimesEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] *= ee; return *this; }
    template <class EE> Vec& scalarTimesEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee * d[i*STRIDE]; return *this; }
    template <class EE> Vec& scalarDivideEq(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] /= ee; return *this; }
    template <class EE> Vec& scalarDivideEqFromLeft(const EE& ee)
      { for(int i=0;i<M;++i) d[i*STRIDE] = ee / d[i*STRIDE]; return *this; }

    // Specialize for int to avoid warnings and ambiguities.
    Vec& scalarEq(int ee)       {return scalarEq(Precision(ee));}
    Vec& scalarPlusEq(int ee)   {return scalarPlusEq(Precision(ee));}
    Vec& scalarMinusEq(int ee)  {return scalarMinusEq(Precision(ee));}
    Vec& scalarTimesEq(int ee)  {return scalarTimesEq(Precision(ee));}
    Vec& scalarDivideEq(int ee) {return scalarDivideEq(Precision(ee));}
    Vec& scalarMinusEqFromLeft(int ee)  {return scalarMinusEqFromLeft(Precision(ee));}
    Vec& scalarTimesEqFromLeft(int ee)  {return scalarTimesEqFromLeft(Precision(ee));}
    Vec& scalarDivideEqFromLeft(int ee) {return scalarDivideEqFromLeft(Precision(ee));}

    /** Set every scalar in this %Vec to NaN; this is the default initial
    value in Debug builds, but not in Release. **/
    void setToNaN() {
        (*this) = CNT<ELT>::getNaN();
    }

    /** Set every scalar in this %Vec to zero. **/
    void setToZero() {
        (*this) = ELT(0);
    }

    /** Extract a const reference to a sub-Vec with size known at compile time. 
    This must be called with an explicit template argument for the size, for
    example, getSubVec<3>(i). This is only a recast; no copying or computation
    is performed. The size and index are range checked in Debug builds but
    not in Release builds. **/
    template <int MM>
    const Vec<MM,ELT,STRIDE>& getSubVec(int i) const {
        assert(0 <= i && i + MM <= M);
        return Vec<MM,ELT,STRIDE>::getAs(&(*this)[i]);
    }
    /** Extract a writable reference to a sub-Vec with size known at compile time. 
    This must be called with an explicit template argument for the size, for
    example, updSubVec<3>(i). This is only a recast; no copying or computation
    is performed. The size and index are range checked in Debug builds but
    not in Release builds. **/
    template <int MM>
    Vec<MM,ELT,STRIDE>& updSubVec(int i) {
        assert(0 <= i && i + MM <= M);
        return Vec<MM,ELT,STRIDE>::updAs(&(*this)[i]);
    }


    /** Extract a subvector of type %Vec from a longer one that has the same
    element type and stride, and return a const reference to the selected 
    subsequence. **/
    template <int MM>
    static const Vec& getSubVec(const Vec<MM,ELT,STRIDE>& v, int i) {
        assert(0 <= i && i + M <= MM);
        return getAs(&v[i]);
    }
    /** Extract a subvector of type %Vec from a longer one that has the same
    element type and stride, and return a writable reference to the selected 
    subsequence. **/
    template <int MM>
    static Vec& updSubVec(Vec<MM,ELT,STRIDE>& v, int i) {
        assert(0 <= i && i + M <= MM);
        return updAs(&v[i]);
    }

    /** Return a vector one smaller than this one by dropping the element
    at the indicated position p. The result is a packed copy with the same
    element type as this one. **/
    Vec<M-1,ELT,1> drop1(int p) const {
        assert(0 <= p && p < M);
        Vec<M-1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<M-1; ++i, ++nxt) {
            if (nxt==p) ++nxt;  // skip the loser
            out[i] = (*this)[nxt];
        }
        return out;
    }

    /** Return a vector one larger than this one by adding an element
    to the end. The result is a packed copy with the same element type as
    this one. Works for any assignment compatible element. **/
    template <class EE> Vec<M+1,ELT,1> append1(const EE& v) const {
        Vec<M+1,ELT,1> out;
        Vec<M,ELT,1>::updAs(&out[0]) = (*this);
        out[M] = v;
        return out;
    }


    /** Return a vector one larger than this one by inserting an element
    \e before the indicated one. The result is a packed copy with the same 
    element type as this one. Works for any assignment compatible element. The 
    index can be one greater than normally allowed in which case the element
    is appended (but use append1() if you know you're appending). **/
    template <class EE> Vec<M+1,ELT,1> insert1(int p, const EE& v) const {
        assert(0 <= p && p <= M);
        if (p==M) return append1(v);
        Vec<M+1,ELT,1> out;
        int nxt=0;
        for (int i=0; i<M; ++i, ++nxt) {
            if (i==p) out[nxt++] = v;
            out[nxt] = (*this)[i];
        }
        return out;
    }
            
    /** Recast an ordinary C++ array E[] to a const %Vec<M,E,S>; assumes 
    compatible length, stride, and packing. **/
    static const Vec& getAs(const ELT* p)  
    {   return *reinterpret_cast<const Vec*>(p); }
    /** Recast a writable ordinary C++ array E[] to a writable %Vec<M,E,S>; 
    assumes compatible length, stride, and packing. **/
    static Vec&       updAs(ELT* p)
    {   return *reinterpret_cast<Vec*>(p); }


    /** Return a %Vec of the same length and element type as this one but
    with all elements set to NaN. The result is packed (stride==1) regardless
    of the stride of this %Vec. **/
    static Vec<M,ELT,1> getNaN() { return Vec<M,ELT,1>(CNT<ELT>::getNaN()); }

    /** Return true if any element of this Vec contains a NaN anywhere. **/
    bool isNaN() const {
        for (int i=0; i<M; ++i)
            if (CNT<ELT>::isNaN((*this)[i]))
                return true;
        return false;
    }

    /** Return true if any element of this Vec contains a +Infinity
    or -Infinity somewhere but no element contains a NaN anywhere. **/
    bool isInf() const {
        bool seenInf = false;
        for (int i=0; i<M; ++i) {
            const ELT& e = (*this)[i];
            if (!CNT<ELT>::isFinite(e)) {
                if (!CNT<ELT>::isInf(e)) 
                    return false; // something bad was found
                seenInf = true; 
            }
        }
        return seenInf;
    }

    /** Return true if no element of this %Vec contains an Infinity or a NaN 
    anywhere. **/
    bool isFinite() const {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isFinite((*this)[i]))
                return false;
        return true;
    }

    /** For approximate comparisons, the default tolerance to use for a vector is
    the same as its elements' default tolerance. **/
    static double getDefaultTolerance() {return CNT<ELT>::getDefaultTolerance();}

    /** %Test whether this vector is numerically equal to some other vector with
    the same shape, using a specified tolerance. **/
    template <class E2, int RS2>
    bool isNumericallyEqual(const Vec<M,E2,RS2>& v, double tol) const {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isNumericallyEqual((*this)[i], v[i], tol))
                return false;
        return true;
    }

    /** %Test whether this vector is numerically equal to some other vector with
    the same shape, using a default tolerance which is the looser of the
    default tolerances of the two objects being compared. **/
    template <class E2, int RS2>
    bool isNumericallyEqual(const Vec<M,E2,RS2>& v) const {
        const double tol = std::max(getDefaultTolerance(),v.getDefaultTolerance());
        return isNumericallyEqual(v, tol);
    }

    /** %Test whether every element of this vector is numerically equal to the given
    element, using either a specified tolerance or the vector's 
    default tolerance (which is always the same or looser than the default
    tolerance for one of its elements). **/
    bool isNumericallyEqual
       (const ELT& e,
        double     tol = getDefaultTolerance()) const 
    {
        for (int i=0; i<M; ++i)
            if (!CNT<ELT>::isNumericallyEqual((*this)[i], e, tol))
                return false;
        return true;
    }

    // Functions to be used for Scripting in MATLAB and languages that do not support operator overloading
    /** Print Vec into a string and return it.  Please refer to operator<< for details. **/
    std::string toString() const {
        std::stringstream stream;
        stream <<  (*this);
        return stream.str(); 
    }

    /** Variant of operator[] that's scripting friendly to set ith entry **/
    void set(int i, const E& value)  
    {   (*this)[i] = value; }

    /** Variant of operator[] that's scripting friendly to get const reference to ith entry **/
    const E& get(int i) const 
    {   return operator[](i); }

private:
    // TODO: should be an array of scalars rather than elements to control
    // packing more carefully.
    ELT d[NActualElements];    // data
};

/////////////////////////////////////////////
// Global operators involving two vectors. //
//   v+v, v-v, v==v, v!=v                  //
/////////////////////////////////////////////

// v3 = v1 + v2 where all v's have the same length M. 
template <int M, class E1, int S1, class E2, int S2> inline
typename Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >::Add
operator+(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) {
    return Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >
        ::AddOp::perform(l,r);
}

// v3 = v1 - v2, similar to +
template <int M, class E1, int S1, class E2, int S2> inline
typename Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >::Sub
operator-(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) { 
    return Vec<M,E1,S1>::template Result< Vec<M,E2,S2> >
        ::SubOp::perform(l,r);
}

/// bool = v1[i] == v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator==(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] != r[i]) return false;
    return true; }
/// bool = v1[i] != v2[i], for any element i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator!=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) {return !(l==r);} 

/// bool = v[i] == e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator==(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] != e) return false;
    return true; }
/// bool = v[i] != e, for any element v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator!=(const Vec<M,E1,S1>& v, const E2& e) {return !(v==e);} 

/// bool = v1[i] < v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator<(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] >= r[i]) return false;
    return true; }
/// bool = v[i] < e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator<(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] >= e) return false;
    return true; }

/// bool = v1[i] > v2[i], for all elements i
template <int M, class E1, int S1, class E2, int S2> inline bool
operator>(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] <= r[i]) return false;
    return true; }
/// bool = v[i] > e, for all elements v[i] and element e
template <int M, class E1, int S1, class E2> inline bool
operator>(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] <= e) return false;
    return true; }

/// bool = v1[i] <= v2[i], for all elements i.
/// This is not the same as !(v1>v2).
template <int M, class E1, int S1, class E2, int S2> inline bool
operator<=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] > r[i]) return false;
    return true; }
/// bool = v[i] <= e, for all elements v[i] and element e.
/// This is not the same as !(v1>e).
template <int M, class E1, int S1, class E2> inline bool
operator<=(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] > e) return false;
    return true; }

/// bool = v1[i] >= v2[i], for all elements i
/// This is not the same as !(v1<v2).
template <int M, class E1, int S1, class E2, int S2> inline bool
operator>=(const Vec<M,E1,S1>& l, const Vec<M,E2,S2>& r) 
{   for (int i=0; i < M; ++i) if (l[i] < r[i]) return false;
    return true; }
/// bool = v[i] >= e, for all elements v[i] and element e.
/// This is not the same as !(v1<e).
template <int M, class E1, int S1, class E2> inline bool
operator>=(const Vec<M,E1,S1>& v, const E2& e) 
{   for (int i=0; i < M; ++i) if (v[i] < e) return false;
    return true; }

///////////////////////////////////////////////////////
// Global operators involving a vector and a scalar. //
///////////////////////////////////////////////////////

// I haven't been able to figure out a nice way to templatize for the
// built-in reals without introducing a lot of unwanted type matches
// as well. So we'll just grind them out explicitly here.

// SCALAR MULTIPLY

// v = v*real, real*v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Mul
operator*(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Mul
operator*(const float& l, const Vec<M,E,S>& r) {return r*l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Mul
operator*(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::MulOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Mul
operator*(const double& l, const Vec<M,E,S>& r) {return r*l;}

// v = v*int, int*v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(const Vec<M,E,S>& l, int r) {return l * (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Mul
operator*(int l, const Vec<M,E,S>& r) {return r * (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v*complex, complex*v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::MulOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const std::complex<R>& l, const Vec<M,E,S>& r) {return r*l;}

// v = v*conjugate, conjugate*v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const Vec<M,E,S>& l, const conjugate<R>& r) {return l*(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Mul
operator*(const conjugate<R>& l, const Vec<M,E,S>& r) {return r*(std::complex<R>)l;}

// v = v*negator, negator*v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const Vec<M,E,S>& l, const negator<R>& r) {return l * (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Mul
operator*(const negator<R>& l, const Vec<M,E,S>& r) {return r * (typename negator<R>::StdNumber)(R)l;}


// SCALAR DIVIDE. This is a scalar operation when the scalar is on the right,
// but when it is on the left it means scalar * pseudoInverse(vec), which is 
// a row.

// v = v/real, real/v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Dvd
operator/(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<Vec<M,E,S> >::Dvd
operator/(const float& l, const Vec<M,E,S>& r)
  { return CNT<float>::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Dvd
operator/(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::DvdOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<Vec<M,E,S> >::Dvd
operator/(const double& l, const Vec<M,E,S>& r)
  { return CNT<double>::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

// v = v/int, int/v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Dvd
operator/(const Vec<M,E,S>& l, int r) {return l / (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Vec<M,E,S> >::Dvd
operator/(int l, const Vec<M,E,S>& r) {return (typename CNT<E>::Precision)l / r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v/complex, complex/v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::DvdOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Dvd
operator/(const std::complex<R>& l, const Vec<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Vec<M,E,S> >::DvdOp::perform(l,r); }

// v = v/conjugate, conjugate/v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Dvd
operator/(const Vec<M,E,S>& l, const conjugate<R>& r) {return l/(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Dvd
operator/(const conjugate<R>& l, const Vec<M,E,S>& r) {return (std::complex<R>)l/r;}

// v = v/negator, negator/v: convert negator to number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Dvd
operator/(const Vec<M,E,S>& l, const negator<R>& r) {return l/(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<Vec<M,E,S> >::Dvd
operator/(const negator<R>& l, const Vec<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l/r;}


// Add and subtract are odd as scalar ops. They behave as though the
// scalar stands for a vector each of whose elements is that scalar,
// and then a normal vector add or subtract is done.

// SCALAR ADD

// v = v+real, real+v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Add
operator+(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Add
operator+(const float& l, const Vec<M,E,S>& r) {return r+l;}

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Add
operator+(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::AddOp::perform(l,r); }
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Add
operator+(const double& l, const Vec<M,E,S>& r) {return r+l;}

// v = v+int, int+v -- just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(const Vec<M,E,S>& l, int r) {return l + (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Add
operator+(int l, const Vec<M,E,S>& r) {return r + (typename CNT<E>::Precision)l;}

// Complex, conjugate, and negator are all easy to templatize.

// v = v+complex, complex+v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::AddOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const std::complex<R>& l, const Vec<M,E,S>& r) {return r+l;}

// v = v+conjugate, conjugate+v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const Vec<M,E,S>& l, const conjugate<R>& r) {return l+(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Add
operator+(const conjugate<R>& l, const Vec<M,E,S>& r) {return r+(std::complex<R>)l;}

// v = v+negator, negator+v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const Vec<M,E,S>& l, const negator<R>& r) {return l + (typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Add
operator+(const negator<R>& l, const Vec<M,E,S>& r) {return r + (typename negator<R>::StdNumber)(R)l;}

// SCALAR SUBTRACT -- careful, not commutative.

// v = v-real, real-v 
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<float>::Sub
operator-(const Vec<M,E,S>& l, const float& r)
  { return Vec<M,E,S>::template Result<float>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<float>::template Result<Vec<M,E,S> >::Sub
operator-(const float& l, const Vec<M,E,S>& r)
  { return CNT<float>::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<double>::Sub
operator-(const Vec<M,E,S>& l, const double& r)
  { return Vec<M,E,S>::template Result<double>::SubOp::perform(l,r); }
template <int M, class E, int S> inline
typename CNT<double>::template Result<Vec<M,E,S> >::Sub
operator-(const double& l, const Vec<M,E,S>& r)
  { return CNT<double>::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

// v = v-int, int-v // just convert int to v's precision float
template <int M, class E, int S> inline
typename Vec<M,E,S>::template Result<typename CNT<E>::Precision>::Sub
operator-(const Vec<M,E,S>& l, int r) {return l - (typename CNT<E>::Precision)r;}
template <int M, class E, int S> inline
typename CNT<typename CNT<E>::Precision>::template Result<Vec<M,E,S> >::Sub
operator-(int l, const Vec<M,E,S>& r) {return (typename CNT<E>::Precision)l - r;}


// Complex, conjugate, and negator are all easy to templatize.

// v = v-complex, complex-v
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const Vec<M,E,S>& l, const std::complex<R>& r)
  { return Vec<M,E,S>::template Result<std::complex<R> >::SubOp::perform(l,r); }
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Sub
operator-(const std::complex<R>& l, const Vec<M,E,S>& r)
  { return CNT<std::complex<R> >::template Result<Vec<M,E,S> >::SubOp::perform(l,r); }

// v = v-conjugate, conjugate-v (convert conjugate->complex)
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<std::complex<R> >::Sub
operator-(const Vec<M,E,S>& l, const conjugate<R>& r) {return l-(std::complex<R>)r;}
template <int M, class E, int S, class R> inline
typename CNT<std::complex<R> >::template Result<Vec<M,E,S> >::Sub
operator-(const conjugate<R>& l, const Vec<M,E,S>& r) {return (std::complex<R>)l-r;}

// v = v-negator, negator-v: convert negator to standard number
template <int M, class E, int S, class R> inline
typename Vec<M,E,S>::template Result<typename negator<R>::StdNumber>::Sub
operator-(const Vec<M,E,S>& l, const negator<R>& r) {return l-(typename negator<R>::StdNumber)(R)r;}
template <int M, class E, int S, class R> inline
typename CNT<R>::template Result<Vec<M,E,S> >::Sub
operator-(const negator<R>& l, const Vec<M,E,S>& r) {return (typename negator<R>::StdNumber)(R)l-r;}

// Vec I/O
template <int M, class E, int S, class CHAR, class TRAITS> inline
std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& o, const Vec<M,E,S>& v) {
    o << "~[" << v[0]; for(int i=1;i<M;++i) o<<','<<v[i]; o<<']'; return o;
}

/** Read a Vec from a stream as M elements separated by white space or
by commas, optionally enclosed in () [] ~() or ~[]. **/
template <int M, class E, int S, class CHAR, class TRAITS> inline
std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, Vec<M,E,S>& v) {
    CHAR tilde;
    is >> tilde; if (is.fail()) return is;
    if (tilde != CHAR('~')) {
        tilde = CHAR(0);
        is.unget(); if (is.fail()) return is;
    }

    CHAR openBracket, closeBracket;
    is >> openBracket; if (is.fail()) return is;
    if (openBracket==CHAR('('))
        closeBracket = CHAR(')');
    else if (openBracket==CHAR('['))
        closeBracket = CHAR(']');
    else {
        closeBracket = CHAR(0);
        is.unget(); if (is.fail()) return is;
    }

    // If we saw a "~" but then we didn't see any brackets, that's an
    // error. Set the fail bit and return.
    if (tilde != CHAR(0) && closeBracket == CHAR(0)) {
        is.setstate( std::ios::failbit );
        return is;
    }

    for (int i=0; i < M; ++i) {
        is >> v[i];
        if (is.fail()) return is;
        if (i != M-1) {
            CHAR c; is >> c; if (is.fail()) return is;
            if (c != ',') is.unget();
            if (is.fail()) return is;
        }
    }

    // Get the closing bracket if there was an opening one. If we don't
    // see the expected character we'll set the fail bit in the istream.
    if (closeBracket != CHAR(0)) {
        CHAR closer; is >> closer; if (is.fail()) return is;
        if (closer != closeBracket) {
            is.unget(); if (is.fail()) return is;
            is.setstate( std::ios::failbit );
        }
    }

    return is;
}

} //namespace SimTK


#endif //SimTK_SIMMATRIX_SMALLMATRIX_VEC_H_
