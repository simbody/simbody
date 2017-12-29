#ifndef SimTK_SIMMATRIX_CONJUGATE_H_
#define SimTK_SIMMATRIX_CONJUGATE_H_

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
 * This file defines the conjugate<R> template class, where R is one of
 * the three built-in real types. It is exactly like the C++ complex<R>
 * template except that the represented value is the conjugate of the value
 * represented by a complex number containing the same bit pattern. That is,
 * complex and conjugate both contain two real numbers, re and im, with
 * complex(re,im) meaning (re + im*i) while conjugate(re,im) means (re - im*i).
 * It is guaranteed that our conjugate type has the identical size and
 * representation as complex. Together, these definitions and
 * guarantee permit conjugation to be done by reinterpretation (i.e. type
 * casting) rather than by computation.  
 *
 * A caution on implementation of complex numbers: it's not as simple
 * as you might think. A mildly nasty issue is handling mixed precision
 * arguments. More substantial is getting the right complex answer using
 * real arithmetic. For example, if c=a+bi then even a seemingly simple
 * calculation like |c|=sqrt(a*a+b*b) may fail if calculated that way,
 * due to unnecessary overflow. Take a look at a good complex<>
 * implementation or see Press, et al., Numerical Recipes in C++, 2nd ed.,
 * 2003, section 5.4, pp 183-4 for a discussion. Consequently, wherever
 * it doesn't matter for performance, conjugate<R> objects are simply
 * converted to complex<R> (which costs one negation) and then the
 * system class is used to do the hard stuff.
 *
 * References for the complex<> class, which I used as a model
 * for SimTK::conjugate<>:
 *   (1) Stroustrup, B., The C++ Programming Language, 3rd ed., 1997, section
 *       22.5, pg 679ff.
 *   (2) The C++ Standard Incorporating Technical Corrigendum 1, British
 *       Standards Institute ISO/IEC 14882:2003 2nd ed., 2003, section 26.2,
 *       pg 592ff.
 */

#include <complex>
#include <iostream>
#include <limits>

using std::complex;


// These mixed-precision real/complex operators should not be needed
// and may have to be removed on some systems. However neither
// gcc 3.4.4 nor VC++ 7 provide them.
// Define the symbol below if these conflict with standard library operators.

#ifndef SimTK_MIXED_PRECISION_REAL_COMPLEX_ALREADY_DEFINED
namespace SimTK {
// complex<float> with int, double
inline complex<float> operator*(const complex<float>& c,int r) {return c*(float)r;}
inline complex<float> operator*(int r,const complex<float>& c) {return (float)r*c;}
inline complex<double> operator*(const complex<float>& c,const double& r)           {return complex<double>(c)*r;}
inline complex<double> operator*(const double& r,const complex<float>& c)           {return r*complex<double>(c);}

inline complex<float> operator/(const complex<float>& c,int r) {return c/(float)r;}
inline complex<float> operator/(int r,const complex<float>& c) {return (float)r/c;}
inline complex<double> operator/(const complex<float>& c,const double& r)           {return complex<double>(c)/r;}
inline complex<double> operator/(const double& r,const complex<float>& c)           {return r/complex<double>(c);}

inline complex<float> operator+(const complex<float>& c,int r) {return c+(float)r;}
inline complex<float> operator+(int r,const complex<float>& c) {return (float)r+c;}
inline complex<double> operator+(const complex<float>& c,const double& r)           {return complex<double>(c)+r;}
inline complex<double> operator+(const double& r,const complex<float>& c)           {return r+complex<double>(c);}

inline complex<float> operator-(const complex<float>& c,int r) {return c-(float)r;}
inline complex<float> operator-(int r,const complex<float>& c) {return (float)r-c;}
inline complex<double> operator-(const complex<float>& c,const double& r)           {return complex<double>(c)-r;}
inline complex<double> operator-(const double& r,const complex<float>& c)           {return r-complex<double>(c);}

// complex<double> with int, float
inline complex<double> operator*(const complex<double>& c,int r) {return c*(double)r;}
inline complex<double> operator*(int r,const complex<double>& c) {return (double)r*c;}
inline complex<double> operator*(const complex<double>& c,const float& r)           {return c*(double)r;}
inline complex<double> operator*(const float& r,const complex<double>& c)           {return (double)r*c;}

inline complex<double> operator/(const complex<double>& c,int r) {return c/(double)r;}
inline complex<double> operator/(int r,const complex<double>& c) {return (double)r/c;}
inline complex<double> operator/(const complex<double>& c,const float& r)           {return c/(double)r;}
inline complex<double> operator/(const float& r,const complex<double>& c)           {return (double)r/c;}

inline complex<double> operator+(const complex<double>& c,int r) {return c+(double)r;}
inline complex<double> operator+(int r,const complex<double>& c) {return (double)r+c;}
inline complex<double> operator+(const complex<double>& c,const float& r)           {return c+(double)r;}
inline complex<double> operator+(const float& r,const complex<double>& c)           {return (double)r+c;}

inline complex<double> operator-(const complex<double>& c,int r) {return c-(double)r;}
inline complex<double> operator-(int r,const complex<double>& c) {return (double)r-c;}
inline complex<double> operator-(const complex<double>& c,const float& r)           {return c-(double)r;}
inline complex<double> operator-(const float& r,const complex<double>& c)           {return (double)r-c;}
} // namespace SimTK
#endif
    
namespace SimTK {

template <class R> class conjugate;    // Only defined for float, double
template <> class conjugate<float>;
template <> class conjugate<double>;

// This is an adaptor for number types which negates the apparent values. A
// negator<N> has exactly the same internal representation as a number
// type N, but it is to be interpreted has having the negative of the value
// it would have if interpreted as an N. This permits negation to be done
// by reinterpretation rather than computation. A full set of arithmetic operators
// are provided involving negator<N>'s and N's. Sometimes we can save an op or
// two this way. For example negator<N>*negator<N> can be performed as an N*N
// since the negations cancel, and we saved two floating point negations.
template <class N> class negator;      // Only defined for numbers

/*
 * This class is specialized for all 9 combinations of built-in real types
 * and contains a typedef for the appropriate "widened" real, complex, or
 * conjugate type for use when R1 & R2 appear in an operation together.
 */
template <class R1, class R2> struct Wider {/* Only defined for built-ins. */};
template <> struct Wider<float,float> {
    typedef float               WReal;
    typedef complex<float>      WCplx;
    typedef conjugate<float>    WConj;
};
template <> struct Wider<float,double> {
    typedef double              WReal;
    typedef complex<double>     WCplx;
    typedef conjugate<double>   WConj;
};
template <> struct Wider<double,float> {
    typedef double              WReal;
    typedef complex<double>     WCplx;
    typedef conjugate<double>   WConj;
};
template <> struct Wider<double,double> {
    typedef double              WReal;
    typedef complex<double>     WCplx;
    typedef conjugate<double>   WConj;
};


/**
 * SimTK::conjugate<R> should be instantiated only for float, double.
 * This should behave just like std::complex<R> and in most cases we'll just
 * convert and punt to the std class.
 *
 * The three specializations are almost identical, but differ in the rules
 * for implicit conversion. Implicit conversions are allowed from a narrow
 * type to a wider one, but the other direction must be explicit.
 *
 * Some of this could be done more compactly with member templates, helper
 * classes and so on, but this can make debugging of client programs very
 * awkward, both by producing indecipherable error messages when compiling
 * and by referencing unexpanded template code when debugging. Since this
 * class need only be done once, and has a grand total of three specializations,
 * I felt it made more sense to do each of them explicitly here (sherm 051006).
 */
template <class R> class conjugate {/*Only defined for float, double*/};

/////////////////////////////////////////
// Specialization for conjugate<float> //
/////////////////////////////////////////

template <>  class conjugate<float> {
public:
    conjugate() {
    #ifndef NDEBUG
        re = negIm = std::numeric_limits<float>::quiet_NaN();
    #endif  
    }
    // default copy constructor, copy assignment, destructor

    /// Construction from reals. Note that the numeric result is (real-imag*i).
    conjugate(const float& real, const float& imag) { re = real; negIm = imag; }
    conjugate(const float& real, int i) { re = real; negIm = float(i); }
    conjugate(int r, const float& imag) { re = float(r); negIm = imag; }
    conjugate(int r, int i) { re = float(r); negIm = float(i); }

    /// Implicit conversion from float to conjugate<float>.
    conjugate(const float& real) { re = real; negIm = 0.f; }
    conjugate(int r) { re = float(r); negIm = 0.f; }

    // No implicit conversions from double because precision will be lost. Some
    // definitions must be deferred until conjugate<double> is defined below.
    inline explicit conjugate(const conjugate<double>& cd);

    explicit conjugate(const double& rd)
      { re = float(rd); negIm = 0.f; }

    // Conversions from complex are always explicit. Note that the value
    // represented by the conjugate must be identical to that represented by
    // the complex, which means we must negate the imaginary part.
    explicit conjugate(const complex<float>& x)
      { re = x.real(); negIm = -x.imag(); }
    explicit conjugate(const complex<double>& x)
      { re = float(x.real()); negIm = float(-x.imag()); }

    /// Implicit conversion to complex<float> when necessary
    /// (costs an actual negation -- yuck!).
    operator complex<float>() const
      { return complex<float>(re,-negIm); } 
    
    // Can't defer here by casting to negator<conjugate> -- this must act
    // like a built-in. But ... we can use this as a chance to convert
    // to complex and save one negation.
    complex<float> operator-() const { return complex<float>(-re,negIm); }

    // Useless.
    const conjugate& operator+() const { return *this; }

    // Computed assignment operators. We don't depend on implicit conversions
    // from reals to conjugates here because we can save a few flops by handling
    // the reals explicitly. Note that we only provide operators for implicitly
    // convertible precisions, though, which in this case means only floats.
    conjugate& operator=(const float& r)
      { re = r; negIm = 0.f; return *this; }
    conjugate& operator+=(const float& r)
      { re += r; return *this; }
    conjugate& operator-=(const float& r)
      { re -= r; return *this; }
    conjugate& operator*=(const float& r)
      { re *= r; negIm *= r; return *this; }
    conjugate& operator/=(const float& r)
      { re /= r; negIm /= r; return *this; }

    conjugate& operator+=(const conjugate<float>& c)
      { re += c.re; negIm += c.negIm; return *this; }
    conjugate& operator-=(const conjugate<float>& c)
      { re -= c.re; negIm -= c.negIm; return *this; }

    conjugate& operator=(const complex<float>& c)
      { re =  c.real(); negIm = -c.imag(); return *this; }
    conjugate& operator+=(const complex<float>& c)
      { re += c.real(); negIm -= c.imag(); return *this; }
    conjugate& operator-=(const complex<float>& c)
      { re -= c.real(); negIm += c.imag(); return *this; }

    // It is pleasant to note that we can self-multiply by either a complex or
    // a conjugate (leaving a conjugate result) in six flops which is the same
    // cost as an ordinary complex multiply:
    //    cplx=cplx*cplx: (a+bi)(r+si) = (ar-bs)+(as+br)i
    //    conj=conj*conj: (a-bi)(r-si) = (ar-bs)-(as+br)i
    //    conj=conj*cplx: (a-bi)(r+si) = (ar+bs)-(br-as)i
    conjugate& operator*=(const conjugate<float>& c) {
        const float r=(re*c.re - negIm*c.negIm);
        negIm=(re*c.negIm + negIm*c.re); re=r; return *this;
    }
    conjugate& operator*=(const complex<float>& t) {
        const float r=(re*t.real() + negIm*t.imag()); 
        negIm=(negIm*t.real() - re*t.imag()); re=r; return *this;
    }

    // Complex divide is messy and slow anyway so we'll convert to complex and back here,
    // making use of the fact that for complex c and d, c/d=conj(conj(c)/conj(d)).
    conjugate& operator/=(const conjugate<float>& d) {
        const complex<float> t = conj()/d.conj();
        re = t.real(); negIm = t.imag(); // conjugating!
        return *this;
    }
    conjugate& operator/=(const complex<float>& d) {
        const complex<float> t = conj()/std::conj(d);
        re = t.real(); negIm = t.imag(); // conjugating!
        return *this;
    }

    const float&               real() const { return re; }
    float&                     real()       { return re; }

    const negator<float>&      imag() const { return reinterpret_cast<const negator<float>&>(negIm); }
    negator<float>&            imag()       { return reinterpret_cast<negator<float>&>(negIm); }

    const complex<float>& conj() const { return reinterpret_cast<const complex<float>&>(*this); }
    complex<float>&       conj()       { return reinterpret_cast<complex<float>&>(*this); }

    // Special conjugate methods of use primarily in operator implementations.
    const float& negImag() const { return negIm; }
    float&       negImag()       { return negIm; }
    bool         isReal()  const { return negIm==0.f; }

private:
    float re;   // The value represented here is re - negIm*i.
    float negIm;
};




//////////////////////////////////////////
// Specialization for conjugate<double> //
//////////////////////////////////////////

template <>  class conjugate<double> {
public:
    conjugate() {
    #ifndef NDEBUG
        re = negIm = std::numeric_limits<double>::quiet_NaN();
    #endif  
    }
    // default copy constructor, copy assignment, destructor

    /// Construction from reals. Note that the numeric result is (real-imag*i).
    conjugate(const double& real, const double& imag) { re = real; negIm = imag; }
    conjugate(const double& real, int i) { re = real; negIm = double(i); }
    conjugate(int r, const double& imag) { re = double(r); negIm = imag; }
    conjugate(int r, int i) { re = double(r); negIm = double(i); }

    /// Implicit conversion from double to conjugate<double>.
    conjugate(const double& real) { re = real; negIm = 0.; }
    conjugate(int r) { re = double(r); negIm = 0.; }

    // Implicit conversions from float are allowed since
    // there is no loss in going to double.
    conjugate(const conjugate<float>& cf)
      { re = double(cf.real()); negIm = double(cf.negImag()); }
    conjugate(const float& rf)
      { re = double(rf); negIm = 0.; }

    // Conversions from complex are always explicit. Note that the value
    // represented by the conjugate must be identical to that represented by
    // the complex, which means we must negate the imaginary part.
    explicit conjugate(const complex<float>& x)
      { re = double(x.real()); negIm = double(-x.imag()); }
    explicit conjugate(const complex<double>& x)
      { re = x.real(); negIm = -x.imag(); }

    /// Implicit conversion to complex<double> when necessary
    /// (costs an actual negation -- yuck!).
    operator complex<double>() const
      { return complex<double>(re,-negIm); } 
    
    // Can't defer here by casting to negator<conjugate> -- this must act
    // like a built-in. But ... we can use this as a chance to convert
    // to complex and save one negation.
    complex<double> operator-() const { return complex<double>(-re,negIm); }

    // Useless.
    const conjugate& operator+() const { return *this; }

    // Computed assignment operators. We don't depend on implicit conversions
    // from reals to conjugates here because we can save a few flops by handling
    // the reals explicitly. Note that we only provide operators for implicitly
    // convertible precisions, though, which in this case means floats and doubles.
    conjugate& operator=(const double& r)
      { re = r; negIm = 0.; return *this; }
    conjugate& operator+=(const double& r)
      { re += r; return *this; }
    conjugate& operator-=(const double& r)
      { re -= r; return *this; }
    conjugate& operator*=(const double& r)
      { re *= r; negIm *= r; return *this; }
    conjugate& operator/=(const double& r)
      { re /= r; negIm /= r; return *this; }

    conjugate& operator=(const float& r)
      { re = r; negIm = 0.; return *this; }
    conjugate& operator+=(const float& r)
      { re += r; return *this; }
    conjugate& operator-=(const float& r)
      { re -= r; return *this; }
    conjugate& operator*=(const float& r)
      { re *= r; negIm *= r; return *this; }
    conjugate& operator/=(const float& r)
      { re /= r; negIm /= r; return *this; }

    // Disambiguate int to be a double.
    conjugate& operator =(int i) {*this =(double)i; return *this;}
    conjugate& operator+=(int i) {*this+=(double)i; return *this;}
    conjugate& operator-=(int i) {*this-=(double)i; return *this;}
    conjugate& operator*=(int i) {*this*=(double)i; return *this;}
    conjugate& operator/=(int i) {*this/=(double)i; return *this;}

    conjugate& operator+=(const conjugate<double>& c)
      { re += c.re; negIm += c.negIm; return *this; }
    conjugate& operator-=(const conjugate<double>& c)
      { re -= c.re; negIm -= c.negIm; return *this; }

    conjugate& operator+=(const conjugate<float>& c)
      { re += c.real(); negIm += c.negImag(); return *this; }
    conjugate& operator-=(const conjugate<float>& c)
      { re -= c.real(); negIm -= c.negImag(); return *this; }

    conjugate& operator=(const complex<double>& c)
      { re =  c.real(); negIm = -c.imag(); return *this; }
    conjugate& operator+=(const complex<double>& c)
      { re += c.real(); negIm -= c.imag(); return *this; }
    conjugate& operator-=(const complex<double>& c)
      { re -= c.real(); negIm += c.imag(); return *this; }

    conjugate& operator=(const complex<float>& c)
      { re =  c.real(); negIm = -c.imag(); return *this; }
    conjugate& operator+=(const complex<float>& c)
      { re += c.real(); negIm -= c.imag(); return *this; }
    conjugate& operator-=(const complex<float>& c)
      { re -= c.real(); negIm += c.imag(); return *this; }

    // It is pleasant to note that we can self-multiply by either a complex or
    // a conjugate (leaving a conjugate result) in six flops which is the same
    // cost as an ordinary complex multiply:
    //    cplx=cplx*cplx: (a+bi)(r+si) = (ar-bs)+(as+br)i
    //    conj=conj*conj: (a-bi)(r-si) = (ar-bs)-(as+br)i
    //    conj=conj*cplx: (a-bi)(r+si) = (ar+bs)-(br-as)i
    conjugate& operator*=(const conjugate<double>& c) {
        const double r=(re*c.re - negIm*c.negIm);
        negIm=(re*c.negIm + negIm*c.re); re=r; return *this;
    }
    conjugate& operator*=(const complex<double>& t) {
        const double r=(re*t.real() + negIm*t.imag()); 
        negIm=(negIm*t.real() - re*t.imag()); re=r; return *this;
    }

    conjugate& operator*=(const conjugate<float>& c)     { return operator*=(conjugate<double>(c)); }
    conjugate& operator*=(const complex<float>& c)  { return operator*=(complex<double>(c)); }

    // Complex divide is messy and slow anyway so we'll convert to complex and back here,
    // making use of the fact that for complex c and d, c/d=conj(conj(c)/conj(d)).
    conjugate& operator/=(const conjugate<double>& d) {
        const complex<double> t = conj()/d.conj();
        re = t.real(); negIm = t.imag(); // conjugating!
        return *this;
    }
    conjugate& operator/=(const complex<double>& d) {
        const complex<double> t = conj()/std::conj(d);
        re = t.real(); negIm = t.imag(); // conjugating!
        return *this;
    }

    conjugate& operator/=(const conjugate<float>& c)     { return operator/=(conjugate<double>(c)); }
    conjugate& operator/=(const complex<float>& c)  { return operator/=(complex<double>(c)); }

    const double&               real() const { return re; }
    double&                     real()       { return re; }

    const negator<double>&      imag() const { return reinterpret_cast<const negator<double>&>(negIm); }
    negator<double>&            imag()       { return reinterpret_cast<negator<double>&>(negIm); }

    const complex<double>& conj() const { return reinterpret_cast<const complex<double>&>(*this); }
    complex<double>&       conj()       { return reinterpret_cast<complex<double>&>(*this); }

    // Special conjugate methods of use primarily in operator implementations.
    const double& negImag() const { return negIm; }
    double&       negImag()       { return negIm; }
    bool          isReal()  const { return negIm==0.; }

private:
    double re;   // The value represented here is re - negIm*i.
    double negIm;
};



// These definitions had to be deferred until all the specializations have been declared.
conjugate<float>::conjugate(const conjugate<double>& cd) { 
    re = float(cd.real()); negIm = float(cd.negImag());
}

// Global functions real(),imag(), conj(), abs(), and norm() are overloaded here
// for efficiency (e.g., abs(c)==abs(conj(c))). The others still work through
// the implicit conversion from conjugate<T> to complex<T>, which costs
// one negation. The non-overloaded functions defined for complex are
// arg(), which returns a real, and a set of functions which return complex:
// polar,cos,cosh,exp,log,log10,pow,sin,sinh,sqrt,tan,tanh.
inline const float&               real(const conjugate<float>& c) { return c.real(); }
inline const negator<float>&      imag(const conjugate<float>& c) { return c.imag(); }
inline const complex<float>& conj(const conjugate<float>& c) { return c.conj(); }
inline float abs (const conjugate<float>& c) { return std::abs(c.conj()); }
inline float norm(const conjugate<float>& c) { return std::norm(c.conj()); }

inline const double&               real(const conjugate<double>& c) { return c.real(); }
inline const negator<double>&      imag(const conjugate<double>& c) { return c.imag(); }
inline const complex<double>& conj(const conjugate<double>& c) { return c.conj(); }
inline double abs (const conjugate<double>& c) { return std::abs(c.conj()); }
inline double norm(const conjugate<double>& c) { return std::norm(c.conj()); }






// Binary operators with conjugate as one of the operands, and the other any
// numerical type (real, complex, conjugate) but NOT a negator type. Each operator
// will silently work with operands of mixed precision, widening the result as
// necessary. We try to return complex rather than conjugate whenever possible.

template <class R, class CHAR, class TRAITS> inline std::basic_istream<CHAR,TRAITS>&
operator>>(std::basic_istream<CHAR,TRAITS>& is, conjugate<R>& c) {
    complex<R> z; is >> z; c=z;
    return is;
}
template <class R, class CHAR, class TRAITS> inline std::basic_ostream<CHAR,TRAITS>&
operator<<(std::basic_ostream<CHAR,TRAITS>& os, const conjugate<R>& c) {
    return os << complex<R>(c);
}

// Operators involving only conjugate and complex can be templatized reliably, as can
// operators which do not mix precision. But we have to deal explicitly with mixes
// of conjugate<R> and some other real type S, because the 'class S' template
// argument can match anything and create ambiguities.

// conjugate<R> with float, double. With 'float' we can be sure that R
// is the right width for the return value. 'double' is trickier and we have
// to use the Wider<R,...> helper class to give us the right return type.

// Commutative ops need be done only once: +, *, ==, and != is defined in terms of ==.

// conjugate = conjugate + real
template <class R> inline conjugate<R>                    operator+(const conjugate<R>& a, const float&       b)
  { return conjugate<R>(a) += b; }
template <class R> inline typename Wider<R,double>::WConj operator+(const conjugate<R>& a, const double&      b)
  { return typename Wider<R,double>::WConj(a) += b; }

// conjugate = real + conjugate
template <class R> inline conjugate<R>                    operator+(const float&       a, const conjugate<R>& b) {return b+a;}
template <class R> inline typename Wider<R,double>::WConj operator+(const double&      a, const conjugate<R>& b) {return b+a;}

// conjugate = conjugate * real
template <class R> inline conjugate<R>                    operator*(const conjugate<R>& a, const float&       b)
  { return conjugate<R>(a) *= b; }
template <class R> inline typename Wider<R,double>::WConj operator*(const conjugate<R>& a, const double&      b)
  { return typename Wider<R,double>::WConj(a) *= b; }

// conjugate = real * conjugate
template <class R> inline conjugate<R>                    operator*(const float&       a, const conjugate<R>& b) {return b*a;}
template <class R> inline typename Wider<R,double>::WConj operator*(const double&      a, const conjugate<R>& b) {return b*a;}

// bool = conjugate==real
template <class R> inline bool                            operator==(const conjugate<R>& a, const float&       b)
  { return a.isReal() && a.real()==b; }
template <class R> inline bool                            operator==(const conjugate<R>& a, const double&      b)
  { return a.isReal() && a.real()==b; }

// bool = real==conjugate, bool = conjugate!=real, bool = real!=conjugate 
template <class R> inline bool operator==(const float&        a, const conjugate<R>& b) {return b==a;}
template <class R> inline bool operator==(const double&       a, const conjugate<R>& b) {return b==a;}
template <class R> inline bool operator!=(const conjugate<R>& a, const float&        b) {return !(a==b);}
template <class R> inline bool operator!=(const conjugate<R>& a, const double&       b) {return !(a==b);}
template <class R> inline bool operator!=(const float&        a, const conjugate<R>& b) {return !(a==b);}
template <class R> inline bool operator!=(const double&       a, const conjugate<R>& b) {return !(a==b);}

// Non-commutative ops are a little messier.

// conjugate = conjugate - real
template <class R> inline conjugate<R>                    operator-(const conjugate<R>& a, const float&       b)
  { return conjugate<R>(a) -= b; }
template <class R> inline typename Wider<R,double>::WConj operator-(const conjugate<R>& a, const double&      b)
  { return typename Wider<R,double>::WConj(a) -= b; }

// complex = real - conjugate 
// This is nice because -conjugate.imag() is free.
template <class R> inline complex<R>                      operator-(const float&       a, const conjugate<R>& b)
  { return complex<R>(a-b.real(), -b.imag()); }
template <class R> inline typename Wider<R,double>::WCplx operator-(const double&      a, const conjugate<R>& b)
  { return typename Wider<R,double>::WCplx(a-b.real(), -b.imag()); }

// conjugate = conjugate / real
template <class R> inline conjugate<R>                    operator/(const conjugate<R>& a, const float&       b)
  { return conjugate<R>(a) /= b; }
template <class R> inline typename Wider<R,double>::WConj operator/(const conjugate<R>& a, const double&      b)
  { return typename Wider<R,double>::WConj(a) /= b; }

// complex = real / conjugate 
// Division by complex is tricky and slow anyway so we'll just convert to complex
// at the cost of one negation.
template <class R> inline complex<R>                      operator/(const float&       a, const conjugate<R>& b)
  { return (R)a/complex<R>(b); }
template <class R> inline typename Wider<R,double>::WCplx operator/(const double&      a, const conjugate<R>& b)
  { return (typename Wider<R,double>::WReal)a/(typename Wider<R,double>::WCplx(b)); }


// That's it for (conjugate, real) combinations. Now we need to do all the (conjugate, conjugate) and
// (conjugate, complex) combinations which are safer to templatize fully. There are many more opportunities
// to return complex rather than conjugate here. Keep in mind here that for a conjugate number c=a-bi,
// a=c.real() and b=-c.imag() are available for free. conjugate provides negImag() to return b directly.
// Below we'll capitalize components that should be accessed with negImag().

// conjugate = conjugate + conjugate: (a-Bi)+(r-Si) = (a+r)-(B+S)i
template <class R, class S> inline typename Wider<R,S>::WConj
operator+(const conjugate<R>& a, const conjugate<S>& r) {
    return typename Wider<R,S>::WConj(a.real()+r.real(), 
                                      a.negImag()+r.negImag());
}

// complex = conjugate + complex = complex + conjugate: (a-Bi)+(r+si) = (a+r)+(s-B)i
template <class R, class S> inline typename Wider<R,S>::WCplx
operator+(const conjugate<R>& a, const complex<S>& r) {
    return typename Wider<R,S>::WCplx(a.real()+r.real(), 
                                      r.imag()-a.negImag());
}
template <class R, class S> inline typename Wider<R,S>::WCplx
operator+(const complex<R>& a, const conjugate<S>& r) { return r+a; }

// complex = conjugate - conjugate: (a-Bi)-(r-Si) = (a-r)+(S-B)i
template <class R, class S> inline typename Wider<R,S>::WCplx
operator-(const conjugate<R>& a, const conjugate<S>& r) {
    return typename Wider<R,S>::WCplx(a.real()-r.real(),
                                      r.negImag()-a.negImag());
}

// Neg<complex> = conjugate - complex:   (a-Bi)-(r+si) = (a-r)+(-B-s)i = -[(r-a)+(B+s)i]
template <class R, class S> inline negator<typename Wider<R,S>::WCplx>
operator-(const conjugate<R>& a, const complex<S>& r) {
    return negator<typename Wider<R,S>::WCplx>::recast
             (typename Wider<R,S>::WCplx(r.real()-a.real(), 
                                         a.negImag()+r.imag()));
}

// complex = complex - conjugate: (a+bi)-(r-Si) = (a-r)+(b+S)i
template <class R, class S> inline typename Wider<R,S>::WCplx
operator-(const complex<R>& a, const conjugate<S>& r) {
    return typename Wider<R,S>::WCplx(a.real()-r.real(),
                                      a.imag()+r.negImag());
}

// We can multiply by either a complex or a conjugate (leaving a complex 
// or negated complex result) in six flops which is the same cost as an
// ordinary complex multiply.
//        (cplx=cplx*cplx: (a+bi)(r+si) =   (ar-bs)+(as+br)i)

// Neg<cplx>=conj*conj: (a-Bi)(r-Si) = -[(BS-ar)+(aS+Br)i]
template <class R, class S> inline negator<typename Wider<R,S>::WCplx>
operator*(const conjugate<R>& a, const conjugate<S>& r) {
    return negator<typename Wider<R,S>::WCplx>::recast
           (typename Wider<R,S>::WCplx(a.negImag()*r.negImag() - a.real()*r.real(),
                                       a.real()*r.negImag()    + a.negImag()*r.real()));
}

// cplx=conj*cplx: (a-Bi)(r+si) = (ar+Bs)+(as-Br)i
template <class R, class S> inline typename Wider<R,S>::WCplx
operator*(const conjugate<R>& a, const complex<S>& r) {
    return typename Wider<R,S>::WCplx(a.real()*r.real() + a.negImag()*r.imag(),
                                      a.real()*r.imag() - a.negImag()*r.real());
}

template <class R, class S> inline typename Wider<R,S>::WCplx
operator*(const complex<R>& a, const conjugate<S>& r)
  { return r*a; }

// If there's a negator on the complex number, move it to the conjugate
// one which will change to complex with just one flop; then this
// is an ordinary complex mutiply.
template <class R, class S> inline typename Wider<R,S>::WCplx
operator*(const negator< complex<R> >& a, const conjugate<S>& r)
  { return (-a)*(-r); } // -a is free here
template <class R, class S> inline typename Wider<R,S>::WCplx
operator*(const conjugate<R>& a, const negator< complex<S> >& r)
  { return (-a)*(-r); } // -r is free here

// Division is tricky and there is little to gain by trying to exploit
// the conjugate class, so we'll just convert to complex. (Remember that
// conj() is free for Conjugates.)

template <class R, class S> inline typename Wider<R,S>::WCplx
operator/(const conjugate<R>& a, const conjugate<S>& r) {
    return std::conj(a.conj()/r.conj());
}

template <class R, class S> inline typename Wider<R,S>::WCplx
operator/(const conjugate<R>& a, const complex<S>& r) {
    return std::conj(a.conj()/std::conj(r));
}

template <class R, class S> inline typename Wider<R,S>::WCplx
operator/(const complex<R>& a, const conjugate<S>& r) {
    return std::conj(std::conj(a)/r.conj());
}


template <class R, class S> inline bool
operator==(const conjugate<R>& a, const conjugate<S>& r) {
    return a.real() == r.real() && a.negImag() == r.negImag();
}

template <class R, class S> inline bool
operator==(const conjugate<R>& a, const complex<S>& r) {
    return a.real() == r.real() && -a.negImag() == r.imag();
}

template <class R, class S> inline bool
operator==(const complex<R>& a, const conjugate<S>& r) {return r==a;}

template <class R, class S> inline bool
operator!=(const conjugate<R>& a, const conjugate<S>& r)    {return !(a==r);}

template <class R, class S> inline bool
operator!=(const conjugate<R>& a, const complex<S>& r) {return !(a==r);}

template <class R, class S> inline bool
operator!=(const complex<R>& a, const conjugate<S>& r) {return !(a==r);}


} // namespace SimTK

#endif //SimTK_SIMMATRIX_CONJUGATE_H_
