#ifndef SimTK_SimTKCOMMON_SERIALIZE_H_
#define SimTK_SimTKCOMMON_SERIALIZE_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
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
This file contains definitions of templatized serialize-to-stream methods 
specialized for the built-in C++ and SimTK low-level classes. These are used
to provide the concrete level of recursively-defined methods for serialization
of container classes (Array_, Vector_, Vec, Mat, etc.) each of which should
define templatized methods of the same names that serialize using the same
serialization method applied to their elements. Methods are provided for 
different kinds of formatting; at the concrete level these are likely to be
identical but they differ for containers. Unformatted output will write just
space separated elements, while formatted output may write surrounding
parentheses or brackets, comma separators, or whatever else might be 
appropriate.
**/

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/String.h"

#include <iostream>

namespace SimTK {

//------------------------------------------------------------------------------
//                           WRITE UNFORMATTED
//------------------------------------------------------------------------------
/**
@defgroup writeUnformatted        writeUnformatted()
@ingroup Serialization

Namespace-scope utility method SimTK::writeUnformatted\<T>() writes a value of
type T to an output stream as a space-separated series of tokens with no 
brackets, commas, semicolons or other formatting characters. If there is only 
one token it is serialized without any added leading or trailing whitespace.
Output of bool values is "true" or "false"; output of non-finite floating point
values is "NaN", "Inf", or "-Inf" as in Matlab. For matrix output, we will 
write a newline between each row, meaning there won't be one after the last 
row. This method is specialized for many SimTK types, and you can specialize
it for your own types in which case containers like Array_\<YourType> will 
work correctly too. **/
/**@{**/

/** The default implementation of writeUnformatted\<T> converts the object to
a String using the templatized String constructor, and then writes that
string to the stream using String::operator<<(). This is suitable for use
with any of the built-in types. Note that bool will be output "true" or
"false" and non-finite floating point values are written as NaN, Inf,
or -Inf as appropriate. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const T& v) {
    o << String(v);
}

/** Specialize for float to help some compilers with template matching. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const float& v)  
{   writeUnformatted<float>(o,v); }
/** Specialize for double to help some compilers with template matching. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const double& v) 
{   writeUnformatted<double>(o,v); }
/** Specialize for long double to help some compilers with template 
matching. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const long double& v) 
{   writeUnformatted<long double>(o,v); }

/** Specialize for SimTK::negator\<T>: convert to T and write. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const negator<T>& v) 
{   writeUnformatted(o, T(v)); }

/** Specialize for std::complex\<T>: just write two T's separated by a space;
no parentheses or comma. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const std::complex<T>& v) 
{   writeUnformatted(o, v.real()); o << " "; writeUnformatted(o, v.imag()); }

/** Specialize for SimTK::conjugate\<T>: same as std::complex\<T>. **/
template <class T> inline void
writeUnformatted(std::ostream& o, const conjugate<T>& v) 
{   writeUnformatted(o, std::complex<T>(v)); }


/** Specialize for Vec<M,E,S> to delegate to element type E, with spaces
separating the elements. 
@relates SimTK::Vec **/
template <int M, class E, int S> inline void
writeUnformatted(std::ostream& o, const Vec<M,E,S>& v) {
    for (int i=0; i < M; ++i) {
        if (i != 0) o << " ";
        writeUnformatted(o, v[i]);
    }
}   
/** Specialize for Row<N,E,S> to delegate to element type E, with spaces
separating the elements; raw output is same as Vec. 
@relates SimTK::Row **/
template <int N, class E, int S> inline void
writeUnformatted(std::ostream& o, const Row<N,E,S>& v) 
{   writeUnformatted(o, ~v); }

/** Specialize for Mat<M,N,E,CS,RS> delegating to Row<N,E,RS> with newlines
separating the rows, but no final newline.
@relates SimTK::Mat **/
template <int M, int N, class E, int CS, int RS> inline void
writeUnformatted(std::ostream& o, const Mat<M,N,E,CS,RS>& v) {
    for (int i=0; i < M; ++i) {
        if (i != 0) o << std::endl;
        writeUnformatted(o, v[i]);
    }
} 

/** Specialize for SymMat<M,E,RS> delegating to Row<M,E,RS> with newlines
separating the rows, but no final newline.
@relates SimTK::SymMat **/
template <int M, class E, int RS> inline void
writeUnformatted(std::ostream& o, const SymMat<M,E,RS>& v) {
    for (int i=0; i < M; ++i) {
        if (i != 0) o << std::endl;
        writeUnformatted(o, v[i]);
    }
} 
/**@}**/


//------------------------------------------------------------------------------
//                             READ UNFORMATTED
//------------------------------------------------------------------------------
/**
@defgroup readFromStream         readUnformatted()
@ingroup Serialization

Namespace-scope utility function SimTK::readUnformatted\<T>() reads the next
whitespace-separated token(s) from a given stream and attempts to interpret 
them as a value of type T. **/
/**@{**/


/** Read in the next whitespace-delimited token as a String, ignoring
leading whitespace. The token is terminated by whitespace or eof and the
terminating character is left in the stream. Failure to find any non-whitespace
characters is treated as a failure. If we're successful at getting a token the 
function returns true; the stream may be good() or eof(). Otherwise this 
function returns false and the fail bit in the stream is also set. **/
inline bool
readOneTokenUnformatted(std::istream& in, String& token) {
    // If the stream is already bad or at eof, we fail.
    if (!in.good()) {in.setstate(std::ios::failbit); return false;}
    // Skip whitespace. Failure or eof here means no token.
    std::ws(in);
    if (!in.good()) {in.setstate(std::ios::failbit); return false;}
    in >> token; if (in.fail()) return false;
    if (token.empty()) {in.setstate(std::ios_base::failbit); return false;}
    return true;
} 

/** The default implementation of readUnformatted\<T> reads in the next
whitespace-separated token and then attempts to convert the whole thing into
one value of type T. If that fails, or if the whole token isn't consumed, we
return false and the fail bit is set in the stream.
@returns true if we are successful in reading the value **/
template <class T> inline bool
readUnformatted(std::istream& in, T& v) {
    String token;
    if (!readOneTokenUnformatted(in, token)) return false;
    if (!token.tryConvertTo<T>(v)) 
    {   in.setstate(std::ios::failbit); return false; }
    return true;
}

/** Specialize for float to help some compilers with template matching. **/
template <class T> inline bool
readUnformatted(std::istream& in, float& v) 
{   return readUnformatted<float>(in,v); }
/** Specialize for double to help some compilers with template matching. **/
template <class T> inline bool
readUnformatted(std::istream& in, double& v) 
{   return readUnformatted<double>(in,v); }
/** Specialize for long double to help some compilers with template 
matching. **/
template <class T> inline bool
readUnformatted(std::istream& in, long double& v) 
{   return readUnformatted<long double>(in,v); }

/** Specialization for negator<T>: read as type T and convert. **/
template <class T> inline bool
readUnformatted(std::istream& in, negator<T>& v) {
    T nonneg; if (!readUnformatted(in, nonneg)) return false;
    v = nonneg; // Changes representation, not value.
    return true;
}

/** Specialization for std::complex<T> (two space-separated T's). **/
template <class T> inline bool
readUnformatted(std::istream& in, std::complex<T>& v) {
    T real, imag;
    if (!readUnformatted(in, real)) return false;
    if (!readUnformatted(in, imag)) return false;
    v = std::complex<T>(real,imag);
    return true;
}

/** Specialization for SimTK::conjugate<T> (same as std::complex<T>). **/
template <class T> inline bool
readUnformatted(std::istream& in, conjugate<T>& v) {
    std::complex<T> cmplx; if (!readUnformatted(in, cmplx)) return false;
    v = cmplx; // Changes representation, not value.
    return true;
}

/** Specialization for SimTK::String (just read token). **/
template <> inline bool 
readUnformatted<String>(std::istream& in, String& v)
{   return readOneTokenUnformatted(in,v); }


/** Specialize for Vec<M,E,S> to delegate to element type E, with spaces
separating the elements. 
@relates SimTK::Vec **/
template <int M, class E, int S> inline bool
readUnformatted(std::istream& in, Vec<M,E,S>& v) {
    for (int i=0; i < M; ++i)
        if (!readUnformatted(in, v[i])) return false;
    return true;
}   

/** Specialize for Row<N,E,S> to delegate to element type E, with spaces
separating the elements; format is same as Vec. 
@relates SimTK::Row **/
template <int N, class E, int S> inline bool
readUnformatted(std::istream& in, Row<N,E,S>& v) 
{   return readUnformatted(in, ~v); }

/** Specialize for Mat<M,N,E,CS,RS> delegating to Row<N,E,RS>.
@relates SimTK::Mat **/
template <int M, int N, class E, int CS, int RS> inline bool
readUnformatted(std::istream& in, Mat<M,N,E,CS,RS>& v) {
    for (int i=0; i < M; ++i)
        if (!readUnformatted(in, v[i])) return false;
    return true;
} 

/** Specialize for SymMat<M,E,RS>. We require the whole matrix, then
verify symmetry and assign to the output argument. We'll return false with
the fail bit set in the stream if we get a Mat\<M,M> but it isn't 
symmetric to a tight numerical tolerance.
@see SimTK::SymMat::setFromSymmetric(), SimTK::Mat::isNumericallySymmetric()
@relates SimTK::SymMat **/
template <int M, class E, int RS> inline bool
readUnformatted(std::istream& in, SymMat<M,E,RS>& v) {
    Mat<M,M,E> m; if (!readUnformatted(in, m)) return false;
    if (!m.isNumericallySymmetric()) {
        in.setstate(std::ios::failbit); 
        return false;
    }
    v.setFromSymmetric(m);
    return true;
} 
/**@}**/

//------------------------------------------------------------------------------
//                             WRITE FORMATTED
//------------------------------------------------------------------------------
/**
@defgroup writeFormatted        writeFormatted()
@ingroup Serialization

Namespace-scope utility method SimTK::writeFormatted\<T>() writes a value of
type T to an output stream in a mildly formatted way depending on the type
of container represented by T. For the basic types this is the same as
what you get from the templatized String<T>() constructor when given a value
of type T. For containers like Array_ or Vector_, this method is specialized as
part of the definition of that container. Usually the containers will surround
their contents with some kind of brackets, and separate elements by commas.

Anything written by writeFormatted() can be read back in with readFormatted().
**/
/**@{**/

/** The default implementation of writeFormatted\<T> converts the object to
a String using the templatized String constructor, and then writes that
string to the stream using String::operator<<(). This is suitable for use
with any of the built-in and scalar SimTK types. Note that bool will be output 
"true" or "false" and non-finite floating point values are written as NaN, 
Inf, or -Inf as appropriate. Complex numbers will serialize as (real,imag). **/
template <class T> inline void
writeFormatted(std::ostream& o, const T& v) {
    o << String(v);
}
/**@}**/

//------------------------------------------------------------------------------
//                             READ FORMATTED
//------------------------------------------------------------------------------
/**
@defgroup readFormatted        readFormatted()
@ingroup Serialization

Namespace-scope utility method SimTK::readFormatted\<T>() reads a value of
type T from an input stream, recognizing a variety of formats such as the
format produced by writeFormatted() or simpler unformatted data such as is
produced by writeUnformatted(). When T is a container, the recognized 
formats depend on the readFormatted() specialization provided with as part of
the definition of that container. **/
/**@{**/

/** The default implementation of readFormatted\<T>() uses 
readUnformatted\<T>(). **/
template <class T> inline bool
readFormatted(std::istream& in, T& v) {
    return readUnformatted(in, v);
}

// TODO: need specializations for complex to support (real,imag) where the
// numbers can be NaN, Inf, -Inf.
/**@}**/



} // namespace SimTK
#endif // SimTK_SimTKCOMMON_SERIALIZE_H_
