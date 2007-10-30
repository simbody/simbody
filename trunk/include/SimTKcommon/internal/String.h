#ifndef SimTK_SimTKCOMMON_STRING_H_
#define SimTK_SimTKCOMMON_STRING_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif


#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/List.h"

#include <cstdio>
#include <string>
#include <limits>
#include <complex>

namespace SimTK {
	
/**
 * SimTK::String is just an std::string with some additional methods 
 * defined. It may be freely interspersed with std:strings.
 *
 * Binary compatibility note: while the std::string implementation, like
 * the other std classes, is fully exposed in its header file, we have 
 * determined that these classes are stable enough, and that compiler
 * implementators work hard to maintain cross-release stability, that we
 * will assume they are "sufficiently" binary compatible to use them in
 * the SimTK API. So std::string and SimTK::String can be passed through
 * the API interface, and will be binary compatible from release to release
 * of the SimTK Core, provided that the compilers used to compile client
 * and library code are binary compatible with respect to the std::string
 * implementation in the C++ Standard Template Library.
 */
class String : public std::string {
public:
    /// default constructor produces an empty string
	String() { }

    // uses default copy constructor, copy assignment, and destructor

    /// This is an implicit conversion from const char* to String.
	String(const char* s) : std::string(s) { }

    /// This is an implicit conversion from std::string to String
	String(const std::string& s) : std::string(s) { }

    /// Construct a String as a copy of a substring begining at 
    /// position \a start with length \a len.
    String(const String& s, int start, int len) : std::string(s,start,len) { }

    /// This is an implicit conversion from String to null-terminated 
    /// C-style string (array of chars).
    operator const char*() const { return c_str(); }

    /// Add operator[] that takes int index instead of size_type.
    char& operator[](int i) {
        assert(i >= 0);
        return std::string::operator[]((std::string::size_type)i);
    }

    /// Add operator[] that takes int index instead of size_type.
    char operator[](int i) const {
        assert(i >= 0);
        return std::string::operator[]((std::string::size_type)i);
    }

    /// Pass through to string::operator[].
    char& operator[](std::string::size_type i) {return std::string::operator[](i);}
    /// Pass through to string::operator[].
    char operator[](std::string::size_type i) const {return std::string::operator[](i);}

    /// @name Formatting constructors
    /// These contructors format the supplied argument into a String.
    //@{
	explicit String(int i) { char buf[32]; sprintf(buf,"%d",i); (*this)=buf; }
	explicit String(long i) { char buf[32]; sprintf(buf,"%ld",i); (*this)=buf; }
    explicit String(unsigned int s)  { char buf[32]; sprintf(buf,"%u",s); (*this)=buf; }
    explicit String(unsigned long s) { char buf[32]; sprintf(buf,"%lu",s); (*this)=buf; }
    explicit String(float r)	{ char buf[64]; sprintf(buf,"%.8g",r); (*this)=buf; }
	explicit String(double r)	{ char buf[64]; sprintf(buf,"%.16g",r); (*this)=buf; }
	explicit String(long double r)	{ char buf[64]; sprintf(buf,"%.20Lg",r); (*this)=buf; }
	explicit String(std::complex<float> r)	
		{ char buf[128]; sprintf(buf,"(%.8g,%.8g)",r.real(),r.imag()); (*this)=buf; }
	explicit String(std::complex<double> r)	
		{ char buf[128]; sprintf(buf,"(%.16g,%.16g)",r.real(),r.imag()); (*this)=buf; }
	explicit String(std::complex<long double> r)	
		{ char buf[128]; sprintf(buf,"(%.20Lg,%.20Lg)",r.real(),r.imag()); (*this)=buf; }
    explicit String(bool b) : std::string(b?"true":"false") { }
    //@}
};	

SimTK_LIST_SPECIALIZE(String);

}
#endif // SimTK_SimTKCOMMON_STRING_H_
