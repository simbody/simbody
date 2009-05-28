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
 * Portions copyright (c) 2005-9 Stanford University and the Authors.         *
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
    /// Default constructor produces an empty string.
	String() { }

    // uses default copy constructor, copy assignment, and destructor

    /// This is an implicit conversion from const char* to String.
	String(const char* s) : std::string(s) { }

    /// We allow creating a String from a char but you have to do
    /// it explicitly.
    explicit String(char c) {push_back(c);}

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

    /// Override std::string size() method to return an int instead of the 
    /// inconvenient unsigned type size_type.
    int size() const {return (int)std::string::size();}

    /// Override std::string length() method to return an int instead of the 
    /// inconvenient unsigned type size_type.
    int length() const {return (int)std::string::length();}

    /// @name Formatting constructors
    /// These contructors format the supplied argument into a String which is 
    /// suitable to use with the "<<" stream operator.
    //@{
    /// Format an int as a printable String.
	explicit String(int i) { char buf[32]; sprintf(buf,"%d",i); (*this)=buf; }
    /// Format a long as a printable String.
	explicit String(long i) { char buf[32]; sprintf(buf,"%ld",i); (*this)=buf; }
    /// Format a long long as a printable String.
	explicit String(long long i) { char buf[64]; sprintf(buf,"%lld",i); (*this)=buf; }
    /// Format an unsigned int as a printable String.
    explicit String(unsigned int s)  { char buf[32]; sprintf(buf,"%u",s); (*this)=buf; }
    /// Format an unsigned long as a printable String.
    explicit String(unsigned long s) { char buf[32]; sprintf(buf,"%lu",s); (*this)=buf; }
    /// Format an unsigned long long as a printable String.
    explicit String(unsigned long long s) { char buf[64]; sprintf(buf,"%llu",s); (*this)=buf; }
    /// Format a float as a printable String.
    explicit String(float r)	{ char buf[64]; sprintf(buf,"%.8g",r); (*this)=buf; }
    /// Format a double as a printable String.
	explicit String(double r)	{ char buf[64]; sprintf(buf,"%.16g",r); (*this)=buf; }
    /// Format a long double as a printable String.
	explicit String(long double r)	{ char buf[64]; sprintf(buf,"%.20Lg",r); (*this)=buf; }
    /// Format a complex<float> as a printable String (real,imag).
	explicit String(std::complex<float> r)	
		{ char buf[128]; sprintf(buf,"(%.8g,%.8g)",r.real(),r.imag()); (*this)=buf; }
    /// Format a complex<double> as a printable String (real,imag).
	explicit String(std::complex<double> r)	
		{ char buf[128]; sprintf(buf,"(%.16g,%.16g)",r.real(),r.imag()); (*this)=buf; }
    /// Format a complex<long double> as a printable String (real,imag).
	explicit String(std::complex<long double> r)	
		{ char buf[128]; sprintf(buf,"(%.20Lg,%.20Lg)",r.real(),r.imag()); (*this)=buf; }
    /// Format a bool as a printable String "true" or "false".
    explicit String(bool b) : std::string(b?"true":"false") { }
    //@}

    /// @name In-place modifications
    /// These are member functions which add to the existing std::string functionality.
    /// These methods return a reference to "this" String, so may be chained like 
    /// assignment statements. If you would like to use these on an std::string,
    /// use the String::updAs() method to recast the std::string to a String.
    /// Note that there is also an equivalent set of static methods which return 
    /// a new String rather than changing the original.
    //@{
    /// Upshift the given String in place, so that lowercase letters are replaced
    /// with their uppercase equivalents as defined by std::toupper().
    String& toUpper();
    /// Downshift the given String in place, so that uppercase letters are replaced
    /// with their lowercase equivalents as defined by std::tolower().
    String& toLower();
    /// Trim this String in place, removing all the initial leading and trailing 
    /// white space, as defined by std::isspace() which typically includes space, 
    /// tab (\\t), newline (\\n), return (\\r),  and form feed (\\f).
    String& trimWhiteSpace();
    /// Substitute in place \a newChar for \a oldChar wherever \a oldChar appears
    /// in this String.
    String& replaceAllChar(char oldChar, char newChar);
    //@}


    /// @name Utility methods
    /// These static methods operate on SimTK::String or std::string objects and return
    /// SimTK::String objects (which are also std::string objects).
    //@{
    /// Cast an std::string to a SimTK::String without copying; subsequent changes
    /// to the std::string will affect the SimTK::String too since it is just a 
    /// reference to the original std::string.
    static const String& getAs(const std::string& s) 
    {   return reinterpret_cast<const String&>(s); }
    /// Cast a non-const std::string to a non-const SimTK::String without copying;
    /// changes made to the SimTK::String will affect the original std::string and
    /// vice versa.
    static String& updAs(std::string& s) 
    {   return reinterpret_cast<String&>(s); }
    /// Upshift the given std::string returning a new SimTK::String in which all
    /// the letters have been made upper case with toupper().
    static String toUpper(const std::string& in)
    {   return String(in).toUpper(); }
    /// Downshift the given std::string returning a new SimTK::String in which all
    /// the letters have be made lower case with tolower().
    static String toLower(const std::string& in)
    {   return String(in).toLower(); }
    /// Copy the input std::string to a new SimTK::String leaving off all the
    /// initial leading and trailing white space, as defined by isspace() which
    /// typically includes space, tab (\\t), newline (\\n), return (\\r), 
    /// and form feed (\\f).
    static String trimWhiteSpace(const std::string& in);
    /// Copy the input std::string to a new SimTK::String while substituting
    /// \a newChar for \a oldChar wherever \a oldChar appears in the input.
    String& replaceAllChar(const std::string& in, char oldChar, char newChar)
    {   return String(in).replaceAllChar(oldChar, newChar); }
    //@}
};	

}
#endif // SimTK_SimTKCOMMON_STRING_H_
