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
 * Portions copyright (c) 2005-10 Stanford University and the Authors.        *
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
#include <sstream>

namespace SimTK {
	
/** SimTK::String is a plug-compatible std::string replacement (plus some
additional functionality) intended to be suitable for passing through the 
SimTK API without introducing binary compatibility problems the way 
std::string does, especially on Windows. You can work in your own code with 
std::strings which will be quietly converted to and from SimTK::Strings when 
invoking SimTK API methods. Or, you can use SimTK::Strings and still pass them 
to standard library or other methods that are expecting std::strings, usually 
transparently. The SimTK::Array_<T> class is used similarly to avoid binary
compatibility problems that arise with std::vector<T>.

@todo Currently this is just derived from std::string and inherits all the
binary compatibility issues. Use it now anyway and you'll pick up the 
compatibility benefits later, when we get the time ...

@see SimTK::Array_  **/
class String : public std::string {
public:

/** Default constructor produces an empty string. **/
String() { }

// uses default copy constructor, copy assignment, and destructor

/** This is an implicit conversion from const char* to String. **/
String(const char* s) : std::string(s) { }

/** We allow creating a String from a char but you have to do it explicitly. **/
explicit String(char c) {push_back(c);}

/** This is an implicit conversion from std::string to String **/
String(const std::string& s) : std::string(s) { }

/** Construct a String as a copy of a substring begining at position \a start 
with length \a len. **/
String(const String& s, int start, int len) : std::string(s,start,len) { }

/** This is an implicit conversion from String to null-terminated C-style 
string (array of chars). **/
operator const char*() const { return c_str(); }

/** Add operator[] that takes int index instead of size_type. **/
char& operator[](int i) {
    assert(i >= 0);
    return std::string::operator[]((std::string::size_type)i);
}

/** Add operator[] that takes int index instead of size_type. **/
char operator[](int i) const {
    assert(i >= 0);
    return std::string::operator[]((std::string::size_type)i);
}

/** Pass through to string::operator[]. **/
char& operator[](std::string::size_type i) {return std::string::operator[](i);}
/** Pass through to string::operator[]. **/
char operator[](std::string::size_type i) const {return std::string::operator[](i);}

/** Override std::string size() method to return an int instead of the 
inconvenient unsigned type size_type. **/
int size() const {return (int)std::string::size();}

/** Override std::string length() method to return an int instead of the 
inconvenient unsigned type size_type. **/
int length() const {return (int)std::string::length();}

/** @name             Formatted output constructors
These contructors format the supplied argument into a human-readable String. **/
/*@{*/
/** Format an int as a printable String. **/
explicit String(int i) { char buf[32]; sprintf(buf,"%d",i); (*this)=buf; }
/** Format a long as a printable String. **/
explicit String(long i) { char buf[32]; sprintf(buf,"%ld",i); (*this)=buf; }
/** Format a long long as a printable String. **/
explicit String(long long i) { char buf[64]; sprintf(buf,"%lld",i); (*this)=buf; }
/** Format an unsigned int as a printable String. **/
explicit String(unsigned int s)  { char buf[32]; sprintf(buf,"%u",s); (*this)=buf; }
/** Format an unsigned long as a printable String. **/
explicit String(unsigned long s) { char buf[32]; sprintf(buf,"%lu",s); (*this)=buf; }
/** Format an unsigned long long as a printable String. **/
explicit String(unsigned long long s) { char buf[64]; sprintf(buf,"%llu",s); (*this)=buf; }
/** Format a float as a printable String. **/
explicit String(float r)	{ char buf[64]; sprintf(buf,"%.8g",r); (*this)=buf; }
/** Format a double as a printable String. **/
explicit String(double r)	{ char buf[64]; sprintf(buf,"%.16g",r); (*this)=buf; }
/** Format a long double as a printable String. **/
explicit String(long double r)	{ char buf[64]; sprintf(buf,"%.20Lg",r); (*this)=buf; }
/** Format a complex<float> as a printable String (real,imag). **/
explicit String(std::complex<float> r)	
	{ char buf[128]; sprintf(buf,"(%.8g,%.8g)",r.real(),r.imag()); (*this)=buf; }
/** Format a complex<double> as a printable String (real,imag). **/
explicit String(std::complex<double> r)	
	{ char buf[128]; sprintf(buf,"(%.16g,%.16g)",r.real(),r.imag()); (*this)=buf; }
/** Format a complex<long double> as a printable String (real,imag). **/
explicit String(std::complex<long double> r)	
	{ char buf[128]; sprintf(buf,"(%.20Lg,%.20Lg)",r.real(),r.imag()); (*this)=buf; }
/** Format a bool as a printable String "true" or "false"; if you want "1"
or "0" cast the bool to an int first. **/
explicit String(bool b) : std::string(b?"true":"false") { }
/** For any type T for which there is no matching constructor, this templatized
constructor will format an object of type T into a String provided that there
is a stream insertion operator<<() available for type T. **/
template <class T> inline explicit String(const T& t); // see below
/*@}*/

/** @name             Formatted input from string
These templatized methods attempt to interpret the entire contents of
this String as a single object of type T (although T may itself be a 
container like an Array or Vector). It is an error if the String has
the wrong format for an object of this type, or if the entire String is
not consumed. The acceptable formatting is defined by type T based on what
it thinks is acceptable stream foramtting. Leading and trailing white space 
are ignored except when type T is itself a String or std::string in which case 
the white space is included in the result. It is not acceptable for type T
to be a pointer type. In particular if you want to convert a String to a null-
terminated C-style char*, use the standard c_str() method rather than any of 
these.
@see Related namespace-level static methods converStringTo<T>().
**/
/*@{*/
/** Attempt to convert this String to an object of type T, returning a status
value to indicate success or failure. We require that the whole string is 
consumed except possibly for some trailing white space. 
@tparam         T   
    A non-pointer type that supports extraction operator>>() from an istream. 
    You will get a compilation failure if you try to use this method for a 
    type T for which no extraction operator is available and a runtime error
    if T is a pointer type.
@param[out]     out
    The converted value if we were able to parse the string successfully
    (i.e., function return is true), otherwise the output value is 
    undefined.
@return true if we got what we're looking for, false if anything went
wrong including failure to consume the entire string. 
@see convertTo<T>() **/
template <class T> inline bool tryConvertTo(T& out) const; // see below

/** Convert this String to an object of type T using the tryConvertTo<T>()
method but throwing an error on failure rather than returning status. Using
this routine can save you a lot of error-checking code if you were going to
have to throw an error anyway.
@param[out]     out     The converted value.
@see tryConvertTo<T>(out), SimTK::convertStringTo<T>() **/
template <class T> inline void convertTo(T& out) const; // see below

/** A more convenient form of convertTo<T>() that returns the result as its 
function argument, although this may involve an extra copy operation. For very
large objects you may want to use the other form where the output is written
to an already-constructed object you provide. 
@return The converted value as an object of type T.
@see convertTo<T>(out), tryConvertTo<T>(out), SimTK::convertStringTo<T>() **/
template <class T> T convertTo() const 
{   T temp; convertTo<T>(temp); return temp; }

/*@}*/

/** @name In-place modifications
These are member functions which add to the existing std::string functionality.
These methods return a reference to "this" String, so may be chained like 
assignment statements. If you would like to use these on an std::string, use 
the String::updAs() method to recast the std::string to a String. Note that 
there is also an equivalent set of static methods which return a new String 
rather than changing the original. **/
/*@{*/
/** Upshift the given String in place, so that lowercase letters are replaced 
with their uppercase equivalents as defined by std::toupper(). **/
SimTK_SimTKCOMMON_EXPORT String& toUpper();
/** Downshift the given String in place, so that uppercase letters are replaced
with their lowercase equivalents as defined by std::tolower(). **/
SimTK_SimTKCOMMON_EXPORT String& toLower();
/** Trim this String in place, removing all the initial leading and trailing 
white space, as defined by std::isspace() which typically includes space, 
tab (\\t), newline (\\n), return (\\r),  and form feed (\\f). **/
SimTK_SimTKCOMMON_EXPORT String& trimWhiteSpace();
/** Substitute in place \a newChar for \a oldChar wherever \a oldChar appears
in this String. **/
SimTK_SimTKCOMMON_EXPORT String& replaceAllChar(char oldChar, char newChar);
/*@}*/


/** @name Utility methods
These static methods operate on SimTK::String or std::string objects and return
SimTK::String objects (which are also std::string objects). **/
/*@{*/
/** Cast an std::string to a SimTK::String without copying; subsequent changes
to the std::string will affect the SimTK::String too since it is just a 
reference to the original std::string. **/
static const String& getAs(const std::string& s) 
{   return reinterpret_cast<const String&>(s); }
/** Cast a non-const std::string to a non-const SimTK::String without copying;
changes made to the SimTK::String will affect the original std::string and
vice versa. **/
static String& updAs(std::string& s) 
{   return reinterpret_cast<String&>(s); }
/** Upshift the given std::string returning a new SimTK::String in which all
the letters have been made upper case with toupper(). **/
static String toUpper(const std::string& in)
{   return String(in).toUpper(); }
/** Downshift the given std::string returning a new SimTK::String in which all
the letters have be made lower case with tolower(). **/
static String toLower(const std::string& in)
{   return String(in).toLower(); }
/** Copy the input std::string to a new SimTK::String leaving off all the
initial leading and trailing white space, as defined by isspace() which
typically includes space, tab (\\t), newline (\\n), return (\\r), 
and form feed (\\f). **/
static SimTK_SimTKCOMMON_EXPORT String trimWhiteSpace(const std::string& in);
/** Copy the input std::string to a new SimTK::String while substituting
\a newChar for \a oldChar wherever \a oldChar appears in the input. **/
String& replaceAllChar(const std::string& in, char oldChar, char newChar)
{   return String(in).replaceAllChar(oldChar, newChar); }
/*@}*/

};	

// All std::stream activity should be dealt with inline so that we don't have
// to worry about binary compatibiliy issues that can arise when passing 
// streams through the API.

// Templatized String constructor.
template <class T> inline
String::String(const T& t) {
    std::ostringstream stream;
    stream << t;
    *this = stream.str();
}

// This namespace-level static method should not be necessary but gcc 4.1
// still has trouble with template specialization for template member
// functions. So rather than specializing the tryConvertTo() member, I'm 
// specializing this helper function instead.
template <class T> inline static
bool tryConvertStringTo(const String& value, T& out) {
	std::istringstream sstream(value);
	sstream >> out; if (sstream.fail()) return false;
    if (sstream.eof()) return true;
    // Successful conversion but didn't use all the characters. Maybe the
    // rest is just whitespace?
    std::ws(sstream);       // Skip trailing whitespace if any.
    return sstream.eof();   // We must have used up the whole string now.
}

// This specialization ensures that we get the whole String including
// leading and trailing white space. Of course this is not useful for 
// anything but may occur as a result of some higher-level templatized 
// method that doesn't know what type it is converting here.
template<> inline static
bool tryConvertStringTo(const String& value, String& out)
{   out = value; return true; }

// Same as above but for std::string output rather than String.
template<> inline static
bool tryConvertStringTo(const String& value, std::string& out)
{   out = value; return true; }

// This partial specialization ensures that you can't interpret
// a String as a pointer.
template<class T> inline static
bool tryConvertStringTo(const String& value, T*& out) {
    SimTK_ERRCHK1_ALWAYS(false, "SimTK::convertStringTo(value,T*)",
        "Can't interpret a string as a pointer (%s*).",
        NiceTypeName<T>::name());
    return false; 
}

template <class T> inline bool 
String::tryConvertTo(T& out) const 
{   return tryConvertStringTo(*this, out); }

template <class T> inline void 
String::convertTo(T& out) const {
    const int MaxStr = 50;
    const bool convertOK = tryConvertTo<T>(out);
    if (convertOK) return;

    // Make sure we don't try to output more than MaxStr characters of
    // the bad string in the error message.
    String shorter = this->substr(0, MaxStr);
    if (shorter.size() < this->size()) shorter += " ...";
    SimTK_ERRCHK2_ALWAYS(convertOK, "Xml::getAs()",
        "Couldn't interpret string '%s' as type T=%s.",
        shorter.c_str(), NiceTypeName<T>::name());
}

/** Alternate namespace-level static convertStringTo() functions are sometimes
more convenient; see the equivalent String::convertTo() method for information.
**/
template <class T> inline static
void convertStringTo(const String& in, T& out)
{   in.convertTo<T>(out); }

template <class T> inline static
T convertStringTo(const String& in)
{   return in.convertTo<T>(); }

}
#endif // SimTK_SimTKCOMMON_STRING_H_
