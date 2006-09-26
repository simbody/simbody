#ifndef SimTK_STRING_H_
#define SimTK_STRING_H_

// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/List.h"

#include <cstdio>
#include <string>
#include <limits>
#include <complex>

namespace SimTK {
	
/**
 * Temporary implementation -- must hide implementation
 */
class String : public std::string {
public:
	String() { }
    String(const String& s, size_t start, size_t len) : std::string(s,start,len) { }
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
	
	String(const char* s) : std::string(s) { }
	String(const String& s) : std::string(s) { }
	String(const std::string& s) : std::string(s) { }
    
    operator const char*() const { return c_str(); }
};	

SimTK_LIST_SPECIALIZE(String);

}
#endif //SimTK_STRING_H_
