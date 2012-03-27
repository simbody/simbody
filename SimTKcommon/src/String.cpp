/* -------------------------------------------------------------------------- *
 *                      SimTK Simbody: SimTKcommon                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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

/** @file
 * This file contains the non-inline implementations of the SimTK::String
 * class.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/NTraits.h"

#include <string>
#include <cctype>

using SimTK::String;

String::String(float r, const char* fmt) {
    if (!isFinite(r)) {
        if (isNaN(r)) {(*this)="NaN"; return;}
        if (isInf(r)) {(*this)=(r<0?"-Infinity":"Infinity"); return;}
        SimTK_ERRCHK2_ALWAYS(false, "SimTK::String(float)",
            "Unrecognized non-finite value %g (0x%x).", 
            (double)r, *reinterpret_cast<const unsigned*>(&r));
        return;
    }
    char buf[64]; sprintf(buf,fmt,r); (*this)=buf; 
}

String::String(double r, const char* fmt) {
    if (!isFinite(r)) {
        if (isNaN(r)) {(*this)="NaN"; return;}
        if (isInf(r)) {(*this)=(r<0?"-Infinity":"Infinity"); return;}
        SimTK_ERRCHK2_ALWAYS(false, "SimTK::String(double)",
            "Unrecognized non-finite value %g (0x%llx).", 
            r, *reinterpret_cast<const unsigned long long*>(&r));
        return;
    }
    char buf[64]; sprintf(buf,fmt,r); (*this)=buf; 
}

String::String(long double r, const char* fmt) {
    if (!isFinite(r)) {
        if (isNaN(r)) {(*this)="NaN"; return;}
        if (isInf(r)) {(*this)=(r<0?"-Infinity":"Infinity"); return;}
        SimTK_ERRCHK1_ALWAYS(false, "SimTK::String(long double)",
            "Unrecognized non-finite value %lg.", r);
        return;
    }
    char buf[128]; sprintf(buf,fmt,r); (*this)=buf; 
}

static String cleanUp(const String& in) {
    return String(in).trimWhiteSpace().toLower();
}

bool String::tryConvertToBool(bool& out) const {
    const String adjusted = cleanUp(*this);
    if (adjusted=="true")  {out=true;  return true;}
    if (adjusted=="false") {out=false; return true;}
	std::istringstream sstream(adjusted);
	sstream >> out;
    return !sstream.fail();
}

bool String::tryConvertToFloat(float& out) const {
    const String adjusted = cleanUp(*this);
    if (adjusted=="nan")  {out=NTraits<float>::getNaN();  return true;}
    if (   adjusted=="inf" || adjusted=="infinity"
        || adjusted=="+inf" || adjusted=="+infinity") 
    {   out = NTraits<float>::getInfinity(); return true;}
    if (adjusted=="-inf" || adjusted=="-infinity") 
    {   out = -NTraits<float>::getInfinity(); return true;}
	std::istringstream sstream(adjusted);
	sstream >> out;
    return !sstream.fail();
}

bool String::tryConvertToDouble(double& out) const {
    const String adjusted = cleanUp(*this);
    if (adjusted=="nan")  {out=NTraits<double>::getNaN();  return true;}
    if (   adjusted=="inf" || adjusted=="infinity"
        || adjusted=="+inf" || adjusted=="+infinity") 
    {   out = NTraits<double>::getInfinity(); return true;}
    if (adjusted=="-inf" || adjusted=="-infinity") 
    {   out = -NTraits<double>::getInfinity(); return true;}
	std::istringstream sstream(adjusted);
	sstream >> out;
    return !sstream.fail();
}

bool String::tryConvertToLongDouble(long double& out) const {
    const String adjusted = cleanUp(*this);
    if (adjusted=="nan")  {out=NTraits<long double>::getNaN();  return true;}
    if (   adjusted=="inf" || adjusted=="infinity"
        || adjusted=="+inf" || adjusted=="+infinity") 
    {   out = NTraits<long double>::getInfinity(); return true;}
    if (adjusted=="-inf" || adjusted=="-infinity") 
    {   out = -NTraits<long double>::getInfinity(); return true;}
	std::istringstream sstream(adjusted);
	sstream >> out;
    return !sstream.fail();
}



String& String::toUpper() {
    for (int i=0; i < size(); ++i)
        (*this)[i] = (char)std::toupper((*this)[i]);
    return *this;
}

String& String::toLower() {
    for (int i=0; i < size(); ++i)
        (*this)[i] = (char)std::tolower((*this)[i]);
    return *this;
}

String& String::replaceAllChar(char oldChar, char newChar) {
    for (int i=0; i < size(); ++i)
        if ((*this)[i] == oldChar)
            (*this)[i] = newChar;
    return *this;
}

String& String::trimWhiteSpace() {
    *this = trimWhiteSpace(*this);
    return *this;
}

String String::trimWhiteSpace(const std::string& in) {
    const int inz = (int)in.size();

    // Find first non-white character position of "in".
    int firstNonWhite = 0;
    for ( ; firstNonWhite < inz; ++firstNonWhite)
        if (!std::isspace(in[firstNonWhite])) break;

    if (firstNonWhite == inz)
        return String();    // "in" was all white space

    // Find last non-white character position of "in".
    int lastNonWhite = inz-1;
    for ( ; lastNonWhite >= 0; --lastNonWhite)
        if (!std::isspace(in[lastNonWhite])) break;

    return String(in, firstNonWhite, (lastNonWhite+1) - firstNonWhite);
}


