/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-12 Stanford University and the Authors.        *
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
 * This file contains the non-inline implementations of the SimTK::String
 * class.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/NTraits.h"
#include "SimTKcommon/internal/Serialize.h"

#include <string>
#include <cctype>

using SimTK::String;

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about sprintf, etc.
#endif

String::String(float r, const char* fmt) {
    if (!isFinite(r)) {
        if (isNaN(r)) {(*this)="NaN"; return;}
        if (isInf(r)) {(*this)=(r<0?"-Inf":"Inf"); return;}
        SimTK_ERRCHK1_ALWAYS(false, "SimTK::String(float)",
            "Unrecognized non-finite value %g.", (double)r);
        return;
    }
    char buf[64]; sprintf(buf,fmt,r); (*this)=buf; 
}

String::String(double r, const char* fmt) {
    if (!isFinite(r)) {
        if (isNaN(r)) {(*this)="NaN"; return;}
        if (isInf(r)) {(*this)=(r<0?"-Inf":"Inf"); return;}
        SimTK_ERRCHK1_ALWAYS(false, "SimTK::String(double)",
            "Unrecognized non-finite value %g.", r);
        return;
    }
    char buf[64]; sprintf(buf,fmt,r); (*this)=buf; 
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
        if (!std::isspace((unsigned char)in[firstNonWhite])) break;

    if (firstNonWhite == inz)
        return String();    // "in" was all white space

    // Find last non-white character position of "in".
    int lastNonWhite = inz-1;
    for ( ; lastNonWhite >= 0; --lastNonWhite)
        if (!std::isspace((unsigned char)in[lastNonWhite])) break;

    return String(in, firstNonWhite, (lastNonWhite+1) - firstNonWhite);
}


