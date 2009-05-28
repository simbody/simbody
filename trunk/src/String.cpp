/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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
#include "SimTKcommon/internal/String.h"

#include <string>
#include <cctype>

using SimTK::String;

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