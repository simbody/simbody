/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Chris Dembia, Michael Sherman                                     *
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

#include "SimTKcommon/internal/common.h"
#include <string>
#include <regex>
#include <array>
#include <algorithm>
#include <iostream>

#if defined(__GNUG__)
// https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.3/a01696.html
    #include <cxxabi.h>
    #include <cstdlib>
#endif

namespace SimTK {

std::string demangle(const char* name) {
    #if defined(__GNUG__)
        int status=-100; // just in case it doesn't get set
        char* ret = abi::__cxa_demangle(name, NULL, NULL, &status);
        const char* const demangled_name = (status == 0) ? ret : name;
        std::string demangled_string(demangled_name);
        if (ret) std::free(ret);
        return demangled_string;
    #else
        // On other platforms, we hope the typeid name is not mangled.
        return name;
    #endif
}

// Given a demangled string, attempt to canonicalize it for platform
// indpendence. We'll remove Microsoft's "class ", "struct ", etc.
// designations, and get rid of all unnecessary spaces.
std::string canonicalizeTypeName(std::string&& demangled) {
    using SPair = std::pair<std::regex,std::string>;
    // These are applied in this order.
    static const std::array<SPair,7> subs{
        // Remove unwanted keywords and following space.
        SPair(std::regex("\\b(class|struct|enum|union) "),      ""),
        // Tidy up anonymous namespace.
        SPair(std::regex("[`(]anonymous namespace[')]"),        "(anonymous)"),
        // Standardize "unsigned int" -> "unsigned".
        SPair(std::regex("\\bunsigned int\\b"),                 "unsigned"),
        // Temporarily replace spaces we want to keep with "!". (\w is 
        // alphanumeric or underscore.)
        SPair(std::regex("(\\w) (\\w)"),    "$1!$2"),
        SPair(std::regex(" "),              ""), // Delete unwanted spaces.
        // OSX clang throws in extra namespaces like "__1". Delete them.
        SPair(std::regex("\\b__[0-9]+::"),  ""),
        SPair(std::regex("!"),              " ") // Restore wanted spaces.
    };
    std::string canonical(std::move(demangled));
    for (const auto& sp : subs) {
        canonical = std::regex_replace(canonical, sp.first, sp.second);
    }
    return canonical;
}

std::string encodeTypeNameForXML(std::string&& nicestr) {
    std::string xmlstr(std::move(nicestr));
    std::replace(xmlstr.begin(), xmlstr.end(), '<', '{');
    std::replace(xmlstr.begin(), xmlstr.end(), '>', '}');
    return xmlstr;
}

std::string decodeXMLTypeName(std::string&& xmlstr) {
    std::string nicestr(std::move(xmlstr));
    std::replace(nicestr.begin(), nicestr.end(), '{', '<');
    std::replace(nicestr.begin(), nicestr.end(), '}', '>');
    return nicestr;
}


}
