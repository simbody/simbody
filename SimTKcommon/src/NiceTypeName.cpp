/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-14 Stanford University and the Authors.        *
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

#if defined(__GNUG__)
// https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.3/a01696.html
    #include <cxxabi.h>
    #include <string>
#endif

namespace SimTK {

std::string demangle(const char* name) {
    #if defined(__GNUG__)
        int status;
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

}
