#ifndef SimTK_SIMBODY_EXAMPLE_HELPER_H_
#define SimTK_SIMBODY_EXAMPLE_HELPER_H_

/* -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

/* This class contains utility methods that can be helpful for Simbody
example programs. In particular it knows the search path to be used to find
auxiliary and shared files needed by particular examples.

We depend on two macros being defined at compile time:
    SIMBODY_EXAMPLE_NAME        
        the example name is also the subdirectory name for auxilary files 
        for that example
    SIMBODY_EXAMPLES_INSTALL_SRC
        this is the directory into which Simbody example source was installed

The Simbody example CMakeList.txt will define these macros when building
Simbody.
*/
#include "Simbody.h"

#include <string>
#include <vector>

class SimbodyExampleHelper {
public:
    // Return the name of the currently-executing example.
    static std::string getExampleName() 
    {   return SIMBODY_EXAMPLE_NAME; }

    // Return the name of the directory to which example source files were
    // installed. This is where the "shared/" subdirectory can be found as
    // well as the subdirectories specific to particular examples.
    static std::string getExamplesSourceInstallDir()
    {   return SIMBODY_EXAMPLES_INSTALL_SRC; }

    // To find the appropriate directory containing auxiliary files for this
    // example, call this method with the name of one of the files that the
    // example depends on. Then use that same directory for all the files.
    // This will first assume this is a freshly-built example and try to find
    // a local copy of the auxiliary directory; if that doesn't work it will 
    // look in the installed location.
    // The returned directory name will end in a "/"; an exception will be
    // thrown with a helpful error message if we can't find what you asked for.
    static std::string 
        findAuxiliaryDirectoryContaining(const std::string& auxFile)
    {
        std::vector<std::string> searchPath;
        searchPath.push_back("./");                         // Windows build
        searchPath.push_back(std::string("./examples/")     // Linux build
                             + SIMBODY_EXAMPLE_NAME + "/");
        searchPath.push_back(                               // Installed
            std::string(SIMBODY_EXAMPLES_INSTALL_SRC) 
                        + SIMBODY_EXAMPLE_NAME + "/");

        for (unsigned i=0; i < searchPath.size(); ++i)
            if (SimTK::Pathname::fileExists(searchPath[i] + auxFile))
                return searchPath[i];

        // Create a chronicle of our failures above.
        std::string failedPath;
        for (unsigned i=0; i < searchPath.size(); ++i)
            failedPath += ("    " + searchPath[i] + "\n");

        SimTK_ERRCHK3_ALWAYS(!"can't find auxiliary directory", 
            "SimbodyExampleHelper::findAuxiliaryDirectoryContaining()",
            "\n  Example %s requires auxiliary files to run, but was"
            "\n  unable to locate an auxiliary directory containing file '%s'."
            "\n  Searched directories:\n%s",
            SIMBODY_EXAMPLE_NAME, auxFile.c_str(), failedPath.c_str());

        return "";
    }

    // To find the directory containing shared files for all the examples,
    // call this method with the name of one of the shared files that the
    // example depends on. Then use that same directory for all other shared
    // files. This will first assume this is a freshly-built example and try to
    // find a local copy of the shared directory; if that doesn't work it will 
    // look in the installed location.
    // The returned directory name will end in a "/"; an exception will be
    // thrown with a helpful error message if we can't find what you asked for.
    static std::string 
        findSharedDirectoryContaining(const std::string& auxFile)
    {
        std::vector<std::string> searchPath;
        searchPath.push_back("../shared/");         // Windows build
        searchPath.push_back("./examples/shared/"); // Linux build
        searchPath.push_back(                       // Installed
            std::string(SIMBODY_EXAMPLES_INSTALL_SRC) + "shared/");

        for (unsigned i=0; i < searchPath.size(); ++i)
            if (SimTK::Pathname::fileExists(searchPath[i] + auxFile))
                return searchPath[i];

        // Create a chronicle of our failures above.
        std::string failedPath;
        for (unsigned i=0; i < searchPath.size(); ++i)
            failedPath += ("    " + searchPath[i] + "\n");

        SimTK_ERRCHK3_ALWAYS(!"can't find shared directory", 
            "SimbodyExampleHelper::findSharedDirectoryContaining()",
            "\n  Example %s requires shared files to run, but was"
            "\n  unable to locate a shared directory containing file '%s'."
            "\n  Searched directories:\n%s",
            SIMBODY_EXAMPLE_NAME, auxFile.c_str(), failedPath.c_str());

        return "";
    }
};


#endif // SimTK_SIMBODY_EXAMPLE_HELPER_H_
