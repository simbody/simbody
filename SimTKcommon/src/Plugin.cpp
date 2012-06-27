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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
#include "SimTKcommon/internal/Pathname.h"
#include "SimTKcommon/internal/Plugin.h"

#include <string>
using std::string;

#include <cctype>
using std::tolower;

#include <iostream>

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include <windows.h>
    #include <direct.h>
    #pragma warning(disable:4996) // getenv() is apparently unsafe
#else
    #ifdef __APPLE__
        #include <mach-o/dyld.h>
    #endif
    #include <dlfcn.h>
    #include <dirent.h>
    #include <unistd.h>
#endif


using namespace SimTK;

Plugin::Plugin(const string& name) : m_handle(0) {
    m_defaultName = name; // might be empty
}

Plugin::~Plugin() {
    unload();
}

bool Plugin::load(const string& name) {
    if (m_handle) {
        m_lastMessage = "Plugin::load(): already loaded " + m_loadedPathname;
        return false;
    }

    const string nameToUse = name.empty() ? m_defaultName : name;

    if (nameToUse.empty()) {
        m_lastMessage = "Plugin::load(): no pathname was supplied and there was no default.";
        return false;
    }

    bool isAbsolutePath;
    string directory, libPrefix, baseName, debugSuffix, extension;
    if (!deconstructLibraryName(nameToUse, isAbsolutePath,
        directory, libPrefix, baseName, debugSuffix, extension)) 
    {
        m_lastMessage = "Plugin::load(): illegal library name '" + nameToUse + "'.";
        return false;
    }

    const std::string base = directory + getDynamicLibPrefix() + baseName;
    if (isAbsolutePath) {
        m_handle = loadDebugOrReleasePlugin(base, extension, m_loadedPathname, m_lastMessage);
    } else {
        for (unsigned i=0; !m_handle && i < m_searchPath.size(); ++i)
            m_handle = loadDebugOrReleasePlugin(
                            m_searchPath[i] + base, extension, m_loadedPathname, m_lastMessage);
    }

    return m_handle != 0;
}

void Plugin::unload() {
    if (m_handle) {
        std::cout << "UNLOADING '" << m_loadedPathname << "'.\n";
        unloadPlugin(m_handle);
        m_handle = 0;
        m_loadedPathname.clear();
    }
}

// This is like deconstructPathname() but adds the assumption that
// the pathname represents a library file. We assume further structure
// for the "fileName" part: it may have a "lib" prefix and "_d"
// debug suffix which are stripped off to reveal the baseName.
// Also, this method returns false if there is no baseName
// left either because the original path was a directory name
// or because there was nothing left after stripping the prefix
// or suffix.
bool Plugin::deconstructLibraryName(const string& name,
                                    bool&   isAbsolutePath,
                                    string& directory,
                                    string& libPrefix,
                                    string& baseName,
                                    string& debugSuffix,
                                    string& extension)
{
    libPrefix.erase(); baseName.erase(); debugSuffix.erase();

    Pathname::deconstructPathname(name, isAbsolutePath, directory, baseName, extension);

    // If there is a "lib" prefix, strip it off.
    if (baseName.substr(0,3) == "lib") {
        libPrefix = "lib";
        baseName.erase(0,3);
    }

    // If there is a trailing "_d" or "_D", strip it off. Preserve
    // the case in the debugSuffix.
    if (baseName.size() >= 2) {
        const char last2 = baseName[baseName.size()-2];
        const char last  = baseName[baseName.size()-1];
        if (last2 == '_' && (last=='d' || last=='D')) {
            debugSuffix = string("_") + last;
            baseName.erase(baseName.size()-2);
        }
    }

    // That must have left us with something in the baseName or
    // this is not a legal library name.

    return !baseName.empty();
}


#ifdef _WIN32
static string getWindowsSystemMessage() {
    LPVOID lpMsgBuf;
    FormatMessage( 
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM | 
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        GetLastError(),
        0, // Default language
        (LPTSTR) &lpMsgBuf,
        0,
        NULL 
    );
    const string retStr((const char*)lpMsgBuf);

    LocalFree( lpMsgBuf );
    return retStr;
}
#endif


void* Plugin::loadPluginByFileName(const string& name, string& errorMessage) {
    errorMessage.clear();

    std::cout << "LOAD ATTEMPT: '" << name << "' ... ";

    #ifdef _WIN32
        // Tell Windows not to bother the user with ugly error boxes.
        const UINT oldErrorMode = SetErrorMode(SEM_FAILCRITICALERRORS);
        HMODULE handle = LoadLibrary(name.c_str());
        SetErrorMode(oldErrorMode); // Restore previous error mode.
        if (!handle) errorMessage = getWindowsSystemMessage();
    #else
        dlerror(); // clear pre-existing message if any
        void* handle = dlopen(name.c_str(), RTLD_LAZY);
        if (!handle) {
            const char* msg = dlerror();
            if (msg) errorMessage = string(msg);
        }
    #endif

    std::cout << (handle ? "SUCCEEDED." : "FAILED!") << std::endl;

    return handle;
}

// This attempts first to load the debug version of the library if we're
// in Debug mode, but falls back to the Release if no Debug library is
// available. If we're running in Release mode then no attempt is made
// to load the Debug library.
void* Plugin::loadDebugOrReleasePlugin(const string& base, const string& extension,
                                       string& loadedFileName, string& errorMessage)
{
    void* handle = 0;
    if (!getDynamicLibDebugSuffix().empty()) {
        // Attempt to load the Debug library if it exists.
        loadedFileName = base + getDynamicLibDebugSuffix() + getDynamicLibExtension();
        handle = loadPluginByFileName(loadedFileName, errorMessage);
    }
    if (!handle) {
        // Attempt to load the Release library.
        loadedFileName = base + getDynamicLibExtension();
        handle = loadPluginByFileName(loadedFileName, errorMessage);
    }
    if (!handle)
        loadedFileName.clear();

    return handle;
}

void Plugin::unloadPlugin(void* handle) {
    // We're ignoring the return values here.
    #ifdef _WIN32
        (void)FreeLibrary((HMODULE)handle);
    #else
        (void)dlclose(handle);
    #endif
}

void* Plugin::getSymbolAddress(void* handle, const string& name, string& errorMessage) {
    void* address;
    errorMessage.clear();

    #ifdef _WIN32
        address = GetProcAddress((HMODULE)handle, name.c_str());
        if (!address)
            errorMessage = getWindowsSystemMessage();
    #else
        dlerror(); // clear pre-existing message if any
        address = dlsym(handle, name.c_str());
        if (!address) {
            const char* msg = dlerror();
            if (msg) errorMessage = string(msg);
        }
    #endif

    return address;
}

std::string Plugin::getDynamicLibPrefix() {
    #ifdef _WIN32
        return "";
    #else
        return "lib";
    #endif
}

std::string Plugin::getDynamicLibExtension() {
    #if defined(_WIN32)
        return ".dll";
    #elif defined(__APPLE__)
        return ".dylib";
    #else   // linux
        return ".so";
    #endif
}

std::string Plugin::getDynamicLibDebugSuffix() {
    #ifdef NDEBUG
        return "";      // running in Release mode
    #else
        return "_d";    // running in Debug mode
    #endif
}

