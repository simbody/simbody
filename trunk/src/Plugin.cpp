/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-10 Stanford University and the Authors.        *
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

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/String.h"
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


#ifdef _WIN32
    static const bool IsWindows          = true;
    static const char MyPathSeparator    = '\\';
    static const char OtherPathSeparator = '/';
#else
    static const bool IsWindows          = false;
    static const char MyPathSeparator    = '/';
    static const char OtherPathSeparator = '\\';
#endif

char Pathname::getPathSeparatorChar() {
    return MyPathSeparator;
}

string Pathname::getPathSeparator() {
    return String(getPathSeparatorChar());
}

static void makeNativeSlashesInPlace(string& inout) {
    String::updAs(inout).replaceAllChar(OtherPathSeparator, MyPathSeparator);
}

static void addFinalSeparatorInPlace(string& inout) {
    if (!inout.empty() && !Pathname::isPathSeparator(inout[inout.size()-1]))
        inout += MyPathSeparator;
}

static bool beginsWithPathSeparator(const string& in) {
    return !in.empty() && Pathname::isPathSeparator(in[0]);
}

// Remove the last segment of a path name and the separator. The 
// returned component does not include a separator.
static void removeLastPathComponentInPlace(string& inout, string& component) {
    component.clear();
    if (inout.empty()) return;

    const string::size_type lastSeparator = inout.find_last_of("/\\");
    if (lastSeparator == string::npos) {
        component = inout;
        inout.clear();
    } else {
        component = inout.substr(lastSeparator+1);
        inout.erase(lastSeparator);
    }
}

static void removeDriveInPlace(string& inout, string& drive) {
    drive.clear();
    if (IsWindows && inout.size() >= 2 && inout[1]==':') {
        drive = (char)tolower(inout[0]);
        inout.erase(0,2);
    }
}


// We assume a path name structure like this:
//   (1) Everything up to and including the final directory separator
//       character is the directory; the rest is the file name. On return
//       we fix the slashes in the directory name to suit the current
//       system ('\' for Windows, '/' otherwise).
//   (2) If the file name contains a ".", characters after the last
//       "." are the extension and the last "." is removed.
//   (5) What's left is the fileName.
// We accept both "/" and "\" as separator characters. We leave the
// case as it was supplied.
// Leading and trailing white space is removed; embedded white space
// remains.
// Leading "X:" for some drive letter X is recognized on Windows as
// a drive specification, otherwise the drive is the current drive.
// That is then removed for further processing.
// Absolute paths are designated like this:
//      Leading "/" means root relative (on the drive).
//      Leading "./" means current working directory relative (on the drive).
//      Leading "../" is interpreted as "./..".
//      Leading "@/" means relative to executable location.
// Above leading characters are removed and replaced with the
// full path name they represent, then the entire path is divided
// into components. If a component is "." or "" (empty) it is
// removed. If a component is ".." it and the previous component
// if any are removed. (".." as a first component will report as
// ill formed.)
// If there is something ill-formed about the file name we'll return
// false.
void Pathname::deconstructPathname( const string&   name,
                                    bool&           isAbsolutePath,
                                    string&         directory,
                                    string&         fileName,
                                    string&         extension)
{
    isAbsolutePath = false;
    directory.erase(); fileName.erase(); extension.erase();

    // Remove all the white space and make all the slashes be forward ones.
    // (For Windows they'll be changed to backslashes later.)
    String processed = String::trimWhiteSpace(name).replaceAllChar('\\', '/');
    if (processed.empty())
        return; // name consisted only of white space

    string drive;
    removeDriveInPlace(processed, drive);

    // Now the drive if any has been removed and we're looking at
    // the beginning of the path name. 

    // If the pathname in its entirety is just one of these, append 
    // a slash to avoid special cases below.
    if (processed=="." || processed==".." || processed=="@")
        processed += "/";

    // If the path begins with "../" we'll make it ./../ to simplify handling.
    if (processed.substr(0,3) == "../")
        processed.insert(0, "./");

    if (processed.substr(0,1) == "/") {
        isAbsolutePath = true;
        processed.erase(0,1);
        if (drive.empty()) drive = getCurrentDriveLetter();
    } else if (processed.substr(0,2) == "./") {
        isAbsolutePath = true;
        processed.replace(0,2,getCurrentWorkingDirectory(drive));
        removeDriveInPlace(processed, drive);
    } else if (processed.substr(0,2) == "@/") {
        isAbsolutePath = true;
        processed.replace(0,2,getThisExecutableDirectory());
        removeDriveInPlace(processed, drive);
    } else if (!drive.empty()) {
        // Looks like a relative path name. But if it had an initial
        // drive specification, e.g. X:something.txt, that is supposed
        // to be interpreted relative to the current working directory
        // on drive X, just as though it were X:./something.txt.
        isAbsolutePath = true;
        processed.insert(0, getCurrentWorkingDirectory(drive));
        removeDriveInPlace(processed, drive);
    }

    // We may have picked up a new batch of backslashes above.
    processed.replaceAllChar('\\', '/');

    // Now we have the full path name if this is absolute, otherwise
    // we're looking at a relative path name. In any case the last
    // component is the file name if it isn't empty, ".", or "..".

    // Process the ".." segments and eliminate meaningless ones
    // as we go through.
    Array_<string> segmentsInReverse;
    bool isFinalSegment = true; // first time around might be the fileName
    int numDotDotsSeen = 0;
    while (!processed.empty()) {
        string component;
        removeLastPathComponentInPlace(processed, component);
        if (component == "..")
            ++numDotDotsSeen;
        else if (!component.empty() && component != ".") {
            if (numDotDotsSeen)
                --numDotDotsSeen;   // skip component
            else if (isFinalSegment) fileName = component;
            else segmentsInReverse.push_back(component);
        }
        isFinalSegment = false;
    }

    // Now we can put together the canonicalized directory.
    if (isAbsolutePath) {
        if (!drive.empty())
            directory = drive + ":";
        directory += "/";
    }

    for (int i = (int)segmentsInReverse.size()-1; i >= 0; --i)
        directory += segmentsInReverse[i] + "/";

    // Fix the slashes.
    makeNativeSlashesInPlace(directory);

    // If there is a .extension, strip it off.
    string::size_type lastDot = fileName.rfind('.');
    if (lastDot != string::npos) {
        extension = fileName.substr(lastDot);
        fileName.erase(lastDot);
    }
}

string Pathname::getDefaultInstallDir() {
    string installDir;
    #ifdef _WIN32
        installDir = getEnvironmentVariable("ProgramFiles");
        if (!installDir.empty()) {
            makeNativeSlashesInPlace(installDir);
            addFinalSeparatorInPlace(installDir);
        } else
            installDir = "c:\\Program Files\\";
    #else
        installDir = "/usr/local/";
    #endif
    return installDir;
}

string Pathname::addDirectoryOffset(const string& base, const string& offset) {
    string cleanOffset = beginsWithPathSeparator(offset)
        ? offset.substr(1) : offset; // remove leading path separator
    addFinalSeparatorInPlace(cleanOffset); // add trailing / if non-empty
    string result = base;
    addFinalSeparatorInPlace(result);
    result += cleanOffset;
    return result;
}

string Pathname::getInstallDir(const std::string& envInstallDir,
                                 const std::string& offsetFromDefaultInstallDir) 
{   std::string installDir = getEnvironmentVariable(envInstallDir);
    if (!installDir.empty()) {
        makeNativeSlashesInPlace(installDir);
        addFinalSeparatorInPlace(installDir);
    } else
        installDir = addDirectoryOffset(getDefaultInstallDir(), offsetFromDefaultInstallDir);
    return installDir;
}

string Pathname::getThisExecutablePath() {
    char buf[2048];
    #if defined(_WIN32)
        const DWORD nBytes = GetModuleFileName((HMODULE)0, (LPTSTR)buf, sizeof(buf));
        buf[0] = (char)tolower(buf[0]); // drive name
    #elif defined(__APPLE__)
        uint32_t sz = (uint32_t)sizeof(buf);
        int status = _NSGetExecutablePath(buf, &sz);
        assert(status==0); // non-zero means path longer than buf, sz says what's needed
    #else // Linux
        // This isn't automatically null terminated.
        const size_t nBytes = readlink("/proc/self/exe", buf, sizeof(buf));
        buf[nBytes] = '\0';
    #endif
    return string(buf);
}

string Pathname::getThisExecutableDirectory() {
    string path = getThisExecutablePath();
    string component;
    removeLastPathComponentInPlace(path, component);
    path += MyPathSeparator;
    return path;
}

string Pathname::getCurrentDriveLetter() {
    #ifdef _WIN32
        const int which = _getdrive();
        return string() + (char)('a' + which-1);
    #else
        return string();
    #endif
}

string Pathname::getCurrentDrive() {
    #ifdef _WIN32
        return getCurrentDriveLetter() + ":";
    #else
        return string();
    #endif
}


string Pathname::getCurrentWorkingDirectory(const string& drive) {
    char buf[1024];

    #ifdef _WIN32
        const int which = drive.empty() ? 0 : (tolower(drive[0]) - 'a') + 1;
        assert(which >= 0);
        if (which != 0) {
            // Make sure this drive exists.
            const ULONG mask = _getdrives();
            if (!(mask & (1<<(which-1))))
                return getRootDirectory(drive);
        }
        _getdcwd(which, buf, sizeof(buf));
        buf[0] = (char)tolower(buf[0]); // drive letter
    #else
        getcwd(buf, sizeof(buf));
    #endif

    string cwd(buf);
    if (cwd.size() && cwd[cwd.size()-1] != MyPathSeparator)
        cwd += MyPathSeparator;
    return cwd;
}

string Pathname::getRootDirectory(const string& drive) {
    #ifdef _WIN32
        if (drive.empty()) 
            return getCurrentDrive() + getPathSeparator();
        return String(drive[0]).toLower() + ":" + getPathSeparator();
    #else
        return getPathSeparator();
    #endif
}

bool Pathname::environmentVariableExists(const string& name) {
    return getenv(name.c_str()) != 0;
}

string Pathname::getEnvironmentVariable(const string& name) {
    char* value;
    value = getenv(name.c_str());
    return value ? string(value) : string();
}


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

