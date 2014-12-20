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

#include <string>
using std::string;

#include <cctype>
using std::tolower;

#include <iostream>
#include <fstream>

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
    for (unsigned i=0; i < inout.size(); ++i)
        if (inout[i] == OtherPathSeparator)
            inout[i] = MyPathSeparator;
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

void Pathname::deconstructPathRelativeToSWD(const std::string& swd,
	const std::string& path,
	std::string& directory,
	std::string& fileName,
	std::string& extension)
{
	directory.erase(); fileName.erase(); extension.erase();

	// Remove all the white space and make all the slashes be forward ones.
	// (For Windows they'll be changed to backslashes later.)
	String pathcleaned = String::trimWhiteSpace(path).replaceAllChar('\\', '/');
	if (pathcleaned.empty())
		return; // path consisted only of white space

	// Remove all the white space and make all the slashes be forward ones.
	// (For Windows they'll be changed to backslashes later.)
	// Then add a path seprator if needed.
	String swdcleaned = String::trimWhiteSpace(swd).replaceAllChar('\\', '/');
	if (!swdcleaned.empty() && swdcleaned[swdcleaned.size() - 1] != getPathSeparatorChar())
		swdcleaned += getPathSeparatorChar();

	string pathdrive, swddrive, finaldrive;
	removeDriveInPlace(pathcleaned, pathdrive);
	removeDriveInPlace(swdcleaned, swddrive);


	// If the pathname in its entirety is just one of these, append 
	// a slash to avoid special cases below.
	if (pathcleaned == "." || pathcleaned == "..")
		pathcleaned += "/";

	// If the path begins with "../" we'll make it ./../ to simplify handling.
	if (pathcleaned.substr(0, 3) == "../")
		pathcleaned.insert(0, "./");
	
	// If swd is empty, then process path as usual.
	if (swd.empty()) {
		
		if (pathcleaned.substr(0, 1) == "/") {
			pathcleaned.erase(0, 1);
			if (pathdrive.empty()) pathdrive = getCurrentDriveLetter();
		}
		else if (pathcleaned.substr(0, 2) == "./") {
			pathcleaned.replace(0, 2, getCurrentWorkingDirectory(pathdrive));
			removeDriveInPlace(pathcleaned, pathdrive);
		}
		// Otherwise, the path was "X:*" or "*", so find the working directory.
		else {
			pathcleaned.insert(0, getCurrentWorkingDirectory(pathdrive));
			removeDriveInPlace(pathcleaned, pathdrive);
		}
		finaldrive = pathdrive;
	}

	// If swd was not empty, find the correct absolute path to the intended directory.
	else if (!swd.empty()) {
		// Path was of the form "X:/*" or "/*"
		if (pathcleaned.substr(0, 1) == "/") {
			pathcleaned.erase(0, 1);

			// Path is of the form "X:/*"
			if (!pathdrive.empty())
				finaldrive = pathdrive;

			// If path does not provide a drive, but swd provides a drive, use it.
			// If neither provides a drive, use the current drive.
			else if (pathdrive.empty()) {
				if (!swddrive.empty())
					finaldrive = swddrive;
				else if (swddrive.empty())
					finaldrive = getCurrentDriveLetter();
			}
		}
		// Path was of the form "X:./*" or "./*". 
		else if (pathcleaned.substr(0, 2) == "./") {
			// If swd provided a drive, use it
			if (!swddrive.empty()) {
				finaldrive = swddrive;
				pathcleaned.replace(0, 2, swdcleaned);
			}

			// Otherwise, check if path provided a drive and use it. If neither
			// has a path, use the absolute path of the swd.
			else if (swddrive.empty()) {
				string swdfullpath = getAbsoluteDirectoryPathname(swdcleaned);
				removeDriveInPlace(swdfullpath, swddrive);
				if (!pathdrive.empty())
					finaldrive = pathdrive;
				else if (pathdrive.empty())
					finaldrive = swddrive;
				pathcleaned.replace(0, 2, swdfullpath);
			}
		}

		// Path was of the form "X:*" or "*"
		else {
			// If swd provided a drive, use it
			if (!swddrive.empty()) {
				finaldrive = swddrive;
				pathcleaned.insert(0, swdcleaned);
			}

			// Otherwise, check if path provided a drive and use it. If neither
			// has a path, use the absolute path of the swd.
			else if (swddrive.empty()) {
				string swdfullpath = getAbsoluteDirectoryPathname(swdcleaned);
				
				removeDriveInPlace(swdfullpath, swddrive);
				
				if (!pathdrive.empty())
					finaldrive = pathdrive;
				else if (pathdrive.empty())
					finaldrive = swddrive;
				pathcleaned.insert(0, swdfullpath);
			}
		}
	}

	// We may have picked up a new batch of backslashes above.
	pathcleaned.replaceAllChar('\\', '/');

	// Process the ".." segments and eliminate meaningless ones
	// as we go through.
	Array_<string> segmentsInReverse;
	bool isFinalSegment = true; // first time around might be the fileName
	int numDotDotsSeen = 0;
	while (!pathcleaned.empty()) {
		string component;
		removeLastPathComponentInPlace(pathcleaned, component);
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
	if (!finaldrive.empty())
		directory = finaldrive + ":";
	directory += "/";

	for (int i = (int)segmentsInReverse.size() - 1; i >= 0; --i)
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

bool Pathname::fileExists(const std::string& fileName) {
    std::ifstream f(fileName.c_str());
    return f.good();
    // File is closed by destruction of f.
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
        const int status = _NSGetExecutablePath(buf, &sz);
        SimTK_ERRCHK_ALWAYS(status==0,
                "Pathname::getThisExecutablePath()",
                "2048-byte buffer is not big enough to store executable path.");
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
    char buf[2048];

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
        const char* bufp = getcwd(buf, sizeof(buf));
        SimTK_ERRCHK_ALWAYS(bufp != 0,
            "Pathname::getCurrentWorkingDirectory()",
            "2048-byte buffer not big enough for current working directory.");
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


