#ifndef SimTK_SimTKCOMMON_PATHNAME_H_
#define SimTK_SimTKCOMMON_PATHNAME_H_

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
Declaration of the SimTK::Pathname class providing platform-independent
manipulation of file pathnames. **/

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Array.h"
#include <string>
#include <stdexcept>

namespace SimTK {

/**
 * This class encapsulates the handling of file and directory pathnames
 * in a platform-independent manner. We consider a pathname to consist
 * of three components:
 * <pre> [directory] [filename [extension]] </pre>
 * where the directory may be an absolute location or relative to a
 * current working directory. 
 *
 * Several special directory names are supported here:
 *  - root (/)
 *  - current working directory (.)
 *  - current executable directory (@)
 *  - platform default installation directory
 *  - parent directory (..)
 *
 * On Windows root and current working directory are drive-specific, referring
 * to the current drive if none is specified. The current executable 
 * directory is the absolute directory name containing the executable
 * program which is currently running.
 *
 * A pathname has "segments" which are separated by either forward 
 * slashes or backslashes. We are relaxed about the slashes and will
 * accept either one and pathnames which use both. However, each operating
 * system platform has a preferred separator character, backslash on
 * Windows and forward slash everywhere else and we will clean up
 * returned pathnames to use exclusively the preferred separator for
 * the current platform.
 *
 * Pathnames that end in an empty segment, or a segment consisting of
 * just "." or ".." are directory path names, meaning that the "filename"
 * and "extension" components are empty. Other pathnames may be directories
 * or filenames depending on context. Whenever we generate a pathname that
 * we know to be a directory, it will end in a final slash.
 *
 * There is also the concept of a "drive" which on Windows is a drive 
 * letter followed by a colon (e.g. "c:") but is always the empty string 
 * on non-Windows platforms. The drive, if present is considered part of 
 * the directory and does not affect whether the directory is considered 
 * relative or absolute. Drive designators are recognized only on Windows;
 * they are just considered ordinary pathname characters on other platforms.
 *
 * This class is useful for generating "canonicalized" pathnames from
 * names that have been pieced together from environment variables and
 * user entry. Canonicalized names are always absolute pathnames; they
 * contain no empty, ".", or ".." segments. On Windows a canonicalized
 * name is always prefixed by an explicit disk designator followed by
 * a backslash; on other platforms the canonicalized name will always
 * begin with a forward slash.
 *
 */
class SimTK_SimTKCOMMON_EXPORT Pathname {
public:


    /// Dismantle a supplied pathname into its component
    /// parts. This can take pathnames like <pre>   
    ///     /usr/local/libMyDll_d.so
    ///     e:\\Program Files\\Something\\myLibrary_d.dll
    /// </pre> and chop them into <pre>
    /// directory                       fileName       extension
    /// ------------------------------- -------------- ---------
    /// /usr/local/                     libMyDll_d     .so 
    /// e:\\Program Files\\Something\\     myLibrary_d    .dll
    /// </pre>
    /// as well as tell you whether the given pathname is absolute or relative 
    /// (and thus subject to search rules). At the beginning of the pathname
    /// (or right after the drive specification on Windows) we recognize
    /// three special symbols:
    /// - "/" means root; i.e., this is an absolute path name starting from
    ///   the root directory (this drive's root for Windows).
    /// - "." starts an absolute path name which is relative to the
    ///   current working directory (or drive's cwd on Windows).
    /// - "@" starts an absolute path name which is relative to the
    ///   directory in which the currently running executable is located.
    ///
    /// Anywhere in the pathname, the name ".." means "go up one level from
    /// the prior directory". ".." at the start is interpreted as "./..".
    /// A '.' appearing anywhere in the path name except the beginning is
    /// ignored. An '@' appearing anywhere in the pathname other than the
    /// beginning is treated as an ordinary file character.
    ///
    /// The pathname components are returned
    /// as separate strings with separators included such that concatenating 
    /// all the strings reproduces the pathname in a canonicalized form.
    /// The "drive" letter prefix is recognized only when running on Windows; 
    /// otherwise a prefix like "C:" is treated as ordinary file name 
    /// characters. Note that we include the drive letter as part of the 
    /// absolute directory.
    /// White space is removed, and path separator characters
    /// in the directory are changed to the appropriate slash
    /// for the currently running platform (i.e. backslash for
    /// Windows and forward slash everywhere else).
    static void deconstructPathname(    const std::string& name,
                                        bool&        dontApplySearchPath,
                                        std::string& directory,
                                        std::string& fileName,
                                        std::string& extension);

    /// An extension of deconstructPathname(). Given a specified working directory 
    /// (swd) and path, this function evaluates the absolute path of a given path
    /// relative to the swd and returns the directory, fileName, and extension of 
    /// the canonicalized path with respect to a swd, if needed. This means, for 
    /// the path, that instead of evaluating "." as the current working directory 
    /// (cwd), the swd is used. Unlike deconstructPathname(), this function will
    /// always return an absolute path, and no bool dontApplySearchPath is returned.
    /// Rules:
    /// - If path is empty, directory, fileName and extension will be returned empty.
    ///   This case probably should just use getCurrentworkingDirectory().
    /// - If the swd is empty (after removing whitespace), deconstructPathname()
    ///   is called, and cwd is prepended if needed to make it an absolute path.
    /// - Otherwise, we evaluate path relative to the swd. These steps are as follows:
    /// 1) If path is a root-relative path name (and on Windows this includes a drive) 
    ///    (e.g. /usr/file.ext or c:/documents/file.ext), then swd is ignored, and the
    ///    absolute path is returned.
    /// 2) Preprocess the swd. This means that if the swd is of any form that denotes
    ///    an absolute path (i.e. "C:/file.ext", "C:file.ext", "./file.ext", "/file.ext")
    ///    we change the swd to reflect the absolute path (e.g. "./file.ext" may change
    ///    to "/cwd/file.ext" or "C:/cwdOnC/file.ext").
    /// 3) Otherwise, if a path is given relative to a directory that is not the root 
    ///    (e.g. "./dir/file.ext" or "dir/file.ext"), then the swd is prepended to path.
    /// 4) To resolve drive ambiguities, if swd provides a drive, it is used. If not, 
    ///    then the path drive is used. If neither provides a drive, then the current 
    ///    drive is used.
    static void deconstructPathnameUsingSpecifiedWorkingDirectory(const std::string& swd,
                                                                  const std::string& path,
                                                                  std::string& directory,
                                                                  std::string& fileName,
                                                                  std::string& extension);

    /// Get canonicalized absolute pathname from a given pathname which 
    /// can be relative or absolute. Canonicalizing means
    ///   - drive designator is recognized if we're on Windows;
    ///   - leading "." and "@" are replaced with the current working
    ///     directory or the executable directory, resp.
    ///   - each ".." segment is processed, removing it and its
    ///     previous segment; initial ".." is treated as "./..".
    ///   - empty segments and interior "." segments are removed
    ///   - if the input pathname ends in a slash after above processing,
    ///     then the returned pathname will also end in a slash.
    ///   - separators are made all-forward slash or all-backslash
    ///   - on Windows, the returned pathname begins with an explicit
    ///     disk designator in lower case, e.g. "c:".
    ///
    /// The result here is what you get by reassembling the components
    /// from deconstructPathname(), plus inserting the current working
    /// directory in front if the path name was relative.
    static std::string getAbsolutePathname(const std::string& pathname) {
        bool isAbsolutePath;
        std::string directory, fileName, extension;
        deconstructPathname(pathname, isAbsolutePath, directory, fileName, extension);
        if (!isAbsolutePath)
            directory = getCurrentWorkingDirectory() + directory;
        return directory + fileName + extension;
    }

    /// This is the same as getAbsolutePathname() except that the final
    /// segment is interpreted as a directory name rather than a file name,
    /// meaning that we append a slash if necessary.
    static std::string getAbsoluteDirectoryPathname(const std::string& dirPathname) {
        std::string absPath = getAbsolutePathname(dirPathname);
        if (!absPath.empty() && absPath[absPath.size()-1] != getPathSeparatorChar())
            absPath += getPathSeparatorChar();
        return absPath;
    }

    static std::string findAbsolutePathnameUsingSpecifiedWorkingDirectory(const std::string& swd,
                                                                          const std::string& path) {
        std::string directory, fileName, extension;
        deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
        return directory + fileName + extension;
    }

    /// Return true if the given pathname names a file that exists and is
    /// readable.
    static bool fileExists(const std::string& fileName);

    /// Get the default installation directory for this platform. This will
    /// be /usr/local/ for Linux and Apple, and the value of the \%ProgramFiles\%
    /// registry entry on Windows (typically c:\\Program Files\\).
    static std::string getDefaultInstallDir();

    /// Append a subdirectory offset to an existing pathname (relative or absolute).
    /// A leading "/" in the offset is ignored, and the result ends in "/".
    static std::string addDirectoryOffset(const std::string& base, const std::string& offset);

    /// Find the installation directory for something, using the named
    /// installation directory environment variable if it exists, otherwise
    /// by appending the supplied path offset to the default install directory.
    static std::string getInstallDir(const std::string& envInstallDir,
                                     const std::string& offsetFromDefaultInstallDir);

    /// Get the absolute pathname of the currently executing program.
    static std::string getThisExecutablePath();
    /// Get the absolute pathname of the directory which contains the 
    /// currently executing program.
    static std::string getThisExecutableDirectory();
    /// Get the absolute pathname of the current working directory
    /// including a trailing separator character. Windows keeps a current
    /// working directory for each drive which can be optionally specified
    /// (otherwise we use the current drive). If the specified drive 
    /// doesn't exist we'll behave as though root were its current
    /// working directory. The drive argument is ignored
    /// on non-Windows platforms.
    static std::string getCurrentWorkingDirectory(const std::string& drive="");
    /// Get the canonicalized name of the root directory. This is "x:\" on Windows
    /// with "x" replaced by the current drive letter or the specified drive
    /// (in lowercase), and just "/" on non-Windows systems.
    static std::string getRootDirectory(const std::string& drive="");
    /// On Windows, return the current drive letter in lowercase, with no
    /// trailing ":"; on other platforms return an empty string.
    static std::string getCurrentDriveLetter();
    /// On Windows, return the current drive letter in lowercase, 
    /// followed by ":"; on other platforms just return an empty string.
    static std::string getCurrentDrive();
    /// Return true if the named environment variable is present
    /// in the environment.
    static bool environmentVariableExists(const std::string& name);
    /// Return the value of the named environment variable or 
    /// the empty string if the variable is not found. Note that that
    /// is indistinguishable from a variable that is present but with
    /// a null value -- use environmentVariableExists() if you really
    /// need to know the difference.
    static std::string getEnvironmentVariable(const std::string& name);
    /// Return this platform's pathname separator character as
    /// a string. This is backslash on Windows and forward slash
    /// everywhere else.
    static std::string getPathSeparator();
    /// Return this platform's pathname separator character as
    /// a char. This is backslash on Windows and forward slash
    /// everywhere else.
    static char getPathSeparatorChar();
    /// Returns true if the character is slash or backslash.
    static bool isPathSeparator(char c) {
        return c=='/' || c=='\\';
    }
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PATHNAME_H_


