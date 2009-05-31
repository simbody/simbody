#ifndef SimTK_SimTKCOMMON_PLUGIN_H_
#define SimTK_SimTKCOMMON_PLUGIN_H_

/* -------------------------------------------------------------------------- *
 *                          SimTK Core: SimTKcommon                           *
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

#include "SimTKcommon/internal/common.h"
#include <string>
#include <vector>

namespace SimTK {

/**
 * This class encapsulates the handling of file and directory pathnames
 * in a platform-indpendent manner. We consider a pathname to consist
 * of three components:
 *      [directory] [filename [extension]]
 * where the directory may be an absolute location or relative to a
 * current working directory. 
 *
 * Several special directory names are supported here:
 *  - root (/)
 *  - current working directory (.)
 *  - current executable directory (@)
 *  - platform default installation directory
 *  - parent directory (..)
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
 * we know to e a directory, it will end in a final slash.
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
    /// e:\\Program Files\\Something\\  myLibrary_d    .dll
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
    /// Anywhere in the pathname, the name ".." means "go up one level from
    /// the prior directory. ".." at the start is interpreted as "./..".
    /// A "." appearing anywhere in the path name except the begining is
    /// ignored. An "@" appearing anywhere in the pathname other than the
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
                                        bool&        isAbsolutePath,
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

    /// Get the default installation directory for this platform. This will
    /// be /usr/local/ for Linux and Apple, and the value of the %ProgramFiles%
    /// registry entry on Windows (typically c:\Program Files\).
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
    /// the empty string if the variable is not found. (Note that that
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

/**
 * This is the base class for representing a runtime-linked dynamic library,
 * also known as a "plugin", in a platform-independent manner. For any particular 
 * kind of plugin, derive a concrete class from this one and use the macros 
 * to describe the functions or data you expect to find there.
 */
class SimTK_SimTKCOMMON_EXPORT Plugin {
public:
    explicit Plugin(const std::string& defaultPathname="");
    ~Plugin();  // unloads the plugin if it was loaded

    /// Attempt to load a plugin of the given name if any, otherwise the
    /// default name (if any). If the name is relative it will be subject to
    /// the search rule associated with this Plugin object. In any case
    /// the filename will be deconstructed to find the library baseName
    /// to which we will add the "lib" prefix if necessary, the "_d"
    /// suffix if appropriate, and the platform-specific suffix.
    ///
    /// load() returns true if it succeeds, otherwise you
    /// can call getLastErrorMessage() to find out what's wrong. If a library
    /// is already loaded, this returns false without checking to see whether
    /// you are trying to reload the same library. Call isLoaded() to check
    /// whether a plugin has already been loaded, getLoadedPathname() to
    /// get a canonicalized absolute pathname for the loaded library.
    bool load(const std::string& name="");

    /// If a plugin is loaded, unload it now otherwise do nothing. Note that
    /// whether the corresponding DLL is actually unloaded from the address
    /// space is up to the operating system, however, this Plugin object will
    /// be disassociated from it regardless.
    void unload();

    /// Is there currently a DLL associated with this Plugin object? If so you
    /// can call getLoadedPathname() to find out which one.
    bool isLoaded()                   const {return m_handle != 0;}

    const std::string& getLoadedPathname() const {
        return m_loadedPathname; // empty if nothing loaded
    }

    /// If anything goes wrong the last error message is stored so you can
    /// retrieve it with this method.
    std::string getLastErrorMessage() const {return m_lastMessage;}

    /// Provide a list of directories in the order you want them
    /// searched to find a plugin that was provided as a relative path name.
    /// This is ignored if an absolute path name is given.
    /// Each of the directory names will be immediately expanded to an
    /// absolute path name if it isn't already.
    void setSearchPath(const std::vector<std::string>& pathIn) {
        m_searchPath.clear();
        for (unsigned i=0; i < pathIn.size(); ++i) 
            addSearchDirectory(pathIn[i]);
    }

    const std::vector<std::string>& getSearchPath() const {return m_searchPath;}

    /// Add a directory to the end of the search path for this kind of
    /// plugin. This will be expanded immediately to an absolute path
    /// name if it isn't already. If the directory is blank it is ignored.
    void addSearchDirectory(const std::string& directory) {
        const std::string absDir = Pathname::getAbsoluteDirectoryPathname(directory);
        if (!absDir.empty())
            m_searchPath.push_back(absDir);
    }

    /// Put a directory on the front of the search path for this kind of
    /// plugin. This will be expanded immediately to an absolute path
    /// name if it isn't already. If the directory name is blank it is
    /// ignored.
    void prependSearchDirectory(const std::string& directory) {
        const std::string absDir = Pathname::getAbsoluteDirectoryPathname(directory);
        if (!absDir.empty())
            m_searchPath.insert(m_searchPath.begin(), absDir);
    }


    /// If this fails the return value will be null and the system's human-readable
    /// error message is in errMsg.
    static void* loadPluginByFileName(const std::string& name, std::string& errMsg);
    /// If we're in Debug mode then this method attempts first to load the Debug
    /// version of the indicated library which it constructs as base+"_d"+extension.
    /// If that fails (or if we're not in Debug mode) it will try to load the 
    /// Release version (base+extension) instead.
    static void* loadDebugOrReleasePlugin(const std::string& base, const std::string& extension,
                                          std::string& loadedFileName, std::string& errMsg);
    /// If this fails the return value will be null and the system's human-readable
    /// error message is in errMsg.
    static void* getSymbolAddress(void* handle, const std::string& name, std::string& errMsg);
    /// Given a handle returned by loadPluginByFileName(), unload the plugin. Nothing
    /// happens until the last reference is unloaded.
    static void unloadPlugin(void* handle);

    /// Dismantle a supplied library's file or pathname into its component
    /// parts. This can take pathnames like <pre>   
    ///     /usr/local/libMyDll_d.so
    ///     e:\\Program Files\\Something\\myLibrary_d.dll
    /// </pre> and chop them into <pre>
    /// directory                       libPrefix baseName    debug extension
    /// ------------------------------- --------- ----------- ----- ---------
    /// /usr/local/                     lib       MyDll         _d    .so 
    /// e:\\Program Files\\Something\\     (none)    myLibrary     _d    .dll
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
    /// Anywhere in the pathname, the name ".." means "go up one level from
    /// the prior directory. ".." at the start is interpreted as "./..".
    /// A "." appearing anywhere in the path name except the begining is
    /// ignored. An "@" appearing anywhere in the pathname other than the
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
    /// The return value is false if the input path name is
    /// ill-formed; in that case it still tries to parse as
    /// much as it can.
    static bool deconstructLibraryName( const std::string& name,
                                        bool&        isAbsolutePath,
                                        std::string& directory,
                                        std::string& libPrefix,
                                        std::string& baseName,
                                        std::string& debugSuffix,
                                        std::string& extension);

    /// This is the platform dependent string which gets prepended to
    /// a dynamic library baseName to form the fileName -- "lib" on Unix-based
    /// systems and "" (empty string) for Windows.
    static std::string getDynamicLibPrefix();

    /// This is the platform dependent extension (including the ".") used 
    /// by default to identify dynamically linked libraries -- ".so" on 
    /// Linux, ".dylib" on Apple, and ".dll" on Windows.
    static std::string getDynamicLibExtension();

    /// Obtain the appropriate debug suffix to use. This is not platform
    /// dependent but rather depends on whether this compilation unit
    /// was built in Debug or Release modes. If Debug, then the string
    /// "_d" is returned, otherwise the empty string.
    static std::string getDynamicLibDebugSuffix();

protected:
    std::string                 m_defaultName; // if any
    std::vector<std::string>    m_searchPath;

    std::string                 m_loadedPathname; // absolute
    void*                       m_handle;
    mutable std::string         m_lastMessage;

};


#define SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName)          \
    struct FuncName##__Holder__ {                       \
        FuncName##__Holder__() : fp(0) {}               \
        bool loadSym(void* h, std::string& msg) const { \
            if(!fp) fp =(FuncName##__Type__)            \
                Plugin::getSymbolAddress(h, #FuncName, msg);   \
            return (fp!=0);                             \
        }                                               \
        mutable FuncName##__Type__ fp;                  \
    } FuncName##__Ref__
#define SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)            \
    if (!FuncName##__Ref__.loadSym(m_handle,m_lastMessage)) \
    throw std::runtime_error                            \
      ("Plugin function " #FuncName " not found: " + m_lastMessage); \
    return FuncName##__Ref__.fp
#define SimTK_PLUGIN_XXX_MAKE_SYMTEST(Symbol)           \
    bool has_##Symbol() const {                          \
        return Symbol##__Ref__.loadSym(m_handle,m_lastMessage);   \
    }

#define SimTK_PLUGIN_DEFINE_SYMBOL(Type, SymName)   \
    typedef Type SymName##__Type__;                 \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(SymName);          \
    const Type& SymName() const {                   \
        if (!SymName##__Ref__.loadSym(m_handle,m_lastMessage))  \
        throw std::runtime_error                                \
          ("Plugin symbol " #SymName " not found: " + m_lastMessage); \
        return *(SymName##__Ref__.fp);              \
    }

#define SimTK_PLUGIN_DEFINE_FUNCTION(RetType, FuncName) \
    typedef RetType (*FuncName##__Type__)();            \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName() const {                          \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)();         \
    }                                                   \
    SimTK_PLUGIN_XXX_MAKE_SYMTEST(FuncName)

#define SimTK_PLUGIN_DEFINE_FUNCTION1(RetType, FuncName, Arg1) \
    typedef RetType (*FuncName##__Type__)(Arg1);        \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName(Arg1 a1) const {                   \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)(a1);       \
    }                                                   \
    SimTK_PLUGIN_XXX_MAKE_SYMTEST(FuncName)

#define SimTK_PLUGIN_DEFINE_FUNCTION2(RetType, FuncName, Arg1, Arg2) \
    typedef RetType (*FuncName##__Type__)(Arg1,Arg2);   \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName(Arg1 a1, Arg2 a2) const {          \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)(a1,a2);    \
    }                                                   \
    SimTK_PLUGIN_XXX_MAKE_SYMTEST(FuncName)


} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PLUGIN_H_
