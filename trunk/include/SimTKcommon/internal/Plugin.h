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
 * This is the base class for representing a runtime-linked dynamic library,
 * also known as a "plugin". For any particular kind of plugin, derive a
 * concrete class from this one and use the macros to describe the functions
 * or data you expect to find there.
 */

class SimTK_SimTKCOMMON_EXPORT Plugin {
public:
    Plugin(const std::string& name);
    ~Plugin();  // unloads the plugin if it was loaded

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
        const std::string absDir = getAbsoluteDirectoryPathname(directory);
        if (!absDir.empty())
            m_searchPath.push_back(absDir);
    }

    /// Put a directory on the front of the search path for this kind of
    /// plugin. This will be expanded immediately to an absolute path
    /// name if it isn't already. If the directory name is blank it is
    /// ignored.
    void prependSearchDirectory(const std::string& directory) {
        const std::string absDir = getAbsoluteDirectoryPathname(directory);
        if (!absDir.empty())
            m_searchPath.insert(m_searchPath.begin(), absDir);
    }


    bool isPluginLoaded()             const {return m_handle != 0;}
    std::string getLastErrorMessage() const {return m_lastMessage;}

    /// If this fails the return value will be null and the system's human-readable
    /// error message is in errMsg.
    static void* loadPluginByFileName(const std::string& name, std::string& errMsg);
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
    ///     disk designator, e.g. "c:".
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
protected:
    void*               m_handle;
    mutable std::string m_lastMessage;

    std::vector<std::string> m_searchPath;
};


#define SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName)          \
    struct FuncName##__Holder__ {                       \
        FuncName##__Holder__() : fp(0) {}               \
        mutable FuncName##__Type__ fp;                  \
    } FuncName##__Ref__
#define SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)            \
    if (FuncName##__Ref__.fp==0)                        \
        FuncName##__Ref__.fp = (FuncName##__Type__)     \
        Plugin::getSymbolAddress(m_handle, #FuncName, m_lastMessage);    \
    return FuncName##__Ref__.fp

#define SimTK_DEFINE_PLUGIN_FUNCTION(RetType, FuncName) \
    typedef RetType (*FuncName##__Type__)();            \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName() const {                          \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)();         \
    }

#define SimTK_DEFINE_PLUGIN_FUNCTION1(RetType, FuncName, Arg1) \
    typedef RetType (*FuncName##__Type__)(Arg1);        \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName(Arg1 a1) const {                   \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)(a1);       \
    }

#define SimTK_DEFINE_PLUGIN_FUNCTION2(RetType, FuncName, Arg1, Arg2) \
    typedef RetType (*FuncName##__Type__)(Arg1,Arg2);   \
    SimTK_PLUGIN_XXX_MAKE_HOLDER(FuncName);             \
    RetType FuncName(Arg1 a1, Arg2 a2) const {          \
        SimTK_PLUGIN_XXX_MAKE_BODY(FuncName)(a1,a2);    \
    }


/**
 * The plugin library's SimTK_registerPlugin() function allocates one
 * of these and passes it to the supplied PluginManager.
 */
class ManagedPluginInterface {
public:
    virtual ~ManagedPluginInterface() {}
    virtual std::string type() const = 0;
    virtual std::string version() const = 0;
private:
};

/**
 * 
 */
class PluginManager {
public:

};

/**
 * 
 */
template <class Interface>
class PluginManager_ : public PluginManager {
public:
    PluginManager_();

    bool registerInterface(Interface*) {
    }

    class ManagedPlugin : public Plugin {
    public:
        SimTK_DEFINE_PLUGIN_FUNCTION1(bool,SimTK_registerPlugin,PluginManager_&);
    };

private:
    static std::vector<Interface*> loadedPlugins;
};




} // namespace SimTK

#endif // SimTK_SimTKCOMMON_PLUGIN_H_
