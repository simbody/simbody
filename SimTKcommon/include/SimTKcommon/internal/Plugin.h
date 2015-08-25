#ifndef SimTK_SimTKCOMMON_PLUGIN_H_
#define SimTK_SimTKCOMMON_PLUGIN_H_

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
Declaration of the SimTK::Plugin class providing platform-independent 
handling of dynamically-loaded libraries. **/

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Array.h"
#include "SimTKcommon/internal/Pathname.h"
#include <string>
#include <stdexcept>

namespace SimTK {

/**
 * This is the base class for representing a runtime-linked dynamic library,
 * also known as a "plugin", in a platform-independent manner. For any particular 
 * kind of plugin, derive a concrete class from this one and use the macros 
 * to describe the functions or data you expect to find there. Then each 
 * plugin library you load of that type is an object of the concrete class
 * you defined. The derived class constructor sets the search policy for
 * plugin libraries of that type. For example:
 * <pre>
 * class MyPlugin : public Plugin {
 * public:
 *     explicit MyPlugin() : Plugin() {
 *         addSearchDirectory(Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK") 
 *                             + "/lib/plugins/");
 *     }
 * 
 *     SimTK_PLUGIN_DEFINE_FUNCTION1(int,exportedFunctionWithOneArg,const std::string&);
 *     SimTK_PLUGIN_DEFINE_FUNCTION(SomeObjectType*,anExportedFunctionNoArgs);
 *     SimTK_PLUGIN_DEFINE_FUNCTION(void, someFunction);
 *     SimTK_PLUGIN_DEFINE_SYMBOL(SymbolType, nameOfExportedSymbol);
 * };
 * </pre>
 *
 * Then this class is used as follows:
 * <pre>
 *     MyPlugin lib1;
 *     if (!lib1.load("libraryName")) {
 *         std::cerr << lib1.getLastErrorMessage() << std::endl;
 *         // exit, retry or whatever
 *     }
 *     // At this point you can call the functions and access
 *     // the symbols, e.g.
 *     int i = lib1.exportedFunctionWithOneArg("hi");
 * </pre>
 *
 * The currently available macros are:
 *  - SimTK_PLUGIN_DEFINE_SYMBOL(type, name)
 *  - SimTK_PLUGIN_DEFINE_FUNCTION(returnType, name)
 *  - SimTK_PLUGIN_DEFINE_FUNCTION1(returnType, name, typeOfArg)
 *  - SimTK_PLUGIN_DEFINE_FUNCTION2(returnType, name, typeOfArg1, typeOfArg2)
 *
 * These define methods in the derived class which can obtain the symbol
 * address from the plugin library as needed. The methods will have the
 * same signature as the library's exported function. In the case of the
 * symbol, a no-argument method named the same as the symbol will be
 * created. For example,
 * <pre> SimTK_PLUGIN_DEFINE_SYMBOL(int, myCounter); </pre>
 * says the library contains an exported external symbol of type int
 * named "myCounter". If the Plugin object is "plugin" then you can
 * access this symbol (read only) as
 * <pre> const int& counter = plugin.myCounter(); </pre>
 * The macros also define another method of the same name as the symbol
 * or function, but prepended with "has_". This returns true if the 
 * named object is exported by the library, otherwise false. If you try
 * to access a method or symbol without checking first, an exception
 * will be thrown if the library doesn't export the symbol.
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
    void setSearchPath(const Array_<std::string>& pathIn) {
        m_searchPath.clear();
        for (unsigned i=0; i < pathIn.size(); ++i) 
            addSearchDirectory(pathIn[i]);
    }

    const Array_<std::string>& getSearchPath() const {return m_searchPath;}

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
    Array_<std::string>         m_searchPath;

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


