/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include "SimTKcommon/internal/Plugin.h"

#include <iostream>
#include <string>
#include <cstdio>
using std::cout;
using std::endl;
using std::string;
using std::printf;

using namespace SimTK;

class ExportedClass;
class MyPlugin : public Plugin {
public:
    explicit MyPlugin(const std::string& name)
    :   Plugin(name) {
        std::string defaultLocation = getEnvironmentVariable("SimTK_INSTALL_DIR");
        if (defaultLocation.empty()) {
            std::string pfiles = getEnvironmentVariable("PROGRAMFILES");
            if (pfiles.empty())
                pfiles = "c:/Program Files";
            defaultLocation = pfiles + "/SimTK";
        }

        defaultLocation += "/core/lib/plugins/";
        addSearchDirectory(defaultLocation);
    }

    SimTK_DEFINE_PLUGIN_FUNCTION1(int,sayHi,const std::string&);
    SimTK_DEFINE_PLUGIN_FUNCTION(ExportedClass*,TestRuntimeDLL_makeExportedClass);

private:
};

class MyThingInterface : public ManagedPluginInterface {
public:
    std::string type() const {return "MyThing";}
    std::string version() const {return "0.0";}

    class Thing {
    public:
        Thing(const std::string& s) : thingString(s) {}
        std::string thingString;
    };

    virtual Thing* createThing(const std::string&) const = 0;
private:
};

//class MyPlugin2 : public Plugin {
//public:
//    explicit MyPlugin(const std::string& name)
//    :   Plugin() {
//        setInstallDirFromEnvVar("SimTK_INSTALL_DIR", "core");
//        setInstallDirFromDefault("SimTK/core");
//        setOffsetFromInstallDir("lib/plugins");
//
//        std::string defaultLocation = getEnvironmentVariable("SimTK_INSTALL_DIR");
//        if (defaultLocation.empty()) {
//            std::string pfiles = getEnvironmentVariable("PROGRAMFILES");
//            if (pfiles.empty())
//                pfiles = "c:/Program Files";
//            defaultLocation = pfiles + "/SimTK";
//        }
//
//        defaultLocation += "/core/lib/plugins/";
//        addSearchDirectory(defaultLocation);
//    }
//
//    SimTK_DEFINE_PLUGIN_FUNCTION1(int,sayHi,const std::string&);
//    SimTK_DEFINE_PLUGIN_FUNCTION(ExportedClass*,TestRuntimeDLL_makeExportedClass);
//
//private:
//};

void testDeconstructFileName() {
    string name, directory, libPrefix, baseName, debugSuffix, extension;
    bool isAbsPath;
    const std::string s = Plugin::getPathSeparator();
    const std::string d = Plugin::getCurrentDriveLetter();
    const std::string dd = d.empty() ? std::string() : d + ":";
    const std::string cwd = Plugin::getCurrentWorkingDirectory();
    const std::string xd = Plugin::getThisExecutableDirectory();


    //printf("'%s': %s %s|%s|%s|%s|%s\n", name.c_str(),
    //    isAbsPath?"ABS":"REL", directory.c_str(), libPrefix.c_str(), baseName.c_str(), 
    //    debugSuffix.c_str(), extension.c_str());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "  c:\\Program Files\\lib\\libMyPlugIn_d.dll \n ";   // OK
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath && 
        directory==("c:"+s+"Program Files"+s+"lib"+s) && libPrefix=="lib"
        && baseName=="MyPlugIn" && debugSuffix=="_d" && extension==".dll");

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "  \t \n \r ";   // Illegal because nothing but white space
    SimTK_TEST(!Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(!isAbsPath && directory.empty()&&libPrefix.empty()&&baseName.empty()
               &&debugSuffix.empty()&&extension.empty());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "   This.Is.Not.A.Suffix_A.dylib  "; // OK
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(!isAbsPath && directory.empty() && libPrefix.empty() 
        && baseName=="This.Is.Not.A.Suffix_A"
        && debugSuffix.empty() && extension==".dylib");

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "/usr/local/lib/"; // illegal because no file name
    SimTK_TEST(!Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath 
        && directory==(dd+s+"usr"+s+"local"+s+"lib"+s) 
        && libPrefix.empty()&&baseName.empty()
        && debugSuffix.empty()&&extension.empty());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "lib_D"; // illegal because nothing left after lib prefix and debug suffix
    SimTK_TEST(!Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(!isAbsPath && directory.empty()&&libPrefix=="lib"&&baseName.empty()
               &&debugSuffix=="_D"&&extension.empty());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "/mylibrary_D."; // OK (empty file type suffix)
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath 
        && directory==(dd+s) && libPrefix.empty()&&baseName=="mylibrary"
        && debugSuffix=="_D"&&extension==".");

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "./first/more/../../filename"; 
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath 
        && directory==cwd
        && libPrefix.empty()&&baseName=="filename"
        && debugSuffix.empty()&&extension.empty());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "@/../../filename"; // executable path-relative name
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));

    SimTK_TEST(isAbsPath 
        //&& directory==xd
        && libPrefix.empty()&&baseName=="filename"
        && debugSuffix.empty()&&extension.empty());

    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "a"; // OK (very short file name)
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(!isAbsPath && directory.empty()&&libPrefix.empty()&&baseName=="a"
               &&debugSuffix.empty()&&extension.empty());

#ifdef _WIN32
    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = dd + "mydir/relative.txt"; // cwd relative to specified disk
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath
        && directory==(cwd + "mydir" + s)
        && libPrefix.empty()&&baseName=="relative"
        && debugSuffix.empty()&&extension==".txt");
#endif
}

void testPathname() {
#ifdef _WIN32
    SimTK_TEST(Plugin::getPathSeparatorChar()=='\\');
    SimTK_TEST(Plugin::getPathSeparator()=="\\");
    SimTK_TEST(Plugin::getCurrentDriveLetter().size() == 1);
    SimTK_TEST(Plugin::getCurrentDrive()==Plugin::getCurrentDriveLetter()+":");
    SimTK_TEST(Plugin::getRootDirectory() == Plugin::getCurrentDrive() + "\\");
#else
    SimTK_TEST(Plugin::getPathSeparatorChar()=='/');
    SimTK_TEST(Plugin::getPathSeparator()=="/");
    SimTK_TEST(Plugin::getCurrentDriveLetter().size() == 0);
    SimTK_TEST(Plugin::getCurrentDrive()=="");
    SimTK_TEST(Plugin::getRootDirectory() == "/");
#endif

    std::string name;
    bool isAbsPath;
    std::string directory, fileName, extension;
    const std::string curDrive = 
        Plugin::getCurrentDriveLetter().empty() 
            ? std::string() 
            : Plugin::getCurrentDriveLetter()+":";
    const std::string sep = Plugin::getPathSeparator();

    directory=fileName=extension="junk"; isAbsPath=false;
    name = "/topdir/seconddir/myFileName.ext";
    Plugin::deconstructPathname(name, isAbsPath, directory, fileName, extension);
    SimTK_TEST(isAbsPath
        && directory==curDrive+sep+"topdir"+sep+"seconddir"+sep
        && fileName=="myFileName" && extension==".ext");

    directory=fileName=extension="junk"; isAbsPath=true;
    name = "topdir/seconddir/myFileName.ext";
    Plugin::deconstructPathname(name, isAbsPath, directory, fileName, extension);
    SimTK_TEST(!isAbsPath
        && directory=="topdir"+sep+"seconddir"+sep
        && fileName=="myFileName" && extension==".ext");
}

void testPlugin() {
    MyPlugin myPlug("c:\\temp\\TestRuntimeDLL_d.dll");

    std::cout << "fellas length is " << myPlug.sayHi("fellas") << std::endl;


    myPlug.TestRuntimeDLL_makeExportedClass();
}

int main() {
    cout << "Path of this executable: '" << Plugin::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Plugin::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Plugin::getCurrentWorkingDirectory() << "'\n";
    cout << "Current drive letter: '" << Plugin::getCurrentDriveLetter() << "'\n";
    cout << "Current drive: '" << Plugin::getCurrentDrive() << "'\n";
    cout << "Path separator: '" << Plugin::getPathSeparator() << "'\n";
    SimTK_START_TEST("TestPlugin");

        SimTK_SUBTEST(testPathname);
        SimTK_SUBTEST(testDeconstructFileName);
        SimTK_SUBTEST(testPlugin);

    SimTK_END_TEST();
}

