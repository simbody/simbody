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
        addSearchDirectory(Pathname::getInstallDir("SimTK_INSTALL_DIR", "SimTK") 
                            + "/lib/plugins/");
    }

    SimTK_PLUGIN_DEFINE_FUNCTION1(int,sayHi,const std::string&);
    SimTK_PLUGIN_DEFINE_FUNCTION(ExportedClass*,TestRuntimeDLL_makeExportedClass);
    SimTK_PLUGIN_DEFINE_FUNCTION(void, NotPresent);

private:
};

void testDeconstructFileName() {
    string name, directory, libPrefix, baseName, debugSuffix, extension;
    bool isAbsPath;
    const std::string s = Pathname::getPathSeparator();
    const std::string d = Pathname::getCurrentDriveLetter();
    const std::string dd = d.empty() ? std::string() : d + ":";
    const std::string cwd = Pathname::getCurrentWorkingDirectory();
    const std::string xd = Pathname::getThisExecutableDirectory();


    //printf("'%s': %s %s|%s|%s|%s|%s\n", name.c_str(),
    //    isAbsPath?"ABS":"REL", directory.c_str(), libPrefix.c_str(), baseName.c_str(), 
    //    debugSuffix.c_str(), extension.c_str());


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
    // Disk specifiers are only recognized on Windows.
    directory=libPrefix=baseName=debugSuffix=extension="junk";
    name = "  c:\\Program Files\\lib\\libMyPlugIn_d.dll \n ";   // OK
    SimTK_TEST(Plugin::deconstructLibraryName(name,
        isAbsPath, directory, libPrefix, baseName, debugSuffix, extension));
    SimTK_TEST(isAbsPath && 
        directory==("c:"+s+"Program Files"+s+"lib"+s) && libPrefix=="lib"
        && baseName=="MyPlugIn" && debugSuffix=="_d" && extension==".dll");

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
    SimTK_TEST(Pathname::getPathSeparatorChar()=='\\');
    SimTK_TEST(Pathname::getPathSeparator()=="\\");
    SimTK_TEST(Pathname::getCurrentDriveLetter().size() == 1);
    SimTK_TEST(Pathname::getCurrentDrive()==Pathname::getCurrentDriveLetter()+":");
    SimTK_TEST(Pathname::getRootDirectory() == Pathname::getCurrentDrive() + "\\");
#else
    SimTK_TEST(Pathname::getPathSeparatorChar()=='/');
    SimTK_TEST(Pathname::getPathSeparator()=="/");
    SimTK_TEST(Pathname::getCurrentDriveLetter().size() == 0);
    SimTK_TEST(Pathname::getCurrentDrive()=="");
    SimTK_TEST(Pathname::getRootDirectory() == "/");
#endif

    std::string name;
    bool isAbsPath;
    std::string directory, fileName, extension;
    const std::string curDrive = 
        Pathname::getCurrentDriveLetter().empty() 
            ? std::string() 
            : Pathname::getCurrentDriveLetter()+":";
    const std::string sep = Pathname::getPathSeparator();

    directory=fileName=extension="junk"; isAbsPath=false;
    name = "/topdir/seconddir/myFileName.ext";
    Pathname::deconstructPathname(name, isAbsPath, directory, fileName, extension);
    SimTK_TEST(isAbsPath
        && directory==curDrive+sep+"topdir"+sep+"seconddir"+sep
        && fileName=="myFileName" && extension==".ext");

    directory=fileName=extension="junk"; isAbsPath=true;
    name = "topdir/seconddir/myFileName.ext";
    Pathname::deconstructPathname(name, isAbsPath, directory, fileName, extension);
    SimTK_TEST(!isAbsPath
        && directory=="topdir"+sep+"seconddir"+sep
        && fileName=="myFileName" && extension==".ext");
}

void testPlugin() {
    MyPlugin myPlug("c:\\temp\\TestRuntimeDLL_d.dll");

    if (myPlug.load()) {
        std::cout << "fellas length is " << myPlug.sayHi("fellas") << std::endl;
        myPlug.TestRuntimeDLL_makeExportedClass();
    } else
        std::cout << "Plugin load failed, err=" << myPlug.getLastErrorMessage() << std::endl;

    SimTK_TEST_MUST_THROW(myPlug.NotPresent());

    if (myPlug.isLoaded()) {
        SimTK_TEST(myPlug.has_TestRuntimeDLL_makeExportedClass());
        SimTK_TEST(!myPlug.has_NotPresent());
    }

    myPlug.unload();
    myPlug.addSearchDirectory("@");
    myPlug.addSearchDirectory(".");
    myPlug.addSearchDirectory("..");
    myPlug.addSearchDirectory("/temp");
    myPlug.addSearchDirectory("c:/temp");
    myPlug.load("TestRuntimeDLL.so");
}

int main() {
    cout << "Path of this executable: '" << Pathname::getThisExecutablePath() << "'\n";
    cout << "Executable directory: '" << Pathname::getThisExecutableDirectory() << "'\n";
    cout << "Current working directory: '" << Pathname::getCurrentWorkingDirectory() << "'\n";
    cout << "Current drive letter: '" << Pathname::getCurrentDriveLetter() << "'\n";
    cout << "Current drive: '" << Pathname::getCurrentDrive() << "'\n";
    cout << "Path separator: '" << Pathname::getPathSeparator() << "'\n";
    cout << "Default install dir: '" << Pathname::getDefaultInstallDir() << "'\n";
    cout << "installDir(SimTK_INSTALL_DIR,/LunchTime/snacks): '"
        << Pathname::getInstallDir("SimTK_INSTALL_DIR", "/LunchTime/snacks") << "'\n";
    cout << "installDir(LocalAppData,SimTK): '"
        << Pathname::getInstallDir("LocalAppData", "SimTK") << "'\n";
    SimTK_START_TEST("TestPlugin");

        SimTK_SUBTEST(testPathname);
        SimTK_SUBTEST(testDeconstructFileName);
        SimTK_SUBTEST(testPlugin);

    SimTK_END_TEST();
}

