/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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

#ifdef _WIN32
    std::string swd, path, dir, pathname;
    const std::string cwd = Pathname::getCurrentWorkingDirectory();
    std::string cwd_nodrive = cwd; cwd_nodrive.erase(0,3);
    const std::string cwdC = Pathname::getCurrentWorkingDirectory("c");
    const std::string cwdD = Pathname::getCurrentWorkingDirectory("d");
    const std::string thisExeDir = Pathname::getThisExecutableDirectory();

    swd = "   "; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "/topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "  "; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwdC + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = " "; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "D:/specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:/specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = "d:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:/specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = "d:" + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:/specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = "d:" + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:/specified"; path = "topdir/seconddir/myFileName.ext";
    dir = "d:" + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:/specified"; path = "../seconddir/myFileName.ext";
    dir = "d:" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "/specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "D:specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = "d:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwdD + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwdD + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwdD + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "./specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = "c:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = curDrive + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwdC + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "specified"; path = "topdir/../../seconddir/myFileName.ext";
    dir = cwd + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "";
    dir = "";
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = "."; path = "";
    dir = "";
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = ""; path = ".";
    dir = cwd;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = "@"; path = "topdir/seconddir/myFileName.ext";
    dir = thisExeDir + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

#else
    std::string swd, path, dir, pathname;
    const std::string cwd = Pathname::getCurrentWorkingDirectory();
    const std::string thisExeDir = Pathname::getThisExecutableDirectory();

    swd = ""; path = "C:/topdir/seconddir/myFileName.ext";
    dir = cwd + "C:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "/topdir/seconddir/myFileName.ext";
    dir = sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwd + "C:topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = ""; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "/specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = sep + "specified" + sep + "C:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = sep + "specified" + sep + "C:topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "topdir/seconddir/myFileName.ext";
    dir = sep + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "/specified"; path = "../seconddir/myFileName.ext";
    dir = sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "D:specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = cwd + "D:specified" + sep + "C:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwd + "D:specified" + sep + "C:topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "D:specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "D:specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "D:specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "./specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "C:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "X:topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "X:topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "./specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "specified"; path = "C:/topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "C:" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "/topdir/seconddir/myFileName.ext";
    dir = sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "C:topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "C:topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "./topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "topdir/seconddir/myFileName.ext";
    dir = cwd + "specified" + sep + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    ///

    swd = "specified"; path = "topdir/../../seconddir/myFileName.ext";
    dir = cwd + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");

    swd = "specified"; path = "";
    dir = "";
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = "."; path = "";
    dir = "";
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = ""; path = ".";
    dir = cwd;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "" && extension == "");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir);

    swd = "@"; path = "topdir/seconddir/myFileName.ext";
    dir = thisExeDir + "topdir" + sep + "seconddir" + sep;
    directory = fileName = extension = "junk";
    Pathname::deconstructPathnameUsingSpecifiedWorkingDirectory(swd, path, directory, fileName, extension);
    SimTK_TEST(directory == dir && fileName == "myFileName" && extension == ".ext");
    pathname = Pathname::getAbsolutePathnameUsingSpecifiedWorkingDirectory(swd, path);
    SimTK_TEST(pathname == dir + "myFileName.ext");
#endif

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

