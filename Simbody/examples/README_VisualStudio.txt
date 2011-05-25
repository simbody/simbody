Instructions for creating a Visual Studio project by hand to run a 
Simbody example program.
Simbody 2.2, June 2011

----------------------------------------------------------------------------
NOTE: By far the easiest way to create a Visual Studio project is to use
CMake (www.cmake.org). We have supplied a CMakeLists.txt file here which 
will build a single Visual Studio project containing all the examples.
That works on all platforms also, using whatever build system is appropriate.
----------------------------------------------------------------------------

But, if you insist on building a Visual Studio project by hand, here are
the instructions.



Below we will use $(SimTK_INSTALL_DIR) to mean the SimTK installation 
directory, which is typically %ProgramFiles%\SimTK but can be anything.
(On English systems this will be C:\Program Files\SimTK or 
C:\Program Files (x86)\SimTK for 64 bit platforms.). You
should create a SimTK_INSTALL_DIR environment variable and
set it to the appropriate Simbody installation directory.

You must already have the Simbody binaries installed from SimTK.org; go 
to https://simtk.org/home/simbody, Downloads tab.

To use the Visual Studio solution file provided in this directory proceed 
to step 5)

To create a new Visual Studio project for a Simbody example from scratch
follow all steps:


1) Create new Visual Studio project based on an example source file:

  * In Visual Studio, select File->New->Project...->Visual C++->Win32->Win32 Console Application

  * Select a name for the new project, select a location, and select "Create new Solution"

  * Click "Next"; select "Console Application" and "Empty project"

  * Drag the desired example .cpp program into the "Source files" folder on the left.


2) Set up include directory

  * Far click the Project name on the left bar, select "Properties"

  * Set the "Configuration:" pulldown to "All Configurations"

  * Under Configuration Properties->C++->General select "Additional Include Directories"

  * Select the "..." button, then the new directory (folder) button

  * Browse to $(SimTK_INSTALL_DIR)/include, and add it to the "Additional Include Directories" section

  * Click "OK"


3) Set up library directories

  * Far click the Project name on the left bar, select "Properties"

  * Set the "Configuration:" pulldown to "All Configurations"

  * Under Configuration Properties->Linker->General select "Additional Library Directories"

  * Select the "..." button, then the new directory (folder) button

  * Add $(SimTK_INSTALL_DIR)/lib to the "Additional Library Directories" section.

  * Click "OK"


4) Set up library names

  * Far click the Project name on the left bar, select "Properties"

  a) Debug versions
     NOTE: Debug libraries are very sensitive to compiler version; if you can't
     use the prebuilt ones follow the Release instructions below unless you 
     have built your own Debug libraries from source.
     CAUTION: You may not be able to use Release libraries with a main 
     program compiled with Debug.

    * Set the "Configuration:" pulldown to "Debug"

    * Under Configuration Properties->Linker->Input select "Additional Dependencies"

    * Select the "..." button
  
    * Add the following library names:
       SimTKsimbody_d.lib
       SimTKmath_d.lib
       SimTKcommon_d.lib
       SimTKlapack.lib
     Note that there is no "_d" lapack library; that's OK.

  * Click "Apply"

  a) Release versions

    * Set the "Configuration:" pulldown to "Release"

    * Under Configuration Properties->Linker->Input select "Additional Dependencies"

    * Select the "..." button
  
    * Add the following library names:
       SimTKsimbody.lib
       SimTKmath.lib
       SimTKcommon.lib
       SimTKlapack.lib

  * Click "Apply"


5) Run the program

* select "Release" configuration in the Visual Studio main toolbar if you 
  want speed; Debug is typically 10X slower.

* Type Ctrl-F5, or click Debug->"Start Without Debugging"

* Your program should run

* Type Ctrl-C in the console (Command Prompt) window to stop the program


6) To run another example program:
  
  * Delete the previous example .cpp file from the "Source Files" folder

  * Drag the new example .cpp file into the "Source Files" folder
  
  * Select "Release" configuration, Ctrl-F5, etc.

