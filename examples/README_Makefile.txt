Instructions for using the Makefile provided.

----------------------------------------------------------------------------
NOTE: By far the easiest way to create a build project is to use
CMake (www.cmake.org). We have supplied a CMakeLists.txt file here which 
will build a single Makefile containing all the examples.
That works on all platforms also, using whatever build system is appropriate,
even Visual Studio.  See CMakeNotes.txt.
----------------------------------------------------------------------------


You must already have the SimTK Core binaries installed from SimTK.org; go 
to https://simtk.org/home/simtkcore, Downloads tab.

If you want to use OpenMM acceleration, you must have installed that separately;
go to https://simtk.org/home/openmm for information and downloads.

You may need to slightly edit the Makefile to make it run on your system.

Open the Makefile and make the appropriate changes (debug libraries or not), and set the Default install directory to the correct location.

To compile all example programs type "make all"

To compile just ExampleAdenylateMobilitiesVTK (for example) type make ExampleAdenylateMobilitiesVTK.

Before you run the executables, remember to add the SimTK Core library directory to your library path. The simplest way to do this is to type the following commands:

For Linux (for the bash shell, assuming installation was done in the default location: /usr/local/SimTK): 
$ export LD_LIBRARY_PATH=/usr/local/SimTK/lib 

For Mac (for the bash shell, assuming installation was done in the default location: /usr/local/SimTK):
$ export DYLD_LIBRARY_PATH=/usr/local/SimTK/lib
