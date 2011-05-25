                     Simbody 2.2 Examples README
                             June, 2011

This is an eclectic set of examples covering a range of systems that can
be modeled with Simbody: mechanical systems, biomechanical systems, and
molecular systems. These were collected rather than built specifically
for their use here, so you will find a variety of styles.

Precompiled binaries of these programs are installed in the examples/bin
directory. In particular, you can run SimbodyInstallTest to check installation,
or SimbodyInstallTestNoViz to check that everything other than the Visualizer
is working. You should be able to re-create these binaries using the source
provided here and instructions below.

We have supplied an example Makefile that will work on Linux and Mac, and
a CMakeLists.txt file that will work on all platforms. Although we have not
provided a pre-built Visual Studio solution file, CMake will quickly 
create one for you from the CMakeLists.txt and we strongly recommend that
approach. If you don't have it already, go to cmake.org and download 
CMake 2.8 or later. Then run the CMake GUI, tell it the source directory 
(where it will look for CMakeLists.txt) hit "Configure", tell it what 
compiler you want to use, then hit "Generate" and exit. You will then 
have a nice Visual Studio solution (.sln) file you can click on. This 
should also work for building an XCode project on the Mac although we 
haven't tested that yet. And CMake also creates ordinary make files like 
the one supplied here.

The CMakeLists.txt file here creates a single build environment containing 
separate projects for each of the examples. That works on all platforms, 
using whatever build system is appropriate.

If you want to make your own Visual Studio solution by hand, see 
README_VisualStudio.txt.

If you want to use the supplied Makefile: see README_Makefile.txt.


General Instructions
--------------------
Download CMake 2.8 or later from http://www.cmake.org. Be sure to tell 
CMake to use the version of the compiler that you intend to use. Note 
that you must use a version that is compatible with the installed 
Simbody binaries. If you can't find a suitable pre-built binary for your
compiler you'll have to build the Simbody libraries from the source
distribution. That's not very hard -- see the step-by-step instructions.

There are platform dependencies regarding running CMake; basically either
click on the icon (Windows) or run the ccmake UI. Make sure it understands
where Simbody was installed (it should see the SimTK_INSTALL_DIR environment
variable).

For help, post to the Simbody project help Forum, at 
https://simtk.org/home/simbody, Advanced tab, Public Forums.


Running the Examples
--------------------
These programs were created by different people and operate differently.
In particular, some will just start running and others want some input
before they do anything. Of the input-wanting ones, some ask for input
at the console (terminal) window, while others expect something from
the Visualizer window. Try one or the other if nothing seems to be
happening.
