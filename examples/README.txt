                     Simbody Examples README

This is an eclectic set of examples covering a range of systems that can
be modeled with Simbody: mechanical systems, biomechanical systems, and
molecular systems. These were collected rather than built specifically
for their use here, so you will find a variety of styles.

Precompiled binaries of these programs are installed in the appropriate
bin directory. In particular, you can run SimbodyInstallTest to check
installation, or SimbodyInstallTestNoViz to check that everything other
than the Visualizer is working. You should be able to re-create these
binaries using the source provided here and instructions below.

We have supplied a CMakeLists.txt file that will work on all platforms.
Although we have not provided a pre-built Visual Studio solution file,
CMake will quickly create one for you from the CMakeLists.txt and we
strongly recommend that approach. If you don't have it already, go to
cmake.org and download the current version of CMake. Then run the CMake
GUI, tell it the source directory (where it will look for CMakeLists.txt)
hit "Configure", tell it what compiler you want to use, then hit
"Generate" and exit. You will then have a nice Visual Studio solution
(.sln) file you can click on. This should also work for building an Xcode
project on Mac/OSX although that has not been well tested. And CMake also
creates ordinary make files for either Linux or OSX.

The CMakeLists.txt file here creates a single build environment containing
separate projects for each of the examples. That works on all platforms,
using whatever build system is appropriate.

General Instructions
--------------------
Download CMake from http://www.cmake.org. Be sure to tell
CMake to use the version of the compiler that you intend to use. Note
that you must use a version that is compatible with the installed
Simbody binaries, which should have been built with the same compiler.

There are platform dependencies regarding running CMake; basically either
click on the icon (Windows) or run the ccmake UI. Make sure it understands
where Simbody was installed (if not in the obvious place, then use the
SIMBODY_HOME environment variable to give it a hint).

For help, post to the Simbody project help Forum, at
https://simtk.org/home/simbody, Public Forums. Simbody can be obtained
from GitHub at https://github.com/simbody.

Running the Examples
--------------------
These programs were created by different people and operate differently.
In particular, some will just start running and others want some input
before they do anything. Of the input-wanting ones, some ask for input
at the console (terminal) window, while others expect something from
the Visualizer window. Try one or the other if nothing seems to be
happening. If you continue to have problems, post to the forum.
