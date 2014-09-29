Instructions for using the Makefile provided.

----------------------------------------------------------------------------
NOTE: The easiest way to create a build project is to use
CMake (www.cmake.org). We have supplied a CMakeLists.txt file here which 
will build a single Makefile containing all the examples.
That works on all platforms also, using whatever build system is appropriate,
even Visual Studio.  
----------------------------------------------------------------------------


You must already have the Simbody binaries installed from SimTK.org; go 
to https://simtk.org/home/simbody, Downloads tab.

You may need to slightly edit the Makefile to make it run on your system.

Open the Makefile and make the appropriate changes -- debug libraries or not,
remove -m32 if you are running 64 bits on Linux, and set the default 
install directory to the correct location.

Type "make" to build just one simple example, SimbodyInstallTestNoViz and
SimbodyInstallTest, and try running them if they succeed: 
./SimbodyInstallTestNoViz, ./SimbodyInstallTest.

If that works, then try to build all the examples with "make all".

To compile just ExampleChain (for example) type "make ExampleChain".

Before you run the executables, remember to add the Simbody library 
directory to your library path. The simplest way to do this is to 
type the following commands:

First, if you didn't install Simbody in the default location, 
/usr/local/SimTK, tell the Makefile where to look:

    $ export SimTK_INSTALL_DIR=/mysimbody/installdir  

Then you have to tell the linker where to find the shared libraries:

    For Linux (for the bash shell):
    $ export LD_LIBRARY_PATH=$SimTK_INSTALL_DIR/lib 

    For Mac (for the bash shell):
    $ export DYLD_LIBRARY_PATH=$SimTK_INSTALL_DIR/lib 

If you didn't set SimTK_INSTALL_DIR, then you would type 
"/usr/local/SimTK/lib" instead above. On 64-bit Linux, the "lib64" directory
is used rather than "lib".
