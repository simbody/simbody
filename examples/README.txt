                      SimTK 2.1 Examples README
                           August, 2010

This is an eclectic set of examples covering a range of systems that can
be modeled with SimTK: mechanical systems, biomechanical systems, and
molecular systems. These were collected rather than built specifically
for their use here, so you will find a variety of styles.

We have supplied an example Visual Studio project and makefile, but
by far the easiest way to create a build project for these examples
is to use CMake (www.cmake.org). We have supplied a CMakeLists.txt file here 
which will build a single build project containing all the examples.
That works on all platforms, using whatever build system is appropriate,
and even knows how to make Visual Studio solution files.

If you want to use the supplied VisualStudio 8 (2005) or 9 (2008) solutions
or make your own Visual Studio solution: see README_VisualStudio.txt.

If you want to use the supplied Makefile: see README_Makefile.txt.


General Instructions
--------------------

Download CMake 2.6 or later from http://www.cmake.org. All the examples
were tested successfully on Windows Vista using CMake 2.8.1, but later
CMake releases *should* work. Be sure to tell CMake to use the version of
the compiler that you intend to use.

There are platform dependencies regarding running CMake; basically either
click on the icon (Windows) or run the ccmake UI. Make sure it understands
where SimTK was installed (it should see the SimTK_INSTALL_DIR environment
variable).

For help, post to the SimTK Core project help Forum, at 
https://simtk.org/home/simtkcore, Advanced tab, Public Forums.


Special Instructions
--------------------

ExampleAdenylateMobilitiesVMD

There are special instructions for running ExampleAdenylateMobilitiesVMD
because you have to install and run VMD for displaying the results. See
README_ExampleAdenylateMobilitiesVMD.txt.

ExampleLoadPdb[2]

The two ExampleLoadPdb examples expect to find their pdb files in the
current directory. If you are running them from somewhere other than
the Examples source directory, copy the pdb files there.

ExampleDynamicWalker

This program has some options so asks some questions when it starts up.
Look at the console or shell window to see the prompts and make your
choices; then it will start running.
