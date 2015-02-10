Simbody [![Travis][buildstatus_image_travis]][travisci] [![Appveyor][buildstatus_image_appveyor]][appveyorci]
=======

[**WARNING: This is the Simbody 4.0 development branch -- work in progress!**]
=======

Simbody is a high-performance, open-source toolkit for science- and
engineering-quality simulation of articulated mechanisms, including
biomechanical structures such as human and animal skeletons,
mechanical systems like robots, vehicles, and machines, and anything
else that can be described as a set of rigid bodies interconnected
by joints, influenced by forces and motions, and restricted by
constraints. Simbody includes a multibody dynamics library for
modeling motion in [generalized/internal coordinates in O(n) time][thy].
This is sometimes called a Featherstone-style physics engine.

Simbody provides a C++ API that is used to build domain-specific applications;
it is not a standalone application itself. For example, it is used by
biomechanists in [OpenSim](http://opensim.stanford.edu), by roboticists in
[Gazebo](http://gazebosim.org), and for biomolecular research in
[MacroMoleculeBuilder (MMB)](https://simtk.org/home/rnatoolbox). Here's an
artful simulation of several RNA molecules containing thousands of bodies,
performed with MMB by [Samuel Flores][flores]:

[![Sam Flores' Simbody RNA simulation][rna]][simbios]

Read more about Simbody at the [Simbody homepage](https://simtk.org/home/simbody).


Simple example: a double pendulum
---------------------------------
Here's some code to simulate and visualize a 2-link chain:

```cpp
#include "Simbody.h"
using namespace SimTK;
int main() {
    // Define the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, 9.8);

    // Describe mass and visualization properties for a generic body.
    Body::Rigid bodyInfo(MassProperties(1.0, Vec3(0), UnitInertia(1)));
    bodyInfo.addDecoration(Transform(), DecorativeSphere(0.1));

    // Create the moving (mobilized) bodies of the pendulum.
    MobilizedBody::Pin pendulum1(matter.Ground(), Transform(Vec3(0)),
            bodyInfo, Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum2(pendulum1, Transform(Vec3(0)),
            bodyInfo, Transform(Vec3(0, 1, 0)));

    // Set up visualization.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));

    // Initialize the system and state.
    State state = system.realizeTopology();
    pendulum2.setRate(state, 5.0);

    // Simulate for 20 seconds.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(20.0);
}
```

![Double-pendulum simulation in Simbody][doublePendulum]

See [Simbody's User Guide][user] for a step-by-step explanation of this
example.


Features
--------
- Wide variety of joint, constraint, and force types; easily user-extended.
- Forward, inverse, and mixed dynamics. Motion driven by forces or
  prescribed motion.
- Contact (Hertz, Hunt and Crossley models).
- Gradient descent, interior point, and global (CMA) optimizers.
- A variety of numerical integrators with error control.
- Visualizer, using OpenGL


You want to...
--------------
* **[install Simbody](#installing)**.
* [use Simbody in your own program][user].
* [view API documentation](https://simtk.org/api_docs/simbody/latest/index.html).
* [learn the theory behind Simbody](https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf).
* [extend Simbody](https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAdvancedProgrammingGuide.pdf).
* [**get support** at the Simbody Forum](https://simtk.org/forums/viewforum.php?f=47).
* [report a bug or suggest a feature](https://github.com/simbody/simbody/issues/new).

---


Dependencies
------------

Simbody depends on the following:

* cross-platform building: [CMake](http://www.cmake.org/cmake/resources/software.html) 2.8.8 or later
* compiler: [Visual Studio](http://www.visualstudio.com) 2013 or later (Windows only), [gcc](http://gcc.gnu.org/) 4.8.1 or later (typically on Linux), or [Clang](http://clang.llvm.org/) 3.4 or later (typically on Mac, possibly through Xcode)
* linear algebra: [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/)
* visualization (optional): [FreeGLUT](http://freeglut.sourceforge.net/), [Xi and Xmu](http://www.x.org/wiki/)
* API documentation (optional): [Doxygen](http://www.stack.nl/~dimitri/doxygen/) 1.8.6 or later; we recommend at least 1.8.8.


Using Simbody
-------------

* **Creating your own Simbody-using project with CMake** To get started with
  your own Simbody-using project, check out the
  [cmake/SampleCMakeLists.txt](cmake/SampleCMakeLists.txt) file.


Installing
----------

Simbody works on Windows, Mac, and Linux. For Windows, you must build from source. For Mac and Linux, you can use a package manager or build from source. In this file, we provide instructions for 4 different ways of installing Simbody:

1. [**Windows**](#windows-using-visual-studio): build from source using Microsoft Visual Studio.
3. [**Linux or Mac (make)**](#linux-or-mac-using-make): build from source using gcc or Clang with make.
2. [**Mac (Homebrew)**](#mac-and-homebrew): automated build/install with Homebrew.
4. [**Ubuntu/Debian**](#ubuntu-and-apt-get): install pre-built binaries with apt-get.

These are not the only ways to install Simbody, however. For example, on a Mac, you could use CMake and Xcode.

#### C++11 and gcc/Clang

Simbody 3.6 and later uses C++11 features (the `-std=c++11` flag). Simbody 3.3
and earlier use only C++03 features, and Simbody 3.4 and 3.5 can use either
C++03 or C++11; see the `SIMBODY_STANDARD_11` CMake variable in these versions.
Note that if you want to use the C++11 flag in your own project, Simbody must
have been compiled with the C++11 flag as well.


Windows using Visual Studio
---------------------------

#### Get the dependencies

All needed library dependencies are provided with the Simbody installation on Windows, including linear algebra and visualization dependencies. 

1. Download and install [Microsoft Visual Studio](http://www.visualstudio.com), version 2013 or higher. The "Community Edition" is free for "non-enterprise" use. The "Express" edition is another free option, in that case use *Visual Studio Express for Windows Desktop*.
2. Download and install [CMake](http://www.cmake.org/download), version 2.8.8 or higher.
3. (optional) If you want to build API documentation, download and install Doxygen, version 1.8.8 or higher.

#### Download the Simbody source code

* Method 1: Download the source code from https://github.com/simbody/simbody/releases. Look for the highest-numbered release, click on the .zip button, and unzip it on your computer. We'll assume you unzipped the source code into `C:/Simbody-source`.
* Method 2: Clone the git repository.
    1. Get git. There are many options: [Git for Windows](http://msysgit.github.io/) (most advanced), [TortoiseGit](https://code.google.com/p/tortoisegit/wiki/Download) (intermediate; good for TortoiseSVN users), [GitHub for Windows](https://windows.github.com/) (easiest).
    2. Clone the github repository into `C:/Simbody-source`. Run the following in a Git Bash / Git Shell, or find a way to run the equivalent commands in a GUI client:

            $ git clone https://github.com/simbody/simbody.git C:/Simbody-source
            $ git checkout Simbody-3.5.1

    3. In the last line above, we assumed you want to build a released version. Feel free to change the version you want to build. If you want to build the latest development version ("bleeding edge") of Simbody off the master branch, you can omit the `checkout` line.

#### Configure and generate project files

1. Open CMake.
2. In the field **Where is the source code**, specify `C:/Simbody-source`.
3. In the field **Where to build the binaries**, specify something like `C:/Simbody-build`, just not inside your source directory. This is *not* where we will install Simbody; see below.
4. Click the **Configure** button.
    1. Choose a "generator" that corresponds to the Visual Studio you're using. For *Visual Studio 2013*, select **Visual Studio 12**. To build as 64-bit (recommended), select an option that ends with **Win64**.
    2. Click **Finish**.
5. Where do you want to install Simbody on your computer? Set this by changing the `CMAKE_INSTALL_PREFIX` variable. We'll assume you set it to `C:/Simbody`. If you choose a different installation location, make sure to use *yours* where we use `C:/Simbody` below.
6. Play around with the other build options:
    * `BUILD_EXAMPLES` to see what Simbody can do. On by default.
    * `BUILD_TESTING` to ensure your Simbody works correctly. On by default.
    * `BUILD_VISUALIZER` to be able to watch your system move about! If building remotely, you could turn this off. On by default.
    * `BUILD_STATIC_LIBRARIES` builds the three libraries as static libraries, whose names will end with `_static`. Off by default.
    * `BUILD_TESTS_AND_EXAMPLES_STATIC` if static libraries, and tests or examples are being built, creates statically-linked tests/examples. Can take a while to build, and it is unlikely you'll use the statically-linked libraries.
    * `BUILD_TESTS_AND_EXAMPLES_SHARED` if tests or examples are being built, creates dynamically-linked tests/examples. Unless you know what you're doing, leave this one on.
7. Click the **Configure** button again. Then, click **Generate** to make Visual Studio project files.

#### Build and install

1. Open `C:/Simbody-build/Simbody.sln` in Visual Studio.
2. Select your desired *Solution configuration* from the drop-down at the top.
    * **Debug**: debugger symbols; no optimizations (more than 10x slower). Library and visualizer names end with `_d`.
    * **RelWithDebInfo**: debugger symbols; optimized. This is the configuration we recommend.
    * **Release**: no debugger symbols; optimized. Generated libraries and executables are smaller but not faster than RelWithDebInfo.
    * **MinSizeRel**: minimum size; optimized. May be slower than RelWithDebInfo or Release.

    You at least want optimized libraries (all configurations but Debug are optimized), but you
    can have Debug libraries coexist with them. To do this, go through the full
    installation process twice, once for each configuration.
3. Build the project **ALL_BUILD** by right-clicking it and selecting **Build**.
4. Run the tests by right-clicking **RUN_TESTS** and selecting **Build**. Make sure all tests pass. You can use **RUN_TESTS_PARALLEL** for a faster test run if you have multiple cores.
5. (Optional) Build the project **doxygen** to get API documentation generated from your Simbody source. You will get some warnings if your doxygen version is earlier than Doxygen 1.8.8; upgrade if you can.
6. Install Simbody by right-clicking **INSTALL** and selecting **Build**.

#### Play around with examples

Within your build in Visual Studio (not the installation):

1. Make sure your configuration is set to a release configuration (e.g., RelWithDebInfo).
2. Right click on one of the targets whose name begins with `Example -` and select **Select as Startup Project**.
3. Type **Ctrl-F5** to start the program.

#### Set environment variables and test the installation

If you are only building Simbody to use it with OpenSim, you can skip this section.

1. Allow executables to find Simbody libraries (.dll's) by adding the Simbody `bin/` directory to your `PATH` environment variable.
    1. In the Start menu (Windows 7) or screen (Windows 8), search `environment`.
    2. Select **Edit the system environment variables**.
    3. Click **Environment Variables...**.
    4. Under **System variables**, click **Path**, then click **Edit**.
    5. Add `C:/Simbody/bin;` to the front of the text field. Don't forget the semicolon!
2. Allow Simbody and other projects (e.g., OpenSim) to find Simbody. In the same Environment Variables window:
    1. Under **User variables for...**, click **New...**.
    2. For **Variable name**, type `SIMBODY_HOME`.
    3. For **Variable value**, type `C:/Simbody`.
3. Changes only take effect in newly-opened windows. Close any Windows Explorer or Command Prompt windows.
4. Test your installation by navigating to `C:/Simbody/examples/bin` and running `SimbodyInstallTest.exe` or `SimbodyInstallTestNoViz.exe`.

**Note**: Example binaries are *not* installed for Debug configurations. They are present in the build environment, however, so you can run them from there. They will run *very* slowly!

#### Layout of installation

How is your Simbody installation organized?

* `bin/` the visualizer and shared libraries (.dll's,  used at runtime).
* `doc/` a few manuals, as well as API docs (`SimbodyAPI.html`).
* `examples/`
    * `src/` the source code for the examples.
    * `bin/` the examples, compiled into executables; run them! (Not installed for Debug builds.)
* `include/` the header (.h) files; necessary for projects that use Simbody.
* `lib/` "import" libraries, used during linking.
* `cmake/` CMake files that are useful for projects that use Simbody.


Linux or Mac using make
-----------------------

These instructions are for building Simbody from source on either a Mac or on
Ubuntu.

#### Get dependencies

On a Mac, the Xcode developer package gives LAPACK and BLAS to you via the Accelerate
framework. Mac's come with the visualization dependencies.

On Ubuntu, we need to get the dependencies ourselves. Open a terminal and run the following commands.

1. Get the necessary dependencies: `$ sudo apt-get install cmake liblapack-dev`. The cmake on Ubuntu 12.04 is not new enough; you could instead download it from [cmake.org](http://www.cmake.org/download/) or use [this third party PPA](https://launchpad.net/~robotology/+archive/ubuntu/ppa).
2. If you want to use the CMake GUI, install `cmake-qt-gui`.
3. For visualization (optional): `$ sudo apt-get install freeglut3-dev libxi-dev libxmu-dev`
4. For API documentation (optional): `$ sudo apt-get install doxygen`

#### Get the Simbody source code

There are two ways to get the source code.

* Method 1: Download the source code from https://github.com/simbody/simbody/releases. Look for the highest-numbered release, click on the .zip button, and unzip it on your computer. We'll assume you unzipped the source code into `~/simbody-source`.
* Method 2: Clone the git repository.
    1. Get git.
        * Mac: You might have it already, especially if you have Xcode, which is free in the App Store. If not, one method is to install [Homebrew](http://brew.sh/) and run `brew install git` in a terminal.
        * Ubuntu: run `sudo apt-get install git` in a terminal.
    2. Clone the github repository into `~/simbody-source`.

            $ git clone https://github.com/simbody/simbody.git ~/simbody-source
            $ git checkout Simbody-3.5.1

    3. In the last line above, we assumed you want to build a released version. Feel free to change the version you want to build. If you want to build the latest development version ("bleeding edge") of Simbody off the master branch, you can omit the `checkout` line.

#### Configure and generate Makefiles

1. Create a directory in which we'll build Simbody. We'll assume you choose `~/simbody-build`. Don't choose a location inside `~/simbody-source`.

        $ mkdir ~/simbody-build
        $ cd ~/simbody-build

2. Configure your Simbody build with CMake. We'll use the `cmake` command but you could also use the interactive tools `ccmake` or `cmake-gui`. You have a few configuration options to play with here.

    * If you don't want to fuss with any options, run:

            $ cmake ~/simbody-source

    * Where do you want to install Simbody? By default, it is installed to `/usr/local/`. That's a great default option, especially if you think you'll only use one version of Simbody at a time. You can change this via the `CMAKE_INSTALL_PREFIX` variable. Let's choose `~/simbody`:

            $ cmake ~/simbody-source -DCMAKE_INSTALL_PREFIX=~/simbody

    * Do you want the libraries to be optimized for speed, or to contain debugger symbols? You can change this via the `CMAKE_BUILD_TYPE` variable. There are 4 options:
        - **Debug**: debugger symbols; no optimizations (more than 10x slower). Library and visualizer names end with `_d`.
        - **RelWithDebInfo**: debugger symbols; optimized. This is the configuration we recommend.
        - **Release**: no debugger symbols; optimized. Generated libraries and executables are smaller but not faster than RelWithDebInfo.
        - **MinSizeRel**: minimum size; optimized. May be slower than RelWithDebInfo or Release.

        You at least want optimized libraries (all configurations but Debug are optimized),
        but you can have Debug libraries coexist with them. To do this, go through
        the full installation process twice, once for each configuration. It is
        typical to use a different build directory for each build type (e.g.,
        `~/simbody-build-debug` and `~/simbody-build-release`).

    * There are a few other variables you might want to play with:
        * `BUILD_EXAMPLES` to see what Simbody can do. On by default.
        * `BUILD_TESTING` to ensure your Simbody works
          correctly. On by default.
        * `BUILD_VISUALIZER` to be able to watch your system
          move about! If building on a cluster, you could turn this off. On by
          default.
        * `BUILD_STATIC_LIBRARIES` builds the three libraries as static libraries, whose names will end with `_static`.
        * `BUILD_TESTS_AND_EXAMPLES_STATIC` if tests or examples are being built, creates statically-linked tests/examples. Can take a while to build, and it is unlikely you'll use the statically-linked libraries.
        * `BUILD_TESTS_AND_EXAMPLES_SHARED` if tests or examples are being built, creates dynamically-linked tests/examples. Unless you know what you're doing, leave this one on.

        You can combine all these options. Here's another example:

            $ cmake ~/simbody-source -DCMAKE_INSTALL_PREFIX=~/simbody -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_VISUALIZER=off

#### Build and install

1. Build the API documentation. This is optional, and you can only do this if
   you have Doxygen. You will get warnings if your doxygen installation is a version older than Doxygen 1.8.8.

        $ make doxygen

2. Compile. Use the `-jn` flag to build using `n` processor cores. For example:

        $ make -j8

3. Run the tests.

        $ ctest -j8

4. Install. If you chose `CMAKE_INSTALL_PREFIX` to be a location which requires sudo access to write to (like `/usr/local/`, prepend this command with a `sudo `.

        $ make -j8 install

Just so you know, you can also uninstall (delete all files that CMake placed into `CMAKE_INSTALL_PREFIX`) if you're in `~/simbody-build`.

    $ make uninstall


#### Play around with examples

From your build directory, you can run Simbody's example programs. For instance, try:

        $ ./ExamplePendulum
        

#### Set environment variables and test the installation

If you are only building Simbody to use it with OpenSim, you can skip this section.

1. Allow executables to find Simbody libraries (.dylib's or so's) by adding the Simbody lib directory to your linker path.
    * If your `CMAKE_INSTALL_PREFIX` is `/usr/local/`, run:

            $ sudo ldconfig

    * If your `CMAKE_INSTALL_PREFIX` is neither `/usr/` nor `/usr/local/` (e.g., `~/simbody`'):
        * Mac:

                $ echo 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/simbody/lib' >> ~/.bash_profile
        * Ubuntu:

                $ echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/simbody/lib/x86_64-linux-gnu' >> ~/.bashrc
        These commands add a line to a configuration file that is loaded every time you open a new terminal. If using Ubuntu, you may need to replace `x86_64-linux-gnu` with the appropriate directory on your computer.
2. Allow Simbody and other projects (e.g., OpenSim) to find Simbody. Make sure to replace `~/simbody` with your `CMAKE_INSTALL_PREFIX`.
    * Mac:

            $ echo 'export SIMBODY_HOME=~/simbody' >> ~/.bash_profile
    * Ubuntu:
            
            $ echo 'export SIMBODY_HOME=~/simbody' >> ~/.bashrc
3. Open a new terminal.
4. Test your installation:

        $ cd ~/simbody/share/doc/simbody/examples/bin
        $ ./SimbodyInstallTest # or ./SimbodyInstallTestNoViz

#### Layout of installation

The installation creates the following directories in `CMAKE_INSTALL_PREFIX`. The directory `[x86_64-linux-gnu]` only exists if you did NOT install to `/usr/local/` and varies by platform. Even in that case, the name of your directory may be different.

* `include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `lib/[x86_64-linux-gnu]/` shared libraries (.dylib's or .so's).
    * `cmake/simbody/` CMake files that are useful for projects that use Simbody.
    * `pkgconfig/` pkg-config files useful for projects that use Simbody.
    * `simbody/examples/` the examples, compiled into executables; run them! (Not installed for Debug builds.)
* `libexec/simbody/` the `simbody-visualizer` executable.
* `share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).
    * `examples/src` source code for the examples.
    * `examples/bin` symbolic link to the runnable examples.


Mac and Homebrew
----------------

If using a Mac and Homebrew, the dependencies are taken care of for you.

#### Install

1. Install [Homebrew](http://brew.sh/).
2. Open a terminal.
3. Add the Open Source Robotics Foundation's list of repositories to Homebrew:
    ```
    $ brew tap osrf/simulation
    ```

2. Install the latest release of Simbody.
    ```
    $ brew install simbody
    ```
    To install from the master branch instead, append ` --HEAD` to the command above.

#### Where is Simbody installed?

Simbody is now installed to `/usr/local/Cellar/simbody/<version>/`, where `<version>` is either the version number (e.g., `3.5.1`), or `HEAD` if you specified `--HEAD` above.

Some directories are symlinked (symbolically linked) to `/usr/local/`, which is where your system typically expects to find executables, shared libraries (.dylib's), headers (.h's), etc. The following directories from the Simbody installation are symlinked:

* `include/simbody   -> /usr/local/include/simbody`
* `lib               -> /usr/local/lib`
* `share/doc/simbody -> /usr/local/share/doc/simbody`

#### Layout of installation

What's in the `/usr/local/Cellar/simbody/<version>/` directory?

* `include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `lib/` shared libraries (.dylib's), used at runtime.
    * `cmake/simbody/` CMake files that are useful for projects that use Simbody.
    * `pkgconfig/` pkg-config files useful for projects that use Simbody.
    * `simbody/examples/` the examples, compiled into executables; run them! (Not installed for Debug builds.)
* `libexec/simbody/` the `simbody-visualizer` executable.
* `share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).
    * `examples/src` source code for the examples.
    * `examples/bin` symbolic link to executable examples.

Ubuntu and apt-get
------------------

You can currently get Simbody via the Open Source Robotics Foundation's Debian repositories. We are currently working on getting Simbody directly into the Debian repositories. `apt-get` will take care of getting the necessary dependencies.

**Caution**: this installation method is still a work in progress. If you try it, please let us know on the [Simbody Forum](https://simtk.org/forums/viewforum.php?f=47) if it worked or if not, what problems you encountered.

#### Install

1. Setup your computer to accept software from packages.osrfoundation.org. For more detailed instructions, see [OSRF's installation instructions](http://gazebosim.org/tutorials?tut=install_ubuntu&ver=4.0&cat=install).
        
        $ sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu `lsb_release -cs` main" > /etc/apt/sources.list.d/gazebo-latest.list'
        $ wget http://packages.osrfoundation.org/gazebo.key -O - | sudo apt-key add -
        $ sudo apt-get update

2. Install Simbody.

        $ sudo apt-get update
        $ sudo apt-get install libsimbody-dev libsimbody-doc

#### Layout of installation

Simbody is installed into the `usr/` directory.  The directory
`[x86_64-linux-gnu]` varies by platform. 

* `usr/include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `usr/lib/[x86_64-linux-gnu]` shared libraries (.so's).
    * `cmake/simbody/` CMake files that are useful for projects that use Simbody.
    * `pkgconfig/` pkg-config files useful for projects that use Simbody.
* `usr/libexec/simbody/` the `simbody-visualizer` executable.
* `usr/share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).
    * `examples/src` source code for the examples.
    * `examples/bin` symbolic link to executable examples.


Acknowledgments
---------------
We are grateful for past and continuing support for Simbody's development in Stanford's Bioengineering department through the following grants:

- NIH U54 GM072970 (Simulation of Biological Structures)
- NIH U54 EB020405 (Mobilize Center)
- NIH R24 HD065690 (Simulation in Rehabilitation Research)
- OSRF subcontract 12-006 to DARPA HR0011-12-C-0111 (Robotics Challenge)

Prof. Scott Delp is the Principal Investigator on these grants and Simbody is used extensively in Scott's [Neuromuscular Biomechanics Lab](https://nmbl.stanford.edu) as the basis for the [OpenSim](http://opensim.stanford.edu) biomechanical simulation software application for medical research.



[buildstatus_image_travis]: https://travis-ci.org/simbody/simbody.png?branch=master
[travisci]: https://travis-ci.org/simbody/simbody
[buildstatus_image_appveyor]: https://ci.appveyor.com/api/projects/status/2dua0qna2m85fts2/branch/master?svg=true
[appveyorci]: https://ci.appveyor.com/project/opensim-org/simbody/branch/master
[user]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAndMolmodelUserGuide.pdf
[rna]: doc/images/simbios_11000_body_RNA.gif
[simbios]: http://simbios.stanford.edu/
[doublePendulum]: doc/images/doublePendulum.gif
[thy]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf
[flores]: http://xray.bmc.uu.se/flores/Home.html
[buildwin]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_Windows.pdf
[buildunix]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_MacLinux.pdf
