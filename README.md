Simbody [![Build Status][buildstatus_image]][travisci]
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

    // Simulate for 50 seconds.
    RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(50.0);
}
```

See [Simbody's User Guide][user] for a step-by-step explanation of this
example.


Features
--------
- Wide variety of joint, constraint, and force types; easily user-extended.
- Forward, inverse, and mixed dynamics. Motion driven by forces or
  prescribed motion.
- Contact (Hertz, Hunt and Crossley models).
- Gradient descent and interior point optimizers.
- A variety of numerical integrators with error control.
- Visualizer, using [OpenGL](http://www.opengl.org/).


You want to...
--------------
* **[install Simbody](#installing)** (see also: old [Windows][buildwin] and [Mac/Linux][buildunix] build instructions, [old install instructions][user]).
* [use Simbody in your own program][user].
* [view API documentation](https://simtk.org/api_docs/simbody/api_docs33/Simbody/html/index.html).
* [learn the theory behind Simbody](https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf).
* [extend Simbody](https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAdvancedProgrammingGuide.pdf).
* [**get support** at the Simbody forum](https://simtk.org/forums/viewforum.php?f=47).
* [report a bug or suggest a feature](https://github.com/simbody/simbody/issues/new).

---


Dependencies
------------

Simbody depends on the following:

* cross-platform building: [CMake](http://www.cmake.org/cmake/resources/software.html) 2.8 or greater
* compiler: [Visual Studio](http://www.visualstudio.com) 2010 or 2013 (Windows only), [gcc](http://gcc.gnu.org/) (typically on Linux), or [Clang](http://clang.llvm.org/) (typically on Mac)
* linear algebra: [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/)
* visualization (optional): [FreeGLUT](http://freeglut.sourceforge.net/), [Xi and Xmu](http://www.x.org/wiki/)
* API documentation (optional): [Doxygen](http://www.stack.nl/~dimitri/doxygen/) 1.7.2 or later


Installing
----------

Simbody works on Windows, Mac, and Linux. For Windows, you must build from source. For Mac and Linux, you can use a package manager or build from source. In this file, we provide instructions for 4 different ways of installing Simbody:

1. [**Windows**](#windows-and-visual-studio): build from source using Microsoft Visual Studio
2. [**Mac**](#mac-and-homebrew): install with Homebrew
3. [**Ubuntu**](#ubuntu-and-apt-get): install with apt-get
4. [**UNIX (Mac, Linux)**](#unix-and-makefiles): build from source using gcc or Clang with Makefile's

These are not the only ways to install Simbody, however. For example, on a Mac, you could use CMake and Xcode.


Windows and Visual Studio
-------------------------

#### Get the dependencies

We give the linear algebra dependencies to you, and Windows comes with the visualization dependencies.

1. Download and install Microsoft Visual Studio. If using an Express (free) version, use *Visual Studio Express 2013 for Windows Desktop* or *Visual C++ 2010 Express*.
2. Download and install CMake.
3. If you want to build API documentation, download and install Doxygen as well.

#### Download the Simbody source code

* Method 1: Download the source code from https://github.com/simbody/simbody/releases. Look for the highest-numbered release, click on the .zip button, and unzip it on your computer. We'll assume you unzipped the source code into `~/simbody-source`.
* Method 2: Clone the git repository.
    1. Get git. There are many options: [Git for Windows](http://msysgit.github.io/) (most advanced), [TortoiseGit](https://code.google.com/p/tortoisegit/wiki/Download) (intermediate; good for TortoiseSVN users), [GitHub for Windows](https://windows.github.com/) (easiest).
    2. Clone the github repository into `~/simbody-source`. Run the following in a Git Bash / Git Shell, or find a way to run the equivalent commands in a GUI client:

            $ git clone https://github.com/simbody/simbody.git ~/simbody-source
            $ git checkout Simbody-3.4

    3. In the last line above, we assumed you want to build a released version. Feel free to change the version you want to build. If you want to build the latest development version ("bleeding edge") of Simbody off the master branch, you can omit the `checkout` line.

#### Configure and generate project files

1. Open CMake.
2. In the field **Where is the source code**, specify `C:/Simbody-source`.
3. In the field **Where to build the binaries**, specify something like `C:/Simbody-build`, just not inside your source directory. This is *not* where we will install Simbody; see below.
4. Click the **Configure** button.
    1. Choose a "generator" that corresponds to the Visual Studio you're using. For *Visual Studio 2013*, select **Visual Studio 12**. To build as 64-bit, select an option that ends with **Win64**.
    2. Click **Finish**.
5. Where do you want to install Simbody on your computer? Set this by changing the `CMAKE_INSTALL_PREFIX` variable. We'll assume you set it to `C:/Simbody`. If you choose a different installation location, make sure to use *yours* where we use `C:/Simbody` below.
6. Play around with the other build options:
    * `BUILD_EXAMPLES` to see what Simbody can do. On by default.
    * `BUILD_TESTING` to ensure your Simbody works correctly. The tests take a long time to build, though. If you need to build Simbody quickly, maybe turn this off. On by default.
    * `BUILD_VISUALIZER` to be able to watch your system move about! If building on a cluster, you could turn this off. On by default.
    * `BUILD_STATIC_LIBRARIES` builds the three libraries as static libraries, whose names will end with `_static`.
    * `BUILD_TESTS_AND_EXAMPLES_STATIC` if tests or examples are being built, creates statically-linked tests/examples. Can take a while to build, and it is unlikely you'll use the statically-linked libraries.
    * `BUILD_TESTS_AND_EXAMPLES_SHARED` if tests or examples are being built, creates dynamically-linked tests/examples. Unless you know what you're doing, leave this one on.
7. Click the **Configure** button again. Then, click **Generate** to make Visual Studio project files.

#### Build and install

1. Open `C:/Simbody-build/Simbody.sln` in Visual Studio.
2. Select your desired *Solution configuration* from the drop-down at the top.
    * **Debug**: debugger symbols; no optimizations (more than 10x slower). Library names end with `_d`.
    * **Release**: no debugger symbols; optimized.
    * **RelWithDebInfo**: debugger symbols; optimized. Bigger but not slower than Release; choose this if unsure.
    * **MinSizeRel**: minimum size; optimized.

    You at least want release libraries (the last 3 count as release), but you can have debug libraries coexist with them. To do this, go through the full installation process twice, once for each configuration. You should install the release configuration *last* to ensure that you use the release version of the `simbody-visualizer` instead of the slow debug version.
3. Build the project **ALL_BUILD** by right-clicking it and selecting **Build**.
4. Run the tests by right-clicking **RUN_TESTS** and selecting **Build**.
5. Install Simbody by right-clicking **INSTALL** and selecting **Build**.

#### Play around with examples

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
    6. Under **User variables for...**, click **New...**.
    7. For **Variable name**, type `SIMBODY_HOME`.
    8. For **Variable value**, type `C:/Simbody`.
3. Changes only take effect in newly-opened windows. Close any Windows Explorer or Command Prompt windows.
4. Test your installation by navigating to `C:/Simbody/examples/bin` and running `SimbodyInstallTest.exe` or `SimbodyInstallTestNoViz.exe`.

#### Layout of installation

How is your Simbody installation organized?

* `bin/` the visualizer and shared libraries (.dll's,  used at runtime).
* `doc/` a few manuals, as well as API docs (`SimbodyAPI.html`).
* `examples/`
    * `src/` the source code for the examples.
    * `bin/` the examples, compiled into executables; run them!
    * `simmath/` source code for examples of Simbody's SimTKmath library.
* `include/` the header (.h) files; necessary for projects that use Simbody.
* `lib/` "import" libraries, used during linking.
* `share/` CMake files that are useful for projects that use Simbody.


Mac and Homebrew
----------------

If using a Mac and Homebrew, the dependencies are taken care of for you.

With this method, Simbody is built without C++11 (the `-std=c++11` compiler flag). Thus, any projects you build on top of Simbody must also NOT use C++11. If you do try to use C++11, you'll run into mysterious errors. See issue #125.

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

Simbody is now installed to `/usr/local/Cellar/simbody/<version>/`, where `<version>` is either the version number (e.g., `3.4`), or `HEAD` if you specified `--HEAD` above.

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
    * `simbody/examples/` the examples, compiled into executables; run them!
* `libexec/simbody/` the `simbody-visualizer` executable.
* `share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).
    * `examples/` source code for the examples.


Ubuntu and apt-get
------------------

You can currently get Simbody via the Open Source Robotics Foundation's Debian repositories. We are currently working on getting Simbody directly into the Debian repositories. `apt-get` will take care of getting the necessary dependencies.

With this method, Simbody is built without C++11 (the `-std=c++11` compiler flag). Thus, any projects you build on top of Simbody must also NOT use C++11. If you do try to use C++11, you'll run into mysterious errors. See issue #125.

#### Install

1. Setup your computer to accept software from packages.osrfoundation.org. This step depends on your version of Ubuntu. For more detailed instructions, see [OSRF's installation instructions](http://gazebosim.org/wiki/3.0/install#Ubuntu_Debians).
    * 12.04:

            sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu precise main" > /etc/apt/sources.list.d/gazebo-latest.list'
    * 13.10:

            sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu saucy main" > /etc/apt/sources.list.d/gazebo-latest.list'

2. Install Simbody.

        $ sudo apt-get update
        $ sudo apt-get install libsimbody-dev

#### Layout of installation

Simbody is installed into the `usr/` directory.

* `usr/include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `usr/lib/` shared libraries (.so's), used at runtime.
    * `cmake/simbody/` CMake files that are useful for projects that use Simbody.
    * `pkgconfig/` pkg-config files useful for projects that use Simbody.
* `usr/libexec/simbody/` the `simbody-visualizer` executable.
* `usr/share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).


UNIX and Makefiles
------------------

These instructions are for building Simbody from source on either a Mac or on Ubuntu.

#### Get dependencies

On a Mac, the Xcode developer package gives LAPACK and BLAS to you via the Accelerate
framework. Mac's come with the visualization dependencies.

On Ubuntu, we need to get the dependencies ourselves. Open a terminal and run the following commands.

1. Get the necessary dependencies: `$ sudo apt-get install cmake liblapack-dev`
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
            $ git checkout Simbody-3.4

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

    * Do you want to use C++11? By default, Simbody assumes no. If you plan to use Simbody in a project that DOES use C++11, then you must build Simbody with C++11 as well. You can change this via the `SIMBODY_STANDARD_11` variable:

            $ cmake ~/simbody-source -DSIMBODY_STANDARD_11=on

    * Do you want the libraries to be optimized for speed, or to contain debugger symbols? You can change this via the `CMAKE_BUILD_TYPE` variable. There are 4 options:
        - **Debug**: debugger symbols; no optimizations (more than 10x slower). Library names end with `_d`.
        - **Release**: no debugger symbols; optimized.
        - **RelWithDebInfo**: debugger symbols; optimized. Bigger but not slower than Release; choose this if unsure.
        - **MinSizeRel**: minimum size; optimized.

        You at least want release libraries (the last 3 count as release), but you can have debug libraries coexist with them. To do this, go through the full installation process twice, once for each configuration. It is typical to use a different build directory for each build type (e.g., `~/simbody-build-debug` and `~/simbody-build-release`). You should install the release configuration *last* to ensure that you use the release version of the `simbody-visualizer` instead of the slow debug version.

    * There are a few other variables you might want to play with:
        * `BUILD_EXAMPLES` to see what Simbody can do. On by default.
        * `BUILD_TESTING` to ensure your Simbody works
          correctly. The tests take a long time to build, though. If you need to
          build Simbody quickly, maybe turn this off. On by default.
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
   you have Doxygen.

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

                $ sudo echo 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/simbody/lib' > /etc/profile
        * Ubuntu:

                $ sudo echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/simbody/lib/x86_64-linux-gnu' > ~/.bashrc
        These commands add a line to a configuration file that is loaded every time you open a new terminal. If using Ubuntu, you may need to replace `x86_64-linux-gnu` with the appropriate directory on your computer.
2. Allow Simbody and other projects (e.g., OpenSim) to find Simbody. Make sure to replace `~/simbody` with your `CMAKE_INSTALL_PREFIX`.
    * Mac:

            $ sudo echo 'export SIMBODY_HOME=~/simbody' > /etc/profile
    * Ubuntu:
            
            $ sudo echo 'export SIMBODY_HOME=~/simbody' > ~/.bashrc
3. Open a new terminal.
4. Test your installation:

        $ cd ~/simbody/share/doc/simbody/examples/bin
        $ ./SimbodyInstallTest # or ./SimbodyInstallTestNoViz

#### Layout of installation

The installation creates the following directories in `CMAKE_INSTALL_PREFIX`. The directory `[x86_64-linux-gnu]` only exists if you're using a recent version of Ubuntu (e.g., 13.10) and did NOT install to `/usr/local/`. Even in that case, the name of your directory may be different.

* `include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `lib/[x86_64-linux-gnu]/` shared libraries (.dylib's or .so's), used at runtime.
    * `cmake/simbody/` CMake files that are useful for projects that use Simbody.
    * `pkgconfig/` pkg-config files useful for projects that use Simbody.
    * `simbody/examples/` the examples, compiled into executables; run them!
* `libexec/simbody/` the `simbody-visualizer` executable.
* `share/doc/simbody/` a few manuals, as well as API docs (`SimbodyAPI.html`).
    * `examples/` source code for the examples.


[buildstatus_image]: https://travis-ci.org/simbody/simbody.png?branch=master
[travisci]: https://travis-ci.org/simbody/simbody
[user]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAndMolmodelUserGuide.pdf
[rna]: doc/images/simbios_11000_body_RNA.gif
[simbios]: http://simbios.stanford.edu/
[thy]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf
[flores]: https://simtk.org/forums/memberlist.php?mode=viewprofile&u=482
[buildwin]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_Windows.pdf
[buildunix]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_MacLinux.pdf
