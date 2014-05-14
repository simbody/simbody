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
* [view the Simbody forum](https://simtk.org/forums/viewforum.php?f=47).
* [report a bug or suggest a feature](https://github.com/simbody/simbody/issues/new).
* [get the source code, or contribute!](https://github.com/simbody/simbody/)

---


Dependencies
------------

Simbody depends on the following:

* cross-platform building: [CMake](http://www.cmake.org/cmake/resources/software.html) 2.8 or greater.
* compiler: [Visual Studio](http://www.visualstudio.com) 2010 or 2013 (Windows only), [gcc](http://gcc.gnu.org/) (typically on Linux), or [Clang](http://clang.llvm.org/) (typically on Mac)
* linear algebra: [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/)
* visualization (optional): [FreeGLUT](http://freeglut.sourceforge.net/), [Xi and Xmu](http://www.x.org/wiki/)
* API documentation (optional): [Doxygen](http://www.stack.nl/~dimitri/doxygen/)


Installing
----------

Simbody works on Windows, Mac, and Linux. For Windows, you must build from source. For Mac and Linux, you can use a package manager or build from source. In this file, we provide instructions for 4 different ways of installing Simbody:

1. [**Windows**](#windows-and-visual-studio): build from source using Microsoft Visual Studio
2. [**Mac**](#mac-and-homebrew): install with Homebrew
3. [**Ubuntu**](#ubuntu-and-apt-get): install with apt-get
4. [**UNIX (Mac, Linux)**](#unix-and-makefiles): build from source using gcc or Clang with Makefile's.

These are not the only ways to install Simbody, however. For example, on a Mac, you could use CMake and Xcode.


Windows and Visual Studio
-------------------------

#### Get the dependencies

We give the linear algebra dependencies to you, and Windows comes with the visualization dependencies.

1. Download and install Microsoft Visual Studio. If using an Express (free) version, use *Visual Studio Express 2013 for Windows Desktop* or *Visual C++ 2010 Express*.
2. Download and install CMake.
3. If you want to build API documentation, download and install Doxygen as well.

#### Download the Simbody source code

Download the source code from https://github.com/simbody/simbody/releases. Look for the highest-numbered release, click on the .zip button, and unzip it on your computer. We'll assume you unzipped the source code into `C:/simbody-source`.

#### Configure and generate project files

1. Open CMake.
2. In the field **Where is the source code**, specify `C:/simbody-source`.
3. In the field **Where to build the binaries**, specify something like `C:/simbody-build`.
4. Click the **Configure** button.
    1. Choose a "generator" that corresponds to the Visual Studio you're using. For *Visual Studio 2013*, select **Visual Studio 12**. To build as 64-bit, select an option that ends with **Win64**.
    2. Click **Finish**.
5. Where do you want to install Simbody on your computer? Set this by changing the `CMAKE_INSTALL_PREFIX` variable. We'll assume you set it to `C:/simbody`.
6. Click the **Configure** button again. Then, click **Generate** to make Visual Studio project files.

#### Build and install

1. Open `C:/simbody-build/Simbody.sln` in Visual Studio.
2. Select your desired *Solution configuration* from the drop-down at the top.
    * **Debug**: debugger symbols; no optimizations (*very* slow).
    * **Release**: no debugger symbols; optimized.
    * **RelWithDebInfo**: debugger symbols; optimized. Select this if you don't know.
    * **MinSizeRel**: minimum size; optimized.
3. Build the project **ALL_BUILD** by right-clicking it and selecting **Build**.
4. Run the tests by right-clicking **RUN_TESTS** and selecting **Build**.
5. Install Simbody by right-clicking **INSTALL** and selecting **Build**.

#### Set environment variables and test the installation

If you are only building Simbody to use it with OpenSim, you can skip this section.

1. Allow executables to find Simbody libraries (.dll's) by adding the Simbody `bin/` directory to your `PATH` environment variable.
    1. In the Start menu (Windows 7) or screen (Windows 8), search `path`.
    2. Select **Edit the system environment variables**.
    3. Click **Environment Variables...**.
    4. Under **System variables**, click **Path**, then click **Edit**.
    5. Add `C:/simbody/bin;` to the front of the text field. Don't forget the semicolon!
    6. Changes only take effect in newly-opened windows.
3. Test your installation by navigating to `C:/simbody/examples/bin` and running `SimbodyInstallTest.exe` or `SimbodyInstallTestNoViz.exe`.

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

* Setup your computer to accept software from packages.osrfoundation.org. This step depends on your version of Ubuntu. For more detailed instructions, see [OSRF's installation instructions](http://gazebosim.org/wiki/3.0/install#Ubuntu_Debians).
    * 12.04:
    
        ```
        sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu precise main" > /etc/apt/sources.list.d/gazebo-latest.list'
        ```

    * 13.10:
     
        ```
        sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu saucy main" > /etc/apt/sources.list.d/gazebo-latest.list'
        ```
        
* Install Simbody.

    ```
    $ sudo apt-get update
    $ sudo apt-get install libsimbody-dev
    ```

#### Layout of installation

Simbody is installed into the `usr/` directory.

* `usr/include/simbody/` the header (.h) files; necessary for projects that use Simbody.
* `usr/lib/` shared libraries (.dylib's), used at runtime.
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
    
        ```
        $ git clone https://github.com/simbody/simbody.git ~/simbody-source
        $ git checkout Simbody-3.4
        ```
        
    3. In the last line above, we assumed you want to build a released version. Feel free to change the version you want to build. If you want to build the latest development version ("bleeding edge") of Simbody off the master branch, you can omit the `checkout` line.

#### Configure and generate Makefiles

1. Create a directory in which we'll build Simbody.

    ```
    $ mkdir ~/simbody-build
    $ cd ~/simbody-build
    ```

2. Configure your Simbody build with CMake. We'll use the `cmake` command but you could also use the interactive tools `ccmake` or `cmake-gui`. You have a few configuration options to play with here.

    * If you don't want to fuss with any options, run:

        ```
        $ cmake ~/simbody-source
        ```
    
    * Where do you want to install Simbody? By default, it is installed to `/usr/local/`. You can change this via the `CMAKE_INSTALL_PREFIX` variable. Let's choose `~/simbody`:
    
        ```
        $ cmake ~/simbody-source -DCMAKE_INSTALL_PREFIX=~/simbody
        ```
    
    * Do you want to use C++11? By default, Simbody assumes not. If you plan to use Simbody in a project that DOES use C++11, then you must build Simbody with C++11 as well. You can change this via the `SIMBODY_STANDARD_11` flag:
    
        ```
        $ cmake ~/simbody-source -DSIMBODY_STANDARD_11=on
        ```
    
    * There are a few other variables you might want to play with:
        * `BUILD_EXAMPLES` on by default
        * `BUILD_TESTING` on by default
        * `BUILD_VISUALIZER` on by default
        
        You can combine all these options. Here's another example:
        
        ```
        $ cmake ~/simbody-source -DCMAKE_INSTALL_PREFIX=~/simbody -DBUILD_VISUALIZER=off 
        ```

#### Build and install

3. Compile. Use the `-jn` flag to build using `n` processor cores. For example:
    ```
    $ make -j8
    ```

4. Run the tests.
    ```
    $ ctest -j8
    ```

5. Install. The `sudo` is there in case your `CMAKE_INSTALL_PREFIX` is something like `/usr/local/` (the default).
    ```
    $ sudo make -j8 install
    ```

Just so you know, you can also uninstall (delete all files that CMake placed into `CMAKE_INSTALL_PREFIX`) if you're in `~/simbody-build`.
```
$ make uninstall
```

#### Set environment variables and test the installation

If you are only building Simbody to use it with OpenSim, you can skip this section.

1. Allow executables to find Simbody libraries (.dylib's or so's) by adding the Simbody lib directory to your linker path. There are two cases in which this is unnecessary:
    1. If you chose your `CMAKE_INSTALL_PREFIX` to be `/usr/`
    2. If you chose your `CMAKE_INSTALL_PREFIX` to be `/usr/local/` (the default), AND your libraries are in `/usr/local/lib/`. Go check! On some platforms, the libraries are in an additional subdirectory (on Ubuntu 13.10: `/usr/local/lib/x86_64-linux-gnu`).

        * Mac:
    
            ```
            $ sudo echo 'export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:~/simbody/lib' > /etc/profile
            ```
        
        * Ubuntu:
        
            ```
            $ sudo echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/simbody/lib/x86_64-linux-gnu' > ~/.bashrc
            ```
        
        These commands add a line to a configuration file that is loaded every time you open a new terminal. If using Ubuntu, you may need to replace `x86_64-linux-gnu` with the appropriate directory on your computer.

2. Open a new terminal.
3. Test your installation:

    ```
    $ cd ~/simbody/share/doc/simbody/examples/bin
    $ ./SimbodyInstallTest # or ./SimbodyInstallTestNoViz
    ```

[buildstatus_image]: https://travis-ci.org/simbody/simbody.png?branch=master
[travisci]: https://travis-ci.org/simbody/simbody
[user]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAndMolmodelUserGuide.pdf
[rna]:https://raw2.github.com/simbody/simbody/master/doc/images/simbios_11000_body_RNA.gif
[simbios]: http://simbios.stanford.edu/
[thy]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf
[flores]: https://simtk.org/forums/memberlist.php?mode=viewprofile&u=482
[buildwin]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_Windows.pdf
[buildunix]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_MacLinux.pdf
