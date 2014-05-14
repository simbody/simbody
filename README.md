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

#### Build and Install

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
6. Test your installation by navigating to `C:/simbody/examples/bin` and running `SimbodyInstallTest.exe` or `SimbodyInstallTestNoViz.exe`.

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


UNIX and Makefiles
------------------

The Xcode developer package gives LAPACK and BLAS to you via the Accelerate
framework. Mac's come with the visualization dependencies.

Building and Installing
-----------------------
We currently do not provide a pre-built binary distribution; you must build the
source code yourself. We're working on getting Simbody into the Debian and
Ubuntu repositories, though, so you'll be able to `apt-get` them.


#### Install CMake

[Download](http://www.cmake.org/cmake/resources/software.html) and install
CMake, which is a program we use to manage the build process. On Ubuntu, you
can obtain CMake via `$ apt-get install cmake` in a terminal.

### Configure your build

Open CMake.

#### Build and install Simbody

Say you've placed the source code into `~/simbody-source/` and that you want to
install Simbody to `<install-dir>` (perhaps, `~/simbody-install`).

```sh
$ cd ~/
$ mkdir simbody-build
$ cd simbody-build
$ cmake -DCMAKE_INSTALL_PREFIX=<install-dir> ../simbody-source
$ make install
```

You may need to run the last line as a superuser (`sudo apt-get ...`), if you
chose an `<install-dir>` that you can't otherwise write to.

If you have the dependencies, you should have no issues with the installation.

#### Configure your system to find Simbody

Set the environment variable SimTK_INSTALL_DIR `<install-dir>`. This is not
strictly required but helps with examples and locating the Visualizer, and we
use it in the next step.

* Mac: `$ echo 'export SimTK_INSTALL_DIR=<install-dir>' > ~/.bash_profile`
* Ubuntu: `$ echo 'export SimTK_INSTALL_DIR=<install-dir>' > ~/.bashrc`


Set the appropriate environment variable so the libraries can be found:

* Mac: `echo 'export DYLD_LIBRARY_PATH:$SimTK_INSTALL_DIR/lib' > `/.bash_profile`
* Ubuntu: `echo 'export
  LD_LIBRARY_PATH:$SimTK_INSTALL_DIR/lib/x86_64-linux-gnu' > ~/.bashrc`.
  Actually, you may need to replace `x86_64-linux-gnu` with the approriate
  directory on your platform.

Close and open a new terminal window.

#### Test your installation

In your new terminal window, run:

* Mac: `$ $SimTK_INSTALL_DIR/examples/bin/SimbodyInstallTest`
* Linux: `$ $SimTK_INSTALL_DIR/lib/x86_64-linux-gnu/simbody/examples/SimbodyInstallTest`



### Linux/Ubuntu
Your system's package manager surely has the dependencies. On Ubuntu you can
enter the following in a terminal:
```sh
$ apt-get install liblapack-dev
```

Optionally, for visualization:

```sh
$ apt-get install freeglut3-dev libxi-dev libxmu-dev
```
You may need to run these lines as a superuser (`sudo apt-get ...`).


[buildstatus_image]: https://travis-ci.org/simbody/simbody.png?branch=master
[travisci]: https://travis-ci.org/simbody/simbody
[user]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAndMolmodelUserGuide.pdf
[rna]:https://raw2.github.com/simbody/simbody/master/doc/images/simbios_11000_body_RNA.gif
[simbios]: http://simbios.stanford.edu/
[thy]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf
[flores]: https://simtk.org/forums/memberlist.php?mode=viewprofile&u=482
[buildwin]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_Windows.pdf
[buildunix]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_MacLinux.pdf
