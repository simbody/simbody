Simbody [![Build Status][buildstatus_image]][travisci]
=======

**[Caution: this README.md intro is under development and is buggy]**
--

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

Read more about Simbody [here](https://simtk.org/home/simbody).


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
* build Simbody from the source code: [Windows][buildwin], [Mac or
  Linux][buildunix] (somewhat outdated).
* [install Simbody from a binary distribution][user] (outdated).
* [use Simbody in your own program][user].
* [learn the theory behind Simbody][thy].
* [extend Simbody][extend].
* [ask questions or share information in the Simbody forum][forum]
* [report a bug or suggest a feature][issue]
* [contribute or get the source code](https://github.com/simbody/simbody/)

Dependencies
------------
Simbody depends on LAPACK and BLAS for math, and on FreeGLUT, Xi and Xmu for
visualization.

### Windows
We give the dependencies to you!

### Mac
The XCode developer package gives LAPACK and BLAS to you via the Accelerate
framework. Mac's come with the visualization dependencies.

### Linux/Ubuntu
Your system's package manager surely has the dependencies. On Ubuntu you can
enter the following in a terminal:
```sh
$ apt-get install liblapack-dev
```

Optionally, for visualization:

```sh
$ apt-get install freeglut-dev libxi-dev libxmu-dev
```
You may need to run these lines as a superuser (`sudo apt-get ...`).


Building and Installing
-----------------------
We currently do not provide a pre-built binary distribution; you must build the
source code yourself. We're working on getting Simbody into the Debian and
Ubuntu repositories, though, so you'll be able to `apt-get` them.

#### Download the source code
Download the source code from https://github.com/simbody/simbody/releases.
Alternatively, if you have [Git](http://git-scm.com) installed, you can clone
the Simbody repository. From a command line, this would look like:

```sh
$ git clone git@github.com/simbody/simbody.git
```

From this point, if you're on Windows or you don't want to do anything on the
command line, head straight for either the [Windows][buildwin] or the [Unix
(Mac/Linux)][buildunix] instructions; we can't do much for you without showing
you screenshots.

#### Install CMake

[Download](http://www.cmake.org/cmake/resources/software.html) and install
CMake, which is a program we use to manage the build process. On Ubuntu, you
can obtain CMake via `$ apt-get install cmake` in a terminal.

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


[buildstatus_image]: https://travis-ci.org/simbody/simbody.png?branch=master
[travisci]: https://travis-ci.org/simbody/simbody
[user]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAndMolmodelUserGuide.pdf
[rna]:https://raw2.github.com/simbody/simbody/master/doc/images/simbios_11000_body_RNA.gif
[simbios]: http://simbios.stanford.edu/
[thy]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyTheoryManual.pdf
[extend]: https://github.com/simbody/simbody/raw/master/Simbody/doc/SimbodyAdvancedProgrammingGuide.pdf
[flores]: https://simtk.org/forums/memberlist.php?mode=viewprofile&u=482
[forum]: https://simtk.org/forums/viewforum.php?f=47
[issue]: https://github.com/simbody/simbody/issues/new
[buildwin]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_Windows.pdf
[buildunix]: https://github.com/simbody/simbody/raw/master/doc/HowToBuildSimbodyFromSource_MacLinux.pdf
