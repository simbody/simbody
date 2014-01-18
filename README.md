Simbody [![Build Status][buildstatus_image]][travisci]
=======

Simbody is a high-performance, open-source toolkit for simulation of
articulated mechanisms, including biological structures, robots, vehicles,
machines, etc. Simbody includes a multibody dynamics library for modeling
motion in [generalized/internal coordinates in O(n) time][thy].

It is used by biomechanists in [OpenSim](http://opensim.stanford.edu) and by
roboticists in [Gazebo](http://gazebosim.org). Here's an 11,000-body simulation
of an RNA molecule performed with Simbody by [Sam Flores][flores]:

[![Sam Flores' Simbody RNA simulation][rna]][simbios]

Read more about it here: https://simtk.org/home/simbody.


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
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));

    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));

    // Create the moving bodies of the pendulum.
    MobilizedBody::Pin pendulum1(matter.Ground(), Transform(Vec3(0)),
            pendulumBody, Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum2(pendulum1, Transform(Vec3(0)),
            pendulumBody, Transform(Vec3(0, 1, 0)));

    // Set up visualization.
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));

    // Initialize the system and state.
    system.realizeTopology();
    State state = system.getDefaultState();
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
- Visualizer, using [OpenGL](http://www.opengl.org/).
- Contact (Hertz, Hunt and Crossley models).
- Gradient descent optimizers (IPOPT, LBFGS, LBFGSB).
- Path wrapping around objects, used to define forces that act along complex
  paths, like muscles do.

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


Other multibody dynamics engines
--------------------------------
* [Open Dynamics Engine](http://www.ode.org/)
* [Bullet](http://bulletphysics.org)
* [Moby](http://physsim.sourceforge.net/)


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
