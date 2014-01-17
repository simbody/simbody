Simbody
=======
[![Build Status](https://travis-ci.org/simbody/simbody.png?branch=master)](https://travis-ci.org/simbody/simbody)

Source code for the Simbody libraries, examples, and build system.

Simbody is a high-performance, open-source toolkit for simulation of 
articulated mechanisms, including biological structures, robots, vehicles,
machines, etc. Simbody includes a multibody dynamics library for modeling
motion in internal coordinates in O(n) time.

Read more about it here: https://simtk.org/home/simbody


Simple example: a double pendulum
---------------------------------
Want to know the basics of what Simbody can do? Look no further! This example
comes from Simbody's User Guide.

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

Installing
----------
Detailed installation instructions for Linux, Mac, and Windows installation
of precompiled binaries are provided in the Simbody User's Guide, which can 
be found in the installation's doc directory or posted online at 
https://simtk.org/home/simbody, Documents tab. If you are building from 
source, use the How To Build From Source document for your platform instead;
you can find that in the same places. 

We highly recommend that you read the detailed installation instructions in
the User's Guide. However, if you are one of those guys who just won't ask 
for directions, here at least is an abbreviated summary of the installation 
procedure for the precompiled binaries:

1. Install prerequisites for your platform (see below)
2. Download an appropriate .zip file from Simbody's Download page.
3. Unzip it into the installation directory, or unzip it somewhere and then
   move the whole directory structure to where you want it to reside. That can
   be anywhere you want, but we suggest 
    - Linux or Mac:   /usr/local/SimTK 
    - 32-bit Windows: C:\Program Files\SimTK
    - 64-bit Windows: C:\Program Files (x86)\SimTK
4. Set the environment variable SimTK_INSTALL_DIR to the installation directory
   you chose in step 3. (Not strictly required but helps with examples and
   locating the Visualizer.)
5. Set the appropriate PATH variable so the libraries can be found:
    - Linux:   export LD_LIBRARY_PATH=$SimTK_INSTALL_DIR/lib[64]
    - Mac:     export DYLD_LIBRARY_PATH=$SimTK_INSTALL_DIR/lib
    - Windows: set    path=%SimTK_INSTALL_DIR%\bin;%path%
6. Test installation using:
    $SimTK_INSTALL_DIR/examples/bin/SimbodyInstallTest

Prerequisites:
-  Windows: None
-  Mac: Must have installed developer package (XCode), to get Lapack and Blas
       (they are part of Apple's Accelerate framework).
-  Linux:
    - Lapack and blas version 3: install package liblapack-dev. (Must end up
      with liblapack.so.3 and libblas.so.3 in /usr/lib[64].)
    - These are needed if you want to run the Simbody Visualizer (you do):
        - Freeglut: install freeglut-dev package
        - Xi, Xmu: install libxi-dev, libxmu-dev

### Getting help

1. Read the installation instructions in the User's Guide (see above).
2. Post to the help forum at https://simtk.org/home/simbody, Advanced tab,
   Public forums.


Building from source
--------------------

### Unix (Ubuntu, Mac OSX, etc.)
Assuming we start in a directory containing `simbody/`, which contains the
source code.

```
mkdir simbody-build
cd simbody-build
cmake -DCMAKE_INSTALL_PREFIX=<desired-installation-path> ../simbody-build
make install
```

You may need to run that last line as a superuser (`sudo make install`),
depending on how you choose `<desired-installation-path>`.



