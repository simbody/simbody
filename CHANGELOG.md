Simbody Changelog and Release Notes
===================================
**(work in progress)**

This is not a comprehensive list of changes but rather
a hand-curated collection of the more notable ones. For
a comprehensive history, see the [Simbody GitHub repo]
(https://github.com/simbody/simbody). You can use the release dates below to find all the PRs and issues that were 
included in a particular release. 

3.5 (in progress)
-----------------

This release focused primarily on infrastructure for and prototyping of rigid contact and impact, and the development of examples showing how to perform task space control using Simbody. These two projects were supported by our DARPA research subcontract with Open Source Robotics Foundation, and were integrated with Gazebo. Further development for rigid contact is required for smooth integration into Simbody; this is planned for Simbody 4.0 and only the bravest among you should attempt to use rigid contact before then. The task space control examples `TaskSpaceControl-UR10` and `TaskSpaceControl-Atlas` can be found in the Simbody examples directory.

Chris Dembia integrated Nikolaus Hansen's [Covariant Matrix Adaptation Evolution Strategy](https://www.lri.fr/~hansen/cmaesintro.html) (CMA-ES) global optimizer into Simbody's existing Optimizer class framework, and implemented a multithreading capability for it. This is ready to use and we would like feedback. 

There were numerous smaller improvements to Simbody since the previous release, in build and installation, documentation, performance, bug fixes, and small enhancements. There are no known incompatibilities with previous release 3.4.1 and we highly recommend that you upgrade.

### New features
* Added Task Space control examples (pr #237 #238) (work supported by OSRF)
* Added IMU (orientation) tracking assembly condition for IK (pr #239)
* Added STL file reader to PolygonalMesh class (issue #57, pr #176)
* Added CMA-ES global optimizer with multithreading support (pr #249 #253 #267 #268)
* Initial rigid contact support (pr #137). Includes a collection of new unilateral constraints and a specialized solver. Not yet well integrated into Simbody overall. (work supported by OSRF)
* Implemented PLUS method impact handling (issue #115 #205, pr #226) (work supported by OSRF)

### Bugs fixed
* Fixed bug in orientation of non-camera-facing text in visualizer (issue #214, pr #215)
* Fixed bug in mesh triangulation (issue #198, pr #199)
* Fixed Assembler bugs; could sometimes make initial solution worse or report spurious failure (issue #164 #167, pr #165 #168)
* Fixed Debug visualizer name to have "_d" (issue #5)

### Misc. improvements
* Improved regression test timing framework to make it useful for multithreaded tests
* Much nicer Doxygen for Rotation class (issue #183)
* Reorganized Simbody example programs (pr #238)
* Added HalfSpace-ConvexImplicitSurface contact tracker (issue #232, pr #233)
* Added methods for principal curvatures and directions on arbitrary implicit surface (pr #231)
* Added missing calcContactPatch() functionality for Hertz contacts (pr #224)
* Added Brick ContactGeometry shape (pr #218)
* Added several new bilateral constraints useful as the basis for unilateral constraints in rigid contact: Constraint::PointOnPlaneContact, SphereOnPlaneContact, SphereOnSphereContact,  LineOnLineContact (edge/edge) (pr #137 #169)
* Improved constraint performance for several basic constraints (pr #137)
* Moved install instructions to README.md where they can be kept up to date (pr #144 and others)
* Replaced distance constraint equations to use length rather than length^2 (issue #3). This improves scaling when distance constraint is combined with other constraints.
* Numerous improvements to build, install, documentation, and performance.

3.4.1 (31 Mar 2014)
-------------------

This was primarily a release for improving our build and install process to comply with Debian's requirements. Thanks to José Rivero and Steve Peters at Open Source Robotics Foundation, and Chris Dembia at Stanford for the bulk of this effort.

* Fixed `SimbodyMatterSubsystem::getTotalCentrifugalForces()` (issue #112, pr #116).
* Many changes to build and install, mostly affecting Linux and OSX. Should now conform better to standards on those platforms and in general be better and finding its dependencies. (pr #114 #118 #120 #121 #122 #127 #130 #131 #133, issue #48 #68 #92 #113) 
* Started using the -Werror flag (treat warnings as errors) in Travis builds to ensure that we are warning-free. (issue #128, pr #129 #131)
* Compile with C++11 enabled by default. However, the code will still build in C++03. (pr #123, #126)


3.3.1 (21 Jan 2014)
-------------------

This is the first release built for use in Open Source Robotic Foundation's Gazebo robot simulator and is also the version of Simbody that ships with OpenSim 3.2. It incorporates many fixes and enhancements prompted by the integration effort with OSRF, and a new Debian builder for smooth incorporation into the Gazebo build.




3.1 (15 Aug 2013)
-----------------

This is the Simbody release that shipped with OpenSim 3.1. The source was managed on Subversion, although its code and history were transferred to GitHub, from which a determined code archeologist could presumably reconstruct the changelog. But why? The source zip is also available on the GitHub release page [here](https://github.com/simbody/simbody/releases/tag/Simbody-3.1), and does include a changelog named `ChangeLog.txt`.