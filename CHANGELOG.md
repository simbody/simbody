Simbody Changelog and Release Notes
===================================

This is not a comprehensive list of changes but rather a hand-curated collection
of the more notable ones. For a comprehensive history, see the
[Simbody GitHub repo](https://github.com/simbody/simbody). You can use the release dates below to find all the
PRs and issues that were included in a particular release.

3.9 (November 2026)
-------------------
This release requires C++20 and includes several changes that may require
updates to code or build scripts that use Simbody; see "Breaking changes" below.
It also includes the CableSpan improvements needed by OpenSim 4.6.

**Breaking changes**
* The required C++ level was increased from C++11 to C++20 (PRs #840 and #849).
  Code that uses Simbody must also be compiled as C++20 or later. Building
  Simbody requires at least GCC 11, Clang 15, Apple Clang (Xcode) 15, or
  Visual Studio 2019 version 16.11.
* Constants in `Scalar.h` (`NaN`, `Infinity`, `Eps`, `Pi`, etc.) are now
  `inline constexpr` (defined in the header) rather than exported library
  symbols, and the `NTraits<T>` getters are now `constexpr`. This eliminates
  static initialization order problems when linking Simbody statically. The
  predefined colors (`Black`, `Red`, ...) and the `XAxis`, `YAxis`, and `ZAxis`
  coordinate axis constants are likewise now `constexpr`. Code that relied on
  these constants being a single exported symbol (e.g., comparing addresses
  across shared-library boundaries) may be affected. (PRs #852 and #860)
* `MobilizedBody::Translation` "lone particle" support has been removed. (PR #859)
* `SimTK::Value<T>` now always value-initializes its contained `T`, so a
  default-constructed `Value<double>` (for example) holds zero rather than
  uninitialized memory. (PR #853)
* Removed the CMake options deprecated in 3.8: `BUILD_DYNAMIC_LIBRARIES`,
  `BUILD_STATIC_LIBRARIES`, `BUILD_UNVERSIONED_LIBRARIES`,
  `BUILD_TESTS_AND_EXAMPLES_SHARED` and `BUILD_TESTS_AND_EXAMPLES_STATIC`. Use
  `SIMBODY_BUILD_SHARED_LIBS` and `BUILD_VERSIONED_LIBRARIES` instead; the old
  options are now ignored. (PR #874)
* When Simbody is built against user-supplied BLAS/LAPACK libraries (e.g.,
  OpenBLAS via `BUILD_USING_OTHER_LAPACK`), the install step no longer copies
  those libraries into the Simbody installation; only the vendored BLAS/LAPACK
  are installed. (PR #861)

**CableSpan** (used for muscle and ligament paths in OpenSim 4.6)
* Added support for via points in a `CableSpan` path; see
  `CableSpan::addViaPoint()`. Via points divide the cable into segments, each
  of which is optimized independently. (PR #828)
* `CableSpan` now implements the Subsystem position and velocity realization
  methods, fixing erratic path behavior (jumping to invalid configurations)
  during time-stepping simulations. (PR #832)
* The `CableSpan` unit force calculation methods (`calcOriginUnitForce()`,
  `calcTerminationUnitForce()`, `calcCurveSegmentUnitForce()`, and
  `calcViaPointUnitForce()`) now ensure that the position-level cache is valid.
  (PR #846)
* Added `CableSpan::calcResampledDecorativePathPoints()`, which generates
  decorative points along the whole path with curve segments resampled at equal
  length intervals. (PR #855)
* Added `CableSpan` methods for calculating path tangent directions at the
  origin, termination, and via points (e.g.,
  `calcOriginTangentDirection()`), to make it easier to apply cable tension as
  a set of point forces. (PR #857)
* Added `CableSpan::getUseWarmStart()` and `setUseWarmStart()` to control
  whether the previous path solution is used as the starting point for the next
  path solve. Disabling warm starts makes the path (and moment arms) depend only
  on the current configuration and the contact point hints. (PR #858)
* `CableSpan::calcResampledDecorativePathPoints()` and
  `CableSpan::calcCurveSegmentResampledPoints()` now always emit decorative
  points whether or not a curve segment is in contact with an obstacle, so that
  the number of emitted points is constant. Previously they emitted NaN points,
  or read out of bounds, for obstacles that were not in contact. (PR #872)
* Fixed the `CableSpan` solver's per-iteration step size limit, which ignored
  the sign of the angular displacement of the curve segment boundary frames, so
  that large rotations in the negative direction were not limited. (PR #881)

**Other new features and improvements**
* Added support for variable mobilizer inboard and outboard frames. (PR #859)
* Added `MobilizedBody::CantileverFreeBeam`, a mobilizer for modeling spinal
  joints and other flexible structures. (PR #844)
* `MobilizedBody::Ellipsoid` can now be mirrored about the plane normal to one
  axis of its inboard frame F, which can be useful for building left/right joint
  pairs (e.g., shoulder joints) driven by the same coordinates. (PR #875) The
  mirror axis can also be set or cleared after construction with
  `setMirrorAxis()` and `clearMirrorAxis()`. (PR #880)
* Added `Rotation::mirror()`, a method for mirroring a Rotation about a
  coordinate axis. (PR #879)
* Added LLT (Cholesky) factorization, `FactorLLT`. (PR #833)
* Added quaternion multiplication (Hamilton product): `Quaternion::multiply()`
  and `operator*`. (PR #830)
* `Vector_<T>` can now be constructed from an initializer list, e.g.
  `Vector v{1, 2, 3};`. (PR #871)
* Added the `SimTK_NODISCARD` macro and applied it to non-mutating Simmatrix
  operators (e.g., `normalize()`, `transpose()`) in `Vec`, `Mat`, `Row`,
  `SymMat`, `Vector`, `Matrix`, and `RowVector`. Discarding the result of one of
  these methods now produces a compiler warning. (PR #865)

**Bug fixes**
* `MobilizedBody::Ellipsoid::setQToFitTranslation()` (and therefore
  `setQToFitTransform()`) now places the M frame's origin exactly on the ray
  towards the requested translation, for any semi-axes. Similarly,
  `MobilizedBody::Ellipsoid::setUToFitLinearVelocity()` (and
  `setUToFitVelocity()`) now reproduces any velocity the joint can produce, for
  any semi-axes. Previously, both methods relied on calculations assuming the
  choice of semi-axes represented a sphere. (PR #875)
* Fixed implementations of `setUToFitAngularVelocityImpl()` for the Gimbal and
  Bushing mobilizers. (PR #844)
* `CPodesIntegrator` now correctly returns `EndOfSimulation` after repeated
  calls to `stepTo()` (issue #824). (PR #825)
* Fixed a race condition in `simbody-visualizer` that could cause it to crash
  on startup when menus or sliders were added while being drawn. (PR #838)
* Fixed out-of-bounds accesses (and the resulting compiler warnings) in `Mat`
  constructors and `Row`/`Vec` element layout. (PR #840)
* Fixed a memory leak on copy assignment of `Factor` classes (`FactorLU`,
  `FactorQTZ`, `FactorSVD`, `Eigen`, etc.); these now also support move
  construction and assignment. (PR #841)
* `ExponentialSpringForce` no longer retains a reference to the caller's
  `MobilizedBody` handle, which could dangle or change the body the spring acts
  on (issue #869). (PR #870)
* The OBJ mesh file parser now accepts any whitespace (tabs, repeated spaces,
  etc.) between the vertices of a face (`f`) record (issue #747). (PR #863)
* Fixed an ambiguous comparison in a `CablePath` assert that failed to compile
  with some compilers in Debug builds. (PR #851)

**Build, packaging, and documentation**
* Fixed the generated `simbody.pc` pkg-config file: its `Version:` field was
  empty (also in 3.8), its `Cflags` no longer forces `-std=c++11`, and `libdir`
  and `includedir` are now correct when the install directories are given as
  absolute paths. The visualizer install location is also correct when
  `CMAKE_INSTALL_LIBEXECDIR` is absolute. Added a Nix flake and Nix CI builds.
  (PR #878)
* Simbody can now be built as a CMake subproject (e.g., via
  `add_subdirectory()`). (PR #721)
* When building static libraries (`BUILD_SHARED_LIBS=OFF`), the vendored
  BLAS/LAPACK libraries are now installed (previously they were omitted).
  (PR #861)
* `CMakePresets.json` no longer hard-codes a Visual Studio version, so the MSVC
  presets work with any installed Visual Studio (e.g., 2022 or 2026). (PR #862)
* Simbody now configures without warnings under current CMake (tested with 4.4). The
  installed examples' `CMakeLists.txt` now requires CMake 3.21 (it previously declared 2.8,
  which CMake 4 rejects) and no longer forces C++11. (PR #874)
* `cmake --install --prefix` and `DESTDIR` installs now place `simbody-visualizer` and the
  `examples/bin` symlink under the requested location (both previously used the configure-time
  prefix), and `uninstall` now removes the symlink. (PR #874)
* Added a "For New Contributors" section to CONTRIBUTING.md. (PR #755)


3.8 (May 2025)
--------------------
* Added CableSpan as a superior replacement for CablePath. (PRs #791 and #814)
* Many improvements/modernization of the CMake build system. (PR #801)
* Minor changes for C++20 compatibility.
* Added calcForceContributionsSum() to GeneralForceSubsystem. (PR #788)
* Added methods for dealing with the position and velocity (PV) subset of
  the full constraint Jacobian G: multiplyByPV(), etc. (PR #783)
* Improved performance of reading mesh files. (PR #773)
* Added ExponentialSpringForce, an efficient point-ground contact model.
  (PRs #720, #734, #746 and others)
* Many minor improvements and bug fixes, see Simbody GitHub for details.

3.7 (December 2019)
-------------------
* The new SmoothSphereHalfSpaceForce provides a continuous and differentiable
  contact model, ideal for use with gradient-based optimization algorithms (PR
  #667).
* Fixed a memory issue with CPodes (PR #642).
* Created a CMake variable to control whether docs are installed (PR #655).
* Fixed a bug with calculating constraint acceleration errors (PR #670).
* Fixed Pathname::getThisExecutablePath() for FreeBSD (PR #672).
* Fixed scaling issue of simbody-visualizer on macOS 10.15 Catalina. Now,
  simbody-visualizer is an app bundle (simbody-visualizer.app) on Mac (PR #676).

3.6.1 (11 June 2018)
--------------------
* Fixed bug wherein a program may crash when using the visualizer if the
  visualizer window was closed manually (PR #633).

3.6 (21 February 2018)
----------------------
* Forced GCC to be at least 4.9.0, so that new C++11 features like regular
  expressions can be used (PR #485).
* Minimum Ubuntu version supported 14.04 LTS (Trusty), with a [manual update of GCC](http://askubuntu.com/questions/466651/how-do-i-use-the-latest-gcc-on-ubuntu-14-04).
* Recommended minimum Ubuntu version : 15.04 (Vivid), that is shipped with GCC 4.9.2.
* Added mixin classes `ResetOnCopy<T>` and `ReinitOnCopy<T>` to force default
  construction or reinitialization of data members on copy construction or copy
  assignment, without requiring a user written copy constructor and copy
  assignment operator just to get those reinitializations done.
* Added clone() method to `SimTK::Function_` base class and implemented it for
  Simbody-defined concrete Function classes. Made concrete Function members
  non-const to permit assignment, and modified Function_<T>::Step to allow
  changing its parameters after construction.
* Added C++11 features to the `SimTK::Array_` container including `std::initializer_list` construction, move construction, move assignment, and `emplace` methods.
* Prevented copy construction of `Array_<T>` from `Array_<T2>` unless T2 is *implicitly*
  convertible to T. Previously this was allowed if there was any conversion possible
  even if it was explicit. Array_ was being too relaxed about this, causing hidden 
  copies to occur. 
* Added CloneOnWritePtr smart pointer (acts like ClonePtr but with deferred cloning).
* Updated ClonePtr and ReferencePtr APIs to follow C++11 standard smart pointer
  terminology. This required deprecating some existing methods and operators, so
  you can expect to get annoying warnings until you switch to the new API. 
* Possible BREAKING CHANGE: ClonePtr's operator==() previously delegated
  to the managed object; now it just operates on the managed pointer as is done 
  in other smart pointers. Consequently now only a clone() method is required for a type
  to be contained in a ClonePtr; previously it had to support comparison also.
* Make doxygen run silently so errors will be easier to see.
* Added new methods to `Pathname` class for interpreting pathnames against a specified working directory instead
of the current working directory (thanks to Carmichael Ong). See [Issue #264](https://github.com/simbody/simbody/issues/264) and [PR #307](https://github.com/simbody/simbody/pull/307). 
* Added the ability to parallelize forces, robustness/performance improvements
to ParallelExecutor, mutex state lock.
[PR #414](https://github.com/simbody/simbody/pull/414).
* Added move constructor and move assignment to State (very fast).
* Added "stage version" counters for time, q, u, and z that are incremented
  whenever one of these changes.
* Separated time-independent position and kinematics calculations so that 
  they are not invalidated by a time change. These can also be initiated
  explicitly with new methods `realizePositionKinematics()` and
  `realizeVelocityKinematics()`. They are invalidated by a change to q or
  to u, respectively. 
* Modified floating point-to-String conversions to use lossless number of 
  digits by default. An attempt to use the default type-to-String conversion 
  when no stream insertion operator is available is now a runtime rather than
  compile time error. [PR #459](https://github.com/simbody/simbody/pull/459).
* Upgraded `SimTK::Xml` from class to namespace; local classes like `Xml::Element`
  were promoted to namespace level within Xml. This was necessary to untangle
  Array_ and Xml classes which can be mutually dependent. This change is backwards
  compatible for those following the recommended use of the previously-available
  `Xml::Document` typedef instead of Xml directly. Otherwise you will have to change
  uses of `SimTK::Xml` to `SimTK::Xml::Document`. Also any `using SimTK::Xml;`
  must be removed or replaced with `using namespace SimTK::Xml` or
  `using SimTK::Xml::Document` depending on the intent. 
  [PR #460](https://github.com/simbody/simbody/pull/460)
* Improved `NiceTypeName<T>::namestr()` to produce a canonicalized name that is
  the same on all platforms (with a few exceptions). Added `xmlstr()` method to
  make an XML-friendly modification of `namestr()` that replaces angle brackets
  with curly braces. Added a new regression test to verify that the names come
  out right. [PR #461](https://github.com/simbody/simbody/pull/461)
* Added helper class IteratorRange to use range-based for loops with a pair of
  iterators. [PR#467](https://github.com/simbody/simbody/pull/467)
* Added the method `State::isConsistent()` to compare two states
  [PR #469](https://github.com/simbody/simbody/pull/469).
* Started using RPATH on OSX so that users need not set `DYLD_LIBRARY_PATH` to
  run `simbody-visualizer` or the example executables, regardless of where you
  install Simbody.
* Improved the ability to find the simbody-visualizer executable when Simbody
  is installed in non-standard locations or if the Simbody installation is
  relocated (even to different computers). This enhancement is only for
  non-Windows operating systems.
* Fixed a bug when compiling on macOS (OSX) with SDK MacOSX10.12.sdk, related
  to the POSIX function `clock_gettime()`.
  [Issue #523](https://github.com/simbody/simbody/issues/523),
  [PR #524](https://github.com/simbody/simbody/pull/524)
* Replaced usages of pthreads with C++11 classes. Removed the following classes:
    - ThreadLocal
    - AtomicInteger
* Deprecated some badly-named methods in MobilizedBody::Translation
  [Issue #604](https://github.com/simbody/simbody/issues/604)
* Removed support for `long double`. 
  [PR #597](https://github.com/simbody/simbody/pull/597)
* Fixed a bug in `MultibodyGraphMaker` where massless bodies were handled
  incorrectly.
  [PR #592](https://github.com/simbody/simbody/pull/592) and 
  [PR #594](https://github.com/simbody/simbody/pull/594).
* Moved Ipopt into the namespace `SimTKIpopt` to avoid duplicate symbols when
  combining Simbody with an independent Ipopt.
* Fixed a bug where RowVectors could not be constructed with 0 elements.
* CMake: Downstream projects no longer need to use `include_directories()`.
* CMAES: You can now specify a different initial step size (standard deviation)
  for each variable. [PR #540](https://github.com/simbody/simbody/pull/540)
* SimbodyMatterSubsystem now provides a function for the error in
  acceleration-level constraint equations.
  [PR #517](https://github.com/simbody/simbody/pull/517)
* The visualizer now supports lazily loading mesh geometry from a file.
  [PR #475](https://github.com/simbody/simbody/pull/475)
* Simbody can now be built with MinGW on Windows.
  [PR #441](https://github.com/simbody/simbody/pull/441)


3.5.3 (15 June 2015)
-------------------
This is the release that shipped with OpenSim 3.3.
* Small changes to allow compilation with Visual Studio 2015 (PRs [#395](https://github.com/simbody/simbody/pull/395) and [#396](https://github.com/simbody/simbody/pull/396)).
* Fixed a problem with SpatialInertia::shift() with non-zero COM offset, see issue [#334](https://github.com/simbody/simbody/issues/334). This also affected calcCompositeBodyInertias(). These are not commonly used.
* Fixed a problem with VectorIterator which could unnecessary copying, possibly affecting mesh handling performance. See issue [#349](https://github.com/simbody/simbody/issues/349). 


3.5.2 (15 May 2015)
-------------------
Same as 3.5.1 except on 64 bit Windows which has a patched version of Lapack that
addresses an error handling problem that caused trouble for some OpenSim users.
This is a patch to Lapack 3.4.2 (64 bit) to fix the bug discussed in [Issue #177](https://github.com/simbody/simbody/issues/177) and [PR #342](https://github.com/simbody/simbody/pull/342).
There were two functions where convergence failures incorrectly caused an abort (XERBLA in Lapack-speak). See discussion on Lapack forum:
http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=4586
This Lapack DLL is binary compatible with the previous one, same functions and ordinals.


3.5.1 (31 Dec 2014)
-------------------
This patch release fixed an installation problem but is otherwise identical to 3.5.
### Bugs fixed
* Fixed a Mac installation problem reported [here](https://github.com/osrf/homebrew-simulation/issues/33).
See [PR #321](https://github.com/simbody/simbody/pull/321). Probably affects Linux too; we think it is due to a change in CMake behavior, noticed with CMake 3.1 but could have happened earlier.
* Updated Travis-CI script to attempt installation to catch these problems earlier. See [PR #322](https://github.com/simbody/simbody/pull/322).


3.5 (18 Dec 2014)
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
* Improved regression test timing framework to make it useful for multithreaded tests (pr #265)
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
* Added CONTRIBUTING.md file explaining ways to contribute to the Simbody project.


3.4.1 (31 Mar 2014)
-------------------

This is primarily a release for improving our build and install process to comply with Debian's requirements. Thanks to José Rivero and Steve Peters at Open Source Robotics Foundation, and Chris Dembia at Stanford for the bulk of this effort.

* Fixed `SimbodyMatterSubsystem::getTotalCentrifugalForces()` (issue #112, pr #116).
* Many changes to build and install, mostly affecting Linux and OSX. Should now conform better to standards on those platforms and in general be better and finding its dependencies. (pr #114 #118 #120 #121 #122 #127 #130 #131 #133, issue #48 #68 #92 #113) 
* Started using the -Werror flag (treat warnings as errors) in Travis builds to ensure that we are warning-free. (issue #128, pr #129 #131)
* Compile with C++11 enabled by default. However, the code will still build in C++03. (pr #123, #126)


3.3.1 (21 Jan 2014)
-------------------
This is the first release built for use in Open Source Robotic Foundation's Gazebo robot simulator and is also the version of Simbody that ships with OpenSim 3.2. It incorporates many fixes and enhancements prompted by the integration effort with OSRF, and a new Debian package builder for smooth incorporation into the Gazebo build.

* Improved matrix/vector documentation and reorganized source to break up large files (pr #70 #76 #77 #78 #87)
* Improved Force::Gravity to make it more flexible and avoid NaNs (pr #33)
* Fix multiplyByMInv() to be usable in forces and controllers; now only requires state to be realized to position stage. (issue #29, pr #31)
* Make MobilizedBody::lock() and lockAt() method immediately modify the state (issue #20, pr #23)
* Added Debian/Ubuntu package building (pr #24 #32 #52 #58 #59 #64)
* Improved TestMultibodyPerformance timing information (pr #22)
* Fixed bug in MobilizedBody::isSameMobilizedBody() (issue #15, pr #16)
* Started using Travis-CI for continuous integration testing (pr #25)
* Improved DecorativeGeometry classes (issue #34, pr #35 #61)
* Improved installation, esp. Linux and OSX (issue #38 #65 #101, pr #64 #91 #102 #107)
* Added uninstall (issue #104, pr #106)
* Visualizer name changed to "simbody-visualizer" (issue #27, pr #53)
* Added torque-limited motor examples
* Added ability to lock/unlock mobilizers and disble/enable Motion objects on the fly, for flexible mixed forward/inverse dynamics
* Re-engineered Force::Gravity for speed and run time flexibility (for Gazebo)
* Extended Force::Gravity to support gravity compensation controllers
* Allow runtime changes to MobilityLinear{Spring,Damper} (for Gazebo)
* Added Theo Jansen Strandbeest example
* Ported some Gazebo regression tests to Simbody's regression suite.
* Added dirent.h support on Windows (for directory searching)
* Many bug fixes, doxygen improvements, and small performance improvements
* Added Semi-Explicit Euler integrators, with and without error control (for Gazebo)
* Added O(n) methods for task space Jacobians


3.1 (22 Apr 2013)
-----------------
This is the Simbody release that shipped with OpenSim 3.1 and contains the initial
work done with Open Source Robotics Foundation for the Gazebo simulator. The source was managed on Subversion, although its code and history were transferred to GitHub. The source zip is available on the GitHub release page [here](https://github.com/simbody/simbody/releases/tag/Simbody-3.1).

* Added MobilityLinearStop for joint stops (for Gazebo)
* Added Simbody example for reading Gazebo sdf files 
* Added MultibodyGraphMaker to map bodies and joints to mobilizers and constraints (for Gazebo)
* Added a variety of prototype implementations of unilateral contact as adhoc tests
* Added some pre-built mesh types to PolygonalMesh (for Gazebo)
* Modified Force::Gravity to allow excluded bodies (needed by Gazebo)
* Added DiscreteForces and MobilityDiscreteForces elements for Gazebo
* Added Measure::Delay for time-delay of input signal
* Added scissor-lift example to show behavior on highly constrained mechanisms
* Added RK2 integrator
* Modified Gimbal mobilizer to use Euler angle derivatives as generalized speeds
* Added Bushing (6 dof) mobilizer parameterized like Gimbal
* Added methods for calculating geodesics over smooth contact geometry
* Added infrastructure for fast cable wrapping and some prototypes
* Many small performance improvements and bug fixes
