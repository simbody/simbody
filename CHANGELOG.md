# Simbody Changelog and Release Notes
**(work in progress)**

This is not a comprehensive list of changes but rather
a hand-curated collection of the more notable ones. For
a comprehensive history, see the [Simbody GitHub repo]
(https://github.com/simbody/simbody). You can use the release dates below to find all the PRs and issues that were 
included in a particular release. 

## 3.5 (in progress)

* Added CHANGELOG.md to track changes for versions.
* More


## 3.4.1 (31 Mar 2014)

This was primarily a release for improving our build and install process to comply with Debian's requirements. Thanks to José Rivero and Steve Peters at Open Source Robotics Foundation, and Chris Dembia at Stanford for the bulk of this effort.

* Fixed `SimbodyMatterSubsystem::getTotalCentrifugalForces()` (issue #112, pr #116).
* Many changes to build and install, mostly affecting Linux and OSX. Should now conform better to standards on those platforms and in general be better and finding its dependencies. (pr #114 #118 #120 #121 #122 #127 #130 #131 #133, issue #48 #68 #92 #113) 
* Started using the -Werror flag (treat warnings as errors) in Travis builds to ensure that we are warning-free. (issue #128, pr #129 #131)
* Compile with C++11 enabled by default. However, the code will still build in C++03. (pr #123, #126)


## 3.3.1 (21 Jan 2014)

This is the first release built for use in Open Source Robotic Foundation's Gazebo robot simulator and is also the version of Simbody that ships with OpenSim 3.2. It incorporates many fixes and enhancements prompted by the integration effort with OSRF, and a new Debian builder for smooth incorporation into the Gazebo build.


## 3.1 (15 Aug 2013)

This is the Simbody release that shipped with OpenSim 3.1. The source was managed on Subversion, although its code and history were transferred to GitHub, from which a determined code archeologist could presumably reconstruct the changelog. But why? The source zip is also available on the GitHub release page [here](https://github.com/simbody/simbody/releases/tag/Simbody-3.1), and does include a changelog named `ChangeLog.txt`.