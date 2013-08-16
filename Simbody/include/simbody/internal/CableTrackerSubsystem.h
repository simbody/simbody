#ifndef SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_H_
#define SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Ian Stavness, Andreas Scholz                     *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKmath.h"
#include "simbody/internal/common.h"

namespace SimTK {

/** @class SimTK::CablePathIndex
This is a unique integer type for quickly identifying specific cables for 
fast lookup purposes. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CablePathIndex);

class MultibodySystem;
class CablePath;

//==============================================================================
//                          CABLE TRACKER SUBSYSTEM
//==============================================================================
/** This subsystem tracks the paths of massless, frictionless cables that take
the shortest route between two distant points of a multibody system, passing 
smoothly over geometric obstacles that are attached to intermediate bodies. 
The calculated path will consist of a series of straight line segments 
between obstacles, and geodesics over the obstacles.

Force elements defined elsewhere may make use of cable paths to apply forces
to the system, by calculating a uniform tension in the cable that may depend 
on the cable kinematics calculated here. Cable kinematics includes the path,
the cable length, and the cable "rate", defined as the time derivative of 
length. The path and length are available at Position stage, the rate is
available at Velocity stage.

During construction, one or more CablePath objects are defined by giving for
each CablePath an origin and end point and an ordered set of geometric 
obstacles represented either by surfaces or "via" points. Via points are like 
frictionless eyelets that the cable must pass through and can generate forces
in any direction perpendicular to the cable; surfaces are one-sided and can
only apply positive forces to the cable. Thus the cable path does not 
necessarily touch all the obstacles; the obstacles that it does touch are
called the "active" obstacles. Via points are always active.

Every obstacle and point is rigidly fixed to Ground or some moving body of
the multibody system, with its pose or station point provided. Any number of
obstacles may be placed on one body. **/
class SimTK_SIMBODY_EXPORT CableTrackerSubsystem : public Subsystem {
public:
CableTrackerSubsystem();
explicit CableTrackerSubsystem(MultibodySystem&);

/** Get the number of cable paths being managed by this cable tracker subsystem.
These are identified by CablePathIndex values from 0 to getNumCablePaths()-1. 
This is available after realizeTopology() and does not change subsequently. **/
int getNumCablePaths() const;

/** Get const access to a particular cable path. **/
const CablePath& getCablePath(CablePathIndex cableIx) const;
/** Get writable access to a particular cable path. **/
CablePath& updCablePath(CablePathIndex cableIx);

/** @cond **/ // Hide from Doxygen.
SimTK_PIMPL_DOWNCAST(CableTrackerSubsystem, Subsystem);
class Impl;
Impl& updImpl();
const Impl& getImpl() const;
/** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_H_
