#ifndef SimTK_SIMBODY_CABLE_PATH_H_
#define SimTK_SIMBODY_CABLE_PATH_H_

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

/** @file
This file declares the CablePath and CableObstacle classes. **/

#include "SimTKmath.h"
#include "simbody/internal/common.h"

namespace SimTK {

/** @class SimTK::CableObstacleIndex
This is a unique integer type for identifying obstacles comprising a particular
cable path. These begin at zero for each cable path. **/
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableObstacleIndex);

class CableTrackerSubsystem;
class MobilizedBody;
class CableObstacle;

//==============================================================================
//                                CABLE PATH
//==============================================================================
/** This class represents the path of one cable from its origin point, through
via points, around obstacles, to its termination point. The cable ends are
fixed to the origin and termination points. The cable passes through each via
point, and takes the shortest allowable path over the surfaces of obstacles. 

During initialization, if there is more than one possible geodesic over the
surface, or if a straight line path would miss the surface altogether, we'll
take the shortest route unless the user has provided a "preference point" on the
surface. In that case whichever path segment runs closer to the preferred point
is chosen. During continuation, only local path movement is allowed so the
path segment will not flip from one side to the other once it has been
initialized. The preferred point is ignored during continuation.

<h3>Notation</h3>

For convenience, we include the origin and termination points as obstacles,
with the origin being obstacle 0, followed by m via point and surface 
obstacles numbered 1 to m, followed by the termination point as obstacle t=m+1.
Every obstacle is represented by two "contact points", P and Q, which we'll
number Pi and Qi for obstacle i. The obstacles are
separated by straight-line segments running from Qi-1 to Pi. Starting at 
the origin Q0, the path first touches the surface of obstacle 1 at P1, 
travels over the surface to Q1, and then leaves the surface towards the 
termination point Pt. That is, Pi's path coordinate must be less than Qi's. 
The segment from P to Q is a geodesic over the surface. For via points,
points P and Q are in the same location but there are two different tangents
associated with them in the incoming and outgoing straight-line directions.
For the origin obstacle, only point Q is relevant and for the termination
obstacle only point P is relevant.

**/
class SimTK_SIMBODY_EXPORT CablePath {
public:

/** Create a straight-line cable path connected a point fixed on one body with
one fixed on another body. You can add additional obstacles later. **/
CablePath(CableTrackerSubsystem&    cables,
          const MobilizedBody&      originBody,
          const Vec3&               defaultOriginPoint,
          const MobilizedBody&      terminationBody,
          const Vec3&               defaultTerminationPoint);

/** Constructor taking only the origin and terminal bodies with the expectation
that you'll set the end point locations later. The default locations are set
to the body frame origins here. **/
CablePath(CableTrackerSubsystem&    cables,
          const MobilizedBody&      originBody,
          const MobilizedBody&      terminationBody);

/** Copy constructor is shallow and reference counted. **/
CablePath(const CablePath& source);

/** Copy assignment is shallow and reference counted. **/
CablePath& operator=(const CablePath& source);

/** Delete the cable path if this handle was the last reference to it. **/
~CablePath() {clear();}

/** Return the total number of obstacles (origin point, surfaces and via 
points, and termination point) that were provided for this cable path, 
regardless of whether they are currently in use. **/
int getNumObstacles() const;
/** Return a const reference to one of the obstacles along this path, given
by its index starting at zero for the origin point. **/
const CableObstacle& getObstacle(CableObstacleIndex obstacleIx) const;

/** Return the total length of the cable that was calculated for the
configuration supplied in \a state. State must have been realized through
Position stage. **/
Real getCableLength(const State& state) const;

/** Return the cable rate (time derivative of cable length) that was 
calculated for the configuration and velocities supplied in \a state. State 
must have been realized through Velocity stage. Calculation of the cable rate
may be initiated if necessary the first time this is called at this state
but will be saved in the cache for subsequent accesses. **/
Real getCableLengthDot(const State& state) const;

/** Given a tension > 0 acting uniformly along this cable, apply the resulting
forces to the bodies it touches. The body forces are added into the appropriate
slots in the supplied Array which has one entry per body in the same format
as is supplied to the calcForce() method of force elements. If the supplied
tension is <= 0, signifying a slack cable, this method does nothing. **/
void applyBodyForces(const State& state, Real tension, 
                     Vector_<SpatialVec>& bodyForcesInG) const;

/** Calculate the power this cable would apply or dissipate at the given
\a tension (>0) value, using the velocities provided in the given \a state,
which must already have been realized to Velocity stage. Power is positive if
the cable is adding energy to the system; negative when dissipating. If the 
supplied \a tension is <= 0, signifying a slack cable, this method returns 
zero. **/
Real calcCablePower(const State& state, Real tension) const;


/** Default constructor creates an empty cable path not associated with any
subsystem; don't use this. **/
CablePath() : impl(0) {}
class Impl;
const Impl& getImpl() const {assert(impl); return *impl;}
Impl& updImpl() {assert(impl); return *impl;}
//--------------------------------------------------------------------------
private:
void clear();
Impl*   impl;
};


//==============================================================================
//                                CABLE OBSTACLE
//==============================================================================
/** An obstacle is any significant object along the cable path -- one of the
end points, a via point, or a surface. This is the base class that can refer
to any kind of obstacles; specific types are derived from this one. **/
class SimTK_SIMBODY_EXPORT CableObstacle {
public:
class ViaPoint; // also used for end points
class Surface;

/** Create an empty obstacle handle that can refer to any type of obstacle. **/
CableObstacle() : impl(0) {}

/** Insert this obstacle into the given cable path. **/
explicit CableObstacle(CablePath& path);
/** Copy constructor is shallow and reference-counted; this handle will 
point to the same object as does the \a source. **/
CableObstacle(const CableObstacle& source);
/** Copy assignment is shallow and reference-counted; this handle will 
point to the same object as does the \a source. **/
CableObstacle& operator=(const CableObstacle& source);
/** Destructor clears the handle, deleting the referenced object if this
was the last reference. **/
~CableObstacle() {clear();}

/** Return the default pose X_BS of the obstacle S on its body B. For a via
point, the point is located at the S frame origin whose position vector in B
is X_BS.p(). **/
const Transform& getDefaultTransform() const;
/** Get a reference to the Mobilized body to which this obstacle is fixed.
There can be multiple objects on a singe body, so this is not necessarily 
unique within a path. **/
const MobilizedBody& getMobilizedBody() const;
/** Return a reference to the CablePath in which this obstacle resides. **/
const CablePath& getCablePath() const;
/** Return the obstacle index within this obstacle's CablePath. **/
CableObstacleIndex getObstacleIndex() const;
/** Return decorative geometry that can be used to display this obstacle. The
decorative geometry's coordinate frame is coincident with the obstacle's 
coordinate frame S. **/
const DecorativeGeometry& getDecorativeGeometry() const;

CableObstacle& setDefaultTransform(const Transform& X_BS);
CableObstacle& setDecorativeGeometry(const DecorativeGeometry& viz);

/** Clear this handle, deleting the referenced object if this
was the last reference.  **/
void clear();
/** See if this handle is empty. **/
bool isEmpty() const {return impl==0;}

//--------------------------------------------------------------------------
class Impl;
const Impl& getImpl() const {assert(impl); return *impl;}
Impl&       updImpl()       {assert(impl); return *impl;}

protected:
explicit CableObstacle(Impl* impl);

private:
Impl*   impl; // opaque pointer to reference-counted implementation object
};


//==============================================================================
//                      CABLE OBSTACLE :: VIA POINT
//==============================================================================
/** This is a point through which the cable must pass. **/
class SimTK_SIMBODY_EXPORT CableObstacle::ViaPoint : public CableObstacle {
public:
/** Default constructor creates an empty handle. **/
ViaPoint() : CableObstacle() {}
/** Insert a via point obstacle to the given cable path. **/
ViaPoint(CablePath& path, const MobilizedBody& viaMobod, 
         const Vec3& defaultStation);

/** Return true if the given CableObstacle is a ViaPoint. **/
static bool isInstance(const CableObstacle&);
/** Cast the given CableObstacle to a const ViaPoint; will throw an exception
if the obstacle is not a via point. **/
static const ViaPoint& downcast(const CableObstacle&);
/** Cast the given CableObstacle to a writable ViaPoint; will throw an 
exception if the obstacle is not a via point. **/
static ViaPoint& updDowncast(CableObstacle&);

class Impl;
};

//==============================================================================
//                      CABLE OBSTACLE :: SURFACE
//==============================================================================
/** This obstacle is a solid object represented by a ContactGeometry surface.
The cable cannot penetrate the surface and will instead take the shortest
path in an allowed direction on the surface. **/
class SimTK_SIMBODY_EXPORT CableObstacle::Surface : public CableObstacle {
public:
/** Default constructor creates an empty handle. **/
Surface() : CableObstacle() {}

Surface(CablePath& path, const MobilizedBody& body);
Surface& setContactSurface(const Transform& X_BS, 
                           const ContactGeometry& surface);
Surface& setDecorativeGeometry(const DecorativeGeometry& viz)
{   CableObstacle::setDecorativeGeometry(viz); return *this; }
/** Optionally provide a "preference point" that can be used during
path initialization to disambiguate when there is more than one geodesic
that can connect the contact points. The geodesic that passes nearest
the preference point will be selected. Without this point the shortest
geodesic will be selected. This point is ignored during path continuation
calculations. The point location is given in the local frame S of the
contact surface. **/
Surface& setPathPreferencePoint(const Vec3& point);
    
/** Optionally provide some hints for the initialization algorithm to
use as starting guesses for the contact point locations. These are
ignored during path continuation calculations. These locations are given
in the local frame S of the contact surface. **/
Surface& setContactPointHints(const Vec3& startHint, 
                              const Vec3& endHint);

/** Return true if the given CableObstacle is a Surface. **/
static bool isInstance(const CableObstacle&);
/** Cast the given CableObstacle to a const Surface; will throw an exception
if the obstacle is not a via point. **/
static const Surface& downcast(const CableObstacle&);
/** Cast the given CableObstacle to a writable Surface; will throw an 
exception if the obstacle is not a via point. **/
static Surface& updDowncast(CableObstacle&);
class Impl;
};


} // namespace SimTK

#endif // SimTK_SIMBODY_CABLE_PATH_H_
