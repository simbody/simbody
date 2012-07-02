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
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/CableTrackerSubsystem.h"
#include "simbody/internal/CablePath.h"

#include "CablePath_Impl.h"
#include "CableTrackerSubsystem_Impl.h"

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

std::ostream& operator<<(std::ostream& o, const PathInstanceInfo& info) {
    o << "  disabled="; 
    for (CableObstacleIndex ox(0); ox < info.obstacleDisabled.size(); ++ox)
        o << " " << String(info.obstacleDisabled[ox]);
    o << "\n  pose: ";
    for (CableObstacleIndex ox(0); ox < info.obstaclePose.size(); ++ox)
        o << " " << info.obstaclePose[ox].p();
    return o << endl;
}

std::ostream& operator<<(std::ostream& o, const PathPosEntry& entry) {
    cout << "PathPosEntry: length=" << entry.length << " straight dirs:\n  ";
    for (CableObstacleIndex i(0); i < entry.straightDirections.size(); ++i)
        cout << entry.straightDirections[i];
    cout << "\n  Unit forces:\n";
    for (CableObstacleIndex ox(0); ox < entry.obstacleUnitForces.size(); ++ox)
        cout << "  " << entry.obstacleUnitForces[ox] << endl;
    return o;
}

std::ostream& operator<<(std::ostream& o, const PathVelEntry& entry) {
    return o;
}
}

using namespace SimTK;

//==============================================================================
//                              CABLE PATH 
//==============================================================================

CablePath::CablePath
   (CableTrackerSubsystem&    cables,
    const MobilizedBody&      originBody,
    const Vec3&               defaultOriginPoint,
    const MobilizedBody&      terminationBody,
    const Vec3&               defaultTerminationPoint) : impl(0)
{
    impl = new Impl(cables);
    impl->referenceCount = 1;
    impl->index = cables.updImpl().adoptCablePath(*this);
    CableObstacle::ViaPoint(*this, originBody, defaultOriginPoint);
    CableObstacle::ViaPoint(*this, terminationBody, defaultTerminationPoint);
}

CablePath::CablePath
   (CableTrackerSubsystem&    cables,
    const MobilizedBody&      originBody,
    const MobilizedBody&      terminationBody) : impl(0)
{   // Invoke above constructor.
    new(this) CablePath(cables, originBody, Vec3(0), terminationBody, Vec3(0));
}

// Copy constructor is shallow and reference counted.
CablePath::CablePath(const CablePath& src) : impl(src.impl) 
{   if (impl) ++impl->referenceCount; }

// Copy assignment is shallow and reference counted.
CablePath& CablePath::operator=(const CablePath& src) {
    if (&src != this) {
        clear();
        if (src.impl) {
            impl = src.impl;
            ++impl->referenceCount;
        }
    }
    return *this;
}

// This is used to implement the destructor and copy assignment.
void CablePath::clear() {
    if (impl) {
        if (--impl->referenceCount == 0)
            delete impl;
        impl = 0;
    }
}

int CablePath::getNumObstacles() const 
{   return getImpl().obstacles.size(); }

const CableObstacle& CablePath::
getObstacle(CableObstacleIndex obstacleIx) const
{   return getImpl().obstacles[obstacleIx]; }

Real CablePath::getCableLength(const State& state) const 
{   return getImpl().getCableLength(state); }

Real CablePath::getCableLengthDot(const State& state) const 
{   return getImpl().getCableLengthDot(state); }

void CablePath::applyBodyForces(const State& state, Real tension, 
                     Vector_<SpatialVec>& bodyForcesInG) const
{   getImpl().applyBodyForces(state,tension,bodyForcesInG); }

Real CablePath::calcCablePower(const State& state, Real tension) const
{   return getImpl().calcCablePower(state,tension); }


//==============================================================================
//                             CABLE OBSTACLE
//==============================================================================


CableObstacle::CableObstacle(CableObstacle::Impl* guts) : impl(guts)
{
    if (impl) ++impl->referenceCount;
}

// Copy constructor is shallow and reference counted.
CableObstacle::CableObstacle(const CableObstacle& src) : impl(src.impl) 
{   if (impl) ++impl->referenceCount; }

// Copy assignment is shallow and reference counted.
CableObstacle& CableObstacle::operator=(const CableObstacle& src) {
    if (&src != this) {
        clear();
        if (src.impl) {
            impl = src.impl;
            ++impl->referenceCount;
        }
    }
    return *this;
}

// This is used to implement the destructor and copy assignment.
void CableObstacle::clear() {
    if (impl) {
        if (--impl->referenceCount == 0)
            delete impl;
        impl = 0;
    }
}

const Transform& CableObstacle::
getDefaultTransform() const {return impl->getDefaultPoseOnBody();}
const MobilizedBody& CableObstacle::
getMobilizedBody() const {return impl->getMobilizedBody();}
const CablePath& CableObstacle::
getCablePath() const {assert(impl->cablePath); return *impl->cablePath;}
CableObstacleIndex CableObstacle::
getObstacleIndex() const {return impl->index;}
const DecorativeGeometry& CableObstacle::
getDecorativeGeometry() const {return impl->decoration;}

CableObstacle& CableObstacle::
setDefaultTransform(const Transform& X_BS) 
{   impl->defaultX_BS=X_BS; return *this;}
CableObstacle& CableObstacle::
setDecorativeGeometry(const DecorativeGeometry& viz)
{   impl->decoration = viz; return *this; }

inline void CableObstacle::Impl::invalidateTopology() 
{   if (cablePath) cablePath->updImpl().invalidateTopology(); }   


//==============================================================================
//                         CABLE OBSTACLE :: VIA POINT
//==============================================================================

CableObstacle::ViaPoint::
ViaPoint(CablePath& path, const MobilizedBody& viaMobod,
         const Vec3& defaultStation)
:   CableObstacle(new ViaPoint::Impl(path, viaMobod, defaultStation))
{
    impl->index = path.updImpl().adoptObstacle(*this); 
}



//==============================================================================
//                         CABLE OBSTACLE :: SURFACE
//==============================================================================

CableObstacle::Surface::
Surface(CablePath& path, const MobilizedBody& viaMobod)
:   CableObstacle(new Surface::Impl(path, viaMobod, Transform()))
{
    impl->index = path.updImpl().adoptObstacle(*this); 
}

CableObstacle::Surface& CableObstacle::Surface::
setPathPreferencePoint(const Vec3& stationInS) { 
    Impl& surfImpl = dynamic_cast<Impl&>(*impl);
    surfImpl.pathPreferencePointInS = stationInS; 
    return *this; 
}

CableObstacle::Surface& CableObstacle::Surface::
setContactSurface(const Transform& X_BS, const ContactGeometry& surface) {
    impl->defaultX_BS = X_BS;
    Impl& surfImpl = dynamic_cast<Impl&>(*impl);
    surfImpl.surface = surface; 
    return *this; 
}

