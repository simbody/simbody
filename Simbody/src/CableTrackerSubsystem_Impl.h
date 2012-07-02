#ifndef SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_IMPL_H_
#define SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_IMPL_H_

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

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

//==============================================================================
//                    CABLE TRACKER SUBSYSTEM :: IMPL
//==============================================================================
// This is the private implementation of CableTrackerSubsystem.   
    
class CableTrackerSubsystem::Impl : public Subsystem::Guts {
public:
// Constructor registers a default set of Trackers to use with geometry
// we know about. These can be overridden later.
Impl() {}

~Impl() {}

Impl* cloneImpl() const OVERRIDE_11 
{   return new Impl(*this); }

int getNumCablePaths() const {return cablePaths.size();}

const CablePath& getCablePath(CablePathIndex index) const 
{   return cablePaths[index]; }

CablePath& updCablePath(CablePathIndex index) 
{   return cablePaths[index]; }

// Add a cable path to the list, bumping the reference count.
CablePathIndex adoptCablePath(CablePath& path) {
    cablePaths.push_back(path);
    return CablePathIndex(cablePaths.size()-1);
}

// Return the MultibodySystem which owns this ContactTrackerSubsystem.
const MultibodySystem& getMultibodySystem() const 
{   return MultibodySystem::downcast(getSystem()); }

// Return the SimbodyMatterSubsystem from which this ContactTrackerSubsystem
// gets the bodies to track.
const SimbodyMatterSubsystem& getMatterSubsystem() const 
{   return getMultibodySystem().getMatterSubsystem(); }

// Get access to state variables and cache entries.

// TODO

// Allocate state variables.
int realizeSubsystemTopologyImpl(State& state) const OVERRIDE_11 {
    // Briefly allow writing into the Topology cache; after this the
    // Topology cache is const.
    Impl* wThis = const_cast<Impl*>(this);

    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        CablePath& path = wThis->updCablePath(ix);
        path.updImpl().realizeTopology(state);
    }

    return 0;
}

int realizeSubsystemPositionImpl(const State& state) const OVERRIDE_11 {
    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        const CablePath& path = getCablePath(ix);
        path.getImpl().realizePosition(state);
    }
    return 0;
}

int realizeSubsystemVelocityImpl(const State& state) const OVERRIDE_11 {
    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        const CablePath& path = getCablePath(ix);
        path.getImpl().realizeVelocity(state);
    }
    return 0;
}

int calcDecorativeGeometryAndAppendImpl
   (const State&                state, 
    Stage                       stage, 
    Array_<DecorativeGeometry>& decorations) const OVERRIDE_11 
{
    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        const CablePath::Impl&  path     = getCablePath(ix).getImpl();
        const PathInstanceInfo& instInfo = path.getInstanceInfo(state);

        Vec3 prevPoint;
        for (CableObstacleIndex ox(0); ox < path.getNumObstacles(); ++ox) {
            const CableObstacle::Impl& obs = path.getObstacleImpl(ox);
            const Transform& X_BS = obs.getObstaclePoseOnBody(state,instInfo);
            const MobilizedBody& bodyB = obs.getMobilizedBody();
            const Vec3 So_G = bodyB.findStationLocationInGround(state,
                                                                X_BS.p());
            const Vec3 color = ox==0 ? Green 
                             : ox==path.getNumObstacles()-1 ? Red : Purple;
            decorations.push_back(DecorativePoint(So_G).setColor(color));

            if (ox!=0)
                decorations.push_back(DecorativeLine(prevPoint,So_G)
                    .setColor(Purple).setLineThickness(3));
            prevPoint = So_G;
        }
    }
    return 0;
}


SimTK_DOWNCAST(Impl, Subsystem::Guts);

private:
// TOPOLOGY STATE
Array_<CablePath, CablePathIndex> cablePaths;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_CABLE_TRACKER_SUBSYSTEM_IMPL_H_

