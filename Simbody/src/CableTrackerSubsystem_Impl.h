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

int realizeSubsystemInstanceImpl(const State& state) const OVERRIDE_11 {
    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        const CablePath& path = getCablePath(ix);
        path.getImpl().realizeInstance(state);
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


int realizeSubsystemAccelerationImpl(const State& state) const OVERRIDE_11 {
    for (CablePathIndex ix(0); ix < cablePaths.size(); ++ix) {
        const CablePath& path = getCablePath(ix);
        path.getImpl().realizeAcceleration(state);
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
        const PathPosEntry&     ppe      = path.getPosEntry(state);

        Vec3 prevPoint;
        for (CableObstacleIndex ox(0); ox < path.getNumObstacles(); ++ox) {
            const CableObstacle::Impl& obs = path.getObstacleImpl(ox);
            Vec3 P_B, Q_B;
            obs.getContactStationsOnBody(state, instInfo, ppe, P_B, Q_B);
            const Transform& X_GB = obs.getBodyTransform(state);
            Vec3 P_G = X_GB*P_B, Q_G = X_GB*Q_B;

            const Vec3 color = ox==0 ? Green 
                             : ox==path.getNumObstacles()-1 ? Red : Purple;
            // Draw point at Pi.
            decorations.push_back(DecorativePoint(P_G).setColor(color));

            if (ox!=0) // draw straight line from Qi-1 to Pi
                decorations.push_back(DecorativeLine(prevPoint,P_G)
                    .setColor(Purple).setLineThickness(3));

            // If there are two distinct points, draw the second one and a
            // red curve connecting them.
            if ((Q_G-P_G).normSqr() >= square(SignificantReal)) {
                decorations.push_back(DecorativePoint(Q_G).setColor(color));

                ActiveSurfaceIndex asx = ppe.mapToActiveSurface[ox];
                assert(asx.isValid());
                const Geodesic& g = ppe.geodesics[asx];
                const Transform& X_BS = obs.getObstaclePoseOnBody(state,instInfo);
                const Transform X_GS = X_GB*X_BS;
                const Array_<Transform>& Kf = g.getFrenetFrames();
                Vec3 prevP_G = X_GS*Kf.front().p();
                for (unsigned i=0; i < Kf.size(); ++i) {
                    Vec3 Q_G = X_GS*Kf[i].p();
                    decorations.push_back(DecorativeLine(prevP_G, Q_G)
                        .setColor(Red).setLineThickness(2));
                    if (i < Kf.size()-1)
                        decorations.push_back(DecorativePoint(Q_G)
                            .setColor(Blue).setScale(.5));

                    prevP_G = Q_G;
                }
            }
            prevPoint = Q_G;
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

