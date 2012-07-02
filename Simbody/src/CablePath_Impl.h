#ifndef SimTK_SIMBODY_CABLE_PATH_IMPL_H_
#define SimTK_SIMBODY_CABLE_PATH_IMPL_H_

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

#include <cassert>
#include <iostream>
using std::cout; using std::endl;

namespace SimTK {

// This is a discrete state variable for holding a path's instance-level 
// information, including placement of points and surfaces on their bodies
// and obstacle enable/disable settings.
class PathInstanceInfo {
public:
    PathInstanceInfo() {}
    // One entry per obstacle but you can't disable the origin or termination
    // points.
    Array_<bool,CableObstacleIndex>         obstacleDisabled;
    Array_<Transform,CableObstacleIndex>    obstaclePose;
};


// This is a cache entry for holding a path's calculated position-level 
// information.
class PathPosEntry {
public:
    PathPosEntry() : length(NaN) {}
    Vector x; // unknowns

    // This is the total length of the path.
    Real length;

    // These are unit vectors along the straight segments following the Q
    // point of each obstacle except the termination point.
    Array_<UnitVec3,CableObstacleIndex> straightDirections;

    // Multiply these by tension to get body forces in Ground.
    Array_<SpatialVec,CableObstacleIndex> obstacleUnitForces;
};


// This is a cache entry for holding a path's calculated velocity-level 
// information.
class PathVelEntry {
public:
    PathVelEntry() : lengthDot(NaN), unitPower(NaN) {}
    Vector xdot;    // calculated time derivatives of x
    Real lengthDot;
    Real unitPower; // unitForces*body velocities
};

std::ostream& operator<<(std::ostream& o, const SimTK::PathInstanceInfo& info);
std::ostream& operator<<(std::ostream& o, const SimTK::PathPosEntry& entry);
std::ostream& operator<<(std::ostream& o, const SimTK::PathVelEntry& entry);

//==============================================================================
//                         CABLE OBSTACLE :: IMPL
//==============================================================================
// This is the internal implementation abstract base class for CableObstacle. 
// It is reference counted so that CableObstacle objects may be multiply 
// referenced.
class CableObstacle::Impl {
public:
    explicit Impl(CablePath& cablePath, const MobilizedBody& mobod,
                  const Transform& defaultPose)
    :   referenceCount(0), cablePath(&cablePath), mobod(mobod),
        defaultX_BS(defaultPose) {}

    virtual ~Impl() {assert(referenceCount == 0);}

    // Concrete CableObstacle classes must implement these.
    virtual void getContactPointsOnObstacle(const State&, 
                                        const PathInstanceInfo&,
                                        const PathPosEntry&,
                                        Vec3& P_S, Vec3& Q_S) const = 0;
    virtual Real getSegmentLength(const State&, 
                                  const PathInstanceInfo&,
                                  const PathPosEntry&) const = 0;
    virtual Real getSegmentLengthDot(const State&, 
                                     const PathInstanceInfo&,
                                     const PathPosEntry&,
                                     const PathVelEntry&) const = 0;
    // Return the contact point stations in the body frame, that is, vectors
    // from body origin Bo to point P and Q, expressed in B.
    void getContactStations(const State&            state, 
                            const PathInstanceInfo& instInfo,
                            const PathPosEntry&     posEntry,
                            Vec3& P_B, Vec3& Q_B) const 
    {   Vec3 P_S, Q_S;
        getContactPointsOnObstacle(state, instInfo, posEntry, P_S, Q_S);
        const Transform& X_BS = getObstaclePoseOnBody(state, instInfo);
        P_B = X_BS*P_S; Q_B = X_BS*Q_S; // 36 flops
    }

    // Return the contact point stations in the body frame, but re-expressed
    // in Ground. This is still the vector from the body origin Bo to the
    // contact points P and G; it is not measured from the ground origin.
    // 66 flops.
    void expressContactStationsInGround(const State&            state, 
                                const PathInstanceInfo& instInfo,
                                const PathPosEntry&     posEntry,
                                Vec3& PB_G, Vec3& QB_G) const 
    {   Vec3 P_B, Q_B;
        getContactStations(state, instInfo, posEntry, P_B, Q_B);
        const Rotation& R_GB = getBodyTransform(state).R();
        PB_G = R_GB*P_B; QB_G = R_GB*Q_B;   // 30 flops
    }

    const MobilizedBody& getMobilizedBody() const {return mobod;}
    const Transform& getDefaultPoseOnBody() const {return defaultX_BS;}
    const Transform& getObstaclePoseOnBody(const State& state,
                                           const PathInstanceInfo& instInfo) const
    {   return instInfo.obstaclePose[index]; }
    const Transform& getBodyTransform(const State& state) const
    {   return mobod.getBodyTransform(state); }
    const SpatialVec& getBodyVelocity(const State& state) const
    {   return mobod.getBodyVelocity(state); }
    void setCablePath(CablePath& path) { cablePath = &path; }
    void setCableObstacleIndex(int ix) { index = CableObstacleIndex(ix); }
    void invalidateTopology(); // see below   

protected:
friend class CableObstacle;

    // CablePath to which this obstacle belongs, and the index within that path.
    CablePath*              cablePath;
    CableObstacleIndex      index;

    // MobilizedBody to which this obstacle is fixed.
    const MobilizedBody     mobod;
    Transform               defaultX_BS;
    DecorativeGeometry      decoration;

    mutable int             referenceCount;
};



//==============================================================================
//                         CABLE PATH :: IMPL
//==============================================================================


// This is the internal implementation class for CablePath. It is reference
// counted so that CablePath objects can be multiply referenced.
class CablePath::Impl {
public:
    Impl() : referenceCount(0), cables(0) {}
    Impl(CableTrackerSubsystem& cables)
    :   referenceCount(0), cables(&cables) {}
    ~Impl() {assert(referenceCount == 0);}

    int getNumObstacles() const {return obstacles.size();}
    const CableObstacle& getObstacle(CableObstacleIndex ix) const
    {   return obstacles[ix]; }
    const CableObstacle::Impl& getObstacleImpl(CableObstacleIndex ix) const
    {   return obstacles[ix].getImpl(); }

    // Insert an obstacle into the path, just before the termination point
    // which is always the last entry. The first two are the origin point
    // and then the termination point. 
    CableObstacleIndex adoptObstacle(CableObstacle& obstacle) 
    {   if (obstacles.size() < 2) {
            obstacles.push_back(obstacle);
            return CableObstacleIndex(obstacles.size()-1);
        }          
        obstacles.insert(obstacles.end()-1, obstacle);
        // Update the relocated termination point's obstacle index.
        obstacles.back().updImpl().setCableObstacleIndex(obstacles.size()-1);
        return CableObstacleIndex(obstacles.size()-2); 
    }

    Real getCableLength(const State& state) const {
        const PathPosEntry& posEntry = getPosEntry(state);
        return posEntry.length;
    }

    Real getCableLengthDot(const State& state) const {
        const PathVelEntry& velEntry = getVelEntry(state);
        return velEntry.lengthDot;
    }

    void applyBodyForces(const State& state, Real tension, 
                         Vector_<SpatialVec>& bodyForcesInG) const
    {
        if (tension <= 0) return;

        const PathPosEntry& posEntry = getPosEntry(state);
        for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
            const CableObstacle::Impl& obs = obstacles[ox].getImpl();
            const MobilizedBody& B = obs.getMobilizedBody();
            const SpatialVec& F_GB = posEntry.obstacleUnitForces[ox];
            B.applyBodyForce(state, tension*F_GB, bodyForcesInG);
        }
    }

    Real calcCablePower(const State& state, Real tension) const {
        if (tension <= 0) return 0;
        const PathVelEntry& velEntry = getVelEntry(state);
        return tension * velEntry.unitPower;
    }

    void realizeTopology(State& state) {
        // Initialize instance info from defaults.
        PathInstanceInfo init;
        for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
            const CableObstacle::Impl& obs = obstacles[ox].getImpl();
            init.obstacleDisabled.push_back(false);
            init.obstaclePose.push_back(obs.getDefaultPoseOnBody());
        }
        cout << init << endl;
        // Allocate and initialize instance state variable.
        instanceInfoIx = cables->allocateDiscreteVariable(state,
            Stage::Instance, new Value<PathInstanceInfo>(init));

        // Allocate cache entries for position and velocity calculations.
        posEntryIx = cables->allocateLazyCacheEntry(state, Stage::Position,
            new Value<PathPosEntry>());
        velEntryIx = cables->allocateLazyCacheEntry(state, Stage::Velocity,
            new Value<PathVelEntry>());
    }

    void realizePosition(const State& state) const {
        ensurePositionKinematicsCalculated(state);
    }

    void realizeVelocity(const State& state) const {
        ensureVelocityKinematicsCalculated(state);
    }

    void ensurePositionKinematicsCalculated(const State& state) const {
        if (cables->isCacheValueRealized(state, posEntryIx))
            return;

        const PathInstanceInfo& instInfo = getInstanceInfo(state);
        PathPosEntry&           posEntry = updPosEntry(state);

        posEntry.length = 0;
        posEntry.straightDirections.clear();
        posEntry.obstacleUnitForces.clear();

        Vec3 prevQB_G, prevQ_G;
        for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
            const CableObstacle::Impl& obs = obstacles[ox].getImpl();
            Vec3 PB_G, QB_G;
            obs.expressContactStationsInGround(state,instInfo,posEntry,
                                               PB_G, QB_G);
            const Vec3& Bo_G = obs.getBodyTransform(state).p();
            const Vec3 P_G = Bo_G + PB_G;
            if (ox == 0) {
                prevQB_G = QB_G;
                prevQ_G  = Bo_G + QB_G;
                posEntry.obstacleUnitForces.push_back(SpatialVec(Vec3(0)));
                continue;
            }

            const Vec3 lineSeg(P_G - prevQ_G);
            const Real lineSegLen = lineSeg.norm();
            const UnitVec3 lineDir(lineSeg/lineSegLen, true);
            posEntry.straightDirections.push_back(lineDir);
            posEntry.length += lineSegLen;
            const SpatialVec unitForcePrevQ(prevQB_G % lineDir, lineDir);
            const SpatialVec unitForceP(lineDir % PB_G, -lineDir);
            posEntry.obstacleUnitForces.back() += unitForcePrevQ;
            posEntry.obstacleUnitForces.push_back(unitForceP);

            const Real segLen = obs.getSegmentLength(state, instInfo, posEntry);
            posEntry.length += segLen;
            prevQB_G = QB_G;
            prevQ_G  = Bo_G + QB_G;
        }

        cables->markCacheValueRealized(state, posEntryIx);
    }

    void ensureVelocityKinematicsCalculated(const State& state) const {
        if (cables->isCacheValueRealized(state, velEntryIx))
            return;

        const PathInstanceInfo& instInfo = getInstanceInfo(state);
        const PathPosEntry& posEntry = getPosEntry(state);
        PathVelEntry&       velEntry = updVelEntry(state);

        velEntry.lengthDot = 0;
        velEntry.unitPower = 0;

        //TODO: calc length dot
        Vec3 vprev_GQ;
        for (CableObstacleIndex ox(0); ox < obstacles.size(); ++ox) {
            const CableObstacle::Impl& obs = obstacles[ox].getImpl();
            const MobilizedBody& B = obs.getMobilizedBody();
            Vec3 P_B, Q_B, v_GP, v_GQ;
            obs.getContactStations(state, instInfo, posEntry, P_B, Q_B);
            v_GP = B.findStationVelocityInGround(state, P_B);
            v_GQ = B.findStationVelocityInGround(state, Q_B);

            if (ox > 0) {
                CableObstacleIndex prevx(ox-1);
                velEntry.lengthDot += dot(v_GP - vprev_GQ,
                                          posEntry.straightDirections[prevx]);
            }
            vprev_GQ = v_GQ;

            const SpatialVec& V_GB = obs.getBodyVelocity(state);
            velEntry.unitPower += 
                ~posEntry.obstacleUnitForces[ox] * V_GB;
        }

        cables->markCacheValueRealized(state, velEntryIx);
    }

    const PathInstanceInfo& getInstanceInfo(const State& state) const 
    {   return Value<PathInstanceInfo>::downcast
           (cables->getDiscreteVariable(state, instanceInfoIx)); }
    PathInstanceInfo& updInstanceInfo(State& state) const 
    {   return Value<PathInstanceInfo>::updDowncast
           (cables->updDiscreteVariable(state, instanceInfoIx)); }

    // Lazy-evaluate position kinematics and return the result.
    const PathPosEntry& getPosEntry(const State& state) const 
    {   ensurePositionKinematicsCalculated(state);
        return Value<PathPosEntry>::downcast
           (cables->getCacheEntry(state, posEntryIx)); }

    PathPosEntry& updPosEntry(const State& state) const 
    {   return Value<PathPosEntry>::updDowncast
           (cables->updCacheEntry(state, posEntryIx)); }

    // Lazy-evaluate velocity kinematics and return the result.
    const PathVelEntry& getVelEntry(const State& state) const 
    {   ensureVelocityKinematicsCalculated(state);
        return Value<PathVelEntry>::downcast
           (cables->getCacheEntry(state, velEntryIx)); }

    PathVelEntry& updVelEntry(const State& state) const 
    {   return Value<PathVelEntry>::updDowncast
           (cables->updCacheEntry(state, velEntryIx)); }

    void invalidateTopology()
    {   if (cables) cables->invalidateSubsystemTopologyCache(); }
private:
friend class CablePath;

    // Subsystem to which this path belongs, and the index within that
    // subsystem.
    CableTrackerSubsystem*  cables;
    CablePathIndex          index;

    // The list of via points and surfaces, ordered by expected path 
    // coordinate. The first obstacle is the origin point; the last is the
    // termination point.
    Array_<CableObstacle>   obstacles;

    // TOPOLOGY CACHE (set during realizeTopology())
    DiscreteVariableIndex   instanceInfoIx;
    CacheEntryIndex         posEntryIx;
    CacheEntryIndex         velEntryIx;

    mutable int             referenceCount;
};



//==============================================================================
//                     CABLE OBSTACLE :: VIA POINT :: IMPL
//==============================================================================
class CableObstacle::ViaPoint::Impl : public CableObstacle::Impl {
    typedef CableObstacle::Impl Super;
public:
    Impl(CablePath& path, const MobilizedBody& viaMobod,
         const Vec3& station) 
    :   Super(path, viaMobod, station) 
    {   decoration = DecorativePoint().setColor(Red); }

    void getContactPointsOnObstacle(const State& state, 
                                    const PathInstanceInfo& instInfo,
                                    const PathPosEntry& posEntry,
                                    Vec3& P_S, Vec3& Q_S) const OVERRIDE_11
    {
        const Transform& X_BS = getObstaclePoseOnBody(state, instInfo);
        P_S = Q_S = X_BS.p();
    }

    Real getSegmentLength(const State& state,
                          const PathInstanceInfo& instInfo,
                          const PathPosEntry& posInfo) const OVERRIDE_11
    {   return 0; } // via point has no segment

    Real getSegmentLengthDot( const State& state,
                              const PathInstanceInfo& instInfo,
                              const PathPosEntry& posInfo,
                              const PathVelEntry& velInfo) const OVERRIDE_11
    {   return 0; } // via point has no segment
};

//==============================================================================
//                     CABLE OBSTACLE :: SURFACE :: IMPL
//==============================================================================
class CableObstacle::Surface::Impl : public CableObstacle::Impl {
    typedef CableObstacle::Impl Super;
public:
    Impl(CablePath& path, const MobilizedBody& viaMobod,
         const Transform& pose) 
    :   Super(path, viaMobod, pose) {}

    void getContactPointsOnObstacle(const State& state, 
                                    const PathInstanceInfo& instInfo,
                                    const PathPosEntry& posEntry,
                                    Vec3& P_S, Vec3& Q_S) const OVERRIDE_11
    {
        // TODO: get real points from posEntry
        const Transform& X_BS = getObstaclePoseOnBody(state, instInfo);
        P_S = Q_S = X_BS.p();
    }

    Real getSegmentLength(const State& state,
                          const PathInstanceInfo& instInfo,
                          const PathPosEntry& posInfo) const OVERRIDE_11
    {   return 0; } // TODO

    Real getSegmentLengthDot( const State& state,
                              const PathInstanceInfo& instInfo,
                              const PathPosEntry& posInfo,
                              const PathVelEntry& velInfo) const OVERRIDE_11
    {   return 0; } // TODO
 
private:
friend class CableObstacle::Surface;

    ContactGeometry     surface;
    Vec3                pathPreferencePointInS; 
};


} // namespace SimTK


#endif // SimTK_SIMBODY_CABLE_PATH_IMPL_H_

