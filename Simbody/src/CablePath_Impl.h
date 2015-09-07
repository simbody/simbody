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
 * Authors: Michael Sherman, Ian Stavness                                     *
 * Contributors: Andreas Scholz                                               *
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

// This file declares the private implementation classes for CablePath and
// CableObstacle classes. This is internal source, not part of the Simbody API.

#include "SimTKmath.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/CableTrackerSubsystem.h"
#include "simbody/internal/CablePath.h"


#include <cassert>
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

// This is a unique integer type for identifying active obstacles in a
// particular path, including both via points and surfaces.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ActiveObstacleIndex);

// This is a unique integer type for identifying those obstacles that 
// are surfaces rather than via points, regardless of whether they are 
// enabled or active.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(SurfaceObstacleIndex);

// These are indices assigned to the subset of active obstacles that are 
// surfaces and thus have geodesics and contribute unknowns to the path
// finding problem. This set is the intersection of active & surface.
SimTK_DEFINE_UNIQUE_INDEX_TYPE(ActiveSurfaceIndex);


// This is a discrete state variable for holding a path's instance-level 
// information, including placement of points and surfaces on their bodies
// and obstacle enable/disable settings.
class PathInstanceInfo {
public:
    // Initialize instance info from the defaults built into the obstacles.
    explicit PathInstanceInfo
       (const Array_<CableObstacle,CableObstacleIndex>& obstacles);

    int getNumObstacles()        const {return mapObstacleToSurface.size();}
    int getNumSurfaceObstacles() const {return mapSurfaceToObstacle.size();}

    // Each surface obstacle gets one of these. Via point obstacles will
    // have a surface obstacle index that tests !isValid(). This is filled
    // in by the constructor.
    Array_<SurfaceObstacleIndex,CableObstacleIndex>  mapObstacleToSurface;
    // This is indexed by surface obstacle index to give back the corresponding
    // obstacle index.
    Array_<CableObstacleIndex,SurfaceObstacleIndex>  mapSurfaceToObstacle;

    // One entry per obstacle but you can't disable the origin or termination
    // points.
    Array_<bool,CableObstacleIndex>         obstacleDisabled;
    Array_<Transform,CableObstacleIndex>    obstaclePose; // X_BS[i]
};


// This is a cache entry for holding a path's calculated position-level 
// information. At the time it is created we know the total number of
// obstacles n (including the end points), but not which ones are active. We 
// don't know how many 
// unknowns there are until we find out later based on the subset of active
// obstacles that have unknown P's and Q's.
// Note that disabled obstacles still get slots here but we never look at 
// them. And disabled != inactive.
//
// Individual obstacles can only know their CableObstacleIndex because all
// possible obstacles for a cable are provided during construction (i.e., that
// is topological information). Then when we determine the set of active 
// obstacles, we assign ActiveObstacleIndex and ActiveSurfaceIndex indices
// which must be maintained in the state here. An obstacle can use its
// CableObstacleIndex to obtain the other indices for access to related
// quantities.
class PathPosEntry {
public:
    PathPosEntry() : length(NaN) {}

    // Set the number of obstacles to n. If there is any information already
    // in this object, it is lost.
    void setNumObstacles(int n, int nSurfaces) {
        mapToActive.clear();        mapToActive.resize(n);
        mapToActiveSurface.clear(); mapToActiveSurface.resize(n);
        mapToCoords.clear();        mapToCoords.resize(n);
        witnesses.clear(); witnesses.resize(nSurfaces, Real(NaN));
        closestSurfacePoint.clear(); closestSurfacePoint.resize(nSurfaces, Vec3(NaN));
        closestPathPoint.clear(); closestPathPoint.resize(nSurfaces, Vec3(NaN));
        initialize(0,0,0); // not yet
    }

    // mapToActive, mapToActiveSurface and mapToCoords are already allocated.
    void initialize(int na, int nas, int nx) {
        length = NaN;
        // Active obstacles
        eIn_G.clear(); eIn_G.resize(na); // all NaN
        Fu_GB.clear(); Fu_GB.resize(na, SpatialVec(Vec3(NaN)));
        // Active surfaces
        geodesics.clear(); geodesics.resize(nas);
        x.clear(); x.resize(nx);
        err.clear(); err.resize(nx);
        J.clear(); J.resize(nx,nx);
    }

    // Return the obstacle index of the first active obstacle following the
    // given obstacle. Returns an invalid index if there are no more.
    CableObstacleIndex findNextActiveObstacle(CableObstacleIndex thisOx) const {
        for (CableObstacleIndex ox(thisOx+1); ox < mapToActive.size(); ++ox)
            if (mapToActive[ox].isValid())
                return ox;
        return CableObstacleIndex();
    }

    // Return the obstacle index of the first active obstacle preceding the
    // given obstacle. Returns an invalid index if there are no more in that
    // direction.
    CableObstacleIndex findPrevActiveObstacle(CableObstacleIndex thisOx) const {
        for (CableObstacleIndex ox(thisOx); ox > 0; --ox)
            if (mapToActive[ox.prev()].isValid())
                return ox.prev();
        return CableObstacleIndex();
    }

    // This is the total length of the path corresponding to the current
    // configuration and values for contact point coordinates x stored here.
    Real length;

    //                         ALL OBSTACLES

    // Map each cable obstacle to its ActiveObstacleIndex if it is currently
    // active, otherwise the entry tests !isValid(). Disabled obstacles are
    // always inactive, enabled via points are always active. The origin and 
    // termination points are always active, so mapToActive[0]=0 and 
    // mapToActive[n-1]=nActive-1.
    Array_<ActiveObstacleIndex,CableObstacleIndex>  mapToActive;

    // Each active, surface obstacle gets one of these. Points and inactive
    // surfaces will have indices that test !isValid().
    Array_<ActiveSurfaceIndex,CableObstacleIndex>   mapToActiveSurface;

    // Each active, surface obstacle is assigned a slot in the vector of 
    // unknowns x. The first index of that slot is saved here, with a -1
    // for any obstacle that is not a surface or is inactive.
    Array_<int,CableObstacleIndex>                  mapToCoords;

    //                       ACTIVE OBSTACLES

    // These are unit vectors along the straight segments following the Q
    // point of each active obstacle except the termination point.
    Array_<UnitVec3,ActiveObstacleIndex>            eIn_G;

    // Multiply these unit spatial forces by tension to get body spatial
    // forces in Ground.
    Array_<SpatialVec,ActiveObstacleIndex>          Fu_GB;


    //                       SURFACE OBSTACLES
    // Each surface obstacle has an event witness (trigger) function that
    // passes through zero when the cable lifts off or touches down on this
    // obstacle.
    Array_<Real,SurfaceObstacleIndex>               witnesses;
    // Each inactive surface tracks the closest point on the surface to the
    // path line segment that it would contact if it were active. We also
    // keep the closest point on the path for display purposes. Both points
    // are in the surface's local frame S. Only the
    // slots for inactive surfaces have meaningful content here.
    Array_<Vec3,SurfaceObstacleIndex>               closestSurfacePoint;
    Array_<Vec3,SurfaceObstacleIndex>               closestPathPoint;


    //                       ACTIVE SURFACES

    // We calculate a geodesic for each active surface obstacle.
    Array_<Geodesic,ActiveSurfaceIndex>             geodesics;

    // This is a dense vector with all the active surface obstacles' unknowns 
    // packed together.
    Vector      x;      // nx unknowns

    // This is the patherr corresponding to x and is always the same length.
    Vector      err;    // patherr (nx of these)

    // This is J(x) where J=partial(patherr)/partial(x), and its LU 
    // factorization for use in solving for length dot at Velocity stage.
    // TODO: this is banded but we're treating it as full.
    Matrix      J;      // nx X nx
    FactorLU    JInv;
};


// This is a cache entry for holding a path's calculated velocity-level 
// information.
class PathVelEntry {
public:
    PathVelEntry() : lengthDot(NaN), unitPower(NaN) {}

    // Set the number of obstacles n (including end points). We don't yet know
    // the number of active obstacles. Any information previously here is
    // lost.
    void setNumObstacles(int n, int nSurfaces) {
        initialize(0,0,0); // not yet
    }

    void initialize(int na, int nas, int nx) {
        xdot.clear(); xdot.resize(nx);
        nerrdotK.clear(); nerrdotK.resize(nx);
        lengthDot = NaN;
        unitPower = NaN;
    }

    Vector nerrdotK; // negative of frozen-contact point portion of patherrdot
    Vector xdot;    // calculated time derivatives of x; xdot=-J\errdotK
    Real lengthDot;
    Real unitPower; // unitForces*body velocities
};

std::ostream& operator<<(std::ostream& o, const PathInstanceInfo& info);
std::ostream& operator<<(std::ostream& o, const PathPosEntry& entry);
std::ostream& operator<<(std::ostream& o, const PathVelEntry& entry);



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
        defaultX_BS(defaultPose), defaultDisabled(false) {}

    virtual ~Impl() {assert(referenceCount == 0);}

    // Return the number of scalar coordinates needed to represent the 
    // location of one of this obstacle's contact points. This will be zero for
    // via points, 2 for parametric surfaces, and 3 for
    // implicit surfaces. This is also necessarily half the number of 
    // coordinates, their time derivatives, and the number of path error terms 
    // associated with this obstacle.
    virtual int getNumCoordsPerContactPoint() const = 0;

    // Concrete CableObstacle classes must implement these.

    // Given coordinates xi=(xP,xQ) for this obstacle (taken from the supplied
    // PathPosEntry), return the Cartesian positions of P and Q measured and
    // expressed in the obstacle's S frame. For via points these are always
    // at (0,0,0)_S. For implicit surfaces, the returned points will *not*
    // necessarily be on the surface; that is, we are not doing any projection
    // of xP and xQ here.
    virtual void getContactPointsOnObstacle(const State&, 
                                        const PathInstanceInfo&,
                                        const PathPosEntry&,
                                        Vec3& P_S, Vec3& Q_S) const = 0;
    // After geodesics have been calculated, this will return the length
    // of the geodesic from P' to Q' for this obstacle, where P' and Q' are
    // the nearest points on the surface corresponding to given points P and
    // Q which for an implicit surface may be off the surface. For via points
    // the length will be returned zero since we consider P and Q to be in the 
    // same place.
    virtual Real getSegmentLength(const State&, 
                                  const PathInstanceInfo&,
                                  const PathPosEntry&) const = 0;

    // Once xdots have been calculated, this will return the rate of length
    // change for the geodesic on this obstacle, or zero for via points. Note
    // that the geodesic connects P' to Q'; see above for what that means.
    virtual Real getSegmentLengthDot(const State&, 
                                     const PathInstanceInfo&,
                                     const PathPosEntry&,
                                     const PathVelEntry&) const = 0;

    // Return the contact point stations in the body frame, that is, vectors
    // from body origin Bo to point P and Q, expressed in B. 
    // Cost is 36 flops.
    // TODO: should this be projected onto the surface? 
    void getContactStationsOnBody(const State&            state, 
                                  const PathInstanceInfo& instInfo,
                                  const PathPosEntry&     ppe,
                                  Vec3& P_B, Vec3& Q_B) const 
    {   Vec3 P_S, Q_S;
        getContactPointsOnObstacle(state, instInfo, ppe, P_S, Q_S);
        const Transform& X_BS = getObstaclePoseOnBody(state, instInfo);
        P_B = X_BS*P_S; Q_B = X_BS*Q_S; // 36 flops
    }

    // Return the contact point stations measured in the body frame, but 
    // re-expressed in Ground. This is still the vector from the body origin 
    // Bo to the contact points P and G; it is not measured from the ground 
    // origin. Cost is 66 flops.
    // TODO: should this be projected onto the surface? 
    void expressContactStationsInGround(const State&            state, 
                                const PathInstanceInfo& instInfo,
                                const PathPosEntry&     ppe,
                                Vec3& PB_G, Vec3& QB_G) const 
    {   Vec3 P_B, Q_B;
        getContactStationsOnBody(state, instInfo, ppe, P_B, Q_B); // 36 flops
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
    void setDisabledByDefault(bool shouldBeDisabled) 
    {   defaultDisabled=shouldBeDisabled; }
    bool isDisabledByDefault() const {return defaultDisabled;}
    void invalidateTopology(); // see below 

    const DecorativeGeometry& getDecoration() const {return decoration;}
    DecorativeGeometry& updDecoration() {return decoration;}
    void setDecoration(const DecorativeGeometry& newGeo)
    {   decoration = newGeo; }

protected:
friend class CableObstacle;

    // CablePath to which this obstacle belongs, and the index within that path.
    CablePath*              cablePath;
    CableObstacleIndex      index;

    // MobilizedBody to which this obstacle is fixed.
    const MobilizedBody     mobod;
    // Pose of surface frame (or point) on its body.
    Transform               defaultX_BS;
    // How to draw this obstacle in the visualizer.
    DecorativeGeometry      decoration;

    // Whether this obstacle should start out disabled.
    bool                    defaultDisabled;

    // The number of handles that reference this implementation object.
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

    void solveForInitialCablePath(State& state) const;

    Real getCableLength(const State& state) const {
        const PathPosEntry& posEntry = getPosEntry(state);
        return posEntry.length;
    }

    Real getCableLengthDot(const State& state) const {
        const PathVelEntry& velEntry = getVelEntry(state);
        return velEntry.lengthDot;
    }

    void applyBodyForces(const State& state, Real tension, 
                         Vector_<SpatialVec>& bodyForcesInG) const;

    // TODO: this isn't right when the tension goes to zero during rapid
    // shrinking due to dissipation cancelling stiffness.
    Real calcCablePower(const State& state, Real tension) const {
        if (tension <= 0) return 0;
        const PathVelEntry& velEntry = getVelEntry(state);
        return tension * velEntry.unitPower;
    }

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);

    // We know which obstacles are disabled. Assume the rest are active
    void realizeInstance(const State& state) const;

    // Solve for the new path and its length L, evaluate event witnesses to
    // monitor for liftoff or touchdown.
    void realizePosition(const State& state) const;

    // Evaluate the path length rate of change (Ldot).
    void realizeVelocity(const State& state) const;

    void realizeAcceleration(const State& state) const;

    void calcEventTriggerInfo
       (const State&, Array_<EventTriggerInfo>&) const;

    void handleEvents
       (State&, Event::Cause, const Array_<EventId>& eventIds,
        const HandleEventsOptions& options, HandleEventsResults& results) const;

    // Update the cable path and its length in the state cache. This is the 
    // expensive part of the cable path computation. State must already have
    // been realized to position stage. This is invoked automatically when
    // you request access to the position cache entry.
    void ensurePositionKinematicsCalculated(const State& state) const;

    // Calculate the cable's rate of length change and put the result in
    // the state cache. State must already have been realized to Velocity stage.
    // This is invoked automatically when you request access to the velocity
    // cache entry.
    void ensureVelocityKinematicsCalculated(const State& state) const;

    // Methods for convenient access to state variables and cache entries.

    const PathInstanceInfo& getInstanceInfo(const State& state) const 
    {   return Value<PathInstanceInfo>::downcast
           (cables->getDiscreteVariable(state, instanceInfoIx)); }
    PathInstanceInfo& updInstanceInfo(State& state) const 
    {   return Value<PathInstanceInfo>::updDowncast
           (cables->updDiscreteVariable(state, instanceInfoIx)); }

    Real getIntegratedLengthDot(const State& state) const
    {   return cables->getZ(state)[integratedLengthDotIx]; }

    void setIntegratedLengthDot(State& state, Real value) const 
    {   cables->updZ(state)[integratedLengthDotIx] = value; }

    // Return the *previous* path pos entry.
    const PathPosEntry& getPrevPosEntry(const State& state) const 
    {   return Value<PathPosEntry>::downcast
           (cables->getDiscreteVariable(state, posEntryIx)); }
    PathPosEntry& updPrevPosEntry(State& state) const 
    {   return Value<PathPosEntry>::updDowncast
           (cables->updDiscreteVariable(state, posEntryIx)); }
    // Return the *previous* path vel entry.
    const PathVelEntry& getPrevVelEntry(const State& state) const 
    {   return Value<PathVelEntry>::downcast
           (cables->getDiscreteVariable(state, velEntryIx)); }
    PathVelEntry& updPrevVelEntry(State& state) const 
    {   return Value<PathVelEntry>::updDowncast
           (cables->updDiscreteVariable(state, velEntryIx)); }

    // Lazy-evaluate position kinematics and return the result.
    const PathPosEntry& getPosEntry(const State& state) const 
    {   ensurePositionKinematicsCalculated(state);
        return Value<PathPosEntry>::downcast
           (cables->getDiscreteVarUpdateValue(state, posEntryIx)); }

    PathPosEntry& updPosEntry(const State& state) const 
    {   return Value<PathPosEntry>::updDowncast
           (cables->updDiscreteVarUpdateValue(state, posEntryIx)); }

    // Lazy-evaluate velocity kinematics and return the result.
    const PathVelEntry& getVelEntry(const State& state) const 
    {   ensureVelocityKinematicsCalculated(state);
        return Value<PathVelEntry>::downcast
           (cables->getDiscreteVarUpdateValue(state, velEntryIx)); }

    PathVelEntry& updVelEntry(const State& state) const 
    {   return Value<PathVelEntry>::updDowncast
           (cables->updDiscreteVarUpdateValue(state, velEntryIx)); }

    // Be sure to call this whenever you make a topology-level change to
    // the cable definition, like adding an obstacle or modifying one in
    // a significant way.
    void invalidateTopology()
    {   if (cables) cables->invalidateSubsystemTopologyCache(); }

    // Given kinematics and a set of contact point coordinates x (already in
    // PathPosEntry), calculate the resulting path errors and related 
    // quantities with the result going back into PathPosEntry.
    void calcPathError
       (const State&, const PathInstanceInfo&, PathPosEntry&) const;

    // Given kinematics K and a set of contact point coordinates x (already in
    // PathPosEntry), calculate the Jacobian D patherr(K;x)/Dx, with the
    // result going back into PathPosEntry. This is a banded matrix assembled
    // from individual blocks provided by the active surface obstacles.
    void calcPathErrorJacobian
       (const State&, const PathInstanceInfo&, PathPosEntry&) const;

private:
friend class CablePath;

    // Starting with the current guess in PathPosEntry, project each contact
    // point to the nearest point on the surface.
    void projectOntoSurface(const PathInstanceInfo&, PathPosEntry&) const;

    // Starting with an initial guess for the contact coordinates x, make
    // repeated calls to calcPathError to drive the path errors to zero by
    // modifying x in PathPosEntry.
    void solveForPathPoints
       (const State&, const PathInstanceInfo&, PathPosEntry&) const;

    // Given a solved path, compute the term of the path error time derivative
    // due only to kinematics, with path points frozen on their surfaces.
    // Results go into the errdotK member of the velocity cache entry. This
    // is used to solve the for xdot such that the total path error time
    // derivative is zero.
    void findKinematicVelocityErrors
        (const State&, const PathInstanceInfo&, const PathPosEntry&,
         PathVelEntry&) const;

    // Call this for any obstacle except the end points to obtain the 
    // segment between the previous and next active obstacles which is the
    // part with which this obstacle may interact. The points are returned
    // in the local frame S of the given obstacle.
    void findPathSegmentForObstacle
       (const State&, const PathInstanceInfo&, const PathPosEntry&,
        CableObstacleIndex ox, Vec3& Qprev_S, Vec3& Pnext_S) const;


    // Subsystem to which this path belongs, and the index within that
    // subsystem.
    CableTrackerSubsystem*  cables;
    CablePathIndex          index;

    // The list of via points and surfaces, ordered by expected path 
    // coordinate. The first obstacle is the origin point; the last is the
    // termination point.
    Array_<CableObstacle,CableObstacleIndex>   obstacles;

    // TOPOLOGY CACHE (set during realizeTopology())
    DiscreteVariableIndex       instanceInfoIx;
    ZIndex                      integratedLengthDotIx;
    DiscreteVariableIndex       posEntryIx; // these are auto-update
    DiscreteVariableIndex       velEntryIx;
    EventTriggerByStageIndex    eventIx; // 1st index; one for each surface
    std::map<EventId, CableObstacleIndex> mapEventIdToObstacle;

    mutable int                 referenceCount;
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

    int getNumCoordsPerContactPoint() const override {return 0;}

    void getContactPointsOnObstacle(const State& state, 
                                    const PathInstanceInfo& instInfo,
                                    const PathPosEntry& posEntry,
                                    Vec3& P_S, Vec3& Q_S) const override
    {
        P_S = Q_S = Vec3(0);
    }

    Real getSegmentLength(const State& state,
                          const PathInstanceInfo& instInfo,
                          const PathPosEntry& posInfo) const override
    {   return 0; } // via point has no segment

    Real getSegmentLengthDot( const State& state,
                              const PathInstanceInfo& instInfo,
                              const PathPosEntry& posInfo,
                              const PathVelEntry& velInfo) const override
    {   return 0; } // via point has no segment
};



//==============================================================================
//                     CABLE OBSTACLE :: SURFACE :: IMPL
//==============================================================================
class CableObstacle::Surface::Impl : public CableObstacle::Impl {
    typedef CableObstacle::Impl Super;
public:
    Impl(CablePath& path, const MobilizedBody& mobod,
         const Transform& pose, const ContactGeometry& geom) 
    :   Super(path, mobod, pose), surface(geom), 
        nearPointInS(NaN), xPhint(NaN), xQhint(NaN) 
    {   decoration = geom.createDecorativeGeometry()
            .setColor(Orange).setOpacity(.75).setResolution(3); }

    const ContactGeometry& getContactGeometry() const {return surface;}

    // Hardcoded for implicit surfaces -- would be 2 for parametric.
    int getNumCoordsPerContactPoint() const override {return 3;}

    void getContactPointsOnObstacle(const State& state, 
                                    const PathInstanceInfo& instInfo,
                                    const PathPosEntry& ppe,
                                    Vec3& P_S, Vec3& Q_S) const override
    {
        const ActiveSurfaceIndex asx   = ppe.mapToActiveSurface[index];
        const int                xSlot = ppe.mapToCoords[index];
        assert(asx.isValid() && xSlot >= 0);
        P_S = Vec3::getAs(&ppe.x[xSlot]);  // implicit coords are point in S
        Q_S = Vec3::getAs(&ppe.x[xSlot+3]);
    }

    Real getSegmentLength(const State& state,
                          const PathInstanceInfo& instInfo,
                          const PathPosEntry& ppe) const override
    {   const ActiveSurfaceIndex asx = ppe.mapToActiveSurface[index];
        assert(asx.isValid());
        const Geodesic& geod = ppe.geodesics[asx];
        return geod.getLength();
    }


    Real getSegmentLengthDot( const State& state,
                              const PathInstanceInfo& instInfo,
                              const PathPosEntry& ppe,
                              const PathVelEntry& pve) const override
    {   const ActiveSurfaceIndex asx   = ppe.mapToActiveSurface[index];
        const int                xSlot = ppe.mapToCoords[index];
        assert(asx.isValid() && xSlot >= 0);
        // contact point velocities
        const Vec3 xdotP = Vec3::getAs(&pve.xdot[xSlot]); 
        const Vec3 xdotQ = Vec3::getAs(&pve.xdot[xSlot+3]);
        const Geodesic& geod = ppe.geodesics[asx];
        return geod.calcLengthDot(xdotP,xdotQ);
    }

    // During initialization, call this to calculate an appropriate geodesic
    // connecting points P' and Q', the nearest surface points to points P and
    // Q whose surface coordinates are supplied.
    // There may be multiple geodesics between P' and Q'. If a near point
    // N has been provided, we'll return the geodesic that passes closest to it,
    // otherwise we'll return the shortest one. If entry and exit hints are 
    // provided they may be used to bias the algorithm. 
    void initSurfacePath(const Vec3&        xP, 
                         const Vec3&        xQ,
                         const UnitVec3&    entryHint_S, // optional (use NaN)
                         const UnitVec3&    exitHint_S,  // optional (use NaN)
                         Geodesic&          path) const;

    // Calculate the surface-local path error due to the given entry/exit
    // points and directions being inconsistent with a geodesic connecting
    // those points.
    Vec6 calcSurfacePathError(  const Geodesic& previous,
                                const UnitVec3& entryDir_S,
                                const Vec3&     xP,
                                const Vec3&     xQ,
                                const UnitVec3& exitDir_S,
                                Geodesic&       next) const;

    // Calculate the 6-vector surface-local path error due to the given 
    // entry/exit points and directions being inconsistent with a geodesic 
    // connecting those points, and the 6x12 Jacobian of the errors with respect 
    // to each of the 12 measure numbers of the 4 input vectors.
    void calcSurfacePathErrorJacobianAnalytically
       (const UnitVec3& entryDir_S,
        const Vec3&     xP_S,
        const Vec3&     xQ_S,
        const UnitVec3& exitDir_S,
        const Geodesic& current,    // Must have used above parameters.
        Mat63&          DerrDentry, // 4x3 for parametric
        Mat63&          DerrDxP,    // 4x2       "
        Mat63&          DerrDxQ,    // 4x2       "
        Mat63&          DerrDexit)  // 4x3       "
        const;


    // Calculate partial derivatives of the surface path error with respect to
    // each of its four arguments.
    void calcSurfacePathErrorJacobianNumerically
       (const Geodesic& previous,
        const UnitVec3& entryDir_S,
        const Vec3&     xP,
        const Vec3&     xQ,
        const UnitVec3& exitDir_S,
        const Geodesic& current,    // Must have used above parameters.
        Mat63&          DerrDentry, // 4x3 for parametric
        Mat63&          DerrDxP,    // 4x2       "
        Mat63&          DerrDxQ,    // 4x2       "
        Mat63&          DerrDexit)  // 4x3       "
        const;

    // Calculate the *negated* surface-local contribution to the path error time 
    // derivative due only to material point velocities, with the surface
    // coordinates xP and xQ held constant at the given values.
    // This is nerrdotK = -(Dpatherr(K;x)/DK)*(dK/dt).
    Vec6 calcSurfaceNegKinematicVelocityError
       (  const Geodesic& geodesic,
          const UnitVec3& entryDir_S,
          const Vec3&     entryDirDot_S,    // time derivative in S
          const Vec3&     xP,
          const Vec3&     xQ,
          const UnitVec3& exitDir_S,
          const Vec3&     exitDirDot_S) const;

    bool hasNearPoint() const {return !nearPointInS.isNaN(); }
    const Vec3& getNearPoint() const {return nearPointInS;}

    bool hasContactPointHints() const 
    {   return !(xPhint.isNaN() || xQhint.isNaN()); }
    void getContactPointHints(Vec3& xPhint, Vec3& xQhint) const
    {   xPhint = this->xPhint; xQhint = this->xQhint; }
 
private:
friend class CableObstacle::Surface;

    ContactGeometry      surface;
    Vec3                nearPointInS; // Cartesian location of N, in S frame
    Vec3                xPhint, xQhint;
};


#endif // SimTK_SIMBODY_CABLE_PATH_IMPL_H_

